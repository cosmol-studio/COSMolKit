use crate::source::base::ichi_bns::{
    AddRemoveIsoProtonsRestr, AddRemoveProtonsRestr, AllocateAndInitBnData, DeAllocateBnData,
    DeAllocateBnStruct, SetForbiddenEdges,
};
use crate::source::base::ichi_io::{inchi_strbuf_close, inchi_strbuf_init};
use crate::source::base::ichimake::CompareReversedINChI2;
use crate::source::base::ichinorm::{FreeInpAtomData, MarkRingSystemsInp};
use crate::source::base::ichiprt1::MergeZzInHillFormula;
use crate::source::base::ichiread::bRevInchiComponentExists;
use crate::source::base::ichiring::is_bond_in_Nmax_memb_ring;
use crate::source::base::ichirvr1::{
    AddCGroups2TCGBnStruct, AddTGroups2TCGBnStruct, AddToEdgeList, AllocBfsQueue, AllocEdgeList,
    AllocateAndInitTCGBnStruct, CN_LIST, ConnectDisconnectedH, DisconnectedConnectedH,
    FillOutpStructEndpointFromInChI, GetAtomRestoreInfo, GetChargeFlowerUpperEdge,
    GetTgroupInfoFromInChI, MakeInChIOutOfStrFromINChI2, MakeOneInChIOutOfStrFromINChI,
    MakeOneInChIOutOfStrFromINChI2, RemoveForbiddenBondFlowBits, RemoveForbiddenEdgeMask,
    RunBnsRestoreOnce, RunBnsTestOnce, SetForbiddenEdgeMask, get_sp_element_type, nAddSuperCGroups,
    nCountBnsSizes,
};
use crate::source::base::ichirvr2::{
    CheckBnsConsistency, Convert_SIV_to_SVI, CopyBnsToAtom, EliminateChargeSeparationOnHeteroatoms,
    EliminateNitrogen5Val3Bonds, FixMetal_Nminus_Ominus, IncrementZeroOrderBondsToHeteroat,
    MoveChargeToMakeCenerpoints, MoveMobileHToAvoidFixedBonds, MovePlusFromS2DiaminoCarbon,
    MoveRadToAtomsAddCharges, PlusFromDB_N_DB_O_to_Metal, RearrangePlusMinusEdgesFlow,
    RemoveRadFromMobileHEndpoint, RemoveRadFromMobileHEndpointFixH,
    RestoreAtomConnectionsSetStereo, RestoreCyanoGroup, RestoreIsoCyanoGroup, RestoreNNNgroup,
    SetStereoBondTypesFrom0DStereo,
};
use crate::source::base::ichirvr3::FixFixedHRestoredStructure;
use crate::source::base::ichirvr5::{FixMobileHRestoredStructure, GetPlusMinusVertex};
use crate::source::base::ichirvr6::FixRestoredStructureStereo;
use crate::source::base::ichister::ReconcileAllCmlBondParities;
use crate::source::base::ichitaut::free_t_group_info;
use crate::source::base::ichitaut::is_centerpoint_elem;
use crate::source::base::runichi4::FreeAllINChIArrays;
use crate::source::base::strutil::{Free_INChI, Free_INChI_Aux, bIsMetalSalt};
use crate::source::base::util::{
    get_endpoint_valence, inchi_calloc, inchi_free, inchi_malloc, inchi_realloc, is_el_a_metal,
};
use crate::source_types::{
    ALL_TC_GROUPS, AT_NUMB, BFS_Q, BFS_Q_CLEAR, BFS_Q_FREE, BN_DATA, BN_MAX_ALTP, BN_STRUCT,
    BNS_EDGE_FORBIDDEN_MASK, BNS_EDGE_FORBIDDEN_TEMP, BNS_EDGE_FORBIDDEN_TEST, BNS_OUT_OF_RAM,
    BNS_VERT_TYPE_C_GROUP, BNS_VT_C_POS_ALL, BOND_ALT12NS, BOND_TYPE_DOUBLE, BOND_TYPE_MASK,
    BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, CANON_GLOBALS, CMP2FHINCHI, CMP2MHINCHI, EDGE_LIST,
    EDGE_LIST_CLEAR, EDGE_LIST_FREE, ICR, INCHI_BAS, INCHI_CLOCK, INCHI_IOS_STRING, INCHI_MODE,
    INCHI_REC, INChI, INPUT_PARMS, MAX_DIFF_FIXH, MAX_NUM_CN_BITS, MAXVAL, NO_VERTEX, NUM_H,
    NUM_H_ISOTOPES, REQ_MODE_BASIC, RI_ERR_ALLOC, RI_ERR_PROGR, RI_ERR_SYNTAX, STRUCT_DATA,
    SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, SourceTGroupInfoPointer,
    StrFromINChI, TAUT_NON, TAUT_NUM, TAUT_YES, VAL_AT, clock_t, cn_bits_M, cn_bits_N, cn_bits_P,
    cn_bits_shift, inp_ATOM, tagTCGroupTypes_TCG_Minus_C0 as TCG_Minus_C0,
    tagTCGroupTypes_TCG_None as TCG_None, tagTCGroupTypes_TCG_Plus_C0 as TCG_Plus_C0,
};

const IDIF_PROBLEM: INCHI_MODE = 0x0000_0001;
const IDIF_MORE_H: INCHI_MODE = 0x0000_0010;
const IDIF_LESS_H: INCHI_MODE = 0x0000_0020;
const IDIF_EXTRA_TG_ENDP: INCHI_MODE = 0x0000_0800;

#[allow(non_snake_case)]
pub(crate) fn ForbidCarbonChargeEdges(
    heap: &mut SourceHeap,
    pBNS: &BN_STRUCT,
    pTCGroups: &ALL_TC_GROUPS,
    pCarbonChargeEdges: &mut EDGE_LIST,
    forbidden_edge_mask: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:57 ForbidCarbonChargeEdges
    /*
    int ForbidCarbonChargeEdges(BN_STRUCT* pBNS,
        ALL_TC_GROUPS* pTCGroups,
        EDGE_LIST* pCarbonChargeEdges,
        int forbidden_edge_mask)
    {
    #define MAX_NUM_CARBON_CHARGE_EDGES 2
        int nType, i, k, ret;
        BNS_EDGE* pEdge;
        if ((ret = AllocEdgeList(pCarbonChargeEdges, MAX_NUM_CARBON_CHARGE_EDGES))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }
        pCarbonChargeEdges->num_edges = 0;
        for (i = 0; i < MAX_NUM_CARBON_CHARGE_EDGES; i++)
        {
            switch (i)
            {
            case 0:
                nType = TCG_Plus_C0;
                break;
            case 1:
                nType = TCG_Minus_C0;
                break;
            default:
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            if ((k = pTCGroups->nGroup[nType]) >= 0)
            {
                k = pTCGroups->pTCG[k].nForwardEdge;
                if (k > 0)
                {
                    pEdge = pBNS->edge + k;
                    if (!(pEdge->forbidden & forbidden_edge_mask))
                    {
                        pEdge->forbidden |= forbidden_edge_mask;
                        if ((ret = AddToEdgeList(pCarbonChargeEdges, k, 0))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
                else
                {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
            }
        }
        ret = pCarbonChargeEdges->num_edges;

    exit_function:

        return ret;
    #undef MAX_NUM_CARBON_CHARGE_EDGES
    }
    */
    // END INCHI C FUNCTION: ForbidCarbonChargeEdges
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ForbidCarbonChargeEdges
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: #define MAX_NUM_CARBON_CHARGE_EDGES 2
    // INCHI✔️❌: SourceHeap checked allocation lookup adds overhead versus direct C pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: ForbidCarbonChargeEdges

    const MAX_NUM_CARBON_CHARGE_EDGES: i32 = 2;
    let mut ret = AllocEdgeList(heap, pCarbonChargeEdges, MAX_NUM_CARBON_CHARGE_EDGES)?;
    if ret != 0 {
        return Ok(ret);
    }
    pCarbonChargeEdges.num_edges = 0;

    for group_type in [TCG_Plus_C0, TCG_Minus_C0] {
        let mut index = pTCGroups.nGroup[group_type as usize];
        if index >= 0 {
            index = heap
                .slice(pTCGroups.pTCG.as_const())?
                .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .nForwardEdge;
            if index > 0 {
                let edge = heap
                    .slice_mut(pBNS.edge)?
                    .get_mut(
                        usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if edge.forbidden as i32 & forbidden_edge_mask == 0 {
                    edge.forbidden = (edge.forbidden as i32 | forbidden_edge_mask) as i8;
                    ret = AddToEdgeList(heap, pCarbonChargeEdges, index, 0)?;
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            } else {
                return Ok(RI_ERR_PROGR);
            }
        }
    }
    Ok(pCarbonChargeEdges.num_edges)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ForbidNintrogenPlus2BondsInSmallRings(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    at: &[inp_ATOM],
    num_at: i32,
    pVA: &[VAL_AT],
    min_ring_size: i32,
    _pTCGroups: &ALL_TC_GROUPS,
    pNplus2BondsEdges: &mut EDGE_LIST,
    forbidden_edge_mask: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:116 ForbidNintrogenPlus2BondsInSmallRings
    /*
    int ForbidNintrogenPlus2BondsInSmallRings(BN_STRUCT* pBNS,
        inp_ATOM* at,
        int num_at,
        VAL_AT* pVA,
        int min_ring_size,
        ALL_TC_GROUPS* pTCGroups,
        EDGE_LIST* pNplus2BondsEdges,
        int forbidden_edge_mask)
    {
        int i, j, ret;
        BNS_EDGE* e;

        /* djb-rwth: removing redundant code */
            /* --- forbid edges that allow to make =N(+)= or #N(+)- in small ring */
        for (i = 0; i < num_at; i++)
        {
            if (at[i].valence == 2 &&
                !at[i].num_H && !at[i].endpoint &&
                pVA[i].cNumValenceElectrons == 5 &&
                pVA[i].cPeriodicRowNumber == 1 &&
                !pVA[i].cMaxFlowToMetal && pVA[i].nCPlusGroupEdge > 0 &&
                pVA[i].cnListIndex > 0 && cnList[pVA[i].cnListIndex - 1].bits == cn_bits_MNP &&
                pVA[i].cMinRingSize && pVA[i].cMinRingSize <= min_ring_size)
            {

                e = pBNS->edge + (j = pVA[i].nCPlusGroupEdge - 1);
                if (!(e->forbidden & forbidden_edge_mask))
                {
                    e->forbidden |= forbidden_edge_mask;
                    if ((ret = AddToEdgeList(pNplus2BondsEdges, j, 128))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
        }
        ret = 0;

    exit_function:

        return ret;
    }
    */
    // END INCHI C FUNCTION: ForbidNintrogenPlus2BondsInSmallRings
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ForbidNintrogenPlus2BondsInSmallRings
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: cn_bits_MNP is the configured MAKE_CN_BITS value 140 (CN_LIST entry 4).
    // INCHI✔️❌: SourceHeap checked edge/list access adds overhead versus direct C pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: ForbidNintrogenPlus2BondsInSmallRings

    let mut i = 0_i32;
    while i < num_at {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = at.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let cn_index = i32::from(valence.cnListIndex).wrapping_sub(1);
        if atom.valence == 2
            && atom.num_H == 0
            && atom.endpoint == 0
            && valence.cNumValenceElectrons == 5
            && valence.cPeriodicRowNumber == 1
            && valence.cMaxFlowToMetal == 0
            && valence.nCPlusGroupEdge > 0
            && valence.cnListIndex > 0
            && CN_LIST
                .get(usize::try_from(cn_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bits
                == CN_LIST[4].bits
            && valence.cMinRingSize != 0
            && i32::from(valence.cMinRingSize) <= min_ring_size
        {
            let edge_index = valence.nCPlusGroupEdge.wrapping_sub(1);
            let edge = heap
                .slice_mut(pBNS.edge)?
                .get_mut(
                    usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(edge.forbidden) & forbidden_edge_mask == 0 {
                edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                let ret = AddToEdgeList(heap, pNplus2BondsEdges, edge_index, 128)?;
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }
        i = i.wrapping_add(1);
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixLessHydrogenInFormula(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &StrFromINChI,
    at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    atf: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:177 FixLessHydrogenInFormula
    // INCHI✔️❌: complete active source frame follows verbatim; SourceHeap pointer checks add known overhead.
    /*
    int FixLessHydrogenInFormula(BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        inp_ATOM* at,
        inp_ATOM* at2,
        inp_ATOM* atf,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups,
        int* pnNumRunBNS,
        int* pnTotalDelta, int forbidden_edge_mask)
    {
        int iBPlus = NO_VERTEX, iNV = NO_VERTEX, iNH = NO_VERTEX, neigh;
        EDGE_LIST NewlyFixedEdges;
        int ret, i, j;
        int num_at = pStruct->num_atoms;
        int inv_forbidden_edge_mask = ~forbidden_edge_mask;
        /* for RunBnsTestOnce */
        Vertex     vPathStart, vPathEnd;
        int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

        AllocEdgeList(&NewlyFixedEdges, EDGE_LIST_CLEAR);
        if ((ret = AllocEdgeList(&NewlyFixedEdges, 2 * num_at))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }
        for (i = 0; i < num_at; i++)
        {
            if ((j = pVA[i].nCMinusGroupEdge - 1) >= 0)
            {
                if ((ret = AddToEdgeList(&NewlyFixedEdges, j, 0))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                pBNS->edge[j].forbidden |= forbidden_edge_mask;
            }
            if ((j = pVA[i].nCPlusGroupEdge - 1) >= 0)
            {
                if ((ret = AddToEdgeList(&NewlyFixedEdges, j, 0))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                pBNS->edge[j].forbidden |= forbidden_edge_mask;
            }
        }
        /* extra H has been removed; check non-tautomeric atoms */
        for (i = 0; i < num_at; i++)
        {
            if (!at2[i].endpoint && !pVA[i].cMetal &&
                pVA[i].cNumValenceElectrons == 5 && pVA[i].cPeriodicRowNumber == 1 &&
                at2[i].num_H == atf[i].num_H + 1)
            {
                /* H was removed from N */
                iNH = i;
                break;
            }
        }
        if (0 <= iNH && iNH < num_at)
        {
            /* check neighbors for  |                 |
                              (a)  -B(+)-  or  (b)   =N-
                                    |                 |
            */
            for (j = 0; j < at2[i].valence; j++)
            {
                neigh = at2[iNH].neighbor[j];
                if (at2[neigh].valence == 4)
                {
                    if (at2[neigh].charge == -1 && at2[neigh].chem_bonds_valence == 4 &&
                        !at2[neigh].radical && !at[neigh].num_H)
                    {
                        iBPlus = neigh;
                    }
                }
            }
        }
        if (0 <= iNH && iNH < num_at)
        {
            int bond_type_at2;
            int bond_type_atf;
            /* djb-rwth: removing redundant variables */
            int delta = -1, nxt = iNH, prv = NO_VERTEX, nxt_is_NPlus;
            /* the changed bond to the dehydrogenated atom H should have greater order */
            /* delta = (new bond order in atf[]) - (restored bond order in at2[]) */
            nxt_is_NPlus = 0;
            do
            {
                i = nxt;
                nxt = NO_VERTEX;
                delta = -delta;
                for (j = 0; j < at2[i].valence; j++)
                {
                    bond_type_at2 = at2[i].bond_type[j] & BOND_TYPE_MASK; /* restored bond */
                    bond_type_atf = atf[i].bond_type[j] & BOND_TYPE_MASK; /* normalized bond */
                    nxt_is_NPlus = 0;
                    if ((bond_type_atf - bond_type_at2 == delta || bond_type_atf == BOND_ALT12NS) &&
                        BOND_TYPE_SINGLE <= bond_type_at2 + delta && bond_type_at2 + delta <= BOND_TYPE_TRIPLE &&
                        !at2[(int)at2[i].neighbor[j]].cFlags)
                    {
                        prv = i;
                        nxt = at2[i].neighbor[j];
                        nxt_is_NPlus = at2[nxt].charge == 1 && atf[nxt].charge == 0 &&
                            pVA[nxt].cNumValenceElectrons == 5 && pVA[nxt].cPeriodicRowNumber == 1;
                        at2[i].cFlags |= 1;  /* avoid cycling */
                        /* djb-rwth: removing redundant code */
                        if (delta == -1 && at2[prv].valence == 4 && at2[prv].chem_bonds_valence == 5 &&
                            !at2[prv].charge && !at2[prv].radical && pVA[prv].cNumValenceElectrons == 5 &&
                            pVA[prv].nCPlusGroupEdge > 0)
                        {
                            iNV = prv;
                        }
                        if (at2[nxt].charge != atf[nxt].charge)
                        {
                            if ((at2[nxt].charge == 1 || atf[nxt].charge == 1) &&
                                pVA[nxt].nCPlusGroupEdge > 0)
                            {
                                pBNS->edge[pVA[nxt].nCPlusGroupEdge - 1].forbidden &= inv_forbidden_edge_mask;
                            }
                            if ((at2[nxt].charge == -1 || atf[nxt].charge == -1) &&
                                pVA[nxt].nCMinusGroupEdge > 0)
                            {
                                pBNS->edge[pVA[nxt].nCMinusGroupEdge - 1].forbidden &= inv_forbidden_edge_mask;
                            }
                        }
                        break; /* found */
                    }
                }
            } while (nxt >= 0 && !(nxt_is_NPlus && delta == -1));
            for (i = 0; i < num_at; i++)
            {
                at2[i].cFlags = 0;
            }
            if (nxt >= 0 && nxt_is_NPlus && delta == -1)
            {
                /* a simple alt path from NH-= to =N(+) has been found */
                if (iBPlus || iNV)
                {
                    /* move (+) charge from N(+) to iNV or, if iBPlus, then to iNH */
                    if ((iNV >= 0 && (j = pVA[iNV].nCPlusGroupEdge - 1) > 0 && pBNS->edge[j].flow > 0) ||
                        (iNH >= 0 && (j = pVA[iNH].nCPlusGroupEdge - 1) > 0 && pBNS->edge[j].flow > 0)) /* djb-rwth: addressing LLVM warnings */
                    {
                        int          ieFlower;
                        BNS_EDGE* pe = pBNS->edge + j, * peFlower = NULL;
                        Vertex      v1 = pe->neighbor1;
                        Vertex      v2 = v1 ^ pe->neighbor12;
                        BNS_VERTEX* pv1 = pBNS->vert + v1;
                        BNS_VERTEX* pv2 = pBNS->vert + v2;

                        delta = 1;
                        /* prevent conversion of >N(+)= into N(V) neutral */
                        ieFlower = GetChargeFlowerUpperEdge(pBNS, pVA, pVA[nxt].nCPlusGroupEdge - 1);
                        if (ieFlower >= 0)
                        {
                            peFlower = pBNS->edge + ieFlower;
                            if (peFlower->flow == delta)
                            {
                                peFlower->forbidden |= forbidden_edge_mask;
                                if ((ret = AddToEdgeList(&NewlyFixedEdges, ieFlower, 0))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                            }
                        }
                        pe->forbidden |= forbidden_edge_mask;
                        pe->flow -= delta;
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;
                        ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                            &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);
                        if (ret < 0)
                        {
                            goto exit_function;
                        }
                        if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1)) &&
                            nDeltaCharge <= 0  /* charge moving to this atom disappers*/) /* djb-rwth: addressing LLVM warnings */
                        {
                            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                            (*pnNumRunBNS)++;
                            if (ret < 0)
                            {
                                goto exit_function;
                            }
                            else
                                if (ret == 1)
                                {
                                    *pnTotalDelta += ret;
                                }
                                else
                                {
                                    ret = RI_ERR_PROGR;
                                    goto exit_function;
                                }
                        }
                        else
                        {
                            ret = 0;
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                    }
                }
            }
        }

    exit_function:
        /* remove bond fixation */
        RemoveForbiddenEdgeMask(pBNS, &NewlyFixedEdges, forbidden_edge_mask);
        AllocEdgeList(&NewlyFixedEdges, EDGE_LIST_FREE);

        return ret;
    }
    */
    // END INCHI C FUNCTION: FixLessHydrogenInFormula
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FixLessHydrogenInFormula
    // INCHI✔️❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: BOND_TYPE_MASK=15, BOND_ALT12NS=9, NO_VERTEX=-2; signed-char narrowing is explicit.
    // END INCHI ACTIVE MACRO CONFIGURATION: FixLessHydrogenInFormula

    let num_at = pStruct.num_atoms;
    let inverse_mask = !forbidden_edge_mask;
    let mut newly_fixed_edges = EDGE_LIST::default();
    let _ = AllocEdgeList(heap, &mut newly_fixed_edges, EDGE_LIST_CLEAR)?;

    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut ret = AllocEdgeList(heap, &mut newly_fixed_edges, num_at.wrapping_mul(2))?;
        if ret != 0 {
            return Ok(ret);
        }

        let mut atom_number = 0_i32;
        while atom_number < num_at {
            let index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            for edge_number in [
                valence.nCMinusGroupEdge.wrapping_sub(1),
                valence.nCPlusGroupEdge.wrapping_sub(1),
            ] {
                if edge_number >= 0 {
                    ret = AddToEdgeList(heap, &mut newly_fixed_edges, edge_number, 0)?;
                    if ret != 0 {
                        return Ok(ret);
                    }
                    let edge = heap
                        .slice_mut(pBNS.edge)?
                        .get_mut(
                            usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                }
            }
            atom_number = atom_number.wrapping_add(1);
        }

        let original_atoms = if num_at <= 0 {
            Vec::new()
        } else {
            heap.slice(at.as_const())?
                .get(..usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        };
        let normalized_atoms = if num_at <= 0 {
            Vec::new()
        } else {
            heap.slice(atf.as_const())?
                .get(..usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        };

        let mut hydrogen_atom = NO_VERTEX;
        atom_number = 0;
        while atom_number < num_at {
            let index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let restored = heap
                .slice(at2.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let normalized = normalized_atoms
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if restored.endpoint == 0
                && valence.cMetal == 0
                && valence.cNumValenceElectrons == 5
                && valence.cPeriodicRowNumber == 1
                && i32::from(restored.num_H) == i32::from(normalized.num_H).wrapping_add(1)
            {
                hydrogen_atom = atom_number;
                break;
            }
            atom_number = atom_number.wrapping_add(1);
        }

        let mut boron_plus = NO_VERTEX;
        if hydrogen_atom >= 0 && hydrogen_atom < num_at {
            let index =
                usize::try_from(hydrogen_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let restored = heap
                .slice(at2.as_const())?
                .get(index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut neighbor_order = 0_i32;
            while neighbor_order < i32::from(restored.valence) {
                let neighbor = usize::from(
                    *restored
                        .neighbor
                        .get(
                            usize::try_from(neighbor_order)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let neighbor_atom = heap
                    .slice(at2.as_const())?
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let original = original_atoms
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if neighbor_atom.valence == 4
                    && neighbor_atom.charge == -1
                    && neighbor_atom.chem_bonds_valence == 4
                    && neighbor_atom.radical == 0
                    && original.num_H == 0
                {
                    boron_plus = neighbor as i32;
                }
                neighbor_order = neighbor_order.wrapping_add(1);
            }
        }

        let mut nitrogen_v = NO_VERTEX;
        if hydrogen_atom >= 0 && hydrogen_atom < num_at {
            let mut delta = -1_i32;
            let mut next = hydrogen_atom;
            let mut next_is_nitrogen_plus = false;
            loop {
                let current = next;
                next = NO_VERTEX;
                delta = delta.wrapping_neg();
                let current_index =
                    usize::try_from(current).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let restored = heap
                    .slice(at2.as_const())?
                    .get(current_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let normalized = normalized_atoms
                    .get(current_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let mut neighbor_order = 0_i32;
                while neighbor_order < i32::from(restored.valence) {
                    let order = usize::try_from(neighbor_order)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let restored_bond =
                        i32::from(restored.bond_type[order]) & BOND_TYPE_MASK as i32;
                    let normalized_bond =
                        i32::from(normalized.bond_type[order]) & BOND_TYPE_MASK as i32;
                    next_is_nitrogen_plus = false;
                    let neighbor = usize::from(restored.neighbor[order]);
                    let neighbor_flags = heap
                        .slice(at2.as_const())?
                        .get(neighbor)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .cFlags;
                    let adjusted_bond = restored_bond.wrapping_add(delta);
                    if (normalized_bond.wrapping_sub(restored_bond) == delta
                        || normalized_bond == BOND_ALT12NS as i32)
                        && BOND_TYPE_SINGLE as i32 <= adjusted_bond
                        && adjusted_bond <= BOND_TYPE_TRIPLE as i32
                        && neighbor_flags == 0
                    {
                        next = neighbor as i32;
                        let next_restored = heap
                            .slice(at2.as_const())?
                            .get(neighbor)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let next_normalized = normalized_atoms
                            .get(neighbor)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let next_valence = pVA
                            .get(neighbor)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        next_is_nitrogen_plus = next_restored.charge == 1
                            && next_normalized.charge == 0
                            && next_valence.cNumValenceElectrons == 5
                            && next_valence.cPeriodicRowNumber == 1;
                        heap.slice_mut(at2)?[current_index].cFlags |= 1;
                        let current_valence = pVA
                            .get(current_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if delta == -1
                            && restored.valence == 4
                            && restored.chem_bonds_valence == 5
                            && restored.charge == 0
                            && restored.radical == 0
                            && current_valence.cNumValenceElectrons == 5
                            && current_valence.nCPlusGroupEdge > 0
                        {
                            nitrogen_v = current;
                        }
                        if next_restored.charge != next_normalized.charge {
                            for edge_number in [
                                if (next_restored.charge == 1 || next_normalized.charge == 1)
                                    && next_valence.nCPlusGroupEdge > 0
                                {
                                    Some(next_valence.nCPlusGroupEdge.wrapping_sub(1))
                                } else {
                                    None
                                },
                                if (next_restored.charge == -1 || next_normalized.charge == -1)
                                    && next_valence.nCMinusGroupEdge > 0
                                {
                                    Some(next_valence.nCMinusGroupEdge.wrapping_sub(1))
                                } else {
                                    None
                                },
                            ]
                            .into_iter()
                            .flatten()
                            {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(
                                        usize::try_from(edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden = (i32::from(edge.forbidden) & inverse_mask) as i8;
                            }
                        }
                        break;
                    }
                    neighbor_order = neighbor_order.wrapping_add(1);
                }
                if next < 0 || (next_is_nitrogen_plus && delta == -1) {
                    break;
                }
            }

            atom_number = 0;
            while atom_number < num_at {
                let index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                heap.slice_mut(at2)?[index].cFlags = 0;
                atom_number = atom_number.wrapping_add(1);
            }

            if next >= 0
                && next_is_nitrogen_plus
                && delta == -1
                && (boron_plus != 0 || nitrogen_v != 0)
            {
                let mut selected_edge = None;
                for atom in [nitrogen_v, hydrogen_atom] {
                    if atom >= 0 {
                        let edge_number = pVA
                            .get(
                                usize::try_from(atom)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCPlusGroupEdge
                            .wrapping_sub(1);
                        if edge_number > 0
                            && heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(edge_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .flow
                                > 0
                        {
                            selected_edge = Some(edge_number);
                            break;
                        }
                    }
                }

                if let Some(edge_number) = selected_edge {
                    let next_plus_edge = pVA
                        .get(
                            usize::try_from(next)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let flower = GetChargeFlowerUpperEdge(heap, pBNS, pVA, next_plus_edge)?;
                    if flower >= 0 {
                        let flower_index = usize::try_from(flower)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if heap
                            .slice(pBNS.edge.as_const())?
                            .get(flower_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            == 1
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(flower_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                            ret = AddToEdgeList(heap, &mut newly_fixed_edges, flower, 0)?;
                            if ret != 0 {
                                return Ok(ret);
                            }
                        }
                    }

                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let selected = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let first_vertex = i32::from(selected.neighbor1);
                    let second_vertex = i32::from(selected.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(1);
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
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                        }
                    }
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
                    if ret < 0 {
                        return Ok(ret);
                    }
                    if ret == 1
                        && ((path_end == first_vertex && path_start == second_vertex)
                            || (path_end == second_vertex && path_start == first_vertex))
                        && delta_charge <= 0
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
                        if ret < 0 {
                            return Ok(ret);
                        }
                        if ret == 1 {
                            *pnTotalDelta = pnTotalDelta.wrapping_add(ret);
                        } else {
                            return Ok(RI_ERR_PROGR);
                        }
                    } else {
                        ret = 0;
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = heap
                            .slice(pBNS.edge.as_const())?[edge_index]
                            .flow
                            .wrapping_add(1);
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
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                }
            }
        }
        Ok(ret)
    })();

    let remove_result =
        RemoveForbiddenEdgeMask(heap, pBNS, &newly_fixed_edges, forbidden_edge_mask);
    let free_result = AllocEdgeList(heap, &mut newly_fixed_edges, EDGE_LIST_FREE);
    match (execution, remove_result, free_result) {
        (Err(error), _, _) => Err(error),
        (Ok(_), Err(error), _) => Err(error),
        (Ok(_), Ok(_), Err(error)) => Err(error),
        (Ok(ret), Ok(()), Ok(_)) => Ok(ret),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixMoreHydrogenInFormula(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &StrFromINChI,
    _at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    atf: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:398 FixMoreHydrogenInFormula
    // INCHI✔️❌: complete active source frame follows verbatim; SourceHeap pointer checks add known overhead.
    /*
    int FixMoreHydrogenInFormula(BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        inp_ATOM* at,
        inp_ATOM* at2,
        inp_ATOM* atf,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups,
        int* pnNumRunBNS,
        int* pnTotalDelta,
        int forbidden_edge_mask)
    {
        int iNH = NO_VERTEX, neigh, neigh2;
        EDGE_LIST NewlyFixedEdges;
        int ret, i, j, k, k2 = 0, delta;
        int num_at = pStruct->num_atoms;
        int inv_forbidden_edge_mask = ~forbidden_edge_mask;
        Vertex v1, v2;
        /* for RunBnsTestOnce */
        Vertex     vPathStart, vPathEnd;
        int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
        BNS_EDGE* pe, * pe2;

        AllocEdgeList(&NewlyFixedEdges, EDGE_LIST_CLEAR);
        if ((ret = AllocEdgeList(&NewlyFixedEdges, 2 * num_at))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }
        /* fix all charges */
        for (i = 0; i < num_at; i++)
        {
            if ((j = pVA[i].nCMinusGroupEdge - 1) >= 0)
            {
                if ((ret = AddToEdgeList(&NewlyFixedEdges, j, 0))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                pBNS->edge[j].forbidden |= forbidden_edge_mask;
            }
            if ((j = pVA[i].nCPlusGroupEdge - 1) >= 0)
            {
                if ((ret = AddToEdgeList(&NewlyFixedEdges, j, 0))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                pBNS->edge[j].forbidden |= forbidden_edge_mask;
            }
        }

        /* H(+) has been added to -O(-); check non-tautomeric atoms */
        for (i = 0; i < num_at; i++)
        {
            neigh = at2[i].neighbor[0]; /* djb-rwth: avoiding unsequenced modification and access to neigh */
            if (!(pStruct->bMobileH ? at2[i].endpoint : pStruct->endpoint[i]) && !pVA[i].cMetal &&
                at2[i].num_H + 1 == atf[i].num_H &&      /* normalization added H ??? What would happen in Fixed-H case?*/
                (k = pVA[i].nCMinusGroupEdge - 1) >= 0 &&
                pBNS->edge[k].flow == 1 &&               /* atom had (-) charge before preprocessing */
                at2[i].charge == -1 && atf[i].charge == 0 && /* and has no charge after preprocessing */
                at2[i].valence == 1 && at2[i].chem_bonds_valence == 1 && /* connected by a single bond */
                pVA[i].cNumValenceElectrons == 6 &&     /* atom is O, S, Se, Te */
                at2[neigh].chem_bonds_valence > at2[neigh].valence
                /* atom's single neighbor has multiple bond(s)*/
                )
            {
                /* H(+) was added to O in Y=X-O(-), where X is the only neighbor of O, X=neigh, Y=neigh2 */
                iNH = i;
                for (j = 0; j < at2[neigh].valence; j++)
                {
                    neigh2 = at2[neigh].neighbor[j];
                    if (neigh2 != iNH && !at2[neigh2].endpoint &&
                        !pBNS->edge[(int)pBNS->vert[neigh].iedge[j]].forbidden &&
                        4 <= pVA[neigh2].cNumValenceElectrons &&
                        pVA[neigh2].cNumValenceElectrons <= 5 && /* neig2 is C or N */
                        (k2 = pVA[neigh2].nCMinusGroupEdge - 1) >= 0 &&
                        !pBNS->edge[k2].flow /* negative charge may be moved to neigh2 */)
                    {
                        break;
                    }
                }
                if (j < at2[neigh].valence)
                {
                    delta = 1;
                    pe = pBNS->edge + k;  /* -O(-) negative charge edge; flow = 1 */
                    pe2 = pBNS->edge + k2; /* X charge edge; flow = 0 */
                    v1 = pe->neighbor1;
                    v2 = pe->neighbor12 ^ v1;
                    pe->flow -= delta;
                    pBNS->vert[v1].st_edge.flow -= delta;
                    pBNS->vert[v2].st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;
                    pe2->forbidden &= inv_forbidden_edge_mask; /* allow the charge to move */

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);
                    if (ret < 0)
                    {
                        goto exit_function;
                    }
                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) &&
                        nDeltaCharge <= 1) /* djb-rwth: addressing LLVM warnings */
                    {
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        (*pnNumRunBNS)++;
                        if (ret < 0)
                        {
                            goto exit_function;
                        }
                        else
                            if (ret)
                            {
                                *pnTotalDelta += ret;
                            }
                            else
                            {
                                ret = RI_ERR_PROGR;
                            }
                        break;
                    }
                    else
                    {
                        /* the attempt has failed; restore the flow */
                        ret = 0;
                        pe->flow += delta;
                        pBNS->vert[v1].st_edge.flow += delta;
                        pBNS->vert[v2].st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                }
            }
        }

    exit_function:
        /* remove bond fixation */
        RemoveForbiddenEdgeMask(pBNS, &NewlyFixedEdges, forbidden_edge_mask);
        AllocEdgeList(&NewlyFixedEdges, EDGE_LIST_FREE);

        return ret;
    }
    */
    // END INCHI C FUNCTION: FixMoreHydrogenInFormula
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FixMoreHydrogenInFormula
    // INCHI✔️❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: BNS edge flow and signed-char forbidden updates preserve source ordering and narrowing.
    // END INCHI ACTIVE MACRO CONFIGURATION: FixMoreHydrogenInFormula

    let num_at = pStruct.num_atoms;
    let inverse_mask = !forbidden_edge_mask;
    let mut newly_fixed_edges = EDGE_LIST::default();
    let _ = AllocEdgeList(heap, &mut newly_fixed_edges, EDGE_LIST_CLEAR)?;
    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut ret = AllocEdgeList(heap, &mut newly_fixed_edges, num_at.wrapping_mul(2))?;
        if ret != 0 {
            return Ok(ret);
        }

        let mut atom_number = 0_i32;
        while atom_number < num_at {
            let index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            for edge_number in [
                valence.nCMinusGroupEdge.wrapping_sub(1),
                valence.nCPlusGroupEdge.wrapping_sub(1),
            ] {
                if edge_number >= 0 {
                    ret = AddToEdgeList(heap, &mut newly_fixed_edges, edge_number, 0)?;
                    if ret != 0 {
                        return Ok(ret);
                    }
                    let edge = heap
                        .slice_mut(pBNS.edge)?
                        .get_mut(
                            usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                }
            }
            atom_number = atom_number.wrapping_add(1);
        }

        let normalized_atoms = if num_at <= 0 {
            Vec::new()
        } else {
            heap.slice(atf.as_const())?
                .get(..usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        };
        let fixed_endpoints = if pStruct.bMobileH == 0 && num_at > 0 {
            Some(
                heap.slice(pStruct.endpoint.as_const())?
                    .get(
                        ..usize::try_from(num_at)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec(),
            )
        } else {
            None
        };

        atom_number = 0;
        while atom_number < num_at {
            let index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let restored = heap
                .slice(at2.as_const())?
                .get(index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let normalized = normalized_atoms
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor = usize::from(restored.neighbor[0]);
            let endpoint = if pStruct.bMobileH != 0 {
                restored.endpoint
            } else {
                *fixed_endpoints
                    .as_ref()
                    .and_then(|values| values.get(index))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            };
            let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
            let eligible = endpoint == 0
                && valence.cMetal == 0
                && i32::from(restored.num_H).wrapping_add(1) == i32::from(normalized.num_H)
                && minus_edge >= 0
                && heap
                    .slice(pBNS.edge.as_const())?
                    .get(
                        usize::try_from(minus_edge)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .flow
                    == 1
                && restored.charge == -1
                && normalized.charge == 0
                && restored.valence == 1
                && restored.chem_bonds_valence == 1
                && valence.cNumValenceElectrons == 6
                && {
                    let center = heap
                        .slice(at2.as_const())?
                        .get(neighbor)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    center.chem_bonds_valence > center.valence
                };
            if eligible {
                let center = heap
                    .slice(at2.as_const())?
                    .get(neighbor)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let mut candidate_minus_edge = 0_i32;
                let mut neighbor_order = 0_i32;
                while neighbor_order < i32::from(center.valence) {
                    let order = usize::try_from(neighbor_order)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let second_neighbor = usize::from(center.neighbor[order]);
                    if second_neighbor != index {
                        let second_atom = heap
                            .slice(at2.as_const())?
                            .get(second_neighbor)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if second_atom.endpoint == 0 {
                            let center_vertex = heap
                                .slice(pBNS.vert.as_const())?
                                .get(neighbor)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let bond_edge = *heap
                                .slice(center_vertex.iedge.as_const())?
                                .get(order)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(bond_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .forbidden
                                == 0
                            {
                                let second_valence = pVA
                                    .get(second_neighbor)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if 4 <= second_valence.cNumValenceElectrons
                                    && second_valence.cNumValenceElectrons <= 5
                                {
                                    candidate_minus_edge =
                                        second_valence.nCMinusGroupEdge.wrapping_sub(1);
                                    if candidate_minus_edge >= 0
                                        && heap
                                            .slice(pBNS.edge.as_const())?
                                            .get(
                                                usize::try_from(candidate_minus_edge).map_err(
                                                    |_| SourceHeapError::PointerOutOfBounds,
                                                )?,
                                            )
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                                            .flow
                                            == 0
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    neighbor_order = neighbor_order.wrapping_add(1);
                }
                if neighbor_order < i32::from(center.valence) {
                    let source_edge_index = usize::try_from(minus_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let source_edge = heap
                        .slice(pBNS.edge.as_const())?
                        .get(source_edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let first_vertex = i32::from(source_edge.neighbor1);
                    let second_vertex = i32::from(source_edge.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[source_edge_index].flow =
                        source_edge.flow.wrapping_sub(1);
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
                    let target_edge_index = usize::try_from(candidate_minus_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let target = heap
                        .slice_mut(pBNS.edge)?
                        .get_mut(target_edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    target.forbidden = (i32::from(target.forbidden) & inverse_mask) as i8;

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
                    if ret < 0 {
                        return Ok(ret);
                    }
                    if ret == 1
                        && ((path_end == first_vertex && path_start == second_vertex)
                            || (path_end == second_vertex && path_start == first_vertex))
                        && delta_charge <= 1
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
                        if ret < 0 {
                            return Ok(ret);
                        }
                        if ret != 0 {
                            *pnTotalDelta = pnTotalDelta.wrapping_add(ret);
                        } else {
                            ret = RI_ERR_PROGR;
                        }
                        break;
                    } else {
                        ret = 0;
                        heap.slice_mut(pBNS.edge)?[source_edge_index].flow = source_edge.flow;
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
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                }
            }
            atom_number = atom_number.wrapping_add(1);
        }
        Ok(ret)
    })();

    let remove_result =
        RemoveForbiddenEdgeMask(heap, pBNS, &newly_fixed_edges, forbidden_edge_mask);
    let free_result = AllocEdgeList(heap, &mut newly_fixed_edges, EDGE_LIST_FREE);
    match (execution, remove_result, free_result) {
        (Err(error), _, _) => Err(error),
        (Ok(_), Err(error), _) => Err(error),
        (Ok(_), Ok(_), Err(error)) => Err(error),
        (Ok(ret), Ok(()), Ok(_)) => Ok(ret),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixRemoveExtraTautEndpoints(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &StrFromINChI,
    _at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    _atf: SourceMutPointer<inp_ATOM>,
    atn: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    picr: &ICR,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:589 FixRemoveExtraTautEndpoints
    // INCHI✔️❌: complete active source frame follows verbatim; SourceHeap pointer checks add known overhead.
    /*
    int FixRemoveExtraTautEndpoints(BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        inp_ATOM* at,
        inp_ATOM* at2,
        inp_ATOM* atf,
        inp_ATOM* atn,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups, ICR* picr,
        int* pnNumRunBNS,
        int* pnTotalDelta,
        int forbidden_edge_mask)
    {
        EDGE_LIST NewlyFixedEdges;
        int ret, i, j, k, delta, centerpoint, endpoint1, endpoint2;
        int num_at = pStruct->num_atoms;
        int inv_forbidden_edge_mask = ~forbidden_edge_mask;
        Vertex v1, v2;
        /* for RunBnsTestOnce */
        Vertex     vPathStart, vPathEnd;
        int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
        BNS_EDGE* pe, * pe2;

        ret = 0; /* djb-rwth: ignoring LLVM warning: value might be returned */

        AllocEdgeList(&NewlyFixedEdges, EDGE_LIST_CLEAR);
        if ((ret = AllocEdgeList(&NewlyFixedEdges, 2 * num_at))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }
        /* fix all charges */
        for (i = 0; i < num_at; i++)
        {
            if ((j = pVA[i].nCMinusGroupEdge - 1) >= 0)
            {
                if ((ret = AddToEdgeList(&NewlyFixedEdges, j, 0))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                pBNS->edge[j].forbidden |= forbidden_edge_mask;
            }
            if ((j = pVA[i].nCPlusGroupEdge - 1) >= 0)
            {
                if ((ret = AddToEdgeList(&NewlyFixedEdges, j, 0))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                pBNS->edge[j].forbidden |= forbidden_edge_mask;
            }
        }

        for (i = 0; i < picr->num_endp_in1_only; i++)
        {
            endpoint1 = picr->endp_in1_only[i] - 1;
            if (at2[endpoint1].valence == at2[endpoint1].chem_bonds_valence ||
                pVA[endpoint1].nCMinusGroupEdge <= 0)
            {
                continue;
            }
            /* find centerpoint */
            for (j = 0; j < at2[endpoint1].valence; j++)
            {
                if (BOND_TYPE_DOUBLE == (BOND_TYPE_MASK & at2[endpoint1].bond_type[j]))
                {
                    centerpoint = at2[endpoint1].neighbor[j];
                    if (at2[centerpoint].charge || pVA[centerpoint].nCPlusGroupEdge <= 0 ||
                        !is_centerpoint_elem(at2[centerpoint].el_number))
                    {
                        continue;
                    }
                    /* -- the centerpoint as depicted has no ChargeStruct flower ---
                    m = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[centerpoint].nCPlusGroupEdge-1 );
                    if ( m < 0 || pBNS->edge[m].flow ) {
                        continue;
                    }
                    */
                    /* find 2nd endpoint */
                    for (k = 0; k < at2[centerpoint].valence; k++)
                    {
                        if (BOND_TYPE_SINGLE != (BOND_TYPE_MASK & at2[centerpoint].bond_type[k]))
                        {
                            continue;
                        }
                        endpoint2 = at2[centerpoint].neighbor[k];
                        if (!at2[endpoint2].endpoint && atn[endpoint2].endpoint)
                        {
                            break;
                        }
                    }
                    if (k == at2[centerpoint].valence)
                    {
                        continue;
                    }
                    /* the centerpoint and two extra endpoints have been found */
                    pe = pBNS->edge + pVA[centerpoint].nCPlusGroupEdge - 1;
                    if (!pe->flow)
                    {
                        continue;
                    }
                    pe2 = pBNS->edge + pVA[endpoint1].nCMinusGroupEdge - 1;
                    if (pe2->flow)
                    {
                        continue;
                    }
                    delta = 1;
                    v1 = pe->neighbor1;
                    v2 = pe->neighbor12 ^ v1;
                    pe->flow -= delta;
                    pBNS->vert[v1].st_edge.flow -= delta;
                    pBNS->vert[v2].st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;
                    pe2->forbidden &= inv_forbidden_edge_mask; /* allow the charge to move */

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);
                    if (ret < 0)
                    {
                        goto exit_function;
                    }
                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) &&
                        nDeltaCharge <= 1) /* djb-rwth: addressing LLVM warnings */
                    {
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        (*pnNumRunBNS)++;
                        if (ret < 0)
                        {
                            goto exit_function;
                        }
                        else
                        {
                            if (ret)
                            {
                                *pnTotalDelta += ret;
                            }
                            else
                            {
                                ret = RI_ERR_PROGR;
                            }
                        }
                        goto exit_function;
                    }
                    else
                    {
                        ret = 0;
                        pe->flow += delta;
                        pBNS->vert[v1].st_edge.flow += delta;
                        pBNS->vert[v2].st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                        pe2->forbidden |= forbidden_edge_mask;
                    }
                }
            }
        }

    exit_function:
        /* remove bond fixation */
        RemoveForbiddenEdgeMask(pBNS, &NewlyFixedEdges, forbidden_edge_mask);
        AllocEdgeList(&NewlyFixedEdges, EDGE_LIST_FREE);

        return ret;
    }
    */
    // END INCHI C FUNCTION: FixRemoveExtraTautEndpoints
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FixRemoveExtraTautEndpoints
    // INCHI✔️❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: ICR endpoint numbers are 1-based AT_NUMB values and are promoted before subtraction.
    // END INCHI ACTIVE MACRO CONFIGURATION: FixRemoveExtraTautEndpoints

    let num_at = pStruct.num_atoms;
    let inverse_mask = !forbidden_edge_mask;
    let mut newly_fixed_edges = EDGE_LIST::default();
    let _ = AllocEdgeList(heap, &mut newly_fixed_edges, EDGE_LIST_CLEAR)?;
    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut ret = AllocEdgeList(heap, &mut newly_fixed_edges, num_at.wrapping_mul(2))?;
        if ret != 0 {
            return Ok(ret);
        }
        let mut atom_number = 0_i32;
        while atom_number < num_at {
            let index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            for edge_number in [
                valence.nCMinusGroupEdge.wrapping_sub(1),
                valence.nCPlusGroupEdge.wrapping_sub(1),
            ] {
                if edge_number >= 0 {
                    ret = AddToEdgeList(heap, &mut newly_fixed_edges, edge_number, 0)?;
                    if ret != 0 {
                        return Ok(ret);
                    }
                    let edge = heap
                        .slice_mut(pBNS.edge)?
                        .get_mut(
                            usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                }
            }
            atom_number = atom_number.wrapping_add(1);
        }

        let mut difference_index = 0_i32;
        while difference_index < picr.num_endp_in1_only {
            let endpoint_number = *picr
                .endp_in1_only
                .get(
                    usize::try_from(difference_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let endpoint1 = i32::from(endpoint_number).wrapping_sub(1);
            let endpoint_index =
                usize::try_from(endpoint1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let endpoint_atom = heap
                .slice(at2.as_const())?
                .get(endpoint_index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let endpoint_valence = pVA
                .get(endpoint_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if endpoint_atom.valence == endpoint_atom.chem_bonds_valence
                || endpoint_valence.nCMinusGroupEdge <= 0
            {
                difference_index = difference_index.wrapping_add(1);
                continue;
            }
            let mut neighbor_order = 0_i32;
            while neighbor_order < i32::from(endpoint_atom.valence) {
                let order = usize::try_from(neighbor_order)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if i32::from(endpoint_atom.bond_type[order]) & BOND_TYPE_MASK as i32
                    == crate::source_types::BOND_TYPE_DOUBLE as i32
                {
                    let centerpoint = usize::from(endpoint_atom.neighbor[order]);
                    let center_atom = heap
                        .slice(at2.as_const())?
                        .get(centerpoint)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let center_valence = pVA
                        .get(centerpoint)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if center_atom.charge != 0
                        || center_valence.nCPlusGroupEdge <= 0
                        || is_centerpoint_elem(center_atom.el_number) == 0
                    {
                        neighbor_order = neighbor_order.wrapping_add(1);
                        continue;
                    }
                    let mut center_order = 0_i32;
                    while center_order < i32::from(center_atom.valence) {
                        let center_slot = usize::try_from(center_order)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if i32::from(center_atom.bond_type[center_slot]) & BOND_TYPE_MASK as i32
                            != BOND_TYPE_SINGLE as i32
                        {
                            center_order = center_order.wrapping_add(1);
                            continue;
                        }
                        let endpoint2 = usize::from(center_atom.neighbor[center_slot]);
                        let endpoint2_restored = heap
                            .slice(at2.as_const())?
                            .get(endpoint2)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if endpoint2_restored.endpoint == 0
                            && heap
                                .slice(atn.as_const())?
                                .get(endpoint2)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .endpoint
                                != 0
                        {
                            break;
                        }
                        center_order = center_order.wrapping_add(1);
                    }
                    if center_order == i32::from(center_atom.valence) {
                        neighbor_order = neighbor_order.wrapping_add(1);
                        continue;
                    }
                    let source_edge_number = center_valence.nCPlusGroupEdge.wrapping_sub(1);
                    let source_edge_index = usize::try_from(source_edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let source_edge = heap
                        .slice(pBNS.edge.as_const())?
                        .get(source_edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if source_edge.flow == 0 {
                        neighbor_order = neighbor_order.wrapping_add(1);
                        continue;
                    }
                    let target_edge_number = endpoint_valence.nCMinusGroupEdge.wrapping_sub(1);
                    let target_edge_index = usize::try_from(target_edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    if heap
                        .slice(pBNS.edge.as_const())?
                        .get(target_edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .flow
                        != 0
                    {
                        neighbor_order = neighbor_order.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(source_edge.neighbor1);
                    let second_vertex = i32::from(source_edge.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[source_edge_index].flow =
                        source_edge.flow.wrapping_sub(1);
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
                    let target = heap
                        .slice_mut(pBNS.edge)?
                        .get_mut(target_edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    target.forbidden = (i32::from(target.forbidden) & inverse_mask) as i8;

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
                    if ret < 0 {
                        return Ok(ret);
                    }
                    if ret == 1
                        && ((path_end == first_vertex && path_start == second_vertex)
                            || (path_end == second_vertex && path_start == first_vertex))
                        && delta_charge <= 1
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
                        if ret < 0 {
                            return Ok(ret);
                        }
                        if ret != 0 {
                            *pnTotalDelta = pnTotalDelta.wrapping_add(ret);
                        } else {
                            ret = RI_ERR_PROGR;
                        }
                        return Ok(ret);
                    }
                    ret = 0;
                    heap.slice_mut(pBNS.edge)?[source_edge_index].flow = source_edge.flow;
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
                    let target = heap
                        .slice_mut(pBNS.edge)?
                        .get_mut(target_edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    target.forbidden = (i32::from(target.forbidden) | forbidden_edge_mask) as i8;
                }
                neighbor_order = neighbor_order.wrapping_add(1);
            }
            difference_index = difference_index.wrapping_add(1);
        }
        Ok(ret)
    })();

    let remove_result =
        RemoveForbiddenEdgeMask(heap, pBNS, &newly_fixed_edges, forbidden_edge_mask);
    let free_result = AllocEdgeList(heap, &mut newly_fixed_edges, EDGE_LIST_FREE);
    match (execution, remove_result, free_result) {
        (Err(error), _, _) => Err(error),
        (Ok(_), Err(error), _) => Err(error),
        (Ok(_), Ok(_), Err(error)) => Err(error),
        (Ok(ret), Ok(()), Ok(_)) => Ok(ret),
    }
}

#[allow(non_snake_case)]
pub(crate) fn FillOutExtraFixedHDataRestr(
    heap: &mut SourceHeap,
    pStruct: &mut StrFromINChI,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:754 FillOutExtraFixedHDataRestr
    /*
    int  FillOutExtraFixedHDataRestr(StrFromINChI* pStruct)
    {
        int i, j, k, len, ret = 0;
        AT_NUMB* pNum;
        for (i = 0; i < TAUT_NUM; i++)
        {
            if (pStruct->pOneINChI_Aux[i])
            {
                pNum = (pStruct->pOneINChI_Aux[i]->nIsotopicOrigAtNosInCanonOrd &&
                    pStruct->pOneINChI_Aux[i]->nIsotopicOrigAtNosInCanonOrd[0]) ?
                    pStruct->pOneINChI_Aux[i]->nIsotopicOrigAtNosInCanonOrd :
                    (pStruct->pOneINChI_Aux[i]->nOrigAtNosInCanonOrd &&
                        pStruct->pOneINChI_Aux[i]->nOrigAtNosInCanonOrd[0]) ?
                    pStruct->pOneINChI_Aux[i]->nOrigAtNosInCanonOrd : NULL;
            }
            else
            {
                pNum = NULL;
            }
            if (pNum)
            {
                len = pStruct->num_atoms * sizeof(pStruct->nCanon2Atno[0][0]);
                if ((!pStruct->nCanon2Atno[i] &&
                    !(pStruct->nCanon2Atno[i] = (AT_NUMB*)inchi_malloc(len))) ||
                    (!pStruct->nAtno2Canon[i] &&
                        !(pStruct->nAtno2Canon[i] = (AT_NUMB*)inchi_malloc(len)))) /* djb-rwth: addressing LLVM warnings */
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }

                INCHI_HEAPCHK

                    memcpy(pStruct->nCanon2Atno[i], pNum, len); /* ??? the next for(...) fills it out */

                INCHI_HEAPCHK

                    for (j = 0; j < pStruct->num_atoms; j++)
                    {
                        k = pNum[j] - 1; /* atom number */
                        pStruct->nCanon2Atno[i][j] = (AT_NUMB)k;
                        pStruct->nAtno2Canon[i][k] = (AT_NUMB)j;
                        INCHI_HEAPCHK
                    }
            }
            else
            {
                if (!i)
                {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                else
                {
                    if (pStruct->nCanon2Atno[i])
                    {
                        inchi_free(pStruct->nCanon2Atno[i]);
                        pStruct->nCanon2Atno[i] = NULL;
                    }
                    INCHI_HEAPCHK
                        if (pStruct->nAtno2Canon[i])
                        {
                            inchi_free(pStruct->nAtno2Canon[i]);
                            pStruct->nAtno2Canon[i] = NULL;
                        }
                    INCHI_HEAPCHK
                }
            }
        }

    exit_function:

        return ret;
    }
    */
    // END INCHI C FUNCTION: FillOutExtraFixedHDataRestr
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillOutExtraFixedHDataRestr
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: TAUT_NUM is 2; AT_NUMB is unsigned short; INCHI_HEAPCHK expands to no code.
    // INCHI✔️❌: inchi_malloc is malloc and inchi_free(X) conditionally calls free(X).
    // INCHI✔️❌: SourceHeap checked pointer lookup adds overhead versus direct C pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: FillOutExtraFixedHDataRestr

    let atom_count = usize::try_from(pStruct.num_atoms)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    pStruct
        .num_atoms
        .checked_mul(std::mem::size_of::<AT_NUMB>() as i32)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;

    for i in 0..TAUT_NUM as usize {
        let p_num = if pStruct.pOneINChI_Aux[i].is_null() {
            SourceMutPointer::null()
        } else {
            let aux = heap
                .slice(pStruct.pOneINChI_Aux[i].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if !aux.nIsotopicOrigAtNosInCanonOrd.is_null()
                && heap
                    .slice(aux.nIsotopicOrigAtNosInCanonOrd.as_const())?
                    .first()
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
            {
                aux.nIsotopicOrigAtNosInCanonOrd
            } else if !aux.nOrigAtNosInCanonOrd.is_null()
                && heap
                    .slice(aux.nOrigAtNosInCanonOrd.as_const())?
                    .first()
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
            {
                aux.nOrigAtNosInCanonOrd
            } else {
                SourceMutPointer::null()
            }
        };

        if !p_num.is_null() {
            if pStruct.nCanon2Atno[i].is_null() {
                pStruct.nCanon2Atno[i] = match heap.allocate(vec![AT_NUMB::default(); atom_count]) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => return Ok(RI_ERR_ALLOC),
                    Err(error) => return Err(error),
                };
            }
            if pStruct.nAtno2Canon[i].is_null() {
                pStruct.nAtno2Canon[i] = match heap.allocate(vec![AT_NUMB::default(); atom_count]) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => return Ok(RI_ERR_ALLOC),
                    Err(error) => return Err(error),
                };
            }

            let canonical_numbers = heap
                .slice(p_num.as_const())?
                .get(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(pStruct.nCanon2Atno[i])?
                .get_mut(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(&canonical_numbers);

            for (j, atom_number) in canonical_numbers.into_iter().enumerate() {
                let k = i32::from(atom_number).wrapping_sub(1);
                heap.slice_mut(pStruct.nCanon2Atno[i])?[j] = k as AT_NUMB;
                let inverse_index =
                    usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                *heap
                    .slice_mut(pStruct.nAtno2Canon[i])?
                    .get_mut(inverse_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = j as AT_NUMB;
            }
        } else if i == 0 {
            return Ok(RI_ERR_PROGR);
        } else {
            if !pStruct.nCanon2Atno[i].is_null() {
                inchi_free(heap, pStruct.nCanon2Atno[i])?;
                pStruct.nCanon2Atno[i] = SourceMutPointer::null();
            }
            if !pStruct.nAtno2Canon[i].is_null() {
                inchi_free(heap, pStruct.nAtno2Canon[i])?;
                pStruct.nAtno2Canon[i] = SourceMutPointer::null();
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn FillOutExtraFixedHDataInChI(
    heap: &mut SourceHeap,
    pStruct: &mut StrFromINChI,
    pInChI: [SourceMutPointer<INChI>; 2],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:831 FillOutExtraFixedHDataInChI
    /*
    int  FillOutExtraFixedHDataInChI(StrFromINChI* pStruct, INChI* pInChI[])
    {
        int ret = 0;
        /*--- allocate memory for Mobile/Fixed-H data from the input InChI ---*/
        if (NULL == pStruct->endpoint)
        {
            pStruct->endpoint = (AT_NUMB*)inchi_calloc(pStruct->num_atoms, sizeof(pStruct->endpoint[0]));
        }
        else
        {
            memset(pStruct->endpoint, 0, pStruct->num_atoms * sizeof(pStruct->endpoint[0])); /* djb-rwth: memset_s C11/Annex K variant? */
        }
        if (NULL == pStruct->fixed_H)
        {
            pStruct->fixed_H = (S_CHAR*)inchi_malloc(pStruct->num_atoms * sizeof(pStruct->fixed_H[0]));
        }
        if (!pStruct->endpoint || !pStruct->fixed_H)
        {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
        /*--- fill out Mobile/Fixed-H data from the input InChI ---*/
        GetTgroupInfoFromInChI(&pStruct->ti, NULL, pStruct->endpoint, pInChI[1]);
        if (pInChI[0]->nNum_H_fixed)
        {
            memcpy(pStruct->fixed_H, pInChI[0]->nNum_H_fixed, pStruct->num_atoms * sizeof(pStruct->fixed_H[0]));
        }
        else
        {
            memset(pStruct->fixed_H, 0, pStruct->num_atoms * sizeof(pStruct->fixed_H[0])); /* djb-rwth: memset_s C11/Annex K variant? */
        }

    exit_function:

        return ret;
    }
    */
    // END INCHI C FUNCTION: FillOutExtraFixedHDataInChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillOutExtraFixedHDataInChI
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔❌: inchi_calloc(X,Y)=calloc(X,Y); inchi_malloc(X)=malloc(X).
    // INCHI✔❌: SourceHeap checked ownership and allocation-map lookup add overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: FillOutExtraFixedHDataInChI

    let atom_count = u64::try_from(pStruct.num_atoms)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let atom_count_usize = usize::try_from(atom_count)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;

    if pStruct.endpoint.is_null() {
        pStruct.endpoint = match inchi_calloc::<AT_NUMB>(
            heap,
            atom_count,
            std::mem::size_of::<AT_NUMB>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    } else {
        heap.slice_mut(pStruct.endpoint)?
            .get_mut(..atom_count_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(0);
    }

    if pStruct.fixed_H.is_null() {
        pStruct.fixed_H = match inchi_malloc(heap, atom_count) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }
    if pStruct.endpoint.is_null() || pStruct.fixed_H.is_null() {
        return Ok(RI_ERR_ALLOC);
    }

    let mobile_inchi = if pInChI[1].is_null() {
        None
    } else {
        Some(
            heap.slice(pInChI[1].as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let _ = GetTgroupInfoFromInChI(
        heap,
        &mut pStruct.ti,
        SourceMutPointer::null(),
        pStruct.endpoint,
        mobile_inchi.as_ref(),
    )?;

    let fixed_inchi = heap
        .slice(pInChI[0].as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if fixed_inchi.nNum_H_fixed.is_null() {
        heap.slice_mut(pStruct.fixed_H)?
            .get_mut(..atom_count_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(0);
    } else {
        let fixed_h = heap
            .slice(fixed_inchi.nNum_H_fixed.as_const())?
            .get(..atom_count_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        heap.slice_mut(pStruct.fixed_H)?
            .get_mut(..atom_count_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .copy_from_slice(&fixed_h);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn FillOutCMP2FHINCHI(
    heap: &SourceHeap,
    pStruct: &StrFromINChI,
    at2: &[inp_ATOM],
    pVA: &[VAL_AT],
    pInChI: [SourceMutPointer<INChI>; 2],
    pc2i: &mut CMP2FHINCHI,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:870 FillOutCMP2FHINCHI
    /*
    int FillOutCMP2FHINCHI(StrFromINChI* pStruct,
        inp_ATOM* at2,
        VAL_AT* pVA,
        INChI* pInChI[],
        CMP2FHINCHI* pc2i)
    {
        int       ret = 0, i, j;
        int       bFixHRevrsExists = pInChI[1] && pInChI[1]->nNumberOfAtoms > 0 && !pInChI[1]->bDeleted;
        inp_ATOM* at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at) ? pStruct->pOne_norm_data[1]->at : NULL;
        S_CHAR* num_Fixed_H_Revrs = pStruct->pOneINChI[0]->nNum_H_fixed ? pStruct->pOneINChI[0]->nNum_H_fixed : NULL;
        /* atom number in structure that produced original InChI is atom number in all inp_ATOM *atoms */
        /* atom number in structure that produced restored InChI is in nAtomRevrs[]: */
        AT_NUMB* nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
        S_CHAR* pnMobHInChI = (pInChI[1] && pInChI[1]->nNum_H) ? pInChI[1]->nNum_H :
            (pInChI[0] && pInChI[0]->nNum_H) ? pInChI[0]->nNum_H : NULL;
        S_CHAR* pnMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H) ?
            pStruct->pOneINChI[1]->nNum_H :
            (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H) ?
            pStruct->pOneINChI[0]->nNum_H : NULL;
        int     nNumTgHInChI, nNumTgMInChI, nNumTgHRevrs, nNumTgMRevrs;
        memset(pc2i, 0, sizeof(*pc2i)); /* djb-rwth: memset_s C11/Annex K variant? */
        pc2i->nNumTgInChI = pStruct->ti.num_t_groups;
        pc2i->nNumTgRevrs = pStruct->One_ti.num_t_groups;
        pc2i->bHasDifference |= pc2i->nNumTgInChI != pc2i->nNumTgRevrs;

        pc2i->nNumRemHInChI = pStruct->nNumRemovedProtonsMobHInChI;
        pc2i->nNumRemHRevrs = pStruct->One_ti.tni.nNumRemovedProtons;
        pc2i->bHasDifference |= pc2i->nNumRemHInChI != pc2i->nNumRemHRevrs;

        pc2i->bFixedHLayerExistsRevrs = bFixHRevrsExists;
        pc2i->bHasDifference |= !bFixHRevrsExists;

        for (i = 0; i < pStruct->ti.num_t_groups && i < pStruct->One_ti.num_t_groups; i++)
        {
            nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
            nNumTgMInChI = pStruct->ti.t_group[i].num[1];
            nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
            nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];

            pc2i->bHasDifference |= nNumTgHInChI != nNumTgHRevrs;
            pc2i->bHasDifference |= nNumTgMInChI != nNumTgMRevrs;

            if (pStruct->ti.t_group[i].nNumEndpoints ==
                pStruct->One_ti.t_group[i].nNumEndpoints)
            {

                if (nNumTgHInChI != nNumTgHRevrs)
                {
                    pc2i->nNumTgDiffH++;
                }
                if (nNumTgMInChI != nNumTgMRevrs)
                {
                    pc2i->nNumTgDiffMinus++;
                }
            }
            pc2i->bHasDifference |= pStruct->ti.t_group[i].nNumEndpoints !=
                pStruct->One_ti.t_group[i].nNumEndpoints;

            pc2i->nNumTgHInChI += nNumTgHInChI;
            pc2i->nNumTgMInChI += nNumTgMInChI;
            pc2i->nNumTgHRevrs += nNumTgHRevrs;
            pc2i->nNumTgMRevrs += nNumTgMRevrs;
        }
        for (; i < pStruct->ti.num_t_groups; i++)
        {
            nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
            nNumTgMInChI = pStruct->ti.t_group[i].num[1];
            pc2i->nNumTgHInChI += nNumTgHInChI;
            pc2i->nNumTgMInChI += nNumTgMInChI;
            pc2i->bHasDifference |= 1;
        }
        for (; i < pStruct->One_ti.num_t_groups; i++)
        {
            nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
            nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];
            pc2i->nNumTgHRevrs += nNumTgHRevrs;
            pc2i->nNumTgMRevrs += nNumTgMRevrs;
            pc2i->bHasDifference |= 1;
        }
        for (i = j = 0; i < pStruct->num_atoms; i++)
        {
            /* i = original InChI canonical number - 1 */
            /* k = atom number from InChI created out of restored Fixed-H structure */
            int iCanonRevrs = nAtno2CanonRevrs[i];
            int endptInChI = pStruct->endpoint[i]; /* endpoint in InChI */
            int endptRevrs = at_Mobile_H_Revrs ? at_Mobile_H_Revrs[i].endpoint : 0;
            int nFixHInChI = pStruct->fixed_H[i];
            int nFixHRevrs = num_Fixed_H_Revrs ? num_Fixed_H_Revrs[iCanonRevrs] : 0;
            int nMobHInChI = pnMobHInChI ? pnMobHInChI[i] : 0;
            int nMobHRevrs = pnMobHRevrs ? pnMobHRevrs[iCanonRevrs] : 0;
            if ( /*(!endptInChI || !endptRevrs) &&*/ (nFixHInChI != nFixHRevrs) ||
                (!endptInChI != !endptRevrs) || nMobHInChI != nMobHRevrs)
            {
                /* in InChI or reversed InChI atom[i] is not tautomeric */
                /* and number of fixed-H on the atom[i] differs */
                if (j >= MAX_DIFF_FIXH)
                {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                pc2i->c2at[j].endptInChI = endptInChI;
                pc2i->c2at[j].endptRevrs = endptRevrs;
                pc2i->bHasDifference |= !endptInChI != !endptRevrs;
                pc2i->c2at[j].atomNumber = i;
                pc2i->c2at[j].nValElectr = pVA[i].cNumValenceElectrons;
                pc2i->c2at[j].nPeriodNum = pVA[i].cPeriodicRowNumber;
                pc2i->c2at[j].nFixHInChI = nFixHInChI;
                pc2i->c2at[j].nFixHRevrs = nFixHRevrs;
                pc2i->bHasDifference |= nFixHInChI != nFixHRevrs;
                pc2i->c2at[j].nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H[i] :
                    pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H[i] : 0;
                pc2i->c2at[j].nMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H) ?
                    pStruct->pOneINChI[1]->nNum_H[iCanonRevrs] :
                    (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H) ?
                    pStruct->pOneINChI[0]->nNum_H[iCanonRevrs] : 0;
                pc2i->nNumDiffMobH += (nMobHInChI != nMobHRevrs && !endptRevrs && !endptInChI);
                pc2i->bHasDifference |= nMobHInChI != nMobHRevrs;
                pc2i->c2at[j].nNumHRevrs = at2[i].num_H;
                pc2i->c2at[j].nAtChargeRevrs = at2[i].charge;
                j++;
            }
            pc2i->nNumEndpInChI += (endptInChI != 0);
            pc2i->nNumEndpRevrs += (endptRevrs != 0);

            if (!pVA[i].cMetal)
            {
                pc2i->nChargeFixHRevrsNonMetal += at2[i].charge;
                pc2i->nChargeMobHRevrsNonMetal += at_Mobile_H_Revrs ? at_Mobile_H_Revrs[i].charge : 0;
            }

            /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
        }
        pc2i->nChargeFixHInChI = pInChI[0] ? pInChI[0]->nTotalCharge : 0;
        pc2i->nChargeMobHInChI = pInChI[1] ? pInChI[1]->nTotalCharge : 0;

        pc2i->nChargeMobHRevrs = pStruct->pOneINChI[1] ? pStruct->pOneINChI[1]->nTotalCharge :
            pStruct->pOneINChI[0] ? pStruct->pOneINChI[0]->nTotalCharge : 0;
        pc2i->nChargeFixHRevrs = pStruct->pOneINChI[0] ? pStruct->pOneINChI[0]->nTotalCharge : 0;

        pc2i->bHasDifference |= pc2i->nChargeFixHInChI != pc2i->nChargeFixHRevrs;
        pc2i->bHasDifference |= pc2i->nChargeMobHInChI != pc2i->nChargeMobHRevrs;

    exit_function:
        pc2i->len_c2at = j;

        return ret;
    }
    */
    // END INCHI C FUNCTION: FillOutCMP2FHINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillOutCMP2FHINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: MAX_DIFF_FIXH is 256; memset zeroes the complete output structure.
    // INCHI✔️❌: GCC narrowing to short, signed char, and unsigned char is reproduced with wrapping casts.
    // INCHI✔️❌: SourceHeap checked pointer lookup adds overhead versus direct C pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: FillOutCMP2FHINCHI

    let input = [
        if pInChI[0].is_null() {
            None
        } else {
            Some(
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        },
        if pInChI[1].is_null() {
            None
        } else {
            Some(
                heap.slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        },
    ];
    let reversed_fixed = heap
        .slice(pStruct.pOneINChI[0].as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let reversed_mobile = if pStruct.pOneINChI[1].is_null() {
        None
    } else {
        Some(
            heap.slice(pStruct.pOneINChI[1].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        )
    };
    let mobile_atoms = if pStruct.pOne_norm_data[1].is_null() {
        SourceMutPointer::null()
    } else {
        heap.slice(pStruct.pOne_norm_data[1].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at
    };
    let fixed_h_reversed = reversed_fixed.nNum_H_fixed;
    let mobile_h_input = input[1]
        .as_ref()
        .filter(|inchi| !inchi.nNum_H.is_null())
        .map_or_else(
            || {
                input[0]
                    .as_ref()
                    .map_or(SourceMutPointer::null(), |inchi| inchi.nNum_H)
            },
            |inchi| inchi.nNum_H,
        );
    let mobile_h_reversed = reversed_mobile
        .as_ref()
        .filter(|inchi| !inchi.nNum_H.is_null())
        .map_or(reversed_fixed.nNum_H, |inchi| inchi.nNum_H);
    let fixed_layer_exists = input[1]
        .as_ref()
        .is_some_and(|inchi| inchi.nNumberOfAtoms > 0 && inchi.bDeleted == 0);

    *pc2i = CMP2FHINCHI::default();
    pc2i.nNumTgInChI = pStruct.ti.num_t_groups as i16;
    pc2i.nNumTgRevrs = pStruct.One_ti.num_t_groups as i16;
    pc2i.bHasDifference |= (pc2i.nNumTgInChI != pc2i.nNumTgRevrs) as i8;
    pc2i.nNumRemHInChI = pStruct.nNumRemovedProtonsMobHInChI as i16;
    pc2i.nNumRemHRevrs = pStruct.One_ti.tni.nNumRemovedProtons as i16;
    pc2i.bHasDifference |= (pc2i.nNumRemHInChI != pc2i.nNumRemHRevrs) as i8;
    pc2i.bFixedHLayerExistsRevrs = fixed_layer_exists as i8;
    pc2i.bHasDifference |= (!fixed_layer_exists) as i8;

    let mut i = 0_i32;
    while i < pStruct.ti.num_t_groups && i < pStruct.One_ti.num_t_groups {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let original = heap
            .slice(pStruct.ti.t_group.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let reversed = heap
            .slice(pStruct.One_ti.t_group.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let h_input = i32::from(original.num[0]).wrapping_sub(i32::from(original.num[1]));
        let minus_input = i32::from(original.num[1]);
        let h_reversed = i32::from(reversed.num[0]).wrapping_sub(i32::from(reversed.num[1]));
        let minus_reversed = i32::from(reversed.num[1]);
        pc2i.bHasDifference |= (h_input != h_reversed) as i8;
        pc2i.bHasDifference |= (minus_input != minus_reversed) as i8;
        if original.nNumEndpoints == reversed.nNumEndpoints {
            if h_input != h_reversed {
                pc2i.nNumTgDiffH = pc2i.nNumTgDiffH.wrapping_add(1);
            }
            if minus_input != minus_reversed {
                pc2i.nNumTgDiffMinus = pc2i.nNumTgDiffMinus.wrapping_add(1);
            }
        }
        pc2i.bHasDifference |= (original.nNumEndpoints != reversed.nNumEndpoints) as i8;
        pc2i.nNumTgHInChI = pc2i.nNumTgHInChI.wrapping_add(h_input as i16);
        pc2i.nNumTgMInChI = pc2i.nNumTgMInChI.wrapping_add(minus_input as i16);
        pc2i.nNumTgHRevrs = pc2i.nNumTgHRevrs.wrapping_add(h_reversed as i16);
        pc2i.nNumTgMRevrs = pc2i.nNumTgMRevrs.wrapping_add(minus_reversed as i16);
        i = i.wrapping_add(1);
    }
    while i < pStruct.ti.num_t_groups {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let original = heap
            .slice(pStruct.ti.t_group.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let h_input = i32::from(original.num[0]).wrapping_sub(i32::from(original.num[1]));
        pc2i.nNumTgHInChI = pc2i.nNumTgHInChI.wrapping_add(h_input as i16);
        pc2i.nNumTgMInChI = pc2i.nNumTgMInChI.wrapping_add(original.num[1] as i16);
        pc2i.bHasDifference |= 1;
        i = i.wrapping_add(1);
    }
    while i < pStruct.One_ti.num_t_groups {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let reversed = heap
            .slice(pStruct.One_ti.t_group.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let h_reversed = i32::from(reversed.num[0]).wrapping_sub(i32::from(reversed.num[1]));
        pc2i.nNumTgHRevrs = pc2i.nNumTgHRevrs.wrapping_add(h_reversed as i16);
        pc2i.nNumTgMRevrs = pc2i.nNumTgMRevrs.wrapping_add(reversed.num[1] as i16);
        pc2i.bHasDifference |= 1;
        i = i.wrapping_add(1);
    }

    let mut j = 0_i32;
    i = 0;
    while i < pStruct.num_atoms {
        let atom_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let canonical_index = usize::from(
            *heap
                .slice(pStruct.nAtno2Canon[0].as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let endpoint_input = i32::from(
            *heap
                .slice(pStruct.endpoint.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let endpoint_reversed = if mobile_atoms.is_null() {
            0
        } else {
            i32::from(
                heap.slice(mobile_atoms.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .endpoint,
            )
        };
        let fixed_input = i32::from(
            *heap
                .slice(pStruct.fixed_H.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let fixed_reversed = if fixed_h_reversed.is_null() {
            0
        } else {
            i32::from(
                *heap
                    .slice(fixed_h_reversed.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let mobile_input = if mobile_h_input.is_null() {
            0
        } else {
            i32::from(
                *heap
                    .slice(mobile_h_input.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let mobile_reversed_value = if mobile_h_reversed.is_null() {
            0
        } else {
            i32::from(
                *heap
                    .slice(mobile_h_reversed.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };

        if fixed_input != fixed_reversed
            || (endpoint_input == 0) != (endpoint_reversed == 0)
            || mobile_input != mobile_reversed_value
        {
            if j >= MAX_DIFF_FIXH as i32 {
                pc2i.len_c2at = j as i16;
                return Ok(RI_ERR_PROGR);
            }
            let difference = &mut pc2i.c2at[j as usize];
            difference.endptInChI = endpoint_input as AT_NUMB;
            difference.endptRevrs = endpoint_reversed as AT_NUMB;
            pc2i.bHasDifference |= ((endpoint_input == 0) != (endpoint_reversed == 0)) as i8;
            difference.atomNumber = i as AT_NUMB;
            let valence = pVA
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            difference.nValElectr = valence.cNumValenceElectrons as u8;
            difference.nPeriodNum = valence.cPeriodicRowNumber as u8;
            difference.nFixHInChI = fixed_input as i8;
            difference.nFixHRevrs = fixed_reversed as i8;
            pc2i.bHasDifference |= (fixed_input != fixed_reversed) as i8;
            difference.nMobHInChI = mobile_input as i8;
            difference.nMobHRevrs = mobile_reversed_value as i8;
            pc2i.nNumDiffMobH = pc2i.nNumDiffMobH.wrapping_add(
                (mobile_input != mobile_reversed_value
                    && endpoint_reversed == 0
                    && endpoint_input == 0) as u8,
            );
            pc2i.bHasDifference |= (mobile_input != mobile_reversed_value) as i8;
            let atom = at2
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            difference.nNumHRevrs = atom.num_H;
            difference.nAtChargeRevrs = atom.charge;
            j = j.wrapping_add(1);
        }
        pc2i.nNumEndpInChI = pc2i
            .nNumEndpInChI
            .wrapping_add((endpoint_input != 0) as i16);
        pc2i.nNumEndpRevrs = pc2i
            .nNumEndpRevrs
            .wrapping_add((endpoint_reversed != 0) as i16);
        let valence = pVA
            .get(atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if valence.cMetal == 0 {
            let atom = at2
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            pc2i.nChargeFixHRevrsNonMetal = pc2i.nChargeFixHRevrsNonMetal.wrapping_add(atom.charge);
            let mobile_charge = if mobile_atoms.is_null() {
                0
            } else {
                heap.slice(mobile_atoms.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .charge
            };
            pc2i.nChargeMobHRevrsNonMetal =
                pc2i.nChargeMobHRevrsNonMetal.wrapping_add(mobile_charge);
        }
        i = i.wrapping_add(1);
    }

    pc2i.nChargeFixHInChI = input[0]
        .as_ref()
        .map_or(0, |inchi| inchi.nTotalCharge as i8);
    pc2i.nChargeMobHInChI = input[1]
        .as_ref()
        .map_or(0, |inchi| inchi.nTotalCharge as i8);
    pc2i.nChargeMobHRevrs = reversed_mobile
        .as_ref()
        .map_or(reversed_fixed.nTotalCharge as i8, |inchi| {
            inchi.nTotalCharge as i8
        });
    pc2i.nChargeFixHRevrs = reversed_fixed.nTotalCharge as i8;
    pc2i.bHasDifference |= (pc2i.nChargeFixHInChI != pc2i.nChargeFixHRevrs) as i8;
    pc2i.bHasDifference |= (pc2i.nChargeMobHInChI != pc2i.nChargeMobHRevrs) as i8;
    pc2i.len_c2at = j as i16;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn FillOutCMP2MHINCHI(
    heap: &SourceHeap,
    pStruct: &StrFromINChI,
    pTCGroups: &ALL_TC_GROUPS,
    at2: &[inp_ATOM],
    pVA: &[VAL_AT],
    pInChI: [SourceMutPointer<INChI>; 2],
    pc2i: &mut CMP2MHINCHI,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1021 FillOutCMP2MHINCHI
    // INCHI✔️❌: complete active source frame follows verbatim; SourceHeap pointer checks add known overhead.
    /*
    int FillOutCMP2MHINCHI(StrFromINChI* pStruct,
        ALL_TC_GROUPS* pTCGroups,
        inp_ATOM* at2,
        VAL_AT* pVA,
        INChI* pInChI[],
        CMP2MHINCHI* pc2i)
    {
        int       ret = 0, i, j, iat;
        int       bFixHRevrsExists = pInChI[1] && pInChI[1]->nNumberOfAtoms > 0 && !pInChI[1]->bDeleted;
        inp_ATOM* at_Mobile_H_Revrs = (pStruct->pOne_norm_data[0] &&
            pStruct->pOne_norm_data[0]->at) ? pStruct->pOne_norm_data[0]->at : NULL;
        /* atom number in structure that produced original InChI is atom number in all inp_ATOM *atoms */
        /* atom number in structure that produced restored InChI is in nAtomRevrs[]: */
        AT_NUMB* nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        AT_NUMB* nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
        S_CHAR* pnMobHInChI = (pInChI[0] && pInChI[0]->nNum_H) ? pInChI[0]->nNum_H : NULL;
        S_CHAR* pnMobHRevrs = (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H) ?
            pStruct->pOneINChI[0]->nNum_H : NULL;
        int     nNumTgHInChI, nNumTgMInChI, nNumTgHRevrs, nNumTgMRevrs;
        memset(pc2i, 0, sizeof(*pc2i)); /* djb-rwth: memset_s C11/Annex K variant? */
        pc2i->nNumTgInChI = pStruct->ti.num_t_groups;
        pc2i->nNumTgRevrs = pStruct->One_ti.num_t_groups;
        pc2i->bHasDifference |= pc2i->nNumTgInChI != pc2i->nNumTgRevrs;

        pc2i->nNumRemHInChI = pStruct->nNumRemovedProtonsMobHInChI;
        pc2i->nNumRemHRevrs = pStruct->One_ti.tni.nNumRemovedProtons;
        /*pc2i->bHasDifference |= pc2i->nNumRemHInChI != pc2i->nNumRemHRevrs;*/

        pc2i->bFixedHLayerExistsRevrs = bFixHRevrsExists;
        /*pc2i->bHasDifference |= !bFixHRevrsExists;*/

        for (i = 0; i < pStruct->ti.num_t_groups; i++)
        {
            int jFst = pStruct->ti.t_group[i].nFirstEndpointAtNoPos;
            int jNum = pStruct->ti.t_group[i].nNumEndpoints;
            int is_N, is_O;
            for (j = 0; j < jNum; j++)
            {
                iat = pStruct->ti.nEndpointAtomNumber[jFst + j];
                is_N = pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1;
                is_O = pVA[iat].cNumValenceElectrons == 6;
                if (is_N + is_O != 1)
                {
                    return RI_ERR_SYNTAX;
                }
                pc2i->nNumTgNInChI += is_N;
                pc2i->nNumTgOInChI += is_O;
                if (at2[iat].chem_bonds_valence == at2[iat].valence)
                {
                    /* donor */
                    if (is_N)
                    {
                        /* N */
                        pc2i->nNumTgNHInChI += at2[iat].charge == 0 && at2[iat].num_H == 1;
                        pc2i->nNumTgNH2InChI += at2[iat].charge == 0 && at2[iat].num_H == 2;
                        pc2i->nNumTgNMinusInChI += at2[iat].charge == -1 && at2[iat].num_H == 0;
                        pc2i->nNumTgNHMinusInChI += at2[iat].charge == -1 && at2[iat].num_H == 1;
                    }
                    else
                    {
                        /* O, S, Se, Te */
                        pc2i->nNumTgOHInChI += at2[iat].charge == 0 && at2[iat].num_H == 1;
                        pc2i->nNumTgOMinusInChI += at2[iat].charge == -1 && at2[iat].num_H == 0;
                    }
                }
                else
                {
                    if (at2[iat].chem_bonds_valence == at2[iat].valence + 1)
                    {
                        /* donor */
                        if (is_N)
                        {
                            /* N */
                            pc2i->nNumTgDBNHInChI += at2[iat].charge == 0 && at2[iat].num_H == 1;
                            pc2i->nNumTgDBNMinusInChI += at2[iat].charge == -1 && at2[iat].num_H == 0;
                            pc2i->nNumTgDBNInChI += at2[iat].charge == 0 && at2[iat].num_H == 0;
                        }
                        else
                        {
                            /* O, S, Se, Te */
                            pc2i->nNumTgDBOInChI += at2[iat].charge == 0 && at2[iat].num_H == 0;
                        }
                    }
                }
            }
        }
        for (i = 0; i < pStruct->One_ti.num_t_groups; i++)
        {
            int jFst = pStruct->One_ti.t_group[i].nFirstEndpointAtNoPos;
            int jNum = pStruct->One_ti.t_group[i].nNumEndpoints;
            int is_N, is_O;
            for (j = 0; j < jNum; j++)
            {
                iat = nCanon2AtnoRevrs[(int)pStruct->One_ti.nEndpointAtomNumber[jFst + j]];
                is_N = pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1;
                is_O = pVA[iat].cNumValenceElectrons == 6;
                if (is_N + is_O != 1)
                {
                    return RI_ERR_PROGR;
                }
                pc2i->nNumTgNRevrs += is_N;
                pc2i->nNumTgORevrs += is_O;
                if (at2[iat].chem_bonds_valence == at2[iat].valence)
                {
                    /* donor */
                    if (is_N)
                    {
                        /* N */
                        pc2i->nNumTgNHRevrs += at2[iat].charge == 0 && at2[iat].num_H == 1;
                        pc2i->nNumTgNH2Revrs += at2[iat].charge == 0 && at2[iat].num_H == 2;
                        pc2i->nNumTgNMinusRevrs += at2[iat].charge == -1 && at2[iat].num_H == 0;
                        pc2i->nNumTgNHMinusRevrs += at2[iat].charge == -1 && at2[iat].num_H == 1;
                    }
                    else
                    {
                        /* O, S, Se, Te */
                        pc2i->nNumTgOHRevrs += at2[iat].charge == 0 && at2[iat].num_H == 1;
                        pc2i->nNumTgOMinusRevrs += at2[iat].charge == -1 && at2[iat].num_H == 0;
                    }
                }
                else
                {
                    if (at2[iat].chem_bonds_valence == at2[iat].valence + 1)
                    {
                        /* donor */
                        if (is_N)
                        {
                            /* N */
                            pc2i->nNumTgDBNHRevrs += at2[iat].charge == 0 && at2[iat].num_H == 1;
                            pc2i->nNumTgDBNMinusRevrs += at2[iat].charge == -1 && at2[iat].num_H == 0;
                            pc2i->nNumTgDBNRevrs += at2[iat].charge == 0 && at2[iat].num_H == 0;
                        }
                        else
                        {
                            /* O, S, Se, Te */
                            pc2i->nNumTgDBORevrs += at2[iat].charge == 0 && at2[iat].num_H == 0;
                        }
                    }
                }
            }
        }

        for (i = 0; i < pStruct->ti.num_t_groups && i < pStruct->One_ti.num_t_groups; i++)
        {
            nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
            nNumTgMInChI = pStruct->ti.t_group[i].num[1];
            nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
            nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];

            pc2i->bHasDifference |= nNumTgHInChI != nNumTgHRevrs;
            pc2i->bHasDifference |= nNumTgMInChI != nNumTgMRevrs;

            if (pStruct->ti.t_group[i].nNumEndpoints ==
                pStruct->One_ti.t_group[i].nNumEndpoints)
            {

                if (nNumTgHInChI != nNumTgHRevrs)
                {
                    pc2i->nNumTgDiffH++;
                }
                if (nNumTgMInChI != nNumTgMRevrs)
                {
                    pc2i->nNumTgDiffMinus++;
                }
            }
            pc2i->bHasDifference |= pStruct->ti.t_group[i].nNumEndpoints !=
                pStruct->One_ti.t_group[i].nNumEndpoints;

            pc2i->nNumTgHInChI += nNumTgHInChI;
            pc2i->nNumTgMInChI += nNumTgMInChI;
            pc2i->nNumTgHRevrs += nNumTgHRevrs;
            pc2i->nNumTgMRevrs += nNumTgMRevrs;
        }
        for (; i < pStruct->ti.num_t_groups; i++)
        {
            nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
            nNumTgMInChI = pStruct->ti.t_group[i].num[1];
            pc2i->nNumTgHInChI += nNumTgHInChI;
            pc2i->nNumTgMInChI += nNumTgMInChI;
            pc2i->bHasDifference |= 1;
        }
        for (; i < pStruct->One_ti.num_t_groups; i++)
        {
            nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
            nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];
            pc2i->nNumTgHRevrs += nNumTgHRevrs;
            pc2i->nNumTgMRevrs += nNumTgMRevrs;
            pc2i->bHasDifference |= 1;
        }
        for (i = j = 0; i < pStruct->num_atoms; i++)
        {
            /* i = original InChI canonical number - 1 */
            /* k = atom number from InChI created out of restored Fixed-H structure */
            int iCanonRevrs = nAtno2CanonRevrs[i];
            int endptInChI = at2[i].endpoint; /* endpoint in InChI */
            int endptRevrs = at_Mobile_H_Revrs ? at_Mobile_H_Revrs[i].endpoint : 0;
            int nMobHInChI = pnMobHInChI ? pnMobHInChI[i] : 0;
            int nMobHRevrs = pnMobHRevrs ? pnMobHRevrs[iCanonRevrs] : 0;
            if ((!endptInChI != !endptRevrs) || nMobHInChI != nMobHRevrs)
            {
                /* in InChI or reversed InChI atom[i] is not tautomeric */
                /* and number of fixed-H on the atom[i] differs */
                if (j >= MAX_DIFF_FIXH)
                {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                pc2i->c2at[j].endptInChI = endptInChI;
                pc2i->c2at[j].endptRevrs = endptRevrs;
                pc2i->bHasDifference |= !endptInChI != !endptRevrs;
                pc2i->c2at[j].atomNumber = i;
                pc2i->c2at[j].nValElectr = pVA[i].cNumValenceElectrons;
                pc2i->c2at[j].nPeriodNum = pVA[i].cPeriodicRowNumber;
                pc2i->c2at[j].nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H[i] :
                    pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H[i] : 0;
                pc2i->c2at[j].nMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H) ?
                    pStruct->pOneINChI[1]->nNum_H[iCanonRevrs] :
                    (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H) ?
                    pStruct->pOneINChI[0]->nNum_H[iCanonRevrs] : 0;

                pc2i->nNumDiffMobH += (nMobHInChI != nMobHRevrs && !endptRevrs && !endptInChI);
                pc2i->bHasDifference |= (nMobHInChI != nMobHRevrs);
                pc2i->c2at[j].nNumHRevrs = at2[i].num_H;
                pc2i->c2at[j].nAtChargeRevrs = at2[i].charge;
                j++;
            }
            pc2i->nNumEndpInChI += (endptInChI != 0);
            pc2i->nNumEndpRevrs += (endptRevrs != 0);

            if (!pVA[i].cMetal)
            {
                pc2i->nChargeMobHRevrsNonMetal += (at_Mobile_H_Revrs && !at_Mobile_H_Revrs[i].endpoint) ? at_Mobile_H_Revrs[i].charge : 0;
            }


            /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
        }
        pc2i->nChargeMobHRevrsNonMetal += pTCGroups->tgroup_charge;

        pc2i->nChargeMobHInChI = pInChI[0] ? pInChI[0]->nTotalCharge : 0;

        pc2i->nChargeMobHRevrs = pStruct->pOneINChI[0] ? pStruct->pOneINChI[0]->nTotalCharge : 0;

        pc2i->bHasDifference |= pc2i->nChargeMobHInChI != pc2i->nChargeMobHRevrs;

    exit_function:
        pc2i->len_c2at = j;

        return ret;
    }
    */
    // END INCHI C FUNCTION: FillOutCMP2MHINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillOutCMP2MHINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: MAX_DIFF_FIXH=256; memset zeroes every output field; short/char compound assignments narrow after promotion.
    // END INCHI ACTIVE MACRO CONFIGURATION: FillOutCMP2MHINCHI

    fn classify(output: &mut CMP2MHINCHI, atom: &inp_ATOM, is_n: bool, reversed: bool) {
        let is_o = !is_n;
        macro_rules! add {
            ($field:ident, $value:expr) => {
                output.$field = output.$field.wrapping_add(i16::from($value));
            };
        }
        if reversed {
            add!(nNumTgNRevrs, is_n);
            add!(nNumTgORevrs, is_o);
            if atom.chem_bonds_valence == atom.valence {
                if is_n {
                    add!(nNumTgNHRevrs, atom.charge == 0 && atom.num_H == 1);
                    add!(nNumTgNH2Revrs, atom.charge == 0 && atom.num_H == 2);
                    add!(nNumTgNMinusRevrs, atom.charge == -1 && atom.num_H == 0);
                    add!(nNumTgNHMinusRevrs, atom.charge == -1 && atom.num_H == 1);
                } else {
                    add!(nNumTgOHRevrs, atom.charge == 0 && atom.num_H == 1);
                    add!(nNumTgOMinusRevrs, atom.charge == -1 && atom.num_H == 0);
                }
            } else if i32::from(atom.chem_bonds_valence) == i32::from(atom.valence) + 1 {
                if is_n {
                    add!(nNumTgDBNHRevrs, atom.charge == 0 && atom.num_H == 1);
                    add!(nNumTgDBNMinusRevrs, atom.charge == -1 && atom.num_H == 0);
                    add!(nNumTgDBNRevrs, atom.charge == 0 && atom.num_H == 0);
                } else {
                    add!(nNumTgDBORevrs, atom.charge == 0 && atom.num_H == 0);
                }
            }
        } else {
            add!(nNumTgNInChI, is_n);
            add!(nNumTgOInChI, is_o);
            if atom.chem_bonds_valence == atom.valence {
                if is_n {
                    add!(nNumTgNHInChI, atom.charge == 0 && atom.num_H == 1);
                    add!(nNumTgNH2InChI, atom.charge == 0 && atom.num_H == 2);
                    add!(nNumTgNMinusInChI, atom.charge == -1 && atom.num_H == 0);
                    add!(nNumTgNHMinusInChI, atom.charge == -1 && atom.num_H == 1);
                } else {
                    add!(nNumTgOHInChI, atom.charge == 0 && atom.num_H == 1);
                    add!(nNumTgOMinusInChI, atom.charge == -1 && atom.num_H == 0);
                }
            } else if i32::from(atom.chem_bonds_valence) == i32::from(atom.valence) + 1 {
                if is_n {
                    add!(nNumTgDBNHInChI, atom.charge == 0 && atom.num_H == 1);
                    add!(nNumTgDBNMinusInChI, atom.charge == -1 && atom.num_H == 0);
                    add!(nNumTgDBNInChI, atom.charge == 0 && atom.num_H == 0);
                } else {
                    add!(nNumTgDBOInChI, atom.charge == 0 && atom.num_H == 0);
                }
            }
        }
    }

    let input = [
        if pInChI[0].is_null() {
            None
        } else {
            Some(
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        },
        if pInChI[1].is_null() {
            None
        } else {
            Some(
                heap.slice(pInChI[1].as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        },
    ];
    let reversed = if pStruct.pOneINChI[0].is_null() {
        None
    } else {
        Some(
            heap.slice(pStruct.pOneINChI[0].as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let reversed_mobile_inchi = if pStruct.pOneINChI[1].is_null() {
        None
    } else {
        Some(
            heap.slice(pStruct.pOneINChI[1].as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let mobile_atoms = if pStruct.pOne_norm_data[0].is_null() {
        SourceMutPointer::null()
    } else {
        heap.slice(pStruct.pOne_norm_data[0].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at
    };
    let input_mobile_h = input[0]
        .as_ref()
        .map_or(SourceMutPointer::null(), |inchi| inchi.nNum_H);
    let reversed_mobile_h = reversed
        .as_ref()
        .map_or(SourceMutPointer::null(), |inchi| inchi.nNum_H);
    *pc2i = CMP2MHINCHI::default();
    pc2i.nNumTgInChI = pStruct.ti.num_t_groups as i16;
    pc2i.nNumTgRevrs = pStruct.One_ti.num_t_groups as i16;
    pc2i.bHasDifference |= (pc2i.nNumTgInChI != pc2i.nNumTgRevrs) as i8;
    pc2i.nNumRemHInChI = pStruct.nNumRemovedProtonsMobHInChI as i16;
    pc2i.nNumRemHRevrs = pStruct.One_ti.tni.nNumRemovedProtons as i16;
    pc2i.bFixedHLayerExistsRevrs = input[1]
        .as_ref()
        .is_some_and(|inchi| inchi.nNumberOfAtoms > 0 && inchi.bDeleted == 0)
        as i8;

    for (info, reversed_side) in [(&pStruct.ti, false), (&pStruct.One_ti, true)] {
        let mut group_index = 0_i32;
        while group_index < info.num_t_groups {
            let group = heap
                .slice(info.t_group.as_const())?
                .get(
                    usize::try_from(group_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut endpoint_offset = 0_i32;
            while endpoint_offset < i32::from(group.nNumEndpoints) {
                let endpoint_position =
                    i32::from(group.nFirstEndpointAtNoPos).wrapping_add(endpoint_offset);
                let endpoint_number = *heap
                    .slice(info.nEndpointAtomNumber.as_const())?
                    .get(
                        usize::try_from(endpoint_position)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atom_index = if reversed_side {
                    usize::from(
                        *heap
                            .slice(pStruct.nCanon2Atno[0].as_const())?
                            .get(usize::from(endpoint_number))
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    )
                } else {
                    usize::from(endpoint_number)
                };
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let is_n = valence.cNumValenceElectrons == 5 && valence.cPeriodicRowNumber == 1;
                let is_o = valence.cNumValenceElectrons == 6;
                if i32::from(is_n) + i32::from(is_o) != 1 {
                    return Ok(if reversed_side {
                        RI_ERR_PROGR
                    } else {
                        RI_ERR_SYNTAX
                    });
                }
                classify(
                    pc2i,
                    at2.get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    is_n,
                    reversed_side,
                );
                endpoint_offset = endpoint_offset.wrapping_add(1);
            }
            group_index = group_index.wrapping_add(1);
        }
    }

    let mut group_index = 0_i32;
    while group_index < pStruct.ti.num_t_groups && group_index < pStruct.One_ti.num_t_groups {
        let original = heap
            .slice(pStruct.ti.t_group.as_const())?
            .get(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let reverse = heap
            .slice(pStruct.One_ti.t_group.as_const())?
            .get(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let original_h = i32::from(original.num[0]).wrapping_sub(i32::from(original.num[1]));
        let original_minus = i32::from(original.num[1]);
        let reverse_h = i32::from(reverse.num[0]).wrapping_sub(i32::from(reverse.num[1]));
        let reverse_minus = i32::from(reverse.num[1]);
        pc2i.bHasDifference |= (original_h != reverse_h) as i8;
        pc2i.bHasDifference |= (original_minus != reverse_minus) as i8;
        if original.nNumEndpoints == reverse.nNumEndpoints {
            pc2i.nNumTgDiffH = pc2i
                .nNumTgDiffH
                .wrapping_add(i16::from(original_h != reverse_h));
            pc2i.nNumTgDiffMinus = pc2i
                .nNumTgDiffMinus
                .wrapping_add(i16::from(original_minus != reverse_minus));
        }
        pc2i.bHasDifference |= (original.nNumEndpoints != reverse.nNumEndpoints) as i8;
        pc2i.nNumTgHInChI = pc2i.nNumTgHInChI.wrapping_add(original_h as i16);
        pc2i.nNumTgMInChI = pc2i.nNumTgMInChI.wrapping_add(original_minus as i16);
        pc2i.nNumTgHRevrs = pc2i.nNumTgHRevrs.wrapping_add(reverse_h as i16);
        pc2i.nNumTgMRevrs = pc2i.nNumTgMRevrs.wrapping_add(reverse_minus as i16);
        group_index = group_index.wrapping_add(1);
    }
    while group_index < pStruct.ti.num_t_groups {
        let group = heap
            .slice(pStruct.ti.t_group.as_const())?
            .get(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let h = i32::from(group.num[0]).wrapping_sub(i32::from(group.num[1]));
        pc2i.nNumTgHInChI = pc2i.nNumTgHInChI.wrapping_add(h as i16);
        pc2i.nNumTgMInChI = pc2i.nNumTgMInChI.wrapping_add(group.num[1] as i16);
        pc2i.bHasDifference |= 1;
        group_index = group_index.wrapping_add(1);
    }
    while group_index < pStruct.One_ti.num_t_groups {
        let group = heap
            .slice(pStruct.One_ti.t_group.as_const())?
            .get(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let h = i32::from(group.num[0]).wrapping_sub(i32::from(group.num[1]));
        pc2i.nNumTgHRevrs = pc2i.nNumTgHRevrs.wrapping_add(h as i16);
        pc2i.nNumTgMRevrs = pc2i.nNumTgMRevrs.wrapping_add(group.num[1] as i16);
        pc2i.bHasDifference |= 1;
        group_index = group_index.wrapping_add(1);
    }

    let atom_count = if pStruct.num_atoms > 0 {
        pStruct.num_atoms as usize
    } else {
        0
    };
    let mut differences = 0_usize;
    for atom_index in 0..atom_count {
        let canonical = usize::from(
            *heap
                .slice(pStruct.nAtno2Canon[0].as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let endpoint_input = at2
            .get(atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint;
        let endpoint_reverse = if mobile_atoms.is_null() {
            0
        } else {
            heap.slice(mobile_atoms.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .endpoint
        };
        let mobile_h_input = if input_mobile_h.is_null() {
            0
        } else {
            *heap
                .slice(input_mobile_h.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let mobile_h_reverse = if reversed_mobile_h.is_null() {
            0
        } else {
            *heap
                .slice(reversed_mobile_h.as_const())?
                .get(canonical)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        if (endpoint_input == 0) != (endpoint_reverse == 0) || mobile_h_input != mobile_h_reverse {
            if differences >= MAX_DIFF_FIXH as usize {
                pc2i.len_c2at = differences as i16;
                return Ok(RI_ERR_PROGR);
            }
            let difference = &mut pc2i.c2at[differences];
            difference.endptInChI = endpoint_input;
            difference.endptRevrs = endpoint_reverse;
            pc2i.bHasDifference |= ((endpoint_input == 0) != (endpoint_reverse == 0)) as i8;
            difference.atomNumber = atom_index as u16;
            let valence = pVA
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            difference.nValElectr = valence.cNumValenceElectrons as u8;
            difference.nPeriodNum = valence.cPeriodicRowNumber as u8;
            difference.nMobHInChI = if let Some(mobile) =
                input[1].as_ref().filter(|inchi| !inchi.nNum_H.is_null())
            {
                *heap
                    .slice(mobile.nNum_H.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            } else if let Some(fixed) = input[0].as_ref().filter(|inchi| !inchi.nNum_H.is_null()) {
                *heap
                    .slice(fixed.nNum_H.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            } else {
                0
            };
            difference.nMobHRevrs = if let Some(mobile) = reversed_mobile_inchi
                .as_ref()
                .filter(|inchi| !inchi.nNum_H.is_null())
            {
                *heap
                    .slice(mobile.nNum_H.as_const())?
                    .get(canonical)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            } else if let Some(fixed) = reversed.as_ref().filter(|inchi| !inchi.nNum_H.is_null()) {
                *heap
                    .slice(fixed.nNum_H.as_const())?
                    .get(canonical)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            } else {
                0
            };
            pc2i.nNumDiffMobH = pc2i.nNumDiffMobH.wrapping_add(
                (mobile_h_input != mobile_h_reverse && endpoint_reverse == 0 && endpoint_input == 0)
                    as u8,
            );
            pc2i.bHasDifference |= (mobile_h_input != mobile_h_reverse) as i8;
            difference.nNumHRevrs = at2[atom_index].num_H;
            difference.nAtChargeRevrs = at2[atom_index].charge;
            differences += 1;
        }
        pc2i.nNumEndpInChI = pc2i
            .nNumEndpInChI
            .wrapping_add(i16::from(endpoint_input != 0));
        pc2i.nNumEndpRevrs = pc2i
            .nNumEndpRevrs
            .wrapping_add(i16::from(endpoint_reverse != 0));
        if pVA[atom_index].cMetal == 0 && !mobile_atoms.is_null() {
            let atom = heap
                .slice(mobile_atoms.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom.endpoint == 0 {
                pc2i.nChargeMobHRevrsNonMetal =
                    pc2i.nChargeMobHRevrsNonMetal.wrapping_add(atom.charge);
            }
        }
    }
    pc2i.nChargeMobHRevrsNonMetal = pc2i
        .nChargeMobHRevrsNonMetal
        .wrapping_add(pTCGroups.tgroup_charge as i8);
    pc2i.nChargeMobHInChI = input[0]
        .as_ref()
        .map_or(0, |inchi| inchi.nTotalCharge as i8);
    pc2i.nChargeMobHRevrs = reversed
        .as_ref()
        .map_or(0, |inchi| inchi.nTotalCharge as i8);
    pc2i.bHasDifference |= (pc2i.nChargeMobHInChI != pc2i.nChargeMobHRevrs) as i8;
    pc2i.len_c2at = differences as i16;
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CheckAndRefixStereobonds(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &mut StrFromINChI,
    at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1678 CheckAndRefixStereobonds
    /*
    int CheckAndRefixStereobonds(BN_STRUCT* pBNS, BN_DATA* pBD, StrFromINChI* pStruct,
        inp_ATOM* at, inp_ATOM* at2, VAL_AT* pVA, ALL_TC_GROUPS* pTCGroups,
        int* pnNumRunBNS, int* pnTotalDelta, int forbidden_edge_mask)
    {
        int forbidden_edge_stereo = BNS_EDGE_FORBIDDEN_MASK;
        int inv_forbidden_edge_stereo = ~forbidden_edge_stereo;

        int i, k, ne, j1, j2, num_wrong, num_fixed;
        int ret2, retBNS, ret;
        int num_at = pStruct->num_atoms;
        int num_deleted_H = pStruct->num_deleted_H;
        int len_at = num_at + num_deleted_H;
        EDGE_LIST FixedEdges, WrongEdges, CarbonChargeEdges;

        BNS_EDGE* pEdge;
        Vertex      v1, v2;
        BNS_VERTEX* pv1, * pv2;

        ret = 0;

        /* to simplify, prepare new at[] from pBNS */
        memcpy(at2, at, len_at * sizeof(at2[0]));
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
        if (ret2 < 0)
        {
            return ret;
        }

        num_wrong = 0;
        /* find wrong double bonds */
        for (i = 0; i < num_at; i++)
        {
            if (at2[i].valence == 3 &&
                at2[i].chem_bonds_valence - at2[i].valence == 1 &&
                at2[i].sb_parity[0] && at2[i].sb_parity[1] && !at2[i].sb_parity[2] &&
                (at2[i].bond_type[j1 = (int)at2[i].sb_ord[0]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE &&
                (at2[i].bond_type[j2 = (int)at2[i].sb_ord[1]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE &&
                j1 != j2)
            {

                num_wrong++;
            }
        }
        if (!num_wrong)
        {
            return 0;
        }
        num_fixed = 0;
        for (i = 0; i < pBNS->num_bonds; i++)
        {
            pEdge = pBNS->edge + i;
            if (pEdge->forbidden & forbidden_edge_stereo)
            {
                num_fixed++;
            }
        }

        /* there may be no fixed stereo bonds at all, see #87607 */
        AllocEdgeList(&CarbonChargeEdges, EDGE_LIST_CLEAR);
        AllocEdgeList(&FixedEdges, EDGE_LIST_CLEAR);
        AllocEdgeList(&WrongEdges, EDGE_LIST_CLEAR);

        /* do not goto exit_function before reaching this point: EdgeLists have not been initiated */

        if (0 > (ret = ForbidCarbonChargeEdges(pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask)))
        {
            goto exit_function;
        }
        if ((ret = AllocEdgeList(&FixedEdges, num_fixed)) ||
            (ret = AllocEdgeList(&WrongEdges, num_wrong)))
        {
            goto exit_function;
        }
        /* collect wrong double bonds and set flow=0 */
        for (i = 0; i < num_at && WrongEdges.num_edges < num_wrong; i++)
        {
            if (at2[i].valence == 3 &&
                at2[i].chem_bonds_valence - at2[i].valence == 1 &&
                at2[i].sb_parity[0] && at2[i].sb_parity[1] && !at2[i].sb_parity[2] &&
                (at2[i].bond_type[j1 = (int)at2[i].sb_ord[0]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE &&
                (at2[i].bond_type[j2 = (int)at2[i].sb_ord[1]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE &&
                j1 != j2)
            {
                switch (j1 + j2)
                {
                case 1: /* 0, 1 */
                    k = 2;
                    break;
                case 2: /* 0, 2 */
                    k = 1;
                    break;
                case 3: /* 1, 2 */
                    k = 0;
                    break;
                default:
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                ne = pBNS->vert[i].iedge[k];
                pEdge = pBNS->edge + ne;
                v1 = pEdge->neighbor1;
                v2 = pEdge->neighbor12 ^ v1;
                pv1 = pBNS->vert + v1;
                pv2 = pBNS->vert + v2;

                if (!pEdge->flow)
                {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                pEdge->flow--;
                pEdge->forbidden |= forbidden_edge_mask;
                pv1->st_edge.flow--;
                pv2->st_edge.flow--;
                pBNS->tot_st_flow -= 2;
                if ((ret = AddToEdgeList(&WrongEdges, ne, 0))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
            }
        }
        /* remove forbidden mark from stereo bonds (unfix stereo bonds) */
        for (i = 0; i < pBNS->num_bonds && FixedEdges.num_edges < num_fixed; i++)
        {
            pEdge = pBNS->edge + i;
            if (pEdge->forbidden & forbidden_edge_stereo)
            {
                pEdge->forbidden &= inv_forbidden_edge_stereo;
                FixedEdges.pnEdges[FixedEdges.num_edges++] = i;
            }
        }
        /* Run BNS to move charges and rearrange bond orders */
        retBNS = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
        (*pnNumRunBNS)++;
        if (retBNS < 0)
        {
            goto exit_function;
        }
        else
        {
            if (retBNS > 0)
            {
                *pnTotalDelta += retBNS;
            }
        }
        /* remove forbidden_edge_mask and set forbidden_edge_stereo */
        RemoveForbiddenEdgeMask(pBNS, &WrongEdges, forbidden_edge_mask);
        /* allow carbon charges to change */
        RemoveForbiddenEdgeMask(pBNS, &CarbonChargeEdges, forbidden_edge_mask);
        /* fix previously unfixed stereo bonds */
        SetForbiddenEdgeMask(pBNS, &FixedEdges, forbidden_edge_stereo);
        /* Run BNS again in case not all edge flows are maximal */
        ret2 = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
        (*pnNumRunBNS)++;
        if (ret2 < 0)
        {
            goto exit_function;
        }
        else
        {
            if (ret2 > 0)
            {
                *pnTotalDelta += retBNS;
            }
        }
        ret = retBNS;

    exit_function:

        AllocEdgeList(&CarbonChargeEdges, EDGE_LIST_FREE);
        AllocEdgeList(&FixedEdges, EDGE_LIST_FREE);
        AllocEdgeList(&WrongEdges, EDGE_LIST_FREE);

        return ret;
    }
    */
    // END INCHI C FUNCTION: CheckAndRefixStereobonds
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CheckAndRefixStereobonds
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // INCHI✔❌: BNS_EDGE_FORBIDDEN_MASK=1; EdgeIndex=int; Vertex=AT_NUMB.
    // INCHI✔❌: SourceHeap snapshots and checked pointer lookups add overhead to direct C traversal.
    // END INCHI ACTIVE MACRO CONFIGURATION: CheckAndRefixStereobonds

    let forbidden_edge_stereo = BNS_EDGE_FORBIDDEN_MASK as i32;
    let inverse_forbidden_edge_stereo = !forbidden_edge_stereo;
    let num_at = pStruct.num_atoms;
    let len_at = num_at.wrapping_add(pStruct.num_deleted_H);
    let copied_len = usize::try_from(len_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if copied_len != 0 {
        let copied = heap
            .slice(at.as_const())?
            .get(..copied_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        heap.slice_mut(at2)?
            .get_mut(..copied_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&copied);
    }
    pStruct.at = at2;
    let copy_result = CopyBnsToAtom(heap, pStruct, pBNS, pVA, pTCGroups, 1);
    pStruct.at = at;
    let copy_result = copy_result?;
    if copy_result < 0 {
        return Ok(0);
    }

    let atom_count = usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let copied_atoms = heap
        .slice(at2.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let wrong_orders = |atom: &inp_ATOM| -> Option<(usize, usize)> {
        if atom.valence != 3
            || i32::from(atom.chem_bonds_valence).wrapping_sub(i32::from(atom.valence)) != 1
            || atom.sb_parity[0] == 0
            || atom.sb_parity[1] == 0
            || atom.sb_parity[2] != 0
        {
            return None;
        }
        let first = usize::try_from(i32::from(atom.sb_ord[0])).ok()?;
        let second = usize::try_from(i32::from(atom.sb_ord[1])).ok()?;
        if atom.bond_type.get(first).copied()? & BOND_TYPE_MASK as u8 != BOND_TYPE_SINGLE as u8
            || atom.bond_type.get(second).copied()? & BOND_TYPE_MASK as u8 != BOND_TYPE_SINGLE as u8
            || first == second
        {
            return None;
        }
        Some((first, second))
    };
    let num_wrong = copied_atoms
        .iter()
        .filter(|atom| wrong_orders(atom).is_some())
        .count();
    if num_wrong == 0 {
        return Ok(0);
    }
    let num_wrong = i32::try_from(num_wrong).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let edge_count =
        usize::try_from(pBNS.num_bonds).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let num_fixed = heap
        .slice(pBNS.edge.as_const())?
        .get(..edge_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .iter()
        .filter(|edge| i32::from(edge.forbidden) & forbidden_edge_stereo != 0)
        .count();
    let num_fixed = i32::try_from(num_fixed).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;

    let mut carbon_charge_edges = EDGE_LIST::default();
    let mut fixed_edges = EDGE_LIST::default();
    let mut wrong_edges = EDGE_LIST::default();
    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut ret = ForbidCarbonChargeEdges(
            heap,
            pBNS,
            pTCGroups,
            &mut carbon_charge_edges,
            forbidden_edge_mask,
        )?;
        if ret < 0 {
            return Ok(ret);
        }
        ret = AllocEdgeList(heap, &mut fixed_edges, num_fixed)?;
        if ret == 0 {
            ret = AllocEdgeList(heap, &mut wrong_edges, num_wrong)?;
        }
        if ret != 0 {
            return Ok(ret);
        }

        let mut atom_index = 0_i32;
        while atom_index < num_at && wrong_edges.num_edges < num_wrong {
            let index =
                usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if let Some((first, second)) = wrong_orders(
                copied_atoms
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            ) {
                let third = match first.wrapping_add(second) {
                    1 => 2,
                    2 => 1,
                    3 => 0,
                    _ => return Ok(RI_ERR_PROGR),
                };
                let vertex = heap
                    .slice(pBNS.vert.as_const())?
                    .get(index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let edge_number = i32::from(
                    *heap
                        .slice(vertex.iedge.as_const())?
                        .get(third)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let edge_index = usize::try_from(edge_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let edge_before = heap
                    .slice(pBNS.edge.as_const())?
                    .get(edge_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if edge_before.flow == 0 {
                    return Ok(RI_ERR_PROGR);
                }
                let first_vertex = i32::from(edge_before.neighbor1);
                let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                {
                    let edge = heap
                        .slice_mut(pBNS.edge)?
                        .get_mut(edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.flow = edge.flow.wrapping_sub(1);
                    edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                }
                for vertex_number in [first_vertex, second_vertex] {
                    let vertex = heap
                        .slice_mut(pBNS.vert)?
                        .get_mut(
                            usize::try_from(vertex_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                }
                pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
                ret = AddToEdgeList(heap, &mut wrong_edges, edge_number, 0)?;
                if ret != 0 {
                    return Ok(ret);
                }
            }
            atom_index = atom_index.wrapping_add(1);
        }

        let mut edge_index = 0_i32;
        while edge_index < pBNS.num_bonds && fixed_edges.num_edges < num_fixed {
            let index =
                usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let edge = heap
                .slice_mut(pBNS.edge)?
                .get_mut(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(edge.forbidden) & forbidden_edge_stereo != 0 {
                edge.forbidden = (i32::from(edge.forbidden) & inverse_forbidden_edge_stereo) as i8;
                let list_index = usize::try_from(fixed_edges.num_edges)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                *heap
                    .slice_mut(fixed_edges.pnEdges)?
                    .get_mut(list_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = edge_index;
                fixed_edges.num_edges = fixed_edges.num_edges.wrapping_add(1);
            }
            edge_index = edge_index.wrapping_add(1);
        }

        let ret_bns = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
        *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
        if ret_bns < 0 {
            return Ok(ret);
        }
        if ret_bns > 0 {
            *pnTotalDelta = pnTotalDelta.wrapping_add(ret_bns);
        }
        RemoveForbiddenEdgeMask(heap, pBNS, &wrong_edges, forbidden_edge_mask)?;
        RemoveForbiddenEdgeMask(heap, pBNS, &carbon_charge_edges, forbidden_edge_mask)?;
        SetForbiddenEdgeMask(heap, pBNS, &fixed_edges, forbidden_edge_stereo)?;
        let ret2 = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
        *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
        if ret2 < 0 {
            return Ok(ret);
        }
        if ret2 > 0 {
            *pnTotalDelta = pnTotalDelta.wrapping_add(ret_bns);
        }
        Ok(ret_bns)
    })();

    let cleanup = (|| -> Result<(), SourceHeapError> {
        let _ = AllocEdgeList(heap, &mut carbon_charge_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut fixed_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut wrong_edges, EDGE_LIST_FREE)?;
        Ok(())
    })();
    let result = execution?;
    cleanup?;
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MoveChargeToRemoveCenerpoints(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &mut StrFromINChI,
    at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1860 MoveChargeToRemoveCenerpoints
    /*
    int MoveChargeToRemoveCenerpoints(BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        inp_ATOM* at,
        inp_ATOM* at2,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups,
        int* pnNumRunBNS,
        int* pnTotalDelta,
        int forbidden_edge_mask)
    {
        int i, j, neigh, num_success; /* djb-rwth: removing redundant variables */
        int num_donors, num_acceptors, bond_type, num_donors_O, num_acceptors_O, is_centerpoint_N, num_known_endpoints, num_wrong_neigh;
        int ret2, ret_forbid_edges, ret, delta;
        int num_at = pStruct->num_atoms;
        int num_deleted_H = pStruct->num_deleted_H;
        int len_at = num_at + num_deleted_H;
        int forbidden_edge_test = BNS_EDGE_FORBIDDEN_TEST;
        int bPossiblyIgnore = pStruct->charge >= 0 && (!pTCGroups->num_tgroups || (pStruct->iMobileH == TAUT_NON && pStruct->ti.num_t_groups)); /* djb-rwth: addressing LLVM warning */
        S_CHAR MobileChargeNeigh[MAXVAL], DoubleBondAcceptors[MAXVAL], DoubleBondNotONeigh[MAXVAL];
        int    numMobileChargeNeigh, numDoubleBondAcceptors, numDoubleBondNotONeigh; /* djb-rwth: removing redundant variables */
        EDGE_LIST ChargeListAllExcept_DB_O;


        BNS_EDGE* pEdgeMinus, * pe;
        Vertex      v1m, v2m;
        BNS_VERTEX* pv1m, * pv2m;
        /* djb-rwth: removing redundant code */
        num_success = 0;

        /* count O(+)H, N(+)H */

        /*
        if ( pStruct->charge >= 0 && (!pTCGroups->num_tgroups || pStruct->iMobileH == TAUT_NON && pStruct->ti.num_t_groups) ) {
            goto exit_function;
        }
        */
        if ((ret = AllocEdgeList(&ChargeListAllExcept_DB_O, EDGE_LIST_CLEAR))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }


        /* to simplify, prepare new at[] from pBNS */
        memcpy(at2, at, len_at * sizeof(at2[0]));
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
        if (ret2 < 0)
        {
            ret = ret2;
            goto exit_function;
        }
    #if ( FIND_RING_SYSTEMS == 1 )
        ret2 = MarkRingSystemsInp(at2, num_at, 0);
        if (ret2 < 0)
        {
            ret = ret2;
            goto exit_function;
        }
    #endif
        /* mark bonds that cannot be tautomeric; do not forget to remove the marks later */
        ret_forbid_edges = SetForbiddenEdges(pBNS, at2, num_at, forbidden_edge_test, 0, NULL);
        if (ret_forbid_edges < 0)
        {
            ret = ret_forbid_edges;
            goto exit_function;
        }

        for (i = 0; i < num_at; i++)
        {
            if (pVA[i].cNumValenceElectrons != 4 && /* not C, Si, Ge */
                !(pVA[i].nTautGroupEdge || (pStruct->iMobileH == TAUT_NON && pStruct->endpoint && pStruct->endpoint[i])) &&
                !at2[i].num_H && !at2[i].charge && at2[i].valence >= 2 &&
                at2[i].valence < at2[i].chem_bonds_valence &&
                is_centerpoint_elem(at2[i].el_number)) /* djb-rwth: addressing LLVM warning */
            {

                is_centerpoint_N = (pVA[i].cNumValenceElectrons == 5 && (pVA[i].cPeriodicRowNumber == 1 || pVA[i].cMetal));
                /* look at the neighbors */
                numMobileChargeNeigh = numDoubleBondAcceptors = numDoubleBondNotONeigh = num_donors = num_acceptors = 0;
                num_donors_O = num_acceptors_O = 0;
                num_known_endpoints = num_wrong_neigh = 0;
                for (j = 0; j < at2[i].valence; j++) /* djb-rwth: removing redundant code */
                {
                    neigh = at2[i].neighbor[j];
                    if ((at2[neigh].endpoint || (pStruct->iMobileH == TAUT_NON && pStruct->endpoint && pStruct->endpoint[neigh])) || at2[neigh].charge > 0) /* djb-rwth: addressing LLVM warning */
                    {
                        num_known_endpoints++;
                        continue;
                    }
                    if (pBNS->edge[pBNS->vert[i].iedge[j]].forbidden & forbidden_edge_test)
                    {
                        continue;
                    }
                    bond_type = at2[i].bond_type[j] & BOND_TYPE_MASK;
                    if (bond_type > BOND_TYPE_DOUBLE)
                    {
                        num_wrong_neigh++;
                        continue;
                    }
                    if (at2[neigh].num_H && bond_type == BOND_TYPE_SINGLE)
                    {
                        break;  /* not this case */
                    }
                    if (at2[neigh].chem_bonds_valence - at2[neigh].charge
                        != get_endpoint_valence(at2[neigh].el_number))
                    {
                        if (bond_type == BOND_TYPE_DOUBLE && pVA[neigh].cNumValenceElectrons != 6)
                        {
                            DoubleBondNotONeigh[numDoubleBondNotONeigh++] = j;
                        }
                        continue;
                    }
                    if (at2[neigh].charge == -1 && bond_type == BOND_TYPE_SINGLE &&
                        (pVA[neigh].nCMinusGroupEdge < 1 || pBNS->edge[pVA[neigh].nCMinusGroupEdge - 1].flow != 1))
                    {
                        break;
                    }
                    switch (bond_type)
                    {
                    case BOND_TYPE_SINGLE:
                        if (at2[neigh].charge != -1 || pVA[neigh].nCMinusGroupEdge <= 0)
                        {
                            num_wrong_neigh++;
                            continue;
                        }
                        num_donors++;
                        num_donors_O += (pVA[neigh].cNumValenceElectrons == 6 && pVA[neigh].cPeriodicRowNumber <= 4);
                        MobileChargeNeigh[numMobileChargeNeigh++] = j;
                        break;
                    case BOND_TYPE_DOUBLE:
                        if (at2[neigh].charge)
                        {
                            num_wrong_neigh++;
                            continue;
                        }
                        DoubleBondAcceptors[numDoubleBondAcceptors++] = j;
                        num_acceptors++;
                        num_acceptors_O += (pVA[neigh].cNumValenceElectrons == 6 && pVA[neigh].cPeriodicRowNumber <= 4);
                    }
                }
                if (j != at2[i].valence || !num_donors || !num_acceptors)
                {
                    continue;
                }
                /* special case NOn(-) */
                if (is_centerpoint_N && (num_donors == num_donors_O) && (num_acceptors == num_acceptors_O))
                {
                    continue;
                }
                if (pStruct->iMobileH == TAUT_NON && num_donors == numDoubleBondNotONeigh)
                {
                    /* fix all charges except on =O */
                    Vertex     vPathStart, vPathEnd;
                    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
                    int k, e, num_MovedCharges = 0;

                    if (!ChargeListAllExcept_DB_O.num_edges)
                    {
                        /* djb-rwth: removing redundant code */
                        for (k = 0; k < num_at; k++)
                        {
                            if (!((1 == at2[k].valence && pBNS->edge[pBNS->vert[k].iedge[0]].flow &&
                                !pBNS->edge[pBNS->vert[k].iedge[0]].forbidden &&
                                !((e = pVA[k].nCMinusGroupEdge - 1) >= 0 && pBNS->edge[e].flow) &&
                                !((e = pVA[k].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].flow) &&
                                /* 0 == at2[k].charge && */
                                pVA[k].cNumValenceElectrons == 6 && !pVA[k].cMetal &&
                                (pStruct->endpoint && pStruct->endpoint[k])) ||
                                (pStruct->fixed_H && pStruct->fixed_H[k]))) /* djb-rwth: addressing LLVM warnings */
                                /* djb-rwth: removing redundant code */
                                if ((e = pVA[k].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].flow &&
                                    !pBNS->edge[e].forbidden &&
                                    (ret = AddToEdgeList(&ChargeListAllExcept_DB_O, e, 64)))
                                {
                                    goto exit_function;
                                }
                            if ((e = pVA[k].nCPlusGroupEdge - 1) >= 0 &&
                                !pBNS->edge[e].forbidden &&
                                (ret = AddToEdgeList(&ChargeListAllExcept_DB_O, e, 64)))
                            {
                                goto exit_function;
                            }
                        }
                    }
                    /* fix double bonds to non-O neighbors connected by double bonds;
                       we will try to make these bons single */
                    for (k = 0; k < numDoubleBondNotONeigh; k++)
                    {
                        e = pBNS->vert[i].iedge[(int)DoubleBondNotONeigh[k]];
                        if (!pBNS->edge[e].forbidden &&
                            (ret = AddToEdgeList(&ChargeListAllExcept_DB_O, e, 64)))
                        {
                            goto exit_function;
                        }
                    }
                    /* attempt to make DoubleBondNotONeigh[] single */
                    SetForbiddenEdgeMask(pBNS, &ChargeListAllExcept_DB_O, forbidden_edge_mask);
                    for (k = 0; k < numDoubleBondNotONeigh && num_MovedCharges < numMobileChargeNeigh; k++)
                    {
                        pe = pBNS->edge + pBNS->vert[i].iedge[(int)DoubleBondNotONeigh[k]];
                        delta = 1;
                        if (pe->flow != delta)
                            continue;
                        pv1m = pBNS->vert + (v1m = pe->neighbor1);
                        pv2m = pBNS->vert + (v2m = pe->neighbor12 ^ v1m);
                        pv1m->st_edge.flow -= delta;
                        pv2m->st_edge.flow -= delta;
                        pe->flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;
                        ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                            &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);
                        if (ret < 0)
                        {
                            goto exit_function;
                        }
                        if (ret == 1 && ((vPathEnd == v1m && vPathStart == v2m) ||
                            (vPathEnd == v2m && vPathStart == v1m)) &&
                            nDeltaCharge == 0  /* (-) moving from one to another atom*/) /* djb-rwth: addressing LLVM warnings */
                        {
                            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                            (*pnNumRunBNS)++;
                            if (ret < 0)
                            {
                                goto exit_function;
                            }
                            else
                            {
                                if (ret == 1)
                                {
                                    *pnTotalDelta += ret;
                                    num_MovedCharges++;
                                }
                                else
                                {
                                    ret = RI_ERR_PROGR;
                                    goto exit_function;
                                }
                            }
                        }
                        else
                        {
                            /* djb-rwth: removing redundant code */
                            pv1m->st_edge.flow += delta;
                            pv2m->st_edge.flow += delta;
                            pe->flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                    }
                    RemoveForbiddenEdgeMask(pBNS, &ChargeListAllExcept_DB_O, forbidden_edge_mask);
                }
                else
                {
                    if (!bPossiblyIgnore || (!num_known_endpoints && !num_wrong_neigh && (num_acceptors_O + num_donors_O >= 3))) /* djb-rwth: addressing LLVM warning */
                    {
                        /* remove negative charges from the neighbors */
                        pBNS->vert[i].st_edge.cap += num_donors; /* enough to make all bonds to donors double */
                        pBNS->tot_st_cap += num_donors;
                        pVA[i].cInitCharge -= num_donors; /* work no matter what are known charge/valence */
                        for (j = 0; j < numMobileChargeNeigh; j++)
                        {
                            neigh = at2[i].neighbor[(int)MobileChargeNeigh[j]];
                            pEdgeMinus = pBNS->edge + ((long long)pVA[neigh].nCMinusGroupEdge - 1); /* djb-rwth: cast operator added */
                            v1m = pEdgeMinus->neighbor1;
                            v2m = pEdgeMinus->neighbor12 ^ v1m;
                            pv1m = pBNS->vert + v1m;
                            pv2m = pBNS->vert + v2m;
                            delta = pEdgeMinus->flow;
                            pv1m->st_edge.flow -= delta;
                            pv2m->st_edge.flow -= delta;
                            if (IS_BNS_VT_C_GR(pv1m->type))
                            {
                                /* irreversible change to ChargeStruct */
                                pv1m->st_edge.cap -= delta;
                            }
                            else
                            {
                                if (IS_BNS_VT_C_GR(pv2m->type))
                                {
                                    /* irreversible change to ChargeStruct */
                                    pv2m->st_edge.cap -= delta;
                                }
                                else
                                {
                                    ret = RI_ERR_PROGR;
                                    goto exit_function;
                                }
                            }
                            pBNS->tot_st_cap -= delta;
                            pBNS->tot_st_flow -= 2 * delta;
                            pEdgeMinus->flow -= delta;
                        }
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        (*pnNumRunBNS)++;
                        if (ret < 0)
                        {
                            goto exit_function;
                        }
                        else
                            if (ret == num_donors)
                            {
                                *pnTotalDelta += ret;
                                num_success++;
                                /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
                            }
                            else
                            {
                                ret = RI_ERR_PROGR;
                                goto exit_function;
                            }
                    }
                }
            }
        }
        if (ret_forbid_edges)
        {
            /* remove the marks */
            RemoveForbiddenBondFlowBits(pBNS, forbidden_edge_test);
        }
        ret = num_success;

    exit_function:

        AllocEdgeList(&ChargeListAllExcept_DB_O, EDGE_LIST_FREE);

        return ret;
    }
    */
    // END INCHI C FUNCTION: MoveChargeToRemoveCenerpoints
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MoveChargeToRemoveCenerpoints
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; FIND_RING_SYSTEMS=1.
    // INCHI✔❌: BNS_EDGE_FORBIDDEN_TEST=4; IS_BNS_VT_C_GR uses BNS_VT_C_POS_ALL/C_GROUP.
    // INCHI✔❌: SourceHeap snapshots and checked pointer lookups add overhead to direct C traversal.
    // END INCHI ACTIVE MACRO CONFIGURATION: MoveChargeToRemoveCenerpoints

    let num_at = pStruct.num_atoms;
    let len_at = num_at.wrapping_add(pStruct.num_deleted_H);
    let forbidden_edge_test = BNS_EDGE_FORBIDDEN_TEST as i32;
    let possibly_ignore = pStruct.charge >= 0
        && (pTCGroups.num_tgroups == 0
            || (i32::from(pStruct.iMobileH) == TAUT_NON as i32 && pStruct.ti.num_t_groups != 0));
    let mut charge_list = EDGE_LIST::default();
    let clear = AllocEdgeList(heap, &mut charge_list, EDGE_LIST_CLEAR)?;
    if clear != 0 {
        return Ok(clear);
    }

    let execution = (|| -> Result<i32, SourceHeapError> {
        let copied_len =
            usize::try_from(len_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if copied_len != 0 {
            let copied = heap
                .slice(at.as_const())?
                .get(..copied_len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(at2)?
                .get_mut(..copied_len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&copied);
        }
        pStruct.at = at2;
        let copy_result = CopyBnsToAtom(heap, pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct.at = at;
        let copy_result = copy_result?;
        if copy_result < 0 {
            return Ok(copy_result);
        }
        let ring_result = MarkRingSystemsInp(heap, at2, num_at, 0)?;
        if ring_result < 0 {
            return Ok(ring_result);
        }
        let ret_forbid_edges = SetForbiddenEdges(
            heap,
            pBNS,
            at2,
            num_at,
            forbidden_edge_test,
            0,
            SourceMutPointer::null(),
        )?;
        if ret_forbid_edges < 0 {
            return Ok(ret_forbid_edges);
        }

        let atom_count =
            usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atoms = heap
            .slice(at2.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        let endpoints = if pStruct.endpoint.is_null() {
            None
        } else {
            Some(
                heap.slice(pStruct.endpoint.as_const())?
                    .get(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec(),
            )
        };
        let fixed_h = if pStruct.fixed_H.is_null() {
            None
        } else {
            Some(
                heap.slice(pStruct.fixed_H.as_const())?
                    .get(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec(),
            )
        };
        let is_charge_group =
            |type_: u16| (u32::from(type_) & BNS_VT_C_POS_ALL) == BNS_VERT_TYPE_C_GROUP;
        let mut num_success = 0_i32;

        let mut atom_number = 0_i32;
        while atom_number < num_at {
            let center =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = atoms
                .get(center)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA.get(center).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let structure_endpoint = endpoints.as_ref().is_some_and(|values| values[center] != 0);
            if valence.cNumValenceElectrons == 4
                || valence.nTautGroupEdge != 0
                || (i32::from(pStruct.iMobileH) == TAUT_NON as i32 && structure_endpoint)
                || atom.num_H != 0
                || atom.charge != 0
                || atom.valence < 2
                || atom.valence >= atom.chem_bonds_valence
                || is_centerpoint_elem(atom.el_number) == 0
            {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }

            let centerpoint_n = valence.cNumValenceElectrons == 5
                && (valence.cPeriodicRowNumber == 1 || valence.cMetal != 0);
            let mut mobile_charge_neighbors = [0_i8; MAXVAL as usize];
            let mut double_bond_acceptors = [0_i8; MAXVAL as usize];
            let mut double_bond_not_o = [0_i8; MAXVAL as usize];
            let mut num_mobile = 0_usize;
            let mut num_double_acceptors = 0_usize;
            let mut num_double_not_o = 0_usize;
            let mut num_donors = 0_i32;
            let mut num_acceptors = 0_i32;
            let mut num_donors_o = 0_i32;
            let mut num_acceptors_o = 0_i32;
            let mut num_known_endpoints = 0_i32;
            let mut num_wrong_neighbors = 0_i32;
            let mut order = 0_i32;
            while order < i32::from(atom.valence) {
                let order_index =
                    usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let neighbor = usize::from(atom.neighbor[order_index]);
                let neighbor_atom = atoms
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let neighbor_endpoint = neighbor_atom.endpoint != 0
                    || (i32::from(pStruct.iMobileH) == TAUT_NON as i32
                        && endpoints
                            .as_ref()
                            .is_some_and(|values| values[neighbor] != 0));
                if neighbor_endpoint || neighbor_atom.charge > 0 {
                    num_known_endpoints = num_known_endpoints.wrapping_add(1);
                    order = order.wrapping_add(1);
                    continue;
                }
                let vertex = heap
                    .slice(pBNS.vert.as_const())?
                    .get(center)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let edge_number = i32::from(
                    *heap
                        .slice(vertex.iedge.as_const())?
                        .get(order_index)
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
                if i32::from(edge.forbidden) & forbidden_edge_test != 0 {
                    order = order.wrapping_add(1);
                    continue;
                }
                let bond_type = i32::from(atom.bond_type[order_index]) & BOND_TYPE_MASK as i32;
                if bond_type > BOND_TYPE_DOUBLE as i32 {
                    num_wrong_neighbors = num_wrong_neighbors.wrapping_add(1);
                    order = order.wrapping_add(1);
                    continue;
                }
                if neighbor_atom.num_H != 0 && bond_type == BOND_TYPE_SINGLE as i32 {
                    break;
                }
                if i32::from(neighbor_atom.chem_bonds_valence)
                    .wrapping_sub(i32::from(neighbor_atom.charge))
                    != get_endpoint_valence(neighbor_atom.el_number)
                {
                    if bond_type == BOND_TYPE_DOUBLE as i32
                        && pVA[neighbor].cNumValenceElectrons != 6
                    {
                        double_bond_not_o[num_double_not_o] = order as i8;
                        num_double_not_o += 1;
                    }
                    order = order.wrapping_add(1);
                    continue;
                }
                if neighbor_atom.charge == -1 && bond_type == BOND_TYPE_SINGLE as i32 {
                    let minus = pVA[neighbor].nCMinusGroupEdge.wrapping_sub(1);
                    if minus < 0
                        || heap.slice(pBNS.edge.as_const())?[usize::try_from(minus)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                        .flow
                            != 1
                    {
                        break;
                    }
                }
                match bond_type {
                    value if value == BOND_TYPE_SINGLE as i32 => {
                        if neighbor_atom.charge != -1 || pVA[neighbor].nCMinusGroupEdge <= 0 {
                            num_wrong_neighbors = num_wrong_neighbors.wrapping_add(1);
                            order = order.wrapping_add(1);
                            continue;
                        }
                        num_donors = num_donors.wrapping_add(1);
                        num_donors_o = num_donors_o.wrapping_add(i32::from(
                            pVA[neighbor].cNumValenceElectrons == 6
                                && pVA[neighbor].cPeriodicRowNumber <= 4,
                        ));
                        mobile_charge_neighbors[num_mobile] = order as i8;
                        num_mobile += 1;
                    }
                    value if value == BOND_TYPE_DOUBLE as i32 => {
                        if neighbor_atom.charge != 0 {
                            num_wrong_neighbors = num_wrong_neighbors.wrapping_add(1);
                            order = order.wrapping_add(1);
                            continue;
                        }
                        double_bond_acceptors[num_double_acceptors] = order as i8;
                        num_double_acceptors += 1;
                        num_acceptors = num_acceptors.wrapping_add(1);
                        num_acceptors_o = num_acceptors_o.wrapping_add(i32::from(
                            pVA[neighbor].cNumValenceElectrons == 6
                                && pVA[neighbor].cPeriodicRowNumber <= 4,
                        ));
                    }
                    _ => {}
                }
                order = order.wrapping_add(1);
            }
            let _ = double_bond_acceptors;
            if order != i32::from(atom.valence) || num_donors == 0 || num_acceptors == 0 {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }
            if centerpoint_n && num_donors == num_donors_o && num_acceptors == num_acceptors_o {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }

            if i32::from(pStruct.iMobileH) == TAUT_NON as i32
                && num_donors == num_double_not_o as i32
            {
                if charge_list.num_edges == 0 {
                    let mut candidate = 0_usize;
                    while candidate < atom_count {
                        let candidate_atom = &atoms[candidate];
                        let mut keep_unfixed = false;
                        if candidate_atom.valence == 1 {
                            let vertex = heap.slice(pBNS.vert.as_const())?[candidate].clone();
                            let bond_edge = i32::from(heap.slice(vertex.iedge.as_const())?[0]);
                            let bond =
                                &heap.slice(pBNS.edge.as_const())?[usize::try_from(bond_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                            let minus = pVA[candidate].nCMinusGroupEdge.wrapping_sub(1);
                            let plus = pVA[candidate].nCPlusGroupEdge.wrapping_sub(1);
                            keep_unfixed = bond.flow != 0
                                && bond.forbidden == 0
                                && !(minus >= 0
                                    && heap.slice(pBNS.edge.as_const())?[usize::try_from(minus)
                                        .map_err(|_| {
                                        SourceHeapError::PointerOutOfBounds
                                    })?]
                                    .flow
                                        != 0)
                                && !(plus >= 0
                                    && heap.slice(pBNS.edge.as_const())?[usize::try_from(plus)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                                    .flow
                                        == 0)
                                && pVA[candidate].cNumValenceElectrons == 6
                                && pVA[candidate].cMetal == 0
                                && endpoints
                                    .as_ref()
                                    .is_some_and(|values| values[candidate] != 0);
                        }
                        let excluded = keep_unfixed
                            || fixed_h
                                .as_ref()
                                .is_some_and(|values| values[candidate] != 0);
                        if !excluded {
                            let minus = pVA[candidate].nCMinusGroupEdge.wrapping_sub(1);
                            if minus >= 0 {
                                let edge =
                                    &heap.slice(pBNS.edge.as_const())?[usize::try_from(minus)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                                if edge.flow == 0 && edge.forbidden == 0 {
                                    let ret = AddToEdgeList(heap, &mut charge_list, minus, 64)?;
                                    if ret != 0 {
                                        return Ok(ret);
                                    }
                                }
                            }
                        }
                        let plus = pVA[candidate].nCPlusGroupEdge.wrapping_sub(1);
                        if plus >= 0
                            && heap.slice(pBNS.edge.as_const())?[usize::try_from(plus)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                            .forbidden
                                == 0
                        {
                            let ret = AddToEdgeList(heap, &mut charge_list, plus, 64)?;
                            if ret != 0 {
                                return Ok(ret);
                            }
                        }
                        candidate += 1;
                    }
                }
                let center_vertex = heap.slice(pBNS.vert.as_const())?[center].clone();
                for &stored_order in double_bond_not_o.iter().take(num_double_not_o) {
                    let edge_number = i32::from(
                        heap.slice(center_vertex.iedge.as_const())?[usize::try_from(i32::from(
                            stored_order,
                        ))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                    );
                    if heap.slice(pBNS.edge.as_const())?[usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                    .forbidden
                        == 0
                    {
                        let ret = AddToEdgeList(heap, &mut charge_list, edge_number, 64)?;
                        if ret != 0 {
                            return Ok(ret);
                        }
                    }
                }
                SetForbiddenEdgeMask(heap, pBNS, &charge_list, forbidden_edge_mask)?;
                let mut moved = 0_i32;
                for &stored_order in double_bond_not_o.iter().take(num_double_not_o) {
                    if moved >= num_mobile as i32 {
                        break;
                    }
                    let edge_number = i32::from(
                        heap.slice(center_vertex.iedge.as_const())?[usize::try_from(i32::from(
                            stored_order,
                        ))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                    );
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap.slice(pBNS.edge.as_const())?[edge_index].clone();
                    if edge_before.flow != 1 {
                        continue;
                    }
                    let first = i32::from(edge_before.neighbor1);
                    let second = i32::from(edge_before.neighbor12) ^ first;
                    for vertex_number in [first, second] {
                        let vertex =
                            &mut heap.slice_mut(pBNS.vert)?[usize::try_from(vertex_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                        vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                    }
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = heap
                        .slice(pBNS.edge.as_const())?[edge_index]
                        .flow
                        .wrapping_sub(1);
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
                    let mut path_start = 0;
                    let mut path_end = 0;
                    let mut path_len = 0;
                    let mut delta_h = 0;
                    let mut delta_charge = 0;
                    let mut visited = 0;
                    let ret = RunBnsTestOnce(
                        heap,
                        pBNS,
                        pBD,
                        pVA,
                        &mut path_start,
                        &mut path_end,
                        &mut path_len,
                        &mut delta_h,
                        &mut delta_charge,
                        &mut visited,
                    )?;
                    if ret < 0 {
                        return Ok(ret);
                    }
                    if ret == 1
                        && ((path_end == first && path_start == second)
                            || (path_end == second && path_start == first))
                        && delta_charge == 0
                    {
                        let ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
                        if ret < 0 {
                            return Ok(ret);
                        }
                        if ret == 1 {
                            *pnTotalDelta = pnTotalDelta.wrapping_add(ret);
                            moved = moved.wrapping_add(1);
                        } else {
                            return Ok(RI_ERR_PROGR);
                        }
                    } else {
                        for vertex_number in [first, second] {
                            let vertex =
                                &mut heap.slice_mut(pBNS.vert)?[usize::try_from(vertex_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(1);
                        }
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = heap
                            .slice(pBNS.edge.as_const())?[edge_index]
                            .flow
                            .wrapping_add(1);
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &charge_list, forbidden_edge_mask)?;
            } else if !possibly_ignore
                || (num_known_endpoints == 0
                    && num_wrong_neighbors == 0
                    && num_acceptors_o.wrapping_add(num_donors_o) >= 3)
            {
                heap.slice_mut(pBNS.vert)?[center].st_edge.cap = heap
                    .slice(pBNS.vert.as_const())?[center]
                    .st_edge
                    .cap
                    .wrapping_add(num_donors);
                pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_add(num_donors);
                pVA[center].cInitCharge =
                    i32::from(pVA[center].cInitCharge).wrapping_sub(num_donors) as i8;
                for &stored_order in mobile_charge_neighbors.iter().take(num_mobile) {
                    let neighbor = usize::from(
                        atom.neighbor[usize::try_from(i32::from(stored_order))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                    );
                    let minus = pVA[neighbor].nCMinusGroupEdge.wrapping_sub(1);
                    let edge_index =
                        usize::try_from(minus).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap.slice(pBNS.edge.as_const())?[edge_index].clone();
                    let first = i32::from(edge_before.neighbor1);
                    let second = i32::from(edge_before.neighbor12) ^ first;
                    let delta = edge_before.flow;
                    for vertex_number in [first, second] {
                        let vertex =
                            &mut heap.slice_mut(pBNS.vert)?[usize::try_from(vertex_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                        vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                    }
                    let first_type = heap.slice(pBNS.vert.as_const())?[usize::try_from(first)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                    .type_;
                    let second_type = heap.slice(pBNS.vert.as_const())?[usize::try_from(second)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                    .type_;
                    if is_charge_group(first_type) {
                        let vertex = &mut heap.slice_mut(pBNS.vert)?[usize::try_from(first)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                        vertex.st_edge.cap = vertex.st_edge.cap.wrapping_sub(delta);
                    } else if is_charge_group(second_type) {
                        let vertex = &mut heap.slice_mut(pBNS.vert)?[usize::try_from(second)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                        vertex.st_edge.cap = vertex.st_edge.cap.wrapping_sub(delta);
                    } else {
                        return Ok(RI_ERR_PROGR);
                    }
                    pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_sub(delta);
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(delta.wrapping_mul(2));
                    heap.slice_mut(pBNS.edge)?[edge_index].flow =
                        edge_before.flow.wrapping_sub(delta);
                }
                let ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
                if ret < 0 {
                    return Ok(ret);
                }
                if ret == num_donors {
                    *pnTotalDelta = pnTotalDelta.wrapping_add(ret);
                    num_success = num_success.wrapping_add(1);
                } else {
                    return Ok(RI_ERR_PROGR);
                }
            }
            atom_number = atom_number.wrapping_add(1);
        }
        if ret_forbid_edges != 0 {
            RemoveForbiddenBondFlowBits(heap, pBNS, forbidden_edge_test)?;
        }
        Ok(num_success)
    })();

    let cleanup = AllocEdgeList(heap, &mut charge_list, EDGE_LIST_FREE);
    let result = execution?;
    let _ = cleanup?;
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeSingleBondsMetal2ChargedHeteroat(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &mut StrFromINChI,
    at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2194 MakeSingleBondsMetal2ChargedHeteroat
    /*
    int MakeSingleBondsMetal2ChargedHeteroat(BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        inp_ATOM* at,
        inp_ATOM* at2,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups,
        int* pnNumRunBNS,
        int* pnTotalDelta,
        int forbidden_edge_mask)
    {
        int i;

        int ret2, ret, pass;
        int num_at = pStruct->num_atoms;
        int num_deleted_H = pStruct->num_deleted_H;
        int len_at = num_at + num_deleted_H;
        int inv_forbidden_edge_mask = ~forbidden_edge_mask;

        int         j, k;
        int        cur_num_edges;
        BNS_EDGE* e;
        Vertex     v1, v2;

        EdgeIndex* pFixedEdges;
        int        nNumEdgesToFix;

        ret = 0;

        /* to simplify, prepare new at[] from pBNS */
        memcpy(at2, at, len_at * sizeof(at2[0]));
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
        if (ret2 < 0)
        {
            ret = ret2;
            goto exit_function;
        }

        pFixedEdges = NULL;

        nNumEdgesToFix = 0; /* cpunt nNumEdgesToFix only when pass==0 */
        cur_num_edges = 0; /* count cur_num_edges  only when pass==1; at the end they must be equal */
        for (pass = 0; pass < 2; pass++)
        {
            if (pass)
            {
                /* 2nd pass: allocate edge storage */
                if (!nNumEdgesToFix)
                {
                    break; /* nothing to do */
                }
                pFixedEdges = (EdgeIndex*)inchi_malloc(nNumEdgesToFix * sizeof(pFixedEdges[0]));
                if (!pFixedEdges)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }
            for (i = 0; i < num_at; i++)
            {
                int neigh;
                if (pVA[i].cMetal)
                {
                    for (j = 0; j < at2[i].valence; j++)
                    {
                        neigh = at2[i].neighbor[j];
                        if (pVA[neigh].cNumValenceElectrons == 4 &&
                            pVA[neigh].cPeriodicRowNumber == 1)
                        {
                            continue; /* ignore carbon */
                        }
                        if (at2[i].bond_type[j] > BOND_TYPE_SINGLE && at2[neigh].charge &&
                            !pVA[neigh].cMetal && pVA[neigh].cnListIndex > 0)
                        {
                            int cnBits = at2[neigh].charge > 0 ? MAKE_CN_BITS(cn_bits_N, cn_bits_P, 0, 0) :
                                MAKE_CN_BITS(cn_bits_N, cn_bits_M, 0, 0);
                            int atBits = cnList[pVA[neigh].cnListIndex - 1].bits;
                            for (k = 0; k < MAX_NUM_CN_BITS - 1; k++, atBits >>= cn_bits_shift)
                            {
                                /* ??? */
                                if ((atBits & cnBits) == cnBits)
                                {
                                    break;
                                }
                            }
                            if (k == MAX_NUM_CN_BITS - 1)
                            {
                                continue;
                            }
                            if (pass == 0)
                            {
                                nNumEdgesToFix++;
                            }
                            else
                            {
                                pFixedEdges[cur_num_edges++] = pBNS->vert[i].iedge[j];
                            }
                        }
                    }
                }
            }
        }

        /* restore the initial structures */
        memcpy(at2, at, ((long long)num_at + (long long)num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operators added */

        if (nNumEdgesToFix && pFixedEdges)
        {
            if (nNumEdgesToFix != cur_num_edges)
            {
                ret = RI_ERR_PROGR;
                goto pre_exit_function; /* djb-rwth: fixing coverity ID #499637 */
            }
            /* change edge flow, fix the edges, and run BNS */
            for (i = 0; i < nNumEdgesToFix; i++)
            {
                e = pBNS->edge + pFixedEdges[i];
                v1 = e->neighbor1;
                v2 = e->neighbor12 ^ v1;
                e->flow--;
                e->forbidden |= forbidden_edge_mask;
                pBNS->vert[v1].st_edge.flow--;
                pBNS->vert[v2].st_edge.flow--;
                pBNS->tot_st_flow -= 2;
                (*pnTotalDelta) -= 2;
            }
            /* Run BNS allowing to change any charges */
            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
            (*pnNumRunBNS)++;
            if (ret < 0)
            {
                goto pre_exit_function; /* djb-rwth: fixing coverity ID #499637 */
            }
            else
            {
                (*pnTotalDelta) += ret;
            }
            /* unfix the edges */
            for (i = 0; i < nNumEdgesToFix; i++)
            {
                e = pBNS->edge + pFixedEdges[i];
                e->forbidden &= inv_forbidden_edge_mask;
            }
            if (ret < 2 * nNumEdgesToFix)
            {
                /* not all fixes succeeded */
                ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                (*pnNumRunBNS)++;
                if (ret < 0)
                {
                    goto pre_exit_function; /* djb-rwth: fixing coverity ID #499637 */
                }
                else
                {
                    (*pnTotalDelta) += ret;
                }
            }
        }

    pre_exit_function:
        if (pFixedEdges)
        {
            inchi_free(pFixedEdges);
            pFixedEdges = NULL;
        }

    exit_function:
        return ret;
    }
    */
    // END INCHI C FUNCTION: MakeSingleBondsMetal2ChargedHeteroat
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MakeSingleBondsMetal2ChargedHeteroat
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // INCHI✔❌: EdgeIndex=int; MAX_NUM_CN_BITS=4; MAKE_CN_BITS uses cn_bits_shift=3.
    // INCHI✔❌: SourceHeap checked pointer lookups and temporary edge vector add overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: MakeSingleBondsMetal2ChargedHeteroat

    let num_at = pStruct.num_atoms;
    let len_at = num_at.wrapping_add(pStruct.num_deleted_H);
    let copied_len = usize::try_from(len_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if copied_len != 0 {
        let copied = heap
            .slice(at.as_const())?
            .get(..copied_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        heap.slice_mut(at2)?
            .get_mut(..copied_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&copied);
    }
    pStruct.at = at2;
    let copy_result = CopyBnsToAtom(heap, pStruct, pBNS, pVA, pTCGroups, 1);
    pStruct.at = at;
    let copy_result = copy_result?;
    if copy_result < 0 {
        return Ok(copy_result);
    }

    let atom_count = usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let copied_atoms = heap
        .slice(at2.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let scan_edges = |heap: &SourceHeap| -> Result<Vec<i32>, SourceHeapError> {
        let mut found = Vec::new();
        for (atom_index, atom) in copied_atoms.iter().enumerate() {
            if pVA
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .cMetal
                == 0
            {
                continue;
            }
            let vertex = heap
                .slice(pBNS.vert.as_const())?
                .get(atom_index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let valence = usize::try_from(i32::from(atom.valence))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            for order in 0..valence {
                let neighbor = usize::from(atom.neighbor[order]);
                let neighbor_valence = pVA
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if neighbor_valence.cNumValenceElectrons == 4
                    && neighbor_valence.cPeriodicRowNumber == 1
                {
                    continue;
                }
                let neighbor_atom = copied_atoms
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if atom.bond_type[order] <= BOND_TYPE_SINGLE as u8
                    || neighbor_atom.charge == 0
                    || neighbor_valence.cMetal != 0
                    || neighbor_valence.cnListIndex <= 0
                {
                    continue;
                }
                let charge_bit = if neighbor_atom.charge > 0 {
                    cn_bits_P
                } else {
                    cn_bits_M
                };
                let required = cn_bits_N as i32 | ((charge_bit as i32) << cn_bits_shift);
                let mut bits = CN_LIST
                    .get(
                        usize::try_from(neighbor_valence.cnListIndex.wrapping_sub(1))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .bits;
                let mut slot = 0_i32;
                while slot < MAX_NUM_CN_BITS as i32 - 1 {
                    if bits & required == required {
                        break;
                    }
                    slot = slot.wrapping_add(1);
                    bits >>= cn_bits_shift;
                }
                if slot == MAX_NUM_CN_BITS as i32 - 1 {
                    continue;
                }
                found.push(i32::from(
                    *heap
                        .slice(vertex.iedge.as_const())?
                        .get(order)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                ));
            }
        }
        Ok(found)
    };

    let counted_edges = scan_edges(heap)?;
    let edge_count =
        i32::try_from(counted_edges.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut fixed_owner = SourceMutPointer::<i8>::null();
    if edge_count != 0 {
        let bytes = u64::try_from(edge_count)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
            .checked_mul(std::mem::size_of::<i32>() as u64)
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        fixed_owner = match inchi_malloc(heap, bytes) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => return Ok(RI_ERR_ALLOC),
            Err(error) => return Err(error),
        };
    }
    let fixed_edges = if edge_count == 0 {
        Vec::new()
    } else {
        scan_edges(heap)?
    };

    if copied_len != 0 {
        let original = heap
            .slice(at.as_const())?
            .get(..copied_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        heap.slice_mut(at2)?
            .get_mut(..copied_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&original);
    }

    let execution = (|| -> Result<i32, SourceHeapError> {
        if edge_count
            != i32::try_from(fixed_edges.len())
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        {
            return Ok(RI_ERR_PROGR);
        }
        if edge_count == 0 {
            return Ok(0);
        }
        for &edge_number in &fixed_edges {
            let edge_index =
                usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let edge_before = heap
                .slice(pBNS.edge.as_const())?
                .get(edge_index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let first = i32::from(edge_before.neighbor1);
            let second = i32::from(edge_before.neighbor12) ^ first;
            {
                let edge = &mut heap.slice_mut(pBNS.edge)?[edge_index];
                edge.flow = edge.flow.wrapping_sub(1);
                edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
            }
            for vertex_number in [first, second] {
                let vertex = &mut heap.slice_mut(pBNS.vert)?[usize::try_from(vertex_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
            }
            pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
            *pnTotalDelta = pnTotalDelta.wrapping_sub(2);
        }
        let mut ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
        *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
        if ret < 0 {
            return Ok(ret);
        }
        *pnTotalDelta = pnTotalDelta.wrapping_add(ret);
        let inverse_mask = !forbidden_edge_mask;
        for &edge_number in &fixed_edges {
            let edge = &mut heap.slice_mut(pBNS.edge)?
                [usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?];
            edge.forbidden = (i32::from(edge.forbidden) & inverse_mask) as i8;
        }
        if ret < edge_count.wrapping_mul(2) {
            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
            *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
            if ret < 0 {
                return Ok(ret);
            }
            *pnTotalDelta = pnTotalDelta.wrapping_add(ret);
        }
        Ok(ret)
    })();
    let cleanup = if fixed_owner.is_null() {
        Ok(())
    } else {
        inchi_free(heap, fixed_owner)
    };
    let result = execution?;
    cleanup?;
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn SaltBondsToCoordBonds(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &mut StrFromINChI,
    at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2373 SaltBondsToCoordBonds
    /*
    int SaltBondsToCoordBonds(BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        inp_ATOM* at,
        inp_ATOM* at2,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups,
        int* pnNumRunBNS,
        int* pnTotalDelta,
        int forbidden_edge_mask)
    {
        int i;

        int ret2, ret, cur_success;
        int num_at = pStruct->num_atoms;
        int num_edges = pBNS->num_bonds + 2 * pBNS->num_atoms;
        int num_deleted_H = pStruct->num_deleted_H;
        int len_at = num_at + num_deleted_H;
        int inv_forbidden_edge_mask = ~forbidden_edge_mask;
        EDGE_LIST AllChargeEdges;

        int         j, k, n;
        BNS_EDGE* pe, * pePlusMetal, * peMinusO;
        BNS_VERTEX* pv1, * pv2, * pvO, * pvM;
        Vertex     v1, v2, vPlusMinus;

        EdgeIndex  ie, iePlusMetal, ieMinusO;

        Vertex     vPathStart, vPathEnd;
        int        delta, nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

        ret = 0;
        cur_success = 0;
        AllocEdgeList(&AllChargeEdges, EDGE_LIST_CLEAR);

        if (pStruct->iInchiRec == INCHI_BAS || !pStruct->pSrm->bMetalAddFlower || pStruct->pSrm->nMetalMinBondOrder)
        {
            goto exit_function;
        }

        /* to simplify, prepare new at[] from pBNS */
        memcpy(at2, at, len_at * sizeof(at2[0]));
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
        if (ret2 < 0)
        {
            ret = ret2;
            goto exit_function;
        }
        for (i = 0; i < num_at; i++)
        {
            if (bIsMetalSalt(at2, i))
            {
                if (!AllChargeEdges.num_edges)
                {
                    /*--------- one-time action: fix all bonds, charges, taut. group edges ------------*/
                    for (j = 0; j < num_at; j++)
                    {
                        /* all bonds */
                        for (k = 0; k < at2[j].valence; k++)
                        {
                            n = at2[j].neighbor[k];
                            if (n < j && !pBNS->edge[ie = pBNS->vert[j].iedge[k]].forbidden &&
                                (ret = AddToEdgeList(&AllChargeEdges, ie, num_edges)))
                            {
                                goto exit_function;
                            }
                        }
                        /* charge edges */
                        if ((ie = pVA[j].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[ie].forbidden &&
                            (ret = AddToEdgeList(&AllChargeEdges, ie, num_edges)))
                        {
                            goto exit_function;
                        }
                        if ((ie = pVA[j].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[ie].forbidden &&
                            (ret = AddToEdgeList(&AllChargeEdges, ie, num_edges)))
                        {
                            goto exit_function;
                        }
                    }
                    /* taut group edges */
                    for (j = 0; j < pTCGroups->num_tgroups; j++)
                    {
                        pv1 = pBNS->vert + (v1 = pTCGroups->pTCG[j].nVertexNumber); /* t-group vertex */ /* djb-rwth: ignoring LLVM warning: see comment below */
                        for (k = 0; k < pv1->num_adj_edges; k++)
                        {
                            /* ie, pe - tautomeric atom edge; pv2 - endpoint vertex */
                            /* Note: pe, pv2, v1 are not used here; they are to show how to traverse t-group */
                            pv2 = pBNS->vert + (pe = pBNS->edge + (ie = pv1->iedge[k]))->neighbor1; /* djb-rwth: ignoring LLVM warning: see comment above */
                            if ((ret = AddToEdgeList(&AllChargeEdges, ie, num_edges))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                    }
                    /*---------------------------------------------------------------*/
                }
                /* replace all single bonds to neutral neighbors with zero-order bonds
                   allow neighbor charge change to (-1) and metal atom charge increment +1 */
                for (k = 0; k < at2[i].valence; k++)
                {
                    n = at2[i].neighbor[k];
                    pe = pBNS->edge + pBNS->vert[i].iedge[k];
                    if (at2[n].charge || at2[i].bond_type[k] != BOND_TYPE_SINGLE)
                    {
                        continue;
                    }
                    iePlusMetal = pVA[i].nCPlusGroupEdge - 1;
                    ieMinusO = pVA[n].nCMinusGroupEdge - 1;

                    if (pe->flow != 1 || pe->forbidden || iePlusMetal < 0)
                    {
                        continue;
                    }
                    pePlusMetal = pBNS->edge + iePlusMetal;
                    if (pePlusMetal->flow <= 0)
                    {
                        continue; /* to add (+) to metal this flow must be decremented */
                    }
                    if (ieMinusO >= 0)
                    {
                        /* usually does not happen */
                        peMinusO = pBNS->edge + ieMinusO;

                        if (peMinusO->flow || pePlusMetal->forbidden || peMinusO->forbidden)
                        {
                            continue;
                        }

                        /* decrement bond order to 0 */
                        delta = 1;
                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                        pe->flow -= delta;
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;

                        SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                        pePlusMetal->forbidden &= inv_forbidden_edge_mask;
                        peMinusO->forbidden &= inv_forbidden_edge_mask;

                        ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                            &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                        if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1)) /*&& nDeltaCharge > 0*/) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* (+)charge was just moved, no change in number of charges */
                            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                            if (ret > 0)
                            {
                                (*pnNumRunBNS)++;
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
                        RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                    }
                    else
                    {
                        if (NO_VERTEX != (vPlusMinus = GetPlusMinusVertex(pBNS, pTCGroups, 1, 1)))
                        {
                            /* manually add (-) charge to O and (+) charge to metal */
                            /* decrement bond order to 0 */
                            /*---------------------------------------------------------------------------*/
                            /*                                                                           */
                            /*      (+/-)*               (+/-)           Result:                         */
                            /*        |                    ||                                            */
                            /*        |                    ||            - Added (+) to M                */
                            /*       (+)super             (+)super       - Incremented bond M-O          */
                            /*        ||                   |                                             */
                            /*        ||          =>       |             To make this attachment H,      */
                            /*       (Y)                  (Y)            increment                       */
                            /*        |                    ||            pTCGroups->pTCG[itg].tg_num_H   */
                            /*        |                    ||                                            */
                            /*       (+)metal             (+)hetero      Technical details:              */
                            /*         \\                   \            increase capacities of          */
                            /*           M                    M(+)       edges to (+/-) otherwise        */
                            /*           |                    ||         flow may not be able to         */
                            /*          -O*                -O-O          increase                        */
                            /*                                                                           */
                            /*   After that change M=O bond order from 2 to 0                            */
                            /*---------------------------------------------------------------------------*/
                            int i1, j1, k1;
                            delta = 1;
                            pvO = pBNS->vert + n;
                            pvM = pBNS->vert + i;
                            /* Increment st_edge.cap on (+/-) vertex */
                            pBNS->vert[vPlusMinus].st_edge.cap += delta;
                            /* Increment st_edge.cap on O */
                            pvO->st_edge.cap += delta;
                            /* increment cap on M-O edge */
                            pe->cap += delta;
                            /* total cap count */
                            pBNS->tot_st_cap += 2 * delta;

                            v1 = vPlusMinus;
                            v2 = n; /* atom O */

                            /* increase capacities of edges to Y  */
                            for (i1 = 0; i1 < pBNS->vert[vPlusMinus].num_adj_edges; i1++)
                            {
                                j1 = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i1]].neighbor12 ^ vPlusMinus;
                                for (k1 = 0; k1 < pBNS->vert[j1].num_adj_edges; k1++)
                                {
                                    pBNS->edge[pBNS->vert[j1].iedge[k1]].cap += delta;
                                }
                            }
                            SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                            pePlusMetal->forbidden &= inv_forbidden_edge_mask;
                            pe->forbidden &= inv_forbidden_edge_mask;

                            ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);
                            cur_success = 0;
                            if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                                (vPathEnd == v2 && vPathStart == v1)) /*&& nDeltaCharge == 1*/) /* djb-rwth: addressing LLVM warnings */
                            {
                                /* Added (+)charge to -N< => nDeltaCharge == 1 */
                                /* Flow change on pe (-)charge edge (atom B-O(-)) is not known to RunBnsTestOnce()) */
                                ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                                if (ret > 0)
                                {
                                    (*pnNumRunBNS)++;
                                    cur_success++; /* 01 */
                                }
                            }
                            if (cur_success)
                            {
                                /* set bond M=O order = 0 */
                                if (pe->flow != 2 * delta)
                                {
                                    ret = RI_ERR_PROGR;
                                    goto exit_function;
                                }
                                /* reduce pe bond order by 2*delta */
                                pe->flow -= 2 * delta;
                                pvO->st_edge.cap -= 2 * delta;
                                pvO->st_edge.flow -= 2 * delta;
                                pvM->st_edge.flow -= 2 * delta;
                                pvM->st_edge.cap -= 2 * delta;
                                pBNS->tot_st_cap -= 3 * delta;
                                pBNS->tot_st_flow -= 4 * delta;
                                /* fix M-O bond order to zero */
                                pe->cap -= 2 * delta;
                                /* add fixed (-) charge to O */
                                pVA[n].cInitCharge -= delta;
                            }
                            else
                            {
                                /* failed */
                                pBNS->vert[vPlusMinus].st_edge.cap -= delta;
                                pvO->st_edge.cap -= delta;
                                /*pTCGroups->pTCG[itg].edges_cap     -= delta;*/ /* ???bug??? - commented out 2006-03-22 */
                                pBNS->tot_st_cap -= 2 * delta;
                                /* decrease capacities of edges to Y  */
                                for (i1 = 0; i1 < pBNS->vert[vPlusMinus].num_adj_edges; i1++)
                                {
                                    j1 = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i1]].neighbor12 ^ vPlusMinus;
                                    for (k1 = 0; k1 < pBNS->vert[j1].num_adj_edges; k1++)
                                    {
                                        pBNS->edge[pBNS->vert[j1].iedge[k1]].cap -= delta;
                                    }
                                }
                            }
                            RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                        }
                    }
                }
            }
        }

    exit_function:

        AllocEdgeList(&AllChargeEdges, EDGE_LIST_FREE);

        return ret;
    }
        */
    // END INCHI C FUNCTION: SaltBondsToCoordBonds
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: SaltBondsToCoordBonds
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; KEEP_METAL_EDGE_FLOW=0.
    // INCHI✔❌: EdgeIndex=Vertex=int; forbidden is signed char; BOND_TYPE_SINGLE=1.
    // INCHI✔❌: SourceHeap checked pointer lookups and edge-list ownership add overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: SaltBondsToCoordBonds

    let mut all_charge_edges = EDGE_LIST::default();
    let _ = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_CLEAR)?;
    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut ret = 0_i32;
        let num_at = pStruct.num_atoms;
        let num_edges = pBNS.num_bonds.wrapping_add(pBNS.num_atoms.wrapping_mul(2));
        let len_at = num_at.wrapping_add(pStruct.num_deleted_H);
        let inverse_mask = !forbidden_edge_mask;
        let _ = pnTotalDelta;

        let restore_mode = heap
            .slice(pStruct.pSrm)?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if pStruct.iInchiRec == INCHI_BAS as i8
            || restore_mode.bMetalAddFlower == 0
            || restore_mode.nMetalMinBondOrder != 0
        {
            return Ok(ret);
        }

        let copied_len =
            usize::try_from(len_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if copied_len != 0 {
            let original = heap
                .slice(at.as_const())?
                .get(..copied_len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(at2)?
                .get_mut(..copied_len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&original);
        }
        pStruct.at = at2;
        let copy_result = CopyBnsToAtom(heap, pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct.at = at;
        let copy_result = copy_result?;
        if copy_result < 0 {
            return Ok(copy_result);
        }

        let atom_count =
            usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atoms = heap
            .slice(at2.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();

        for atom_number in 0..num_at {
            if bIsMetalSalt(&atoms, atom_number)? == 0 {
                continue;
            }
            if all_charge_edges.num_edges == 0 {
                for current_atom in 0..num_at {
                    let current_index = usize::try_from(current_atom)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = atoms
                        .get(current_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(current_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut bond_position = 0_i32;
                    while bond_position < i32::from(atom.valence) {
                        let position = usize::try_from(bond_position)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let neighbor = i32::from(atom.neighbor[position]);
                        let edge_number = i32::from(
                            *heap
                                .slice(vertex.iedge.as_const())?
                                .get(position)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if neighbor < current_atom && edge.forbidden == 0 {
                            ret =
                                AddToEdgeList(heap, &mut all_charge_edges, edge_number, num_edges)?;
                            if ret != 0 {
                                return Ok(ret);
                            }
                        }
                        bond_position = bond_position.wrapping_add(1);
                    }

                    for edge_number in [
                        pVA.get(current_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCMinusGroupEdge
                            .wrapping_sub(1),
                        pVA.get(current_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCPlusGroupEdge
                            .wrapping_sub(1),
                    ] {
                        if edge_number >= 0 {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(edge_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge.forbidden == 0 {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut all_charge_edges,
                                    edge_number,
                                    num_edges,
                                )?;
                                if ret != 0 {
                                    return Ok(ret);
                                }
                            }
                        }
                    }
                }

                let mut group_number = 0_i32;
                while group_number < pTCGroups.num_tgroups {
                    let group = heap
                        .slice(pTCGroups.pTCG.as_const())?
                        .get(
                            usize::try_from(group_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let vertex_number = group.nVertexNumber;
                    let vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(
                            usize::try_from(vertex_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut edge_position = 0_i32;
                    while edge_position < i32::from(vertex.num_adj_edges) {
                        let edge_number = i32::from(
                            *heap
                                .slice(vertex.iedge.as_const())?
                                .get(
                                    usize::try_from(edge_position)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        ret = AddToEdgeList(heap, &mut all_charge_edges, edge_number, num_edges)?;
                        if ret != 0 {
                            return Ok(ret);
                        }
                        edge_position = edge_position.wrapping_add(1);
                    }
                    group_number = group_number.wrapping_add(1);
                }
            }

            let metal_index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let metal_atom = atoms
                .get(metal_index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let metal_vertex = heap
                .slice(pBNS.vert.as_const())?
                .get(metal_index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut bond_position = 0_i32;
            while bond_position < i32::from(metal_atom.valence) {
                let position = usize::try_from(bond_position)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let neighbor = usize::from(metal_atom.neighbor[position]);
                let neighbor_atom = atoms
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if neighbor_atom.charge != 0
                    || metal_atom.bond_type[position] != BOND_TYPE_SINGLE as u8
                {
                    bond_position = bond_position.wrapping_add(1);
                    continue;
                }

                let edge_number = i32::from(
                    *heap
                        .slice(metal_vertex.iedge.as_const())?
                        .get(position)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let edge_index = usize::try_from(edge_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let edge_before = heap
                    .slice(pBNS.edge.as_const())?
                    .get(edge_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge_number = pVA
                    .get(metal_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCPlusGroupEdge
                    .wrapping_sub(1);
                let minus_edge_number = pVA
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                if edge_before.flow != 1 || edge_before.forbidden != 0 || plus_edge_number < 0 {
                    bond_position = bond_position.wrapping_add(1);
                    continue;
                }
                let plus_index = usize::try_from(plus_edge_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let plus_before = heap
                    .slice(pBNS.edge.as_const())?
                    .get(plus_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if plus_before.flow <= 0 {
                    bond_position = bond_position.wrapping_add(1);
                    continue;
                }

                if minus_edge_number >= 0 {
                    let minus_index = usize::try_from(minus_edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let minus_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(minus_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if minus_before.flow != 0
                        || plus_before.forbidden != 0
                        || minus_before.forbidden != 0
                    {
                        bond_position = bond_position.wrapping_add(1);
                        continue;
                    }

                    let delta = 1_i32;
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[edge_index].flow =
                        edge_before.flow.wrapping_sub(delta);
                    for vertex_number in [first_vertex, second_vertex] {
                        let vertex = heap
                            .slice_mut(pBNS.vert)?
                            .get_mut(
                                usize::try_from(vertex_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(delta.wrapping_mul(2));

                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    {
                        let plus = &mut heap.slice_mut(pBNS.edge)?[plus_index];
                        plus.forbidden = (i32::from(plus.forbidden) & inverse_mask) as i8;
                    }
                    {
                        let minus = &mut heap.slice_mut(pBNS.edge)?[minus_index];
                        minus.forbidden = (i32::from(minus.forbidden) & inverse_mask) as i8;
                    }

                    let mut path_start = 0_i32;
                    let mut path_end = 0_i32;
                    let mut path_len = 0_i32;
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
                        &mut path_len,
                        &mut delta_h,
                        &mut delta_charge,
                        &mut visited_atoms,
                    )?;
                    if ret == 1
                        && ((path_end == first_vertex && path_start == second_vertex)
                            || (path_end == second_vertex && path_start == first_vertex))
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = heap
                            .slice(pBNS.edge.as_const())?
                            .get(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            .wrapping_add(delta);
                        for vertex_number in [first_vertex, second_vertex] {
                            let vertex = heap
                                .slice_mut(pBNS.vert)?
                                .get_mut(
                                    usize::try_from(vertex_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(delta.wrapping_mul(2));
                    }
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                } else {
                    let plus_minus = GetPlusMinusVertex(heap, pBNS, pTCGroups, 1, 1)?;
                    if plus_minus != NO_VERTEX {
                        let delta = 1_i32;
                        let plus_minus_index = usize::try_from(plus_minus)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        heap.slice_mut(pBNS.vert)?[plus_minus_index].st_edge.cap = heap
                            .slice(pBNS.vert.as_const())?[plus_minus_index]
                            .st_edge
                            .cap
                            .wrapping_add(delta);
                        heap.slice_mut(pBNS.vert)?[neighbor].st_edge.cap = heap
                            .slice(pBNS.vert.as_const())?[neighbor]
                            .st_edge
                            .cap
                            .wrapping_add(delta);
                        heap.slice_mut(pBNS.edge)?[edge_index].cap =
                            edge_before.cap.wrapping_add(delta);
                        pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_add(delta.wrapping_mul(2));

                        let plus_minus_vertex = heap
                            .slice(pBNS.vert.as_const())?
                            .get(plus_minus_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let mut first_position = 0_i32;
                        while first_position < i32::from(plus_minus_vertex.num_adj_edges) {
                            let connector_edge = i32::from(
                                *heap
                                    .slice(plus_minus_vertex.iedge.as_const())?
                                    .get(
                                        usize::try_from(first_position)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            );
                            let connector = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(connector_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let adjacent_vertex = i32::from(connector.neighbor12) ^ plus_minus;
                            let adjacent = heap
                                .slice(pBNS.vert.as_const())?
                                .get(
                                    usize::try_from(adjacent_vertex)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut second_position = 0_i32;
                            while second_position < i32::from(adjacent.num_adj_edges) {
                                let adjusted_edge = i32::from(
                                    *heap
                                        .slice(adjacent.iedge.as_const())?
                                        .get(
                                            usize::try_from(second_position)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                );
                                let adjusted = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(
                                        usize::try_from(adjusted_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                adjusted.cap = adjusted.cap.wrapping_add(delta);
                                second_position = second_position.wrapping_add(1);
                            }
                            first_position = first_position.wrapping_add(1);
                        }

                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        {
                            let plus = &mut heap.slice_mut(pBNS.edge)?[plus_index];
                            plus.forbidden = (i32::from(plus.forbidden) & inverse_mask) as i8;
                        }
                        {
                            let metal_bond = &mut heap.slice_mut(pBNS.edge)?[edge_index];
                            metal_bond.forbidden =
                                (i32::from(metal_bond.forbidden) & inverse_mask) as i8;
                        }

                        let mut path_start = 0_i32;
                        let mut path_end = 0_i32;
                        let mut path_len = 0_i32;
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
                            &mut path_len,
                            &mut delta_h,
                            &mut delta_charge,
                            &mut visited_atoms,
                        )?;
                        let mut current_success = 0_i32;
                        if ret == 1
                            && ((path_end == plus_minus && path_start == neighbor as i32)
                                || (path_end == neighbor as i32 && path_start == plus_minus))
                        {
                            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                            if ret > 0 {
                                *pnNumRunBNS = pnNumRunBNS.wrapping_add(1);
                                current_success = current_success.wrapping_add(1);
                            }
                        }
                        if current_success != 0 {
                            let current_flow = heap
                                .slice(pBNS.edge.as_const())?
                                .get(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .flow;
                            if current_flow != delta.wrapping_mul(2) {
                                return Ok(RI_ERR_PROGR);
                            }
                            heap.slice_mut(pBNS.edge)?[edge_index].flow =
                                current_flow.wrapping_sub(delta.wrapping_mul(2));
                            {
                                let oxygen = &mut heap.slice_mut(pBNS.vert)?[neighbor];
                                oxygen.st_edge.cap =
                                    oxygen.st_edge.cap.wrapping_sub(delta.wrapping_mul(2));
                                oxygen.st_edge.flow =
                                    oxygen.st_edge.flow.wrapping_sub(delta.wrapping_mul(2));
                            }
                            {
                                let metal = &mut heap.slice_mut(pBNS.vert)?[metal_index];
                                metal.st_edge.flow =
                                    metal.st_edge.flow.wrapping_sub(delta.wrapping_mul(2));
                                metal.st_edge.cap =
                                    metal.st_edge.cap.wrapping_sub(delta.wrapping_mul(2));
                            }
                            pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_sub(delta.wrapping_mul(3));
                            pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(delta.wrapping_mul(4));
                            heap.slice_mut(pBNS.edge)?[edge_index].cap = heap
                                .slice(pBNS.edge.as_const())?
                                .get(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .cap
                                .wrapping_sub(delta.wrapping_mul(2));
                            pVA.get_mut(neighbor)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .cInitCharge = pVA[neighbor].cInitCharge.wrapping_sub(delta as i8);
                        } else {
                            heap.slice_mut(pBNS.vert)?[plus_minus_index].st_edge.cap = heap
                                .slice(pBNS.vert.as_const())?[plus_minus_index]
                                .st_edge
                                .cap
                                .wrapping_sub(delta);
                            heap.slice_mut(pBNS.vert)?[neighbor].st_edge.cap = heap
                                .slice(pBNS.vert.as_const())?[neighbor]
                                .st_edge
                                .cap
                                .wrapping_sub(delta);
                            pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_sub(delta.wrapping_mul(2));

                            let plus_minus_vertex = heap
                                .slice(pBNS.vert.as_const())?
                                .get(plus_minus_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut first_position = 0_i32;
                            while first_position < i32::from(plus_minus_vertex.num_adj_edges) {
                                let connector_edge = i32::from(
                                    *heap
                                        .slice(plus_minus_vertex.iedge.as_const())?
                                        .get(
                                            usize::try_from(first_position)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                );
                                let connector = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(connector_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let adjacent_vertex = i32::from(connector.neighbor12) ^ plus_minus;
                                let adjacent = heap
                                    .slice(pBNS.vert.as_const())?
                                    .get(
                                        usize::try_from(adjacent_vertex)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let mut second_position = 0_i32;
                                while second_position < i32::from(adjacent.num_adj_edges) {
                                    let adjusted_edge =
                                        i32::from(
                                            *heap
                                                .slice(adjacent.iedge.as_const())?
                                                .get(usize::try_from(second_position).map_err(
                                                    |_| SourceHeapError::PointerOutOfBounds,
                                                )?)
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                        );
                                    let adjusted = heap
                                        .slice_mut(pBNS.edge)?
                                        .get_mut(
                                            usize::try_from(adjusted_edge)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    adjusted.cap = adjusted.cap.wrapping_sub(delta);
                                    second_position = second_position.wrapping_add(1);
                                }
                                first_position = first_position.wrapping_add(1);
                            }
                        }
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &all_charge_edges,
                            forbidden_edge_mask,
                        )?;
                    }
                }
                bond_position = bond_position.wrapping_add(1);
            }
        }
        Ok(ret)
    })();

    let cleanup = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_FREE);
    let result = execution?;
    let _ = cleanup?;
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RunBnsRestore1(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
    ip: &INPUT_PARMS,
    sd: &STRUCT_DATA,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &mut StrFromINChI,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    pInChI: [SourceMutPointer<INChI>; 2],
    num_inp: i64,
    bHasSomeFixedH: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2740 RunBnsRestore1
    /*
    int RunBnsRestore1(CANON_GLOBALS* pCG,
        INCHI_CLOCK* ic,
        ICHICONST INPUT_PARMS* ip,
        STRUCT_DATA* sd,
        BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups,
        INChI* pInChI[],
        long num_inp,
        int bHasSomeFixedH)
    {
        int        nNumRunBNS = 0;

        EDGE_LIST CarbonChargeEdges, MetalCarbonEdges, Nplus2BondsEdges;

        int nTotalDelta = 0, ret = 0; /* djb-rwth: removing redundant variables */
        inp_ATOM* at = pStruct->at;
        inp_ATOM* at2 = NULL; /* restored structure */
        inp_ATOM* at3 = NULL; /* structure for calculating one InChI */
        int     num_at = pStruct->num_atoms;
        int     num_deleted_H = pStruct->num_deleted_H;
    #ifdef _DEBUG
        int ret2;
    #endif

    #if ( KEEP_METAL_EDGE_FLOW == 1 )
        BNS_VERTEX* pVert, * pNeigh;
        int         j, neigh;
    #endif

        /* Edge lists initialization */
        AllocEdgeList(&CarbonChargeEdges, EDGE_LIST_CLEAR);
        AllocEdgeList(&MetalCarbonEdges, EDGE_LIST_CLEAR);
        AllocEdgeList(&Nplus2BondsEdges, EDGE_LIST_CLEAR);

        if (pStruct->iMobileH == TAUT_NON &&
            (ret = FillOutExtraFixedHDataInChI(pStruct, pInChI)))
        {
            goto exit_function;
        }

        if ((!at2 && !(at2 = (inp_ATOM*)inchi_malloc(((long long)num_at + (long long)num_deleted_H) * sizeof(at2[0])))) ||
            (!at3 && !(at3 = (inp_ATOM*)inchi_malloc(((long long)num_at + (long long)num_deleted_H) * sizeof(at3[0]))))) /* djb-rwth: cast operators added; addressing LLVM warning */
        {
            inchi_free(at2);
            inchi_free(at3);
            return RI_ERR_ALLOC;
        }

        if (0 > (ret = ForbidCarbonChargeEdges(pBNS, pTCGroups, &CarbonChargeEdges, BNS_EDGE_FORBIDDEN_TEMP)))
        {
            goto exit_function;
        }

    #if ( KEEP_METAL_EDGE_FLOW == 1 )
        /* count edges of -C(IV)< carbons connected to metals */
        if (0 > (ret = ForbidMetalCarbonEdges(pBNS, at, num_at, pVA, pTCGroups, &MetalCarbonEdges, BNS_EDGE_FORBIDDEN_TEMP)))
        {
            goto exit_function;
        }
    #endif
        if (0 > (ret = ForbidNintrogenPlus2BondsInSmallRings(pBNS, at, num_at, pVA, 6,
            pTCGroups, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP)))
        {
            goto exit_function;
        }

        /*********** Run BNS #1: no charge on carbons and =N= ***************/
        if (Nplus2BondsEdges.num_edges)
        {
            /* Run BNS leaving carbon charges unchanged */
            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
            nNumRunBNS++;
            if (ret < 0)
            {
                goto exit_function;
            }
            else
            {
                nTotalDelta += ret;
            }
            RemoveForbiddenEdgeMask(pBNS, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP);
            AllocEdgeList(&Nplus2BondsEdges, EDGE_LIST_FREE);
        }
    #ifdef _DEBUG
        /* debug only */
        memcpy(at2, at, ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operators added */
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
    #endif
        /*************************** extend min ring size to 8 ****************************/
        if (0 > (ret = ForbidNintrogenPlus2BondsInSmallRings(pBNS, at, num_at, pVA, 8,
            pTCGroups, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP)))
        {
            goto exit_function;
        }
        if (Nplus2BondsEdges.num_edges)
        {
            /* Run BNS leaving carbon charges unchanged */
            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
            nNumRunBNS++;
            if (ret < 0)
            {
                goto exit_function;
            }
            else
            {
                nTotalDelta += ret;
            }
            RemoveForbiddenEdgeMask(pBNS, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP);
            AllocEdgeList(&Nplus2BondsEdges, EDGE_LIST_FREE);
        }
    #ifdef _DEBUG
        /* debug only */
        memcpy(at2, at, ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operators added */
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
    #endif
        /*******************************************************************/
        if (CarbonChargeEdges.num_edges > 0)
        {
            /* Run BNS leaving carbon charges unchanged */
            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
            nNumRunBNS++;
            if (ret < 0)
            {
                goto exit_function;
            }
            else
            {
                nTotalDelta += ret;
            }
            RemoveForbiddenEdgeMask(pBNS, &CarbonChargeEdges, BNS_EDGE_FORBIDDEN_TEMP);
            AllocEdgeList(&CarbonChargeEdges, EDGE_LIST_FREE);
        }
    #ifdef _DEBUG
        /* debug only */
        memcpy(at2, at, ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operators added */
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
    #endif
        /*******************************************************************/
        if (MetalCarbonEdges.num_edges > 0)
        {
            /* Run BNS leaving carbon charges unchanged */
            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
            nNumRunBNS++;
            if (ret < 0)
            {
                goto exit_function;
            }
            else
            {
                nTotalDelta += ret;
            }
            RemoveForbiddenEdgeMask(pBNS, &MetalCarbonEdges, BNS_EDGE_FORBIDDEN_TEMP);
            AllocEdgeList(&MetalCarbonEdges, EDGE_LIST_FREE);
        }
        /*******************************************************************/
        /* Run BNS allowing to change any charges */
        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
        nNumRunBNS++;
        if (ret < 0)
        {
            goto exit_function;
        }
        else
        {
            nTotalDelta += ret;
        }
    #ifdef _DEBUG
        /* debug only */
        memcpy(at2, at, ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operators added */
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
    #endif

    #if ( BNS_RAD_SEARCH == 1 )
        /****************************************************************************/
        /* move unfulfilled 'radicals' from ChargeStruct to atoms         */
        /* and set change charges of affected atoms to fit total charge   */
        ret = MoveRadToAtomsAddCharges(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
    #endif
        /**************************************************************/
        /**************************************************************/
        /*****           fix restore inconsistencies              *****/
        /**************************************************************/
        /**************************************************************/
    #ifdef _DEBUG
        /* debug only */
        memcpy(at2, at, ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operators added */
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
    #endif

        /* rearrange (+) and (-) edges flow so that there is no (+)flow=0 and (-)flow=1 */
        ret = RearrangePlusMinusEdgesFlow(pBNS, pBD, pVA, pTCGroups, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }

        /*****************************************************************/
        /*       Increment zero order metal bonds to heteroatoms         */
        /*****************************************************************/
        ret = IncrementZeroOrderBondsToHeteroat(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }

    #ifdef _DEBUG
        /* debug only */
        memcpy(at2, at, ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operators added */
        pStruct->at = at2;
        ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
        pStruct->at = at;
    #endif

    #if (MOVE_CHARGES_FROM_HETEREO_TO_METAL == 1 )
        /*****************************************************************/
        /* move charges from heteroatoms to metal atoms                  */
        /*****************************************************************/
        ret = MoveChargeFromHeteroatomsToMetals(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
    #endif
        /***********************************************************************
                NH2                NH2
                   \                  \
                    C==S(+)-   =>      C(+)-S-   where NH2 are not tautomeric
                   /                  /
                NH2                NH2
        ************************************************************************/
        ret = MovePlusFromS2DiaminoCarbon(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /*****************************************************************/
        /*       Avoid charge separation on heteroatoms                  */
        /*****************************************************************/
        ret = EliminateChargeSeparationOnHeteroatoms(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP, 0);
        if (ret < 0)
        {
            goto exit_function;
        }
        if (ret)
        {
            /*charge separation remains; allow changes of stereobonds in a ring and try again */
            ret = EliminateChargeSeparationOnHeteroatoms(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP,
                BNS_EDGE_FORBIDDEN_MASK);
            if (ret < 0)
            {
                goto exit_function;
            }
        }
        /*****************************************************************/
        /*         convert N#N(+)-N= into N(-)=N(+)=N-                   */
        /*****************************************************************/
        ret = RestoreNNNgroup(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /*****************************************************************/
        /*     convert Metal(q)-N(-)-O(-) Metal(q-2)-N=O (local change)  */
        /*****************************************************************/
        ret = FixMetal_Nminus_Ominus(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /*****************************************************************/
        /*         convert N(-)=C= into N#C-         -                   */
        /*****************************************************************/
        ret = RestoreCyanoGroup(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /*****************************************************************/
        /*         convert C(+)#N(+)- into C(-)#N(+)-                    */
        /*****************************************************************/
        ret = RestoreIsoCyanoGroup(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /*****************************************************************/
        /*         eliminate =N(V)= if possible                          */
        /*                    |                                          */
        /*****************************************************************/
        ret = EliminateNitrogen5Val3Bonds(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }

        /*****************************************************************/
        /*                    |      |                                   */
        /*         convert   -S- to =S= if possible                      */
        /*                    |      |                                   */
        /*****************************************************************/
        ret = Convert_SIV_to_SVI(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }

        /*****************************************************************/
        /*                  =N(+)=O     =N-O(-)                          */
        /*         convert           => if possible                      */
        /*                  Metal(q)    Metal(q+2)                       */
        /*****************************************************************/
        ret = PlusFromDB_N_DB_O_to_Metal(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }

        /*****************************************************************/
        /*  forbidden edges prevents required in InChI tautomerism       */
        /*  incorrectly restored mobile H mix separate tautomeric groups */
        /*  because an edge may not become forbidden                     */
        /* note: removes this 'forbidden_edge' bit from ALL edges        */
        /*****************************************************************/
        ret = MoveMobileHToAvoidFixedBonds(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);

        if (ret < 0)
        {
            goto exit_function;
        }
        /**************************************************************************/
        /* 2. Mobile H endpoint has radical on it (typical for wrong P(VI)(=O)3OH */
        /* djb-rwth: removing redundant code */
        if (pStruct->iMobileH == TAUT_NON)
        {
            ret = RemoveRadFromMobileHEndpointFixH(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        }
        else
        {
            ret = RemoveRadFromMobileHEndpoint(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        }
        if (ret < 0)
        {
            goto exit_function;
        }
        /* djb-rwth: removing redundant code */
        /**************************************************************/
        /* make bonds between a charged heteroatom and a metal single */
        ret = MakeSingleBondsMetal2ChargedHeteroat(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /**************************************************************/
        /* move (+) charges to >N- and other centerpoints             */
        ret = MoveChargeToMakeCenerpoints(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }

        /**************************************************************************/
        /* Find and eliminate false Mobile-H groups: Cl(=O)3(-O(-)) => Cl(-)(=O)4 */
        ret = MoveChargeToRemoveCenerpoints(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /**************************************************************************/
        /* Find A=X< where all bonds to X except A=X are marked as stereogenic    */
        /* make bonds A=X single                                                  */
        ret = CheckAndRefixStereobonds(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /**************************************************************************/
        /* In Reconnected structure change 'salt bonds' to 'coordination bonds    */
        /* for example, M-O-C=  ->  M(+)-O(-)-C=                                  */
        /* Defect: instead of NH2-C=O(+)-M it will restore NH2(+)=C-O(-)-M(+)     */
        /* However, in this release metal-organic compounds do not get much care  */
        ret = SaltBondsToCoordBonds(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
        if (ret < 0)
        {
            goto exit_function;
        }
        /**************************************************************************/
        /* Normalize the structure and compare t-groups and stereobonds           */
        ret = NormalizeAndCompare(pCG, ic, ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups, pInChI, num_inp, bHasSomeFixedH,
            &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP, BNS_EDGE_FORBIDDEN_MASK);
        if (ret < 0)
        {
            goto exit_function;
        }
        /**************************************************************************/
        /* Create InChI out of the restored structure                             */


        /*ret = nTotalDelta;*/

    exit_function:

        pStruct->at = at;
        pStruct->at2 = at2;
        at2 = NULL;
        AllocEdgeList(&CarbonChargeEdges, EDGE_LIST_FREE);
        AllocEdgeList(&MetalCarbonEdges, EDGE_LIST_FREE);
        AllocEdgeList(&Nplus2BondsEdges, EDGE_LIST_FREE);
        if (at2)
        {
            inchi_free(at2);
        }
        if (at3)
        {
            inchi_free(at3);
        }

        return ret;
    }
        */
    // END INCHI C FUNCTION: RunBnsRestore1
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RunBnsRestore1
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; _DEBUG is undefined.
    // INCHI✔❌: BNS_RAD_SEARCH=1; KEEP_METAL_EDGE_FLOW=0; MOVE_CHARGES_FROM_HETEREO_TO_METAL=0.
    // INCHI✔❌: SourceHeap checked ownership and typed atom-buffer allocation add overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: RunBnsRestore1

    let at = pStruct.at;
    let num_at = pStruct.num_atoms;
    let num_deleted_h = pStruct.num_deleted_H;
    let mut carbon_charge_edges = EDGE_LIST::default();
    let mut metal_carbon_edges = EDGE_LIST::default();
    let mut nplus_two_bonds_edges = EDGE_LIST::default();
    let _ = AllocEdgeList(heap, &mut carbon_charge_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut metal_carbon_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut nplus_two_bonds_edges, EDGE_LIST_CLEAR)?;

    if pStruct.iMobileH == TAUT_NON as i8 {
        let ret = FillOutExtraFixedHDataInChI(heap, pStruct, pInChI)?;
        if ret != 0 {
            pStruct.at = at;
            pStruct.at2 = SourceMutPointer::null();
            let _ = AllocEdgeList(heap, &mut carbon_charge_edges, EDGE_LIST_FREE)?;
            let _ = AllocEdgeList(heap, &mut metal_carbon_edges, EDGE_LIST_FREE)?;
            let _ = AllocEdgeList(heap, &mut nplus_two_bonds_edges, EDGE_LIST_FREE)?;
            return Ok(ret);
        }
    }

    let atom_count_i64 = i64::from(num_at).wrapping_add(i64::from(num_deleted_h));
    let atom_count = usize::try_from(atom_count_i64)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    u64::try_from(atom_count)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
        .checked_mul(std::mem::size_of::<inp_ATOM>() as u64)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;

    let at2 = match heap.allocate(vec![inp_ATOM::default(); atom_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(RI_ERR_ALLOC),
        Err(error) => return Err(error),
    };
    let at3 = match heap.allocate(vec![inp_ATOM::default(); atom_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_free(heap, at2)?;
            return Ok(RI_ERR_ALLOC);
        }
        Err(error) => {
            inchi_free(heap, at2)?;
            return Err(error);
        }
    };

    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut ret: i32;
        let mut number_bns_runs = 0_i32;
        let mut total_delta = 0_i32;
        let atoms = heap
            .slice(at.as_const())?
            .get(..usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();

        let forbid_carbon_result = ForbidCarbonChargeEdges(
            heap,
            pBNS,
            pTCGroups,
            &mut carbon_charge_edges,
            BNS_EDGE_FORBIDDEN_TEMP as i32,
        );
        ret = forbid_carbon_result?;
        if ret < 0 {
            return Ok(ret);
        }

        let forbid_n6_result = ForbidNintrogenPlus2BondsInSmallRings(
            heap,
            pBNS,
            &atoms,
            num_at,
            pVA,
            6,
            pTCGroups,
            &mut nplus_two_bonds_edges,
            BNS_EDGE_FORBIDDEN_TEMP as i32,
        );
        ret = forbid_n6_result?;
        if ret < 0 {
            return Ok(ret);
        }
        if nplus_two_bonds_edges.num_edges != 0 {
            let result = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result);
            ret = result?;
            number_bns_runs = number_bns_runs.wrapping_add(1);
            if ret < 0 {
                return Ok(ret);
            }
            total_delta = total_delta.wrapping_add(ret);
            RemoveForbiddenEdgeMask(
                heap,
                pBNS,
                &nplus_two_bonds_edges,
                BNS_EDGE_FORBIDDEN_TEMP as i32,
            )?;
            let _ = AllocEdgeList(heap, &mut nplus_two_bonds_edges, EDGE_LIST_FREE)?;
        }

        let forbid_n8_result = ForbidNintrogenPlus2BondsInSmallRings(
            heap,
            pBNS,
            &atoms,
            num_at,
            pVA,
            8,
            pTCGroups,
            &mut nplus_two_bonds_edges,
            BNS_EDGE_FORBIDDEN_TEMP as i32,
        );
        ret = forbid_n8_result?;
        if ret < 0 {
            return Ok(ret);
        }
        if nplus_two_bonds_edges.num_edges != 0 {
            let result = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result);
            ret = result?;
            number_bns_runs = number_bns_runs.wrapping_add(1);
            if ret < 0 {
                return Ok(ret);
            }
            total_delta = total_delta.wrapping_add(ret);
            RemoveForbiddenEdgeMask(
                heap,
                pBNS,
                &nplus_two_bonds_edges,
                BNS_EDGE_FORBIDDEN_TEMP as i32,
            )?;
            let _ = AllocEdgeList(heap, &mut nplus_two_bonds_edges, EDGE_LIST_FREE)?;
        }

        if carbon_charge_edges.num_edges > 0 {
            let result = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result);
            ret = result?;
            number_bns_runs = number_bns_runs.wrapping_add(1);
            if ret < 0 {
                return Ok(ret);
            }
            total_delta = total_delta.wrapping_add(ret);
            RemoveForbiddenEdgeMask(
                heap,
                pBNS,
                &carbon_charge_edges,
                BNS_EDGE_FORBIDDEN_TEMP as i32,
            )?;
            let _ = AllocEdgeList(heap, &mut carbon_charge_edges, EDGE_LIST_FREE)?;
        }

        if metal_carbon_edges.num_edges > 0 {
            let result = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result);
            ret = result?;
            number_bns_runs = number_bns_runs.wrapping_add(1);
            if ret < 0 {
                return Ok(ret);
            }
            total_delta = total_delta.wrapping_add(ret);
            RemoveForbiddenEdgeMask(
                heap,
                pBNS,
                &metal_carbon_edges,
                BNS_EDGE_FORBIDDEN_TEMP as i32,
            )?;
            let _ = AllocEdgeList(heap, &mut metal_carbon_edges, EDGE_LIST_FREE)?;
        }

        let result = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result);
        ret = result?;
        number_bns_runs = number_bns_runs.wrapping_add(1);
        if ret < 0 {
            return Ok(ret);
        }
        total_delta = total_delta.wrapping_add(ret);

        let move_rad_result = MoveRadToAtomsAddCharges(
            heap,
            pBNS,
            pBD,
            pStruct,
            at,
            at2,
            pVA,
            pTCGroups,
            BNS_EDGE_FORBIDDEN_TEMP as i32,
            clock_result,
        );
        ret = move_rad_result?;
        if ret < 0 {
            return Ok(ret);
        }

        let rearrange_result = RearrangePlusMinusEdgesFlow(
            heap,
            pBNS,
            pBD,
            pVA,
            pTCGroups,
            BNS_EDGE_FORBIDDEN_TEMP as i32,
            clock_result,
        );
        ret = rearrange_result?;
        if ret < 0 {
            return Ok(ret);
        }

        macro_rules! run_structure_helper {
            ($function:path) => {{
                let helper_result = $function(
                    heap,
                    pBNS,
                    pBD,
                    pStruct,
                    at,
                    at2,
                    pVA,
                    pTCGroups,
                    &mut number_bns_runs,
                    &mut total_delta,
                    BNS_EDGE_FORBIDDEN_TEMP as i32,
                    clock_result,
                );
                ret = helper_result?;
                if ret < 0 {
                    return Ok(ret);
                }
            }};
        }

        run_structure_helper!(IncrementZeroOrderBondsToHeteroat);
        run_structure_helper!(MovePlusFromS2DiaminoCarbon);

        ret = EliminateChargeSeparationOnHeteroatoms(
            heap,
            pBNS,
            pBD,
            pStruct,
            at,
            at2,
            pVA,
            pTCGroups,
            &mut number_bns_runs,
            &mut total_delta,
            BNS_EDGE_FORBIDDEN_TEMP as i32,
            0,
            clock_result,
        )?;
        if ret < 0 {
            return Ok(ret);
        }
        if ret != 0 {
            ret = EliminateChargeSeparationOnHeteroatoms(
                heap,
                pBNS,
                pBD,
                pStruct,
                at,
                at2,
                pVA,
                pTCGroups,
                &mut number_bns_runs,
                &mut total_delta,
                BNS_EDGE_FORBIDDEN_TEMP as i32,
                BNS_EDGE_FORBIDDEN_MASK as i32,
                clock_result,
            )?;
            if ret < 0 {
                return Ok(ret);
            }
        }

        run_structure_helper!(RestoreNNNgroup);
        run_structure_helper!(FixMetal_Nminus_Ominus);
        run_structure_helper!(RestoreCyanoGroup);
        run_structure_helper!(RestoreIsoCyanoGroup);
        run_structure_helper!(EliminateNitrogen5Val3Bonds);
        run_structure_helper!(Convert_SIV_to_SVI);
        run_structure_helper!(PlusFromDB_N_DB_O_to_Metal);
        run_structure_helper!(MoveMobileHToAvoidFixedBonds);

        ret = if pStruct.iMobileH == TAUT_NON as i8 {
            RemoveRadFromMobileHEndpointFixH(
                heap,
                pBNS,
                pBD,
                pStruct,
                at,
                at2,
                pVA,
                pTCGroups,
                &mut number_bns_runs,
                &mut total_delta,
                BNS_EDGE_FORBIDDEN_TEMP as i32,
                clock_result,
            )?
        } else {
            RemoveRadFromMobileHEndpoint(
                heap,
                pBNS,
                pBD,
                pStruct,
                at,
                at2,
                pVA,
                pTCGroups,
                &mut number_bns_runs,
                &mut total_delta,
                BNS_EDGE_FORBIDDEN_TEMP as i32,
                clock_result,
            )?
        };
        if ret < 0 {
            return Ok(ret);
        }

        run_structure_helper!(MakeSingleBondsMetal2ChargedHeteroat);
        run_structure_helper!(MoveChargeToMakeCenerpoints);
        run_structure_helper!(MoveChargeToRemoveCenerpoints);
        run_structure_helper!(CheckAndRefixStereobonds);
        run_structure_helper!(SaltBondsToCoordBonds);

        let normalize_result = NormalizeAndCompare(
            heap,
            pCG,
            ic,
            ip,
            sd,
            pBNS,
            pBD,
            pStruct,
            at,
            at2,
            at3,
            pVA,
            pTCGroups,
            pInChI,
            num_inp,
            bHasSomeFixedH,
            &mut number_bns_runs,
            &mut total_delta,
            BNS_EDGE_FORBIDDEN_TEMP as i32,
            BNS_EDGE_FORBIDDEN_MASK as i32,
            clock_result,
        );
        ret = normalize_result?;
        Ok(ret)
    })();

    pStruct.at = at;
    pStruct.at2 = at2;
    let free_carbon = AllocEdgeList(heap, &mut carbon_charge_edges, EDGE_LIST_FREE);
    let free_metal = AllocEdgeList(heap, &mut metal_carbon_edges, EDGE_LIST_FREE);
    let free_nplus = AllocEdgeList(heap, &mut nplus_two_bonds_edges, EDGE_LIST_FREE);
    let free_at3 = inchi_free(heap, at3);
    let ret = execution?;
    let _ = free_carbon?;
    let _ = free_metal?;
    let _ = free_nplus?;
    free_at3?;
    Ok(ret)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RestoreAtomMakeBNS(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: &mut CANON_GLOBALS,
    ip: &INPUT_PARMS,
    sd: &mut STRUCT_DATA,
    pStruct: &mut StrFromINChI,
    iComponent: i32,
    iAtNoOffset: i32,
    pInChI: [SourceMutPointer<INChI>; 2],
    _szCurHdr: SourceConstPointer<i8>,
    num_inp: i64,
    bHasSomeFixedH: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3198 RestoreAtomMakeBNS
    // INCHI✔❌: complete source frame follows verbatim.
    /*
    int RestoreAtomMakeBNS(INCHI_CLOCK* ic, CANON_GLOBALS* pCG,
        ICHICONST INPUT_PARMS* ip,
        STRUCT_DATA* sd,
        StrFromINChI* pStruct,
        int iComponent,
        int iAtNoOffset,
        INChI* pInChI[],
        const char* szCurHdr,
        long num_inp,
        int bHasSomeFixedH)
    {
        int i, j, ret = 0, ret2;
        /*int nDelta, nTotalDelta;*/
        VAL_AT* pVA = NULL;
        VAL_AT    va1;
        int    num_at = pStruct->num_atoms;
        inp_ATOM* at = pStruct->at;
        ALL_TC_GROUPS   TCGroups;
        ALL_TC_GROUPS* pTCGroups = &TCGroups;
        int            nAddEdges2eachAtom = 2, nAddVertices = 0;

        BFS_Q bfsq;

        /* BNS creation */
        BN_STRUCT* pBNS = NULL;
        BN_DATA* pBD = NULL;
        int            nNum_changed_bonds = 0;
        int            bTreatMoreAtomsAsMetals = 0, bSecondPassNewMetals = 0;
        int            nMaxAddAtoms = 2, nMaxAddEdges = 2, max_altp = BN_MAX_ALTP;

        memset(pTCGroups, 0, sizeof(pTCGroups[0])); /* djb-rwth: memset_s C11/Annex K variant? */
        for (i = 0; i < NUM_TCGROUP_TYPES; i++)
        {
            pTCGroups->nGroup[i] = TCG_None; /* unassigned */
        }
        pTCGroups->iComponent = iComponent;
        pTCGroups->iAtNoOffset = iAtNoOffset;

        if (num_at == 1)
        {
            /* single atom -- no bonds to restore */
            inp_ATOM* at2 = (inp_ATOM*)inchi_malloc(sizeof(at2[0]) * ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H)); /* djb-rwth: cast operators added */
            inp_ATOM* at3 = (inp_ATOM*)inchi_malloc(sizeof(at3[0]) * ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H)); /* djb-rwth: cast operators added */
            pStruct->at2 = at2;
            at[0].charge = pInChI[0]->nTotalCharge;
            if (at2)
            {
                memcpy(at2, at, sizeof(at2[0]) * ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H)); /* djb-rwth: cast operators added */
            }
            if (!at2 || !at3)
            {
                if (at3) inchi_free(at3);
                return RI_ERR_ALLOC;
            }
            ret = MakeOneInChIOutOfStrFromINChI(pCG, ic, ip, sd, pStruct, pStruct->at2, at3, pTCGroups);
            /* clean up */
            for (i = 0; i < TAUT_NUM; i++)
            {
                Free_INChI(&pStruct->pOneINChI[i]);
                Free_INChI_Aux(&pStruct->pOneINChI_Aux[i]);
                FreeInpAtomData(pStruct->pOne_norm_data[i]);
                if (pStruct->pOne_norm_data[i])
                {
                    inchi_free(pStruct->pOne_norm_data[i]);
                    pStruct->pOne_norm_data[i] = NULL;
                }
            }
            /* djb-rwth: fixing oss-fuzz issue #69602 */
            /* free_t_group_info(&pStruct->One_ti); */
            inchi_free(at3);

            return ret;
        }

        AllocBfsQueue(&bfsq, BFS_Q_CLEAR, 0);
        if (!(pVA = (VAL_AT*)inchi_calloc(num_at, sizeof(pVA[0]))))
        {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
        pStruct->pVA = pVA;
        memset(&va1, 0, sizeof(va1)); /* djb-rwth: memset_s C11/Annex K variant? */
        pTCGroups->total_charge = pInChI[0]->nTotalCharge;
        if (0 > (ret = AllocBfsQueue(&bfsq, num_at, 0 /* min ring size undefined */)))
        {
            goto exit_function;
        }
        pStruct->pbfsq = &bfsq;

        if (pStruct->iMobileH == TAUT_NON && pInChI[1] && pInChI[1]->nNumberOfAtoms > 1 &&
            (ret = FillOutpStructEndpointFromInChI(pInChI[1], &pStruct->endpoint)))
        {
            goto exit_function;
        }

        /* mark metal atoms; find min ring sizes for atoms that have 2 bonds */
        for (i = 0; i < num_at; i++)
        {
            pVA[i].cNumValenceElectrons = get_sp_element_type(at[i].el_number, &j);
            pVA[i].cPeriodicRowNumber = j;
            pVA[i].cPeriodicNumber = at[i].el_number;
            pVA[i].cNumValenceElectrons--; /* = -1 d- and f- metals, 0 for H, 1 for Na, 2 for Mg,.. = (ATYPE_Xx-1)  */

            if (is_el_a_metal(at[i].el_number))
            {
                if (pStruct->pSrm->bStereoRemovesMetalFlag)
                {
                    /* treat metal as non-metal if it is stereogenic or has a stereobond */
                    pVA[i].cMetal = !(at[i].p_parity || at[i].sb_parity[0]);
                }
                else
                {
                    pVA[i].cMetal = 1;
                }
            }
            if (at[i].valence == 2 && !at[i].num_H)
            {
                pVA[i].cMinRingSize = is_bond_in_Nmax_memb_ring(at, i, 0, bfsq.q, bfsq.nAtomLevel,
                    bfsq.cSource, 99 /* max ring size */);
            }
            else
            {
                pVA[i].cMinRingSize = 0;
            }
        }
        /* AllocBfsQueue( &bfsq, BFS_Q_FREE, 0 ); */

    repeat_for_new_metals:
        /* set valences for the first time; find ChargeValence structures for each atom */
        for (i = 0; i < num_at; i++)
        {
            /* get additional fictitious atoms information */
            pVA[i].cInitFreeValences = 0;

            ret = GetAtomRestoreInfo(pCG, at, i, pVA, pStruct->pSrm, pStruct->bMobileH, pStruct->endpoint);

            if (ret < 0)
            {
                goto exit_function;
            }
            if (ret == TREAT_ATOM_AS_METAL && !bSecondPassNewMetals && !pVA[i].cMetal)
            {
                if (pStruct->pSrm->bStereoRemovesMetalFlag)
                {
                    /* treat metal as non-metal if it is stereogenic or has a stereobond */
                    pVA[i].cMetal = !(at[i].p_parity || at[i].sb_parity[0]);
                }
                else
                {
                    pVA[i].cMetal = 1;
                }
                if (pVA[i].cMetal)
                {
                    bTreatMoreAtomsAsMetals++;
                }
            }
            pTCGroups->charge_on_atoms += pVA[i].cInitCharge;
        }
        if (bTreatMoreAtomsAsMetals && !bSecondPassNewMetals)
        {
            for (i = 0; i < num_at; i++)
            {
                /* clear all members of pVA[i] except two */
                pTCGroups->charge_on_atoms -= pVA[i].cInitCharge;
                va1.cMetal = pVA[i].cMetal;
                va1.cMinRingSize = pVA[i].cMinRingSize;
                va1.cNumValenceElectrons = pVA[i].cNumValenceElectrons;
                va1.cPeriodicRowNumber = pVA[i].cPeriodicRowNumber;
                va1.cPeriodicNumber = pVA[i].cPeriodicNumber;
                pVA[i] = va1;
            }
            bSecondPassNewMetals = 1;
            goto repeat_for_new_metals;
        }

        /* count atoms, bonds, additional edges and vertices in ChargeValence structures and t-groups */
        ret = nCountBnsSizes(at, num_at, nAddEdges2eachAtom, nAddVertices, &pStruct->ti,
            pVA, pStruct->pSrm, pTCGroups);
        if (ret < 0)
        {
            goto exit_function;
        }

        /* find and count groups; add counts of all other vertices to be created */
        ret = nAddSuperCGroups(pTCGroups);
        if (ret < 0)
        {
            goto exit_function;
        }

        /* create the BNS and fill it with all real atoms */
        pBNS = AllocateAndInitTCGBnStruct(pStruct, pVA, pTCGroups,
            nMaxAddAtoms, nMaxAddEdges, max_altp, &nNum_changed_bonds);
        if (!pBNS)
        {
            ret = BNS_OUT_OF_RAM;
            goto exit_function;
        }
        /* add t-groups to the BNS */
        ret = AddTGroups2TCGBnStruct(pBNS, pStruct, pVA, pTCGroups, nMaxAddEdges);
        if (ret < 0)
        {
            goto exit_function;
        }

        /* add c-groups to the BNS; adjust charges */
        ret = AddCGroups2TCGBnStruct(pBNS, pStruct, pVA, pTCGroups, nMaxAddEdges);
        if (ret < 0)
        {
            goto exit_function;
        }

        pBNS->ulTimeOutTime = NULL;             /* v. 1.05 */
        pBNS->ic = ic;                          /* v. 1.05 */

        /* allocate BNData */
        pBD = AllocateAndInitBnData(pBNS->max_vertices + pBNS->max_vertices / 2);
        if (!pBD)
        {
            ret = BNS_OUT_OF_RAM;
            goto exit_function;
        }
        CheckBnsConsistency(pStruct, pBNS, pVA, pTCGroups, 0);

        /* restore bonds & charges */
        ret = RunBnsRestore1(pCG, ic, ip, sd, pBNS, pBD, pStruct, pVA, pTCGroups, pInChI, num_inp, bHasSomeFixedH);
        if (ret < 0)
        {
            goto exit_function;
        }

        ret = CheckBnsConsistency(pStruct, pBNS, pVA, pTCGroups, 1);
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
        if (ret)
        {
            fprintf(stdout, "Msg for: %ld %s comp=%d %c%c\n", num_inp, (szCurHdr && szCurHdr[0]) ? szCurHdr : "", iComponent, pStruct->iInchiRec ? 'R' : 'D', pStruct->iMobileH ? 'M' : 'F');
        }
        if (pStruct->iMobileH == TAUT_YES && pStruct->nNumRemovedProtons)
        {
            fprintf(stdout, "REMOVED_PROTONS%+d %ld %s\n", pStruct->nNumRemovedProtons, num_inp, (szCurHdr && szCurHdr[0]) ? szCurHdr : "");
            /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
        }
        if (pStruct->bExtract & EXTRACT_STRUCT_NUMBER)
        {
            fprintf(stdout, "EXTRACT: %ld: %s\n", num_inp, (szCurHdr && szCurHdr[0]) ? szCurHdr : "");
        }
    #endif
    #endif
        {  /* create the final structure in pStruct->at2 */
            inp_ATOM* at_tmp = pStruct->at;
            pStruct->at = pStruct->at2;
            memcpy(pStruct->at, at_tmp, sizeof(pStruct->at[0]) * ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H)); /* djb-rwth: cast operators added */
            ret2 = CopyBnsToAtom(pStruct, pBNS, pVA, pTCGroups, 1);
            pStruct->at2 = pStruct->at;
            pStruct->at = at_tmp;
            if (ret2 < 0)
            {
                ret = ret2;
            }
        }

    exit_function:

        pStruct->pbfsq = NULL;
        AllocBfsQueue(&bfsq, BFS_Q_FREE, 0);

        pBD = DeAllocateBnData(pBD); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        pBNS = DeAllocateBnStruct(pBNS); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        /*
        if ( pVA ) inchi_free( pVA );
        */
        if (pTCGroups->pTCG) inchi_free(pTCGroups->pTCG);

        return ret;
    }
    */
    // END INCHI C FUNCTION: RestoreAtomMakeBNS
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RestoreAtomMakeBNS
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; bRELEASE_VERSION=1.
    // INCHI✔❌: TARGET_API_LIB excludes the diagnostic fprintf block; BN_MAX_ALTP=16.
    // INCHI✔❌: SourceHeap checked pointer access and model storage for the stack BFS_Q add overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: RestoreAtomMakeBNS

    fn atom_count(structure: &StrFromINChI) -> Result<usize, SourceHeapError> {
        let count = i64::from(structure.num_atoms).wrapping_add(i64::from(structure.num_deleted_H));
        let count = usize::try_from(count)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        (count as u64)
            .checked_mul(std::mem::size_of::<inp_ATOM>() as u64)
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        Ok(count)
    }

    fn source_atom_allocation(
        heap: &mut SourceHeap,
        count: usize,
    ) -> Result<SourceMutPointer<inp_ATOM>, SourceHeapError> {
        match heap.allocate(vec![inp_ATOM::default(); count]) {
            Ok(pointer) => Ok(pointer),
            Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
            Err(error) => Err(error),
        }
    }

    let num_at = pStruct.num_atoms;
    let at = pStruct.at;
    let mut groups = ALL_TC_GROUPS::default();
    groups.nGroup.fill(TCG_None);
    groups.iComponent = iComponent;
    groups.iAtNoOffset = iAtNoOffset;

    if num_at == 1 {
        let count = match atom_count(pStruct) {
            Ok(count) => count,
            Err(
                SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationSizeOverflow,
            ) => {
                return Ok(RI_ERR_ALLOC);
            }
            Err(error) => return Err(error),
        };
        let at2 = source_atom_allocation(heap, count)?;
        let at3 = source_atom_allocation(heap, count)?;
        pStruct.at2 = at2;

        let total_charge = heap
            .slice(pInChI[0].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nTotalCharge;
        heap.slice_mut(at)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .charge = total_charge as i8;
        if !at2.is_null() {
            let source = heap
                .slice(at.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(at2)?
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&source);
        }
        if at2.is_null() || at3.is_null() {
            if !at3.is_null() {
                inchi_free(heap, at3)?;
            }
            return Ok(RI_ERR_ALLOC);
        }

        let result = MakeOneInChIOutOfStrFromINChI(
            heap,
            pCG,
            ic,
            ip,
            sd,
            pStruct,
            at2,
            at3,
            &groups,
            clock_result,
        );
        let cleanup = (|| -> Result<(), SourceHeapError> {
            for index in 0..TAUT_NUM as usize {
                Free_INChI(heap, &mut pStruct.pOneINChI[index])?;
                Free_INChI_Aux(heap, &mut pStruct.pOneINChI_Aux[index])?;
                let normalized = pStruct.pOne_norm_data[index];
                if !normalized.is_null() {
                    let mut data = heap
                        .slice(normalized.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    FreeInpAtomData(heap, &mut data)?;
                    inchi_free(heap, normalized)?;
                    pStruct.pOne_norm_data[index] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, at3)?;
            Ok(())
        })();
        return match (result, cleanup) {
            (Err(error), _) => Err(error),
            (Ok(_), Err(error)) => Err(error),
            (Ok(ret), Ok(())) => Ok(ret),
        };
    }

    let mut bfsq = BFS_Q::default();
    let _ = AllocBfsQueue(heap, &mut bfsq, BFS_Q_CLEAR, 0)?;
    let valence_pointer = match inchi_calloc::<VAL_AT>(
        heap,
        u64::try_from(num_at).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
        std::mem::size_of::<VAL_AT>() as u64,
    ) {
        Ok(pointer) => pointer,
        Err(
            SourceHeapError::AllocationFailed
            | SourceHeapError::AllocationSizeOverflow
            | SourceHeapError::AllocationElementCountOutOfRange,
        ) => return Ok(RI_ERR_ALLOC),
        Err(error) => return Err(error),
    };
    pStruct.pVA = valence_pointer;
    let count = usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut valence_atoms = vec![VAL_AT::default(); count];
    groups.total_charge = heap
        .slice(pInChI[0].as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .nTotalCharge;

    let queue_ret = AllocBfsQueue(heap, &mut bfsq, num_at, 0)?;
    if queue_ret < 0 {
        heap.slice_mut(valence_pointer)?[..count].clone_from_slice(&valence_atoms);
        let _ = AllocBfsQueue(heap, &mut bfsq, BFS_Q_FREE, 0)?;
        return Ok(queue_ret);
    }
    let bfs_storage = heap.allocate_model_storage(vec![bfsq.clone()])?;
    pStruct.pbfsq = bfs_storage;

    let restore_mode = heap
        .slice(pStruct.pSrm)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let mut atoms = heap
        .slice(at.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let mut bns_pointer = SourceMutPointer::<BN_STRUCT>::null();
    let mut data_pointer = SourceMutPointer::<BN_DATA>::null();

    let execution = (|| -> Result<i32, SourceHeapError> {
        if pStruct.iMobileH == TAUT_NON as i8 && !pInChI[1].is_null() {
            let second = heap
                .slice(pInChI[1].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if second.nNumberOfAtoms > 1 {
                let ret = FillOutpStructEndpointFromInChI(heap, &second, &mut pStruct.endpoint)?;
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }

        for (index, atom) in atoms.iter().enumerate() {
            let mut row = 0_i32;
            valence_atoms[index].cNumValenceElectrons =
                (get_sp_element_type(i32::from(atom.el_number), &mut row) as i8).wrapping_sub(1);
            valence_atoms[index].cPeriodicRowNumber = row as i8;
            valence_atoms[index].cPeriodicNumber = atom.el_number;
            if is_el_a_metal(i32::from(atom.el_number))? != 0 {
                valence_atoms[index].cMetal = if restore_mode.bStereoRemovesMetalFlag != 0 {
                    i8::from(atom.p_parity == 0 && atom.sb_parity[0] == 0)
                } else {
                    1
                };
            }
            valence_atoms[index].cMinRingSize = if atom.valence == 2 && atom.num_H == 0 {
                is_bond_in_Nmax_memb_ring(
                    heap,
                    at,
                    index as i32,
                    0,
                    bfsq.q,
                    bfsq.nAtomLevel,
                    bfsq.cSource,
                    99,
                )? as i8
            } else {
                0
            };
        }

        let endpoint = if pStruct.endpoint.is_null() {
            None
        } else {
            Some(
                heap.slice(pStruct.endpoint.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec(),
            )
        };
        let mut treat_more_as_metals = 0_i32;
        let mut second_pass = false;
        loop {
            for index in 0..count {
                valence_atoms[index].cInitFreeValences = 0;
                let ret = GetAtomRestoreInfo(
                    pCG,
                    &mut atoms,
                    index as i32,
                    &mut valence_atoms,
                    &restore_mode,
                    i32::from(pStruct.bMobileH),
                    endpoint.as_deref(),
                )?;
                if ret < 0 {
                    return Ok(ret);
                }
                if ret == crate::source_types::TREAT_ATOM_AS_METAL as i32
                    && !second_pass
                    && valence_atoms[index].cMetal == 0
                {
                    valence_atoms[index].cMetal = if restore_mode.bStereoRemovesMetalFlag != 0 {
                        i8::from(atoms[index].p_parity == 0 && atoms[index].sb_parity[0] == 0)
                    } else {
                        1
                    };
                    if valence_atoms[index].cMetal != 0 {
                        treat_more_as_metals = treat_more_as_metals.wrapping_add(1);
                    }
                }
                groups.charge_on_atoms = groups
                    .charge_on_atoms
                    .wrapping_add(i32::from(valence_atoms[index].cInitCharge));
            }
            if treat_more_as_metals == 0 || second_pass {
                break;
            }
            for value in &mut valence_atoms {
                groups.charge_on_atoms = groups
                    .charge_on_atoms
                    .wrapping_sub(i32::from(value.cInitCharge));
                *value = VAL_AT {
                    cMetal: value.cMetal,
                    cMinRingSize: value.cMinRingSize,
                    cNumValenceElectrons: value.cNumValenceElectrons,
                    cPeriodicRowNumber: value.cPeriodicRowNumber,
                    cPeriodicNumber: value.cPeriodicNumber,
                    ..VAL_AT::default()
                };
            }
            second_pass = true;
        }
        heap.slice_mut(at)?[..count].clone_from_slice(&atoms);

        let mut ret = nCountBnsSizes(
            heap,
            &atoms,
            num_at,
            2,
            0,
            &pStruct.ti,
            &valence_atoms,
            &restore_mode,
            &mut groups,
        )?;
        if ret < 0 {
            return Ok(ret);
        }
        ret = nAddSuperCGroups(heap, &mut groups)?;
        if ret < 0 {
            return Ok(ret);
        }

        let mut changed_bonds = 0_i32;
        bns_pointer = AllocateAndInitTCGBnStruct(
            heap,
            pStruct,
            &valence_atoms,
            &groups,
            2,
            2,
            BN_MAX_ALTP as i32,
            &mut changed_bonds,
        )?;
        if bns_pointer.is_null() {
            return Ok(BNS_OUT_OF_RAM);
        }
        let mut bns = heap
            .slice(bns_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        ret = AddTGroups2TCGBnStruct(heap, &mut bns, pStruct, &mut valence_atoms, &mut groups, 2)?;
        if ret < 0 {
            heap.slice_mut(bns_pointer)?[0] = bns;
            return Ok(ret);
        }
        ret = AddCGroups2TCGBnStruct(heap, &mut bns, pStruct, &mut valence_atoms, &mut groups, 2)?;
        if ret < 0 {
            heap.slice_mut(bns_pointer)?[0] = bns;
            return Ok(ret);
        }
        bns.ulTimeOutTime = SourceMutPointer::null();
        bns.ic = ic;

        data_pointer =
            AllocateAndInitBnData(heap, bns.max_vertices.wrapping_add(bns.max_vertices / 2))?;
        if data_pointer.is_null() {
            heap.slice_mut(bns_pointer)?[0] = bns;
            return Ok(BNS_OUT_OF_RAM);
        }
        let mut data = heap
            .slice(data_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let _ = CheckBnsConsistency(
            Some(pStruct),
            Some(&bns),
            Some(&valence_atoms),
            Some(&groups),
            0,
        );
        let restore_bns_result = RunBnsRestore1(
            heap,
            pCG,
            ic,
            ip,
            sd,
            &mut bns,
            &mut data,
            pStruct,
            &mut valence_atoms,
            &mut groups,
            pInChI,
            num_inp,
            bHasSomeFixedH,
            clock_result,
        );
        ret = restore_bns_result?;
        if ret >= 0 {
            ret = CheckBnsConsistency(
                Some(pStruct),
                Some(&bns),
                Some(&valence_atoms),
                Some(&groups),
                1,
            );
            let original = pStruct.at;
            let restored = pStruct.at2;
            pStruct.at = restored;
            let source = heap
                .slice(original.as_const())?
                .get(..atom_count(pStruct)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(restored)?
                .get_mut(..source.len())
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&source);
            let ret2 = CopyBnsToAtom(heap, pStruct, &bns, &valence_atoms, &groups, 1)?;
            pStruct.at2 = pStruct.at;
            pStruct.at = original;
            if ret2 < 0 {
                ret = ret2;
            }
        }
        heap.slice_mut(data_pointer)?[0] = data;
        heap.slice_mut(bns_pointer)?[0] = bns;
        Ok(ret)
    })();

    pStruct.pbfsq = SourceMutPointer::null();
    let cleanup = (|| -> Result<(), SourceHeapError> {
        heap.slice_mut(valence_pointer)?[..count].clone_from_slice(&valence_atoms);
        let _ = AllocBfsQueue(heap, &mut bfsq, BFS_Q_FREE, 0)?;
        data_pointer = DeAllocateBnData(heap, data_pointer)?;
        bns_pointer = DeAllocateBnStruct(heap, bns_pointer)?;
        if !groups.pTCG.is_null() {
            inchi_free(heap, groups.pTCG)?;
            groups.pTCG = SourceMutPointer::null();
        }
        heap.free(bfs_storage)?;
        Ok(())
    })();
    match (execution, cleanup) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(ret), Ok(())) => Ok(ret),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OneInChI2Atom(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: &mut CANON_GLOBALS,
    ip_inp: &INPUT_PARMS,
    sd: &mut STRUCT_DATA,
    szCurHdr: SourceConstPointer<i8>,
    num_inp: i64,
    pStruct: &mut StrFromINChI,
    iComponent: i32,
    iAtNoOffset: i32,
    bHasSomeFixedH: i32,
    pInChI: [SourceMutPointer<INChI>; 2],
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3477 OneInChI2Atom
    // INCHI✔❌: complete source frame follows verbatim.
    /*
    int OneInChI2Atom(INCHI_CLOCK* ic,
        CANON_GLOBALS* pCG,
        ICHICONST INPUT_PARMS* ip_inp,
        STRUCT_DATA* sd,
        const char* szCurHdr,
        long num_inp,
        StrFromINChI* pStruct,
        int iComponent,
        int iAtNoOffset,
        int bHasSomeFixedH,
        INChI* pInChI[])
    {
        int ret;
        INPUT_PARMS* ip, ip_loc;

        ip_loc = *ip_inp;
        ip = &ip_loc;

        sd->pStrErrStruct[0] = '\0';
        ret = RestoreAtomConnectionsSetStereo(pStruct, iComponent, iAtNoOffset, pInChI[0], pInChI[1]);
        if (ret < 0)
        {
            goto exit_function;
        }
        ret = SetStereoBondTypesFrom0DStereo(pStruct, pInChI[0]);
        if (ret < 0)
        {
            goto exit_function;
        }
        ret = ReconcileAllCmlBondParities(pStruct->at, pStruct->num_atoms, 0);
        if (ret < 0)
        {
            goto exit_function;
        }

        /* main InChI restore function */
        ret = RestoreAtomMakeBNS(ic, pCG, ip, sd, pStruct, iComponent, iAtNoOffset, pInChI, szCurHdr, num_inp, bHasSomeFixedH);

    #ifndef COMPILE_ANSI_ONLY
        if ((pStruct->num_inp_actual > 0 ? pStruct->num_inp_actual : num_inp) >= ip->first_struct_number &&
            ((/*ret > 0 &&*/ ip->bDisplayIfRestoreWarnings) && pStruct->pXYZ))
        {
            inchiTime     ulTStart;
            InchiTimeGet(&ulTStart);
            DisplayRestoredComponent(pCG, pStruct, iComponent, iAtNoOffset, pInChI[0], szCurHdr);
            sd->ulStructTime -= InchiTimeElapsed(ic, &ulTStart); /* subtract display time */
        }
    #endif

        if (ret < 0)
        {
            goto exit_function;
        }
        if ((pStruct->num_inp_actual ? pStruct->num_inp_actual : num_inp) >= ip->first_struct_number && ret >= 0)
        {
            /* remove t-group markings and increment zero-order bonds,
               otherwise MakeInChIOutOfStrFromINChI2() woild fail */
               /* --- moved to MakeInChIOutOfStrFromINChI2 ---
               IncrZeroBondsAndClearEndpts(pStruct->at2, pStruct->num_atoms, iComponent+1);
               CopySt2At( pStruct->at2, pStruct->st, pStruct->num_atoms );
               */
               /* include all restored structure features in pStruct->at2 */
               /* make full InChI out of pStruct->at2, pStruct->num_atoms */
               /***************************************************************************************/
               /* !!! pStruct->One_InChI etc. were removed at the exit from NormalizeAndCompare() !!! */
               /***************************************************************************************/
            if (bHasSomeFixedH && pStruct->iInchiRec == INCHI_REC && pStruct->iMobileH == TAUT_YES &&
                !pStruct->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC))
            {
                /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
                ip->nMode |= REQ_MODE_BASIC;
            }

            ret = MakeInChIOutOfStrFromINChI2(ic, pCG, ip, sd, pStruct, iComponent, iAtNoOffset, num_inp);

            if (ret >= 0)
            {
                ;
            }
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
            else
            {
                fprintf(stdout, "\nERROR in MakeInChI-1: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                    szCurHdr ? szCurHdr : "???", iComponent, pStruct->iInchiRec ? 'R' : 'D', pStruct->iMobileH ? 'M' : 'F', ret);
            }
    #endif
    #endif
        }

    exit_function:

        return ret;
    }
    */
    // END INCHI C FUNCTION: OneInChI2Atom
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OneInChI2Atom
    // INCHI✔❌: COMPILE_ANSI_ONLY=1 excludes DisplayRestoredComponent and display timing.
    // INCHI✔❌: TARGET_API_LIB=1 excludes the non-release diagnostic fprintf branch.
    // INCHI✔❌: SourceHeap checked pointer access and CANON_GLOBALS model aliasing add overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: OneInChI2Atom

    let mut ip = ip_inp.clone();
    sd.pStrErrStruct[0] = 0;

    let mut primary = heap
        .slice(pInChI[0].as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let mobile = if pInChI[1].is_null() {
        None
    } else {
        Some(
            heap.slice(pInChI[1].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        )
    };
    let mut ret = RestoreAtomConnectionsSetStereo(
        heap,
        pStruct,
        iComponent,
        iAtNoOffset,
        &mut primary,
        mobile.as_ref(),
    )?;
    heap.slice_mut(pInChI[0])?[0] = primary.clone();
    if ret < 0 {
        return Ok(ret);
    }
    ret = SetStereoBondTypesFrom0DStereo(heap, pStruct, &primary)?;
    if ret < 0 {
        return Ok(ret);
    }
    ret = ReconcileAllCmlBondParities(heap, pStruct.at, pStruct.num_atoms, 0)?;
    if ret < 0 {
        return Ok(ret);
    }

    let restore_result = RestoreAtomMakeBNS(
        heap,
        ic,
        pCG,
        &ip,
        sd,
        pStruct,
        iComponent,
        iAtNoOffset,
        pInChI,
        szCurHdr,
        num_inp,
        bHasSomeFixedH,
        clock_result,
    );
    ret = restore_result?;
    if ret < 0 {
        return Ok(ret);
    }

    let actual_input_number = if pStruct.num_inp_actual != 0 {
        pStruct.num_inp_actual
    } else {
        num_inp
    };
    if actual_input_number >= ip.first_struct_number {
        if bHasSomeFixedH != 0
            && pStruct.iInchiRec == INCHI_REC as i8
            && pStruct.iMobileH == TAUT_YES as i8
            && pStruct.bFixedHExists == 0
            && ip.nMode & u64::from(REQ_MODE_BASIC) == 0
        {
            ip.nMode |= u64::from(REQ_MODE_BASIC);
        }

        let globals = heap.allocate_model_storage(vec![pCG.clone()])?;
        let generated = MakeInChIOutOfStrFromINChI2(
            heap,
            ic,
            globals,
            Some(&ip),
            Some(sd),
            Some(pStruct),
            iComponent,
            iAtNoOffset,
            num_inp,
            clock_result,
        );
        *pCG = heap
            .slice(globals.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        heap.free(globals)?;
        ret = generated?;
    }
    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MakeProtonComponent(
    heap: &mut SourceHeap,
    pStruct: &mut StrFromINChI,
    _iComponent: i32,
    num_prot: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3574 MakeProtonComponent
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int MakeProtonComponent(StrFromINChI* pStruct, int iComponent, int num_prot)
    {
        inp_ATOM* at = NULL;
        int        i;

        if (num_prot <= 0)
        {
            return 0;
        }
        /* allocate */
        pStruct->at = (inp_ATOM*)inchi_calloc(num_prot, sizeof(pStruct->at[0]));
        pStruct->at2 = (inp_ATOM*)inchi_calloc(num_prot, sizeof(pStruct->at2[0]));
        if (!pStruct->at || !pStruct->at2)
        {
            return 0;
        }
        /* create protons */
        at = pStruct->at;
        /* fill out proton atom info */
        for (i = 0; i < num_prot; i++)
        {
            strcpy(at[i].elname, "H");
            at[i].el_number = EL_NUMBER_H;
            at[i].orig_at_number = i + 1;
            /*
            at[i].orig_compt_at_numb = i + 1;
            at[i].component = i + 1;
            */
            at[i].charge = 1;
        }
        memcpy(pStruct->at2, at, num_prot * sizeof(pStruct->at2[0]));
        pStruct->bDeleted = 0;
        pStruct->num_atoms = num_prot;
        pStruct->bMobileH = TAUT_YES;
        pStruct->iMobileH = TAUT_YES;

        return num_prot;
    }

    */
    // END INCHI C FUNCTION: MakeProtonComponent
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MakeProtonComponent
    // INCHI✔️❌: #define inchi_calloc calloc
    // INCHI✔️❌: #define EL_NUMBER_H ((U_CHAR)1)
    // INCHI✔️❌: #define TAUT_YES 1
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // INCHI✔️❌: SourceHeap allocation lookup adds overhead versus direct C pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: MakeProtonComponent

    if num_prot <= 0 {
        return Ok(0);
    }
    let count = u64::try_from(num_prot).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let source_size = std::mem::size_of::<inp_ATOM>() as u64;

    let first_error = match inchi_calloc::<inp_ATOM>(heap, count, source_size) {
        Ok(pointer) => {
            pStruct.at = pointer;
            None
        }
        Err(error) => {
            pStruct.at = SourceMutPointer::null();
            Some(error)
        }
    };
    let second_error = match inchi_calloc::<inp_ATOM>(heap, count, source_size) {
        Ok(pointer) => {
            pStruct.at2 = pointer;
            None
        }
        Err(error) => {
            pStruct.at2 = SourceMutPointer::null();
            Some(error)
        }
    };
    for error in [first_error, second_error].into_iter().flatten() {
        match error {
            SourceHeapError::AllocationFailed => return Ok(0),
            other => return Err(other),
        }
    }

    let count = usize::try_from(num_prot).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    for index in 0..count {
        let atom = heap
            .slice_mut(pStruct.at)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.elname[0] = b'H' as i8;
        atom.elname[1] = 0;
        atom.el_number = 1;
        atom.orig_at_number = (index as u32).wrapping_add(1) as AT_NUMB;
        atom.charge = 1;
    }
    for index in 0..count {
        let atom = heap
            .slice(pStruct.at.as_const())?
            .get(index)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(pStruct.at2)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = atom;
    }
    pStruct.bDeleted = 0;
    pStruct.num_atoms = num_prot;
    pStruct.bMobileH = TAUT_YES as i8;
    pStruct.iMobileH = TAUT_YES as i8;
    Ok(num_prot)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AddRemProtonsInRestrStruct(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: SourceMutPointer<CANON_GLOBALS>,
    ip_inp: &INPUT_PARMS,
    sd: &mut STRUCT_DATA,
    num_inp: i64,
    bHasSomeFixedH: i32,
    pStruct: &mut [StrFromINChI],
    num_components: i32,
    pStructR: Option<&[StrFromINChI]>,
    num_componentsR: i32,
    nProtonsToBeRemovedByNormFromRevrs: &mut NUM_H,
    mut recmet_change_balance: Option<&mut i32>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3615 AddRemProtonsInRestrStruct
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int AddRemProtonsInRestrStruct(INCHI_CLOCK* ic,
    CANON_GLOBALS* pCG,
    ICHICONST INPUT_PARMS* ip_inp,
    STRUCT_DATA* sd, long num_inp,
    int bHasSomeFixedH,
    StrFromINChI* pStruct,
    int num_components,
    StrFromINChI* pStructR,
    int num_componentsR,
    NUM_H* nProtonsToBeRemovedByNormFromRevrs,
    int* recmet_change_balance)
{
    /* on entry and exit, all at[i].num_H do not include isotopic H  and explicit terminal H are connected */
    int  iComp, q, ret = 0;
    int      num_atoms, num_deleted_H, num_tg, num_changed, num_deleted_components; /* djb-rwth: removing redundant variables */
    inp_ATOM* at;
    INPUT_PARMS* ip, ip_loc;
    int      num_prot = *nProtonsToBeRemovedByNormFromRevrs;
    int      delta_recmet_prot, num_prot_prev, bAccumulateChanges = 0, nNumProtAddedByRevrs;
    INChI_Aux* pINChI_Aux;
    INCHI_MODE bNormalizationFlags;
    int        nChargeRevrs, nChargeInChI;

    if (!num_prot)
    {
        return 0;
    }
    delta_recmet_prot = 0;
    num_changed = 0;
    num_deleted_components = 0;
    ip_loc = *ip_inp;
    ip = &ip_loc;
    /*----------------------------------------------------------------------------------
    nLink < 0 && num_componentsR > 0 => This is a Disconnected structure component; it is
                                        same as already processed reconnected one
                                        Do no preicess it

    nLink > 0 && num_componentsR > 0 => This is a Disconnected structure component;
    (should not happen)                 It it is a result of (nLink-1)th Reconeected
                                        component disconnection (NOT IMPLEMENTED YET)

    nLink = 0                        => Process this component. It is either a reconnected
                                        component, or a result of a disconnection (for now)

    nLink > 0 && num_componentsR = 0 => This is a Reconnected component that is same as
                                        a disconnected one that will not be processed.
                                        Process and save charge delta.
    -----------------------------------------------------------------------------------*/

    for (iComp = 0; iComp < num_components && num_prot; iComp++)
    {
        bAccumulateChanges = 0;
        if (pStruct[iComp].nLink < 0 && num_componentsR > 0)
        {
            /* check */
            q = -(pStruct[iComp].nLink + 1);
            if (!pStructR || !num_componentsR || q >= num_componentsR || pStructR[q].nLink != (iComp + 1))
            {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            continue; /* Disconnected structure component has already been processed as a Reconnected one */
        }

        at = pStruct[iComp].at2;
        num_atoms = pStruct[iComp].num_atoms;
        num_deleted_H = pStruct[iComp].num_deleted_H; /* djb-rwth: removing redundant code */
        bAccumulateChanges = (pStruct[iComp].nLink > 0 && !num_componentsR);
        nChargeRevrs = pStruct[iComp].nChargeRevrs;
        nChargeInChI = pStruct[iComp].nChargeInChI;
        num_deleted_components += (0 != pStruct[iComp].bDeleted);
        if (!at || !num_atoms)
        {
            continue;
        }
        /* find whether it is a reconnected structure */
        q = bRevInchiComponentExists(pStruct + iComp, INCHI_REC, TAUT_YES, 0) ? INCHI_REC : INCHI_BAS;
        /*
        q = pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES]->nNumberOfAtoms? INCHI_REC : INCHI_BAS;
        */
        pINChI_Aux = pStruct[iComp].RevInChI.pINChI_Aux[q][0][TAUT_YES]; /* 0 = 1st component in RevInChI */
        /*nNumProtAddedByRevrs = pINChI_Aux->nNumRemovedProtons;*/
        nNumProtAddedByRevrs = -pStruct[iComp].nNumRemovedProtonsByRevrs;
        bNormalizationFlags = pINChI_Aux->bNormalizationFlags;
        num_tg = pINChI_Aux->nNumberOfTGroups;


        /* disconnect all explicit H and add the number of implicit iso H and all explicit terminal H to the number of implicit H */
        if (0 > (ret = DisconnectedConnectedH(at, num_atoms, num_deleted_H)))
        {
            goto exit_function;
        }
        num_prot_prev = num_prot;
        ret = AddRemoveProtonsRestr(at, num_atoms, &num_prot, nNumProtAddedByRevrs,
            bNormalizationFlags, num_tg, nChargeRevrs, nChargeInChI);

        pStruct[iComp].bPostProcessed = ret;
        num_changed += (ret > 0);
        if (ret < 0)
        {
            goto exit_function;
        }
        if (ret > 0)
        {
            /* recalculate InChI; it will reconnect at */
            StrFromINChI* pStruct1 = pStruct + iComp;
            INCHI_MODE    nMode = ip->nMode;
            FreeAllINChIArrays(pStruct1->RevInChI.pINChI,
                pStruct1->RevInChI.pINChI_Aux,
                pStruct1->RevInChI.num_components);

            if (bHasSomeFixedH && pStruct1->iInchiRec == INCHI_REC && pStruct1->iMobileH == TAUT_YES &&
                !pStruct1->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC))
            {
                /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
                ip->nMode |= REQ_MODE_BASIC;
            }
            /* calls ConnectDisconnectedH(...): subtracts number of implicit iso H from implicit H */

            ret = MakeInChIOutOfStrFromINChI2(ic, pCG, ip, sd, pStruct1, 0, 0, num_inp);

            ip->nMode = nMode;
            if (ret < 0)
            {
                goto exit_function;
            }
        }
        else
        {
            /* reconnect disconnected terminal H and subtracts number of implicit iso H from implicit H */
            if (0 > (ret = ConnectDisconnectedH(at, num_atoms, num_deleted_H)))
            {
                goto exit_function;
            }
        }
        if (bAccumulateChanges && recmet_change_balance)
        {
            /* processed Reconnected layer component that is also present in Disconnected layer */
            delta_recmet_prot += num_prot - num_prot_prev;
        }
    }

    iComp = num_components - 1;
    if (!bHasSomeFixedH && num_prot > 0 && 1 == num_deleted_components && iComp >= 0 && pStruct[iComp].bDeleted)
    {
        /* add bare protons to the deleted Mobile-H component; undelete the component */
        num_prot_prev = num_prot;
        if (!MakeProtonComponent(pStruct + iComp, iComp, num_prot))
        {
            goto exit_function;
        }
        else
        {
            /* recalculate InChI; it will reconnect at */
            StrFromINChI* pStruct1 = pStruct + iComp;
            INCHI_MODE    nMode = ip->nMode;
            num_changed++;
            num_prot = 0;
            FreeAllINChIArrays(pStruct1->RevInChI.pINChI,
                pStruct1->RevInChI.pINChI_Aux,
                pStruct1->RevInChI.num_components);

            if (bHasSomeFixedH && pStruct1->iInchiRec == INCHI_REC && pStruct1->iMobileH == TAUT_YES &&
                !pStruct1->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC))
            {
                /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
                ip->nMode |= REQ_MODE_BASIC;
            }
            /* Although MakeInChIOutOfStrFromINChI2() calls ConnectDisconnectedH(...) */
            /* to subtracts number of implicit iso H from implicit H */
            /* this CANNOT have any effect on the deleted H component */

            ret = MakeInChIOutOfStrFromINChI2(ic, pCG, ip, sd, pStruct1, 0, 0, num_inp);

            ip->nMode = nMode;
            if (ret < 0)
            {
                goto exit_function;
            }
            if (bAccumulateChanges && recmet_change_balance)
            {
                /* processed Reconnected layer component that is also present in Disconnected layer */
                delta_recmet_prot += num_prot - num_prot_prev;
            }
        }
    }
    *nProtonsToBeRemovedByNormFromRevrs = num_prot;
    if (recmet_change_balance)
    {
        *recmet_change_balance = delta_recmet_prot;
    }

exit_function:

    return ret < 0 ? ret : num_changed;
}
    */
    // END INCHI C FUNCTION: AddRemProtonsInRestrStruct
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddRemProtonsInRestrStruct
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // INCHI✔️❌: NUM_H is signed 16-bit; INCHI_BAS=0, INCHI_REC=1, TAUT_YES=1.
    // INCHI✔️❌: SourceHeap checked pointer lookup adds overhead versus direct C pointer access.
    // END INCHI ACTIVE MACRO CONFIGURATION: AddRemProtonsInRestrStruct

    let mut ret = 0_i32;
    let mut num_prot = i32::from(*nProtonsToBeRemovedByNormFromRevrs);
    if num_prot == 0 {
        return Ok(0);
    }
    let mut delta_recmet_prot = 0_i32;
    let mut num_changed = 0_i32;
    let mut num_deleted_components = 0_i32;
    let mut ip = ip_inp.clone();
    let mut bAccumulateChanges = 0_i32;
    let mut iComp = 0_i32;

    while iComp < num_components && num_prot != 0 {
        bAccumulateChanges = 0;
        let index = usize::try_from(iComp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let structure = pStruct.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if structure.nLink < 0 && num_componentsR > 0 {
            let q = structure.nLink.wrapping_add(1).wrapping_neg();
            let linked = if let (Some(reconnected), Ok(q_index)) =
                (pStructR, usize::try_from(q))
            {
                q < num_componentsR
                    && reconnected
                        .get(q_index)
                        .is_some_and(|candidate| candidate.nLink == iComp.wrapping_add(1))
            } else {
                false
            };
            if !linked {
                ret = RI_ERR_PROGR;
                return Ok(ret);
            }
            iComp = iComp.wrapping_add(1);
            continue;
        }

        let at = structure.at2;
        let num_atoms = structure.num_atoms;
        let num_deleted_H = structure.num_deleted_H;
        bAccumulateChanges = i32::from(structure.nLink > 0 && num_componentsR == 0);
        let nChargeRevrs = structure.nChargeRevrs;
        let nChargeInChI = structure.nChargeInChI;
        num_deleted_components =
            num_deleted_components.wrapping_add(i32::from(structure.bDeleted != 0));
        if at.is_null() || num_atoms == 0 {
            iComp = iComp.wrapping_add(1);
            continue;
        }

        let q = if bRevInchiComponentExists(
            heap,
            pStruct.get(index),
            INCHI_REC as i32,
            TAUT_YES as i32,
            0,
        )? != 0
        {
            INCHI_REC as usize
        } else {
            INCHI_BAS as usize
        };
        let aux_rows = pStruct[index].RevInChI.pINChI_Aux[q];
        let aux = heap
            .slice(aux_rows.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?[TAUT_YES as usize];
        let aux = heap
            .slice(aux.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let nNumProtAddedByRevrs = pStruct[index]
            .nNumRemovedProtonsByRevrs
            .wrapping_neg();
        let bNormalizationFlags = aux.bNormalizationFlags;
        let num_tg = aux.nNumberOfTGroups;

        ret = {
            let atoms = heap.slice_mut(at)?;
            DisconnectedConnectedH(atoms, num_atoms, num_deleted_H)?
        };
        if ret < 0 {
            return Ok(ret);
        }
        let num_prot_prev = num_prot;
        ret = {
            let atoms = heap.slice_mut(at)?;
            AddRemoveProtonsRestr(
                atoms,
                num_atoms,
                &mut num_prot,
                nNumProtAddedByRevrs,
                bNormalizationFlags,
                num_tg,
                nChargeRevrs,
                nChargeInChI,
            )?
        };
        pStruct[index].bPostProcessed = ret;
        num_changed = num_changed.wrapping_add(i32::from(ret > 0));
        if ret < 0 {
            return Ok(ret);
        }
        if ret > 0 {
            let old_mode = ip.nMode;
            {
                let reversed = &mut pStruct[index].RevInChI;
                FreeAllINChIArrays(
                    heap,
                    &mut reversed.pINChI,
                    &mut reversed.pINChI_Aux,
                    &mut reversed.num_components,
                )?;
            }
            if bHasSomeFixedH != 0
                && pStruct[index].iInchiRec == INCHI_REC as i8
                && pStruct[index].iMobileH == TAUT_YES as i8
                && pStruct[index].bFixedHExists == 0
                && ip.nMode & u64::from(REQ_MODE_BASIC) == 0
            {
                ip.nMode |= u64::from(REQ_MODE_BASIC);
            }
            ret = MakeInChIOutOfStrFromINChI2(
                heap,
                ic,
                pCG,
                Some(&ip),
                Some(sd),
                pStruct.get_mut(index),
                0,
                0,
                num_inp,
                clock_result,
            )?;
            ip.nMode = old_mode;
            if ret < 0 {
                return Ok(ret);
            }
        } else {
            ret = {
                let atoms = heap.slice_mut(at)?;
                ConnectDisconnectedH(atoms, num_atoms, num_deleted_H)?
            };
            if ret < 0 {
                return Ok(ret);
            }
        }
        if bAccumulateChanges != 0 && recmet_change_balance.is_some() {
            delta_recmet_prot =
                delta_recmet_prot.wrapping_add(num_prot.wrapping_sub(num_prot_prev));
        }
        iComp = iComp.wrapping_add(1);
    }

    iComp = num_components.wrapping_sub(1);
    if bHasSomeFixedH == 0
        && num_prot > 0
        && num_deleted_components == 1
        && iComp >= 0
        && pStruct
            .get(usize::try_from(iComp).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .bDeleted
            != 0
    {
        let index = usize::try_from(iComp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let num_prot_prev = num_prot;
        ret = MakeProtonComponent(heap, &mut pStruct[index], iComp, num_prot)?;
        if ret != 0 {
            let old_mode = ip.nMode;
            num_changed = num_changed.wrapping_add(1);
            num_prot = 0;
            {
                let reversed = &mut pStruct[index].RevInChI;
                FreeAllINChIArrays(
                    heap,
                    &mut reversed.pINChI,
                    &mut reversed.pINChI_Aux,
                    &mut reversed.num_components,
                )?;
            }
            if bHasSomeFixedH != 0
                && pStruct[index].iInchiRec == INCHI_REC as i8
                && pStruct[index].iMobileH == TAUT_YES as i8
                && pStruct[index].bFixedHExists == 0
                && ip.nMode & u64::from(REQ_MODE_BASIC) == 0
            {
                ip.nMode |= u64::from(REQ_MODE_BASIC);
            }
            ret = MakeInChIOutOfStrFromINChI2(
                heap,
                ic,
                pCG,
                Some(&ip),
                Some(sd),
                pStruct.get_mut(index),
                0,
                0,
                num_inp,
                clock_result,
            )?;
            ip.nMode = old_mode;
            if ret < 0 {
                return Ok(ret);
            }
            if bAccumulateChanges != 0 && recmet_change_balance.is_some() {
                delta_recmet_prot =
                    delta_recmet_prot.wrapping_add(num_prot.wrapping_sub(num_prot_prev));
            }
        }
    }
    *nProtonsToBeRemovedByNormFromRevrs = num_prot as NUM_H;
    if let Some(balance) = recmet_change_balance.as_mut() {
        **balance = delta_recmet_prot;
    }
    Ok(if ret < 0 { ret } else { num_changed })
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AddRemIsoProtonsInRestrStruct(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: SourceMutPointer<CANON_GLOBALS>,
    ip_inp: &INPUT_PARMS,
    sd: &mut STRUCT_DATA,
    num_inp: i64,
    bHasSomeFixedH: i32,
    pStruct: &mut [StrFromINChI],
    num_components: i32,
    pStructR: Option<&[StrFromINChI]>,
    num_componentsR: i32,
    pProtonBalance: &mut [NUM_H; NUM_H_ISOTOPES as usize],
    mut recmet_change_balance: Option<&mut [NUM_H; NUM_H_ISOTOPES as usize]>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3816 AddRemIsoProtonsInRestrStruct
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int AddRemIsoProtonsInRestrStruct(INCHI_CLOCK* ic,
    CANON_GLOBALS* pCG,
    ICHICONST INPUT_PARMS* ip_inp,
    STRUCT_DATA* sd,
    long num_inp,
    int bHasSomeFixedH,
    StrFromINChI* pStruct,
    int num_components,
    StrFromINChI* pStructR,
    int num_componentsR,
    NUM_H pProtonBalance[],
    NUM_H recmet_change_balance[])
{
    /* on entry and exit, all at[i].num_H do not include isotopic H and explicit terminal H are connected */
    int  iComp, q, k, ret = 0, bNotEmpty;
    int      num_atoms, num_deleted_H, num_tg, num_changed; /* djb-rwth: removing redundant variables */
    inp_ATOM* at;
    NUM_H    num_prot[NUM_H_ISOTOPES], delta_recmet_prot[NUM_H_ISOTOPES], num_prot_prev[NUM_H_ISOTOPES];
    int      bAccumulateChanges;
    INChI_Aux* pINChI_Aux;
    /* djb-rwth: removing redundant variables */
    INPUT_PARMS* ip, ip_loc;

    ip_loc = *ip_inp;
    ip = &ip_loc;

    memcpy(num_prot, pProtonBalance, sizeof(num_prot));
    for (bNotEmpty = 0, k = 0; k < NUM_H_ISOTOPES; k++)
    {
        bNotEmpty |= num_prot[k];
    }
    if (!bNotEmpty)
    {
        return 0;
    }
    memset(delta_recmet_prot, 0, sizeof(delta_recmet_prot)); /* djb-rwth: memset_s C11/Annex K variant? */
    num_changed = 0;
    /*----------------------------------------------------------------------------------
    nLink < 0 && num_componentsR > 0 => This is a Disconnected structure component; it is
                                        same as already processed reconnected one
                                        Do no process it

    nLink > 0 && num_componentsR > 0 => This is a Disconnected structure component;
    (should not happen)                 It it is a result of (nLink-1)th Reconeected
                                        component disconnection (NOT IMPLEMENTED YET)

    nLink = 0                        => Process this component. It is either a reconnected
                                        component, or a result of a disconnection (for now)

    nLink > 0 && num_componentsR = 0 => This is a Reconnected component that is same as
                                        a disconnected one that will not be processed.
                                        Process and save charge delta.
    -----------------------------------------------------------------------------------*/

    for (iComp = 0; iComp < num_components && num_prot; iComp++) /* djb-rwth: the condition will always evaluate to true only if pProtonBalance is not NULL */
    {
        /* djb-rwth: removing redundant code */
        if (pStruct[iComp].nLink < 0 && num_componentsR > 0)
        {
            /* check */
            q = -(pStruct[iComp].nLink + 1);
            if (!pStructR || !num_componentsR || q >= num_componentsR || pStructR[q].nLink != (iComp + 1))
            {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            continue; /* Disconnected structure component has already been processed as a Reconnected one */
        }

        at = pStruct[iComp].at2;
        num_atoms = pStruct[iComp].num_atoms;
        num_deleted_H = pStruct[iComp].num_deleted_H; /* djb-rwth: removing redundant code */
        bAccumulateChanges = (pStruct[iComp].nLink > 0 && !num_componentsR);

        if (!at || !num_atoms)
        {
            continue;
        }
        /* find whether it is a reconnected structure */
        q = pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES]->nNumberOfAtoms ? INCHI_REC : INCHI_BAS;

        pINChI_Aux = pStruct[iComp].RevInChI.pINChI_Aux[q][0][TAUT_YES]; /* 0 = 1st component in RevInChI */
        /* djb-rwth: removing redundant code */
        num_tg = pINChI_Aux->nNumberOfTGroups;
        memcpy(num_prot_prev, num_prot, sizeof(num_prot_prev));

        /* pass CONNECTED explicit H to AddRemoveIsoProtonsRestr() for isotopic H addition */
        ret = AddRemoveIsoProtonsRestr(at, num_atoms, num_prot, num_tg);

        pStruct[iComp].bPostProcessed |= ret;
        num_changed += (ret > 0);
        if (ret < 0)
        {
            goto exit_function;
        }
        if (ret > 0)
        {
            StrFromINChI* pStruct1 = pStruct + iComp;
            INCHI_MODE    nMode = ip->nMode;
            /* recalculate InChI; MakeInChIOutOfStrFromINChI2() will reconnect explicit H */
            /* disconnect all explicit H and add the number of implicit iso H and all explicit terminal H to the number of implicit H */
            if (0 > (ret = DisconnectedConnectedH(at, num_atoms, num_deleted_H)))
            {
                goto exit_function;
            }
            FreeAllINChIArrays(pStruct1->RevInChI.pINChI,
                pStruct1->RevInChI.pINChI_Aux,
                pStruct1->RevInChI.num_components);
            if (bHasSomeFixedH && pStruct1->iInchiRec == INCHI_REC && pStruct1->iMobileH == TAUT_YES &&
                !pStruct1->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC))
            {
                /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
                ip->nMode |= REQ_MODE_BASIC;
            }
            /* input: disconnected explicit H, output: connected explicit H */
            ret = MakeInChIOutOfStrFromINChI2(ic, pCG, ip, sd, pStruct1, 0, 0, num_inp);
            ip->nMode = nMode;
            if (ret < 0)
            {
                goto exit_function;
            }
        }
        /* the following was commented out 2007-08-28 by DT. Reason: it's a bug since H must be already connected */
        /* else {
            if ( 0 > ( ret = ConnectDisconnectedH( at, num_atoms, num_deleted_H ) ) ) {
                goto exit_function;
            }
        } */
        if (bAccumulateChanges)
        {
            /* processed Reconnected layer component that is also present in Disconnected layer */
            for (k = 0; k < NUM_H_ISOTOPES; k++)
            {
                delta_recmet_prot[k] += num_prot[k] - num_prot_prev[k];
            }
        }
    }

    memcpy(pProtonBalance, num_prot, sizeof(num_prot));
    if (recmet_change_balance)
    {
        memcpy(recmet_change_balance, delta_recmet_prot, sizeof(delta_recmet_prot));
    }

exit_function:

    return ret < 0 ? ret : num_changed;
}
    */
    // END INCHI C FUNCTION: AddRemIsoProtonsInRestrStruct
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddRemIsoProtonsInRestrStruct
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this function; COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter its body.
    // INCHI✔️❌: NUM_H is signed 16-bit; NUM_H_ISOTOPES=3; INCHI_BAS=0, INCHI_REC=1, TAUT_YES=1.
    // INCHI✔️❌: SourceHeap checked pointer lookup adds overhead versus direct C pointer access.
    // END INCHI ACTIVE MACRO CONFIGURATION: AddRemIsoProtonsInRestrStruct

    let mut ip = ip_inp.clone();
    let mut num_prot = *pProtonBalance;
    let mut bNotEmpty = 0_i32;
    for value in num_prot {
        bNotEmpty |= i32::from(value);
    }
    if bNotEmpty == 0 {
        return Ok(0);
    }

    let mut delta_recmet_prot = [0 as NUM_H; NUM_H_ISOTOPES as usize];
    let mut num_changed = 0_i32;
    let mut ret = 0_i32;
    let mut iComp = 0_i32;
    while iComp < num_components {
        let index = usize::try_from(iComp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let structure = pStruct.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if structure.nLink < 0 && num_componentsR > 0 {
            let q = structure.nLink.wrapping_add(1).wrapping_neg();
            let linked = if let (Some(reconnected), Ok(q_index)) =
                (pStructR, usize::try_from(q))
            {
                q < num_componentsR
                    && reconnected
                        .get(q_index)
                        .is_some_and(|candidate| candidate.nLink == iComp.wrapping_add(1))
            } else {
                false
            };
            if !linked {
                return Ok(RI_ERR_PROGR);
            }
            iComp = iComp.wrapping_add(1);
            continue;
        }

        let at = structure.at2;
        let num_atoms = structure.num_atoms;
        let num_deleted_H = structure.num_deleted_H;
        let bAccumulateChanges = i32::from(structure.nLink > 0 && num_componentsR == 0);
        if at.is_null() || num_atoms == 0 {
            iComp = iComp.wrapping_add(1);
            continue;
        }

        let rec_rows = pStruct[index].RevInChI.pINChI_Aux[INCHI_REC as usize];
        let has_rec_aux = if rec_rows.is_null() {
            false
        } else {
            let rec_row = heap
                .slice(rec_rows.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let rec_aux = rec_row[TAUT_YES as usize];
            if rec_aux.is_null() {
                false
            } else {
                heap.slice(rec_aux.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNumberOfAtoms
                    != 0
            }
        };
        let q = if has_rec_aux {
            INCHI_REC as usize
        } else {
            INCHI_BAS as usize
        };
        let aux_rows = pStruct[index].RevInChI.pINChI_Aux[q];
        let aux_row = heap
            .slice(aux_rows.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let aux = heap
            .slice(aux_row[TAUT_YES as usize].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let num_tg = aux.nNumberOfTGroups;
        let num_prot_prev = num_prot;

        ret = {
            let atoms = heap.slice_mut(at)?;
            AddRemoveIsoProtonsRestr(atoms, num_atoms, &mut num_prot, num_tg)?
        };
        pStruct[index].bPostProcessed |= ret;
        num_changed = num_changed.wrapping_add(i32::from(ret > 0));
        if ret < 0 {
            return Ok(ret);
        }
        if ret > 0 {
            ret = {
                let atoms = heap.slice_mut(at)?;
                DisconnectedConnectedH(atoms, num_atoms, num_deleted_H)?
            };
            if ret < 0 {
                return Ok(ret);
            }
            {
                let reversed = &mut pStruct[index].RevInChI;
                FreeAllINChIArrays(
                    heap,
                    &mut reversed.pINChI,
                    &mut reversed.pINChI_Aux,
                    &mut reversed.num_components,
                )?;
            }
            let old_mode = ip.nMode;
            if bHasSomeFixedH != 0
                && pStruct[index].iInchiRec == INCHI_REC as i8
                && pStruct[index].iMobileH == TAUT_YES as i8
                && pStruct[index].bFixedHExists == 0
                && ip.nMode & u64::from(REQ_MODE_BASIC) == 0
            {
                ip.nMode |= u64::from(REQ_MODE_BASIC);
            }
            ret = MakeInChIOutOfStrFromINChI2(
                heap,
                ic,
                pCG,
                Some(&ip),
                Some(sd),
                pStruct.get_mut(index),
                0,
                0,
                num_inp,
                clock_result,
            )?;
            ip.nMode = old_mode;
            if ret < 0 {
                return Ok(ret);
            }
        }
        if bAccumulateChanges != 0 {
            for isotope in 0..NUM_H_ISOTOPES as usize {
                delta_recmet_prot[isotope] = delta_recmet_prot[isotope]
                    .wrapping_add(num_prot[isotope].wrapping_sub(num_prot_prev[isotope]));
            }
        }
        iComp = iComp.wrapping_add(1);
    }

    *pProtonBalance = num_prot;
    if let Some(balance) = recmet_change_balance.as_mut() {
        **balance = delta_recmet_prot;
    }
    Ok(if ret < 0 { ret } else { num_changed })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn NormalizeAndCompare(
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
    pInChI: [SourceMutPointer<INChI>; 2],
    _num_inp: i64,
    bHasSomeFixedH: i32,
    pnNumRunBNS: &mut i32,
    pnTotalDelta: &mut i32,
    forbidden_edge_mask: i32,
    forbidden_stereo_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1274 NormalizeAndCompare
    // INCHI❗❌: exact active source frame is inserted below; implementation is under focused verification.
    /*
    int NormalizeAndCompare(CANON_GLOBALS* pCG,
        INCHI_CLOCK* ic,
        ICHICONST INPUT_PARMS* ip,
        STRUCT_DATA* sd,
        BN_STRUCT* pBNS,
        BN_DATA* pBD,
        StrFromINChI* pStruct,
        inp_ATOM* at,
        inp_ATOM* at2,
        inp_ATOM* at3,
        VAL_AT* pVA,
        ALL_TC_GROUPS* pTCGroups,
        INChI* pInChI[],
        long num_inp,
        int bHasSomeFixedH,
        int* pnNumRunBNS,
        int* pnTotalDelta,
        int forbidden_edge_mask,
        int forbidden_stereo_edge_mask)
    {
        int i;
        int err;
        ICR icr, icr2;
        int num_norm_endpoints, num_endpoints, num_norm_t_groups, ret = 0; /* djb-rwth: ignoring LLVM warning: variables used; removing redundant variables */
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
        const char* szCurHdr = (ip->pSdfValue && ip->pSdfValue[0]) ? ip->pSdfValue : "???";
        int         iComponent = pTCGroups->iComponent;
    #endif
    #endif
        T_GROUP_INFO* t_group_info = NULL;
        inp_ATOM* at_norm = NULL; /* normalized */
        inp_ATOM* at_prep = NULL; /* preprocessed */
        INCHI_MODE  cmpInChI, cmpInChI2;
        int         nDeltaPrev, nDeltaCur;
        int         iOrigInChI, iRevrInChI;


        /***********************************************************/
        /* normalize and create one component InChI                */
        /***********************************************************/
        ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
            &t_group_info, &at_norm, &at_prep);
        if (ret < 0)
        {
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
            fprintf(stdout, "\nERROR in MakeOneInchi-1: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                szCurHdr ? szCurHdr : "???", iComponent, pStruct->iInchiRec ? 'R' : 'D', pStruct->iMobileH ? 'M' : 'F', ret);
    #endif
    #endif
            goto exit_function;
        }
        if (pStruct->bMobileH == TAUT_NON)
        {
            /* these indexes are used to compare Mobile-H InChI */
            iOrigInChI = (pInChI[1] && pInChI[1]->nNumberOfAtoms && !pInChI[1]->bDeleted) ? 1 : 0;
            iRevrInChI = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted) ? 1 : 0;
        }
        else
        {
            iOrigInChI = 0;
            iRevrInChI = 0;
        }

        /* Intercept and correct non-polymer Zz to Zy if applicable */
        if (pStruct->n_zy && pStruct->n_pzz)
        {
            if (pStruct->pOneINChI[iRevrInChI]->szHillFormula)
            {
                INCHI_IOS_STRING temp_string_container;
                INCHI_IOS_STRING* strbuf = &temp_string_container;
                int len0 = strlen(pStruct->pOneINChI[iRevrInChI]->szHillFormula);
                if (inchi_strbuf_init(strbuf, len0 + 1, len0 + 1) > 0)
                {
                    inchi_strbuf_printf(strbuf, "%-s", pStruct->pOneINChI[iRevrInChI]->szHillFormula);
                }
                MergeZzInHillFormula(strbuf);
                if (strbuf->nUsedLength > len0 + 1)
                {
                    char* ctmp; /* djb-rwth: supplementary variable */
                    ctmp = (char*)inchi_realloc(pStruct->pOneINChI[iRevrInChI]->szHillFormula, (long long)strbuf->nUsedLength + 1); /* djb-rwth: cast operator added */
                    if (ctmp != NULL) /* djb-rwth: NULL pointer must not be assigned to pStruct->pOneINChI[iRevrInChI]->szHillFormula */
                        pStruct->pOneINChI[iRevrInChI]->szHillFormula = ctmp;
                }
                strcpy(pStruct->pOneINChI[iRevrInChI]->szHillFormula, strbuf->pStr);
                inchi_strbuf_close(strbuf);
            }
        }


        /************************************************************/
        /* compare                                                  */
        /************************************************************/
        if (pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr(pStruct)))
        {
            goto exit_function;
        }
        cmpInChI = CompareReversedINChI2(pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *a2*/, &icr, &err);
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
        /********** InChI from restored structure has LESS hydrogen atoms ******************************/
        if ((cmpInChI & IDIF_LESS_H) && at_prep && 0 < (nDeltaCur = icr.tot_num_H2 - icr.tot_num_H1))
        {
            do
            {
                ret = FixLessHydrogenInFormula(pBNS, pBD, pStruct, at, at2, at_prep, pVA, pTCGroups,
                    pnNumRunBNS, pnTotalDelta, forbidden_edge_mask);
                if (ret < 0)
                {
                    goto exit_function;
                }
                if (ret)
                {
                    /* Probably success. The changes are in pBNS. Create new InChI out of the new restored structure */
                    ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                        &t_group_info, &at_norm, &at_prep);
                    if (ret < 0)
                    {
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
                        fprintf(stdout, "\nERROR in MakeOneInchi-2: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                            szCurHdr ? szCurHdr : "???", iComponent, pStruct->iInchiRec ? 'R' : 'D', pStruct->iMobileH ? 'M' : 'F', ret);
    #endif
    #endif
                        goto exit_function;
                    }
                    /* compare new InChI to the original InChI */
                    if (pStruct->bMobileH == TAUT_NON)
                    {
                        iRevrInChI = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted) ? 1 : 0;
                    }
                    else
                    {
                        iRevrInChI = 0;
                    }
                    if (pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr(pStruct)))
                    {
                        goto exit_function;
                    }
                    cmpInChI = CompareReversedINChI2(pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL, &icr, &err);
                    nDeltaPrev = nDeltaCur;
                    nDeltaCur = icr.tot_num_H2 - icr.tot_num_H1;
                }
                else
                {
                    break;
                }
            } while ((cmpInChI & IDIF_LESS_H) && at_prep && nDeltaCur && nDeltaCur < nDeltaPrev);
        }
        /********** InChI from restored structure has MORE hydrogen atoms ******************************/
        if ((cmpInChI & IDIF_MORE_H) && at_prep && 0 < (nDeltaCur = icr.tot_num_H1 - icr.tot_num_H2))
        {
            do
            {
                ret = FixMoreHydrogenInFormula(pBNS, pBD, pStruct, at, at2, at_prep, pVA, pTCGroups,
                    pnNumRunBNS, pnTotalDelta, forbidden_edge_mask);
                if (ret < 0)
                {
                    goto exit_function;
                }
                if (ret)
                {
                    /* Probably success. The changes are in pBNS. Create new InChI out of the new restored structure */
                    ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                        &t_group_info, &at_norm, &at_prep);
                    if (ret < 0)
                    {
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
                        fprintf(stdout, "\nERROR in MakeOneInchi-3: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                            szCurHdr ? szCurHdr : "???", iComponent, pStruct->iInchiRec ? 'R' : 'D', pStruct->iMobileH ? 'M' : 'F', ret);
    #endif
    #endif
                        goto exit_function;
                    }
                    /* compare new InChI to the original InChI */
                    if (pStruct->bMobileH == TAUT_NON)
                    {
                        iRevrInChI = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted) ? 1 : 0;
                    }
                    else
                    {
                        iRevrInChI = 0;
                    }
                    if (pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr(pStruct)))
                    {
                        goto exit_function;
                    }
                    cmpInChI = CompareReversedINChI2(pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL, &icr, &err);
                    nDeltaPrev = nDeltaCur;
                    nDeltaCur = icr.tot_num_H1 - icr.tot_num_H2;
                }
                else
                {
                    break;
                }
            } while ((cmpInChI & IDIF_MORE_H) && at_prep && nDeltaCur && nDeltaCur < nDeltaPrev);
        }
        /***************** Fix non-taut atoms normalized to tautomeric endpoints ***********************/
        if ((cmpInChI & IDIF_EXTRA_TG_ENDP) && at_norm && 0 < (nDeltaCur = icr.num_endp_in1_only))
        {
            do
            {
                ret = FixRemoveExtraTautEndpoints(pBNS, pBD, pStruct, at, at2, at_prep, at_norm, pVA, pTCGroups, &icr,
                    pnNumRunBNS, pnTotalDelta, forbidden_edge_mask);
                if (ret < 0)
                {
                    goto exit_function;
                }
                if (ret)
                {
                    /* Probably success. The changes are in pBNS. Create new InChI out of the new restored structure */
                    ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                        &t_group_info, &at_norm, &at_prep);
                    if (ret < 0)
                    {
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
                        fprintf(stdout, "\nERROR in MakeOneInchi-4: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                            szCurHdr ? szCurHdr : "???", iComponent, pStruct->iInchiRec ? 'R' : 'D', pStruct->iMobileH ? 'M' : 'F', ret);
    #endif
    #endif
                        goto exit_function;
                    }
                    /* compare new InChI to the original InChI */
                    if (pStruct->bMobileH == TAUT_NON)
                    {
                        iRevrInChI = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted) ? 1 : 0;
                    }
                    else
                    {
                        iRevrInChI = 0;
                    }
                    if (pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr(pStruct)))
                    {
                        goto exit_function;
                    }
                    cmpInChI = CompareReversedINChI2(pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL, &icr, &err);
                    nDeltaPrev = nDeltaCur;
                    nDeltaCur = icr.num_endp_in1_only;
                }
                else
                {
                    break;
                }
            } while ((cmpInChI & IDIF_EXTRA_TG_ENDP) && at_norm && nDeltaCur && nDeltaCur < nDeltaPrev);
        }
        /************************ case of Fixed-H ******************************************************/

        if (pStruct->bMobileH == TAUT_NON)
        {
            int num_tries = 0;
            do
            {
                if (0 > (ret = FixFixedHRestoredStructure(pCG, ic, ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups,
                    &t_group_info, &at_norm, &at_prep, pInChI,
                    num_inp, bHasSomeFixedH, pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
                    forbidden_stereo_edge_mask)))
                {
                    goto exit_function;
                }
            } while (num_tries++ < 2 && ret > 0);
        }
        /************************ case of Fixed-H ******************************************************/
        if (pStruct->bMobileH == TAUT_YES)
        {
            if (0 > (ret = FixMobileHRestoredStructure(pCG, ic, ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups,
                &t_group_info, &at_norm, &at_prep, pInChI,
                num_inp, bHasSomeFixedH, pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
                forbidden_stereo_edge_mask)))
            {
                goto exit_function;
            }
        }
        /**********************************************************************************************/
        /* stereo */
        cmpInChI = CompareReversedINChI2(pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *a2*/, &icr, &err);
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
        memset(&icr2, 0, sizeof(icr2)); /* djb-rwth: memset_s C11/Annex K variant? */
        if (iRevrInChI || iOrigInChI)
        {
            /* additional mobile-H compare in case of Fixed-H */
            cmpInChI2 = CompareReversedINChI2(pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *a2*/, &icr2, &err);
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
        ret = FixRestoredStructureStereo(pCG, ic,
            cmpInChI, &icr, cmpInChI2, &icr2,
            ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups,
            &t_group_info, &at_norm, &at_prep, pInChI,
            num_inp, pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
            forbidden_stereo_edge_mask);

        if (ret < 0)
        {
            goto exit_function;
        }
    #if ( FIX_ADD_PROTON_FOR_ADP == 1 )
        /************************ check and fix ADP by adding a proton (dummy) *************************/
        if (cmpInChI && pTCGroups->num_tgroups && pBNS->tot_st_cap > pBNS->tot_st_flow)
        {
            ret = FixAddProtonForADP(pBNS, pBD, pStruct, at, at2, at_prep, pVA, pTCGroups, &icr,
                pnNumRunBNS, pnTotalDelta, forbidden_edge_mask);
            if (ret < 0)
            {
                goto exit_function;
            }
        }
    #endif
        /* moved to MakeOneInChIOutOfStrFromINChI():
          pStruct->nNumRemovedProtons = (pStruct->iMobileH == TAUT_YES)? pStruct->One_ti.tni.nNumRemovedProtons : 0;
        */

        /* count endpoints */
        num_endpoints = 0;
        num_norm_endpoints = 0;
        num_norm_t_groups = 0;
        /* djb-rwth: removing redundant code */
        at_norm = pStruct->pOne_norm_data[0]->at;
        for (i = 0; i < pTCGroups->num_tgroups; i++)
        {
            num_endpoints += pTCGroups->pTCG[i].num_edges;
            /* djb-rwth: removing redundant code */
        }

        if (t_group_info)
        {
            /* after canonicalization, t_group_info->t_group[i].num[0] = number of H   */
            /*                         t_group_info->t_group[i].num[1] = number of (-) */
            for (i = 0; i < t_group_info->num_t_groups; i++)
            {
                if (t_group_info->t_group[i].num[0])
                {
                    num_norm_t_groups++;
                    num_norm_endpoints += t_group_info->t_group[i].nNumEndpoints;
                    /* djb-rwth: removing redundant code */
                }
            }
        }
    #if ( bRELEASE_VERSION == 0 )
    #ifndef TARGET_API_LIB
        if (num_norm_t_groups != pTCGroups->num_tgroups || num_norm_endpoints != num_endpoints)
        {
            /* need aggressive (de)protonation */
            /* pStruct->bExtract |= EXTRACT_STRUCT_NUMBER; */
            fprintf(stdout, "NORMCOMP: %s comp=%d %c%c: InChI/NormRvrs NumTg=%d/%d NumEndp=%d/%d\n",
                (*ip).pSdfValue, (*pTCGroups).iComponent,
                pStruct->iInchiRec ? 'R' : 'D', pStruct->iMobileH ? 'M' : 'F',
                pTCGroups->num_tgroups, num_norm_t_groups,
                num_endpoints, num_norm_endpoints);
        }
    #endif
    #endif

    exit_function:

        for (i = 0; i < TAUT_NUM; i++)
        {
            Free_INChI(&pStruct->pOneINChI[i]);
            Free_INChI_Aux(&pStruct->pOneINChI_Aux[i]);
            FreeInpAtomData(pStruct->pOne_norm_data[i]);
            if (pStruct->pOne_norm_data[i])
            {
                inchi_free(pStruct->pOne_norm_data[i]);
                pStruct->pOne_norm_data[i] = NULL;
            }
        }

        free_t_group_info(&pStruct->One_ti);

        return ret;
    }
    */
    // END INCHI C FUNCTION: NormalizeAndCompare
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: NormalizeAndCompare
    // INCHI✔❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; bRELEASE_VERSION=1.
    // INCHI✔❌: FIX_ADD_PROTON_FOR_ADP=0; debug fprintf and dummy-proton blocks are inactive.
    // END INCHI ACTIVE MACRO CONFIGURATION: NormalizeAndCompare

    let mut ret = 0_i32;
    let mut t_group_info = SourceTGroupInfoPointer::Null;
    let mut at_norm = SourceMutPointer::null();
    let mut at_prep = SourceMutPointer::null();

    let execution = (|| -> Result<(), SourceHeapError> {
        macro_rules! first_model {
            ($pointer:expr) => {{
                heap.slice($pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            }};
        }
        macro_rules! compare {
            ($reverse:expr, $original:expr, $output:expr) => {{
                #[cfg(test)]
                let forced_compare =
                    normalize_and_compare_test_compare($reverse, $original, $output);
                #[cfg(not(test))]
                let forced_compare: Option<(INCHI_MODE, i32)> = None;
                if let Some(forced_compare) = forced_compare {
                    forced_compare
                } else {
                    let reverse = first_model!(pStruct.pOneINChI[$reverse]);
                    let original = first_model!(pInChI[$original]);
                    let auxiliary = if pStruct.pOneINChI_Aux[$reverse].is_null() {
                        None
                    } else {
                        Some(first_model!(pStruct.pOneINChI_Aux[$reverse]))
                    };
                    let mut comparison_error = 0_i32;
                    let flags = CompareReversedINChI2(
                        heap,
                        Some(&reverse),
                        Some(&original),
                        auxiliary.as_ref(),
                        None,
                        $output,
                        &mut comparison_error,
                    )?;
                    (flags, comparison_error)
                }
            }};
        }
        macro_rules! rebuild {
            () => {{
                #[cfg(test)]
                let forced_rebuild = normalize_and_compare_test_rebuild_return();
                #[cfg(not(test))]
                let forced_rebuild: Option<(i32, bool, bool, bool)> = None;
                ret = if let Some((forced_rebuild, prep_present, norm_present, tinfo_present)) =
                    forced_rebuild
                {
                    at_prep = if prep_present {
                        at
                    } else {
                        SourceMutPointer::null()
                    };
                    at_norm = if norm_present {
                        at
                    } else {
                        SourceMutPointer::null()
                    };
                    t_group_info = if tinfo_present {
                        SourceTGroupInfoPointer::StructureOne
                    } else {
                        SourceTGroupInfoPointer::Null
                    };
                    forced_rebuild
                } else {
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
                        Some(&mut t_group_info),
                        Some(&mut at_norm),
                        Some(&mut at_prep),
                        clock_result,
                    );
                    rebuild_result?
                };
                if ret < 0 {
                    return Ok(());
                }
            }};
        }
        macro_rules! fill_fixed_h {
            () => {{
                #[cfg(test)]
                let forced = normalize_and_compare_test_fill();
                #[cfg(not(test))]
                let forced: Option<i32> = None;
                if let Some(forced) = forced {
                    forced
                } else {
                    FillOutExtraFixedHDataRestr(heap, pStruct)?
                }
            }};
        }
        macro_rules! fix_less_h {
            () => {{
                #[cfg(test)]
                let forced = normalize_and_compare_test_fix_less();
                #[cfg(not(test))]
                let forced: Option<i32> = None;
                if let Some(forced) = forced {
                    forced
                } else {
                    FixLessHydrogenInFormula(
                        heap,
                        pBNS,
                        pBD,
                        pStruct,
                        at,
                        at2,
                        at_prep,
                        pVA,
                        pTCGroups,
                        pnNumRunBNS,
                        pnTotalDelta,
                        forbidden_edge_mask,
                        clock_result,
                    )?
                }
            }};
        }
        macro_rules! fix_more_h {
            () => {{
                #[cfg(test)]
                let forced = normalize_and_compare_test_fix_more();
                #[cfg(not(test))]
                let forced: Option<i32> = None;
                if let Some(forced) = forced {
                    forced
                } else {
                    FixMoreHydrogenInFormula(
                        heap,
                        pBNS,
                        pBD,
                        pStruct,
                        at,
                        at2,
                        at_prep,
                        pVA,
                        pTCGroups,
                        pnNumRunBNS,
                        pnTotalDelta,
                        forbidden_edge_mask,
                        clock_result,
                    )?
                }
            }};
        }
        macro_rules! fix_extra_endpoints {
            ($comparison:expr) => {{
                #[cfg(test)]
                let forced = normalize_and_compare_test_fix_extra();
                #[cfg(not(test))]
                let forced: Option<i32> = None;
                if let Some(forced) = forced {
                    forced
                } else {
                    FixRemoveExtraTautEndpoints(
                        heap,
                        pBNS,
                        pBD,
                        pStruct,
                        at,
                        at2,
                        at_prep,
                        at_norm,
                        pVA,
                        pTCGroups,
                        $comparison,
                        pnNumRunBNS,
                        pnTotalDelta,
                        forbidden_edge_mask,
                        clock_result,
                    )?
                }
            }};
        }

        rebuild!();

        let (iOrigInChI, mut iRevrInChI) = if i32::from(pStruct.bMobileH) == TAUT_NON as i32 {
            let original = if !pInChI[1].is_null() {
                let inchi = first_model!(pInChI[1]);
                usize::from(inchi.nNumberOfAtoms != 0 && inchi.bDeleted == 0)
            } else {
                0
            };
            let reversed = if !pStruct.pOneINChI[1].is_null() {
                let inchi = first_model!(pStruct.pOneINChI[1]);
                usize::from(inchi.nNumberOfAtoms != 0 && inchi.bDeleted == 0)
            } else {
                0
            };
            (original, reversed)
        } else {
            (0, 0)
        };

        if pStruct.n_zy != 0 && pStruct.n_pzz != 0 {
            let reversed_pointer = pStruct.pOneINChI[iRevrInChI];
            let formula_pointer = first_model!(reversed_pointer).szHillFormula;
            if !formula_pointer.is_null() {
                let formula = heap.slice(formula_pointer.as_const())?;
                let len0 = formula
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                let len0_i32 =
                    i32::try_from(len0).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let initial_size = len0_i32.wrapping_add(1);
                let mut buffer = INCHI_IOS_STRING::default();
                let string_result = (|| -> Result<(), SourceHeapError> {
                    #[cfg(test)]
                    normalize_and_compare_test_event("inchi_strbuf_init:begin".to_owned());
                    let init_result =
                        inchi_strbuf_init(heap, &mut buffer, initial_size, initial_size)?;
                    #[cfg(test)]
                    normalize_and_compare_test_event(if init_result > 0 {
                        "inchi_strbuf_init:positive".to_owned()
                    } else {
                        "inchi_strbuf_init:non-positive".to_owned()
                    });
                    if init_result > 0 {
                        let source = heap.slice(formula_pointer.as_const())?[..=len0].to_vec();
                        heap.slice_mut(buffer.pStr)?[..=len0].copy_from_slice(&source);
                        buffer.nUsedLength = len0_i32;
                    }
                    #[cfg(test)]
                    normalize_and_compare_test_event("MergeZzInHillFormula:begin".to_owned());
                    let merge_result = MergeZzInHillFormula(heap, &mut buffer)?;
                    #[cfg(test)]
                    {
                        if normalize_and_compare_test_force_zz_growth() && merge_result == 0 {
                            buffer.nUsedLength = buffer.nAllocatedLength.wrapping_add(1);
                            normalize_and_compare_test_event(
                                "MergeZzInHillFormula:forced-growth".to_owned(),
                            );
                        }
                        normalize_and_compare_test_event(if merge_result == 0 {
                            "MergeZzInHillFormula:zero".to_owned()
                        } else {
                            "MergeZzInHillFormula:negative".to_owned()
                        });
                    }
                    if buffer.pStr.is_null() {
                        return Err(SourceHeapError::AllocationFailed);
                    }
                    if buffer.nUsedLength > len0_i32.wrapping_add(1) {
                        let requested = u64::try_from(buffer.nUsedLength.wrapping_add(1))
                            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
                        let replacement = inchi_realloc(heap, formula_pointer, requested)?;
                        #[cfg(test)]
                        normalize_and_compare_test_event("inchi_realloc:begin".to_owned());
                        #[cfg(test)]
                        normalize_and_compare_test_event(if replacement.is_null() {
                            "inchi_realloc:null".to_owned()
                        } else {
                            "inchi_realloc:non-null".to_owned()
                        });
                        if !replacement.is_null() {
                            heap.slice_mut(reversed_pointer)?[0].szHillFormula = replacement;
                        }
                    }
                    let destination = first_model!(reversed_pointer).szHillFormula;
                    let source_buffer = heap.slice(buffer.pStr.as_const())?;
                    let nul = source_buffer
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or(SourceHeapError::MissingNulTerminator)?;
                    let source = source_buffer[..=nul].to_vec();
                    heap.slice_mut(destination)?
                        .get_mut(..source.len())
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .copy_from_slice(&source);
                    Ok(())
                })();
                #[cfg(test)]
                normalize_and_compare_test_event("inchi_strbuf_close:begin".to_owned());
                let close_result = inchi_strbuf_close(heap, Some(&mut buffer));
                #[cfg(test)]
                normalize_and_compare_test_event("inchi_strbuf_close:end".to_owned());
                string_result?;
                close_result?;
            }
        }

        if i32::from(pStruct.iMobileH) == TAUT_NON as i32 {
            ret = fill_fixed_h!();
            if ret != 0 {
                return Ok(());
            }
        }
        let mut icr = ICR::default();
        let compared = compare!(iRevrInChI, iOrigInChI, &mut icr);
        let mut cmpInChI = compared.0;
        if cmpInChI & IDIF_PROBLEM != 0 {
            ret = RI_ERR_PROGR;
            return Ok(());
        }
        if compared.1 != 0 {
            ret = RI_ERR_ALLOC;
            return Ok(());
        }

        let mut nDeltaCur = 0_i32;
        if cmpInChI & IDIF_LESS_H != 0 && !at_prep.is_null() {
            nDeltaCur = icr.tot_num_H2.wrapping_sub(icr.tot_num_H1);
            if nDeltaCur > 0 {
                loop {
                    ret = fix_less_h!();
                    if ret < 0 {
                        return Ok(());
                    }
                    if ret == 0 {
                        break;
                    }
                    rebuild!();
                    iRevrInChI = if i32::from(pStruct.bMobileH) == TAUT_NON as i32
                        && !pStruct.pOneINChI[1].is_null()
                    {
                        let inchi = first_model!(pStruct.pOneINChI[1]);
                        usize::from(inchi.nNumberOfAtoms != 0 && inchi.bDeleted == 0)
                    } else {
                        0
                    };
                    if i32::from(pStruct.iMobileH) == TAUT_NON as i32 {
                        ret = fill_fixed_h!();
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    let compared = compare!(iRevrInChI, iOrigInChI, &mut icr);
                    cmpInChI = compared.0;
                    let nDeltaPrev = nDeltaCur;
                    nDeltaCur = icr.tot_num_H2.wrapping_sub(icr.tot_num_H1);
                    if cmpInChI & IDIF_LESS_H == 0
                        || at_prep.is_null()
                        || nDeltaCur == 0
                        || nDeltaCur >= nDeltaPrev
                    {
                        break;
                    }
                }
            }
        }

        if cmpInChI & IDIF_MORE_H != 0 && !at_prep.is_null() {
            nDeltaCur = icr.tot_num_H1.wrapping_sub(icr.tot_num_H2);
            if nDeltaCur > 0 {
                loop {
                    ret = fix_more_h!();
                    if ret < 0 {
                        return Ok(());
                    }
                    if ret == 0 {
                        break;
                    }
                    rebuild!();
                    iRevrInChI = if i32::from(pStruct.bMobileH) == TAUT_NON as i32
                        && !pStruct.pOneINChI[1].is_null()
                    {
                        let inchi = first_model!(pStruct.pOneINChI[1]);
                        usize::from(inchi.nNumberOfAtoms != 0 && inchi.bDeleted == 0)
                    } else {
                        0
                    };
                    if i32::from(pStruct.iMobileH) == TAUT_NON as i32 {
                        ret = fill_fixed_h!();
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    let compared = compare!(iRevrInChI, iOrigInChI, &mut icr);
                    cmpInChI = compared.0;
                    let nDeltaPrev = nDeltaCur;
                    nDeltaCur = icr.tot_num_H1.wrapping_sub(icr.tot_num_H2);
                    if cmpInChI & IDIF_MORE_H == 0
                        || at_prep.is_null()
                        || nDeltaCur == 0
                        || nDeltaCur >= nDeltaPrev
                    {
                        break;
                    }
                }
            }
        }

        if cmpInChI & IDIF_EXTRA_TG_ENDP != 0 && !at_norm.is_null() {
            nDeltaCur = icr.num_endp_in1_only;
            if nDeltaCur > 0 {
                loop {
                    ret = fix_extra_endpoints!(&icr);
                    if ret < 0 {
                        return Ok(());
                    }
                    if ret == 0 {
                        break;
                    }
                    rebuild!();
                    iRevrInChI = if i32::from(pStruct.bMobileH) == TAUT_NON as i32
                        && !pStruct.pOneINChI[1].is_null()
                    {
                        let inchi = first_model!(pStruct.pOneINChI[1]);
                        usize::from(inchi.nNumberOfAtoms != 0 && inchi.bDeleted == 0)
                    } else {
                        0
                    };
                    if i32::from(pStruct.iMobileH) == TAUT_NON as i32 {
                        ret = fill_fixed_h!();
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    let compared = compare!(iRevrInChI, iOrigInChI, &mut icr);
                    cmpInChI = compared.0;
                    let nDeltaPrev = nDeltaCur;
                    nDeltaCur = icr.num_endp_in1_only;
                    if cmpInChI & IDIF_EXTRA_TG_ENDP == 0
                        || at_norm.is_null()
                        || nDeltaCur == 0
                        || nDeltaCur >= nDeltaPrev
                    {
                        break;
                    }
                }
            }
        }

        if i32::from(pStruct.bMobileH) == TAUT_NON as i32 {
            let mut num_tries = 0_i32;
            loop {
                #[cfg(test)]
                let forced = normalize_and_compare_test_fix_fixed(_num_inp);
                #[cfg(not(test))]
                let forced: Option<i32> = None;
                ret = if let Some(forced) = forced {
                    forced
                } else {
                    FixFixedHRestoredStructure(
                        heap,
                        pCG,
                        ic,
                        ip,
                        sd,
                        pBNS,
                        pBD,
                        pStruct,
                        at,
                        at2,
                        at3,
                        pVA,
                        pTCGroups,
                        Some(&mut t_group_info),
                        Some(&mut at_norm),
                        Some(&mut at_prep),
                        pInChI,
                        _num_inp,
                        bHasSomeFixedH,
                        Some(pnNumRunBNS),
                        Some(pnTotalDelta),
                        forbidden_edge_mask,
                        forbidden_stereo_edge_mask,
                        clock_result,
                    )?
                };
                if ret < 0 {
                    return Ok(());
                }
                let repeat = num_tries < 2 && ret > 0;
                num_tries = num_tries.wrapping_add(1);
                if !repeat {
                    break;
                }
            }
        }
        if i32::from(pStruct.bMobileH) == TAUT_YES as i32 {
            #[cfg(test)]
            let forced = normalize_and_compare_test_fix_mobile(_num_inp);
            #[cfg(not(test))]
            let forced: Option<i32> = None;
            ret = if let Some(forced) = forced {
                forced
            } else {
                let fix_result = FixMobileHRestoredStructure(
                    heap,
                    pCG,
                    ic,
                    ip,
                    sd,
                    pBNS,
                    pBD,
                    pStruct,
                    at,
                    at2,
                    at3,
                    pVA,
                    pTCGroups,
                    Some(&mut t_group_info),
                    Some(&mut at_norm),
                    Some(&mut at_prep),
                    pInChI,
                    _num_inp,
                    bHasSomeFixedH,
                    Some(pnNumRunBNS),
                    Some(pnTotalDelta),
                    forbidden_edge_mask,
                    forbidden_stereo_edge_mask,
                    clock_result,
                );
                fix_result?
            };
            if ret < 0 {
                return Ok(());
            }
        }

        let compared = compare!(0, 0, &mut icr);
        let cmpInChI = compared.0;
        if cmpInChI & IDIF_PROBLEM != 0 {
            ret = RI_ERR_PROGR;
            return Ok(());
        }
        if compared.1 != 0 {
            ret = RI_ERR_ALLOC;
            return Ok(());
        }
        let mut icr2 = ICR::default();
        let mut cmpInChI2 = 0;
        if iRevrInChI != 0 || iOrigInChI != 0 {
            let compared = compare!(iRevrInChI, iOrigInChI, &mut icr2);
            cmpInChI2 = compared.0;
            if cmpInChI & IDIF_PROBLEM != 0 {
                ret = RI_ERR_PROGR;
                return Ok(());
            }
            if compared.1 != 0 {
                ret = RI_ERR_ALLOC;
                return Ok(());
            }
        }
        #[cfg(test)]
        let forced = normalize_and_compare_test_fix_stereo();
        #[cfg(not(test))]
        let forced: Option<i32> = None;
        ret = if let Some(forced) = forced {
            forced
        } else {
            FixRestoredStructureStereo(
                heap,
                pCG,
                ic,
                cmpInChI,
                &mut icr,
                cmpInChI2,
                &mut icr2,
                ip,
                sd,
                pBNS,
                pBD,
                pStruct,
                at,
                at2,
                at3,
                pVA,
                pTCGroups,
                Some(&mut t_group_info),
                Some(&mut at_norm),
                Some(&mut at_prep),
                pInChI,
                _num_inp,
                pnNumRunBNS,
                pnTotalDelta,
                forbidden_edge_mask,
                forbidden_stereo_edge_mask,
                clock_result,
            )?
        };
        if ret < 0 {
            return Ok(());
        }

        at_norm = first_model!(pStruct.pOne_norm_data[0]).at;
        let mut num_endpoints = 0_i32;
        let group_count = usize::try_from(pTCGroups.num_tgroups)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if group_count != 0 {
            for group in heap
                .slice(pTCGroups.pTCG.as_const())?
                .get(..group_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                num_endpoints = num_endpoints.wrapping_add(group.num_edges);
            }
        }
        let mut num_norm_endpoints = 0_i32;
        let mut num_norm_t_groups = 0_i32;
        if t_group_info == SourceTGroupInfoPointer::StructureOne {
            let count = usize::try_from(pStruct.One_ti.num_t_groups)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            if count != 0 {
                for group in heap
                    .slice(pStruct.One_ti.t_group.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    if group.num[0] != 0 {
                        num_norm_t_groups = num_norm_t_groups.wrapping_add(1);
                        num_norm_endpoints =
                            num_norm_endpoints.wrapping_add(i32::from(group.nNumEndpoints));
                    }
                }
            }
        }
        let _ = (
            at_norm,
            num_endpoints,
            num_norm_endpoints,
            num_norm_t_groups,
        );
        Ok(())
    })();

    let cleanup = (|| -> Result<(), SourceHeapError> {
        for index in 0..TAUT_NUM as usize {
            #[cfg(test)]
            normalize_and_compare_test_cleanup_inchi(heap, index, pStruct.pOneINChI[index]);
            let had_inchi = !pStruct.pOneINChI[index].is_null();
            let _ = Free_INChI(heap, &mut pStruct.pOneINChI[index])?;
            #[cfg(test)]
            if had_inchi {
                normalize_and_compare_test_event(format!("free:pOneINChI[{index}]"));
            }
            #[cfg(test)]
            normalize_and_compare_test_cleanup_aux(heap, index, pStruct.pOneINChI_Aux[index]);
            let had_aux = !pStruct.pOneINChI_Aux[index].is_null();
            let _ = Free_INChI_Aux(heap, &mut pStruct.pOneINChI_Aux[index])?;
            #[cfg(test)]
            if had_aux {
                normalize_and_compare_test_event(format!("free:pOneINChI_Aux[{index}]"));
            }
            let holder = pStruct.pOne_norm_data[index];
            #[cfg(test)]
            normalize_and_compare_test_cleanup_norm(heap, index, holder);
            if !holder.is_null() {
                let mut data = heap
                    .slice(holder.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                #[cfg(test)]
                {
                    if !data.at.is_null() {
                        normalize_and_compare_test_event(format!(
                            "free:pOne_norm_data[{index}].at"
                        ));
                    }
                    if !data.at_fixed_bonds.is_null() {
                        normalize_and_compare_test_event(format!(
                            "free:pOne_norm_data[{index}].at_fixed_bonds"
                        ));
                    }
                }
                FreeInpAtomData(heap, &mut data)?;
                inchi_free(heap, holder)?;
                #[cfg(test)]
                normalize_and_compare_test_event(format!("free:pOne_norm_data[{index}]"));
                pStruct.pOne_norm_data[index] = SourceMutPointer::null();
            }
        }
        #[cfg(test)]
        normalize_and_compare_test_cleanup_t_group(heap, &pStruct.One_ti);
        let had_t_group = !pStruct.One_ti.t_group.is_null();
        let _ = free_t_group_info(heap, Some(&mut pStruct.One_ti))?;
        #[cfg(test)]
        if had_t_group {
            normalize_and_compare_test_event("free:One_ti.t_group".to_owned());
        }
        Ok(())
    })();

    execution?;
    cleanup?;
    Ok(ret)
}

#[cfg(test)]
#[derive(Default)]
struct NormalizeAndCompareTestControl {
    forced_rebuild_return: Option<i32>,
    forced_compare: Option<(INCHI_MODE, i32)>,
    expected_inchi_atoms: [i32; 2],
    holder_mask: u8,
    events: Vec<String>,
    prefree_state_exact: bool,
    force_zz_growth: bool,
    formula_before_cleanup: String,
    rebuild_script: Vec<(i32, bool, bool, bool)>,
    compare_script: Vec<(INCHI_MODE, i32, i32, i32, i32)>,
    fill_script: Vec<i32>,
    fix_less_script: Vec<i32>,
    fix_more_script: Vec<i32>,
    fix_extra_script: Vec<i32>,
    fix_fixed_script: Vec<i32>,
    fix_mobile_script: Vec<i32>,
    fix_stereo_script: Vec<i32>,
    final_cleanup: Option<NormalizeAndCompareFinalCleanup>,
}

#[cfg(test)]
#[derive(Default)]
struct NormalizeAndCompareFinalCleanup {
    expected_norm_at: [SourceMutPointer<inp_ATOM>; TAUT_NUM as usize],
    expected_t_groups: Option<Vec<crate::source_types::T_GROUP>>,
}

#[cfg(test)]
std::thread_local! {
    static NORMALIZE_AND_COMPARE_TEST_CONTROL: std::cell::RefCell<Option<NormalizeAndCompareTestControl>> = const { std::cell::RefCell::new(None) };
}

#[cfg(test)]
fn normalize_and_compare_test_begin(holder_mask: u8, forced_rebuild_return: i32) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        *control.borrow_mut() = Some(NormalizeAndCompareTestControl {
            forced_rebuild_return: Some(forced_rebuild_return),
            forced_compare: None,
            expected_inchi_atoms: [10, 11],
            holder_mask,
            events: Vec::new(),
            prefree_state_exact: true,
            force_zz_growth: false,
            formula_before_cleanup: String::new(),
            rebuild_script: Vec::new(),
            compare_script: Vec::new(),
            fill_script: Vec::new(),
            fix_less_script: Vec::new(),
            fix_more_script: Vec::new(),
            fix_extra_script: Vec::new(),
            fix_fixed_script: Vec::new(),
            fix_mobile_script: Vec::new(),
            fix_stereo_script: Vec::new(),
            final_cleanup: None,
        });
    });
}

#[cfg(test)]
fn normalize_and_compare_test_begin_layer_selection(reversed_state: i32) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        *control.borrow_mut() = Some(NormalizeAndCompareTestControl {
            forced_rebuild_return: Some(0),
            forced_compare: Some((IDIF_PROBLEM, 0)),
            expected_inchi_atoms: [1, if reversed_state == 1 { 0 } else { 1 }],
            holder_mask: 0,
            events: Vec::new(),
            prefree_state_exact: true,
            force_zz_growth: false,
            formula_before_cleanup: String::new(),
            rebuild_script: Vec::new(),
            compare_script: Vec::new(),
            fill_script: Vec::new(),
            fix_less_script: Vec::new(),
            fix_more_script: Vec::new(),
            fix_extra_script: Vec::new(),
            fix_fixed_script: Vec::new(),
            fix_mobile_script: Vec::new(),
            fix_stereo_script: Vec::new(),
            final_cleanup: None,
        });
    });
}

#[cfg(test)]
fn normalize_and_compare_test_begin_common_success() {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        *control.borrow_mut() = Some(NormalizeAndCompareTestControl {
            forced_rebuild_return: Some(0),
            forced_compare: None,
            expected_inchi_atoms: [1, 0],
            holder_mask: 0,
            events: Vec::new(),
            prefree_state_exact: true,
            force_zz_growth: false,
            formula_before_cleanup: String::new(),
            rebuild_script: Vec::new(),
            compare_script: Vec::new(),
            fill_script: Vec::new(),
            fix_less_script: Vec::new(),
            fix_more_script: Vec::new(),
            fix_extra_script: Vec::new(),
            fix_fixed_script: Vec::new(),
            fix_mobile_script: Vec::new(),
            fix_stereo_script: Vec::new(),
            final_cleanup: None,
        });
    });
}

#[cfg(test)]
fn normalize_and_compare_test_begin_zz(force_growth: bool) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        *control.borrow_mut() = Some(NormalizeAndCompareTestControl {
            forced_rebuild_return: Some(0),
            forced_compare: Some((IDIF_PROBLEM, 0)),
            expected_inchi_atoms: [1, 0],
            holder_mask: 0,
            events: Vec::new(),
            prefree_state_exact: true,
            force_zz_growth: force_growth,
            formula_before_cleanup: String::new(),
            rebuild_script: Vec::new(),
            compare_script: Vec::new(),
            fill_script: Vec::new(),
            fix_less_script: Vec::new(),
            fix_more_script: Vec::new(),
            fix_extra_script: Vec::new(),
            fix_fixed_script: Vec::new(),
            fix_mobile_script: Vec::new(),
            fix_stereo_script: Vec::new(),
            final_cleanup: None,
        });
    });
}

#[cfg(test)]
fn normalize_and_compare_test_begin_scripted(
    rebuild_script: Vec<(i32, bool, bool, bool)>,
    compare_script: Vec<(INCHI_MODE, i32, i32, i32, i32)>,
    fill_script: Vec<i32>,
    fix_less_script: Vec<i32>,
    fix_more_script: Vec<i32>,
    fix_extra_script: Vec<i32>,
    fix_fixed_script: Vec<i32>,
    fix_mobile_script: Vec<i32>,
    fix_stereo_script: Vec<i32>,
) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        *control.borrow_mut() = Some(NormalizeAndCompareTestControl {
            forced_rebuild_return: None,
            forced_compare: None,
            expected_inchi_atoms: [1, 0],
            holder_mask: 0,
            events: Vec::new(),
            prefree_state_exact: true,
            force_zz_growth: false,
            formula_before_cleanup: String::new(),
            rebuild_script,
            compare_script,
            fill_script,
            fix_less_script,
            fix_more_script,
            fix_extra_script,
            fix_fixed_script,
            fix_mobile_script,
            fix_stereo_script,
            final_cleanup: None,
        });
    });
}

#[cfg(test)]
fn normalize_and_compare_test_set_final_cleanup(
    expected_inchi_atoms: [i32; TAUT_NUM as usize],
    expected_norm_at: [SourceMutPointer<inp_ATOM>; TAUT_NUM as usize],
    expected_t_groups: Option<Vec<crate::source_types::T_GROUP>>,
) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control
            .as_mut()
            .expect("NormalizeAndCompare test control must be active");
        control.expected_inchi_atoms = expected_inchi_atoms;
        control.final_cleanup = Some(NormalizeAndCompareFinalCleanup {
            expected_norm_at,
            expected_t_groups,
        });
    });
}

#[cfg(test)]
fn normalize_and_compare_test_force_zz_growth() -> bool {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        control
            .borrow()
            .as_ref()
            .is_some_and(|control| control.force_zz_growth)
    })
}

#[cfg(test)]
fn normalize_and_compare_test_rebuild_return() -> Option<(i32, bool, bool, bool)> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        control
            .events
            .push("MakeOneInChIOutOfStrFromINChI2".to_owned());
        if !control.rebuild_script.is_empty() {
            Some(control.rebuild_script.remove(0))
        } else {
            control
                .forced_rebuild_return
                .take()
                .map(|value| (value, false, false, false))
        }
    })
}

#[cfg(test)]
fn normalize_and_compare_test_compare(
    reverse_index: usize,
    original_index: usize,
    output: &mut ICR,
) -> Option<(INCHI_MODE, i32)> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        control.events.push(format!(
            "CompareReversedINChI2:r{reverse_index}:o{original_index}"
        ));
        if !control.compare_script.is_empty() {
            let (flags, error, h1, h2, endpoints) = control.compare_script.remove(0);
            output.tot_num_H1 = h1;
            output.tot_num_H2 = h2;
            output.num_endp_in1_only = endpoints;
            Some((flags, error))
        } else {
            control.forced_compare
        }
    })
}

#[cfg(test)]
fn normalize_and_compare_test_fill() -> Option<i32> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        if control.fill_script.is_empty() {
            return None;
        }
        control
            .events
            .push("FillOutExtraFixedHDataRestr".to_owned());
        Some(control.fill_script.remove(0))
    })
}

#[cfg(test)]
fn normalize_and_compare_test_fix_less() -> Option<i32> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        if control.fix_less_script.is_empty() {
            return None;
        }
        control.events.push("FixLessHydrogenInFormula".to_owned());
        Some(control.fix_less_script.remove(0))
    })
}

#[cfg(test)]
fn normalize_and_compare_test_fix_more() -> Option<i32> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        if control.fix_more_script.is_empty() {
            return None;
        }
        control.events.push("FixMoreHydrogenInFormula".to_owned());
        Some(control.fix_more_script.remove(0))
    })
}

#[cfg(test)]
fn normalize_and_compare_test_fix_extra() -> Option<i32> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        if control.fix_extra_script.is_empty() {
            return None;
        }
        control
            .events
            .push("FixRemoveExtraTautEndpoints".to_owned());
        Some(control.fix_extra_script.remove(0))
    })
}

#[cfg(test)]
fn normalize_and_compare_test_fix_fixed(num_inp: i64) -> Option<i32> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        if control.fix_fixed_script.is_empty() {
            return None;
        }
        control
            .events
            .push(format!("FixFixedHRestoredStructure:num_inp={num_inp}"));
        Some(control.fix_fixed_script.remove(0))
    })
}

#[cfg(test)]
fn normalize_and_compare_test_fix_mobile(num_inp: i64) -> Option<i32> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        if control.fix_mobile_script.is_empty() {
            return None;
        }
        control
            .events
            .push(format!("FixMobileHRestoredStructure:num_inp={num_inp}"));
        Some(control.fix_mobile_script.remove(0))
    })
}

#[cfg(test)]
fn normalize_and_compare_test_fix_stereo() -> Option<i32> {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let control = control.as_mut()?;
        if control.fix_stereo_script.is_empty() {
            return None;
        }
        control.events.push("FixRestoredStructureStereo".to_owned());
        Some(control.fix_stereo_script.remove(0))
    })
}

#[cfg(test)]
fn normalize_and_compare_test_event(event: String) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        if let Some(control) = control.borrow_mut().as_mut() {
            control.events.push(event);
        }
    });
}

#[cfg(test)]
fn normalize_and_compare_test_cleanup_inchi(
    heap: &SourceHeap,
    index: usize,
    pointer: SourceMutPointer<INChI>,
) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let Some(control) = control.as_mut() else {
            return;
        };
        let present = !pointer.is_null();
        control.events.push(format!(
            "Free_INChI[{index}]:{}",
            if present { "present" } else { "null" }
        ));
        if present {
            control.prefree_state_exact &= heap.slice(pointer.as_const()).is_ok_and(|values| {
                values.first().is_some_and(|inchi| {
                    inchi.nErrorCode == 100 + index as i32
                        && inchi.nNumberOfAtoms == control.expected_inchi_atoms[index]
                        && inchi.nRefCount == 0
                })
            });
            if let Ok(values) = heap.slice(pointer.as_const())
                && let Some(inchi) = values.first()
                && !inchi.szHillFormula.is_null()
                && let Ok(formula) = heap.slice(inchi.szHillFormula.as_const())
            {
                let length = formula.iter().position(|byte| *byte == 0).unwrap_or(0);
                control.formula_before_cleanup = formula[..length]
                    .iter()
                    .map(|byte| *byte as u8 as char)
                    .collect();
                control.events.push("free:formula".to_owned());
            }
        }
    });
}

#[cfg(test)]
fn normalize_and_compare_test_cleanup_aux(
    heap: &SourceHeap,
    index: usize,
    pointer: SourceMutPointer<crate::source_types::INChI_Aux>,
) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let Some(control) = control.as_mut() else {
            return;
        };
        let present = !pointer.is_null();
        control.events.push(format!(
            "Free_INChI_Aux[{index}]:{}",
            if present { "present" } else { "null" }
        ));
        if present {
            control.prefree_state_exact &= heap.slice(pointer.as_const()).is_ok_and(|values| {
                values.first().is_some_and(|aux| {
                    aux.nErrorCode == 200 + index as i32
                        && aux.nNumberOfAtoms == 20 + index as i32
                        && aux.nRefCount == 0
                })
            });
        }
    });
}

#[cfg(test)]
fn normalize_and_compare_test_cleanup_norm(
    heap: &SourceHeap,
    index: usize,
    pointer: SourceMutPointer<crate::source_types::INP_ATOM_DATA>,
) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let Some(control) = control.as_mut() else {
            return;
        };
        let present = !pointer.is_null();
        control.events.push(format!(
            "FreeInpAtomData[{index}]:{}",
            if present { "present" } else { "null" }
        ));
        if present {
            control.prefree_state_exact &= heap.slice(pointer.as_const()).is_ok_and(|values| {
                values.first().is_some_and(|data| {
                    let expected_at = control
                        .final_cleanup
                        .as_ref()
                        .map_or(SourceMutPointer::null(), |final_cleanup| {
                            final_cleanup.expected_norm_at[index]
                        });
                    data.num_at == 30 + index as i32
                        && data.num_bonds == 40 + index as i32
                        && data.at == expected_at
                        && data.at_fixed_bonds.is_null()
                })
            });
        }
    });
}

#[cfg(test)]
fn normalize_and_compare_test_cleanup_t_group(
    heap: &SourceHeap,
    info: &crate::source_types::T_GROUP_INFO,
) {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        let mut control = control.borrow_mut();
        let Some(control) = control.as_mut() else {
            return;
        };
        let present = !info.t_group.is_null();
        control.events.push(format!(
            "free_t_group_info:{}",
            if present { "present" } else { "null" }
        ));
        if let Some(final_cleanup) = &control.final_cleanup {
            control.prefree_state_exact &= match &final_cleanup.expected_t_groups {
                Some(expected) => {
                    present
                        && info.num_t_groups == expected.len() as i32
                        && heap
                            .slice(info.t_group.as_const())
                            .is_ok_and(|groups| groups == expected)
                }
                None => !present && info.num_t_groups == 0,
            };
        } else {
            control.prefree_state_exact &= present
                && info.num_t_groups == 1
                && heap.slice(info.t_group.as_const()).is_ok_and(|groups| {
                    groups.first().is_some_and(|group| {
                        group.nNumEndpoints == 51 && group.num == [52, 53, 0, 0, 0]
                    })
                });
        }
    });
}

#[cfg(test)]
fn normalize_and_compare_test_finish() -> NormalizeAndCompareTestControl {
    NORMALIZE_AND_COMPARE_TEST_CONTROL.with(|control| {
        control
            .borrow_mut()
            .take()
            .expect("NormalizeAndCompare test control must be active")
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        BNS_EDGE, BNS_VERTEX, INChI_Aux, INP_ATOM_DATA, REQ_MODE_NON_ISO, REQ_MODE_TAUT, SRM,
        T_GROUP, T_GROUP_INFO, TAUT_INI, TC_GROUP, TG_FLAG_FIX_ISO_FIXEDH_BUG,
        TG_FLAG_FIX_TERM_H_CHRG_BUG,
    };

    fn absent_groups() -> ALL_TC_GROUPS {
        ALL_TC_GROUPS {
            nGroup: [-1; 18],
            ..ALL_TC_GROUPS::default()
        }
    }

    #[test]
    fn source_port__ichirvr4__restoreatommakebns__line_3198() {
        struct SingleFixture {
            heap: SourceHeap,
            structure: StrFromINChI,
            atoms: SourceMutPointer<inp_ATOM>,
            inchi: SourceMutPointer<INChI>,
        }

        fn single_fixture() -> SingleFixture {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM {
                    orig_at_number: 1,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: 1,
                    nTotalCharge: 130,
                    ..INChI::default()
                }])
                .unwrap();
            SingleFixture {
                heap,
                structure: StrFromINChI {
                    at: atoms,
                    num_atoms: 1,
                    ..StrFromINChI::default()
                },
                atoms,
                inchi,
            }
        }

        let mut first = single_fixture();
        let first_live = first.heap.live_allocation_count();
        first.heap.fail_after_allocations(0);
        assert_eq!(
            RestoreAtomMakeBNS(
                &mut first.heap,
                SourceMutPointer::null(),
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                &mut first.structure,
                3,
                5,
                [first.inchi, SourceMutPointer::null()],
                SourceConstPointer::null(),
                7,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(first.structure.at2.is_null());
        assert_eq!(first.heap.source_allocation_calls(), 2);
        assert_eq!(first.heap.live_allocation_count(), first_live);
        assert_eq!(
            first.heap.slice(first.atoms.as_const()).unwrap()[0].charge,
            -126
        );

        let mut second = single_fixture();
        let second_live = second.heap.live_allocation_count();
        second.heap.fail_after_allocations(1);
        assert_eq!(
            RestoreAtomMakeBNS(
                &mut second.heap,
                SourceMutPointer::null(),
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                &mut second.structure,
                11,
                13,
                [second.inchi, SourceMutPointer::null()],
                SourceConstPointer::null(),
                17,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(!second.structure.at2.is_null());
        assert_eq!(second.heap.source_allocation_calls(), 2);
        assert_eq!(second.heap.live_allocation_count(), second_live + 1);
        assert_eq!(
            second.heap.slice(second.structure.at2.as_const()).unwrap()[0].charge,
            -126
        );

        struct MultiFixture {
            heap: SourceHeap,
            structure: StrFromINChI,
            inchi: SourceMutPointer<INChI>,
            clock: SourceMutPointer<INCHI_CLOCK>,
        }

        fn multi_fixture() -> MultiFixture {
            let mut heap = SourceHeap::default();
            let mut first = inp_ATOM {
                el_number: 6,
                valence: 1,
                chem_bonds_valence: 1,
                num_H: 3,
                orig_at_number: 1,
                ..inp_ATOM::default()
            };
            first.elname[0] = b'C' as i8;
            first.neighbor[0] = 1;
            first.bond_type[0] = BOND_TYPE_SINGLE as u8;
            let mut second = first.clone();
            second.orig_at_number = 2;
            second.neighbor[0] = 0;
            let atoms = heap.allocate_model_storage(vec![first, second]).unwrap();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let formula = heap
                .allocate_model_storage(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'6' as i8, 0])
                .unwrap();
            let elements = heap.allocate_model_storage(vec![6_u8, 6]).unwrap();
            let hydrogens = heap.allocate_model_storage(vec![3_i8, 3]).unwrap();
            let connections = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: 2,
                    szHillFormula: formula,
                    nAtom: elements,
                    nNum_H: hydrogens,
                    lenConnTable: 2,
                    nConnTable: connections,
                    ..INChI::default()
                }])
                .unwrap();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            MultiFixture {
                heap,
                structure: StrFromINChI {
                    at: atoms,
                    num_atoms: 2,
                    pSrm: restore_mode.as_const(),
                    bMobileH: TAUT_YES as i8,
                    iMobileH: TAUT_YES as i8,
                    ..StrFromINChI::default()
                },
                inchi,
                clock,
            }
        }

        let mut valence_failure = multi_fixture();
        let valence_live = valence_failure.heap.live_allocation_count();
        valence_failure.heap.fail_after_allocations(0);
        assert_eq!(
            RestoreAtomMakeBNS(
                &mut valence_failure.heap,
                valence_failure.clock,
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                &mut valence_failure.structure,
                19,
                23,
                [valence_failure.inchi, SourceMutPointer::null()],
                SourceConstPointer::null(),
                29,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(valence_failure.structure.pVA.is_null());
        assert!(valence_failure.structure.pbfsq.is_null());
        assert_eq!(valence_failure.heap.source_allocation_calls(), 1);
        assert_eq!(valence_failure.heap.live_allocation_count(), valence_live);

        let mut queue_failure = multi_fixture();
        let queue_live = queue_failure.heap.live_allocation_count();
        queue_failure.heap.fail_after_allocations(1);
        assert_eq!(
            RestoreAtomMakeBNS(
                &mut queue_failure.heap,
                queue_failure.clock,
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                &mut queue_failure.structure,
                31,
                37,
                [queue_failure.inchi, SourceMutPointer::null()],
                SourceConstPointer::null(),
                41,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(!queue_failure.structure.pVA.is_null());
        assert!(queue_failure.structure.pbfsq.is_null());
        assert_eq!(queue_failure.heap.source_allocation_calls(), 4);
        assert_eq!(queue_failure.heap.live_allocation_count(), queue_live + 1);
        assert_eq!(
            queue_failure
                .heap
                .slice(queue_failure.structure.pVA.as_const())
                .unwrap(),
            &[VAL_AT::default(), VAL_AT::default()]
        );

        let mut complete = multi_fixture();
        let input_parameters = INPUT_PARMS {
            nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
            bTautFlags: u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG),
            ..INPUT_PARMS::default()
        };
        let source_atoms = complete
            .heap
            .slice(complete.structure.at.as_const())
            .unwrap()
            .to_vec();
        let generated_at = complete
            .heap
            .allocate_model_storage(source_atoms.clone())
            .unwrap();
        let generated_scratch = complete
            .heap
            .allocate_model_storage(vec![inp_ATOM::default(); 2])
            .unwrap();
        let mut generated_structure = StrFromINChI {
            num_atoms: 2,
            bMobileH: TAUT_YES as i8,
            iMobileH: TAUT_YES as i8,
            ..StrFromINChI::default()
        };
        assert_eq!(
            MakeOneInChIOutOfStrFromINChI(
                &mut complete.heap,
                &mut CANON_GLOBALS::default(),
                complete.clock,
                &input_parameters,
                &mut STRUCT_DATA::default(),
                &mut generated_structure,
                generated_at,
                generated_scratch,
                &absent_groups(),
                0,
            ),
            Ok(0)
        );
        let generated_input = generated_structure.pOneINChI[0];
        assert!(!generated_input.is_null());
        generated_structure.pOneINChI[0] = SourceMutPointer::null();
        complete.heap.trace_source_allocations();
        assert_eq!(
            RestoreAtomMakeBNS(
                &mut complete.heap,
                complete.clock,
                &mut CANON_GLOBALS::default(),
                &input_parameters,
                &mut STRUCT_DATA::default(),
                &mut complete.structure,
                0,
                0,
                [generated_input, SourceMutPointer::null()],
                SourceConstPointer::null(),
                43,
                0,
                0,
            ),
            Ok(0)
        );
        assert!(complete.structure.pbfsq.is_null());
        assert!(!complete.structure.pVA.is_null());
        assert!(!complete.structure.at2.is_null());
        assert_eq!(
            complete
                .heap
                .slice(complete.structure.at2.as_const())
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn source_port__ichirvr4__oneinchi2atom__line_3477() {
        struct Fixture {
            heap: SourceHeap,
            clock: SourceMutPointer<INCHI_CLOCK>,
            input: SourceMutPointer<INChI>,
            parameters: INPUT_PARMS,
            structure: StrFromINChI,
        }

        fn fixture(first_struct_number: i64) -> Fixture {
            let mut heap = SourceHeap::default();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let mut carbon = inp_ATOM {
                el_number: 6,
                num_H: 4,
                orig_at_number: 1,
                component: 1,
                ..inp_ATOM::default()
            };
            carbon.elname[0] = b'C' as i8;
            let at = heap.allocate_model_storage(vec![carbon.clone()]).unwrap();
            let scratch = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let parameters = INPUT_PARMS {
                first_struct_number,
                nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                bTautFlags: u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG),
                ..INPUT_PARMS::default()
            };
            let mut generated = StrFromINChI {
                num_atoms: 1,
                bMobileH: TAUT_YES as i8,
                iMobileH: TAUT_YES as i8,
                ..StrFromINChI::default()
            };
            assert_eq!(
                MakeOneInChIOutOfStrFromINChI(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &parameters,
                    &mut STRUCT_DATA::default(),
                    &mut generated,
                    at,
                    scratch,
                    &absent_groups(),
                    0,
                ),
                Ok(0)
            );
            let input = generated.pOneINChI[0];
            assert!(!input.is_null());
            generated.pOneINChI[0] = SourceMutPointer::null();
            Fixture {
                heap,
                clock,
                input,
                parameters,
                structure: StrFromINChI {
                    pSrm: restore_mode.as_const(),
                    bMobileH: TAUT_YES as i8,
                    iMobileH: TAUT_YES as i8,
                    ..StrFromINChI::default()
                },
            }
        }

        let mut allocation = fixture(0);
        let mut allocation_data = STRUCT_DATA::default();
        allocation_data.pStrErrStruct[0] = 77;
        allocation.heap.fail_after_allocations(0);
        assert_eq!(
            OneInChI2Atom(
                &mut allocation.heap,
                allocation.clock,
                &mut CANON_GLOBALS::default(),
                &allocation.parameters,
                &mut allocation_data,
                SourceConstPointer::null(),
                1,
                &mut allocation.structure,
                0,
                0,
                0,
                [allocation.input, SourceMutPointer::null()],
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(allocation_data.pStrErrStruct[0], 0);
        assert_eq!(allocation.heap.source_allocation_calls(), 1);
        assert!(allocation.structure.at.is_null());

        let mut gated = fixture(2);
        gated.structure.RevInChI.nRetVal = 73;
        let mut gated_data = STRUCT_DATA::default();
        gated_data.pStrErrStruct[0] = 79;
        assert_eq!(
            OneInChI2Atom(
                &mut gated.heap,
                gated.clock,
                &mut CANON_GLOBALS::default(),
                &gated.parameters,
                &mut gated_data,
                SourceConstPointer::null(),
                1,
                &mut gated.structure,
                3,
                5,
                0,
                [gated.input, SourceMutPointer::null()],
                0,
            ),
            Ok(0)
        );
        assert_eq!(gated_data.pStrErrStruct[0], 0);
        assert_eq!(gated.structure.RevInChI.nRetVal, 73);
        assert_eq!(
            gated.structure.RevInChI.pINChI,
            [SourceMutPointer::null(); 2]
        );
        assert!(gated.structure.pbfsq.is_null());
        assert!(gated.structure.pVA.is_null());

        let mut generated = fixture(0);
        generated.structure.iInchiRec = INCHI_REC as i8;
        let original_mode = generated.parameters.nMode;
        assert_eq!(
            OneInChI2Atom(
                &mut generated.heap,
                generated.clock,
                &mut CANON_GLOBALS::default(),
                &generated.parameters,
                &mut STRUCT_DATA::default(),
                SourceConstPointer::null(),
                1,
                &mut generated.structure,
                7,
                11,
                1,
                [generated.input, SourceMutPointer::null()],
                0,
            ),
            Ok(0)
        );
        assert_eq!(generated.parameters.nMode, original_mode);
        assert_eq!(generated.structure.RevInChI.nRetVal, 0);
        assert!(
            generated
                .structure
                .RevInChI
                .pINChI
                .iter()
                .any(|pointer| !pointer.is_null())
        );
        assert!(generated.structure.pbfsq.is_null());
        assert!(!generated.structure.at.is_null());
        assert!(!generated.structure.at2.is_null());
    }

    #[test]
    fn source_port__ichirvr4__makeprotoncomponent__line_3574() {
        let mut no_op_heap = SourceHeap::default();
        let retained_at = no_op_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let retained_at2 = no_op_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut no_op = StrFromINChI {
            at: retained_at,
            at2: retained_at2,
            bDeleted: 7,
            num_atoms: 11,
            bMobileH: 12,
            iMobileH: 13,
            ..StrFromINChI::default()
        };
        let unchanged = no_op.clone();
        no_op_heap.trace_source_allocations();
        assert_eq!(
            MakeProtonComponent(&mut no_op_heap, &mut no_op, i32::MIN, 0),
            Ok(0)
        );
        assert_eq!(no_op, unchanged);
        assert_eq!(no_op_heap.source_allocation_calls(), 0);
        assert_eq!(
            MakeProtonComponent(&mut no_op_heap, &mut no_op, i32::MAX, -1),
            Ok(0)
        );
        assert_eq!(no_op, unchanged);

        let mut heap = SourceHeap::default();
        let old_at = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let old_at2 = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut structure = StrFromINChI {
            at: old_at,
            at2: old_at2,
            bDeleted: -1,
            num_atoms: -2,
            bMobileH: -3,
            iMobileH: -4,
            ..StrFromINChI::default()
        };
        heap.trace_source_allocations();
        assert_eq!(
            MakeProtonComponent(&mut heap, &mut structure, i32::MAX, 3),
            Ok(3)
        );
        assert_eq!(heap.source_allocation_calls(), 2);
        assert_ne!(structure.at, old_at);
        assert_ne!(structure.at2, old_at2);
        assert_ne!(structure.at, structure.at2);
        assert_eq!(heap.live_allocation_count(), 4);
        let primary = heap.slice(structure.at.as_const()).unwrap();
        let copy = heap.slice(structure.at2.as_const()).unwrap();
        assert_eq!(primary, copy);
        assert_eq!(primary.len(), 3);
        for (index, atom) in primary.iter().enumerate() {
            assert_eq!(atom.elname, [b'H' as i8, 0, 0, 0, 0, 0]);
            assert_eq!(atom.el_number, 1);
            assert_eq!(atom.orig_at_number, (index + 1) as AT_NUMB);
            assert_eq!(atom.charge, 1);
            let mut expected = inp_ATOM::default();
            expected.elname[0] = b'H' as i8;
            expected.el_number = 1;
            expected.orig_at_number = (index + 1) as AT_NUMB;
            expected.charge = 1;
            assert_eq!(*atom, expected);
        }
        assert_eq!(
            (
                structure.bDeleted,
                structure.num_atoms,
                structure.bMobileH,
                structure.iMobileH,
            ),
            (0, 3, TAUT_YES as i8, TAUT_YES as i8)
        );

        let mut first_failure_heap = SourceHeap::default();
        let first_old_at = first_failure_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let first_old_at2 = first_failure_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut first_failure = StrFromINChI {
            at: first_old_at,
            at2: first_old_at2,
            num_atoms: 91,
            bDeleted: 5,
            ..StrFromINChI::default()
        };
        first_failure_heap.fail_after_allocations(0);
        assert_eq!(
            MakeProtonComponent(&mut first_failure_heap, &mut first_failure, 0, 2),
            Ok(0)
        );
        assert!(first_failure.at.is_null());
        assert!(!first_failure.at2.is_null());
        assert_eq!((first_failure.num_atoms, first_failure.bDeleted), (91, 5));
        assert_eq!(first_failure_heap.source_allocation_calls(), 2);
        assert_eq!(first_failure_heap.live_allocation_count(), 3);
        assert_eq!(
            first_failure_heap
                .slice(first_failure.at2.as_const())
                .unwrap(),
            &[inp_ATOM::default(), inp_ATOM::default()]
        );

        let mut second_failure_heap = SourceHeap::default();
        let second_old_at2 = second_failure_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut second_failure = StrFromINChI {
            at2: second_old_at2,
            num_atoms: 73,
            bMobileH: 6,
            ..StrFromINChI::default()
        };
        second_failure_heap.fail_after_allocations(1);
        assert_eq!(
            MakeProtonComponent(&mut second_failure_heap, &mut second_failure, 0, 2),
            Ok(0)
        );
        assert!(!second_failure.at.is_null());
        assert!(second_failure.at2.is_null());
        assert_eq!((second_failure.num_atoms, second_failure.bMobileH), (73, 6));
        assert_eq!(second_failure_heap.source_allocation_calls(), 2);
        assert_eq!(second_failure_heap.live_allocation_count(), 2);
        assert_eq!(
            second_failure_heap
                .slice(second_failure.at.as_const())
                .unwrap(),
            &[inp_ATOM::default(), inp_ATOM::default()]
        );
    }

    #[test]
    fn source_port__ichirvr4__addremprotonsinrestrstruct__line_3615() {
        fn call(
            heap: &mut SourceHeap,
            structures: &mut [StrFromINChI],
            num_components: i32,
            reconnected: Option<&[StrFromINChI]>,
            num_reconnected: i32,
            fixed_h: i32,
            balance: &mut NUM_H,
            recmet: Option<&mut i32>,
        ) -> Result<i32, SourceHeapError> {
            AddRemProtonsInRestrStruct(
                heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                17,
                fixed_h,
                structures,
                num_components,
                reconnected,
                num_reconnected,
                balance,
                recmet,
                0,
            )
        }

        fn with_aux(
            heap: &mut SourceHeap,
            atoms: Vec<inp_ATOM>,
            num_atoms: i32,
            num_deleted_h: i32,
            representation: usize,
            aux: INChI_Aux,
        ) -> StrFromINChI {
            let at2 = heap.allocate_model_storage(atoms).unwrap();
            let aux = heap.allocate_model_storage(vec![aux]).unwrap();
            let aux_rows = heap
                .allocate_model_storage(vec![[SourceMutPointer::null(), aux]])
                .unwrap();
            let mut structure = StrFromINChI {
                at2,
                num_atoms,
                num_deleted_H: num_deleted_h,
                ..StrFromINChI::default()
            };
            structure.RevInChI.pINChI_Aux[representation] = aux_rows;
            structure.RevInChI.num_components[representation] = 1;
            structure
        }

        let mut no_op_heap = SourceHeap::default();
        let mut no_op = [StrFromINChI {
            nLink: i32::MIN,
            bPostProcessed: 71,
            ..StrFromINChI::default()
        }];
        let unchanged = no_op.clone();
        let mut zero_balance: NUM_H = 0;
        let mut zero_recmet = 93;
        assert_eq!(
            call(
                &mut no_op_heap,
                &mut no_op,
                i32::MAX,
                None,
                i32::MAX,
                0,
                &mut zero_balance,
                Some(&mut zero_recmet),
            ),
            Ok(0)
        );
        assert_eq!(no_op, unchanged);
        assert_eq!((zero_balance, zero_recmet), (0, 93));

        let mut link_heap = SourceHeap::default();
        let mut disconnected = [StrFromINChI {
            nLink: -1,
            ..StrFromINChI::default()
        }];
        let mut link_balance: NUM_H = -1;
        let mut link_recmet = 41;
        assert_eq!(
            call(
                &mut link_heap,
                &mut disconnected,
                1,
                None,
                1,
                0,
                &mut link_balance,
                Some(&mut link_recmet),
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!((link_balance, link_recmet), (-1, 41));
        let valid_reconnected = [StrFromINChI {
            nLink: 1,
            ..StrFromINChI::default()
        }];
        assert_eq!(
            call(
                &mut link_heap,
                &mut disconnected,
                1,
                Some(&valid_reconnected),
                1,
                0,
                &mut link_balance,
                Some(&mut link_recmet),
            ),
            Ok(0)
        );
        assert_eq!((link_balance, link_recmet), (-1, 0));

        let mut reconnect_heap = SourceHeap::default();
        let mut reconnect = [with_aux(
            &mut reconnect_heap,
            vec![inp_ATOM {
                el_number: 6,
                ..inp_ATOM::default()
            }],
            1,
            0,
            INCHI_BAS as usize,
            INChI_Aux {
                bNormalizationFlags: u64::MAX,
                nNumberOfTGroups: 0,
                ..INChI_Aux::default()
            },
        )];
        reconnect[0].nLink = 1;
        let mut reconnect_balance: NUM_H = -1;
        let mut reconnect_recmet = 77;
        assert_eq!(
            call(
                &mut reconnect_heap,
                &mut reconnect,
                1,
                None,
                0,
                1,
                &mut reconnect_balance,
                Some(&mut reconnect_recmet),
            ),
            Ok(0)
        );
        assert_eq!(reconnect[0].bPostProcessed, 0);
        assert_eq!((reconnect_balance, reconnect_recmet), (-1, 0));
        assert_eq!(
            reconnect_heap.slice(reconnect[0].at2.as_const()).unwrap()[0].el_number,
            6
        );

        let mut rec_aux_heap = SourceHeap::default();
        let mut rec_aux = with_aux(
            &mut rec_aux_heap,
            vec![inp_ATOM {
                el_number: 6,
                ..inp_ATOM::default()
            }],
            1,
            0,
            INCHI_REC as usize,
            INChI_Aux::default(),
        );
        let inchi = rec_aux_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        rec_aux.RevInChI.pINChI[INCHI_REC as usize] = rec_aux_heap
            .allocate_model_storage(vec![[SourceMutPointer::null(), inchi]])
            .unwrap();
        let mut rec_aux = [rec_aux];
        let mut rec_aux_balance: NUM_H = -1;
        assert_eq!(
            call(
                &mut rec_aux_heap,
                &mut rec_aux,
                1,
                None,
                0,
                1,
                &mut rec_aux_balance,
                None,
            ),
            Ok(0)
        );
        assert_eq!(rec_aux_balance, -1);

        let mut bad_h_heap = SourceHeap::default();
        let mut bad_h = [with_aux(
            &mut bad_h_heap,
            vec![
                inp_ATOM {
                    el_number: 6,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    el_number: 1,
                    neighbor: [0; 20],
                    ..inp_ATOM::default()
                },
            ],
            1,
            1,
            INCHI_BAS as usize,
            INChI_Aux::default(),
        )];
        let mut bad_h_balance: NUM_H = -1;
        let mut bad_h_recmet = 55;
        assert_eq!(
            call(
                &mut bad_h_heap,
                &mut bad_h,
                1,
                None,
                0,
                1,
                &mut bad_h_balance,
                Some(&mut bad_h_recmet),
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!((bad_h_balance, bad_h_recmet), (-1, 55));
        assert_eq!(bad_h[0].bPostProcessed, 0);

        let mut changed_heap = SourceHeap::default();
        let mut changed_atoms = vec![
            inp_ATOM {
                el_number: 8,
                valence: 1,
                chem_bonds_valence: 1,
                charge: -1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                valence: 1,
                chem_bonds_valence: 4,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 1,
                charge: 1,
                ..inp_ATOM::default()
            },
        ];
        changed_atoms[0].neighbor[0] = 1;
        changed_atoms[1].neighbor[0] = 0;
        let mut changed = [with_aux(
            &mut changed_heap,
            changed_atoms,
            3,
            0,
            INCHI_BAS as usize,
            INChI_Aux::default(),
        )];
        changed[0].iInchiRec = INCHI_REC as i8;
        changed[0].iMobileH = TAUT_YES as i8;
        let mut changed_balance: NUM_H = 1;
        let mut changed_recmet = 66;
        changed_heap.fail_after_allocations(0);
        assert_eq!(
            call(
                &mut changed_heap,
                &mut changed,
                1,
                None,
                0,
                1,
                &mut changed_balance,
                Some(&mut changed_recmet),
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(changed[0].bPostProcessed, 1);
        let changed_atom = &changed_heap.slice(changed[0].at2.as_const()).unwrap()[0];
        assert_eq!((changed_atom.charge, changed_atom.num_H), (0, 1));
        assert_eq!((changed_balance, changed_recmet), (1, 66));
        assert_eq!(changed[0].RevInChI.num_components, [0, 0]);

        let mut deleted_failure_heap = SourceHeap::default();
        let mut deleted_failure = [StrFromINChI {
            bDeleted: 1,
            ..StrFromINChI::default()
        }];
        let mut deleted_failure_balance: NUM_H = 2;
        let mut deleted_failure_recmet = 88;
        deleted_failure_heap.fail_after_allocations(0);
        assert_eq!(
            call(
                &mut deleted_failure_heap,
                &mut deleted_failure,
                1,
                None,
                0,
                0,
                &mut deleted_failure_balance,
                Some(&mut deleted_failure_recmet),
            ),
            Ok(0)
        );
        assert!(deleted_failure[0].at.is_null());
        assert!(!deleted_failure[0].at2.is_null());
        assert_eq!((deleted_failure_balance, deleted_failure_recmet), (2, 0));

        let mut deleted_recalc_heap = SourceHeap::default();
        let mut deleted_recalc = [StrFromINChI {
            bDeleted: 1,
            ..StrFromINChI::default()
        }];
        let mut deleted_recalc_balance: NUM_H = 1;
        let mut deleted_recalc_recmet = 99;
        deleted_recalc_heap.fail_after_allocations(2);
        assert_eq!(
            call(
                &mut deleted_recalc_heap,
                &mut deleted_recalc,
                1,
                None,
                0,
                0,
                &mut deleted_recalc_balance,
                Some(&mut deleted_recalc_recmet),
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(
            (
                deleted_recalc[0].bDeleted,
                deleted_recalc[0].num_atoms,
                deleted_recalc[0].bMobileH,
                deleted_recalc[0].iMobileH,
            ),
            (0, 1, TAUT_YES as i8, TAUT_YES as i8)
        );
        assert_eq!((deleted_recalc_balance, deleted_recalc_recmet), (1, 99));
    }

    #[test]
    fn source_port__ichirvr4__addremisoprotonsinrestrstruct__line_3816() {
        fn call(
            heap: &mut SourceHeap,
            clock: SourceMutPointer<INCHI_CLOCK>,
            canonical_globals: SourceMutPointer<CANON_GLOBALS>,
            structures: &mut [StrFromINChI],
            num_components: i32,
            reconnected: Option<&[StrFromINChI]>,
            num_reconnected: i32,
            fixed_h: i32,
            balance: &mut [NUM_H; NUM_H_ISOTOPES as usize],
            recmet: Option<&mut [NUM_H; NUM_H_ISOTOPES as usize]>,
        ) -> Result<i32, SourceHeapError> {
            AddRemIsoProtonsInRestrStruct(
                heap,
                clock,
                canonical_globals,
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                19,
                fixed_h,
                structures,
                num_components,
                reconnected,
                num_reconnected,
                balance,
                recmet,
                0,
            )
        }

        fn add_aux(
            heap: &mut SourceHeap,
            structure: &mut StrFromINChI,
            representation: usize,
            number_of_atoms: i32,
            number_of_t_groups: i32,
        ) {
            let aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumberOfAtoms: number_of_atoms,
                    nNumberOfTGroups: number_of_t_groups,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let rows = heap
                .allocate_model_storage(vec![[SourceMutPointer::null(), aux]])
                .unwrap();
            structure.RevInChI.pINChI_Aux[representation] = rows;
            structure.RevInChI.num_components[representation] = 1;
        }

        fn structure_with_atom(
            heap: &mut SourceHeap,
            atom: inp_ATOM,
            representation: usize,
            number_of_atoms: i32,
            number_of_t_groups: i32,
        ) -> StrFromINChI {
            let at2 = heap.allocate_model_storage(vec![atom]).unwrap();
            let mut structure = StrFromINChI {
                at2,
                num_atoms: 1,
                ..StrFromINChI::default()
            };
            add_aux(
                heap,
                &mut structure,
                representation,
                number_of_atoms,
                number_of_t_groups,
            );
            structure
        }

        let mut zero_heap = SourceHeap::default();
        let mut zero_structures = [StrFromINChI {
            nLink: i32::MIN,
            bPostProcessed: 73,
            ..StrFromINChI::default()
        }];
        let zero_before = zero_structures.clone();
        let mut zero_balance = [0; NUM_H_ISOTOPES as usize];
        let mut zero_recmet = [11, 13, 17];
        assert_eq!(
            call(
                &mut zero_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut zero_structures,
                i32::MAX,
                None,
                i32::MAX,
                0,
                &mut zero_balance,
                Some(&mut zero_recmet),
            ),
            Ok(0)
        );
        assert_eq!(zero_structures, zero_before);
        assert_eq!(zero_recmet, [11, 13, 17]);

        let mut link_heap = SourceHeap::default();
        let mut disconnected = [StrFromINChI {
            nLink: -1,
            ..StrFromINChI::default()
        }];
        let mut link_balance = [1, 2, 3];
        let mut link_recmet = [19, 23, 29];
        assert_eq!(
            call(
                &mut link_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut disconnected,
                1,
                None,
                1,
                0,
                &mut link_balance,
                Some(&mut link_recmet),
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(link_balance, [1, 2, 3]);
        assert_eq!(link_recmet, [19, 23, 29]);
        let valid_reconnected = [StrFromINChI {
            nLink: 1,
            ..StrFromINChI::default()
        }];
        assert_eq!(
            call(
                &mut link_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut disconnected,
                1,
                Some(&valid_reconnected),
                1,
                0,
                &mut link_balance,
                Some(&mut link_recmet),
            ),
            Ok(0)
        );
        assert_eq!(link_balance, [1, 2, 3]);
        assert_eq!(link_recmet, [0, 0, 0]);

        let mut empty_heap = SourceHeap::default();
        let mut empty = [StrFromINChI::default()];
        let mut empty_balance = [1, 0, 0];
        let mut empty_recmet = [31, 37, 41];
        assert_eq!(
            call(
                &mut empty_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty,
                1,
                None,
                0,
                0,
                &mut empty_balance,
                Some(&mut empty_recmet),
            ),
            Ok(0)
        );
        assert_eq!(empty_balance, [1, 0, 0]);
        assert_eq!(empty_recmet, [0, 0, 0]);
        let mut negative_count_recmet = [43, 47, 53];
        assert_eq!(
            call(
                &mut empty_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty,
                i32::MIN,
                None,
                0,
                0,
                &mut empty_balance,
                Some(&mut negative_count_recmet),
            ),
            Ok(0)
        );
        assert_eq!(negative_count_recmet, [0, 0, 0]);

        let endpoint = inp_ATOM {
            el_number: 7,
            endpoint: 1,
            num_H: 1,
            chem_bonds_valence: 3,
            ..inp_ATOM::default()
        };
        let mut base_heap = SourceHeap::default();
        let mut base = [structure_with_atom(
            &mut base_heap,
            endpoint.clone(),
            INCHI_BAS as usize,
            1,
            0,
        )];
        let mut base_balance = [0, 0, 1];
        assert_eq!(
            call(
                &mut base_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut base,
                1,
                None,
                0,
                0,
                &mut base_balance,
                None,
            ),
            Ok(0)
        );
        assert_eq!(base_balance, [0, 0, 1]);
        assert_eq!(
            base_heap.slice(base[0].at2.as_const()).unwrap()[0].num_iso_H,
            [0, 0, 0]
        );

        let mut rec_zero_heap = SourceHeap::default();
        let mut rec_zero = structure_with_atom(
            &mut rec_zero_heap,
            endpoint.clone(),
            INCHI_BAS as usize,
            1,
            0,
        );
        add_aux(&mut rec_zero_heap, &mut rec_zero, INCHI_REC as usize, 0, 1);
        let mut rec_zero = [rec_zero];
        let mut rec_zero_balance = [0, 1, 0];
        assert_eq!(
            call(
                &mut rec_zero_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut rec_zero,
                1,
                None,
                0,
                0,
                &mut rec_zero_balance,
                None,
            ),
            Ok(0)
        );
        assert_eq!(rec_zero_balance, [0, 1, 0]);

        let mut negative_heap = SourceHeap::default();
        let mut negative = [structure_with_atom(
            &mut negative_heap,
            inp_ATOM::default(),
            INCHI_BAS as usize,
            1,
            0,
        )];
        negative[0].bPostProcessed = 4;
        let mut negative_balance = [0, 0, -1];
        let mut negative_recmet = [59, 61, 67];
        assert_eq!(
            call(
                &mut negative_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut negative,
                1,
                None,
                0,
                0,
                &mut negative_balance,
                Some(&mut negative_recmet),
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(negative[0].bPostProcessed, 4 | RI_ERR_PROGR);
        assert_eq!(negative_balance, [0, 0, -1]);
        assert_eq!(negative_recmet, [59, 61, 67]);

        let mut recalc_heap = SourceHeap::default();
        let mut recalc = [structure_with_atom(
            &mut recalc_heap,
            endpoint,
            INCHI_REC as usize,
            1,
            1,
        )];
        recalc[0].bPostProcessed = 2;
        recalc[0].iInchiRec = INCHI_REC as i8;
        recalc[0].iMobileH = TAUT_YES as i8;
        recalc[0].nLink = 1;
        let recalc_at = recalc[0].at2;
        let mut recalc_balance = [0, 0, 1];
        let mut recalc_recmet = [71, 73, 79];
        recalc_heap.fail_after_allocations(0);
        assert_eq!(
            call(
                &mut recalc_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut recalc,
                1,
                None,
                0,
                1,
                &mut recalc_balance,
                Some(&mut recalc_recmet),
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(recalc[0].bPostProcessed, 3);
        let recalc_atom = &recalc_heap.slice(recalc_at.as_const()).unwrap()[0];
        assert_eq!((recalc_atom.num_H, recalc_atom.num_iso_H), (1, [0, 0, 1]));
        assert_eq!(recalc_balance, [0, 0, 1]);
        assert_eq!(recalc_recmet, [71, 73, 79]);
        assert_eq!(recalc[0].RevInChI.num_components, [0, 0]);
    }

    #[test]
    fn source_port__ichirvr4__runbnsrestore1__line_2740() {
        struct AllocationFixture {
            heap: SourceHeap,
            structure: StrFromINChI,
            atoms: SourceMutPointer<inp_ATOM>,
            old_at2: SourceMutPointer<inp_ATOM>,
        }

        fn fixture(mobile_h: i8) -> AllocationFixture {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM {
                    orig_at_number: 1,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            let old_at2 = heap
                .allocate_model_storage(vec![inp_ATOM {
                    orig_at_number: 77,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            AllocationFixture {
                heap,
                structure: StrFromINChI {
                    at: atoms,
                    at2: old_at2,
                    num_atoms: 1,
                    iMobileH: mobile_h,
                    ..StrFromINChI::default()
                },
                atoms,
                old_at2,
            }
        }

        let mut first = fixture(TAUT_YES as i8);
        let first_live = first.heap.live_allocation_count();
        first.heap.fail_after_allocations(0);
        assert_eq!(
            RunBnsRestore1(
                &mut first.heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                &mut first.structure,
                &mut [VAL_AT::default()],
                &mut absent_groups(),
                [SourceMutPointer::null(); 2],
                3,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(first.structure.at, first.atoms);
        assert_eq!(first.structure.at2, first.old_at2);
        assert_eq!(first.heap.live_allocation_count(), first_live);
        assert_eq!(first.heap.source_allocation_calls(), 1);

        let mut second = fixture(TAUT_YES as i8);
        let second_live = second.heap.live_allocation_count();
        second.heap.fail_after_allocations(1);
        assert_eq!(
            RunBnsRestore1(
                &mut second.heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                &mut second.structure,
                &mut [VAL_AT::default()],
                &mut absent_groups(),
                [SourceMutPointer::null(); 2],
                5,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(second.structure.at, second.atoms);
        assert_eq!(second.structure.at2, second.old_at2);
        assert_eq!(second.heap.live_allocation_count(), second_live);
        assert_eq!(second.heap.source_allocation_calls(), 2);

        let mut fixed_h = fixture(TAUT_NON as i8);
        let fixed_h_live = fixed_h.heap.live_allocation_count();
        fixed_h.heap.fail_after_allocations(0);
        assert_eq!(
            RunBnsRestore1(
                &mut fixed_h.heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                &mut fixed_h.structure,
                &mut [VAL_AT::default()],
                &mut absent_groups(),
                [SourceMutPointer::null(); 2],
                7,
                1,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(fixed_h.structure.at, fixed_h.atoms);
        assert!(fixed_h.structure.at2.is_null());
        assert_eq!(fixed_h.heap.source_allocation_calls(), 2);
        assert_eq!(fixed_h.heap.live_allocation_count(), fixed_h_live + 1);
        assert_eq!(
            fixed_h.heap.slice(fixed_h.old_at2.as_const()).unwrap()[0].orig_at_number,
            77
        );
        assert!(!fixed_h.structure.fixed_H.is_null());
    }

    #[test]
    fn source_port__ichirvr4__saltbondstocoordbonds__line_2373() {
        struct Fixture {
            heap: SourceHeap,
            structure: StrFromINChI,
            atoms: SourceMutPointer<inp_ATOM>,
            copied: SourceMutPointer<inp_ATOM>,
            bns: BN_STRUCT,
            valence: Vec<VAL_AT>,
            groups: ALL_TC_GROUPS,
        }

        fn fixture(restore: SRM, inchi_record: i8) -> Fixture {
            let mut heap = SourceHeap::default();
            let restore_mode = heap.allocate_model_storage(vec![restore]).unwrap();
            let mut sodium = inp_ATOM {
                el_number: 11,
                valence: 1,
                orig_at_number: 1,
                ..inp_ATOM::default()
            };
            sodium.neighbor[0] = 1;
            sodium.bond_type[0] = BOND_TYPE_SINGLE as u8;
            let mut fluorine = inp_ATOM {
                el_number: 9,
                valence: 1,
                orig_at_number: 2,
                ..inp_ATOM::default()
            };
            fluorine.neighbor[0] = 0;
            fluorine.bond_type[0] = BOND_TYPE_SINGLE as u8;
            let atoms = heap.allocate_model_storage(vec![sodium, fluorine]).unwrap();
            let copied = heap
                .allocate_model_storage(vec![
                    inp_ATOM {
                        orig_at_number: 91,
                        ..inp_ATOM::default()
                    },
                    inp_ATOM {
                        orig_at_number: 92,
                        ..inp_ATOM::default()
                    },
                ])
                .unwrap();
            let incident_0 = heap.allocate_model_storage(vec![0]).unwrap();
            let incident_1 = heap.allocate_model_storage(vec![0]).unwrap();
            let vertices = heap
                .allocate_model_storage(vec![
                    BNS_VERTEX {
                        num_adj_edges: 1,
                        max_adj_edges: 1,
                        iedge: incident_0,
                        st_edge: crate::source_types::BNS_ST_EDGE {
                            cap: 1,
                            flow: 1,
                            ..crate::source_types::BNS_ST_EDGE::default()
                        },
                        ..BNS_VERTEX::default()
                    },
                    BNS_VERTEX {
                        num_adj_edges: 1,
                        max_adj_edges: 1,
                        iedge: incident_1,
                        st_edge: crate::source_types::BNS_ST_EDGE {
                            cap: 1,
                            flow: 1,
                            ..crate::source_types::BNS_ST_EDGE::default()
                        },
                        ..BNS_VERTEX::default()
                    },
                ])
                .unwrap();
            let edges = heap
                .allocate_model_storage(vec![BNS_EDGE {
                    neighbor1: 0,
                    neighbor12: 1,
                    neigh_ord: [0, 0],
                    cap: 1,
                    flow: 1,
                    ..BNS_EDGE::default()
                }])
                .unwrap();
            Fixture {
                heap,
                structure: StrFromINChI {
                    at: atoms,
                    num_atoms: 2,
                    pSrm: restore_mode.as_const(),
                    iInchiRec: inchi_record,
                    ..StrFromINChI::default()
                },
                atoms,
                copied,
                bns: BN_STRUCT {
                    num_atoms: 2,
                    num_vertices: 2,
                    num_bonds: 1,
                    num_edges: 1,
                    vert: vertices,
                    edge: edges,
                    tot_st_cap: 2,
                    tot_st_flow: 2,
                    ..BN_STRUCT::default()
                },
                valence: vec![
                    VAL_AT {
                        cMetal: 1,
                        cPeriodicNumber: 11,
                        ..VAL_AT::default()
                    },
                    VAL_AT {
                        cPeriodicNumber: 9,
                        ..VAL_AT::default()
                    },
                ],
                groups: absent_groups(),
            }
        }

        for (restore, record) in [
            (
                SRM {
                    bMetalAddFlower: 1,
                    ..SRM::default()
                },
                INCHI_BAS as i8,
            ),
            (SRM::default(), 1),
            (
                SRM {
                    bMetalAddFlower: 1,
                    nMetalMinBondOrder: 1,
                    ..SRM::default()
                },
                1,
            ),
        ] {
            let mut gated = fixture(restore, record);
            let mut runs = 13;
            let mut delta = -17;
            assert_eq!(
                SaltBondsToCoordBonds(
                    &mut gated.heap,
                    &mut gated.bns,
                    &mut BN_DATA::default(),
                    &mut gated.structure,
                    gated.atoms,
                    gated.copied,
                    &mut gated.valence,
                    &mut gated.groups,
                    &mut runs,
                    &mut delta,
                    0x20,
                    0,
                ),
                Ok(0)
            );
            assert_eq!((runs, delta), (13, -17));
            assert_eq!(gated.structure.at, gated.atoms);
            assert_eq!(
                gated.heap.slice(gated.copied.as_const()).unwrap()[0].orig_at_number,
                91
            );
            assert_eq!(gated.heap.source_allocation_calls(), 0);
        }

        let active_mode = SRM {
            bMetalAddFlower: 1,
            ..SRM::default()
        };
        let mut active = fixture(active_mode.clone(), 1);
        let active_allocations = active.heap.live_allocation_count();
        let mut runs = 19;
        let mut delta = 23;
        assert_eq!(
            SaltBondsToCoordBonds(
                &mut active.heap,
                &mut active.bns,
                &mut BN_DATA::default(),
                &mut active.structure,
                active.atoms,
                active.copied,
                &mut active.valence,
                &mut active.groups,
                &mut runs,
                &mut delta,
                0x40,
                0,
            ),
            Ok(0)
        );
        assert_eq!((runs, delta), (19, 23));
        assert_eq!(active.structure.at, active.atoms);
        assert_eq!(active.heap.live_allocation_count(), active_allocations);
        assert_eq!(active.heap.source_allocation_calls(), 1);
        let copied_atoms = active.heap.slice(active.copied.as_const()).unwrap();
        assert_eq!(copied_atoms[0].bond_type[0], BOND_TYPE_SINGLE as u8);
        assert_eq!(copied_atoms[1].chem_bonds_valence, 1);
        assert_eq!(
            active.heap.slice(active.bns.edge.as_const()).unwrap()[0].forbidden,
            0
        );

        let mut failure = fixture(active_mode, 1);
        let failure_allocations = failure.heap.live_allocation_count();
        failure.heap.fail_after_allocations(0);
        runs = 29;
        delta = -31;
        assert_eq!(
            SaltBondsToCoordBonds(
                &mut failure.heap,
                &mut failure.bns,
                &mut BN_DATA::default(),
                &mut failure.structure,
                failure.atoms,
                failure.copied,
                &mut failure.valence,
                &mut failure.groups,
                &mut runs,
                &mut delta,
                0x40,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!((runs, delta), (29, -31));
        assert_eq!(failure.structure.at, failure.atoms);
        assert_eq!(failure.heap.live_allocation_count(), failure_allocations);
        assert_eq!(failure.heap.source_allocation_calls(), 1);
    }

    #[test]
    fn source_port__ichirvr4__makesinglebondsmetal2chargedheteroat__line_2194() {
        let mut empty_heap = SourceHeap::default();
        let restore_mode = empty_heap
            .allocate_model_storage(vec![SRM::default()])
            .unwrap();
        let atoms = empty_heap
            .allocate_model_storage(vec![inp_ATOM {
                orig_at_number: 7,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let copied = empty_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let vertices = empty_heap
            .allocate_model_storage(vec![BNS_VERTEX::default()])
            .unwrap();
        let mut structure = StrFromINChI {
            at: atoms,
            num_atoms: 1,
            pSrm: restore_mode.as_const(),
            ..StrFromINChI::default()
        };
        let mut bns = BN_STRUCT {
            num_atoms: 1,
            num_vertices: 1,
            vert: vertices,
            ..BN_STRUCT::default()
        };
        let mut runs = 17;
        let mut delta = -19;
        assert_eq!(
            MakeSingleBondsMetal2ChargedHeteroat(
                &mut empty_heap,
                &mut bns,
                &mut BN_DATA::default(),
                &mut structure,
                atoms,
                copied,
                &mut [VAL_AT::default()],
                &mut absent_groups(),
                &mut runs,
                &mut delta,
                0x20,
                0,
            ),
            Ok(0)
        );
        assert_eq!((runs, delta), (17, -19));
        assert_eq!(structure.at, atoms);
        assert_eq!(
            empty_heap.slice(copied.as_const()).unwrap()[0].orig_at_number,
            7
        );

        let mut heap = SourceHeap::default();
        let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
        let mut metal = inp_ATOM {
            el_number: 11,
            valence: 1,
            orig_at_number: 1,
            ..inp_ATOM::default()
        };
        metal.neighbor[0] = 1;
        metal.bond_type[0] = BOND_TYPE_DOUBLE as u8;
        let mut nitrogen = inp_ATOM {
            el_number: 7,
            valence: 1,
            orig_at_number: 2,
            ..inp_ATOM::default()
        };
        nitrogen.neighbor[0] = 0;
        nitrogen.bond_type[0] = BOND_TYPE_DOUBLE as u8;
        let atoms = heap.allocate_model_storage(vec![metal, nitrogen]).unwrap();
        let copied = heap
            .allocate_model_storage(vec![inp_ATOM::default(); 2])
            .unwrap();
        let incident_0 = heap.allocate_model_storage(vec![0]).unwrap();
        let incident_1 = heap.allocate_model_storage(vec![0]).unwrap();
        let vertices = heap
            .allocate_model_storage(vec![
                BNS_VERTEX {
                    num_adj_edges: 1,
                    max_adj_edges: 1,
                    iedge: incident_0,
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    num_adj_edges: 1,
                    max_adj_edges: 1,
                    iedge: incident_1,
                    ..BNS_VERTEX::default()
                },
            ])
            .unwrap();
        let edges = heap
            .allocate_model_storage(vec![BNS_EDGE {
                neighbor1: 0,
                neighbor12: 1,
                neigh_ord: [0, 0],
                cap: 2,
                flow: 1,
                ..BNS_EDGE::default()
            }])
            .unwrap();
        let mut structure = StrFromINChI {
            at: atoms,
            num_atoms: 2,
            pSrm: restore_mode.as_const(),
            ..StrFromINChI::default()
        };
        let mut bns = BN_STRUCT {
            num_atoms: 2,
            num_vertices: 2,
            num_bonds: 1,
            num_edges: 1,
            vert: vertices,
            edge: edges,
            ..BN_STRUCT::default()
        };
        let mut valence = vec![
            VAL_AT {
                cMetal: 1,
                ..VAL_AT::default()
            },
            VAL_AT {
                cInitCharge: 1,
                cNumValenceElectrons: 5,
                cPeriodicRowNumber: 1,
                cnListIndex: 2,
                ..VAL_AT::default()
            },
        ];
        let before_allocations = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        runs = 23;
        delta = 29;
        assert_eq!(
            MakeSingleBondsMetal2ChargedHeteroat(
                &mut heap,
                &mut bns,
                &mut BN_DATA::default(),
                &mut structure,
                atoms,
                copied,
                &mut valence,
                &mut absent_groups(),
                &mut runs,
                &mut delta,
                0x40,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!((runs, delta), (23, 29));
        assert_eq!(structure.at, atoms);
        assert_eq!(heap.live_allocation_count(), before_allocations);
        assert_eq!(heap.source_allocation_calls(), 1);
        let generated = &heap.slice(copied.as_const()).unwrap()[1];
        assert_eq!(generated.charge, 1);
        assert_eq!(generated.bond_type[0], BOND_TYPE_DOUBLE as u8);
    }

    #[test]
    fn source_port__ichirvr4__movechargetoremovecenerpoints__line_1860() {
        fn fixture() -> (
            SourceHeap,
            StrFromINChI,
            SourceMutPointer<inp_ATOM>,
            SourceMutPointer<inp_ATOM>,
            BN_STRUCT,
            Vec<VAL_AT>,
        ) {
            let mut heap = SourceHeap::default();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM {
                    el_number: 6,
                    orig_at_number: 1,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            let copied = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let vertices = heap
                .allocate_model_storage(vec![BNS_VERTEX::default()])
                .unwrap();
            let structure = StrFromINChI {
                at: atoms,
                num_atoms: 1,
                pSrm: restore_mode.as_const(),
                ..StrFromINChI::default()
            };
            let bns = BN_STRUCT {
                num_atoms: 1,
                num_vertices: 1,
                vert: vertices,
                ..BN_STRUCT::default()
            };
            (
                heap,
                structure,
                atoms,
                copied,
                bns,
                vec![VAL_AT {
                    cNumValenceElectrons: 4,
                    ..VAL_AT::default()
                }],
            )
        }

        for successful_allocations in 0_u64..5 {
            let (mut heap, mut structure, atoms, copied, mut bns, mut valence) = fixture();
            let before_allocations = heap.live_allocation_count();
            heap.fail_after_allocations(successful_allocations);
            let mut runs = 17;
            let mut delta = -19;
            assert_eq!(
                MoveChargeToRemoveCenerpoints(
                    &mut heap,
                    &mut bns,
                    &mut BN_DATA::default(),
                    &mut structure,
                    atoms,
                    copied,
                    &mut valence,
                    &mut absent_groups(),
                    &mut runs,
                    &mut delta,
                    0x20,
                    0,
                ),
                Ok(crate::source_types::CT_OUT_OF_RAM),
                "allocation ordinal={successful_allocations}"
            );
            assert_eq!((runs, delta), (17, -19));
            assert_eq!(structure.at, atoms);
            assert_eq!(heap.live_allocation_count(), before_allocations);
            assert_eq!(bns.edge_forbidden_mask, 0);
        }

        let (mut heap, mut structure, atoms, copied, mut bns, mut valence) = fixture();
        let before_allocations = heap.live_allocation_count();
        heap.trace_source_allocations();
        let mut runs = 23;
        let mut delta = 29;
        assert_eq!(
            MoveChargeToRemoveCenerpoints(
                &mut heap,
                &mut bns,
                &mut BN_DATA::default(),
                &mut structure,
                atoms,
                copied,
                &mut valence,
                &mut absent_groups(),
                &mut runs,
                &mut delta,
                0x40,
                0,
            ),
            Ok(0)
        );
        assert_eq!((runs, delta), (23, 29));
        assert_eq!(structure.at, atoms);
        assert_eq!(heap.live_allocation_count(), before_allocations);
        assert_eq!(heap.source_allocation_calls(), 5);
        assert_eq!(bns.edge_forbidden_mask, BNS_EDGE_FORBIDDEN_TEST as i8);
        let copied_atom = &heap.slice(copied.as_const()).unwrap()[0];
        assert_eq!(copied_atom.orig_at_number, 1);
        assert_eq!(copied_atom.nRingSystem, 1);
    }

    #[test]
    fn source_port__ichirvr4__checkandrefixstereobonds__line_1678() {
        fn no_wrong_fixture() -> (
            SourceHeap,
            StrFromINChI,
            SourceMutPointer<inp_ATOM>,
            SourceMutPointer<inp_ATOM>,
            BN_STRUCT,
            Vec<VAL_AT>,
        ) {
            let mut heap = SourceHeap::default();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let atom = inp_ATOM {
                el_number: 6,
                orig_at_number: 1,
                ..inp_ATOM::default()
            };
            let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
            let copied = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let vertices = heap
                .allocate_model_storage(vec![BNS_VERTEX::default()])
                .unwrap();
            let structure = StrFromINChI {
                at: atoms,
                num_atoms: 1,
                pSrm: restore_mode.as_const(),
                ..StrFromINChI::default()
            };
            let bns = BN_STRUCT {
                num_atoms: 1,
                num_vertices: 1,
                vert: vertices,
                ..BN_STRUCT::default()
            };
            (heap, structure, atoms, copied, bns, vec![VAL_AT::default()])
        }

        fn wrong_fixture() -> (
            SourceHeap,
            StrFromINChI,
            SourceMutPointer<inp_ATOM>,
            SourceMutPointer<inp_ATOM>,
            BN_STRUCT,
            Vec<VAL_AT>,
        ) {
            let mut heap = SourceHeap::default();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let mut center = inp_ATOM {
                el_number: 6,
                valence: 3,
                orig_at_number: 1,
                ..inp_ATOM::default()
            };
            center.neighbor[..3].copy_from_slice(&[1, 2, 3]);
            center.bond_type[..4].fill(BOND_TYPE_SINGLE as u8);
            center.sb_parity[..3].copy_from_slice(&[1, 1, 0]);
            center.sb_ord[..2].copy_from_slice(&[1, 3]);
            let mut terminal_1 = inp_ATOM {
                el_number: 6,
                valence: 1,
                orig_at_number: 2,
                ..inp_ATOM::default()
            };
            terminal_1.neighbor[0] = 0;
            terminal_1.bond_type[0] = BOND_TYPE_SINGLE as u8;
            let mut terminal_2 = terminal_1.clone();
            terminal_2.orig_at_number = 3;
            let mut terminal_3 = terminal_1.clone();
            terminal_3.orig_at_number = 4;
            let atoms = heap
                .allocate_model_storage(vec![center, terminal_1, terminal_2, terminal_3])
                .unwrap();
            let copied = heap
                .allocate_model_storage(vec![inp_ATOM::default(); 4])
                .unwrap();
            let center_edges = heap.allocate_model_storage(vec![0, 1, 2]).unwrap();
            let terminal_1_edges = heap.allocate_model_storage(vec![0]).unwrap();
            let terminal_2_edges = heap.allocate_model_storage(vec![1]).unwrap();
            let terminal_3_edges = heap.allocate_model_storage(vec![2]).unwrap();
            let vertices = heap
                .allocate_model_storage(vec![
                    BNS_VERTEX {
                        st_edge: crate::source_types::BNS_ST_EDGE {
                            cap: 1,
                            flow: 1,
                            ..crate::source_types::BNS_ST_EDGE::default()
                        },
                        type_: crate::source_types::BNS_VERT_TYPE_ATOM as u16,
                        num_adj_edges: 3,
                        max_adj_edges: 3,
                        iedge: center_edges,
                        ..BNS_VERTEX::default()
                    },
                    BNS_VERTEX {
                        type_: crate::source_types::BNS_VERT_TYPE_ATOM as u16,
                        num_adj_edges: 1,
                        max_adj_edges: 1,
                        iedge: terminal_1_edges,
                        ..BNS_VERTEX::default()
                    },
                    BNS_VERTEX {
                        type_: crate::source_types::BNS_VERT_TYPE_ATOM as u16,
                        num_adj_edges: 1,
                        max_adj_edges: 1,
                        iedge: terminal_2_edges,
                        ..BNS_VERTEX::default()
                    },
                    BNS_VERTEX {
                        st_edge: crate::source_types::BNS_ST_EDGE {
                            cap: 1,
                            flow: 1,
                            ..crate::source_types::BNS_ST_EDGE::default()
                        },
                        type_: crate::source_types::BNS_VERT_TYPE_ATOM as u16,
                        num_adj_edges: 1,
                        max_adj_edges: 1,
                        iedge: terminal_3_edges,
                        ..BNS_VERTEX::default()
                    },
                ])
                .unwrap();
            let edges = heap
                .allocate_model_storage(vec![
                    BNS_EDGE {
                        neighbor1: 0,
                        neighbor12: 1,
                        neigh_ord: [0, 0],
                        flow: 0,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE {
                        neighbor1: 0,
                        neighbor12: 2,
                        neigh_ord: [1, 0],
                        flow: 0,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE {
                        neighbor1: 0,
                        neighbor12: 3,
                        neigh_ord: [2, 0],
                        cap: 1,
                        flow: 1,
                        ..BNS_EDGE::default()
                    },
                ])
                .unwrap();
            let structure = StrFromINChI {
                at: atoms,
                num_atoms: 4,
                pSrm: restore_mode.as_const(),
                ..StrFromINChI::default()
            };
            let bns = BN_STRUCT {
                num_atoms: 4,
                num_vertices: 4,
                num_bonds: 3,
                num_edges: 3,
                tot_st_cap: 2,
                tot_st_flow: 2,
                vert: vertices,
                edge: edges,
                ..BN_STRUCT::default()
            };
            (
                heap,
                structure,
                atoms,
                copied,
                bns,
                vec![VAL_AT::default(); 4],
            )
        }

        let (mut heap, mut structure, atoms, copied, mut bns, mut valence) = no_wrong_fixture();
        let before_allocations = heap.live_allocation_count();
        let mut runs = 17;
        let mut delta = -19;
        assert_eq!(
            CheckAndRefixStereobonds(
                &mut heap,
                &mut bns,
                &mut BN_DATA::default(),
                &mut structure,
                atoms,
                copied,
                &mut valence,
                &mut absent_groups(),
                &mut runs,
                &mut delta,
                0x20,
                0,
            ),
            Ok(0)
        );
        assert_eq!((runs, delta), (17, -19));
        assert_eq!(structure.at, atoms);
        assert_eq!(heap.slice(copied.as_const()).unwrap()[0].orig_at_number, 1);
        assert_eq!(heap.live_allocation_count(), before_allocations);

        for (successful_allocations, expected) in [
            (0_u64, RI_ERR_ALLOC),
            (1_u64, RI_ERR_ALLOC),
            (u64::MAX, RI_ERR_PROGR),
        ] {
            let (mut heap, mut structure, atoms, copied, mut bns, mut valence) = wrong_fixture();
            let before_allocations = heap.live_allocation_count();
            if successful_allocations != u64::MAX {
                heap.fail_after_allocations(successful_allocations);
            }
            let mut runs = 23;
            let mut delta = 29;
            assert_eq!(
                CheckAndRefixStereobonds(
                    &mut heap,
                    &mut bns,
                    &mut BN_DATA::default(),
                    &mut structure,
                    atoms,
                    copied,
                    &mut valence,
                    &mut absent_groups(),
                    &mut runs,
                    &mut delta,
                    0x40,
                    0,
                ),
                Ok(expected),
                "allocation ordinal={successful_allocations}"
            );
            assert_eq!((runs, delta), (23, 29));
            assert_eq!(structure.at, atoms);
            assert_eq!(heap.live_allocation_count(), before_allocations);
            assert_eq!(
                heap.slice(copied.as_const()).unwrap()[0].chem_bonds_valence,
                4
            );
            assert_eq!(heap.slice(bns.edge.as_const()).unwrap()[2].flow, 1);
            assert_eq!(heap.slice(bns.edge.as_const()).unwrap()[2].forbidden, 0);
        }

        let (mut heap, mut structure, atoms, copied, mut bns, mut valence) = wrong_fixture();
        heap.slice_mut(atoms).unwrap()[0].sb_ord[..2].copy_from_slice(&[0, 1]);
        heap.slice_mut(bns.edge).unwrap()[0].forbidden = BNS_EDGE_FORBIDDEN_MASK as i8;
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        bns.ic = clock;
        bns.max_altp = bns.altp.len() as i32;
        for path in &mut bns.altp {
            *path = heap
                .allocate_model_storage(vec![crate::source_types::BNS_ALT_PATH::default(); 16])
                .unwrap();
        }
        let slots = 16;
        let mut data = BN_DATA {
            BasePtr: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
            SwitchEdge: heap
                .allocate_model_storage(vec![[NO_VERTEX, -1]; slots])
                .unwrap(),
            Tree: heap.allocate_model_storage(vec![0_i8; slots]).unwrap(),
            ScanQ: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
            Pu: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
            Pv: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
            max_num_vertices: slots as i32,
            max_len_Pu_Pv: slots as i32,
            ..BN_DATA::default()
        };
        let before_allocations = heap.live_allocation_count();
        let mut runs = 31;
        let mut delta = 37;
        assert_eq!(
            CheckAndRefixStereobonds(
                &mut heap,
                &mut bns,
                &mut data,
                &mut structure,
                atoms,
                copied,
                &mut valence,
                &mut absent_groups(),
                &mut runs,
                &mut delta,
                0x40,
                0,
            ),
            Ok(0)
        );
        assert_eq!((runs, delta), (33, 37));
        assert_eq!(structure.at, atoms);
        assert_eq!(heap.live_allocation_count(), before_allocations);
        let edges = heap.slice(bns.edge.as_const()).unwrap();
        assert_eq!(edges[0].forbidden, BNS_EDGE_FORBIDDEN_MASK as i8);
        assert_eq!(edges[2].flow, 1);
        assert_eq!(edges[2].forbidden, 0);
        assert_eq!(bns.tot_st_flow, 0);
    }

    #[test]
    fn source_port__ichirvr4__normalizeandcompare__line_1274() {
        struct Fixture {
            heap: SourceHeap,
            clock: SourceMutPointer<INCHI_CLOCK>,
            bns: BN_STRUCT,
            structure: StrFromINChI,
            at: SourceMutPointer<inp_ATOM>,
            at2: SourceMutPointer<inp_ATOM>,
            at3: SourceMutPointer<inp_ATOM>,
            parameters: INPUT_PARMS,
            data: STRUCT_DATA,
            valence: Vec<VAL_AT>,
            groups: ALL_TC_GROUPS,
        }

        fn fixture() -> Fixture {
            let mut heap = SourceHeap::default();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let vertices = heap
                .allocate_model_storage(vec![BNS_VERTEX::default()])
                .unwrap();
            let bns = BN_STRUCT {
                num_atoms: 1,
                num_vertices: 1,
                vert: vertices,
                ..BN_STRUCT::default()
            };
            let mut carbon = inp_ATOM {
                el_number: 6,
                num_H: 4,
                orig_at_number: 1,
                component: 1,
                ..inp_ATOM::default()
            };
            carbon.elname[0] = b'C' as i8;
            let at = heap.allocate_model_storage(vec![carbon]).unwrap();
            let at2 = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let at3 = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            Fixture {
                heap,
                clock,
                bns,
                structure: StrFromINChI {
                    at,
                    num_atoms: 1,
                    bMobileH: TAUT_YES as i8,
                    iMobileH: TAUT_YES as i8,
                    pSrm: restore_mode.as_const(),
                    ..StrFromINChI::default()
                },
                at,
                at2,
                at3,
                parameters: INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                    bTautFlags: u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG),
                    ..INPUT_PARMS::default()
                },
                data: STRUCT_DATA::default(),
                valence: vec![VAL_AT::default()],
                groups: ALL_TC_GROUPS::default(),
            }
        }

        fn generate_input(fixture: &mut Fixture) -> SourceMutPointer<INChI> {
            assert_eq!(
                MakeOneInChIOutOfStrFromINChI2(
                    &mut fixture.heap,
                    &mut CANON_GLOBALS::default(),
                    fixture.clock,
                    &fixture.parameters,
                    &fixture.data,
                    &fixture.bns,
                    &mut fixture.structure,
                    fixture.at,
                    fixture.at2,
                    fixture.at3,
                    &fixture.valence,
                    &fixture.groups,
                    None,
                    None,
                    None,
                    0,
                ),
                Ok(0)
            );
            let input = fixture.structure.pOneINChI[0];
            assert!(!input.is_null());
            fixture.structure.pOneINChI[0] = SourceMutPointer::null();
            input
        }

        fn run(
            fixture: &mut Fixture,
            input: SourceMutPointer<INChI>,
        ) -> Result<i32, SourceHeapError> {
            let mut runs = 17;
            let mut delta = -19;
            let result = NormalizeAndCompare(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.clock,
                &fixture.parameters,
                &fixture.data,
                &mut fixture.bns,
                &mut BN_DATA::default(),
                &mut fixture.structure,
                fixture.at,
                fixture.at2,
                fixture.at3,
                &mut fixture.valence,
                &mut fixture.groups,
                [input, SourceMutPointer::null()],
                i64::MAX,
                0,
                &mut runs,
                &mut delta,
                4,
                0,
                0,
            );
            assert_eq!((runs, delta), (17, -19));
            result
        }

        fn assert_cleanup(fixture: &Fixture) {
            assert_eq!(
                fixture.structure.pOneINChI,
                [SourceMutPointer::null(), SourceMutPointer::null()]
            );
            assert_eq!(
                fixture.structure.pOneINChI_Aux,
                [SourceMutPointer::null(), SourceMutPointer::null()]
            );
            assert_eq!(
                fixture.structure.pOne_norm_data,
                [SourceMutPointer::null(), SourceMutPointer::null()]
            );
            assert_eq!(fixture.structure.One_ti, T_GROUP_INFO::default());
        }

        let mut equal = fixture();
        let equal_input = generate_input(&mut equal);
        assert_eq!(run(&mut equal, equal_input), Ok(0));
        assert_cleanup(&equal);
        assert!(!equal.heap.slice(equal_input.as_const()).unwrap().is_empty());

        let mut zz = fixture();
        let zz_input = generate_input(&mut zz);
        zz.structure.n_zy = 1;
        zz.structure.n_pzz = 1;
        assert_eq!(run(&mut zz, zz_input), Ok(0));
        assert_cleanup(&zz);

        let mut problem = fixture();
        let problem_input = generate_input(&mut problem);
        problem.heap.slice_mut(problem_input).unwrap()[0].nErrorCode = 7;
        assert_eq!(run(&mut problem, problem_input), Ok(RI_ERR_PROGR));
        assert_cleanup(&problem);

        let mut allocation = fixture();
        let allocation_input = generate_input(&mut allocation);
        allocation.heap.fail_after_allocations(0);
        assert_eq!(run(&mut allocation, allocation_input), Ok(RI_ERR_ALLOC));
        assert_cleanup(&allocation);
        assert_eq!(allocation.heap.source_allocation_calls(), 1);
    }

    #[test]
    fn source_policy__normalizeandcompare__undefined_initial_buffer_allocation_returns_structured_error()
     {
        let mut heap = SourceHeap::default();
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
        let vertices = heap
            .allocate_model_storage(vec![BNS_VERTEX::default()])
            .unwrap();
        let mut bns = BN_STRUCT {
            num_atoms: 1,
            num_vertices: 1,
            vert: vertices,
            ..BN_STRUCT::default()
        };
        let mut carbon = inp_ATOM {
            el_number: 6,
            num_H: 4,
            orig_at_number: 1,
            component: 1,
            ..inp_ATOM::default()
        };
        carbon.elname[0] = b'C' as i8;
        let at = heap.allocate_model_storage(vec![carbon]).unwrap();
        let at2 = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let at3 = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let formula = heap
            .allocate(b"CZzZz2\0".iter().copied().map(|byte| byte as i8).collect())
            .unwrap();
        let reversed = heap
            .allocate(vec![INChI {
                nErrorCode: 100,
                nNumberOfAtoms: 1,
                szHillFormula: formula,
                ..INChI::default()
            }])
            .unwrap();
        let original = heap
            .allocate_model_storage(vec![INChI {
                nErrorCode: 70,
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let parameters = INPUT_PARMS {
            nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
            ..INPUT_PARMS::default()
        };
        let data = STRUCT_DATA::default();
        let mut structure = StrFromINChI {
            at,
            num_atoms: 1,
            bMobileH: TAUT_YES as i8,
            iMobileH: TAUT_YES as i8,
            pSrm: restore_mode.as_const(),
            n_zy: 1,
            n_pzz: 1,
            pOneINChI: [reversed, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let mut expected_structure = structure.clone();
        expected_structure.pOneINChI = [SourceMutPointer::null(); TAUT_NUM as usize];
        let mut canonical_globals = CANON_GLOBALS::default();
        let mut bns_data = BN_DATA::default();
        let mut valence = vec![VAL_AT::default()];
        let mut groups = ALL_TC_GROUPS::default();
        let canonical_globals_before = canonical_globals.clone();
        let clock_before = heap.slice(clock.as_const()).unwrap()[0].clone();
        let bns_before = bns.clone();
        let bns_data_before = bns_data.clone();
        let vertices_before = heap.slice(vertices.as_const()).unwrap().to_vec();
        let atoms_before = heap.slice(at.as_const()).unwrap().to_vec();
        let atoms2_before = heap.slice(at2.as_const()).unwrap().to_vec();
        let atoms3_before = heap.slice(at3.as_const()).unwrap().to_vec();
        let source_allocations_before = heap.live_source_allocation_count();

        heap.trace_source_allocations();
        heap.fail_after_allocations(0);
        normalize_and_compare_test_begin_zz(false);
        let mut runs = 17;
        let mut delta = -19;
        let result = NormalizeAndCompare(
            &mut heap,
            &mut canonical_globals,
            clock,
            &parameters,
            &data,
            &mut bns,
            &mut bns_data,
            &mut structure,
            at,
            at2,
            at3,
            &mut valence,
            &mut groups,
            [original, SourceMutPointer::null()],
            i64::MAX,
            0,
            &mut runs,
            &mut delta,
            4,
            0,
            0,
        );
        let control = normalize_and_compare_test_finish();

        assert_eq!(result, Err(SourceHeapError::AllocationFailed));
        assert_eq!((runs, delta), (17, -19));
        assert_eq!(source_allocations_before, 2);
        assert_eq!(heap.source_allocation_calls(), 1);
        assert_eq!(heap.live_source_allocation_count(), 0);
        assert!(heap.slice(reversed.as_const()).is_err());
        assert!(heap.slice(formula.as_const()).is_err());
        assert_eq!(heap.slice(original.as_const()).unwrap()[0].nErrorCode, 70);
        assert_eq!(structure, expected_structure);
        assert_eq!(canonical_globals, canonical_globals_before);
        assert_eq!(heap.slice(clock.as_const()).unwrap()[0], clock_before);
        assert_eq!(bns, bns_before);
        assert_eq!(bns_data, bns_data_before);
        assert_eq!(heap.slice(vertices.as_const()).unwrap(), vertices_before);
        assert_eq!(heap.slice(at.as_const()).unwrap(), atoms_before);
        assert_eq!(heap.slice(at2.as_const()).unwrap(), atoms2_before);
        assert_eq!(heap.slice(at3.as_const()).unwrap(), atoms3_before);
        assert_eq!(valence, vec![VAL_AT::default()]);
        assert_eq!(groups, ALL_TC_GROUPS::default());
        assert_eq!(control.formula_before_cleanup, "CZzZz2");
        assert_eq!(
            control.events,
            [
                "MakeOneInChIOutOfStrFromINChI2",
                "inchi_strbuf_init:begin",
                "inchi_strbuf_init:non-positive",
                "MergeZzInHillFormula:begin",
                "MergeZzInHillFormula:zero",
                "inchi_strbuf_close:begin",
                "inchi_strbuf_close:end",
                "Free_INChI[0]:present",
                "free:formula",
                "free:pOneINChI[0]",
                "Free_INChI_Aux[0]:null",
                "FreeInpAtomData[0]:null",
                "Free_INChI[1]:null",
                "Free_INChI_Aux[1]:null",
                "FreeInpAtomData[1]:null",
                "free_t_group_info:null",
            ]
        );
        assert!(
            !control
                .events
                .iter()
                .any(|event| event.starts_with("CompareReversedINChI2") || event.starts_with("Fix"))
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__normalizeandcompare__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::{Value, json};

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--normalize-and-compare-records")
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
        let mut record_count = 0_usize;

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "NormalizeAndCompare");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            if official["family"] == "zz-zy-undefined" {
                assert_eq!(official["output"]["source_defined"], false, "{case_id}");
                assert!(
                    matches!(
                        official["output"]["termination_kind"].as_str(),
                        Some("signal" | "exit")
                    ),
                    "{case_id}"
                );
                assert!(official["output"].get("result").is_none(), "{case_id}");
                assert!(
                    official["output"].get("cleanup_events").is_none(),
                    "{case_id}"
                );
                record_count += 1;
                continue;
            }
            if official["family"] == "final" {
                let integers = |name: &str| -> Vec<i32> {
                    official["input"][name]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: {name} must be an array"))
                        .iter()
                        .map(|value| {
                            value
                                .as_i64()
                                .and_then(|value| i32::try_from(value).ok())
                                .unwrap_or_else(|| panic!("{case_id}: {name} value must fit i32"))
                        })
                        .collect()
                };
                let b_mobile_h = official["input"]["b_mobile_h"]
                    .as_i64()
                    .and_then(|value| i8::try_from(value).ok())
                    .expect("b_mobile_h must fit i8");
                let i_mobile_h = official["input"]["i_mobile_h"]
                    .as_i64()
                    .and_then(|value| i8::try_from(value).ok())
                    .expect("i_mobile_h must fit i8");
                let original_layer1 = official["input"]["original_layer1"]
                    .as_bool()
                    .expect("original_layer1 must be boolean");
                let reversed_layer1 = official["input"]["reversed_layer1"]
                    .as_bool()
                    .expect("reversed_layer1 must be boolean");
                let num_inp = official["input"]["num_inp"]
                    .as_i64()
                    .expect("num_inp must fit the GCC/Linux long model");
                let fixed_returns = integers("fixed_ret");
                let mobile_returns = integers("mobile_ret");
                let stereo_returns = integers("stereo_ret");
                let tc_edges = integers("tc_edges");
                let tinfo_h: Vec<u16> = integers("tinfo_h")
                    .into_iter()
                    .map(|value| {
                        u16::try_from(value)
                            .unwrap_or_else(|_| panic!("{case_id}: tinfo_h must fit AT_RANK"))
                    })
                    .collect();
                let tinfo_endpoints: Vec<u16> = integers("tinfo_endpoints")
                    .into_iter()
                    .map(|value| {
                        u16::try_from(value).unwrap_or_else(|_| {
                            panic!("{case_id}: tinfo_endpoints must fit AT_NUMB")
                        })
                    })
                    .collect();
                let tinfo_present = official["input"]["tinfo_present"]
                    .as_bool()
                    .expect("tinfo_present must be boolean");
                assert_eq!(tinfo_h.len(), tinfo_endpoints.len(), "{case_id}");
                assert_eq!(tinfo_present, !tinfo_h.is_empty(), "{case_id}");
                let primary_flags = official["input"]["primary_flags"]
                    .as_u64()
                    .expect("primary_flags must be unsigned");
                let primary_error = official["input"]["primary_error"]
                    .as_i64()
                    .and_then(|value| i32::try_from(value).ok())
                    .expect("primary_error must fit i32");
                let optional_flags = official["input"]["optional_flags"]
                    .as_u64()
                    .expect("optional_flags must be unsigned");
                let optional_error = official["input"]["optional_error"]
                    .as_i64()
                    .and_then(|value| i32::try_from(value).ok())
                    .expect("optional_error must fit i32");
                let optional_active =
                    b_mobile_h == TAUT_NON as i8 && (original_layer1 || reversed_layer1);
                let mut compare_script =
                    vec![(0, 0, 0, 0, 0), (primary_flags, primary_error, 0, 0, 0)];
                if optional_active {
                    compare_script.push((optional_flags, optional_error, 0, 0, 0));
                }
                let fill_returns = if i_mobile_h == TAUT_NON as i8 {
                    vec![0]
                } else {
                    Vec::new()
                };

                let mut heap = SourceHeap::default();
                let clock = heap
                    .allocate_model_storage(vec![INCHI_CLOCK::default()])
                    .unwrap();
                let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
                let vertices = heap
                    .allocate_model_storage(vec![BNS_VERTEX::default()])
                    .unwrap();
                let mut bns = BN_STRUCT {
                    num_atoms: 1,
                    num_vertices: 1,
                    vert: vertices,
                    ..BN_STRUCT::default()
                };
                let mut carbon = inp_ATOM {
                    el_number: 6,
                    num_H: 4,
                    orig_at_number: 1,
                    component: 1,
                    ..inp_ATOM::default()
                };
                carbon.elname[0] = b'C' as i8;
                let at = heap.allocate_model_storage(vec![carbon]).unwrap();
                let at2 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let at3 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let reversed_zero = heap
                    .allocate(vec![INChI {
                        nErrorCode: 100,
                        nNumberOfAtoms: 1,
                        ..INChI::default()
                    }])
                    .unwrap();
                let reversed_one = if reversed_layer1 {
                    heap.allocate(vec![INChI {
                        nErrorCode: 101,
                        nNumberOfAtoms: 1,
                        ..INChI::default()
                    }])
                    .unwrap()
                } else {
                    SourceMutPointer::null()
                };
                let norm_at = heap.allocate(vec![inp_ATOM::default()]).unwrap();
                let norm_holder = heap
                    .allocate(vec![INP_ATOM_DATA {
                        at: norm_at,
                        num_at: 30,
                        num_bonds: 40,
                        ..INP_ATOM_DATA::default()
                    }])
                    .unwrap();
                let original = heap
                    .allocate_model_storage(vec![
                        INChI {
                            nErrorCode: 70,
                            nNumberOfAtoms: 1,
                            ..INChI::default()
                        },
                        INChI {
                            nErrorCode: 71,
                            nNumberOfAtoms: usize::from(original_layer1) as i32,
                            ..INChI::default()
                        },
                    ])
                    .unwrap();
                let original_one = if original_layer1 {
                    original.offset(1).unwrap()
                } else {
                    SourceMutPointer::null()
                };
                let tc_group_values: Vec<TC_GROUP> = tc_edges
                    .iter()
                    .map(|num_edges| TC_GROUP {
                        num_edges: *num_edges,
                        ..TC_GROUP::default()
                    })
                    .collect();
                let tc_group_pointer = if tc_group_values.is_empty() {
                    SourceMutPointer::null()
                } else {
                    heap.allocate_model_storage(tc_group_values).unwrap()
                };
                let t_group_values: Vec<T_GROUP> = tinfo_h
                    .iter()
                    .copied()
                    .zip(tinfo_endpoints.iter().copied())
                    .map(|(mobile_h, endpoints)| T_GROUP {
                        num: [mobile_h, 0, 0, 0, 0],
                        nNumEndpoints: endpoints,
                        ..T_GROUP::default()
                    })
                    .collect();
                let t_group_pointer = if tinfo_present {
                    heap.allocate(t_group_values.clone()).unwrap()
                } else {
                    SourceMutPointer::null()
                };
                let parameters = INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                    ..INPUT_PARMS::default()
                };
                let data = STRUCT_DATA::default();
                let mut structure = StrFromINChI {
                    at,
                    num_atoms: 1,
                    bMobileH: b_mobile_h,
                    iMobileH: i_mobile_h,
                    pSrm: restore_mode.as_const(),
                    pOneINChI: [reversed_zero, reversed_one],
                    pOne_norm_data: [norm_holder, SourceMutPointer::null()],
                    One_ti: T_GROUP_INFO {
                        t_group: t_group_pointer,
                        num_t_groups: t_group_values.len() as i32,
                        ..T_GROUP_INFO::default()
                    },
                    ..StrFromINChI::default()
                };
                let mut expected_structure = structure.clone();
                expected_structure.pOneINChI = [SourceMutPointer::null(); TAUT_NUM as usize];
                expected_structure.pOne_norm_data = [SourceMutPointer::null(); TAUT_NUM as usize];
                expected_structure.One_ti = T_GROUP_INFO::default();
                let input = [original, original_one];
                let mut valence = vec![VAL_AT::default()];
                let mut groups = ALL_TC_GROUPS {
                    pTCG: tc_group_pointer,
                    num_tgroups: tc_edges.len() as i32,
                    ..ALL_TC_GROUPS::default()
                };
                let mut bns_data = BN_DATA::default();
                let source_allocations_before = heap.live_source_allocation_count();
                normalize_and_compare_test_begin_scripted(
                    vec![(0, true, true, tinfo_present)],
                    compare_script.clone(),
                    fill_returns.clone(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    fixed_returns.clone(),
                    mobile_returns.clone(),
                    stereo_returns.clone(),
                );
                normalize_and_compare_test_set_final_cleanup(
                    [1, i32::from(reversed_layer1)],
                    [norm_at, SourceMutPointer::null()],
                    tinfo_present.then_some(t_group_values),
                );
                let mut runs = 17;
                let mut delta = -19;
                let result = NormalizeAndCompare(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &parameters,
                    &data,
                    &mut bns,
                    &mut bns_data,
                    &mut structure,
                    at,
                    at2,
                    at3,
                    &mut valence,
                    &mut groups,
                    input,
                    num_inp,
                    0,
                    &mut runs,
                    &mut delta,
                    4,
                    0,
                    0,
                )
                .unwrap_or_else(|error| panic!("{case_id}: Rust source error: {error:?}"));
                let control = normalize_and_compare_test_finish();
                let source_allocations_after = heap.live_source_allocation_count();
                let actual = json!({
                    "result": result,
                    "runs": runs,
                    "delta": delta,
                    "rebuild_used": 1 - control.rebuild_script.len(),
                    "compare_used": compare_script.len() - control.compare_script.len(),
                    "fill_used": fill_returns.len() - control.fill_script.len(),
                    "fixed_used": fixed_returns.len() - control.fix_fixed_script.len(),
                    "mobile_used": mobile_returns.len() - control.fix_mobile_script.len(),
                    "stereo_used": stereo_returns.len() - control.fix_stereo_script.len(),
                    "complete_caller_state_exact": structure == expected_structure
                        && heap.slice(original.as_const()).unwrap()[0].nErrorCode == 70
                        && heap.slice(original.as_const()).unwrap()[1].nErrorCode == 71,
                    "prefree_state_exact": control.prefree_state_exact,
                    "holders_null": structure.pOneINChI
                        == [SourceMutPointer::null(); TAUT_NUM as usize]
                        && structure.pOneINChI_Aux
                            == [SourceMutPointer::null(); TAUT_NUM as usize]
                        && structure.pOne_norm_data
                            == [SourceMutPointer::null(); TAUT_NUM as usize],
                    "one_ti_cleared": structure.One_ti == T_GROUP_INFO::default(),
                    "input_pointers_preserved": input == [original, original_one],
                    "allocation_free_exact": source_allocations_after == 0,
                    "source_allocations_before": source_allocations_before,
                    "source_allocations_after": source_allocations_after,
                    "cleanup_events": control.events,
                });
                assert_eq!(actual, official["output"], "{case_id}");
                record_count += 1;
                continue;
            }
            if official["family"] == "first-less" || official["family"] == "more-extra" {
                let integers = |name: &str| -> Vec<i32> {
                    official["input"][name]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: {name} must be an array"))
                        .iter()
                        .map(|value| {
                            value
                                .as_i64()
                                .and_then(|value| i32::try_from(value).ok())
                                .unwrap_or_else(|| panic!("{case_id}: {name} value must fit i32"))
                        })
                        .collect()
                };
                let modes = |name: &str| -> Vec<INCHI_MODE> {
                    official["input"][name]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: {name} must be an array"))
                        .iter()
                        .map(|value| {
                            value
                                .as_u64()
                                .unwrap_or_else(|| panic!("{case_id}: {name} must be unsigned"))
                        })
                        .collect()
                };
                let mobile_h = official["input"]["mobile_h"]
                    .as_i64()
                    .and_then(|value| i8::try_from(value).ok())
                    .expect("mobile_h must fit i8");
                let rebuild_returns = integers("rebuild_ret");
                let rebuild_prep = integers("rebuild_prep");
                let rebuild_norm = integers("rebuild_norm");
                let compare_flags = modes("compare_flags");
                let compare_errors = integers("compare_error");
                let compare_h1 = integers("compare_h1");
                let compare_h2 = integers("compare_h2");
                let compare_endpoints = integers("compare_endpoints");
                let fill_returns = integers("fill_ret");
                let fix_returns = integers("fix_ret");
                assert_eq!(rebuild_returns.len(), rebuild_prep.len(), "{case_id}");
                assert_eq!(rebuild_returns.len(), rebuild_norm.len(), "{case_id}");
                assert_eq!(compare_flags.len(), compare_errors.len(), "{case_id}");
                assert_eq!(compare_flags.len(), compare_h1.len(), "{case_id}");
                assert_eq!(compare_flags.len(), compare_h2.len(), "{case_id}");
                assert_eq!(compare_flags.len(), compare_endpoints.len(), "{case_id}");
                let rebuild_script: Vec<_> = rebuild_returns
                    .iter()
                    .copied()
                    .zip(rebuild_prep.iter().map(|value| *value != 0))
                    .zip(rebuild_norm.iter().map(|value| *value != 0))
                    .map(|((result, prep), norm)| (result, prep, norm, false))
                    .collect();
                let compare_script: Vec<_> = compare_flags
                    .iter()
                    .copied()
                    .zip(compare_errors.iter().copied())
                    .zip(compare_h1.iter().copied())
                    .zip(compare_h2.iter().copied())
                    .zip(compare_endpoints.iter().copied())
                    .map(|((((flags, error), h1), h2), endpoints)| {
                        (flags, error, h1, h2, endpoints)
                    })
                    .collect();
                let repair_kind = official["input"]["repair_kind"]
                    .as_str()
                    .expect("repair_kind must be text");
                assert!(
                    matches!(repair_kind, "less" | "more" | "extra"),
                    "{case_id}: unknown repair kind {repair_kind}"
                );

                let mut heap = SourceHeap::default();
                let clock = heap
                    .allocate_model_storage(vec![INCHI_CLOCK::default()])
                    .unwrap();
                let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
                let vertices = heap
                    .allocate_model_storage(vec![BNS_VERTEX::default()])
                    .unwrap();
                let mut bns = BN_STRUCT {
                    num_atoms: 1,
                    num_vertices: 1,
                    vert: vertices,
                    ..BN_STRUCT::default()
                };
                let mut carbon = inp_ATOM {
                    el_number: 6,
                    num_H: 4,
                    orig_at_number: 1,
                    component: 1,
                    ..inp_ATOM::default()
                };
                carbon.elname[0] = b'C' as i8;
                let at = heap.allocate_model_storage(vec![carbon]).unwrap();
                let at2 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let at3 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let reversed = heap
                    .allocate(vec![INChI {
                        nErrorCode: 100,
                        nNumberOfAtoms: 1,
                        ..INChI::default()
                    }])
                    .unwrap();
                let original = heap
                    .allocate_model_storage(vec![INChI {
                        nErrorCode: 70,
                        nNumberOfAtoms: 1,
                        ..INChI::default()
                    }])
                    .unwrap();
                let parameters = INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                    ..INPUT_PARMS::default()
                };
                let data = STRUCT_DATA::default();
                let mut structure = StrFromINChI {
                    at,
                    num_atoms: 1,
                    bMobileH: TAUT_INI as i8,
                    iMobileH: mobile_h,
                    pSrm: restore_mode.as_const(),
                    pOneINChI: [reversed, SourceMutPointer::null()],
                    ..StrFromINChI::default()
                };
                let mut expected_structure = structure.clone();
                expected_structure.pOneINChI = [SourceMutPointer::null(); TAUT_NUM as usize];
                let input = [original, SourceMutPointer::null()];
                let mut valence = vec![VAL_AT::default()];
                let mut groups = ALL_TC_GROUPS::default();
                let mut bns_data = BN_DATA::default();
                normalize_and_compare_test_begin_scripted(
                    rebuild_script.clone(),
                    compare_script.clone(),
                    fill_returns.clone(),
                    if repair_kind == "less" {
                        fix_returns.clone()
                    } else {
                        Vec::new()
                    },
                    if repair_kind == "more" {
                        fix_returns.clone()
                    } else {
                        Vec::new()
                    },
                    if repair_kind == "extra" {
                        fix_returns.clone()
                    } else {
                        Vec::new()
                    },
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                );
                let mut runs = 17;
                let mut delta = -19;
                let result = NormalizeAndCompare(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &parameters,
                    &data,
                    &mut bns,
                    &mut bns_data,
                    &mut structure,
                    at,
                    at2,
                    at3,
                    &mut valence,
                    &mut groups,
                    input,
                    i64::MAX,
                    0,
                    &mut runs,
                    &mut delta,
                    4,
                    0,
                    0,
                )
                .unwrap_or_else(|error| panic!("{case_id}: Rust source error: {error:?}"));
                let control = normalize_and_compare_test_finish();
                let holders_null = structure.pOneINChI
                    == [SourceMutPointer::null(); TAUT_NUM as usize]
                    && structure.pOneINChI_Aux == [SourceMutPointer::null(); TAUT_NUM as usize]
                    && structure.pOne_norm_data == [SourceMutPointer::null(); TAUT_NUM as usize];
                let actual = json!({
                    "result": result,
                    "runs": runs,
                    "delta": delta,
                    "rebuild_used": rebuild_script.len() - control.rebuild_script.len(),
                    "compare_used": compare_script.len() - control.compare_script.len(),
                    "fill_used": fill_returns.len() - control.fill_script.len(),
                    "fix_used": match repair_kind {
                        "less" => fix_returns.len() - control.fix_less_script.len(),
                        "more" => fix_returns.len() - control.fix_more_script.len(),
                        "extra" => fix_returns.len() - control.fix_extra_script.len(),
                        _ => unreachable!(),
                    },
                    "complete_caller_state_exact": structure == expected_structure
                        && heap.slice(original.as_const()).unwrap()[0].nErrorCode == 70,
                    "prefree_state_exact": control.prefree_state_exact,
                    "holders_null": holders_null,
                    "allocation_free_exact": heap.live_source_allocation_count() == 0,
                    "cleanup_events": control.events,
                });
                assert_eq!(actual, official["output"], "{case_id}");
                record_count += 1;
                continue;
            }
            if official["family"] == "zz-zy" {
                let n_zy = official["input"]["n_zy"]
                    .as_i64()
                    .and_then(|value| i32::try_from(value).ok())
                    .expect("n_zy must fit i32");
                let n_pzz = official["input"]["n_pzz"]
                    .as_i64()
                    .and_then(|value| i32::try_from(value).ok())
                    .expect("n_pzz must fit i32");
                let formula_present = official["input"]["formula_present"]
                    .as_bool()
                    .expect("formula_present must be boolean");
                let formula = official["input"]["formula"]
                    .as_str()
                    .expect("formula must be text");
                let failure_stage = official["input"]["failure_stage"]
                    .as_u64()
                    .expect("failure_stage must be unsigned");
                let force_growth = official["input"]["force_growth"]
                    .as_bool()
                    .expect("force_growth must be boolean");

                let mut heap = SourceHeap::default();
                let clock = heap
                    .allocate_model_storage(vec![INCHI_CLOCK::default()])
                    .unwrap();
                let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
                let vertices = heap
                    .allocate_model_storage(vec![BNS_VERTEX::default()])
                    .unwrap();
                let mut bns = BN_STRUCT {
                    num_atoms: 1,
                    num_vertices: 1,
                    vert: vertices,
                    ..BN_STRUCT::default()
                };
                let mut carbon = inp_ATOM {
                    el_number: 6,
                    num_H: 4,
                    orig_at_number: 1,
                    component: 1,
                    ..inp_ATOM::default()
                };
                carbon.elname[0] = b'C' as i8;
                let at = heap.allocate_model_storage(vec![carbon]).unwrap();
                let at2 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let at3 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let formula_pointer = if formula_present {
                    heap.allocate(
                        formula
                            .bytes()
                            .chain(std::iter::once(0))
                            .map(|byte| byte as i8)
                            .collect(),
                    )
                    .unwrap()
                } else {
                    SourceMutPointer::null()
                };
                let reversed = heap
                    .allocate(vec![INChI {
                        nErrorCode: 100,
                        nNumberOfAtoms: 1,
                        szHillFormula: formula_pointer,
                        ..INChI::default()
                    }])
                    .unwrap();
                let original = heap
                    .allocate_model_storage(vec![INChI {
                        nErrorCode: 70,
                        nNumberOfAtoms: 1,
                        ..INChI::default()
                    }])
                    .unwrap();
                let parameters = INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                    ..INPUT_PARMS::default()
                };
                let data = STRUCT_DATA::default();
                let mut structure = StrFromINChI {
                    at,
                    num_atoms: 1,
                    bMobileH: TAUT_YES as i8,
                    iMobileH: TAUT_YES as i8,
                    pSrm: restore_mode.as_const(),
                    n_zy,
                    n_pzz,
                    pOneINChI: [reversed, SourceMutPointer::null()],
                    ..StrFromINChI::default()
                };
                let mut expected_structure = structure.clone();
                expected_structure.pOneINChI = [SourceMutPointer::null(); TAUT_NUM as usize];
                let mut valence = vec![VAL_AT::default()];
                let mut groups = ALL_TC_GROUPS::default();
                let mut canonical_globals = CANON_GLOBALS::default();
                let mut bns_data = BN_DATA::default();
                let canonical_globals_before = canonical_globals.clone();
                let clock_before = heap.slice(clock.as_const()).unwrap()[0].clone();
                let parameters_before = parameters.clone();
                let data_before = data.clone();
                let bns_before = bns.clone();
                let bns_data_before = bns_data.clone();
                let vertices_before = heap.slice(vertices.as_const()).unwrap().to_vec();
                let atoms_before = heap.slice(at.as_const()).unwrap().to_vec();
                let atoms2_before = heap.slice(at2.as_const()).unwrap().to_vec();
                let atoms3_before = heap.slice(at3.as_const()).unwrap().to_vec();
                let valence_before = valence.clone();
                let groups_before = groups.clone();
                let input = [original, SourceMutPointer::null()];
                let initial_source_allocations = heap.live_source_allocation_count();
                heap.trace_source_allocations();
                if failure_stage != 0 {
                    heap.fail_after_allocations(failure_stage - 1);
                }
                normalize_and_compare_test_begin_zz(force_growth);
                let mut runs = 17;
                let mut delta = -19;
                let result = NormalizeAndCompare(
                    &mut heap,
                    &mut canonical_globals,
                    clock,
                    &parameters,
                    &data,
                    &mut bns,
                    &mut bns_data,
                    &mut structure,
                    at,
                    at2,
                    at3,
                    &mut valence,
                    &mut groups,
                    input,
                    i64::MAX,
                    0,
                    &mut runs,
                    &mut delta,
                    4,
                    0,
                    0,
                )
                .unwrap_or_else(|error| panic!("{case_id}: Rust source error: {error:?}"));
                let control = normalize_and_compare_test_finish();
                let source_allocation_calls = heap.source_allocation_calls();
                let source_allocations_after = heap.live_source_allocation_count();
                let failure_observed = usize::from(failure_stage != 0);
                let successful_realloc = usize::from(force_growth && failure_stage != 4);
                let successful_allocations = initial_source_allocations
                    + usize::try_from(source_allocation_calls).unwrap()
                    - failure_observed
                    - successful_realloc;
                let complete_caller_state_exact = canonical_globals == canonical_globals_before
                    && heap.slice(clock.as_const()).unwrap()[0] == clock_before
                    && parameters == parameters_before
                    && data == data_before
                    && bns == bns_before
                    && bns_data == bns_data_before
                    && heap.slice(vertices.as_const()).unwrap() == vertices_before
                    && heap.slice(at.as_const()).unwrap() == atoms_before
                    && heap.slice(at2.as_const()).unwrap() == atoms2_before
                    && heap.slice(at3.as_const()).unwrap() == atoms3_before
                    && valence == valence_before
                    && groups == groups_before
                    && structure == expected_structure
                    && heap.slice(original.as_const()).unwrap()[0].nErrorCode == 70;
                let holders_null = structure.pOneINChI
                    == [SourceMutPointer::null(); TAUT_NUM as usize]
                    && structure.pOneINChI_Aux == [SourceMutPointer::null(); TAUT_NUM as usize]
                    && structure.pOne_norm_data == [SourceMutPointer::null(); TAUT_NUM as usize];
                let actual = json!({
                    "result": result,
                    "runs": runs,
                    "delta": delta,
                    "formula_before_cleanup": control.formula_before_cleanup,
                    "complete_caller_state_exact": complete_caller_state_exact,
                    "prefree_state_exact": control.prefree_state_exact,
                    "holders_null": holders_null,
                    "source_allocation_calls": source_allocation_calls,
                    "successful_allocations": successful_allocations,
                    "frees": successful_allocations - source_allocations_after,
                    "allocation_free_exact": source_allocations_after == 0,
                    "cleanup_events": control.events,
                });
                assert_eq!(actual, official["output"], "{case_id}");
                record_count += 1;
                continue;
            }
            if official["family"] == "layer-selection" || official["family"] == "common-success" {
                let common_success = official["family"] == "common-success";
                let mobile_h = official["input"]["mobile_h"]
                    .as_i64()
                    .and_then(|value| i8::try_from(value).ok())
                    .expect("mobile_h must fit i8");
                let original_state = official["input"]["original_state"]
                    .as_i64()
                    .and_then(|value| i32::try_from(value).ok())
                    .expect("original state must fit i32");
                let reversed_state = official["input"]["reversed_state"]
                    .as_i64()
                    .and_then(|value| i32::try_from(value).ok())
                    .expect("reversed state must fit i32");
                let mut heap = SourceHeap::default();
                let clock = heap
                    .allocate_model_storage(vec![INCHI_CLOCK::default()])
                    .unwrap();
                let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
                let vertices = heap
                    .allocate_model_storage(vec![BNS_VERTEX::default()])
                    .unwrap();
                let mut bns = BN_STRUCT {
                    num_atoms: 1,
                    num_vertices: 1,
                    vert: vertices,
                    ..BN_STRUCT::default()
                };
                let mut carbon = inp_ATOM {
                    el_number: 6,
                    num_H: 4,
                    orig_at_number: 1,
                    component: 1,
                    ..inp_ATOM::default()
                };
                carbon.elname[0] = b'C' as i8;
                let at = heap.allocate_model_storage(vec![carbon]).unwrap();
                let at2 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let at3 = heap
                    .allocate_model_storage(vec![inp_ATOM::default()])
                    .unwrap();
                let original = heap
                    .allocate_model_storage(vec![
                        INChI {
                            nErrorCode: 70,
                            nNumberOfAtoms: 1,
                            ..INChI::default()
                        },
                        INChI {
                            nErrorCode: 71,
                            nNumberOfAtoms: if original_state == 1 { 0 } else { 1 },
                            bDeleted: i32::from(original_state == 2),
                            ..INChI::default()
                        },
                    ])
                    .unwrap();
                let original_one = if original_state == 0 {
                    SourceMutPointer::null()
                } else {
                    original.offset(1).unwrap()
                };
                let reversed_zero = heap
                    .allocate(vec![INChI {
                        nErrorCode: 100,
                        nNumberOfAtoms: 1,
                        ..INChI::default()
                    }])
                    .unwrap();
                let reversed_one = if reversed_state == 0 {
                    SourceMutPointer::null()
                } else {
                    heap.allocate(vec![INChI {
                        nErrorCode: 101,
                        nNumberOfAtoms: if reversed_state == 1 { 0 } else { 1 },
                        bDeleted: i32::from(reversed_state == 2),
                        ..INChI::default()
                    }])
                    .unwrap()
                };
                let norm_zero = if common_success {
                    heap.allocate(vec![INP_ATOM_DATA {
                        num_at: 30,
                        num_bonds: 40,
                        ..INP_ATOM_DATA::default()
                    }])
                    .unwrap()
                } else {
                    SourceMutPointer::null()
                };
                let parameters = INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                    ..INPUT_PARMS::default()
                };
                let data = STRUCT_DATA::default();
                let mut structure = StrFromINChI {
                    at,
                    num_atoms: 1,
                    bMobileH: mobile_h,
                    iMobileH: TAUT_YES as i8,
                    pSrm: restore_mode.as_const(),
                    pOneINChI: [reversed_zero, reversed_one],
                    pOne_norm_data: [norm_zero, SourceMutPointer::null()],
                    ..StrFromINChI::default()
                };
                let mut expected_structure = structure.clone();
                expected_structure.pOneINChI = [SourceMutPointer::null(); TAUT_NUM as usize];
                expected_structure.pOne_norm_data = [SourceMutPointer::null(); TAUT_NUM as usize];
                let mut valence = vec![VAL_AT::default()];
                let mut groups = ALL_TC_GROUPS::default();
                let mut bns_data = BN_DATA::default();
                let input = [original, original_one];
                let source_allocations_before = heap.live_source_allocation_count();
                if common_success {
                    normalize_and_compare_test_begin_common_success();
                } else {
                    normalize_and_compare_test_begin_layer_selection(reversed_state);
                }
                let mut runs = 17;
                let mut delta = -19;
                let result = NormalizeAndCompare(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &parameters,
                    &data,
                    &mut bns,
                    &mut bns_data,
                    &mut structure,
                    at,
                    at2,
                    at3,
                    &mut valence,
                    &mut groups,
                    input,
                    i64::MAX,
                    0,
                    &mut runs,
                    &mut delta,
                    4,
                    0,
                    0,
                )
                .unwrap_or_else(|error| panic!("{case_id}: Rust source error: {error:?}"));
                let control = normalize_and_compare_test_finish();
                let source_allocations_after = heap.live_source_allocation_count();
                let selected_original_index =
                    i32::from(mobile_h == TAUT_NON as i8 && original_state == 3);
                let selected_reversed_index =
                    i32::from(mobile_h == TAUT_NON as i8 && reversed_state == 3);
                let actual = json!({
                    "result": result,
                    "runs": runs,
                    "delta": delta,
                    "selected_original_index": selected_original_index,
                    "selected_reversed_index": selected_reversed_index,
                    "complete_caller_state_exact": structure == expected_structure
                        && heap.slice(original.as_const()).unwrap()[0].nErrorCode == 70
                        && heap.slice(original.as_const()).unwrap()[1].nErrorCode == 71,
                    "prefree_state_exact": control.prefree_state_exact,
                    "allocation_free_exact": source_allocations_after == 0,
                    "source_allocations_before": source_allocations_before,
                    "source_allocations_after": source_allocations_after,
                    "holders_null": structure.pOneINChI
                        == [SourceMutPointer::null(); TAUT_NUM as usize]
                        && structure.pOneINChI_Aux
                            == [SourceMutPointer::null(); TAUT_NUM as usize]
                        && structure.pOne_norm_data
                            == [SourceMutPointer::null(); TAUT_NUM as usize],
                    "one_ti_cleared": structure.One_ti == T_GROUP_INFO::default(),
                    "input_pointers_preserved": input == [original, original_one],
                    "cleanup_events": control.events,
                });
                assert_eq!(actual, official["output"], "{case_id}");
                record_count += 1;
                continue;
            }
            let holder_mask = official["input"]["holder_mask"]
                .as_u64()
                .and_then(|value| u8::try_from(value).ok())
                .expect("holder mask must fit u8");
            let forced_return = official["input"]["forced_return"]
                .as_i64()
                .and_then(|value| i32::try_from(value).ok())
                .expect("forced return must fit i32");
            assert!(matches!(forced_return, RI_ERR_ALLOC | RI_ERR_PROGR));

            let mut heap = SourceHeap::default();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let vertices = heap
                .allocate_model_storage(vec![BNS_VERTEX::default()])
                .unwrap();
            let mut bns = BN_STRUCT {
                num_atoms: 1,
                num_vertices: 1,
                vert: vertices,
                ..BN_STRUCT::default()
            };
            let mut carbon = inp_ATOM {
                el_number: 6,
                num_H: 4,
                orig_at_number: 1,
                component: 1,
                ..inp_ATOM::default()
            };
            carbon.elname[0] = b'C' as i8;
            let at = heap.allocate_model_storage(vec![carbon]).unwrap();
            let at2 = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let at3 = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let parameters = INPUT_PARMS {
                nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                bTautFlags: u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG),
                ..INPUT_PARMS::default()
            };
            let data = STRUCT_DATA::default();
            let mut structure = StrFromINChI {
                at,
                num_atoms: 1,
                bMobileH: TAUT_YES as i8,
                iMobileH: TAUT_YES as i8,
                pSrm: restore_mode.as_const(),
                ..StrFromINChI::default()
            };
            for index in 0..TAUT_NUM as usize {
                let shift = index * 3;
                if holder_mask & (1 << shift) != 0 {
                    structure.pOneINChI[index] = heap
                        .allocate(vec![INChI {
                            nErrorCode: 100 + index as i32,
                            nNumberOfAtoms: 10 + index as i32,
                            ..INChI::default()
                        }])
                        .unwrap();
                }
                if holder_mask & (1 << (shift + 1)) != 0 {
                    structure.pOneINChI_Aux[index] = heap
                        .allocate(vec![INChI_Aux {
                            nErrorCode: 200 + index as i32,
                            nNumberOfAtoms: 20 + index as i32,
                            ..INChI_Aux::default()
                        }])
                        .unwrap();
                }
                if holder_mask & (1 << (shift + 2)) != 0 {
                    structure.pOne_norm_data[index] = heap
                        .allocate(vec![INP_ATOM_DATA {
                            num_at: 30 + index as i32,
                            num_bonds: 40 + index as i32,
                            ..INP_ATOM_DATA::default()
                        }])
                        .unwrap();
                }
            }
            structure.One_ti.t_group = heap
                .allocate(vec![T_GROUP {
                    nNumEndpoints: 51,
                    num: [52, 53, 0, 0, 0],
                    ..T_GROUP::default()
                }])
                .unwrap();
            structure.One_ti.num_t_groups = 1;

            let mut valence = vec![VAL_AT::default()];
            let mut groups = ALL_TC_GROUPS::default();
            let mut canonical_globals = CANON_GLOBALS::default();
            let mut bns_data = BN_DATA::default();
            let canonical_globals_before = canonical_globals.clone();
            let clock_before = heap.slice(clock.as_const()).unwrap()[0].clone();
            let parameters_before = parameters.clone();
            let data_before = data.clone();
            let bns_before = bns.clone();
            let bns_data_before = bns_data.clone();
            let vertices_before = heap.slice(vertices.as_const()).unwrap().to_vec();
            let atoms_before = heap.slice(at.as_const()).unwrap().to_vec();
            let atoms2_before = heap.slice(at2.as_const()).unwrap().to_vec();
            let atoms3_before = heap.slice(at3.as_const()).unwrap().to_vec();
            let valence_before = valence.clone();
            let groups_before = groups.clone();
            let mut expected_structure = structure.clone();
            expected_structure.pOneINChI = [SourceMutPointer::null(); TAUT_NUM as usize];
            expected_structure.pOneINChI_Aux = [SourceMutPointer::null(); TAUT_NUM as usize];
            expected_structure.pOne_norm_data = [SourceMutPointer::null(); TAUT_NUM as usize];
            expected_structure.One_ti = T_GROUP_INFO::default();
            let input = [SourceMutPointer::null(); TAUT_NUM as usize];
            let source_allocations_before = heap.live_source_allocation_count();
            normalize_and_compare_test_begin(holder_mask, forced_return);
            let mut runs = 17;
            let mut delta = -19;
            let result = NormalizeAndCompare(
                &mut heap,
                &mut canonical_globals,
                clock,
                &parameters,
                &data,
                &mut bns,
                &mut bns_data,
                &mut structure,
                at,
                at2,
                at3,
                &mut valence,
                &mut groups,
                input,
                i64::MAX,
                0,
                &mut runs,
                &mut delta,
                4,
                0,
                0,
            )
            .unwrap_or_else(|error| panic!("{case_id}: Rust source error: {error:?}"));
            let control = normalize_and_compare_test_finish();
            assert_eq!(control.holder_mask, holder_mask, "{case_id}");
            let complete_caller_state_exact = canonical_globals == canonical_globals_before
                && heap.slice(clock.as_const()).unwrap()[0] == clock_before
                && parameters == parameters_before
                && data == data_before
                && bns == bns_before
                && bns_data == bns_data_before
                && heap.slice(vertices.as_const()).unwrap() == vertices_before
                && heap.slice(at.as_const()).unwrap() == atoms_before
                && heap.slice(at2.as_const()).unwrap() == atoms2_before
                && heap.slice(at3.as_const()).unwrap() == atoms3_before
                && valence == valence_before
                && groups == groups_before
                && structure == expected_structure
                && input == [SourceMutPointer::null(); TAUT_NUM as usize];
            let holders_null = structure.pOneINChI == [SourceMutPointer::null(); TAUT_NUM as usize]
                && structure.pOneINChI_Aux == [SourceMutPointer::null(); TAUT_NUM as usize]
                && structure.pOne_norm_data == [SourceMutPointer::null(); TAUT_NUM as usize];
            let source_allocations_after = heap.live_source_allocation_count();
            let actual = json!({
                "result": result,
                "runs": runs,
                "delta": delta,
                "complete_caller_state_exact": complete_caller_state_exact,
                "prefree_state_exact": control.prefree_state_exact,
                "allocation_free_exact": source_allocations_after == 0,
                "source_allocations_before": source_allocations_before,
                "source_allocations_after": source_allocations_after,
                "holders_null": holders_null,
                "one_ti_cleared": structure.One_ti.t_group.is_null()
                    && structure.One_ti.num_t_groups == 0,
                "input_pointers_preserved": input
                    == [SourceMutPointer::null(); TAUT_NUM as usize],
                "cleanup_events": control.events,
            });
            assert_eq!(actual, official["output"], "{case_id}");
            record_count += 1;
        }

        assert_eq!(record_count, 241, "NormalizeAndCompare oracle case count");
    }

    #[test]
    fn source_port__ichirvr4__fixlesshydrogeninformula__line_177() {
        fn run(
            heap: &mut SourceHeap,
            atoms: SourceMutPointer<inp_ATOM>,
            restored: SourceMutPointer<inp_ATOM>,
            normalized: SourceMutPointer<inp_ATOM>,
            atom_count: i32,
            bns: &mut BN_STRUCT,
            valence: &mut [VAL_AT],
            runs: &mut i32,
            total_delta: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            FixLessHydrogenInFormula(
                heap,
                bns,
                &mut BN_DATA::default(),
                &StrFromINChI {
                    num_atoms: atom_count,
                    ..StrFromINChI::default()
                },
                atoms,
                restored,
                normalized,
                valence,
                &mut ALL_TC_GROUPS::default(),
                runs,
                total_delta,
                0x20,
                clock_t::default(),
            )
        }

        for atom_count in [-1, 0] {
            let mut heap = SourceHeap::default();
            heap.trace_source_allocations();
            let mut runs = i32::MIN;
            let mut total_delta = i32::MAX;
            assert_eq!(
                run(
                    &mut heap,
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    atom_count,
                    &mut BN_STRUCT::default(),
                    &mut [],
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(0),
                "atom count={atom_count}"
            );
            assert_eq!(heap.source_allocation_calls(), 0);
            assert_eq!(heap.live_allocation_count(), 0);
            assert_eq!((runs, total_delta), (i32::MIN, i32::MAX));
        }

        let mut failure_heap = SourceHeap::default();
        failure_heap.fail_after_allocations(0);
        let mut failure_runs = 5;
        let mut failure_delta = 7;
        assert_eq!(
            run(
                &mut failure_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                1,
                &mut BN_STRUCT::default(),
                &mut [VAL_AT::default()],
                &mut failure_runs,
                &mut failure_delta,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert_eq!(failure_heap.live_allocation_count(), 0);
        assert_eq!((failure_runs, failure_delta), (5, 7));

        let mut mask_heap = SourceHeap::default();
        let mask_atoms = mask_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mask_restored = mask_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mask_normalized = mask_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mask_edges = mask_heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    forbidden: 0x40,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    forbidden: 0x01,
                    ..BNS_EDGE::default()
                },
            ])
            .unwrap();
        let mut mask_bns = BN_STRUCT {
            edge: mask_edges,
            num_edges: 2,
            ..BN_STRUCT::default()
        };
        let mut mask_valence = vec![VAL_AT {
            nCMinusGroupEdge: 1,
            nCPlusGroupEdge: 2,
            ..VAL_AT::default()
        }];
        let fixture_allocations = mask_heap.live_allocation_count();
        mask_heap.trace_source_allocations();
        let mut mask_runs = 11;
        let mut mask_delta = 13;
        assert_eq!(
            run(
                &mut mask_heap,
                mask_atoms,
                mask_restored,
                mask_normalized,
                1,
                &mut mask_bns,
                &mut mask_valence,
                &mut mask_runs,
                &mut mask_delta,
            ),
            Ok(0)
        );
        assert_eq!(mask_heap.source_allocation_calls(), 1);
        assert_eq!(mask_heap.live_allocation_count(), fixture_allocations);
        assert_eq!(
            mask_heap
                .slice(mask_edges.as_const())
                .unwrap()
                .iter()
                .map(|edge| edge.forbidden)
                .collect::<Vec<_>>(),
            vec![0x40, 0x01]
        );
        assert_eq!((mask_runs, mask_delta), (11, 13));

        let mut clear_heap = SourceHeap::default();
        let original = clear_heap
            .allocate_model_storage(vec![inp_ATOM::default(); 2])
            .unwrap();
        let restored = clear_heap
            .allocate_model_storage(vec![
                inp_ATOM {
                    num_H: 1,
                    valence: 1,
                    neighbor: {
                        let mut neighbors = [0_u16; 20];
                        neighbors[0] = 1;
                        neighbors
                    },
                    cFlags: 7,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    cFlags: 9,
                    ..inp_ATOM::default()
                },
            ])
            .unwrap();
        let normalized = clear_heap
            .allocate_model_storage(vec![inp_ATOM::default(); 2])
            .unwrap();
        let mut clear_valence = vec![VAL_AT::default(); 2];
        clear_valence[0].cNumValenceElectrons = 5;
        clear_valence[0].cPeriodicRowNumber = 1;
        let fixture_allocations = clear_heap.live_allocation_count();
        let mut clear_runs = 17;
        let mut clear_delta = 19;
        assert_eq!(
            run(
                &mut clear_heap,
                original,
                restored,
                normalized,
                2,
                &mut BN_STRUCT::default(),
                &mut clear_valence,
                &mut clear_runs,
                &mut clear_delta,
            ),
            Ok(0)
        );
        assert_eq!(clear_heap.live_allocation_count(), fixture_allocations);
        assert_eq!(
            clear_heap
                .slice(restored.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.cFlags)
                .collect::<Vec<_>>(),
            vec![0, 0]
        );
        assert_eq!((clear_runs, clear_delta), (17, 19));

        let mut path_heap = SourceHeap::default();
        let original = path_heap
            .allocate_model_storage(vec![inp_ATOM::default(); 3])
            .unwrap();
        let mut restored_atoms = vec![inp_ATOM::default(); 3];
        restored_atoms[0].num_H = 1;
        restored_atoms[0].valence = 1;
        restored_atoms[0].neighbor[0] = 1;
        restored_atoms[0].bond_type[0] = BOND_TYPE_SINGLE as u8;
        restored_atoms[1].valence = 2;
        restored_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        restored_atoms[1].bond_type[..2].copy_from_slice(&[BOND_TYPE_SINGLE as u8, 2]);
        restored_atoms[2].charge = 1;
        restored_atoms[2].valence = 1;
        restored_atoms[2].neighbor[0] = 1;
        restored_atoms[2].bond_type[0] = 2;
        let restored = path_heap.allocate_model_storage(restored_atoms).unwrap();
        let mut normalized_atoms = vec![inp_ATOM::default(); 3];
        normalized_atoms[0].valence = 1;
        normalized_atoms[0].neighbor[0] = 1;
        normalized_atoms[0].bond_type[0] = 2;
        normalized_atoms[1].valence = 2;
        normalized_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        normalized_atoms[1].bond_type[..2].copy_from_slice(&[2, BOND_TYPE_SINGLE as u8]);
        normalized_atoms[2].valence = 1;
        normalized_atoms[2].neighbor[0] = 1;
        normalized_atoms[2].bond_type[0] = BOND_TYPE_SINGLE as u8;
        let normalized = path_heap.allocate_model_storage(normalized_atoms).unwrap();
        let path_edges = path_heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    flow: 1,
                    forbidden: 0x40,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    forbidden: 0x01,
                    ..BNS_EDGE::default()
                },
            ])
            .unwrap();
        let mut path_bns = BN_STRUCT {
            edge: path_edges,
            num_edges: 2,
            ..BN_STRUCT::default()
        };
        let mut path_valence = vec![VAL_AT::default(); 3];
        path_valence[0].cNumValenceElectrons = 5;
        path_valence[0].cPeriodicRowNumber = 1;
        path_valence[0].nCPlusGroupEdge = 1;
        path_valence[2].cNumValenceElectrons = 5;
        path_valence[2].cPeriodicRowNumber = 1;
        path_valence[2].nCPlusGroupEdge = 2;
        let fixture_allocations = path_heap.live_allocation_count();
        let before_edges = path_heap.slice(path_edges.as_const()).unwrap().to_vec();
        let mut path_runs = 23;
        let mut path_delta = 29;
        assert_eq!(
            run(
                &mut path_heap,
                original,
                restored,
                normalized,
                3,
                &mut path_bns,
                &mut path_valence,
                &mut path_runs,
                &mut path_delta,
            ),
            Ok(0)
        );
        assert_eq!(path_heap.live_allocation_count(), fixture_allocations);
        assert_eq!(
            path_heap.slice(path_edges.as_const()).unwrap(),
            before_edges
        );
        assert_eq!(
            path_heap
                .slice(restored.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.cFlags)
                .collect::<Vec<_>>(),
            vec![0, 0, 0]
        );
        assert_eq!((path_runs, path_delta), (23, 29));
    }

    #[test]
    fn source_port__ichirvr4__fixmorehydrogeninformula__line_398() {
        fn run(
            heap: &mut SourceHeap,
            restored: SourceMutPointer<inp_ATOM>,
            normalized: SourceMutPointer<inp_ATOM>,
            structure: &StrFromINChI,
            bns: &mut BN_STRUCT,
            valence: &mut [VAL_AT],
            runs: &mut i32,
            total_delta: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            FixMoreHydrogenInFormula(
                heap,
                bns,
                &mut BN_DATA::default(),
                structure,
                SourceMutPointer::null(),
                restored,
                normalized,
                valence,
                &mut ALL_TC_GROUPS::default(),
                runs,
                total_delta,
                0x20,
                clock_t::default(),
            )
        }

        for atom_count in [-1, 0] {
            let mut heap = SourceHeap::default();
            heap.trace_source_allocations();
            let mut runs = i32::MIN;
            let mut total_delta = i32::MAX;
            assert_eq!(
                run(
                    &mut heap,
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    &StrFromINChI {
                        num_atoms: atom_count,
                        ..StrFromINChI::default()
                    },
                    &mut BN_STRUCT::default(),
                    &mut [],
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(0)
            );
            assert_eq!(heap.source_allocation_calls(), 0);
            assert_eq!(heap.live_allocation_count(), 0);
            assert_eq!((runs, total_delta), (i32::MIN, i32::MAX));
        }

        let mut failure_heap = SourceHeap::default();
        failure_heap.fail_after_allocations(0);
        let mut failure_runs = 3;
        let mut failure_delta = 5;
        assert_eq!(
            run(
                &mut failure_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &StrFromINChI {
                    num_atoms: 1,
                    ..StrFromINChI::default()
                },
                &mut BN_STRUCT::default(),
                &mut [VAL_AT::default()],
                &mut failure_runs,
                &mut failure_delta,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert_eq!(failure_heap.live_allocation_count(), 0);
        assert_eq!((failure_runs, failure_delta), (3, 5));

        for endpoint_mode in [0_i8, 1, -1] {
            let mut heap = SourceHeap::default();
            let mut restored_atoms = vec![inp_ATOM::default(); 2];
            restored_atoms[0].neighbor[0] = 1;
            restored_atoms[0].charge = -1;
            restored_atoms[0].valence = 1;
            restored_atoms[0].chem_bonds_valence = 1;
            restored_atoms[0].endpoint = u16::from(endpoint_mode == 1);
            restored_atoms[1].neighbor[0] = 0;
            restored_atoms[1].valence = 1;
            restored_atoms[1].chem_bonds_valence = 2;
            let restored = heap.allocate_model_storage(restored_atoms).unwrap();
            let mut normalized_atoms = vec![inp_ATOM::default(); 2];
            normalized_atoms[0].num_H = 1;
            let normalized = heap.allocate_model_storage(normalized_atoms).unwrap();
            let fixed_endpoints = heap
                .allocate_model_storage(vec![u16::from(endpoint_mode == 0), 0])
                .unwrap();
            let edges = heap
                .allocate_model_storage(vec![BNS_EDGE {
                    flow: 1,
                    forbidden: 0x40,
                    ..BNS_EDGE::default()
                }])
                .unwrap();
            let mut bns = BN_STRUCT {
                edge: edges,
                num_edges: 1,
                ..BN_STRUCT::default()
            };
            let mut valence = vec![VAL_AT::default(); 2];
            valence[0].nCMinusGroupEdge = 1;
            valence[0].cNumValenceElectrons = 6;
            let structure = StrFromINChI {
                num_atoms: 2,
                bMobileH: endpoint_mode,
                endpoint: fixed_endpoints,
                ..StrFromINChI::default()
            };
            let fixture_allocations = heap.live_allocation_count();
            heap.trace_source_allocations();
            let mut runs = 7;
            let mut total_delta = 11;
            assert_eq!(
                run(
                    &mut heap,
                    restored,
                    normalized,
                    &structure,
                    &mut bns,
                    &mut valence,
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(0),
                "endpoint mode={endpoint_mode}"
            );
            assert_eq!(heap.source_allocation_calls(), 1);
            assert_eq!(heap.live_allocation_count(), fixture_allocations);
            assert_eq!(heap.slice(edges.as_const()).unwrap()[0].flow, 1);
            assert_eq!(heap.slice(edges.as_const()).unwrap()[0].forbidden, 0x40);
            assert_eq!((runs, total_delta), (7, 11));
        }

        let mut masks_heap = SourceHeap::default();
        let restored = masks_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let normalized = masks_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let endpoint = masks_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let edges = masks_heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    forbidden: 0x40,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    forbidden: 0x01,
                    ..BNS_EDGE::default()
                },
            ])
            .unwrap();
        let mut bns = BN_STRUCT {
            edge: edges,
            num_edges: 2,
            ..BN_STRUCT::default()
        };
        let mut valence = vec![VAL_AT {
            nCMinusGroupEdge: 1,
            nCPlusGroupEdge: 2,
            ..VAL_AT::default()
        }];
        let fixture_allocations = masks_heap.live_allocation_count();
        let mut runs = 13;
        let mut total_delta = 17;
        assert_eq!(
            run(
                &mut masks_heap,
                restored,
                normalized,
                &StrFromINChI {
                    num_atoms: 1,
                    endpoint,
                    ..StrFromINChI::default()
                },
                &mut bns,
                &mut valence,
                &mut runs,
                &mut total_delta,
            ),
            Ok(0)
        );
        assert_eq!(masks_heap.live_allocation_count(), fixture_allocations);
        assert_eq!(
            masks_heap
                .slice(edges.as_const())
                .unwrap()
                .iter()
                .map(|edge| edge.forbidden)
                .collect::<Vec<_>>(),
            vec![0x40, 0x01]
        );
        assert_eq!((runs, total_delta), (13, 17));
    }

    #[test]
    fn source_port__ichirvr4__fixremoveextratautendpoints__line_589() {
        fn run(
            heap: &mut SourceHeap,
            restored: SourceMutPointer<inp_ATOM>,
            normalized: SourceMutPointer<inp_ATOM>,
            atom_count: i32,
            bns: &mut BN_STRUCT,
            valence: &mut [VAL_AT],
            comparison: &ICR,
            runs: &mut i32,
            total_delta: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            FixRemoveExtraTautEndpoints(
                heap,
                bns,
                &mut BN_DATA::default(),
                &StrFromINChI {
                    num_atoms: atom_count,
                    ..StrFromINChI::default()
                },
                SourceMutPointer::null(),
                restored,
                SourceMutPointer::null(),
                normalized,
                valence,
                &mut ALL_TC_GROUPS::default(),
                comparison,
                runs,
                total_delta,
                0x20,
                clock_t::default(),
            )
        }

        for atom_count in [-1, 0] {
            let mut heap = SourceHeap::default();
            heap.trace_source_allocations();
            let mut runs = 2;
            let mut total_delta = 3;
            assert_eq!(
                run(
                    &mut heap,
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    atom_count,
                    &mut BN_STRUCT::default(),
                    &mut [],
                    &ICR::default(),
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(0)
            );
            assert_eq!(heap.source_allocation_calls(), 0);
            assert_eq!(heap.live_allocation_count(), 0);
            assert_eq!((runs, total_delta), (2, 3));
        }

        let mut failure_heap = SourceHeap::default();
        failure_heap.fail_after_allocations(0);
        let mut runs = 5;
        let mut total_delta = 7;
        assert_eq!(
            run(
                &mut failure_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                1,
                &mut BN_STRUCT::default(),
                &mut [VAL_AT::default()],
                &ICR::default(),
                &mut runs,
                &mut total_delta,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert_eq!(failure_heap.live_allocation_count(), 0);
        assert_eq!((runs, total_delta), (5, 7));

        for center_mode in [0_u8, 1, 2] {
            let mut heap = SourceHeap::default();
            let mut atoms = vec![inp_ATOM::default(); 2];
            atoms[0].valence = i8::from(center_mode != 0);
            atoms[0].chem_bonds_valence = if center_mode == 0 { 0 } else { 2 };
            atoms[0].neighbor[0] = 1;
            atoms[0].bond_type[0] = crate::source_types::BOND_TYPE_DOUBLE as u8;
            atoms[1].charge = i8::from(center_mode == 1);
            atoms[1].el_number = if center_mode == 2 { 6 } else { 1 };
            atoms[1].valence = 1;
            atoms[1].chem_bonds_valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].bond_type[0] = crate::source_types::BOND_TYPE_DOUBLE as u8;
            let restored = heap.allocate_model_storage(atoms).unwrap();
            let edges = heap
                .allocate_model_storage(vec![
                    BNS_EDGE {
                        forbidden: 0x40,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE {
                        forbidden: 0x01,
                        ..BNS_EDGE::default()
                    },
                ])
                .unwrap();
            let mut bns = BN_STRUCT {
                edge: edges,
                num_edges: 2,
                ..BN_STRUCT::default()
            };
            let mut valence = vec![VAL_AT::default(); 2];
            valence[0].nCMinusGroupEdge = 1;
            valence[1].nCPlusGroupEdge = 2;
            let mut comparison = ICR::default();
            comparison.num_endp_in1_only = 1;
            comparison.endp_in1_only[0] = 1;
            let fixture_allocations = heap.live_allocation_count();
            heap.trace_source_allocations();
            let mut runs = 11;
            let mut total_delta = 13;
            assert_eq!(
                run(
                    &mut heap,
                    restored,
                    SourceMutPointer::null(),
                    2,
                    &mut bns,
                    &mut valence,
                    &comparison,
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(0),
                "center mode={center_mode}"
            );
            assert_eq!(heap.source_allocation_calls(), 1);
            assert_eq!(heap.live_allocation_count(), fixture_allocations);
            assert_eq!(
                heap.slice(edges.as_const())
                    .unwrap()
                    .iter()
                    .map(|edge| edge.forbidden)
                    .collect::<Vec<_>>(),
                vec![0x40, 0x01]
            );
            assert_eq!((runs, total_delta), (11, 13));
        }
    }

    #[test]
    fn source_port__ichirvr4__forbidcarbonchargeedges__line_57() {
        let mut empty_heap = SourceHeap::default();
        let empty_edges = empty_heap
            .allocate_model_storage(vec![BNS_EDGE::default()])
            .unwrap();
        let empty_groups = absent_groups();
        let mut empty_list = EDGE_LIST::default();
        assert_eq!(
            ForbidCarbonChargeEdges(
                &mut empty_heap,
                &BN_STRUCT {
                    edge: empty_edges,
                    ..BN_STRUCT::default()
                },
                &empty_groups,
                &mut empty_list,
                0x04,
            ),
            Ok(0)
        );
        assert_eq!((empty_list.num_alloc, empty_list.num_edges), (2, 0));

        let mut heap = SourceHeap::default();
        let edges = heap
            .allocate_model_storage(vec![BNS_EDGE::default(); 4])
            .unwrap();
        let groups_pointer = heap
            .allocate_model_storage(vec![
                TC_GROUP {
                    nForwardEdge: 1,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    nForwardEdge: 2,
                    ..TC_GROUP::default()
                },
            ])
            .unwrap();
        let mut groups = absent_groups();
        groups.pTCG = groups_pointer;
        groups.nGroup[TCG_Plus_C0 as usize] = 0;
        groups.nGroup[TCG_Minus_C0 as usize] = 1;
        let network = BN_STRUCT {
            edge: edges,
            ..BN_STRUCT::default()
        };
        let mut list = EDGE_LIST::default();
        assert_eq!(
            ForbidCarbonChargeEdges(&mut heap, &network, &groups, &mut list, 0x05),
            Ok(2)
        );
        assert_eq!((list.num_alloc, list.num_edges), (2, 2));
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap(), [1, 2]);
        let after = heap.slice(edges.as_const()).unwrap();
        assert_eq!(
            (after[0].forbidden, after[1].forbidden, after[2].forbidden),
            (0, 5, 5)
        );

        heap.slice_mut(edges).unwrap()[1].forbidden = 7;
        assert_eq!(
            ForbidCarbonChargeEdges(&mut heap, &network, &groups, &mut list, 0x02),
            Ok(1)
        );
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap()[0], 2);
        assert_eq!(heap.slice(edges.as_const()).unwrap()[1].forbidden, 7);

        heap.slice_mut(groups_pointer).unwrap()[1].nForwardEdge = 0;
        assert_eq!(
            ForbidCarbonChargeEdges(&mut heap, &network, &groups, &mut list, 0x08),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(list.num_edges, 1);
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap()[0], 1);
        assert_eq!(heap.slice(edges.as_const()).unwrap()[1].forbidden, 15);

        let mut failure_heap = SourceHeap::default();
        let failure_edges = failure_heap
            .allocate_model_storage(vec![BNS_EDGE::default()])
            .unwrap();
        failure_heap.fail_after_allocations(0);
        let mut failure_list = EDGE_LIST::default();
        assert_eq!(
            ForbidCarbonChargeEdges(
                &mut failure_heap,
                &BN_STRUCT {
                    edge: failure_edges,
                    ..BN_STRUCT::default()
                },
                &absent_groups(),
                &mut failure_list,
                i32::MIN,
            ),
            Ok(crate::source_types::RI_ERR_ALLOC)
        );
        assert!(failure_list.pnEdges.is_null());
        assert_eq!((failure_list.num_alloc, failure_list.num_edges), (0, 0));
    }

    #[test]
    fn source_port__ichirvr4__forbidnintrogenplus2bondsinsmallrings__line_116() {
        fn run_case(
            atom: inp_ATOM,
            valence: VAL_AT,
            min_ring_size: i32,
            initial_forbidden: i8,
            fail_allocation: bool,
        ) -> (i32, i8, i32, i32, Option<i32>) {
            let mut heap = SourceHeap::default();
            let edge_pointer = heap
                .allocate_model_storage(vec![BNS_EDGE {
                    forbidden: initial_forbidden,
                    ..BNS_EDGE::default()
                }])
                .unwrap();
            if fail_allocation {
                heap.fail_after_allocations(0);
            }
            let mut network = BN_STRUCT {
                edge: edge_pointer,
                ..BN_STRUCT::default()
            };
            let mut list = EDGE_LIST::default();
            let result = ForbidNintrogenPlus2BondsInSmallRings(
                &mut heap,
                &mut network,
                &[atom],
                1,
                &[valence],
                min_ring_size,
                &ALL_TC_GROUPS::default(),
                &mut list,
                0x20,
            )
            .unwrap();
            let first = if list.num_edges > 0 {
                Some(heap.slice(list.pnEdges.as_const()).unwrap()[0])
            } else {
                None
            };
            (
                result,
                heap.slice(edge_pointer.as_const()).unwrap()[0].forbidden,
                list.num_alloc,
                list.num_edges,
                first,
            )
        }

        let atom = inp_ATOM {
            valence: 2,
            ..inp_ATOM::default()
        };
        let valence = VAL_AT {
            cNumValenceElectrons: 5,
            cPeriodicRowNumber: 1,
            nCPlusGroupEdge: 1,
            cnListIndex: 5,
            cMinRingSize: 6,
            ..VAL_AT::default()
        };

        assert_eq!(
            run_case(atom.clone(), valence.clone(), 6, 2, false),
            (0, 34, 128, 1, Some(0))
        );
        assert_eq!(
            run_case(atom.clone(), valence.clone(), 6, 34, false),
            (0, 34, 0, 0, None)
        );
        assert_eq!(
            run_case(atom.clone(), valence.clone(), 6, 2, true),
            (crate::source_types::RI_ERR_ALLOC, 34, 0, 0, None)
        );

        let mut rejected = Vec::new();
        let mut changed_atom = atom.clone();
        changed_atom.valence = 1;
        rejected.push((changed_atom, valence.clone(), 6));
        let mut changed_atom = atom.clone();
        changed_atom.num_H = 1;
        rejected.push((changed_atom, valence.clone(), 6));
        let mut changed_atom = atom.clone();
        changed_atom.endpoint = 1;
        rejected.push((changed_atom, valence.clone(), 6));

        let mut changed_valence = valence.clone();
        changed_valence.cNumValenceElectrons = 4;
        rejected.push((atom.clone(), changed_valence, 6));
        let mut changed_valence = valence.clone();
        changed_valence.cPeriodicRowNumber = 2;
        rejected.push((atom.clone(), changed_valence, 6));
        let mut changed_valence = valence.clone();
        changed_valence.cMaxFlowToMetal = 1;
        rejected.push((atom.clone(), changed_valence, 6));
        let mut changed_valence = valence.clone();
        changed_valence.nCPlusGroupEdge = 0;
        rejected.push((atom.clone(), changed_valence, 6));
        let mut changed_valence = valence.clone();
        changed_valence.cnListIndex = 0;
        rejected.push((atom.clone(), changed_valence, 6));
        let mut changed_valence = valence.clone();
        changed_valence.cnListIndex = 1;
        rejected.push((atom.clone(), changed_valence, 6));
        let mut changed_valence = valence.clone();
        changed_valence.cMinRingSize = 0;
        rejected.push((atom.clone(), changed_valence, 6));
        let mut changed_valence = valence.clone();
        changed_valence.cMinRingSize = 7;
        rejected.push((atom, changed_valence, 6));

        for (rejected_atom, rejected_valence, threshold) in rejected {
            assert_eq!(
                run_case(rejected_atom, rejected_valence, threshold, 2, false),
                (0, 2, 0, 0, None)
            );
        }
    }

    #[test]
    fn source_port__ichirvr4__filloutextrafixedhdatainchi__line_831() {
        fn inchi_pair(
            heap: &mut SourceHeap,
            fixed_h: Option<Vec<i8>>,
            tautomer: Option<Vec<AT_NUMB>>,
            atom_count: i32,
        ) -> [SourceMutPointer<INChI>; 2] {
            let fixed_h = fixed_h
                .map(|values| heap.allocate_model_storage(values).unwrap())
                .unwrap_or_default();
            let tautomer_length = tautomer.as_ref().map_or(0, |values| values.len() as i32);
            let tautomer = tautomer
                .map(|values| heap.allocate_model_storage(values).unwrap())
                .unwrap_or_default();
            let fixed = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: atom_count,
                    nNum_H_fixed: fixed_h,
                    ..INChI::default()
                }])
                .unwrap();
            let mobile = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: atom_count,
                    lenTautomer: tautomer_length,
                    nTautomer: tautomer,
                    ..INChI::default()
                }])
                .unwrap();
            [fixed, mobile]
        }

        let mut allocated_heap = SourceHeap::default();
        let allocated_inchi = inchi_pair(&mut allocated_heap, Some(vec![3, -2, 1]), None, 3);
        let mut allocated = StrFromINChI {
            num_atoms: 3,
            ..StrFromINChI::default()
        };
        allocated_heap.trace_source_allocations();
        assert_eq!(
            FillOutExtraFixedHDataInChI(&mut allocated_heap, &mut allocated, allocated_inchi,),
            Ok(0)
        );
        assert_eq!(allocated_heap.source_allocation_calls(), 2);
        assert_eq!(
            allocated_heap.slice(allocated.endpoint.as_const()).unwrap(),
            [0, 0, 0]
        );
        assert_eq!(
            allocated_heap.slice(allocated.fixed_H.as_const()).unwrap(),
            [3, -2, 1]
        );

        let mut reused_heap = SourceHeap::default();
        let reused_inchi = inchi_pair(&mut reused_heap, None, None, 3);
        let old_endpoint = reused_heap
            .allocate_model_storage(vec![7_u16, 8, 9])
            .unwrap();
        let old_fixed = reused_heap
            .allocate_model_storage(vec![4_i8, 5, 6])
            .unwrap();
        let mut reused = StrFromINChI {
            endpoint: old_endpoint,
            fixed_H: old_fixed,
            num_atoms: 3,
            ..StrFromINChI::default()
        };
        reused_heap.trace_source_allocations();
        assert_eq!(
            FillOutExtraFixedHDataInChI(&mut reused_heap, &mut reused, reused_inchi),
            Ok(0)
        );
        assert_eq!(reused.endpoint, old_endpoint);
        assert_eq!(reused.fixed_H, old_fixed);
        assert_eq!(reused_heap.source_allocation_calls(), 0);
        assert_eq!(
            reused_heap.slice(old_endpoint.as_const()).unwrap(),
            [0, 0, 0]
        );
        assert_eq!(reused_heap.slice(old_fixed.as_const()).unwrap(), [0, 0, 0]);

        let mut first_failure_heap = SourceHeap::default();
        let first_failure_inchi = inchi_pair(&mut first_failure_heap, None, None, 2);
        let mut first_failure = StrFromINChI {
            num_atoms: 2,
            ..StrFromINChI::default()
        };
        first_failure_heap.fail_after_allocations(0);
        assert_eq!(
            FillOutExtraFixedHDataInChI(
                &mut first_failure_heap,
                &mut first_failure,
                first_failure_inchi,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(first_failure.endpoint.is_null());
        assert!(!first_failure.fixed_H.is_null());
        assert_eq!(first_failure_heap.source_allocation_calls(), 2);
        assert_eq!(
            first_failure_heap
                .slice(first_failure.fixed_H.as_const())
                .unwrap(),
            [0, 0]
        );

        let mut second_failure_heap = SourceHeap::default();
        let second_failure_inchi = inchi_pair(&mut second_failure_heap, None, None, 2);
        let mut second_failure = StrFromINChI {
            num_atoms: 2,
            ..StrFromINChI::default()
        };
        second_failure_heap.fail_after_allocations(1);
        assert_eq!(
            FillOutExtraFixedHDataInChI(
                &mut second_failure_heap,
                &mut second_failure,
                second_failure_inchi,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(!second_failure.endpoint.is_null());
        assert!(second_failure.fixed_H.is_null());
        assert_eq!(second_failure_heap.source_allocation_calls(), 2);
        assert_eq!(
            second_failure_heap
                .slice(second_failure.endpoint.as_const())
                .unwrap(),
            [0, 0]
        );

        let mut tautomer_heap = SourceHeap::default();
        let tautomer_inchi = inchi_pair(
            &mut tautomer_heap,
            Some(vec![-3, 4]),
            Some(vec![1, 3, 2, 1, 2]),
            2,
        );
        let endpoint = tautomer_heap
            .allocate_model_storage(vec![11_u16, 12])
            .unwrap();
        let fixed = tautomer_heap
            .allocate_model_storage(vec![13_i8, 14])
            .unwrap();
        let mut tautomer = StrFromINChI {
            endpoint,
            fixed_H: fixed,
            num_atoms: 2,
            ..StrFromINChI::default()
        };
        tautomer_heap.trace_source_allocations();
        assert_eq!(
            FillOutExtraFixedHDataInChI(&mut tautomer_heap, &mut tautomer, tautomer_inchi),
            Ok(0)
        );
        assert_eq!(tautomer_heap.source_allocation_calls(), 3);
        assert_eq!(tautomer_heap.slice(endpoint.as_const()).unwrap(), [0, 1]);
        assert_eq!(tautomer_heap.slice(fixed.as_const()).unwrap(), [-3, 4]);
        assert_eq!(tautomer.ti.num_t_groups, 1);
        assert_eq!(tautomer.ti.nNumEndpoints, 1);
    }

    #[test]
    fn source_port__ichirvr4__filloutextrafixedhdatarestr__line_754() {
        fn aux(
            heap: &mut SourceHeap,
            ordinary: Option<Vec<AT_NUMB>>,
            isotopic: Option<Vec<AT_NUMB>>,
        ) -> SourceMutPointer<INChI_Aux> {
            let ordinary = ordinary
                .map(|values| heap.allocate_model_storage(values).unwrap())
                .unwrap_or_default();
            let isotopic = isotopic
                .map(|values| heap.allocate_model_storage(values).unwrap())
                .unwrap_or_default();
            heap.allocate_model_storage(vec![INChI_Aux {
                nOrigAtNosInCanonOrd: ordinary,
                nIsotopicOrigAtNosInCanonOrd: isotopic,
                ..INChI_Aux::default()
            }])
            .unwrap()
        }

        let mut missing_heap = SourceHeap::default();
        let old_forward = missing_heap
            .allocate_model_storage(vec![41_u16, 42])
            .unwrap();
        let old_inverse = missing_heap
            .allocate_model_storage(vec![43_u16, 44])
            .unwrap();
        let mut missing = StrFromINChI {
            num_atoms: 2,
            nCanon2Atno: [old_forward, SourceMutPointer::null()],
            nAtno2Canon: [old_inverse, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        assert_eq!(
            FillOutExtraFixedHDataRestr(&mut missing_heap, &mut missing),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(missing.nCanon2Atno[0], old_forward);
        assert_eq!(missing.nAtno2Canon[0], old_inverse);
        assert_eq!(
            missing_heap.slice(old_forward.as_const()).unwrap(),
            [41, 42]
        );
        assert_eq!(
            missing_heap.slice(old_inverse.as_const()).unwrap(),
            [43, 44]
        );

        let mut preferred_heap = SourceHeap::default();
        let preferred_aux = aux(
            &mut preferred_heap,
            Some(vec![1, 2, 3]),
            Some(vec![3, 1, 2]),
        );
        let old_second_forward = preferred_heap
            .allocate_model_storage(vec![71_u16, 72, 73])
            .unwrap();
        let old_second_inverse = preferred_heap
            .allocate_model_storage(vec![81_u16, 82, 83])
            .unwrap();
        let mut preferred = StrFromINChI {
            num_atoms: 3,
            pOneINChI_Aux: [preferred_aux, SourceMutPointer::null()],
            nCanon2Atno: [SourceMutPointer::null(), old_second_forward],
            nAtno2Canon: [SourceMutPointer::null(), old_second_inverse],
            ..StrFromINChI::default()
        };
        assert_eq!(
            FillOutExtraFixedHDataRestr(&mut preferred_heap, &mut preferred),
            Ok(0)
        );
        assert_eq!(
            preferred_heap
                .slice(preferred.nCanon2Atno[0].as_const())
                .unwrap(),
            [2, 0, 1]
        );
        assert_eq!(
            preferred_heap
                .slice(preferred.nAtno2Canon[0].as_const())
                .unwrap(),
            [1, 2, 0]
        );
        assert!(preferred.nCanon2Atno[1].is_null());
        assert!(preferred.nAtno2Canon[1].is_null());
        assert_eq!(
            preferred_heap.slice(old_second_forward.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            preferred_heap.slice(old_second_inverse.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut fallback_heap = SourceHeap::default();
        let fallback_aux_0 = aux(&mut fallback_heap, Some(vec![2, 3, 1]), Some(vec![0, 1, 2]));
        let fallback_aux_1 = aux(&mut fallback_heap, Some(vec![3, 1, 2]), None);
        let reused_forward = fallback_heap
            .allocate_model_storage(vec![90_u16, 91, 92])
            .unwrap();
        let reused_inverse = fallback_heap
            .allocate_model_storage(vec![93_u16, 94, 95])
            .unwrap();
        let mut fallback = StrFromINChI {
            num_atoms: 3,
            pOneINChI_Aux: [fallback_aux_0, fallback_aux_1],
            nCanon2Atno: [reused_forward, SourceMutPointer::null()],
            nAtno2Canon: [reused_inverse, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        fallback_heap.trace_source_allocations();
        assert_eq!(
            FillOutExtraFixedHDataRestr(&mut fallback_heap, &mut fallback),
            Ok(0)
        );
        assert_eq!(fallback.nCanon2Atno[0], reused_forward);
        assert_eq!(fallback.nAtno2Canon[0], reused_inverse);
        assert_eq!(fallback_heap.source_allocation_calls(), 2);
        assert_eq!(
            fallback_heap.slice(reused_forward.as_const()).unwrap(),
            [1, 2, 0]
        );
        assert_eq!(
            fallback_heap.slice(reused_inverse.as_const()).unwrap(),
            [2, 0, 1]
        );
        assert_eq!(
            fallback_heap
                .slice(fallback.nCanon2Atno[1].as_const())
                .unwrap(),
            [2, 0, 1]
        );
        assert_eq!(
            fallback_heap
                .slice(fallback.nAtno2Canon[1].as_const())
                .unwrap(),
            [1, 2, 0]
        );

        let mut zero_heap = SourceHeap::default();
        let zero_aux = aux(&mut zero_heap, Some(vec![0, 1]), None);
        let mut zero_list = StrFromINChI {
            num_atoms: 2,
            pOneINChI_Aux: [zero_aux, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        assert_eq!(
            FillOutExtraFixedHDataRestr(&mut zero_heap, &mut zero_list),
            Ok(RI_ERR_PROGR)
        );
        assert!(zero_list.nCanon2Atno[0].is_null());

        let mut first_failure_heap = SourceHeap::default();
        let first_failure_aux = aux(&mut first_failure_heap, Some(vec![1, 2]), None);
        let mut first_failure = StrFromINChI {
            num_atoms: 2,
            pOneINChI_Aux: [first_failure_aux, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        first_failure_heap.fail_after_allocations(0);
        assert_eq!(
            FillOutExtraFixedHDataRestr(&mut first_failure_heap, &mut first_failure),
            Ok(RI_ERR_ALLOC)
        );
        assert!(first_failure.nCanon2Atno[0].is_null());
        assert!(first_failure.nAtno2Canon[0].is_null());

        let mut second_failure_heap = SourceHeap::default();
        let second_failure_aux = aux(&mut second_failure_heap, Some(vec![2, 1]), None);
        let mut second_failure = StrFromINChI {
            num_atoms: 2,
            pOneINChI_Aux: [second_failure_aux, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        second_failure_heap.fail_after_allocations(1);
        assert_eq!(
            FillOutExtraFixedHDataRestr(&mut second_failure_heap, &mut second_failure),
            Ok(RI_ERR_ALLOC)
        );
        assert!(!second_failure.nCanon2Atno[0].is_null());
        assert!(second_failure.nAtno2Canon[0].is_null());
        assert_eq!(
            second_failure_heap
                .slice(second_failure.nCanon2Atno[0].as_const())
                .unwrap(),
            [0, 0]
        );

        let mut invalid_heap = SourceHeap::default();
        let invalid_aux = aux(&mut invalid_heap, Some(vec![0, 1]), Some(vec![1, 0]));
        let valid_second_aux = aux(&mut invalid_heap, Some(vec![1, 2]), None);
        let mut invalid = StrFromINChI {
            num_atoms: 2,
            pOneINChI_Aux: [invalid_aux, valid_second_aux],
            ..StrFromINChI::default()
        };
        assert_eq!(
            FillOutExtraFixedHDataRestr(&mut invalid_heap, &mut invalid),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            invalid_heap
                .slice(invalid.nCanon2Atno[0].as_const())
                .unwrap(),
            [0, AT_NUMB::MAX]
        );
        assert_eq!(
            invalid_heap
                .slice(invalid.nAtno2Canon[0].as_const())
                .unwrap(),
            [0, 0]
        );
        assert!(invalid.nCanon2Atno[1].is_null());
    }

    #[test]
    fn source_port__ichirvr4__filloutcmp2fhinchi__line_870() {
        fn inchi(
            heap: &mut SourceHeap,
            mobile_h: Option<Vec<i8>>,
            fixed_h: Option<Vec<i8>>,
            charge: i32,
            number_of_atoms: i32,
            deleted: i32,
        ) -> SourceMutPointer<INChI> {
            let mobile_h = mobile_h
                .map(|values| heap.allocate_model_storage(values).unwrap())
                .unwrap_or_default();
            let fixed_h = fixed_h
                .map(|values| heap.allocate_model_storage(values).unwrap())
                .unwrap_or_default();
            heap.allocate_model_storage(vec![INChI {
                nNum_H: mobile_h,
                nNum_H_fixed: fixed_h,
                nTotalCharge: charge,
                nNumberOfAtoms: number_of_atoms,
                bDeleted: deleted,
                ..INChI::default()
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let input_fixed = inchi(&mut heap, Some(vec![11, 12, 13]), None, 5, 3, 0);
        let input_mobile = inchi(&mut heap, Some(vec![5, 6, 7]), None, -3, 3, 0);
        let reversed_fixed = inchi(
            &mut heap,
            Some(vec![21, 22, 23]),
            Some(vec![2, 3, 4]),
            6,
            3,
            0,
        );
        let reversed_mobile = inchi(&mut heap, Some(vec![6, 8, 9]), None, -4, 3, 0);
        let mobile_atoms = heap
            .allocate_model_storage(vec![
                inp_ATOM {
                    endpoint: 0,
                    charge: 4,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    endpoint: 0,
                    charge: 5,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    endpoint: 2,
                    charge: 6,
                    ..inp_ATOM::default()
                },
            ])
            .unwrap();
        let norm_data = heap
            .allocate_model_storage(vec![INP_ATOM_DATA {
                at: mobile_atoms,
                ..INP_ATOM_DATA::default()
            }])
            .unwrap();
        let original_groups = heap
            .allocate_model_storage(vec![
                T_GROUP {
                    num: [5, 2, 0, 0, 0],
                    nNumEndpoints: 1,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    num: [7, 2, 0, 0, 0],
                    nNumEndpoints: 2,
                    ..T_GROUP::default()
                },
            ])
            .unwrap();
        let reversed_groups = heap
            .allocate_model_storage(vec![
                T_GROUP {
                    num: [4, 1, 0, 0, 0],
                    nNumEndpoints: 1,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    num: [8, 1, 0, 0, 0],
                    nNumEndpoints: 3,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    num: [6, 2, 0, 0, 0],
                    nNumEndpoints: 4,
                    ..T_GROUP::default()
                },
            ])
            .unwrap();
        let mut structure = StrFromINChI {
            num_atoms: 3,
            endpoint: heap.allocate_model_storage(vec![0_u16, 1, 0]).unwrap(),
            fixed_H: heap.allocate_model_storage(vec![1_i8, 2, 3]).unwrap(),
            nAtno2Canon: [
                heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap(),
                SourceMutPointer::null(),
            ],
            pOneINChI: [reversed_fixed, reversed_mobile],
            pOne_norm_data: [SourceMutPointer::null(), norm_data],
            nNumRemovedProtonsMobHInChI: 4,
            ..StrFromINChI::default()
        };
        structure.ti.num_t_groups = 2;
        structure.ti.t_group = original_groups;
        structure.One_ti.num_t_groups = 3;
        structure.One_ti.t_group = reversed_groups;
        structure.One_ti.tni.nNumRemovedProtons = -2;
        let atoms = vec![
            inp_ATOM {
                num_H: 1,
                charge: -2,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                num_H: 2,
                charge: 3,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                num_H: 3,
                charge: 4,
                ..inp_ATOM::default()
            },
        ];
        let valence = vec![
            VAL_AT {
                cNumValenceElectrons: -1,
                cPeriodicRowNumber: 1,
                ..VAL_AT::default()
            },
            VAL_AT {
                cNumValenceElectrons: 5,
                cPeriodicRowNumber: 2,
                cMetal: 1,
                ..VAL_AT::default()
            },
            VAL_AT {
                cNumValenceElectrons: 6,
                cPeriodicRowNumber: 3,
                ..VAL_AT::default()
            },
        ];
        let mut comparison = CMP2FHINCHI::default();
        comparison.c2at[0].nValue = 99;
        comparison.nNumTgHInChI = 123;
        assert_eq!(
            FillOutCMP2FHINCHI(
                &heap,
                &structure,
                &atoms,
                &valence,
                [input_fixed, input_mobile],
                &mut comparison,
            ),
            Ok(0)
        );
        assert_eq!(comparison.len_c2at, 3);
        assert_eq!(comparison.nNumTgInChI, 2);
        assert_eq!(comparison.nNumTgRevrs, 3);
        assert_eq!(
            (comparison.nNumRemHInChI, comparison.nNumRemHRevrs),
            (4, -2)
        );
        assert_eq!((comparison.nNumTgDiffH, comparison.nNumTgDiffMinus), (0, 1));
        assert_eq!((comparison.nNumTgHInChI, comparison.nNumTgMInChI), (8, 4));
        assert_eq!((comparison.nNumTgHRevrs, comparison.nNumTgMRevrs), (14, 4));
        assert_eq!((comparison.nNumEndpInChI, comparison.nNumEndpRevrs), (1, 1));
        assert_eq!(comparison.bFixedHLayerExistsRevrs, 1);
        assert_eq!(comparison.bHasDifference, 1);
        assert_eq!(comparison.nNumDiffMobH, 1);
        assert_eq!(
            (
                comparison.nChargeFixHInChI,
                comparison.nChargeMobHInChI,
                comparison.nChargeFixHRevrs,
                comparison.nChargeMobHRevrs,
            ),
            (5, -3, 6, -4)
        );
        assert_eq!(
            (
                comparison.nChargeFixHRevrsNonMetal,
                comparison.nChargeMobHRevrsNonMetal,
            ),
            (2, 10)
        );
        assert_eq!(
            (
                comparison.c2at[0].atomNumber,
                comparison.c2at[0].endptInChI,
                comparison.c2at[0].endptRevrs,
                comparison.c2at[0].nValElectr,
                comparison.c2at[0].nPeriodNum,
                comparison.c2at[0].nFixHInChI,
                comparison.c2at[0].nFixHRevrs,
                comparison.c2at[0].nMobHInChI,
                comparison.c2at[0].nMobHRevrs,
                comparison.c2at[0].nNumHRevrs,
                comparison.c2at[0].nAtChargeRevrs,
                comparison.c2at[0].nValue,
            ),
            (0, 0, 0, u8::MAX, 1, 1, 4, 5, 9, 1, -2, 0)
        );
        assert_eq!(
            (
                comparison.c2at[1].atomNumber,
                comparison.c2at[1].endptInChI,
                comparison.c2at[1].endptRevrs,
                comparison.c2at[1].nFixHRevrs,
                comparison.c2at[1].nMobHRevrs,
            ),
            (1, 1, 0, 2, 6)
        );
        assert_eq!(
            (
                comparison.c2at[2].atomNumber,
                comparison.c2at[2].endptInChI,
                comparison.c2at[2].endptRevrs,
                comparison.c2at[2].nFixHRevrs,
                comparison.c2at[2].nMobHRevrs,
            ),
            (2, 0, 2, 3, 8)
        );

        let mut fallback_heap = SourceHeap::default();
        let fallback_mobile_h = fallback_heap.allocate_model_storage(vec![0_i8]).unwrap();
        let fallback_reversed = fallback_heap
            .allocate_model_storage(vec![INChI {
                nNum_H: fallback_mobile_h,
                ..INChI::default()
            }])
            .unwrap();
        let fallback_input = inchi(&mut fallback_heap, Some(vec![0]), None, 0, 1, 0);
        let fallback_mapping = fallback_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let fallback_endpoint = fallback_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let fallback_fixed = fallback_heap.allocate_model_storage(vec![0_i8]).unwrap();
        let fallback_structure = StrFromINChI {
            num_atoms: 1,
            endpoint: fallback_endpoint,
            fixed_H: fallback_fixed,
            nAtno2Canon: [fallback_mapping, SourceMutPointer::null()],
            pOneINChI: [fallback_reversed, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let mut fallback_output = CMP2FHINCHI::default();
        assert_eq!(
            FillOutCMP2FHINCHI(
                &fallback_heap,
                &fallback_structure,
                &[inp_ATOM::default()],
                &[VAL_AT::default()],
                [fallback_input, SourceMutPointer::null()],
                &mut fallback_output,
            ),
            Ok(0)
        );
        assert_eq!(fallback_output.len_c2at, 0);
        assert_eq!(fallback_output.bFixedHLayerExistsRevrs, 0);
        assert_eq!(fallback_output.bHasDifference, 1);

        let mut limit_heap = SourceHeap::default();
        let limit_count = MAX_DIFF_FIXH as usize + 1;
        let limit_reversed = inchi(
            &mut limit_heap,
            None,
            Some(vec![0_i8; limit_count]),
            0,
            limit_count as i32,
            0,
        );
        let limit_input = inchi(&mut limit_heap, None, None, 0, limit_count as i32, 0);
        let limit_mapping = limit_heap
            .allocate_model_storage((0..limit_count as u16).collect())
            .unwrap();
        let limit_endpoint = limit_heap
            .allocate_model_storage(vec![0_u16; limit_count])
            .unwrap();
        let limit_fixed = limit_heap
            .allocate_model_storage(vec![1_i8; limit_count])
            .unwrap();
        let limit_structure = StrFromINChI {
            num_atoms: limit_count as i32,
            endpoint: limit_endpoint,
            fixed_H: limit_fixed,
            nAtno2Canon: [limit_mapping, SourceMutPointer::null()],
            pOneINChI: [limit_reversed, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let limit_atoms = vec![
            inp_ATOM {
                charge: 1,
                ..inp_ATOM::default()
            };
            limit_count
        ];
        let mut limit_output = CMP2FHINCHI::default();
        assert_eq!(
            FillOutCMP2FHINCHI(
                &limit_heap,
                &limit_structure,
                &limit_atoms,
                &vec![VAL_AT::default(); limit_count],
                [limit_input, SourceMutPointer::null()],
                &mut limit_output,
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(limit_output.len_c2at, MAX_DIFF_FIXH as i16);
        assert_eq!(
            limit_output.c2at[MAX_DIFF_FIXH as usize - 1].atomNumber,
            255
        );
        assert_eq!(limit_output.nChargeFixHRevrsNonMetal, 0);
        assert_eq!(limit_output.nChargeFixHInChI, 0);
        assert_eq!(limit_output.nChargeFixHRevrs, 0);
    }

    #[test]
    fn source_port__ichirvr4__filloutcmp2mhinchi__line_1021() {
        fn inchi(
            heap: &mut SourceHeap,
            mobile_h: Option<Vec<i8>>,
            charge: i32,
            number_of_atoms: i32,
            deleted: i32,
        ) -> SourceMutPointer<INChI> {
            let mobile_h = mobile_h
                .map(|values| heap.allocate_model_storage(values).unwrap())
                .unwrap_or_default();
            heap.allocate_model_storage(vec![INChI {
                nNum_H: mobile_h,
                nTotalCharge: charge,
                nNumberOfAtoms: number_of_atoms,
                bDeleted: deleted,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn group_info(
            heap: &mut SourceHeap,
            groups: Vec<T_GROUP>,
            endpoints: Vec<u16>,
        ) -> T_GROUP_INFO {
            let count = groups.len() as i32;
            T_GROUP_INFO {
                t_group: if groups.is_empty() {
                    SourceMutPointer::null()
                } else {
                    heap.allocate_model_storage(groups).unwrap()
                },
                nEndpointAtomNumber: if endpoints.is_empty() {
                    SourceMutPointer::null()
                } else {
                    heap.allocate_model_storage(endpoints).unwrap()
                },
                num_t_groups: count,
                ..T_GROUP_INFO::default()
            }
        }

        for atom_count in [-1, 0] {
            let heap = SourceHeap::default();
            let mut output = CMP2MHINCHI::default();
            output.c2at[0].nValue = 99;
            output.bHasDifference = 1;
            assert_eq!(
                FillOutCMP2MHINCHI(
                    &heap,
                    &StrFromINChI {
                        num_atoms: atom_count,
                        ..StrFromINChI::default()
                    },
                    &ALL_TC_GROUPS::default(),
                    &[],
                    &[],
                    [SourceMutPointer::null(), SourceMutPointer::null()],
                    &mut output,
                ),
                Ok(0),
                "atom count={atom_count}"
            );
            assert_eq!(output, CMP2MHINCHI::default());
        }

        let mut heap = SourceHeap::default();
        let input_base = inchi(&mut heap, Some(vec![0; 10]), 130, 10, 0);
        let input_fixed = inchi(&mut heap, Some((20_i8..30).collect()), 77, 10, 0);
        let atom_to_canonical: Vec<u16> = (0_u16..10).rev().collect();
        let mut reversed_base_h = vec![0_i8; 10];
        reversed_base_h[9] = 7;
        let reversed_base = inchi(&mut heap, Some(reversed_base_h), -130, 10, 0);
        let reversed_fixed = inchi(&mut heap, Some((40_i8..50).collect()), -77, 10, 0);

        let mut atoms = vec![inp_ATOM::default(); 10];
        for atom in &mut atoms {
            atom.valence = 1;
            atom.chem_bonds_valence = 1;
        }
        atoms[0].num_H = 1;
        atoms[1].num_H = 2;
        atoms[1].endpoint = 1;
        atoms[2].charge = -1;
        atoms[3].charge = -1;
        atoms[3].num_H = 1;
        atoms[4].chem_bonds_valence = 2;
        atoms[4].num_H = 1;
        atoms[5].chem_bonds_valence = 2;
        atoms[5].charge = -1;
        atoms[6].chem_bonds_valence = 2;
        atoms[7].num_H = 1;
        atoms[8].charge = -1;
        atoms[9].chem_bonds_valence = 2;

        let mut valence = vec![VAL_AT::default(); 10];
        for item in &mut valence[..7] {
            item.cNumValenceElectrons = 5;
            item.cPeriodicRowNumber = 1;
        }
        for item in &mut valence[7..] {
            item.cNumValenceElectrons = 6;
            item.cPeriodicRowNumber = 2;
        }
        valence[1].cMetal = 1;

        let mut mobile_atoms = atoms.clone();
        for (index, atom) in mobile_atoms.iter_mut().enumerate() {
            atom.endpoint = 0;
            atom.charge = (index + 1) as i8;
        }
        mobile_atoms[2].endpoint = 2;
        let mobile_atoms = heap.allocate_model_storage(mobile_atoms).unwrap();
        let norm_data = heap
            .allocate_model_storage(vec![INP_ATOM_DATA {
                at: mobile_atoms,
                ..INP_ATOM_DATA::default()
            }])
            .unwrap();

        let original_groups = vec![
            T_GROUP {
                num: [5, 2, 0, 0, 0],
                nNumEndpoints: 10,
                ..T_GROUP::default()
            },
            T_GROUP {
                num: [7, 3, 0, 0, 0],
                ..T_GROUP::default()
            },
        ];
        let reversed_groups = vec![
            T_GROUP {
                num: [4, 1, 0, 0, 0],
                nNumEndpoints: 10,
                ..T_GROUP::default()
            },
            T_GROUP {
                num: [7, 3, 0, 0, 0],
                ..T_GROUP::default()
            },
            T_GROUP {
                num: [6, 2, 0, 0, 0],
                ..T_GROUP::default()
            },
        ];
        let original_info = group_info(&mut heap, original_groups, (0_u16..10).collect());
        let reversed_info = group_info(&mut heap, reversed_groups, (0_u16..10).rev().collect());
        let canonical_to_atom = heap
            .allocate_model_storage((0_u16..10).rev().collect())
            .unwrap();
        let atom_to_canonical = heap.allocate_model_storage(atom_to_canonical).unwrap();
        let structure = StrFromINChI {
            num_atoms: 10,
            ti: original_info,
            One_ti: T_GROUP_INFO {
                tni: crate::source_types::TNI {
                    nNumRemovedProtons: -4,
                    ..crate::source_types::TNI::default()
                },
                ..reversed_info
            },
            nNumRemovedProtonsMobHInChI: 6,
            pOneINChI: [reversed_base, reversed_fixed],
            pOne_norm_data: [norm_data, SourceMutPointer::null()],
            nAtno2Canon: [atom_to_canonical, SourceMutPointer::null()],
            nCanon2Atno: [canonical_to_atom, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let tc_groups = ALL_TC_GROUPS {
            tgroup_charge: 3,
            ..ALL_TC_GROUPS::default()
        };
        let mut output = CMP2MHINCHI::default();
        output.c2at[0].nValue = 99;
        output.nNumTgHInChI = 123;
        assert_eq!(
            FillOutCMP2MHINCHI(
                &heap,
                &structure,
                &tc_groups,
                &atoms,
                &valence,
                [input_base, input_fixed],
                &mut output,
            ),
            Ok(0)
        );
        assert_eq!(output.len_c2at, 3);
        assert_eq!((output.nNumTgInChI, output.nNumTgRevrs), (2, 3));
        assert_eq!((output.nNumRemHInChI, output.nNumRemHRevrs), (6, -4));
        assert_eq!((output.nNumTgDiffH, output.nNumTgDiffMinus), (0, 1));
        assert_eq!((output.nNumTgHInChI, output.nNumTgMInChI), (7, 5));
        assert_eq!((output.nNumTgHRevrs, output.nNumTgMRevrs), (11, 6));
        assert_eq!((output.nNumTgNInChI, output.nNumTgOInChI), (7, 3));
        assert_eq!((output.nNumTgNRevrs, output.nNumTgORevrs), (7, 3));
        assert_eq!(
            (
                output.nNumTgNHInChI,
                output.nNumTgNH2InChI,
                output.nNumTgNMinusInChI,
                output.nNumTgNHMinusInChI,
                output.nNumTgDBNHInChI,
                output.nNumTgDBNMinusInChI,
                output.nNumTgDBNInChI,
                output.nNumTgOHInChI,
                output.nNumTgOMinusInChI,
                output.nNumTgDBOInChI,
            ),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        );
        assert_eq!(
            (
                output.nNumTgNHRevrs,
                output.nNumTgNH2Revrs,
                output.nNumTgNMinusRevrs,
                output.nNumTgNHMinusRevrs,
                output.nNumTgDBNHRevrs,
                output.nNumTgDBNMinusRevrs,
                output.nNumTgDBNRevrs,
                output.nNumTgOHRevrs,
                output.nNumTgOMinusRevrs,
                output.nNumTgDBORevrs,
            ),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        );
        assert_eq!((output.nNumEndpInChI, output.nNumEndpRevrs), (1, 1));
        assert_eq!(output.bFixedHLayerExistsRevrs, 1);
        assert_eq!(output.bHasDifference, 1);
        assert_eq!(output.nNumDiffMobH, 1);
        assert_eq!(
            (
                output.nChargeMobHInChI,
                output.nChargeMobHRevrs,
                output.nChargeMobHRevrsNonMetal,
            ),
            (-126, 126, 53)
        );
        assert_eq!(
            (
                output.c2at[0].atomNumber,
                output.c2at[0].endptInChI,
                output.c2at[0].endptRevrs,
                output.c2at[0].nValElectr,
                output.c2at[0].nPeriodNum,
                output.c2at[0].nMobHInChI,
                output.c2at[0].nMobHRevrs,
                output.c2at[0].nNumHRevrs,
                output.c2at[0].nAtChargeRevrs,
                output.c2at[0].nValue,
            ),
            (0, 0, 0, 5, 1, 20, 49, 1, 0, 0)
        );
        assert_eq!(
            (
                output.c2at[1].atomNumber,
                output.c2at[1].endptInChI,
                output.c2at[1].endptRevrs,
                output.c2at[1].nMobHInChI,
                output.c2at[1].nMobHRevrs,
            ),
            (1, 1, 0, 21, 48)
        );
        assert_eq!(
            (
                output.c2at[2].atomNumber,
                output.c2at[2].endptInChI,
                output.c2at[2].endptRevrs,
                output.c2at[2].nMobHInChI,
                output.c2at[2].nMobHRevrs,
            ),
            (2, 0, 2, 22, 47)
        );

        let mut syntax_heap = SourceHeap::default();
        let syntax_info = group_info(
            &mut syntax_heap,
            vec![T_GROUP {
                nNumEndpoints: 1,
                ..T_GROUP::default()
            }],
            vec![0],
        );
        let mut syntax_output = CMP2MHINCHI::default();
        syntax_output.c2at[0].nValue = 99;
        assert_eq!(
            FillOutCMP2MHINCHI(
                &syntax_heap,
                &StrFromINChI {
                    ti: syntax_info,
                    ..StrFromINChI::default()
                },
                &ALL_TC_GROUPS::default(),
                &[inp_ATOM::default()],
                &[VAL_AT::default()],
                [SourceMutPointer::null(), SourceMutPointer::null()],
                &mut syntax_output,
            ),
            Ok(RI_ERR_SYNTAX)
        );
        assert_eq!((syntax_output.nNumTgInChI, syntax_output.len_c2at), (1, 0));
        assert_eq!(syntax_output.c2at[0].nValue, 0);

        let mut program_heap = SourceHeap::default();
        let program_info = group_info(
            &mut program_heap,
            vec![T_GROUP {
                nNumEndpoints: 1,
                ..T_GROUP::default()
            }],
            vec![0],
        );
        let program_map = program_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut program_output = CMP2MHINCHI::default();
        assert_eq!(
            FillOutCMP2MHINCHI(
                &program_heap,
                &StrFromINChI {
                    One_ti: program_info,
                    nCanon2Atno: [program_map, SourceMutPointer::null()],
                    ..StrFromINChI::default()
                },
                &ALL_TC_GROUPS::default(),
                &[inp_ATOM::default()],
                &[VAL_AT::default()],
                [SourceMutPointer::null(), SourceMutPointer::null()],
                &mut program_output,
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(
            (program_output.nNumTgRevrs, program_output.len_c2at),
            (1, 0)
        );

        let mut limit_heap = SourceHeap::default();
        let limit_count = MAX_DIFF_FIXH as usize + 1;
        let limit_input = inchi(
            &mut limit_heap,
            Some(vec![0; limit_count]),
            8,
            limit_count as i32,
            0,
        );
        let limit_reversed = inchi(
            &mut limit_heap,
            Some(vec![1; limit_count]),
            7,
            limit_count as i32,
            0,
        );
        let limit_mapping = limit_heap
            .allocate_model_storage((0..limit_count as u16).collect())
            .unwrap();
        let limit_structure = StrFromINChI {
            num_atoms: limit_count as i32,
            pOneINChI: [limit_reversed, SourceMutPointer::null()],
            nAtno2Canon: [limit_mapping, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let mut limit_output = CMP2MHINCHI::default();
        assert_eq!(
            FillOutCMP2MHINCHI(
                &limit_heap,
                &limit_structure,
                &ALL_TC_GROUPS::default(),
                &vec![inp_ATOM::default(); limit_count],
                &vec![VAL_AT::default(); limit_count],
                [limit_input, SourceMutPointer::null()],
                &mut limit_output,
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(limit_output.len_c2at, MAX_DIFF_FIXH as i16);
        assert_eq!(limit_output.c2at[255].atomNumber, 255);
        assert_eq!(limit_output.nNumDiffMobH, 0);
        assert_eq!(limit_output.nChargeMobHInChI, 0);
        assert_eq!(limit_output.nChargeMobHRevrs, 0);
    }
}
