use crate::source::base::ichi_bns::{
    BalancedNetworkSearch, BnsSearchWorkspace, DeAllocateBnStruct, ReInitBnData, ReInitBnStructAltPaths,
    RunBalancedNetworkSearch, run_balanced_network_search_with_workspace,
};
use crate::source::base::ichi_io::{inchi_strbuf_close, inchi_strbuf_init};
use crate::source::base::ichimake::Create_INChI;
use crate::source::base::ichinorm::FreeInpAtomData;
use crate::source::base::ichiread::bRevInchiComponentExists;
use crate::source::base::ichiring::{QueueCreate, QueueDelete};
use crate::source::base::ichirvr2::{CopyBnsToAtom, CopySt2At};
use crate::source::base::ichister::{FixUnkn0DStereoBonds, ReconcileAllCmlBondParities};
use crate::source::base::ichitaut::free_t_group_info;
use crate::source::base::mol2atom::{CreateInpAtomData, FreeOrigAtData};
use crate::source::base::runichi::ProcessOneStructure;
use crate::source::base::strutil::{
    Alloc_INChI, Alloc_INChI_Aux, Free_INChI, Free_INChI_Aux, SetConnectedComponentNumber, bNumHeterAtomHasIsotopicH,
    fix_odd_things,
};
use crate::source::base::util::{
    get_endpoint_valence, inchi_calloc, inchi_free, is_el_a_metal, n_no_metal_bonds_valence, n_no_metal_num_bonds,
};
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, _IS_OKAY, _IS_UNKNOWN, _IS_WARNING, ALL_TC_GROUPS, AT_NUMB, AT_RANK, BFS_Q, BFS_Q_CLEAR,
    BFS_Q_FREE, BN_DATA, BN_MAX_ALTP, BN_STRUCT, BNS_ADD_EDGES, BNS_ADD_SUPER_TGROUP, BNS_ALT_PATH, BNS_BOND_ERR,
    BNS_CPOINT_ERR, BNS_EDGE, BNS_EDGE_FORBIDDEN_MASK, BNS_EF_CHNG_FLOW, BNS_ERR, BNS_IEDGE, BNS_MAX_ERR_VALUE,
    BNS_PROGRAM_ERR, BNS_VERT_EDGE_OVFL, BNS_VERT_TYPE__AUX, BNS_VERT_TYPE_ATOM, BNS_VERT_TYPE_C_GROUP,
    BNS_VERT_TYPE_C_NEGATIVE, BNS_VERT_TYPE_ENDPOINT, BNS_VERT_TYPE_METAL_GR, BNS_VERT_TYPE_SUPER_TGROUP,
    BNS_VERT_TYPE_TEMP, BNS_VERT_TYPE_TGROUP, BNS_VERTEX, BNS_VT_C_NEG, BNS_VT_C_NEG_ALL, BNS_VT_C_NEG_C,
    BNS_VT_C_NEG_M, BNS_VT_C_POS, BNS_VT_C_POS_ALL, BNS_VT_C_POS_C, BNS_VT_C_POS_M, BNS_VT_M_GROUP, BNS_VT_YVCONNECTOR,
    BOND_TYPE_DOUBLE, BOND_TYPE_MASK, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, CANON_GLOBALS, CC_CAND, EDGE_FLOW_MASK,
    EDGE_FLOW_ST_MASK, EDGE_LIST, EDGE_LIST_CLEAR, EDGE_LIST_FREE, EL_TYPE_C, EL_TYPE_N, EL_TYPE_O, EL_TYPE_OSt,
    EL_TYPE_P, EL_TYPE_PT, EL_TYPE_S, EL_TYPE_X, INC_NUM_TCGROUPS, INCHI_BAS, INCHI_CLOCK, INCHI_IOS_STRING, INCHI_NUM,
    INCHI_OUT_NO_AUX_INFO, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT, INChI, INP_ATOM_DATA, INPUT_PARMS,
    MAX_BOND_EDGE_CAP, MAX_NUM_STEREO_BONDS, MAX_TGROUP_EDGE_CAP, MAXVAL, NO_VERTEX, NUM_H_ISOTOPES,
    NUM_KINDS_OF_GROUPS, ORIG_ATOM_DATA, REQ_MODE_BASIC, REQ_MODE_ISO, REQ_MODE_TAUT, RI_ERR_ALLOC, RI_ERR_PROGR,
    RI_ERR_SYNTAX, SRM, STRUCT_DATA, SourceHeap, SourceHeapError, SourceMutPointer, SourceTGroupInfoPointer,
    StrFromINChI, T_GROUP_HDR_LEN, T_GROUP_INFO, TAUT_NON, TAUT_NUM, TAUT_YES, TC_GROUP, TG_FLAG_FIX_ODD_THINGS_DONE,
    TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE, TG_FLAG_H_ALREADY_REMOVED, TGSO_SYMM_IORDER,
    TGSO_TOTAL_LEN, VAL_AT, clock_t, inp_ATOM,
    local_ichirvr1::{NUM_VF, VF, VF_USED_ALL, VF_USED_IN, VF_USED_OUT},
    tagAltPathConst_iALTP_END_ATOM, tagAltPathConst_iALTP_FLOW, tagAltPathConst_iALTP_HDR_LEN,
    tagAltPathConst_iALTP_MAX_LEN, tagAltPathConst_iALTP_PATH_LEN, tagAltPathConst_iALTP_START_ATOM,
    tagTCGroupTypes_TCG_MeFlower0 as TCG_MeFlower0, tagTCGroupTypes_TCG_MeFlower1 as TCG_MeFlower1,
    tagTCGroupTypes_TCG_MeFlower2 as TCG_MeFlower2, tagTCGroupTypes_TCG_MeFlower3 as TCG_MeFlower3,
    tagTCGroupTypes_TCG_Minus as TCG_Minus, tagTCGroupTypes_TCG_Minus_C0 as TCG_Minus_C0,
    tagTCGroupTypes_TCG_Minus_M0 as TCG_Minus_M0, tagTCGroupTypes_TCG_Minus0 as TCG_Minus0,
    tagTCGroupTypes_TCG_Plus as TCG_Plus, tagTCGroupTypes_TCG_Plus_C0 as TCG_Plus_C0,
    tagTCGroupTypes_TCG_Plus_M0 as TCG_Plus_M0, tagTCGroupTypes_TCG_Plus0 as TCG_Plus0,
};

fn copy_inp_atom_prefix(
    heap: &mut SourceHeap,
    destination: SourceMutPointer<inp_ATOM>,
    source: SourceMutPointer<inp_ATOM>,
    count: usize,
) -> Result<(), SourceHeapError> {
    if destination == source {
        heap.slice(source.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        return Ok(());
    }
    if destination.allocation_identity() != source.allocation_identity() {
        // SAFETY: distinct allocation identities prove that source and
        // destination do not alias. Neither allocation is resized or freed
        // during this source-equivalent memcpy.
        let source_view = unsafe { heap.stable_slice(source.as_const())? };
        let source = source_view.prefix(count)?;
        heap.slice_mut(destination)?
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(source);
        return Ok(());
    }

    // Preserve the modeled overlapping-pointer behavior even though the
    // Official production call sites use distinct allocations.
    let source = heap
        .slice(source.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(destination)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone_from_slice(&source);
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn clear_t_group_info(
    heap: &mut SourceHeap,
    info: Option<&mut T_GROUP_INFO>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:509 clear_t_group_info
    // INCHI✔️✔️: void clear_t_group_info( T_GROUP_INFO *ti )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (!ti)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         T_GROUP   *t_group = ti->t_group;
    // INCHI✔️✔️:         int       max_num_t_groups = ti->max_num_t_groups;
    // INCHI✔️✔️:         AT_NUMB   *tGroupNumber = ti->tGroupNumber;
    // INCHI✔️✔️:         int       num_t_groups = ti->num_t_groups;
    // INCHI✔️✔️:         AT_NUMB   *nEndpointAtomNumber = ti->nEndpointAtomNumber;
    // INCHI✔️✔️:         int       nNumEndpoints = ti->nNumEndpoints;
    // INCHI✔️✔️:         AT_NUMB   *nIsotopicEndpointAtomNumber = ti->nIsotopicEndpointAtomNumber;
    // INCHI✔️✔️:         int       nNumIsotopicEndpoints = ti->nNumIsotopicEndpoints;
    // INCHI✔️✔️:         memset( ti, 0, sizeof( *ti ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         if (t_group)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             memset( t_group, 0, sizeof( t_group[0] )*max_num_t_groups ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             max_num_t_groups = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (tGroupNumber)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             memset( tGroupNumber, 0, sizeof( tGroupNumber[0] )*num_t_groups ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             num_t_groups = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nEndpointAtomNumber)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             memset( nEndpointAtomNumber, 0, sizeof( nEndpointAtomNumber[0] )*nNumEndpoints ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nNumEndpoints = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nIsotopicEndpointAtomNumber)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             memset( nIsotopicEndpointAtomNumber, 0, sizeof( nIsotopicEndpointAtomNumber[0] )*nNumIsotopicEndpoints ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nNumIsotopicEndpoints = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         ti->t_group = t_group;
    // INCHI✔️✔️:         ti->max_num_t_groups = max_num_t_groups;
    // INCHI✔️✔️:         ti->tGroupNumber = tGroupNumber;
    // INCHI✔️✔️:         ti->num_t_groups = num_t_groups;
    // INCHI✔️✔️:         ti->nEndpointAtomNumber = nEndpointAtomNumber;
    // INCHI✔️✔️:         ti->nNumEndpoints = nNumEndpoints;
    // INCHI✔️✔️:         ti->nIsotopicEndpointAtomNumber = nIsotopicEndpointAtomNumber;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         ti->nNumIsotopicEndpoints = nNumIsotopicEndpoints;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: clear_t_group_info

    let Some(info) = info else {
        return Ok(());
    };
    let t_group = info.t_group;
    let mut max_num_t_groups = info.max_num_t_groups;
    let t_group_number = info.tGroupNumber;
    let mut num_t_groups = info.num_t_groups;
    let endpoint_atom_number = info.nEndpointAtomNumber;
    let mut num_endpoints = info.nNumEndpoints;
    let isotopic_endpoint_atom_number = info.nIsotopicEndpointAtomNumber;
    let mut num_isotopic_endpoints = info.nNumIsotopicEndpoints;
    *info = T_GROUP_INFO::default();

    if t_group.is_null() {
        max_num_t_groups = 0;
    } else {
        let count = usize::try_from(max_num_t_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice_mut(t_group)?
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(Default::default());
    }
    if t_group_number.is_null() {
        num_t_groups = 0;
    } else {
        let count = usize::try_from(num_t_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice_mut(t_group_number)?
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(0);
    }
    if endpoint_atom_number.is_null() {
        num_endpoints = 0;
    } else {
        let count = usize::try_from(num_endpoints).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice_mut(endpoint_atom_number)?
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(0);
    }
    if isotopic_endpoint_atom_number.is_null() {
        num_isotopic_endpoints = 0;
    } else {
        let count = usize::try_from(num_isotopic_endpoints).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice_mut(isotopic_endpoint_atom_number)?
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(0);
    }

    info.t_group = t_group;
    info.max_num_t_groups = max_num_t_groups;
    info.tGroupNumber = t_group_number;
    info.num_t_groups = num_t_groups;
    info.nEndpointAtomNumber = endpoint_atom_number;
    info.nNumEndpoints = num_endpoints;
    info.nIsotopicEndpointAtomNumber = isotopic_endpoint_atom_number;
    info.nNumIsotopicEndpoints = num_isotopic_endpoints;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn FillOutpStructEndpointFromInChI(
    heap: &mut SourceHeap,
    inchi: &INChI,
    endpoint_output: &mut SourceMutPointer<AT_NUMB>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:692 FillOutpStructEndpointFromInChI
    // INCHI✔️❌: int FillOutpStructEndpointFromInChI( INChI *pInChI, AT_NUMB **pEndpoint )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int num_at = pInChI->nNumberOfAtoms;
    // INCHI✔️❌:     AT_NUMB *endpoint = *pEndpoint;
    // INCHI✔️❌:     int     itg, i, j, k, len_tg;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!endpoint && !( endpoint = (AT_NUMB*) inchi_malloc( num_at * sizeof( endpoint[0] ) ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return RI_ERR_ALLOC;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( endpoint, 0, num_at * sizeof( endpoint[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     if (pInChI->lenTautomer <= 1 || !pInChI->nTautomer)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     j = 1; /* index in pInChI->nTautomer[] */
    // INCHI✔️❌:     i = 0; /* index in ti->nEndpointAtomNumber[] */
    // INCHI✔️❌:     for (itg = 0; itg < pInChI->nTautomer[0]; itg++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len_tg = pInChI->nTautomer[j]; /* t-group length not including pInChI->nTautomer[j] */
    // INCHI✔️❌:         j += T_GROUP_HDR_LEN;   /* skip t-group header */
    // INCHI✔️❌:         len_tg -= T_GROUP_HDR_LEN - 1;
    // INCHI✔️❌:         /* ti->t_group[itg].nNumEndpoints = len_tg; */
    // INCHI✔️❌:         for (; 0 < len_tg--; j++, i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             k = pInChI->nTautomer[j] - 1;
    // INCHI✔️❌:             endpoint[k] = itg + 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     *pEndpoint = endpoint;
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FillOutpStructEndpointFromInChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillOutpStructEndpointFromInChI
    // INCHI✔️❌: #define inchi_malloc malloc
    // INCHI✔️❌: typedef unsigned short AT_NUMB
    // INCHI✔️❌: #define T_GROUP_HDR_LEN (1 + INCHI_T_NUM_MOVABLE)
    // END INCHI ACTIVE MACRO CONFIGURATION: FillOutpStructEndpointFromInChI

    let num_atoms = inchi.nNumberOfAtoms;
    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut endpoint = *endpoint_output;
    if endpoint.is_null() {
        endpoint = match inchi_calloc::<AT_NUMB>(heap, num_atoms as i64 as u64, std::mem::size_of::<AT_NUMB>() as u64) {
            Ok(pointer) => pointer,
            Err(
                SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange,
            ) => return Ok(RI_ERR_ALLOC),
            Err(error) => return Err(error),
        };
    }
    heap.slice_mut(endpoint)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(0);

    if inchi.lenTautomer > 1 && !inchi.nTautomer.is_null() {
        let tautomer = heap.slice(inchi.nTautomer.as_const())?.to_vec();
        let number_of_groups = i32::from(*tautomer.first().ok_or(SourceHeapError::PointerOutOfBounds)?);
        let mut group_index = 0_i32;
        let mut source_index = 1_i32;
        while group_index < number_of_groups {
            let mut group_length = i32::from(
                *tautomer
                    .get(usize::try_from(source_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            source_index = source_index.wrapping_add(T_GROUP_HDR_LEN as i32);
            group_length = group_length.wrapping_sub(T_GROUP_HDR_LEN as i32 - 1);
            while group_length > 0 {
                let atom_number = i32::from(
                    *tautomer
                        .get(usize::try_from(source_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
                .wrapping_sub(1);
                heap.slice_mut(endpoint)?
                    .get_mut(usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone_from(&(group_index.wrapping_add(1) as AT_NUMB));
                group_length = group_length.wrapping_sub(1);
                source_index = source_index.wrapping_add(1);
            }
            group_index = group_index.wrapping_add(1);
        }
    }
    *endpoint_output = endpoint;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn cmp_charge_val(first: &crate::source_types::CHARGE_VAL, second: &crate::source_types::CHARGE_VAL) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:732 cmp_charge_val
    // INCHI✔️✔️: int cmp_charge_val( const void *a1, const void *a2, void *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     const CHARGE_VAL *p1 = (const CHARGE_VAL *) a1;
    // INCHI✔️✔️:     const CHARGE_VAL *p2 = (const CHARGE_VAL *) a2;
    // INCHI✔️✔️:     int    diff;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if ((diff = (int) p1->nValence - (int) p2->nValence))  /* smaller valence first */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return diff;
    // INCHI✔️✔️:     if ((diff = abs( (int) p1->nCharge ) - abs( (int) p2->nCharge ))) /* smaller abs charge first */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return diff;
    // INCHI✔️✔️:     if ((diff = (int) p2->nCharge - (int) p1->nCharge)) /* (+) first, (-) second */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return diff;
    // INCHI✔️✔️:     return (int) p1->nValenceOrderingNumber - (int) p2->nValenceOrderingNumber;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: cmp_charge_val
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: cmp_charge_val
    // INCHI✔️✔️: typedef struct tagChargeValence {
    // INCHI✔️✔️:     int nValence;
    // INCHI✔️✔️:     int nCharge;
    // INCHI✔️✔️:     int nValenceOrderingNumber;
    // INCHI✔️✔️: } CHARGE_VAL;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: cmp_charge_val

    let mut difference = first.nValence.wrapping_sub(second.nValence);
    if difference != 0 {
        return difference;
    }
    difference = first.nCharge.wrapping_abs().wrapping_sub(second.nCharge.wrapping_abs());
    if difference != 0 {
        return difference;
    }
    difference = second.nCharge.wrapping_sub(first.nCharge);
    if difference != 0 {
        return difference;
    }
    first.nValenceOrderingNumber.wrapping_sub(second.nValenceOrderingNumber)
}

#[allow(non_snake_case)]
pub(crate) fn bMayBeACationInMobileHLayer(
    atoms: &[inp_ATOM],
    valence_atoms: &[crate::source_types::VAL_AT],
    atom_index: i32,
    mobile_h: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:749 bMayBeACationInMobileHLayer
    // INCHI✔️✔️: int bMayBeACationInMobileHLayer( inp_ATOM *at, VAL_AT *pVA, int iat, int bMobileH )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int    j, neigh;
    // INCHI✔️✔️:     U_CHAR cVal;  /* moved from below 2024-09-01 DT */
    // INCHI✔️✔️:     if (!bMobileH || !at[iat].num_H)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* cVal, cation valence */
    // INCHI✔️✔️:     switch ( at[iat].el_number ) {
    // INCHI✔️✔️:         case EL_NUMBER_N: /* fallthrough */
    // INCHI✔️✔️:         case EL_NUMBER_P:
    // INCHI✔️✔️:             cVal = 4;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case EL_NUMBER_O: /* fallthrough */
    // INCHI✔️✔️:         case EL_NUMBER_S:
    // INCHI✔️✔️:         case EL_NUMBER_SE:
    // INCHI✔️✔️:         case EL_NUMBER_TE:
    // INCHI✔️✔️:             cVal = 3;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at[iat].valence + at[iat].num_H <= cVal)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (j = 0; j < at[iat].valence; j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             neigh = at[iat].neighbor[j];
    // INCHI✔️✔️:             if (at[neigh].valence == 4 && at[neigh].chem_bonds_valence == 4 && !at[neigh].num_H &&
    // INCHI✔️✔️:                     pVA[neigh].cNumValenceElectrons == 3 && pVA[neigh].cPeriodicRowNumber == 1)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: bMayBeACationInMobileHLayer

    let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    if mobile_h == 0 || atom.num_H == 0 {
        return Ok(1);
    }
    let cation_valence = match atom.el_number {
        7 | 15 => 4_i32,
        8 | 16 | 34 | 52 => 3_i32,
        _ => return Ok(1),
    };
    if i32::from(atom.valence) + i32::from(atom.num_H) <= cation_valence {
        let mut ordinal = 0_i32;
        while ordinal < i32::from(atom.valence) {
            let neighbor = usize::from(
                *atom
                    .neighbor
                    .get(usize::try_from(ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let neighbor_atom = atoms.get(neighbor).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor_valence = valence_atoms.get(neighbor).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor_atom.valence == 4
                && neighbor_atom.chem_bonds_valence == 4
                && neighbor_atom.num_H == 0
                && neighbor_valence.cNumValenceElectrons == 3
                && neighbor_valence.cPeriodicRowNumber == 1
            {
                return Ok(1);
            }
            ordinal = ordinal.wrapping_add(1);
        }
        return Ok(0);
    }
    Ok(1)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn clean_charge_val(
    canonical_globals: &mut CANON_GLOBALS,
    charge_values: &mut [crate::source_types::CHARGE_VAL],
    len: i32,
    atoms: &[inp_ATOM],
    valence_atoms: &[crate::source_types::VAL_AT],
    atom_index: i32,
    is_metal: i32,
    mobile_h: i32,
    endpoint: Option<&[AT_NUMB]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:792 clean_charge_val
    // INCHI✔️✔️: complete active source frame follows verbatim.
    // INCHI✔️✔️: /************************************************************************************/
    // INCHI✔️✔️: int clean_charge_val( struct tagCANON_GLOBALS *pCG, CHARGE_VAL *pChargeVal, int len, inp_ATOM *atom, VAL_AT *pVA,
    // INCHI✔️✔️:                       int iat, int bIsMetal, int bMobileH, AT_NUMB *endpoint )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     inp_ATOM *at = atom + iat;
    // INCHI✔️✔️:     int nPeriodicNum = at->el_number;
    // INCHI✔️✔️:     int num_bonds = at->valence;
    // INCHI✔️✔️:     int min_valence = at->valence + at->num_H;
    // INCHI✔️✔️:     /* in fixed-H case treat tautomeric -O as tautomeric to avoid #O(+) */
    // INCHI✔️✔️:     int bTautomeric = ( at->endpoint != 0 );
    // INCHI✔️✔️:     int bFixedHTautomeric = !bMobileH && ( endpoint && endpoint[iat] &&
    // INCHI✔️✔️:                             pVA[iat].cNumValenceElectrons == 6 && 1 == num_bonds &&
    // INCHI✔️✔️:                             !at->num_H && !bIsMetal );
    // INCHI✔️✔️:     /* int bIsMetal      = is_el_a_metal( nPeriodicNum );*/
    // INCHI✔️✔️:     int bDoNotAddH = if_skip_add_H( nPeriodicNum );
    // INCHI✔️✔️:     int nPeriod, nNumEqAbsCharges;
    // INCHI✔️✔️:     int nNumValenceEl = get_sp_element_type( nPeriodicNum, &nPeriod ) - 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int i, j;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!len)
    // INCHI✔️✔️:         return len;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     insertions_sort( pCG, pChargeVal, len, sizeof( pChargeVal[0] ), cmp_charge_val );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* metals -- very preliminary code */
    // INCHI✔️✔️:     if (bIsMetal && bDoNotAddH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️: /* keep the 1st found */
    // INCHI✔️✔️:         return inchi_min( 1, len );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     /* Mobile-H layer cannot have H on positively charged N, P (all IV), O, S, Se, Te (all III) */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     if ( abs( pChargeVal[0].nCharge ) > 1 && pChargeVal[0].nValence >= min_valence ) {
    // INCHI✔️✔️:         return inchi_min( 1, len );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:     nNumEqAbsCharges = 0;
    // INCHI✔️✔️:     for (i = j = 0; i < len && j < ( nNumEqAbsCharges ? 3 + nNumEqAbsCharges : 4 ); i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* for now accept only charge = 0, -1, +1 */
    // INCHI✔️✔️:         if (abs( pChargeVal[i].nCharge ) > 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * ( min_valence - 1 ) < pChargeVal[i].nValence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue; /* not more than one triple and the rest - double bonds per atom */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (( bTautomeric || (j && bFixedHTautomeric) ) && pChargeVal[i].nCharge < 0) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue; /* negative charge must be included in the tautomeric group */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (( bTautomeric || bFixedHTautomeric ) && pChargeVal[i].nCharge > 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue; /* positive charge for now cannot reach a tautomeric group */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (j && !bMayBeACationInMobileHLayer( atom, pVA, iat, bMobileH ) && pChargeVal[i].nCharge > 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (i + 1 < len &&
    // INCHI✔️✔️:                  pChargeVal[i].nValence == pChargeVal[i + 1].nValence &&
    // INCHI✔️✔️:                  pChargeVal[i].nCharge == -pChargeVal[i + 1].nCharge)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* (-) if exists is always after (+) */
    // INCHI✔️✔️:                 i += 1; /* also skip the next element */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             continue; /* in case of Mobile-H, a hydrogen cannot be on a (+)-charged heteroatom */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* accept same valence opposite charges only for C and its group in Periodic Table */
    // INCHI✔️✔️:         if (j && !bTautomeric &&
    // INCHI✔️✔️:              pChargeVal[i].nValence == pChargeVal[j - 1].nValence &&
    // INCHI✔️✔️:              pChargeVal[i].nCharge == -pChargeVal[j - 1].nCharge)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (nNumValenceEl == VALUE_OCTET / 2 && pChargeVal[i].nCharge && !nNumEqAbsCharges)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 pChargeVal[j++] = pChargeVal[i];
    // INCHI✔️✔️:                 nNumEqAbsCharges++;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* do not accept valence=5 for neutral NHn in case of not Mobile-H 2005-01-26 ???? */
    // INCHI✔️✔️:         if (nNumValenceEl == 5 && nPeriod == 1 && at->num_H &&
    // INCHI✔️✔️:              j && !bMobileH &&
    // INCHI✔️✔️:              pChargeVal[i].nValence == 5 && !pChargeVal[i].nCharge)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* do not accept gaps in allowed valences */
    // INCHI✔️✔️:         if (j && pChargeVal[i].nValence > pChargeVal[j - 1].nValence + 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         pChargeVal[j++] = pChargeVal[i];
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     len = j;
    // INCHI✔️✔️:     if (!nNumEqAbsCharges && num_bonds < 3 && len == 4)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         len--; /* prohibit =S#  where # is a triple bond */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return len;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: clean_charge_val

    let _ = canonical_globals;
    if len == 0 {
        return Ok(len);
    }
    let count = usize::try_from(len).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if count > charge_values.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut sorted_index = 1_usize;
    while sorted_index < count {
        let inserted = charge_values[sorted_index].clone();
        let mut destination = sorted_index;
        while destination > 0 && cmp_charge_val(&inserted, &charge_values[destination - 1]) < 0 {
            charge_values[destination] = charge_values[destination - 1].clone();
            destination -= 1;
        }
        charge_values[destination] = inserted;
        sorted_index += 1;
    }

    let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let num_bonds = i32::from(atom.valence);
    let min_valence = num_bonds + i32::from(atom.num_H);
    let tautomeric = atom.endpoint != 0;
    let fixed_h_tautomeric = mobile_h == 0
        && endpoint.and_then(|values| values.get(index)).copied().unwrap_or(0) != 0
        && valence_atoms
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .cNumValenceElectrons
            == 6
        && num_bonds == 1
        && atom.num_H == 0
        && is_metal == 0;
    let do_not_add_h = crate::source::base::util::if_skip_add_H(i32::from(atom.el_number))?;
    let mut period = 0_i32;
    let num_valence_electrons = get_sp_element_type(i32::from(atom.el_number), &mut period).wrapping_sub(1);

    if is_metal != 0 && do_not_add_h != 0 {
        return Ok(len.min(1));
    }

    let mut equal_absolute_charges = 0_i32;
    let mut source_index = 0_i32;
    let mut accepted = 0_i32;
    while source_index < len
        && accepted
            < if equal_absolute_charges != 0 {
                3_i32.wrapping_add(equal_absolute_charges)
            } else {
                4
            }
    {
        let source = usize::try_from(source_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let candidate = charge_values[source].clone();
        if candidate.nCharge.wrapping_abs() > 1 {
            source_index = source_index.wrapping_add(1);
            continue;
        }
        if (BOND_TYPE_TRIPLE as i32).wrapping_add((BOND_TYPE_DOUBLE as i32).wrapping_mul(min_valence.wrapping_sub(1)))
            < candidate.nValence
        {
            source_index = source_index.wrapping_add(1);
            continue;
        }
        if (tautomeric || (accepted != 0 && fixed_h_tautomeric)) && candidate.nCharge < 0 {
            source_index = source_index.wrapping_add(1);
            continue;
        }
        if (tautomeric || fixed_h_tautomeric) && candidate.nCharge > 0 {
            source_index = source_index.wrapping_add(1);
            continue;
        }
        if accepted != 0
            && bMayBeACationInMobileHLayer(atoms, valence_atoms, atom_index, mobile_h)? == 0
            && candidate.nCharge > 0
        {
            if source_index.wrapping_add(1) < len {
                let next = &charge_values
                    [usize::try_from(source_index.wrapping_add(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                if candidate.nValence == next.nValence && candidate.nCharge == next.nCharge.wrapping_neg() {
                    source_index = source_index.wrapping_add(1);
                }
            }
            source_index = source_index.wrapping_add(1);
            continue;
        }
        if accepted != 0
            && !tautomeric
            && candidate.nValence
                == charge_values
                    [usize::try_from(accepted.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .nValence
            && candidate.nCharge
                == charge_values
                    [usize::try_from(accepted.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .nCharge
                .wrapping_neg()
        {
            if num_valence_electrons == crate::source_types::VALUE_OCTET as i32 / 2
                && candidate.nCharge != 0
                && equal_absolute_charges == 0
            {
                charge_values[usize::try_from(accepted).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = candidate;
                accepted = accepted.wrapping_add(1);
                equal_absolute_charges = equal_absolute_charges.wrapping_add(1);
            }
            source_index = source_index.wrapping_add(1);
            continue;
        }
        if num_valence_electrons == 5
            && period == 1
            && atom.num_H != 0
            && accepted != 0
            && mobile_h == 0
            && candidate.nValence == 5
            && candidate.nCharge == 0
        {
            source_index = source_index.wrapping_add(1);
            continue;
        }
        if accepted != 0
            && candidate.nValence
                > charge_values
                    [usize::try_from(accepted.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .nValence
                .wrapping_add(1)
        {
            break;
        }
        charge_values[usize::try_from(accepted).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = candidate;
        accepted = accepted.wrapping_add(1);
        source_index = source_index.wrapping_add(1);
    }
    let mut result = accepted;
    if equal_absolute_charges == 0 && num_bonds < 3 && result == 4 {
        result = result.wrapping_sub(1);
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetAtomRestoreInfo(
    canonical_globals: &mut CANON_GLOBALS,
    atoms: &mut [inp_ATOM],
    atom_index: i32,
    valence_atoms: &mut [crate::source_types::VAL_AT],
    restore_mode: &SRM,
    mobile_h: i32,
    endpoint: Option<&[AT_NUMB]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:912 GetAtomRestoreInfo
    // INCHI✔️✔️: complete active source frame follows verbatim.
    /*
    int GetAtomRestoreInfo(struct tagCANON_GLOBALS* pCG,
        inp_ATOM* atom, int iat, VAL_AT* pVArray,
        ICHICONST SRM* pSrm, int bMobileH, AT_NUMB* endpoint)
    {
        /* #defines from util.c */
    #define MIN_ATOM_CHARGE        (-2)
    #define MAX_ATOM_CHARGE         2
    #define NEUTRAL_STATE          (-MIN_ATOM_CHARGE)
    #define NUM_ATOM_CHARGES       (MAX_ATOM_CHARGE - MIN_ATOM_CHARGE + 1)
    #define MAX_NUM_VALENCES        5  /* max. number + 1 to provide zero termination */

        int i, j, j2, k, k2, charge, cur_charge, num_non_bonding_electrons;
        int nNumStates, nNumSelectedStates, num_H, num_bonds;
        int nOctetNeutralValenceExcess, nFirstNeutralValenceExcess;
        int nFoundNeutralValenceExcess, nFoundNeutralValenceOrdNumber;
        int nLastFoundValenceOrdNumber, nLastFoundValenceState;
        int cn_bits, cn_bits_array[5], len_cn_bits_array;
        inp_ATOM* at = atom + iat;
        VAL_AT* pVA = pVArray + iat;
        int nPeriodicNum = at->el_number;
        int cur_chem_valence, cur_chem_valence_fixed, min_chem_valence, known_chem_valence;
        int metal_bonds_chem_valence, not_metal_bonds_chem_valence, alt_bonds_delta_valence, bonds_chem_valence, bond_type;
        CHARGE_VAL  ChargeVal[NUM_ATOM_CHARGES * MAX_NUM_VALENCES];

        memset(ChargeVal, 0, sizeof(ChargeVal));

        pVA->cDoNotAddH = if_skip_add_H(nPeriodicNum);  /* InChI never adds H to this atom */
        /*pVA->cMetal     = is_el_a_metal( nPeriodicNum );*/ /* the atom is a metal */

        /* count bonds to metal atoms; metals have already been marked */
        metal_bonds_chem_valence = not_metal_bonds_chem_valence = alt_bonds_delta_valence = 0;

        if (pVA->cMetal)
        {
            j = at->valence; /* all bonds to metal */
            for (i = k = j2 = k2 = 0; i < at->valence; i++)
            {
                bond_type = (at->bond_type[i] & BOND_TYPE_MASK);
                if (bond_type <= BOND_TYPE_TRIPLE)
                {
                    metal_bonds_chem_valence += inchi_max(BOND_TYPE_SINGLE, bond_type);
                }
                else
                {
                    metal_bonds_chem_valence += BOND_TYPE_SINGLE;
                    k++; /* count alternating bonds */
                }
            }
        }
        else
        {
            for (i = j = j2 = k = k2 = 0; i < at->valence; i++)
            {
                bond_type = (at->bond_type[i] & BOND_TYPE_MASK);
                if (pVArray[(int)at->neighbor[i]].cMetal)
                {
                    j++;  /* number of bonds to metal atoms */
                    if (bond_type <= BOND_TYPE_TRIPLE)
                    {
                        metal_bonds_chem_valence += inchi_max(BOND_TYPE_SINGLE, bond_type);
                    }
                    else
                    {
                        metal_bonds_chem_valence += BOND_TYPE_SINGLE;
                        k++; /* count alternating bonds */
                    }
                }
                else
                {
                    j2++;
                    if (bond_type <= BOND_TYPE_TRIPLE)
                    {
                        not_metal_bonds_chem_valence += inchi_max(BOND_TYPE_SINGLE, bond_type);
                    }
                    else
                    {
                        not_metal_bonds_chem_valence += BOND_TYPE_SINGLE;
                        k2++; /* count alternating bonds */
                    }
                }
            }
        }

        bonds_chem_valence = metal_bonds_chem_valence + not_metal_bonds_chem_valence;

        if (at->chem_bonds_valence > bonds_chem_valence)
        {
            if (at->chem_bonds_valence - bonds_chem_valence > 1)
            {
                at->chem_bonds_valence = bonds_chem_valence + 1;  /* should not happen */
            }
            alt_bonds_delta_valence = at->chem_bonds_valence - bonds_chem_valence;
        }

        pVA->cNumBondsToMetal = j;

        if (nPeriodicNum == EL_NUMBER_H)
        {
            /* ignore bridging H; ??? later add ??? */
            return 0;
        }

        num_H = at->num_H;
        num_bonds = at->valence;

        if (!num_bonds && !num_H)
        {
            return 0; /* do not know the answer: isolated atom */
        }
        /* at the beginning all bonds are single */
        min_chem_valence = num_bonds + num_H;
        cur_chem_valence = bonds_chem_valence + alt_bonds_delta_valence + num_H; /* includes double & alternating bond contribution */

        /* number of non-bonding electrons in case of all single bonds */
        num_non_bonding_electrons = (int)pVA->cNumValenceElectrons - min_chem_valence;
        /* Octet rule: charge = bonds_valence + NumValenceElectrons - 8 */
        charge = min_chem_valence + (int)pVA->cNumValenceElectrons - VALUE_OCTET; /* wrong */

        /* typical (ad hoc) minimal neutral valence */
        known_chem_valence = (pVA->cNumValenceElectrons > VALUE_OCTET / 2) ?
            VALUE_OCTET - pVA->cNumValenceElectrons :
            pVA->cNumValenceElectrons;
        /* excess of typical valence over all-single-bonds valence */
        nOctetNeutralValenceExcess = known_chem_valence - min_chem_valence;
        /*  (NB=num.bonds, NV=neutral valence, NVX=neutral valence excess, LFVS=last found valence state, val.=valence)

           element NB  knownFst octet Last  octetNVX  firstNVX foundNVX chargeLFVS  LFVS
                       valence  val.  NV>=

           -B       1    3        3    3        2         2   =     2        +2
           >B       2    3        3    3        1         1   =     1        +1
           >B-      3    3        3    3        0         0   =     0         0
           >B<      4    3        3    3       -1        -1   <>   N/A       -1

           -C       1    4        4    4        3         3   =     3        N/A
           >C       2    4        4    4        2         2   =     2        +2  (-2)
           >C-      3    4        4    4        1         1   =     1        +1  (-1)
           >C<      4    4        4    4        0         0   =     0         0
           C(V)     5    4        4    N/A     -1        -1   <>   N/A       N/A

           -Si      1    4        4    4        3         3   =     3        N/A
           >Si      2    4        4    4        2         2   =     2        +2  (-2)
           >Si-     3    4        4    4        1         1   =     1        +1  (-1)
           >Si-     4    4        4    4        0         0   =     0         0
           Si(V)    5    4        4    N/A     -1        -1   <>   N/A       -1

           -N       1    3        3    3        2         2   =     2        -2
           >N       2    3        3    3        1         1   =     1        -1
           >N-      3    3        3    3        0         0   =     0         0  (+2)
           >N<      4    3        3    5       -1        -1   <>    1        +1
           N(V)     5    3        3    5       -2        -2   <>    0         0
           N(VI)    6    3        3    N/A     -3        -3   <>   N/A       N/A
           N(VII)   7    3        3    N/A     -4        -4   <>   N/A       N/A

           -P       1    3        3    3        2         2   =     2        -2
           >P       2    3        3    3        1         1   =     1        -1
           >P-      3    3        3    3        0         0   =     0         0  (-2, +2)
           >P<      4    3        3    5       -1        -1   <>    1        +1  (-1)
           P(V)     5    3        3    5       -2        -2   <>    0         0  (-2)
           P(VI)    6    3        3    N/A     -3        -3   <>   N/A       -1
           P(VII)   7    3        3    N/A     -4        -4   <>   N/A       -2
           P(VIII)  8    3        3    N/A     -5        -5   <>   N/A       N/A

           -O       1    2        2    2        1         1   =     1        -1
           >O       2    2        2    2        0         0   =     0         0
           >O-      3    2        2    N/A     -1        -1   <>   N/A       +1
           >O<      4    2        2    N/A     -2        -2   <>   N/A       +2
           O(V)     5    2        2    N/A     -3        -3   <>   N/A       +1
           O(VI)    6    2        2    N/A     -4        -4   <>   N/A       N/A

           -S       1    2        2    2        1         1   =     1        -1
           >S       2    2        2    2        0         0   =     0         0       NPNP - prohibit
           >S-      3    2        2    4       -1        -1   <>    1        +1  (-1) PNPN
           >S<      4    2        2    4       -2        -2   <>    0         0  (+2)
           S(V)     5    2        2    6       -3        -3   <>    1        +1  (-1)
           S(VI)    6    2        2    6       -4        -4   <>    0         0
           S(VII)   7    2        2    N/A     -5        -5   <>    0        -1
           S(VIII)  8    2        2    N/A     -6        -6   <>   N/A       N/A

           -F       1    1        1    1        0         0   =     0         0
           >F       2    1        1    1       -1        -1   <>   N/A       +1
           >F-      3    1        1    1       -2        -2   <>   N/A       +2
           >F<      4    1        1    1       -3        -3   <>   N/A       N/A
           F(V)     5    1        1    1       -4        -4   <>   N/A       +2
           F(VI)    6    1        1    1       -5        -5   <>   N/A       N/A

           -Cl      1    1        1    1        0         0   =     0         0      NPNP - prohibit
           >Cl      2    1        1    3       -1        -1   <>    1        +1      PNPN - prohibit
           >Cl-     3    1        1    3       -2        -2   <>    0         0 (+2) NPNP
           >Cl<     4    1        1    5       -3        -3   <>    1        +1      PNPN
           Cl(V)    5    1        1    5       -4        -4   <>    0         0
           Cl(VI)   6    1        1    7       -5        -5   <>    1        +1
           Cl(VII)  7    1        1    7       -6        -6   <>    0         0
           Cl(VIII) 8    1        1    N/A     -7        -7   <>   N/A       N/A


          NB                 = num_bonds+num_H

          knownFst valence   = nFirstNeutralValenceExcess + min_chem_valence
          octet val.         = nOctetNeutralValenceExcess + min_chem_valence
          Last NV>=          = nFoundNeutralValenceExcess + min_chem_valence

          octetNVX           = nOctetNeutralValenceExcess
          firstNVX           = nFirstNeutralValenceExcess
          foundNVX           = nFoundNeutralValenceExcess

          chargeLFVS         = ChargeVal[nLastFoundValenceState].nCharge

        */
        /* minimal known neutral atom valence; different for Sn(2/4), Tl(1/3), Pb(2/4): (known/typical ad hoc) */
        known_chem_valence = get_el_valence( nPeriodicNum, 0, 0 );

        if (pSrm->bMetalAddFlower)
        {
            /* bond orders of bonds to metal may be as they are (pSrm->nMetalInitBondOrder==1)
                                            or decreased by one (pSrm->nMetalInitBondOrder==0)
               nMetalInitBondOrder == nMetalMinBondOrder + nMetalInitEdgeFlow
            */
            cur_chem_valence_fixed = cur_chem_valence - pVA->cNumBondsToMetal * ( 1 - pSrm->nMetalInitBondOrder );
            pVA->cInitOrigValenceToMetal = metal_bonds_chem_valence;
            pVA->cInitValenceToMetal = metal_bonds_chem_valence - pVA->cNumBondsToMetal * ( 1 - pSrm->nMetalInitBondOrder );
            pVA->cInitFlowToMetal = pVA->cInitValenceToMetal - pVA->cNumBondsToMetal * pSrm->nMetalMinBondOrder;
            if (pVA->cMetal)
            {
                pVA->cInitFreeValences += alt_bonds_delta_valence;
            }
            if (pSrm->nMetalInitEdgeFlow < pSrm->nMetalInitBondOrder - pSrm->nMetalMinBondOrder)
            {
                /* single bond has zero initial flow + 2 radicals at incident atoms */
                if (pVA->cInitFlowToMetal <= pVA->cNumBondsToMetal)
                {
                    if (pVA->cMetal)
                    {
                        pVA->cInitFreeValences += pVA->cInitFlowToMetal;
                    }
                    pVA->cInitFlowToMetal = 0;
                }
                else
                {
                    if (pVA->cMetal)
                    {
                        pVA->cInitFreeValences += pVA->cNumBondsToMetal * ( 1 - pSrm->nMetalInitEdgeFlow );
                    }
                    pVA->cInitFlowToMetal -= pVA->cNumBondsToMetal * ( 1 - pSrm->nMetalInitEdgeFlow );
                }
            }
        }
        else
        {
            /* treat metal atoms as ordinary non-metal atoms */
            cur_chem_valence_fixed = cur_chem_valence;
            pVA->cInitFlowToMetal = metal_bonds_chem_valence - pVA->cNumBondsToMetal;
            pVA->cInitValenceToMetal = metal_bonds_chem_valence;
            pVA->cInitOrigValenceToMetal = metal_bonds_chem_valence;
        }


        if (pVA->cMetal && pSrm->bMetalAddFlower)
        {
            pVA->cnListIndex = cnListIndexMe + 1;
            /*
            pVA->cInitOrigValenceToMetal += alt_bonds_delta_valence;
            pVA->cInitValenceToMetal     += alt_bonds_delta_valence;
            pVA->cInitFreeValences = (pSrm->nMetalInitBondOrder + alt_bonds_delta_valence
                                     - (pSrm->nMetalMinBondOrder + pSrm->nMetalInitEdgeFlow)) * pVA->cNumBondsToMetal;
            */
            return 0;  /* metal */
        }

        if (!known_chem_valence)
        {
            /* a noble gas like He, Ne, ... */
            pVA->cInitFreeValences = at->chem_bonds_valence - at->valence;
            return TREAT_ATOM_AS_METAL; /* do not know anything about this atom; needs 2nd pass */
        }

        nFirstNeutralValenceExcess = known_chem_valence - min_chem_valence;

        nFoundNeutralValenceExcess = NO_VALUE_INT;
        nFoundNeutralValenceOrdNumber = NO_VALUE_INT;
        nLastFoundValenceOrdNumber = NO_VALUE_INT;
        nLastFoundValenceState = NO_VALUE_INT;

        /* find the lowest known valence >= all-single-bonds valence */
        for (cur_charge = MIN_ATOM_CHARGE, nNumStates = 0; cur_charge <= MAX_ATOM_CHARGE; cur_charge++)
        {
            for (i = 0; i < MAX_NUM_VALENCES; i++)
            {
                known_chem_valence = get_el_valence(nPeriodicNum, cur_charge, i);
                if (cur_chem_valence_fixed > known_chem_valence || !known_chem_valence)
                {
                    continue; /* known valence < all-single-bonds valence */
                }
                if (BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * (num_bonds - 1) + num_H < known_chem_valence)
                {
                    continue; /* not more than one triple and the rest - double bonds per atom */
                }

                /* keep all found */
                ChargeVal[nNumStates].nValence = known_chem_valence;
                ChargeVal[nNumStates].nCharge = cur_charge;
                ChargeVal[nNumStates].nValenceOrderingNumber = i;
                if (!cur_charge && nFoundNeutralValenceExcess == NO_VALUE_INT)
                {
                    /* neutral state; compare to the lowest typical valence */
                    nFoundNeutralValenceExcess = known_chem_valence - min_chem_valence;
                    nFoundNeutralValenceOrdNumber = i;
                }
                if (min_chem_valence == known_chem_valence)
                {
                    if (nLastFoundValenceState == NO_VALUE_INT)
                    {
                        /* accept the first found */
                        nLastFoundValenceState = nNumStates;
                    }
                    else
                        if (abs(ChargeVal[nLastFoundValenceState].nCharge) >= abs(cur_charge))
                        {
                            /* accept smaller abs(charge); if abs(charges) are same, accept (+) */
                            nLastFoundValenceState = nNumStates;
                        }
                }

                nNumStates++;
            }
        }

        /***********************************************************************************/
        /* select only appropriate charge & valence so that a suitable ChargeStruct exists */
        /***********************************************************************************/

        nNumSelectedStates = clean_charge_val( pCG, ChargeVal, nNumStates, atom,
                                               pVArray, iat, pVA->cMetal,
                                               bMobileH, endpoint );

        if (!nNumSelectedStates)
        {
            return TREAT_ATOM_AS_METAL; /* nothing to do */
        }

        /***********************************************************************************/
        /*       Find an appropriate ChargeStruct index for the ChargeVal found            */
        /***********************************************************************************/
        cn_bits = 0;
        memset( cn_bits_array, 0, sizeof( cn_bits_array ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        /***** set bits identifying a suitable ChargeStruct ******/
        for (i = len_cn_bits_array = 0; i < nNumSelectedStates && len_cn_bits_array < 4; i++)
        {
            switch (ChargeVal[i].nCharge)
            {
                case -1:
                    cn_bits_array[len_cn_bits_array] |= cn_bits_M; /* Minus 1 */
                    break;
                case 0:
                    cn_bits_array[len_cn_bits_array] |= cn_bits_N; /* Neutral */
                    break;
                case 1:
                    cn_bits_array[len_cn_bits_array] |= cn_bits_P; /* Plus 1 */
                    break;
                default:
                    return RI_ERR_PROGR; /* program error */
            }
            if (i + 1 < nNumSelectedStates &&
                 ChargeVal[i].nValence == ChargeVal[i + 1].nValence &&
                 ChargeVal[i].nCharge &&
                 ChargeVal[i].nCharge == -ChargeVal[i + 1].nCharge)
            {
                ; /* add opposite charge to the same element of cn_bits_array[] */
            }
            else
            {
                len_cn_bits_array++;
            }
        }
        if (!len_cn_bits_array || len_cn_bits_array > 4)
        {
            return RI_ERR_PROGR; /* program error */
        }
        /* accommodate added 4-state ChargeStruct: +/- cannot be in case of 4 states */
        if (len_cn_bits_array + 1 == nNumSelectedStates && nNumSelectedStates == 4)
        {
            len_cn_bits_array--;
            nNumSelectedStates--;
            cn_bits_array[len_cn_bits_array] = 0;
        }
        /* fix for terminal hydrogenless -C as in isocyano or CO: there is no just cnE_[] ChargeStruct */
        if (len_cn_bits_array == 1 &&
             cn_bits_array[0] == ( cn_bits_P | cn_bits_M ) &&
             ChargeVal[0].nValence + 1 > BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * ( num_bonds - 1 ) + num_H)
        {
            cn_bits_array[len_cn_bits_array++] = cn_bits_N;
            ChargeVal[nNumSelectedStates].nValence = ChargeVal[nNumSelectedStates - 1].nValence;
            ChargeVal[nNumSelectedStates].nCharge = 0;
            ChargeVal[nNumSelectedStates].nValenceOrderingNumber = 0;
        }

    make_cn_bits:
        cn_bits = MAKE_CN_BITS( cn_bits_array[0], cn_bits_array[1], cn_bits_array[2], cn_bits_array[3] );

        /*********** find ChargeStructure **************/
        for (i = 0, j = -1; i < cnListNumEl; i++)
        {
            if (cnList[i].bits == cn_bits)
            {
                j = i;
                break; /* found */
            }
        }

        if (j < 0)
        {
            /* ChargeStructure was not found */
            if (1 < len_cn_bits_array && len_cn_bits_array + 1 == nNumSelectedStates)
            {
    /* a pair of opposite charges was combined */
                len_cn_bits_array--;
                cn_bits_array[len_cn_bits_array] = 0;
                goto make_cn_bits;
            }
            else if (nNumSelectedStates == 4)
            {
                /* reduce number of states */
                len_cn_bits_array--;
                cn_bits_array[len_cn_bits_array] = 0;
                nNumSelectedStates--;
                goto make_cn_bits;
            }
            return RI_ERR_PROGR; /* charge structure not found */
        }

        /********** ChargeStructure has been found **********/
        pVA->cnListIndex = j + 1; /* charge structure index + 1 */
        pVA->cInitCharge = cnList[j].nInitialCharge;
        /********** Calculate "Free Valence" ****************/
    #if ( ALLOW_METAL_BOND_ZERO == 1 )

    #if ( INIT_METAL_BOND_ZERO == 1 )
        if (pVA->cMetal)
        {
            j = 0;
        }
        else
        {
            j = ChargeVal[0].nValence - cur_chem_valence_fixed;
        }
    #else
        j = ChargeVal[0].nValence - cur_chem_valence_fixed;
    #endif

    #else
        j = ChargeVal[0].nValence - cur_chem_valence_fixed;
    #endif
        if (j < 0)
        {
            return RI_ERR_PROGR; /* program error */
        }
        pVA->cInitFreeValences = j; /* number of initial unsatisfied valences; should be combined with */
                                    /* (cap - flow) of vertex=0 in the charge structure[pVA->cnListIndex-1] */
        return 1; /* success */

    #undef MIN_ATOM_CHARGE
    #undef MAX_ATOM_CHARGE
    #undef NEUTRAL_STATE
    #undef NUM_ATOM_CHARGES
    #undef MAX_NUM_VALENCES
    }
        */
    // END INCHI C FUNCTION: GetAtomRestoreInfo
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetAtomRestoreInfo
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️✔️: ALLOW_METAL_BOND_ZERO=1; INIT_METAL_BOND_ZERO=0.
    // INCHI✔️✔️: MIN_ATOM_CHARGE=-2; MAX_ATOM_CHARGE=2; MAX_NUM_VALENCES=5.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetAtomRestoreInfo

    const MIN_ATOM_CHARGE: i32 = -2;
    const MAX_ATOM_CHARGE: i32 = 2;
    const MAX_NUM_VALENCES: i32 = 5;

    let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if index >= atoms.len() || index >= valence_atoms.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let periodic_number = i32::from(atoms[index].el_number);
    valence_atoms[index].cDoNotAddH = crate::source::base::util::if_skip_add_H(periodic_number)? as i8;

    let atom_valence = i32::from(atoms[index].valence);
    let current_is_metal = valence_atoms[index].cMetal != 0;
    let mut metal_bonds_chemical_valence = 0_i32;
    let mut nonmetal_bonds_chemical_valence = 0_i32;
    let mut alternating_bonds_delta_valence = 0_i32;
    let mut bonds_to_metal;
    let mut alternating_metal_bonds = 0_i32;
    let mut alternating_nonmetal_bonds = 0_i32;

    if current_is_metal {
        bonds_to_metal = atom_valence;
        let mut ordinal = 0_i32;
        while ordinal < atom_valence {
            let bond_type = i32::from(
                atoms[index].bond_type[usize::try_from(ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?],
            ) & crate::source_types::BOND_TYPE_MASK as i32;
            if bond_type <= crate::source_types::BOND_TYPE_TRIPLE as i32 {
                metal_bonds_chemical_valence = metal_bonds_chemical_valence
                    .wrapping_add((crate::source_types::BOND_TYPE_SINGLE as i32).max(bond_type));
            } else {
                metal_bonds_chemical_valence =
                    metal_bonds_chemical_valence.wrapping_add(crate::source_types::BOND_TYPE_SINGLE as i32);
                alternating_metal_bonds = alternating_metal_bonds.wrapping_add(1);
            }
            ordinal = ordinal.wrapping_add(1);
        }
    } else {
        bonds_to_metal = 0;
        let mut nonmetal_bond_count = 0_i32;
        let mut ordinal = 0_i32;
        while ordinal < atom_valence {
            let ordinal_index = usize::try_from(ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let bond_type =
                i32::from(atoms[index].bond_type[ordinal_index]) & crate::source_types::BOND_TYPE_MASK as i32;
            let neighbor = usize::from(atoms[index].neighbor[ordinal_index]);
            let neighbor_is_metal = valence_atoms
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .cMetal
                != 0;
            if neighbor_is_metal {
                bonds_to_metal = bonds_to_metal.wrapping_add(1);
                if bond_type <= crate::source_types::BOND_TYPE_TRIPLE as i32 {
                    metal_bonds_chemical_valence = metal_bonds_chemical_valence
                        .wrapping_add((crate::source_types::BOND_TYPE_SINGLE as i32).max(bond_type));
                } else {
                    metal_bonds_chemical_valence =
                        metal_bonds_chemical_valence.wrapping_add(crate::source_types::BOND_TYPE_SINGLE as i32);
                    alternating_metal_bonds = alternating_metal_bonds.wrapping_add(1);
                }
            } else {
                nonmetal_bond_count = nonmetal_bond_count.wrapping_add(1);
                if bond_type <= crate::source_types::BOND_TYPE_TRIPLE as i32 {
                    nonmetal_bonds_chemical_valence = nonmetal_bonds_chemical_valence
                        .wrapping_add((crate::source_types::BOND_TYPE_SINGLE as i32).max(bond_type));
                } else {
                    nonmetal_bonds_chemical_valence =
                        nonmetal_bonds_chemical_valence.wrapping_add(crate::source_types::BOND_TYPE_SINGLE as i32);
                    alternating_nonmetal_bonds = alternating_nonmetal_bonds.wrapping_add(1);
                }
            }
            ordinal = ordinal.wrapping_add(1);
        }
        let _ = nonmetal_bond_count;
    }
    let _ = (alternating_metal_bonds, alternating_nonmetal_bonds);

    let bonds_chemical_valence = metal_bonds_chemical_valence.wrapping_add(nonmetal_bonds_chemical_valence);
    if i32::from(atoms[index].chem_bonds_valence) > bonds_chemical_valence {
        if i32::from(atoms[index].chem_bonds_valence).wrapping_sub(bonds_chemical_valence) > 1 {
            atoms[index].chem_bonds_valence = bonds_chemical_valence.wrapping_add(1) as i8;
        }
        alternating_bonds_delta_valence =
            i32::from(atoms[index].chem_bonds_valence).wrapping_sub(bonds_chemical_valence);
    }
    valence_atoms[index].cNumBondsToMetal = bonds_to_metal as i8;

    if periodic_number == 1 {
        return Ok(0);
    }

    let number_hydrogens = i32::from(atoms[index].num_H);
    let number_bonds = atom_valence;
    if number_bonds == 0 && number_hydrogens == 0 {
        return Ok(0);
    }

    let minimum_chemical_valence = number_bonds.wrapping_add(number_hydrogens);
    let current_chemical_valence = bonds_chemical_valence
        .wrapping_add(alternating_bonds_delta_valence)
        .wrapping_add(number_hydrogens);
    let _number_nonbonding_electrons =
        i32::from(valence_atoms[index].cNumValenceElectrons).wrapping_sub(minimum_chemical_valence);
    let _charge = minimum_chemical_valence
        .wrapping_add(i32::from(valence_atoms[index].cNumValenceElectrons))
        .wrapping_sub(crate::source_types::VALUE_OCTET as i32);
    let _octet_neutral_valence_excess =
        if i32::from(valence_atoms[index].cNumValenceElectrons) > crate::source_types::VALUE_OCTET as i32 / 2 {
            (crate::source_types::VALUE_OCTET as i32).wrapping_sub(i32::from(valence_atoms[index].cNumValenceElectrons))
        } else {
            i32::from(valence_atoms[index].cNumValenceElectrons)
        }
        .wrapping_sub(minimum_chemical_valence);

    let mut known_chemical_valence = crate::source::base::util::get_el_valence(periodic_number, 0, 0)?;
    let current_chemical_valence_fixed;
    if restore_mode.bMetalAddFlower != 0 {
        current_chemical_valence_fixed = current_chemical_valence.wrapping_sub(
            i32::from(valence_atoms[index].cNumBondsToMetal)
                .wrapping_mul(1_i32.wrapping_sub(restore_mode.nMetalInitBondOrder)),
        );
        valence_atoms[index].cInitOrigValenceToMetal = metal_bonds_chemical_valence as i8;
        valence_atoms[index].cInitValenceToMetal = metal_bonds_chemical_valence.wrapping_sub(
            i32::from(valence_atoms[index].cNumBondsToMetal)
                .wrapping_mul(1_i32.wrapping_sub(restore_mode.nMetalInitBondOrder)),
        ) as i8;
        valence_atoms[index].cInitFlowToMetal = i32::from(valence_atoms[index].cInitValenceToMetal).wrapping_sub(
            i32::from(valence_atoms[index].cNumBondsToMetal).wrapping_mul(restore_mode.nMetalMinBondOrder),
        ) as i8;
        if current_is_metal {
            valence_atoms[index].cInitFreeValences = valence_atoms[index]
                .cInitFreeValences
                .wrapping_add(alternating_bonds_delta_valence as i8);
        }
        if restore_mode.nMetalInitEdgeFlow
            < restore_mode
                .nMetalInitBondOrder
                .wrapping_sub(restore_mode.nMetalMinBondOrder)
        {
            if i32::from(valence_atoms[index].cInitFlowToMetal) <= i32::from(valence_atoms[index].cNumBondsToMetal) {
                if current_is_metal {
                    valence_atoms[index].cInitFreeValences = valence_atoms[index]
                        .cInitFreeValences
                        .wrapping_add(valence_atoms[index].cInitFlowToMetal);
                }
                valence_atoms[index].cInitFlowToMetal = 0;
            } else {
                if current_is_metal {
                    let delta = i32::from(valence_atoms[index].cNumBondsToMetal)
                        .wrapping_mul(1_i32.wrapping_sub(restore_mode.nMetalInitEdgeFlow));
                    valence_atoms[index].cInitFreeValences =
                        valence_atoms[index].cInitFreeValences.wrapping_add(delta as i8);
                }
                let delta = i32::from(valence_atoms[index].cNumBondsToMetal)
                    .wrapping_mul(1_i32.wrapping_sub(restore_mode.nMetalInitEdgeFlow));
                valence_atoms[index].cInitFlowToMetal =
                    i32::from(valence_atoms[index].cInitFlowToMetal).wrapping_sub(delta) as i8;
            }
        }
    } else {
        current_chemical_valence_fixed = current_chemical_valence;
        valence_atoms[index].cInitFlowToMetal =
            metal_bonds_chemical_valence.wrapping_sub(i32::from(valence_atoms[index].cNumBondsToMetal)) as i8;
        valence_atoms[index].cInitValenceToMetal = metal_bonds_chemical_valence as i8;
        valence_atoms[index].cInitOrigValenceToMetal = metal_bonds_chemical_valence as i8;
    }

    if current_is_metal && restore_mode.bMetalAddFlower != 0 {
        valence_atoms[index].cnListIndex = (crate::source_types::cnListIndexMe as i32).wrapping_add(1) as i8;
        return Ok(0);
    }

    if known_chemical_valence == 0 {
        valence_atoms[index].cInitFreeValences =
            i32::from(atoms[index].chem_bonds_valence).wrapping_sub(atom_valence) as i8;
        return Ok(crate::source_types::TREAT_ATOM_AS_METAL as i32);
    }

    let _first_neutral_valence_excess = known_chemical_valence.wrapping_sub(minimum_chemical_valence);
    let mut charge_values: [crate::source_types::CHARGE_VAL; 25] =
        std::array::from_fn(|_| crate::source_types::CHARGE_VAL::default());
    let mut number_states = 0_i32;
    let mut found_neutral_valence_excess = crate::source_types::NO_VALUE_INT as i32;
    let mut _found_neutral_valence_order = crate::source_types::NO_VALUE_INT as i32;
    let mut last_found_valence_state = crate::source_types::NO_VALUE_INT as i32;

    let mut current_charge = MIN_ATOM_CHARGE;
    while current_charge <= MAX_ATOM_CHARGE {
        let mut valence_number = 0_i32;
        while valence_number < MAX_NUM_VALENCES {
            known_chemical_valence =
                crate::source::base::util::get_el_valence(periodic_number, current_charge, valence_number)?;
            if current_chemical_valence_fixed > known_chemical_valence || known_chemical_valence == 0 {
                valence_number = valence_number.wrapping_add(1);
                continue;
            }
            if (crate::source_types::BOND_TYPE_TRIPLE as i32)
                .wrapping_add((crate::source_types::BOND_TYPE_DOUBLE as i32).wrapping_mul(number_bonds.wrapping_sub(1)))
                .wrapping_add(number_hydrogens)
                < known_chemical_valence
            {
                valence_number = valence_number.wrapping_add(1);
                continue;
            }

            let state_index = usize::try_from(number_states).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let state = charge_values
                .get_mut(state_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            state.nValence = known_chemical_valence;
            state.nCharge = current_charge;
            state.nValenceOrderingNumber = valence_number;
            if current_charge == 0 && found_neutral_valence_excess == crate::source_types::NO_VALUE_INT as i32 {
                found_neutral_valence_excess = known_chemical_valence.wrapping_sub(minimum_chemical_valence);
                _found_neutral_valence_order = valence_number;
            }
            if minimum_chemical_valence == known_chemical_valence {
                if last_found_valence_state == crate::source_types::NO_VALUE_INT as i32 {
                    last_found_valence_state = number_states;
                } else {
                    let last_index =
                        usize::try_from(last_found_valence_state).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    if charge_values[last_index].nCharge.wrapping_abs() >= current_charge.wrapping_abs() {
                        last_found_valence_state = number_states;
                    }
                }
            }
            number_states = number_states.wrapping_add(1);
            valence_number = valence_number.wrapping_add(1);
        }
        current_charge = current_charge.wrapping_add(1);
    }

    let mut number_selected_states = clean_charge_val(
        canonical_globals,
        &mut charge_values,
        number_states,
        atoms,
        valence_atoms,
        atom_index,
        i32::from(valence_atoms[index].cMetal),
        mobile_h,
        endpoint,
    )?;
    if number_selected_states == 0 {
        return Ok(crate::source_types::TREAT_ATOM_AS_METAL as i32);
    }

    let mut charge_bits = [0_i32; 5];
    let mut bits_length = 0_i32;
    let mut selected_index = 0_i32;
    while selected_index < number_selected_states && bits_length < 4 {
        let state_index = usize::try_from(selected_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let bits_index = usize::try_from(bits_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        charge_bits[bits_index] |= match charge_values[state_index].nCharge {
            -1 => crate::source_types::cn_bits_M as i32,
            0 => crate::source_types::cn_bits_N as i32,
            1 => crate::source_types::cn_bits_P as i32,
            _ => return Ok(RI_ERR_PROGR),
        };
        if selected_index.wrapping_add(1) < number_selected_states {
            let next = &charge_values
                [usize::try_from(selected_index.wrapping_add(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?];
            if charge_values[state_index].nValence == next.nValence
                && charge_values[state_index].nCharge != 0
                && charge_values[state_index].nCharge == next.nCharge.wrapping_neg()
            {
                selected_index = selected_index.wrapping_add(1);
                continue;
            }
        }
        bits_length = bits_length.wrapping_add(1);
        selected_index = selected_index.wrapping_add(1);
    }
    if bits_length == 0 || bits_length > 4 {
        return Ok(RI_ERR_PROGR);
    }
    if bits_length.wrapping_add(1) == number_selected_states && number_selected_states == 4 {
        bits_length = bits_length.wrapping_sub(1);
        number_selected_states = number_selected_states.wrapping_sub(1);
        charge_bits[usize::try_from(bits_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = 0;
    }
    if bits_length == 1
        && charge_bits[0] == (crate::source_types::cn_bits_P | crate::source_types::cn_bits_M) as i32
        && charge_values[0].nValence.wrapping_add(1)
            > (crate::source_types::BOND_TYPE_TRIPLE as i32)
                .wrapping_add((crate::source_types::BOND_TYPE_DOUBLE as i32).wrapping_mul(number_bonds.wrapping_sub(1)))
                .wrapping_add(number_hydrogens)
    {
        charge_bits[usize::try_from(bits_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?] =
            crate::source_types::cn_bits_N as i32;
        bits_length = bits_length.wrapping_add(1);
        let write_index = usize::try_from(number_selected_states).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        charge_values[write_index].nValence = charge_values[usize::try_from(number_selected_states.wrapping_sub(1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
        .nValence;
        charge_values[write_index].nCharge = 0;
        charge_values[write_index].nValenceOrderingNumber = 0;
    }

    let make_charge_bits = |values: &[i32; 5]| {
        (((values[3] << crate::source_types::cn_bits_shift) | values[2]) << crate::source_types::cn_bits_shift
            | values[1])
            << crate::source_types::cn_bits_shift
            | values[0]
    };
    let list_index = loop {
        let combined_bits = make_charge_bits(&charge_bits);
        if let Some(found) = CN_LIST.iter().position(|entry| entry.bits == combined_bits) {
            break found;
        }
        if bits_length > 1 && bits_length.wrapping_add(1) == number_selected_states {
            bits_length = bits_length.wrapping_sub(1);
            charge_bits[usize::try_from(bits_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = 0;
        } else if number_selected_states == 4 {
            bits_length = bits_length.wrapping_sub(1);
            charge_bits[usize::try_from(bits_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = 0;
            number_selected_states = number_selected_states.wrapping_sub(1);
        } else {
            return Ok(RI_ERR_PROGR);
        }
    };

    valence_atoms[index].cnListIndex = list_index.wrapping_add(1) as i8;
    valence_atoms[index].cInitCharge = CN_LIST_INITIAL_CHARGES[list_index];
    let free_valences = charge_values[0].nValence.wrapping_sub(current_chemical_valence_fixed);
    if free_valences < 0 {
        return Ok(RI_ERR_PROGR);
    }
    valence_atoms[index].cInitFreeValences = free_valences as i8;
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn GetTgroupInfoFromInChI(
    heap: &mut SourceHeap,
    ti: &mut T_GROUP_INFO,
    at: SourceMutPointer<inp_ATOM>,
    endpoint: SourceMutPointer<AT_NUMB>,
    p_inchi: Option<&INChI>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:575 GetTgroupInfoFromInChI
    // INCHI✔️❌: int GetTgroupInfoFromInChI( T_GROUP_INFO *ti,
    // INCHI✔️❌:                             inp_ATOM *at,
    // INCHI✔️❌:                             AT_NUMB *endpoint,
    // INCHI✔️❌:                             INChI *pInChI )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret, i, j, k, itg, num_atoms, len_tg, num_t_groups; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     AT_NUMB   *tGroupNumber = NULL;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     AT_NUMB   *tiGroupNumber = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     clear_t_group_info( ti );
    // INCHI✔️❌:     if (pInChI && pInChI->lenTautomer > 1 &&
    // INCHI✔️❌:          pInChI->nTautomer && pInChI->nTautomer[0] > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_atoms = pInChI->nNumberOfAtoms;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         num_t_groups = pInChI->nTautomer[0];
    // INCHI✔️❌:         len_tg = pInChI->lenTautomer - T_GROUP_HDR_LEN*pInChI->nTautomer[0] - 1; /* number of endpoints */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* allocation ti->t_group */
    // INCHI✔️❌:         if (ti->max_num_t_groups != num_atoms / 2 + 1 || !ti->t_group)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ti->max_num_t_groups = num_atoms / 2 + 1;
    // INCHI✔️❌:             if (ti->t_group)
    // INCHI✔️❌:                 inchi_free( ti->t_group );
    // INCHI✔️❌:             ti->t_group = (T_GROUP *) inchi_calloc( (long long)ti->max_num_t_groups + 1, sizeof( ti->t_group[0] ) ); /* djb-rwth: correcting 0 bytes allocation */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* allocation ti->tGroupNumber */
    // INCHI✔️❌:         if (ti->num_t_groups != num_t_groups || !ti->tGroupNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ti->num_t_groups = num_t_groups;
    // INCHI✔️❌:             if (ti->tGroupNumber)
    // INCHI✔️❌:                 inchi_free( ti->tGroupNumber );
    // INCHI✔️❌:             ti->tGroupNumber = (AT_NUMB *) inchi_calloc( ( (long long)ti->num_t_groups + 1 )*TGSO_TOTAL_LEN, sizeof( ti->tGroupNumber[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* allocation ti->tGroupNumber */
    // INCHI✔️❌:         if (len_tg != ti->nNumEndpoints || !ti->nEndpointAtomNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ti->nNumEndpoints = len_tg;
    // INCHI✔️❌:             if (ti->nEndpointAtomNumber)
    // INCHI✔️❌:                 inchi_free( ti->nEndpointAtomNumber );
    // INCHI✔️❌:             ti->nEndpointAtomNumber = (AT_NUMB *) inchi_calloc( (long long)len_tg + 1, sizeof( ti->nEndpointAtomNumber[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* check */
    // INCHI✔️❌:         if (!ti->t_group || !ti->tGroupNumber || !ti->nEndpointAtomNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = RI_ERR_ALLOC;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         tGroupNumber = ti->tGroupNumber;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         tiGroupNumber = tGroupNumber + TGSO_SYMM_IORDER * (long long)ti->num_t_groups; /* djb-rwth: cast operator added */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:             j = 1; /* index in pInChI->nTautomer[] */
    // INCHI✔️❌:         i = 0; /* index in ti->nEndpointAtomNumber[] */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* djb-rwth: fixing oss-fuzz issues #67681, #67641 */
    // INCHI✔️❌:         for (itg = 0; itg < pInChI->nTautomer[0] && itg <= ti->max_num_t_groups; itg++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len_tg = pInChI->nTautomer[j]; /* t-group length not including pInChI->nTautomer[j] */
    // INCHI✔️❌:             ti->t_group[itg].num[0] = pInChI->nTautomer[j + 1] + pInChI->nTautomer[j + 2]; /* num mobile H & (-) */
    // INCHI✔️❌:             ti->t_group[itg].num[1] = pInChI->nTautomer[j + 2]; /* num mobile (-) */
    // INCHI✔️❌:             tGroupNumber[itg] = tiGroupNumber[itg] = itg; /* index */
    // INCHI✔️❌:             ti->t_group[itg].nGroupNumber = /*tSymmRank[itg] = tiSymmRank[itg] =*/ itg + 1; /* t-group number */
    // INCHI✔️❌:             j += T_GROUP_HDR_LEN;   /* skip t-group header */
    // INCHI✔️❌:             len_tg -= T_GROUP_HDR_LEN - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:             ti->t_group[itg].nNumEndpoints = len_tg;
    // INCHI✔️❌:             ti->t_group[itg].nFirstEndpointAtNoPos = i;
    // INCHI✔️❌:
    // INCHI✔️❌:             while (len_tg > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k = ti->nEndpointAtomNumber[i] = pInChI->nTautomer[j] - 1; /* djb-rwth: buffer overrun avoided implicitly in loop condition */
    // INCHI✔️❌: #if ( FIX_GAF_2019_1==1 )
    // INCHI✔️❌:                 if (k<0 || k>num_atoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = RI_ERR_PROGR;
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 if (at)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at[k].endpoint = itg + 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (endpoint)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     endpoint[k] = itg + 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 len_tg--;
    // INCHI✔️❌:                 j++;
    // INCHI✔️❌:                 i++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (i != ti->nNumEndpoints)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetTgroupInfoFromInChI
    // BEGIN ACTIVE INCHI MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:396
    // INCHI✔️❌: #define FIX_GAF_2019_1 1
    // END ACTIVE INCHI MACRO: FIX_GAF_2019_1

    clear_t_group_info(heap, Some(ti))?;
    let Some(p_inchi) = p_inchi else {
        return Ok(0);
    };
    if p_inchi.lenTautomer <= 1 || p_inchi.nTautomer.is_null() {
        return Ok(0);
    }
    let tautomer = heap.slice(p_inchi.nTautomer.as_const())?.to_vec();
    let num_t_groups = i32::from(*tautomer.first().ok_or(SourceHeapError::PointerOutOfBounds)?);
    if num_t_groups <= 0 {
        return Ok(0);
    }
    let num_atoms = p_inchi.nNumberOfAtoms;
    let mut len_tg = p_inchi
        .lenTautomer
        .wrapping_sub((T_GROUP_HDR_LEN as i32).wrapping_mul(num_t_groups))
        .wrapping_sub(1);

    let required_max_groups = num_atoms / 2 + 1;
    if ti.max_num_t_groups != required_max_groups || ti.t_group.is_null() {
        ti.max_num_t_groups = required_max_groups;
        if !ti.t_group.is_null() {
            inchi_free(heap, ti.t_group)?;
            ti.t_group = SourceMutPointer::null();
        }
        let count = u64::try_from(i64::from(ti.max_num_t_groups) + 1)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        ti.t_group = match inchi_calloc(heap, count, 40) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }
    if ti.num_t_groups != num_t_groups || ti.tGroupNumber.is_null() {
        ti.num_t_groups = num_t_groups;
        if !ti.tGroupNumber.is_null() {
            inchi_free(heap, ti.tGroupNumber)?;
            ti.tGroupNumber = SourceMutPointer::null();
        }
        let count = (i64::from(ti.num_t_groups) + 1)
            .checked_mul(i64::from(TGSO_TOTAL_LEN))
            .and_then(|value| u64::try_from(value).ok())
            .ok_or(SourceHeapError::AllocationElementCountOutOfRange)?;
        ti.tGroupNumber = match inchi_calloc(heap, count, 2) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }
    if len_tg != ti.nNumEndpoints || ti.nEndpointAtomNumber.is_null() {
        ti.nNumEndpoints = len_tg;
        if !ti.nEndpointAtomNumber.is_null() {
            inchi_free(heap, ti.nEndpointAtomNumber)?;
            ti.nEndpointAtomNumber = SourceMutPointer::null();
        }
        let count =
            u64::try_from(i64::from(len_tg) + 1).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        ti.nEndpointAtomNumber = match inchi_calloc(heap, count, 2) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }
    if ti.t_group.is_null() || ti.tGroupNumber.is_null() || ti.nEndpointAtomNumber.is_null() {
        return Ok(RI_ERR_ALLOC);
    }

    let mut j = 1_i32;
    let mut i = 0_i32;
    let mut itg = 0_i32;
    while itg < num_t_groups && itg <= ti.max_num_t_groups {
        let j_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        len_tg = i32::from(*tautomer.get(j_index).ok_or(SourceHeapError::PointerOutOfBounds)?);
        let mobile_h = *tautomer.get(j_index + 1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mobile_minus = *tautomer.get(j_index + 2).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let group_index = usize::try_from(itg).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        {
            let group = heap
                .slice_mut(ti.t_group)?
                .get_mut(group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            group.num[0] = mobile_h.wrapping_add(mobile_minus);
            group.num[1] = mobile_minus;
            group.nGroupNumber = itg.wrapping_add(1) as AT_NUMB;
        }
        let inverse_base = (TGSO_SYMM_IORDER as i64)
            .checked_mul(i64::from(ti.num_t_groups))
            .and_then(|value| usize::try_from(value).ok())
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        {
            let numbers = heap.slice_mut(ti.tGroupNumber)?;
            *numbers
                .get_mut(group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = itg as AT_NUMB;
            *numbers
                .get_mut(inverse_base + group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = itg as AT_NUMB;
        }
        j = j.wrapping_add(T_GROUP_HDR_LEN as i32);
        len_tg = len_tg.wrapping_sub(T_GROUP_HDR_LEN as i32 - 1);
        {
            let group = heap
                .slice_mut(ti.t_group)?
                .get_mut(group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            group.nNumEndpoints = len_tg as AT_NUMB;
            group.nFirstEndpointAtNoPos = i as AT_NUMB;
        }
        while len_tg > 0 {
            let source_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let k = i32::from(*tautomer.get(source_index).ok_or(SourceHeapError::PointerOutOfBounds)?) - 1;
            let endpoint_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(ti.nEndpointAtomNumber)?
                .get_mut(endpoint_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = k as AT_NUMB;
            if k < 0 || k > num_atoms {
                return Ok(RI_ERR_PROGR);
            }
            let k_index = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let group_number = itg.wrapping_add(1) as AT_NUMB;
            if !at.is_null() {
                heap.slice_mut(at)?
                    .get_mut(k_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .endpoint = group_number;
            }
            if !endpoint.is_null() {
                *heap
                    .slice_mut(endpoint)?
                    .get_mut(k_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = group_number;
            }
            len_tg -= 1;
            j += 1;
            i += 1;
        }
        itg += 1;
    }
    Ok(if i != ti.nNumEndpoints { RI_ERR_PROGR } else { 0 })
}

#[allow(non_snake_case)]
pub(crate) fn IncrZeroBonds(atoms: &mut [inp_ATOM], num_atoms: i32, component: i32) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4172 IncrZeroBonds
    // INCHI✔️❌: void IncrZeroBonds( inp_ATOM *at, int num_at, int iComponent )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j;
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[i].component = iComponent;
    // INCHI✔️❌:         for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!at[i].bond_type[j])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[i].bond_type[j] = BOND_TYPE_SINGLE;
    // INCHI✔️❌:                 at[i].chem_bonds_valence += BOND_TYPE_SINGLE;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: IncrZeroBonds
    // BEGIN ACTIVE INCHI MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/incomdef.h:67
    // INCHI✔️❌: #define BOND_TYPE_SINGLE    1
    // END ACTIVE INCHI MACRO

    let count = usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    for index in 0..count {
        let atom = atoms.get_mut(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.component = component as u16;
        let valence =
            usize::try_from(i32::from(atom.valence).max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        for bond in 0..valence {
            let bond_type = atom
                .bond_type
                .get_mut(bond)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if *bond_type == 0 {
                *bond_type = BOND_TYPE_SINGLE as u8;
                atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(BOND_TYPE_SINGLE as i8);
            }
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn ClearEndpts(atoms: &mut [inp_ATOM], num_atoms: i32) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4191 ClearEndpts
    // INCHI✔️❌: void ClearEndpts( inp_ATOM *at, int num_at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[i].endpoint = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ClearEndpts

    let count = usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    for index in 0..count {
        atoms
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint = 0;
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn ConnectDisconnectedH(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
    num_deleted_h: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5480 ConnectDisconnectedH
    // INCHI✔️❌: int ConnectDisconnectedH( inp_ATOM *at, int num_atoms, int num_deleted_H )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, n, m, num_H;
    // INCHI✔️❌:     int tot_atoms = num_atoms + num_deleted_H;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = num_atoms; i < tot_atoms; i = j)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         k = at[i].neighbor[0]; /* a[k] is the atom connected to the explicit hydrogen at[i] */
    // INCHI✔️❌:
    // INCHI✔️❌:         for (j = i; j < tot_atoms && at[j].neighbor[0] == k; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         num_H = j - i; /* number of explicit H for at[k] */
    // INCHI✔️❌:         if (num_H > at[k].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (num_H + at[k].valence > MAXVAL)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_SYNTAX;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* insert links to explicit H before all other links in the connection list */
    // INCHI✔️❌:         n = at[k].valence;
    // INCHI✔️❌:         memmove(at[k].neighbor + num_H, at[k].neighbor, sizeof(at[k].neighbor[0]) * n);
    // INCHI✔️❌:         memmove(at[k].bond_stereo + num_H, at[k].bond_stereo, sizeof(at[k].bond_stereo[0]) * n);
    // INCHI✔️❌:         memmove(at[k].bond_type + num_H, at[k].bond_type, sizeof(at[k].bond_type[0]) * n);
    // INCHI✔️❌:         for (n = 0; n < num_H; n++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[k].neighbor[n] = i + n;
    // INCHI✔️❌:             at[k].bond_stereo[n] = 0;
    // INCHI✔️❌:             at[k].bond_type[n] = BOND_TYPE_SINGLE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[k].sb_ord[m] += num_H;
    // INCHI✔️❌:             if (at[k].sn_ord[m] < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (n = i; n < j; n++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[n].orig_at_number == at[k].sn_orig_at_num[m])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at[k].sn_ord[m] = n - i;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (n == j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return RI_ERR_PROGR;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[k].sn_ord[m] += num_H;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         at[k].valence += num_H;
    // INCHI✔️❌:         at[k].chem_bonds_valence += num_H;
    // INCHI✔️❌:         at[k].num_H -= num_H; /* cannot be negative */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*memset( at[k].num_iso_H, 0, sizeof(at[0].num_iso_H) );*/ /* attached H must carry all isotopic shifts */
    // INCHI✔️❌:         for (n = i; n < j; n++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[n].chem_bonds_valence = BOND_TYPE_SINGLE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* isotopic H */
    // INCHI✔️❌:         for (m = j - 1; i <= m && at[m].iso_atw_diff > 0; m--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[m].iso_atw_diff > NUM_H_ISOTOPES)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return RI_ERR_PROGR;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (0 >= at[k].num_iso_H[(int) at[m].iso_atw_diff - 1] --)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return RI_ERR_PROGR;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* subtract isotopic H */
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (m = 0; m < NUM_H_ISOTOPES; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[i].num_H -= at[i].num_iso_H[m];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (0 > at[i].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return tot_atoms;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ConnectDisconnectedH
    // BEGIN ACTIVE INCHI MACROS
    // INCHI✔️❌: #define MAXVAL                20 /* max number of bonds per atom */
    // INCHI✔️❌: #define NUM_H_ISOTOPES        3  /* number of hydrogen isotopes: protium, deuterium, tritium */
    // INCHI✔️❌: #define BOND_TYPE_SINGLE       1
    // INCHI✔️❌: #define MAX_NUM_STEREO_BONDS   3
    // INCHI✔️❌: #define RI_ERR_SYNTAX          (-2)
    // INCHI✔️❌: #define RI_ERR_PROGR           (-3)
    // END ACTIVE INCHI MACROS

    let total_atoms = num_atoms
        .checked_add(num_deleted_h)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let atom_index = |value: i32| -> Result<usize, SourceHeapError> {
        usize::try_from(value).map_err(|_| SourceHeapError::PointerOutOfBounds)
    };

    let mut i = num_atoms;
    while i < total_atoms {
        let explicit_index = atom_index(i)?;
        let k = i32::from(
            atoms
                .get(explicit_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .neighbor[0],
        );
        let base_index = atom_index(k)?;

        let mut j = i;
        while j < total_atoms {
            let current = atoms.get(atom_index(j)?).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(current.neighbor[0]) != k {
                break;
            }
            j += 1;
        }

        let num_h = j - i;
        let base = atoms.get(base_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if num_h > i32::from(base.num_H) {
            return Ok(RI_ERR_PROGR);
        }
        if num_h + i32::from(base.valence) > MAXVAL as i32 {
            return Ok(RI_ERR_SYNTAX);
        }

        let old_valence = usize::try_from(i32::from(base.valence)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let insert_count = usize::try_from(num_h).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let base = atoms.get_mut(base_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        base.neighbor.copy_within(0..old_valence, insert_count);
        base.bond_stereo.copy_within(0..old_valence, insert_count);
        base.bond_type.copy_within(0..old_valence, insert_count);
        for n in 0..insert_count {
            base.neighbor[n] = (i + n as i32) as u16;
            base.bond_stereo[n] = 0;
            base.bond_type[n] = BOND_TYPE_SINGLE as u8;
        }

        for m in 0..MAX_NUM_STEREO_BONDS as usize {
            if atoms[base_index].sb_parity[m] == 0 {
                break;
            }
            atoms[base_index].sb_ord[m] = atoms[base_index].sb_ord[m].wrapping_add(num_h as i8);
            if atoms[base_index].sn_ord[m] < 0 {
                let target = atoms[base_index].sn_orig_at_num[m];
                let mut n = i;
                while n < j {
                    if atoms
                        .get(atom_index(n)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .orig_at_number
                        == target
                    {
                        atoms[base_index].sn_ord[m] = (n - i) as i8;
                        break;
                    }
                    n += 1;
                }
                if n == j {
                    return Ok(RI_ERR_PROGR);
                }
            } else {
                atoms[base_index].sn_ord[m] = atoms[base_index].sn_ord[m].wrapping_add(num_h as i8);
            }
        }

        let base = atoms.get_mut(base_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        base.valence = base.valence.wrapping_add(num_h as i8);
        base.chem_bonds_valence = base.chem_bonds_valence.wrapping_add(num_h as i8);
        base.num_H = base.num_H.wrapping_sub(num_h as i8);

        let mut n = i;
        while n < j {
            atoms
                .get_mut(atom_index(n)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .chem_bonds_valence = BOND_TYPE_SINGLE as i8;
            n += 1;
        }

        let mut m = j - 1;
        while i <= m {
            let isotope_difference = atoms
                .get(atom_index(m)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iso_atw_diff;
            if isotope_difference <= 0 {
                break;
            }
            if i32::from(isotope_difference) > NUM_H_ISOTOPES as i32 {
                return Ok(RI_ERR_PROGR);
            }
            let isotope_index =
                usize::try_from(i32::from(isotope_difference) - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let count = &mut atoms
                .get_mut(base_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .num_iso_H[isotope_index];
            let old_count = *count;
            *count = count.wrapping_sub(1);
            if old_count <= 0 {
                return Ok(RI_ERR_PROGR);
            }
            m -= 1;
        }
        i = j;
    }

    let mut i = 0;
    while i < num_atoms {
        let atom = atoms
            .get_mut(atom_index(i)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for m in 0..NUM_H_ISOTOPES as usize {
            atom.num_H = atom.num_H.wrapping_sub(atom.num_iso_H[m]);
        }
        if atom.num_H < 0 {
            return Ok(RI_ERR_PROGR);
        }
        i += 1;
    }

    Ok(total_atoms)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn DisconnectedConnectedH(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
    num_deleted_h: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5594 DisconnectedConnectedH
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int DisconnectedConnectedH( inp_ATOM *at, int num_atoms, int num_deleted_H )
    {
        int i, j, k, n, m, num_H, num_iso_H;
        int tot_atoms = num_atoms + num_deleted_H;
        S_CHAR ctmp = '0';

        /* add implicit isotopic H to total implicit H */
        for (i = 0; i < num_atoms; i++)
        {
            for (m = 0; m < NUM_H_ISOTOPES; m++)
            {
                at[i].num_H += at[i].num_iso_H[m];
            }
        }

        for (i = num_atoms; i < tot_atoms; i = j)
        {
            k = at[i].neighbor[0]; /* a[k] is the atom connected to the explicit hydrogen at[i] */

            for (j = i; j < tot_atoms && at[j].neighbor[0] == k; j++)
            {
                at[j].chem_bonds_valence = 0;
            }
            num_H = j - i; /* number of explicit H for at[k] */

            /* verify correct number of explicit H */
            for (n = 0; n < at[k].valence && at[k].neighbor[n] >= num_atoms; n++)
            {
                ;
            }
            if (n != num_H)
            {
                return RI_ERR_PROGR;
            }

            /* remove bonds to explicit H located in front of all other bonds in the connection list */
            at[k].valence -= num_H; /* djb-rwth: use of cast operators avoided */
            n = at[k].valence; /* new number of bonds */ 
            at[k].chem_bonds_valence -= num_H; /* new no-H valence */
            if (n)
            {
                memmove(at[k].neighbor, at[k].neighbor + num_H, sizeof(at[k].neighbor[0]) * n);
                memmove(at[k].bond_stereo, at[k].bond_stereo + num_H, sizeof(at[k].bond_stereo[0]) * n);
                memmove(at[k].bond_type, at[k].bond_type + num_H, sizeof(at[k].bond_type[0]) * n);
            }
            /* clear the 'tails' */
            memset( at[k].neighbor + n, 0, sizeof( at[k].neighbor[0] ) * num_H ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( at[k].bond_stereo + n, 0, sizeof( at[k].bond_stereo[0] ) * num_H ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( at[k].bond_type + n, 0, sizeof( at[k].bond_type[0] ) * num_H ); /* djb-rwth: memset_s C11/Annex K variant? */

            for (m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m++)
            {
                at[k].sb_ord[m] -= num_H;
                if (0 <= at[k].sn_ord[m] && at[k].sn_ord[m] < num_H)
                {
                    at[k].sn_ord[m] = -1; /* disconnected explicit H */
                }
            }
            /* add explicit isotopic H (already included in num_H) */
            for (num_iso_H = 0, m = j - 1; i <= m && at[m].iso_atw_diff > 0; m--)
            {
                if (at[m].iso_atw_diff > NUM_H_ISOTOPES)
                {
                    return RI_ERR_PROGR;
                }
                at[k].num_iso_H[(int) at[m].iso_atw_diff - 1] ++;
            }
            at[k].num_H += num_H; /* add all explicit H including isotopic */
        }

        return tot_atoms;
    }

    */
    // END INCHI C FUNCTION: DisconnectedConnectedH
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: DisconnectedConnectedH
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // INCHI✔️✔️: #define MAX_NUM_STEREO_BONDS 3
    // INCHI✔️✔️: #define RI_ERR_PROGR (-3)
    // INCHI✔️✔️: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: DisconnectedConnectedH

    let total_atoms = num_atoms
        .checked_add(num_deleted_h)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let atom_index = |value: i32| -> Result<usize, SourceHeapError> {
        usize::try_from(value).map_err(|_| SourceHeapError::PointerOutOfBounds)
    };

    for atom_number in 0..num_atoms {
        let atom = atoms
            .get_mut(atom_index(atom_number)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for isotope in 0..NUM_H_ISOTOPES as usize {
            atom.num_H = atom.num_H.wrapping_add(atom.num_iso_H[isotope]);
        }
    }

    let mut i = num_atoms;
    while i < total_atoms {
        let explicit_index = atom_index(i)?;
        let attached_atom = i32::from(
            atoms
                .get(explicit_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .neighbor[0],
        );
        let attached_index = atom_index(attached_atom)?;

        let mut j = i;
        while j < total_atoms {
            let hydrogen = atoms
                .get(atom_index(j)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(hydrogen.neighbor[0]) != attached_atom {
                break;
            }
            atoms
                .get_mut(atom_index(j)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .chem_bonds_valence = 0;
            j = j.wrapping_add(1);
        }
        let num_h = j.wrapping_sub(i);

        let attached = atoms
            .get(attached_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let old_valence = i32::from(attached.valence);
        let mut leading_h = 0_i32;
        while leading_h < old_valence {
            let order = usize::try_from(leading_h)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let neighbor = attached
                .neighbor
                .get(order)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(*neighbor) < num_atoms {
                break;
            }
            leading_h = leading_h.wrapping_add(1);
        }
        if leading_h != num_h {
            return Ok(RI_ERR_PROGR);
        }

        let new_valence = old_valence.wrapping_sub(num_h);
        let old_count = usize::try_from(old_valence)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let remove_count = usize::try_from(num_h)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let new_count = usize::try_from(new_valence)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if old_count > atoms[attached_index].neighbor.len()
            || remove_count > old_count
            || new_count.wrapping_add(remove_count) > atoms[attached_index].neighbor.len()
        {
            return Err(SourceHeapError::PointerOutOfBounds);
        }

        {
            let attached = atoms
                .get_mut(attached_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            attached.valence = attached.valence.wrapping_sub(num_h as i8);
            attached.chem_bonds_valence =
                attached.chem_bonds_valence.wrapping_sub(num_h as i8);
            if new_count != 0 {
                attached.neighbor.copy_within(remove_count..old_count, 0);
                attached
                    .bond_stereo
                    .copy_within(remove_count..old_count, 0);
                attached.bond_type.copy_within(remove_count..old_count, 0);
            }
            attached.neighbor[new_count..new_count + remove_count].fill(0);
            attached.bond_stereo[new_count..new_count + remove_count].fill(0);
            attached.bond_type[new_count..new_count + remove_count].fill(0);

            for stereo_index in 0..MAX_NUM_STEREO_BONDS as usize {
                if attached.sb_parity[stereo_index] == 0 {
                    break;
                }
                attached.sb_ord[stereo_index] =
                    attached.sb_ord[stereo_index].wrapping_sub(num_h as i8);
                if attached.sn_ord[stereo_index] >= 0
                    && i32::from(attached.sn_ord[stereo_index]) < num_h
                {
                    attached.sn_ord[stereo_index] = -1;
                }
            }
        }

        let mut m = j.wrapping_sub(1);
        while i <= m {
            let isotope_difference = atoms
                .get(atom_index(m)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iso_atw_diff;
            if isotope_difference <= 0 {
                break;
            }
            if i32::from(isotope_difference) > NUM_H_ISOTOPES as i32 {
                return Ok(RI_ERR_PROGR);
            }
            let isotope_index = usize::try_from(i32::from(isotope_difference).wrapping_sub(1))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let attached = atoms
                .get_mut(attached_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            attached.num_iso_H[isotope_index] =
                attached.num_iso_H[isotope_index].wrapping_add(1);
            m = m.wrapping_sub(1);
        }
        let attached = atoms
            .get_mut(attached_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        attached.num_H = attached.num_H.wrapping_add(num_h as i8);
        i = j;
    }

    Ok(total_atoms)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeInChIOutOfStrFromINChI2(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    input_parameters: Option<&INPUT_PARMS>,
    input_structure_data: Option<&mut STRUCT_DATA>,
    mut structure: Option<&mut StrFromINChI>,
    component: i32,
    _atom_number_offset: i32,
    input_number: i64,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5669 MakeInChIOutOfStrFromINChI2
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: MakeInChIOutOfStrFromINChI2
    // INCHI✔️❌: int MakeInChIOutOfStrFromINChI2( INCHI_CLOCK *ic,
    // INCHI✔️❌:                                  CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                  ICHICONST INPUT_PARMS *ip_inp,
    // INCHI✔️❌:                                  STRUCT_DATA *sd_inp,
    // INCHI✔️❌:                                  StrFromINChI *pStruct,
    // INCHI✔️❌:                                  int iComponent,
    // INCHI✔️❌:                                  int iAtNoOffset,
    // INCHI✔️❌:                                  long num_inp )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szTitle[MAX_SDF_HEADER + MAX_SDF_VALUE + 256];
    // INCHI✔️❌:
    // INCHI✔️❌:     int   len, ret;
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     PINChI2     *pINChI[INCHI_NUM];
    // INCHI✔️❌:     PINChI_Aux2 *pINChI_Aux[INCHI_NUM];
    // INCHI✔️❌:     */
    // INCHI✔️❌:     INPUT_PARMS local_ip;
    // INCHI✔️❌:     STRUCT_DATA local_sd;
    // INCHI✔️❌:     INPUT_PARMS *ip = &local_ip;
    // INCHI✔️❌:     STRUCT_DATA *sd = &local_sd;
    // INCHI✔️❌:
    // INCHI✔️❌:     ORIG_ATOM_DATA OrigAtData; /* 0=> disconnected, 1=> original */
    // INCHI✔️❌:     ORIG_ATOM_DATA *orig_inp_data = &OrigAtData;
    // INCHI✔️❌:     ORIG_ATOM_DATA PrepAtData[2]; /* 0=> disconnected, 1=> original */
    // INCHI✔️❌:     ORIG_ATOM_DATA *prep_inp_data = PrepAtData;
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_IOS_STRING temp_string_container;
    // INCHI✔️❌:     INCHI_IOS_STRING *strbuf = &temp_string_container;
    // INCHI✔️❌:     memset( strbuf, 0, sizeof( *strbuf ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 >= inchi_strbuf_init(strbuf, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_ALLOC;
    // INCHI✔️❌:         goto early_exit_error; /* djb-rwth: avoiding garbage value */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #70552 */
    // INCHI✔️❌:     if (!ip_inp || !sd_inp || !pStruct)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_ALLOC;
    // INCHI✔️❌:         goto early_exit_error; /* djb-rwth: avoiding garbage value */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *ip = *ip_inp;
    // INCHI✔️❌:     ip->bDisplay = 0;
    // INCHI✔️❌:     ip->bDisplayCompositeResults = 0;
    // INCHI✔️❌:     ip->bDisplayEachComponentINChI = 0;
    // INCHI✔️❌:     ip->bDisplayIfRestoreWarnings = 0;
    // INCHI✔️❌:     ip->bINChIOutputOptions = INCHI_OUT_NO_AUX_INFO;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( pStruct->bMobileH ) {
    // INCHI✔️❌:         ip->nMode &= ~REQ_MODE_BASIC;
    // INCHI✔️❌:         ip->nMode |= REQ_MODE_TAUT;
    // INCHI✔️❌:     } else {
    // INCHI✔️❌:         ip->nMode |= (REQ_MODE_TAUT | REQ_MODE_BASIC);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     sd->fPtrStart = -1;
    // INCHI✔️❌:     sd->fPtrEnd = -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( ip->nMode & REQ_MODE_STEREO ) {
    // INCHI✔️❌:         if ( ip->nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO) ) {
    // INCHI✔️❌:             sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL;
    // INCHI✔️❌:         } else {
    // INCHI✔️❌:             sd->bChiralFlag |= FLAG_INP_AT_CHIRAL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( orig_inp_data, 0, sizeof( *orig_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( prep_inp_data, 0, 2 * sizeof( *prep_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( pStruct->RevInChI.pINChI, 0, sizeof( pStruct->RevInChI.pINChI ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( pStruct->RevInChI.pINChI_Aux, 0, sizeof( pStruct->RevInChI.pINChI_Aux ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( szTitle, 0, sizeof( szTitle ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     len = sizeof( orig_inp_data->at[0] )*( (long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H ); /* djb-rwth: cast operators added */
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_inp_data->at = (inp_ATOM *) inchi_malloc( len );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_inp_data->at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*memcpy( orig_inp_data->at, pStruct->at2, len );*/
    // INCHI✔️❌:         /*ret = ConnectDisconnectedH( orig_inp_data->at, pStruct->num_atoms, pStruct->num_deleted_H );*/
    // INCHI✔️❌:
    // INCHI✔️❌:         CopySt2At( pStruct->at2, pStruct->st, pStruct->num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = ConnectDisconnectedH( pStruct->at2, pStruct->num_atoms, pStruct->num_deleted_H );
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         orig_inp_data->num_inp_atoms = ret;
    // INCHI✔️❌:         /* connections changed => reconcile parities even if they were reconciled before */
    // INCHI✔️❌:         /* remove t-group markings and increment zero-order bonds,
    // INCHI✔️❌:            otherwise MakeInChIOutOfStrFromINChI2() woild fail */
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         IncrZeroBondsAndClearEndpts(pStruct->at2, pStruct->num_atoms, iComponent+1);
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         IncrZeroBonds( pStruct->at2, pStruct->num_atoms, iComponent + 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* CopySt2At() moved to the position before ConnectDisconnectedH() because
    // INCHI✔️❌:            in case stereo exists only in Mobile-H layer and the processd here
    // INCHI✔️❌:            component is restored in Fixed-H layer the parities needed by
    // INCHI✔️❌:            ConnectDisconnectedH() must be there before calling
    // INCHI✔️❌:            ConnectDisconnectedH()
    // INCHI✔️❌:         */
    // INCHI✔️❌:         /*CopySt2At( pStruct->at2, pStruct->st, pStruct->num_atoms );*/
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = ReconcileAllCmlBondParities( pStruct->at2, orig_inp_data->num_inp_atoms, 0 );
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if USE_BCF
    // INCHI✔️❌:         memcpy_s( orig_inp_data->at, (long long)len, pStruct->at2, len ); /* djb-rwth: function replaced with its safe C11 variant */
    // INCHI✔️❌: #else
    // INCHI✔️❌:         memcpy(orig_inp_data->at, pStruct->at2, len);
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         ClearEndpts( orig_inp_data->at, pStruct->num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (FixUnkn0DStereoBonds( orig_inp_data->at, pStruct->num_atoms ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = ReconcileAllCmlBondParities( pStruct->at2, orig_inp_data->num_inp_atoms, 0 );
    // INCHI✔️❌:             if (ret < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* keep endpoint[] markings in at2[] for subsequent add/remove protons */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_ALLOC;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( sd->num_components, 0, sizeof( sd->num_components ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( sd->num_taut, 0, sizeof( sd->num_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( sd->num_non_taut, 0, sizeof( sd->num_non_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( sd->bTautFlagsDone, 0, sizeof( sd->bTautFlagsDone ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( sd->bTautFlags, 0, sizeof( sd->bTautFlags ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = ProcessOneStructure( ic, pCG, sd, ip, szTitle,
    // INCHI✔️❌:                                pStruct->RevInChI.pINChI,
    // INCHI✔️❌:                                pStruct->RevInChI.pINChI_Aux,
    // INCHI✔️❌:                                NULL /*inp_file*/,
    // INCHI✔️❌:                                NULL /*log_file*/,
    // INCHI✔️❌:                                NULL /*out_file*/,
    // INCHI✔️❌:                                NULL /*prb_file*/,
    // INCHI✔️❌:                                orig_inp_data, prep_inp_data,
    // INCHI✔️❌:                                num_inp, strbuf,
    // INCHI✔️❌:                                0 /* save_opt_bits */ );
    // INCHI✔️❌:
    // INCHI✔️❌:     memcpy(pStruct->RevInChI.num_components, sd->num_components, sizeof(pStruct->RevInChI.num_components));
    // INCHI✔️❌:     memcpy(sd_inp->pStrErrStruct, sd->pStrErrStruct, sizeof(sd_inp->pStrErrStruct));
    // INCHI✔️❌:     pStruct->RevInChI.nRetVal = ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* translate returned value */
    // INCHI✔️❌:     if (ret == _IS_ERROR || ret == _IS_FATAL || ret == _IS_UNKNOWN)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_PROGR;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (ret == _IS_OKAY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (ret == _IS_WARNING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_PROGR;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* save total charge from Mobile-H layer */
    // INCHI✔️❌:
    // INCHI✔️❌:     pStruct->nChargeRevrs = 0;
    // INCHI✔️❌:     if (ret >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bRevInchiComponentExists( pStruct, INCHI_REC, TAUT_YES, 0 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pStruct->nChargeRevrs = pStruct->RevInChI.pINChI[INCHI_REC][0][TAUT_YES]->nTotalCharge;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (bRevInchiComponentExists( pStruct, INCHI_BAS, TAUT_YES, 0 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pStruct->nChargeRevrs = pStruct->RevInChI.pINChI[INCHI_BAS][0][TAUT_YES]->nTotalCharge;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* free structure data */
    // INCHI✔️❌:     FreeOrigAtData( orig_inp_data );
    // INCHI✔️❌:     FreeOrigAtData( prep_inp_data );
    // INCHI✔️❌:     FreeOrigAtData( prep_inp_data + 1 );
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌: #if ( FIX_GAF_2019_2==1 )
    // INCHI✔️❌:     if (orig_inp_data->at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         FreeInpAtom(&orig_inp_data->at);
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: early_exit_error:
    // INCHI✔️❌:     inchi_strbuf_close( strbuf );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: MakeInChIOutOfStrFromINChI2
    // END INCHI C FUNCTION: MakeInChIOutOfStrFromINChI2
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MakeInChIOutOfStrFromINChI2
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; USE_BCF=0; FIX_GAF_2019_2=1.
    // INCHI✔️❌: SourceHeap model-storage slots represent C stack structs passed by pointer.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MakeInChIOutOfStrFromINChI2

    let mut string_buffer = INCHI_IOS_STRING::default();
    if inchi_strbuf_init(
        heap,
        &mut string_buffer,
        INCHI_STRBUF_INITIAL_SIZE as i32,
        INCHI_STRBUF_SIZE_INCREMENT as i32,
    )? <= 0
    {
        inchi_strbuf_close(heap, Some(&mut string_buffer))?;
        return Ok(RI_ERR_ALLOC);
    }

    let (Some(input_parameters), Some(input_structure_data), Some(structure)) =
        (input_parameters, input_structure_data, structure.as_deref_mut())
    else {
        inchi_strbuf_close(heap, Some(&mut string_buffer))?;
        return Ok(RI_ERR_ALLOC);
    };

    let mut local_parameters = input_parameters.clone();
    local_parameters.bDisplay = 0;
    local_parameters.bDisplayCompositeResults = 0;
    local_parameters.bDisplayEachComponentINChI = 0;
    local_parameters.bDisplayIfRestoreWarnings = 0;
    local_parameters.bINChIOutputOptions = INCHI_OUT_NO_AUX_INFO as i32;
    let mut local_structure_data = STRUCT_DATA {
        fPtrStart: -1,
        fPtrEnd: -1,
        ..STRUCT_DATA::default()
    };

    structure.RevInChI.pINChI = [SourceMutPointer::null(); INCHI_NUM as usize];
    structure.RevInChI.pINChI_Aux = [SourceMutPointer::null(); INCHI_NUM as usize];
    let mut title = vec![0_i8; 1024];

    let total_atoms_i32 = structure.num_atoms.wrapping_add(structure.num_deleted_H);
    let total_atoms = match usize::try_from(total_atoms_i32) {
        Ok(value) => value,
        Err(_) => {
            inchi_strbuf_close(heap, Some(&mut string_buffer))?;
            return Ok(RI_ERR_ALLOC);
        }
    };
    let original_atoms = match heap.allocate(vec![inp_ATOM::default(); total_atoms]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_strbuf_close(heap, Some(&mut string_buffer))?;
            return Ok(RI_ERR_ALLOC);
        }
        Err(error) => return Err(error),
    };
    let mut original_atom_data = ORIG_ATOM_DATA {
        at: original_atoms,
        ..ORIG_ATOM_DATA::default()
    };
    let mut prepared_atom_data = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];

    let operation_result = (|| -> Result<(i32, bool), SourceHeapError> {
        let stereo_count =
            usize::try_from(structure.num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let connected = if structure.at2.allocation_identity() != structure.st.allocation_identity() {
            heap.with_slice_mut_and_heap(structure.at2, |atoms, heap| {
                let atoms = atoms
                    .get_mut(..total_atoms)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let stereo = if structure.st.is_null() {
                    None
                } else {
                    Some(
                        heap.slice(structure.st.as_const())?
                            .get(..stereo_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    )
                };
                CopySt2At(atoms, stereo, structure.num_atoms)?;
                let connected = ConnectDisconnectedH(atoms, structure.num_atoms, structure.num_deleted_H)?;
                if connected >= 0 {
                    IncrZeroBonds(atoms, structure.num_atoms, component.wrapping_add(1))?;
                }
                Ok(connected)
            })?
        } else {
            // Preserve the checked error behavior of malformed modeled alias
            // pointers; Official production uses distinct typed allocations.
            let mut atoms = heap
                .slice(structure.at2.as_const())?
                .get(..total_atoms)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let stereo = if structure.st.is_null() {
                None
            } else {
                Some(
                    heap.slice(structure.st.as_const())?
                        .get(..stereo_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec(),
                )
            };
            CopySt2At(&mut atoms, stereo.as_deref(), structure.num_atoms)?;
            let connected = ConnectDisconnectedH(&mut atoms, structure.num_atoms, structure.num_deleted_H)?;
            if connected >= 0 {
                IncrZeroBonds(&mut atoms, structure.num_atoms, component.wrapping_add(1))?;
            }
            heap.slice_mut(structure.at2)?[..total_atoms].clone_from_slice(&atoms);
            connected
        };
        if connected < 0 {
            return Ok((connected, false));
        }
        original_atom_data.num_inp_atoms = connected;

        let reconciled = ReconcileAllCmlBondParities(heap, structure.at2, original_atom_data.num_inp_atoms, 0)?;
        if reconciled < 0 {
            return Ok((reconciled, false));
        }

        copy_inp_atom_prefix(heap, original_atom_data.at, structure.at2, total_atoms)?;
        let original_atom_count =
            usize::try_from(structure.num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let copied_atoms = heap
            .slice_mut(original_atom_data.at)?
            .get_mut(..original_atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        ClearEndpts(copied_atoms, structure.num_atoms)?;
        if FixUnkn0DStereoBonds(heap, original_atom_data.at, structure.num_atoms)? != 0 {
            let reconciled = ReconcileAllCmlBondParities(heap, structure.at2, original_atom_data.num_inp_atoms, 0)?;
            if reconciled < 0 {
                return Ok((reconciled, false));
            }
        }

        let original_slot = heap.allocate_model_storage(vec![original_atom_data.clone()])?;
        let prepared_slot = heap.allocate_model_storage(prepared_atom_data.to_vec())?;
        let process_result = ProcessOneStructure(
            heap,
            clock,
            canonical_globals,
            &mut local_structure_data,
            &mut local_parameters,
            Some(&mut title),
            &mut structure.RevInChI.pINChI,
            &mut structure.RevInChI.pINChI_Aux,
            None,
            None,
            None,
            None,
            original_slot,
            prepared_slot,
            input_number,
            Some(&mut string_buffer),
            0,
            SourceMutPointer::null(),
            clock_result,
        );
        original_atom_data = heap.slice(original_slot.as_const())?[0].clone();
        prepared_atom_data.clone_from_slice(&heap.slice(prepared_slot.as_const())?[..2]);
        heap.free(original_slot)?;
        heap.free(prepared_slot)?;
        let process_result = process_result?;

        structure.RevInChI.num_components = local_structure_data.num_components;
        input_structure_data.pStrErrStruct = local_structure_data.pStrErrStruct;
        structure.RevInChI.nRetVal = process_result;
        let translated = if process_result == _IS_ERROR as i32
            || process_result == _IS_FATAL as i32
            || process_result == _IS_UNKNOWN as i32
        {
            RI_ERR_PROGR
        } else if process_result == _IS_OKAY as i32 {
            0
        } else if process_result == _IS_WARNING as i32 {
            1
        } else {
            RI_ERR_PROGR
        };

        structure.nChargeRevrs = 0;
        if translated >= 0 {
            let (representation, exists) = if bRevInchiComponentExists(
                heap,
                Some(structure),
                crate::source_types::INCHI_REC as i32,
                TAUT_YES as i32,
                0,
            )? != 0
            {
                (crate::source_types::INCHI_REC as usize, true)
            } else if bRevInchiComponentExists(heap, Some(structure), INCHI_BAS as i32, TAUT_YES as i32, 0)? != 0 {
                (INCHI_BAS as usize, true)
            } else {
                (0, false)
            };
            if exists {
                let components = heap.slice(structure.RevInChI.pINChI[representation].as_const())?;
                let inchi = heap
                    .slice(components[0][TAUT_YES as usize].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                structure.nChargeRevrs = inchi.nTotalCharge;
            }
        }
        Ok((translated, true))
    })();

    let (result, reached_normal_cleanup) = match operation_result {
        Ok(value) => value,
        Err(error) => {
            inchi_strbuf_close(heap, Some(&mut string_buffer))?;
            return Err(error);
        }
    };
    if !reached_normal_cleanup {
        if !original_atom_data.at.is_null() {
            inchi_free(heap, original_atom_data.at)?;
            original_atom_data.at = SourceMutPointer::null();
        }
    } else {
        FreeOrigAtData(heap, Some(&mut original_atom_data))?;
        FreeOrigAtData(heap, Some(&mut prepared_atom_data[0]))?;
        FreeOrigAtData(heap, Some(&mut prepared_atom_data[1]))?;
    }
    inchi_strbuf_close(heap, Some(&mut string_buffer))?;
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn get_sp_element_type(n_periodic_number: i32, row: &mut i32) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1414 get_sp_element_type
    // INCHI✔️✔️: int get_sp_element_type( int nPeriodicNumber, int *nRow )
    // INCHI✔️✔️: /*
    // INCHI✔️✔️:                                                 num  el
    // INCHI✔️✔️:                                                 el   neg
    // INCHI✔️✔️:    1 => H                          ATYPE_H   1    1  21
    // INCHI✔️✔️:    2 => Li, Na, K,  Rb, Cs, Fr     ATYPE_Na  2    1  10 09 08 08 07
    // INCHI✔️✔️:    3 => Be, Mg, Ca, Sr, Ba, Ra     ATYPE_Mg  3    2  15 12 10 10 09
    // INCHI✔️✔️:    4 => B,  Al, Ga, In, Tl         ATYPE_B   4    3  20 15 18 17 18
    // INCHI✔️✔️:    5 => C,  Si, Ge, Sn, Pb         ATYPE_C   5    4  25 18 18 18 18
    // INCHI✔️✔️:    6 => N,  P,  As, Sb, Bi         ATYPE_N   6    5  30 21 20 19 19
    // INCHI✔️✔️:    7 => O,  S,  Se, Te, Po         ATYPE_O   7    6  35 25 24 21 20
    // INCHI✔️✔️:    8 => F,  Cl, Br, I,  At         ATYPE_Cl  8    7  40 30 28 25 22
    // INCHI✔️✔️:
    // INCHI✔️✔️: number of valence electrons = (type>1)? type-1: type
    // INCHI✔️✔️:
    // INCHI✔️✔️:   */
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int row = 0, type = 0;
    // INCHI✔️✔️:     if (nPeriodicNumber == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = 1; /* H: 1 */
    // INCHI✔️✔️:         row = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber == 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = 0; row = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 10)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Li: 2, Be: 3, B: 4, C: 5, N: 6, O: 7, F: 8, Ne: 9; later subtract 1 */
    // INCHI✔️✔️:         type = nPeriodicNumber - 1; row = 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 18)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 9; row = 2;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 20)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 17; row = 3;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 30)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = 0; row = 3;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 36)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 27; row = 3;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 38)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 35; row = 4;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 48)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = 0; row = 4;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 54)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 45; row = 4;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 56)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 53; row = 5;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 80)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = 0; row = 5;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 86)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 77; row = 5;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (nPeriodicNumber <= 88)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = nPeriodicNumber - 85; row = 6;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = 0; row = 6;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     *nRow = row;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return type == 9 ? 0 : type;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: get_sp_element_type

    let (element_type, element_row) = if n_periodic_number == 1 {
        (1, 0)
    } else if n_periodic_number == 2 {
        (0, 0)
    } else if n_periodic_number <= 10 {
        (n_periodic_number.wrapping_sub(1), 1)
    } else if n_periodic_number <= 18 {
        (n_periodic_number.wrapping_sub(9), 2)
    } else if n_periodic_number <= 20 {
        (n_periodic_number.wrapping_sub(17), 3)
    } else if n_periodic_number <= 30 {
        (0, 3)
    } else if n_periodic_number <= 36 {
        (n_periodic_number.wrapping_sub(27), 3)
    } else if n_periodic_number <= 38 {
        (n_periodic_number.wrapping_sub(35), 4)
    } else if n_periodic_number <= 48 {
        (0, 4)
    } else if n_periodic_number <= 54 {
        (n_periodic_number.wrapping_sub(45), 4)
    } else if n_periodic_number <= 56 {
        (n_periodic_number.wrapping_sub(53), 5)
    } else if n_periodic_number <= 80 {
        (0, 5)
    } else if n_periodic_number <= 86 {
        (n_periodic_number.wrapping_sub(77), 5)
    } else if n_periodic_number <= 88 {
        (n_periodic_number.wrapping_sub(85), 6)
    } else {
        (0, 6)
    };
    *row = element_row;
    if element_type == 9 { 0 } else { element_type }
}

#[allow(non_snake_case)]
pub(crate) fn ReallocTCGroups(
    heap: &mut SourceHeap,
    groups: &mut ALL_TC_GROUPS,
    add: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1502 ReallocTCGroups
    // INCHI✔️❌: int ReallocTCGroups( ALL_TC_GROUPS *pTCGroups, int nAdd )
    // INCHI✔️❌: {
    // INCHI✔️❌:     TC_GROUP *pTCGroup = (TC_GROUP *) inchi_malloc( sizeof( pTCGroup[0] )*( (long long)pTCGroups->max_tc_groups + nAdd ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pTCGroup)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pTCGroups->num_tc_groups)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memcpy(pTCGroup, pTCGroups->pTCG, sizeof(pTCGroup[0]) * pTCGroups->num_tc_groups);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         memset( pTCGroup + pTCGroups->max_tc_groups, 0, sizeof( pTCGroup[0] )*nAdd ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (pTCGroups->pTCG)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( pTCGroups->pTCG );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pTCGroups->pTCG = pTCGroup;
    // INCHI✔️❌:         pTCGroups->max_tc_groups += nAdd;
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return RI_ERR_ALLOC;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ReallocTCGroups

    let total = i64::from(groups.max_tc_groups) + i64::from(add);
    let Ok(total) = usize::try_from(total) else {
        return Ok(RI_ERR_ALLOC);
    };
    let mut replacement = Vec::new();
    if replacement.try_reserve_exact(total).is_err() {
        return Ok(RI_ERR_ALLOC);
    }
    replacement.resize_with(total, TC_GROUP::default);
    let copy_count = usize::try_from(groups.num_tc_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if copy_count != 0 {
        let source = heap
            .slice(groups.pTCG.as_const())?
            .get(..copy_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        replacement
            .get_mut(..copy_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(source);
    }
    let replacement_pointer = match heap.allocate(replacement) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(RI_ERR_ALLOC),
        Err(error) => return Err(error),
    };
    if !groups.pTCG.is_null() {
        inchi_free(heap, groups.pTCG)?;
    }
    groups.pTCG = replacement_pointer;
    groups.max_tc_groups = groups.max_tc_groups.wrapping_add(add);
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RegisterTCGroup(
    heap: &mut SourceHeap,
    groups: &mut ALL_TC_GROUPS,
    group_type: i32,
    group_order_number: i32,
    vertex_capacity: i32,
    vertex_flow: i32,
    edge_capacity: i32,
    edge_flow: i32,
    number_of_edges: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1527 RegisterTCGroup
    // INCHI✔️❌: int RegisterTCGroup( ALL_TC_GROUPS *pTCGroups, int nGroupType, int nGroupOrdNum,
    // INCHI✔️❌:                      int nVertexCap, int nVertexFlow, int nEdgeCap, int nEdgeFlow, int nNumEdges )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* search */
    // INCHI✔️❌:     for (i = 0; i < pTCGroups->num_tc_groups; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pTCGroups->pTCG[i].type == nGroupType &&
    // INCHI✔️❌:              pTCGroups->pTCG[i].ord_num == nGroupOrdNum)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i == pTCGroups->num_tc_groups)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* add one more group */
    // INCHI✔️❌:         if (pTCGroups->num_tc_groups == pTCGroups->max_tc_groups)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = ReallocTCGroups( pTCGroups, INC_NUM_TCGROUPS );
    // INCHI✔️❌:             if (ret)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = i + 1; /* added new group */
    // INCHI✔️❌:         pTCGroups->num_tc_groups++;
    // INCHI✔️❌:         pTCGroups->pTCG[i].type = nGroupType;
    // INCHI✔️❌:         pTCGroups->pTCG[i].ord_num = nGroupOrdNum;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     pTCGroups->pTCG[i].num_edges += nNumEdges;
    // INCHI✔️❌:
    // INCHI✔️❌:     pTCGroups->pTCG[i].st_cap += nVertexCap;
    // INCHI✔️❌:     pTCGroups->pTCG[i].st_flow += nVertexFlow;
    // INCHI✔️❌:
    // INCHI✔️❌:     pTCGroups->pTCG[i].edges_cap += nEdgeCap;
    // INCHI✔️❌:     pTCGroups->pTCG[i].edges_flow += nEdgeFlow;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: RegisterTCGroup

    let group_count = usize::try_from(groups.num_tc_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut index = 0_usize;
    if group_count != 0 {
        let values = heap
            .slice(groups.pTCG.as_const())?
            .get(..group_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        while index < group_count {
            if values[index].type_ == group_type && values[index].ord_num == group_order_number {
                break;
            }
            index += 1;
        }
    }
    let mut result = 0_i32;
    if index == group_count {
        if groups.num_tc_groups == groups.max_tc_groups {
            result = ReallocTCGroups(heap, groups, INC_NUM_TCGROUPS as i32)?;
            if result != 0 {
                return Ok(result);
            }
        }
        result = (index as i32).wrapping_add(1);
        groups.num_tc_groups = groups.num_tc_groups.wrapping_add(1);
        let group = heap
            .slice_mut(groups.pTCG)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        group.type_ = group_type;
        group.ord_num = group_order_number;
    }
    let group = heap
        .slice_mut(groups.pTCG)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    group.num_edges = group.num_edges.wrapping_add(number_of_edges);
    group.st_cap = group.st_cap.wrapping_add(vertex_capacity);
    group.st_flow = group.st_flow.wrapping_add(vertex_flow);
    group.edges_cap = group.edges_cap.wrapping_add(edge_capacity);
    group.edges_flow = group.edges_flow.wrapping_add(edge_flow);
    Ok(result)
}

// Exact active `cnList` node data from ichirvr1.c:102-504 under TARGET_API_LIB.
// ECF edge tuples preserve the source field order:
// (neighbor, capacity, bForbiddenEdge, flow).
pub(crate) type CnNodeData = (i32, i32, i32, [(i32, i32, i32, i32); 3]);
pub(crate) struct CnListData {
    pub(crate) bits: i32,
    nodes: &'static [CnNodeData],
}
pub(crate) const CN_LIST: [CnListData; 18] = [
    CnListData {
        bits: 650,
        nodes: &[
            (1, 3, 3, [(2, 3, 0, 3), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 3, 3, [(3, 1, 0, 0), (4, 2, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(5, 1, 0, 0), (4, 2, 0, 2), (0, 0, 0, 0)]),
            (192, 2, 2, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 2, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 1105,
        nodes: &[
            (1, 3, 3, [(2, 3, 0, 3), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 3, 3, [(3, 2, 0, 0), (4, 1, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(5, 1, 0, 1), (4, 1, 0, 1), (0, 0, 0, 0)]),
            (192, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 81,
        nodes: &[
            (1, 2, 2, [(2, 2, 0, 2), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(3, 1, 0, 0), (4, 1, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(5, 1, 0, 1), (4, 1, 0, 1), (0, 0, 0, 0)]),
            (192, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 138,
        nodes: &[
            (1, 2, 2, [(2, 2, 0, 2), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(3, 1, 0, 0), (4, 1, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(5, 1, 0, 0), (4, 1, 0, 1), (0, 0, 0, 0)]),
            (192, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 2, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 140,
        nodes: &[
            (1, 2, 1, [(2, 1, 0, 1), (3, 1, 0, 0), (0, 0, 0, 0)]),
            (16, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (272, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 266,
        nodes: &[
            (1, 2, 2, [(2, 2, 0, 2), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(3, 1, 0, 0), (4, 1, 0, 0), (0, 0, 0, 0)]),
            (16, 2, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (272, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 14,
        nodes: &[
            (1, 1, 0, [(2, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(3, 1, 0, 1), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(4, 1, 0, 0), (6, 1, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(5, 1, 0, 1), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (1040, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (1296, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 97,
        nodes: &[
            (1, 2, 2, [(2, 2, 0, 2), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(3, 1, 0, 0), (4, 1, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(5, 1, 0, 0), (4, 1, 0, 1), (0, 0, 0, 0)]),
            (192, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (272, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 49,
        nodes: &[
            (1, 1, 1, [(2, 1, 0, 1), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(3, 1, 0, 0), (5, 1, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(4, 1, 0, 1), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (272, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 105,
        nodes: &[
            (1, 2, 2, [(2, 2, 0, 2), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(3, 1, 0, 0), (4, 1, 0, 0), (0, 0, 0, 0)]),
            (192, 2, 2, [(5, 1, 0, 1), (4, 1, 0, 1), (0, 0, 0, 0)]),
            (192, 1, 1, [(6, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (272, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 17,
        nodes: &[
            (1, 1, 1, [(2, 1, 0, 1), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 1, 1, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 10,
        nodes: &[
            (1, 1, 1, [(2, 1, 0, 1), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(3, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 2, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 33,
        nodes: &[
            (1, 1, 1, [(2, 1, 0, 1), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (192, 1, 1, [(3, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (272, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 12,
        nodes: &[
            (1, 1, 0, [(2, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (272, 1, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 2,
        nodes: &[
            (1, 0, 0, [(2, 1, 1, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (16, 2, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
    CnListData {
        bits: 4,
        nodes: &[(1, 0, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)])],
    },
    CnListData {
        bits: 1,
        nodes: &[(1, 0, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)])],
    },
    CnListData {
        bits: -1,
        nodes: &[
            (1, 0, 0, [(2, 16, 0, 0), (4, 16, 0, 0), (0, 0, 0, 0)]),
            (192, 16, 16, [(3, 16, 0, 16), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (2064, 16, 16, [(0, 0, 0, 0), (0, 0, 0, 2), (0, 0, 0, 0)]),
            (2320, 16, 0, [(0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
            (2048, 0, 0, [(1, 3, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]),
        ],
    },
];
const CN_LIST_INITIAL_CHARGES: [i8; 18] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0];

// Exact `C_NODE.v.valence` values for the active cnList definitions. They are
// kept separately because the pre-existing compact node tuple stores type,
// capacity, flow, and edge records.
const CN_LIST_VALENCES: [&[i32]; 18] = [
    &[1, 3, 3, 2, 1],
    &[1, 3, 3, 2, 1],
    &[1, 3, 3, 2, 1],
    &[1, 3, 3, 2, 1],
    &[2, 1, 1],
    &[1, 3, 1, 1],
    &[1, 2, 3, 2, 1, 1],
    &[1, 3, 3, 2, 1],
    &[1, 3, 2, 1, 1],
    &[1, 3, 3, 3, 1, 1],
    &[1, 1],
    &[1, 2, 1],
    &[1, 2, 1],
    &[1, 1],
    &[1, 1],
    &[0],
    &[0],
    &[3, 2, 1, 1, 3],
];

// `nTautEndpointEdgeCap` reads only the first C_NODE of each active cnList entry.
// These are the exact `v.cap - v.flow` values from ichirvr1.c:102-479.
const CN_LIST_FIRST_NODE_FREE_VALENCE: [i32; 18] = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];

#[allow(non_snake_case)]
pub(crate) fn nTautEndpointEdgeCap(
    atoms: &[inp_ATOM],
    valence_atoms: &[crate::source_types::VAL_AT],
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1574 nTautEndpointEdgeCap
    // INCHI✔️✔️: int nTautEndpointEdgeCap( inp_ATOM *at, VAL_AT *pVA, int i )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* There are 3 sources of cap-flow = number of unsatisfied valences:
    // INCHI✔️✔️:        -----------------------------------------------------------------
    // INCHI✔️✔️:        1. pVA[i].cInitFreeValences
    // INCHI✔️✔️:        2. pCN[0].v.cap - pCN[0].v.flow
    // INCHI✔️✔️:        3. st[i].chem_bonds_valence - SUM(SINGLE, DOUBLE, TRIPLE bond orders)
    // INCHI✔️✔️:           Reasons: (a) This sum will not include 'ALTERN' bonds
    // INCHI✔️✔️:                    (b) until now at[i].chem_bonds_valence was used as a
    // INCHI✔️✔️:                        number of satisfied valences. In case of adjacent
    // INCHI✔️✔️:                        stereobonds marked as BOND_TYPE_ALTERN the value of
    // INCHI✔️✔️:                        at[i].chem_bonds_valence may be = at[i].valence+1.
    // INCHI✔️✔️:        4. Since tautomerism is defined for a neutral atom, do not add
    // INCHI✔️✔️:           initial flows from the atom to the ChargeStruct
    // INCHI✔️✔️:           CORRECTION: tautomeric endpoints do not have ChargeStruct.
    // INCHI✔️✔️:
    // INCHI✔️✔️:      */
    // INCHI✔️✔️:     int j, k, nEdgeCap, bonds_valence, stereo_bond_excess_valence;
    // INCHI✔️✔️:     MY_CONST C_NODE *pCN = pVA[i].cnListIndex > 0 ? cnList[pVA[i].cnListIndex - 1].pCN : NULL;
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* 1: free valences to reach the minimum known atom valence */
    // INCHI✔️✔️:     nEdgeCap = pVA[i].cInitFreeValences;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* 2: atom free valence in the ChargeStruct */
    // INCHI✔️✔️:     if (pCN)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nEdgeCap += pCN[0].v.cap - pCN[0].v.flow; /* normally should not happen */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* 3: atom free valence due to known from stereochemistry stereogenic bond types */
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     for ( j = 0, bonds_valence = 0; j < at[i].valence; j ++ )
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if ( at[i].bond_type[j] <= BOND_TYPE_TRIPLE )
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bonds_valence += at[i].bond_type[j];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* bonds > SINGLE are assumed fixed stereobonds; fixed bond cannot increase t-group edge flow */
    // INCHI✔️✔️:     for (stereo_bond_excess_valence = 0, j = 0; j < MAX_NUM_STEREO_BONDS && at[i].sb_parity[j]; j++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         k = at[i].sb_ord[j];
    // INCHI✔️✔️:         if (at[i].bond_type[k] < BOND_TYPE_TRIPLE)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             stereo_bond_excess_valence += at[i].bond_type[k] - BOND_TYPE_SINGLE;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     bonds_valence = (at[i].chem_bonds_valence - bonds_valence) + (bonds_valence -at[i].valence - stereo_bond_excess_valence);
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:     bonds_valence = ( at[i].chem_bonds_valence - at[i].valence ) - stereo_bond_excess_valence;
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*---- add 1, 2, 3 ----*/
    // INCHI✔️✔️:     if (bonds_valence >= 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nEdgeCap += bonds_valence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nEdgeCap = RI_ERR_PROGR;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nEdgeCap;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nTautEndpointEdgeCap

    let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence_atom = valence_atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;

    let mut edge_capacity = i32::from(valence_atom.cInitFreeValences);
    if valence_atom.cnListIndex > 0 {
        let charge_node_index = usize::try_from(i32::from(valence_atom.cnListIndex) - 1)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        edge_capacity = edge_capacity.wrapping_add(
            *CN_LIST_FIRST_NODE_FREE_VALENCE
                .get(charge_node_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
    }

    let mut stereo_bond_excess_valence = 0_i32;
    for stereo_index in 0..MAX_NUM_STEREO_BONDS as usize {
        if atom.sb_parity[stereo_index] == 0 {
            break;
        }
        let bond_index =
            usize::try_from(i32::from(atom.sb_ord[stereo_index])).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let bond_type = i32::from(
            *atom
                .bond_type
                .get(bond_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if bond_type < BOND_TYPE_TRIPLE as i32 {
            stereo_bond_excess_valence =
                stereo_bond_excess_valence.wrapping_add(bond_type.wrapping_sub(BOND_TYPE_SINGLE as i32));
        }
    }

    let bonds_valence = i32::from(atom.chem_bonds_valence)
        .wrapping_sub(i32::from(atom.valence))
        .wrapping_sub(stereo_bond_excess_valence);
    if bonds_valence >= 0 {
        Ok(edge_capacity.wrapping_add(bonds_valence))
    } else {
        Ok(RI_ERR_PROGR)
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
#[inline(always)]
pub(crate) fn BondFlowMaxcapMinorder(
    atoms: &[inp_ATOM],
    valence_atoms: &[crate::source_types::VAL_AT],
    restore_mode: &SRM,
    atom_index: i32,
    neighbor_order: i32,
    max_capacity: Option<&mut i32>,
    minimum_order: Option<&mut i32>,
    needs_flower: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1690 BondFlowMaxcapMinorder
    // INCHI✔️✔️: int BondFlowMaxcapMinorder( inp_ATOM *atom,
    // INCHI✔️✔️:                             VAL_AT *pVA,
    // INCHI✔️✔️:                             ICHICONST SRM *pSrm,
    // INCHI✔️✔️:                             int iat,
    // INCHI✔️✔️:                             int ineigh,
    // INCHI✔️✔️:                             int *pnMaxcap,
    // INCHI✔️✔️:                             int *pnMinorder,
    // INCHI✔️✔️:                             int *pbNeedsFlower )
    // INCHI✔️✔️: {
    let iat = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let ineigh = usize::try_from(neighbor_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if iat >= atoms.len() || iat >= valence_atoms.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    // SAFETY: the source-array prefix containing `iat` was validated once
    // above. Both slices remain immutably borrowed for this complete call.
    let atom = unsafe { atoms.get_unchecked(iat) };
    if ineigh >= atom.neighbor.len() || ineigh >= atom.bond_type.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    // SAFETY: the two fixed source arrays were validated for `ineigh` above.
    let neigh = usize::from(unsafe { *atom.neighbor.get_unchecked(ineigh) });
    if neigh >= atoms.len() || neigh >= valence_atoms.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    // SAFETY: `iat` and `neigh` are within both complete source arrays, which
    // are not mutated during this call.
    let neighbor = unsafe { atoms.get_unchecked(neigh) };
    let atom_va = unsafe { valence_atoms.get_unchecked(iat) };
    let neigh_va = unsafe { valence_atoms.get_unchecked(neigh) };

    // SAFETY: all four source-array references and `ineigh` were validated
    // above and remain immutable for the duration of the source calculation.
    Ok(unsafe {
        bond_flow_maxcap_minorder_source(
            atom,
            neighbor,
            atom_va,
            neigh_va,
            restore_mode,
            ineigh,
            max_capacity,
            minimum_order,
            needs_flower,
        )
    })
}

/// Executes the Official InChI calculation after its array operands have been
/// resolved by the owning source loop.
///
/// # Safety
///
/// `neighbor_order` must address both fixed bond arrays in `atom`. The four
/// references must describe the same atom and neighbor pair from the live
/// source arrays.
#[allow(non_snake_case, clippy::too_many_arguments)]
#[inline(always)]
pub(crate) unsafe fn bond_flow_maxcap_minorder_source(
    atom: &inp_ATOM,
    neighbor: &inp_ATOM,
    atom_va: &crate::source_types::VAL_AT,
    neigh_va: &crate::source_types::VAL_AT,
    restore_mode: &SRM,
    neighbor_order: usize,
    max_capacity: Option<&mut i32>,
    minimum_order: Option<&mut i32>,
    needs_flower: Option<&mut i32>,
) -> i32 {
    // INCHI✔️✔️:     int nFlow, nMaxcap, nMinorder, nInitorder, bNeedsFlower = 0;
    let mut flow;
    let max_capacity_value;
    let min_order;
    let initial_order;
    let mut flower = 0;
    // INCHI✔️✔️:     inp_ATOM *at = atom + iat;
    // INCHI✔️✔️:     int       neigh = at->neighbor[ineigh];
    // INCHI✔️✔️:     int       bond_type = at->bond_type[ineigh] & BOND_TYPE_MASK;
    let mut bond_type = i32::from(unsafe { *atom.bond_type.get_unchecked(neighbor_order) }) & BOND_TYPE_MASK as i32;
    // INCHI✔️✔️:     int       nMetal = ( 0 != pVA[iat].cMetal ) + ( 0 != pVA[neigh].cMetal );
    let n_metal = i32::from(atom_va.cMetal != 0) + i32::from(neigh_va.cMetal != 0);
    // INCHI✔️✔️:     int       nEndpoint = ( 0 != at->endpoint ) + ( 0 != atom[neigh].endpoint );
    let n_endpoint = i32::from(atom.endpoint != 0) + i32::from(neighbor.endpoint != 0);
    // INCHI✔️✔️:     int       nStereo = ( at->p_parity || at->sb_parity[0] ) + ( atom[neigh].p_parity || atom[neigh].sb_parity[0] );
    let n_stereo = i32::from(atom.p_parity != 0 || atom.sb_parity[0] != 0)
        + i32::from(neighbor.p_parity != 0 || neighbor.sb_parity[0] != 0);

    // INCHI✔️✔️:     if (bond_type > BOND_TYPE_TRIPLE)
    if bond_type > BOND_TYPE_TRIPLE as i32 {
        // INCHI✔️✔️:         bond_type = BOND_TYPE_SINGLE;
        bond_type = BOND_TYPE_SINGLE as i32;
    }

    // INCHI✔️✔️:     /* M=metal, A=non-metal atom, e=endpoint */
    // INCHI✔️✔️:     if ((nStereo && pSrm->bFixStereoBonds) || !nMetal || !pSrm->bMetalAddFlower)
    if (n_stereo != 0 && restore_mode.bFixStereoBonds != 0) || n_metal == 0 || restore_mode.bMetalAddFlower == 0 {
        // INCHI✔️✔️:         /* atom-atom rules, no metal atoms involved (1: A-A, A-Ae, Ae-Ae) */
        // INCHI✔️✔️:         nMinorder = BOND_TYPE_SINGLE;
        min_order = BOND_TYPE_SINGLE as i32;
        // INCHI✔️✔️:         nInitorder = bond_type;
        initial_order = bond_type;
        // INCHI✔️✔️:         nFlow = nInitorder - nMinorder;
        flow = initial_order.wrapping_sub(min_order);
    // INCHI✔️✔️:     else if (nMetal && !nEndpoint)
    } else if n_metal != 0 && n_endpoint == 0 {
        // INCHI✔️✔️:         /* M-a, M-M */
        // INCHI✔️✔️:         nMinorder = pSrm->nMetalMinBondOrder;
        min_order = restore_mode.nMetalMinBondOrder;
        // INCHI✔️✔️:         nInitorder = pSrm->nMetalInitBondOrder + bond_type - BOND_TYPE_SINGLE;
        initial_order = restore_mode
            .nMetalInitBondOrder
            .wrapping_add(bond_type)
            .wrapping_sub(BOND_TYPE_SINGLE as i32);
        // INCHI✔️✔️:         nFlow = nInitorder - nMinorder;
        flow = initial_order.wrapping_sub(min_order);
        // INCHI✔️✔️:         if (!pSrm->nMetalInitEdgeFlow &&
        // INCHI✔️✔️:             pSrm->nMetalInitBondOrder > pSrm->nMetalMinBondOrder && nFlow > 0)
        if restore_mode.nMetalInitEdgeFlow == 0
            && restore_mode.nMetalInitBondOrder > restore_mode.nMetalMinBondOrder
            && flow > 0
        {
            // INCHI✔️✔️:             nFlow--;
            flow = flow.wrapping_sub(1);
        }
        // INCHI✔️✔️:         bNeedsFlower = ( 0 != pVA[iat].cMetal );
        flower = i32::from(atom_va.cMetal != 0);
    // INCHI✔️✔️:     else if ((pVA[iat].cMetal && !at->endpoint && !pVA[neigh].cMetal && atom[neigh].endpoint) ||
    // INCHI✔️✔️:              (pVA[neigh].cMetal && !atom[neigh].endpoint && !pVA[iat].cMetal && at->endpoint))
    } else if (atom_va.cMetal != 0 && atom.endpoint == 0 && neigh_va.cMetal == 0 && neighbor.endpoint != 0)
        || (neigh_va.cMetal != 0 && neighbor.endpoint == 0 && atom_va.cMetal == 0 && atom.endpoint != 0)
    {
        // INCHI✔️✔️:             /* M-ae */
        // INCHI✔️✔️:             nMinorder = pSrm->nMetal2EndpointMinBondOrder;
        min_order = restore_mode.nMetal2EndpointMinBondOrder;
        // INCHI✔️✔️:             nInitorder = pSrm->nMetal2EndpointInitBondOrder + bond_type - BOND_TYPE_SINGLE;
        initial_order = restore_mode
            .nMetal2EndpointInitBondOrder
            .wrapping_add(bond_type)
            .wrapping_sub(BOND_TYPE_SINGLE as i32);
        // INCHI✔️✔️:             nFlow = nInitorder - nMinorder;
        flow = initial_order.wrapping_sub(min_order);
        // INCHI✔️✔️:             if (!pSrm->nMetal2EndpointInitEdgeFlow &&
        // INCHI✔️✔️:                 pSrm->nMetal2EndpointInitBondOrder > pSrm->nMetal2EndpointMinBondOrder && nFlow > 0)
        if restore_mode.nMetal2EndpointInitEdgeFlow == 0
            && restore_mode.nMetal2EndpointInitBondOrder > restore_mode.nMetal2EndpointMinBondOrder
            && flow > 0
        {
            // INCHI✔️✔️:                 nFlow--;
            flow = flow.wrapping_sub(1);
        }
        // INCHI✔️✔️:             bNeedsFlower = ( 0 != pVA[iat].cMetal );
        flower = i32::from(atom_va.cMetal != 0);
    // INCHI✔️✔️:         else
    } else {
        // INCHI✔️✔️:             /* endpoint is metal => no flower (4: M-Me, Me-Me, Me-A, Me-Ae) */
        // INCHI✔️✔️:             nMinorder = pSrm->nMetal2EndpointMinBondOrder;
        min_order = restore_mode.nMetal2EndpointMinBondOrder;
        // INCHI✔️✔️:             nInitorder = pSrm->nMetal2EndpointInitBondOrder + bond_type - BOND_TYPE_SINGLE;
        initial_order = restore_mode
            .nMetal2EndpointInitBondOrder
            .wrapping_add(bond_type)
            .wrapping_sub(BOND_TYPE_SINGLE as i32);
        // INCHI✔️✔️:             nFlow = nInitorder - nMinorder;
        flow = initial_order.wrapping_sub(min_order);
        // INCHI✔️✔️:             if (!pSrm->nMetal2EndpointInitEdgeFlow &&
        // INCHI✔️✔️:                 pSrm->nMetal2EndpointInitBondOrder > pSrm->nMetal2EndpointMinBondOrder && nFlow > 0)
        if restore_mode.nMetal2EndpointInitEdgeFlow == 0
            && restore_mode.nMetal2EndpointInitBondOrder > restore_mode.nMetal2EndpointMinBondOrder
            && flow > 0
        {
            // INCHI✔️✔️:                 nFlow--;
            flow = flow.wrapping_sub(1);
        }
        // INCHI✔️✔️:             bNeedsFlower = ( pVA[iat].cMetal && !at->endpoint );
        flower = i32::from(atom_va.cMetal != 0 && atom.endpoint == 0);
    }

    // INCHI✔️✔️:     nMaxcap = BOND_TYPE_TRIPLE - nMinorder;
    max_capacity_value = (BOND_TYPE_TRIPLE as i32).wrapping_sub(min_order);
    // INCHI✔️✔️:     if (pnMaxcap) { *pnMaxcap = nMaxcap; }
    if let Some(value) = max_capacity {
        *value = max_capacity_value;
    }
    // INCHI✔️✔️:     if (pnMinorder) { *pnMinorder = nMinorder; }
    if let Some(value) = minimum_order {
        *value = min_order;
    }
    // INCHI✔️✔️:     if (pbNeedsFlower) { *pbNeedsFlower = bNeedsFlower; }
    if let Some(value) = needs_flower {
        *value = flower;
    }
    // INCHI✔️✔️:     return nFlow;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: BondFlowMaxcapMinorder
    flow
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AtomStcapStflow(
    atoms: &[inp_ATOM],
    valence_atoms: &[crate::source_types::VAL_AT],
    restore_mode: &SRM,
    atom_index: i32,
    stationary_capacity: Option<&mut i32>,
    stationary_flow: Option<&mut i32>,
    metal_group_edge_capacity: Option<&mut i32>,
    metal_group_edge_flow: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1790 AtomStcapStflow
    // INCHI✔️✔️: int AtomStcapStflow( inp_ATOM *atom, VAL_AT *pVA, ICHICONST SRM *pSrm, int iat, int *pnStcap, int *pnStflow,
    // INCHI✔️✔️:                      EdgeFlow *pnMGroupEdgeCap, EdgeFlow *pnMGroupEdgeFlow )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int ineigh, bFlower;
    // INCHI✔️✔️:     int nStflow = 0, nMaxBondCap, nMinBondOrder, bNeedsFlower = 0;
    // INCHI✔️✔️:     int valence = atom[iat].valence;
    // INCHI✔️✔️:     int nStcap = atom[iat].chem_bonds_valence;
    // INCHI✔️✔️:     int nMGroupEdgeCap = 0, nMGroupEdgeFlow = 0, nFlow;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (pSrm->bMetalAddFlower)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nStcap -= pVA[iat].cInitOrigValenceToMetal - pVA[iat].cInitValenceToMetal;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (ineigh = 0; ineigh < valence; ineigh++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nFlow = BondFlowMaxcapMinorder( atom, pVA, pSrm, iat, ineigh, &nMaxBondCap, &nMinBondOrder, &bFlower );
    // INCHI✔️✔️:         nStflow += nFlow;
    // INCHI✔️✔️:         nStcap -= nMinBondOrder;
    // INCHI✔️✔️:         if (bFlower)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bNeedsFlower++;
    // INCHI✔️✔️:             nMGroupEdgeFlow += nFlow;
    // INCHI✔️✔️:             nMGroupEdgeCap += BOND_TYPE_TRIPLE - nMinBondOrder + pSrm->nMetalMaxCharge_D;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (pnStcap)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *pnStcap = bNeedsFlower ? nStflow : nStcap; /* initially, metal atoms are not radicals */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (pnStflow)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *pnStflow = nStflow;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (pnMGroupEdgeFlow)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *pnMGroupEdgeFlow = nMGroupEdgeCap - nMGroupEdgeFlow;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (pnMGroupEdgeCap)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *pnMGroupEdgeCap = nMGroupEdgeCap;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return bNeedsFlower; /* number of variable bonds to metal */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AtomStcapStflow

    let iat = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms.get(iat).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence_atom = valence_atoms.get(iat).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = i32::from(atom.valence);
    let mut st_capacity = i32::from(atom.chem_bonds_valence);
    let mut st_flow = 0_i32;
    let mut number_needing_flower = 0_i32;
    let mut metal_edge_capacity = 0_i32;
    let mut metal_edge_flow = 0_i32;

    if restore_mode.bMetalAddFlower != 0 {
        st_capacity = st_capacity.wrapping_sub(
            i32::from(valence_atom.cInitOrigValenceToMetal).wrapping_sub(i32::from(valence_atom.cInitValenceToMetal)),
        );
    }
    let mut ineigh = 0_i32;
    while ineigh < valence {
        let mut maximum_bond_capacity = 0_i32;
        let mut minimum_bond_order = 0_i32;
        let mut flower = 0_i32;
        let flow = BondFlowMaxcapMinorder(
            atoms,
            valence_atoms,
            restore_mode,
            atom_index,
            ineigh,
            Some(&mut maximum_bond_capacity),
            Some(&mut minimum_bond_order),
            Some(&mut flower),
        )?;
        st_flow = st_flow.wrapping_add(flow);
        st_capacity = st_capacity.wrapping_sub(minimum_bond_order);
        if flower != 0 {
            number_needing_flower = number_needing_flower.wrapping_add(1);
            metal_edge_flow = metal_edge_flow.wrapping_add(flow);
            metal_edge_capacity = metal_edge_capacity
                .wrapping_add((BOND_TYPE_TRIPLE as i32).wrapping_sub(minimum_bond_order))
                .wrapping_add(restore_mode.nMetalMaxCharge_D);
        }
        ineigh = ineigh.wrapping_add(1);
    }

    if let Some(value) = stationary_capacity {
        *value = if number_needing_flower != 0 {
            st_flow
        } else {
            st_capacity
        };
    }
    if let Some(value) = stationary_flow {
        *value = st_flow;
    }
    if let Some(value) = metal_group_edge_flow {
        *value = metal_edge_capacity.wrapping_sub(metal_edge_flow);
    }
    if let Some(value) = metal_group_edge_capacity {
        *value = metal_edge_capacity;
    }
    Ok(number_needing_flower)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn nCountBnsSizes(
    heap: &mut SourceHeap,
    atoms: &[inp_ATOM],
    number_of_atoms: i32,
    _add_edges_to_each_atom: i32,
    _add_vertices: i32,
    tautomer_info: &T_GROUP_INFO,
    valence_atoms: &[crate::source_types::VAL_AT],
    restore_mode: &SRM,
    groups: &mut ALL_TC_GROUPS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1870 nCountBnsSizes
    // INCHI✔️❌: int nCountBnsSizes( inp_ATOM *at, int num_at, int nAddEdges2eachAtom, int nAddVertices,
    // INCHI✔️❌:                     T_GROUP_INFO *ti, VAL_AT *pVA, ICHICONST SRM *pSrm, ALL_TC_GROUPS *pTCGroups )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, n, k, ret = 0, nBonds, nOtherEdges, nVertices, bMetalAtoms, bNeedsFlower;
    // INCHI✔️❌:     int nTgroupEdges, nTgroupEdgesFromTg, nTotNegChargInTgroups, cap, flow;
    // INCHI✔️❌:     MY_CONST C_NODE *pCN = NULL;
    // INCHI✔️❌:     nVertices = nBonds = nOtherEdges = nTgroupEdges = nTgroupEdgesFromTg = nTotNegChargInTgroups = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count metal atoms and electrons */
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pTCGroups->num_metal_atoms += ( pVA[i].cMetal != 0 );
    // INCHI✔️❌:         pTCGroups->num_metal_bonds += pVA[i].cNumBondsToMetal;
    // INCHI✔️❌:         pTCGroups->total_electrons += at[i].el_number;
    // INCHI✔️❌:         pTCGroups->total_electrons_metals += pVA[i].cMetal ? at[i].el_number : 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     pTCGroups->total_electrons -= pTCGroups->total_charge;
    // INCHI✔️❌:     pTCGroups->num_metal_bonds /= 2;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* register tautomeric groups */
    // INCHI✔️❌:     for (i = 0; i < ti->num_t_groups; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                BNS_VERT_TYPE_TGROUP,
    // INCHI✔️❌:                                ti->t_group[i].nGroupNumber,
    // INCHI✔️❌:                                ti->t_group[i].num[0] /* st_cap */,
    // INCHI✔️❌:                                0 /* st_flow */,
    // INCHI✔️❌:                                0 /* edge cap */,
    // INCHI✔️❌:                                0 /* edge flow */,
    // INCHI✔️❌:                                ti->t_group[i].nNumEndpoints /* num Edges */
    // INCHI✔️❌:         );
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* edges to tautomeric groups */
    // INCHI✔️❌:         nOtherEdges += ti->t_group[i].nNumEndpoints;
    // INCHI✔️❌:         nTgroupEdgesFromTg += ti->t_group[i].nNumEndpoints;
    // INCHI✔️❌:         /* total negative charge in t-groups */
    // INCHI✔️❌:         nTotNegChargInTgroups += ti->t_group[i].num[1];
    // INCHI✔️❌:         if (ret > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* should always happen since this is the first time this t-group is added */
    // INCHI✔️❌:             j = ret - 1;
    // INCHI✔️❌:             pTCGroups->pTCG[j].tg_num_H = ti->t_group[i].num[0] - ti->t_group[i].num[1];
    // INCHI✔️❌:             pTCGroups->pTCG[j].tg_num_Minus = ti->t_group[i].num[1];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     bMetalAtoms = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: repeat_for_metals:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count vertices and register ChargeValence groups */
    // INCHI✔️❌:     /* for now an atom may belong either to a t-group or to a ChargeValence group, but not to both */
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* number of bonds */
    // INCHI✔️❌:         nBonds += at[i].valence;
    // INCHI✔️❌:         /* Process ChargeStruct vertices and edges */
    // INCHI✔️❌:         if (pVA[i].cnListIndex)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* count vertices & edges in the ChargeValence Substructure attached to an atom */
    // INCHI✔️❌:             /* Important: unlike inp_ATOM, each edge e appears in pCN[*].e[*] only ONE time */
    // INCHI✔️❌:             int     len = cnList[j = pVA[i].cnListIndex - 1].len;
    // INCHI✔️❌:             int     bits = cnList[j].bits;
    // INCHI✔️❌:             int     type, neigh_type; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:             pCN = cnList[j].pCN;
    // INCHI✔️❌:
    // INCHI✔️❌:             /* first process all non-metals, after that -- all metals */
    // INCHI✔️❌:             if (( bits != cn_bits_Me ) != !bMetalAtoms)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             for (j = 0; j < len; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 type = pCN[j].v.type; /* ChargeStruct vertex type: atom is the first, c-groups are last */
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* process all pCN[j] neighbors */
    // INCHI✔️❌:                 for (k = 0; k < MAX_CN_VAL && ( n = pCN[j].e[k].neigh ); k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nOtherEdges++;  /* edges inside ChargeStruct */
    // INCHI✔️❌:                     n--; /* neighbor vertex position inside cnList[j].pCN */
    // INCHI✔️❌:                     neigh_type = pCN[n].v.type; /* type of the neighboring atom */
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (IS_BNS_VT_C_GR( neigh_type ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* register this edge to a CN-group vertex */
    // INCHI✔️❌:                         cap = !bMetalAtoms ? pCN[j].e[k].cap : pCN[j].e[k].cap ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                         flow = !bMetalAtoms ? pCN[j].e[k].flow : pCN[j].e[k].flow ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                                neigh_type,
    // INCHI✔️❌:                                                0 /* ord_num*/,
    // INCHI✔️❌:                                                0 /* st_cap */,
    // INCHI✔️❌:                                                0 /* st_flow */,
    // INCHI✔️❌:                                                cap /* edge cap*/,
    // INCHI✔️❌:                                                flow /* edge flow */,
    // INCHI✔️❌:                                                1 /* nNumEdges*/ );
    // INCHI✔️❌:                         if (ret < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto exit_function;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (ret > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* the group has just been created; add one more edge to (+/-) or supergroup */
    // INCHI✔️❌:                             ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                                    neigh_type,
    // INCHI✔️❌:                                                    0 /* ord_num*/,
    // INCHI✔️❌:                                                    0 /* st_cap */,
    // INCHI✔️❌:                                                    0 /* st_flow */,
    // INCHI✔️❌:                                                    0 /* edge cap*/,
    // INCHI✔️❌:                                                    0/* edge flow*/,
    // INCHI✔️❌:                                                    1 /* nNumEdges*/ );
    // INCHI✔️❌:                             if (ret < 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 goto exit_function;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nOtherEdges++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (IS_BNS_VT_C_GR( type ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* register this edge to a CN-group vertex; normally this does not happen */
    // INCHI✔️❌:
    // INCHI✔️❌:                         cap = !bMetalAtoms ? pCN[j].e[k].cap : pCN[j].e[k].cap ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                         flow = !bMetalAtoms ? pCN[j].e[k].flow : pCN[j].e[k].flow ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                         ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
    // INCHI✔️❌:                                                0 /* st_cap */, 0 /* st_flow */,
    // INCHI✔️❌:                                                cap /* edge cap*/, flow /* edge flow */, 1 /* nNumEdges*/ );
    // INCHI✔️❌:                         if (ret < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto exit_function;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (ret > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* the group has just been created; add one more edge to (+/-) or supergroup */
    // INCHI✔️❌:                             ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
    // INCHI✔️❌:                                                    0 /* st_cap */, 0 /* st_flow */,
    // INCHI✔️❌:                                                    0 /* edge cap*/, 0/* edge flow*/, 1 /* nNumEdges*/ );
    // INCHI✔️❌:                             if (ret < 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 goto exit_function;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nOtherEdges++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 } /* end of the current vertex pCN[j] neighbors */
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* process  pCN[j] vertex */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (type & BNS_VERT_TYPE_ATOM)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;  /* do not count regular atoms here */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (IS_BNS_VT_CHRG_STRUCT( type ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nVertices++;
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (pSrm->bMetalAddFlower && IS_BNS_VT_M_GR( type ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* special treatment: flow and cap are known as well as structure */
    // INCHI✔️❌:                     /* initial bond valence to metal is either 0 or 1 */
    // INCHI✔️❌:                     EdgeFlow nEdgeFlow, nEdgeCap;
    // INCHI✔️❌:
    // INCHI✔️❌:                     bNeedsFlower = AtomStcapStflow( at, pVA, pSrm, i, NULL /*pnStcap*/, NULL /*pnStflow*/,
    // INCHI✔️❌:                                                     &nEdgeCap, &nEdgeFlow );
    // INCHI✔️❌:                     if (!bNeedsFlower)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = RI_ERR_PROGR;
    // INCHI✔️❌:                         goto exit_function;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /*
    // INCHI✔️❌:                     GetAtomToMCGroupInitEdgeCapFlow( &nEdgeCap, &nEdgeFlow, pSrm, at,  pVA, i );
    // INCHI✔️❌:                     GetAtomToMCGroupInitEdgeCapFlow( &nEdgeCap, &nEdgeFlow, pSrm );
    // INCHI✔️❌:                     */
    // INCHI✔️❌:                     /* the 1st is the flower base */
    // INCHI✔️❌:                     /* atom - G0 edge and G0 vertex */
    // INCHI✔️❌:                     ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                            type,
    // INCHI✔️❌:                                            0 /* ord_num*/,
    // INCHI✔️❌:                                            /*pVA[i].cInitFreeValences*/
    // INCHI✔️❌:                                            0 /* st_cap */,
    // INCHI✔️❌:                                            0 /* st_flow */,
    // INCHI✔️❌:                                            (int) nEdgeCap,
    // INCHI✔️❌:                                            (int) nEdgeFlow,
    // INCHI✔️❌:                                            1 /* nNumEdges*/ );
    // INCHI✔️❌:                     if (ret < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto exit_function;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* count edge atom-G0 */
    // INCHI✔️❌:                     nOtherEdges++;
    // INCHI✔️❌:                     if (ret > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* first time registration: add G0-G1 and G0-G2 edges to G0 */
    // INCHI✔️❌:
    // INCHI✔️❌:                         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                                type,
    // INCHI✔️❌:                                                0 /* ord_num*/,
    // INCHI✔️❌:                                                0 /* st_cap */,
    // INCHI✔️❌:                                                0 /* st_flow */,
    // INCHI✔️❌:                                                0,/* edge cap*/
    // INCHI✔️❌:                                                0 /*edge flow*/,
    // INCHI✔️❌:                                                2 /* nNumEdges*/ );
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (ret < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto exit_function;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         /* first time registration: add G1; it has 3 edges */
    // INCHI✔️❌:
    // INCHI✔️❌:                         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                                type,
    // INCHI✔️❌:                                                1 /* ord_num*/,
    // INCHI✔️❌:                                                0 /* st_cap */,
    // INCHI✔️❌:                                                0 /* st_flow */,
    // INCHI✔️❌:                                                0,/* edge cap*/
    // INCHI✔️❌:                                                0 /*edge flow*/,
    // INCHI✔️❌:                                                3 /* nNumEdges*/ );
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (ret <= 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ret = !ret ? RI_ERR_PROGR : ret;
    // INCHI✔️❌:                             goto exit_function;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         /* first time registration: add G2; it has 3 edges */
    // INCHI✔️❌:
    // INCHI✔️❌:                         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                                type,
    // INCHI✔️❌:                                                2 /* ord_num*/,
    // INCHI✔️❌:                                                0 /* st_cap */,
    // INCHI✔️❌:                                                0 /* st_flow */,
    // INCHI✔️❌:                                                0,/* edge cap*/
    // INCHI✔️❌:                                                0 /*edge flow*/,
    // INCHI✔️❌:                                                3 /* nNumEdges*/ );
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (ret <= 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ret = !ret ? RI_ERR_PROGR : ret;
    // INCHI✔️❌:                             goto exit_function;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         /* first time registration: add G3; it has 2 edges */
    // INCHI✔️❌:
    // INCHI✔️❌:                         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                                type,
    // INCHI✔️❌:                                                3 /* ord_num*/,
    // INCHI✔️❌:                                                0 /* st_cap */,
    // INCHI✔️❌:                                                0 /* st_flow */,
    // INCHI✔️❌:                                                0,/* edge cap*/
    // INCHI✔️❌:                                                0 /*edge flow*/,
    // INCHI✔️❌:                                                2 /* nNumEdges*/ );
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (ret <= 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ret = !ret ? RI_ERR_PROGR : ret;
    // INCHI✔️❌:                             goto exit_function;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /* count added metal flower vertices: G0, G1, G2, G3 */
    // INCHI✔️❌:                         nVertices += 4;
    // INCHI✔️❌:                         /* count added metal flower edges: C0-C1, C0-C2, C1-C2, C1-C3, C2-C3 */
    // INCHI✔️❌:                         nOtherEdges += 5;
    // INCHI✔️❌:                         /* add connections of G0 to G1 and G2 */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 nVertices++; /* count BNS_VT_C_POS* types; all contain BNS_VERT_TYPE_C_GROUP bit */
    // INCHI✔️❌:                 if (!IS_BNS_VT_C_GR( type ))
    // INCHI✔️❌:                 {  /* check */
    // INCHI✔️❌:                     ret = RI_ERR_PROGR;
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* add st_cap and st_flow for a charge group */
    // INCHI✔️❌:                 cap = !bMetalAtoms ? pCN[j].v.cap : pCN[j].v.cap ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                 flow = !bMetalAtoms ? pCN[j].v.flow : pCN[j].v.flow ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                 ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
    // INCHI✔️❌:                                        cap /* st-cap*/, flow /* st-flow */,
    // INCHI✔️❌:                                        0 /* edge cap */, 0 /* edge flow */, 0 /* edges already counted */ );
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pCN = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* count edge caps to t-groups */
    // INCHI✔️❌:         if (at[i].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int nEdgeCap = nTautEndpointEdgeCap( at, pVA, i );
    // INCHI✔️❌:             nTgroupEdges++;
    // INCHI✔️❌:             if (nEdgeCap < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = nEdgeCap;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* add number of unsatisfied valences for a t-group; the unknown flow = 0 */
    // INCHI✔️❌:             ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                    BNS_VERT_TYPE_TGROUP,
    // INCHI✔️❌:                                    at[i].endpoint,
    // INCHI✔️❌:                                    0 /* st_cap */,
    // INCHI✔️❌:                                    0 /* st_flow */,
    // INCHI✔️❌:                                    nEdgeCap /* edge cap */,
    // INCHI✔️❌:                                    0 /* edge flow */,
    // INCHI✔️❌:                                    0 /* t-group edges have already been counted */ );
    // INCHI✔️❌:             if (ret < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!bMetalAtoms && pTCGroups->num_metal_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bMetalAtoms = 1;
    // INCHI✔️❌:         nBonds = 0; /* added 2006-05-15 */
    // INCHI✔️❌:         goto repeat_for_metals;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count real atoms and bonds */
    // INCHI✔️❌:     nBonds /= 2;
    // INCHI✔️❌:     pTCGroups->num_atoms = num_at;
    // INCHI✔️❌:     pTCGroups->num_bonds = nBonds;
    // INCHI✔️❌:
    // INCHI✔️❌:     pTCGroups->num_tgroups = ti->num_t_groups;
    // INCHI✔️❌:     pTCGroups->num_tgroup_edges = nTgroupEdges;
    // INCHI✔️❌:     pTCGroups->tgroup_charge = -nTotNegChargInTgroups;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 <= ret && nTgroupEdgesFromTg != nTgroupEdges)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = BNS_PROGRAM_ERR;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nVertices += num_at;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count other vertices */
    // INCHI✔️❌:     nVertices += ti->num_t_groups;
    // INCHI✔️❌:     nBonds += nOtherEdges;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* return edges and vertices */
    // INCHI✔️❌:     pTCGroups->nVertices = nVertices;
    // INCHI✔️❌:     pTCGroups->nEdges = nBonds;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: nCountBnsSizes
    let is_charge_group = |type_: i32| type_ & BNS_VT_C_POS_ALL as i32 == BNS_VERT_TYPE_C_GROUP as i32;
    let is_charge_structure =
        |type_: i32| type_ & BNS_VERT_TYPE__AUX as i32 != 0 && type_ & BNS_VERT_TYPE_TEMP as i32 != 0;
    let is_metal_group = |type_: i32| type_ == BNS_VERT_TYPE_METAL_GR as i32;

    let mut bonds = 0_i32;
    let mut other_edges = 0_i32;
    let mut vertices = 0_i32;
    let mut tgroup_edges = 0_i32;
    let mut tgroup_edges_from_groups = 0_i32;
    let mut negative_charge_in_tgroups = 0_i32;

    let mut i = 0_i32;
    while i < number_of_atoms {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let va = valence_atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        groups.num_metal_atoms = groups.num_metal_atoms.wrapping_add(i32::from(va.cMetal != 0));
        groups.num_metal_bonds = groups.num_metal_bonds.wrapping_add(i32::from(va.cNumBondsToMetal));
        groups.total_electrons = groups.total_electrons.wrapping_add(i32::from(atom.el_number));
        if va.cMetal != 0 {
            groups.total_electrons_metals = groups.total_electrons_metals.wrapping_add(i32::from(atom.el_number));
        }
        i = i.wrapping_add(1);
    }
    groups.total_electrons = groups.total_electrons.wrapping_sub(groups.total_charge);
    groups.num_metal_bonds /= 2;

    let tgroup_count =
        usize::try_from(tautomer_info.num_t_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let tautomer_groups = if tgroup_count == 0 {
        Vec::new()
    } else {
        heap.slice(tautomer_info.t_group.as_const())?
            .get(..tgroup_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    for tautomer_group in &tautomer_groups {
        let result = RegisterTCGroup(
            heap,
            groups,
            BNS_VERT_TYPE_TGROUP as i32,
            i32::from(tautomer_group.nGroupNumber),
            i32::from(tautomer_group.num[0]),
            0,
            0,
            0,
            i32::from(tautomer_group.nNumEndpoints),
        )?;
        if result < 0 {
            return Ok(result);
        }
        other_edges = other_edges.wrapping_add(i32::from(tautomer_group.nNumEndpoints));
        tgroup_edges_from_groups = tgroup_edges_from_groups.wrapping_add(i32::from(tautomer_group.nNumEndpoints));
        negative_charge_in_tgroups = negative_charge_in_tgroups.wrapping_add(i32::from(tautomer_group.num[1]));
        if result > 0 {
            let group_index = usize::try_from(result - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let group = heap
                .slice_mut(groups.pTCG)?
                .get_mut(group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            group.tg_num_H = tautomer_group.num[0].wrapping_sub(tautomer_group.num[1]) as i16;
            group.tg_num_Minus = tautomer_group.num[1] as i16;
        }
    }

    let mut metal_pass = false;
    loop {
        let mut atom_number = 0_i32;
        while atom_number < number_of_atoms {
            let index = usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let va = valence_atoms.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            bonds = bonds.wrapping_add(i32::from(atom.valence));

            if va.cnListIndex != 0 {
                let list_index =
                    usize::try_from(i32::from(va.cnListIndex) - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let list = CN_LIST.get(list_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                if (list.bits != -1) != !metal_pass {
                    atom_number = atom_number.wrapping_add(1);
                    continue;
                }
                for (node_index, node) in list.nodes.iter().enumerate() {
                    let (type_, node_cap, node_flow, edges) = *node;
                    for (neighbor_number, edge_cap, _forbidden, edge_flow) in edges {
                        if neighbor_number == 0 {
                            break;
                        }
                        other_edges = other_edges.wrapping_add(1);
                        let neighbor_index =
                            usize::try_from(neighbor_number - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let neighbor_type = list
                            .nodes
                            .get(neighbor_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .0;
                        for group_type in [neighbor_type, type_] {
                            if is_charge_group(group_type) {
                                let cap = if !metal_pass {
                                    edge_cap
                                } else if edge_cap != 0 {
                                    restore_mode.nMetalMaxCharge_D
                                } else {
                                    0
                                };
                                let flow = if !metal_pass {
                                    edge_flow
                                } else if edge_flow != 0 {
                                    restore_mode.nMetalMaxCharge_D
                                } else {
                                    0
                                };
                                let mut result = RegisterTCGroup(heap, groups, group_type, 0, 0, 0, cap, flow, 1)?;
                                if result < 0 {
                                    return Ok(result);
                                }
                                if result > 0 {
                                    result = RegisterTCGroup(heap, groups, group_type, 0, 0, 0, 0, 0, 1)?;
                                    if result < 0 {
                                        return Ok(result);
                                    }
                                    other_edges = other_edges.wrapping_add(1);
                                }
                            }
                        }
                    }

                    if type_ & BNS_VERT_TYPE_ATOM as i32 != 0 {
                        continue;
                    }
                    if is_charge_structure(type_) {
                        vertices = vertices.wrapping_add(1);
                        continue;
                    }
                    if restore_mode.bMetalAddFlower != 0 && is_metal_group(type_) {
                        let mut edge_capacity = 0_i32;
                        let mut edge_flow = 0_i32;
                        let needs_flower = AtomStcapStflow(
                            atoms,
                            valence_atoms,
                            restore_mode,
                            atom_number,
                            None,
                            None,
                            Some(&mut edge_capacity),
                            Some(&mut edge_flow),
                        )?;
                        if needs_flower == 0 {
                            return Ok(RI_ERR_PROGR);
                        }
                        let result = RegisterTCGroup(heap, groups, type_, 0, 0, 0, edge_capacity, edge_flow, 1)?;
                        if result < 0 {
                            return Ok(result);
                        }
                        other_edges = other_edges.wrapping_add(1);
                        if result > 0 {
                            let result = RegisterTCGroup(heap, groups, type_, 0, 0, 0, 0, 0, 2)?;
                            if result < 0 {
                                return Ok(result);
                            }
                            for (order, edges) in [(1, 3), (2, 3), (3, 2)] {
                                let result = RegisterTCGroup(heap, groups, type_, order, 0, 0, 0, 0, edges)?;
                                if result <= 0 {
                                    return Ok(if result == 0 { RI_ERR_PROGR } else { result });
                                }
                            }
                            vertices = vertices.wrapping_add(4);
                            other_edges = other_edges.wrapping_add(5);
                        }
                        continue;
                    }

                    vertices = vertices.wrapping_add(1);
                    if !is_charge_group(type_) {
                        return Ok(RI_ERR_PROGR);
                    }
                    let cap = if !metal_pass {
                        node_cap
                    } else if node_cap != 0 {
                        restore_mode.nMetalMaxCharge_D
                    } else {
                        0
                    };
                    let flow = if !metal_pass {
                        node_flow
                    } else if node_flow != 0 {
                        restore_mode.nMetalMaxCharge_D
                    } else {
                        0
                    };
                    let result = RegisterTCGroup(heap, groups, type_, 0, cap, flow, 0, 0, 0)?;
                    if result < 0 {
                        return Ok(result);
                    }
                    let _ = node_index;
                }
            }

            if atom.endpoint != 0 {
                let edge_capacity = nTautEndpointEdgeCap(atoms, valence_atoms, atom_number)?;
                tgroup_edges = tgroup_edges.wrapping_add(1);
                if edge_capacity < 0 {
                    return Ok(edge_capacity);
                }
                let result = RegisterTCGroup(
                    heap,
                    groups,
                    BNS_VERT_TYPE_TGROUP as i32,
                    i32::from(atom.endpoint),
                    0,
                    0,
                    edge_capacity,
                    0,
                    0,
                )?;
                if result < 0 {
                    return Ok(result);
                }
            }
            atom_number = atom_number.wrapping_add(1);
        }
        if !metal_pass && groups.num_metal_atoms != 0 {
            metal_pass = true;
            bonds = 0;
            continue;
        }
        break;
    }

    bonds /= 2;
    groups.num_atoms = number_of_atoms;
    groups.num_bonds = bonds;
    groups.num_tgroups = tautomer_info.num_t_groups;
    groups.num_tgroup_edges = tgroup_edges;
    groups.tgroup_charge = negative_charge_in_tgroups.wrapping_neg();
    let mut result = 0;
    if tgroup_edges_from_groups != tgroup_edges {
        result = BNS_PROGRAM_ERR;
    }
    vertices = vertices.wrapping_add(number_of_atoms);
    vertices = vertices.wrapping_add(tautomer_info.num_t_groups);
    bonds = bonds.wrapping_add(other_edges);
    groups.nVertices = vertices;
    groups.nEdges = bonds;
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn nAddSuperCGroups(heap: &mut SourceHeap, groups: &mut ALL_TC_GROUPS) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2279 nAddSuperCGroups
    // INCHI✔️❌: int nAddSuperCGroups( ALL_TC_GROUPS *pTCGroups )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, k, n, n1, n2, n3, ret = 0, nNumToConnect; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < pTCGroups->num_tc_groups; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pTCGroups->pTCG[i].type & BNS_VERT_TYPE_TGROUP)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             continue; /* t-group */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (IS_BNS_VT_C_GR( pTCGroups->pTCG[i].type ) ||
    // INCHI✔️❌:              IS_BNS_VT_M_GR( pTCGroups->pTCG[i].type ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* ChargeValence (cn) group */
    // INCHI✔️❌:             switch (pTCGroups->pTCG[i].type)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case BNS_VT_C_POS:
    // INCHI✔️❌:                     k = TCG_Plus0;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case BNS_VT_C_NEG:
    // INCHI✔️❌:                     k = TCG_Minus0;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case BNS_VT_C_POS_C:
    // INCHI✔️❌:                     k = TCG_Plus_C0;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case BNS_VT_C_NEG_C:
    // INCHI✔️❌:                     k = TCG_Minus_C0;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case BNS_VT_C_POS_M:
    // INCHI✔️❌:                     k = TCG_Plus_M0;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case BNS_VT_C_NEG_M:
    // INCHI✔️❌:                     k = TCG_Minus_M0;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case BNS_VT_M_GROUP:
    // INCHI✔️❌:                     switch (pTCGroups->pTCG[i].ord_num)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         case 0:
    // INCHI✔️❌:                             k = TCG_MeFlower0;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case 1:
    // INCHI✔️❌:                             k = TCG_MeFlower1;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case 2:
    // INCHI✔️❌:                             k = TCG_MeFlower2;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case 3:
    // INCHI✔️❌:                             k = TCG_MeFlower3;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         default:
    // INCHI✔️❌:                             ret = RI_ERR_PROGR; /* unexpected group type */
    // INCHI✔️❌:                             goto exit_function;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:
    // INCHI✔️❌:                 default:
    // INCHI✔️❌:                     ret = RI_ERR_PROGR; /* unexpected group type */
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (pTCGroups->nGroup[k] >= 0 || (pTCGroups->pTCG[i].ord_num && !IS_BNS_VT_M_GR( pTCGroups->pTCG[i].type ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = RI_ERR_PROGR;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             pTCGroups->nGroup[k] = i; /* ordering number of the Charge group, starting from 0 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* add (+) supergroup */
    // INCHI✔️❌:     n1 = pTCGroups->nGroup[TCG_Plus0];
    // INCHI✔️❌:     n2 = pTCGroups->nGroup[TCG_Plus_C0];
    // INCHI✔️❌:     n3 = pTCGroups->nGroup[TCG_Plus_M0];
    // INCHI✔️❌:     nNumToConnect = ( n1 >= 0 ) + ( n2 >= 0 ) + ( n3 >= 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumToConnect)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* if both groups are present then add a supergroup */
    // INCHI✔️❌:         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                BNS_VT_C_POS_ALL,
    // INCHI✔️❌:                                0,
    // INCHI✔️❌:                                0 /* st_cap */,
    // INCHI✔️❌:                                0 /* st_flow */,
    // INCHI✔️❌:                                0 /* edge cap */,
    // INCHI✔️❌:                                0 /* edge flow */,
    // INCHI✔️❌:                                1 + nNumToConnect
    // INCHI✔️❌:                                 /* one more edge to connect to */
    // INCHI✔️❌:         /* an additional (+/-) vertex */ );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ret <= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = !ret ? RI_ERR_PROGR : ret;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         pTCGroups->nGroup[TCG_Plus] = ret - 1; /* newly added group number */
    // INCHI✔️❌:         pTCGroups->nVertices += 2; /* two vertices including itself */
    // INCHI✔️❌:         pTCGroups->nEdges += 1 + nNumToConnect; /* one more edge to connect to an additional (+/-) vertex */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* add (-) supergroup */
    // INCHI✔️❌:     n1 = pTCGroups->nGroup[TCG_Minus0];
    // INCHI✔️❌:     n2 = pTCGroups->nGroup[TCG_Minus_C0];
    // INCHI✔️❌:     n3 = pTCGroups->nGroup[TCG_Minus_M0];
    // INCHI✔️❌:     nNumToConnect = ( n1 >= 0 ) + ( n2 >= 0 ) + ( n3 >= 0 );
    // INCHI✔️❌:     if (nNumToConnect)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* if both groups are present then add a supergroup */
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = RegisterTCGroup( pTCGroups,
    // INCHI✔️❌:                                BNS_VT_C_NEG_ALL,
    // INCHI✔️❌:                                0,
    // INCHI✔️❌:                                0 /* st_cap */,
    // INCHI✔️❌:                                0 /* st_flow */,
    // INCHI✔️❌:                                0 /* edge cap */,
    // INCHI✔️❌:                                0 /* edge flow */,
    // INCHI✔️❌:                                1 + nNumToConnect /* one more edge to connect to an additional (+/-) vertex */ );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pTCGroups->nGroup[TCG_Minus] = ret - 1; /* newly added group number */
    // INCHI✔️❌:         pTCGroups->nVertices += 2; /* needs two vertices including itself */
    // INCHI✔️❌:         pTCGroups->nEdges += 1 + nNumToConnect; /* one more edge to connect to an additional (+/-) vertex */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* add neutralization vertex: (+)-()=(-) connection */
    // INCHI✔️❌:     k = pTCGroups->nGroup[TCG_Minus];
    // INCHI✔️❌:     n = pTCGroups->nGroup[TCG_Plus];
    // INCHI✔️❌:     nNumToConnect = ( k >= 0 ) + ( n >= 0 );
    // INCHI✔️❌:     if (nNumToConnect)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pTCGroups->nVertices += 1;
    // INCHI✔️❌:         pTCGroups->nEdges += nNumToConnect; /* one edge per super-c-group */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: nAddSuperCGroups
    let count = usize::try_from(groups.num_tc_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let existing = if count == 0 {
        Vec::new()
    } else {
        heap.slice(groups.pTCG.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    for (index, group) in existing.iter().enumerate() {
        if group.type_ & BNS_VERT_TYPE_TGROUP as i32 != 0 {
            continue;
        }
        let is_charge_group = group.type_ & BNS_VT_C_POS_ALL as i32 == BNS_VERT_TYPE_C_GROUP as i32;
        let is_metal_group = group.type_ == BNS_VERT_TYPE_METAL_GR as i32;
        if is_charge_group || is_metal_group {
            let slot = match group.type_ as u32 {
                BNS_VT_C_POS => TCG_Plus0,
                BNS_VT_C_NEG => TCG_Minus0,
                BNS_VT_C_POS_C => TCG_Plus_C0,
                BNS_VT_C_NEG_C => TCG_Minus_C0,
                BNS_VT_C_POS_M => TCG_Plus_M0,
                BNS_VT_C_NEG_M => TCG_Minus_M0,
                BNS_VT_M_GROUP => match group.ord_num {
                    0 => TCG_MeFlower0,
                    1 => TCG_MeFlower1,
                    2 => TCG_MeFlower2,
                    3 => TCG_MeFlower3,
                    _ => return Ok(RI_ERR_PROGR),
                },
                _ => return Ok(RI_ERR_PROGR),
            } as usize;
            if groups.nGroup[slot] >= 0 || (group.ord_num != 0 && !is_metal_group) {
                return Ok(RI_ERR_PROGR);
            }
            groups.nGroup[slot] = index as i32;
        }
    }

    for (slots, super_type, super_slot, strict_new) in [
        ([TCG_Plus0, TCG_Plus_C0, TCG_Plus_M0], BNS_VT_C_POS_ALL, TCG_Plus, true),
        (
            [TCG_Minus0, TCG_Minus_C0, TCG_Minus_M0],
            BNS_VT_C_NEG_ALL,
            TCG_Minus,
            false,
        ),
    ] {
        let number_to_connect = slots
            .into_iter()
            .map(|slot| i32::from(groups.nGroup[slot as usize] >= 0))
            .sum::<i32>();
        if number_to_connect != 0 {
            let result = RegisterTCGroup(
                heap,
                groups,
                super_type as i32,
                0,
                0,
                0,
                0,
                0,
                1_i32.wrapping_add(number_to_connect),
            )?;
            if (strict_new && result <= 0) || (!strict_new && result < 0) {
                return Ok(if result == 0 { RI_ERR_PROGR } else { result });
            }
            groups.nGroup[super_slot as usize] = result.wrapping_sub(1);
            groups.nVertices = groups.nVertices.wrapping_add(2);
            groups.nEdges = groups.nEdges.wrapping_add(1_i32.wrapping_add(number_to_connect));
        }
    }

    let number_to_connect = i32::from(groups.nGroup[TCG_Minus as usize] >= 0)
        .wrapping_add(i32::from(groups.nGroup[TCG_Plus as usize] >= 0));
    if number_to_connect != 0 {
        groups.nVertices = groups.nVertices.wrapping_add(1);
        groups.nEdges = groups.nEdges.wrapping_add(number_to_connect);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn GetDeltaChargeFromVF(
    heap: &SourceHeap,
    pBNS: &BN_STRUCT,
    pVA: &[VAL_AT],
    vf: &mut VF,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4229 GetDeltaChargeFromVF
    /*
    int GetDeltaChargeFromVF( BN_STRUCT *pBNS, VAL_AT *pVA, VF *vf )
    {
        int i, v = NO_VERTEX;
        int ieIn1 = ( !( vf->bUsed & VF_USED_IN ) && vf->e_In >= 0 && vf->delta_In ) ? vf->e_In + 1 : NO_VERTEX;
        int ieOut1 = ( !( vf->bUsed & VF_USED_OUT ) && vf->e_Out >= 0 && vf->delta_Out ) ? vf->e_Out + 1 : NO_VERTEX;
        int nInitCharge, nPlusFlow, nMinusFlow, nDeltaCharge, nNumDeltaCharge, eCPlus, eCMinus;

        if (!( vf->type & BNS_VERT_TYPE_C_GROUP ) ||
            ( vf->type & BNS_VERT_TYPE_SUPER_TGROUP ) ||
            ( ieIn1 == NO_VERTEX && ieOut1 == NO_VERTEX ))
        {
            return 0;
        }
        if (vf->type & BNS_VERT_TYPE_C_NEGATIVE)
        {
            /* negative charge edge */
            for (i = 0; i < pBNS->num_atoms; i++)
            {
                if (pVA[i].nCMinusGroupEdge == ieIn1 || pVA[i].nCMinusGroupEdge == ieOut1)
                {
                    v = i;
                    break;
                }
            }
        }
        else
        {
            /* positive charge edge */
            for (i = 0; i < pBNS->num_atoms; i++)
            {
                if (pVA[i].nCPlusGroupEdge == ieIn1 || pVA[i].nCPlusGroupEdge == ieOut1)
                {
                    v = i;
                    break;
                }
            }
        }

        if (v == NO_VERTEX)
            return 0;

        nInitCharge = pVA[v].cInitCharge;
        nPlusFlow = nMinusFlow = 0;
        nNumDeltaCharge = 0;

        if (( eCPlus = pVA[v].nCPlusGroupEdge - 1 ) >= 0)
        {
            nPlusFlow = pBNS->edge[eCPlus].cap
                - pBNS->edge[eCPlus].flow;
        }

        if (( eCMinus = pVA[v].nCMinusGroupEdge - 1 ) >= 0)
        {
            nMinusFlow = -pBNS->edge[eCMinus].flow;
        }
        nInitCharge += nPlusFlow + nMinusFlow;

        nDeltaCharge = 0;

        if (!( vf[0].bUsed & VF_USED_OUT ))
        {
            if (vf[0].e_Out == eCPlus || vf[0].e_Out == eCMinus)
            {
                nDeltaCharge -= vf[0].delta_Out;
                vf[0].bUsed |= VF_USED_OUT;
            }
        }

        if (!( vf[0].bUsed & VF_USED_IN ))
        {
            if (vf[0].e_In == eCPlus || vf[0].e_In == eCMinus)
            {
                nDeltaCharge -= vf[0].delta_In;
                vf[0].bUsed |= VF_USED_IN;
            }
        }

        if (!nInitCharge && nDeltaCharge)
        {
            nNumDeltaCharge++;
        }
        else if (nInitCharge && 0 == nInitCharge + nDeltaCharge)
        {
            nNumDeltaCharge--;
        }

        return nNumDeltaCharge;
    }
    */
    // END INCHI C FUNCTION: GetDeltaChargeFromVF
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetDeltaChargeFromVF
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: #define VF_USED_IN 1
    // INCHI✔️❌: #define VF_USED_OUT 2
    // INCHI✔️❌: #define NO_VERTEX (-2)
    // INCHI✔️❌: SourceHeap allocation lookup adds overhead versus direct C edge indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetDeltaChargeFromVF

    let input_edge = if vf.bUsed & VF_USED_IN as i32 == 0 && vf.e_In >= 0 && vf.delta_In != 0 {
        vf.e_In.wrapping_add(1)
    } else {
        NO_VERTEX
    };
    let output_edge = if vf.bUsed & VF_USED_OUT as i32 == 0 && vf.e_Out >= 0 && vf.delta_Out != 0 {
        vf.e_Out.wrapping_add(1)
    } else {
        NO_VERTEX
    };

    if vf.type_ & BNS_VERT_TYPE_C_GROUP as i32 == 0
        || vf.type_ & BNS_VERT_TYPE_SUPER_TGROUP as i32 != 0
        || (input_edge == NO_VERTEX && output_edge == NO_VERTEX)
    {
        return Ok(0);
    }

    let negative_group = vf.type_ & BNS_VERT_TYPE_C_NEGATIVE as i32 != 0;
    let mut vertex = NO_VERTEX;
    let mut atom_index = 0_i32;
    while atom_index < pBNS.num_atoms {
        let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let charge_edge = if negative_group {
            atom.nCMinusGroupEdge
        } else {
            atom.nCPlusGroupEdge
        };
        if charge_edge == input_edge || charge_edge == output_edge {
            vertex = atom_index;
            break;
        }
        atom_index = atom_index.wrapping_add(1);
    }
    if vertex == NO_VERTEX {
        return Ok(0);
    }

    let atom = pVA
        .get(usize::try_from(vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let plus_edge = atom.nCPlusGroupEdge.wrapping_sub(1);
    let minus_edge = atom.nCMinusGroupEdge.wrapping_sub(1);
    let plus_flow = if plus_edge >= 0 {
        let edge = heap
            .slice(pBNS.edge.as_const())?
            .get(usize::try_from(plus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        edge.cap.wrapping_sub(edge.flow)
    } else {
        0
    };
    let minus_flow = if minus_edge >= 0 {
        let edge = heap
            .slice(pBNS.edge.as_const())?
            .get(usize::try_from(minus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        edge.flow.wrapping_neg()
    } else {
        0
    };
    let initial_charge = i32::from(atom.cInitCharge)
        .wrapping_add(plus_flow)
        .wrapping_add(minus_flow);

    let mut delta_charge = 0_i32;
    if vf.bUsed & VF_USED_OUT as i32 == 0 && (vf.e_Out == plus_edge || vf.e_Out == minus_edge) {
        delta_charge = delta_charge.wrapping_sub(vf.delta_Out);
        vf.bUsed |= VF_USED_OUT as i32;
    }
    if vf.bUsed & VF_USED_IN as i32 == 0 && (vf.e_In == plus_edge || vf.e_In == minus_edge) {
        delta_charge = delta_charge.wrapping_sub(vf.delta_In);
        vf.bUsed |= VF_USED_IN as i32;
    }

    if initial_charge == 0 && delta_charge != 0 {
        Ok(1)
    } else if initial_charge != 0 && initial_charge.wrapping_add(delta_charge) == 0 {
        Ok(-1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case)]
pub(crate) fn EvaluateChargeChanges(
    heap: &SourceHeap,
    pBNS: &mut BN_STRUCT,
    pVA: &[VAL_AT],
    pnDeltaH: &mut i32,
    pnDeltaCharge: &mut i32,
    pnNumVisitedAtoms: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4320 EvaluateChargeChanges
    /*
    int EvaluateChargeChanges( BN_STRUCT *pBNS, VAL_AT *pVA, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms )
    {
        int       pass, i, j, v0, v1, v2, v, ineigh1, /*ineigh2,*/ vLast, n, delta, ret, ie, err = 0;
        BNS_EDGE *edge;
        int       nDeltaH, nDeltaCharge, iPrev, nInitCharge, nPlusFlow, nMinusFlow;
        int       nNumDeltaH = 0;
        int       nNumDeltaCharge = 0;
        int       nNumVisitedAtoms = 0;
        VF   vf[NUM_VF + 1];

        *pnDeltaH = 0;
        *pnDeltaCharge = 0;
        *pnNumVisitedAtoms = 0;

        for (pass = pBNS->num_altp - 1, ret = 0; 0 <= pass; pass--)
        {

            pBNS->alt_path = pBNS->altp[pass];
            v1 = ALTP_START_ATOM( pBNS->alt_path );
            n = ALTP_PATH_LEN( pBNS->alt_path );
            delta = ALTP_DELTA( pBNS->alt_path );
            vLast = ALTP_END_ATOM( pBNS->alt_path );
            v0 = v2 = NO_VERTEX;

            memset( vf, 0, sizeof( vf ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            for (i = 0; i < (int) ( sizeof( vf ) / sizeof( vf[0] ) ); i++)
            {
                vf[i].v = NO_VERTEX; /* = -2 */
                vf[i].e_In = NO_VERTEX;
                vf[i].e_Out = NO_VERTEX;
            }
            iPrev = 0;
            /* add to the queue */
            if (ANY_VERT_TYPE( pBNS->vert[v1].type ))
            {
                if (pBNS->vert[v1].type & BNS_VERT_TYPE_ATOM)
                {
                    nNumVisitedAtoms++;
                }
                vf[2].type = pBNS->vert[v1].type;
                vf[2].v = v1;
                iPrev = 2;
            }

            nNumDeltaH = 0;
            nNumDeltaCharge = 0;
            nNumVisitedAtoms = 0;

            for (i = 0; i < n; i++, delta = -delta, v0 = v1, v1 = v2) /* djb-rwth: removing redundant code */
            {
                ineigh1 = ALTP_THIS_ATOM_NEIGHBOR( pBNS->alt_path, i );  /* v1->v2 neighbor */
                /*ineigh2 = ALTP_NEXT_ATOM_NEIGHBOR(pBNS->alt_path, i);*/  /* v2->v1 neighbor */
                edge = pBNS->edge + ( ie = pBNS->vert[v1].iedge[ineigh1] );
                /* follow the BN Structure, not the inp_ATOM, to take care of swithching to
                   t-groups, c-groups or other fictitious edges/vertices
                */

                if (iPrev)
                {
                    /* add exit delta and edge */
                    vf[2].e_Out = ie;
                    vf[2].delta_Out = delta;
                }

                v2 = edge->neighbor12 ^ v1;  /* next vertex */
                if (pBNS->vert[v2].type & BNS_VERT_TYPE_ATOM)
                {
                    nNumVisitedAtoms++;
                }

                if (( ANY_VERT_TYPE( pBNS->vert[v2].type ) || i == n - 1 ) &&
                    ( vf[0].type & BNS_VERT_TYPE_C_GROUP ) && vf[0].bUsed != VF_USED_ALL)
                {
                    /* unused vertex is about to be discarded */
                    nNumDeltaCharge += GetDeltaChargeFromVF( pBNS, pVA, &vf[0] );
                }

                if (ANY_VERT_TYPE( pBNS->vert[v2].type ))
                {
                    /* shift the queue */
                    vf[0] = vf[1];
                    vf[1] = vf[2];
                    vf[2] = vf[3]; /* make vf[2] empty */
                    /* add next vertex */
                    vf[2].v = v2;
                    vf[2].type = pBNS->vert[v2].type;
                    vf[2].e_In = ie;
                    vf[2].delta_In = delta;
                    iPrev = 2; /* indicates a newly added vertex */
                }
                else if (i == n - 1)
                {
                    /* shift the queue */
                    vf[0] = vf[1];
                    vf[1] = vf[2];
                    vf[2] = vf[3]; /* make vf[2] empty */
                    iPrev = 1; /* indicates the last vertex */
                }
                else
                {
                    iPrev = 0; /* no new vertex has been added */
                }

                if (iPrev && ( vf[1].type & BNS_VERT_TYPE_ATOM ))
                {
                    /* a new vertex has just been added and  */
                    /* an atom is in the middle of the queue */
                    EdgeIndex eCPlus, eCMinus;
                    v = vf[1].v;
                    nInitCharge = pVA[v].cInitCharge;
                    nPlusFlow = nMinusFlow = 0;
                    if (( eCPlus = pVA[v].nCPlusGroupEdge - 1 ) >= 0)
                    {
                        nPlusFlow = pBNS->edge[eCPlus].cap
                            - pBNS->edge[eCPlus].flow;
                    }
                    if (( eCMinus = pVA[v].nCMinusGroupEdge - 1 ) >= 0)
                    {
                        nMinusFlow = -pBNS->edge[eCMinus].flow;
                    }
                    nInitCharge += nPlusFlow + nMinusFlow;

                    nDeltaH = nDeltaCharge = 0;

                    if (vf[0].type & BNS_VERT_TYPE_TGROUP)
                    {
                        nDeltaH -= delta;
                    }
                    else if (( vf[0].type & BNS_VERT_TYPE_C_GROUP ) &&
                             !( vf[0].bUsed & VF_USED_OUT ))
                    {
                        if (vf[0].e_Out == eCPlus || vf[0].e_Out == eCMinus)
                        {
                            nDeltaCharge -= vf[0].delta_Out;
                            vf[0].bUsed |= VF_USED_OUT;
                        }
                    }

                    if (vf[2].type & BNS_VERT_TYPE_TGROUP)
                    {
                        nDeltaH += delta;
                    }
                    else if (( vf[2].type & BNS_VERT_TYPE_C_GROUP ) &&
                             !( vf[2].bUsed & VF_USED_IN ))
                    {
                        if (vf[2].e_In == eCPlus || vf[2].e_In == eCMinus)
                        {
                            nDeltaCharge -= vf[2].delta_In;
                            vf[2].bUsed |= VF_USED_IN;
                        }
                    }

                    if (!nInitCharge && nDeltaCharge)
                    {
                        nNumDeltaCharge++;
                    }
                    else if (nInitCharge && 0 == nInitCharge + nDeltaCharge)
                    {
                        nNumDeltaCharge--;
                    }

                    nNumDeltaH += abs( nDeltaH );
                    /* nNumDeltaCharge += abs(nDeltaCharge); */
                    vf[1].bUsed = VF_USED_ALL;
                }
            }

            for (j = 0; j < 3; j++)
            {
                nNumDeltaCharge += GetDeltaChargeFromVF( pBNS, pVA, &vf[j] );
            }

            *pnDeltaH += nNumDeltaH;
            *pnDeltaCharge += nNumDeltaCharge;
            *pnNumVisitedAtoms += nNumVisitedAtoms;


            if (v2 != vLast)
            {
                err = BNS_PROGRAM_ERR;
            }
        }
        return err ? err : ret;
    }
    */
    // END INCHI C FUNCTION: EvaluateChargeChanges
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: EvaluateChargeChanges
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: #define ANY_VERT_TYPE(X) (((X) & (BNS_VERT_TYPE_ATOM | BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP)) && !((X) & BNS_VERT_TYPE_SUPER_TGROUP))
    // INCHI✔️❌: #define ALTP_DELTA(altp) (altp)[iALTP_FLOW].flow[0]
    // INCHI✔️❌: #define ALTP_PATH_LEN(altp) (altp)[iALTP_PATH_LEN].number
    // INCHI✔️❌: #define ALTP_START_ATOM(altp) (altp)[iALTP_START_ATOM].number
    // INCHI✔️❌: #define ALTP_END_ATOM(altp) (altp)[iALTP_END_ATOM].number
    // INCHI✔️❌: #define ALTP_THIS_ATOM_NEIGHBOR(altp,i) (altp)[iALTP_NEIGHBOR+(i)].ineigh[0]
    // INCHI✔️❌: SourceHeap checked allocation lookup adds overhead versus direct C pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: EvaluateChargeChanges

    fn any_vertex_type(type_: i32) -> bool {
        type_ & (BNS_VERT_TYPE_ATOM | BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP) as i32 != 0
            && type_ & BNS_VERT_TYPE_SUPER_TGROUP as i32 == 0
    }
    fn path_entry(
        heap: &SourceHeap,
        path: SourceMutPointer<BNS_ALT_PATH>,
        index: i32,
    ) -> Result<BNS_ALT_PATH, SourceHeapError> {
        heap.slice(path.as_const())?
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn vertex(heap: &SourceHeap, network: &BN_STRUCT, index: i32) -> Result<BNS_VERTEX, SourceHeapError> {
        heap.slice(network.vert.as_const())?
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn edge(heap: &SourceHeap, network: &BN_STRUCT, index: i32) -> Result<BNS_EDGE, SourceHeapError> {
        heap.slice(network.edge.as_const())?
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    *pnDeltaH = 0;
    *pnDeltaCharge = 0;
    *pnNumVisitedAtoms = 0;

    let mut err = 0_i32;
    let mut pass = pBNS.num_altp.wrapping_sub(1);
    while pass >= 0 {
        let pass_index = usize::try_from(pass).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        pBNS.alt_path = *pBNS.altp.get(pass_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let path = pBNS.alt_path;
        let mut v1 = path_entry(heap, path, tagAltPathConst_iALTP_START_ATOM as i32)?.number();
        let path_len = path_entry(heap, path, tagAltPathConst_iALTP_PATH_LEN as i32)?.number();
        let mut delta = path_entry(heap, path, tagAltPathConst_iALTP_FLOW as i32)?.flow(0);
        let v_last = path_entry(heap, path, tagAltPathConst_iALTP_END_ATOM as i32)?.number();
        let mut v0 = NO_VERTEX;
        let mut v2 = NO_VERTEX;

        let mut vf: [VF; NUM_VF as usize + 1] = std::array::from_fn(|_| VF {
            v: NO_VERTEX,
            e_In: NO_VERTEX,
            e_Out: NO_VERTEX,
            ..VF::default()
        });
        let mut i_prev = 0_i32;
        let first_type = i32::from(vertex(heap, pBNS, v1)?.type_);
        let mut number_visited_atoms = 0_i32;
        if any_vertex_type(first_type) {
            if first_type & BNS_VERT_TYPE_ATOM as i32 != 0 {
                number_visited_atoms = number_visited_atoms.wrapping_add(1);
            }
            vf[2].type_ = first_type;
            vf[2].v = v1;
            i_prev = 2;
        }

        let mut number_delta_h = 0_i32;
        let mut number_delta_charge = 0_i32;
        number_visited_atoms = 0;

        let mut i = 0_i32;
        while i < path_len {
            let ineigh1 =
                i32::from(path_entry(heap, path, (tagAltPathConst_iALTP_HDR_LEN as i32).wrapping_add(i))?.ineigh(0));
            let current_vertex = vertex(heap, pBNS, v1)?;
            let edge_index = *heap
                .slice(current_vertex.iedge.as_const())?
                .get(usize::try_from(ineigh1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let current_edge = edge(heap, pBNS, edge_index)?;

            if i_prev != 0 {
                vf[2].e_Out = edge_index;
                vf[2].delta_Out = delta;
            }

            v2 = i32::from(current_edge.neighbor12) ^ v1;
            let next_type = i32::from(vertex(heap, pBNS, v2)?.type_);
            if next_type & BNS_VERT_TYPE_ATOM as i32 != 0 {
                number_visited_atoms = number_visited_atoms.wrapping_add(1);
            }

            if (any_vertex_type(next_type) || i == path_len.wrapping_sub(1))
                && vf[0].type_ & BNS_VERT_TYPE_C_GROUP as i32 != 0
                && vf[0].bUsed != VF_USED_ALL as i32
            {
                number_delta_charge =
                    number_delta_charge.wrapping_add(GetDeltaChargeFromVF(heap, pBNS, pVA, &mut vf[0])?);
            }

            if any_vertex_type(next_type) {
                vf[0] = vf[1].clone();
                vf[1] = vf[2].clone();
                vf[2] = vf[3].clone();
                vf[2].v = v2;
                vf[2].type_ = next_type;
                vf[2].e_In = edge_index;
                vf[2].delta_In = delta;
                i_prev = 2;
            } else if i == path_len.wrapping_sub(1) {
                vf[0] = vf[1].clone();
                vf[1] = vf[2].clone();
                vf[2] = vf[3].clone();
                i_prev = 1;
            } else {
                i_prev = 0;
            }

            if i_prev != 0 && vf[1].type_ & BNS_VERT_TYPE_ATOM as i32 != 0 {
                let atom_index = usize::try_from(vf[1].v).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = pVA.get(atom_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = atom.nCPlusGroupEdge.wrapping_sub(1);
                let minus_edge = atom.nCMinusGroupEdge.wrapping_sub(1);
                let plus_flow = if plus_edge >= 0 {
                    let charge_edge = edge(heap, pBNS, plus_edge)?;
                    charge_edge.cap.wrapping_sub(charge_edge.flow)
                } else {
                    0
                };
                let minus_flow = if minus_edge >= 0 {
                    edge(heap, pBNS, minus_edge)?.flow.wrapping_neg()
                } else {
                    0
                };
                let initial_charge = i32::from(atom.cInitCharge)
                    .wrapping_add(plus_flow)
                    .wrapping_add(minus_flow);

                let mut delta_h = 0_i32;
                let mut delta_charge = 0_i32;
                if vf[0].type_ & BNS_VERT_TYPE_TGROUP as i32 != 0 {
                    delta_h = delta_h.wrapping_sub(delta);
                } else if vf[0].type_ & BNS_VERT_TYPE_C_GROUP as i32 != 0
                    && vf[0].bUsed & VF_USED_OUT as i32 == 0
                    && (vf[0].e_Out == plus_edge || vf[0].e_Out == minus_edge)
                {
                    delta_charge = delta_charge.wrapping_sub(vf[0].delta_Out);
                    vf[0].bUsed |= VF_USED_OUT as i32;
                }

                if vf[2].type_ & BNS_VERT_TYPE_TGROUP as i32 != 0 {
                    delta_h = delta_h.wrapping_add(delta);
                } else if vf[2].type_ & BNS_VERT_TYPE_C_GROUP as i32 != 0
                    && vf[2].bUsed & VF_USED_IN as i32 == 0
                    && (vf[2].e_In == plus_edge || vf[2].e_In == minus_edge)
                {
                    delta_charge = delta_charge.wrapping_sub(vf[2].delta_In);
                    vf[2].bUsed |= VF_USED_IN as i32;
                }

                if initial_charge == 0 && delta_charge != 0 {
                    number_delta_charge = number_delta_charge.wrapping_add(1);
                } else if initial_charge != 0 && initial_charge.wrapping_add(delta_charge) == 0 {
                    number_delta_charge = number_delta_charge.wrapping_sub(1);
                }
                number_delta_h = number_delta_h.wrapping_add(delta_h.wrapping_abs());
                vf[1].bUsed = VF_USED_ALL as i32;
            }

            i = i.wrapping_add(1);
            delta = delta.wrapping_neg();
            v0 = v1;
            v1 = v2;
        }

        for queued in vf.iter_mut().take(NUM_VF as usize) {
            number_delta_charge = number_delta_charge.wrapping_add(GetDeltaChargeFromVF(heap, pBNS, pVA, queued)?);
        }
        *pnDeltaH = pnDeltaH.wrapping_add(number_delta_h);
        *pnDeltaCharge = pnDeltaCharge.wrapping_add(number_delta_charge);
        *pnNumVisitedAtoms = pnNumVisitedAtoms.wrapping_add(number_visited_atoms);

        if v2 != v_last {
            err = BNS_PROGRAM_ERR;
        }
        pass = pass.wrapping_sub(1);
        let _ = v0;
    }
    Ok(if err != 0 { err } else { 0 })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RunBnsTestOnce(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pVA: &[VAL_AT],
    pvFirst: &mut i32,
    pvLast: &mut i32,
    pPathLen: &mut i32,
    pnDeltaH: &mut i32,
    pnDeltaCharge: &mut i32,
    pnNumVisitedAtoms: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4507 RunBnsTestOnce
    /*
    int RunBnsTestOnce( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, Vertex *pvFirst, Vertex *pvLast,
                        int *pPathLen, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms )
    {
        int bChangeFlow = 0; /* do not change flow */
        int delta, ret, ret2, pass; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

        ReInitBnStructAltPaths( pBNS );
        pass = 0;
        pBNS->alt_path = pBNS->altp[pass];
        pBNS->num_altp = 0;
        pBNS->bChangeFlow = 0;
        delta = BalancedNetworkSearch( pBNS, pBD, bChangeFlow );
        if (delta > 0)
        {
            pBNS->alt_path = pBNS->altp[pass];
            *pvFirst = ALTP_START_ATOM( pBNS->alt_path );
            *pPathLen = ALTP_PATH_LEN( pBNS->alt_path );
            *pvLast = ALTP_END_ATOM( pBNS->alt_path );
            pBNS->num_altp++;
            ret2 = EvaluateChargeChanges( pBNS, pVA, pnDeltaH, pnDeltaCharge, pnNumVisitedAtoms ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        }
        else
        {
            *pvFirst = NO_VERTEX;
            *pPathLen = 0;
            *pvLast = NO_VERTEX;
            ret2 = 0;
        }

        ReInitBnStructAltPaths( pBNS );

        ret = ReInitBnData( pBD );

        return ( delta >= 0 && ret > 0 ) ? -ret
            : delta;
    }
    */
    // END INCHI C FUNCTION: RunBnsTestOnce
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RunBnsTestOnce
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: #define ALTP_START_ATOM(altp) (altp)[iALTP_START_ATOM].number
    // INCHI✔️❌: #define ALTP_PATH_LEN(altp) (altp)[iALTP_PATH_LEN].number
    // INCHI✔️❌: #define ALTP_END_ATOM(altp) (altp)[iALTP_END_ATOM].number
    // INCHI✔️❌: SourceHeap checked allocation lookup adds overhead versus direct C pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: RunBnsTestOnce

    let _ = ReInitBnStructAltPaths(heap, pBNS)?;
    let path = pBNS.altp[0];
    pBNS.alt_path = path;
    pBNS.num_altp = 0;
    pBNS.bChangeFlow = 0;
    let delta = BalancedNetworkSearch(heap, pBNS, pBD, 0)?;
    if delta > 0 {
        pBNS.alt_path = path;
        let path_entries = heap.slice(path.as_const())?;
        *pvFirst = path_entries
            .get(tagAltPathConst_iALTP_START_ATOM as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .number();
        *pPathLen = path_entries
            .get(tagAltPathConst_iALTP_PATH_LEN as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .number();
        *pvLast = path_entries
            .get(tagAltPathConst_iALTP_END_ATOM as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .number();
        pBNS.num_altp = pBNS.num_altp.wrapping_add(1);
        let _ret2 = EvaluateChargeChanges(heap, pBNS, pVA, pnDeltaH, pnDeltaCharge, pnNumVisitedAtoms)?;
    } else {
        *pvFirst = NO_VERTEX;
        *pPathLen = 0;
        *pvLast = NO_VERTEX;
    }

    let _ = ReInitBnStructAltPaths(heap, pBNS)?;
    let ret = ReInitBnData(heap, Some(pBD))?;
    Ok(if delta >= 0 && ret > 0 {
        ret.wrapping_neg()
    } else {
        delta
    })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RunBnsRestoreOnce(
    heap: &mut SourceHeap,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    _pVA: &mut [crate::source_types::VAL_AT],
    _pTCGroups: &mut ALL_TC_GROUPS,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4546 RunBnsRestoreOnce
    /*
    int RunBnsRestoreOnce( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups )
    {
        /* run BNS for the first time */
        int nTotalDelta = 0, ret = 0;
        int nDelta;

        ReInitBnStructAltPaths( pBNS );

        do
        {
            nDelta = RunBalancedNetworkSearch( pBNS, pBD, BNS_EF_CHNG_FLOW );
            if (IS_BNS_ERROR( nDelta ))
            {
                ret = nDelta;
                goto exit_function;
            }
            nTotalDelta += nDelta;
            ReInitBnStructAltPaths( pBNS );
            ret = ReInitBnData( pBD );
            if (ret > 0)
            {
                ret = -ret;
                goto exit_function;
            }
        }
        while (nDelta > 0 && ret == 0);

        pBNS->tot_st_flow += 2 * nTotalDelta;

        ret = nTotalDelta;

    exit_function:
        return ret;
    }
    */
    // END INCHI C FUNCTION: RunBnsRestoreOnce
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RunBnsRestoreOnce
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: #define IS_BNS_ERROR(X) (BNS_ERR <= (X) && (X) <= BNS_MAX_ERR_VALUE)
    // INCHI✔️❌: The control flow adds no extra asymptotic work, but its SourceHeap-backed
    // INCHI✔️❌: callees retain their documented snapshot and allocation overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: RunBnsRestoreOnce

    let mut total_delta = 0_i32;
    let _ = ReInitBnStructAltPaths(heap, pBNS)?;
    // Some source-level error fixtures intentionally omit BNS work arrays.
    // Keep the checked public path for those inputs so its original error code
    // is preserved; valid BNS instances use the stable workspace below.
    let mut workspace = BnsSearchWorkspace::new(heap, pBNS, pBD).ok();

    loop {
        let delta_result = if let Some(workspace) = workspace.as_mut() {
            run_balanced_network_search_with_workspace(
                heap,
                pBNS,
                pBD,
                BNS_EF_CHNG_FLOW as i32,
                clock_result,
                workspace,
            )
        } else {
            RunBalancedNetworkSearch(heap, pBNS, pBD, BNS_EF_CHNG_FLOW as i32, clock_result)
        };
        let delta = delta_result?;
        if BNS_ERR <= delta && delta <= BNS_MAX_ERR_VALUE {
            return Ok(delta);
        }
        total_delta = total_delta.wrapping_add(delta);
        let _ = ReInitBnStructAltPaths(heap, pBNS)?;
        let ret = ReInitBnData(heap, Some(pBD))?;
        if ret > 0 {
            return Ok(ret.wrapping_neg());
        }
        if delta <= 0 || ret != 0 {
            break;
        }
    }

    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(total_delta));
    Ok(total_delta)
}

pub(crate) fn comp_cc_cand(a1: &CC_CAND, a2: &CC_CAND) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4584 comp_cc_cand
    // INCHI✔️✔️: int comp_cc_cand( const void *a1, const void *a2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     const CC_CAND *p1 = (const CC_CAND *) a1;
    // INCHI✔️✔️:     const CC_CAND *p2 = (const CC_CAND *) a2;
    // INCHI✔️✔️:     int            ret;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if ((ret = (int) p2->cMetal - (int) p1->cMetal)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return ret; /* metal first */
    // INCHI✔️✔️:     if ((ret = (int) p2->cNumBondsToMetal - (int) p1->cNumBondsToMetal)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return ret; /* connected to metal first */
    // INCHI✔️✔️:     if ((ret = (int) p2->cPeriodicRowNumber - (int) p1->cPeriodicRowNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return ret; /* heaviest first */
    // INCHI✔️✔️:     if ((ret = (int) p2->num_bonds - (int) p1->num_bonds)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return ret; /* more bonds first */
    // INCHI✔️✔️:     if ((ret = (int) p1->chem_valence - (int) p2->chem_valence)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return ret; /* less bond order first */
    // INCHI✔️✔️:     if ((!p1->cNumValenceElectrons && p2->cNumValenceElectrons)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return -1; /* no valence electrons first */
    // INCHI✔️✔️:     if ((!p2->cNumValenceElectrons && p1->cNumValenceElectrons)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return -1; /* no valence electrons first */
    // INCHI✔️✔️:     if (((int) p2->cNumValenceElectrons - (int) p1->cNumValenceElectrons)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         return ret; /* more valence electrons first */
    // INCHI✔️✔️:     ret = (int) p2->iat - (int) p1->iat; /* greater canon number first */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: comp_cc_cand
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: comp_cc_cand
    // INCHI✔️✔️: typedef struct tagChargeChangeCandidate {
    // INCHI✔️✔️:     Vertex iat;
    // INCHI✔️✔️:     char   num_bonds;
    // INCHI✔️✔️:     char   chem_valence;
    // INCHI✔️✔️:     char   cMetal;
    // INCHI✔️✔️:     char   cNumBondsToMetal;
    // INCHI✔️✔️:     char   cNumValenceElectrons;
    // INCHI✔️✔️:     char   cPeriodicRowNumber;
    // INCHI✔️✔️:     char   cNumChargeStates;
    // INCHI✔️✔️:     U_CHAR el_number;
    // INCHI✔️✔️: } CC_CAND;
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux signed-char ABI.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: comp_cc_cand

    let mut ret = i32::from(a2.cMetal).wrapping_sub(i32::from(a1.cMetal));
    if ret != 0 {
        return ret;
    }
    ret = i32::from(a2.cNumBondsToMetal).wrapping_sub(i32::from(a1.cNumBondsToMetal));
    if ret != 0 {
        return ret;
    }
    ret = i32::from(a2.cPeriodicRowNumber).wrapping_sub(i32::from(a1.cPeriodicRowNumber));
    if ret != 0 {
        return ret;
    }
    ret = i32::from(a2.num_bonds).wrapping_sub(i32::from(a1.num_bonds));
    if ret != 0 {
        return ret;
    }
    ret = i32::from(a1.chem_valence).wrapping_sub(i32::from(a2.chem_valence));
    if ret != 0 {
        return ret;
    }
    if a1.cNumValenceElectrons == 0 && a2.cNumValenceElectrons != 0 {
        return -1;
    }
    if a2.cNumValenceElectrons == 0 && a1.cNumValenceElectrons != 0 {
        return -1;
    }
    if i32::from(a2.cNumValenceElectrons).wrapping_sub(i32::from(a1.cNumValenceElectrons)) != 0 {
        return ret;
    }
    i32::from(a2.iat).wrapping_sub(i32::from(a1.iat))
}

#[allow(non_snake_case)]
pub(crate) fn get_pVA_atom_type(
    pVA: &[VAL_AT],
    at: &[inp_ATOM],
    iat: i32,
    bond_type: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4628 get_pVA_atom_type
    /*
    int get_pVA_atom_type( VAL_AT *pVA, inp_ATOM *at, int iat, int bond_type )
    {
        int type = 0, val;
        if (pVA[iat].cNumValenceElectrons == 4)
        {
            if (pVA[iat].cPeriodicRowNumber == 1)
            {
                type |= EL_TYPE_C;
            }
        }

        else if (pVA[iat].cNumValenceElectrons == 6)
        {
            if (pVA[iat].cPeriodicRowNumber == 1)
            {
                type |= EL_TYPE_O;
            }
            else if (pVA[iat].cPeriodicRowNumber < 5)
            {
                type |= EL_TYPE_S;
            }
            if (bond_type == BOND_TYPE_SINGLE &&
                ( type & ( EL_TYPE_O | EL_TYPE_S ) ) &&
                 1 == nNoMetalBondsValence( at, iat ) &&
                 1 == nNoMetalNumBonds( at, iat ))
            {
                type |= EL_TYPE_OSt;
            }
        }

        else if (pVA[iat].cNumValenceElectrons == 5)
        {
            if (pVA[iat].cPeriodicRowNumber == 1)
            {
                type |= EL_TYPE_N;
            }
            else
            {
                type |= EL_TYPE_P;
            }
        }

        else if (!is_el_a_metal( pVA[iat].cPeriodicNumber ))
        {
            type |= EL_TYPE_X;
        }

        /* check for possibility to be a tautomeric endpoint (that is, be a Mobile H site) */
        val = get_endpoint_valence( at[iat].el_number );

        if (val && val > at[iat].valence && !at[iat].radical &&
             -1 <= at[iat].charge && at[iat].charge <= 0 &&
             val == at[iat].chem_bonds_valence - at[iat].charge + at[iat].num_H)
        {
            type |= EL_TYPE_PT;
        }

        return type;
    }
    */
    // END INCHI C FUNCTION: get_pVA_atom_type
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: get_pVA_atom_type
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; S_VI_O_PLUS_METAL_FIX_BOND=1.
    // INCHI✔️✔️: EL_TYPE_O=1, S=2, N=4, P=8, C=16, X=32, OSt=256, PT=512.
    // INCHI✔️✔️: Direct indexed slices and source-backed helpers preserve linear neighbor scans.
    // END INCHI ACTIVE MACRO CONFIGURATION: get_pVA_atom_type

    let atom_index = usize::try_from(iat).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let valence_atom = pVA.get(atom_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atom = at.get(atom_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut type_ = 0_i32;
    if valence_atom.cNumValenceElectrons == 4 {
        if valence_atom.cPeriodicRowNumber == 1 {
            type_ |= EL_TYPE_C as i32;
        }
    } else if valence_atom.cNumValenceElectrons == 6 {
        if valence_atom.cPeriodicRowNumber == 1 {
            type_ |= EL_TYPE_O as i32;
        } else if valence_atom.cPeriodicRowNumber < 5 {
            type_ |= EL_TYPE_S as i32;
        }
        if bond_type == BOND_TYPE_SINGLE as i32
            && type_ & (EL_TYPE_O | EL_TYPE_S) as i32 != 0
            && n_no_metal_bonds_valence(Some(at), iat)? == 1
            && n_no_metal_num_bonds(Some(at), iat)? == 1
        {
            type_ |= EL_TYPE_OSt as i32;
        }
    } else if valence_atom.cNumValenceElectrons == 5 {
        if valence_atom.cPeriodicRowNumber == 1 {
            type_ |= EL_TYPE_N as i32;
        } else {
            type_ |= EL_TYPE_P as i32;
        }
    } else if is_el_a_metal(i32::from(valence_atom.cPeriodicNumber))? == 0 {
        type_ |= EL_TYPE_X as i32;
    }

    let endpoint_valence = get_endpoint_valence(atom.el_number);
    if endpoint_valence != 0
        && endpoint_valence > i32::from(atom.valence)
        && atom.radical == 0
        && atom.charge >= -1
        && atom.charge <= 0
        && endpoint_valence
            == i32::from(atom.chem_bonds_valence)
                .wrapping_sub(i32::from(atom.charge))
                .wrapping_add(i32::from(atom.num_H))
    {
        type_ |= EL_TYPE_PT as i32;
    }
    Ok(type_)
}

#[allow(non_snake_case)]
pub(crate) fn AllocEdgeList(heap: &mut SourceHeap, pEdges: &mut EDGE_LIST, nLen: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4690 AllocEdgeList
    /*
    int AllocEdgeList( EDGE_LIST *pEdges, int nLen )
    {
        switch (nLen)
        {
            case EDGE_LIST_FREE:
                if (NULL != pEdges->pnEdges)
                {
                    inchi_free( pEdges->pnEdges );
                }
                /* fall through */
            case EDGE_LIST_CLEAR:
                memset( pEdges, 0, sizeof( *pEdges ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                break;
            default:
                if (nLen > 0 && nLen != pEdges->num_alloc)
                {
                    EdgeIndex *tmp_edges = pEdges->pnEdges;
                    int        tmp_num = pEdges->num_edges;
                    pEdges->pnEdges = (EdgeIndex *) inchi_calloc( nLen, sizeof( pEdges->pnEdges[0] ) );
                    if (!pEdges->pnEdges)
                    {
                        return RI_ERR_ALLOC;
                    }
                    tmp_num = inchi_min( tmp_num, nLen );
                    if (tmp_edges && tmp_num > 0)
                    {
                        memcpy(pEdges->pnEdges, tmp_edges, tmp_num * sizeof(pEdges->pnEdges[0]));
                        pEdges->num_edges = tmp_num;
                    }
                    else
                    {
                        pEdges->num_edges = 0;
                    }
                    if (tmp_edges)
                    {
                        inchi_free( tmp_edges );
                    }
                    pEdges->num_alloc = nLen;
                    return 0;
                }
                break;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: AllocEdgeList
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AllocEdgeList
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: EDGE_LIST_CLEAR=-1; EDGE_LIST_FREE=-2.
    // INCHI✔️❌: SourceHeap preserves ownership behavior but adds allocation-map overhead
    // INCHI✔️❌: compared with the active GCC/Linux calloc/free macros.
    // END INCHI ACTIVE MACRO CONFIGURATION: AllocEdgeList

    match nLen {
        EDGE_LIST_FREE => {
            if !pEdges.pnEdges.is_null() {
                inchi_free(heap, pEdges.pnEdges)?;
            }
            *pEdges = EDGE_LIST::default();
        }
        EDGE_LIST_CLEAR => {
            *pEdges = EDGE_LIST::default();
        }
        requested if requested > 0 && requested != pEdges.num_alloc => {
            let old_edges = pEdges.pnEdges;
            let old_count = pEdges.num_edges;
            let new_edges = match inchi_calloc::<i32>(heap, requested as u64, std::mem::size_of::<i32>() as u64) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    pEdges.pnEdges = SourceMutPointer::null();
                    return Ok(RI_ERR_ALLOC);
                }
                Err(error) => return Err(error),
            };
            pEdges.pnEdges = new_edges;
            let copied_count = old_count.min(requested);
            if !old_edges.is_null() && copied_count > 0 {
                let copied_count = usize::try_from(copied_count).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let old_values = heap
                    .slice(old_edges.as_const())?
                    .get(..copied_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                heap.slice_mut(pEdges.pnEdges)?
                    .get_mut(..copied_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(&old_values);
                pEdges.num_edges = copied_count as i32;
            } else {
                pEdges.num_edges = 0;
            }
            if !old_edges.is_null() {
                inchi_free(heap, old_edges)?;
            }
            pEdges.num_alloc = requested;
        }
        _ => {}
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn AddToEdgeList(
    heap: &mut SourceHeap,
    pEdges: &mut EDGE_LIST,
    iedge: i32,
    nAddLen: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4738 AddToEdgeList
    // INCHI✔️❌: int AddToEdgeList( EDGE_LIST *pEdges, int iedge, int nAddLen )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (pEdges->num_alloc == pEdges->num_edges)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int ret;
    // INCHI✔️❌:         if (nAddLen <= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if ((ret = AllocEdgeList( pEdges, pEdges->num_alloc + nAddLen ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     pEdges->pnEdges[pEdges->num_edges++] = (EdgeIndex) iedge;
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AddToEdgeList
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddToEdgeList
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: typedef int EdgeIndex;
    // END INCHI ACTIVE MACRO CONFIGURATION: AddToEdgeList

    if pEdges.num_alloc == pEdges.num_edges {
        if nAddLen <= 0 {
            return Ok(RI_ERR_PROGR);
        }
        let ret = AllocEdgeList(heap, pEdges, pEdges.num_alloc.wrapping_add(nAddLen))?;
        if ret != 0 {
            return Ok(ret);
        }
    }
    let index = usize::try_from(pEdges.num_edges).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(pEdges.pnEdges)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = iedge;
    pEdges.num_edges = pEdges.num_edges.wrapping_add(1);
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveFromEdgeListByIndex(
    heap: &mut SourceHeap,
    pEdges: &mut EDGE_LIST,
    index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4759 RemoveFromEdgeListByIndex
    /*
    int RemoveFromEdgeListByIndex( EDGE_LIST *pEdges, int index )
    {
        int len;
        if (0 <= ( len = pEdges->num_edges - index - 1 ))
        {
            if (len)
            {
                memmove(pEdges->pnEdges + index, pEdges->pnEdges + index + 1, len * sizeof(pEdges->pnEdges[0]));
            }
            pEdges->num_edges--;
            pEdges->pnEdges[pEdges->num_edges] = 0;
            return 0;
        }
        return -1;
    }
    */
    // END INCHI C FUNCTION: RemoveFromEdgeListByIndex
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RemoveFromEdgeListByIndex
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: EdgeIndex is int; memmove permits the overlapping left shift.
    // INCHI✔️❌: SourceHeap checked allocation lookup adds overhead versus direct C pointer access.
    // END INCHI ACTIVE MACRO CONFIGURATION: RemoveFromEdgeListByIndex

    let len = pEdges.num_edges.wrapping_sub(index).wrapping_sub(1);
    if len < 0 {
        return Ok(-1);
    }
    let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let len = usize::try_from(len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let values = heap.slice_mut(pEdges.pnEdges)?;
    if len != 0 {
        let source_start = index.checked_add(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let source_end = source_start
            .checked_add(len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if source_end > values.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        values.copy_within(source_start..source_end, index);
    }
    pEdges.num_edges = pEdges.num_edges.wrapping_sub(1);
    let old_tail = usize::try_from(pEdges.num_edges).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *values.get_mut(old_tail).ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn FindInEdgeList(heap: &SourceHeap, pEdges: &EDGE_LIST, iedge: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4777 FindInEdgeList
    /*
    int FindInEdgeList( EDGE_LIST *pEdges, int iedge )
    {
        int i;
        EdgeIndex ie = iedge;
        for (i = pEdges->num_edges - 1; 0 <= i; i--)
        {
            if (ie == pEdges->pnEdges[i])
            {
                return i;
            }
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: FindInEdgeList
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FindInEdgeList
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: EdgeIndex is int; the reverse scan returns the highest matching index.
    // INCHI✔️❌: SourceHeap checked allocation lookup adds overhead versus direct C pointer access.
    // END INCHI ACTIVE MACRO CONFIGURATION: FindInEdgeList

    let mut index = pEdges.num_edges.wrapping_sub(1);
    while index >= 0 {
        let value = *heap
            .slice(pEdges.pnEdges.as_const())?
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if value == iedge {
            return Ok(index);
        }
        index = index.wrapping_sub(1);
    }
    Ok(-1)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveFromEdgeListByValue(
    heap: &mut SourceHeap,
    pEdges: &mut EDGE_LIST,
    iedge: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4794 RemoveFromEdgeListByValue
    /*
    int RemoveFromEdgeListByValue( EDGE_LIST *pEdges, int iedge )
    {
        int i, ret, n = 0;
        EdgeIndex ie = iedge;
        for (i = pEdges->num_edges - 1; 0 <= i; i--)
        {
            if (ie == pEdges->pnEdges[i])
            {
                if ((ret = RemoveFromEdgeListByIndex( pEdges, i ))) /* djb-rwth: addressing LLVM warning */
                {
                    return ret;
                }
                n++;
            }
        }

        return n;
    }
    */
    // END INCHI C FUNCTION: RemoveFromEdgeListByValue
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RemoveFromEdgeListByValue
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: EdgeIndex is int; the reverse scan preserves duplicate-removal behavior.
    // INCHI✔️❌: SourceHeap checked allocation lookup adds overhead versus direct C pointer access.
    // END INCHI ACTIVE MACRO CONFIGURATION: RemoveFromEdgeListByValue

    let mut number_removed = 0_i32;
    let mut index = pEdges.num_edges.wrapping_sub(1);
    while index >= 0 {
        let value = *heap
            .slice(pEdges.pnEdges.as_const())?
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if value == iedge {
            let ret = RemoveFromEdgeListByIndex(heap, pEdges, index)?;
            if ret != 0 {
                return Ok(ret);
            }
            number_removed = number_removed.wrapping_add(1);
        }
        index = index.wrapping_sub(1);
    }
    Ok(number_removed)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveForbiddenEdgeMask(
    heap: &mut SourceHeap,
    pBNS: &BN_STRUCT,
    pEdges: &EDGE_LIST,
    forbidden_edge_mask: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4869 RemoveForbiddenEdgeMask
    /*
    void RemoveForbiddenEdgeMask( BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask )
    {
        int i, mask = ~forbidden_edge_mask;
        for (i = 0; i < pEdges->num_edges; i++)
        {
            pBNS->edge[pEdges->pnEdges[i]].forbidden &= mask;
        }
    }
    */
    // END INCHI C FUNCTION: RemoveForbiddenEdgeMask
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RemoveForbiddenEdgeMask
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: C signed-char promotion and narrowing are reproduced explicitly.
    // INCHI✔️❌: SourceHeap allocation lookups add overhead versus direct C pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: RemoveForbiddenEdgeMask

    let inverse_mask = !forbidden_edge_mask;
    let mut index = 0_i32;
    while index < pEdges.num_edges {
        let list_index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let edge_index = *heap
            .slice(pEdges.pnEdges.as_const())?
            .get(list_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let edge_index = usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let edge = heap
            .slice_mut(pBNS.edge)?
            .get_mut(edge_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        edge.forbidden = (i32::from(edge.forbidden) & inverse_mask) as i8;
        index = index.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn SetForbiddenEdgeMask(
    heap: &mut SourceHeap,
    pBNS: &BN_STRUCT,
    pEdges: &EDGE_LIST,
    forbidden_edge_mask: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4880 SetForbiddenEdgeMask
    /*
    void SetForbiddenEdgeMask( BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask )
    {
        int i;
        for (i = 0; i < pEdges->num_edges; i++)
        {
            pBNS->edge[pEdges->pnEdges[i]].forbidden |= forbidden_edge_mask;
        }
    }
    */
    // END INCHI C FUNCTION: SetForbiddenEdgeMask
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: SetForbiddenEdgeMask
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: C signed-char promotion and narrowing are reproduced explicitly.
    // INCHI✔️❌: SourceHeap allocation lookups add overhead versus direct C pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: SetForbiddenEdgeMask

    let mut index = 0_i32;
    while index < pEdges.num_edges {
        let edge_index = *heap
            .slice(pEdges.pnEdges.as_const())?
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let edge = heap
            .slice_mut(pBNS.edge)?
            .get_mut(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
        index = index.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn RemoveForbiddenBondFlowBits(
    heap: &mut SourceHeap,
    pBNS: &BN_STRUCT,
    forbidden_edge_mask_int: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4891 RemoveForbiddenBondFlowBits
    /*
    void RemoveForbiddenBondFlowBits( BN_STRUCT *pBNS, int forbidden_edge_mask_int )
    {
        BNS_EDGE   *e;
        int         i;
        int         inv_forbidden_edge_mask = ~forbidden_edge_mask_int;
        for (i = 0, e = pBNS->edge; i < pBNS->num_bonds; i++, e++)
        {
            e->forbidden &= inv_forbidden_edge_mask;
        }
    }
    */
    // END INCHI C FUNCTION: RemoveForbiddenBondFlowBits
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RemoveForbiddenBondFlowBits
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: BNS_EDGE.forbidden is signed char; integer promotion and narrowing are explicit.
    // INCHI✔️❌: SourceHeap checked allocation lookup adds overhead versus direct pointer traversal.
    // END INCHI ACTIVE MACRO CONFIGURATION: RemoveForbiddenBondFlowBits

    let inverse_mask = !forbidden_edge_mask_int;
    let mut edge_number = 0_i32;
    while edge_number < pBNS.num_bonds {
        let edge = heap
            .slice_mut(pBNS.edge)?
            .get_mut(usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        edge.forbidden = (i32::from(edge.forbidden) & inverse_mask) as i8;
        edge_number = edge_number.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn GetChargeFlowerUpperEdge(
    heap: &SourceHeap,
    pBNS: &BN_STRUCT,
    _pVA: &[VAL_AT],
    nChargeEdge: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4915 GetChargeFlowerUpperEdge
    /*
    int GetChargeFlowerUpperEdge( BN_STRUCT *pBNS, VAL_AT *pVA, int nChargeEdge )
    {
        int ret = NO_VERTEX, i, j, k, i0, i1;
        Vertex  v0, v1[3], vc, v_t, v;
        BNS_EDGE   *pe, *pe1[3], *pe_t;
        BNS_VERTEX *pv0, *pv1[3], *pv_t;

        if (nChargeEdge < 0)
        {
            goto exit_function;
        }

        pe = pBNS->edge + nChargeEdge;
        vc = pe->neighbor1; /* charge vertex */
        if (!IS_BNS_VT_C_GR( pBNS->vert[vc].type ))
        {
            vc = vc ^ pe->neighbor12;
        }

        v0 = vc ^ pe->neighbor12; /* ChargeStruct vertex ? */
        pv0 = pBNS->vert + v0;
        if (IS_BNS_VT_ATOM( pv0->type ))
        {
            goto exit_function; /* no charge flower exists */
        }

        /* 2 edges from v0 */
        for (i = j = 0; i < pv0->num_adj_edges && j < 3; i++)
        {
            pe1[j] = pBNS->edge + pv0->iedge[i];
            if (vc != ( v1[j] = pe1[j]->neighbor12 ^ v0 ) &&
                ( pv1[j] = pBNS->vert + v1[j],
                    !IS_BNS_VT_ATOM( pv1[j]->type ) && !IS_BNS_VT_C_GR( pv1[j]->type ) ))
            {
                j++;
            }
        }

        if (j != 2 || i != pv0->num_adj_edges)
        {
            goto exit_function;
        }

        if (pv1[1]->num_adj_edges == 2 &&
             pv1[0]->num_adj_edges == 3)
        {
            i0 = 1;
            i1 = 0;
        }
        else if (pv1[0]->num_adj_edges == 2 &&
             pv1[1]->num_adj_edges == 3)
        {
            i0 = 0;
            i1 = 1;
        }
        else
        {
            goto exit_function;
        }

        /* additional check: traverse edges around v1[i1] */
        pv_t = pv1[i1];
        v_t = v1[i1];
        for (i = k = 0; i < pv_t->num_adj_edges; i++)
        {
            pe_t = pBNS->edge + pv_t->iedge[i];
            v = pe_t->neighbor12 ^ v_t; /* v1[i1] neighbor */
            if (v == v0)
            {
                k += 1;
            }
            if (v == v1[i0])
            {
                k += 2;
            }
            if (IS_BNS_VT_ATOM( pBNS->vert[v].type ))
            {
                k += 4;
            }
        }
        if (k != 7)
        {
            goto exit_function;
        }

        ret = (int) ( pe1[i0] - pBNS->edge );

    exit_function:
        return ret;
    }
    */
    // END INCHI C FUNCTION: GetChargeFlowerUpperEdge
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetChargeFlowerUpperEdge
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: IS_BNS_VT_C_GR(X) is (((X) & BNS_VT_C_POS_ALL) == BNS_VERT_TYPE_C_GROUP).
    // INCHI✔️❌: IS_BNS_VT_ATOM(X) is ((X) & BNS_VERT_TYPE_ATOM).
    // INCHI✔️❌: SourceHeap checked access adds overhead versus direct C pointer traversal.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetChargeFlowerUpperEdge

    if nChargeEdge < 0 {
        return Ok(NO_VERTEX);
    }
    let edges = heap.slice(pBNS.edge.as_const())?;
    let vertices = heap.slice(pBNS.vert.as_const())?;
    let charge_edge_index = usize::try_from(nChargeEdge).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let charge_edge = edges
        .get(charge_edge_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let is_charge_group = |type_: u16| (type_ & BNS_VT_C_POS_ALL as u16) == BNS_VERT_TYPE_C_GROUP as u16;
    let is_atom = |type_: u16| type_ & BNS_VERT_TYPE_ATOM as u16 != 0;

    let mut charge_vertex = i32::from(charge_edge.neighbor1);
    if !is_charge_group(
        vertices
            .get(usize::try_from(charge_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .type_,
    ) {
        charge_vertex ^= i32::from(charge_edge.neighbor12);
    }
    let center_vertex = charge_vertex ^ i32::from(charge_edge.neighbor12);
    let center = vertices
        .get(usize::try_from(center_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if is_atom(center.type_) {
        return Ok(NO_VERTEX);
    }

    let center_adjacency = heap.slice(center.iedge.as_const())?;
    let mut edge_indices = [0_i32; 3];
    let mut neighbor_vertices = [0_i32; 3];
    let mut neighbor_count = 0_i32;
    let mut adjacency_index = 0_i32;
    while adjacency_index < i32::from(center.num_adj_edges) && neighbor_count < 3 {
        let edge_index = *center_adjacency
            .get(usize::try_from(adjacency_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let edge = edges
            .get(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let neighbor = i32::from(edge.neighbor12) ^ center_vertex;
        let neighbor_vertex = vertices
            .get(usize::try_from(neighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if charge_vertex != neighbor && !is_atom(neighbor_vertex.type_) && !is_charge_group(neighbor_vertex.type_) {
            let output_index = usize::try_from(neighbor_count).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            edge_indices[output_index] = edge_index;
            neighbor_vertices[output_index] = neighbor;
            neighbor_count = neighbor_count.wrapping_add(1);
        }
        adjacency_index = adjacency_index.wrapping_add(1);
    }
    if neighbor_count != 2 || adjacency_index != i32::from(center.num_adj_edges) {
        return Ok(NO_VERTEX);
    }

    let first_degree = vertices[usize::try_from(neighbor_vertices[0]).unwrap()].num_adj_edges;
    let second_degree = vertices[usize::try_from(neighbor_vertices[1]).unwrap()].num_adj_edges;
    let (upper_index, three_edge_index) = if second_degree == 2 && first_degree == 3 {
        (1_usize, 0_usize)
    } else if first_degree == 2 && second_degree == 3 {
        (0_usize, 1_usize)
    } else {
        return Ok(NO_VERTEX);
    };

    let three_edge_vertex = neighbor_vertices[three_edge_index];
    let three_edge_data = &vertices[usize::try_from(three_edge_vertex).unwrap()];
    let three_edge_adjacency = heap.slice(three_edge_data.iedge.as_const())?;
    let mut check = 0_i32;
    let mut index = 0_i32;
    while index < i32::from(three_edge_data.num_adj_edges) {
        let edge_index = *three_edge_adjacency
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let edge = edges
            .get(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let neighbor = i32::from(edge.neighbor12) ^ three_edge_vertex;
        if neighbor == center_vertex {
            check = check.wrapping_add(1);
        }
        if neighbor == neighbor_vertices[upper_index] {
            check = check.wrapping_add(2);
        }
        if is_atom(
            vertices
                .get(usize::try_from(neighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .type_,
        ) {
            check = check.wrapping_add(4);
        }
        index = index.wrapping_add(1);
    }
    Ok(if check == 7 {
        edge_indices[upper_index]
    } else {
        NO_VERTEX
    })
}

#[allow(non_snake_case)]
pub(crate) fn ConnectTwoVertices(
    heap: &mut SourceHeap,
    network: &BN_STRUCT,
    first_vertex_index: i32,
    second_vertex_index: i32,
    edge_index: i32,
    clear_edge: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2627 ConnectTwoVertices
    // INCHI✔️❌: int ConnectTwoVertices( BNS_VERTEX *p1, BNS_VERTEX *p2, BNS_EDGE *e, BN_STRUCT *pBNS, int bClearEdge )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ip1 = (int) ( p1 - pBNS->vert );
    // INCHI✔️❌:     int ip2 = (int) ( p2 - pBNS->vert );
    // INCHI✔️❌:     int ie = (int) ( e - pBNS->edge );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* debug: check bounds */
    // INCHI✔️❌:     if (ip1 >= pBNS->max_vertices || ip1 < 0 ||
    // INCHI✔️❌:          ip2 >= pBNS->max_vertices || ip2 < 0 ||
    // INCHI✔️❌:          ie >= pBNS->max_edges || ie < 0 ||
    // INCHI✔️❌:          ( p1->iedge - pBNS->iedge ) < 0 ||
    // INCHI✔️❌:          ( p1->iedge - pBNS->iedge ) + p1->max_adj_edges > pBNS->max_iedges ||
    // INCHI✔️❌:          ( p2->iedge - pBNS->iedge ) < 0 ||
    // INCHI✔️❌:          ( p2->iedge - pBNS->iedge ) + p2->max_adj_edges > pBNS->max_iedges ||
    // INCHI✔️❌:          p1->num_adj_edges >= p1->max_adj_edges ||
    // INCHI✔️❌:          p2->num_adj_edges >= p2->max_adj_edges)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* clear the edge */
    // INCHI✔️❌:     if (bClearEdge)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memset( e, 0, sizeof( *e ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:         if (e->neighbor1 || e->neighbor12)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return BNS_PROGRAM_ERR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* connect */
    // INCHI✔️❌:     e->neighbor1 = inchi_min( ip1, ip2 );
    // INCHI✔️❌:     e->neighbor12 = ip1 ^ ip2;
    // INCHI✔️❌:     p1->iedge[p1->num_adj_edges] = ie;
    // INCHI✔️❌:     p2->iedge[p2->num_adj_edges] = ie;
    // INCHI✔️❌:     e->neigh_ord[ip1 > ip2] = p1->num_adj_edges++;
    // INCHI✔️❌:
    // INCHI✔️❌:     e->neigh_ord[ip1 < ip2] = p2->num_adj_edges++;
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ConnectTwoVertices

    if first_vertex_index < 0
        || first_vertex_index >= network.max_vertices
        || second_vertex_index < 0
        || second_vertex_index >= network.max_vertices
        || edge_index < 0
        || edge_index >= network.max_edges
    {
        return Ok(BNS_VERT_EDGE_OVFL);
    }
    let first_index = usize::try_from(first_vertex_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let second_index = usize::try_from(second_vertex_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let edge_index_usize = usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let neighbor1 = u16::try_from(first_vertex_index.min(second_vertex_index))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let neighbor12 =
        u16::try_from(first_vertex_index ^ second_vertex_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    heap.with_two_slices_mut_and_optional_third(
        network.vert,
        network.edge,
        Some(network.iedge),
        |vertices, edges, incident_edges| {
            let incident_edges = incident_edges.expect("the incident-edge allocation was supplied");
            let first = vertices.get(first_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let second = vertices.get(second_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let first_start = first.iedge.difference(network.iedge)?;
            let second_start = second.iedge.difference(network.iedge)?;
            if first_start < 0
                || first_start + i64::from(first.max_adj_edges) > i64::from(network.max_iedges)
                || second_start < 0
                || second_start + i64::from(second.max_adj_edges) > i64::from(network.max_iedges)
                || first.num_adj_edges >= first.max_adj_edges
                || second.num_adj_edges >= second.max_adj_edges
            {
                return Ok(BNS_VERT_EDGE_OVFL);
            }

            let edge = edges.get(edge_index_usize).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if clear_edge == 0 && (edge.neighbor1 != 0 || edge.neighbor12 != 0) {
                return Ok(BNS_PROGRAM_ERR);
            }
            let first_slot = usize::try_from(first_start + i64::from(vertices[first_index].num_adj_edges))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let second_slot = usize::try_from(second_start + i64::from(vertices[second_index].num_adj_edges))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            incident_edges
                .get(first_slot)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            incident_edges
                .get(second_slot)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;

            let edge = &mut edges[edge_index_usize];
            if clear_edge != 0 {
                *edge = BNS_EDGE::default();
            }
            edge.neighbor1 = neighbor1;
            edge.neighbor12 = neighbor12;
            incident_edges[first_slot] = edge_index;
            incident_edges[second_slot] = edge_index;

            let first_order_slot = usize::from(first_vertex_index > second_vertex_index);
            edge.neigh_ord[first_order_slot] = vertices[first_index].num_adj_edges;
            vertices[first_index].num_adj_edges = vertices[first_index].num_adj_edges.wrapping_add(1);
            let second_order_slot = usize::from(first_vertex_index < second_vertex_index);
            edge.neigh_ord[second_order_slot] = vertices[second_index].num_adj_edges;
            vertices[second_index].num_adj_edges = vertices[second_index].num_adj_edges.wrapping_add(1);
            Ok(0)
        },
    )
}

#[allow(non_snake_case)]
pub(crate) fn AddTGroups2TCGBnStruct(
    heap: &mut SourceHeap,
    network: &mut BN_STRUCT,
    structure: &StrFromINChI,
    valence_atoms: &mut [crate::source_types::VAL_AT],
    groups: &mut ALL_TC_GROUPS,
    max_add_edges: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2425 AddTGroups2TCGBnStruct
    // INCHI✔️❌: int AddTGroups2TCGBnStruct( BN_STRUCT *pBNS,
    // INCHI✔️❌:                             StrFromINChI *pStruct,
    // INCHI✔️❌:                             VAL_AT *pVA,
    // INCHI✔️❌:                             ALL_TC_GROUPS *pTCGroups,
    // INCHI✔️❌:                             int nMaxAddEdges )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     inp_ATOM *at = pStruct->at;
    // INCHI✔️❌:     int       num_atoms = pStruct->num_atoms;
    // INCHI✔️❌:     int tot_st_cap, tot_st_flow;
    // INCHI✔️❌:     /* ret = ReInitBnStruct( pBNS ); */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pTCGroups->num_tgroups /* tgi && tgi->num_t_groups && tgi->t_group*/)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int         i, k, endpoint, /*centerpoint,*/ fictpoint;
    // INCHI✔️❌:         int         num_tg = pTCGroups->num_tgroups;
    // INCHI✔️❌:         int         num_edges = pBNS->num_edges;
    // INCHI✔️❌:         int         num_vertices = pBNS->num_vertices;
    // INCHI✔️❌:         BNS_VERTEX *vert_ficpoint, *vert_ficpoint_prev;  /* fictitious vertex describing t-group */
    // INCHI✔️❌:         BNS_VERTEX *vert_endpoint;
    // INCHI✔️❌:         BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric endpoint */
    // INCHI✔️❌:         int        nMaxTGroupNumber = 0;
    // INCHI✔️❌:         /*ENDPOINT_INFO eif;*/
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Debug: check overflow */
    // INCHI✔️❌:         if (num_vertices + num_tg >= pBNS->max_vertices)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (num_edges + pTCGroups->num_tgroup_edges >= pBNS->max_edges)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* find the largest t-group ID */
    // INCHI✔️❌:         for (i = 0; i < pTCGroups->num_tc_groups; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pTCGroups->pTCG[i].type & BNS_VERT_TYPE_TGROUP)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k = pTCGroups->pTCG[i].ord_num;
    // INCHI✔️❌:                 if (k <= 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_CPOINT_ERR; /* t-group does not have a number or has a wrong number */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (k > pTCGroups->num_tc_groups)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_CPOINT_ERR; /* t-group has a wrong number */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (k != nMaxTGroupNumber + 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_CPOINT_ERR; /* t-group numbers are not contiguously ascending */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 nMaxTGroupNumber = k;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break; /* t-groups are contiguous and first in the list */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (i != num_tg)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return BNS_CPOINT_ERR; /* number of t-groups is wrong */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* since t-group IDs may be not contiguous, clear all vertices that will be added.
    // INCHI✔️❌:            all-zeroes-vertex will be ignored by the BNS
    // INCHI✔️❌:         */
    // INCHI✔️❌:         memset( pBNS->vert + num_vertices, 0, nMaxTGroupNumber * sizeof( pBNS->vert[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* initialize new fictitious vertices */
    // INCHI✔️❌:         vert_ficpoint_prev = pBNS->vert + num_vertices - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:         tot_st_cap = tot_st_flow = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < num_tg; i++, vert_ficpoint_prev = vert_ficpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*
    // INCHI✔️❌:               vert_ficpoint-1 is the last vertex;
    // INCHI✔️❌:               vert_ficpoint   is the vertex that is being added
    // INCHI✔️❌:               Note: nGroupNumber are not contiguous
    // INCHI✔️❌:             */
    // INCHI✔️❌:             vert_ficpoint = pBNS->vert + num_vertices + pTCGroups->pTCG[i].ord_num - 1;
    // INCHI✔️❌:             vert_ficpoint->iedge = vert_ficpoint_prev->iedge + vert_ficpoint_prev->max_adj_edges;
    // INCHI✔️❌:             vert_ficpoint->max_adj_edges = pTCGroups->pTCG[i].num_edges + nMaxAddEdges + BNS_ADD_SUPER_TGROUP;
    // INCHI✔️❌:             vert_ficpoint->num_adj_edges = 0;
    // INCHI✔️❌:             vert_ficpoint->st_edge.flow = vert_ficpoint->st_edge.flow0 = 0;
    // INCHI✔️❌:             vert_ficpoint->st_edge.cap = vert_ficpoint->st_edge.cap0 = pTCGroups->pTCG[i].st_cap;
    // INCHI✔️❌:             tot_st_cap += pTCGroups->pTCG[i].st_cap;
    // INCHI✔️❌:             vert_ficpoint->type = pTCGroups->pTCG[i].type;
    // INCHI✔️❌:             pTCGroups->pTCG[i].nVertexNumber = (int) ( vert_ficpoint - pBNS->vert );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (endpoint = 0; endpoint < num_atoms; endpoint++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!at[endpoint].endpoint)
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             fictpoint = at[endpoint].endpoint + num_vertices - 1;
    // INCHI✔️❌:             vert_ficpoint = pBNS->vert + fictpoint;  /* t-group vertex */
    // INCHI✔️❌:             vert_endpoint = pBNS->vert + endpoint;   /* endpoint vertex */
    // INCHI✔️❌:             /* Debug: check overflow */
    // INCHI✔️❌:             if (fictpoint >= pBNS->max_vertices ||
    // INCHI✔️❌:                  num_edges >= pBNS->max_edges ||
    // INCHI✔️❌:                  vert_ficpoint->num_adj_edges >= vert_ficpoint->max_adj_edges ||
    // INCHI✔️❌:                  vert_endpoint->num_adj_edges >= vert_endpoint->max_adj_edges)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef NEVER
    // INCHI✔️❌:             /* obtain donor/acceptor info */
    // INCHI✔️❌:             if (!nGetEndpointInfo( at, endpoint, &eif ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = BNS_BOND_ERR;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             vert_endpoint->type |= BNS_VERT_TYPE_ENDPOINT;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef NEVER
    // INCHI✔️❌:             /* set capacity = 1 to the edges from the endpoint to the centerpoint(s) */
    // INCHI✔️❌:             for (k = 0; k < vert_endpoint->num_adj_edges; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int iedge = vert_endpoint->iedge[k];
    // INCHI✔️❌:                 if (!pBNS->edge[iedge].cap)
    // INCHI✔️❌:                 {
    // INCHI✔️❌: /* single bond, possibly between endpoint and centerpoint */
    // INCHI✔️❌:                     centerpoint = ( pBNS->edge[iedge].neighbor12 ^ endpoint );
    // INCHI✔️❌:                     if (centerpoint < pBNS->num_atoms &&
    // INCHI✔️❌:                          pBNS->vert[centerpoint].st_edge.cap >= 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int bond_type = ( at[endpoint].bond_type[k] & BOND_TYPE_MASK );
    // INCHI✔️❌:                         if (bond_type == BOND_TAUTOM ||
    // INCHI✔️❌:                             bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                             bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                             bond_type == BOND_SINGLE)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             pBNS->edge[iedge].cap = 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             /* create a new edge connecting endpoint to the new fictitious t-group vertex vert_ficpoint */
    // INCHI✔️❌:             edge = pBNS->edge + num_edges;
    // INCHI✔️❌:             edge->cap = vert_endpoint->st_edge.cap - vert_endpoint->st_edge.flow;
    // INCHI✔️❌:             edge->cap = inchi_min( edge->cap, MAX_TGROUP_EDGE_CAP );
    // INCHI✔️❌:             edge->cap = inchi_max( edge->cap, 0 );
    // INCHI✔️❌:             edge->flow = 0;
    // INCHI✔️❌:             edge->pass = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( RESET_EDGE_FORBIDDEN_MASK == 1 )
    // INCHI✔️❌:             edge->forbidden &= pBNS->edge_forbidden_mask;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef NEVER
    // INCHI✔️❌:             /* later include case when the charge change allows the endpoint to become tautomeric */
    // INCHI✔️❌:             /* mark endoint having moveable H atom with flow=1 */
    // INCHI✔️❌:
    // INCHI✔️❌:             /* -- old "no charges" version -- */
    // INCHI✔️❌:             /* if (at[endpoint].chem_bonds_valence == at[endpoint].valence) */
    // INCHI✔️❌:             /* -- the following line takes charges into account -- */
    // INCHI✔️❌:             if (eif.cDonor) /* means the endpoint has an H-atom to donate */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* increment edge flow */
    // INCHI✔️❌:                 edge->flow++;
    // INCHI✔️❌:                 /* increment one vertex st-flow & cap */
    // INCHI✔️❌:                 vert_ficpoint->st_edge.flow++;
    // INCHI✔️❌:                 vert_ficpoint->st_edge.cap++;
    // INCHI✔️❌:                 /* increment another vertex st-flow & cap */
    // INCHI✔️❌:                 vert_endpoint->st_edge.flow++;
    // INCHI✔️❌:                 vert_endpoint->st_edge.cap++;
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             /* connect edge to endpoint and fictpoint and increment the counters of neighbors and edges */
    // INCHI✔️❌:             ret = ConnectTwoVertices( vert_endpoint, vert_ficpoint, edge, pBNS, 0 );
    // INCHI✔️❌:             if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             num_edges++;
    // INCHI✔️❌:             edge->cap0 = edge->cap;
    // INCHI✔️❌:             edge->flow0 = edge->flow;
    // INCHI✔️❌:             pVA[endpoint].nTautGroupEdge = num_edges; /* edge index + 1 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         pBNS->num_edges = num_edges;
    // INCHI✔️❌:         pBNS->num_vertices += nMaxTGroupNumber;
    // INCHI✔️❌:         pBNS->num_t_groups = num_tg;
    // INCHI✔️❌:         pBNS->tot_st_cap += tot_st_cap;
    // INCHI✔️❌:         pBNS->tot_st_flow += tot_st_flow;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AddTGroups2TCGBnStruct
    let mut result = 0_i32;
    if groups.num_tgroups != 0 {
        let number_of_atoms = structure.num_atoms;
        let number_of_tgroups = groups.num_tgroups;
        let mut number_of_edges = network.num_edges;
        let number_of_vertices = network.num_vertices;
        if number_of_vertices.wrapping_add(number_of_tgroups) >= network.max_vertices
            || number_of_edges.wrapping_add(groups.num_tgroup_edges) >= network.max_edges
        {
            return Ok(BNS_VERT_EDGE_OVFL);
        }
        let group_count = usize::try_from(groups.num_tc_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let mut group_values = heap
            .slice(groups.pTCG.as_const())?
            .get(..group_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        let mut max_tgroup_number = 0_i32;
        let mut scanned = 0_usize;
        for group in &group_values {
            if group.type_ & BNS_VERT_TYPE_TGROUP as i32 == 0 {
                break;
            }
            let order = group.ord_num;
            if order <= 0 || order > groups.num_tc_groups || order != max_tgroup_number.wrapping_add(1) {
                return Ok(BNS_CPOINT_ERR);
            }
            max_tgroup_number = order;
            scanned += 1;
        }
        if scanned != usize::try_from(number_of_tgroups).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
            return Ok(BNS_CPOINT_ERR);
        }
        if number_of_vertices <= 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }

        let mut vertices = heap.slice(network.vert.as_const())?.to_vec();
        let clear_start = usize::try_from(number_of_vertices).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let clear_count = usize::try_from(max_tgroup_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        vertices
            .get_mut(clear_start..clear_start + clear_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(BNS_VERTEX::default());
        let mut previous_index = clear_start - 1;
        let mut total_stationary_capacity = 0_i32;
        let total_stationary_flow = 0_i32;
        for (index, group) in group_values
            .iter_mut()
            .take(usize::try_from(number_of_tgroups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?)
            .enumerate()
        {
            let vertex_index = usize::try_from(number_of_vertices.wrapping_add(group.ord_num).wrapping_sub(1))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let previous = vertices
                .get(previous_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let incident_pointer = previous.iedge.offset(i64::from(previous.max_adj_edges))?;
            let vertex = vertices
                .get_mut(vertex_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            vertex.iedge = incident_pointer;
            vertex.max_adj_edges = group
                .num_edges
                .wrapping_add(max_add_edges)
                .wrapping_add(BNS_ADD_SUPER_TGROUP as i32) as u16;
            vertex.num_adj_edges = 0;
            vertex.st_edge.flow = 0;
            vertex.st_edge.flow0 = 0;
            vertex.st_edge.cap = group.st_cap;
            vertex.st_edge.cap0 = group.st_cap;
            total_stationary_capacity = total_stationary_capacity.wrapping_add(group.st_cap);
            vertex.type_ = group.type_ as u16;
            group.nVertexNumber = vertex_index as i32;
            previous_index = vertex_index;
            let _ = index;
        }
        heap.slice_mut(network.vert)?.clone_from_slice(&vertices);
        heap.slice_mut(groups.pTCG)?
            .get_mut(..group_values.len())
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&group_values);

        let atoms = if number_of_atoms == 0 {
            Vec::new()
        } else {
            heap.slice(structure.at.as_const())?.to_vec()
        };
        let mut endpoint = 0_i32;
        while endpoint < number_of_atoms {
            let endpoint_index = usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = atoms.get(endpoint_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom.endpoint == 0 {
                endpoint = endpoint.wrapping_add(1);
                continue;
            }
            let fictitious_index = i32::from(atom.endpoint)
                .wrapping_add(number_of_vertices)
                .wrapping_sub(1);
            let fictitious_usize =
                usize::try_from(fictitious_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let edge_usize = usize::try_from(number_of_edges).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let vertex_snapshot = heap.slice(network.vert.as_const())?;
            if fictitious_index >= network.max_vertices
                || number_of_edges >= network.max_edges
                || vertex_snapshot
                    .get(fictitious_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .num_adj_edges
                    >= vertex_snapshot[fictitious_usize].max_adj_edges
                || vertex_snapshot
                    .get(endpoint_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .num_adj_edges
                    >= vertex_snapshot[endpoint_index].max_adj_edges
            {
                result = BNS_VERT_EDGE_OVFL;
                break;
            }
            {
                let endpoint_vertex = heap
                    .slice_mut(network.vert)?
                    .get_mut(endpoint_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                endpoint_vertex.type_ |= BNS_VERT_TYPE_ENDPOINT as u16;
                let capacity = endpoint_vertex
                    .st_edge
                    .cap
                    .wrapping_sub(endpoint_vertex.st_edge.flow)
                    .clamp(0, MAX_TGROUP_EDGE_CAP as i32);
                let edge = heap
                    .slice_mut(network.edge)?
                    .get_mut(edge_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                edge.cap = capacity;
                edge.flow = 0;
                edge.pass = 0;
            }
            result = ConnectTwoVertices(heap, network, endpoint, fictitious_index, number_of_edges, 0)?;
            if result >= BNS_ERR && result <= BNS_MAX_ERR_VALUE {
                break;
            }
            number_of_edges = number_of_edges.wrapping_add(1);
            let edge = heap
                .slice_mut(network.edge)?
                .get_mut(edge_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            edge.cap0 = edge.cap;
            edge.flow0 = edge.flow;
            valence_atoms
                .get_mut(endpoint_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .nTautGroupEdge = number_of_edges;
            endpoint = endpoint.wrapping_add(1);
        }
        network.num_edges = number_of_edges;
        network.num_vertices = network.num_vertices.wrapping_add(max_tgroup_number);
        network.num_t_groups = number_of_tgroups;
        network.tot_st_cap = network.tot_st_cap.wrapping_add(total_stationary_capacity);
        network.tot_st_flow = network.tot_st_flow.wrapping_add(total_stationary_flow);
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AddRadicalToMetal(
    heap: &mut SourceHeap,
    total_stationary_capacity: &mut i32,
    _total_stationary_flow: &mut i32,
    restore_mode: &SRM,
    network: &BN_STRUCT,
    groups: &ALL_TC_GROUPS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2706 AddRadicalToMetal
    // INCHI✔️✔️: int AddRadicalToMetal(int* tot_st_cap, int* tot_st_flow, ICHICONST SRM* pSrm, BN_STRUCT* pBNS, ALL_TC_GROUPS* pTCGroups)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int iG0 = pTCGroups->nGroup[TCG_MeFlower0]; /* index in pTCGroups->pTCG[] */
    // INCHI✔️✔️:     int iG1 = pTCGroups->nGroup[TCG_MeFlower1];
    // INCHI✔️✔️:     int iG2 = pTCGroups->nGroup[TCG_MeFlower2];
    // INCHI✔️✔️:     int iG3 = pTCGroups->nGroup[TCG_MeFlower3];
    // INCHI✔️✔️:     int n = (iG0 >= 0) + (iG1 >= 0) + (iG2 >= 0) + (iG3 >= 0);
    // INCHI✔️✔️:     int vG0, vG1, vG2, vG3;  /* M-vertex number */
    // INCHI✔️✔️:     BNS_VERTEX* pG0 = NULL, * pG1 = NULL, * pG2 = NULL, * pG3 = NULL;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (pTCGroups->num_metal_atoms &&
    // INCHI✔️✔️:         pSrm->bMetalAddFlower &&
    // INCHI✔️✔️:         *tot_st_cap % 2 &&
    // INCHI✔️✔️:         n == 4)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         vG0 = pTCGroups->pTCG[iG0].nVertexNumber;
    // INCHI✔️✔️:         vG1 = pTCGroups->pTCG[iG1].nVertexNumber;
    // INCHI✔️✔️:         vG2 = pTCGroups->pTCG[iG2].nVertexNumber;
    // INCHI✔️✔️:         vG3 = pTCGroups->pTCG[iG3].nVertexNumber;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         pG0 = pBNS->vert + vG0;
    // INCHI✔️✔️:         pG1 = pBNS->vert + vG1;
    // INCHI✔️✔️:         pG2 = pBNS->vert + vG2;
    // INCHI✔️✔️:         pG3 = pBNS->vert + vG3;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* add 1 unit to metal flower st_cap */
    // INCHI✔️✔️:         pG3->st_edge.cap++;
    // INCHI✔️✔️:         pG3->st_edge.cap0++;
    // INCHI✔️✔️:         (*tot_st_cap)++;
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddRadicalToMetal

    let group_indices = [
        groups.nGroup[TCG_MeFlower0 as usize],
        groups.nGroup[TCG_MeFlower1 as usize],
        groups.nGroup[TCG_MeFlower2 as usize],
        groups.nGroup[TCG_MeFlower3 as usize],
    ];
    let number_present = group_indices
        .into_iter()
        .map(|index| i32::from(index >= 0))
        .sum::<i32>();
    if groups.num_metal_atoms != 0
        && restore_mode.bMetalAddFlower != 0
        && *total_stationary_capacity % 2 != 0
        && number_present == 4
    {
        let group_values = heap.slice(groups.pTCG.as_const())?;
        let mut vertex_indices = [0_usize; 4];
        for (target, group_index) in vertex_indices.iter_mut().zip(group_indices) {
            let group_index = usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *target = usize::try_from(
                group_values
                    .get(group_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nVertexNumber,
            )
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        }
        let vertices = heap.slice_mut(network.vert)?;
        for vertex_index in vertex_indices {
            vertices.get(vertex_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        }
        let vertex = vertices
            .get_mut(vertex_indices[3])
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        vertex.st_edge.cap = vertex.st_edge.cap.wrapping_add(1);
        vertex.st_edge.cap0 = vertex.st_edge.cap0.wrapping_add(1);
        *total_stationary_capacity = total_stationary_capacity.wrapping_add(1);
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ConnectMetalFlower(
    heap: &mut SourceHeap,
    current_vertices: &mut i32,
    current_edges: &mut i32,
    total_stationary_capacity: &mut i32,
    total_stationary_flow: &mut i32,
    restore_mode: &SRM,
    network: &BN_STRUCT,
    groups: &ALL_TC_GROUPS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2743 ConnectMetalFlower
    // INCHI✔️❌: int ConnectMetalFlower( int *pcur_num_vertices, int *pcur_num_edges,
    // INCHI✔️❌:                         int *tot_st_cap, int *tot_st_flow, ICHICONST SRM *pSrm,
    // INCHI✔️❌:                         BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int iG0 = pTCGroups->nGroup[TCG_MeFlower0]; /* index in pTCGroups->pTCG[] */
    // INCHI✔️❌:     int iG1 = pTCGroups->nGroup[TCG_MeFlower1];
    // INCHI✔️❌:     int iG2 = pTCGroups->nGroup[TCG_MeFlower2];
    // INCHI✔️❌:     int iG3 = pTCGroups->nGroup[TCG_MeFlower3];
    // INCHI✔️❌:     int n = ( iG0 >= 0 ) + ( iG1 >= 0 ) + ( iG2 >= 0 ) + ( iG3 >= 0 );
    // INCHI✔️❌:     int vG0, vG1, vG2, vG3;  /* M-vertex number */
    // INCHI✔️❌:     int cur_num_edges = *pcur_num_edges;
    // INCHI✔️❌:     int cur_num_vertices = *pcur_num_vertices;
    // INCHI✔️❌:     BNS_VERTEX *pG0 = NULL, *pG1 = NULL, *pG2 = NULL, *pG3 = NULL;
    // INCHI✔️❌:     BNS_EDGE   *ea = NULL, *eb = NULL, *ed = NULL, *ex = NULL, *ey = NULL, *e;
    // INCHI✔️❌:     int         ia, ib, id, ix, iy;
    // INCHI✔️❌:     int         c, f, dc, df, ca, fa, cb, fb, cd, fd, cx, fx, cy, fy;
    // INCHI✔️❌:     int         C0, F0, C1, F1, C2, F2, C3, F3, D;
    // INCHI✔️❌:     int         ret = 0, i;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 == n)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (4 != n)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_PROGR;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     vG0 = pTCGroups->pTCG[iG0].nVertexNumber;
    // INCHI✔️❌:     vG1 = pTCGroups->pTCG[iG1].nVertexNumber;
    // INCHI✔️❌:     vG2 = pTCGroups->pTCG[iG2].nVertexNumber;
    // INCHI✔️❌:     vG3 = pTCGroups->pTCG[iG3].nVertexNumber;
    // INCHI✔️❌:
    // INCHI✔️❌:     pG0 = pBNS->vert + vG0;
    // INCHI✔️❌:     pG1 = pBNS->vert + vG1;
    // INCHI✔️❌:     pG2 = pBNS->vert + vG2;
    // INCHI✔️❌:     pG3 = pBNS->vert + vG3;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count G0 edges cap and flow (currently only atoms are connected to G0) */
    // INCHI✔️❌:     for (i = 0, c = 0, f = 0; i < pG0->num_adj_edges; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         e = pBNS->edge + pG0->iedge[i];
    // INCHI✔️❌:         c += e->cap;
    // INCHI✔️❌:         f += e->flow;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* consistency checks */
    // INCHI✔️❌:     if (!IS_BNS_VT_M_GR( pTCGroups->pTCG[iG0].type ) &&
    // INCHI✔️❌:         ( pTCGroups->pTCG[iG0].edges_cap != pG0->st_edge.cap ||
    // INCHI✔️❌:             pTCGroups->pTCG[iG0].edges_flow != pG0->st_edge.flow ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_PROGR;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pTCGroups->pTCG[iG0].edges_cap != c ||
    // INCHI✔️❌:          pTCGroups->pTCG[iG0].edges_flow != f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_PROGR;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* get new edges */
    // INCHI✔️❌:
    // INCHI✔️❌:     ea = pBNS->edge + (ia = cur_num_edges++);
    // INCHI✔️❌:     eb = pBNS->edge + (ib = cur_num_edges++);
    // INCHI✔️❌:     ed = pBNS->edge + (id = cur_num_edges++);
    // INCHI✔️❌:     ex = pBNS->edge + (ix = cur_num_edges++);
    // INCHI✔️❌:     ey = pBNS->edge + (iy = cur_num_edges++);
    // INCHI✔️❌:
    // INCHI✔️❌:     /* connect vertices with edges */
    // INCHI✔️❌:     ret = ConnectTwoVertices( pG0, pG1, eb, pBNS, 1 );
    // INCHI✔️❌:     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = ConnectTwoVertices( pG0, pG2, ea, pBNS, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = ConnectTwoVertices( pG1, pG2, ed, pBNS, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = ConnectTwoVertices( pG1, pG3, ey, pBNS, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = ConnectTwoVertices( pG2, pG3, ex, pBNS, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* calculate caps and flows */
    // INCHI✔️❌:
    // INCHI✔️❌:     dc = c % 2;
    // INCHI✔️❌:     c /= 2;
    // INCHI✔️❌:     df = f % 2;
    // INCHI✔️❌:     f /= 2;
    // INCHI✔️❌:
    // INCHI✔️❌:     D = pSrm->nMetalFlowerParam_D;
    // INCHI✔️❌:
    // INCHI✔️❌:     C0 = F0 = 2 * c + 2 * D + dc;
    // INCHI✔️❌:     C1 = F1 = c + 2 * D + dc - df;
    // INCHI✔️❌:     C2 = F2 = c + 2 * D;
    // INCHI✔️❌:     C3 = F3 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     ca = c + 2 * D;
    // INCHI✔️❌:     fa = c + D - f;
    // INCHI✔️❌:
    // INCHI✔️❌:     cb = c + 2 * D + dc;
    // INCHI✔️❌:     fb = c + D + dc - ( f + df );
    // INCHI✔️❌:
    // INCHI✔️❌:     cd = c + 2 * D;
    // INCHI✔️❌:     fd = f + D;
    // INCHI✔️❌:
    // INCHI✔️❌:     cx = cy = D;
    // INCHI✔️❌:     fx = fy = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check overflow */
    // INCHI✔️❌:     if (C0 >= EDGE_FLOW_ST_MASK || F0 >= EDGE_FLOW_ST_MASK ||
    // INCHI✔️❌:          C1 >= EDGE_FLOW_ST_MASK || F1 >= EDGE_FLOW_ST_MASK ||
    // INCHI✔️❌:          C2 >= EDGE_FLOW_ST_MASK || F2 >= EDGE_FLOW_ST_MASK ||
    // INCHI✔️❌:          C3 >= EDGE_FLOW_ST_MASK || F3 >= EDGE_FLOW_ST_MASK)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return BNS_PROGRAM_ERR; /* cannot handle too large st-cap or st-flow */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* set st caps and flows */
    // INCHI✔️❌:
    // INCHI✔️❌:     SetStCapFlow( pG0, tot_st_flow, tot_st_cap, C0, F0 );
    // INCHI✔️❌:     SetStCapFlow( pG1, tot_st_flow, tot_st_cap, C1, F1 );
    // INCHI✔️❌:     SetStCapFlow( pG2, tot_st_flow, tot_st_cap, C2, F2 );
    // INCHI✔️❌:     SetStCapFlow( pG3, tot_st_flow, tot_st_cap, C3, F3 );
    // INCHI✔️❌:
    // INCHI✔️❌:     SetEdgeCapFlow( ea, ca, fa );
    // INCHI✔️❌:     SetEdgeCapFlow( eb, cb, fb );
    // INCHI✔️❌:     SetEdgeCapFlow( ed, cd, fd );
    // INCHI✔️❌:     SetEdgeCapFlow( ex, cx, fx );
    // INCHI✔️❌:     SetEdgeCapFlow( ey, cy, fy );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     *pcur_num_edges = cur_num_edges;
    // INCHI✔️❌:     *pcur_num_vertices = cur_num_vertices;
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ConnectMetalFlower

    let current_vertex_count = *current_vertices;
    let group_indices = [
        groups.nGroup[TCG_MeFlower0 as usize],
        groups.nGroup[TCG_MeFlower1 as usize],
        groups.nGroup[TCG_MeFlower2 as usize],
        groups.nGroup[TCG_MeFlower3 as usize],
    ];
    let number_present = group_indices
        .into_iter()
        .map(|index| i32::from(index >= 0))
        .sum::<i32>();
    if number_present == 0 {
        return Ok(0);
    }
    if number_present != 4 {
        return Ok(RI_ERR_PROGR);
    }

    let group_values = heap.slice(groups.pTCG.as_const())?;
    let mut flower_groups: [TC_GROUP; 4] = std::array::from_fn(|_| TC_GROUP::default());
    for (target, group_index) in flower_groups.iter_mut().zip(group_indices) {
        *target = group_values
            .get(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
    }
    let vertex_indices: [i32; 4] = std::array::from_fn(|index| flower_groups[index].nVertexNumber);
    let vertices = heap.slice(network.vert.as_const())?;
    let g0 = vertices
        .get(usize::try_from(vertex_indices[0]).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for vertex_index in vertex_indices {
        vertices
            .get(usize::try_from(vertex_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
    }
    let incident_start = g0.iedge.difference(network.iedge)?;
    let incident_edges = heap.slice(network.iedge.as_const())?;
    let edges = heap.slice(network.edge.as_const())?;
    let mut capacity = 0_i32;
    let mut flow = 0_i32;
    for offset in 0..usize::from(g0.num_adj_edges) {
        let incident_index = usize::try_from(
            incident_start
                .checked_add(i64::try_from(offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let edge_index = *incident_edges
            .get(incident_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let edge = edges
            .get(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        capacity = capacity.wrapping_add(edge.cap);
        flow = flow.wrapping_add(edge.flow);
    }
    if flower_groups[0].type_ != BNS_VERT_TYPE_METAL_GR as i32
        && (flower_groups[0].edges_cap != g0.st_edge.cap || flower_groups[0].edges_flow != g0.st_edge.flow)
    {
        return Ok(RI_ERR_PROGR);
    }
    if flower_groups[0].edges_cap != capacity || flower_groups[0].edges_flow != flow {
        return Ok(RI_ERR_PROGR);
    }

    let mut next_edge = *current_edges;
    let edge_a = next_edge;
    next_edge = next_edge.wrapping_add(1);
    let edge_b = next_edge;
    next_edge = next_edge.wrapping_add(1);
    let edge_d = next_edge;
    next_edge = next_edge.wrapping_add(1);
    let edge_x = next_edge;
    next_edge = next_edge.wrapping_add(1);
    let edge_y = next_edge;
    next_edge = next_edge.wrapping_add(1);

    for (first, second, edge) in [
        (vertex_indices[0], vertex_indices[1], edge_b),
        (vertex_indices[0], vertex_indices[2], edge_a),
        (vertex_indices[1], vertex_indices[2], edge_d),
        (vertex_indices[1], vertex_indices[3], edge_y),
        (vertex_indices[2], vertex_indices[3], edge_x),
    ] {
        let result = ConnectTwoVertices(heap, network, first, second, edge, 1)?;
        if result >= BNS_ERR && result <= BNS_MAX_ERR_VALUE {
            return Ok(result);
        }
    }

    let capacity_remainder = capacity % 2;
    capacity /= 2;
    let flow_remainder = flow % 2;
    flow /= 2;
    let d = restore_mode.nMetalFlowerParam_D;
    let twice_d = d.wrapping_mul(2);

    let c0 = capacity
        .wrapping_mul(2)
        .wrapping_add(twice_d)
        .wrapping_add(capacity_remainder);
    let f0 = c0;
    let c1 = capacity
        .wrapping_add(twice_d)
        .wrapping_add(capacity_remainder)
        .wrapping_sub(flow_remainder);
    let f1 = c1;
    let c2 = capacity.wrapping_add(twice_d);
    let f2 = c2;
    let c3 = 0_i32;
    let f3 = 0_i32;

    let cap_a = capacity.wrapping_add(twice_d);
    let flow_a = capacity.wrapping_add(d).wrapping_sub(flow);
    let cap_b = capacity.wrapping_add(twice_d).wrapping_add(capacity_remainder);
    let flow_b = capacity
        .wrapping_add(d)
        .wrapping_add(capacity_remainder)
        .wrapping_sub(flow.wrapping_add(flow_remainder));
    let cap_d = capacity.wrapping_add(twice_d);
    let flow_d = flow.wrapping_add(d);
    let cap_x = d;
    let flow_x = 0_i32;
    let cap_y = d;
    let flow_y = 0_i32;

    let stationary_mask = EDGE_FLOW_ST_MASK as i32;
    if c0 >= stationary_mask
        || f0 >= stationary_mask
        || c1 >= stationary_mask
        || f1 >= stationary_mask
        || c2 >= stationary_mask
        || f2 >= stationary_mask
        || c3 >= stationary_mask
        || f3 >= stationary_mask
    {
        return Ok(BNS_PROGRAM_ERR);
    }

    for (vertex_index, cap, vertex_flow) in [
        (vertex_indices[0], c0, f0),
        (vertex_indices[1], c1, f1),
        (vertex_indices[2], c2, f2),
        (vertex_indices[3], c3, f3),
    ] {
        let vertex = heap
            .slice_mut(network.vert)?
            .get_mut(usize::try_from(vertex_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        SetStCapFlow(
            vertex,
            total_stationary_flow,
            total_stationary_capacity,
            cap,
            vertex_flow,
        );
    }
    for (edge_index, cap, edge_flow) in [
        (edge_a, cap_a, flow_a),
        (edge_b, cap_b, flow_b),
        (edge_d, cap_d, flow_d),
        (edge_x, cap_x, flow_x),
        (edge_y, cap_y, flow_y),
    ] {
        let edge = heap
            .slice_mut(network.edge)?
            .get_mut(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        SetEdgeCapFlow(edge, cap, edge_flow);
    }
    *current_edges = next_edge;
    *current_vertices = current_vertex_count;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn SetEdgeCapFlow(edge: &mut BNS_EDGE, edge_capacity: i32, edge_flow: i32) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2911 SetEdgeCapFlow
    // INCHI✔️✔️: void SetEdgeCapFlow( BNS_EDGE *e, int edge_cap, int edge_flow )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     e->cap = e->cap0 = edge_cap;
    // INCHI✔️✔️:     e->flow = e->flow0 = edge_flow;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SetEdgeCapFlow

    edge.cap0 = edge_capacity;
    edge.cap = edge.cap0;
    edge.flow0 = edge_flow;
    edge.flow = edge.flow0;
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AddEdgeFlow(
    edge_capacity: i32,
    edge_flow: i32,
    edge: &mut BNS_EDGE,
    source: &mut BNS_VERTEX,
    destination: &mut BNS_VERTEX,
    total_stationary_capacity: &mut i32,
    total_stationary_flow: &mut i32,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2923 AddEdgeFlow
    // INCHI✔️✔️: int AddEdgeFlow( int edge_cap, int edge_flow, BNS_EDGE *e01, BNS_VERTEX *pSrc /*src*/,
    // INCHI✔️✔️:                   BNS_VERTEX *pDst/*dest*/, int *tot_st_cap, int *tot_st_flow )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* overflow check */
    // INCHI✔️✔️:     if (e01->cap < 0 || edge_cap < 0 || (int) e01->cap + edge_cap >= EDGE_FLOW_MASK)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return BNS_PROGRAM_ERR;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (pDst->st_edge.cap < 0 || (int) pDst->st_edge.cap + edge_cap >= EDGE_FLOW_ST_MASK ||
    // INCHI✔️✔️:          pDst->st_edge.flow < 0 || (int) pDst->st_edge.flow + edge_flow >= EDGE_FLOW_ST_MASK ||
    // INCHI✔️✔️:          pSrc->st_edge.cap < 0 || pSrc->st_edge.flow < 0 ||
    // INCHI✔️✔️:          (int) pSrc->st_edge.flow + edge_flow >= EDGE_FLOW_ST_MASK)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return BNS_PROGRAM_ERR;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* add flow */
    // INCHI✔️✔️:     e01->cap += edge_cap;
    // INCHI✔️✔️:     e01->flow += edge_flow;
    // INCHI✔️✔️:     e01->cap0 = e01->cap;
    // INCHI✔️✔️:     e01->flow0 = e01->flow;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pDst->st_edge.cap += edge_cap;
    // INCHI✔️✔️:     pDst->st_edge.cap0 = pDst->st_edge.cap;
    // INCHI✔️✔️:     *tot_st_cap += edge_cap;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pDst->st_edge.flow += edge_flow;
    // INCHI✔️✔️:     pDst->st_edge.flow0 = pDst->st_edge.flow;
    // INCHI✔️✔️:     *tot_st_flow += edge_flow;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pSrc->st_edge.flow += edge_flow;
    // INCHI✔️✔️:     pSrc->st_edge.flow0 = pSrc->st_edge.flow;
    // INCHI✔️✔️:     *tot_st_flow += edge_flow;
    // INCHI✔️✔️:
    // INCHI✔️✔️: /*
    // INCHI✔️✔️:     pDst->st_edge.cap  += e01->cap;
    // INCHI✔️✔️:     pDst->st_edge.cap0  = pDst->st_edge.cap;
    // INCHI✔️✔️:     *tot_st_cap       += e01->cap;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pDst->st_edge.flow += e01->flow;
    // INCHI✔️✔️:     pDst->st_edge.flow0 = pDst->st_edge.flow;
    // INCHI✔️✔️:     *tot_st_flow      += e01->flow;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pSrc->st_edge.flow += e01->flow;
    // INCHI✔️✔️:     pSrc->st_edge.flow0 = pSrc->st_edge.flow;
    // INCHI✔️✔️:     *tot_st_flow      += e01->flow;
    // INCHI✔️✔️: */
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddEdgeFlow

    if edge.cap < 0 || edge_capacity < 0 || edge.cap.wrapping_add(edge_capacity) >= EDGE_FLOW_MASK as i32 {
        return BNS_PROGRAM_ERR;
    }
    if destination.st_edge.cap < 0
        || destination.st_edge.cap.wrapping_add(edge_capacity) >= EDGE_FLOW_ST_MASK as i32
        || destination.st_edge.flow < 0
        || destination.st_edge.flow.wrapping_add(edge_flow) >= EDGE_FLOW_ST_MASK as i32
        || source.st_edge.cap < 0
        || source.st_edge.flow < 0
        || source.st_edge.flow.wrapping_add(edge_flow) >= EDGE_FLOW_ST_MASK as i32
    {
        return BNS_PROGRAM_ERR;
    }

    edge.cap = edge.cap.wrapping_add(edge_capacity);
    edge.flow = edge.flow.wrapping_add(edge_flow);
    edge.cap0 = edge.cap;
    edge.flow0 = edge.flow;

    destination.st_edge.cap = destination.st_edge.cap.wrapping_add(edge_capacity);
    destination.st_edge.cap0 = destination.st_edge.cap;
    *total_stationary_capacity = total_stationary_capacity.wrapping_add(edge_capacity);

    destination.st_edge.flow = destination.st_edge.flow.wrapping_add(edge_flow);
    destination.st_edge.flow0 = destination.st_edge.flow;
    *total_stationary_flow = total_stationary_flow.wrapping_add(edge_flow);

    source.st_edge.flow = source.st_edge.flow.wrapping_add(edge_flow);
    source.st_edge.flow0 = source.st_edge.flow;
    *total_stationary_flow = total_stationary_flow.wrapping_add(edge_flow);
    0
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ConnectSuperCGroup(
    heap: &mut SourceHeap,
    super_group: i32,
    add_groups: &[i32],
    number_to_add: i32,
    current_vertices: &mut i32,
    current_edges: &mut i32,
    total_stationary_capacity: &mut i32,
    total_stationary_flow: &mut i32,
    network: &BN_STRUCT,
    groups: &mut ALL_TC_GROUPS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3058 ConnectSuperCGroup
    // INCHI✔️❌: int ConnectSuperCGroup( int nSuperCGroup, int nAddGroups[], int num_add,
    // INCHI✔️❌:                         int *pcur_num_vertices, int *pcur_num_edges,
    // INCHI✔️❌:                         int *tot_st_cap, int *tot_st_flow,
    // INCHI✔️❌:                         BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups )
    // INCHI✔️❌: {
    // INCHI✔️❌:     BNS_EDGE   **e0X = NULL, *e;
    // INCHI✔️❌:     BNS_VERTEX **pvX = NULL, *pv0 = NULL, *pv1 = NULL, *pv = NULL;
    // INCHI✔️❌:     int         *jX = NULL, *iX = NULL;
    // INCHI✔️❌:     int        i, j, num_groups, j0, i1, j1, iXX, ret = 0, fst = 0;
    // INCHI✔️❌:     int        cur_num_vertices = *pcur_num_vertices;
    // INCHI✔️❌:     int        cur_num_edges = *pcur_num_edges;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nSuperCGroup >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         i1 = pTCGroups->nGroup[nSuperCGroup]; /* the supergroup */
    // INCHI✔️❌:         if (i1 < 0)
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         i1 = -1;
    // INCHI✔️❌:         fst = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = num_groups = 0; i < num_add; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         iXX = pTCGroups->nGroup[nAddGroups[i]];
    // INCHI✔️❌:         num_groups += ( iXX >= 0 && iXX != i1 );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (num_groups < 1)
    // INCHI✔️❌:     {  /* Y connect only 2 or more groups; V connects even 1 group */
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     e0X = (BNS_EDGE   **) inchi_calloc( (long long)num_groups + 1, sizeof( e0X[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     pvX = (BNS_VERTEX **) inchi_calloc( (long long)num_groups + 1, sizeof( pvX[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     jX = (int         *) inchi_calloc( (long long)num_groups + 1, sizeof( jX[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     iX = (int         *) inchi_calloc( (long long)num_groups + 1, sizeof( iX[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!e0X || !pvX || !jX || !iX)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = RI_ERR_ALLOC;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* create vert_ficpoint -- central Y-connection vertex */
    // INCHI✔️❌:
    // INCHI✔️❌:     j0 = cur_num_vertices;
    // INCHI✔️❌:     pv0 = pBNS->vert + j0; /* center of the Y-connection; has number j0 */
    // INCHI✔️❌:     pv0->iedge = ( pv0 - 1 )->iedge + ( pv0 - 1 )->max_adj_edges;
    // INCHI✔️❌:     pv0->max_adj_edges = num_groups + 1 + BNS_ADD_EDGES; /* Y-connection num. edges */
    // INCHI✔️❌:     pv0->num_adj_edges = 0; /* nothing connected yet */
    // INCHI✔️❌:     pv0->type = BNS_VT_YVCONNECTOR;
    // INCHI✔️❌:     cur_num_vertices++;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (fst == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* find super c-group vertex pv1, number j1 */
    // INCHI✔️❌:         jX[0] = j1 = pTCGroups->pTCG[i1].nVertexNumber;
    // INCHI✔️❌:         iX[0] = i1;
    // INCHI✔️❌:         pvX[0] = pv1 = pBNS->vert + j1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* find other c-group vertices */
    // INCHI✔️❌:     for (i = 0, j = 1; i < num_add; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         iXX = pTCGroups->nGroup[nAddGroups[i]];
    // INCHI✔️❌:         if (( iXX >= 0 ) && ( iXX != i1 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iX[j] = iXX;
    // INCHI✔️❌:             jX[j] = pTCGroups->pTCG[iXX].nVertexNumber;
    // INCHI✔️❌:             pvX[j] = pBNS->vert + jX[j];
    // INCHI✔️❌:             j++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* grab (num_groups+1) free edges */
    // INCHI✔️❌:     for (i = fst; i <= num_groups; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         e = e0X[i] = pBNS->edge + cur_num_edges;
    // INCHI✔️❌:         pv = pvX[i];
    // INCHI✔️❌:         j = jX[i];
    // INCHI✔️❌:         iXX = iX[i];
    // INCHI✔️❌:         /* connect all to pv0 */
    // INCHI✔️❌:         ret = ConnectTwoVertices( pv0, pv, e, pBNS, 1 );
    // INCHI✔️❌:         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* from c-group to central Y-connecting vertex of from supergroup to (+/-) vertex */
    // INCHI✔️❌:             pTCGroups->pTCG[iX[i]].nForwardEdge = cur_num_edges;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* from central Y-connecting vertex to supergroup */
    // INCHI✔️❌:             pTCGroups->pTCG[iX[i]].nBackwardEdge = cur_num_edges;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         cur_num_edges++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* set flow and cap for incoming into pv0 edges */
    // INCHI✔️❌:     for (i = 1; i <= num_groups; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int nDelta = pTCGroups->pTCG[iX[i]].st_cap - pTCGroups->pTCG[iX[i]].edges_cap;
    // INCHI✔️❌:         int edge_cap = pTCGroups->pTCG[iX[i]].edges_cap + nDelta; /* added nDelta */
    // INCHI✔️❌:         int edge_flow = pTCGroups->pTCG[iX[i]].edges_cap - pTCGroups->pTCG[iX[i]].edges_flow /*-nDelta*/;
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = AddEdgeFlow( edge_cap, edge_flow,
    // INCHI✔️❌:                      e0X[i], pvX[i]/*src*/, pv0 /* dest*/, tot_st_cap, tot_st_flow );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (fst == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* set flow and cap for going out of pv0 and into pv1 edge */
    // INCHI✔️❌:         int edge_cap = pv0->st_edge.cap;
    // INCHI✔️❌:         int edge_flow = pv0->st_edge.cap - pv0->st_edge.flow;
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = AddEdgeFlow( pv0->st_edge.cap, pv0->st_edge.cap - pv0->st_edge.flow,
    // INCHI✔️❌:                      e0X[0], pv0/*src*/, pv1 /* dest*/, tot_st_cap, tot_st_flow );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         pTCGroups->pTCG[iX[0]].edges_cap += edge_cap;
    // INCHI✔️❌:         pTCGroups->pTCG[iX[0]].edges_flow += edge_flow;
    // INCHI✔️❌:         pTCGroups->pTCG[iX[0]].st_cap += edge_cap;
    // INCHI✔️❌:         pTCGroups->pTCG[iX[0]].st_flow += edge_flow;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* no supergroup => change cap to flow */
    // INCHI✔️❌:         *tot_st_cap += pv0->st_edge.flow - pv0->st_edge.cap;
    // INCHI✔️❌:         pv0->st_edge.cap += pv0->st_edge.flow - pv0->st_edge.cap;
    // INCHI✔️❌:         pv0->st_edge.cap0 = pv0->st_edge.cap;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *pcur_num_vertices = cur_num_vertices;
    // INCHI✔️❌:     *pcur_num_edges = cur_num_edges;
    // INCHI✔️❌:     ret = num_groups;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (e0X) inchi_free( e0X );
    // INCHI✔️❌:     if (pvX) inchi_free( pvX );
    // INCHI✔️❌:     if (jX) inchi_free( jX );
    // INCHI✔️❌:     if (iX) inchi_free( iX );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ConnectSuperCGroup

    let mut next_vertex = *current_vertices;
    let mut next_edge = *current_edges;
    let (super_group_index, first_connection) = if super_group >= 0 {
        let slot = usize::try_from(super_group).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let group_index = *groups.nGroup.get(slot).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if group_index < 0 {
            return Ok(0);
        }
        (group_index, 0_i32)
    } else {
        (-1_i32, 1_i32)
    };

    let mut number_of_groups = 0_i32;
    let mut index = 0_i32;
    while index < number_to_add {
        let add_slot = *add_groups
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let group_index = *groups
            .nGroup
            .get(usize::try_from(add_slot).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        number_of_groups =
            number_of_groups.wrapping_add(i32::from(group_index >= 0 && group_index != super_group_index));
        index = index.wrapping_add(1);
    }
    if number_of_groups < 1 {
        return Ok(0);
    }

    let allocation_count = u64::try_from(i64::from(number_of_groups) + 1)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let edge_pointers: SourceMutPointer<SourceMutPointer<BNS_EDGE>> =
        inchi_calloc(heap, allocation_count, 8).unwrap_or_else(|_| SourceMutPointer::null());
    let vertex_pointers: SourceMutPointer<SourceMutPointer<BNS_VERTEX>> =
        inchi_calloc(heap, allocation_count, 8).unwrap_or_else(|_| SourceMutPointer::null());
    let vertex_indices: SourceMutPointer<i32> =
        inchi_calloc(heap, allocation_count, 4).unwrap_or_else(|_| SourceMutPointer::null());
    let group_indices: SourceMutPointer<i32> =
        inchi_calloc(heap, allocation_count, 4).unwrap_or_else(|_| SourceMutPointer::null());

    if edge_pointers.is_null() || vertex_pointers.is_null() || vertex_indices.is_null() || group_indices.is_null() {
        inchi_free(heap, edge_pointers)?;
        inchi_free(heap, vertex_pointers)?;
        inchi_free(heap, vertex_indices)?;
        inchi_free(heap, group_indices)?;
        return Ok(RI_ERR_ALLOC);
    }

    let result = (|| -> Result<i32, SourceHeapError> {
        let central_index = next_vertex;
        let central_usize = usize::try_from(central_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if central_usize == 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let mut vertices = heap.slice(network.vert.as_const())?.to_vec();
        let previous = vertices
            .get(central_usize - 1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let central_iedge = previous.iedge.offset(i64::from(previous.max_adj_edges))?;
        let central = vertices
            .get_mut(central_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        central.iedge = central_iedge;
        central.max_adj_edges = number_of_groups.wrapping_add(1).wrapping_add(BNS_ADD_EDGES as i32) as u16;
        central.num_adj_edges = 0;
        central.type_ = BNS_VT_YVCONNECTOR as u16;
        heap.slice_mut(network.vert)?.clone_from_slice(&vertices);
        next_vertex = next_vertex.wrapping_add(1);

        if first_connection == 0 {
            let super_usize = usize::try_from(super_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let super_vertex = heap
                .slice(groups.pTCG.as_const())?
                .get(super_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .nVertexNumber;
            heap.slice_mut(vertex_indices)?[0] = super_vertex;
            heap.slice_mut(group_indices)?[0] = super_group_index;
            heap.slice_mut(vertex_pointers)?[0] = network.vert.offset(i64::from(super_vertex))?;
        }

        let mut add_index = 0_i32;
        let mut connection_index = 1_i32;
        while add_index < number_to_add {
            let add_slot = *add_groups
                .get(usize::try_from(add_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let group_index = *groups
                .nGroup
                .get(usize::try_from(add_slot).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if group_index >= 0 && group_index != super_group_index {
                let group_usize = usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let vertex_index = heap
                    .slice(groups.pTCG.as_const())?
                    .get(group_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nVertexNumber;
                let target = usize::try_from(connection_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                heap.slice_mut(group_indices)?[target] = group_index;
                heap.slice_mut(vertex_indices)?[target] = vertex_index;
                heap.slice_mut(vertex_pointers)?[target] = network.vert.offset(i64::from(vertex_index))?;
                connection_index = connection_index.wrapping_add(1);
            }
            add_index = add_index.wrapping_add(1);
        }

        let mut connection = first_connection;
        while connection <= number_of_groups {
            let slot = usize::try_from(connection).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let vertex_index = heap.slice(vertex_indices.as_const())?[slot];
            let group_index = heap.slice(group_indices.as_const())?[slot];
            heap.slice_mut(edge_pointers)?[slot] = network.edge.offset(i64::from(next_edge))?;
            let connect_result = ConnectTwoVertices(heap, network, central_index, vertex_index, next_edge, 1)?;
            if connect_result >= BNS_ERR && connect_result <= BNS_MAX_ERR_VALUE {
                return Ok(connect_result);
            }
            let group = heap
                .slice_mut(groups.pTCG)?
                .get_mut(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if connection != 0 {
                group.nForwardEdge = next_edge;
            } else {
                group.nBackwardEdge = next_edge;
            }
            next_edge = next_edge.wrapping_add(1);
            connection = connection.wrapping_add(1);
        }

        let mut incoming = 1_i32;
        while incoming <= number_of_groups {
            let slot = usize::try_from(incoming).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let group_index = heap.slice(group_indices.as_const())?[slot];
            let source_vertex = heap.slice(vertex_indices.as_const())?[slot];
            let group = heap
                .slice(groups.pTCG.as_const())?
                .get(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let delta = group.st_cap.wrapping_sub(group.edges_cap);
            let edge_capacity = group.edges_cap.wrapping_add(delta);
            let edge_flow = group.edges_cap.wrapping_sub(group.edges_flow);
            let edge_index = heap.slice(edge_pointers.as_const())?[slot].difference(network.edge)?;
            let edge_index = usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let source_usize = usize::try_from(source_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let central_usize = usize::try_from(central_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let add_result = heap.with_two_slices_mut_and_optional_third::<BNS_EDGE, BNS_VERTEX, i32, _>(
                network.edge,
                network.vert,
                None,
                |edges, vertices, _| {
                    let edge = edges.get_mut(edge_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let (source, destination) = if source_usize < central_usize {
                        let (left, right) = vertices.split_at_mut(central_usize);
                        (
                            left.get_mut(source_usize).ok_or(SourceHeapError::PointerOutOfBounds)?,
                            right.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                    } else {
                        let (left, right) = vertices.split_at_mut(source_usize);
                        (
                            right.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?,
                            left.get_mut(central_usize).ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                    };
                    Ok(AddEdgeFlow(
                        edge_capacity,
                        edge_flow,
                        edge,
                        source,
                        destination,
                        total_stationary_capacity,
                        total_stationary_flow,
                    ))
                },
            )?;
            if add_result >= BNS_ERR && add_result <= BNS_MAX_ERR_VALUE {
                return Ok(add_result);
            }
            incoming = incoming.wrapping_add(1);
        }

        if first_connection == 0 {
            let super_vertex = heap.slice(vertex_indices.as_const())?[0];
            let super_group = heap.slice(group_indices.as_const())?[0];
            let edge_index = usize::try_from(heap.slice(edge_pointers.as_const())?[0].difference(network.edge)?)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let central_usize = usize::try_from(central_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let super_usize = usize::try_from(super_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let central_snapshot = heap
                .slice(network.vert.as_const())?
                .get(central_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let edge_capacity = central_snapshot.st_edge.cap;
            let edge_flow = central_snapshot.st_edge.cap.wrapping_sub(central_snapshot.st_edge.flow);
            let add_result = heap.with_two_slices_mut_and_optional_third::<BNS_EDGE, BNS_VERTEX, i32, _>(
                network.edge,
                network.vert,
                None,
                |edges, vertices, _| {
                    let edge = edges.get_mut(edge_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let (source, destination) = if central_usize < super_usize {
                        let (left, right) = vertices.split_at_mut(super_usize);
                        (
                            left.get_mut(central_usize).ok_or(SourceHeapError::PointerOutOfBounds)?,
                            right.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                    } else {
                        let (left, right) = vertices.split_at_mut(central_usize);
                        (
                            right.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?,
                            left.get_mut(super_usize).ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                    };
                    Ok(AddEdgeFlow(
                        edge_capacity,
                        edge_flow,
                        edge,
                        source,
                        destination,
                        total_stationary_capacity,
                        total_stationary_flow,
                    ))
                },
            )?;
            if add_result >= BNS_ERR && add_result <= BNS_MAX_ERR_VALUE {
                return Ok(add_result);
            }
            let group = heap
                .slice_mut(groups.pTCG)?
                .get_mut(usize::try_from(super_group).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            group.edges_cap = group.edges_cap.wrapping_add(edge_capacity);
            group.edges_flow = group.edges_flow.wrapping_add(edge_flow);
            group.st_cap = group.st_cap.wrapping_add(edge_capacity);
            group.st_flow = group.st_flow.wrapping_add(edge_flow);
        } else {
            let central = heap
                .slice_mut(network.vert)?
                .get_mut(usize::try_from(central_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let difference = central.st_edge.flow.wrapping_sub(central.st_edge.cap);
            *total_stationary_capacity = total_stationary_capacity.wrapping_add(difference);
            central.st_edge.cap = central.st_edge.cap.wrapping_add(difference);
            central.st_edge.cap0 = central.st_edge.cap;
        }

        *current_vertices = next_vertex;
        *current_edges = next_edge;
        Ok(number_of_groups)
    })();

    inchi_free(heap, edge_pointers)?;
    inchi_free(heap, vertex_pointers)?;
    inchi_free(heap, vertex_indices)?;
    inchi_free(heap, group_indices)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn AddStCapFlow(
    vertex: &mut BNS_VERTEX,
    total_stationary_flow: &mut i32,
    total_stationary_capacity: &mut i32,
    capacity: i32,
    flow: i32,
) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3218 AddStCapFlow
    // INCHI✔️✔️: void AddStCapFlow( BNS_VERTEX *vert_ficpoint, int *tot_st_flow, int *tot_st_cap, int cap, int flow )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     vert_ficpoint->st_edge.flow += flow;
    // INCHI✔️✔️:     *tot_st_flow += flow;
    // INCHI✔️✔️:     vert_ficpoint->st_edge.cap += cap;
    // INCHI✔️✔️:     *tot_st_cap += cap;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     vert_ficpoint->st_edge.flow0 = vert_ficpoint->st_edge.flow;
    // INCHI✔️✔️:     vert_ficpoint->st_edge.cap0 = vert_ficpoint->st_edge.cap;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddStCapFlow

    vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(flow);
    *total_stationary_flow = total_stationary_flow.wrapping_add(flow);
    vertex.st_edge.cap = vertex.st_edge.cap.wrapping_add(capacity);
    *total_stationary_capacity = total_stationary_capacity.wrapping_add(capacity);
    vertex.st_edge.flow0 = vertex.st_edge.flow;
    vertex.st_edge.cap0 = vertex.st_edge.cap;
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AddCGroups2TCGBnStruct(
    heap: &mut SourceHeap,
    network: &mut BN_STRUCT,
    structure: &StrFromINChI,
    valence_atoms: &mut [crate::source_types::VAL_AT],
    groups: &mut ALL_TC_GROUPS,
    max_add_edges: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3249 AddCGroups2TCGBnStruct
    // INCHI✔️❌: int AddCGroups2TCGBnStruct( BN_STRUCT *pBNS, StrFromINChI *pStruct, VAL_AT *pVA,
    // INCHI✔️❌:                             ALL_TC_GROUPS *pTCGroups, int nMaxAddEdges )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, ret1, ret2, ret3, bNeedsFlower; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     inp_ATOM *at = pStruct->at;
    // INCHI✔️❌:     int       num_atoms = pStruct->num_atoms;
    // INCHI✔️❌:     /*int       num_tg           = pTCGroups->num_tgroups;*/
    // INCHI✔️❌:     int       num_cg = pTCGroups->num_tc_groups - pTCGroups->num_tgroups;
    // INCHI✔️❌:     int       fst_cg_vertex = pBNS->num_vertices;
    // INCHI✔️❌:     int       fst_cg_group = pTCGroups->num_tgroups;
    // INCHI✔️❌:     int       num_vertices = pBNS->num_vertices;
    // INCHI✔️❌:     int       num_edges = pBNS->num_edges;
    // INCHI✔️❌:     int       cg_charge = 0;
    // INCHI✔️❌:     ICHICONST SRM *pSrm = pStruct->pSrm;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* ret = ReInitBnStruct( pBNS ); */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_cg > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:     /* if ( cgi && cgi->num_c_groups && cgi->c_group ) */
    // INCHI✔️❌:         int         i, i1, i2, j, j1, j2, k, k1, k2, n, c_point, c_neigh, cap, flow;
    // INCHI✔️❌:         int         cur_num_vertices, cur_num_edges;
    // INCHI✔️❌:         BNS_VERTEX *vert_ficpoint, *vert_ficpoint_prev, *vert_ficpoint_base;  /* fictitious vertex describing charge c-group */
    // INCHI✔️❌:         BNS_VERTEX *pv1, *pv2;
    // INCHI✔️❌:         BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric c_point */
    // INCHI✔️❌:         int        nMaxCGroupNumber = 0;
    // INCHI✔️❌:         MY_CONST C_NODE *pCN;
    // INCHI✔️❌:         int              cn_len, cn_bits, bMetalAtoms;
    // INCHI✔️❌:         int              type;
    // INCHI✔️❌:         int        tot_st_cap, tot_st_flow;
    // INCHI✔️❌:         int        nAddGroups[16];
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Debug: check overflow */
    // INCHI✔️❌:         if (num_vertices >= pBNS->max_vertices)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         nMaxCGroupNumber = num_cg;
    // INCHI✔️❌:         /* clear all vertices not used until now */
    // INCHI✔️❌:         memset( pBNS->vert + num_vertices, 0, ( (long long)pBNS->max_vertices - num_vertices ) * sizeof( pBNS->vert[0] ) ); /* djb-rwth: cast operator added; memset_s C11/Annex K variant? */
    // INCHI✔️❌:         tot_st_cap = pBNS->tot_st_cap;
    // INCHI✔️❌:         tot_st_flow = pBNS->tot_st_flow;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*****************************************/
    // INCHI✔️❌:         /* initialize new fictitious vertices    */
    // INCHI✔️❌:         /* representing c-point groups, c-groups */
    // INCHI✔️❌:         /*****************************************/
    // INCHI✔️❌:         vert_ficpoint_prev = pBNS->vert + fst_cg_vertex - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < num_cg; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*
    // INCHI✔️❌:               vert_ficpoint-1 is the last vertex;
    // INCHI✔️❌:               vert_ficpoint   is the being added vertex
    // INCHI✔️❌:               Note: nGroupNumber are not contiguous
    // INCHI✔️❌:             */
    // INCHI✔️❌:             vert_ficpoint = vert_ficpoint_prev + 1;
    // INCHI✔️❌:             vert_ficpoint->iedge = vert_ficpoint_prev->iedge + vert_ficpoint_prev->max_adj_edges;
    // INCHI✔️❌:             vert_ficpoint->max_adj_edges = pTCGroups->pTCG[i + fst_cg_group].num_edges + nMaxAddEdges;
    // INCHI✔️❌:             vert_ficpoint->num_adj_edges = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:             vert_ficpoint->st_edge.flow += pTCGroups->pTCG[i + fst_cg_group].st_flow;
    // INCHI✔️❌:             tot_st_flow += pTCGroups->pTCG[i + fst_cg_group].st_flow;
    // INCHI✔️❌:             vert_ficpoint->st_edge.cap += pTCGroups->pTCG[i + fst_cg_group].st_cap;
    // INCHI✔️❌:             tot_st_cap += pTCGroups->pTCG[i + fst_cg_group].st_cap;
    // INCHI✔️❌:
    // INCHI✔️❌:             vert_ficpoint->st_edge.flow0 = vert_ficpoint->st_edge.flow;
    // INCHI✔️❌:             vert_ficpoint->st_edge.cap0 = vert_ficpoint->st_edge.cap;
    // INCHI✔️❌:
    // INCHI✔️❌:             vert_ficpoint->type = pTCGroups->pTCG[i + fst_cg_group].type;
    // INCHI✔️❌:             /* save the vertex number */
    // INCHI✔️❌:             pTCGroups->pTCG[i + fst_cg_group].nVertexNumber = (int) ( vert_ficpoint - pBNS->vert );
    // INCHI✔️❌:
    // INCHI✔️❌:             vert_ficpoint_prev = vert_ficpoint; /* keep track of iedges */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         cur_num_vertices = (int) ( vert_ficpoint_prev - pBNS->vert ) + 1;
    // INCHI✔️❌:         cur_num_edges = num_edges;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*************************************************************/
    // INCHI✔️❌:         /* pass 1:                                                   */
    // INCHI✔️❌:         /* create ChargeStruct for c-points and connect them to      */
    // INCHI✔️❌:         /* the vertices representing c-point groups;                 */
    // INCHI✔️❌:         /* set final atom st_cap, st_flow                            */
    // INCHI✔️❌:         /*************************************************************/
    // INCHI✔️❌:         for (c_point = 0; c_point < num_atoms; c_point++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!( k = pVA[c_point].cnListIndex ))
    // INCHI✔️❌:                 continue;  /* not a c-point */
    // INCHI✔️❌:
    // INCHI✔️❌:             k--;
    // INCHI✔️❌:             pCN = cnList[k].pCN;   /* pointer to the ChargeStruct */
    // INCHI✔️❌:             cn_len = cnList[k].len;   /* length of the ChargeStruct  */
    // INCHI✔️❌:             cn_bits = cnList[k].bits;  /* bits: for M-recognition */
    // INCHI✔️❌:             /* cn_bits = cnList[k].bits; */ /* ChargeStruct type */
    // INCHI✔️❌:             bMetalAtoms = ( cn_bits == cn_bits_Me );
    // INCHI✔️❌:             vert_ficpoint_base = vert_ficpoint_prev; /* add aux vertices after this */
    // INCHI✔️❌:
    // INCHI✔️❌:             /* create disconnected auxiliary vertices of the at[c_point] ChargeStruct; add to them st_flow & st_cap */
    // INCHI✔️❌:             for (i1 = 0; i1 < cn_len; i1++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!IS_BNS_VT_CHRG_STRUCT( pCN[i1].v.type ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* the atom is always the first; the attached c-points are always the last */
    // INCHI✔️❌:                 vert_ficpoint = vert_ficpoint_base + i1; /* i1 = 1, 2,.. less number of attached c-points */
    // INCHI✔️❌:                 vert_ficpoint->iedge = vert_ficpoint_prev->iedge + vert_ficpoint_prev->max_adj_edges;
    // INCHI✔️❌:                 vert_ficpoint->max_adj_edges = pCN[i1].v.valence; /* do not add additional edges to aux vertices */
    // INCHI✔️❌:                 vert_ficpoint->num_adj_edges = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                 cap = !bMetalAtoms ? pCN[i1].v.cap : pCN[i1].v.cap ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                 flow = !bMetalAtoms ? pCN[i1].v.flow : pCN[i1].v.flow ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                 AddStCapFlow( vert_ficpoint, &tot_st_flow, &tot_st_cap, cap, flow );
    // INCHI✔️❌:                 vert_ficpoint->type = pCN[i1].v.type; /* =BNS_VERT_TYPE__AUX */
    // INCHI✔️❌:
    // INCHI✔️❌:                 vert_ficpoint_prev = vert_ficpoint; /* the last one will be vert_ficpoint for the next c-point */
    // INCHI✔️❌:                 cur_num_vertices = (int) ( vert_ficpoint - pBNS->vert ) + 1;
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (vert_ficpoint->iedge + vert_ficpoint->max_adj_edges - pBNS->iedge >= pBNS->max_iedges)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (cur_num_vertices >= pBNS->max_vertices)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* connect the vertices with new edges, add edge flow and cap */
    // INCHI✔️❌:             for (i1 = 0; i1 < cn_len; i1++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pv1 = NULL;
    // INCHI✔️❌:                 k1 = -1;
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* find vertex corresponding to i1 */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (pCN[i1].v.type & BNS_VERT_TYPE_ATOM)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pv1 = pBNS->vert + c_point; /* may be only one atom -- the current c_point at i1==0 */
    // INCHI✔️❌:                     /* add atom vertex st_cap and st_flow */
    // INCHI✔️❌:                     cap = !bMetalAtoms ? pCN[i1].v.cap : pCN[i1].v.cap ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                     flow = !bMetalAtoms ? pCN[i1].v.flow : pCN[i1].v.flow ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                     AddStCapFlow( pv1, &tot_st_flow, &tot_st_cap, cap, flow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 else if (IS_BNS_VT_C_GR( pCN[i1].v.type ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* find c-group vertex by looking for its type */
    // INCHI✔️❌:                     for (j = 0; j < num_cg; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (pCN[i1].v.type == pBNS->vert[fst_cg_vertex + j].type)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             pv1 = pBNS->vert + fst_cg_vertex + j;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* index of the pTCGroups->pTCG[] */
    // INCHI✔️❌:                     if (pv1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k1 = j + fst_cg_group;
    // INCHI✔️❌:                         if (pTCGroups->pTCG[k1].type != pCN[i1].v.type ||
    // INCHI✔️❌:                              pTCGroups->pTCG[k1].ord_num)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return RI_ERR_PROGR;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 else if (IS_BNS_VT_M_GR( pCN[i1].v.type ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     k1 = pTCGroups->nGroup[TCG_MeFlower0];
    // INCHI✔️❌:                     if (k1 < 0 ||
    // INCHI✔️❌:                          pTCGroups->pTCG[k1].type != pCN[i1].v.type ||
    // INCHI✔️❌:                              pTCGroups->pTCG[k1].ord_num ||
    // INCHI✔️❌:                              !pSrm->bMetalAddFlower)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return RI_ERR_PROGR;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     pv1 = pBNS->vert + pTCGroups->pTCG[k1].nVertexNumber;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 else if (IS_BNS_VT_CHRG_STRUCT( pCN[i1].v.type ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* aux vertex */
    // INCHI✔️❌:                     pv1 = vert_ficpoint_base + i1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (!pv1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_BOND_ERR;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* connect pairs of vertices with new edges */
    // INCHI✔️❌:                 for (k = 0; k < MAX_CN_VAL && ( i2 = pCN[i1].e[k].neigh ); k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pv2 = NULL;
    // INCHI✔️❌:                     k2 = -1;
    // INCHI✔️❌:                     i2--; /* neighbor */
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* find vertex corresponding to i2 */
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (pCN[i2].v.type & BNS_VERT_TYPE_ATOM)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         pv2 = pBNS->vert + c_point;
    // INCHI✔️❌:                         cap = !bMetalAtoms ? pCN[i2].v.cap : pCN[i2].v.cap ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                         flow = !bMetalAtoms ? pCN[i2].v.flow : pCN[i2].v.flow ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                         /* add atom vertex st_cap and st_flow; this normally should not happen */
    // INCHI✔️❌:                         AddStCapFlow( pv2, &tot_st_flow, &tot_st_cap, cap, flow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     else if (IS_BNS_VT_C_GR( pCN[i2].v.type ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* find c-group vertex by looking for its type */
    // INCHI✔️❌:                         for (j = 0; j < num_cg; j++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (pCN[i2].v.type == pBNS->vert[fst_cg_vertex + j].type)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 pv2 = pBNS->vert + fst_cg_vertex + j;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (pv2)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k2 = j + fst_cg_group;
    // INCHI✔️❌:                             if (pTCGroups->pTCG[k2].type != pCN[i2].v.type ||
    // INCHI✔️❌:                                  pTCGroups->pTCG[k2].ord_num)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 return RI_ERR_PROGR;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     else if (IS_BNS_VT_M_GR( pCN[i2].v.type ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k2 = pTCGroups->nGroup[TCG_MeFlower0];
    // INCHI✔️❌:                         if (k2 < 0 ||
    // INCHI✔️❌:                              pTCGroups->pTCG[k2].type != pCN[i2].v.type ||
    // INCHI✔️❌:                                  pTCGroups->pTCG[k2].ord_num ||
    // INCHI✔️❌:                                  !pSrm->bMetalAddFlower)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return RI_ERR_PROGR;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         pv2 = pBNS->vert + pTCGroups->pTCG[k2].nVertexNumber;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     else if (IS_BNS_VT_CHRG_STRUCT( pCN[i2].v.type ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         pv2 = vert_ficpoint_base + i2;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* connect pv1 and pv2 */
    // INCHI✔️❌:                     if (!pv1 || !pv2 || pv1 == pv2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return BNS_BOND_ERR;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     j1 = (int)(pv1 - pBNS->vert);
    // INCHI✔️❌:                     j2 = (int)(pv2 - pBNS->vert);
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* create a new edge connecting pv1 and pv2 */
    // INCHI✔️❌:                     edge = pBNS->edge + cur_num_edges;
    // INCHI✔️❌:                     if ((IS_BNS_VT_M_GR( pCN[i1].v.type ) && IS_BNS_VT_ATOM( pCN[i2].v.type )) ||
    // INCHI✔️❌:                          (IS_BNS_VT_M_GR( pCN[i2].v.type ) && IS_BNS_VT_ATOM( pCN[i1].v.type ))) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* at[c_point] is a metal or is treated as a metal; connect it to M-group */
    // INCHI✔️❌:                         /* metal - M-group (i.e. Metal-Flower) edge */
    // INCHI✔️❌:                         int nStCap, nStFlow;
    // INCHI✔️❌:                         bNeedsFlower = AtomStcapStflow( at, pVA, pSrm, c_point, &nStCap, &nStFlow, &edge->cap, &edge->flow );
    // INCHI✔️❌:                         /* GetAtomToMCGroupInitEdgeCapFlow( &edge->cap, &edge->flow, pSrm, at, pVA, c_point ); */
    // INCHI✔️❌:                         if (!bNeedsFlower)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return RI_ERR_PROGR;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         pVA[c_point].nMetalGroupEdge = cur_num_edges + 1;
    // INCHI✔️❌:                         /* pBNS->vert[c_point].st_edge.cap  += edge->flow;*/ /* where was this done ???*/
    // INCHI✔️❌:                         pBNS->vert[c_point].st_edge.flow += edge->flow;
    // INCHI✔️❌:                         pBNS->vert[c_point].st_edge.cap += edge->flow + pVA[c_point].cInitFreeValences;
    // INCHI✔️❌:                         pBNS->vert[c_point].st_edge.flow0 = pBNS->vert[c_point].st_edge.flow;
    // INCHI✔️❌:                         pBNS->vert[c_point].st_edge.cap0 = pBNS->vert[c_point].st_edge.cap;
    // INCHI✔️❌:                         tot_st_flow += edge->flow;
    // INCHI✔️❌:                         tot_st_cap += edge->flow + pVA[c_point].cInitFreeValences;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         edge->cap = !bMetalAtoms ? pCN[i1].e[k].cap : pCN[i1].e[k].cap ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                         edge->flow = !bMetalAtoms ? pCN[i1].e[k].flow : pCN[i1].e[k].flow ? pSrm->nMetalMaxCharge_D : 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     edge->forbidden = pCN[i1].e[k].bForbiddenEdge ? BNS_EDGE_FORBIDDEN_MASK : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* c-group incoming edges cap and flow needed in ConnectSuperCGroup() */
    // INCHI✔️❌:                     /*
    // INCHI✔️❌:                     if ( k1 >= 0 ) {
    // INCHI✔️❌:                         pTCGroups->pTCG[k1].edges_cap  += pCN[i1].e[k].cap;
    // INCHI✔️❌:                         pTCGroups->pTCG[k1].edges_flow += pCN[i1].e[k].flow;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ( k2 >= 0 ) {
    // INCHI✔️❌:                         pTCGroups->pTCG[k2].edges_cap  += pCN[i1].e[k].cap;
    // INCHI✔️❌:                         pTCGroups->pTCG[k2].edges_flow += pCN[i1].e[k].flow;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     */
    // INCHI✔️❌:
    // INCHI✔️❌:                     edge->pass = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( RESET_EDGE_FORBIDDEN_MASK == 1 )
    // INCHI✔️❌:                     edge->forbidden &= pBNS->edge_forbidden_mask;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* check edge overflow */
    // INCHI✔️❌:                     if (pv1->num_adj_edges >= pv1->max_adj_edges ||
    // INCHI✔️❌:                          pv2->num_adj_edges >= pv2->max_adj_edges ||
    // INCHI✔️❌:                          cur_num_edges >= pBNS->max_edges)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* connect edge to the incident vertices and increment the counters of neighbors and edges */
    // INCHI✔️❌:
    // INCHI✔️❌:                     ret = ConnectTwoVertices( pv1, pv2, edge, pBNS, 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return ret;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     edge->cap0 = edge->cap;
    // INCHI✔️❌:                     edge->flow0 = edge->flow;
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* save the edge index */
    // INCHI✔️❌:                     type = IS_BNS_VT_C_GR( pv1->type ) ? pv1->type :
    // INCHI✔️❌:                         IS_BNS_VT_C_GR( pv2->type ) ? pv2->type : 0;
    // INCHI✔️❌:                     if (type)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* the edge connects to a c-group */
    // INCHI✔️❌:                         if (type & BNS_VERT_TYPE_C_NEGATIVE)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             pVA[c_point].nCMinusGroupEdge = cur_num_edges + 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             pVA[c_point].nCPlusGroupEdge = cur_num_edges + 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     cur_num_edges++; /* end of new edge creation */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /*************************************************************/
    // INCHI✔️❌:         /* pass 2:                                                   */
    // INCHI✔️❌:         /* adjust bond cap, flow from the final atom st_cap, st_flow */
    // INCHI✔️❌:         /*************************************************************/
    // INCHI✔️❌:         for (c_point = 0; c_point < num_atoms; c_point++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int st_cap, st_cap2, max_edge_flow;
    // INCHI✔️❌:             pv1 = pBNS->vert + c_point;  /* atom vertex */
    // INCHI✔️❌:             st_cap = pv1->st_edge.cap;
    // INCHI✔️❌:             for (k = 0; k < pv1->num_adj_edges; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 edge = pBNS->edge + pv1->iedge[k];      /* incident edge */
    // INCHI✔️❌:                 c_neigh = edge->neighbor12 ^ c_point;   /* adjacent vertex */
    // INCHI✔️❌:                 pv2 = pBNS->vert + c_neigh;
    // INCHI✔️❌:                 if (c_neigh > c_point || !( pv2->type & BNS_VERT_TYPE_ATOM ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* adjacent vertex is an atom; the edge is a bond; process each bond only once */
    // INCHI✔️❌:                 st_cap2 = pv2->st_edge.cap;
    // INCHI✔️❌:                 /* the edge flow <= min( incident atom st_caps) */
    // INCHI✔️❌:                 max_edge_flow = inchi_min( st_cap, st_cap2 );
    // INCHI✔️❌:                 /* bond order <= triple bond (flow=2) */
    // INCHI✔️❌:                 if (((pSrm->bMetalAddFlower && !pSrm->nMetalMinBondOrder &&
    // INCHI✔️❌:                     ( pVA[c_point].cMetal && pVA[c_point].cNumBondsToMetal)) ||
    // INCHI✔️❌:                         (pVA[c_neigh].cMetal && pVA[c_neigh].cNumBondsToMetal) )) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     max_edge_flow = inchi_min( max_edge_flow, MAX_BOND_EDGE_CAP + 1 );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     max_edge_flow = inchi_min( max_edge_flow, MAX_BOND_EDGE_CAP );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[c_point].bond_type[k] == BOND_TYPE_SINGLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* the bond has not been changed due to stereo */
    // INCHI✔️❌:                     edge->cap = edge->cap0 = max_edge_flow;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /***********************************************************/
    // INCHI✔️❌:         /**************                                 ************/
    // INCHI✔️❌:         /************** connect M-flower with new edges ************/
    // INCHI✔️❌:         /**************                                 ************/
    // INCHI✔️❌:         /***********************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = ConnectMetalFlower( &cur_num_vertices, &cur_num_edges, &tot_st_cap, &tot_st_flow, pSrm, pBNS, pTCGroups );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /***********************************************************/
    // INCHI✔️❌:         /**************                                 ************/
    // INCHI✔️❌:         /************** add additional vertices & edges ************/
    // INCHI✔️❌:         /************** to connect c-groups             ************/
    // INCHI✔️❌:         /**************                                 ************/
    // INCHI✔️❌:         /***********************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:         /* (+) supergroup, Y-connection */
    // INCHI✔️❌:
    // INCHI✔️❌:         k = 0;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Plus0;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Plus_C0;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Plus_M0;
    // INCHI✔️❌:
    // INCHI✔️❌:         ret1 = ConnectSuperCGroup( TCG_Plus,
    // INCHI✔️❌:                                    nAddGroups,
    // INCHI✔️❌:                                    k,
    // INCHI✔️❌:                                    &cur_num_vertices,
    // INCHI✔️❌:                                    &cur_num_edges,
    // INCHI✔️❌:                                    &tot_st_cap,
    // INCHI✔️❌:                                    &tot_st_flow,
    // INCHI✔️❌:                                    pBNS,
    // INCHI✔️❌:                                    pTCGroups ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* (-) supergroup, Y-connection */
    // INCHI✔️❌:
    // INCHI✔️❌:         k = 0;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Minus0;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Minus_C0;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Minus_M0;
    // INCHI✔️❌:
    // INCHI✔️❌:         ret2 = ConnectSuperCGroup( TCG_Minus,
    // INCHI✔️❌:                                    nAddGroups,
    // INCHI✔️❌:                                    k,
    // INCHI✔️❌:                                    &cur_num_vertices,
    // INCHI✔️❌:                                    &cur_num_edges,
    // INCHI✔️❌:                                    &tot_st_cap,
    // INCHI✔️❌:                                    &tot_st_flow,
    // INCHI✔️❌:                                    pBNS,
    // INCHI✔️❌:                                    pTCGroups ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:
    // INCHI✔️❌:         /******** connect (+) and (-) ***************/
    // INCHI✔️❌:
    // INCHI✔️❌:         k = 0;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Plus;
    // INCHI✔️❌:         nAddGroups[k++] = TCG_Minus;
    // INCHI✔️❌:
    // INCHI✔️❌:         ret3 = ConnectSuperCGroup( -1,
    // INCHI✔️❌:                                     nAddGroups,
    // INCHI✔️❌:                                     k,
    // INCHI✔️❌:                                     &cur_num_vertices,
    // INCHI✔️❌:                                     &cur_num_edges,
    // INCHI✔️❌:                                     &tot_st_cap,
    // INCHI✔️❌:                                     &tot_st_flow,
    // INCHI✔️❌:                                     pBNS,
    // INCHI✔️❌:                                     pTCGroups );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Take care of the full charge */
    // INCHI✔️❌:
    // INCHI✔️❌:         cg_charge = pTCGroups->total_charge - pTCGroups->tgroup_charge - pTCGroups->charge_on_atoms;
    // INCHI✔️❌:         ret = 1;
    // INCHI✔️❌:         if (ret3 > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* (+) and (-) or at least one of them have been connected */
    // INCHI✔️❌:             int nVertPlusMinus = cur_num_vertices - 1;
    // INCHI✔️❌:             BNS_VERTEX *pVertPlusMinus = pBNS->vert + nVertPlusMinus;
    // INCHI✔️❌:             BNS_VERTEX *pVertPlus = NULL, *pVertMinus = NULL, *pVert = NULL;
    // INCHI✔️❌:             BNS_EDGE   *pEdgePlus = NULL, *pEdgeMinus = NULL, *pEdge = NULL;
    // INCHI✔️❌:             n = pTCGroups->nGroup[TCG_Plus] >= 0;   /* (+)-supergroup exists */
    // INCHI✔️❌:             k = pTCGroups->nGroup[TCG_Minus] >= 0;  /* (-)-supergroup exists */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (pVertPlusMinus->num_adj_edges == 2 && k + n == 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pEdgePlus = pBNS->edge + pVertPlusMinus->iedge[0];  /* TCG_Plus was the 1st */
    // INCHI✔️❌:                 pEdgeMinus = pBNS->edge + pVertPlusMinus->iedge[1];  /* TCG_Minus was the 2nd */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             else if (pVertPlusMinus->num_adj_edges == 1 && k + n == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pTCGroups->nGroup[TCG_Plus] >= 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pEdgePlus = pBNS->edge + pVertPlusMinus->iedge[0];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (pTCGroups->nGroup[TCG_Minus] >= 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pEdgeMinus = pBNS->edge + pVertPlusMinus->iedge[0];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             else if (k + n)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* program error check */
    // INCHI✔️❌:                 ret = BNS_BOND_ERR;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (pEdgePlus)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pVertPlus = pBNS->vert + ( pEdgePlus->neighbor12 ^ nVertPlusMinus );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pEdgeMinus)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pVertMinus = pBNS->vert + ( pEdgeMinus->neighbor12 ^ nVertPlusMinus );
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             pVert = pVertPlus ? pVertPlus : pVertMinus ? pVertMinus : NULL;
    // INCHI✔️❌:             pEdge = pEdgePlus ? pEdgePlus : pEdgeMinus ? pEdgeMinus : NULL;
    // INCHI✔️❌:             if (pEdgeMinus)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pTCGroups->nEdgeMinus = (int) ( pEdgeMinus - pBNS->edge );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pEdgePlus)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pTCGroups->nEdgePlus = (int) ( pEdgePlus - pBNS->edge );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pEdge)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pTCGroups->nEdge4charge = (int) ( pEdge - pBNS->edge );
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* set total charge */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (pVert && pEdge)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* do not check rescaps for now */
    // INCHI✔️❌:                 if (cg_charge > 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pVertPlusMinus->st_edge.cap += cg_charge;
    // INCHI✔️❌:                     tot_st_cap += cg_charge;
    // INCHI✔️❌:                     pVertPlusMinus->st_edge.cap0 = pVertPlusMinus->st_edge.cap;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (cg_charge < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pVert->st_edge.cap -= cg_charge;
    // INCHI✔️❌:                     tot_st_cap -= cg_charge;
    // INCHI✔️❌:                     pVert->st_edge.cap0 = pVert->st_edge.cap;
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (pEdge->cap - pEdge->flow + cg_charge < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* 2006-02-06: increase edge capacity to avoid clogging */
    // INCHI✔️❌:                         pEdge->cap = pEdge->flow - cg_charge;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pTCGroups->added_charge = cg_charge;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!cg_charge || ( pVert && pEdge ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = 2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         AddRadicalToMetal( &tot_st_cap, &tot_st_flow, pSrm, pBNS, pTCGroups );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         pBNS->num_edges = cur_num_edges;
    // INCHI✔️❌:         pBNS->num_vertices = cur_num_vertices;
    // INCHI✔️❌:         pBNS->num_c_groups = num_cg;
    // INCHI✔️❌:         pBNS->tot_st_cap = tot_st_cap;
    // INCHI✔️❌:         pBNS->tot_st_flow = tot_st_flow;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // END INCHI C FUNCTION: AddCGroups2TCGBnStruct

    let number_of_cgroups = groups.num_tc_groups.wrapping_sub(groups.num_tgroups);
    if number_of_cgroups <= 0 {
        return Ok(0);
    }

    let first_cgroup_vertex = network.num_vertices;
    let first_cgroup = groups.num_tgroups;
    let mut current_vertices = network.num_vertices;
    let mut current_edges = network.num_edges;
    if current_vertices >= network.max_vertices {
        return Ok(BNS_VERT_EDGE_OVFL);
    }
    if current_vertices <= 0 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let number_of_atoms = structure.num_atoms;
    let atom_count = usize::try_from(number_of_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let atoms = if atom_count == 0 {
        Vec::new()
    } else {
        heap.slice(structure.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let restore_mode = heap
        .slice(structure.pSrm)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if valence_atoms.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let maximum_vertices = usize::try_from(network.max_vertices).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let clear_start = usize::try_from(first_cgroup_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut vertices = heap.slice(network.vert.as_const())?.to_vec();
    vertices
        .get_mut(clear_start..maximum_vertices)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(BNS_VERTEX::default());
    let group_count = usize::try_from(groups.num_tc_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut group_values = heap
        .slice(groups.pTCG.as_const())?
        .get(..group_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let mut total_stationary_capacity = network.tot_st_cap;
    let mut total_stationary_flow = network.tot_st_flow;
    let mut previous_vertex = clear_start - 1;
    let cgroup_count = usize::try_from(number_of_cgroups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let first_cgroup_usize = usize::try_from(first_cgroup).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for offset in 0..cgroup_count {
        let vertex_index = clear_start
            .checked_add(offset)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        let group_index = first_cgroup_usize
            .checked_add(offset)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        let previous = vertices
            .get(previous_vertex)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let next_iedge = previous.iedge.offset(i64::from(previous.max_adj_edges))?;
        let group = group_values
            .get_mut(group_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let vertex = vertices
            .get_mut(vertex_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        vertex.iedge = next_iedge;
        vertex.max_adj_edges = group.num_edges.wrapping_add(max_add_edges) as AT_NUMB;
        vertex.num_adj_edges = 0;
        vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(group.st_flow);
        total_stationary_flow = total_stationary_flow.wrapping_add(group.st_flow);
        vertex.st_edge.cap = vertex.st_edge.cap.wrapping_add(group.st_cap);
        total_stationary_capacity = total_stationary_capacity.wrapping_add(group.st_cap);
        vertex.st_edge.flow0 = vertex.st_edge.flow;
        vertex.st_edge.cap0 = vertex.st_edge.cap;
        vertex.type_ = group.type_ as AT_NUMB;
        group.nVertexNumber = vertex_index as i32;
        previous_vertex = vertex_index;
    }
    heap.slice_mut(network.vert)?.clone_from_slice(&vertices);
    heap.slice_mut(groups.pTCG)?
        .get_mut(..group_values.len())
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone_from_slice(&group_values);
    current_vertices = i32::try_from(previous_vertex + 1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;

    let mut c_point = 0_i32;
    while c_point < number_of_atoms {
        let cpoint_index = usize::try_from(c_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let list_number = i32::from(valence_atoms[cpoint_index].cnListIndex);
        if list_number == 0 {
            c_point = c_point.wrapping_add(1);
            continue;
        }
        let list_index =
            usize::try_from(list_number.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let charge_structure = CN_LIST.get(list_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let node_valences = CN_LIST_VALENCES
            .get(list_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if node_valences.len() != charge_structure.nodes.len() {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        let metal_atoms = charge_structure.bits == -1;
        let base_vertex = previous_vertex;

        for (node_index, node) in charge_structure.nodes.iter().enumerate() {
            let node_type = node.0;
            let is_auxiliary = node_type & BNS_VERT_TYPE__AUX as i32 != 0 && node_type & BNS_VERT_TYPE_TEMP as i32 != 0;
            if !is_auxiliary {
                continue;
            }
            let vertex_index = base_vertex
                .checked_add(node_index)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            let capacity = if !metal_atoms {
                node.1
            } else if node.1 != 0 {
                restore_mode.nMetalMaxCharge_D
            } else {
                0
            };
            let flow = if !metal_atoms {
                node.2
            } else if node.2 != 0 {
                restore_mode.nMetalMaxCharge_D
            } else {
                0
            };
            {
                let vertices = heap.slice_mut(network.vert)?;
                let previous = vertices
                    .get(previous_vertex)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let next_iedge = previous.iedge.offset(i64::from(previous.max_adj_edges))?;
                let vertex = vertices
                    .get_mut(vertex_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                vertex.iedge = next_iedge;
                vertex.max_adj_edges = node_valences[node_index] as AT_NUMB;
                vertex.num_adj_edges = 0;
                AddStCapFlow(
                    vertex,
                    &mut total_stationary_flow,
                    &mut total_stationary_capacity,
                    capacity,
                    flow,
                );
                vertex.type_ = node_type as AT_NUMB;
            }
            previous_vertex = vertex_index;
            current_vertices = i32::try_from(vertex_index + 1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let vertex = &heap.slice(network.vert.as_const())?[vertex_index];
            let incident_offset = vertex.iedge.difference(network.iedge)?;
            if incident_offset + i64::from(vertex.max_adj_edges) >= i64::from(network.max_iedges)
                || current_vertices >= network.max_vertices
            {
                return Ok(BNS_VERT_EDGE_OVFL);
            }
        }

        for (first_node_index, first_node) in charge_structure.nodes.iter().enumerate() {
            let first_type = first_node.0;
            let mut first_vertex_index = None;
            let mut first_group_index = -1_i32;
            if first_type & BNS_VERT_TYPE_ATOM as i32 != 0 {
                first_vertex_index = Some(c_point);
                let capacity = if !metal_atoms {
                    first_node.1
                } else if first_node.1 != 0 {
                    restore_mode.nMetalMaxCharge_D
                } else {
                    0
                };
                let flow = if !metal_atoms {
                    first_node.2
                } else if first_node.2 != 0 {
                    restore_mode.nMetalMaxCharge_D
                } else {
                    0
                };
                let vertex = heap
                    .slice_mut(network.vert)?
                    .get_mut(cpoint_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                AddStCapFlow(
                    vertex,
                    &mut total_stationary_flow,
                    &mut total_stationary_capacity,
                    capacity,
                    flow,
                );
            } else if first_type & BNS_VT_C_POS_ALL as i32 == BNS_VERT_TYPE_C_GROUP as i32 {
                let vertex_values = heap.slice(network.vert.as_const())?;
                for offset in 0..cgroup_count {
                    let candidate = clear_start + offset;
                    if i32::from(vertex_values[candidate].type_) == first_type {
                        first_vertex_index = Some(candidate as i32);
                        first_group_index = (first_cgroup_usize + offset) as i32;
                        break;
                    }
                }
                if let Some(_) = first_vertex_index {
                    let group = heap
                        .slice(groups.pTCG.as_const())?
                        .get(usize::try_from(first_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if group.type_ != first_type || group.ord_num != 0 {
                        return Ok(RI_ERR_PROGR);
                    }
                }
            } else if first_type == BNS_VERT_TYPE_METAL_GR as i32 {
                first_group_index = groups.nGroup[TCG_MeFlower0 as usize];
                if first_group_index < 0 {
                    return Ok(RI_ERR_PROGR);
                }
                let group = heap
                    .slice(groups.pTCG.as_const())?
                    .get(usize::try_from(first_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if group.type_ != first_type || group.ord_num != 0 || restore_mode.bMetalAddFlower == 0 {
                    return Ok(RI_ERR_PROGR);
                }
                first_vertex_index = Some(group.nVertexNumber);
            } else if first_type & BNS_VERT_TYPE__AUX as i32 != 0 && first_type & BNS_VERT_TYPE_TEMP as i32 != 0 {
                first_vertex_index = Some(
                    i32::try_from(base_vertex + first_node_index)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                );
            }
            let Some(first_vertex) = first_vertex_index else {
                return Ok(BNS_BOND_ERR);
            };

            for &(neighbor_number, source_edge_capacity, source_forbidden, source_edge_flow) in &first_node.3 {
                if neighbor_number == 0 {
                    break;
                }
                let second_node_index = usize::try_from(neighbor_number.wrapping_sub(1))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let second_node = charge_structure
                    .nodes
                    .get(second_node_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let second_type = second_node.0;
                let mut second_vertex_index = None;
                let mut second_group_index = -1_i32;
                if second_type & BNS_VERT_TYPE_ATOM as i32 != 0 {
                    second_vertex_index = Some(c_point);
                    let capacity = if !metal_atoms {
                        second_node.1
                    } else if second_node.1 != 0 {
                        restore_mode.nMetalMaxCharge_D
                    } else {
                        0
                    };
                    let flow = if !metal_atoms {
                        second_node.2
                    } else if second_node.2 != 0 {
                        restore_mode.nMetalMaxCharge_D
                    } else {
                        0
                    };
                    let vertex = heap
                        .slice_mut(network.vert)?
                        .get_mut(cpoint_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    AddStCapFlow(
                        vertex,
                        &mut total_stationary_flow,
                        &mut total_stationary_capacity,
                        capacity,
                        flow,
                    );
                } else if second_type & BNS_VT_C_POS_ALL as i32 == BNS_VERT_TYPE_C_GROUP as i32 {
                    let vertex_values = heap.slice(network.vert.as_const())?;
                    for offset in 0..cgroup_count {
                        let candidate = clear_start + offset;
                        if i32::from(vertex_values[candidate].type_) == second_type {
                            second_vertex_index = Some(candidate as i32);
                            second_group_index = (first_cgroup_usize + offset) as i32;
                            break;
                        }
                    }
                    if second_vertex_index.is_some() {
                        let group = heap
                            .slice(groups.pTCG.as_const())?
                            .get(usize::try_from(second_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if group.type_ != second_type || group.ord_num != 0 {
                            return Ok(RI_ERR_PROGR);
                        }
                    }
                } else if second_type == BNS_VERT_TYPE_METAL_GR as i32 {
                    second_group_index = groups.nGroup[TCG_MeFlower0 as usize];
                    if second_group_index < 0 {
                        return Ok(RI_ERR_PROGR);
                    }
                    let group = heap
                        .slice(groups.pTCG.as_const())?
                        .get(usize::try_from(second_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if group.type_ != second_type || group.ord_num != 0 || restore_mode.bMetalAddFlower == 0 {
                        return Ok(RI_ERR_PROGR);
                    }
                    second_vertex_index = Some(group.nVertexNumber);
                } else if second_type & BNS_VERT_TYPE__AUX as i32 != 0 && second_type & BNS_VERT_TYPE_TEMP as i32 != 0 {
                    second_vertex_index = Some(
                        i32::try_from(base_vertex + second_node_index)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                    );
                }
                let Some(second_vertex) = second_vertex_index else {
                    return Ok(BNS_BOND_ERR);
                };
                if first_vertex == second_vertex {
                    return Ok(BNS_BOND_ERR);
                }

                let edge_index = usize::try_from(current_edges).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let metal_atom_edge = (first_type == BNS_VERT_TYPE_METAL_GR as i32
                    && second_type & BNS_VERT_TYPE_ATOM as i32 != 0)
                    || (second_type == BNS_VERT_TYPE_METAL_GR as i32 && first_type & BNS_VERT_TYPE_ATOM as i32 != 0);
                if metal_atom_edge {
                    let mut stationary_capacity = 0_i32;
                    let mut stationary_flow = 0_i32;
                    let mut edge_capacity = 0_i32;
                    let mut edge_flow = 0_i32;
                    let needs_flower = AtomStcapStflow(
                        &atoms,
                        valence_atoms,
                        &restore_mode,
                        c_point,
                        Some(&mut stationary_capacity),
                        Some(&mut stationary_flow),
                        Some(&mut edge_capacity),
                        Some(&mut edge_flow),
                    )?;
                    {
                        let edge = heap
                            .slice_mut(network.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.cap = edge_capacity;
                        edge.flow = edge_flow;
                    }
                    if needs_flower == 0 {
                        return Ok(RI_ERR_PROGR);
                    }
                    valence_atoms[cpoint_index].nMetalGroupEdge = current_edges.wrapping_add(1);
                    let free_valences = i32::from(valence_atoms[cpoint_index].cInitFreeValences);
                    let atom_vertex = heap
                        .slice_mut(network.vert)?
                        .get_mut(cpoint_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    atom_vertex.st_edge.flow = atom_vertex.st_edge.flow.wrapping_add(edge_flow);
                    atom_vertex.st_edge.cap = atom_vertex
                        .st_edge
                        .cap
                        .wrapping_add(edge_flow.wrapping_add(free_valences));
                    atom_vertex.st_edge.flow0 = atom_vertex.st_edge.flow;
                    atom_vertex.st_edge.cap0 = atom_vertex.st_edge.cap;
                    total_stationary_flow = total_stationary_flow.wrapping_add(edge_flow);
                    total_stationary_capacity =
                        total_stationary_capacity.wrapping_add(edge_flow.wrapping_add(free_valences));
                } else {
                    let edge_capacity = if !metal_atoms {
                        source_edge_capacity
                    } else if source_edge_capacity != 0 {
                        restore_mode.nMetalMaxCharge_D
                    } else {
                        0
                    };
                    let edge_flow = if !metal_atoms {
                        source_edge_flow
                    } else if source_edge_flow != 0 {
                        restore_mode.nMetalMaxCharge_D
                    } else {
                        0
                    };
                    let edge = heap
                        .slice_mut(network.edge)?
                        .get_mut(edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.cap = edge_capacity;
                    edge.flow = edge_flow;
                }
                {
                    let edge = heap
                        .slice_mut(network.edge)?
                        .get_mut(edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.forbidden = if source_forbidden != 0 {
                        BNS_EDGE_FORBIDDEN_MASK as i8
                    } else {
                        0
                    };
                    edge.pass = 0;
                }

                let vertex_values = heap.slice(network.vert.as_const())?;
                let first_usize = usize::try_from(first_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let second_usize = usize::try_from(second_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let first_vertex_value = vertex_values
                    .get(first_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let second_vertex_value = vertex_values
                    .get(second_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if first_vertex_value.num_adj_edges >= first_vertex_value.max_adj_edges
                    || second_vertex_value.num_adj_edges >= second_vertex_value.max_adj_edges
                    || current_edges >= network.max_edges
                {
                    return Ok(BNS_VERT_EDGE_OVFL);
                }
                let connect_result = ConnectTwoVertices(heap, network, first_vertex, second_vertex, current_edges, 0)?;
                if connect_result >= BNS_ERR && connect_result <= BNS_MAX_ERR_VALUE {
                    return Ok(connect_result);
                }
                let edge = heap
                    .slice_mut(network.edge)?
                    .get_mut(edge_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                edge.cap0 = edge.cap;
                edge.flow0 = edge.flow;

                let first_vertex_type = i32::from(heap.slice(network.vert.as_const())?[first_usize].type_);
                let second_vertex_type = i32::from(heap.slice(network.vert.as_const())?[second_usize].type_);
                let cgroup_type = if first_vertex_type & BNS_VT_C_POS_ALL as i32 == BNS_VERT_TYPE_C_GROUP as i32 {
                    first_vertex_type
                } else if second_vertex_type & BNS_VT_C_POS_ALL as i32 == BNS_VERT_TYPE_C_GROUP as i32 {
                    second_vertex_type
                } else {
                    0
                };
                if cgroup_type != 0 {
                    if cgroup_type & BNS_VERT_TYPE_C_NEGATIVE as i32 != 0 {
                        valence_atoms[cpoint_index].nCMinusGroupEdge = current_edges.wrapping_add(1);
                    } else {
                        valence_atoms[cpoint_index].nCPlusGroupEdge = current_edges.wrapping_add(1);
                    }
                }
                current_edges = current_edges.wrapping_add(1);
                let _ = (first_group_index, second_group_index);
            }
        }
        c_point = c_point.wrapping_add(1);
    }

    c_point = 0;
    while c_point < number_of_atoms {
        let cpoint_index = usize::try_from(c_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let vertex = heap
            .slice(network.vert.as_const())?
            .get(cpoint_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let stationary_capacity = vertex.st_edge.cap;
        let mut adjacency = 0_i32;
        while adjacency < i32::from(vertex.num_adj_edges) {
            let adjacency_usize = usize::try_from(adjacency).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let edge_number = *heap
                .slice(vertex.iedge.as_const())?
                .get(adjacency_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let edge_index = usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let edge_snapshot = heap
                .slice(network.edge.as_const())?
                .get(edge_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let neighbor = i32::from(edge_snapshot.neighbor12) ^ c_point;
            let neighbor_usize = usize::try_from(neighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let neighbor_vertex = heap
                .slice(network.vert.as_const())?
                .get(neighbor_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor > c_point || neighbor_vertex.type_ & BNS_VERT_TYPE_ATOM as AT_NUMB == 0 {
                adjacency = adjacency.wrapping_add(1);
                continue;
            }
            let mut maximum_edge_flow = stationary_capacity.min(neighbor_vertex.st_edge.cap);
            let metal_capacity = (restore_mode.bMetalAddFlower != 0
                && restore_mode.nMetalMinBondOrder == 0
                && valence_atoms[cpoint_index].cMetal != 0
                && valence_atoms[cpoint_index].cNumBondsToMetal != 0)
                || (valence_atoms
                    .get(neighbor_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .cMetal
                    != 0
                    && valence_atoms[neighbor_usize].cNumBondsToMetal != 0);
            maximum_edge_flow = maximum_edge_flow.min(if metal_capacity {
                MAX_BOND_EDGE_CAP as i32 + 1
            } else {
                MAX_BOND_EDGE_CAP as i32
            });
            if u32::from(
                *atoms[cpoint_index]
                    .bond_type
                    .get(adjacency_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            ) == BOND_TYPE_SINGLE
            {
                let edge = heap
                    .slice_mut(network.edge)?
                    .get_mut(edge_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                edge.cap = maximum_edge_flow;
                edge.cap0 = maximum_edge_flow;
            }
            adjacency = adjacency.wrapping_add(1);
        }
        c_point = c_point.wrapping_add(1);
    }

    let mut result = ConnectMetalFlower(
        heap,
        &mut current_vertices,
        &mut current_edges,
        &mut total_stationary_capacity,
        &mut total_stationary_flow,
        &restore_mode,
        network,
        groups,
    )?;
    if result < 0 {
        return Ok(result);
    }

    let _plus_result = ConnectSuperCGroup(
        heap,
        TCG_Plus,
        &[TCG_Plus0, TCG_Plus_C0, TCG_Plus_M0],
        3,
        &mut current_vertices,
        &mut current_edges,
        &mut total_stationary_capacity,
        &mut total_stationary_flow,
        network,
        groups,
    )?;
    let _minus_result = ConnectSuperCGroup(
        heap,
        TCG_Minus,
        &[TCG_Minus0, TCG_Minus_C0, TCG_Minus_M0],
        3,
        &mut current_vertices,
        &mut current_edges,
        &mut total_stationary_capacity,
        &mut total_stationary_flow,
        network,
        groups,
    )?;
    let plus_minus_result = ConnectSuperCGroup(
        heap,
        -1,
        &[TCG_Plus, TCG_Minus],
        2,
        &mut current_vertices,
        &mut current_edges,
        &mut total_stationary_capacity,
        &mut total_stationary_flow,
        network,
        groups,
    )?;

    let cgroup_charge = groups
        .total_charge
        .wrapping_sub(groups.tgroup_charge)
        .wrapping_sub(groups.charge_on_atoms);
    result = 1;
    if plus_minus_result > 0 {
        let plus_minus_vertex = current_vertices.wrapping_sub(1);
        let plus_minus_usize = usize::try_from(plus_minus_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let vertex_snapshot = heap
            .slice(network.vert.as_const())?
            .get(plus_minus_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let plus_exists = groups.nGroup[TCG_Plus as usize] >= 0;
        let minus_exists = groups.nGroup[TCG_Minus as usize] >= 0;
        let existing_groups = i32::from(plus_exists).wrapping_add(i32::from(minus_exists));
        let mut plus_edge = None;
        let mut minus_edge = None;
        if vertex_snapshot.num_adj_edges == 2 && existing_groups == 2 {
            plus_edge = Some(
                *heap
                    .slice(vertex_snapshot.iedge.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            minus_edge = Some(
                *heap
                    .slice(vertex_snapshot.iedge.as_const())?
                    .get(1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
        } else if vertex_snapshot.num_adj_edges == 1 && existing_groups == 1 {
            let edge = *heap
                .slice(vertex_snapshot.iedge.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if plus_exists {
                plus_edge = Some(edge);
            } else if minus_exists {
                minus_edge = Some(edge);
            }
        } else if existing_groups != 0 {
            return Ok(BNS_BOND_ERR);
        }

        let plus_vertex = if let Some(edge_index) = plus_edge {
            let edge = heap
                .slice(network.edge.as_const())?
                .get(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            Some(i32::from(edge.neighbor12) ^ plus_minus_vertex)
        } else {
            None
        };
        let minus_vertex = if let Some(edge_index) = minus_edge {
            let edge = heap
                .slice(network.edge.as_const())?
                .get(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            Some(i32::from(edge.neighbor12) ^ plus_minus_vertex)
        } else {
            None
        };
        let selected_vertex = plus_vertex.or(minus_vertex);
        let selected_edge = plus_edge.or(minus_edge);
        if let Some(edge) = minus_edge {
            groups.nEdgeMinus = edge;
        }
        if let Some(edge) = plus_edge {
            groups.nEdgePlus = edge;
        }
        if let Some(edge) = selected_edge {
            groups.nEdge4charge = edge;
        }

        if let (Some(vertex_number), Some(edge_number)) = (selected_vertex, selected_edge) {
            if cgroup_charge > 0 {
                let vertex = heap
                    .slice_mut(network.vert)?
                    .get_mut(plus_minus_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                vertex.st_edge.cap = vertex.st_edge.cap.wrapping_add(cgroup_charge);
                total_stationary_capacity = total_stationary_capacity.wrapping_add(cgroup_charge);
                vertex.st_edge.cap0 = vertex.st_edge.cap;
            }
            if cgroup_charge < 0 {
                let vertex = heap
                    .slice_mut(network.vert)?
                    .get_mut(usize::try_from(vertex_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                vertex.st_edge.cap = vertex.st_edge.cap.wrapping_sub(cgroup_charge);
                total_stationary_capacity = total_stationary_capacity.wrapping_sub(cgroup_charge);
                vertex.st_edge.cap0 = vertex.st_edge.cap;
                let edge = heap
                    .slice_mut(network.edge)?
                    .get_mut(usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if edge.cap.wrapping_sub(edge.flow).wrapping_add(cgroup_charge) < 0 {
                    edge.cap = edge.flow.wrapping_sub(cgroup_charge);
                }
            }
            groups.added_charge = cgroup_charge;
        }
        if cgroup_charge == 0 || (selected_vertex.is_some() && selected_edge.is_some()) {
            result = 2;
        }
    }

    let _ = AddRadicalToMetal(
        heap,
        &mut total_stationary_capacity,
        &mut total_stationary_flow,
        &restore_mode,
        network,
        groups,
    )?;
    network.num_edges = current_edges;
    network.num_vertices = current_vertices;
    network.num_c_groups = number_of_cgroups;
    network.tot_st_cap = total_stationary_capacity;
    network.tot_st_flow = total_stationary_flow;
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn nNumEdgesToCnVertex(
    charge_nodes: &[CnNodeData],
    length: i32,
    vertex: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3824 nNumEdgesToCnVertex
    // INCHI✔️✔️: int nNumEdgesToCnVertex( MY_CONST C_NODE *pCN, int len, int v )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j, n, num_edges, v1 = v + 1;
    // INCHI✔️✔️:     for (i = 0, num_edges = 0; i < len; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (j = 0; j < MAX_CN_VAL && ( n = pCN[i].e[j].neigh ); j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             num_edges += ( i == v || n == v1 );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return num_edges;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nNumEdgesToCnVertex

    if length <= 0 {
        return Ok(0);
    }
    let count = usize::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let nodes = charge_nodes.get(..count).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let vertex_number = vertex.wrapping_add(1);
    let mut number_of_edges = 0_i32;
    for (node_index, node) in nodes.iter().enumerate() {
        let node_index = i32::try_from(node_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        for &(neighbor, _, _, _) in &node.3 {
            if neighbor == 0 {
                break;
            }
            number_of_edges =
                number_of_edges.wrapping_add(i32::from(node_index == vertex || neighbor == vertex_number));
        }
    }
    Ok(number_of_edges)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AllocateAndInitTCGBnStruct(
    heap: &mut SourceHeap,
    structure: &StrFromINChI,
    valence_atoms: &[crate::source_types::VAL_AT],
    groups: &ALL_TC_GROUPS,
    nMaxAddAtoms: i32,
    nMaxAddEdges: i32,
    max_altp: i32,
    pNum_changed_bonds: &mut i32,
) -> Result<SourceMutPointer<BN_STRUCT>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3879 AllocateAndInitTCGBnStruct
    // INCHI✔️❌: BN_STRUCT* AllocateAndInitTCGBnStruct( StrFromINChI *pStruct, VAL_AT *pVA,
    // INCHI✔️❌:                                        ALL_TC_GROUPS *pTCGroups,
    // INCHI✔️❌:                                        int nMaxAddAtoms, int nMaxAddEdges,
    // INCHI✔️❌:                                        int max_altp, int *pNum_changed_bonds )
    // INCHI✔️❌: {
    // INCHI✔️❌:     inp_ATOM *at = pStruct->at;
    // INCHI✔️❌:     int          num_atoms = pStruct->num_atoms;
    // INCHI✔️❌:     ICHICONST SRM *pSrm = pStruct->pSrm;
    // INCHI✔️❌:
    // INCHI✔️❌:     BN_STRUCT   *pBNS = NULL;
    // INCHI✔️❌:     BNS_VERTEX  *vert;
    // INCHI✔️❌:     BNS_IEDGE   *iedge;
    // INCHI✔️❌:
    // INCHI✔️❌:     int    neigh, num_changed_bonds = 0;
    // INCHI✔️❌:     U_CHAR bond_type, bond_mark;
    // INCHI✔️❌:     int bNeedsFlower1, bNeedsFlower2, min_order;
    // INCHI✔️❌:
    // INCHI✔️❌:     int i, j, k, m, n_edges, num_bonds, num_edges;
    // INCHI✔️❌:     int f1, f2, c1, c2, edge_cap, edge_flow, st_cap, st_flow, flag_alt_bond;
    // INCHI✔️❌:     int tot_st_cap, tot_st_flow;
    // INCHI✔️❌:     int max_tg, max_edges, max_vertices, len_alt_path, max_iedges, num_iedges, num_altp;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count vertices */
    // INCHI✔️❌:     max_tg = pTCGroups->num_tgroups;
    // INCHI✔️❌:     /* +1 for a super-tautomeric group */
    // INCHI✔️❌:     /* max_vertices = num_atoms + nMaxAddAtoms + max_tg + 1; */
    // INCHI✔️❌:     max_vertices = pTCGroups->nVertices + nMaxAddAtoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count edges */
    // INCHI✔️❌:     num_changed_bonds = 0;
    // INCHI✔️❌:     num_bonds = pTCGroups->num_bonds;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* each atom has enough edges to belong to a tautomeric group + nMaxAddEdges */
    // INCHI✔️❌:     /* number of atoms is large enough to accommodate max. possible number of t-groups + nMaxAddAtoms */
    // INCHI✔️❌:     /* max_altp cannot be larger than BN_MAX_ALTP = 16 */
    // INCHI✔️❌:     num_edges = pTCGroups->nEdges;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* +max_tg for edges between t-groups and super-tautomeric group */
    // INCHI✔️❌:     max_edges = num_edges + ( nMaxAddEdges + NUM_KINDS_OF_GROUPS )*max_vertices;
    // INCHI✔️❌:     max_iedges = 2 * max_edges + pTCGroups->nAddIedges;
    // INCHI✔️❌:     len_alt_path = max_vertices + iALTP_HDR_LEN + 1; /* may overflow if an edge is traversed in 2 directions */
    // INCHI✔️❌:     len_alt_path += inchi_max( max_vertices / 2, 16 ); /* to avoid the overflow */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!(pBNS = (BN_STRUCT*)inchi_calloc(1, sizeof(BN_STRUCT))) ||
    // INCHI✔️❌:         !(pBNS->edge = (BNS_EDGE*)inchi_calloc(max_edges, sizeof(BNS_EDGE))) ||
    // INCHI✔️❌:         !(pBNS->vert = (BNS_VERTEX*)inchi_calloc(max_vertices, sizeof(BNS_VERTEX))) ||
    // INCHI✔️❌:         !(pBNS->iedge = (BNS_IEDGE*)inchi_calloc(max_iedges, sizeof(BNS_IEDGE))))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return DeAllocateBnStruct( pBNS );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* alt path init (standard spell) */
    // INCHI✔️❌:     for (num_altp = 0; num_altp < max_altp && num_altp < BN_MAX_ALTP; num_altp++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!( pBNS->altp[num_altp] = (BNS_ALT_PATH*) inchi_calloc( len_alt_path, sizeof( BNS_ALT_PATH ) ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return DeAllocateBnStruct( pBNS );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ALTP_ALLOCATED_LEN( pBNS->altp[num_altp] ) = len_alt_path;
    // INCHI✔️❌:         pBNS->len_alt_path = len_alt_path;  /* ??? duplication ??? */
    // INCHI✔️❌:         /* re-init */
    // INCHI✔️❌:         ALTP_DELTA( pBNS->altp[num_altp] ) = 0;
    // INCHI✔️❌:         ALTP_START_ATOM( pBNS->altp[num_altp] ) = NO_VERTEX;
    // INCHI✔️❌:         ALTP_END_ATOM( pBNS->altp[num_altp] ) = NO_VERTEX;
    // INCHI✔️❌:         ALTP_PATH_LEN( pBNS->altp[num_altp] ) = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     pBNS->alt_path = NULL;
    // INCHI✔️❌:     pBNS->num_altp = 0;
    // INCHI✔️❌:     pBNS->max_altp = num_altp;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* fill vertices (no connectivity) */
    // INCHI✔️❌:     iedge = pBNS->iedge;
    // INCHI✔️❌:     num_iedges = 0;
    // INCHI✔️❌:     tot_st_cap = tot_st_flow = 0;
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* count edges incident to pBNS->vert[i] */
    // INCHI✔️❌:         k = at[i].valence + ( at[i].endpoint != 0 ) + ( nMaxAddEdges /*+ NUM_KINDS_OF_GROUPS*/ );
    // INCHI✔️❌:         if (( j = pVA[i].cnListIndex - 1 ) >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* add number of neighbors in the ChargeStruct */
    // INCHI✔️❌:             k += nNumEdgesToCnVertex( cnList[j].pCN, cnList[j].len, 0 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* set max number of edges for the vertex */
    // INCHI✔️❌:         pBNS->vert[i].max_adj_edges = k;
    // INCHI✔️❌:         pBNS->vert[i].iedge = iedge;
    // INCHI✔️❌:         iedge += k;
    // INCHI✔️❌:         /* add atom vertex cap */
    // INCHI✔️❌:         st_cap = 0;
    // INCHI✔️❌:         st_flow = 0;
    // INCHI✔️❌:         bNeedsFlower1 = AtomStcapStflow( at, pVA, pSrm, i, &c1, &f1, NULL, NULL );
    // INCHI✔️❌:         /* pVA[i].cNumBondsToMetal = bNeedsFlower1; */
    // INCHI✔️❌:         /* GetAtomStCapFlow( at, pVA, pSrm, i, &c1, &f1 ); */
    // INCHI✔️❌:         st_cap += c1;
    // INCHI✔️❌:         st_cap += bNeedsFlower1 ? 0 : pVA[i].cInitFreeValences;
    // INCHI✔️❌:         pBNS->vert[i].st_edge.cap = st_cap; /* the 1st time st_cap is set */
    // INCHI✔️❌:         pBNS->vert[i].st_edge.cap0 = pBNS->vert[i].st_edge.cap;
    // INCHI✔️❌:         tot_st_cap += st_cap;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     num_iedges = (int) ( iedge - pBNS->iedge );
    // INCHI✔️❌:     if (max_iedges - num_iedges < ( nMaxAddEdges + NUM_KINDS_OF_GROUPS )*max_vertices)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return DeAllocateBnStruct( pBNS );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     pBNS->num_atoms = num_atoms;      /* number of real atoms */
    // INCHI✔️❌:     pBNS->num_added_atoms = 0;
    // INCHI✔️❌:     pBNS->num_t_groups = 0;              /* number of added t-groups */
    // INCHI✔️❌:     pBNS->num_c_groups = 0;
    // INCHI✔️❌:     pBNS->nMaxAddAtoms = nMaxAddAtoms;
    // INCHI✔️❌:     pBNS->nMaxAddEdges = nMaxAddEdges;
    // INCHI✔️❌:
    // INCHI✔️❌:     pBNS->num_vertices = num_atoms;      /* current number of vertices, in general a sum of
    // INCHI✔️❌:                                                pBNS->num_atoms
    // INCHI✔️❌:                                                pBNS->num_t_groups
    // INCHI✔️❌:                                                number of c-groups
    // INCHI✔️❌:                                                number of auxiliary vertices
    // INCHI✔️❌:                                                pBNS->num_added_atoms
    // INCHI✔️❌:                                             */
    // INCHI✔️❌:     pBNS->max_vertices = max_vertices;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     pBNS->num_bonds = num_bonds;      /* number of real edges (bonds) */
    // INCHI✔️❌:     pBNS->max_edges = max_edges;
    // INCHI✔️❌:     pBNS->max_iedges = max_iedges;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:        To remove t-groups and added atoms:
    // INCHI✔️❌:
    // INCHI✔️❌:         for ( i = 0; i < pBNS->num_atoms; i ++ )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for ( j = pBNS->vert[i].num_adj_edges-1; 0 <= j; j -- )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k = pBNS->edge[pBNS->vert[i].iedge[j]].neighbor12 ^ i;
    // INCHI✔️❌:                 if ( pBNS->vert[k].type & BNS_VERT_TYPE_ATOM )
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pBNS->vert[i].num_adj_edges = j+1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:        pBNS->num_vertices    = pBNS->num_atoms;
    // INCHI✔️❌:        pBNS->num_edges       = pBNS->num_bonds;
    // INCHI✔️❌:        pBNS->num_added_atoms = 0;
    // INCHI✔️❌:        pBNS->num_t_groups    = 0;
    // INCHI✔️❌:        pBNS->num_added_edges = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:        ALTP_DELTA(pBNS->alt_path)      = 0;
    // INCHI✔️❌:        ALTP_START_ATOM(pBNS->alt_path) = NO_VERTEX;
    // INCHI✔️❌:        ALTP_END_ATOM(pBNS->alt_path)   = NO_VERTEX;
    // INCHI✔️❌:        ALTP_PATH_LEN(pBNS->alt_path)   = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* add and fill edges and connectivity */
    // INCHI✔️❌:     for (i = 0, n_edges = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         vert = pBNS->vert + i; /* pointer to the ith vertex */
    // INCHI✔️❌:         st_cap = 0;
    // INCHI✔️❌:         st_flow = 0;
    // INCHI✔️❌:         flag_alt_bond = 0;
    // INCHI✔️❌:         for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[i].neighbor[j];
    // INCHI✔️❌:             /* find this bond at the neighbor */
    // INCHI✔️❌:             for (k = 0; k < at[neigh].valence; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].neighbor[k] == i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             bond_type = ( at[i].bond_type[j] & BOND_TYPE_MASK );
    // INCHI✔️❌:             bond_mark = ( at[i].bond_type[j] & ~BOND_TYPE_MASK );
    // INCHI✔️❌:             if (bond_type != BOND_SINGLE && bond_type != BOND_DOUBLE && bond_type != BOND_TRIPLE)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* make unknown bonds single */
    // INCHI✔️❌:                 bond_type = BOND_SINGLE;
    // INCHI✔️❌:                 at[i].bond_type[j] = bond_mark | bond_type;
    // INCHI✔️❌:                 num_changed_bonds++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (neigh > i)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* this is the first time we encounter this bond */
    // INCHI✔️❌:                 bNeedsFlower1 = AtomStcapStflow( at, pVA, pSrm, i, &c1, &f1, NULL, NULL );
    // INCHI✔️❌:                 /* GetAtomStCapFlow( at, pVA, pSrm, i, &c1, &f1 ); */
    // INCHI✔️❌:                 c1 += bNeedsFlower1 ? 0 : pVA[i].cInitFreeValences;  /* elevate cap to the lowest valence in ChargeStruct */
    // INCHI✔️❌:                 bNeedsFlower2 = AtomStcapStflow( at, pVA, pSrm, neigh, &c2, &f2, NULL, NULL );
    // INCHI✔️❌:                 /* GetAtomStCapFlow( at, pVA, pSrm, neigh, &c2, &f2 ); */
    // INCHI✔️❌:                 c2 += bNeedsFlower2 ? 0 : pVA[neigh].cInitFreeValences; /* elevate cap to the lowest valence in ChargeStruct */
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* at this point -O would have st_cap=st_flow=0 because the lowest valence=1 for charge=-1 */
    // INCHI✔️❌:                 /* however, if -O belongs to a t-group its cap would be 1, flow = 0 */
    // INCHI✔️❌:                 /*f1 = MAX_AT_FLOW(at[i]);*/
    // INCHI✔️❌:                 /*f2 = MAX_AT_FLOW(at[neigh]);*/
    // INCHI✔️❌:
    // INCHI✔️❌:                 edge_flow = BondFlowMaxcapMinorder( at, pVA, pSrm, i, j, &edge_cap, &min_order, NULL );
    // INCHI✔️❌:
    // INCHI✔️❌:                 pBNS->edge[n_edges].neighbor1 = (AT_NUMB) i;
    // INCHI✔️❌:                 pBNS->edge[n_edges].neighbor12 = (AT_NUMB) ( i ^ neigh );
    // INCHI✔️❌:                 pBNS->edge[n_edges].flow =
    // INCHI✔️❌:                     pBNS->edge[n_edges].flow0 = edge_flow;
    // INCHI✔️❌:                 pBNS->edge[n_edges].cap =
    // INCHI✔️❌:                     pBNS->edge[n_edges].cap0 = edge_cap;
    // INCHI✔️❌:                 pBNS->edge[n_edges].neigh_ord[0] = j;  /* iedge to neigh index at vertex[i],     i < neigh */
    // INCHI✔️❌:                 pBNS->edge[n_edges].neigh_ord[1] = k;  /* iedge to i     index at vertex[neigh], i < neigh */
    // INCHI✔️❌:                 pBNS->edge[n_edges].pass = 0;
    // INCHI✔️❌:                 pBNS->edge[n_edges].forbidden = 0; /* may be forbidden if edge_flow = 1: stereogenic fixed double bond */
    // INCHI✔️❌:                 if (bond_type == BOND_TYPE_DOUBLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* forbid changing stereogenic double bonds */
    // INCHI✔️❌:                     for (m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (at[i].sb_ord[m] == j)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             pBNS->edge[n_edges].forbidden |= BNS_EDGE_FORBIDDEN_MASK;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (pBNS->vert[neigh].iedge) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     vert->iedge[j] = pBNS->vert[neigh].iedge[k] = n_edges++; /* same iedge index as neighbor index in at[] */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* this is the second time we encounter this bond. It was stored at */
    // INCHI✔️❌:                 int  iedge2 = pBNS->vert[neigh].iedge[k];
    // INCHI✔️❌:                 edge_cap = pBNS->edge[iedge2].cap;
    // INCHI✔️❌:                 edge_flow = pBNS->edge[iedge2].flow;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             st_flow += edge_flow;
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             st_cap  += edge_cap;
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         vert->num_adj_edges = j;
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         vert->st_edge.cap   =
    // INCHI✔️❌:         vert->st_edge.cap0  = st_cap;
    // INCHI✔️❌:         */
    // INCHI✔️❌:         vert->st_edge.flow =
    // INCHI✔️❌:             vert->st_edge.flow0 = st_flow;
    // INCHI✔️❌:         vert->type = BNS_VERT_TYPE_ATOM;
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         tot_st_cap  += vert->st_edge.cap;
    // INCHI✔️❌:         */
    // INCHI✔️❌:         tot_st_flow += vert->st_edge.flow;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *pNum_changed_bonds = num_changed_bonds / 2;
    // INCHI✔️❌:
    // INCHI✔️❌:     pBNS->num_edges = n_edges;   /* number of edges */
    // INCHI✔️❌:     pBNS->num_iedges = num_iedges;
    // INCHI✔️❌:     pBNS->num_added_edges = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     pBNS->tot_st_cap = tot_st_cap;
    // INCHI✔️❌:     pBNS->tot_st_flow = tot_st_flow;
    // INCHI✔️❌:
    // INCHI✔️❌: /* exit_function: */
    // INCHI✔️❌:
    // INCHI✔️❌:     return pBNS;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AllocateAndInitTCGBnStruct
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AllocateAndInitTCGBnStruct
    // INCHI✔️❌: #define inchi_calloc calloc
    // INCHI✔️❌: #define inchi_free(X) do{ if(X) free(X); }while(0)
    // INCHI✔️❌: #define BN_MAX_ALTP 16
    // INCHI✔️❌: #define NUM_KINDS_OF_GROUPS 2
    // INCHI✔️❌: #define BOND_TYPE_MASK 0x0f
    // INCHI✔️❌: #define MAX_NUM_STEREO_BONDS 3
    // INCHI✔️❌: #define BNS_EDGE_FORBIDDEN_MASK 0x01
    // END INCHI ACTIVE MACRO CONFIGURATION: AllocateAndInitTCGBnStruct

    fn source_count(value: i32) -> Result<u64, SourceHeapError> {
        u64::try_from(value).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)
    }

    fn cleanup(
        heap: &mut SourceHeap,
        network: SourceMutPointer<BN_STRUCT>,
        allocated_alt_paths: i32,
    ) -> Result<SourceMutPointer<BN_STRUCT>, SourceHeapError> {
        if !network.is_null() {
            let (vertices, incident_edges) = {
                let stored = heap
                    .slice_mut(network)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                stored.max_altp = allocated_alt_paths;
                (stored.vert, stored.iedge)
            };
            if !vertices.is_null() && !incident_edges.is_null() {
                if let Some(first_vertex) = heap.slice_mut(vertices)?.first_mut() {
                    first_vertex.iedge = incident_edges;
                }
            }
        }
        DeAllocateBnStruct(heap, network)
    }

    let num_atoms = structure.num_atoms;
    let atom_count = if num_atoms > 0 {
        usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    let mut atoms = if atom_count == 0 {
        Vec::new()
    } else {
        heap.slice(structure.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    if valence_atoms.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let restore_mode = if atom_count == 0 {
        SRM::default()
    } else {
        heap.slice(structure.pSrm)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone()
    };

    let _max_tg = groups.num_tgroups;
    let max_vertices = groups.nVertices.wrapping_add(nMaxAddAtoms);
    let num_bonds = groups.num_bonds;
    let num_edges = groups.nEdges;
    let reserve_per_vertex = nMaxAddEdges.wrapping_add(NUM_KINDS_OF_GROUPS as i32);
    let max_edges = num_edges.wrapping_add(reserve_per_vertex.wrapping_mul(max_vertices));
    let max_iedges = 2_i32.wrapping_mul(max_edges).wrapping_add(groups.nAddIedges);
    let mut len_alt_path = max_vertices
        .wrapping_add(tagAltPathConst_iALTP_HDR_LEN as i32)
        .wrapping_add(1);
    len_alt_path = len_alt_path.wrapping_add((max_vertices / 2).max(16));

    let network = match inchi_calloc::<BN_STRUCT>(heap, 1, std::mem::size_of::<BN_STRUCT>() as u64) {
        Ok(pointer) => pointer,
        Err(_) => return Ok(SourceMutPointer::null()),
    };

    let edges = match source_count(max_edges)
        .and_then(|count| inchi_calloc::<BNS_EDGE>(heap, count, std::mem::size_of::<BNS_EDGE>() as u64))
    {
        Ok(pointer) => pointer,
        Err(_) => return cleanup(heap, network, 0),
    };
    heap.slice_mut(network)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .edge = edges;

    let vertices = match source_count(max_vertices)
        .and_then(|count| inchi_calloc::<BNS_VERTEX>(heap, count, std::mem::size_of::<BNS_VERTEX>() as u64))
    {
        Ok(pointer) => pointer,
        Err(_) => return cleanup(heap, network, 0),
    };
    heap.slice_mut(network)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .vert = vertices;

    let incident_edges = match source_count(max_iedges)
        .and_then(|count| inchi_calloc::<BNS_IEDGE>(heap, count, std::mem::size_of::<BNS_IEDGE>() as u64))
    {
        Ok(pointer) => pointer,
        Err(_) => return cleanup(heap, network, 0),
    };
    heap.slice_mut(network)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .iedge = incident_edges;

    let mut num_altp = 0_i32;
    while num_altp < max_altp && num_altp < BN_MAX_ALTP as i32 {
        let alt_path = match source_count(len_alt_path)
            .and_then(|count| inchi_calloc::<BNS_ALT_PATH>(heap, count, std::mem::size_of::<BNS_ALT_PATH>() as u64))
        {
            Ok(pointer) => pointer,
            Err(_) => return cleanup(heap, network, num_altp),
        };
        {
            let path = heap.slice_mut(alt_path)?;
            path[tagAltPathConst_iALTP_MAX_LEN as usize].set_number(len_alt_path);
            path[tagAltPathConst_iALTP_FLOW as usize].set_flow(0, 0);
            path[tagAltPathConst_iALTP_START_ATOM as usize].set_number(NO_VERTEX);
            path[tagAltPathConst_iALTP_END_ATOM as usize].set_number(NO_VERTEX);
            path[tagAltPathConst_iALTP_PATH_LEN as usize].set_number(0);
        }
        let network_value = heap
            .slice_mut(network)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        network_value.altp[usize::try_from(num_altp).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = alt_path;
        network_value.len_alt_path = len_alt_path;
        num_altp = num_altp.wrapping_add(1);
    }

    let mut network_value = heap
        .slice(network.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    network_value.alt_path = SourceMutPointer::null();
    network_value.num_altp = 0;
    network_value.max_altp = num_altp;

    let mut current_incident = network_value.iedge;
    let mut total_stationary_capacity = 0_i32;
    for i in 0..atom_count {
        let mut capacity = i32::from(atoms[i].valence)
            .wrapping_add(i32::from(atoms[i].endpoint != 0))
            .wrapping_add(nMaxAddEdges);
        let charge_list_index = i32::from(valence_atoms[i].cnListIndex).wrapping_sub(1);
        if charge_list_index >= 0 {
            let list = CN_LIST
                .get(usize::try_from(charge_list_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            capacity = capacity.wrapping_add(nNumEdgesToCnVertex(
                list.nodes,
                i32::try_from(list.nodes.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                0,
            )?);
        }

        {
            let vertex = heap
                .slice_mut(network_value.vert)?
                .get_mut(i)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            vertex.max_adj_edges = capacity as AT_NUMB;
            vertex.iedge = current_incident;
        }
        current_incident = current_incident.offset(i64::from(capacity))?;

        let mut atom_capacity = 0_i32;
        let mut atom_flow = 0_i32;
        let needs_flower = AtomStcapStflow(
            &atoms,
            valence_atoms,
            &restore_mode,
            i32::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
            Some(&mut atom_capacity),
            Some(&mut atom_flow),
            None,
            None,
        )?;
        atom_capacity = atom_capacity.wrapping_add(if needs_flower != 0 {
            0
        } else {
            i32::from(valence_atoms[i].cInitFreeValences)
        });
        let vertex = heap
            .slice_mut(network_value.vert)?
            .get_mut(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        vertex.st_edge.cap = atom_capacity;
        vertex.st_edge.cap0 = atom_capacity;
        total_stationary_capacity = total_stationary_capacity.wrapping_add(atom_capacity);
    }

    let num_iedges = current_incident.difference(network_value.iedge)? as i32;
    if max_iedges.wrapping_sub(num_iedges) < reserve_per_vertex.wrapping_mul(max_vertices) {
        return cleanup(heap, network, num_altp);
    }

    network_value.num_atoms = num_atoms;
    network_value.num_added_atoms = 0;
    network_value.num_t_groups = 0;
    network_value.num_c_groups = 0;
    network_value.nMaxAddAtoms = nMaxAddAtoms;
    network_value.nMaxAddEdges = nMaxAddEdges;
    network_value.num_vertices = num_atoms;
    network_value.max_vertices = max_vertices;
    network_value.num_bonds = num_bonds;
    network_value.max_edges = max_edges;
    network_value.max_iedges = max_iedges;

    let mut total_stationary_flow = 0_i32;
    let mut n_edges = 0_i32;
    let mut num_changed_bonds = 0_i32;
    for i in 0..atom_count {
        let mut stationary_flow = 0_i32;
        let mut j = 0_i32;
        let atom_valence = i32::from(atoms[i].valence);
        while j < atom_valence {
            let j_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let neighbor = i32::from(atoms[i].neighbor[j_index]);
            let neighbor_index = usize::try_from(neighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let neighbor_atom = atoms.get(neighbor_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor_valence = i32::from(neighbor_atom.valence);
            let mut k = 0_i32;
            while k < neighbor_valence {
                let k_index = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if neighbor_atom.neighbor[k_index]
                    == i32::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)? as AT_NUMB
                {
                    break;
                }
                k = k.wrapping_add(1);
            }

            let mut bond_type = atoms[i].bond_type[j_index] & BOND_TYPE_MASK as u8;
            let bond_mark = atoms[i].bond_type[j_index] & !(BOND_TYPE_MASK as u8);
            if bond_type != BOND_TYPE_SINGLE as u8
                && bond_type != BOND_TYPE_DOUBLE as u8
                && bond_type != BOND_TYPE_TRIPLE as u8
            {
                bond_type = BOND_TYPE_SINGLE as u8;
                atoms[i].bond_type[j_index] = bond_mark | bond_type;
                num_changed_bonds = num_changed_bonds.wrapping_add(1);
            }

            let edge_flow;
            if neighbor > i32::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
                let mut atom_capacity = 0_i32;
                let mut atom_flow = 0_i32;
                let needs_flower_1 = AtomStcapStflow(
                    &atoms,
                    valence_atoms,
                    &restore_mode,
                    i as i32,
                    Some(&mut atom_capacity),
                    Some(&mut atom_flow),
                    None,
                    None,
                )?;
                atom_capacity = atom_capacity.wrapping_add(if needs_flower_1 != 0 {
                    0
                } else {
                    i32::from(valence_atoms[i].cInitFreeValences)
                });

                let mut neighbor_capacity = 0_i32;
                let mut neighbor_flow = 0_i32;
                let needs_flower_2 = AtomStcapStflow(
                    &atoms,
                    valence_atoms,
                    &restore_mode,
                    neighbor,
                    Some(&mut neighbor_capacity),
                    Some(&mut neighbor_flow),
                    None,
                    None,
                )?;
                neighbor_capacity = neighbor_capacity.wrapping_add(if needs_flower_2 != 0 {
                    0
                } else {
                    i32::from(valence_atoms[neighbor_index].cInitFreeValences)
                });

                let mut edge_capacity = 0_i32;
                let mut minimum_order = 0_i32;
                edge_flow = BondFlowMaxcapMinorder(
                    &atoms,
                    valence_atoms,
                    &restore_mode,
                    i as i32,
                    j,
                    Some(&mut edge_capacity),
                    Some(&mut minimum_order),
                    None,
                )?;

                let edge_index = usize::try_from(n_edges).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                {
                    let edge = heap
                        .slice_mut(network_value.edge)?
                        .get_mut(edge_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.neighbor1 = i as AT_NUMB;
                    edge.neighbor12 = (i as i32 ^ neighbor) as AT_NUMB;
                    edge.flow = edge_flow;
                    edge.flow0 = edge_flow;
                    edge.cap = edge_capacity;
                    edge.cap0 = edge_capacity;
                    edge.neigh_ord[0] = j as AT_NUMB;
                    edge.neigh_ord[1] = k as AT_NUMB;
                    edge.pass = 0;
                    edge.forbidden = 0;
                    if bond_type == BOND_TYPE_DOUBLE as u8 {
                        let mut m = 0_usize;
                        while m < MAX_NUM_STEREO_BONDS as usize && atoms[i].sb_parity[m] != 0 {
                            if i32::from(atoms[i].sb_ord[m]) == j {
                                edge.forbidden |= BNS_EDGE_FORBIDDEN_MASK as i8;
                                break;
                            }
                            m += 1;
                        }
                    }
                }

                let atom_incident = heap
                    .slice(network_value.vert.as_const())?
                    .get(i)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iedge;
                let neighbor_incident = heap
                    .slice(network_value.vert.as_const())?
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iedge;
                if !neighbor_incident.is_null() {
                    *heap
                        .slice_mut(atom_incident)?
                        .get_mut(j_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = n_edges;
                    *heap
                        .slice_mut(neighbor_incident)?
                        .get_mut(usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = n_edges;
                    n_edges = n_edges.wrapping_add(1);
                }
            } else {
                let neighbor_incident = heap
                    .slice(network_value.vert.as_const())?
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iedge;
                let edge_index = *heap
                    .slice(neighbor_incident.as_const())?
                    .get(usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                edge_flow = heap
                    .slice(network_value.edge.as_const())?
                    .get(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .flow;
            }
            stationary_flow = stationary_flow.wrapping_add(edge_flow);
            j = j.wrapping_add(1);
        }

        let vertex = heap
            .slice_mut(network_value.vert)?
            .get_mut(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        vertex.num_adj_edges = j as AT_NUMB;
        vertex.st_edge.flow = stationary_flow;
        vertex.st_edge.flow0 = stationary_flow;
        vertex.type_ = BNS_VERT_TYPE_ATOM as AT_NUMB;
        total_stationary_flow = total_stationary_flow.wrapping_add(vertex.st_edge.flow);
    }

    *pNum_changed_bonds = num_changed_bonds / 2;
    network_value.num_edges = n_edges;
    network_value.num_iedges = num_iedges;
    network_value.num_added_edges = 0;
    network_value.tot_st_cap = total_stationary_capacity;
    network_value.tot_st_flow = total_stationary_flow;

    if atom_count != 0 {
        heap.slice_mut(structure.at)?
            .get_mut(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&atoms);
    }
    *heap
        .slice_mut(network)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = network_value;
    Ok(network)
}

#[allow(non_snake_case)]
pub(crate) fn IncrZeroBondsAndClearEndpts(
    atoms: &mut [inp_ATOM],
    num_at: i32,
    component: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4152 IncrZeroBondsAndClearEndpts
    // INCHI✔️✔️: void IncrZeroBondsAndClearEndpts( inp_ATOM *at, int num_at, int iComponent )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j;
    // INCHI✔️✔️:     for (i = 0; i < num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         at[i].endpoint = 0;
    // INCHI✔️✔️:         at[i].component = iComponent;
    // INCHI✔️✔️:         for (j = 0; j < at[i].valence; j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (!at[i].bond_type[j])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[i].bond_type[j] = BOND_TYPE_SINGLE;
    // INCHI✔️✔️:                 at[i].chem_bonds_valence += BOND_TYPE_SINGLE;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: IncrZeroBondsAndClearEndpts

    if num_at <= 0 {
        return Ok(());
    }
    let count = usize::try_from(num_at).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let selected = atoms.get_mut(..count).ok_or(SourceHeapError::PointerOutOfBounds)?;
    for atom in selected {
        atom.endpoint = 0;
        atom.component = component as AT_NUMB;
        let mut neighbor_order = 0_i32;
        let valence = i32::from(atom.valence);
        while neighbor_order < valence {
            let index = usize::try_from(neighbor_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if atom.bond_type[index] == 0 {
                atom.bond_type[index] = BOND_TYPE_SINGLE as u8;
                atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(BOND_TYPE_SINGLE as i8);
            }
            neighbor_order = neighbor_order.wrapping_add(1);
        }
    }
    Ok(())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeOneInChIOutOfStrFromINChI2(
    heap: &mut SourceHeap,
    canonical_globals: &mut CANON_GLOBALS,
    clock: SourceMutPointer<INCHI_CLOCK>,
    input_parameters: &INPUT_PARMS,
    structure_data: &STRUCT_DATA,
    bns: &BN_STRUCT,
    structure: &mut StrFromINChI,
    at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    at3: SourceMutPointer<inp_ATOM>,
    valence_atoms: &[VAL_AT],
    groups: &ALL_TC_GROUPS,
    t_group_info: Option<&mut SourceTGroupInfoPointer>,
    at_norm: Option<&mut SourceMutPointer<inp_ATOM>>,
    at_prep: Option<&mut SourceMutPointer<inp_ATOM>>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5087 MakeOneInChIOutOfStrFromINChI2
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: MakeOneInChIOutOfStrFromINChI2
    // INCHI✔️✔️: int MakeOneInChIOutOfStrFromINChI2( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️✔️:                                     INCHI_CLOCK *ic,
    // INCHI✔️✔️:                                     ICHICONST INPUT_PARMS *ip_inp,
    // INCHI✔️✔️:                                     STRUCT_DATA *sd_inp,
    // INCHI✔️✔️:                                     BN_STRUCT *pBNS, StrFromINChI *pStruct,
    // INCHI✔️✔️:                                     inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
    // INCHI✔️✔️:                                     VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
    // INCHI✔️✔️:                                     T_GROUP_INFO **t_group_info,
    // INCHI✔️✔️:                                     inp_ATOM **at_norm, inp_ATOM **at_prep )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int ret;
    // INCHI✔️✔️:     INPUT_PARMS ip_loc, *ip;
    // INCHI✔️✔️:     STRUCT_DATA sd_loc, *sd;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     ip_loc = *ip_inp;
    // INCHI✔️✔️:     sd_loc = *sd_inp;
    // INCHI✔️✔️:     ip = &ip_loc;
    // INCHI✔️✔️:     sd = &sd_loc;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* create structure out of BNS */
    // INCHI✔️✔️:     memcpy(at2, at, ((long long)pStruct->num_atoms + (long long)pStruct->num_deleted_H) * sizeof(at2[0])); /* djb-rwth: cast operator added */
    // INCHI✔️✔️:     pStruct->at = at2;
    // INCHI✔️✔️:     ret = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    // INCHI✔️✔️:     pStruct->at = at;
    // INCHI✔️✔️:     if (ret < 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto exit_function;/*  out of RAM or other normalization problem */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pStruct->at = at;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     ret = MakeOneInChIOutOfStrFromINChI( pCG, ic, ip, sd, pStruct,
    // INCHI✔️✔️:                                          at2, at3, pTCGroups );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (ret < 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto exit_function;/*  out of RAM or other normalization problem */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (at_norm)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *at_norm = pStruct->pOne_norm_data[0]->at;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (at_prep)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (pStruct->pOne_norm_data[0]->bTautPreprocessed && pStruct->pOne_norm_data[0]->at_fixed_bonds)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *at_prep = pStruct->pOne_norm_data[0]->at_fixed_bonds;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         /* get preprocessed structure in case of Fixed-H */
    // INCHI✔️✔️:             if (pStruct->iMobileH == TAUT_NON && pStruct->pOne_norm_data[1] && pStruct->pOne_norm_data[1]->bTautPreprocessed)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *at_prep = pStruct->pOne_norm_data[1]->at_fixed_bonds;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *at_prep = NULL;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (t_group_info)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (pStruct->iMobileH == TAUT_YES &&
    // INCHI✔️✔️:              pStruct->One_ti.num_t_groups &&
    // INCHI✔️✔️:              pStruct->One_ti.t_group && pStruct->One_ti.nEndpointAtomNumber)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *t_group_info = &pStruct->One_ti;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *t_group_info = NULL;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_function:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: MakeOneInChIOutOfStrFromINChI2
    // END INCHI C FUNCTION: MakeOneInChIOutOfStrFromINChI2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MakeOneInChIOutOfStrFromINChI2
    // INCHI✔️✔️: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // INCHI✔️✔️: ICHICONST is const and the active memcpy is the GCC LP64 inp_ATOM layout.
    // END INCHI ACTIVE MACRO CONFIGURATION: MakeOneInChIOutOfStrFromINChI2

    let input_parameters_local = input_parameters.clone();
    let mut structure_data_local = structure_data.clone();
    structure_data_local = STRUCT_DATA::default();

    let count_i64 = i64::from(structure.num_atoms) + i64::from(structure.num_deleted_H);
    let count = usize::try_from(count_i64).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if count != 0 {
        copy_inp_atom_prefix(heap, at2, at, count)?;
    }

    structure.at = at2;
    let copy_result = CopyBnsToAtom(heap, structure, bns, valence_atoms, groups, 1);
    structure.at = at;
    let ret = copy_result?;
    if ret < 0 {
        return Ok(ret);
    }

    structure.at = at;
    let ret = MakeOneInChIOutOfStrFromINChI(
        heap,
        canonical_globals,
        clock,
        &input_parameters_local,
        &mut structure_data_local,
        structure,
        at2,
        at3,
        groups,
        clock_result,
    )?;
    if ret < 0 {
        return Ok(ret);
    }

    let normalized_holder = structure.pOne_norm_data[0];
    if let Some(output) = at_norm {
        *output = heap
            .slice(normalized_holder.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at;
    }

    if let Some(output) = at_prep {
        let normalized = heap
            .slice(normalized_holder.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if normalized.bTautPreprocessed != 0 && !normalized.at_fixed_bonds.is_null() {
            *output = normalized.at_fixed_bonds;
        } else if structure.iMobileH == TAUT_NON as i8 && !structure.pOne_norm_data[1].is_null() {
            let fixed_holder = heap
                .slice(structure.pOne_norm_data[1].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if fixed_holder.bTautPreprocessed != 0 {
                *output = fixed_holder.at_fixed_bonds;
            } else {
                *output = SourceMutPointer::null();
            }
        } else {
            *output = SourceMutPointer::null();
        }
    }

    if let Some(output) = t_group_info {
        if structure.iMobileH == TAUT_YES as i8
            && structure.One_ti.num_t_groups != 0
            && !structure.One_ti.t_group.is_null()
            && !structure.One_ti.nEndpointAtomNumber.is_null()
        {
            *output = SourceTGroupInfoPointer::StructureOne;
        } else {
            *output = SourceTGroupInfoPointer::Null;
        }
    }

    Ok(ret)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeOneInChIOutOfStrFromINChI(
    heap: &mut SourceHeap,
    canonical_globals: &mut CANON_GLOBALS,
    clock: SourceMutPointer<INCHI_CLOCK>,
    input_parameters: &INPUT_PARMS,
    structure_data: &mut STRUCT_DATA,
    structure: &mut StrFromINChI,
    at2: SourceMutPointer<inp_ATOM>,
    at3: SourceMutPointer<inp_ATOM>,
    groups: &ALL_TC_GROUPS,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5168 MakeOneInChIOutOfStrFromINChI
    // INCHI✔️❌: complete active source frame follows verbatim; typed SourceHeap pointer validation,
    // INCHI✔️❌: temporary cloning, and checked ownership transfers add overhead.
    /*
    int MakeOneInChIOutOfStrFromINChI( struct tagCANON_GLOBALS *pCG,
                                       INCHI_CLOCK *ic,
                                       ICHICONST INPUT_PARMS *ip,
                                       STRUCT_DATA *sd,
                                       StrFromINChI *pStruct,
                                       inp_ATOM *at2,
                                       inp_ATOM *at3,
                                       ALL_TC_GROUPS *pTCGroups )
    {

        INCHI_MODE     bTautFlags = ip->bTautFlags | TG_FLAG_H_ALREADY_REMOVED;
        INCHI_MODE     bTautFlagsDone = 0; /*(ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS]);*/
        INChI       *cur_INChI[TAUT_NUM];
        INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
        int           i, j, k;
        int           iComponent = pTCGroups->iComponent;
        int           len_at = pStruct->num_atoms + pStruct->num_deleted_H;
        int           num_atoms = pStruct->num_atoms;
        long          ulStructTime;

        INP_ATOM_DATA InpCurAtData;
        INP_ATOM_DATA *inp_cur_data;

        INP_ATOM_DATA InpNormAtData, InpNormTautData;
        INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */

        int            bOrigCoord = 0;
        int            num_at, ret = RI_ERR_PROGR;
        struct tagInchiTime ulMaxTime;

        T_GROUP_INFO *t_group_info = NULL;
        /* initialization */
        inp_cur_data = &InpCurAtData;
        inp_norm_data[TAUT_NON] = &InpNormAtData;
        inp_norm_data[TAUT_YES] = &InpNormTautData;

        memset( inp_cur_data, 0, sizeof( *inp_cur_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( inp_norm_data[TAUT_NON], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( inp_norm_data[TAUT_YES], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        ulStructTime = sd->ulStructTime;
        memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        /* deallocate old results */
        free_t_group_info( &pStruct->One_ti );
        for (k = 0; k < TAUT_NUM; k++)
        {
            Free_INChI( &pStruct->pOneINChI[k] );
            Free_INChI_Aux( &pStruct->pOneINChI_Aux[k] );
            if (pStruct->pOne_norm_data[k])
            {
                FreeInpAtomData( pStruct->pOne_norm_data[k] );
                inchi_free( pStruct->pOne_norm_data[k] );
                pStruct->pOne_norm_data[k] = NULL;
            }
            cur_INChI[k] = NULL;
            cur_INChI_Aux[k] = NULL;
        }

        memcpy(at3, at2, sizeof(at3[0]) * len_at);

        /* prepare the structure */
        IncrZeroBondsAndClearEndpts( at3, num_atoms, iComponent + 1 );

        CopySt2At( at3, pStruct->st, pStruct->num_atoms );

        FixUnkn0DStereoBonds( at3, pStruct->num_atoms );

        ret = ReconcileAllCmlBondParities( at3, pStruct->num_atoms, 0 );

        if (ret < 0)
        {
            goto exit_function;
        }

        if (0 < fix_odd_things( num_atoms, at3, 1, ip->bFixNonUniformDraw ))
        {
            if (sd->nErrorType < _IS_WARNING)
            {
                sd->nErrorType = _IS_WARNING;
            }
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
        }

        /* allocate and set parameters */
        inp_cur_data->at = at3;
        inp_cur_data->num_at = num_atoms;
        inp_cur_data->num_removed_H = pStruct->num_deleted_H;

        bTautFlagsDone &= ~( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE );

        if ((i = bNumHeterAtomHasIsotopicH( at3, num_atoms ))) /* djb-rwth: addressing LLVM warning */
        {
            if (i & 1)
            {
                bTautFlagsDone |= TG_FLAG_FOUND_ISOTOPIC_H_DONE;
            }
            if (i & 2)
            {
                bTautFlagsDone |= TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE;
            }
        }

        memset( &ulMaxTime, 0, sizeof( ulMaxTime ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        /*  allocate memory for non-tautimeric (k=0) and tautomeric (k=1) results */
        for (k = 0; k < TAUT_NUM; k++)
        {

            if (!pStruct->bMobileH || k == pStruct->bMobileH)
            {
                /* pStruct->bMobileH=0: k = 0, 1   => allow allocation of both Fixed-H and Mobile-H InChI
                   pStruct->bMobileH=1: k = 1 only => allow allocation of only Mobile-H InChI              */

                /* djb-rwth: introducing variables for correct nAllocMode expression */
                int nAM1 = 0, nAM2 = 0;
                int nAllocMode = 0;   /* copied from below 2024-09-01 DT */

                if (k == TAUT_YES)
                    nAM1 = REQ_MODE_TAUT;

                if (bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))
                    nAM2 = ip->nMode & REQ_MODE_ISO;

                nAllocMode = nAM1 | nAM2; /* djb-rwth: original sequence of bit-wise operations had to be rewritten */

                if ((k == TAUT_NON && ( ip->nMode & REQ_MODE_BASIC )) ||
                     (k == TAUT_YES && ( ip->nMode & REQ_MODE_TAUT ))) /* djb-rwth: addressing LLVM warnings */
                {
                    /*  alloc INChI and INChI_Aux only if ip->nMode allows this */
                    cur_INChI[k] = Alloc_INChI( inp_cur_data->at, inp_cur_data->num_at, &inp_cur_data->num_bonds,
                                                  &inp_cur_data->num_isotopic, nAllocMode );
                    cur_INChI_Aux[k] = Alloc_INChI_Aux( inp_cur_data->num_at,
                                                  inp_cur_data->num_isotopic, nAllocMode, bOrigCoord );
                    if (cur_INChI_Aux[k])
                    {
                        cur_INChI_Aux[k]->bIsIsotopic = inp_cur_data->num_isotopic;
                    }

                    /*  alloc memory for the output structure: non-tautomeric and tautomeric (for displaying) */
                    CreateInpAtomData( inp_norm_data[k], inp_cur_data->num_at + inp_cur_data->num_removed_H, k );

                    inp_norm_data[k]->num_removed_H = inp_cur_data->num_removed_H;
                }
                else
                {
                    FreeInpAtomData( inp_norm_data[k] );
                }
            }
            else
            {
                FreeInpAtomData( inp_norm_data[k] );
            }
        }

        k = pStruct->bMobileH;

        /* In case of Fixed-H we have to create InChI for both Fixed-H and Mobile-H */

        num_at = Create_INChI( pCG,
                               ic,
                               (INPUT_PARMS *) ip,
                               cur_INChI,
                               cur_INChI_Aux,
                               NULL,
                               inp_cur_data->at,
                               inp_norm_data,
                               inp_cur_data->num_at + inp_cur_data->num_removed_H,
                               ip->nMode,
                               &bTautFlags,
                               &bTautFlagsDone,
                               NULL /* &ulMaxTime*/,
                               &pStruct->One_ti,
                               sd->pStrErrStruct );

        SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, iComponent + 1 ); /*  normalization alters structure component number */

        /* Detect InChI errors */

        if (num_at < 0)
        {
            ret = num_at;
        }
        else if (cur_INChI[k] && cur_INChI[k]->nErrorCode)
        {
            ret = cur_INChI[k]->nErrorCode;
        }
        else if (cur_INChI_Aux[k] && cur_INChI_Aux[k]->nErrorCode)
        {
            ret = cur_INChI_Aux[k]->nErrorCode;
        }
        else
        {
            ret = 0;
        }

        /* Fill out the output */

        if (!ret) /* djb-rwth: fixing a NULL pointer dereference */
        {
            int bMobileH = pStruct->bMobileH;
            if (bMobileH == TAUT_NON &&
                 0 == cur_INChI[TAUT_NON]->nNumberOfAtoms &&
                 0 < cur_INChI[TAUT_YES]->nNumberOfAtoms)
            {
                /* tautomerism or H(+) removal/addition was not discovered */
                bMobileH = TAUT_YES;
            }

            if (cur_INChI[1])
                pStruct->nChargeRevrs = cur_INChI[TAUT_YES]->nTotalCharge; /* djb-rwth: fixing a NULL pointer dereference */

            pStruct->pOneINChI[0] = cur_INChI[bMobileH];
            pStruct->pOneINChI_Aux[0] = cur_INChI_Aux[bMobileH];
            pStruct->nOneINChI_bMobileH = bMobileH;
            cur_INChI[bMobileH] = NULL;  /* remove pointer to avoid deallocation at exit_function */
            cur_INChI_Aux[bMobileH] = NULL;  /* remove pointer to avoid deallocation at exit_function */

            pStruct->nNumRemovedProtons = ( pStruct->iMobileH == TAUT_YES ) ? pStruct->One_ti.tni.nNumRemovedProtons : 0;


            /* set correct t-group numbers to endpoints */

            t_group_info = &pStruct->One_ti;

            if (t_group_info->num_t_groups && t_group_info->t_group && t_group_info->nEndpointAtomNumber)
            {
                inp_ATOM     *at_norm = inp_norm_data[TAUT_YES]->at;
                int          num_at_norm = inp_norm_data[TAUT_YES]->num_at;
                for (i = 0; i < num_at_norm; i++)
                {
                    at_norm[i].endpoint = 0;
                }
                for (i = 0; i < t_group_info->num_t_groups; i++)
                {
                    k = t_group_info->t_group[i].nFirstEndpointAtNoPos;
                    /* add number of mobile (-) to the number of mobile H */
                    t_group_info->t_group[i].num[0] += t_group_info->t_group[i].num[1];
                    for (j = 0; j < t_group_info->t_group[i].nNumEndpoints; j++, k++)
                    {
                        at_norm[t_group_info->nEndpointAtomNumber[k]].endpoint = t_group_info->t_group[i].nGroupNumber;
                    }
                }
            }

            pStruct->pOne_norm_data[0] = (INP_ATOM_DATA *) inchi_malloc( sizeof( pStruct->pOne_norm_data[0][0] ) );

            if (pStruct->pOne_norm_data[0])
            {
                memcpy(pStruct->pOne_norm_data[0], inp_norm_data[bMobileH], sizeof(pStruct->pOne_norm_data[0][0]));
                memset( inp_norm_data[bMobileH], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            }
            else
            {
                ret = RI_ERR_ALLOC;
            }

            if (bMobileH == TAUT_NON && cur_INChI[TAUT_YES]->nNumberOfAtoms > 0)
            {
                int bMobileHalt = ALT_TAUT( bMobileH ); /* = TAUT_YES */
                pStruct->pOneINChI[1] = cur_INChI[bMobileHalt];
                pStruct->pOneINChI_Aux[1] = cur_INChI_Aux[bMobileHalt];
                cur_INChI[bMobileHalt] = NULL;
                cur_INChI_Aux[bMobileHalt] = NULL;
                pStruct->pOne_norm_data[1] = (INP_ATOM_DATA *) inchi_malloc( sizeof( pStruct->pOne_norm_data[0][0] ) );
                if (pStruct->pOne_norm_data[1])
                {
                    memcpy(pStruct->pOne_norm_data[1], inp_norm_data[bMobileHalt], sizeof(pStruct->pOne_norm_data[0][0]));
                    memset( inp_norm_data[bMobileHalt], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                }
                else
                {
                    ret = RI_ERR_ALLOC;
                }
            }
        }
        else
        {
    #if ( bRELEASE_VERSION != 1 )
    #ifndef TARGET_API_LIB
            fprintf( stdout, "ERROR: Create_INChI returned %d\n", ret );
    #endif
    #endif
        }

    exit_function:
        /* deallocate unused */
        for (k = 0; k < TAUT_NUM; k++)
        {
            Free_INChI( &cur_INChI[k] );
            Free_INChI_Aux( &cur_INChI_Aux[k] );
            FreeInpAtomData( inp_norm_data[k] );
        }
        sd->ulStructTime = ulStructTime;

        return ret;
    }
        */
    // END INCHI C FUNCTION: MakeOneInChIOutOfStrFromINChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MakeOneInChIOutOfStrFromINChI
    // INCHI✔️❌: COMPILE_ANSI_ONLY=1
    // INCHI✔️❌: TARGET_API_LIB=1
    // INCHI✔️❌: GCC/Linux branch
    // INCHI✔️❌: #define inchi_malloc malloc
    // INCHI✔️❌: #define inchi_free(X) do{ if(X) free(X); }while(0)
    // INCHI✔️❌: The non-TARGET_API_LIB diagnostic fprintf branch is inactive.
    // END INCHI ACTIVE MACRO CONFIGURATION: MakeOneInChIOutOfStrFromINChI

    fn allocate_atom_data_holder(heap: &mut SourceHeap) -> Result<SourceMutPointer<INP_ATOM_DATA>, SourceHeapError> {
        match inchi_calloc::<INP_ATOM_DATA>(heap, 1, std::mem::size_of::<INP_ATOM_DATA>() as u64) {
            Ok(pointer) => Ok(pointer),
            Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
            Err(error) => Err(error),
        }
    }

    let saved_structure_time = structure_data.ulStructTime;
    *structure_data = STRUCT_DATA::default();

    let mut current_inchi = [SourceMutPointer::null(); TAUT_NUM as usize];
    let mut current_aux = [SourceMutPointer::null(); TAUT_NUM as usize];
    let mut current_data = INP_ATOM_DATA::default();
    let mut normalized_data: [INP_ATOM_DATA; TAUT_NUM as usize] = std::array::from_fn(|_| INP_ATOM_DATA::default());
    let mut ret = RI_ERR_PROGR;

    let processing = (|| -> Result<(), SourceHeapError> {
        free_t_group_info(heap, Some(&mut structure.One_ti))?;
        for tautomer in 0..TAUT_NUM as usize {
            Free_INChI(heap, &mut structure.pOneINChI[tautomer])?;
            Free_INChI_Aux(heap, &mut structure.pOneINChI_Aux[tautomer])?;
            let old_normalized = structure.pOne_norm_data[tautomer];
            if !old_normalized.is_null() {
                let mut owned_data = heap
                    .slice(old_normalized.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                FreeInpAtomData(heap, &mut owned_data)?;
                inchi_free(heap, old_normalized)?;
                structure.pOne_norm_data[tautomer] = SourceMutPointer::null();
            }
        }

        let len_at = structure.num_atoms.wrapping_add(structure.num_deleted_H);
        let copied_count = usize::try_from(len_at).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let copied_atoms = heap
            .slice(at2.as_const())?
            .get(..copied_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        heap.slice_mut(at3)?
            .get_mut(..copied_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&copied_atoms);

        let component_number = groups.iComponent.wrapping_add(1);
        IncrZeroBondsAndClearEndpts(heap.slice_mut(at3)?, structure.num_atoms, component_number)?;

        let stereo = if structure.st.is_null() {
            None
        } else {
            let count = usize::try_from(structure.num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            Some(
                heap.slice(structure.st.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec(),
            )
        };
        CopySt2At(heap.slice_mut(at3)?, stereo.as_deref(), structure.num_atoms)?;
        let _ = FixUnkn0DStereoBonds(heap, at3, structure.num_atoms)?;
        ret = ReconcileAllCmlBondParities(heap, at3, structure.num_atoms, 0)?;
        if ret < 0 {
            return Ok(());
        }

        if fix_odd_things(
            structure.num_atoms,
            heap.slice_mut(at3)?,
            1,
            input_parameters.bFixNonUniformDraw,
        )? > 0
        {
            if structure_data.nErrorType < _IS_WARNING as i32 {
                structure_data.nErrorType = _IS_WARNING as i32;
            }
            structure_data.bTautFlagsDone[INCHI_BAS as usize] |= u64::from(TG_FLAG_FIX_ODD_THINGS_DONE);
        }

        current_data.at = at3;
        current_data.num_at = structure.num_atoms;
        current_data.num_removed_H = structure.num_deleted_H;

        let mut taut_flags = input_parameters.bTautFlags | u64::from(TG_FLAG_H_ALREADY_REMOVED);
        let mut taut_flags_done = 0_u64;
        taut_flags_done &= !(u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE) | u64::from(TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE));
        let isotope_state = bNumHeterAtomHasIsotopicH(heap.slice(at3.as_const())?, structure.num_atoms)?;
        if isotope_state != 0 {
            if isotope_state & 1 != 0 {
                taut_flags_done |= u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE);
            }
            if isotope_state & 2 != 0 {
                taut_flags_done |= u64::from(TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
            }
        }

        for tautomer in 0..TAUT_NUM as usize {
            let tautomer_i32 = tautomer as i32;
            if structure.bMobileH == 0 || tautomer_i32 == i32::from(structure.bMobileH) {
                let taut_mode = if tautomer_i32 == TAUT_YES as i32 {
                    REQ_MODE_TAUT
                } else {
                    0
                };
                let isotope_mode = if taut_flags_done
                    & (u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE) | u64::from(TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))
                    != 0
                {
                    (input_parameters.nMode & u64::from(REQ_MODE_ISO)) as u32
                } else {
                    0
                };
                let allocation_mode = (taut_mode | isotope_mode) as i32;
                let requested = (tautomer_i32 == TAUT_NON as i32
                    && input_parameters.nMode & u64::from(REQ_MODE_BASIC) != 0)
                    || (tautomer_i32 == TAUT_YES as i32 && input_parameters.nMode & u64::from(REQ_MODE_TAUT) != 0);
                if requested {
                    let atom_count =
                        usize::try_from(current_data.num_at).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    let atoms = heap
                        .slice(current_data.at.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec();
                    current_inchi[tautomer] = Alloc_INChI(
                        heap,
                        &atoms,
                        current_data.num_at,
                        &mut current_data.num_bonds,
                        &mut current_data.num_isotopic,
                        allocation_mode,
                    )?;
                    current_aux[tautomer] =
                        Alloc_INChI_Aux(heap, current_data.num_at, current_data.num_isotopic, allocation_mode, 0)?;
                    if !current_aux[tautomer].is_null() {
                        heap.slice_mut(current_aux[tautomer])?[0].bIsIsotopic = current_data.num_isotopic;
                    }
                    let _ = CreateInpAtomData(
                        heap,
                        &mut normalized_data[tautomer],
                        current_data.num_at.wrapping_add(current_data.num_removed_H),
                        tautomer_i32,
                    )?;
                    normalized_data[tautomer].num_removed_H = current_data.num_removed_H;
                } else {
                    FreeInpAtomData(heap, &mut normalized_data[tautomer])?;
                }
            } else {
                FreeInpAtomData(heap, &mut normalized_data[tautomer])?;
            }
        }

        let selected = usize::try_from(structure.bMobileH).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let create_result = Create_INChI(
            heap,
            canonical_globals,
            clock,
            input_parameters,
            current_inchi,
            current_aux,
            None,
            current_data.at,
            &mut normalized_data,
            current_data.num_at.wrapping_add(current_data.num_removed_H),
            input_parameters.nMode,
            &mut taut_flags,
            &mut taut_flags_done,
            SourceMutPointer::null(),
            Some(&mut structure.One_ti),
            Some(&mut structure_data.pStrErrStruct),
            clock_result,
        );
        let num_at = create_result?;
        SetConnectedComponentNumber(heap.slice_mut(current_data.at)?, current_data.num_at, component_number)?;

        ret = if num_at < 0 {
            num_at
        } else {
            let inchi_error = if current_inchi
                .get(selected)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .is_null()
            {
                0
            } else {
                heap.slice(current_inchi[selected].as_const())?[0].nErrorCode
            };
            let aux_error = if current_aux
                .get(selected)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .is_null()
            {
                0
            } else {
                heap.slice(current_aux[selected].as_const())?[0].nErrorCode
            };
            if inchi_error != 0 {
                inchi_error
            } else if aux_error != 0 {
                aux_error
            } else {
                0
            }
        };
        if ret != 0 {
            return Ok(());
        }

        let mut mobile_h = i32::from(structure.bMobileH);
        if mobile_h == TAUT_NON as i32 {
            let fixed = current_inchi[TAUT_NON as usize];
            let taut = current_inchi[TAUT_YES as usize];
            if heap.slice(fixed.as_const())?[0].nNumberOfAtoms == 0
                && heap.slice(taut.as_const())?[0].nNumberOfAtoms > 0
            {
                mobile_h = TAUT_YES as i32;
            }
        }

        let taut_inchi = current_inchi[TAUT_YES as usize];
        if !taut_inchi.is_null() {
            structure.nChargeRevrs = heap.slice(taut_inchi.as_const())?[0].nTotalCharge;
        }

        let selected = usize::try_from(mobile_h).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        structure.pOneINChI[0] = current_inchi[selected];
        structure.pOneINChI_Aux[0] = current_aux[selected];
        structure.nOneINChI_bMobileH = mobile_h;
        current_inchi[selected] = SourceMutPointer::null();
        current_aux[selected] = SourceMutPointer::null();
        structure.nNumRemovedProtons = if i32::from(structure.iMobileH) == TAUT_YES as i32 {
            i32::from(structure.One_ti.tni.nNumRemovedProtons)
        } else {
            0
        };

        if structure.One_ti.num_t_groups != 0
            && !structure.One_ti.t_group.is_null()
            && !structure.One_ti.nEndpointAtomNumber.is_null()
        {
            let normalized_atoms = normalized_data[TAUT_YES as usize].at;
            let normalized_count = usize::try_from(normalized_data[TAUT_YES as usize].num_at)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice_mut(normalized_atoms)?
                .get_mut(..normalized_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iter_mut()
                .for_each(|atom| atom.endpoint = 0);

            let group_count =
                usize::try_from(structure.One_ti.num_t_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            for group_index in 0..group_count {
                let (first_endpoint, endpoint_count, group_number) = {
                    let group = heap
                        .slice_mut(structure.One_ti.t_group)?
                        .get_mut(group_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    group.num[0] = group.num[0].wrapping_add(group.num[1]);
                    (
                        usize::from(group.nFirstEndpointAtNoPos),
                        usize::from(group.nNumEndpoints),
                        group.nGroupNumber,
                    )
                };
                let endpoint_indices = heap
                    .slice(structure.One_ti.nEndpointAtomNumber.as_const())?
                    .get(
                        first_endpoint
                            ..first_endpoint
                                .checked_add(endpoint_count)
                                .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let atoms = heap.slice_mut(normalized_atoms)?;
                for atom_index in endpoint_indices {
                    atoms
                        .get_mut(usize::from(atom_index))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint = group_number;
                }
            }
        }

        structure.pOne_norm_data[0] = allocate_atom_data_holder(heap)?;
        if !structure.pOne_norm_data[0].is_null() {
            heap.slice_mut(structure.pOne_norm_data[0])?[0] = normalized_data[selected].clone();
            normalized_data[selected] = INP_ATOM_DATA::default();
        } else {
            ret = RI_ERR_ALLOC;
        }

        if mobile_h == TAUT_NON as i32 {
            let alternative_inchi = current_inchi[TAUT_YES as usize];
            if heap.slice(alternative_inchi.as_const())?[0].nNumberOfAtoms > 0 {
                let alternative = TAUT_YES as usize;
                structure.pOneINChI[1] = current_inchi[alternative];
                structure.pOneINChI_Aux[1] = current_aux[alternative];
                current_inchi[alternative] = SourceMutPointer::null();
                current_aux[alternative] = SourceMutPointer::null();
                structure.pOne_norm_data[1] = allocate_atom_data_holder(heap)?;
                if !structure.pOne_norm_data[1].is_null() {
                    heap.slice_mut(structure.pOne_norm_data[1])?[0] = normalized_data[alternative].clone();
                    normalized_data[alternative] = INP_ATOM_DATA::default();
                } else {
                    ret = RI_ERR_ALLOC;
                }
            }
        }
        Ok(())
    })();

    let cleanup = (|| -> Result<(), SourceHeapError> {
        for tautomer in 0..TAUT_NUM as usize {
            Free_INChI(heap, &mut current_inchi[tautomer])?;
            Free_INChI_Aux(heap, &mut current_aux[tautomer])?;
            FreeInpAtomData(heap, &mut normalized_data[tautomer])?;
        }
        Ok(())
    })();
    structure_data.ulStructTime = saved_structure_time;
    processing?;
    cleanup?;
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn AllocBfsQueue(
    heap: &mut SourceHeap,
    queue: &mut BFS_Q,
    num_at: i32,
    min_ring_size: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4815 AllocBfsQueue
    // INCHI✔️❌: int AllocBfsQueue( BFS_Q *pQ, int num_at, int min_ring_size )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     switch (num_at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         case BFS_Q_FREE:
    // INCHI✔️❌:             if (pQ->q)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pQ->q = QueueDelete( pQ->q );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pQ->nAtomLevel)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( pQ->nAtomLevel );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pQ->cSource)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( pQ->cSource );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* fall through */
    // INCHI✔️❌:         case BFS_Q_CLEAR:
    // INCHI✔️❌:             memset( pQ, 0, sizeof( *pQ ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         default:
    // INCHI✔️❌:             if (num_at <= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = RI_ERR_PROGR;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num_at > pQ->num_at)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pQ->num_at)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     AllocBfsQueue( pQ, BFS_Q_FREE, 0 );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pQ->q = QueueCreate( num_at + 1, sizeof( qInt ) );
    // INCHI✔️❌:                 pQ->nAtomLevel = (AT_RANK*) inchi_calloc( num_at, sizeof( pQ->nAtomLevel[0] ) );
    // INCHI✔️❌:                 pQ->cSource = (S_CHAR *) inchi_calloc( num_at, sizeof( pQ->cSource[0] ) );
    // INCHI✔️❌:                 if (!pQ->q || !pQ->cSource || !pQ->nAtomLevel)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = RI_ERR_ALLOC;
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pQ->num_at = num_at;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             pQ->min_ring_size = min_ring_size;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AllocBfsQueue
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AllocBfsQueue
    // INCHI✔️❌: #define BFS_Q_CLEAR (-1)
    // INCHI✔️❌: #define BFS_Q_FREE (-2)
    // INCHI✔️❌: typedef AT_RANK qInt;
    // INCHI✔️❌: #define inchi_calloc calloc
    // INCHI✔️❌: #define inchi_free(X) do{ if(X) free(X); }while(0)
    // END INCHI ACTIVE MACRO CONFIGURATION: AllocBfsQueue

    match num_at {
        BFS_Q_FREE => {
            if !queue.q.is_null() {
                queue.q = QueueDelete(heap, queue.q)?;
            }
            if !queue.nAtomLevel.is_null() {
                inchi_free(heap, queue.nAtomLevel)?;
            }
            if !queue.cSource.is_null() {
                inchi_free(heap, queue.cSource)?;
            }
            *queue = BFS_Q::default();
            Ok(0)
        }
        BFS_Q_CLEAR => {
            *queue = BFS_Q::default();
            Ok(0)
        }
        _ => {
            if num_at <= 0 {
                return Ok(RI_ERR_PROGR);
            }
            if num_at > queue.num_at {
                if queue.num_at != 0 {
                    let _ = AllocBfsQueue(heap, queue, BFS_Q_FREE, 0)?;
                }
                queue.q = QueueCreate(
                    heap,
                    num_at.wrapping_add(1),
                    std::mem::size_of::<crate::source_types::qInt>() as i32,
                )?;
                queue.nAtomLevel = match inchi_calloc::<AT_RANK>(heap, num_at as u64, 2) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
                queue.cSource = match inchi_calloc::<i8>(heap, num_at as u64, 1) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
                if queue.q.is_null() || queue.cSource.is_null() || queue.nAtomLevel.is_null() {
                    return Ok(RI_ERR_ALLOC);
                }
                queue.num_at = num_at;
            }
            queue.min_ring_size = min_ring_size as AT_RANK;
            Ok(0)
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn SetStCapFlow(
    vertex: &mut BNS_VERTEX,
    total_stationary_flow: &mut i32,
    total_stationary_capacity: &mut i32,
    capacity: i32,
    flow: i32,
) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3231 SetStCapFlow
    // INCHI✔️✔️: void SetStCapFlow( BNS_VERTEX *vert_ficpoint, int *tot_st_flow, int *tot_st_cap, int cap, int flow )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     *tot_st_flow += flow - vert_ficpoint->st_edge.flow;
    // INCHI✔️✔️:     vert_ficpoint->st_edge.flow = flow;
    // INCHI✔️✔️:     *tot_st_cap += cap - vert_ficpoint->st_edge.cap;
    // INCHI✔️✔️:     vert_ficpoint->st_edge.cap = cap;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     vert_ficpoint->st_edge.flow0 = vert_ficpoint->st_edge.flow;
    // INCHI✔️✔️:     vert_ficpoint->st_edge.cap0 = vert_ficpoint->st_edge.cap;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SetStCapFlow

    *total_stationary_flow = total_stationary_flow.wrapping_add(flow.wrapping_sub(vertex.st_edge.flow));
    vertex.st_edge.flow = flow;
    *total_stationary_capacity = total_stationary_capacity.wrapping_add(capacity.wrapping_sub(vertex.st_edge.cap));
    vertex.st_edge.cap = capacity;
    vertex.st_edge.flow0 = vertex.st_edge.flow;
    vertex.st_edge.cap0 = vertex.st_edge.cap;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{BOND_TYPE_DOUBLE, T_GROUP};

    #[test]
    fn copy_inp_atom_prefix_preserves_source_memcpy_boundaries() {
        let atom = |number| inp_ATOM {
            orig_at_number: number,
            ..inp_ATOM::default()
        };

        let mut heap = SourceHeap::default();
        let source = heap.allocate_model_storage(vec![atom(1), atom(2), atom(3)]).unwrap();
        let destination = heap.allocate_model_storage(vec![atom(10), atom(20), atom(30)]).unwrap();
        copy_inp_atom_prefix(&mut heap, destination, source, 2).unwrap();
        assert_eq!(
            heap.slice(source.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.orig_at_number)
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );
        assert_eq!(
            heap.slice(destination.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.orig_at_number)
                .collect::<Vec<_>>(),
            [1, 2, 30]
        );

        copy_inp_atom_prefix(&mut heap, source, source, 3).unwrap();
        assert_eq!(
            heap.slice(source.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.orig_at_number)
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );
        assert_eq!(
            copy_inp_atom_prefix(&mut heap, source, source, 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn copy_inp_atom_prefix_preserves_modeled_overlap_and_error_order() {
        let atom = |number| inp_ATOM {
            orig_at_number: number,
            ..inp_ATOM::default()
        };

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![atom(1), atom(2), atom(3), atom(4)])
            .unwrap();
        let shifted = atoms.offset(1).unwrap();
        copy_inp_atom_prefix(&mut heap, shifted, atoms, 3).unwrap();
        assert_eq!(
            heap.slice(atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.orig_at_number)
                .collect::<Vec<_>>(),
            [1, 1, 2, 3]
        );

        let short_source = heap.allocate_model_storage(vec![atom(5)]).unwrap();
        assert_eq!(
            copy_inp_atom_prefix(&mut heap, SourceMutPointer::null(), short_source, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            copy_inp_atom_prefix(&mut heap, SourceMutPointer::null(), short_source, 1),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            copy_inp_atom_prefix(&mut heap, shifted, atoms, 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr1__getdeltachargefromvf__line_4229() {
        let mut heap = SourceHeap::default();
        let edges = heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    cap: 3,
                    flow: 1,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    flow: 1,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    cap: 1,
                    flow: 1,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    flow: 2,
                    ..BNS_EDGE::default()
                },
            ])
            .unwrap();
        let bns = BN_STRUCT {
            num_atoms: 3,
            edge: edges,
            ..BN_STRUCT::default()
        };
        let valence_atoms = [
            VAL_AT {
                cInitCharge: -1,
                nCPlusGroupEdge: 1,
                nCMinusGroupEdge: 2,
                ..VAL_AT::default()
            },
            VAL_AT {
                cInitCharge: 1,
                nCPlusGroupEdge: 3,
                ..VAL_AT::default()
            },
            VAL_AT {
                cInitCharge: 2,
                nCMinusGroupEdge: 4,
                ..VAL_AT::default()
            },
        ];

        let base = VF {
            type_: BNS_VERT_TYPE_C_GROUP as i32,
            e_In: -1,
            e_Out: 0,
            delta_Out: 1,
            ..VF::default()
        };
        let mut wrong_type = VF {
            type_: BNS_VERT_TYPE_TGROUP as i32,
            ..base.clone()
        };
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut wrong_type),
            Ok(0)
        );
        assert_eq!(wrong_type.bUsed, 0);

        let mut super_group = VF {
            type_: (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_SUPER_TGROUP) as i32,
            ..base.clone()
        };
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut super_group),
            Ok(0)
        );
        assert_eq!(super_group.bUsed, 0);

        for mut unavailable in [
            VF {
                bUsed: VF_USED_OUT as i32,
                ..base.clone()
            },
            VF {
                delta_Out: 0,
                ..base.clone()
            },
            VF {
                e_Out: -1,
                ..base.clone()
            },
        ] {
            let before = unavailable.clone();
            assert_eq!(
                GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut unavailable),
                Ok(0)
            );
            assert_eq!(unavailable, before);
        }

        let mut no_match = VF {
            e_Out: 19,
            ..base.clone()
        };
        assert_eq!(GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut no_match), Ok(0));
        assert_eq!(no_match.bUsed, 0);

        let empty_bns = BN_STRUCT {
            num_atoms: 0,
            edge: edges,
            ..BN_STRUCT::default()
        };
        let mut empty = base.clone();
        assert_eq!(GetDeltaChargeFromVF(&heap, &empty_bns, &[], &mut empty), Ok(0));
        assert_eq!(empty.bUsed, 0);

        let mut zero_to_charged = base.clone();
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut zero_to_charged),
            Ok(1)
        );
        assert_eq!(zero_to_charged.bUsed, VF_USED_OUT as i32);

        let mut both_edges = VF {
            type_: BNS_VERT_TYPE_C_GROUP as i32,
            e_In: 1,
            e_Out: 0,
            delta_In: -2,
            delta_Out: 1,
            ..VF::default()
        };
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut both_edges),
            Ok(1)
        );
        assert_eq!(both_edges.bUsed, (VF_USED_IN | VF_USED_OUT) as i32);

        let mut cancelling_edges = VF {
            delta_In: -1,
            ..both_edges.clone()
        };
        cancelling_edges.bUsed = 0;
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut cancelling_edges),
            Ok(0)
        );
        assert_eq!(cancelling_edges.bUsed, (VF_USED_IN | VF_USED_OUT) as i32);

        let mut charged_to_zero = VF {
            type_: BNS_VERT_TYPE_C_GROUP as i32,
            e_In: -1,
            e_Out: 2,
            delta_Out: 1,
            ..VF::default()
        };
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut charged_to_zero),
            Ok(-1)
        );
        assert_eq!(charged_to_zero.bUsed, (VF_USED_IN | VF_USED_OUT) as i32);

        let mut remains_charged = VF {
            delta_Out: 2,
            ..charged_to_zero.clone()
        };
        remains_charged.bUsed = 0;
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut remains_charged),
            Ok(0)
        );
        assert_eq!(remains_charged.bUsed, (VF_USED_IN | VF_USED_OUT) as i32);

        let mut negative_group = VF {
            type_: (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE) as i32,
            e_In: -1,
            e_Out: 3,
            delta_Out: 1,
            ..VF::default()
        };
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut negative_group),
            Ok(1)
        );
        assert_eq!(negative_group.bUsed, (VF_USED_IN | VF_USED_OUT) as i32);

        let mut output_already_used = VF {
            type_: (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE) as i32,
            bUsed: VF_USED_OUT as i32,
            e_In: 1,
            delta_In: 1,
            ..base
        };
        assert_eq!(
            GetDeltaChargeFromVF(&heap, &bns, &valence_atoms, &mut output_already_used),
            Ok(1)
        );
        assert_eq!(output_already_used.bUsed, (VF_USED_IN | VF_USED_OUT) as i32);
    }

    #[test]
    fn source_port__ichirvr1__evaluatechargechanges__line_4320() {
        fn fixture(
            heap: &mut SourceHeap,
            types: &[u16],
            path_vertices: &[usize],
            delta: i32,
            end_override: Option<i32>,
            num_atoms: i32,
        ) -> BN_STRUCT {
            let mut adjacency = vec![Vec::<i32>::new(); types.len()];
            let mut edges = Vec::new();
            let mut path_neighbor_orders = Vec::new();
            for pair in path_vertices.windows(2) {
                let first = pair[0];
                let second = pair[1];
                let edge_index = edges.len() as i32;
                path_neighbor_orders.push(adjacency[first].len() as u16);
                adjacency[first].push(edge_index);
                adjacency[second].push(edge_index);
                edges.push(BNS_EDGE {
                    neighbor1: first as u16,
                    neighbor12: (first as u16) ^ (second as u16),
                    ..BNS_EDGE::default()
                });
            }
            let adjacency_lengths = adjacency.iter().map(|indices| indices.len() as u16).collect::<Vec<_>>();
            let incident = adjacency
                .into_iter()
                .map(|indices| heap.allocate_model_storage(indices).unwrap())
                .collect::<Vec<_>>();
            let vertices = heap
                .allocate_model_storage(
                    types
                        .iter()
                        .enumerate()
                        .map(|(index, type_)| BNS_VERTEX {
                            type_: *type_,
                            num_adj_edges: adjacency_lengths[index],
                            max_adj_edges: adjacency_lengths[index],
                            iedge: incident[index],
                            ..BNS_VERTEX::default()
                        })
                        .collect(),
                )
                .unwrap();
            let edges = heap.allocate_model_storage(edges).unwrap();
            let mut path = vec![BNS_ALT_PATH::default(); 5 + path_neighbor_orders.len()];
            path[tagAltPathConst_iALTP_FLOW as usize].set_flow(0, delta);
            path[tagAltPathConst_iALTP_PATH_LEN as usize].set_number(path_neighbor_orders.len() as i32);
            path[tagAltPathConst_iALTP_START_ATOM as usize].set_number(path_vertices[0] as i32);
            path[tagAltPathConst_iALTP_END_ATOM as usize]
                .set_number(end_override.unwrap_or_else(|| *path_vertices.last().unwrap() as i32));
            for (index, neighbor_order) in path_neighbor_orders.into_iter().enumerate() {
                path[tagAltPathConst_iALTP_HDR_LEN as usize + index].set_ineigh(0, neighbor_order);
            }
            let path = heap.allocate_model_storage(path).unwrap();
            let mut network = BN_STRUCT {
                num_atoms,
                num_vertices: types.len() as i32,
                num_edges: path_vertices.len().saturating_sub(1) as i32,
                vert: vertices,
                edge: edges,
                num_altp: 1,
                ..BN_STRUCT::default()
            };
            network.altp[0] = path;
            network
        }

        let mut empty_h = 7;
        let mut empty_charge = -9;
        let mut empty_visited = 11;
        let mut empty_heap = SourceHeap::default();
        let retained_alt_path = empty_heap
            .allocate_model_storage(vec![BNS_ALT_PATH::default()])
            .unwrap();
        let mut empty = BN_STRUCT {
            alt_path: retained_alt_path,
            ..BN_STRUCT::default()
        };
        assert_eq!(
            EvaluateChargeChanges(
                &empty_heap,
                &mut empty,
                &[],
                &mut empty_h,
                &mut empty_charge,
                &mut empty_visited,
            ),
            Ok(0)
        );
        assert_eq!((empty_h, empty_charge, empty_visited), (0, 0, 0));
        assert_eq!(empty.alt_path, retained_alt_path);

        let mut atom_heap = SourceHeap::default();
        let mut atoms = fixture(
            &mut atom_heap,
            &[BNS_VERT_TYPE_ATOM as u16, BNS_VERT_TYPE_ATOM as u16],
            &[0, 1],
            1,
            None,
            2,
        );
        let mut delta_h = -1;
        let mut delta_charge = -1;
        let mut visited = -1;
        assert_eq!(
            EvaluateChargeChanges(
                &atom_heap,
                &mut atoms,
                &[VAL_AT::default(), VAL_AT::default()],
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(0)
        );
        assert_eq!((delta_h, delta_charge, visited), (0, 0, 1));
        assert_eq!(atoms.alt_path, atoms.altp[0]);

        let mut wrong_end_heap = SourceHeap::default();
        let mut wrong_end = fixture(
            &mut wrong_end_heap,
            &[BNS_VERT_TYPE_ATOM as u16, BNS_VERT_TYPE_ATOM as u16],
            &[0, 1],
            1,
            Some(0),
            2,
        );
        assert_eq!(
            EvaluateChargeChanges(
                &wrong_end_heap,
                &mut wrong_end,
                &[VAL_AT::default(), VAL_AT::default()],
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(BNS_PROGRAM_ERR)
        );
        assert_eq!((delta_h, delta_charge, visited), (0, 0, 1));

        let mut hydrogen_heap = SourceHeap::default();
        let mut hydrogen = fixture(
            &mut hydrogen_heap,
            &[BNS_VERT_TYPE_ATOM as u16, BNS_VERT_TYPE_TGROUP as u16],
            &[0, 1],
            -3,
            None,
            1,
        );
        assert_eq!(
            EvaluateChargeChanges(
                &hydrogen_heap,
                &mut hydrogen,
                &[VAL_AT::default()],
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(0)
        );
        assert_eq!((delta_h, delta_charge, visited), (3, 0, 0));

        let second_path = hydrogen_heap
            .allocate_model_storage(hydrogen_heap.slice(hydrogen.altp[0].as_const()).unwrap().to_vec())
            .unwrap();
        hydrogen.altp[1] = second_path;
        hydrogen.num_altp = 2;
        assert_eq!(
            EvaluateChargeChanges(
                &hydrogen_heap,
                &mut hydrogen,
                &[VAL_AT::default()],
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(0)
        );
        assert_eq!((delta_h, delta_charge, visited), (6, 0, 0));
        assert_eq!(hydrogen.alt_path, hydrogen.altp[0]);

        for (initial_charge, expected) in [(0_i8, 1), (1_i8, -1)] {
            let mut charge_heap = SourceHeap::default();
            let mut charge = fixture(
                &mut charge_heap,
                &[BNS_VERT_TYPE_ATOM as u16, BNS_VERT_TYPE_C_GROUP as u16],
                &[0, 1],
                1,
                None,
                1,
            );
            let valence = [VAL_AT {
                cInitCharge: initial_charge,
                nCPlusGroupEdge: 1,
                ..VAL_AT::default()
            }];
            assert_eq!(
                EvaluateChargeChanges(
                    &charge_heap,
                    &mut charge,
                    &valence,
                    &mut delta_h,
                    &mut delta_charge,
                    &mut visited,
                ),
                Ok(0)
            );
            assert_eq!((delta_h, delta_charge, visited), (0, expected, 0));
        }

        let mut discarded_heap = SourceHeap::default();
        let mut discarded = fixture(
            &mut discarded_heap,
            &[BNS_VERT_TYPE_C_GROUP as u16, BNS_VERT_TYPE_ATOM as u16],
            &[0, 1],
            1,
            None,
            1,
        );
        let discarded_valence = [VAL_AT {
            nCPlusGroupEdge: 1,
            ..VAL_AT::default()
        }];
        assert_eq!(
            EvaluateChargeChanges(
                &discarded_heap,
                &mut discarded,
                &discarded_valence,
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(0)
        );
        assert_eq!((delta_h, delta_charge, visited), (0, 1, 1));
    }

    #[test]
    fn source_port__ichirvr1__runbnstestonce__line_4507() {
        fn data(heap: &mut SourceHeap, slots: usize) -> BN_DATA {
            BN_DATA {
                BasePtr: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
                SwitchEdge: heap.allocate_model_storage(vec![[NO_VERTEX, -1]; slots]).unwrap(),
                Tree: heap.allocate_model_storage(vec![0_i8; slots]).unwrap(),
                ScanQ: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
                Pu: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
                Pv: heap.allocate_model_storage(vec![NO_VERTEX; slots]).unwrap(),
                max_num_vertices: slots as i32,
                max_len_Pu_Pv: slots as i32,
                ..BN_DATA::default()
            }
        }

        fn no_path_network(heap: &mut SourceHeap) -> BN_STRUCT {
            let incident = heap.allocate_model_storage(Vec::<i32>::new()).unwrap();
            let vertices = heap
                .allocate_model_storage(vec![BNS_VERTEX {
                    type_: BNS_VERT_TYPE_ATOM as u16,
                    iedge: incident,
                    ..BNS_VERTEX::default()
                }])
                .unwrap();
            let edges = heap.allocate_model_storage(Vec::<BNS_EDGE>::new()).unwrap();
            let path = heap.allocate_model_storage(vec![BNS_ALT_PATH::default(); 5]).unwrap();
            let mut network = BN_STRUCT {
                num_atoms: 1,
                num_vertices: 1,
                vert: vertices,
                edge: edges,
                max_altp: 1,
                ..BN_STRUCT::default()
            };
            network.altp[0] = path;
            network
        }

        let mut heap = SourceHeap::default();
        let mut network = no_path_network(&mut heap);
        {
            let path = heap.slice_mut(network.altp[0]).unwrap();
            path[tagAltPathConst_iALTP_FLOW as usize].set_flow(0, 12);
            path[tagAltPathConst_iALTP_PATH_LEN as usize].set_number(3);
            path[tagAltPathConst_iALTP_START_ATOM as usize].set_number(4);
            path[tagAltPathConst_iALTP_END_ATOM as usize].set_number(5);
        }
        network.num_altp = 1;
        network.bChangeFlow = 27;
        let mut search = data(&mut heap, 6);
        let mut first = 10;
        let mut last = 11;
        let mut path_len = 12;
        let mut delta_h = 13;
        let mut delta_charge = 14;
        let mut visited = 15;
        assert_eq!(
            RunBnsTestOnce(
                &mut heap,
                &mut network,
                &mut search,
                &[VAL_AT::default()],
                &mut first,
                &mut last,
                &mut path_len,
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(0)
        );
        assert_eq!((first, last, path_len), (NO_VERTEX, NO_VERTEX, 0));
        assert_eq!((delta_h, delta_charge, visited), (13, 14, 15));
        assert_eq!(network.alt_path, SourceMutPointer::null());
        assert_eq!((network.num_altp, network.bChangeFlow), (0, 0));
        let cleared = heap.slice(network.altp[0].as_const()).unwrap();
        assert_eq!(
            (
                cleared[tagAltPathConst_iALTP_FLOW as usize].flow(0),
                cleared[tagAltPathConst_iALTP_PATH_LEN as usize].number(),
                cleared[tagAltPathConst_iALTP_START_ATOM as usize].number(),
                cleared[tagAltPathConst_iALTP_END_ATOM as usize].number(),
            ),
            (0, 0, NO_VERTEX, NO_VERTEX)
        );
        assert_eq!(search.QSize, -1);

        let mut invalid_data_heap = SourceHeap::default();
        let mut invalid_network = no_path_network(&mut invalid_data_heap);
        let mut invalid_search = data(&mut invalid_data_heap, 6);
        invalid_search.Pu = SourceMutPointer::null();
        invalid_search.Pv = SourceMutPointer::null();
        assert_eq!(
            RunBnsTestOnce(
                &mut invalid_data_heap,
                &mut invalid_network,
                &mut invalid_search,
                &[VAL_AT::default()],
                &mut first,
                &mut last,
                &mut path_len,
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(-96)
        );
        assert_eq!((first, last, path_len), (NO_VERTEX, NO_VERTEX, 0));
        assert_eq!(invalid_search.QSize, -1);

        let mut augment_heap = SourceHeap::default();
        let incident0 = augment_heap.allocate_model_storage(vec![0_i32]).unwrap();
        let incident1 = augment_heap.allocate_model_storage(vec![0_i32]).unwrap();
        let vertices = augment_heap
            .allocate_model_storage(vec![
                BNS_VERTEX {
                    st_edge: crate::source_types::BNS_ST_EDGE {
                        cap: 1,
                        ..crate::source_types::BNS_ST_EDGE::default()
                    },
                    type_: BNS_VERT_TYPE_ATOM as u16,
                    num_adj_edges: 1,
                    max_adj_edges: 1,
                    iedge: incident0,
                },
                BNS_VERTEX {
                    st_edge: crate::source_types::BNS_ST_EDGE {
                        cap: 1,
                        ..crate::source_types::BNS_ST_EDGE::default()
                    },
                    type_: BNS_VERT_TYPE_ATOM as u16,
                    num_adj_edges: 1,
                    max_adj_edges: 1,
                    iedge: incident1,
                },
            ])
            .unwrap();
        let edges = augment_heap
            .allocate_model_storage(vec![BNS_EDGE {
                neighbor1: 0,
                neighbor12: 1,
                neigh_ord: [0, 0],
                cap: 1,
                ..BNS_EDGE::default()
            }])
            .unwrap();
        let mut augment_path = vec![BNS_ALT_PATH::default(); 6];
        augment_path[tagAltPathConst_iALTP_MAX_LEN as usize].set_number(6);
        let alt_path = augment_heap.allocate_model_storage(augment_path).unwrap();
        let mut augment_network = BN_STRUCT {
            num_atoms: 2,
            num_vertices: 2,
            num_bonds: 1,
            num_edges: 1,
            vert: vertices,
            edge: edges,
            max_altp: 1,
            ..BN_STRUCT::default()
        };
        augment_network.altp[0] = alt_path;
        let mut augment_search = data(&mut augment_heap, 8);
        delta_h = -1;
        delta_charge = -1;
        visited = -1;
        assert_eq!(
            RunBnsTestOnce(
                &mut augment_heap,
                &mut augment_network,
                &mut augment_search,
                &[VAL_AT::default(), VAL_AT::default()],
                &mut first,
                &mut last,
                &mut path_len,
                &mut delta_h,
                &mut delta_charge,
                &mut visited,
            ),
            Ok(1)
        );
        assert_eq!((first, last, path_len), (0, 1, 1));
        assert_eq!((delta_h, delta_charge, visited), (0, 0, 1));
        assert_eq!(augment_network.num_altp, 0);
        assert_eq!(augment_network.alt_path, SourceMutPointer::null());
    }

    #[test]
    fn source_port__ichirvr1__runbnsrestoreonce__line_4546() {
        fn empty_valid_data(heap: &mut SourceHeap) -> BN_DATA {
            BN_DATA {
                BasePtr: heap.allocate_model_storage(Vec::<i32>::new()).unwrap(),
                SwitchEdge: heap
                    .allocate_model_storage(Vec::<crate::source_types::Edge>::new())
                    .unwrap(),
                Tree: heap.allocate_model_storage(Vec::<i8>::new()).unwrap(),
                ScanQ: heap.allocate_model_storage(Vec::<i32>::new()).unwrap(),
                QSize: -1,
                Pu: heap.allocate_model_storage(Vec::<i32>::new()).unwrap(),
                Pv: heap.allocate_model_storage(Vec::<i32>::new()).unwrap(),
                ..BN_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let clock = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_CLOCK::default()])
            .unwrap();
        let mut bns = BN_STRUCT {
            max_altp: 0,
            num_altp: 7,
            tot_st_flow: 11,
            ic: clock,
            ..BN_STRUCT::default()
        };
        let mut data = empty_valid_data(&mut heap);
        let mut valence_atoms = Vec::new();
        let mut groups = ALL_TC_GROUPS::default();
        assert_eq!(
            RunBnsRestoreOnce(&mut heap, &mut bns, &mut data, &mut valence_atoms, &mut groups, 0,),
            Ok(0)
        );
        assert_eq!(bns.num_altp, 0);
        assert!(bns.alt_path.is_null());
        assert_eq!(bns.tot_st_flow, 11);
        assert_eq!(data.QSize, -1);

        let mut missing_heap = SourceHeap::default();
        let missing_clock = missing_heap
            .allocate_model_storage(vec![crate::source_types::INCHI_CLOCK::default()])
            .unwrap();
        let mut missing_bns = BN_STRUCT {
            max_altp: 0,
            tot_st_flow: 13,
            ic: missing_clock,
            ..BN_STRUCT::default()
        };
        let mut missing_data = BN_DATA::default();
        assert_eq!(
            RunBnsRestoreOnce(
                &mut missing_heap,
                &mut missing_bns,
                &mut missing_data,
                &mut [],
                &mut ALL_TC_GROUPS::default(),
                0,
            ),
            Ok(-126)
        );
        assert_eq!(missing_bns.tot_st_flow, 13);
        assert_eq!(missing_data.QSize, -1);

        let mut timeout_heap = SourceHeap::default();
        let timeout_clock = timeout_heap
            .allocate_model_storage(vec![crate::source_types::INCHI_CLOCK::default()])
            .unwrap();
        let timeout_start = timeout_heap
            .allocate_model_storage(vec![crate::source_types::inchiTime { clockTime: -1 }])
            .unwrap();
        let mut timeout_bns = BN_STRUCT {
            max_altp: 0,
            num_altp: 9,
            tot_st_flow: 17,
            ic: timeout_clock,
            ulTimeOutTime: timeout_start,
            ..BN_STRUCT::default()
        };
        let mut timeout_data = BN_DATA::default();
        assert_eq!(
            RunBnsRestoreOnce(
                &mut timeout_heap,
                &mut timeout_bns,
                &mut timeout_data,
                &mut [],
                &mut ALL_TC_GROUPS::default(),
                0,
            ),
            Ok(crate::source_types::BNS_TIMEOUT)
        );
        assert_eq!(timeout_bns.num_altp, 0);
        assert!(timeout_bns.alt_path.is_null());
        assert_eq!(timeout_bns.tot_st_flow, 17);
        assert_eq!(timeout_data.QSize, 0);
    }

    #[test]
    fn source_port__ichirvr1__get_pva_atom_type__line_4628() {
        let classify = |electrons: i8, row: i8, periodic_number: u8, bond_type: i32| -> Result<i32, SourceHeapError> {
            get_pVA_atom_type(
                &[VAL_AT {
                    cNumValenceElectrons: electrons,
                    cPeriodicRowNumber: row,
                    cPeriodicNumber: periodic_number,
                    ..VAL_AT::default()
                }],
                &[inp_ATOM {
                    el_number: 6,
                    ..inp_ATOM::default()
                }],
                0,
                bond_type,
            )
        };
        for (electrons, row, periodic_number, expected) in [
            (4, 1, 5, EL_TYPE_C as i32),
            (4, 2, 5, 0),
            (6, 1, 7, EL_TYPE_O as i32),
            (6, 2, 15, EL_TYPE_S as i32),
            (6, 4, 51, EL_TYPE_S as i32),
            (6, 5, 83, 0),
            (6, -1, 15, EL_TYPE_S as i32),
            (5, 1, 6, EL_TYPE_N as i32),
            (5, 2, 14, EL_TYPE_P as i32),
            (5, -1, 14, EL_TYPE_P as i32),
            (7, 1, 5, EL_TYPE_X as i32),
            (7, 1, 11, 0),
        ] {
            assert_eq!(
                classify(electrons, row, periodic_number, BOND_TYPE_DOUBLE as i32),
                Ok(expected),
                "electrons={electrons}, row={row}, periodic={periodic_number}"
            );
        }

        let terminal_atoms = [
            inp_ATOM {
                el_number: 8,
                neighbor: {
                    let mut neighbors = [0_u16; 20];
                    neighbors[0] = 1;
                    neighbors
                },
                bond_type: {
                    let mut bonds = [0_u8; 20];
                    bonds[0] = BOND_TYPE_SINGLE as u8;
                    bonds
                },
                valence: 1,
                chem_bonds_valence: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                valence: 1,
                chem_bonds_valence: 1,
                neighbor: {
                    let mut neighbors = [0_u16; 20];
                    neighbors[0] = 0;
                    neighbors
                },
                bond_type: {
                    let mut bonds = [0_u8; 20];
                    bonds[0] = BOND_TYPE_SINGLE as u8;
                    bonds
                },
                ..inp_ATOM::default()
            },
        ];
        let oxygen_va = [
            VAL_AT {
                cNumValenceElectrons: 6,
                cPeriodicRowNumber: 1,
                cPeriodicNumber: 7,
                ..VAL_AT::default()
            },
            VAL_AT::default(),
        ];
        assert_eq!(
            get_pVA_atom_type(&oxygen_va, &terminal_atoms, 0, BOND_TYPE_SINGLE as i32),
            Ok((EL_TYPE_O | EL_TYPE_OSt) as i32)
        );
        assert_eq!(
            get_pVA_atom_type(&oxygen_va, &terminal_atoms, 0, BOND_TYPE_DOUBLE as i32),
            Ok(EL_TYPE_O as i32)
        );
        let sulfur_va = [
            VAL_AT {
                cNumValenceElectrons: 6,
                cPeriodicRowNumber: 2,
                cPeriodicNumber: 15,
                ..VAL_AT::default()
            },
            VAL_AT::default(),
        ];
        assert_eq!(
            get_pVA_atom_type(&sulfur_va, &terminal_atoms, 0, BOND_TYPE_SINGLE as i32),
            Ok((EL_TYPE_S | EL_TYPE_OSt) as i32)
        );

        let endpoint_case = |atom: inp_ATOM| {
            get_pVA_atom_type(
                &[VAL_AT {
                    cNumValenceElectrons: 6,
                    cPeriodicRowNumber: 1,
                    cPeriodicNumber: 7,
                    ..VAL_AT::default()
                }],
                &[atom],
                0,
                BOND_TYPE_DOUBLE as i32,
            )
        };
        let endpoint = inp_ATOM {
            el_number: 8,
            valence: 1,
            chem_bonds_valence: 1,
            num_H: 1,
            ..inp_ATOM::default()
        };
        assert_eq!(endpoint_case(endpoint.clone()), Ok((EL_TYPE_O | EL_TYPE_PT) as i32));
        assert_eq!(
            endpoint_case(inp_ATOM {
                charge: -1,
                num_H: 0,
                ..endpoint.clone()
            }),
            Ok((EL_TYPE_O | EL_TYPE_PT) as i32)
        );
        for rejected in [
            inp_ATOM {
                radical: 1,
                ..endpoint.clone()
            },
            inp_ATOM {
                charge: 1,
                ..endpoint.clone()
            },
            inp_ATOM {
                charge: -2,
                ..endpoint.clone()
            },
            inp_ATOM {
                valence: 2,
                ..endpoint.clone()
            },
            inp_ATOM {
                chem_bonds_valence: 0,
                ..endpoint.clone()
            },
            inp_ATOM {
                el_number: 6,
                ..endpoint
            },
        ] {
            assert_eq!(endpoint_case(rejected), Ok(EL_TYPE_O as i32));
        }

        assert_eq!(
            get_pVA_atom_type(&[], &[], -1, BOND_TYPE_SINGLE as i32),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            get_pVA_atom_type(&[VAL_AT::default()], &[], 0, BOND_TYPE_SINGLE as i32),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr1__allocedgelist__line_4690() {
        let mut heap = SourceHeap::default();
        let mut list = EDGE_LIST {
            num_alloc: 9,
            num_edges: 7,
            pnEdges: SourceMutPointer::null(),
        };
        assert_eq!(AllocEdgeList(&mut heap, &mut list, EDGE_LIST_CLEAR), Ok(0));
        assert_eq!(list, EDGE_LIST::default());

        assert_eq!(AllocEdgeList(&mut heap, &mut list, 3), Ok(0));
        assert_eq!((list.num_alloc, list.num_edges), (3, 0));
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap()[..3], [0, 0, 0]);
        assert_eq!(heap.live_allocations_of::<i32>(), 1);
        heap.slice_mut(list.pnEdges).unwrap()[..3].copy_from_slice(&[7, 8, 9]);
        list.num_edges = 3;

        let first_pointer = list.pnEdges;
        assert_eq!(AllocEdgeList(&mut heap, &mut list, 3), Ok(0));
        assert_eq!(list.pnEdges, first_pointer);
        assert_eq!(list.num_edges, 3);

        assert_eq!(AllocEdgeList(&mut heap, &mut list, 2), Ok(0));
        assert_eq!((list.num_alloc, list.num_edges), (2, 2));
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap()[..2], [7, 8]);
        assert_eq!(heap.live_allocations_of::<i32>(), 1);

        assert_eq!(AllocEdgeList(&mut heap, &mut list, 4), Ok(0));
        assert_eq!((list.num_alloc, list.num_edges), (4, 2));
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap()[..4], [7, 8, 0, 0]);
        assert_eq!(heap.live_allocations_of::<i32>(), 1);

        let unchanged = list.clone();
        assert_eq!(AllocEdgeList(&mut heap, &mut list, -3), Ok(0));
        assert_eq!(list, unchanged);

        assert_eq!(AllocEdgeList(&mut heap, &mut list, EDGE_LIST_CLEAR), Ok(0));
        assert_eq!(list, EDGE_LIST::default());
        assert_eq!(heap.live_allocations_of::<i32>(), 1);

        let free_pointer = heap.allocate_model_storage(vec![4_i32, 5]).unwrap();
        let mut free_list = EDGE_LIST {
            num_alloc: 2,
            num_edges: 2,
            pnEdges: free_pointer,
        };
        assert_eq!(heap.live_allocations_of::<i32>(), 2);
        assert_eq!(AllocEdgeList(&mut heap, &mut free_list, EDGE_LIST_FREE), Ok(0));
        assert_eq!(free_list, EDGE_LIST::default());
        assert_eq!(heap.live_allocations_of::<i32>(), 1);

        let mut failure_heap = SourceHeap::default();
        let leaked_pointer = failure_heap.allocate_model_storage(vec![11_i32, 12]).unwrap();
        let mut failure_list = EDGE_LIST {
            num_alloc: 2,
            num_edges: 2,
            pnEdges: leaked_pointer,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(AllocEdgeList(&mut failure_heap, &mut failure_list, 5), Ok(RI_ERR_ALLOC));
        assert!(failure_list.pnEdges.is_null());
        assert_eq!((failure_list.num_alloc, failure_list.num_edges), (2, 2));
        assert_eq!(failure_heap.live_allocations_of::<i32>(), 1);
        assert_eq!(failure_heap.slice(leaked_pointer.as_const()).unwrap(), [11, 12]);
    }

    #[test]
    fn source_port__ichirvr1__addtoedgelist__line_4738() {
        let mut invalid_growth_heap = SourceHeap::default();
        for n_add_len in [0, -1, i32::MIN] {
            let mut list = EDGE_LIST::default();
            assert_eq!(
                AddToEdgeList(&mut invalid_growth_heap, &mut list, 17, n_add_len),
                Ok(RI_ERR_PROGR)
            );
            assert_eq!(list, EDGE_LIST::default());
            assert_eq!(invalid_growth_heap.live_allocations_of::<i32>(), 0);
        }

        let mut heap = SourceHeap::default();
        let mut list = EDGE_LIST::default();
        assert_eq!(AddToEdgeList(&mut heap, &mut list, i32::MIN, 2), Ok(0));
        assert_eq!((list.num_alloc, list.num_edges), (2, 1));
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap()[..2], [i32::MIN, 0]);

        assert_eq!(AddToEdgeList(&mut heap, &mut list, -7, 0), Ok(0));
        assert_eq!((list.num_alloc, list.num_edges), (2, 2));
        assert_eq!(heap.slice(list.pnEdges.as_const()).unwrap()[..2], [i32::MIN, -7]);

        let first_pointer = list.pnEdges;
        assert_eq!(AddToEdgeList(&mut heap, &mut list, i32::MAX, 3), Ok(0));
        assert_ne!(list.pnEdges, first_pointer);
        assert_eq!((list.num_alloc, list.num_edges), (5, 3));
        assert_eq!(
            heap.slice(list.pnEdges.as_const()).unwrap()[..5],
            [i32::MIN, -7, i32::MAX, 0, 0]
        );
        assert_eq!(heap.live_allocations_of::<i32>(), 1);

        assert_eq!(AddToEdgeList(&mut heap, &mut list, 31, -9), Ok(0));
        assert_eq!((list.num_alloc, list.num_edges), (5, 4));
        assert_eq!(
            heap.slice(list.pnEdges.as_const()).unwrap()[..5],
            [i32::MIN, -7, i32::MAX, 31, 0]
        );

        let mut failure_heap = SourceHeap::default();
        let leaked_pointer = failure_heap.allocate_model_storage(vec![101_i32, 202]).unwrap();
        let mut failure_list = EDGE_LIST {
            num_alloc: 2,
            num_edges: 2,
            pnEdges: leaked_pointer,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            AddToEdgeList(&mut failure_heap, &mut failure_list, 303, 4),
            Ok(RI_ERR_ALLOC)
        );
        assert!(failure_list.pnEdges.is_null());
        assert_eq!((failure_list.num_alloc, failure_list.num_edges), (2, 2));
        assert_eq!(failure_heap.live_allocations_of::<i32>(), 1);
        assert_eq!(failure_heap.slice(leaked_pointer.as_const()).unwrap(), [101, 202]);
    }

    #[test]
    fn source_port__ichirvr1__removefromedgelistbyindex__line_4759() {
        let mut empty_heap = SourceHeap::default();
        let mut empty = EDGE_LIST::default();
        assert_eq!(RemoveFromEdgeListByIndex(&mut empty_heap, &mut empty, 0), Ok(-1));
        assert_eq!(empty, EDGE_LIST::default());

        for (index, expected) in [
            (0, [20, 30, 40, 0, 91, 92]),
            (1, [10, 30, 40, 0, 91, 92]),
            (2, [10, 20, 40, 0, 91, 92]),
            (3, [10, 20, 30, 0, 91, 92]),
        ] {
            let mut heap = SourceHeap::default();
            let pointer = heap.allocate_model_storage(vec![10_i32, 20, 30, 40, 91, 92]).unwrap();
            let mut list = EDGE_LIST {
                num_alloc: 6,
                num_edges: 4,
                pnEdges: pointer,
            };
            assert_eq!(RemoveFromEdgeListByIndex(&mut heap, &mut list, index), Ok(0));
            assert_eq!((list.num_alloc, list.num_edges), (6, 3));
            assert_eq!(heap.slice(pointer.as_const()).unwrap(), expected);
        }

        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![i32::MIN, -1, 0, i32::MAX]).unwrap();
        let mut list = EDGE_LIST {
            num_alloc: 4,
            num_edges: 4,
            pnEdges: pointer,
        };
        let unchanged = heap.slice(pointer.as_const()).unwrap().to_vec();
        for index in [4, 5, i32::MAX] {
            assert_eq!(RemoveFromEdgeListByIndex(&mut heap, &mut list, index), Ok(-1));
            assert_eq!((list.num_alloc, list.num_edges), (4, 4));
            assert_eq!(heap.slice(pointer.as_const()).unwrap(), unchanged);
        }
    }

    #[test]
    fn source_port__ichirvr1__removefromedgelistbyvalue__line_4794() {
        let mut empty_heap = SourceHeap::default();
        let mut empty = EDGE_LIST::default();
        assert_eq!(RemoveFromEdgeListByValue(&mut empty_heap, &mut empty, i32::MIN), Ok(0));
        assert_eq!(empty, EDGE_LIST::default());

        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![7_i32, 2, 7, 7, 3, 7]).unwrap();
        let mut list = EDGE_LIST {
            num_alloc: 6,
            num_edges: 6,
            pnEdges: pointer,
        };
        assert_eq!(RemoveFromEdgeListByValue(&mut heap, &mut list, 7), Ok(4));
        assert_eq!((list.num_alloc, list.num_edges), (6, 2));
        assert_eq!(heap.slice(pointer.as_const()).unwrap(), &[2, 3, 0, 0, 0, 0]);

        let unchanged = heap.slice(pointer.as_const()).unwrap().to_vec();
        assert_eq!(RemoveFromEdgeListByValue(&mut heap, &mut list, i32::MAX), Ok(0));
        assert_eq!((list.num_alloc, list.num_edges), (6, 2));
        assert_eq!(heap.slice(pointer.as_const()).unwrap(), unchanged);

        let extremes = heap
            .allocate_model_storage(vec![i32::MIN, 0, i32::MAX, i32::MIN])
            .unwrap();
        let mut extreme_list = EDGE_LIST {
            num_alloc: 4,
            num_edges: 4,
            pnEdges: extremes,
        };
        assert_eq!(RemoveFromEdgeListByValue(&mut heap, &mut extreme_list, i32::MIN), Ok(2));
        assert_eq!(extreme_list.num_edges, 2);
        assert_eq!(heap.slice(extremes.as_const()).unwrap(), &[0, i32::MAX, 0, 0]);
    }

    #[test]
    fn source_port__ichirvr1__findinedgelist__line_4777() {
        let heap = SourceHeap::default();
        let empty = EDGE_LIST::default();
        assert_eq!(FindInEdgeList(&heap, &empty, i32::MIN), Ok(-1));

        let mut heap = SourceHeap::default();
        let pointer = heap
            .allocate_model_storage(vec![i32::MIN, 7_i32, 0, 7, i32::MAX, 7, 91, 92])
            .unwrap();
        let list = EDGE_LIST {
            num_alloc: 8,
            num_edges: 6,
            pnEdges: pointer,
        };

        assert_eq!(FindInEdgeList(&heap, &list, 7), Ok(5));
        assert_eq!(FindInEdgeList(&heap, &list, i32::MIN), Ok(0));
        assert_eq!(FindInEdgeList(&heap, &list, i32::MAX), Ok(4));
        assert_eq!(FindInEdgeList(&heap, &list, -1), Ok(-1));
        assert_eq!(FindInEdgeList(&heap, &list, 91), Ok(-1));

        let one_visible = EDGE_LIST { num_edges: 1, ..list };
        assert_eq!(FindInEdgeList(&heap, &one_visible, i32::MIN), Ok(0));
        assert_eq!(FindInEdgeList(&heap, &one_visible, 7), Ok(-1));
    }

    #[test]
    fn source_port__ichirvr1__removeforbiddenedgemask__line_4869() {
        let mut empty_heap = SourceHeap::default();
        let empty_bns = BN_STRUCT::default();
        let empty_list = EDGE_LIST {
            num_edges: -1,
            ..EDGE_LIST::default()
        };
        assert_eq!(
            RemoveForbiddenEdgeMask(&mut empty_heap, &empty_bns, &empty_list, i32::MIN),
            Ok(())
        );

        let mut heap = SourceHeap::default();
        let edge_pointer = heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    forbidden: -1,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    forbidden: 0x7f,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    forbidden: 0x55,
                    ..BNS_EDGE::default()
                },
            ])
            .unwrap();
        let list_pointer = heap.allocate_model_storage(vec![0_i32, 1, 0]).unwrap();
        let bns = BN_STRUCT {
            edge: edge_pointer,
            num_edges: 3,
            ..BN_STRUCT::default()
        };
        let list = EDGE_LIST {
            num_alloc: 3,
            num_edges: 3,
            pnEdges: list_pointer,
        };
        assert_eq!(RemoveForbiddenEdgeMask(&mut heap, &bns, &list, 0x05), Ok(()));
        let edges = heap.slice(edge_pointer.as_const()).unwrap();
        assert_eq!(edges[0].forbidden, -6);
        assert_eq!(edges[1].forbidden, 0x7a);
        assert_eq!(edges[2].forbidden, 0x55);

        assert_eq!(RemoveForbiddenEdgeMask(&mut heap, &bns, &list, 0), Ok(()));
        assert_eq!(heap.slice(edge_pointer.as_const()).unwrap()[0].forbidden, -6);

        assert_eq!(RemoveForbiddenEdgeMask(&mut heap, &bns, &list, -1), Ok(()));
        let edges = heap.slice(edge_pointer.as_const()).unwrap();
        assert_eq!(
            edges.iter().map(|edge| edge.forbidden).collect::<Vec<_>>(),
            [0, 0, 0x55]
        );
    }

    #[test]
    fn source_port__ichirvr1__setforbiddenedgemask__line_4880() {
        let mut empty_heap = SourceHeap::default();
        assert_eq!(
            SetForbiddenEdgeMask(
                &mut empty_heap,
                &BN_STRUCT::default(),
                &EDGE_LIST {
                    num_edges: -1,
                    ..EDGE_LIST::default()
                },
                i32::MIN,
            ),
            Ok(())
        );

        let mut heap = SourceHeap::default();
        let edges = heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    forbidden: 0x40,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    forbidden: -128,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE::default(),
            ])
            .unwrap();
        let listed = heap.allocate_model_storage(vec![0_i32, 1, 0]).unwrap();
        let network = BN_STRUCT {
            edge: edges,
            ..BN_STRUCT::default()
        };
        let list = EDGE_LIST {
            num_alloc: 3,
            num_edges: 3,
            pnEdges: listed,
        };
        assert_eq!(SetForbiddenEdgeMask(&mut heap, &network, &list, 0x05), Ok(()));
        assert_eq!(
            heap.slice(edges.as_const())
                .unwrap()
                .iter()
                .map(|edge| edge.forbidden)
                .collect::<Vec<_>>(),
            [0x45, -123, 0]
        );
        assert_eq!(SetForbiddenEdgeMask(&mut heap, &network, &list, 0), Ok(()));
        assert_eq!(SetForbiddenEdgeMask(&mut heap, &network, &list, -1), Ok(()));
        assert_eq!(
            heap.slice(edges.as_const())
                .unwrap()
                .iter()
                .map(|edge| edge.forbidden)
                .collect::<Vec<_>>(),
            [-1, -1, 0]
        );
    }

    #[test]
    fn source_port__ichirvr1__removeforbiddenbondflowbits__line_4891() {
        let mut empty_heap = SourceHeap::default();
        for bond_count in [i32::MIN, -1, 0] {
            assert_eq!(
                RemoveForbiddenBondFlowBits(
                    &mut empty_heap,
                    &BN_STRUCT {
                        num_bonds: bond_count,
                        ..BN_STRUCT::default()
                    },
                    i32::MIN,
                ),
                Ok(())
            );
        }

        for (mask, expected) in [
            (0_i32, vec![-1_i8, 0x7f, 0x55, -128]),
            (0x05_i32, vec![-6_i8, 0x7a, 0x50, -128]),
            (-1_i32, vec![0_i8, 0, 0, -128]),
            (i32::MIN, vec![-1_i8, 0x7f, 0x55, -128]),
        ] {
            let mut heap = SourceHeap::default();
            let edges = heap
                .allocate_model_storage(vec![
                    BNS_EDGE {
                        forbidden: -1,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE {
                        forbidden: 0x7f,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE {
                        forbidden: 0x55,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE {
                        forbidden: -128,
                        ..BNS_EDGE::default()
                    },
                ])
                .unwrap();
            let network = BN_STRUCT {
                edge: edges,
                num_bonds: 3,
                num_edges: 4,
                ..BN_STRUCT::default()
            };
            assert_eq!(RemoveForbiddenBondFlowBits(&mut heap, &network, mask), Ok(()));
            assert_eq!(
                heap.slice(edges.as_const())
                    .unwrap()
                    .iter()
                    .map(|edge| edge.forbidden)
                    .collect::<Vec<_>>(),
                expected,
                "mask={mask}"
            );
        }

        let mut all_heap = SourceHeap::default();
        let all_edges = all_heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    forbidden: 0x3f,
                    ..BNS_EDGE::default()
                };
                4
            ])
            .unwrap();
        let all_network = BN_STRUCT {
            edge: all_edges,
            num_bonds: 4,
            num_edges: 4,
            ..BN_STRUCT::default()
        };
        assert_eq!(RemoveForbiddenBondFlowBits(&mut all_heap, &all_network, 0x30), Ok(()));
        assert_eq!(
            all_heap
                .slice(all_edges.as_const())
                .unwrap()
                .iter()
                .map(|edge| edge.forbidden)
                .collect::<Vec<_>>(),
            [0x0f, 0x0f, 0x0f, 0x0f]
        );
    }

    #[test]
    fn source_port__ichirvr1__getchargeflowerupperedge__line_4915() {
        assert_eq!(
            GetChargeFlowerUpperEdge(&SourceHeap::default(), &BN_STRUCT::default(), &[], -1,),
            Ok(NO_VERTEX)
        );

        let mut heap = SourceHeap::default();
        let adjacency = [
            Vec::new(),
            vec![0_i32],
            vec![0_i32, 1, 2],
            vec![1_i32, 3],
            vec![2_i32, 3, 4],
            vec![4_i32],
        ];
        let incident = adjacency
            .iter()
            .map(|edges| {
                if edges.is_empty() {
                    SourceMutPointer::null()
                } else {
                    heap.allocate_model_storage(edges.clone()).unwrap()
                }
            })
            .collect::<Vec<_>>();
        let vertices = heap
            .allocate_model_storage(
                adjacency
                    .iter()
                    .enumerate()
                    .map(|(index, edges)| BNS_VERTEX {
                        type_: if index == 0 || index == 5 {
                            BNS_VERT_TYPE_ATOM as u16
                        } else if index == 1 {
                            BNS_VERT_TYPE_C_GROUP as u16
                        } else {
                            BNS_VERT_TYPE__AUX as u16
                        },
                        num_adj_edges: edges.len() as u16,
                        max_adj_edges: edges.len() as u16,
                        iedge: incident[index],
                        ..BNS_VERTEX::default()
                    })
                    .collect::<Vec<_>>(),
            )
            .unwrap();
        let edges = heap
            .allocate_model_storage(vec![
                BNS_EDGE {
                    neighbor1: 1,
                    neighbor12: 3,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    neighbor1: 2,
                    neighbor12: 1,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    neighbor1: 2,
                    neighbor12: 6,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    neighbor1: 3,
                    neighbor12: 7,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    neighbor1: 4,
                    neighbor12: 1,
                    ..BNS_EDGE::default()
                },
            ])
            .unwrap();
        let network = BN_STRUCT {
            num_vertices: 6,
            num_edges: 5,
            vert: vertices,
            edge: edges,
            ..BN_STRUCT::default()
        };

        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(1));

        heap.slice_mut(edges).unwrap()[0].neighbor1 = 2;
        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(1));
        heap.slice_mut(edges).unwrap()[0].neighbor1 = 1;

        heap.slice_mut(incident[2]).unwrap()[1..3].copy_from_slice(&[2, 1]);
        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(1));
        heap.slice_mut(incident[2]).unwrap()[1..3].copy_from_slice(&[1, 2]);

        heap.slice_mut(vertices).unwrap()[2].type_ = BNS_VERT_TYPE_ATOM as u16;
        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(NO_VERTEX));
        heap.slice_mut(vertices).unwrap()[2].type_ = BNS_VERT_TYPE__AUX as u16;

        heap.slice_mut(vertices).unwrap()[3].type_ = BNS_VERT_TYPE_ATOM as u16;
        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(NO_VERTEX));
        heap.slice_mut(vertices).unwrap()[3].type_ = BNS_VERT_TYPE__AUX as u16;

        heap.slice_mut(vertices).unwrap()[4].num_adj_edges = 2;
        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(NO_VERTEX));
        heap.slice_mut(vertices).unwrap()[4].num_adj_edges = 3;

        heap.slice_mut(vertices).unwrap()[5].type_ = BNS_VERT_TYPE__AUX as u16;
        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(NO_VERTEX));
        heap.slice_mut(vertices).unwrap()[5].type_ = BNS_VERT_TYPE_ATOM as u16;

        heap.slice_mut(vertices).unwrap()[1].type_ = BNS_VERT_TYPE__AUX as u16;
        assert_eq!(GetChargeFlowerUpperEdge(&heap, &network, &[], 0), Ok(NO_VERTEX));
    }

    #[test]
    fn source_port__ichirvr1__cmp_charge_val__line_732() {
        let value = |valence, charge, order| crate::source_types::CHARGE_VAL {
            nValence: valence,
            nCharge: charge,
            nValenceOrderingNumber: order,
        };
        let cases = [
            (value(1, 0, 0), value(2, 0, 0), -1),
            (value(2, 0, 0), value(1, 0, 0), 1),
            (value(2, 0, 0), value(2, 1, 0), -1),
            (value(2, 2, 0), value(2, 1, 0), 1),
            (value(2, 1, 0), value(2, -1, 0), -2),
            (value(2, -1, 0), value(2, 1, 0), 2),
            (value(2, 1, 3), value(2, 1, 5), -2),
            (value(2, 1, 5), value(2, 1, 3), 2),
            (value(2, 1, 3), value(2, 1, 3), 0),
            (value(i32::MIN, 0, 0), value(0, 0, 0), i32::MIN),
            (value(i32::MAX, 0, 0), value(0, 0, 0), i32::MAX),
            (value(0, i32::MAX, 0), value(0, 0, 0), i32::MAX),
            (value(0, 0, i32::MIN), value(0, 0, 0), i32::MIN),
            (value(0, 0, i32::MAX), value(0, 0, 0), i32::MAX),
        ];
        for (first, second, expected) in cases {
            let first_before = first.clone();
            let second_before = second.clone();
            assert_eq!(cmp_charge_val(&first, &second), expected);
            assert_eq!(first, first_before);
            assert_eq!(second, second_before);
        }
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__cmp_charge_val__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--cmp-charge-val-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "cmp_charge_val");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_value = |side: &str, field: &str| {
                i32::try_from(
                    official["input"][side][field]
                        .as_i64()
                        .unwrap_or_else(|| panic!("{case_id}: {side}.{field} must be i32")),
                )
                .unwrap_or_else(|_| panic!("{case_id}: {side}.{field} exceeds i32"))
            };
            let first = crate::source_types::CHARGE_VAL {
                nValence: parse_value("first", "valence"),
                nCharge: parse_value("first", "charge"),
                nValenceOrderingNumber: parse_value("first", "order"),
            };
            let second = crate::source_types::CHARGE_VAL {
                nValence: parse_value("second", "valence"),
                nCharge: parse_value("second", "charge"),
                nValenceOrderingNumber: parse_value("second", "order"),
            };
            let first_before = first.clone();
            let second_before = second.clone();
            let rust = cmp_charge_val(&first, &second);
            let expected = i32::try_from(official["output"]["result"].as_i64().expect("result must be i32"))
                .expect("official result exceeds i32");
            assert_eq!(rust, expected, "{case_id}");
            assert_eq!(first, first_before, "{case_id}");
            assert_eq!(second, second_before, "{case_id}");
            assert_eq!(official["output"]["first_unchanged"], true, "{case_id}");
            assert_eq!(official["output"]["second_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 14);
    }

    #[test]
    fn source_port__ichirvr1__comp_cc_cand__line_4584() {
        let candidate = |iat,
                         num_bonds,
                         chem_valence,
                         metal,
                         bonds_to_metal,
                         valence_electrons,
                         periodic_row,
                         charge_states,
                         element| CC_CAND {
            iat,
            num_bonds,
            chem_valence,
            cMetal: metal,
            cNumBondsToMetal: bonds_to_metal,
            cNumValenceElectrons: valence_electrons,
            cPeriodicRowNumber: periodic_row,
            cNumChargeStates: charge_states,
            el_number: element,
        };
        let base = candidate(7, 2, 3, 0, 0, 4, 2, 1, 6);
        let cases = [
            (
                candidate(7, 2, 3, i8::MIN, 0, 4, 2, 1, 6),
                candidate(7, 2, 3, i8::MAX, 0, 4, 2, 1, 6),
                255,
            ),
            (
                candidate(7, 2, 3, i8::MAX, 0, 4, 2, 1, 6),
                candidate(7, 2, 3, i8::MIN, 0, 4, 2, 1, 6),
                -255,
            ),
            (base.clone(), candidate(7, 2, 3, 0, 5, 4, 2, 1, 6), 5),
            (base.clone(), candidate(7, 2, 3, 0, 0, 4, 5, 1, 6), 3),
            (base.clone(), candidate(7, 5, 3, 0, 0, 4, 2, 1, 6), 3),
            (base.clone(), candidate(7, 2, 5, 0, 0, 4, 2, 1, 6), -2),
            (candidate(7, 2, 3, 0, 0, 0, 2, 1, 6), base.clone(), -1),
            (base.clone(), candidate(7, 2, 3, 0, 0, 0, 2, 1, 6), -1),
            (base.clone(), candidate(7, 2, 3, 0, 0, -4, 2, 1, 6), 0),
            (
                candidate(i32::MIN, 2, 3, 0, 0, 4, 2, 1, 6),
                candidate(i32::MIN + 1, 2, 3, 0, 0, 4, 2, 1, 6),
                1,
            ),
            (
                candidate(i32::MAX, 2, 3, 0, 0, 4, 2, 1, 6),
                candidate(i32::MAX - 1, 2, 3, 0, 0, 4, 2, 1, 6),
                -1,
            ),
            (base.clone(), candidate(7, 2, 3, 0, 0, 4, 2, i8::MAX, u8::MAX), 0),
        ];
        for (first, second, expected) in cases {
            let first_before = first.clone();
            let second_before = second.clone();
            assert_eq!(comp_cc_cand(&first, &second), expected);
            assert_eq!(first, first_before);
            assert_eq!(second, second_before);
        }
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__comp_cc_cand__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--comp-cc-cand-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "comp_cc_cand");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_i32 = |side: &str, field: &str| {
                i32::try_from(
                    official["input"][side][field]
                        .as_i64()
                        .unwrap_or_else(|| panic!("{case_id}: {side}.{field} must be i32")),
                )
                .unwrap_or_else(|_| panic!("{case_id}: {side}.{field} exceeds i32"))
            };
            let parse_i8 = |side: &str, field: &str| {
                i8::try_from(parse_i32(side, field)).unwrap_or_else(|_| panic!("{case_id}: {side}.{field} exceeds i8"))
            };
            let parse_candidate = |side: &str| CC_CAND {
                iat: parse_i32(side, "iat"),
                num_bonds: parse_i8(side, "num_bonds"),
                chem_valence: parse_i8(side, "chem_valence"),
                cMetal: parse_i8(side, "metal"),
                cNumBondsToMetal: parse_i8(side, "bonds_to_metal"),
                cNumValenceElectrons: parse_i8(side, "valence_electrons"),
                cPeriodicRowNumber: parse_i8(side, "periodic_row"),
                cNumChargeStates: parse_i8(side, "charge_states"),
                el_number: u8::try_from(parse_i32(side, "element"))
                    .unwrap_or_else(|_| panic!("{case_id}: {side}.element exceeds u8")),
            };
            let first = parse_candidate("first");
            let second = parse_candidate("second");
            let first_before = first.clone();
            let second_before = second.clone();
            let rust = comp_cc_cand(&first, &second);
            let expected = i32::try_from(official["output"]["result"].as_i64().expect("result must be i32"))
                .expect("official result exceeds i32");
            assert_eq!(rust, expected, "{case_id}");
            assert_eq!(first, first_before, "{case_id}");
            assert_eq!(second, second_before, "{case_id}");
            assert_eq!(official["output"]["first_unchanged"], true, "{case_id}");
            assert_eq!(official["output"]["second_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 12);
    }

    #[test]
    fn source_port__ichirvr1__clean_charge_val__line_792() {
        let value = |valence, charge, order| crate::source_types::CHARGE_VAL {
            nValence: valence,
            nCharge: charge,
            nValenceOrderingNumber: order,
        };
        let run = |values: &mut [crate::source_types::CHARGE_VAL],
                   atoms: &[inp_ATOM],
                   valence_atoms: &[crate::source_types::VAL_AT],
                   is_metal: i32,
                   mobile_h: i32,
                   endpoint: Option<&[AT_NUMB]>| {
            clean_charge_val(
                &mut CANON_GLOBALS::default(),
                values,
                values.len() as i32,
                atoms,
                valence_atoms,
                0,
                is_metal,
                mobile_h,
                endpoint,
            )
        };

        assert_eq!(
            clean_charge_val(
                &mut CANON_GLOBALS::default(),
                &mut [],
                0,
                &[],
                &[],
                i32::MIN,
                0,
                0,
                None,
            ),
            Ok(0)
        );

        let carbon = inp_ATOM {
            el_number: 6,
            valence: 2,
            ..inp_ATOM::default()
        };
        let default_valence = crate::source_types::VAL_AT::default();
        let mut sorted = [value(2, 0, 9), value(1, -1, 5), value(1, 1, 7), value(1, 0, 8)];
        assert_eq!(
            run(
                &mut sorted,
                std::slice::from_ref(&carbon),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(4)
        );
        assert_eq!(
            sorted,
            [value(1, 0, 8), value(1, 1, 7), value(1, -1, 5), value(2, 0, 9),]
        );

        let scandium = inp_ATOM {
            el_number: 21,
            valence: 1,
            ..inp_ATOM::default()
        };
        let mut metal = [value(2, 0, 1), value(1, 0, 2)];
        assert_eq!(
            run(
                &mut metal,
                std::slice::from_ref(&scandium),
                std::slice::from_ref(&default_valence),
                1,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!(metal[0], value(1, 0, 2));

        let low_valence_carbon = inp_ATOM {
            el_number: 6,
            valence: 1,
            ..inp_ATOM::default()
        };
        let mut rejected_limits = [value(4, 0, 3), value(2, 2, 2), value(1, 0, 1)];
        assert_eq!(
            run(
                &mut rejected_limits,
                std::slice::from_ref(&low_valence_carbon),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!(rejected_limits[0], value(1, 0, 1));

        let tautomeric_carbon = inp_ATOM {
            endpoint: 1,
            ..low_valence_carbon.clone()
        };
        let mut tautomeric = [value(1, -1, 2), value(1, 1, 1), value(1, 0, 0)];
        assert_eq!(
            run(
                &mut tautomeric,
                std::slice::from_ref(&tautomeric_carbon),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!(tautomeric[0], value(1, 0, 0));

        let oxygen = inp_ATOM {
            el_number: 8,
            valence: 1,
            ..inp_ATOM::default()
        };
        let oxygen_valence = crate::source_types::VAL_AT {
            cNumValenceElectrons: 6,
            ..crate::source_types::VAL_AT::default()
        };
        let endpoints = [1_u16];
        let mut first_fixed_h_negative = [value(1, -1, 0)];
        assert_eq!(
            run(
                &mut first_fixed_h_negative,
                std::slice::from_ref(&oxygen),
                std::slice::from_ref(&oxygen_valence),
                0,
                0,
                Some(&endpoints),
            ),
            Ok(1)
        );
        let mut later_fixed_h_charges = [value(1, -1, 2), value(1, 1, 1), value(1, 0, 0)];
        assert_eq!(
            run(
                &mut later_fixed_h_charges,
                std::slice::from_ref(&oxygen),
                std::slice::from_ref(&oxygen_valence),
                0,
                0,
                Some(&endpoints),
            ),
            Ok(1)
        );
        assert_eq!(later_fixed_h_charges[0], value(1, 0, 0));

        let mobile_nitrogen = inp_ATOM {
            el_number: 7,
            valence: 1,
            num_H: 1,
            ..inp_ATOM::default()
        };
        let mut mobile_h_pair = [value(2, -1, 2), value(2, 1, 1), value(1, 0, 0)];
        assert_eq!(
            run(
                &mut mobile_h_pair,
                std::slice::from_ref(&mobile_nitrogen),
                std::slice::from_ref(&default_valence),
                0,
                1,
                None,
            ),
            Ok(1)
        );
        assert_eq!(mobile_h_pair[0], value(1, 0, 0));

        let mut carbon_pairs = [value(2, -1, 3), value(2, 1, 2), value(1, -1, 1), value(1, 1, 0)];
        assert_eq!(
            run(
                &mut carbon_pairs,
                std::slice::from_ref(&carbon),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(3)
        );
        assert_eq!(&carbon_pairs[..3], &[value(1, 1, 0), value(1, -1, 1), value(2, 1, 2)]);
        let noncarbon = inp_ATOM {
            el_number: 8,
            valence: 2,
            ..inp_ATOM::default()
        };
        let mut noncarbon_pair = [value(1, -1, 1), value(1, 1, 0)];
        assert_eq!(
            run(
                &mut noncarbon_pair,
                std::slice::from_ref(&noncarbon),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!(noncarbon_pair[0], value(1, 1, 0));

        let fixed_h_nitrogen = inp_ATOM {
            el_number: 7,
            valence: 2,
            num_H: 1,
            ..inp_ATOM::default()
        };
        let mut nitrogen_valence_five = [value(5, 0, 1), value(4, 0, 0)];
        assert_eq!(
            run(
                &mut nitrogen_valence_five,
                std::slice::from_ref(&fixed_h_nitrogen),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!(nitrogen_valence_five[0], value(4, 0, 0));

        let mut valence_gap = [value(3, 0, 1), value(1, 0, 0)];
        assert_eq!(
            run(
                &mut valence_gap,
                std::slice::from_ref(&carbon),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!(valence_gap[0], value(1, 0, 0));

        let mut four_without_opposite_charges = [value(3, 0, 3), value(2, 0, 2), value(1, 0, 1), value(0, 0, 0)];
        assert_eq!(
            run(
                &mut four_without_opposite_charges,
                std::slice::from_ref(&carbon),
                std::slice::from_ref(&default_valence),
                0,
                0,
                None,
            ),
            Ok(3)
        );
        assert_eq!(
            &four_without_opposite_charges[..4],
            &[value(0, 0, 0), value(1, 0, 1), value(2, 0, 2), value(3, 0, 3),]
        );
    }

    #[test]
    fn source_port__ichirvr1__getatomrestoreinfo__line_912() {
        let run = |atoms: &mut [inp_ATOM],
                   valence_atoms: &mut [crate::source_types::VAL_AT],
                   restore_mode: &SRM,
                   mobile_h: i32,
                   endpoint: Option<&[AT_NUMB]>| {
            GetAtomRestoreInfo(
                &mut CANON_GLOBALS::default(),
                atoms,
                0,
                valence_atoms,
                restore_mode,
                mobile_h,
                endpoint,
            )
        };

        assert_eq!(
            GetAtomRestoreInfo(
                &mut CANON_GLOBALS::default(),
                &mut [],
                0,
                &mut [],
                &SRM::default(),
                0,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut bridging_hydrogen = inp_ATOM {
            el_number: 1,
            valence: 1,
            chem_bonds_valence: 5,
            ..inp_ATOM::default()
        };
        bridging_hydrogen.neighbor[0] = 1;
        bridging_hydrogen.bond_type[0] = BOND_TYPE_TRIPLE as u8;
        let mut hydrogen_atoms = [bridging_hydrogen, inp_ATOM::default()];
        let mut hydrogen_valence = [
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT {
                cMetal: 1,
                ..crate::source_types::VAL_AT::default()
            },
        ];
        assert_eq!(
            run(&mut hydrogen_atoms, &mut hydrogen_valence, &SRM::default(), 0, None,),
            Ok(0)
        );
        assert_eq!(hydrogen_atoms[0].chem_bonds_valence, 4);
        assert_eq!(hydrogen_valence[0].cNumBondsToMetal, 1);
        assert_eq!(hydrogen_valence[0].cDoNotAddH, 0);

        let mut isolated = [inp_ATOM {
            el_number: 6,
            ..inp_ATOM::default()
        }];
        let mut isolated_valence = [crate::source_types::VAL_AT {
            cNumValenceElectrons: 4,
            ..crate::source_types::VAL_AT::default()
        }];
        assert_eq!(
            run(&mut isolated, &mut isolated_valence, &SRM::default(), 0, None,),
            Ok(0)
        );
        assert_eq!(isolated_valence[0].cDoNotAddH, 0);
        assert_eq!(isolated_valence[0].cnListIndex, 0);

        let mut mixed = inp_ATOM {
            el_number: 6,
            valence: 2,
            chem_bonds_valence: 5,
            ..inp_ATOM::default()
        };
        mixed.neighbor[0] = 1;
        mixed.neighbor[1] = 2;
        mixed.bond_type[0] = BOND_TYPE_DOUBLE as u8;
        mixed.bond_type[1] = crate::source_types::BOND_TYPE_ALTERN as u8;
        let mut mixed_atoms = [mixed, inp_ATOM::default(), inp_ATOM::default()];
        let mut mixed_valence = [
            crate::source_types::VAL_AT {
                cNumValenceElectrons: 4,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT {
                cMetal: 1,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT::default(),
        ];
        assert_eq!(
            run(&mut mixed_atoms, &mut mixed_valence, &SRM::default(), 0, None,),
            Ok(1)
        );
        assert_eq!(mixed_atoms[0].chem_bonds_valence, 4);
        assert_eq!(mixed_valence[0].cNumBondsToMetal, 1);
        assert_eq!(mixed_valence[0].cInitOrigValenceToMetal, 2);
        assert_eq!(mixed_valence[0].cInitValenceToMetal, 2);
        assert_eq!(mixed_valence[0].cInitFlowToMetal, 1);
        assert_eq!(mixed_valence[0].cInitFreeValences, 0);
        assert_eq!(mixed_valence[0].cInitCharge, 0);
        assert_eq!(mixed_valence[0].cnListIndex, 17);

        let flower_mode = SRM {
            bMetalAddFlower: 1,
            nMetalInitBondOrder: 1,
            nMetalMinBondOrder: 0,
            nMetalInitEdgeFlow: 0,
            ..SRM::default()
        };
        let mut low_flow_metal = inp_ATOM {
            el_number: 26,
            valence: 2,
            chem_bonds_valence: 3,
            ..inp_ATOM::default()
        };
        low_flow_metal.bond_type[0] = BOND_TYPE_SINGLE as u8;
        low_flow_metal.bond_type[1] = crate::source_types::BOND_TYPE_ALTERN as u8;
        let mut low_flow_atoms = [low_flow_metal];
        let mut low_flow_valence = [crate::source_types::VAL_AT {
            cMetal: 1,
            cInitFreeValences: 5,
            ..crate::source_types::VAL_AT::default()
        }];
        assert_eq!(
            run(&mut low_flow_atoms, &mut low_flow_valence, &flower_mode, 0, None,),
            Ok(0)
        );
        assert_eq!(low_flow_valence[0].cNumBondsToMetal, 2);
        assert_eq!(low_flow_valence[0].cInitOrigValenceToMetal, 2);
        assert_eq!(low_flow_valence[0].cInitValenceToMetal, 2);
        assert_eq!(low_flow_valence[0].cInitFlowToMetal, 0);
        assert_eq!(low_flow_valence[0].cInitFreeValences, 8);
        assert_eq!(low_flow_valence[0].cnListIndex, 18);

        let mut high_flow_metal = inp_ATOM {
            el_number: 26,
            valence: 2,
            chem_bonds_valence: 6,
            ..inp_ATOM::default()
        };
        high_flow_metal.bond_type[0] = BOND_TYPE_TRIPLE as u8;
        high_flow_metal.bond_type[1] = BOND_TYPE_TRIPLE as u8;
        let mut high_flow_atoms = [high_flow_metal];
        let mut high_flow_valence = [crate::source_types::VAL_AT {
            cMetal: 1,
            cInitFreeValences: 5,
            ..crate::source_types::VAL_AT::default()
        }];
        assert_eq!(
            run(&mut high_flow_atoms, &mut high_flow_valence, &flower_mode, 0, None,),
            Ok(0)
        );
        assert_eq!(high_flow_valence[0].cInitOrigValenceToMetal, 6);
        assert_eq!(high_flow_valence[0].cInitValenceToMetal, 6);
        assert_eq!(high_flow_valence[0].cInitFlowToMetal, 4);
        assert_eq!(high_flow_valence[0].cInitFreeValences, 7);
        assert_eq!(high_flow_valence[0].cnListIndex, 18);

        let mut neon = inp_ATOM {
            el_number: 10,
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        neon.bond_type[0] = BOND_TYPE_SINGLE as u8;
        let mut neon_atoms = [neon];
        let mut neon_valence = [crate::source_types::VAL_AT::default()];
        assert_eq!(
            run(&mut neon_atoms, &mut neon_valence, &SRM::default(), 0, None,),
            Ok(crate::source_types::TREAT_ATOM_AS_METAL as i32)
        );
        assert_eq!(neon_valence[0].cInitFreeValences, 0);
        assert_eq!(neon_valence[0].cInitFlowToMetal, 0);
        assert_eq!(neon_valence[0].cnListIndex, 0);

        let mut methane_carbon = inp_ATOM {
            el_number: 6,
            valence: 4,
            chem_bonds_valence: 4,
            ..inp_ATOM::default()
        };
        methane_carbon.bond_type[..4].fill(BOND_TYPE_SINGLE as u8);
        let mut methane_atoms = [
            methane_carbon,
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        methane_atoms[0].neighbor[..4].copy_from_slice(&[1, 2, 3, 4]);
        let mut methane_valence = [
            crate::source_types::VAL_AT {
                cNumValenceElectrons: 4,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
        ];
        assert_eq!(
            run(&mut methane_atoms, &mut methane_valence, &SRM::default(), 0, None,),
            Ok(1)
        );
        assert_eq!(methane_valence[0].cnListIndex, 17);
        assert_eq!(methane_valence[0].cInitCharge, 0);
        assert_eq!(methane_valence[0].cInitFreeValences, 0);

        let mut terminal_carbon = inp_ATOM {
            el_number: 6,
            valence: 3,
            chem_bonds_valence: 3,
            ..inp_ATOM::default()
        };
        terminal_carbon.bond_type[..3].fill(BOND_TYPE_SINGLE as u8);
        terminal_carbon.neighbor[..3].copy_from_slice(&[1, 2, 3]);
        let mut terminal_atoms = [
            terminal_carbon,
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        let mut terminal_valence = [
            crate::source_types::VAL_AT {
                cNumValenceElectrons: 4,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
        ];
        assert_eq!(
            run(&mut terminal_atoms, &mut terminal_valence, &SRM::default(), 0, None,),
            Ok(1)
        );
        assert_eq!(terminal_valence[0].cnListIndex, 7);
        assert_eq!(terminal_valence[0].cInitCharge, 0);
        assert_eq!(terminal_valence[0].cInitFreeValences, 0);

        let mut isocyano_carbon = inp_ATOM {
            el_number: 6,
            valence: 1,
            chem_bonds_valence: 3,
            ..inp_ATOM::default()
        };
        isocyano_carbon.bond_type[0] = BOND_TYPE_TRIPLE as u8;
        isocyano_carbon.neighbor[0] = 1;
        let mut isocyano_atoms = [isocyano_carbon, inp_ATOM::default()];
        let mut isocyano_valence = [
            crate::source_types::VAL_AT {
                cNumValenceElectrons: 4,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT::default(),
        ];
        assert_eq!(
            run(&mut isocyano_atoms, &mut isocyano_valence, &SRM::default(), 0, None,),
            Ok(1)
        );
        assert_eq!(isocyano_valence[0].cnListIndex, 7);
        assert_eq!(isocyano_valence[0].cInitCharge, 0);
        assert_eq!(isocyano_valence[0].cInitFreeValences, 0);

        let mut overbonded_carbon = inp_ATOM {
            el_number: 6,
            valence: 6,
            chem_bonds_valence: 6,
            ..inp_ATOM::default()
        };
        overbonded_carbon.bond_type[..6].fill(BOND_TYPE_SINGLE as u8);
        overbonded_carbon.neighbor[..6].copy_from_slice(&[1, 2, 3, 4, 5, 6]);
        let mut overbonded_atoms = [
            overbonded_carbon,
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        let mut overbonded_valence = [
            crate::source_types::VAL_AT {
                cNumValenceElectrons: 4,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
            crate::source_types::VAL_AT::default(),
        ];
        assert_eq!(
            run(&mut overbonded_atoms, &mut overbonded_valence, &SRM::default(), 0, None,),
            Ok(crate::source_types::TREAT_ATOM_AS_METAL as i32)
        );
        assert_eq!(overbonded_valence[0].cnListIndex, 0);
    }

    #[test]
    fn source_port__ichirvr1__bmaybeacationinmobilehlayer__line_749() {
        let valence_atoms = vec![crate::source_types::VAL_AT::default(); 2];
        assert_eq!(
            bMayBeACationInMobileHLayer(
                &[inp_ATOM {
                    el_number: 7,
                    num_H: 1,
                    ..inp_ATOM::default()
                }],
                &valence_atoms,
                0,
                0,
            ),
            Ok(1)
        );
        assert_eq!(
            bMayBeACationInMobileHLayer(
                &[inp_ATOM {
                    el_number: 7,
                    ..inp_ATOM::default()
                }],
                &valence_atoms,
                0,
                1,
            ),
            Ok(1)
        );
        assert_eq!(
            bMayBeACationInMobileHLayer(
                &[inp_ATOM {
                    el_number: 6,
                    num_H: 1,
                    ..inp_ATOM::default()
                }],
                &valence_atoms,
                0,
                1,
            ),
            Ok(1)
        );

        for (element, maximum) in [(7, 4), (15, 4), (8, 3), (16, 3), (34, 3), (52, 3)] {
            assert_eq!(
                bMayBeACationInMobileHLayer(
                    &[inp_ATOM {
                        el_number: element,
                        valence: (maximum - 1) as i8,
                        num_H: 1,
                        ..inp_ATOM::default()
                    }],
                    &valence_atoms,
                    0,
                    1,
                ),
                Ok(0)
            );
            assert_eq!(
                bMayBeACationInMobileHLayer(
                    &[inp_ATOM {
                        el_number: element,
                        valence: maximum as i8,
                        num_H: 1,
                        ..inp_ATOM::default()
                    }],
                    &valence_atoms,
                    0,
                    1,
                ),
                Ok(1)
            );
        }

        let mut center = inp_ATOM {
            el_number: 7,
            valence: 1,
            num_H: 1,
            ..inp_ATOM::default()
        };
        center.neighbor[0] = 1;
        let special_neighbor = inp_ATOM {
            valence: 4,
            chem_bonds_valence: 4,
            ..inp_ATOM::default()
        };
        let special_valence = crate::source_types::VAL_AT {
            cNumValenceElectrons: 3,
            cPeriodicRowNumber: 1,
            ..crate::source_types::VAL_AT::default()
        };
        assert_eq!(
            bMayBeACationInMobileHLayer(
                &[center.clone(), special_neighbor.clone()],
                &[crate::source_types::VAL_AT::default(), special_valence.clone()],
                0,
                1,
            ),
            Ok(1)
        );

        let mut rejected_neighbors = Vec::new();
        let mut neighbor = special_neighbor.clone();
        neighbor.valence = 3;
        rejected_neighbors.push((neighbor, special_valence.clone()));
        let mut neighbor = special_neighbor.clone();
        neighbor.chem_bonds_valence = 3;
        rejected_neighbors.push((neighbor, special_valence.clone()));
        let mut neighbor = special_neighbor.clone();
        neighbor.num_H = 1;
        rejected_neighbors.push((neighbor, special_valence.clone()));
        let mut valence = special_valence.clone();
        valence.cNumValenceElectrons = 4;
        rejected_neighbors.push((special_neighbor.clone(), valence));
        let mut valence = special_valence;
        valence.cPeriodicRowNumber = 2;
        rejected_neighbors.push((special_neighbor, valence));
        for (neighbor, valence) in rejected_neighbors {
            assert_eq!(
                bMayBeACationInMobileHLayer(
                    &[center.clone(), neighbor],
                    &[crate::source_types::VAL_AT::default(), valence],
                    0,
                    1,
                ),
                Ok(0)
            );
        }

        let negative_valence = inp_ATOM {
            el_number: 7,
            valence: -1,
            num_H: 1,
            ..inp_ATOM::default()
        };
        assert_eq!(bMayBeACationInMobileHLayer(&[negative_valence], &[], 0, 1), Ok(0));
    }

    #[test]
    fn source_port__ichirvr1__filloutpstructendpointfrominchi__line_692() {
        let mut heap = SourceHeap::default();
        let existing = heap.allocate_model_storage(vec![11_u16, 13, 17]).unwrap();
        let mut endpoint = existing;
        assert_eq!(
            FillOutpStructEndpointFromInChI(
                &mut heap,
                &INChI {
                    nNumberOfAtoms: 3,
                    lenTautomer: 1,
                    ..INChI::default()
                },
                &mut endpoint,
            ),
            Ok(0)
        );
        assert_eq!(endpoint, existing);
        assert_eq!(heap.slice(endpoint.as_const()).unwrap(), &[0, 0, 0]);

        let encoded_groups = heap
            .allocate_model_storage(vec![2_u16, 4, 1, 0, 1, 2, 3, 0, 0, 3])
            .unwrap();
        let mut allocated = SourceMutPointer::null();
        assert_eq!(
            FillOutpStructEndpointFromInChI(
                &mut heap,
                &INChI {
                    nNumberOfAtoms: 3,
                    lenTautomer: 10,
                    nTautomer: encoded_groups,
                    ..INChI::default()
                },
                &mut allocated,
            ),
            Ok(0)
        );
        assert!(!allocated.is_null());
        assert_eq!(heap.slice(allocated.as_const()).unwrap(), &[1, 1, 2]);

        let mut failure_heap = SourceHeap::default();
        failure_heap.fail_after_allocations(0);
        let mut failed = SourceMutPointer::null();
        assert_eq!(
            FillOutpStructEndpointFromInChI(
                &mut failure_heap,
                &INChI {
                    nNumberOfAtoms: 3,
                    ..INChI::default()
                },
                &mut failed,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(failed.is_null());
        assert_eq!(failure_heap.source_allocation_calls(), 1);

        inchi_free(&mut heap, allocated).unwrap();
        inchi_free(&mut heap, existing).unwrap();
        inchi_free(&mut heap, encoded_groups).unwrap();
        assert_eq!(heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__ichirvr1__makeoneinchioutofstrfrominchi2__line_5087() {
        let mut heap = SourceHeap::default();
        let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()]).unwrap();
        let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
        let vertex = heap.allocate_model_storage(vec![BNS_VERTEX::default()]).unwrap();
        let bns = BN_STRUCT {
            num_atoms: 1,
            num_vertices: 1,
            vert: vertex,
            ..BN_STRUCT::default()
        };
        let mut carbon = inp_ATOM {
            el_number: 6,
            num_H: 4,
            charge: 7,
            orig_at_number: 1,
            endpoint: 8,
            component: 9,
            ..inp_ATOM::default()
        };
        carbon.elname[0] = b'C' as i8;
        let at = heap.allocate_model_storage(vec![carbon.clone()]).unwrap();
        let at2 = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let at3 = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let mut structure = StrFromINChI {
            at,
            num_atoms: 1,
            bMobileH: TAUT_YES as i8,
            iMobileH: TAUT_YES as i8,
            pSrm: restore_mode.as_const(),
            ..StrFromINChI::default()
        };
        let input_parameters = INPUT_PARMS {
            nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
            bTautFlags: u64::from(
                crate::source_types::TG_FLAG_FIX_ISO_FIXEDH_BUG | crate::source_types::TG_FLAG_FIX_TERM_H_CHRG_BUG,
            ),
            ..INPUT_PARMS::default()
        };
        let structure_data = STRUCT_DATA {
            ulStructTime: 123,
            nErrorCode: 17,
            nErrorType: 19,
            pStrErrStruct: [23; 256],
            ..STRUCT_DATA::default()
        };
        let original_parameters = input_parameters.clone();
        let original_structure_data = structure_data.clone();
        let mut canonical_globals = CANON_GLOBALS::default();
        let mut normalized_output = SourceMutPointer::null();
        let mut prepared_output = at;
        let mut tgroup_output = SourceTGroupInfoPointer::External(99);

        assert_eq!(
            MakeOneInChIOutOfStrFromINChI2(
                &mut heap,
                &mut canonical_globals,
                clock,
                &input_parameters,
                &structure_data,
                &bns,
                &mut structure,
                at,
                at2,
                at3,
                &[VAL_AT::default()],
                &ALL_TC_GROUPS::default(),
                Some(&mut tgroup_output),
                Some(&mut normalized_output),
                Some(&mut prepared_output),
                0,
            ),
            Ok(0)
        );
        assert_eq!(tgroup_output, SourceTGroupInfoPointer::Null);
        assert_eq!(input_parameters, original_parameters);
        assert_eq!(structure_data, original_structure_data);
        assert_eq!(structure.at, at);
        assert_eq!(heap.slice(at.as_const()).unwrap(), &[carbon.clone()]);
        assert_eq!(heap.slice(at2.as_const()).unwrap()[0].charge, 0);
        assert_eq!(heap.slice(at2.as_const()).unwrap()[0].endpoint, 8);
        assert_eq!(heap.slice(at3.as_const()).unwrap()[0].endpoint, 0);
        assert_eq!(heap.slice(at3.as_const()).unwrap()[0].component, 1);
        let normalized_holder = heap.slice(structure.pOne_norm_data[0].as_const()).unwrap()[0].clone();
        assert_eq!(normalized_output, normalized_holder.at);
        assert!(prepared_output.is_null());
        assert!(structure.One_ti.num_t_groups == 0);

        let mut error_heap = SourceHeap::default();
        let group_vertex = error_heap.allocate_model_storage(vec![BNS_VERTEX::default()]).unwrap();
        let group = error_heap.allocate_model_storage(vec![TC_GROUP::default()]).unwrap();
        let groups = ALL_TC_GROUPS {
            pTCG: group,
            num_tgroups: 1,
            ..ALL_TC_GROUPS::default()
        };
        let invalid_bns = BN_STRUCT {
            num_t_groups: 1,
            num_vertices: 1,
            vert: group_vertex,
            ..BN_STRUCT::default()
        };
        let restored_at = error_heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let sentinel = error_heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let mut error_structure = StrFromINChI {
            at: sentinel,
            ..StrFromINChI::default()
        };
        let mut normalized_output = sentinel;
        let mut prepared_output = sentinel;
        let mut tgroup_output = SourceTGroupInfoPointer::External(77);
        assert_eq!(
            MakeOneInChIOutOfStrFromINChI2(
                &mut error_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                &invalid_bns,
                &mut error_structure,
                restored_at,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &[],
                &groups,
                Some(&mut tgroup_output),
                Some(&mut normalized_output),
                Some(&mut prepared_output),
                0,
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(normalized_output, sentinel);
        assert_eq!(prepared_output, sentinel);
        assert_eq!(tgroup_output, SourceTGroupInfoPointer::External(77));
        assert_eq!(error_structure.at, restored_at);

        let mut failure_heap = SourceHeap::default();
        let clock = failure_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let restore_mode = failure_heap.allocate_model_storage(vec![SRM::default()]).unwrap();
        let vertex = failure_heap
            .allocate_model_storage(vec![BNS_VERTEX::default()])
            .unwrap();
        let bns = BN_STRUCT {
            num_atoms: 1,
            num_vertices: 1,
            vert: vertex,
            ..BN_STRUCT::default()
        };
        let at = failure_heap.allocate_model_storage(vec![carbon]).unwrap();
        let at2 = failure_heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let at3 = failure_heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let mut structure = StrFromINChI {
            at,
            num_atoms: 1,
            bMobileH: TAUT_YES as i8,
            iMobileH: TAUT_YES as i8,
            pSrm: restore_mode.as_const(),
            ..StrFromINChI::default()
        };
        let mut normalized_output = at;
        let mut prepared_output = at;
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            MakeOneInChIOutOfStrFromINChI2(
                &mut failure_heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &input_parameters,
                &structure_data,
                &bns,
                &mut structure,
                at,
                at2,
                at3,
                &[VAL_AT::default()],
                &ALL_TC_GROUPS::default(),
                None,
                Some(&mut normalized_output),
                Some(&mut prepared_output),
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(structure.at, at);
        assert_eq!(normalized_output, at);
        assert_eq!(prepared_output, at);
    }

    #[test]
    fn source_port__ichirvr1__makeoneinchioutofstrfrominchi__line_5168() {
        let mut heap = SourceHeap::default();
        let mut source_atom = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 0,
            endpoint: 7,
            component: 9,
            ..inp_ATOM::default()
        };
        source_atom.neighbor[0] = 0;
        source_atom.bond_type[0] = 0;
        let at2 = heap.allocate_model_storage(vec![source_atom.clone()]).unwrap();
        let at3 = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();

        let old_inchi = heap.allocate_model_storage(vec![INChI::default()]).unwrap();
        let old_aux = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Aux::default()])
            .unwrap();
        let old_normalized_atoms = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let old_normalized = heap
            .allocate_model_storage(vec![INP_ATOM_DATA {
                at: old_normalized_atoms,
                num_at: 1,
                ..INP_ATOM_DATA::default()
            }])
            .unwrap();
        let old_groups = heap.allocate_model_storage(vec![T_GROUP::default()]).unwrap();
        let mut structure = StrFromINChI {
            num_atoms: 1,
            pOneINChI: [old_inchi, SourceMutPointer::null()],
            pOneINChI_Aux: [old_aux, SourceMutPointer::null()],
            pOne_norm_data: [old_normalized, SourceMutPointer::null()],
            One_ti: T_GROUP_INFO {
                t_group: old_groups,
                max_num_t_groups: 1,
                ..T_GROUP_INFO::default()
            },
            ..StrFromINChI::default()
        };
        let mut structure_data = STRUCT_DATA {
            ulStructTime: 123,
            nErrorCode: 17,
            nErrorType: 19,
            pStrErrStruct: [23; 256],
            ..STRUCT_DATA::default()
        };
        let groups = ALL_TC_GROUPS {
            iComponent: 3,
            ..ALL_TC_GROUPS::default()
        };
        let mut canonical_globals = CANON_GLOBALS::default();

        heap.fail_after_allocations(0);
        assert_eq!(
            MakeOneInChIOutOfStrFromINChI(
                &mut heap,
                &mut canonical_globals,
                SourceMutPointer::null(),
                &INPUT_PARMS::default(),
                &mut structure_data,
                &mut structure,
                at2,
                at3,
                &groups,
                0,
            ),
            Ok(-1)
        );

        assert_eq!(heap.source_allocation_calls(), 1);
        assert_eq!(heap.live_allocation_count(), 2);
        assert_eq!(heap.slice(at2.as_const()).unwrap(), &[source_atom]);
        let prepared = &heap.slice(at3.as_const()).unwrap()[0];
        assert_eq!(prepared.endpoint, 0);
        assert_eq!(prepared.component, 4);
        assert_eq!(prepared.bond_type[0], BOND_TYPE_SINGLE as u8);
        assert_eq!(prepared.chem_bonds_valence, BOND_TYPE_SINGLE as i8);
        assert!(structure.pOneINChI[0].is_null());
        assert!(structure.pOneINChI_Aux[0].is_null());
        assert!(structure.pOne_norm_data[0].is_null());
        assert_eq!(structure.One_ti, T_GROUP_INFO::default());
        assert_eq!(structure_data.ulStructTime, 123);
        assert_eq!(structure_data.nErrorCode, 0);
        assert_eq!(structure_data.nErrorType, 0);
        assert_eq!(structure_data.pStrErrStruct, [0; 256]);

        let mut heap = SourceHeap::default();
        let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()]).unwrap();
        let mut carbon = inp_ATOM {
            el_number: 6,
            num_H: 4,
            orig_at_number: 1,
            endpoint: 29,
            component: 31,
            ..inp_ATOM::default()
        };
        carbon.elname[0] = b'C' as i8;
        let at2 = heap.allocate_model_storage(vec![carbon.clone()]).unwrap();
        let at3 = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let mut structure = StrFromINChI {
            num_atoms: 1,
            bMobileH: TAUT_YES as i8,
            iMobileH: TAUT_YES as i8,
            ..StrFromINChI::default()
        };
        let input_parameters = INPUT_PARMS {
            nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
            bTautFlags: u64::from(
                crate::source_types::TG_FLAG_FIX_ISO_FIXEDH_BUG | crate::source_types::TG_FLAG_FIX_TERM_H_CHRG_BUG,
            ),
            ..INPUT_PARMS::default()
        };
        let mut structure_data = STRUCT_DATA {
            ulStructTime: 456,
            ..STRUCT_DATA::default()
        };
        let mut canonical_globals = CANON_GLOBALS::default();

        assert_eq!(
            MakeOneInChIOutOfStrFromINChI(
                &mut heap,
                &mut canonical_globals,
                clock,
                &input_parameters,
                &mut structure_data,
                &mut structure,
                at2,
                at3,
                &ALL_TC_GROUPS::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(structure_data.ulStructTime, 456);
        assert_eq!(structure.nOneINChI_bMobileH, TAUT_YES as i32);
        assert_eq!(structure.nNumRemovedProtons, 0);
        assert_eq!(structure.nChargeRevrs, 0);
        assert!(!structure.pOneINChI[0].is_null());
        assert!(!structure.pOneINChI_Aux[0].is_null());
        assert!(!structure.pOne_norm_data[0].is_null());
        assert!(structure.pOneINChI[1].is_null());
        assert!(structure.pOneINChI_Aux[1].is_null());
        assert!(structure.pOne_norm_data[1].is_null());

        let generated = &heap.slice(structure.pOneINChI[0].as_const()).unwrap()[0];
        assert_eq!(generated.nErrorCode, 0);
        assert_eq!(generated.nNumberOfAtoms, 1);
        assert_eq!(heap.slice(generated.nAtom.as_const()).unwrap(), &[6]);
        assert_eq!(heap.slice(generated.nNum_H.as_const()).unwrap(), &[4]);
        let generated_aux = &heap.slice(structure.pOneINChI_Aux[0].as_const()).unwrap()[0];
        assert_eq!(generated_aux.nErrorCode, 0);
        assert_eq!(generated_aux.nNumberOfAtoms, 1);
        let normalized = &heap.slice(structure.pOne_norm_data[0].as_const()).unwrap()[0];
        assert_eq!(normalized.num_at, 1);
        assert_eq!(normalized.num_removed_H, 0);
        assert_eq!(heap.slice(normalized.at.as_const()).unwrap()[0].el_number, 6);
        assert_eq!(heap.slice(at2.as_const()).unwrap(), &[carbon]);
        assert_eq!(heap.slice(at3.as_const()).unwrap()[0].component, 1);

        Free_INChI(&mut heap, &mut structure.pOneINChI[0]).unwrap();
        Free_INChI_Aux(&mut heap, &mut structure.pOneINChI_Aux[0]).unwrap();
        let normalized_pointer = structure.pOne_norm_data[0];
        let mut normalized_owner = heap.slice(normalized_pointer.as_const()).unwrap()[0].clone();
        FreeInpAtomData(&mut heap, &mut normalized_owner).unwrap();
        inchi_free(&mut heap, normalized_pointer).unwrap();
        structure.pOne_norm_data[0] = SourceMutPointer::null();
        free_t_group_info(&mut heap, Some(&mut structure.One_ti)).unwrap();
        assert!(!canonical_globals.m_bBit.is_null());
        assert_eq!(
            crate::source::base::ichican2::SetBitFree(&mut heap, &mut canonical_globals),
            Ok(1)
        );
        assert_eq!(heap.live_allocation_count(), 3);
    }

    fn bond_flow_fixture(
        atom_metal: i8,
        neighbor_metal: i8,
        atom_endpoint: u16,
        neighbor_endpoint: u16,
        bond_type: u8,
    ) -> (Vec<inp_ATOM>, Vec<crate::source_types::VAL_AT>) {
        let mut atom = inp_ATOM {
            endpoint: atom_endpoint,
            valence: 1,
            ..inp_ATOM::default()
        };
        atom.neighbor[0] = 1;
        atom.bond_type[0] = bond_type;
        (
            vec![
                atom,
                inp_ATOM {
                    endpoint: neighbor_endpoint,
                    ..inp_ATOM::default()
                },
            ],
            vec![
                crate::source_types::VAL_AT {
                    cMetal: atom_metal,
                    ..crate::source_types::VAL_AT::default()
                },
                crate::source_types::VAL_AT {
                    cMetal: neighbor_metal,
                    ..crate::source_types::VAL_AT::default()
                },
            ],
        )
    }

    #[test]
    fn source_port__ichirvr1__bondflowmaxcapminorder__line_1690() {
        let mode = SRM {
            bMetalAddFlower: 1,
            nMetalMinBondOrder: 0,
            nMetalInitEdgeFlow: 0,
            nMetalInitBondOrder: 1,
            nMetal2EndpointMinBondOrder: 0,
            nMetal2EndpointInitBondOrder: 2,
            nMetal2EndpointInitEdgeFlow: 0,
            bFixStereoBonds: 1,
            ..SRM::default()
        };

        let (atoms, va) = bond_flow_fixture(0, 0, 0, 0, BOND_TYPE_DOUBLE as u8);
        let (mut max, mut min, mut flower) = (-9, -9, -9);
        assert_eq!(
            BondFlowMaxcapMinorder(
                &atoms,
                &va,
                &mode,
                0,
                0,
                Some(&mut max),
                Some(&mut min),
                Some(&mut flower),
            ),
            Ok(1)
        );
        assert_eq!((max, min, flower), (2, 1, 0));

        let (mut atoms, va) = bond_flow_fixture(1, 0, 0, 0, 0xf4);
        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &va, &mode, 0, 0, None, None, None),
            Ok(0)
        );
        atoms[0].bond_type[0] = 0xf3;
        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &va, &mode, 0, 0, None, None, None),
            Ok(2)
        );

        let (mut atoms, va) = bond_flow_fixture(1, 0, 0, 0, BOND_TYPE_DOUBLE as u8);
        atoms[1].sb_parity[0] = 1;
        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &va, &mode, 0, 0, None, None, None),
            Ok(1)
        );
        let mut no_fix = mode.clone();
        no_fix.bFixStereoBonds = 0;
        let (mut max, mut min, mut flower) = (-1, -1, -1);
        assert_eq!(
            BondFlowMaxcapMinorder(
                &atoms,
                &va,
                &no_fix,
                0,
                0,
                Some(&mut max),
                Some(&mut min),
                Some(&mut flower),
            ),
            Ok(1)
        );
        assert_eq!((max, min, flower), (3, 0, 1));

        let mut disabled = mode.clone();
        disabled.bMetalAddFlower = 0;
        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &va, &disabled, 0, 0, None, None, None),
            Ok(1)
        );

        for (atom_metal, neighbor_metal, atom_endpoint, neighbor_endpoint, expected_flower) in
            [(1, 0, 0, 1, 1), (0, 1, 1, 0, 0), (1, 0, 1, 0, 0), (1, 1, 0, 1, 1)]
        {
            let (atoms, va) = bond_flow_fixture(
                atom_metal,
                neighbor_metal,
                atom_endpoint,
                neighbor_endpoint,
                BOND_TYPE_SINGLE as u8,
            );
            let (mut max, mut min, mut flower) = (-1, -1, -1);
            assert_eq!(
                BondFlowMaxcapMinorder(
                    &atoms,
                    &va,
                    &mode,
                    0,
                    0,
                    Some(&mut max),
                    Some(&mut min),
                    Some(&mut flower),
                ),
                Ok(1)
            );
            assert_eq!((max, min, flower), (3, 0, expected_flower));
        }

        let mut nonzero_edge_flow = mode.clone();
        nonzero_edge_flow.nMetalInitEdgeFlow = 1;
        let (atoms, va) = bond_flow_fixture(1, 1, 0, 0, BOND_TYPE_SINGLE as u8);
        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &va, &nonzero_edge_flow, 0, 0, None, None, None,),
            Ok(1)
        );

        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &va, &mode, -1, 0, None, None, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &va, &mode, 0, -1, None, None, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            BondFlowMaxcapMinorder(&[], &va, &mode, 0, 0, None, None, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            BondFlowMaxcapMinorder(&atoms, &[], &mode, 0, 0, None, None, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr1__atomstcapstflow__line_1790() {
        let mut atom = inp_ATOM {
            valence: 2,
            chem_bonds_valence: 5,
            ..inp_ATOM::default()
        };
        atom.neighbor[..2].copy_from_slice(&[1, 2]);
        atom.bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_TYPE_TRIPLE as u8]);
        let atoms = vec![atom, inp_ATOM::default(), inp_ATOM::default()];
        let valence_atoms = vec![crate::source_types::VAL_AT::default(); 3];
        let disabled = SRM::default();
        let (mut st_cap, mut st_flow, mut metal_cap, mut metal_flow) = (-1, -1, -1, -1);
        assert_eq!(
            AtomStcapStflow(
                &atoms,
                &valence_atoms,
                &disabled,
                0,
                Some(&mut st_cap),
                Some(&mut st_flow),
                Some(&mut metal_cap),
                Some(&mut metal_flow),
            ),
            Ok(0)
        );
        assert_eq!((st_cap, st_flow, metal_cap, metal_flow), (3, 3, 0, 0));

        let mut metal_atom = atoms[0].clone();
        metal_atom.bond_type[..2].copy_from_slice(&[BOND_TYPE_SINGLE as u8, BOND_TYPE_DOUBLE as u8]);
        let metal_atoms = vec![metal_atom, inp_ATOM::default(), inp_ATOM::default()];
        let mut metal_valence_atoms = valence_atoms.clone();
        metal_valence_atoms[0].cMetal = 1;
        metal_valence_atoms[0].cInitOrigValenceToMetal = 4;
        metal_valence_atoms[0].cInitValenceToMetal = 2;
        let mode = SRM {
            bMetalAddFlower: 1,
            nMetalMinBondOrder: 0,
            nMetalInitBondOrder: 1,
            nMetalInitEdgeFlow: 0,
            nMetalMaxCharge_D: 4,
            ..SRM::default()
        };
        let (mut st_cap, mut st_flow, mut metal_cap, mut metal_flow) = (-1, -1, -1, -1);
        assert_eq!(
            AtomStcapStflow(
                &metal_atoms,
                &metal_valence_atoms,
                &mode,
                0,
                Some(&mut st_cap),
                Some(&mut st_flow),
                Some(&mut metal_cap),
                Some(&mut metal_flow),
            ),
            Ok(2)
        );
        assert_eq!((st_cap, st_flow, metal_cap, metal_flow), (1, 1, 14, 13));
        assert_eq!(
            AtomStcapStflow(&metal_atoms, &metal_valence_atoms, &mode, 0, None, None, None, None,),
            Ok(2)
        );

        let mut negative_valence = inp_ATOM {
            valence: -1,
            chem_bonds_valence: -7,
            ..inp_ATOM::default()
        };
        negative_valence.neighbor[0] = u16::MAX;
        let mut st_cap = 0;
        assert_eq!(
            AtomStcapStflow(
                &[negative_valence],
                &[crate::source_types::VAL_AT::default()],
                &disabled,
                0,
                Some(&mut st_cap),
                None,
                None,
                None,
            ),
            Ok(0)
        );
        assert_eq!(st_cap, -7);

        let mut reverse_valence = valence_atoms.clone();
        reverse_valence[1].cMetal = 1;
        assert_eq!(
            AtomStcapStflow(&atoms, &reverse_valence, &mode, 0, None, None, None, None,),
            Ok(0)
        );
        assert_eq!(
            AtomStcapStflow(&atoms, &valence_atoms, &mode, -1, None, None, None, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            AtomStcapStflow(&[], &valence_atoms, &mode, 0, None, None, None, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            AtomStcapStflow(&atoms, &[], &mode, 0, None, None, None, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr1__ncountbnssizes__line_1870() {
        let mut hash = 1_469_598_103_934_665_603_u64;
        let mut hash_i32 = |value: i32| {
            for byte in value.to_le_bytes() {
                hash = (hash ^ u64::from(byte)).wrapping_mul(1_099_511_628_211);
            }
        };
        for list in &CN_LIST {
            hash_i32(list.bits);
            for (type_, cap, flow, edges) in list.nodes {
                hash_i32(*type_);
                hash_i32(*cap);
                hash_i32(*flow);
                for (neighbor, edge_cap, forbidden, edge_flow) in edges {
                    hash_i32(*neighbor);
                    hash_i32(*edge_cap);
                    hash_i32(*forbidden);
                    hash_i32(*edge_flow);
                }
            }
        }
        assert_eq!(hash, 6_410_081_296_829_866_238);

        let mut first = inp_ATOM {
            valence: 1,
            el_number: 6,
            ..inp_ATOM::default()
        };
        first.neighbor[0] = 1;
        let mut second = inp_ATOM {
            valence: 1,
            el_number: 8,
            ..inp_ATOM::default()
        };
        second.neighbor[0] = 0;
        let atoms = vec![first, second];
        let va = vec![crate::source_types::VAL_AT::default(); 2];
        let mut heap = SourceHeap::default();
        let mut groups = ALL_TC_GROUPS {
            total_charge: -1,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(
            nCountBnsSizes(
                &mut heap,
                &atoms,
                2,
                17,
                19,
                &T_GROUP_INFO::default(),
                &va,
                &SRM::default(),
                &mut groups,
            ),
            Ok(0)
        );
        assert_eq!((groups.num_atoms, groups.num_bonds), (2, 1));
        assert_eq!((groups.nVertices, groups.nEdges), (2, 1));
        assert_eq!(groups.total_electrons, 15);

        for list_index in 1..=17_i8 {
            let mut heap = SourceHeap::default();
            let mut groups = ALL_TC_GROUPS::default();
            let atom = inp_ATOM {
                el_number: 7,
                ..inp_ATOM::default()
            };
            let va = crate::source_types::VAL_AT {
                cnListIndex: list_index,
                ..crate::source_types::VAL_AT::default()
            };
            assert_eq!(
                nCountBnsSizes(
                    &mut heap,
                    &[atom],
                    1,
                    0,
                    0,
                    &T_GROUP_INFO::default(),
                    &[va],
                    &SRM::default(),
                    &mut groups,
                ),
                Ok(0),
                "cnListIndex={list_index}"
            );
        }

        let mut heap = SourceHeap::default();
        let mut groups = ALL_TC_GROUPS::default();
        assert_eq!(
            nCountBnsSizes(
                &mut heap,
                &[inp_ATOM {
                    el_number: 8,
                    ..inp_ATOM::default()
                }],
                1,
                0,
                0,
                &T_GROUP_INFO::default(),
                &[crate::source_types::VAL_AT {
                    cnListIndex: 11,
                    ..crate::source_types::VAL_AT::default()
                }],
                &SRM::default(),
                &mut groups,
            ),
            Ok(0)
        );
        let positive_group = heap
            .slice(groups.pTCG.as_const())
            .unwrap()
            .iter()
            .find(|group| group.type_ == BNS_VT_C_POS as i32)
            .unwrap();
        assert_eq!(
            (
                positive_group.num_edges,
                positive_group.st_cap,
                positive_group.st_flow,
                positive_group.edges_cap,
                positive_group.edges_flow,
            ),
            (2, 1, 1, 1, 1)
        );

        let mut heap = SourceHeap::default();
        let tgroup_pointer = heap
            .allocate(vec![T_GROUP {
                num: [3, 1, 0, 0, 0],
                nGroupNumber: 1,
                nNumEndpoints: 1,
                ..T_GROUP::default()
            }])
            .unwrap();
        let tautomer_info = T_GROUP_INFO {
            t_group: tgroup_pointer,
            num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        let endpoint_atom = inp_ATOM {
            endpoint: 1,
            ..inp_ATOM::default()
        };
        let mut groups = ALL_TC_GROUPS::default();
        assert_eq!(
            nCountBnsSizes(
                &mut heap,
                &[endpoint_atom],
                1,
                0,
                0,
                &tautomer_info,
                &[crate::source_types::VAL_AT::default()],
                &SRM::default(),
                &mut groups,
            ),
            Ok(0)
        );
        assert_eq!((groups.num_tgroups, groups.num_tgroup_edges), (1, 1));
        assert_eq!(groups.tgroup_charge, -1);
        let registered = heap.slice(groups.pTCG.as_const()).unwrap();
        assert_eq!((registered[0].tg_num_H, registered[0].tg_num_Minus), (2, 1));

        let mut mismatch_groups = ALL_TC_GROUPS::default();
        assert_eq!(
            nCountBnsSizes(
                &mut heap,
                &[inp_ATOM::default()],
                1,
                0,
                0,
                &tautomer_info,
                &[crate::source_types::VAL_AT::default()],
                &SRM::default(),
                &mut mismatch_groups,
            ),
            Ok(BNS_PROGRAM_ERR)
        );

        let mut metal = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        metal.neighbor[0] = 1;
        metal.bond_type[0] = BOND_TYPE_SINGLE as u8;
        let mut ligand = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        ligand.neighbor[0] = 0;
        ligand.bond_type[0] = BOND_TYPE_SINGLE as u8;
        let metal_va = crate::source_types::VAL_AT {
            cMetal: 1,
            cNumBondsToMetal: 1,
            cnListIndex: 18,
            ..crate::source_types::VAL_AT::default()
        };
        let ligand_va = crate::source_types::VAL_AT {
            cNumBondsToMetal: 1,
            ..crate::source_types::VAL_AT::default()
        };
        let mode = SRM {
            bMetalAddFlower: 1,
            nMetalMinBondOrder: 0,
            nMetalInitBondOrder: 1,
            nMetalMaxCharge_D: 4,
            ..SRM::default()
        };
        let mut metal_heap = SourceHeap::default();
        let mut metal_groups = ALL_TC_GROUPS::default();
        assert_eq!(
            nCountBnsSizes(
                &mut metal_heap,
                &[metal, ligand],
                2,
                0,
                0,
                &T_GROUP_INFO::default(),
                &[metal_va, ligand_va],
                &mode,
                &mut metal_groups,
            ),
            Ok(0)
        );
        assert_eq!(metal_groups.num_metal_atoms, 1);
        assert_eq!(metal_groups.num_metal_bonds, 1);
        assert!(metal_groups.num_tc_groups >= 4);

        assert_eq!(
            nCountBnsSizes(
                &mut SourceHeap::default(),
                &[],
                1,
                0,
                0,
                &T_GROUP_INFO::default(),
                &[],
                &SRM::default(),
                &mut ALL_TC_GROUPS::default(),
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr1__naddsupercgroups__line_2279() {
        let mut heap = SourceHeap::default();
        let existing = vec![
            TC_GROUP {
                type_: BNS_VERT_TYPE_TGROUP as i32,
                ord_num: 7,
                ..TC_GROUP::default()
            },
            TC_GROUP {
                type_: BNS_VT_C_POS as i32,
                ..TC_GROUP::default()
            },
            TC_GROUP {
                type_: BNS_VT_C_NEG_C as i32,
                ..TC_GROUP::default()
            },
            TC_GROUP {
                type_: BNS_VT_C_POS_M as i32,
                ..TC_GROUP::default()
            },
            TC_GROUP {
                type_: BNS_VT_M_GROUP as i32,
                ord_num: 0,
                ..TC_GROUP::default()
            },
            TC_GROUP {
                type_: BNS_VT_M_GROUP as i32,
                ord_num: 1,
                ..TC_GROUP::default()
            },
            TC_GROUP {
                type_: BNS_VT_M_GROUP as i32,
                ord_num: 2,
                ..TC_GROUP::default()
            },
            TC_GROUP {
                type_: BNS_VT_M_GROUP as i32,
                ord_num: 3,
                ..TC_GROUP::default()
            },
        ];
        let pointer = heap.allocate(existing).unwrap();
        let mut groups = ALL_TC_GROUPS {
            pTCG: pointer,
            num_tc_groups: 8,
            max_tc_groups: 8,
            nGroup: [-1; 18],
            nVertices: 10,
            nEdges: 20,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(nAddSuperCGroups(&mut heap, &mut groups), Ok(0));
        assert_eq!(groups.num_tc_groups, 10);
        assert_eq!((groups.nVertices, groups.nEdges), (15, 27));
        assert_eq!(groups.nGroup[TCG_Plus0 as usize], 1);
        assert_eq!(groups.nGroup[TCG_Minus_C0 as usize], 2);
        assert_eq!(groups.nGroup[TCG_Plus_M0 as usize], 3);
        assert_eq!(groups.nGroup[TCG_MeFlower0 as usize], 4);
        assert_eq!(groups.nGroup[TCG_MeFlower3 as usize], 7);
        assert_eq!(groups.nGroup[TCG_Plus as usize], 8);
        assert_eq!(groups.nGroup[TCG_Minus as usize], 9);
        let values = heap.slice(groups.pTCG.as_const()).unwrap();
        assert_eq!((values[8].type_, values[8].num_edges), (BNS_VT_C_POS_ALL as i32, 3));
        assert_eq!((values[9].type_, values[9].num_edges), (BNS_VT_C_NEG_ALL as i32, 2));

        let mut empty = ALL_TC_GROUPS {
            nGroup: [-1; 18],
            nVertices: 3,
            nEdges: 4,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(nAddSuperCGroups(&mut SourceHeap::default(), &mut empty), Ok(0));
        assert_eq!((empty.nVertices, empty.nEdges), (3, 4));

        for invalid in [
            vec![TC_GROUP {
                type_: BNS_VT_C_POS as i32,
                ord_num: 1,
                ..TC_GROUP::default()
            }],
            vec![
                TC_GROUP {
                    type_: BNS_VT_C_POS as i32,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: BNS_VT_C_POS as i32,
                    ..TC_GROUP::default()
                },
            ],
            vec![TC_GROUP {
                type_: BNS_VT_M_GROUP as i32,
                ord_num: 4,
                ..TC_GROUP::default()
            }],
        ] {
            let mut heap = SourceHeap::default();
            let len = invalid.len() as i32;
            let pointer = heap.allocate(invalid).unwrap();
            let mut groups = ALL_TC_GROUPS {
                pTCG: pointer,
                num_tc_groups: len,
                max_tc_groups: len,
                nGroup: [-1; 18],
                ..ALL_TC_GROUPS::default()
            };
            assert_eq!(nAddSuperCGroups(&mut heap, &mut groups), Ok(RI_ERR_PROGR));
        }
    }

    #[test]
    fn source_port__ichirvr1__connecttwovertices__line_2627() {
        let mut heap = SourceHeap::default();
        let iedge = heap.allocate(vec![-1_i32; 6]).unwrap();
        let vertices = heap
            .allocate(vec![
                BNS_VERTEX {
                    num_adj_edges: 1,
                    max_adj_edges: 3,
                    iedge,
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    max_adj_edges: 3,
                    iedge: iedge.offset(3).unwrap(),
                    ..BNS_VERTEX::default()
                },
            ])
            .unwrap();
        let edges = heap
            .allocate(vec![
                BNS_EDGE {
                    cap: 9,
                    flow: 7,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE::default(),
            ])
            .unwrap();
        let network = BN_STRUCT {
            vert: vertices,
            edge: edges,
            iedge,
            max_vertices: 2,
            max_edges: 2,
            max_iedges: 6,
            ..BN_STRUCT::default()
        };
        assert_eq!(ConnectTwoVertices(&mut heap, &network, 0, 1, 0, 1), Ok(0));
        let connected = &heap.slice(edges.as_const()).unwrap()[0];
        assert_eq!((connected.neighbor1, connected.neighbor12), (0, 1));
        assert_eq!(connected.neigh_ord, [1, 0]);
        assert_eq!((connected.cap, connected.flow), (0, 0));
        assert_eq!(heap.slice(iedge.as_const()).unwrap(), &[-1, 0, -1, 0, -1, -1]);
        let vertices_after = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            (vertices_after[0].num_adj_edges, vertices_after[1].num_adj_edges),
            (2, 1)
        );

        let edge_before = heap.slice(edges.as_const()).unwrap()[0].clone();
        let vertices_before = heap.slice(vertices.as_const()).unwrap().to_vec();
        assert_eq!(ConnectTwoVertices(&mut heap, &network, 0, 1, 0, 0), Ok(BNS_PROGRAM_ERR));
        assert_eq!(heap.slice(edges.as_const()).unwrap()[0], edge_before);
        assert_eq!(heap.slice(vertices.as_const()).unwrap(), vertices_before);

        assert_eq!(
            ConnectTwoVertices(&mut heap, &network, -1, 1, 1, 1),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!(
            ConnectTwoVertices(&mut heap, &network, 0, 2, 1, 1),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!(
            ConnectTwoVertices(&mut heap, &network, 0, 1, 2, 1),
            Ok(BNS_VERT_EDGE_OVFL)
        );

        let mut self_heap = SourceHeap::default();
        let self_iedge = self_heap.allocate(vec![-1_i32; 2]).unwrap();
        let self_vertices = self_heap
            .allocate(vec![BNS_VERTEX {
                max_adj_edges: 2,
                iedge: self_iedge,
                ..BNS_VERTEX::default()
            }])
            .unwrap();
        let self_edges = self_heap.allocate(vec![BNS_EDGE::default()]).unwrap();
        let self_network = BN_STRUCT {
            vert: self_vertices,
            edge: self_edges,
            iedge: self_iedge,
            max_vertices: 1,
            max_edges: 1,
            max_iedges: 2,
            ..BN_STRUCT::default()
        };
        assert_eq!(ConnectTwoVertices(&mut self_heap, &self_network, 0, 0, 0, 1), Ok(0));
        assert_eq!(self_heap.slice(self_vertices.as_const()).unwrap()[0].num_adj_edges, 2);
        assert_eq!(self_heap.slice(self_edges.as_const()).unwrap()[0].neigh_ord, [1, 0]);

        let bad_vertices = heap
            .allocate(vec![BNS_VERTEX {
                max_adj_edges: 2,
                iedge: iedge.offset(5).unwrap(),
                ..BNS_VERTEX::default()
            }])
            .unwrap();
        let bad_network = BN_STRUCT {
            vert: bad_vertices,
            max_vertices: 1,
            ..network.clone()
        };
        assert_eq!(
            ConnectTwoVertices(&mut heap, &bad_network, 0, 0, 1, 1),
            Ok(BNS_VERT_EDGE_OVFL)
        );
    }

    #[test]
    fn source_port__ichirvr1__addtgroups2tcgbnstruct__line_2425() {
        let mut heap = SourceHeap::default();
        let atom_pointer = heap
            .allocate(vec![
                inp_ATOM {
                    endpoint: 1,
                    ..inp_ATOM::default()
                },
                inp_ATOM::default(),
            ])
            .unwrap();
        let structure = StrFromINChI {
            at: atom_pointer,
            num_atoms: 2,
            ..StrFromINChI::default()
        };
        let incident = heap.allocate(vec![-1_i32; 10]).unwrap();
        let vertices = heap
            .allocate(vec![
                BNS_VERTEX {
                    st_edge: crate::source_types::BNS_ST_EDGE {
                        cap: 5,
                        flow: 1,
                        ..crate::source_types::BNS_ST_EDGE::default()
                    },
                    max_adj_edges: 2,
                    iedge: incident,
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    max_adj_edges: 2,
                    iedge: incident.offset(2).unwrap(),
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX::default(),
                BNS_VERTEX::default(),
            ])
            .unwrap();
        let edges = heap.allocate(vec![BNS_EDGE::default(); 3]).unwrap();
        let group_pointer = heap
            .allocate(vec![TC_GROUP {
                type_: BNS_VERT_TYPE_TGROUP as i32,
                ord_num: 1,
                num_edges: 1,
                st_cap: 3,
                ..TC_GROUP::default()
            }])
            .unwrap();
        let mut groups = ALL_TC_GROUPS {
            pTCG: group_pointer,
            num_tc_groups: 1,
            max_tc_groups: 1,
            num_tgroups: 1,
            num_tgroup_edges: 1,
            ..ALL_TC_GROUPS::default()
        };
        let mut network = BN_STRUCT {
            num_vertices: 2,
            max_vertices: 4,
            max_edges: 3,
            max_iedges: 10,
            vert: vertices,
            edge: edges,
            iedge: incident,
            ..BN_STRUCT::default()
        };
        let mut valence_atoms = vec![crate::source_types::VAL_AT::default(); 2];
        assert_eq!(
            AddTGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 0,),
            Ok(0)
        );
        assert_eq!(
            (network.num_vertices, network.num_edges, network.num_t_groups),
            (3, 1, 1)
        );
        assert_eq!((network.tot_st_cap, network.tot_st_flow), (3, 0));
        assert_eq!(valence_atoms[0].nTautGroupEdge, 1);
        let vertices_after = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            vertices_after[0].type_ & BNS_VERT_TYPE_ENDPOINT as u16,
            BNS_VERT_TYPE_ENDPOINT as u16
        );
        assert_eq!((vertices_after[2].st_edge.cap, vertices_after[2].max_adj_edges), (3, 2));
        assert_eq!(vertices_after[2].iedge, incident.offset(4).unwrap());
        let edge = &heap.slice(edges.as_const()).unwrap()[0];
        assert_eq!((edge.neighbor1, edge.neighbor12, edge.cap, edge.cap0), (0, 2, 2, 2));
        assert_eq!(heap.slice(group_pointer.as_const()).unwrap()[0].nVertexNumber, 2);

        let mut no_groups_network = network.clone();
        let mut no_groups = ALL_TC_GROUPS::default();
        assert_eq!(
            AddTGroups2TCGBnStruct(
                &mut heap,
                &mut no_groups_network,
                &structure,
                &mut valence_atoms,
                &mut no_groups,
                0,
            ),
            Ok(0)
        );

        let mut invalid_heap = SourceHeap::default();
        let invalid_group = invalid_heap
            .allocate(vec![TC_GROUP {
                type_: BNS_VERT_TYPE_TGROUP as i32,
                ord_num: 2,
                ..TC_GROUP::default()
            }])
            .unwrap();
        let mut invalid_groups = ALL_TC_GROUPS {
            pTCG: invalid_group,
            num_tc_groups: 1,
            num_tgroups: 1,
            ..ALL_TC_GROUPS::default()
        };
        let mut invalid_network = BN_STRUCT {
            num_vertices: 1,
            max_vertices: 3,
            max_edges: 2,
            ..BN_STRUCT::default()
        };
        assert_eq!(
            AddTGroups2TCGBnStruct(
                &mut invalid_heap,
                &mut invalid_network,
                &StrFromINChI::default(),
                &mut [],
                &mut invalid_groups,
                0,
            ),
            Ok(BNS_CPOINT_ERR)
        );

        let mut overflow = network.clone();
        overflow.max_vertices = overflow.num_vertices + groups.num_tgroups;
        assert_eq!(
            AddTGroups2TCGBnStruct(&mut heap, &mut overflow, &structure, &mut valence_atoms, &mut groups, 0,),
            Ok(BNS_VERT_EDGE_OVFL)
        );
    }

    #[test]
    fn source_port__ichirvr1__addradicaltometal__line_2706() {
        let mut heap = SourceHeap::default();
        let vertices = heap
            .allocate(vec![
                BNS_VERTEX::default(),
                BNS_VERTEX::default(),
                BNS_VERTEX::default(),
                BNS_VERTEX {
                    st_edge: crate::source_types::BNS_ST_EDGE {
                        cap: i32::MAX,
                        cap0: -1,
                        ..crate::source_types::BNS_ST_EDGE::default()
                    },
                    ..BNS_VERTEX::default()
                },
            ])
            .unwrap();
        let group_pointer = heap
            .allocate(
                (0..4)
                    .map(|index| TC_GROUP {
                        nVertexNumber: index,
                        ..TC_GROUP::default()
                    })
                    .collect(),
            )
            .unwrap();
        let mut slots = [-1; 18];
        slots[TCG_MeFlower0 as usize] = 0;
        slots[TCG_MeFlower1 as usize] = 1;
        slots[TCG_MeFlower2 as usize] = 2;
        slots[TCG_MeFlower3 as usize] = 3;
        let groups = ALL_TC_GROUPS {
            pTCG: group_pointer,
            nGroup: slots,
            num_metal_atoms: 1,
            ..ALL_TC_GROUPS::default()
        };
        let network = BN_STRUCT {
            vert: vertices,
            ..BN_STRUCT::default()
        };
        let mode = SRM {
            bMetalAddFlower: 1,
            ..SRM::default()
        };
        let (mut cap, mut flow) = (5, 17);
        assert_eq!(
            AddRadicalToMetal(&mut heap, &mut cap, &mut flow, &mode, &network, &groups),
            Ok(1)
        );
        assert_eq!((cap, flow), (6, 17));
        let g3 = &heap.slice(vertices.as_const()).unwrap()[3];
        assert_eq!((g3.st_edge.cap, g3.st_edge.cap0), (i32::MIN, 0));

        for (metal_atoms, flower, capacity, missing) in
            [(0, 1, 5, false), (1, 0, 5, false), (1, 1, 4, false), (1, 1, 5, true)]
        {
            let mut groups = groups.clone();
            groups.num_metal_atoms = metal_atoms;
            if missing {
                groups.nGroup[TCG_MeFlower2 as usize] = -1;
            }
            let mode = SRM {
                bMetalAddFlower: flower,
                ..SRM::default()
            };
            let mut cap = capacity;
            assert_eq!(
                AddRadicalToMetal(&mut heap, &mut cap, &mut flow, &mode, &network, &groups),
                Ok(0)
            );
            assert_eq!(cap, capacity);
        }
    }

    #[test]
    fn source_port__ichirvr1__connectmetalflower__line_2743() {
        let fixture = |max_adj_edges: [u16; 4], edge_capacity: i32, edge_flow: i32, group_type: i32| {
            let mut heap = SourceHeap::default();
            let incident = heap.allocate(vec![-1_i32; 11]).unwrap();
            heap.slice_mut(incident).unwrap()[0] = 0;
            let vertices = heap
                .allocate(
                    (0..4)
                        .map(|index| BNS_VERTEX {
                            num_adj_edges: u16::from(index == 0),
                            max_adj_edges: max_adj_edges[index],
                            iedge: incident.offset([0, 3, 6, 9][index]).unwrap(),
                            ..BNS_VERTEX::default()
                        })
                        .collect(),
                )
                .unwrap();
            let edges = heap
                .allocate(vec![
                    BNS_EDGE {
                        cap: edge_capacity,
                        flow: edge_flow,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE::default(),
                    BNS_EDGE::default(),
                    BNS_EDGE::default(),
                    BNS_EDGE::default(),
                    BNS_EDGE::default(),
                ])
                .unwrap();
            let group_pointer = heap
                .allocate(
                    (0..4)
                        .map(|index| TC_GROUP {
                            type_: if index == 0 { group_type } else { 0 },
                            edges_cap: if index == 0 { edge_capacity } else { 0 },
                            edges_flow: if index == 0 { edge_flow } else { 0 },
                            nVertexNumber: index,
                            ..TC_GROUP::default()
                        })
                        .collect(),
                )
                .unwrap();
            let mut slots = [-1; 18];
            slots[TCG_MeFlower0 as usize] = 0;
            slots[TCG_MeFlower1 as usize] = 1;
            slots[TCG_MeFlower2 as usize] = 2;
            slots[TCG_MeFlower3 as usize] = 3;
            let groups = ALL_TC_GROUPS {
                pTCG: group_pointer,
                nGroup: slots,
                ..ALL_TC_GROUPS::default()
            };
            let network = BN_STRUCT {
                vert: vertices,
                edge: edges,
                iedge: incident,
                max_vertices: 4,
                max_edges: 6,
                max_iedges: 11,
                ..BN_STRUCT::default()
            };
            (heap, network, groups, vertices, edges, incident)
        };

        let (mut heap, network, groups, vertices, edges, incident) = fixture([3, 3, 3, 2], 6, 2, BNS_VT_M_GROUP as i32);
        let (mut current_vertices, mut current_edges) = (4, 1);
        let (mut total_capacity, mut total_flow) = (0, 0);
        assert_eq!(
            ConnectMetalFlower(
                &mut heap,
                &mut current_vertices,
                &mut current_edges,
                &mut total_capacity,
                &mut total_flow,
                &SRM {
                    nMetalFlowerParam_D: 2,
                    ..SRM::default()
                },
                &network,
                &groups,
            ),
            Ok(0)
        );
        assert_eq!((current_vertices, current_edges), (4, 6));
        assert_eq!((total_capacity, total_flow), (24, 24));
        let vertex_values = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            vertex_values
                .iter()
                .map(|vertex| (
                    vertex.st_edge.cap,
                    vertex.st_edge.cap0,
                    vertex.st_edge.flow,
                    vertex.st_edge.flow0,
                    vertex.num_adj_edges,
                ))
                .collect::<Vec<_>>(),
            vec![(10, 10, 10, 10, 3), (7, 7, 7, 7, 3), (7, 7, 7, 7, 3), (0, 0, 0, 0, 2)]
        );
        let edge_values = heap.slice(edges.as_const()).unwrap();
        assert_eq!(
            edge_values[1..]
                .iter()
                .map(|edge| (edge.cap, edge.cap0, edge.flow, edge.flow0))
                .collect::<Vec<_>>(),
            vec![(7, 7, 4, 4), (7, 7, 4, 4), (7, 7, 3, 3), (2, 2, 0, 0), (2, 2, 0, 0)]
        );
        assert_eq!(
            heap.slice(incident.as_const()).unwrap(),
            &[0, 2, 1, 2, 3, 5, 1, 3, 4, 5, 4]
        );

        let mut no_groups = ALL_TC_GROUPS::default();
        no_groups.nGroup = [-1; 18];
        let (mut no_group_vertices, mut no_group_edges) = (17, 19);
        let (mut no_group_capacity, mut no_group_flow) = (23, 29);
        assert_eq!(
            ConnectMetalFlower(
                &mut SourceHeap::default(),
                &mut no_group_vertices,
                &mut no_group_edges,
                &mut no_group_capacity,
                &mut no_group_flow,
                &SRM::default(),
                &BN_STRUCT::default(),
                &no_groups,
            ),
            Ok(0)
        );
        assert_eq!(
            (no_group_vertices, no_group_edges, no_group_capacity, no_group_flow),
            (17, 19, 23, 29)
        );
        no_groups.nGroup[TCG_MeFlower0 as usize] = 0;
        assert_eq!(
            ConnectMetalFlower(
                &mut SourceHeap::default(),
                &mut no_group_vertices,
                &mut no_group_edges,
                &mut no_group_capacity,
                &mut no_group_flow,
                &SRM::default(),
                &BN_STRUCT::default(),
                &no_groups,
            ),
            Ok(RI_ERR_PROGR)
        );

        let (mut heap, network, groups, vertices, _, _) = fixture([3, 3, 3, 2], 6, 2, 0);
        let mut values = (4, 1, 0, 0);
        assert_eq!(
            ConnectMetalFlower(
                &mut heap,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &SRM::default(),
                &network,
                &groups,
            ),
            Ok(RI_ERR_PROGR)
        );
        heap.slice_mut(vertices).unwrap()[0].st_edge.cap = 6;
        heap.slice_mut(vertices).unwrap()[0].st_edge.flow = 2;
        assert_eq!(
            ConnectMetalFlower(
                &mut heap,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &SRM::default(),
                &network,
                &groups,
            ),
            Ok(0)
        );

        let (mut heap, network, groups, _, edges, _) = fixture([3, 3, 3, 2], 6, 2, BNS_VT_M_GROUP as i32);
        heap.slice_mut(groups.pTCG).unwrap()[0].edges_cap = 5;
        let mut values = (4, 1, 31, 37);
        assert_eq!(
            ConnectMetalFlower(
                &mut heap,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &SRM::default(),
                &network,
                &groups,
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(heap.slice(edges.as_const()).unwrap()[1..], vec![BNS_EDGE::default(); 5]);

        for (failure_position, maxima) in [[1, 3, 3, 2], [2, 3, 3, 2], [3, 1, 3, 2], [3, 2, 3, 2], [3, 3, 2, 2]]
            .into_iter()
            .enumerate()
        {
            let (mut heap, network, groups, _, edges, _) = fixture(maxima, 6, 2, BNS_VT_M_GROUP as i32);
            let mut values = (4, 1, 41, 43);
            assert_eq!(
                ConnectMetalFlower(
                    &mut heap,
                    &mut values.0,
                    &mut values.1,
                    &mut values.2,
                    &mut values.3,
                    &SRM::default(),
                    &network,
                    &groups,
                ),
                Ok(BNS_VERT_EDGE_OVFL),
                "connection {failure_position}"
            );
            assert_eq!(values, (4, 1, 41, 43));
            assert_eq!(
                heap.slice(edges.as_const()).unwrap()[1..]
                    .iter()
                    .filter(|edge| edge.neighbor1 != 0 || edge.neighbor12 != 0)
                    .count(),
                failure_position
            );
        }

        let (mut heap, network, groups, vertices, edges, _) = fixture([3, 3, 3, 2], 6, 2, BNS_VT_M_GROUP as i32);
        let mut values = (4, 1, 47, 53);
        assert_eq!(
            ConnectMetalFlower(
                &mut heap,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &SRM {
                    nMetalFlowerParam_D: EDGE_FLOW_ST_MASK as i32,
                    ..SRM::default()
                },
                &network,
                &groups,
            ),
            Ok(BNS_PROGRAM_ERR)
        );
        assert_eq!(values, (4, 1, 47, 53));
        assert_eq!(
            heap.slice(vertices.as_const())
                .unwrap()
                .iter()
                .map(|vertex| vertex.num_adj_edges)
                .collect::<Vec<_>>(),
            vec![3, 3, 3, 2]
        );
        assert!(
            heap.slice(edges.as_const()).unwrap()[1..]
                .iter()
                .all(|edge| edge.cap == 0 && edge.flow == 0)
        );

        let (mut heap, network, groups, vertices, edges, _) = fixture([3, 3, 3, 2], -3, -1, BNS_VT_M_GROUP as i32);
        let mut values = (4, 1, 0, 0);
        assert_eq!(
            ConnectMetalFlower(
                &mut heap,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &SRM::default(),
                &network,
                &groups,
            ),
            Ok(0)
        );
        assert_eq!(values, (4, 6, -5, -5));
        assert_eq!(
            heap.slice(vertices.as_const())
                .unwrap()
                .iter()
                .map(|vertex| (vertex.st_edge.cap, vertex.st_edge.flow))
                .collect::<Vec<_>>(),
            vec![(-3, -3), (-1, -1), (-1, -1), (0, 0)]
        );
        assert_eq!(
            heap.slice(edges.as_const()).unwrap()[1..]
                .iter()
                .map(|edge| (edge.cap, edge.flow))
                .collect::<Vec<_>>(),
            vec![(-1, -1), (-2, -1), (-1, 0), (0, 0), (0, 0)]
        );
    }

    #[test]
    fn source_port__ichirvr1__setedgecapflow__line_2911() {
        for (capacity, flow) in [(0, 0), (1, -1), (-17, 23), (i32::MIN, i32::MAX), (i32::MAX, i32::MIN)] {
            let mut edge = BNS_EDGE {
                neighbor1: 7,
                neighbor12: 11,
                cap: 101,
                cap0: 102,
                flow: 103,
                flow0: 104,
                pass: 5,
                forbidden: 6,
                ..BNS_EDGE::default()
            };
            SetEdgeCapFlow(&mut edge, capacity, flow);
            assert_eq!(
                (edge.cap, edge.cap0, edge.flow, edge.flow0),
                (capacity, capacity, flow, flow)
            );
            assert_eq!(
                (edge.neighbor1, edge.neighbor12, edge.pass, edge.forbidden),
                (7, 11, 5, 6)
            );
        }
    }

    #[test]
    fn source_port__ichirvr1__addedgeflow__line_2923() {
        let edge_fixture = || BNS_EDGE {
            cap: 3,
            cap0: 101,
            flow: 5,
            flow0: 102,
            pass: 7,
            forbidden: 9,
            ..BNS_EDGE::default()
        };
        let source_fixture = || BNS_VERTEX {
            st_edge: crate::source_types::BNS_ST_EDGE {
                cap: 7,
                cap0: 103,
                flow: 11,
                flow0: 104,
                pass: 13,
            },
            type_: 17,
            ..BNS_VERTEX::default()
        };
        let destination_fixture = || BNS_VERTEX {
            st_edge: crate::source_types::BNS_ST_EDGE {
                cap: 13,
                cap0: 105,
                flow: 17,
                flow0: 106,
                pass: 19,
            },
            type_: 23,
            ..BNS_VERTEX::default()
        };
        let assert_rejected = |edge_capacity: i32,
                               edge_flow: i32,
                               mut edge: BNS_EDGE,
                               mut source: BNS_VERTEX,
                               mut destination: BNS_VERTEX| {
            let before = (edge.clone(), source.clone(), destination.clone());
            let (mut total_capacity, mut total_flow) = (29, 31);
            assert_eq!(
                AddEdgeFlow(
                    edge_capacity,
                    edge_flow,
                    &mut edge,
                    &mut source,
                    &mut destination,
                    &mut total_capacity,
                    &mut total_flow,
                ),
                BNS_PROGRAM_ERR
            );
            assert_eq!((edge, source, destination), before);
            assert_eq!((total_capacity, total_flow), (29, 31));
        };

        let mut edge = edge_fixture();
        edge.cap = -1;
        assert_rejected(1, 1, edge, source_fixture(), destination_fixture());
        assert_rejected(-1, 1, edge_fixture(), source_fixture(), destination_fixture());
        edge = edge_fixture();
        edge.cap = EDGE_FLOW_MASK as i32 - 1;
        assert_rejected(1, 1, edge, source_fixture(), destination_fixture());

        let mut destination = destination_fixture();
        destination.st_edge.cap = -1;
        assert_rejected(1, 1, edge_fixture(), source_fixture(), destination);
        destination = destination_fixture();
        destination.st_edge.cap = EDGE_FLOW_ST_MASK as i32 - 1;
        assert_rejected(1, 1, edge_fixture(), source_fixture(), destination);
        destination = destination_fixture();
        destination.st_edge.flow = -1;
        assert_rejected(1, 1, edge_fixture(), source_fixture(), destination);
        destination = destination_fixture();
        destination.st_edge.flow = EDGE_FLOW_ST_MASK as i32 - 1;
        assert_rejected(1, 1, edge_fixture(), source_fixture(), destination);

        let mut source = source_fixture();
        source.st_edge.cap = -1;
        assert_rejected(1, 1, edge_fixture(), source, destination_fixture());
        source = source_fixture();
        source.st_edge.flow = -1;
        assert_rejected(1, 1, edge_fixture(), source, destination_fixture());
        source = source_fixture();
        source.st_edge.flow = EDGE_FLOW_ST_MASK as i32 - 1;
        assert_rejected(1, 1, edge_fixture(), source, destination_fixture());

        let (mut edge, mut source, mut destination) = (edge_fixture(), source_fixture(), destination_fixture());
        let (mut total_capacity, mut total_flow) = (37, 41);
        assert_eq!(
            AddEdgeFlow(
                2,
                3,
                &mut edge,
                &mut source,
                &mut destination,
                &mut total_capacity,
                &mut total_flow,
            ),
            0
        );
        assert_eq!((edge.cap, edge.cap0, edge.flow, edge.flow0), (5, 5, 8, 8));
        assert_eq!(
            (
                destination.st_edge.cap,
                destination.st_edge.cap0,
                destination.st_edge.flow,
                destination.st_edge.flow0,
            ),
            (15, 15, 20, 20)
        );
        assert_eq!(
            (
                source.st_edge.cap,
                source.st_edge.cap0,
                source.st_edge.flow,
                source.st_edge.flow0,
            ),
            (7, 103, 14, 14)
        );
        assert_eq!((total_capacity, total_flow), (39, 47));
        assert_eq!((edge.pass, edge.forbidden), (7, 9));
        assert_eq!((source.st_edge.pass, destination.st_edge.pass), (13, 19));

        let (mut edge, mut source, mut destination) = (edge_fixture(), source_fixture(), destination_fixture());
        source.st_edge.flow = 1;
        destination.st_edge.flow = 1;
        let (mut total_capacity, mut total_flow) = (0, 0);
        assert_eq!(
            AddEdgeFlow(
                0,
                -3,
                &mut edge,
                &mut source,
                &mut destination,
                &mut total_capacity,
                &mut total_flow,
            ),
            0
        );
        assert_eq!((edge.flow, destination.st_edge.flow, source.st_edge.flow), (2, -2, -2));
        assert_eq!((total_capacity, total_flow), (0, -6));

        let mut edge = edge_fixture();
        edge.cap = EDGE_FLOW_MASK as i32 - 2;
        let mut source = source_fixture();
        source.st_edge.flow = EDGE_FLOW_ST_MASK as i32 - 2;
        let mut destination = destination_fixture();
        destination.st_edge.cap = EDGE_FLOW_ST_MASK as i32 - 2;
        destination.st_edge.flow = EDGE_FLOW_ST_MASK as i32 - 2;
        let (mut total_capacity, mut total_flow) = (i32::MAX, i32::MAX);
        assert_eq!(
            AddEdgeFlow(
                1,
                1,
                &mut edge,
                &mut source,
                &mut destination,
                &mut total_capacity,
                &mut total_flow,
            ),
            0
        );
        assert_eq!(edge.cap, EDGE_FLOW_MASK as i32 - 1);
        assert_eq!(destination.st_edge.cap, EDGE_FLOW_ST_MASK as i32 - 1);
        assert_eq!(destination.st_edge.flow, EDGE_FLOW_ST_MASK as i32 - 1);
        assert_eq!(source.st_edge.flow, EDGE_FLOW_ST_MASK as i32 - 1);
        assert_eq!((total_capacity, total_flow), (i32::MIN, i32::MIN + 1));
    }

    #[test]
    fn source_port__ichirvr1__connectsupercgroup__line_3058() {
        let v_fixture = |source_max_adj_edges: u16, source_st_cap: i32| {
            let mut heap = SourceHeap::default();
            let incident = heap.allocate_model_storage(vec![-1_i32; 5]).unwrap();
            let vertices = heap
                .allocate_model_storage(vec![
                    BNS_VERTEX {
                        max_adj_edges: source_max_adj_edges,
                        iedge: incident,
                        ..BNS_VERTEX::default()
                    },
                    BNS_VERTEX {
                        st_edge: crate::source_types::BNS_ST_EDGE {
                            cap: 71,
                            cap0: 73,
                            flow: 79,
                            flow0: 83,
                            pass: 89,
                        },
                        type_: 97,
                        iedge: incident.offset(2).unwrap(),
                        ..BNS_VERTEX::default()
                    },
                ])
                .unwrap();
            let edges = heap.allocate_model_storage(vec![BNS_EDGE::default()]).unwrap();
            let group_pointer = heap
                .allocate_model_storage(vec![TC_GROUP {
                    st_cap: source_st_cap,
                    edges_cap: 2,
                    edges_flow: 1,
                    nVertexNumber: 0,
                    nForwardEdge: -1,
                    nBackwardEdge: -1,
                    ..TC_GROUP::default()
                }])
                .unwrap();
            let mut slots = [-1; 18];
            slots[TCG_Minus0 as usize] = 0;
            let groups = ALL_TC_GROUPS {
                pTCG: group_pointer,
                nGroup: slots,
                ..ALL_TC_GROUPS::default()
            };
            let network = BN_STRUCT {
                vert: vertices,
                edge: edges,
                iedge: incident,
                max_vertices: 2,
                max_edges: 1,
                max_iedges: 5,
                ..BN_STRUCT::default()
            };
            (heap, network, groups, vertices, edges, incident)
        };

        let mut no_groups = ALL_TC_GROUPS::default();
        no_groups.nGroup = [-1; 18];
        let mut empty_heap = SourceHeap::default();
        empty_heap.trace_source_allocations();
        let mut empty_values = (11, 13, 17, 19);
        assert_eq!(
            ConnectSuperCGroup(
                &mut empty_heap,
                TCG_Plus as i32,
                &[TCG_Minus0],
                1,
                &mut empty_values.0,
                &mut empty_values.1,
                &mut empty_values.2,
                &mut empty_values.3,
                &BN_STRUCT::default(),
                &mut no_groups,
            ),
            Ok(0)
        );
        assert_eq!(empty_heap.source_allocation_calls(), 0);
        assert_eq!(empty_values, (11, 13, 17, 19));

        let mut same_groups = no_groups.clone();
        same_groups.nGroup[TCG_Plus as usize] = 0;
        assert_eq!(
            ConnectSuperCGroup(
                &mut empty_heap,
                TCG_Plus as i32,
                &[TCG_Plus],
                1,
                &mut empty_values.0,
                &mut empty_values.1,
                &mut empty_values.2,
                &mut empty_values.3,
                &BN_STRUCT::default(),
                &mut same_groups,
            ),
            Ok(0)
        );
        assert_eq!(empty_heap.source_allocation_calls(), 0);

        let (mut heap, network, mut groups, vertices, edges, incident) = v_fixture(2, 4);
        let baseline_allocations = heap.live_allocation_count();
        heap.trace_source_allocations();
        let mut values = (1, 0, 0, 0);
        assert_eq!(
            ConnectSuperCGroup(
                &mut heap,
                -1,
                &[TCG_Minus0],
                1,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &network,
                &mut groups,
            ),
            Ok(1)
        );
        assert_eq!(values, (2, 1, 9, 2));
        assert_eq!(heap.source_allocation_calls(), 4);
        assert_eq!(heap.live_allocation_count(), baseline_allocations);
        let vertex_values = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            (
                vertex_values[0].st_edge.cap,
                vertex_values[0].st_edge.flow,
                vertex_values[1].st_edge.cap,
                vertex_values[1].st_edge.cap0,
                vertex_values[1].st_edge.flow,
                vertex_values[1].st_edge.flow0,
                vertex_values[1].max_adj_edges,
                vertex_values[1].num_adj_edges,
                vertex_values[1].type_,
                vertex_values[1].st_edge.pass,
            ),
            (0, 1, 80, 80, 80, 80, 3, 1, BNS_VT_YVCONNECTOR as u16, 89)
        );
        assert_eq!(heap.slice(edges.as_const()).unwrap()[0].cap, 4);
        assert_eq!(heap.slice(edges.as_const()).unwrap()[0].flow, 1);
        assert_eq!(heap.slice(incident.as_const()).unwrap(), &[0, -1, 0, -1, -1]);
        assert_eq!(heap.slice(groups.pTCG.as_const()).unwrap()[0].nForwardEdge, 0);

        let mut heap = SourceHeap::default();
        let incident = heap.allocate_model_storage(vec![-1_i32; 7]).unwrap();
        let vertices = heap
            .allocate_model_storage(vec![
                BNS_VERTEX {
                    max_adj_edges: 2,
                    iedge: incident,
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    max_adj_edges: 2,
                    iedge: incident.offset(2).unwrap(),
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    iedge: incident.offset(4).unwrap(),
                    ..BNS_VERTEX::default()
                },
            ])
            .unwrap();
        let edges = heap
            .allocate_model_storage(vec![BNS_EDGE::default(), BNS_EDGE::default()])
            .unwrap();
        let group_pointer = heap
            .allocate_model_storage(vec![
                TC_GROUP {
                    nVertexNumber: 0,
                    nForwardEdge: -1,
                    nBackwardEdge: -1,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    st_cap: 4,
                    edges_cap: 2,
                    edges_flow: 1,
                    nVertexNumber: 1,
                    nForwardEdge: -1,
                    nBackwardEdge: -1,
                    ..TC_GROUP::default()
                },
            ])
            .unwrap();
        let mut slots = [-1; 18];
        slots[TCG_Plus as usize] = 0;
        slots[TCG_Minus0 as usize] = 1;
        let mut groups = ALL_TC_GROUPS {
            pTCG: group_pointer,
            nGroup: slots,
            ..ALL_TC_GROUPS::default()
        };
        let network = BN_STRUCT {
            vert: vertices,
            edge: edges,
            iedge: incident,
            max_vertices: 3,
            max_edges: 2,
            max_iedges: 7,
            ..BN_STRUCT::default()
        };
        let baseline_allocations = heap.live_allocation_count();
        let mut values = (2, 0, 0, 0);
        assert_eq!(
            ConnectSuperCGroup(
                &mut heap,
                TCG_Plus as i32,
                &[TCG_Plus, TCG_Minus0],
                2,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &network,
                &mut groups,
            ),
            Ok(1)
        );
        assert_eq!(values, (3, 2, 8, 8));
        assert_eq!(heap.live_allocation_count(), baseline_allocations);
        let vertex_values = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            vertex_values
                .iter()
                .map(|vertex| (vertex.st_edge.cap, vertex.st_edge.flow, vertex.num_adj_edges))
                .collect::<Vec<_>>(),
            vec![(4, 3, 1), (0, 1, 1), (4, 4, 2)]
        );
        assert_eq!(
            heap.slice(edges.as_const())
                .unwrap()
                .iter()
                .map(|edge| (edge.cap, edge.flow))
                .collect::<Vec<_>>(),
            vec![(4, 3), (4, 1)]
        );
        let group_values = heap.slice(group_pointer.as_const()).unwrap();
        assert_eq!(
            (
                group_values[0].nBackwardEdge,
                group_values[0].edges_cap,
                group_values[0].edges_flow,
                group_values[0].st_cap,
                group_values[0].st_flow,
                group_values[1].nForwardEdge,
            ),
            (0, 4, 3, 4, 3, 1)
        );

        for failure_after in 0..4 {
            let (mut heap, network, mut groups, vertices, _, _) = v_fixture(2, 4);
            let baseline_allocations = heap.live_allocation_count();
            let before_vertices = heap.slice(vertices.as_const()).unwrap().to_vec();
            heap.fail_after_allocations(failure_after);
            let mut values = (1, 0, 5, 7);
            assert_eq!(
                ConnectSuperCGroup(
                    &mut heap,
                    -1,
                    &[TCG_Minus0],
                    1,
                    &mut values.0,
                    &mut values.1,
                    &mut values.2,
                    &mut values.3,
                    &network,
                    &mut groups,
                ),
                Ok(RI_ERR_ALLOC),
                "allocation {failure_after}"
            );
            assert_eq!(heap.source_allocation_calls(), 4);
            assert_eq!(heap.live_allocation_count(), baseline_allocations);
            assert_eq!(heap.slice(vertices.as_const()).unwrap(), before_vertices);
            assert_eq!(values, (1, 0, 5, 7));
        }

        let (mut heap, network, mut groups, vertices, _, _) = v_fixture(0, 4);
        let mut values = (1, 0, 11, 13);
        assert_eq!(
            ConnectSuperCGroup(
                &mut heap,
                -1,
                &[TCG_Minus0],
                1,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &network,
                &mut groups,
            ),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!(values, (1, 0, 11, 13));
        let central = &heap.slice(vertices.as_const()).unwrap()[1];
        assert_eq!((central.max_adj_edges, central.type_), (3, BNS_VT_YVCONNECTOR as u16));

        let (mut heap, network, mut groups, vertices, edges, _) = v_fixture(2, -1);
        let mut values = (1, 0, 17, 19);
        assert_eq!(
            ConnectSuperCGroup(
                &mut heap,
                -1,
                &[TCG_Minus0],
                1,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &network,
                &mut groups,
            ),
            Ok(BNS_PROGRAM_ERR)
        );
        assert_eq!(values, (1, 0, 17, 19));
        assert_eq!(heap.slice(vertices.as_const()).unwrap()[1].num_adj_edges, 1);
        assert_eq!(heap.slice(edges.as_const()).unwrap()[0].neighbor12, 1);
        assert_eq!(heap.slice(groups.pTCG.as_const()).unwrap()[0].nForwardEdge, 0);

        let mut heap = SourceHeap::default();
        let incident = heap.allocate_model_storage(vec![-1_i32; 7]).unwrap();
        let vertices = heap
            .allocate_model_storage(vec![
                BNS_VERTEX {
                    st_edge: crate::source_types::BNS_ST_EDGE {
                        cap: -1,
                        ..crate::source_types::BNS_ST_EDGE::default()
                    },
                    max_adj_edges: 2,
                    iedge: incident,
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    max_adj_edges: 2,
                    iedge: incident.offset(2).unwrap(),
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    iedge: incident.offset(4).unwrap(),
                    ..BNS_VERTEX::default()
                },
            ])
            .unwrap();
        let edges = heap
            .allocate_model_storage(vec![BNS_EDGE::default(), BNS_EDGE::default()])
            .unwrap();
        let group_pointer = heap
            .allocate_model_storage(vec![
                TC_GROUP {
                    nVertexNumber: 0,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    st_cap: 4,
                    edges_cap: 2,
                    edges_flow: 1,
                    nVertexNumber: 1,
                    ..TC_GROUP::default()
                },
            ])
            .unwrap();
        let mut slots = [-1; 18];
        slots[TCG_Plus as usize] = 0;
        slots[TCG_Minus0 as usize] = 1;
        let mut groups = ALL_TC_GROUPS {
            pTCG: group_pointer,
            nGroup: slots,
            ..ALL_TC_GROUPS::default()
        };
        let network = BN_STRUCT {
            vert: vertices,
            edge: edges,
            iedge: incident,
            max_vertices: 3,
            max_edges: 2,
            max_iedges: 7,
            ..BN_STRUCT::default()
        };
        let mut values = (2, 0, 0, 0);
        assert_eq!(
            ConnectSuperCGroup(
                &mut heap,
                TCG_Plus as i32,
                &[TCG_Minus0],
                1,
                &mut values.0,
                &mut values.1,
                &mut values.2,
                &mut values.3,
                &network,
                &mut groups,
            ),
            Ok(BNS_PROGRAM_ERR)
        );
        assert_eq!(values, (2, 0, 4, 2));
        assert_eq!(heap.slice(vertices.as_const()).unwrap()[2].st_edge.cap, 4);
        assert_eq!(heap.slice(vertices.as_const()).unwrap()[2].st_edge.flow, 1);
        assert_eq!(heap.slice(group_pointer.as_const()).unwrap()[0].edges_cap, 0);
    }

    #[test]
    fn source_port__ichirvr1__addstcapflow__line_3218() {
        for (old_cap, old_flow, add_cap, add_flow, total_cap, total_flow) in [
            (3, 2, 7, 5, 11, 13),
            (7, 5, -3, -2, -11, -13),
            (19, 23, 0, 0, 29, 31),
            (i32::MAX, i32::MIN, 1, -1, i32::MAX, i32::MIN),
            (i32::MIN, i32::MAX, -1, 1, i32::MIN, i32::MAX),
        ] {
            let mut vertex = BNS_VERTEX {
                st_edge: crate::source_types::BNS_ST_EDGE {
                    cap: old_cap,
                    cap0: 101,
                    flow: old_flow,
                    flow0: 103,
                    pass: 107,
                },
                type_: 109,
                num_adj_edges: 2,
                max_adj_edges: 3,
                ..BNS_VERTEX::default()
            };
            let (mut capacity_total, mut flow_total) = (total_cap, total_flow);
            AddStCapFlow(&mut vertex, &mut flow_total, &mut capacity_total, add_cap, add_flow);
            let expected_cap = old_cap.wrapping_add(add_cap);
            let expected_flow = old_flow.wrapping_add(add_flow);
            assert_eq!(
                (
                    vertex.st_edge.cap,
                    vertex.st_edge.cap0,
                    vertex.st_edge.flow,
                    vertex.st_edge.flow0,
                ),
                (expected_cap, expected_cap, expected_flow, expected_flow)
            );
            assert_eq!(
                (capacity_total, flow_total),
                (total_cap.wrapping_add(add_cap), total_flow.wrapping_add(add_flow),)
            );
            assert_eq!(
                (
                    vertex.st_edge.pass,
                    vertex.type_,
                    vertex.num_adj_edges,
                    vertex.max_adj_edges,
                ),
                (107, 109, 2, 3)
            );
        }
    }

    #[test]
    fn source_port__ichirvr1__setstcapflow__line_3231() {
        for (old_cap, old_flow, new_cap, new_flow, total_cap, total_flow) in [
            (3, 2, 7, 5, 11, 13),
            (7, 5, 3, 2, -11, -13),
            (i32::MAX, i32::MIN, i32::MIN, i32::MAX, i32::MAX, i32::MIN),
        ] {
            let mut vertex = BNS_VERTEX {
                st_edge: crate::source_types::BNS_ST_EDGE {
                    cap: old_cap,
                    cap0: 101,
                    flow: old_flow,
                    flow0: 102,
                    pass: 7,
                },
                type_: 19,
                ..BNS_VERTEX::default()
            };
            let (mut cap_total, mut flow_total) = (total_cap, total_flow);
            let expected_cap_total = total_cap.wrapping_add(new_cap.wrapping_sub(old_cap));
            let expected_flow_total = total_flow.wrapping_add(new_flow.wrapping_sub(old_flow));
            SetStCapFlow(&mut vertex, &mut flow_total, &mut cap_total, new_cap, new_flow);
            assert_eq!((cap_total, flow_total), (expected_cap_total, expected_flow_total));
            assert_eq!((vertex.st_edge.cap, vertex.st_edge.cap0), (new_cap, new_cap));
            assert_eq!((vertex.st_edge.flow, vertex.st_edge.flow0), (new_flow, new_flow));
            assert_eq!((vertex.st_edge.pass, vertex.type_), (7, 19));
        }
    }

    #[test]
    fn source_port__ichirvr1__ntautendpointedgecap__line_1574() {
        let mut atom = inp_ATOM {
            valence: 2,
            chem_bonds_valence: 4,
            ..inp_ATOM::default()
        };
        let mut valence_atom = crate::source_types::VAL_AT {
            cInitFreeValences: -3,
            ..crate::source_types::VAL_AT::default()
        };
        assert_eq!(
            nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0),
            Ok(-1)
        );

        let expected_deltas = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
        for (index, delta) in expected_deltas.into_iter().enumerate() {
            valence_atom.cnListIndex = (index + 1) as i8;
            assert_eq!(
                nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0),
                Ok(-1 + delta),
                "cnList index {}",
                index + 1
            );
        }

        valence_atom.cnListIndex = 0;
        valence_atom.cInitFreeValences = 5;
        atom.sb_parity = [1, 1, 1];
        atom.sb_ord = [0, 1, 2];
        atom.bond_type[..3].copy_from_slice(&[0, BOND_TYPE_DOUBLE as u8, BOND_TYPE_TRIPLE as u8]);
        assert_eq!(nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0), Ok(7));

        atom.bond_type[0] = BOND_TYPE_SINGLE as u8;
        atom.bond_type[1] = BOND_TYPE_TRIPLE as u8 + 1;
        assert_eq!(nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0), Ok(7));

        atom.sb_parity = [1, 0, 1];
        atom.bond_type[0] = BOND_TYPE_DOUBLE as u8;
        atom.sb_ord[2] = -1;
        assert_eq!(nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0), Ok(6));

        atom.sb_parity = [0, 1, 1];
        atom.sb_ord = [-1, -1, -1];
        assert_eq!(nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0), Ok(7));

        atom.sb_parity = [1, 0, 0];
        assert_eq!(
            nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        atom.sb_parity = [0, 0, 0];
        atom.valence = 5;
        atom.chem_bonds_valence = 4;
        assert_eq!(
            nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0),
            Ok(RI_ERR_PROGR)
        );

        valence_atom.cnListIndex = 19;
        assert_eq!(
            nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            nTautEndpointEdgeCap(&[atom.clone()], &[valence_atom.clone()], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            nTautEndpointEdgeCap(&[], &[valence_atom], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            nTautEndpointEdgeCap(&[atom], &[], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr1__registertcgroup__line_1527() {
        let mut heap = SourceHeap::default();
        let mut groups = ALL_TC_GROUPS::default();
        assert_eq!(
            RegisterTCGroup(&mut heap, &mut groups, 7, 11, i32::MAX, i32::MIN, -3, 5, 2,),
            Ok(1)
        );
        assert_eq!(groups.num_tc_groups, 1);
        assert_eq!(groups.max_tc_groups, INC_NUM_TCGROUPS as i32);
        let first = &heap.slice(groups.pTCG.as_const()).unwrap()[0];
        assert_eq!((first.type_, first.ord_num), (7, 11));
        assert_eq!(first.num_edges, 2);
        assert_eq!(first.st_cap, i32::MAX);
        assert_eq!(first.st_flow, i32::MIN);
        assert_eq!(first.edges_cap, -3);
        assert_eq!(first.edges_flow, 5);

        assert_eq!(
            RegisterTCGroup(&mut heap, &mut groups, 7, 11, 1, -1, 4, -6, i32::MAX),
            Ok(0)
        );
        let first = &heap.slice(groups.pTCG.as_const()).unwrap()[0];
        assert_eq!(first.num_edges, i32::MIN + 1);
        assert_eq!(first.st_cap, i32::MIN);
        assert_eq!(first.st_flow, i32::MAX);
        assert_eq!(first.edges_cap, 1);
        assert_eq!(first.edges_flow, -1);

        assert_eq!(RegisterTCGroup(&mut heap, &mut groups, 7, 12, 2, 3, 4, 5, 6), Ok(2));
        assert_eq!(groups.num_tc_groups, 2);
        let second = &heap.slice(groups.pTCG.as_const()).unwrap()[1];
        assert_eq!((second.type_, second.ord_num), (7, 12));
        assert_eq!(
            (
                second.num_edges,
                second.st_cap,
                second.st_flow,
                second.edges_cap,
                second.edges_flow,
            ),
            (6, 2, 3, 4, 5)
        );

        let duplicate_pointer = heap
            .allocate_model_storage(vec![
                TC_GROUP {
                    type_: 3,
                    ord_num: 4,
                    num_edges: 10,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: 3,
                    ord_num: 4,
                    num_edges: 20,
                    ..TC_GROUP::default()
                },
            ])
            .unwrap();
        let mut duplicates = ALL_TC_GROUPS {
            pTCG: duplicate_pointer,
            num_tc_groups: 2,
            max_tc_groups: 2,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(RegisterTCGroup(&mut heap, &mut duplicates, 3, 4, 0, 0, 0, 0, 1), Ok(0));
        let duplicate_values = heap.slice(duplicate_pointer.as_const()).unwrap();
        assert_eq!(duplicate_values[0].num_edges, 11);
        assert_eq!(duplicate_values[1].num_edges, 20);

        let full_pointer = heap
            .allocate_model_storage(vec![TC_GROUP {
                type_: 1,
                ord_num: 1,
                ..TC_GROUP::default()
            }])
            .unwrap();
        let mut full = ALL_TC_GROUPS {
            pTCG: full_pointer,
            num_tc_groups: 1,
            max_tc_groups: 1,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(RegisterTCGroup(&mut heap, &mut full, 2, 2, 1, 2, 3, 4, 5), Ok(2));
        assert_eq!(full.max_tc_groups, 1 + INC_NUM_TCGROUPS as i32);
        assert_eq!(full.num_tc_groups, 2);
        assert_eq!(
            heap.slice(full_pointer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut failure_heap = SourceHeap::default();
        let retained = failure_heap
            .allocate_model_storage(vec![TC_GROUP {
                type_: 1,
                ord_num: 1,
                ..TC_GROUP::default()
            }])
            .unwrap();
        let mut failure = ALL_TC_GROUPS {
            pTCG: retained,
            num_tc_groups: 1,
            max_tc_groups: 1,
            ..ALL_TC_GROUPS::default()
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            RegisterTCGroup(&mut failure_heap, &mut failure, 2, 2, 1, 2, 3, 4, 5,),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(failure.pTCG, retained);
        assert_eq!(failure.num_tc_groups, 1);
        assert_eq!(failure.max_tc_groups, 1);
        assert_eq!(
            failure_heap.slice(retained.as_const()).unwrap()[0],
            TC_GROUP {
                type_: 1,
                ord_num: 1,
                ..TC_GROUP::default()
            }
        );
    }

    #[test]
    fn source_port__ichirvr1__realloctcgroups__line_1502() {
        let mut heap = SourceHeap::default();
        let old = heap
            .allocate_model_storage(vec![
                TC_GROUP {
                    type_: 11,
                    ord_num: 12,
                    num_edges: 13,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: 21,
                    ord_num: 22,
                    num_edges: 23,
                    ..TC_GROUP::default()
                },
            ])
            .unwrap();
        let mut groups = ALL_TC_GROUPS {
            pTCG: old,
            num_tc_groups: 2,
            max_tc_groups: 2,
            nVertices: 77,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(ReallocTCGroups(&mut heap, &mut groups, 3), Ok(0));
        assert_ne!(groups.pTCG, old);
        assert_eq!(groups.num_tc_groups, 2);
        assert_eq!(groups.max_tc_groups, 5);
        assert_eq!(groups.nVertices, 77);
        let values = heap.slice(groups.pTCG.as_const()).unwrap();
        assert_eq!(values.len(), 5);
        assert_eq!((values[0].type_, values[0].ord_num, values[0].num_edges), (11, 12, 13));
        assert_eq!((values[1].type_, values[1].ord_num, values[1].num_edges), (21, 22, 23));
        assert_eq!(&values[2..], vec![TC_GROUP::default(); 3]);
        assert_eq!(heap.slice(old.as_const()), Err(SourceHeapError::MissingAllocation));

        let mut null_groups = ALL_TC_GROUPS::default();
        assert_eq!(ReallocTCGroups(&mut heap, &mut null_groups, 2), Ok(0));
        assert_eq!(null_groups.max_tc_groups, 2);
        assert_eq!(
            heap.slice(null_groups.pTCG.as_const()).unwrap(),
            &[TC_GROUP::default(), TC_GROUP::default()]
        );

        let empty_old = heap
            .allocate_model_storage(vec![TC_GROUP {
                type_: 91,
                ..TC_GROUP::default()
            }])
            .unwrap();
        let mut empty_copy = ALL_TC_GROUPS {
            pTCG: empty_old,
            num_tc_groups: 0,
            max_tc_groups: 1,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(ReallocTCGroups(&mut heap, &mut empty_copy, 1), Ok(0));
        assert_eq!(
            heap.slice(empty_copy.pTCG.as_const()).unwrap(),
            &[TC_GROUP::default(), TC_GROUP::default()]
        );
        assert_eq!(
            heap.slice(empty_old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let zero_old = empty_copy.pTCG;
        assert_eq!(ReallocTCGroups(&mut heap, &mut empty_copy, 0), Ok(0));
        assert_eq!(empty_copy.max_tc_groups, 2);
        assert_ne!(empty_copy.pTCG, zero_old);
        assert_eq!(heap.slice(zero_old.as_const()), Err(SourceHeapError::MissingAllocation));

        let mut failure_heap = SourceHeap::default();
        let retained = failure_heap
            .allocate_model_storage(vec![TC_GROUP {
                type_: 31,
                ..TC_GROUP::default()
            }])
            .unwrap();
        let mut retained_groups = ALL_TC_GROUPS {
            pTCG: retained,
            num_tc_groups: 1,
            max_tc_groups: 1,
            ..ALL_TC_GROUPS::default()
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            ReallocTCGroups(&mut failure_heap, &mut retained_groups, 4),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(retained_groups.pTCG, retained);
        assert_eq!(retained_groups.max_tc_groups, 1);
        assert_eq!(failure_heap.slice(retained.as_const()).unwrap()[0].type_, 31);

        assert_eq!(
            ReallocTCGroups(
                &mut heap,
                &mut ALL_TC_GROUPS {
                    max_tc_groups: 0,
                    ..ALL_TC_GROUPS::default()
                },
                -1,
            ),
            Ok(RI_ERR_ALLOC)
        );
    }

    #[test]
    fn source_port__ichirvr1__get_sp_element_type__line_1414() {
        let cases = [
            (i32::MIN, i32::MAX, 1),
            (-1, -2, 1),
            (0, -1, 1),
            (1, 1, 0),
            (2, 0, 0),
            (3, 2, 1),
            (9, 8, 1),
            (10, 0, 1),
            (11, 2, 2),
            (17, 8, 2),
            (18, 0, 2),
            (19, 2, 3),
            (20, 3, 3),
            (21, 0, 3),
            (30, 0, 3),
            (31, 4, 3),
            (35, 8, 3),
            (36, 0, 3),
            (37, 2, 4),
            (38, 3, 4),
            (39, 0, 4),
            (48, 0, 4),
            (49, 4, 4),
            (53, 8, 4),
            (54, 0, 4),
            (55, 2, 5),
            (56, 3, 5),
            (57, 0, 5),
            (80, 0, 5),
            (81, 4, 5),
            (85, 8, 5),
            (86, 0, 5),
            (87, 2, 6),
            (88, 3, 6),
            (89, 0, 6),
            (i32::MAX, 0, 6),
        ];
        for (periodic_number, expected_type, expected_row) in cases {
            let mut row = -77;
            assert_eq!(
                get_sp_element_type(periodic_number, &mut row),
                expected_type,
                "periodic number {periodic_number}"
            );
            assert_eq!(row, expected_row, "periodic number {periodic_number}");
        }
        let zero_types = [
            10, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 36, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 54, 57, 58, 59,
            60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 86,
        ];
        for periodic_number in 3..=88 {
            let mut row = -1;
            let element_type = get_sp_element_type(periodic_number, &mut row);
            assert!((0..=8).contains(&element_type));
            assert!((1..=6).contains(&row));
            assert_eq!(element_type == 0, zero_types.contains(&periodic_number));
        }
    }

    #[test]
    fn source_port__ichirvr1__clear_t_group_info__line_509() {
        let mut heap = SourceHeap::default();
        assert_eq!(clear_t_group_info(&mut heap, None), Ok(()));

        let groups = heap
            .allocate_model_storage(vec![
                T_GROUP {
                    nGroupNumber: 11,
                    nNumEndpoints: 12,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    nGroupNumber: 21,
                    nNumEndpoints: 22,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    nGroupNumber: 31,
                    nNumEndpoints: 32,
                    ..T_GROUP::default()
                },
            ])
            .unwrap();
        let group_numbers = heap.allocate_model_storage(vec![41_u16, 42, 43]).unwrap();
        let endpoints = heap.allocate_model_storage(vec![51_u16, 52, 53]).unwrap();
        let isotope_endpoints = heap.allocate_model_storage(vec![61_u16, 62, 63, 64]).unwrap();
        let mut info = T_GROUP_INFO {
            t_group: groups,
            nEndpointAtomNumber: endpoints,
            tGroupNumber: group_numbers,
            nNumEndpoints: 1,
            num_t_groups: 2,
            max_num_t_groups: 2,
            bIgnoreIsotopic: 71,
            nIsotopicEndpointAtomNumber: isotope_endpoints,
            nNumIsotopicEndpoints: 3,
            num_iso_H: [72, 73, 74],
            bTautFlags: 75,
            bTautFlagsDone: 76,
            ..T_GROUP_INFO::default()
        };
        info.tni.nNumRemovedExplicitH = 77;
        assert_eq!(clear_t_group_info(&mut heap, Some(&mut info)), Ok(()));
        assert_eq!(info.t_group, groups);
        assert_eq!(info.max_num_t_groups, 2);
        assert_eq!(info.tGroupNumber, group_numbers);
        assert_eq!(info.num_t_groups, 2);
        assert_eq!(info.nEndpointAtomNumber, endpoints);
        assert_eq!(info.nNumEndpoints, 1);
        assert_eq!(info.nIsotopicEndpointAtomNumber, isotope_endpoints);
        assert_eq!(info.nNumIsotopicEndpoints, 3);
        assert_eq!(info.bIgnoreIsotopic, 0);
        assert_eq!(info.num_iso_H, [0; 3]);
        assert_eq!(info.tni, Default::default());
        assert_eq!(info.bTautFlags, 0);
        assert_eq!(info.bTautFlagsDone, 0);
        assert_eq!(
            heap.slice(groups.as_const()).unwrap(),
            &[
                T_GROUP::default(),
                T_GROUP::default(),
                T_GROUP {
                    nGroupNumber: 31,
                    nNumEndpoints: 32,
                    ..T_GROUP::default()
                }
            ]
        );
        assert_eq!(heap.slice(group_numbers.as_const()).unwrap(), &[0, 0, 43]);
        assert_eq!(heap.slice(endpoints.as_const()).unwrap(), &[0, 52, 53]);
        assert_eq!(heap.slice(isotope_endpoints.as_const()).unwrap(), &[0, 0, 0, 64]);

        let mut all_null = T_GROUP_INFO {
            max_num_t_groups: 81,
            num_t_groups: 82,
            nNumEndpoints: 83,
            nNumIsotopicEndpoints: 84,
            bIgnoreIsotopic: 85,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(clear_t_group_info(&mut heap, Some(&mut all_null)), Ok(()));
        assert_eq!(all_null, T_GROUP_INFO::default());

        let zero_groups = heap
            .allocate_model_storage(vec![T_GROUP {
                nGroupNumber: 91,
                ..T_GROUP::default()
            }])
            .unwrap();
        let zero_group_numbers = heap.allocate_model_storage(vec![92_u16]).unwrap();
        let zero_endpoints = heap.allocate_model_storage(vec![93_u16]).unwrap();
        let zero_isotope_endpoints = heap.allocate_model_storage(vec![94_u16]).unwrap();
        let mut zero_lengths = T_GROUP_INFO {
            t_group: zero_groups,
            tGroupNumber: zero_group_numbers,
            nEndpointAtomNumber: zero_endpoints,
            nIsotopicEndpointAtomNumber: zero_isotope_endpoints,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(clear_t_group_info(&mut heap, Some(&mut zero_lengths)), Ok(()));
        assert_eq!(zero_lengths.t_group, zero_groups);
        assert_eq!(zero_lengths.tGroupNumber, zero_group_numbers);
        assert_eq!(zero_lengths.nEndpointAtomNumber, zero_endpoints);
        assert_eq!(zero_lengths.nIsotopicEndpointAtomNumber, zero_isotope_endpoints);
        assert_eq!(heap.slice(zero_groups.as_const()).unwrap()[0].nGroupNumber, 91);
        assert_eq!(heap.slice(zero_group_numbers.as_const()).unwrap(), &[92]);
        assert_eq!(heap.slice(zero_endpoints.as_const()).unwrap(), &[93]);
        assert_eq!(heap.slice(zero_isotope_endpoints.as_const()).unwrap(), &[94]);
    }

    #[test]
    fn source_port__ichirvr1__gettgroupinfofrominchi__line_575() {
        let mut heap = SourceHeap::default();
        let retained_groups = heap
            .allocate_model_storage(vec![T_GROUP {
                nGroupNumber: 99,
                ..T_GROUP::default()
            }])
            .unwrap();
        let mut inactive_info = T_GROUP_INFO {
            t_group: retained_groups,
            max_num_t_groups: 1,
            bIgnoreIsotopic: 7,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            GetTgroupInfoFromInChI(
                &mut heap,
                &mut inactive_info,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                None,
            ),
            Ok(0)
        );
        assert_eq!(inactive_info.t_group, retained_groups);
        assert_eq!(inactive_info.max_num_t_groups, 1);
        assert_eq!(inactive_info.bIgnoreIsotopic, 0);
        assert_eq!(heap.slice(retained_groups.as_const()).unwrap(), &[T_GROUP::default()]);

        let tautomer = heap
            .allocate_model_storage(vec![2_u16, 3, 1, 0, 1, 4, 2, 1, 2, 3])
            .unwrap();
        let atoms = heap.allocate_model_storage(vec![inp_ATOM::default(); 3]).unwrap();
        let endpoints = heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let inchi = INChI {
            nNumberOfAtoms: 3,
            lenTautomer: 10,
            nTautomer: tautomer,
            ..INChI::default()
        };
        let mut info = T_GROUP_INFO::default();
        assert_eq!(
            GetTgroupInfoFromInChI(&mut heap, &mut info, atoms, endpoints, Some(&inchi)),
            Ok(0)
        );
        assert_eq!(info.max_num_t_groups, 2);
        assert_eq!(info.num_t_groups, 2);
        assert_eq!(info.nNumEndpoints, 3);
        let groups = heap.slice(info.t_group.as_const()).unwrap();
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].num[..2], [1, 0]);
        assert_eq!(groups[0].nGroupNumber, 1);
        assert_eq!(groups[0].nNumEndpoints, 1);
        assert_eq!(groups[0].nFirstEndpointAtNoPos, 0);
        assert_eq!(groups[1].num[..2], [3, 1]);
        assert_eq!(groups[1].nGroupNumber, 2);
        assert_eq!(groups[1].nNumEndpoints, 2);
        assert_eq!(groups[1].nFirstEndpointAtNoPos, 1);
        assert_eq!(groups[2], T_GROUP::default());
        let numbers = heap.slice(info.tGroupNumber.as_const()).unwrap();
        assert_eq!(numbers.len(), 12);
        assert_eq!(&numbers[..2], &[0, 1]);
        assert_eq!(&numbers[4..6], &[0, 1]);
        assert_eq!(heap.slice(info.nEndpointAtomNumber.as_const()).unwrap(), &[0, 1, 2, 0]);
        assert_eq!(
            heap.slice(atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.endpoint)
                .collect::<Vec<_>>(),
            vec![1, 2, 2]
        );
        assert_eq!(heap.slice(endpoints.as_const()).unwrap(), &[1, 2, 2]);

        heap.trace_source_allocations();
        assert_eq!(
            GetTgroupInfoFromInChI(&mut heap, &mut info, atoms, endpoints, Some(&inchi)),
            Ok(0)
        );
        assert_eq!(heap.source_allocation_calls(), 0);

        let invalid_tautomer = heap.allocate_model_storage(vec![1_u16, 3, 1, 0, 0]).unwrap();
        let invalid_inchi = INChI {
            nNumberOfAtoms: 1,
            lenTautomer: 5,
            nTautomer: invalid_tautomer,
            ..INChI::default()
        };
        let mut invalid_info = T_GROUP_INFO::default();
        assert_eq!(
            GetTgroupInfoFromInChI(
                &mut heap,
                &mut invalid_info,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                Some(&invalid_inchi),
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(
            heap.slice(invalid_info.nEndpointAtomNumber.as_const()).unwrap()[0],
            u16::MAX
        );

        let boundary_tautomer = heap.allocate_model_storage(vec![1_u16, 3, 0, 0, 2]).unwrap();
        let boundary_inchi = INChI {
            nNumberOfAtoms: 1,
            lenTautomer: 5,
            nTautomer: boundary_tautomer,
            ..INChI::default()
        };
        let mut boundary_info = T_GROUP_INFO::default();
        assert_eq!(
            GetTgroupInfoFromInChI(
                &mut heap,
                &mut boundary_info,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                Some(&boundary_inchi),
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(boundary_info.nEndpointAtomNumber.as_const()).unwrap()[0], 1);

        let mismatch_tautomer = heap.allocate_model_storage(vec![1_u16, 3, 0, 0, 1]).unwrap();
        let mismatch_inchi = INChI {
            nNumberOfAtoms: 1,
            lenTautomer: 6,
            nTautomer: mismatch_tautomer,
            ..INChI::default()
        };
        let mut mismatch_info = T_GROUP_INFO::default();
        assert_eq!(
            GetTgroupInfoFromInChI(
                &mut heap,
                &mut mismatch_info,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                Some(&mismatch_inchi),
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(mismatch_info.nNumEndpoints, 2);

        for failure_after in 0..3 {
            let mut failure_heap = SourceHeap::default();
            let failure_tautomer = failure_heap.allocate_model_storage(vec![1_u16, 3, 0, 0, 1]).unwrap();
            let failure_inchi = INChI {
                nNumberOfAtoms: 1,
                lenTautomer: 5,
                nTautomer: failure_tautomer,
                ..INChI::default()
            };
            failure_heap.fail_after_allocations(failure_after);
            let mut failure_info = T_GROUP_INFO::default();
            assert_eq!(
                GetTgroupInfoFromInChI(
                    &mut failure_heap,
                    &mut failure_info,
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    Some(&failure_inchi),
                ),
                Ok(RI_ERR_ALLOC)
            );
            assert_eq!(failure_heap.source_allocation_calls(), 3);
            assert_eq!(failure_info.t_group.is_null(), failure_after == 0);
            assert_eq!(failure_info.tGroupNumber.is_null(), failure_after == 1);
            assert_eq!(failure_info.nEndpointAtomNumber.is_null(), failure_after == 2);
        }
    }

    #[test]
    fn source_port__ichirvr1__makeinchioutofstrfrominchi2__line_5669() {
        let mut allocation_heap = SourceHeap::default();
        allocation_heap.fail_after_allocations(0);
        assert_eq!(
            MakeInChIOutOfStrFromINChI2(
                &mut allocation_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                None,
                None,
                None,
                i32::MIN,
                i32::MAX,
                i64::MIN,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(allocation_heap.source_allocation_calls(), 1);
        assert_eq!(allocation_heap.live_allocation_count(), 0);

        let mut null_heap = SourceHeap::default();
        assert_eq!(
            MakeInChIOutOfStrFromINChI2(
                &mut null_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                None,
                None,
                None,
                0,
                0,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(null_heap.source_allocation_calls(), 1);
        assert_eq!(null_heap.live_allocation_count(), 0);

        let mut atom_allocation_heap = SourceHeap::default();
        let empty_atoms = atom_allocation_heap
            .allocate_model_storage(Vec::<inp_ATOM>::new())
            .unwrap();
        let mut atom_allocation_structure = StrFromINChI {
            at2: empty_atoms,
            ..StrFromINChI::default()
        };
        atom_allocation_heap.fail_after_allocations(1);
        assert_eq!(
            MakeInChIOutOfStrFromINChI2(
                &mut atom_allocation_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                Some(&INPUT_PARMS::default()),
                Some(&mut STRUCT_DATA::default()),
                Some(&mut atom_allocation_structure),
                0,
                0,
                0,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(atom_allocation_heap.source_allocation_calls(), 2);
        assert_eq!(atom_allocation_heap.live_allocation_count(), 1);

        let mut connect_heap = SourceHeap::default();
        let mut deleted_h = inp_ATOM::default();
        deleted_h.neighbor[0] = 0;
        let connect_atoms = connect_heap
            .allocate_model_storage(vec![inp_ATOM::default(), deleted_h])
            .unwrap();
        let connect_clock = connect_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let connect_canonical = connect_heap
            .allocate_model_storage(vec![CANON_GLOBALS::default()])
            .unwrap();
        let mut connect_structure = StrFromINChI {
            at2: connect_atoms,
            num_atoms: 1,
            num_deleted_H: 1,
            ..StrFromINChI::default()
        };
        connect_structure.RevInChI.nRetVal = 91;
        assert_eq!(
            MakeInChIOutOfStrFromINChI2(
                &mut connect_heap,
                connect_clock,
                connect_canonical,
                Some(&INPUT_PARMS::default()),
                Some(&mut STRUCT_DATA::default()),
                Some(&mut connect_structure),
                2,
                17,
                -3,
                0,
            ),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(
            connect_structure.RevInChI.pINChI,
            [SourceMutPointer::null(); INCHI_NUM as usize]
        );
        assert_eq!(
            connect_structure.RevInChI.pINChI_Aux,
            [SourceMutPointer::null(); INCHI_NUM as usize]
        );
        assert_eq!(connect_structure.RevInChI.nRetVal, 91);

        let mut success_heap = SourceHeap::default();
        let success_atoms = success_heap.allocate_model_storage(Vec::<inp_ATOM>::new()).unwrap();
        let success_clock = success_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let success_canonical = success_heap
            .allocate_model_storage(vec![CANON_GLOBALS::default()])
            .unwrap();
        let mut success_input_data = STRUCT_DATA::default();
        success_input_data.pStrErrStruct[..6]
            .copy_from_slice(&[b'p' as i8, b'r' as i8, b'i' as i8, b'o' as i8, b'r' as i8, 0]);
        let mut success_structure = StrFromINChI {
            at2: success_atoms,
            nChargeRevrs: 37,
            ..StrFromINChI::default()
        };
        assert_eq!(
            MakeInChIOutOfStrFromINChI2(
                &mut success_heap,
                success_clock,
                success_canonical,
                Some(&INPUT_PARMS {
                    bDisplay: 1,
                    bDisplayCompositeResults: 1,
                    bDisplayEachComponentINChI: 1,
                    bDisplayIfRestoreWarnings: 1,
                    bINChIOutputOptions: i32::MAX,
                    ..INPUT_PARMS::default()
                }),
                Some(&mut success_input_data),
                Some(&mut success_structure),
                i32::MAX,
                i32::MIN,
                i64::MAX,
                0,
            ),
            Ok(0)
        );
        assert_eq!(success_structure.RevInChI.nRetVal, _IS_OKAY as i32);
        assert_eq!(success_structure.RevInChI.num_components, [0, 0]);
        assert_eq!(success_structure.nChargeRevrs, 0);
        assert!(success_input_data.pStrErrStruct.iter().all(|byte| *byte == 0));
    }

    #[test]
    fn source_port__ichirvr1__incrzerobonds__line_4172() {
        let untouched = inp_ATOM {
            component: 9,
            valence: 1,
            bond_type: {
                let mut value = [0; 20];
                value[0] = 3;
                value
            },
            chem_bonds_valence: 3,
            ..inp_ATOM::default()
        };
        let mut atoms = vec![untouched.clone()];
        assert_eq!(IncrZeroBonds(&mut atoms, -1, 7), Ok(()));
        assert_eq!(atoms, vec![untouched.clone()]);
        assert_eq!(IncrZeroBonds(&mut atoms, 0, 7), Ok(()));
        assert_eq!(atoms, vec![untouched]);

        let mut first = inp_ATOM {
            component: 8,
            endpoint: 17,
            valence: 4,
            chem_bonds_valence: 126,
            ..inp_ATOM::default()
        };
        first.bond_type[..4].copy_from_slice(&[0, 1, 2, 0]);
        let mut second = inp_ATOM {
            component: 8,
            endpoint: 19,
            valence: -1,
            chem_bonds_valence: -7,
            ..inp_ATOM::default()
        };
        second.bond_type[0] = 0;
        let third = inp_ATOM {
            component: 21,
            valence: 1,
            ..inp_ATOM::default()
        };
        let mut atoms = vec![first, second, third.clone()];
        assert_eq!(IncrZeroBonds(&mut atoms, 2, -1), Ok(()));
        assert_eq!(atoms[0].component, u16::MAX);
        assert_eq!(atoms[0].endpoint, 17);
        assert_eq!(&atoms[0].bond_type[..4], &[1, 1, 2, 1]);
        assert_eq!(atoms[0].chem_bonds_valence, i8::MIN);
        assert_eq!(atoms[1].component, u16::MAX);
        assert_eq!(atoms[1].endpoint, 19);
        assert_eq!(atoms[1].bond_type[0], 0);
        assert_eq!(atoms[1].chem_bonds_valence, -7);
        assert_eq!(atoms[2], third);

        assert_eq!(IncrZeroBonds(&mut atoms, 1, 65_536), Ok(()));
        assert_eq!(atoms[0].component, 0);
        assert_eq!(&atoms[0].bond_type[..4], &[1, 1, 2, 1]);
        assert_eq!(atoms[0].chem_bonds_valence, i8::MIN);
    }

    #[test]
    fn source_port__ichirvr1__clearendpts__line_4191() {
        let original = vec![
            inp_ATOM {
                endpoint: 1,
                component: 11,
                charge: -3,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                endpoint: u16::MAX,
                component: 13,
                charge: 4,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                endpoint: 17,
                component: 19,
                charge: 5,
                ..inp_ATOM::default()
            },
        ];
        let mut atoms = original.clone();
        assert_eq!(ClearEndpts(&mut atoms, -1), Ok(()));
        assert_eq!(atoms, original);
        assert_eq!(ClearEndpts(&mut atoms, 0), Ok(()));
        assert_eq!(atoms, original);

        assert_eq!(ClearEndpts(&mut atoms, 2), Ok(()));
        assert_eq!(atoms[0].endpoint, 0);
        assert_eq!(atoms[1].endpoint, 0);
        assert_eq!(atoms[2].endpoint, 17);
        assert_eq!(atoms[0].component, 11);
        assert_eq!(atoms[1].component, 13);
        assert_eq!(atoms[2].component, 19);
        assert_eq!(atoms[0].charge, -3);
        assert_eq!(atoms[1].charge, 4);
        assert_eq!(atoms[2].charge, 5);

        assert_eq!(ClearEndpts(&mut atoms, 3), Ok(()));
        assert!(atoms.iter().all(|atom| atom.endpoint == 0));
    }

    #[test]
    fn source_port__ichirvr1__connectdisconnectedh__line_5480() {
        let mut base = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 2,
            num_H: 5,
            num_iso_H: [1, 1, 1],
            ..inp_ATOM::default()
        };
        base.neighbor[0] = 9;
        base.bond_stereo[0] = -4;
        base.bond_type[0] = 2;
        base.sb_parity[..2].copy_from_slice(&[1, 2]);
        base.sb_ord[..2].copy_from_slice(&[0, 2]);
        base.sn_ord[..2].copy_from_slice(&[-1, 0]);
        base.sn_orig_at_num[0] = 102;
        let mut atoms = vec![
            base,
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 0;
                    value
                },
                orig_at_number: 101,
                iso_atw_diff: 0,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 0;
                    value
                },
                orig_at_number: 102,
                iso_atw_diff: 2,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 0;
                    value
                },
                orig_at_number: 103,
                iso_atw_diff: 3,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(ConnectDisconnectedH(&mut atoms, 1, 3), Ok(4));
        assert_eq!(&atoms[0].neighbor[..4], &[1, 2, 3, 9]);
        assert_eq!(&atoms[0].bond_stereo[..4], &[0, 0, 0, -4]);
        assert_eq!(&atoms[0].bond_type[..4], &[1, 1, 1, 2]);
        assert_eq!(&atoms[0].sb_ord[..2], &[3, 5]);
        assert_eq!(&atoms[0].sn_ord[..2], &[1, 3]);
        assert_eq!(atoms[0].valence, 4);
        assert_eq!(atoms[0].chem_bonds_valence, 5);
        assert_eq!(atoms[0].num_iso_H, [1, 0, 0]);
        assert_eq!(atoms[0].num_H, 1);
        assert!(atoms[1..].iter().all(|atom| atom.chem_bonds_valence == 1));

        let mut isotope_only = vec![inp_ATOM {
            num_H: 6,
            num_iso_H: [1, 2, 3],
            ..inp_ATOM::default()
        }];
        assert_eq!(ConnectDisconnectedH(&mut isotope_only, 1, 0), Ok(1));
        assert_eq!(isotope_only[0].num_H, 0);

        let mut too_many = vec![
            inp_ATOM {
                num_H: 0,
                ..inp_ATOM::default()
            },
            inp_ATOM::default(),
        ];
        assert_eq!(ConnectDisconnectedH(&mut too_many, 1, 1), Ok(RI_ERR_PROGR));
        assert_eq!(too_many[0].valence, 0);

        let mut maxval = vec![
            inp_ATOM {
                valence: MAXVAL as i8,
                num_H: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM::default(),
        ];
        assert_eq!(ConnectDisconnectedH(&mut maxval, 1, 1), Ok(RI_ERR_SYNTAX));
        assert_eq!(maxval[0].valence, MAXVAL as i8);

        let mut missing_stereo_h = vec![
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 7;
                    value
                },
                bond_stereo: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 6;
                    value
                },
                bond_type: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 2;
                    value
                },
                valence: 1,
                num_H: 1,
                sb_ord: [4, 0, 0],
                sn_ord: [-1, 0, 0],
                sb_parity: [1, 0, 0],
                sn_orig_at_num: [999, 0, 0],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                orig_at_number: 100,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(ConnectDisconnectedH(&mut missing_stereo_h, 1, 1), Ok(RI_ERR_PROGR));
        assert_eq!(&missing_stereo_h[0].neighbor[..2], &[1, 7]);
        assert_eq!(&missing_stereo_h[0].bond_stereo[..2], &[0, 6]);
        assert_eq!(&missing_stereo_h[0].bond_type[..2], &[1, 2]);
        assert_eq!(missing_stereo_h[0].sb_ord[0], 5);
        assert_eq!(missing_stereo_h[0].valence, 1);

        let mut isotope_too_large = vec![
            inp_ATOM {
                num_H: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                iso_atw_diff: (NUM_H_ISOTOPES + 1) as i8,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(ConnectDisconnectedH(&mut isotope_too_large, 1, 1), Ok(RI_ERR_PROGR));
        assert_eq!(isotope_too_large[0].valence, 1);
        assert_eq!(isotope_too_large[1].chem_bonds_valence, 1);

        let mut isotope_underflow = vec![
            inp_ATOM {
                num_H: 1,
                num_iso_H: [0, 0, 0],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                iso_atw_diff: 1,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(ConnectDisconnectedH(&mut isotope_underflow, 1, 1), Ok(RI_ERR_PROGR));
        assert_eq!(isotope_underflow[0].num_iso_H[0], -1);

        let mut final_underflow = vec![inp_ATOM {
            num_H: 1,
            num_iso_H: [1, 1, 0],
            ..inp_ATOM::default()
        }];
        assert_eq!(ConnectDisconnectedH(&mut final_underflow, 1, 0), Ok(RI_ERR_PROGR));
        assert_eq!(final_underflow[0].num_H, -1);
    }

    #[test]
    fn source_port__ichirvr1__disconnectedconnectedh__line_5594() {
        let mut base = inp_ATOM {
            valence: 3,
            chem_bonds_valence: 4,
            num_H: 1,
            num_iso_H: [1, 2, 3],
            sb_parity: [1, 2, 0],
            sb_ord: [2, 5, 0],
            sn_ord: [1, 2, 0],
            ..inp_ATOM::default()
        };
        base.neighbor[..3].copy_from_slice(&[2, 3, 1]);
        base.bond_stereo[..3].copy_from_slice(&[-2, -3, 7]);
        base.bond_type[..3].copy_from_slice(&[1, 1, 2]);
        let heavy = inp_ATOM::default();
        let mut protium = inp_ATOM {
            iso_atw_diff: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        protium.neighbor[0] = 0;
        let mut tritium = inp_ATOM {
            iso_atw_diff: 3,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        tritium.neighbor[0] = 0;
        let mut atoms = vec![base, heavy, protium, tritium];
        assert_eq!(DisconnectedConnectedH(&mut atoms, 2, 2), Ok(4));
        assert_eq!(atoms[0].num_H, 9);
        assert_eq!(atoms[0].num_iso_H, [2, 2, 4]);
        assert_eq!((atoms[0].valence, atoms[0].chem_bonds_valence), (1, 2));
        assert_eq!(&atoms[0].neighbor[..3], &[1, 0, 0]);
        assert_eq!(&atoms[0].bond_stereo[..3], &[7, 0, 0]);
        assert_eq!(&atoms[0].bond_type[..3], &[2, 0, 0]);
        assert_eq!(&atoms[0].sb_ord[..2], &[0, 3]);
        assert_eq!(&atoms[0].sn_ord[..2], &[-1, 2]);
        assert_eq!(atoms[2].chem_bonds_valence, 0);
        assert_eq!(atoms[3].chem_bonds_valence, 0);

        let mut isotope_only = vec![inp_ATOM {
            num_H: 1,
            num_iso_H: [1, 2, 3],
            ..inp_ATOM::default()
        }];
        assert_eq!(DisconnectedConnectedH(&mut isotope_only, 1, 0), Ok(1));
        assert_eq!(isotope_only[0].num_H, 7);

        let mut wrapping = vec![inp_ATOM {
            num_H: i8::MAX,
            num_iso_H: [1, 0, 0],
            ..inp_ATOM::default()
        }];
        assert_eq!(DisconnectedConnectedH(&mut wrapping, 1, 0), Ok(1));
        assert_eq!(wrapping[0].num_H, i8::MIN);

        let mut negative_deleted = vec![
            inp_ATOM {
                num_H: 3,
                num_iso_H: [1, 0, 0],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                num_H: 4,
                num_iso_H: [0, 2, 0],
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(DisconnectedConnectedH(&mut negative_deleted, 2, -1), Ok(1));
        assert_eq!((negative_deleted[0].num_H, negative_deleted[1].num_H), (4, 6));

        let mut mismatch = vec![
            inp_ATOM {
                valence: 1,
                neighbor: [0; MAXVAL as usize],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                neighbor: [0; MAXVAL as usize],
                chem_bonds_valence: 1,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(DisconnectedConnectedH(&mut mismatch, 1, 1), Ok(RI_ERR_PROGR));
        assert_eq!(mismatch[1].chem_bonds_valence, 0);
        assert_eq!(mismatch[0].valence, 1);

        let mut isotope_error = vec![
            inp_ATOM {
                valence: 1,
                chem_bonds_valence: 1,
                num_H: 5,
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 1;
                    value
                },
                ..inp_ATOM::default()
            },
            inp_ATOM {
                iso_atw_diff: (NUM_H_ISOTOPES + 1) as i8,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(DisconnectedConnectedH(&mut isotope_error, 1, 1), Ok(RI_ERR_PROGR));
        assert_eq!(
            (
                isotope_error[0].valence,
                isotope_error[0].chem_bonds_valence,
                isotope_error[0].num_H,
            ),
            (0, 0, 5)
        );
        assert_eq!(isotope_error[1].chem_bonds_valence, 0);

        let mut short = vec![inp_ATOM::default()];
        assert_eq!(
            DisconnectedConnectedH(&mut short, 2, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            DisconnectedConnectedH(&mut [], i32::MAX, 1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );

        let mut oversized = vec![
            inp_ATOM {
                valence: 21,
                neighbor: [1; MAXVAL as usize],
                ..inp_ATOM::default()
            },
            inp_ATOM::default(),
        ];
        assert_eq!(
            DisconnectedConnectedH(&mut oversized, 1, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(oversized[1].chem_bonds_valence, 0);
    }

    #[test]
    fn source_port__ichirvr1__addcgroups2tcgbnstruct__line_3249() {
        let mut untouched_network = BN_STRUCT {
            num_vertices: i32::MAX,
            max_vertices: i32::MIN,
            ..BN_STRUCT::default()
        };
        assert_eq!(
            AddCGroups2TCGBnStruct(
                &mut SourceHeap::default(),
                &mut untouched_network,
                &StrFromINChI::default(),
                &mut [],
                &mut ALL_TC_GROUPS::default(),
                i32::MAX,
            ),
            Ok(0)
        );
        assert_eq!(untouched_network.num_vertices, i32::MAX);

        let mut overflow_network = BN_STRUCT {
            num_vertices: 1,
            max_vertices: 1,
            ..BN_STRUCT::default()
        };
        let mut overflow_groups = ALL_TC_GROUPS {
            num_tc_groups: 1,
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(
            AddCGroups2TCGBnStruct(
                &mut SourceHeap::default(),
                &mut overflow_network,
                &StrFromINChI::default(),
                &mut [],
                &mut overflow_groups,
                0,
            ),
            Ok(BNS_VERT_EDGE_OVFL)
        );

        let fixture = |list_index: i8, group_type: i32, group_edges: i32| {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate(vec![inp_ATOM::default()]).unwrap();
            let restore_mode = heap.allocate(vec![SRM::default()]).unwrap().as_const();
            let incident = heap.allocate(vec![-1_i32; 24]).unwrap();
            let mut vertex_values = vec![BNS_VERTEX::default(); 8];
            vertex_values[0] = BNS_VERTEX {
                type_: BNS_VERT_TYPE_ATOM as AT_NUMB,
                max_adj_edges: 4,
                iedge: incident,
                ..BNS_VERTEX::default()
            };
            let vertices = heap.allocate(vertex_values).unwrap();
            let edges = heap.allocate(vec![BNS_EDGE::default(); 12]).unwrap();
            let group_pointer = heap
                .allocate(vec![TC_GROUP {
                    type_: group_type,
                    num_edges: group_edges,
                    ..TC_GROUP::default()
                }])
                .unwrap();
            let groups = ALL_TC_GROUPS {
                pTCG: group_pointer,
                num_tc_groups: 1,
                max_tc_groups: 1,
                nGroup: [-1; 18],
                ..ALL_TC_GROUPS::default()
            };
            let network = BN_STRUCT {
                num_atoms: 1,
                num_vertices: 1,
                max_vertices: 8,
                max_edges: 12,
                max_iedges: 24,
                vert: vertices,
                edge: edges,
                iedge: incident,
                ..BN_STRUCT::default()
            };
            let structure = StrFromINChI {
                at: atoms,
                num_atoms: 1,
                pSrm: restore_mode,
                ..StrFromINChI::default()
            };
            let valence_atoms = vec![crate::source_types::VAL_AT {
                cnListIndex: list_index,
                ..crate::source_types::VAL_AT::default()
            }];
            (
                heap,
                network,
                structure,
                valence_atoms,
                groups,
                vertices,
                edges,
                incident,
                group_pointer,
            )
        };

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, vertices, edges, incident, group_pointer) =
            fixture(11, BNS_VT_C_POS as i32, 1);
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(1)
        );
        assert_eq!(
            (
                network.num_vertices,
                network.num_edges,
                network.num_c_groups,
                network.tot_st_cap,
                network.tot_st_flow,
            ),
            (2, 1, 1, 1, 1)
        );
        let vertex_values = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            (
                vertex_values[0].st_edge.cap,
                vertex_values[0].st_edge.flow,
                vertex_values[0].num_adj_edges,
            ),
            (1, 1, 1)
        );
        assert_eq!(
            (
                vertex_values[1].type_,
                vertex_values[1].max_adj_edges,
                vertex_values[1].num_adj_edges,
                vertex_values[1].iedge,
            ),
            (BNS_VT_C_POS as AT_NUMB, 3, 1, incident.offset(4).unwrap())
        );
        let edge = &heap.slice(edges.as_const()).unwrap()[0];
        assert_eq!(
            (
                edge.neighbor1,
                edge.neighbor12,
                edge.cap,
                edge.cap0,
                edge.flow,
                edge.flow0,
                edge.forbidden,
            ),
            (0, 1, 1, 1, 1, 1, 0)
        );
        assert_eq!(valence_atoms[0].nCPlusGroupEdge, 1);
        assert_eq!(heap.slice(group_pointer.as_const()).unwrap()[0].nVertexNumber, 1);
        assert_eq!(
            &heap.slice(incident.as_const()).unwrap()[..7],
            &[0, -1, -1, -1, 0, -1, -1]
        );

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, _, edges, _, group_pointer) =
            fixture(15, BNS_VT_C_POS as i32, 1);
        heap.slice_mut(group_pointer).unwrap()[0].st_cap = 2;
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(1)
        );
        let edge = &heap.slice(edges.as_const()).unwrap()[0];
        assert_eq!(
            (edge.cap, edge.cap0, edge.flow, edge.flow0, edge.forbidden),
            (1, 1, 0, 0, 1)
        );
        assert_eq!(valence_atoms[0].nCPlusGroupEdge, 1);
        assert_eq!((network.tot_st_cap, network.tot_st_flow), (2, 0));

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, vertices, edges, incident, _) =
            fixture(13, BNS_VT_C_NEG as i32, 1);
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(1)
        );
        assert_eq!(
            (
                network.num_vertices,
                network.num_edges,
                network.tot_st_cap,
                network.tot_st_flow,
            ),
            (3, 2, 2, 2)
        );
        let vertex_values = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            (
                vertex_values[2].type_,
                vertex_values[2].max_adj_edges,
                vertex_values[2].num_adj_edges,
                vertex_values[2].st_edge.cap,
                vertex_values[2].st_edge.flow,
                vertex_values[2].iedge,
            ),
            (
                BNS_VERT_TYPE__AUX as AT_NUMB | BNS_VERT_TYPE_TEMP as AT_NUMB,
                2,
                2,
                1,
                1,
                incident.offset(7).unwrap(),
            )
        );
        assert_eq!(
            heap.slice(edges.as_const()).unwrap()[..2]
                .iter()
                .map(|edge| { (edge.neighbor1, edge.neighbor12, edge.cap, edge.flow, edge.forbidden,) })
                .collect::<Vec<_>>(),
            vec![(0, 2, 1, 1, 0), (1, 3, 1, 0, 0)]
        );
        assert_eq!(valence_atoms[0].nCMinusGroupEdge, 2);
        assert_eq!(
            &heap.slice(incident.as_const()).unwrap()[..9],
            &[0, -1, -1, -1, 1, -1, -1, 0, 1]
        );

        let bond_fixture = |bond_type: u8, metal_neighbor: bool| {
            let mut heap = SourceHeap::default();
            let mut first_atom = inp_ATOM {
                valence: 1,
                ..inp_ATOM::default()
            };
            first_atom.neighbor[0] = 1;
            first_atom.bond_type[0] = bond_type;
            let mut second_atom = inp_ATOM {
                valence: 1,
                ..inp_ATOM::default()
            };
            second_atom.neighbor[0] = 0;
            second_atom.bond_type[0] = bond_type;
            let atoms = heap.allocate(vec![first_atom, second_atom]).unwrap();
            let restore_mode = heap.allocate(vec![SRM::default()]).unwrap().as_const();
            let incident = heap.allocate(vec![-1_i32; 28]).unwrap();
            {
                let values = heap.slice_mut(incident).unwrap();
                values[0] = 0;
                values[4] = 0;
            }
            let vertices = heap
                .allocate(vec![
                    BNS_VERTEX {
                        type_: BNS_VERT_TYPE_ATOM as AT_NUMB,
                        num_adj_edges: 1,
                        max_adj_edges: 4,
                        iedge: incident,
                        st_edge: crate::source_types::BNS_ST_EDGE {
                            cap: 5,
                            cap0: 5,
                            ..crate::source_types::BNS_ST_EDGE::default()
                        },
                    },
                    BNS_VERTEX {
                        type_: BNS_VERT_TYPE_ATOM as AT_NUMB,
                        num_adj_edges: 1,
                        max_adj_edges: 4,
                        iedge: incident.offset(4).unwrap(),
                        st_edge: crate::source_types::BNS_ST_EDGE {
                            cap: 4,
                            cap0: 4,
                            ..crate::source_types::BNS_ST_EDGE::default()
                        },
                    },
                    BNS_VERTEX::default(),
                    BNS_VERTEX::default(),
                    BNS_VERTEX::default(),
                    BNS_VERTEX::default(),
                    BNS_VERTEX::default(),
                ])
                .unwrap();
            let edges = heap
                .allocate(vec![
                    BNS_EDGE {
                        neighbor1: 0,
                        neighbor12: 1,
                        neigh_ord: [0, 0],
                        cap: 9,
                        cap0: 9,
                        ..BNS_EDGE::default()
                    },
                    BNS_EDGE::default(),
                    BNS_EDGE::default(),
                    BNS_EDGE::default(),
                ])
                .unwrap();
            let group_pointer = heap
                .allocate(vec![TC_GROUP {
                    type_: BNS_VT_C_POS as i32,
                    num_edges: 1,
                    ..TC_GROUP::default()
                }])
                .unwrap();
            let mut groups = ALL_TC_GROUPS {
                pTCG: group_pointer,
                num_tc_groups: 1,
                max_tc_groups: 1,
                nGroup: [-1; 18],
                ..ALL_TC_GROUPS::default()
            };
            let mut network = BN_STRUCT {
                num_atoms: 2,
                num_vertices: 2,
                num_edges: 1,
                max_vertices: 7,
                max_edges: 4,
                max_iedges: 28,
                vert: vertices,
                edge: edges,
                iedge: incident,
                tot_st_cap: 9,
                ..BN_STRUCT::default()
            };
            let structure = StrFromINChI {
                at: atoms,
                num_atoms: 2,
                pSrm: restore_mode,
                ..StrFromINChI::default()
            };
            let mut valence_atoms = vec![
                crate::source_types::VAL_AT {
                    cnListIndex: 11,
                    cMetal: i8::from(metal_neighbor),
                    cNumBondsToMetal: i8::from(metal_neighbor),
                    ..crate::source_types::VAL_AT::default()
                },
                crate::source_types::VAL_AT::default(),
            ];
            assert_eq!(
                AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
                Ok(1)
            );
            (
                heap.slice(edges.as_const()).unwrap()[0].cap,
                network.num_edges,
                network.tot_st_cap,
            )
        };
        assert_eq!(bond_fixture(BOND_TYPE_SINGLE as u8, false), (2, 2, 10));
        assert_eq!(bond_fixture(BOND_TYPE_SINGLE as u8, true), (3, 2, 10));
        assert_eq!(bond_fixture(BOND_TYPE_DOUBLE as u8, false), (9, 2, 10));

        let supergroup_fixture = |charge: i32| {
            let mut heap = SourceHeap::default();
            let restore_mode = heap.allocate(vec![SRM::default()]).unwrap().as_const();
            let incident = heap.allocate(vec![-1_i32; 48]).unwrap();
            let mut vertex_values = vec![BNS_VERTEX::default(); 10];
            vertex_values[0] = BNS_VERTEX {
                type_: BNS_VERT_TYPE_ATOM as AT_NUMB,
                max_adj_edges: 4,
                iedge: incident,
                ..BNS_VERTEX::default()
            };
            let vertices = heap.allocate(vertex_values).unwrap();
            let edges = heap.allocate(vec![BNS_EDGE::default(); 12]).unwrap();
            let group_pointer = heap
                .allocate(vec![
                    TC_GROUP {
                        type_: BNS_VT_C_POS as i32,
                        st_cap: 2,
                        ..TC_GROUP::default()
                    },
                    TC_GROUP {
                        type_: BNS_VT_C_POS_ALL as i32,
                        ..TC_GROUP::default()
                    },
                ])
                .unwrap();
            let mut slots = [-1; 18];
            slots[TCG_Plus0 as usize] = 0;
            slots[TCG_Plus as usize] = 1;
            let mut groups = ALL_TC_GROUPS {
                pTCG: group_pointer,
                num_tc_groups: 2,
                max_tc_groups: 2,
                nGroup: slots,
                total_charge: charge,
                ..ALL_TC_GROUPS::default()
            };
            let mut network = BN_STRUCT {
                num_vertices: 1,
                max_vertices: 10,
                max_edges: 12,
                max_iedges: 48,
                vert: vertices,
                edge: edges,
                iedge: incident,
                ..BN_STRUCT::default()
            };
            let structure = StrFromINChI {
                pSrm: restore_mode,
                ..StrFromINChI::default()
            };
            let result = AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut [], &mut groups, 4);
            (
                result,
                network,
                groups,
                heap.slice(vertices.as_const()).unwrap().to_vec(),
                heap.slice(edges.as_const()).unwrap().to_vec(),
            )
        };
        let (result, network, groups, vertices, edges) = supergroup_fixture(3);
        assert_eq!(result, Ok(2));
        assert_eq!((network.num_vertices, network.num_edges), (5, 3));
        assert_eq!((groups.nEdgePlus, groups.nEdge4charge, groups.added_charge), (2, 2, 3));
        assert_eq!(vertices[4].st_edge.cap, 3);
        assert_eq!(edges[2].neighbor12, 6);

        let (result, network, groups, vertices, edges) = supergroup_fixture(-3);
        assert_eq!(result, Ok(2));
        assert_eq!((network.num_vertices, network.num_edges), (5, 3));
        assert_eq!((groups.nEdgePlus, groups.nEdge4charge, groups.added_charge), (2, 2, -3));
        assert_eq!(vertices[2].st_edge.cap, 5);
        assert_eq!(edges[2].cap, 3);

        let mut heap = SourceHeap::default();
        let mut metal_atom = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        metal_atom.neighbor[0] = 1;
        metal_atom.bond_type[0] = BOND_TYPE_SINGLE as u8;
        let mut ligand_atom = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        ligand_atom.neighbor[0] = 0;
        ligand_atom.bond_type[0] = BOND_TYPE_SINGLE as u8;
        let atoms = heap.allocate(vec![metal_atom, ligand_atom]).unwrap();
        let restore_mode = heap
            .allocate(vec![SRM {
                bMetalAddFlower: 1,
                nMetalMaxCharge_D: 2,
                nMetalFlowerParam_D: 2,
                ..SRM::default()
            }])
            .unwrap()
            .as_const();
        let incident = heap.allocate(vec![-1_i32; 48]).unwrap();
        {
            let values = heap.slice_mut(incident).unwrap();
            values[0] = 0;
            values[4] = 0;
        }
        let mut initial_vertices = vec![BNS_VERTEX::default(); 12];
        initial_vertices[0] = BNS_VERTEX {
            type_: BNS_VERT_TYPE_ATOM as AT_NUMB,
            num_adj_edges: 1,
            max_adj_edges: 4,
            iedge: incident,
            st_edge: crate::source_types::BNS_ST_EDGE {
                cap: 1,
                cap0: 1,
                ..crate::source_types::BNS_ST_EDGE::default()
            },
        };
        initial_vertices[1] = BNS_VERTEX {
            type_: BNS_VERT_TYPE_ATOM as AT_NUMB,
            num_adj_edges: 1,
            max_adj_edges: 2,
            iedge: incident.offset(4).unwrap(),
            st_edge: crate::source_types::BNS_ST_EDGE {
                cap: 1,
                cap0: 1,
                ..crate::source_types::BNS_ST_EDGE::default()
            },
        };
        let vertices = heap.allocate(initial_vertices).unwrap();
        let mut initial_edges = vec![BNS_EDGE::default(); 16];
        initial_edges[0] = BNS_EDGE {
            neighbor1: 0,
            neighbor12: 1,
            neigh_ord: [0, 0],
            cap: 2,
            cap0: 2,
            ..BNS_EDGE::default()
        };
        let edges = heap.allocate(initial_edges).unwrap();
        let group_pointer = heap
            .allocate(vec![
                TC_GROUP {
                    type_: BNS_VT_C_POS_M as i32,
                    num_edges: 1,
                    st_cap: 2,
                    st_flow: 2,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: BNS_VT_C_NEG_M as i32,
                    num_edges: 1,
                    st_cap: 2,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: BNS_VT_M_GROUP as i32,
                    num_edges: 1,
                    edges_cap: 5,
                    edges_flow: 5,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: BNS_VT_M_GROUP as i32,
                    ord_num: 1,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: BNS_VT_M_GROUP as i32,
                    ord_num: 2,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    type_: BNS_VT_M_GROUP as i32,
                    ord_num: 3,
                    ..TC_GROUP::default()
                },
            ])
            .unwrap();
        let mut slots = [-1; 18];
        slots[TCG_MeFlower0 as usize] = 2;
        slots[TCG_MeFlower1 as usize] = 3;
        slots[TCG_MeFlower2 as usize] = 4;
        slots[TCG_MeFlower3 as usize] = 5;
        let mut groups = ALL_TC_GROUPS {
            pTCG: group_pointer,
            num_tc_groups: 6,
            max_tc_groups: 6,
            nGroup: slots,
            num_metal_atoms: 1,
            ..ALL_TC_GROUPS::default()
        };
        let mut network = BN_STRUCT {
            num_atoms: 2,
            num_vertices: 2,
            num_edges: 1,
            max_vertices: 12,
            max_edges: 16,
            max_iedges: 48,
            vert: vertices,
            edge: edges,
            iedge: incident,
            tot_st_cap: 2,
            ..BN_STRUCT::default()
        };
        let structure = StrFromINChI {
            at: atoms,
            num_atoms: 2,
            pSrm: restore_mode,
            ..StrFromINChI::default()
        };
        let mut valence_atoms = vec![
            crate::source_types::VAL_AT {
                cMetal: 1,
                cNumBondsToMetal: 1,
                cnListIndex: 18,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT::default(),
        ];
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 3,),
            Ok(1)
        );
        assert_eq!((network.num_vertices, network.num_edges), (9, 10));
        assert_eq!(valence_atoms[0].nMetalGroupEdge, 5);
        let edge_values = heap.slice(edges.as_const()).unwrap();
        assert_eq!((edge_values[4].cap, edge_values[4].flow), (5, 5));
        assert_eq!(
            edge_values[5..10]
                .iter()
                .map(|edge| (edge.cap, edge.flow))
                .collect::<Vec<_>>(),
            vec![(6, 2), (7, 2), (6, 4), (2, 0), (2, 0)]
        );
        let vertex_values = heap.slice(vertices.as_const()).unwrap();
        assert_eq!(
            vertex_values[4..8]
                .iter()
                .map(|vertex| (vertex.st_edge.cap, vertex.st_edge.flow))
                .collect::<Vec<_>>(),
            vec![(9, 9), (6, 6), (6, 6), (0, 0)]
        );

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, vertices, _, _, _) =
            fixture(11, BNS_VT_C_NEG as i32, 1);
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(BNS_BOND_ERR)
        );
        assert_eq!((network.num_vertices, network.num_edges), (1, 0));
        assert_eq!(
            (
                heap.slice(vertices.as_const()).unwrap()[0].st_edge.cap,
                heap.slice(vertices.as_const()).unwrap()[0].st_edge.flow,
            ),
            (1, 1)
        );

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, _, _, _, group_pointer) =
            fixture(11, BNS_VT_C_POS as i32, 1);
        heap.slice_mut(group_pointer).unwrap()[0].ord_num = 1;
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(RI_ERR_PROGR)
        );

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, vertices, _, _, _) =
            fixture(13, BNS_VT_C_NEG as i32, 1);
        network.max_vertices = 3;
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!((network.num_vertices, network.num_edges), (1, 0));
        assert_eq!(
            heap.slice(vertices.as_const()).unwrap()[2].type_,
            BNS_VERT_TYPE__AUX as AT_NUMB | BNS_VERT_TYPE_TEMP as AT_NUMB
        );

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, _, edges, _, _) =
            fixture(11, BNS_VT_C_POS as i32, 1);
        network.max_edges = 0;
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!((network.num_vertices, network.num_edges), (1, 0));
        assert_eq!(
            (
                heap.slice(edges.as_const()).unwrap()[0].cap,
                heap.slice(edges.as_const()).unwrap()[0].flow,
                heap.slice(edges.as_const()).unwrap()[0].forbidden,
            ),
            (1, 1, 0)
        );

        let (mut heap, mut network, structure, mut valence_atoms, mut groups, _, edges, _, _) =
            fixture(11, BNS_VT_C_POS as i32, 1);
        heap.slice_mut(edges).unwrap()[0].neighbor12 = 7;
        assert_eq!(
            AddCGroups2TCGBnStruct(&mut heap, &mut network, &structure, &mut valence_atoms, &mut groups, 2,),
            Ok(BNS_PROGRAM_ERR)
        );
        assert_eq!((network.num_vertices, network.num_edges), (1, 0));
        assert_eq!(
            (
                heap.slice(edges.as_const()).unwrap()[0].neighbor12,
                heap.slice(edges.as_const()).unwrap()[0].cap,
                heap.slice(edges.as_const()).unwrap()[0].flow,
                heap.slice(edges.as_const()).unwrap()[0].forbidden,
            ),
            (7, 1, 1, 0)
        );
    }

    #[test]
    fn source_port__ichirvr1__nnumedgestocnvertex__line_3824() {
        for (list, expected_valences) in CN_LIST.iter().zip(CN_LIST_VALENCES) {
            assert_eq!(list.nodes.len(), expected_valences.len());
            let expected_incident_edges: &[i32] = if list.bits == -1 {
                &[3, 2, 1, 1, 1]
            } else {
                expected_valences
            };
            for (vertex, &expected) in expected_incident_edges.iter().enumerate() {
                assert_eq!(
                    nNumEdgesToCnVertex(list.nodes, list.nodes.len() as i32, vertex as i32),
                    Ok(expected),
                    "bits={}, vertex={vertex}",
                    list.bits,
                );
            }
            assert_eq!(nNumEdgesToCnVertex(list.nodes, list.nodes.len() as i32, -1), Ok(0));
            assert_eq!(
                nNumEdgesToCnVertex(list.nodes, list.nodes.len() as i32, list.nodes.len() as i32,),
                Ok(0)
            );
        }

        assert_eq!(nNumEdgesToCnVertex(&[], i32::MIN, 0), Ok(0));
        assert_eq!(nNumEdgesToCnVertex(&[], 0, i32::MAX), Ok(0));
        assert_eq!(nNumEdgesToCnVertex(CN_LIST[0].nodes, 1, 0), Ok(1));
        assert_eq!(nNumEdgesToCnVertex(CN_LIST[0].nodes, 2, 1), Ok(3));
        assert_eq!(
            nNumEdgesToCnVertex(CN_LIST[0].nodes, 6, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr1__allocateandinittcgbnstruct__line_3879() {
        fn two_atom_fixture(
            heap: &mut SourceHeap,
            bond_type: u8,
        ) -> (StrFromINChI, Vec<crate::source_types::VAL_AT>, ALL_TC_GROUPS) {
            let mut first = inp_ATOM {
                el_number: 6,
                valence: 1,
                chem_bonds_valence: (bond_type & BOND_TYPE_MASK as u8).min(3) as i8,
                ..inp_ATOM::default()
            };
            first.neighbor[0] = 1;
            first.bond_type[0] = bond_type;
            let mut second = first.clone();
            second.neighbor[0] = 0;
            let atoms = heap.allocate_model_storage(vec![first, second]).unwrap();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            (
                StrFromINChI {
                    at: atoms,
                    num_atoms: 2,
                    pSrm: restore_mode.as_const(),
                    ..StrFromINChI::default()
                },
                vec![
                    crate::source_types::VAL_AT::default(),
                    crate::source_types::VAL_AT::default(),
                ],
                ALL_TC_GROUPS {
                    nVertices: 2,
                    nEdges: 1,
                    num_atoms: 2,
                    num_bonds: 1,
                    ..ALL_TC_GROUPS::default()
                },
            )
        }

        for successful_allocations in 0..6_u64 {
            let mut heap = SourceHeap::default();
            let (structure, valence_atoms, groups) = two_atom_fixture(&mut heap, BOND_TYPE_SINGLE as u8);
            let baseline = heap.live_allocation_count();
            heap.fail_after_allocations(successful_allocations);
            let mut changed = -7;
            assert_eq!(
                AllocateAndInitTCGBnStruct(&mut heap, &structure, &valence_atoms, &groups, 0, 0, 2, &mut changed,),
                Ok(SourceMutPointer::null()),
                "allocation failure after {successful_allocations} successful allocations",
            );
            assert_eq!(heap.source_allocation_calls(), successful_allocations + 1);
            assert_eq!(heap.live_allocation_count(), baseline);
            assert_eq!(changed, -7);
        }

        let mut capacity_heap = SourceHeap::default();
        let (capacity_structure, capacity_va, mut capacity_groups) =
            two_atom_fixture(&mut capacity_heap, BOND_TYPE_SINGLE as u8);
        capacity_groups.nAddIedges = -6;
        let capacity_baseline = capacity_heap.live_allocation_count();
        let mut capacity_changed = 19;
        assert_eq!(
            AllocateAndInitTCGBnStruct(
                &mut capacity_heap,
                &capacity_structure,
                &capacity_va,
                &capacity_groups,
                0,
                0,
                0,
                &mut capacity_changed,
            ),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(capacity_heap.live_allocation_count(), capacity_baseline);
        assert_eq!(capacity_changed, 19);

        for (bond_type, expected_flow) in [
            (BOND_TYPE_SINGLE as u8, 0),
            (BOND_TYPE_DOUBLE as u8, 1),
            (BOND_TYPE_TRIPLE as u8, 2),
        ] {
            let mut heap = SourceHeap::default();
            let (structure, valence_atoms, groups) = two_atom_fixture(&mut heap, bond_type);
            let baseline = heap.live_allocation_count();
            let mut changed = -1;
            let network =
                AllocateAndInitTCGBnStruct(&mut heap, &structure, &valence_atoms, &groups, 0, 0, 1, &mut changed)
                    .unwrap();
            assert!(!network.is_null());
            assert_eq!(changed, 0);
            let value = heap.slice(network.as_const()).unwrap()[0].clone();
            assert_eq!(
                (
                    value.num_atoms,
                    value.num_vertices,
                    value.num_bonds,
                    value.num_edges,
                    value.num_iedges,
                ),
                (2, 2, 1, 1, 2)
            );
            assert_eq!((value.max_vertices, value.max_edges, value.max_iedges), (2, 5, 10));
            assert_eq!(
                (value.tot_st_cap, value.tot_st_flow),
                (2 * expected_flow, 2 * expected_flow)
            );
            let edge = &heap.slice(value.edge.as_const()).unwrap()[0];
            assert_eq!(
                (
                    edge.neighbor1,
                    edge.neighbor12,
                    edge.neigh_ord,
                    edge.cap,
                    edge.flow,
                    edge.cap0,
                    edge.flow0,
                    edge.forbidden,
                ),
                (0, 1, [0, 0], 2, expected_flow, 2, expected_flow, 0)
            );
            let vertices = heap.slice(value.vert.as_const()).unwrap();
            assert_eq!(vertices[0].iedge, value.iedge);
            assert_eq!(vertices[1].iedge, value.iedge.offset(1).unwrap());
            assert_eq!((vertices[0].num_adj_edges, vertices[1].num_adj_edges), (1, 1));
            assert_eq!(heap.slice(value.iedge.as_const()).unwrap()[..2], [0, 0]);
            let path = heap.slice(value.altp[0].as_const()).unwrap();
            assert_eq!(value.len_alt_path, 24);
            assert_eq!(path[tagAltPathConst_iALTP_MAX_LEN as usize].number(), 24);
            assert_eq!(path[tagAltPathConst_iALTP_FLOW as usize].flow(0), 0);
            assert_eq!(path[tagAltPathConst_iALTP_START_ATOM as usize].number(), NO_VERTEX);
            assert_eq!(path[tagAltPathConst_iALTP_END_ATOM as usize].number(), NO_VERTEX);
            assert_eq!(path[tagAltPathConst_iALTP_PATH_LEN as usize].number(), 0);
            DeAllocateBnStruct(&mut heap, network).unwrap();
            assert_eq!(heap.live_allocation_count(), baseline);
        }

        let mut heap = SourceHeap::default();
        let mut atom0 = inp_ATOM {
            el_number: 6,
            valence: 1,
            chem_bonds_valence: 2,
            endpoint: 1,
            ..inp_ATOM::default()
        };
        atom0.neighbor[0] = 1;
        atom0.bond_type[0] = BOND_TYPE_DOUBLE as u8;
        atom0.sb_parity[0] = 1;
        atom0.sb_ord[0] = 0;
        let mut atom1 = inp_ATOM {
            el_number: 7,
            valence: 2,
            chem_bonds_valence: 3,
            ..inp_ATOM::default()
        };
        atom1.neighbor[..2].copy_from_slice(&[0, 2]);
        atom1.bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, 0x89]);
        let mut atom2 = inp_ATOM {
            el_number: 8,
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        atom2.neighbor[0] = 1;
        atom2.bond_type[0] = 0x89;
        let atoms = heap.allocate_model_storage(vec![atom0, atom1, atom2]).unwrap();
        let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
        let structure = StrFromINChI {
            at: atoms,
            num_atoms: 3,
            pSrm: restore_mode.as_const(),
            ..StrFromINChI::default()
        };
        let valence_atoms = vec![
            crate::source_types::VAL_AT {
                cInitFreeValences: 1,
                cnListIndex: 1,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT {
                cInitFreeValences: 2,
                ..crate::source_types::VAL_AT::default()
            },
            crate::source_types::VAL_AT {
                cInitFreeValences: 3,
                ..crate::source_types::VAL_AT::default()
            },
        ];
        let groups = ALL_TC_GROUPS {
            nVertices: 3,
            nEdges: 2,
            nAddIedges: 2,
            num_atoms: 3,
            num_bonds: 2,
            num_tgroups: 4,
            ..ALL_TC_GROUPS::default()
        };
        let baseline = heap.live_allocation_count();
        let mut changed = -1;
        let network =
            AllocateAndInitTCGBnStruct(&mut heap, &structure, &valence_atoms, &groups, 0, 1, 20, &mut changed).unwrap();
        assert!(!network.is_null());
        assert_eq!(changed, 1);
        let value = heap.slice(network.as_const()).unwrap()[0].clone();
        assert_eq!((value.max_altp, value.num_altp), (BN_MAX_ALTP as i32, 0));
        assert_eq!((value.max_vertices, value.max_edges, value.max_iedges), (3, 11, 24));
        assert_eq!((value.num_iedges, value.tot_st_cap, value.tot_st_flow), (9, 8, 2));
        let vertices = heap.slice(value.vert.as_const()).unwrap();
        assert_eq!(
            vertices[..3]
                .iter()
                .map(|vertex| (vertex.max_adj_edges, vertex.num_adj_edges))
                .collect::<Vec<_>>(),
            vec![(4, 1), (3, 2), (2, 1)]
        );
        assert_eq!(vertices[0].iedge, value.iedge);
        assert_eq!(vertices[1].iedge, value.iedge.offset(4).unwrap());
        assert_eq!(vertices[2].iedge, value.iedge.offset(7).unwrap());
        assert_eq!(
            vertices[..3]
                .iter()
                .map(|vertex| (vertex.st_edge.cap, vertex.st_edge.flow, vertex.type_))
                .collect::<Vec<_>>(),
            vec![
                (2, 1, BNS_VERT_TYPE_ATOM as AT_NUMB),
                (3, 1, BNS_VERT_TYPE_ATOM as AT_NUMB),
                (3, 0, BNS_VERT_TYPE_ATOM as AT_NUMB),
            ]
        );
        let edges = heap.slice(value.edge.as_const()).unwrap();
        assert_eq!((edges[0].flow, edges[0].cap, edges[0].forbidden), (1, 2, 1));
        assert_eq!((edges[1].flow, edges[1].cap, edges[1].forbidden), (0, 2, 0));
        assert_eq!(edges[0].neigh_ord, [0, 0]);
        assert_eq!(edges[1].neigh_ord, [1, 0]);
        let incident = heap.slice(value.iedge.as_const()).unwrap();
        assert_eq!((incident[0], incident[4], incident[5], incident[7]), (0, 0, 1, 1));
        let atoms_after = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(atoms_after[1].bond_type[1], 0x81);
        assert_eq!(atoms_after[2].bond_type[0], 0x81);
        DeAllocateBnStruct(&mut heap, network).unwrap();
        assert_eq!(heap.live_allocation_count(), baseline);

        let mut metal_heap = SourceHeap::default();
        let (metal_structure, mut metal_va, metal_groups) = two_atom_fixture(&mut metal_heap, BOND_TYPE_SINGLE as u8);
        metal_va[0].cMetal = 1;
        metal_va[0].cInitOrigValenceToMetal = 1;
        metal_va[0].cInitValenceToMetal = 1;
        metal_va[0].cInitFreeValences = 5;
        metal_va[1].cInitFreeValences = 2;
        let metal_mode = SRM {
            bMetalAddFlower: 1,
            nMetalMinBondOrder: 0,
            nMetalInitBondOrder: 1,
            nMetalInitEdgeFlow: 1,
            ..SRM::default()
        };
        let old_mode = metal_structure.pSrm.as_mut();
        metal_heap.slice_mut(old_mode).unwrap()[0] = metal_mode;
        let mut metal_changed = -1;
        let metal_network = AllocateAndInitTCGBnStruct(
            &mut metal_heap,
            &metal_structure,
            &metal_va,
            &metal_groups,
            0,
            0,
            0,
            &mut metal_changed,
        )
        .unwrap();
        let metal_value = metal_heap.slice(metal_network.as_const()).unwrap()[0].clone();
        assert_eq!((metal_value.tot_st_cap, metal_value.tot_st_flow), (4, 2));
        let metal_edge = &metal_heap.slice(metal_value.edge.as_const()).unwrap()[0];
        assert_eq!((metal_edge.cap, metal_edge.flow), (3, 1));
        let metal_vertices = metal_heap.slice(metal_value.vert.as_const()).unwrap();
        assert_eq!(
            (
                metal_vertices[0].st_edge.cap,
                metal_vertices[0].st_edge.flow,
                metal_vertices[1].st_edge.cap,
                metal_vertices[1].st_edge.flow,
            ),
            (1, 1, 3, 1)
        );
    }

    #[test]
    fn source_port__ichirvr1__incrzerobondsandclearendpts__line_4152() {
        let original = vec![
            inp_ATOM {
                endpoint: 9,
                component: 8,
                valence: 4,
                chem_bonds_valence: i8::MAX,
                bond_type: {
                    let mut values = [7; 20];
                    values[..5].copy_from_slice(&[0, 2, 0, 3, 0]);
                    values
                },
                ..inp_ATOM::default()
            },
            inp_ATOM {
                endpoint: 6,
                component: 5,
                valence: -1,
                chem_bonds_valence: 12,
                bond_type: [0; 20],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                endpoint: 4,
                component: 3,
                valence: 1,
                chem_bonds_valence: -1,
                bond_type: [0; 20],
                ..inp_ATOM::default()
            },
        ];

        let mut unchanged = original.clone();
        assert_eq!(IncrZeroBondsAndClearEndpts(&mut unchanged, i32::MIN, 99), Ok(()));
        assert_eq!(unchanged, original);
        assert_eq!(IncrZeroBondsAndClearEndpts(&mut unchanged, 0, 99), Ok(()));
        assert_eq!(unchanged, original);

        let mut atoms = original.clone();
        assert_eq!(IncrZeroBondsAndClearEndpts(&mut atoms, 2, -1), Ok(()));
        assert_eq!((atoms[0].endpoint, atoms[0].component), (0, AT_NUMB::MAX));
        assert_eq!(&atoms[0].bond_type[..5], &[1, 2, 1, 3, 0]);
        assert_eq!(atoms[0].chem_bonds_valence, i8::MIN.wrapping_add(1));
        assert_eq!((atoms[1].endpoint, atoms[1].component), (0, AT_NUMB::MAX));
        assert_eq!(atoms[1].bond_type, [0; 20]);
        assert_eq!(atoms[1].chem_bonds_valence, 12);
        assert_eq!(atoms[2], original[2]);

        assert_eq!(IncrZeroBondsAndClearEndpts(&mut atoms, 3, 65_537), Ok(()));
        assert_eq!(
            atoms.iter().map(|atom| atom.component).collect::<Vec<_>>(),
            vec![1, 1, 1]
        );
        assert_eq!(atoms[2].endpoint, 0);
        assert_eq!(atoms[2].bond_type[0], BOND_TYPE_SINGLE as u8);
        assert_eq!(atoms[2].chem_bonds_valence, 0);

        let before_error = atoms.clone();
        assert_eq!(
            IncrZeroBondsAndClearEndpts(&mut atoms, 4, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(atoms, before_error);
    }

    #[test]
    fn source_port__ichirvr1__allocbfsqueue__line_4815() {
        let mut heap = SourceHeap::default();
        let mut queue = BFS_Q {
            num_at: 7,
            min_ring_size: 11,
            ..BFS_Q::default()
        };
        let unchanged = queue.clone();
        assert_eq!(AllocBfsQueue(&mut heap, &mut queue, 0, 99), Ok(RI_ERR_PROGR));
        assert_eq!(queue, unchanged);
        assert_eq!(AllocBfsQueue(&mut heap, &mut queue, -3, 99), Ok(RI_ERR_PROGR));
        assert_eq!(queue, unchanged);

        queue = BFS_Q::default();
        assert_eq!(AllocBfsQueue(&mut heap, &mut queue, 3, -1), Ok(0));
        assert_eq!((queue.num_at, queue.min_ring_size), (3, AT_RANK::MAX));
        let first_q = queue.q;
        let first_levels = queue.nAtomLevel;
        let first_sources = queue.cSource;
        let queue_value = &heap.slice(first_q.as_const()).unwrap()[0];
        assert_eq!(queue_value.nTotLength, 4);
        assert_eq!(heap.slice(queue_value.Val.as_const()).unwrap(), &[0; 4]);
        assert_eq!(heap.slice(first_levels.as_const()).unwrap(), &[0; 3]);
        assert_eq!(heap.slice(first_sources.as_const()).unwrap(), &[0; 3]);

        heap.trace_source_allocations();
        assert_eq!(AllocBfsQueue(&mut heap, &mut queue, 3, 17), Ok(0));
        assert_eq!(heap.source_allocation_calls(), 0);
        assert_eq!(
            (queue.q, queue.nAtomLevel, queue.cSource),
            (first_q, first_levels, first_sources)
        );
        assert_eq!(queue.min_ring_size, 17);
        assert_eq!(AllocBfsQueue(&mut heap, &mut queue, 2, 19), Ok(0));
        assert_eq!(heap.source_allocation_calls(), 0);
        assert_eq!(queue.min_ring_size, 19);

        heap.trace_source_allocations();
        assert_eq!(AllocBfsQueue(&mut heap, &mut queue, 5, 23), Ok(0));
        assert_eq!(heap.source_allocation_calls(), 4);
        assert_eq!((queue.num_at, queue.min_ring_size), (5, 23));
        assert_ne!(queue.q, first_q);
        assert!(matches!(
            heap.slice(first_q.as_const()),
            Err(SourceHeapError::MissingAllocation)
        ));
        assert_eq!(heap.slice(queue.nAtomLevel.as_const()).unwrap(), &[0; 5]);
        assert_eq!(heap.slice(queue.cSource.as_const()).unwrap(), &[0; 5]);
        assert_eq!(AllocBfsQueue(&mut heap, &mut queue, BFS_Q_FREE, 999), Ok(0));
        assert_eq!(queue, BFS_Q::default());
        assert_eq!(heap.live_allocation_count(), 0);

        for successful_allocations in 0..4_u64 {
            let mut failure_heap = SourceHeap::default();
            let mut failed_queue = BFS_Q {
                min_ring_size: 31,
                ..BFS_Q::default()
            };
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                AllocBfsQueue(&mut failure_heap, &mut failed_queue, 4, 37),
                Ok(RI_ERR_ALLOC),
                "allocation failure after {successful_allocations} successful calls",
            );
            assert_eq!(failed_queue.num_at, 0);
            assert_eq!(failed_queue.min_ring_size, 31);
            assert_eq!(failed_queue.q.is_null(), successful_allocations < 2);
            assert_eq!(failed_queue.nAtomLevel.is_null(), successful_allocations == 2);
            assert_eq!(failed_queue.cSource.is_null(), successful_allocations == 3);
            assert_eq!(
                AllocBfsQueue(&mut failure_heap, &mut failed_queue, BFS_Q_FREE, 0),
                Ok(0)
            );
            assert_eq!(failed_queue, BFS_Q::default());
            assert_eq!(failure_heap.live_allocation_count(), 0);
        }

        let mut clear_heap = SourceHeap::default();
        let mut clear_queue = BFS_Q::default();
        assert_eq!(AllocBfsQueue(&mut clear_heap, &mut clear_queue, 2, 5), Ok(0));
        let retained_q = clear_queue.q;
        let retained_levels = clear_queue.nAtomLevel;
        let retained_sources = clear_queue.cSource;
        assert_eq!(AllocBfsQueue(&mut clear_heap, &mut clear_queue, BFS_Q_CLEAR, 0), Ok(0));
        assert_eq!(clear_queue, BFS_Q::default());
        assert_eq!(clear_heap.live_allocation_count(), 4);
        QueueDelete(&mut clear_heap, retained_q).unwrap();
        inchi_free(&mut clear_heap, retained_levels).unwrap();
        inchi_free(&mut clear_heap, retained_sources).unwrap();
        assert_eq!(clear_heap.live_allocation_count(), 0);
    }
}
