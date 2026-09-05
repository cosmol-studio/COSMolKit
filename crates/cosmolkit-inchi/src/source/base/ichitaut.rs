use crate::source::base::ichiqueu::{
    bIsCenterPointStrict, nGet12TautIn5MembAltRing, nGet14TautIn5MembAltRing,
    nGet14TautIn7MembAltRing, nGet15TautIn6MembAltRing, nGet15TautInAltPath,
};
use crate::source::base::{
    ichi_bns::{
        AddCGroups2BnStruct, AddTGroups2BnStruct, ReInitBnStruct, bExistsAltPath,
        bExistsAnyAltPath, bHasAcidicHydrogen, bHasAcidicMinus, bHasOtherExchangableH,
        bIsHardRemHCandidate,
    },
    ichisort::insertions_sort,
    util::{
        get_el_valence, get_endpoint_valence, get_endpoint_valence_KET, inchi_calloc, inchi_free,
        inchi_realloc,
    },
};
use crate::source_types::{
    ALT_PATH_MODE_4_SALT, ALT_PATH_MODE_ADD2H_CHG, ALT_PATH_MODE_ADD2H_TST, ALT_PATH_MODE_CHARGE,
    ALT_PATH_MODE_REM2H_CHG, ALT_PATH_MODE_REM2H_TST, ALT_PATH_MODE_TAUTOM,
    ALT_PATH_MODE_TAUTOM_KET, ALT_PATH_MODE_TAUTOM_PT_06_00, ALT_PATH_MODE_TAUTOM_PT_13_00,
    ALT_PATH_MODE_TAUTOM_PT_16_00, ALT_PATH_MODE_TAUTOM_PT_18_00, ALT_PATH_MODE_TAUTOM_PT_22_00,
    ALT_PATH_MODE_TAUTOM_PT_39_00, ALWAYS_ADD_TG_ON_THE_FLY, AT_FLAG_ISO_H_POINT, AT_NUMB, AT_RANK,
    ATT_ACIDIC_CO, BN_DATA, BN_STRUCT, BNS_CPOINT_ERR, BNS_ERR, BNS_MAX_ERR_VALUE, BNS_OUT_OF_RAM,
    BNS_PROGRAM_ERR, BNS_RADICAL_ERR, BNS_VERT_EDGE_OVFL, BOND_ALT_13, BOND_ALT12NS, BOND_ALTERN,
    BOND_DOUBLE, BOND_MARK_ALL, BOND_SINGLE, BOND_TAUTOM, BOND_TRIPLE, BOND_TYPE_MASK, C_CANDIDATE,
    C_GROUP, C_GROUP_INFO, CANON_GLOBALS, CT_OUT_OF_RAM, CT_TAUCOUNT_ERR, DFS_PATH, ENDPOINT_INFO,
    FLAG_FORCE_SALT_TAUT, FLAG_NORM_CONSIDER_TAUT, KETO_ENOL_TAUT, MAX_ATOMS, MAXVAL,
    NUM_H_ISOTOPES, RADICAL_SINGLET, S_CANDIDATE, S_CHAR, S_GROUP_INFO, SALT_ACCEPTOR, SALT_DONOR,
    SALT_DONOR_ALL, SALT_DONOR_H, SALT_DONOR_H2, SALT_DONOR_Neg, SALT_DONOR_Neg2, SALT_SELECTED,
    SALT_p_ACCEPTOR, SALT_p_DONOR, SourceHeap, SourceHeapError, SourceMutPointer, T_BONDPOS,
    T_ENDPOINT, T_GROUP, T_GROUP_HDR_LEN, T_GROUP_INFO, T_NUM_ISOTOPIC, T_NUM_NO_ISOTOPIC,
    TAUT_PT_06_00, TAUT_PT_13_00, TAUT_PT_16_00, TAUT_PT_18_00, TAUT_PT_22_00, TAUT_PT_39_00,
    TG_FLAG_1_5_TAUT, TG_FLAG_ALLOW_NO_NEGTV_O, TG_FLAG_ALLOW_NO_NEGTV_O_DONE,
    TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE,
    TG_FLAG_FOUND_SALT_CHARGES_DONE, TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_MOVE_POS_CHARGES,
    TG_FLAG_PT_06_00, TG_FLAG_PT_13_00, TG_FLAG_PT_16_00, TG_FLAG_PT_18_00, TG_FLAG_PT_22_00,
    TG_FLAG_PT_39_00, TG_FLAG_TEST_TAUT__ATOMS, clock_t, inp_ATOM,
    local_ichitaut::{
        C_SUBTYPE_CHARGED_H_ACCEPT, C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR, C_SUBTYPE_CHARGED_H_DONOR,
        C_SUBTYPE_CHARGED_NON_TAUT, C_SUBTYPE_CHARGED_p_DONOR, C_SUBTYPE_H_ACCEPT,
        C_SUBTYPE_H_DONOR, C_SUBTYPE_NEUTRAL, C_SUBTYPE_NEUTRAL_H_ACCEPT,
        C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT, C_SUBTYPE_NEUTRAL_H_DONOR, C_SUBTYPE_NEUTRAL_NON_TAUT,
        CHARGE_TYPE, MAX_STACK_ARRAY_LEN,
    },
    sp_ATOM, tagTG_NumDA_TG_NUM_DA, tagTG_NumDA_TG_Num_aH, tagTG_NumDA_TG_Num_aM,
    tagTG_NumDA_TG_Num_aO, tagTG_NumDA_TG_Num_dH, tagTG_NumDA_TG_Num_dM, tagTG_NumDA_TG_Num_dO,
};
use std::mem::MaybeUninit;

const CTYPE: [CHARGE_TYPE; 6] = [
    CHARGE_TYPE {
        elname: [b'N' as i8, 0, 0],
        charge: 1,
        neutral_valence: 3,
        neutral_bonds_valence: 3,
        cChangeValence: 1,
        cChargeType: 0,
        num_bonds: 0,
    },
    CHARGE_TYPE {
        elname: [b'P' as i8, 0, 0],
        charge: 1,
        neutral_valence: 3,
        neutral_bonds_valence: 3,
        cChangeValence: 1,
        cChargeType: 1,
        num_bonds: 0,
    },
    CHARGE_TYPE {
        elname: [b'O' as i8, 0, 0],
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        cChangeValence: 1,
        cChargeType: 2,
        num_bonds: 2,
    },
    CHARGE_TYPE {
        elname: [b'S' as i8, 0, 0],
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        cChangeValence: 1,
        cChargeType: 3,
        num_bonds: 2,
    },
    CHARGE_TYPE {
        elname: [b'S' as i8, b'e' as i8, 0],
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        cChangeValence: 1,
        cChargeType: 4,
        num_bonds: 2,
    },
    CHARGE_TYPE {
        elname: [b'T' as i8, b'e' as i8, 0],
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        cChangeValence: 1,
        cChargeType: 5,
        num_bonds: 2,
    },
];

fn source_strcmp_zero(left: &[i8], right: &[i8]) -> i32 {
    let mut index = 0_usize;
    loop {
        let left_byte = left.get(index).copied().unwrap_or(0);
        let right_byte = right.get(index).copied().unwrap_or(0);
        if left_byte != right_byte {
            return i32::from(left_byte) - i32::from(right_byte);
        }
        if left_byte == 0 {
            return 0;
        }
        index += 1;
    }
}

fn add_assign_at_rank(slot: &mut AT_RANK, delta: i32) {
    *slot = (i32::from(*slot) + delta) as AT_RANK;
}

#[allow(non_snake_case)]
pub(crate) fn AddAtom2num(
    num: &mut [AT_RANK],
    atom: &[inp_ATOM],
    at_no: i32,
    bSubtract: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:211 AddAtom2num
    // INCHI✔️✔️: int AddAtom2num( AT_RANK num[], inp_ATOM *atom, int at_no, int bSubtract )
    // INCHI✔️✔️: {  /*  bSubtract: 0=> add, 1=>subtract, 2=> fill */
    // INCHI✔️✔️:     inp_ATOM *at = atom + at_no;
    // INCHI✔️✔️:     int       k;
    // INCHI✔️✔️:     int       nMobile = ( at->charge == -1 );
    // INCHI✔️✔️:     if (bSubtract == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* 1: subtract */
    // INCHI✔️✔️:         num[1] -= nMobile;
    // INCHI✔️✔️:         nMobile += at->num_H;
    // INCHI✔️✔️:         num[0] -= nMobile;
    // INCHI✔️✔️:         for (k = 0; k < T_NUM_ISOTOPIC; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /*  T (3H isotope) first because it has higher weight */
    // INCHI✔️✔️:             num[T_NUM_NO_ISOTOPIC + k] -= at->num_iso_H[NUM_H_ISOTOPES - k - 1];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (bSubtract == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* fill */
    // INCHI✔️✔️:             memset( num, 0, ( T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC ) * sizeof( num[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* else (0): add */
    // INCHI✔️✔️:         num[1] += nMobile;
    // INCHI✔️✔️:         nMobile += at->num_H;
    // INCHI✔️✔️:         num[0] += nMobile;
    // INCHI✔️✔️:         for (k = 0; k < T_NUM_ISOTOPIC; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /*  T (3H isotope) first because it has higher weight */
    // INCHI✔️✔️:             num[T_NUM_NO_ISOTOPIC + k] += at->num_iso_H[NUM_H_ISOTOPES - k - 1];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return nMobile;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddAtom2num
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddAtom2num
    // INCHI✔️✔️: #define T_NUM_NO_ISOTOPIC 2
    // INCHI✔️✔️: #define T_NUM_ISOTOPIC NUM_H_ISOTOPES
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // INCHI✔️✔️: typedef unsigned short AT_RANK;
    // END INCHI ACTIVE MACRO CONFIGURATION: AddAtom2num

    let at_index = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let at = atom
        .get(at_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let needed = (T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC) as usize;
    if num.len() < needed {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let mut n_mobile = i32::from(at.charge == -1);
    if bSubtract == 1 {
        add_assign_at_rank(&mut num[1], -n_mobile);
        n_mobile += i32::from(at.num_H);
        add_assign_at_rank(&mut num[0], -n_mobile);
        for k in 0..T_NUM_ISOTOPIC as usize {
            let isotope_index = NUM_H_ISOTOPES as usize - k - 1;
            add_assign_at_rank(
                &mut num[T_NUM_NO_ISOTOPIC as usize + k],
                -i32::from(at.num_iso_H[isotope_index]),
            );
        }
    } else {
        if bSubtract == 2 {
            num[..needed].fill(0);
        }
        add_assign_at_rank(&mut num[1], n_mobile);
        n_mobile += i32::from(at.num_H);
        add_assign_at_rank(&mut num[0], n_mobile);
        for k in 0..T_NUM_ISOTOPIC as usize {
            let isotope_index = NUM_H_ISOTOPES as usize - k - 1;
            add_assign_at_rank(
                &mut num[T_NUM_NO_ISOTOPIC as usize + k],
                i32::from(at.num_iso_H[isotope_index]),
            );
        }
    }

    Ok(n_mobile)
}

#[allow(non_snake_case)]
pub(crate) fn AddAtom2DA(
    num_da: &mut [AT_RANK],
    atom: &[inp_ATOM],
    at_no: i32,
    bSubtract: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:250 AddAtom2DA
    // INCHI✔️✔️: void AddAtom2DA( AT_RANK num_DA[], inp_ATOM *atom, int at_no, int bSubtract )
    // INCHI✔️✔️: {   /*  bSubtract: 0=> add, 1=>subtract, 2=> fill */
    // INCHI✔️✔️:     inp_ATOM *at = atom + at_no;
    // INCHI✔️✔️:     int       nDelta, nAcidic_O;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at->charge < -1 || (at->charge == 1 && !at->c_point) || at->charge > 1) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nDelta = ( bSubtract == 1 ) ? -1 : 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* "Acidic" O, S, Se, Te recognition */
    // INCHI✔️✔️:     if (at->at_type & ATT_ACIDIC_CO)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nAcidic_O = nDelta;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nAcidic_O = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bSubtract == 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* 2: fill, otherwise add */
    // INCHI✔️✔️:         memset( num_DA, 0, TG_NUM_DA * sizeof( num_DA[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if ((at->charge <= 0 && at->valence == at->chem_bonds_valence) ||
    // INCHI✔️✔️:          /* neutral or negative donor */
    // INCHI✔️✔️:          (at->charge > 0 && at->valence + 1 == at->chem_bonds_valence)
    // INCHI✔️✔️:          /* positively charged donor */
    // INCHI✔️✔️:          ) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (at->charge < 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             num_DA[TG_Num_dM] += nDelta;
    // INCHI✔️✔️:             num_DA[TG_Num_dO] += nAcidic_O;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at->num_H)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_DA[TG_Num_dH] += nDelta;
    // INCHI✔️✔️:                 num_DA[TG_Num_dO] += nAcidic_O;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if ((at->charge <= 0 && at->valence + 1 == at->chem_bonds_valence) ||
    // INCHI✔️✔️:              (at->charge  > 0 && at->valence + 2 == at->chem_bonds_valence)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* acceptor */
    // INCHI✔️✔️:             if (at->charge < 0)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_DA[TG_Num_aM] += nDelta;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (at->num_H)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     num_DA[TG_Num_aH] += nDelta;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     num_DA[TG_Num_aO] += nAcidic_O; /* acidic O-acceptor has no H or charge */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddAtom2DA
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddAtom2DA
    // INCHI✔️✔️: #define ATT_ACIDIC_CO 0x0001
    // INCHI✔️✔️: #define TG_NUM_DA 6
    // INCHI✔️✔️: #define TG_Num_dH 0
    // INCHI✔️✔️: #define TG_Num_dM 1
    // INCHI✔️✔️: #define TG_Num_aH 2
    // INCHI✔️✔️: #define TG_Num_aM 3
    // INCHI✔️✔️: #define TG_Num_dO 4
    // INCHI✔️✔️: #define TG_Num_aO 5
    // END INCHI ACTIVE MACRO CONFIGURATION: AddAtom2DA

    let at_index = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let at = atom
        .get(at_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if at.charge < -1 || (at.charge == 1 && at.c_point == 0) || at.charge > 1 {
        return Ok(());
    }
    let n_delta = if bSubtract == 1 { -1 } else { 1 };
    let n_acidic_o = if (u32::from(at.at_type) & ATT_ACIDIC_CO) != 0 {
        n_delta
    } else {
        0
    };
    if num_da.len() < tagTG_NumDA_TG_NUM_DA as usize {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    if bSubtract == 2 {
        num_da[..tagTG_NumDA_TG_NUM_DA as usize].fill(0);
    }
    if (at.charge <= 0 && at.valence == at.chem_bonds_valence)
        || (at.charge > 0 && at.valence + 1 == at.chem_bonds_valence)
    {
        if at.charge < 0 {
            add_assign_at_rank(&mut num_da[tagTG_NumDA_TG_Num_dM as usize], n_delta);
            add_assign_at_rank(&mut num_da[tagTG_NumDA_TG_Num_dO as usize], n_acidic_o);
        } else if at.num_H != 0 {
            add_assign_at_rank(&mut num_da[tagTG_NumDA_TG_Num_dH as usize], n_delta);
            add_assign_at_rank(&mut num_da[tagTG_NumDA_TG_Num_dO as usize], n_acidic_o);
        }
    } else if (at.charge <= 0 && at.valence + 1 == at.chem_bonds_valence)
        || (at.charge > 0 && at.valence + 2 == at.chem_bonds_valence)
    {
        if at.charge < 0 {
            add_assign_at_rank(&mut num_da[tagTG_NumDA_TG_Num_aM as usize], n_delta);
        } else if at.num_H != 0 {
            add_assign_at_rank(&mut num_da[tagTG_NumDA_TG_Num_aH as usize], n_delta);
        } else {
            add_assign_at_rank(&mut num_da[tagTG_NumDA_TG_Num_aO as usize], n_acidic_o);
        }
    }

    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn AddEndPoint(
    pEndPoint: &mut T_ENDPOINT,
    at: &[inp_ATOM],
    iat: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:330 AddEndPoint
    // INCHI✔️✔️: int AddEndPoint( T_ENDPOINT *pEndPoint, inp_ATOM *at, int iat )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     pEndPoint->nAtomNumber = iat;
    // INCHI✔️✔️:     pEndPoint->nEquNumber = 0;
    // INCHI✔️✔️:     pEndPoint->nGroupNumber = at[iat].endpoint;
    // INCHI✔️✔️:     if (at[iat].endpoint)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* already an endpoint */
    // INCHI✔️✔️:         memset( pEndPoint->num, 0, sizeof( pEndPoint->num ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* not an endpoint yet, make it an endpoint */
    // INCHI✔️✔️:         AddAtom2num( pEndPoint->num, at, iat, 2 );  /* fill */
    // INCHI✔️✔️:         AddAtom2DA( pEndPoint->num_DA, at, iat, 2 );
    // INCHI✔️✔️:         /*
    // INCHI✔️✔️:         nMobile  = pEndPoint->num[1] = (at[iat].charge == -1);
    // INCHI✔️✔️:         nMobile  = pEndPoint->num[0] = at[iat].num_H + nMobile;
    // INCHI✔️✔️:         for ( k = 0; k < T_NUM_ISOTOPIC; k ++ ) {
    // INCHI✔️✔️:         pEndPoint->num[T_NUM_NO_ISOTOPIC+k] = at[iat].num_iso_H[NUM_H_ISOTOPES-k-1];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddEndPoint
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddEndPoint
    // INCHI✔️✔️: typedef unsigned short AT_NUMB;
    // INCHI✔️✔️: #define T_NUM_NO_ISOTOPIC 2
    // INCHI✔️✔️: #define T_NUM_ISOTOPIC NUM_H_ISOTOPES
    // INCHI✔️✔️: #define TG_NUM_DA 6
    // END INCHI ACTIVE MACRO CONFIGURATION: AddEndPoint

    let at_index = usize::try_from(iat).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = at
        .get(at_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    pEndPoint.nAtomNumber = iat as AT_NUMB;
    pEndPoint.nEquNumber = 0;
    pEndPoint.nGroupNumber = atom.endpoint;
    if atom.endpoint != 0 {
        pEndPoint.num.fill(0);
    } else {
        AddAtom2num(&mut pEndPoint.num, at, iat, 2)?;
        AddAtom2DA(&mut pEndPoint.num_DA, at, iat, 2)?;
    }

    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn comp_candidates(s1: &S_CANDIDATE, s2: &S_CANDIDATE) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2930 comp_candidates
    // INCHI✔️✔️: int comp_candidates( const void *a1, const void *a2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     const S_CANDIDATE *s1 = (const S_CANDIDATE *) a1;
    // INCHI✔️✔️:     const S_CANDIDATE *s2 = (const S_CANDIDATE *) a2;
    // INCHI✔️✔️:     int ret;
    // INCHI✔️✔️:     if (s1->type >= 0 /* enabled < */ && s2->type < 0 /* disabled */)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* enabled goes first */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (s1->type < 0 /* disabled > */ && s2->type >= 0 /* enabled */)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (s1->endpoint && !s2->endpoint)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* tautomeric goes first; only tautomeric may be disabled */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (!s1->endpoint && s2->endpoint)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1; /* tautomeric goes first; only tautomeric may be disabled */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (s1->endpoint && s2->endpoint && ( ret = (int) s1->endpoint - (int) s2->endpoint ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ret;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return (int) s1->atnumber - (int) s2->atnumber;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: comp_candidates
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: comp_candidates
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // INCHI✔️✔️: typedef unsigned short AT_NUMB;
    // END INCHI ACTIVE MACRO CONFIGURATION: comp_candidates

    if s1.type_ >= 0 && s2.type_ < 0 {
        return -1;
    }
    if s1.type_ < 0 && s2.type_ >= 0 {
        return 1;
    }
    if s1.endpoint != 0 && s2.endpoint == 0 {
        return -1;
    }
    if s1.endpoint == 0 && s2.endpoint != 0 {
        return 1;
    }
    if s1.endpoint != 0 && s2.endpoint != 0 {
        let ret = i32::from(s1.endpoint) - i32::from(s2.endpoint);
        if ret != 0 {
            return ret;
        }
    }
    i32::from(s1.atnumber) - i32::from(s2.atnumber)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MarkSaltChargeGroups2(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    s_group_info: Option<&mut S_GROUP_INFO>,
    t_group_info: Option<&mut T_GROUP_INFO>,
    mut c_group_info: Option<&mut C_GROUP_INFO>,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2961 MarkSaltChargeGroups2
    // INCHI✔️❌: int MarkSaltChargeGroups2( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                            inp_ATOM *at,
    // INCHI✔️❌:                            int num_atoms,
    // INCHI✔️❌:                            S_GROUP_INFO *s_group_info,
    // INCHI✔️❌:                            T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                            C_GROUP_INFO *c_group_info,
    // INCHI✔️❌:                            struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                            struct BalancedNetworkData *pBD )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* BNS_EDGE_FORBIDDEN_TEMP */
    // INCHI✔️❌: #define ALT_PATH_FOUND    (MAX_ATOMS+1)
    // INCHI✔️❌: #define NO_ENDPOINT       (MAX_ATOMS+2)  /* the two defines must be different */
    // INCHI✔️❌: #define DISABLE_CANDIDATE 10
    // INCHI✔️❌: #define cPAIR(a,b) cPair[a+b*nNumLeftCandidates]
    // INCHI✔️❌: #define ACCEPTOR_PAIR 1
    // INCHI✔️❌: #define DONOR_PAIR    2
    // INCHI✔️❌:
    // INCHI✔️❌:     int nNumOtherChanges = 0, nTotNumChanges = 0; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     S_CHAR      *cPair = NULL;
    // INCHI✔️❌:     T_ENDPOINT  *EndPoint = NULL;
    // INCHI✔️❌:     if (s_group_info && s_group_info->s_candidate && s_group_info->max_num_candidates > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i, j, i1, j1;
    // INCHI✔️❌:         S_CANDIDATE *s_candidate = s_group_info->s_candidate;
    // INCHI✔️❌:         int          nMaxNumCandidates = s_group_info->max_num_candidates;
    // INCHI✔️❌:         int          nNumCandidates = s_group_info->num_candidates;
    // INCHI✔️❌:         int          nNumOtherCandidates; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:         int          nNumLeftCandidates = 0;
    // INCHI✔️❌:         int          nNumMarkedCandidates = 0;
    // INCHI✔️❌:         int          s_type, s_subtype;
    // INCHI✔️❌:         int          ret, nDelta;
    // INCHI✔️❌:         int          bHardAddedRemovedProtons = t_group_info && ( t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT );
    // INCHI✔️❌:
    // INCHI✔️❌:         int          s_subtype_all = 0;
    // INCHI✔️❌:         int          nDonorPairs, nAcceptorPairs, nCurDonorPairs, nCurAcceptorPairs, bAlreadyTested;
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         ENDPOINT_INFO    eif;
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( IGNORE_TGROUP_WITHOUT_H == 1 )
    // INCHI✔️❌:         int          bTGroupHasNegativeChargesOnly = 1;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         /*return 0;*/ /* debug only */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumCandidates <= -2 || !t_group_info || !t_group_info->t_group)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*************************************************************************/
    // INCHI✔️❌:         /* Find all candidates including those with differen s_type (other type) */
    // INCHI✔️❌:         /*************************************************************************/
    // INCHI✔️❌:         for (i = 0, nNumCandidates = nNumOtherCandidates = 0; i < num_atoms; i++) /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (0 == ( s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype ) ) ||
    // INCHI✔️❌:                  /* -C=O or =C-OH, O = S, Se, Te */
    // INCHI✔️❌:                  1 == ( s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1/* bAccept_O*/ ) ) ||
    // INCHI✔️❌:                  /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above */
    // INCHI✔️❌:                  2 == ( s_type = GetOtherSaltType( at, i, &s_subtype ) ) ||
    // INCHI✔️❌:                  ( bHardAddedRemovedProtons && 4 == ( s_type = bIsHardRemHCandidate( at, i, &s_subtype ) ) )
    // INCHI✔️❌:                  /* >C-SH, >C-S(-); S=S,Se,Te */
    // INCHI✔️❌:                  )
    // INCHI✔️❌:             {
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (nNumCandidates >= nMaxNumCandidates)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 s_candidate[nNumCandidates].atnumber = i;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].type = s_type;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].subtype = s_subtype;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].endpoint = at[i].endpoint;
    // INCHI✔️❌:                 nNumCandidates++;
    // INCHI✔️❌:                 nNumOtherCandidates += ( 1 == s_type );
    // INCHI✔️❌:                 s_subtype_all |= s_subtype;
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumCandidates <= 1 ||  /* TG_FLAG_ALLOW_NO_NEGTV_O <=> CHARGED_SALTS_ONLY=0 */
    // INCHI✔️❌:              !( s_subtype_all & SALT_ACCEPTOR ) ||
    // INCHI✔️❌:              ( ( ( t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O ) ||
    // INCHI✔️❌:              ( t_group_info->bTautFlagsDone & TG_FLAG_FOUND_SALT_CHARGES_DONE ) ||
    // INCHI✔️❌:                  ( t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT ) ) ?
    // INCHI✔️❌:                !( s_subtype_all & ( SALT_DONOR ) ) :
    // INCHI✔️❌:                ( !( s_subtype_all & SALT_DONOR_Neg ) || nNumOtherCandidates == nNumCandidates ) )
    // INCHI✔️❌:              )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             s_group_info->num_candidates = 0; /* no candidate exists */
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!( s_subtype_all & ( SALT_DONOR_Neg ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t_group_info->bTautFlagsDone |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /************************************************************************************/
    // INCHI✔️❌:         /* Mark redundant candidates so that only one candidate from one t-group is left in */
    // INCHI✔️❌:         /************************************************************************************/
    // INCHI✔️❌:         for (i = 0; i < nNumCandidates; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #45059 */
    // INCHI✔️❌:             if ((nNumCandidates < nMaxNumCandidates) && (2 == s_candidate[nNumCandidates].type))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 s_candidate[i].type -= DISABLE_CANDIDATE; /* disable >C-SH candidates */
    // INCHI✔️❌:                 nNumLeftCandidates++; /* count rejected */
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (s_candidate[i].endpoint)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (j = i - 1; 0 <= j; j--)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (s_candidate[i].endpoint == s_candidate[j].endpoint)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         s_candidate[i].type -= DISABLE_CANDIDATE; /* disable subsequent redundant */
    // INCHI✔️❌:                         nNumLeftCandidates++; /* count rejected */
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         nNumLeftCandidates = nNumCandidates - nNumLeftCandidates; /* subtract num. rejected from the total */
    // INCHI✔️❌:         s_group_info->num_candidates = 0; /* reinit next time */
    // INCHI✔️❌:                                           /*********************************************************************/
    // INCHI✔️❌:                                           /* reorder so that all disabled are at the end, tautomeric are first */
    // INCHI✔️❌:                                           /*********************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:         qsort( s_candidate, nNumCandidates, sizeof( s_candidate[0] ), comp_candidates );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumLeftCandidates) /* djb-rwth: avoiding zero-length allocation */
    // INCHI✔️❌:             cPair = (S_CHAR *) inchi_calloc( (long long)nNumLeftCandidates * (long long)nNumLeftCandidates, sizeof( cPair[0] ) ); /* djb-rwth: cast operators added */
    // INCHI✔️❌:         if (!cPair)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*printf("BNS_OUT_OF_RAM-6\n");*/
    // INCHI✔️❌:             nTotNumChanges = BNS_OUT_OF_RAM;
    // INCHI✔️❌:             goto quick_exit;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nDonorPairs = nAcceptorPairs = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         /**********************************************************************/
    // INCHI✔️❌:         /* Find whether we have at least one donor pair and one acceptor pair */
    // INCHI✔️❌:         /**********************************************************************/
    // INCHI✔️❌:         for (i = 0; i < nNumLeftCandidates; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nCurDonorPairs = nCurAcceptorPairs = 0;
    // INCHI✔️❌:             for (j = 0; j <= i; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i == j && !s_candidate[i].endpoint)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;  /* same non-taut atom. However, success for i==j means     *
    // INCHI✔️❌:                                * that the whole tautomeric group may donate or accept 2H */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* check for acceptor pair */
    // INCHI✔️❌:                 if (( s_candidate[i].subtype & SALT_ACCEPTOR ) && ( s_candidate[j].subtype & SALT_ACCEPTOR ) &&
    // INCHI✔️❌:                     ( ret = bExistsAltPath( pCG, pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
    // INCHI✔️❌:                                             s_candidate[j].atnumber, ALT_PATH_MODE_ADD2H_TST ) ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nTotNumChanges = ret;
    // INCHI✔️❌:                         goto quick_exit;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (ret & 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nDelta = ( ret & ~3 ) >> 2;
    // INCHI✔️❌:                         /*nNumChanges += (ret & 2);*/
    // INCHI✔️❌:                         if (nDelta)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* alt path unleashed previously localized radicals and they annihilated */
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             nTotNumChanges = BNS_RADICAL_ERR;
    // INCHI✔️❌:                             goto quick_exit;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         cPAIR( i, j ) |= ACCEPTOR_PAIR; /* the result: mark the pair */
    // INCHI✔️❌:                                                         /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* Check for donor pair */
    // INCHI✔️❌:                 if (( s_candidate[i].subtype & SALT_DONOR ) && ( s_candidate[j].subtype & SALT_DONOR ) &&
    // INCHI✔️❌:                     ( ret = bExistsAltPath( pCG, pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
    // INCHI✔️❌:                                             s_candidate[j].atnumber, ALT_PATH_MODE_REM2H_TST ) ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nTotNumChanges = ret;
    // INCHI✔️❌:                         goto quick_exit;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (ret & 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nDelta = ( ret & ~3 ) >> 2;
    // INCHI✔️❌:                         /*nNumChanges += (ret & 2);*/
    // INCHI✔️❌:                         if (nDelta)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* Alt path unleashed previously localized radicals and they annihilated */
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             nTotNumChanges = BNS_RADICAL_ERR;
    // INCHI✔️❌:                             goto quick_exit;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         cPAIR( i, j ) |= DONOR_PAIR; /* the result: mark the pair */
    // INCHI✔️❌:                                                      /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* Since the results will be used later to change bonds, check only now */
    // INCHI✔️❌:                 /* when both results for (i,j) have been obtained. */
    // INCHI✔️❌:                 if (cPAIR( i, j ) & ACCEPTOR_PAIR)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nCurAcceptorPairs++;
    // INCHI✔️❌:                     if (nDonorPairs)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* Find donor pair (i1,j1) such that i!=i1, i!=j1, j!=i1, j!=j1 */
    // INCHI✔️❌:                         for (i1 = 0; i1 < i; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             for (j1 = 0; j1 <= i1; j1++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* Here always j1 < i && i1 < i therefore we do not compare i to i1 or j1 */
    // INCHI✔️❌:                                 if (j1 != j && i1 != j && ( cPAIR( i1, j1 ) & DONOR_PAIR ))
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     /* Both the donor and the acceptor pairs have been found */
    // INCHI✔️❌:                                     goto bFound2Pairs;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (cPAIR( i, j ) & DONOR_PAIR)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nCurDonorPairs++;
    // INCHI✔️❌:                     if (nAcceptorPairs)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* Find acceptor pair (i1,j1) such that i!=i1, i!=j1, j!=i1, j!=j1 */
    // INCHI✔️❌:                         for (i1 = 0; i1 < i; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             for (j1 = 0; j1 <= i1; j1++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* Here always j1 < i && i1 < i therefore we do not compare i to i1 or j1 */
    // INCHI✔️❌:                                 if (j1 != j && i1 != j && ( cPAIR( i1, j1 ) & ACCEPTOR_PAIR ))
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     /* Both the donor and the acceptor pairs have been found */
    // INCHI✔️❌:                                     goto bFound2Pairs;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nDonorPairs += nCurDonorPairs;
    // INCHI✔️❌:             nAcceptorPairs += nCurAcceptorPairs;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Nothing has been found */
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         inchi_free( cPair );
    // INCHI✔️❌:         cPair = NULL;
    // INCHI✔️❌:         goto quick_exit;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Both the donor and the acceptor pairs have been found */
    // INCHI✔️❌:     bFound2Pairs:
    // INCHI✔️❌:         /* first, try already found pairs */
    // INCHI✔️❌:         i1 = i;
    // INCHI✔️❌:         j1 = j;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Find all possible donor and acceptor pairs */
    // INCHI✔️❌:         nNumMarkedCandidates = 0;
    // INCHI✔️❌:         for (i = 0; i < nNumLeftCandidates; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nCurDonorPairs = nCurAcceptorPairs = 0;
    // INCHI✔️❌:             for (j = 0; j <= i; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bAlreadyTested = ( i < i1 || (i == i1 && j <= j1 )); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 if ((bAlreadyTested && ( cPAIR( i, j ) & ACCEPTOR_PAIR )) || !bAlreadyTested) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* Checking for acceptor pair */
    // INCHI✔️❌:                     if (( s_candidate[i].subtype & SALT_ACCEPTOR ) && ( s_candidate[j].subtype & SALT_ACCEPTOR ) &&
    // INCHI✔️❌:                         ( ret = bExistsAltPath( pCG, pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
    // INCHI✔️❌:                                                 s_candidate[j].atnumber, ALT_PATH_MODE_ADD2H_CHG ) ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nTotNumChanges = ret;
    // INCHI✔️❌:                             goto quick_exit;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (ret & 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nDelta = ( ret & ~3 ) >> 2;
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             if (nDelta)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* Alt path unleashed previously localized radicals and they annihilated */
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 nTotNumChanges = BNS_RADICAL_ERR;
    // INCHI✔️❌:                                 goto quick_exit;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             cPAIR( i, j ) |= ACCEPTOR_PAIR;
    // INCHI✔️❌:                             /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
    // INCHI✔️❌:                             nCurAcceptorPairs += !bAlreadyTested;
    // INCHI✔️❌:                             if (!( s_candidate[i].subtype & SALT_SELECTED ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 s_candidate[i].subtype |= SALT_SELECTED;
    // INCHI✔️❌:                                 nNumMarkedCandidates++;
    // INCHI✔️❌:                                 if (!s_candidate[i].endpoint && s_candidate[i].type)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     nNumOtherChanges++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!( s_candidate[j].subtype & SALT_SELECTED ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 s_candidate[j].subtype |= SALT_SELECTED;
    // INCHI✔️❌:                                 nNumMarkedCandidates++;
    // INCHI✔️❌:                                 if (!s_candidate[j].endpoint && s_candidate[j].type)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     nNumOtherChanges++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if ((bAlreadyTested && ( cPAIR( i, j ) & DONOR_PAIR )) || !bAlreadyTested) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* Checking for donor pair */
    // INCHI✔️❌:                     if (( s_candidate[i].subtype & SALT_DONOR ) && ( s_candidate[j].subtype & SALT_DONOR ) &&
    // INCHI✔️❌:                         ( ret = bExistsAltPath( pCG, pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
    // INCHI✔️❌:                                                 s_candidate[j].atnumber, ALT_PATH_MODE_REM2H_CHG ) ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nTotNumChanges = ret;
    // INCHI✔️❌:                             goto quick_exit;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (ret & 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nDelta = ( ret & ~3 ) >> 2;
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             if (nDelta)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* Alt path unleashed previously localized radicals and they annihilated */
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 nTotNumChanges = BNS_RADICAL_ERR;
    // INCHI✔️❌:                                 goto quick_exit;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             cPAIR( i, j ) |= DONOR_PAIR;
    // INCHI✔️❌:                             /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
    // INCHI✔️❌:                             nCurDonorPairs += !bAlreadyTested;
    // INCHI✔️❌:                             if (!( s_candidate[i].subtype & SALT_SELECTED ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 s_candidate[i].subtype |= SALT_SELECTED;
    // INCHI✔️❌:                                 nNumMarkedCandidates++;
    // INCHI✔️❌:                                 if (!s_candidate[i].endpoint && s_candidate[i].type)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     nNumOtherChanges++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!( s_candidate[j].subtype & SALT_SELECTED ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 s_candidate[j].subtype |= SALT_SELECTED;
    // INCHI✔️❌:                                 nNumMarkedCandidates++;
    // INCHI✔️❌:                                 if (!s_candidate[j].endpoint && s_candidate[j].type)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     nNumOtherChanges++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nDonorPairs += nCurDonorPairs;
    // INCHI✔️❌:             nAcceptorPairs += nCurAcceptorPairs;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         inchi_free( cPair );
    // INCHI✔️❌:         cPair = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumMarkedCandidates)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             EndPoint = (T_ENDPOINT *) inchi_calloc( nNumMarkedCandidates, sizeof( EndPoint[0] ) );
    // INCHI✔️❌:             if (!EndPoint)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*printf("BNS_OUT_OF_RAM-7\n");*/
    // INCHI✔️❌:                 nTotNumChanges = BNS_OUT_OF_RAM;
    // INCHI✔️❌:                 goto quick_exit;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (i = 0, j = 0; i < nNumLeftCandidates; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (s_candidate[i].subtype & SALT_SELECTED)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     s_candidate[i].subtype ^= SALT_SELECTED; /* remove the flag */
    // INCHI✔️❌:                     if (j < nNumMarkedCandidates)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         i1 = s_candidate[i].atnumber; /* save a representative of the t-group to be created */
    // INCHI✔️❌:                         AddEndPoint( EndPoint + j, at, i1 );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     j++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (j != nNumMarkedCandidates)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nTotNumChanges = BNS_PROGRAM_ERR;
    // INCHI✔️❌:                 goto quick_exit;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* Merge all marked atoms and their t-groups into one t-group */
    // INCHI✔️❌:             ret = RegisterEndPoints( pCG, t_group_info, EndPoint,
    // INCHI✔️❌:                                      nNumMarkedCandidates, at,
    // INCHI✔️❌:                                      num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (ret == -1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = BNS_PROGRAM_ERR;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (ret < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nTotNumChanges = ret;
    // INCHI✔️❌:                 goto quick_exit;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nTotNumChanges += ( ret > 0 );
    // INCHI✔️❌:             inchi_free( EndPoint );
    // INCHI✔️❌:             EndPoint = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:             if (nNumMarkedCandidates)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = nNumLeftCandidates; i < nNumCandidates; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     s_candidate[i].type += DISABLE_CANDIDATE;
    // INCHI✔️❌:                     j1 = s_candidate[i].atnumber;
    // INCHI✔️❌:                     if (at[j1].endpoint == at[i1].endpoint)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (!s_candidate[i].endpoint && s_candidate[i].type)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nNumOtherChanges++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = nNumLeftCandidates; i < nNumCandidates; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     s_candidate[i].type += DISABLE_CANDIDATE;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* Find whether the new t-group have any movable H */
    // INCHI✔️❌:             for (i = 0, bTGroupHasNegativeChargesOnly = 0; i < t_group_info->num_t_groups; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (t_group_info->t_group[i].nGroupNumber == at[i1].endpoint &&
    // INCHI✔️❌:                      t_group_info->t_group[i].num[0] == t_group_info->t_group[i].num[1])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bTGroupHasNegativeChargesOnly = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nTotNumChanges = ( nTotNumChanges > 0 );
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( IGNORE_TGROUP_WITHOUT_H == 1 )
    // INCHI✔️❌:         if (nTotNumChanges && bTGroupHasNegativeChargesOnly)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nTotNumChanges = 2;  /* means no moveable H has been affected */
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: quick_exit:
    // INCHI✔️❌:     if (nNumOtherChanges && nTotNumChanges == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nTotNumChanges = 5; /* not only acidic atoms merged */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (cPair)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( cPair );
    // INCHI✔️❌:         /*cPair = NULL;*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (EndPoint)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( EndPoint );
    // INCHI✔️❌:         /*EndPoint = NULL;*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nTotNumChanges; /* 0=>no changes, 1=>new salt tautomerism found, 2=>only new charge tautomerism found */
    // INCHI✔️❌: #undef ALT_PATH_FOUND
    // INCHI✔️❌: #undef NO_ENDPOINT
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkSaltChargeGroups2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MarkSaltChargeGroups2
    // INCHI✔️❌: #define IGNORE_TGROUP_WITHOUT_H 1
    // INCHI✔️❌: #define FLAG_FORCE_SALT_TAUT 50
    // INCHI✔️❌: #define TG_FLAG_ALLOW_NO_NEGTV_O 32
    // INCHI✔️❌: #define TG_FLAG_ALLOW_NO_NEGTV_O_DONE 64
    // INCHI✔️❌: #define TG_FLAG_FOUND_SALT_CHARGES_DONE 8192
    // INCHI✔️❌: #define SALT_DONOR 3
    // INCHI✔️❌: #define SALT_ACCEPTOR 4
    // INCHI✔️❌: #define SALT_SELECTED 32
    // INCHI✔️❌: #define ALT_PATH_MODE_REM2H_CHG 5
    // INCHI✔️❌: #define ALT_PATH_MODE_ADD2H_CHG 6
    // INCHI✔️❌: #define ALT_PATH_MODE_REM2H_TST 7
    // INCHI✔️❌: #define ALT_PATH_MODE_ADD2H_TST 8
    // INCHI✔️❌: #define BNS_OUT_OF_RAM -9998
    // INCHI✔️❌: #define BNS_PROGRAM_ERR -9997
    // INCHI✔️❌: #define BNS_VERT_EDGE_OVFL -9993
    // INCHI✔️❌: #define BNS_RADICAL_ERR -9988
    // INCHI✔️❌: #define IS_BNS_ERROR(X) (BNS_ERR <= (X) && (X) <= BNS_MAX_ERR_VALUE)
    // END INCHI ACTIVE MACRO CONFIGURATION: MarkSaltChargeGroups2

    const DISABLE_CANDIDATE: i32 = 10;
    const ACCEPTOR_PAIR: S_CHAR = 1;
    const DONOR_PAIR: S_CHAR = 2;

    fn allocation_failure(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
        )
    }

    let mut c_pair = None;
    let mut end_point = None;
    let result = (|| -> Result<i32, SourceHeapError> {
        let Some(s_group_info) = s_group_info else {
            return Ok(0);
        };
        if s_group_info.s_candidate.is_null() || s_group_info.max_num_candidates <= 0 {
            return Ok(0);
        }

        let s_candidate = s_group_info.s_candidate;
        let n_max_num_candidates = s_group_info.max_num_candidates;
        let mut n_num_candidates = s_group_info.num_candidates;
        let mut n_num_other_changes = 0_i32;
        let mut n_tot_num_changes = 0_i32;
        let mut n_num_left_candidates = 0_i32;
        let mut n_num_marked_candidates = 0_i32;
        let mut b_t_group_has_negative_charges_only = 1_i32;

        let Some(t_group_info) = t_group_info else {
            return Ok(0);
        };
        if n_num_candidates <= -2 || t_group_info.t_group.is_null() {
            return Ok(0);
        }

        let max_count = usize::try_from(n_max_num_candidates)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if heap.slice(s_candidate.as_const())?.len() < max_count {
            return Err(SourceHeapError::PointerOutOfBounds);
        }

        let b_hard_added_removed_protons =
            (t_group_info.tni.bNormalizationFlags & u64::from(FLAG_FORCE_SALT_TAUT)) != 0;
        let mut n_num_other_candidates = 0_i32;
        let mut s_subtype_all = 0_i32;
        n_num_candidates = 0;
        for i in 0..num_atoms {
            let mut s_subtype = 0_i32;
            let mut s_type = {
                let atoms = heap.slice(at.as_const())?;
                GetSaltChargeType(heap, atoms, i, Some(&*t_group_info), &mut s_subtype)?
            };
            let mut accepted = s_type == 0;
            if !accepted {
                s_type = {
                    let atoms = heap.slice(at.as_const())?;
                    GetOtherSaltChargeType(heap, atoms, i, Some(&*t_group_info), &mut s_subtype, 1)?
                };
                accepted = s_type == 1;
            }
            if !accepted {
                s_type = {
                    let atoms = heap.slice(at.as_const())?;
                    GetOtherSaltType(atoms, i, &mut s_subtype)?
                };
                accepted = s_type == 2;
            }
            if !accepted && b_hard_added_removed_protons {
                s_type = {
                    let atoms = heap.slice(at.as_const())?;
                    bIsHardRemHCandidate(atoms, i, &mut s_subtype)?
                };
                accepted = s_type == 4;
            }
            if !accepted {
                continue;
            }
            if n_num_candidates >= n_max_num_candidates {
                return Ok(BNS_VERT_EDGE_OVFL);
            }
            let atom_endpoint = heap
                .slice(at.as_const())?
                .get(i as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .endpoint;
            let candidate = heap
                .slice_mut(s_candidate)?
                .get_mut(n_num_candidates as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            candidate.atnumber = i as AT_NUMB;
            candidate.type_ = s_type as S_CHAR;
            candidate.subtype = s_subtype as S_CHAR;
            candidate.endpoint = atom_endpoint;
            n_num_candidates += 1;
            n_num_other_candidates += i32::from(s_type == 1);
            s_subtype_all |= s_subtype;
        }

        let permissive_donor_rule = (t_group_info.bTautFlags & u64::from(TG_FLAG_ALLOW_NO_NEGTV_O))
            != 0
            || (t_group_info.bTautFlagsDone & u64::from(TG_FLAG_FOUND_SALT_CHARGES_DONE)) != 0
            || (t_group_info.tni.bNormalizationFlags & u64::from(FLAG_FORCE_SALT_TAUT)) != 0;
        if n_num_candidates <= 1
            || (s_subtype_all & SALT_ACCEPTOR as i32) == 0
            || if permissive_donor_rule {
                (s_subtype_all & SALT_DONOR as i32) == 0
            } else {
                (s_subtype_all & SALT_DONOR_Neg as i32) == 0
                    || n_num_other_candidates == n_num_candidates
            }
        {
            s_group_info.num_candidates = 0;
            return Ok(0);
        }

        if (s_subtype_all & SALT_DONOR_Neg as i32) == 0 {
            t_group_info.bTautFlagsDone |= u64::from(TG_FLAG_ALLOW_NO_NEGTV_O_DONE);
        }

        let spare_type_is_two = n_num_candidates < n_max_num_candidates
            && heap.slice(s_candidate.as_const())?[n_num_candidates as usize].type_ == 2;
        for i in 0..n_num_candidates as usize {
            if spare_type_is_two {
                let candidate = &mut heap.slice_mut(s_candidate)?[i];
                candidate.type_ = (i32::from(candidate.type_) - DISABLE_CANDIDATE) as S_CHAR;
                n_num_left_candidates += 1;
                continue;
            }
            let endpoint = heap.slice(s_candidate.as_const())?[i].endpoint;
            if endpoint != 0 {
                for j in (0..i).rev() {
                    if endpoint == heap.slice(s_candidate.as_const())?[j].endpoint {
                        let candidate = &mut heap.slice_mut(s_candidate)?[i];
                        candidate.type_ =
                            (i32::from(candidate.type_) - DISABLE_CANDIDATE) as S_CHAR;
                        n_num_left_candidates += 1;
                        break;
                    }
                }
            }
        }

        n_num_left_candidates = n_num_candidates - n_num_left_candidates;
        s_group_info.num_candidates = 0;
        let candidate_count =
            usize::try_from(n_num_candidates).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut candidates = heap.slice(s_candidate.as_const())?[..candidate_count].to_vec();
        candidates.sort_by(|left, right| comp_candidates(left, right).cmp(&0));
        heap.slice_mut(s_candidate)?[..candidate_count].clone_from_slice(&candidates);

        if n_num_left_candidates != 0 {
            let pair_count = i64::from(n_num_left_candidates)
                .checked_mul(i64::from(n_num_left_candidates))
                .ok_or(SourceHeapError::AllocationSizeOverflow)?;
            match inchi_calloc::<S_CHAR>(heap, pair_count as u64, 1) {
                Ok(pointer) => c_pair = Some(pointer),
                Err(error) if allocation_failure(error) => return Ok(BNS_OUT_OF_RAM),
                Err(error) => return Err(error),
            }
        }
        let Some(c_pair_pointer) = c_pair else {
            return Ok(BNS_OUT_OF_RAM);
        };

        let pair_index = |i: i32, j: i32| -> usize { (i + j * n_num_left_candidates) as usize };
        let mut n_donor_pairs = 0_i32;
        let mut n_acceptor_pairs = 0_i32;
        let mut found_pair = None;
        'find_two_pairs: for i in 0..n_num_left_candidates {
            let mut n_cur_donor_pairs = 0_i32;
            let mut n_cur_acceptor_pairs = 0_i32;
            for j in 0..=i {
                if i == j && candidates[i as usize].endpoint == 0 {
                    continue;
                }
                if (i32::from(candidates[i as usize].subtype) & SALT_ACCEPTOR as i32) != 0
                    && (i32::from(candidates[j as usize].subtype) & SALT_ACCEPTOR as i32) != 0
                {
                    let ret = bExistsAltPath(
                        heap,
                        pCG,
                        pBNS,
                        pBD,
                        None,
                        at,
                        num_atoms,
                        i32::from(candidates[i as usize].atnumber),
                        i32::from(candidates[j as usize].atnumber),
                        ALT_PATH_MODE_ADD2H_TST as i32,
                        clock_result,
                    )?;
                    if ret != 0 {
                        if ichitaut_is_bns_error(ret) {
                            return Ok(ret);
                        }
                        if (ret & 1) != 0 {
                            if ((ret & !3) >> 2) != 0 {
                                return Ok(BNS_RADICAL_ERR);
                            }
                            heap.slice_mut(c_pair_pointer)?[pair_index(i, j)] |= ACCEPTOR_PAIR;
                        }
                    }
                }
                if (i32::from(candidates[i as usize].subtype) & SALT_DONOR as i32) != 0
                    && (i32::from(candidates[j as usize].subtype) & SALT_DONOR as i32) != 0
                {
                    let ret = bExistsAltPath(
                        heap,
                        pCG,
                        pBNS,
                        pBD,
                        None,
                        at,
                        num_atoms,
                        i32::from(candidates[i as usize].atnumber),
                        i32::from(candidates[j as usize].atnumber),
                        ALT_PATH_MODE_REM2H_TST as i32,
                        clock_result,
                    )?;
                    if ret != 0 {
                        if ichitaut_is_bns_error(ret) {
                            return Ok(ret);
                        }
                        if (ret & 1) != 0 {
                            if ((ret & !3) >> 2) != 0 {
                                return Ok(BNS_RADICAL_ERR);
                            }
                            heap.slice_mut(c_pair_pointer)?[pair_index(i, j)] |= DONOR_PAIR;
                        }
                    }
                }
                let pair = heap.slice(c_pair_pointer.as_const())?[pair_index(i, j)];
                if (pair & ACCEPTOR_PAIR) != 0 {
                    n_cur_acceptor_pairs += 1;
                    if n_donor_pairs != 0 {
                        for i1 in 0..i {
                            for j1 in 0..=i1 {
                                if j1 != j
                                    && i1 != j
                                    && (heap.slice(c_pair_pointer.as_const())?[pair_index(i1, j1)]
                                        & DONOR_PAIR)
                                        != 0
                                {
                                    found_pair = Some((i, j));
                                    break 'find_two_pairs;
                                }
                            }
                        }
                    }
                }
                if (pair & DONOR_PAIR) != 0 {
                    n_cur_donor_pairs += 1;
                    if n_acceptor_pairs != 0 {
                        for i1 in 0..i {
                            for j1 in 0..=i1 {
                                if j1 != j
                                    && i1 != j
                                    && (heap.slice(c_pair_pointer.as_const())?[pair_index(i1, j1)]
                                        & ACCEPTOR_PAIR)
                                        != 0
                                {
                                    found_pair = Some((i, j));
                                    break 'find_two_pairs;
                                }
                            }
                        }
                    }
                }
            }
            n_donor_pairs += n_cur_donor_pairs;
            n_acceptor_pairs += n_cur_acceptor_pairs;
        }

        let Some((first_i, first_j)) = found_pair else {
            inchi_free(heap, c_pair_pointer)?;
            c_pair = None;
            return Ok(0);
        };

        for i in 0..n_num_left_candidates {
            let mut n_cur_donor_pairs = 0_i32;
            let mut n_cur_acceptor_pairs = 0_i32;
            for j in 0..=i {
                let already_tested = i < first_i || (i == first_i && j <= first_j);
                let pair = heap.slice(c_pair_pointer.as_const())?[pair_index(i, j)];
                if ((already_tested && (pair & ACCEPTOR_PAIR) != 0) || !already_tested)
                    && (i32::from(candidates[i as usize].subtype) & SALT_ACCEPTOR as i32) != 0
                    && (i32::from(candidates[j as usize].subtype) & SALT_ACCEPTOR as i32) != 0
                {
                    let ret = bExistsAltPath(
                        heap,
                        pCG,
                        pBNS,
                        pBD,
                        None,
                        at,
                        num_atoms,
                        i32::from(candidates[i as usize].atnumber),
                        i32::from(candidates[j as usize].atnumber),
                        ALT_PATH_MODE_ADD2H_CHG as i32,
                        clock_result,
                    )?;
                    if ret != 0 {
                        if ichitaut_is_bns_error(ret) {
                            return Ok(ret);
                        }
                        if (ret & 1) != 0 {
                            if ((ret & !3) >> 2) != 0 {
                                return Ok(BNS_RADICAL_ERR);
                            }
                            heap.slice_mut(c_pair_pointer)?[pair_index(i, j)] |= ACCEPTOR_PAIR;
                            n_cur_acceptor_pairs += i32::from(!already_tested);
                            for candidate_index in [i as usize, j as usize] {
                                if (i32::from(candidates[candidate_index].subtype)
                                    & SALT_SELECTED as i32)
                                    == 0
                                {
                                    candidates[candidate_index].subtype =
                                        (i32::from(candidates[candidate_index].subtype)
                                            | SALT_SELECTED as i32)
                                            as S_CHAR;
                                    heap.slice_mut(s_candidate)?[candidate_index].subtype =
                                        candidates[candidate_index].subtype;
                                    n_num_marked_candidates += 1;
                                    if candidates[candidate_index].endpoint == 0
                                        && candidates[candidate_index].type_ != 0
                                    {
                                        n_num_other_changes += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                let pair = heap.slice(c_pair_pointer.as_const())?[pair_index(i, j)];
                if ((already_tested && (pair & DONOR_PAIR) != 0) || !already_tested)
                    && (i32::from(candidates[i as usize].subtype) & SALT_DONOR as i32) != 0
                    && (i32::from(candidates[j as usize].subtype) & SALT_DONOR as i32) != 0
                {
                    let ret = bExistsAltPath(
                        heap,
                        pCG,
                        pBNS,
                        pBD,
                        None,
                        at,
                        num_atoms,
                        i32::from(candidates[i as usize].atnumber),
                        i32::from(candidates[j as usize].atnumber),
                        ALT_PATH_MODE_REM2H_CHG as i32,
                        clock_result,
                    )?;
                    if ret != 0 {
                        if ichitaut_is_bns_error(ret) {
                            return Ok(ret);
                        }
                        if (ret & 1) != 0 {
                            if ((ret & !3) >> 2) != 0 {
                                return Ok(BNS_RADICAL_ERR);
                            }
                            heap.slice_mut(c_pair_pointer)?[pair_index(i, j)] |= DONOR_PAIR;
                            n_cur_donor_pairs += i32::from(!already_tested);
                            for candidate_index in [i as usize, j as usize] {
                                if (i32::from(candidates[candidate_index].subtype)
                                    & SALT_SELECTED as i32)
                                    == 0
                                {
                                    candidates[candidate_index].subtype =
                                        (i32::from(candidates[candidate_index].subtype)
                                            | SALT_SELECTED as i32)
                                            as S_CHAR;
                                    heap.slice_mut(s_candidate)?[candidate_index].subtype =
                                        candidates[candidate_index].subtype;
                                    n_num_marked_candidates += 1;
                                    if candidates[candidate_index].endpoint == 0
                                        && candidates[candidate_index].type_ != 0
                                    {
                                        n_num_other_changes += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            n_donor_pairs += n_cur_donor_pairs;
            n_acceptor_pairs += n_cur_acceptor_pairs;
        }

        inchi_free(heap, c_pair_pointer)?;
        c_pair = None;

        if n_num_marked_candidates != 0 {
            let endpoint_pointer = match inchi_calloc::<T_ENDPOINT>(
                heap,
                n_num_marked_candidates as u64,
                std::mem::size_of::<T_ENDPOINT>() as u64,
            ) {
                Ok(pointer) => pointer,
                Err(error) if allocation_failure(error) => return Ok(BNS_OUT_OF_RAM),
                Err(error) => return Err(error),
            };
            end_point = Some(endpoint_pointer);
            let atoms = heap.slice(at.as_const())?.to_vec();
            let mut endpoints = heap.slice(endpoint_pointer.as_const())?.to_vec();
            let mut j = 0_i32;
            let mut representative = 0_i32;
            for i in 0..n_num_left_candidates as usize {
                if (i32::from(candidates[i].subtype) & SALT_SELECTED as i32) != 0 {
                    candidates[i].subtype =
                        (i32::from(candidates[i].subtype) ^ SALT_SELECTED as i32) as S_CHAR;
                    heap.slice_mut(s_candidate)?[i].subtype = candidates[i].subtype;
                    if j < n_num_marked_candidates {
                        representative = i32::from(candidates[i].atnumber);
                        AddEndPoint(&mut endpoints[j as usize], &atoms, representative)?;
                    }
                    j += 1;
                }
            }
            if j != n_num_marked_candidates {
                return Ok(BNS_PROGRAM_ERR);
            }

            let mut ret = RegisterEndPoints(
                heap,
                pCG,
                t_group_info,
                &mut endpoints,
                n_num_marked_candidates,
                at,
                num_atoms,
                c_group_info.as_deref_mut(),
                Some(pBNS),
            )?;
            if ret == -1 {
                ret = BNS_PROGRAM_ERR;
            }
            if ret < 0 {
                return Ok(ret);
            }
            n_tot_num_changes += i32::from(ret > 0);
            inchi_free(heap, endpoint_pointer)?;
            end_point = None;

            let representative_endpoint = heap
                .slice(at.as_const())?
                .get(representative as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .endpoint;
            if n_num_marked_candidates != 0 {
                for i in n_num_left_candidates as usize..candidate_count {
                    candidates[i].type_ =
                        (i32::from(candidates[i].type_) + DISABLE_CANDIDATE) as S_CHAR;
                    heap.slice_mut(s_candidate)?[i].type_ = candidates[i].type_;
                    let atom_number = usize::from(candidates[i].atnumber);
                    let atom_endpoint = heap
                        .slice(at.as_const())?
                        .get(atom_number)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint;
                    if atom_endpoint == representative_endpoint
                        && candidates[i].endpoint == 0
                        && candidates[i].type_ != 0
                    {
                        n_num_other_changes += 1;
                    }
                }
            } else {
                for i in n_num_left_candidates as usize..candidate_count {
                    candidates[i].type_ =
                        (i32::from(candidates[i].type_) + DISABLE_CANDIDATE) as S_CHAR;
                    heap.slice_mut(s_candidate)?[i].type_ = candidates[i].type_;
                }
            }

            b_t_group_has_negative_charges_only = 0;
            let groups = heap.slice(t_group_info.t_group.as_const())?;
            let group_count = usize::try_from(t_group_info.num_t_groups)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            for group in groups.iter().take(group_count) {
                if group.nGroupNumber == representative_endpoint && group.num[0] == group.num[1] {
                    b_t_group_has_negative_charges_only = 1;
                    break;
                }
            }
        }

        n_tot_num_changes = i32::from(n_tot_num_changes > 0);
        if n_tot_num_changes != 0 && b_t_group_has_negative_charges_only != 0 {
            n_tot_num_changes = 2;
        }
        if n_num_other_changes != 0 && n_tot_num_changes == 1 {
            n_tot_num_changes = 5;
        }
        Ok(n_tot_num_changes)
    })();

    if let Some(pointer) = c_pair {
        inchi_free(heap, pointer)?;
    }
    if let Some(pointer) = end_point {
        inchi_free(heap, pointer)?;
    }
    result
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MarkSaltChargeGroups(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    s_group_info: Option<&mut S_GROUP_INFO>,
    t_group_info: Option<&mut T_GROUP_INFO>,
    mut c_group_info: Option<&mut C_GROUP_INFO>,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:3483 MarkSaltChargeGroups
    // INCHI✔️❌: int MarkSaltChargeGroups( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                           inp_ATOM *at,
    // INCHI✔️❌:                           int num_atoms,
    // INCHI✔️❌:                           S_GROUP_INFO *s_group_info,
    // INCHI✔️❌:                           T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                           C_GROUP_INFO *c_group_info,
    // INCHI✔️❌:                           struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                           struct BalancedNetworkData *pBD )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     int nNumChanges = 0, nTotNumChanges = 0;
    // INCHI✔️❌:     if (s_group_info && s_group_info->s_candidate && s_group_info->max_num_candidates > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i, i1, i2, j, j1, j2, jj, ii1, ii2, jj1, jj2, /*k,*/ num_tested;
    // INCHI✔️❌:         S_CANDIDATE *s_candidate = s_group_info->s_candidate;
    // INCHI✔️❌:         int          nMaxNumCandidates = s_group_info->max_num_candidates;
    // INCHI✔️❌:         int          nNumCandidates = s_group_info->num_candidates;
    // INCHI✔️❌:         int          nNumOtherCandidates; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         int          nNumPOnlyCandidates; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         int          s_type, s_subtype;
    // INCHI✔️❌:         int          ret, nDelta, /*nMobile,*/ err = 0;
    // INCHI✔️❌:         int          s_subtype_all = 0;
    // INCHI✔️❌:         int          nGroupNumber;
    // INCHI✔️❌:         T_ENDPOINT   EndPoint[2];
    // INCHI✔️❌: #if ( MAX_LOCAL_TGNUM > 0 )
    // INCHI✔️❌:         TGroupData   tgData[MAX_LOCAL_TGNUM];
    // INCHI✔️❌:         TGroupData   *ptgData = tgData;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         int cond1 = 0, cond2a = 0, cond2b = 0, cond2c = 0, cond2 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumCandidates <= -1 || !t_group_info || !t_group_info->t_group)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* count t-groups */
    // INCHI✔️❌:         for (i = 0, nGroupNumber = 0; i < t_group_info->num_t_groups; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nGroupNumber < t_group_info->t_group[i].nGroupNumber)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nGroupNumber = t_group_info->t_group[i].nGroupNumber; /* max. t-group number */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #if ( MAX_LOCAL_TGNUM > 0 )
    // INCHI✔️❌:         /* prepare memory */
    // INCHI✔️❌:         if (nGroupNumber >= MAX_LOCAL_TGNUM)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!( ptgData = (TGroupData*) inchi_calloc( nGroupNumber + 1, sizeof( TGroupData ) ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 err = BNS_OUT_OF_RAM;
    // INCHI✔️❌:                 goto quick_exit;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memset( ptgData, 0, sizeof( tgData ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ptgData[0].nGroupIndex = -1; /* data for non-tautomeric atoms */
    // INCHI✔️❌:         for (i = 0, nGroupNumber = 0; i < t_group_info->num_t_groups; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nGroupNumber = t_group_info->t_group[i].nGroupNumber)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ptgData[nGroupNumber].nGroupIndex = i;
    // INCHI✔️❌:                 ptgData[i].nGroupNumber = nGroupNumber;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         nNumCandidates = 0; /* always recalculate 2004-03-22 */
    // INCHI✔️❌:         num_tested = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumCandidates == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0, nNumCandidates = nNumOtherCandidates = nNumPOnlyCandidates = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (0 == ( s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype ) ) ||
    // INCHI✔️❌:                      /* -C=O or =C-OH, O = S, Se, Te */
    // INCHI✔️❌: #if ( INCL_NON_SALT_CANDIDATATES == 1 )
    // INCHI✔️❌:                      1 == ( s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1 ) ) ||
    // INCHI✔️❌:                      /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                      2 == ( s_type = GetOtherSaltType( at, i, &s_subtype ) )
    // INCHI✔️❌:                      /* >C-SH, >C-S(-); S=S,Se,Te */
    // INCHI✔️❌:                      )
    // INCHI✔️❌:                 {
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (nNumCandidates >= nMaxNumCandidates)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         err = BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                         goto quick_exit;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     s_candidate[nNumCandidates].atnumber = i;
    // INCHI✔️❌:                     s_candidate[nNumCandidates].type = s_type;
    // INCHI✔️❌:                     s_candidate[nNumCandidates].subtype = s_subtype;
    // INCHI✔️❌:                     s_candidate[nNumCandidates].endpoint = at[i].endpoint;
    // INCHI✔️❌:                     nNumCandidates++;
    // INCHI✔️❌:                     nNumOtherCandidates += ( 1 == s_type );
    // INCHI✔️❌:                     nNumPOnlyCandidates += ( 2 == s_type );
    // INCHI✔️❌:                     s_subtype_all |= s_subtype;
    // INCHI✔️❌:                     /*i1 = i;*/ /* save a representative of a tautomeric group */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             } /* for */
    // INCHI✔️❌:
    // INCHI✔️❌:               /* changes: TG_FLAG_ALLOW_NO_NEGTV_O replaced CHARGED_SALTS_ONLY==0 */
    // INCHI✔️❌:             cond1 = s_subtype_all & SALT_ACCEPTOR;
    // INCHI✔️❌:             cond2a = t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O;
    // INCHI✔️❌:             cond2b = t_group_info->bTautFlagsDone & TG_FLAG_FOUND_SALT_CHARGES_DONE;
    // INCHI✔️❌:             cond2c = t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT;
    // INCHI✔️❌:             if (cond2a || cond2b || cond2c)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 cond2 = !( s_subtype_all & ( SALT_DONOR_Neg | SALT_DONOR_H ) );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 cond2 = !( s_subtype_all & SALT_DONOR_Neg ) || nNumOtherCandidates == nNumCandidates;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nNumCandidates <= 1 || !cond1 || cond2
    // INCHI✔️❌:                  /*(
    // INCHI✔️❌:                  ( cond2a || cond2b    || cond2c )
    // INCHI✔️❌:                  ?        !(s_subtype_all & (SALT_DONOR_Neg | SALT_DONOR_H))
    // INCHI✔️❌:                  :        ( !(s_subtype_all & SALT_DONOR_Neg) || nNumOtherCandidates==nNumCandidates) ) */
    // INCHI✔️❌:                  )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 s_group_info->num_candidates = -1; /* no candidate exists */
    // INCHI✔️❌:                 goto quick_exit;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!( s_subtype_all & ( SALT_DONOR_Neg ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 t_group_info->bTautFlagsDone |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < nNumCandidates; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 i1 = s_candidate[i].atnumber;
    // INCHI✔️❌:                 if (0 <= ( s_type = GetSaltChargeType( at, i1, t_group_info, &s_subtype ) )
    // INCHI✔️❌: #if ( INCL_NON_SALT_CANDIDATATES == 1 )
    // INCHI✔️❌:                      || 0 < ( s_type = GetOtherSaltChargeType( at, i1, t_group_info, &s_subtype, 1 /* bAccept_O*/ ) )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                      )
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     s_candidate[nNumCandidates].type = s_type;
    // INCHI✔️❌:                     s_candidate[nNumCandidates].subtype = s_subtype;
    // INCHI✔️❌:                     s_candidate[nNumCandidates].endpoint = at[i1].endpoint;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Look for alt paths connecting:
    // INCHI✔️❌:         SALT_DONOR_Neg to SALT_ACCEPTOR  : long distance migration of negative charges
    // INCHI✔️❌:         SALT_DONOR_H   to SALT_ACCEPTOR  : long distance migration of H-atoms
    // INCHI✔️❌:         */
    // INCHI✔️❌:         do
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nNumChanges = 0;
    // INCHI✔️❌:             for (i1 = 0; i1 < nNumCandidates; i1++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 j1 = s_candidate[i1].atnumber;
    // INCHI✔️❌:                 for (i2 = i1 + 1; i2 < nNumCandidates; i2++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* prev. approach: do not test if both candidates are not "salt-type". Disabled 2004-03-18
    // INCHI✔️❌:                     if ( s_candidate[i1].type && s_candidate[i2].type )
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                     */
    // INCHI✔️❌:                     j2 = s_candidate[i2].atnumber;
    // INCHI✔️❌:                     if (at[j1].endpoint && at[j1].endpoint == at[j2].endpoint)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     for (j = 0; j < 2; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (j)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ii1 = i2; /* candidate 1 (donor)    ordering number */
    // INCHI✔️❌:                             ii2 = i1; /* candidate 2 (acceptor) ordering number */
    // INCHI✔️❌:                             jj1 = j2; /* candidate 1 (donor)    atom number */
    // INCHI✔️❌:                             jj2 = j1; /* candidate 2 (acceptor) atom number */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* transposition */
    // INCHI✔️❌:                             ii1 = i1; /* candidate 1 (donor)    ordering number */
    // INCHI✔️❌:                             ii2 = i2; /* candidate 2 (acceptor) ordering number */
    // INCHI✔️❌:                             jj1 = j1; /* candidate 1 (donor)    atom number     */
    // INCHI✔️❌:                             jj2 = j2; /* candidate 2 (acceptor) atom number     */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (( s_candidate[ii1].subtype & ( SALT_DONOR_Neg | SALT_DONOR_H ) ) &&
    // INCHI✔️❌:                             ( s_candidate[ii2].subtype & SALT_ACCEPTOR ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /****printf("\nChkpt @ %-s:%-d ", __FILE__,__LINE__); fflush(stdout);*/
    // INCHI✔️❌:
    // INCHI✔️❌:                             ret = bExistsAltPath( pCG, pBNS, pBD, NULL, at, num_atoms, jj2, jj1, ALT_PATH_MODE_4_SALT );
    // INCHI✔️❌:
    // INCHI✔️❌:                             num_tested++;
    // INCHI✔️❌:                             if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 err = ret;
    // INCHI✔️❌:                                 goto quick_exit;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (ret & 1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 nDelta = ( ret & ~3 ) >> 2;
    // INCHI✔️❌:                                 nNumChanges += ( ret & 2 );
    // INCHI✔️❌:                                 for (i = 0; i < 2; i++)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     jj = i ? jj2 : jj1;
    // INCHI✔️❌:                                     AddEndPoint( EndPoint + i, at, jj );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                                 /* add/merge taut groups and reinit pBNS in the fly */
    // INCHI✔️❌:                                 ret = RegisterEndPoints( pCG, t_group_info,
    // INCHI✔️❌:                                                          EndPoint, 2, at,
    // INCHI✔️❌:                                                          num_atoms, c_group_info,
    // INCHI✔️❌:                                                          pBNS );
    // INCHI✔️❌:                                 if (ret == -1)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     ret = BNS_PROGRAM_ERR;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (ret < 0)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     err = ret;
    // INCHI✔️❌:                                     goto quick_exit;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nDelta)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     err = BNS_RADICAL_ERR;
    // INCHI✔️❌:                                     goto quick_exit;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 nNumChanges += ( ret > 0 );
    // INCHI✔️❌:                                 break; /* avoid redundant repetition */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nTotNumChanges += nNumChanges;
    // INCHI✔️❌:         } while (num_tested && nNumChanges);
    // INCHI✔️❌:
    // INCHI✔️❌:     quick_exit:
    // INCHI✔️❌:         if (!err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nTotNumChanges += nNumChanges; /* nNumChanges != 0 only in case of 'goto quick_exit' */
    // INCHI✔️❌:             if (s_group_info->num_candidates == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* first time: initialize */
    // INCHI✔️❌:                 s_group_info->num_candidates = num_tested ? nNumCandidates : -1; /* no candidate exists */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nTotNumChanges = err;
    // INCHI✔️❌:         }
    // INCHI✔️❌: #if ( MAX_LOCAL_TGNUM > 0 )
    // INCHI✔️❌:         if (ptgData != tgData)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( ptgData );
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nTotNumChanges;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkSaltChargeGroups
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MarkSaltChargeGroups
    // INCHI✔️❌: #define SALT_WITH_PROTONS 1
    // INCHI✔️❌: #define INCL_NON_SALT_CANDIDATATES 1
    // INCHI✔️❌: #define MAX_LOCAL_TGNUM 0
    // INCHI✔️❌: #define ALT_PATH_MODE_4_SALT 3
    // INCHI✔️❌: #define TG_FLAG_ALLOW_NO_NEGTV_O 32
    // INCHI✔️❌: #define TG_FLAG_ALLOW_NO_NEGTV_O_DONE 64
    // INCHI✔️❌: #define TG_FLAG_FOUND_SALT_CHARGES_DONE 8192
    // INCHI✔️❌: #define FLAG_FORCE_SALT_TAUT 50
    // INCHI✔️❌: #define BNS_PROGRAM_ERR -9997
    // INCHI✔️❌: #define BNS_VERT_EDGE_OVFL -9993
    // INCHI✔️❌: #define BNS_RADICAL_ERR -9988
    // INCHI✔️❌: #define IS_BNS_ERROR(X) (BNS_ERR <= (X) && (X) <= BNS_MAX_ERR_VALUE)
    // END INCHI ACTIVE MACRO CONFIGURATION: MarkSaltChargeGroups

    let Some(s_group_info) = s_group_info else {
        return Ok(0);
    };
    if s_group_info.s_candidate.is_null() || s_group_info.max_num_candidates <= 0 {
        return Ok(0);
    }

    let s_candidate = s_group_info.s_candidate;
    let n_max_num_candidates = s_group_info.max_num_candidates;
    let mut n_num_candidates = s_group_info.num_candidates;
    let Some(t_group_info) = t_group_info else {
        return Ok(0);
    };
    if n_num_candidates <= -1 || t_group_info.t_group.is_null() {
        return Ok(0);
    }

    let max_count =
        usize::try_from(n_max_num_candidates).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if heap.slice(s_candidate.as_const())?.len() < max_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let group_count = usize::try_from(t_group_info.num_t_groups)
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let groups = heap.slice(t_group_info.t_group.as_const())?;
    if groups.len() < group_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut n_group_number = 0_i32;
    for group in groups.iter().take(group_count) {
        n_group_number = n_group_number.max(i32::from(group.nGroupNumber));
    }
    let _ = n_group_number;

    n_num_candidates = 0;
    let mut n_num_other_candidates = 0_i32;
    let mut n_num_p_only_candidates = 0_i32;
    let mut s_subtype_all = 0_i32;
    for i in 0..num_atoms {
        let mut s_subtype = 0_i32;
        let mut s_type = {
            let atoms = heap.slice(at.as_const())?;
            GetSaltChargeType(heap, atoms, i, Some(&*t_group_info), &mut s_subtype)?
        };
        let mut accepted = s_type == 0;
        if !accepted {
            s_type = {
                let atoms = heap.slice(at.as_const())?;
                GetOtherSaltChargeType(heap, atoms, i, Some(&*t_group_info), &mut s_subtype, 1)?
            };
            accepted = s_type == 1;
        }
        if !accepted {
            s_type = {
                let atoms = heap.slice(at.as_const())?;
                GetOtherSaltType(atoms, i, &mut s_subtype)?
            };
            accepted = s_type == 2;
        }
        if !accepted {
            continue;
        }
        if n_num_candidates >= n_max_num_candidates {
            return Ok(BNS_VERT_EDGE_OVFL);
        }
        let endpoint = heap
            .slice(at.as_const())?
            .get(i as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint;
        let candidate = heap
            .slice_mut(s_candidate)?
            .get_mut(n_num_candidates as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        candidate.atnumber = i as AT_NUMB;
        candidate.type_ = s_type as S_CHAR;
        candidate.subtype = s_subtype as S_CHAR;
        candidate.endpoint = endpoint;
        n_num_candidates += 1;
        n_num_other_candidates += i32::from(s_type == 1);
        n_num_p_only_candidates += i32::from(s_type == 2);
        s_subtype_all |= s_subtype;
    }
    let _ = n_num_p_only_candidates;

    let cond1 = (s_subtype_all & SALT_ACCEPTOR as i32) != 0;
    let cond2a = (t_group_info.bTautFlags & u64::from(TG_FLAG_ALLOW_NO_NEGTV_O)) != 0;
    let cond2b = (t_group_info.bTautFlagsDone & u64::from(TG_FLAG_FOUND_SALT_CHARGES_DONE)) != 0;
    let cond2c = (t_group_info.tni.bNormalizationFlags & u64::from(FLAG_FORCE_SALT_TAUT)) != 0;
    let cond2 = if cond2a || cond2b || cond2c {
        (s_subtype_all & (SALT_DONOR_Neg | SALT_DONOR_H) as i32) == 0
    } else {
        (s_subtype_all & SALT_DONOR_Neg as i32) == 0 || n_num_other_candidates == n_num_candidates
    };
    if n_num_candidates <= 1 || !cond1 || cond2 {
        s_group_info.num_candidates = -1;
        return Ok(0);
    }
    if (s_subtype_all & SALT_DONOR_Neg as i32) == 0 {
        t_group_info.bTautFlagsDone |= u64::from(TG_FLAG_ALLOW_NO_NEGTV_O_DONE);
    }

    let mut n_tot_num_changes = 0_i32;
    let mut num_tested = 0_i32;
    loop {
        let mut n_num_changes = 0_i32;
        for i1 in 0..n_num_candidates {
            let j1 = i32::from(heap.slice(s_candidate.as_const())?[i1 as usize].atnumber);
            for i2 in i1 + 1..n_num_candidates {
                let j2 = i32::from(heap.slice(s_candidate.as_const())?[i2 as usize].atnumber);
                let atoms = heap.slice(at.as_const())?;
                if atoms[j1 as usize].endpoint != 0
                    && atoms[j1 as usize].endpoint == atoms[j2 as usize].endpoint
                {
                    continue;
                }
                for orientation in 0..2 {
                    let (donor_index, acceptor_index, donor_atom, acceptor_atom) =
                        if orientation != 0 {
                            (i2, i1, j2, j1)
                        } else {
                            (i1, i2, j1, j2)
                        };
                    let candidates = heap.slice(s_candidate.as_const())?;
                    if (i32::from(candidates[donor_index as usize].subtype)
                        & (SALT_DONOR_Neg | SALT_DONOR_H) as i32)
                        == 0
                        || (i32::from(candidates[acceptor_index as usize].subtype)
                            & SALT_ACCEPTOR as i32)
                            == 0
                    {
                        continue;
                    }

                    let ret = bExistsAltPath(
                        heap,
                        pCG,
                        pBNS,
                        pBD,
                        None,
                        at,
                        num_atoms,
                        acceptor_atom,
                        donor_atom,
                        ALT_PATH_MODE_4_SALT as i32,
                        clock_result,
                    )?;
                    num_tested += 1;
                    if ichitaut_is_bns_error(ret) {
                        return Ok(ret);
                    }
                    if (ret & 1) == 0 {
                        continue;
                    }

                    let n_delta = (ret & !3) >> 2;
                    n_num_changes += ret & 2;
                    let atoms = heap.slice(at.as_const())?.to_vec();
                    let mut endpoints = [T_ENDPOINT::default(), T_ENDPOINT::default()];
                    AddEndPoint(&mut endpoints[0], &atoms, donor_atom)?;
                    AddEndPoint(&mut endpoints[1], &atoms, acceptor_atom)?;
                    let mut register_ret = RegisterEndPoints(
                        heap,
                        pCG,
                        t_group_info,
                        &mut endpoints,
                        2,
                        at,
                        num_atoms,
                        c_group_info.as_deref_mut(),
                        Some(pBNS),
                    )?;
                    if register_ret == -1 {
                        register_ret = BNS_PROGRAM_ERR;
                    }
                    if register_ret < 0 {
                        return Ok(register_ret);
                    }
                    if n_delta != 0 {
                        return Ok(BNS_RADICAL_ERR);
                    }
                    n_num_changes += i32::from(register_ret > 0);
                    break;
                }
            }
        }
        n_tot_num_changes += n_num_changes;
        if num_tested == 0 || n_num_changes == 0 {
            if s_group_info.num_candidates == 0 {
                s_group_info.num_candidates = if num_tested != 0 {
                    n_num_candidates
                } else {
                    -1
                };
            }
            return Ok(n_tot_num_changes);
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn is_centerpoint_elem(el_number: u8) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:157 is_centerpoint_elem
    // INCHI✔️✔️: int is_centerpoint_elem( U_CHAR el_number )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     switch (el_number) {
    // INCHI✔️✔️:         case EL_NUMBER_C:
    // INCHI✔️✔️:         case EL_NUMBER_N:
    // INCHI✔️✔️:         case EL_NUMBER_P:
    // INCHI✔️✔️:         case EL_NUMBER_S:
    // INCHI✔️✔️:         case EL_NUMBER_I:
    // INCHI✔️✔️:         case EL_NUMBER_AS:
    // INCHI✔️✔️:         case EL_NUMBER_SB:
    // INCHI✔️✔️:         case EL_NUMBER_SE:
    // INCHI✔️✔️:         case EL_NUMBER_TE:
    // INCHI✔️✔️:         case EL_NUMBER_CL:
    // INCHI✔️✔️:         case EL_NUMBER_BR:
    // INCHI✔️✔️:             return 1;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: is_centerpoint_elem
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: is_centerpoint_elem
    // INCHI✔️✔️: #define EL_NUMBER_C 6
    // INCHI✔️✔️: #define EL_NUMBER_N 7
    // INCHI✔️✔️: #define EL_NUMBER_P 15
    // INCHI✔️✔️: #define EL_NUMBER_S 16
    // INCHI✔️✔️: #define EL_NUMBER_I 53
    // INCHI✔️✔️: #define EL_NUMBER_AS 33
    // INCHI✔️✔️: #define EL_NUMBER_SB 51
    // INCHI✔️✔️: #define EL_NUMBER_SE 34
    // INCHI✔️✔️: #define EL_NUMBER_TE 52
    // INCHI✔️✔️: #define EL_NUMBER_CL 17
    // INCHI✔️✔️: #define EL_NUMBER_BR 35
    // END INCHI ACTIVE MACRO CONFIGURATION: is_centerpoint_elem

    match el_number {
        6 | 7 | 15 | 16 | 53 | 33 | 51 | 34 | 52 | 17 | 35 => 1,
        _ => 0,
    }
}

#[allow(non_snake_case)]
pub(crate) fn is_centerpoint_elem_KET(el_number: u8) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:182 is_centerpoint_elem_KET
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )  /* post v.1 feature */
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️: /****************************************************************************/
    // INCHI✔️✔️: int is_centerpoint_elem_KET( U_CHAR el_number )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     return el_number == EL_NUMBER_C;
    // INCHI✔️✔️: }
    // INCHI✔️✔️: #endif
    // END INCHI C FUNCTION: is_centerpoint_elem_KET
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: is_centerpoint_elem_KET
    // INCHI✔️✔️: #define KETO_ENOL_TAUT 1
    // INCHI✔️✔️: #define EL_NUMBER_C 6
    // END INCHI ACTIVE MACRO CONFIGURATION: is_centerpoint_elem_KET

    i32::from(el_number == 6)
}

#[allow(non_snake_case)]
pub(crate) fn is_centerpoint_elem_strict(el_number: u8) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:190 is_centerpoint_elem_strict
    // INCHI✔️✔️: int is_centerpoint_elem_strict( U_CHAR el_number )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     switch (el_number) {
    // INCHI✔️✔️:         case EL_NUMBER_C:
    // INCHI✔️✔️:         case EL_NUMBER_N:
    // INCHI✔️✔️:         case EL_NUMBER_P:
    // INCHI✔️✔️:         case EL_NUMBER_AS:
    // INCHI✔️✔️:         case EL_NUMBER_SB:
    // INCHI✔️✔️:             return 1;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: is_centerpoint_elem_strict
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: is_centerpoint_elem_strict
    // INCHI✔️✔️: #define EL_NUMBER_C 6
    // INCHI✔️✔️: #define EL_NUMBER_N 7
    // INCHI✔️✔️: #define EL_NUMBER_P 15
    // INCHI✔️✔️: #define EL_NUMBER_AS 33
    // INCHI✔️✔️: #define EL_NUMBER_SB 51
    // END INCHI ACTIVE MACRO CONFIGURATION: is_centerpoint_elem_strict

    match el_number {
        6 | 7 | 15 | 33 | 51 => 1,
        _ => 0,
    }
}

#[allow(non_snake_case)]
pub(crate) fn bCanBeACPoint(
    at: &inp_ATOM,
    cCharge: S_CHAR,
    cChangeValence: S_CHAR,
    neutral_bonds_valence: S_CHAR,
    neutral_valence: S_CHAR,
    nEndpointValence: S_CHAR,
    cChargeSubtype: &mut S_CHAR,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2047 bCanBeACPoint
    // INCHI✔️✔️: int bCanBeACPoint( inp_ATOM *at,
    // INCHI✔️✔️:                    S_CHAR cCharge,
    // INCHI✔️✔️:                    S_CHAR cChangeValence,
    // INCHI✔️✔️:                    S_CHAR neutral_bonds_valence,
    // INCHI✔️✔️:                    S_CHAR neutral_valence,
    // INCHI✔️✔️:                    S_CHAR nEndpointValence,
    // INCHI✔️✔️:                    S_CHAR *cChargeSubtype )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int nChangeValence;
    // INCHI✔️✔️:     int nNumBonds;
    // INCHI✔️✔️:     int nBondsValence;
    // INCHI✔️✔️:     int bNegCharge = ( at->charge == -1 );  /* add fict. bonds to (-) 2004-02-24*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at->charge == cCharge && at->valence == at->chem_bonds_valence && at->num_H)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* proton donors candidates >NH(+)-, >NH2(+), -NH3(+), >OH(+), -OH2(+) */
    // INCHI✔️✔️:         /* charged, added p-transfer -- 01-28-2004 */
    // INCHI✔️✔️:         nChangeValence = at->charge * cChangeValence; /* +1 or -1; currently only +1 */
    // INCHI✔️✔️:         nBondsValence = at->chem_bonds_valence + at->num_H;
    // INCHI✔️✔️:         if (nBondsValence == neutral_bonds_valence + nChangeValence && nEndpointValence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *cChargeSubtype = C_SUBTYPE_CHARGED_p_DONOR; /* ignore Phosphorus p-donors for now */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (at->charge == cCharge && at->valence < at->chem_bonds_valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* the requirement at->valence < at->chem_bonds_valence rejects
    // INCHI✔️✔️:             candidates >NH(+)-, >NH2(+), -NH3(+), >N(+)<, >OH(+), -OH2(+), >O(+)-
    // INCHI✔️✔️:             Moveable charge requires double bonds; these ions have no double bonds
    // INCHI✔️✔️:             */
    // INCHI✔️✔️:
    // INCHI✔️✔️:             /* charged */
    // INCHI✔️✔️:             nChangeValence = at->charge * cChangeValence; /* +1 or -1; currently only +1 */
    // INCHI✔️✔️:             nBondsValence = at->chem_bonds_valence + at->num_H;
    // INCHI✔️✔️:             nNumBonds = at->valence + at->num_H;
    // INCHI✔️✔️:             if (nBondsValence == neutral_bonds_valence + nChangeValence)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* known valence */
    // INCHI✔️✔️:                 if (nNumBonds == neutral_valence)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* non-tautomeric: >N(+)=, =O(+)-
    // INCHI✔️✔️:                     possibly tautomeric donor: =NH(+)-, =NH2(+), =OH(+) */
    // INCHI✔️✔️:                     if (at->valence == neutral_valence || !nEndpointValence)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* non-tautomeric: >N(+)=, =O(+)-; any suitable P+: >P(+)=, =PH(+)-, =PH2(+) */
    // INCHI✔️✔️:                         *cChargeSubtype = C_SUBTYPE_CHARGED_NON_TAUT;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* possibly tautomeric donor: =NH(+)-, =NH2(+), =OH(+) */
    // INCHI✔️✔️:                         *cChargeSubtype = C_SUBTYPE_CHARGED_H_DONOR;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     return 1;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (nNumBonds == neutral_valence - 1)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* possibly tutomeric acceptor: =N(+)=, #N(+)-, #NH(+), #O(+) */
    // INCHI✔️✔️:                     if (nEndpointValence)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         *cChargeSubtype = at->num_H ? C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR : C_SUBTYPE_CHARGED_H_ACCEPT;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         *cChargeSubtype = C_SUBTYPE_CHARGED_NON_TAUT;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     return 1; /* charge type, charged */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at->charge == 0 || bNegCharge)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* neutral atom or anion, all bonds are single */
    // INCHI✔️✔️:                 nBondsValence = at->chem_bonds_valence + at->num_H + bNegCharge; /* add fict. bonds to (-) 2004-02-24*/
    // INCHI✔️✔️:                 nNumBonds = at->valence + at->num_H + bNegCharge; /* add fict. bonds to (-) 2004-02-24*/
    // INCHI✔️✔️:                 if (nBondsValence == neutral_bonds_valence)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (nNumBonds == neutral_valence)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* only single bonds: >N-, >NH, -NH2, -O-, -OH, >P- >PH -PH2 */
    // INCHI✔️✔️:                         /*                    >N(-), -NH(-), -O(-). >P(-) -PH(-) */
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (at->valence == neutral_valence || !nEndpointValence)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 /* >N-, -O-, any P(3 single bonds): >P- >PH -PH2  */
    // INCHI✔️✔️:                                 *cChargeSubtype = C_SUBTYPE_NEUTRAL_NON_TAUT;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             else
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 if (at->valence < neutral_valence /*&& nEndpointValence */)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     /* num_H > 0: >NH -NH2 -OH */
    // INCHI✔️✔️:                                     /* num_H = 0: none C_SUBTYPE_NEUTRAL_H_ACCEPT for now */
    // INCHI✔️✔️:                                     *cChargeSubtype = at->num_H ? C_SUBTYPE_NEUTRAL_H_DONOR : C_SUBTYPE_NEUTRAL_H_ACCEPT;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                                 else
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     return 0;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         return 1; /* charge type, neutral */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     if (nNumBonds == neutral_valence - 1)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* possibly tautomeric acceptor =N-, =NH, =O or non-taut =P-, =PH */
    // INCHI✔️✔️:                         if (nEndpointValence)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             /* =N-,  =NH, =O  */
    // INCHI✔️✔️:                             *cChargeSubtype = C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         else
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             /* =P-, =PH */
    // INCHI✔️✔️:                             *cChargeSubtype = C_SUBTYPE_NEUTRAL_NON_TAUT;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         return 1; /* charge type, (+) => neutral */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: bCanBeACPoint
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: bCanBeACPoint
    // INCHI✔️✔️: #define C_SUBTYPE_CHARGED_NON_TAUT          (C_SUBTYPE_CHARGED)
    // INCHI✔️✔️: #define C_SUBTYPE_CHARGED_p_DONOR           (C_SUBTYPE_CHARGED|C_SUBTYPE_p_DONOR)
    // INCHI✔️✔️: #define C_SUBTYPE_CHARGED_H_ACCEPT          (C_SUBTYPE_CHARGED|C_SUBTYPE_H_ACCEPT)
    // INCHI✔️✔️: #define C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR  (C_SUBTYPE_CHARGED|C_SUBTYPE_H_ACCEPT|C_SUBTYPE_p_DONOR)
    // INCHI✔️✔️: #define C_SUBTYPE_CHARGED_H_DONOR           (C_SUBTYPE_CHARGED|C_SUBTYPE_H_DONOR |C_SUBTYPE_p_DONOR)
    // INCHI✔️✔️: #define C_SUBTYPE_NEUTRAL_NON_TAUT          (C_SUBTYPE_NEUTRAL)
    // INCHI✔️✔️: #define C_SUBTYPE_NEUTRAL_H_ACCEPT          (C_SUBTYPE_NEUTRAL|C_SUBTYPE_H_ACCEPT)
    // INCHI✔️✔️: #define C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT (C_SUBTYPE_NEUTRAL|C_SUBTYPE_H_ACCEPT|C_SUBTYPE_p_ACCEPT)
    // INCHI✔️✔️: #define C_SUBTYPE_NEUTRAL_H_DONOR           (C_SUBTYPE_NEUTRAL|C_SUBTYPE_H_DONOR)
    // END INCHI ACTIVE MACRO CONFIGURATION: bCanBeACPoint

    let b_neg_charge = i32::from(at.charge == -1);
    if at.charge == cCharge && at.valence == at.chem_bonds_valence && at.num_H != 0 {
        let n_change_valence = i32::from(at.charge) * i32::from(cChangeValence);
        let n_bonds_valence = i32::from(at.chem_bonds_valence) + i32::from(at.num_H);
        if n_bonds_valence == i32::from(neutral_bonds_valence) + n_change_valence
            && nEndpointValence != 0
        {
            *cChargeSubtype = C_SUBTYPE_CHARGED_p_DONOR as S_CHAR;
        }
        return 0;
    }

    if at.charge == cCharge && at.valence < at.chem_bonds_valence {
        let n_change_valence = i32::from(at.charge) * i32::from(cChangeValence);
        let n_bonds_valence = i32::from(at.chem_bonds_valence) + i32::from(at.num_H);
        let n_num_bonds = i32::from(at.valence) + i32::from(at.num_H);
        if n_bonds_valence == i32::from(neutral_bonds_valence) + n_change_valence {
            if n_num_bonds == i32::from(neutral_valence) {
                if at.valence == neutral_valence || nEndpointValence == 0 {
                    *cChargeSubtype = C_SUBTYPE_CHARGED_NON_TAUT as S_CHAR;
                } else {
                    *cChargeSubtype = C_SUBTYPE_CHARGED_H_DONOR as S_CHAR;
                }
                return 1;
            }
            if n_num_bonds == i32::from(neutral_valence) - 1 {
                if nEndpointValence != 0 {
                    *cChargeSubtype = if at.num_H != 0 {
                        C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR as S_CHAR
                    } else {
                        C_SUBTYPE_CHARGED_H_ACCEPT as S_CHAR
                    };
                } else {
                    *cChargeSubtype = C_SUBTYPE_CHARGED_NON_TAUT as S_CHAR;
                }
                return 1;
            }
        }
    } else if at.charge == 0 || b_neg_charge != 0 {
        let n_bonds_valence = i32::from(at.chem_bonds_valence) + i32::from(at.num_H) + b_neg_charge;
        let n_num_bonds = i32::from(at.valence) + i32::from(at.num_H) + b_neg_charge;
        if n_bonds_valence == i32::from(neutral_bonds_valence) {
            if n_num_bonds == i32::from(neutral_valence) {
                if at.valence == neutral_valence || nEndpointValence == 0 {
                    *cChargeSubtype = C_SUBTYPE_NEUTRAL_NON_TAUT as S_CHAR;
                } else if at.valence < neutral_valence {
                    *cChargeSubtype = if at.num_H != 0 {
                        C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR
                    } else {
                        C_SUBTYPE_NEUTRAL_H_ACCEPT as S_CHAR
                    };
                } else {
                    return 0;
                }
                return 1;
            }
            if n_num_bonds == i32::from(neutral_valence) - 1 {
                if nEndpointValence != 0 {
                    *cChargeSubtype = C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT as S_CHAR;
                } else {
                    *cChargeSubtype = C_SUBTYPE_NEUTRAL_NON_TAUT as S_CHAR;
                }
                return 1;
            }
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn GetChargeType(atom: &[inp_ATOM], iat: i32, cChargeSubtype: &mut S_CHAR) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2181 GetChargeType
    // INCHI✔️✔️: int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, n;
    // INCHI✔️✔️:     S_CHAR    nEndpointValence;
    // INCHI✔️✔️:     inp_ATOM *at = atom + iat;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     *cChargeSubtype = 0;
    // INCHI✔️✔️:     /* ignore ion pairs and charges != 1 */
    // INCHI✔️✔️:     if (abs( at->charge ) == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0; i < at->valence; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             n = at->neighbor[i];
    // INCHI✔️✔️:             /* allow negatively charged tautomeric neighbors 2004-02-26 */
    // INCHI✔️✔️:             if (abs( atom[n].charge + at->charge ) < abs( atom[n].charge - at->charge ) && !atom[n].endpoint)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return -1; /* charges have different signs */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (at->charge)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1; /* abs(charge) != 1 */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* find candidates */
    // INCHI✔️✔️:     for (i = 0; i < NUM_C_TYPES; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!strcmp( at->elname, CType[i].elname ) &&
    // INCHI✔️✔️:             ( !CType[i].num_bonds || (CType[i].num_bonds == at->valence && at->nNumAtInRingSystem >= 5) )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nEndpointValence = (S_CHAR) get_endpoint_valence( at->el_number );
    // INCHI✔️✔️:             if (bCanBeACPoint( at, CType[i].charge, CType[i].cChangeValence, CType[i].neutral_bonds_valence,
    // INCHI✔️✔️:                                CType[i].neutral_valence, nEndpointValence, cChargeSubtype ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return CType[i].cChargeType;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return -1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetChargeType
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetChargeType
    // INCHI✔️✔️: const CHARGE_TYPE CType[] =
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     { "N\0",  1, 3, 3, 1, 0, 0 },
    // INCHI✔️✔️:     { "P\0",  1, 3, 3, 1, 1, 0 },
    // INCHI✔️✔️: #if ( ADD_MOVEABLE_O_PLUS == 1 )
    // INCHI✔️✔️:     { "O\0",  1, 2, 2, 1, 2, 2 }, /* added 02-06-2005 */
    // INCHI✔️✔️:     { "S\0",  1, 2, 2, 1, 3, 2 }, /* added 03-18-2005 */
    // INCHI✔️✔️:     { "Se",   1, 2, 2, 1, 4, 2 }, /* added 03-18-2005 */
    // INCHI✔️✔️:     { "Te",   1, 2, 2, 1, 5, 2 }, /* added 03-18-2005 */
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️: };
    // INCHI✔️✔️: #define NUM_C_TYPES  (int)(sizeof( CType )/sizeof(CType[0]))
    // INCHI✔️✔️: ADD_MOVEABLE_O_PLUS == 1 is active in the selected libinchi build.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetChargeType

    let at = &atom[iat as usize];
    *cChargeSubtype = 0;
    if i32::from(at.charge).abs() == 1 {
        for i in 0..usize::try_from(at.valence).unwrap_or(0) {
            let n = at.neighbor[i] as usize;
            if (i32::from(atom[n].charge) + i32::from(at.charge)).abs()
                < (i32::from(atom[n].charge) - i32::from(at.charge)).abs()
                && atom[n].endpoint == 0
            {
                return -1;
            }
        }
    } else if at.charge != 0 {
        return -1;
    }

    for ctype in &CTYPE {
        if source_strcmp_zero(&at.elname, &ctype.elname) == 0
            && (ctype.num_bonds == 0
                || (ctype.num_bonds == at.valence && at.nNumAtInRingSystem >= 5))
        {
            let n_endpoint_valence = get_endpoint_valence(at.el_number) as S_CHAR;
            if bCanBeACPoint(
                at,
                ctype.charge,
                ctype.cChangeValence,
                ctype.neutral_bonds_valence,
                ctype.neutral_valence,
                n_endpoint_valence,
                cChargeSubtype,
            ) != 0
            {
                return i32::from(ctype.cChargeType);
            }
        }
    }

    -1
}

fn charged_cpoint(at: &[inp_ATOM], point: i32) -> u16 {
    u16::from(at[point as usize].charge == 1)
}

#[allow(non_snake_case)]
pub(crate) fn GetSaltChargeType(
    heap: &SourceHeap,
    at: &[inp_ATOM],
    at_no: i32,
    t_group_info: Option<&T_GROUP_INFO>,
    s_subtype: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2565 GetSaltChargeType
    // INCHI✔️✔️: int GetSaltChargeType( inp_ATOM *at,
    // INCHI✔️✔️:                        int at_no,
    // INCHI✔️✔️:                        T_GROUP_INFO *t_group_info,
    // INCHI✔️✔️:                        int *s_subtype )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     type (returned value):
    // INCHI✔️✔️:     -1 => ignore
    // INCHI✔️✔️:     0 => oxygen
    // INCHI✔️✔️:     subtype:
    // INCHI✔️✔️:     1 = SALT_DONOR_H   => has H
    // INCHI✔️✔️:     2 = SALT_DONOR_Neg => has (-) charge
    // INCHI✔️✔️:     4 = SALT_ACCEPTOR  => may be an acceptor of H or (-), but not necessarily
    // INCHI✔️✔️:
    // INCHI✔️✔️:     O-atom should be:
    // INCHI✔️✔️:     - a terminal atom
    // INCHI✔️✔️:     - connected to unsaturated, uncharged, non-radical atom C that has chemical valence 4:
    // INCHI✔️✔️:     H-donors:             =CH-OH, =C(-X)-OH
    // INCHI✔️✔️:     possible H-acceptors: -CH=O, >C=O
    // INCHI✔️✔️:     H-acceptors are true if O is tautomeric
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int iC, tg, i, type;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     *s_subtype = 0; /* initialize the output */
    // INCHI✔️✔️:                     /* check whether it is a candidate */
    // INCHI✔️✔️:     if (at[at_no].valence != 1 ||
    // INCHI✔️✔️:          (at[at_no].radical && at[at_no].radical != RADICAL_SINGLET) ||
    // INCHI✔️✔️:          at[at_no].charge < -1 ||
    // INCHI✔️✔️:          (at[at_no].charge > 0 && !at[at_no].c_point)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at[at_no].el_number == EL_NUMBER_O ||
    // INCHI✔️✔️:          at[at_no].el_number == EL_NUMBER_S ||
    // INCHI✔️✔️:          at[at_no].el_number == EL_NUMBER_SE ||
    // INCHI✔️✔️:          at[at_no].el_number == EL_NUMBER_TE)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = 0;  /* terminal oxygen atom, needs more to be checked... */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         type = -1; /* ignore this atom */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (type < 0 ||
    // INCHI✔️✔️:          at[at_no].chem_bonds_valence + at[at_no].num_H !=
    // INCHI✔️✔️:          get_el_valence( at[at_no].el_number, at[at_no].charge, 0 ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* non-standard valence or not an oxygen */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     iC = at[at_no].neighbor[0];
    // INCHI✔️✔️:
    // INCHI✔️✔️: #if ( SALT_WITH_PROTONS == 1 )
    // INCHI✔️✔️:     if (at[iC].el_number != EL_NUMBER_C ||
    // INCHI✔️✔️:          at[iC].chem_bonds_valence + at[iC].num_H != 4 || /* allow =C(H)-OH or -C(H)=O */
    // INCHI✔️✔️:          at[iC].charge ||
    // INCHI✔️✔️:          (at[iC].radical && at[iC].radical != RADICAL_SINGLET) ||
    // INCHI✔️✔️:          at[iC].valence == at[iC].chem_bonds_valence) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* oxigen is connected to a wrong atom */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     if (( tg = at[at_no].endpoint ) && t_group_info && t_group_info->t_group)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* O-atom is in a tautomeric group */
    // INCHI✔️✔️:         for (i = 0; i < t_group_info->num_t_groups; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (tg == t_group_info->t_group[i].nGroupNumber)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /*
    // INCHI✔️✔️:                 t_group_info->t_group[i].num[0] = number of attached H-atoms and negative charges
    // INCHI✔️✔️:                 t_group_info->t_group[i].num[1] = number of attached negative charges
    // INCHI✔️✔️:                 */
    // INCHI✔️✔️:                 if (t_group_info->t_group[i].num[0] > t_group_info->t_group[i].num[1])
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     *s_subtype |= SALT_DONOR_H; /* has H */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (t_group_info->t_group[i].num[1])
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     *s_subtype |= SALT_DONOR_Neg; /* has (-) */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 *s_subtype |= SALT_ACCEPTOR; /* there is always an acceptor in a t-group */
    // INCHI✔️✔️:                 return type;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return -1; /* error: t-group not found */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* O is not not in a tautomeric group */
    // INCHI✔️✔️:     /* assume valence(O-) < valence(O) < valence(O+) */
    // INCHI✔️✔️:     if (at[at_no].charge == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *s_subtype |= SALT_DONOR_Neg; /* has (-) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (at[at_no].charge <= 0 && at[at_no].num_H)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *s_subtype |= SALT_DONOR_H; /* has H */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (at[at_no].charge == 0 && at[at_no].chem_bonds_valence == 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *s_subtype |= SALT_ACCEPTOR;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     /* since O cannot be a charge point, the following cannot happen: */
    // INCHI✔️✔️:     if (at[at_no].charge == 1 && at[at_no].c_point && at[at_no].chem_bonds_valence == 2 && at[at_no].num_H)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *s_subtype |= SALT_DONOR_H; /* has H */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return type;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetSaltChargeType
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetSaltChargeType
    // INCHI✔️✔️: #define SALT_WITH_PROTONS 1
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define EL_NUMBER_S  ((U_CHAR)16)
    // INCHI✔️✔️: #define EL_NUMBER_SE ((U_CHAR)34)
    // INCHI✔️✔️: #define EL_NUMBER_TE ((U_CHAR)52)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define SALT_DONOR_H      1
    // INCHI✔️✔️: #define SALT_DONOR_Neg    2
    // INCHI✔️✔️: #define SALT_ACCEPTOR     4
    // INCHI✔️✔️: The #else branch for SALT_WITH_PROTONS != 1 is inactive in the selected libinchi target.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetSaltChargeType

    let at_index = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = at
        .get(at_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *s_subtype = 0;
    if atom.valence != 1
        || (atom.radical != 0 && atom.radical != RADICAL_SINGLET as S_CHAR)
        || atom.charge < -1
        || (atom.charge > 0 && atom.c_point == 0)
    {
        return Ok(-1);
    }

    let type_ = match atom.el_number {
        8 | 16 | 34 | 52 => 0,
        _ => -1,
    };

    if type_ < 0
        || i32::from(atom.chem_bonds_valence) + i32::from(atom.num_H)
            != get_el_valence(i32::from(atom.el_number), i32::from(atom.charge), 0)?
    {
        return Ok(-1);
    }

    let carbon_index = usize::from(atom.neighbor[0]);
    let carbon = at
        .get(carbon_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if carbon.el_number != 6
        || i32::from(carbon.chem_bonds_valence) + i32::from(carbon.num_H) != 4
        || carbon.charge != 0
        || (carbon.radical != 0 && carbon.radical != RADICAL_SINGLET as S_CHAR)
        || carbon.valence == carbon.chem_bonds_valence
    {
        return Ok(-1);
    }

    let tg = atom.endpoint;
    if tg != 0
        && let Some(t_group_info) = t_group_info
        && !t_group_info.t_group.is_null()
    {
        let t_groups = heap.slice(t_group_info.t_group.as_const())?;
        for i in 0..usize::try_from(t_group_info.num_t_groups)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?
        {
            let t_group = t_groups.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if tg == t_group.nGroupNumber {
                if t_group.num[0] > t_group.num[1] {
                    *s_subtype |= SALT_DONOR_H as i32;
                }
                if t_group.num[1] != 0 {
                    *s_subtype |= SALT_DONOR_Neg as i32;
                }
                *s_subtype |= SALT_ACCEPTOR as i32;
                return Ok(type_);
            }
        }
        return Ok(-1);
    }

    if atom.charge == -1 {
        *s_subtype |= SALT_DONOR_Neg as i32;
    }
    if atom.charge <= 0 && atom.num_H != 0 {
        *s_subtype |= SALT_DONOR_H as i32;
    }
    if atom.charge == 0 && atom.chem_bonds_valence == 2 {
        *s_subtype |= SALT_ACCEPTOR as i32;
    }
    if atom.charge == 1 && atom.c_point != 0 && atom.chem_bonds_valence == 2 && atom.num_H != 0 {
        *s_subtype |= SALT_DONOR_H as i32;
    }

    Ok(type_)
}

#[allow(non_snake_case)]
pub(crate) fn bDoNotMergeNonTautAtom(at: &[inp_ATOM], at_no: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2691 bDoNotMergeNonTautAtom
    // INCHI✔️✔️: int bDoNotMergeNonTautAtom( inp_ATOM *at, int at_no )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (at[at_no].el_number == EL_NUMBER_N)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: bDoNotMergeNonTautAtom
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: bDoNotMergeNonTautAtom
    // INCHI✔️✔️: #define EL_NUMBER_N ((U_CHAR)7)
    // END INCHI ACTIVE MACRO CONFIGURATION: bDoNotMergeNonTautAtom

    let atom = at
        .get(usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(atom.el_number == 7))
}

#[allow(non_snake_case)]
pub(crate) fn GetOtherSaltChargeType(
    heap: &SourceHeap,
    at: &[inp_ATOM],
    at_no: i32,
    t_group_info: Option<&T_GROUP_INFO>,
    s_subtype: &mut i32,
    bAccept_O: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2702 GetOtherSaltChargeType
    // INCHI✔️✔️: int GetOtherSaltChargeType( inp_ATOM *at,
    // INCHI✔️✔️:                             int at_no,
    // INCHI✔️✔️:                             T_GROUP_INFO *t_group_info,
    // INCHI✔️✔️:                             int *s_subtype,
    // INCHI✔️✔️:                             int bAccept_O )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     type (returned value):
    // INCHI✔️✔️:     -1 => ignore
    // INCHI✔️✔️:     1 => not an oxygen
    // INCHI✔️✔️:     subtype:
    // INCHI✔️✔️:     1 = SALT_DONOR_H   => has H
    // INCHI✔️✔️:     2 = SALT_DONOR_Neg => has (-) charge
    // INCHI✔️✔️:     4 = SALT_ACCEPTOR  => may be an acceptor of H or (-), but not necessarily
    // INCHI✔️✔️:
    // INCHI✔️✔️:     the atom should be:
    // INCHI✔️✔️:     - a tautomeric endpoint atom
    // INCHI✔️✔️:     - connected to possible centerpoint atom
    // INCHI✔️✔️:
    // INCHI✔️✔️:     another description of the atom searched here:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     any possibly tautomeric atom adjacent to a possibly centerpoint
    // INCHI✔️✔️:     that has at least one double bond (possibly if positively charged);
    // INCHI✔️✔️:     if eif.cAcceptor then the bond between the atom and the centerpoint must be possibly double
    // INCHI✔️✔️:     if eif.cAcceptor then the bond must be possibly single
    // INCHI✔️✔️:     Donors that belong to a t-group are also acceptors
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int tg, i, j, type, endpoint_valence, num_centerpoints, bond_type, centerpoint;
    // INCHI✔️✔️:     ENDPOINT_INFO eif;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     *s_subtype = 0; /* initialize the output */
    // INCHI✔️✔️:     if (!bAccept_O /* only N */ &&
    // INCHI✔️✔️:         ( at[at_no].el_number == EL_NUMBER_O ||
    // INCHI✔️✔️:           at[at_no].el_number == EL_NUMBER_S ||
    // INCHI✔️✔️:           at[at_no].el_number == EL_NUMBER_SE ||
    // INCHI✔️✔️:           at[at_no].el_number == EL_NUMBER_TE ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* we are not looking for oxygen here */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     type = 1;
    // INCHI✔️✔️:     if (!( endpoint_valence = nGetEndpointInfo( at, at_no, &eif ) )) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* not a possible endpoint */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* at[at_no] is not not in a tautomeric group; use eif previously filled out by nGetEndpointInfo */
    // INCHI✔️✔️:         /* check whether there is adjacent atom-candidate for a centerpoint */
    // INCHI✔️✔️:         num_centerpoints = 0;
    // INCHI✔️✔️:         for (j = 0; j < at[at_no].valence; j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bond_type = (int) at[at_no].bond_type[j] & BOND_TYPE_MASK;
    // INCHI✔️✔️:             centerpoint = (int) at[at_no].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️✔️:             if (( (eif.cAcceptor && ( bond_type == BOND_DOUBLE ||
    // INCHI✔️✔️:                                      bond_type == BOND_ALTERN || /* possibly double */
    // INCHI✔️✔️:                                      bond_type == BOND_ALT12NS ||
    // INCHI✔️✔️:                                      bond_type == BOND_TAUTOM )) ||
    // INCHI✔️✔️:                   (eif.cDonor && ( bond_type == BOND_SINGLE ||
    // INCHI✔️✔️:                                   bond_type == BOND_ALTERN || /* possibly single */
    // INCHI✔️✔️:                                   bond_type == BOND_ALT12NS ||
    // INCHI✔️✔️:                                   bond_type == BOND_TAUTOM )) ) &&
    // INCHI✔️✔️:                                   ( at[centerpoint].chem_bonds_valence > at[centerpoint].valence ||
    // INCHI✔️✔️:                                     /* check for possible endpoint added 2004-02-24 */
    // INCHI✔️✔️:                                     (at[centerpoint].chem_bonds_valence == at[centerpoint].valence &&
    // INCHI✔️✔️:                                     ( at[centerpoint].endpoint || at[centerpoint].c_point )) /* tautomerism or charge may increment at[centerpoint].chem_bonds_valence*/ ) &&
    // INCHI✔️✔️:                  is_centerpoint_elem( at[centerpoint].el_number )) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_centerpoints++;
    // INCHI✔️✔️:                 break; /* at least one possibly centerpoint neighbor has been found */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (!num_centerpoints)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* moved here from just after "type = 1;" line 2004-02-26 */
    // INCHI✔️✔️:         if (( tg = at[at_no].endpoint ) && t_group_info && t_group_info->t_group)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* atom is in a tautomeric group */
    // INCHI✔️✔️:             for (i = 0; i < t_group_info->num_t_groups; i++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (tg == t_group_info->t_group[i].nGroupNumber)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /*
    // INCHI✔️✔️:                     t_group_info->t_group[i].num[0] = number of attached H-atoms and negative charges
    // INCHI✔️✔️:                     t_group_info->t_group[i].num[1] = number of attached negative charges
    // INCHI✔️✔️:                     */
    // INCHI✔️✔️:                     if (t_group_info->t_group[i].num[0] > t_group_info->t_group[i].num[1])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         *s_subtype |= SALT_DONOR_H; /* has H */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     if (t_group_info->t_group[i].num[1])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         *s_subtype |= SALT_DONOR_Neg; /* has (-) */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     *s_subtype |= SALT_ACCEPTOR; /* there is always an acceptor in a t-group */
    // INCHI✔️✔️:                     return type;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             return -1; /* error: t-group not found */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (eif.cAcceptor)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *s_subtype |= SALT_ACCEPTOR;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (eif.cDonor)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at[at_no].charge == -1)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *s_subtype |= SALT_DONOR_Neg; /* has (-) */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (at[at_no].num_H)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *s_subtype |= SALT_DONOR_H; /* has H */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return type;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetOtherSaltChargeType
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetOtherSaltChargeType
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define EL_NUMBER_S  ((U_CHAR)16)
    // INCHI✔️✔️: #define EL_NUMBER_SE ((U_CHAR)34)
    // INCHI✔️✔️: #define EL_NUMBER_TE ((U_CHAR)52)
    // INCHI✔️✔️: #define BOND_SINGLE 1
    // INCHI✔️✔️: #define BOND_DOUBLE 2
    // INCHI✔️✔️: #define BOND_ALTERN 4
    // INCHI✔️✔️: #define BOND_TAUTOM 8
    // INCHI✔️✔️: #define BOND_ALT12NS 9
    // INCHI✔️✔️: #define BOND_TYPE_MASK 0x0f
    // END INCHI ACTIVE MACRO CONFIGURATION: GetOtherSaltChargeType

    let at_index = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = at
        .get(at_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *s_subtype = 0;
    if bAccept_O == 0 && matches!(atom.el_number, 8 | 16 | 34 | 52) {
        return Ok(-1);
    }

    let type_ = 1;
    let mut eif = ENDPOINT_INFO::default();
    let endpoint_valence = nGetEndpointInfo(at, at_no, &mut eif);
    if endpoint_valence == 0 {
        return Ok(-1);
    }

    let mut num_centerpoints = 0_i32;
    let valence = usize::try_from(atom.valence).unwrap_or(0);
    for j in 0..valence {
        let bond_type = u32::from(atom.bond_type[j]) & BOND_TYPE_MASK;
        let centerpoint_index = usize::from(atom.neighbor[j]);
        let centerpoint = at
            .get(centerpoint_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bond_matches_endpoint_role = (eif.cAcceptor != 0
            && matches!(
                bond_type,
                BOND_DOUBLE | BOND_ALTERN | BOND_ALT12NS | BOND_TAUTOM
            ))
            || (eif.cDonor != 0
                && matches!(
                    bond_type,
                    BOND_SINGLE | BOND_ALTERN | BOND_ALT12NS | BOND_TAUTOM
                ));
        let centerpoint_has_available_valence = centerpoint.chem_bonds_valence
            > centerpoint.valence
            || (centerpoint.chem_bonds_valence == centerpoint.valence
                && (centerpoint.endpoint != 0 || centerpoint.c_point != 0));
        if bond_matches_endpoint_role
            && centerpoint_has_available_valence
            && is_centerpoint_elem(centerpoint.el_number) != 0
        {
            num_centerpoints += 1;
            break;
        }
    }
    if num_centerpoints == 0 {
        return Ok(-1);
    }

    let tg = atom.endpoint;
    if tg != 0
        && let Some(t_group_info) = t_group_info
        && !t_group_info.t_group.is_null()
    {
        let t_groups = heap.slice(t_group_info.t_group.as_const())?;
        let count = if t_group_info.num_t_groups > 0 {
            t_group_info.num_t_groups as usize
        } else {
            0
        };
        for i in 0..count {
            let t_group = t_groups.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if tg == t_group.nGroupNumber {
                if t_group.num[0] > t_group.num[1] {
                    *s_subtype |= SALT_DONOR_H as i32;
                }
                if t_group.num[1] != 0 {
                    *s_subtype |= SALT_DONOR_Neg as i32;
                }
                *s_subtype |= SALT_ACCEPTOR as i32;
                return Ok(type_);
            }
        }
        return Ok(-1);
    }

    if eif.cAcceptor != 0 {
        *s_subtype |= SALT_ACCEPTOR as i32;
    }
    if eif.cDonor != 0 {
        if atom.charge == -1 {
            *s_subtype |= SALT_DONOR_Neg as i32;
        }
        if atom.num_H != 0 {
            *s_subtype |= SALT_DONOR_H as i32;
        }
    }

    Ok(type_)
}

#[allow(non_snake_case)]
pub(crate) fn GetOtherSaltType(
    at: &[inp_ATOM],
    at_no: i32,
    s_subtype: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2828 GetOtherSaltType
    // INCHI✔️✔️: int GetOtherSaltType( inp_ATOM *at, int at_no, int *s_subtype )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     type (returned value):
    // INCHI✔️✔️:     -1 => ignore
    // INCHI✔️✔️:     2 => found:                           SH
    // INCHI✔️✔️:     proton donor     -CH2-SH, >CH-SH, >C<    S(-)
    // INCHI✔️✔️:     proton acceptor  -CH2-S(-), >CH-S(-), >C<
    // INCHI✔️✔️:     subtype:
    // INCHI✔️✔️:     1 = SALT_DONOR_H   => has H
    // INCHI✔️✔️:     2 = SALT_DONOR_Neg => has (-) charge
    // INCHI✔️✔️:     4 = SALT_ACCEPTOR  => may be an acceptor of H or (-), but not necessarily
    // INCHI✔️✔️:
    // INCHI✔️✔️:     non-O-atom should be:
    // INCHI✔️✔️:     - a tautomeric endpoint atom
    // INCHI✔️✔️:     - connected to possible middle point atom
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int type, endpoint_valence, centerpoint; /* djb-rwth: removing redundant variables */
    // INCHI✔️✔️:     ENDPOINT_INFO eif;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at[at_no].valence != 1 || at[at_no].chem_bonds_valence != 1 ||
    // INCHI✔️✔️:          1 != ( at[at_no].num_H == 1 ) + ( at[at_no].charge == -1 ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     *s_subtype = 0; /* initialize the output */
    // INCHI✔️✔️:     if (!( at[at_no].el_number == EL_NUMBER_S ||
    // INCHI✔️✔️:            at[at_no].el_number == EL_NUMBER_SE ||
    // INCHI✔️✔️:            at[at_no].el_number == EL_NUMBER_TE ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* we are not looking for oxygen here */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     type = 2; /* non-tautomeric p-donor or acceptor: C-SH, C-S(-) */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!( endpoint_valence = nGetEndpointInfo( at, at_no, &eif ) ) ||
    // INCHI✔️✔️:          (eif.cMoveableCharge && !at[at_no].c_point) || !eif.cDonor || eif.cAcceptor) /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* not a possible -SH or -S(-) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* at[at_no] is not not in a tautomeric group; use eif previously filled out by nGetEndpointInfo */
    // INCHI✔️✔️:         /* check whether there is adjacent atom-candidate for a centerpoint */
    // INCHI✔️✔️:         centerpoint = (int) at[at_no].neighbor[0];
    // INCHI✔️✔️:         /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:         if (at[centerpoint].el_number != EL_NUMBER_C ||
    // INCHI✔️✔️:              at[centerpoint].charge ||
    // INCHI✔️✔️:              (at[centerpoint].radical && at[centerpoint].radical != RADICAL_SINGLET) ||
    // INCHI✔️✔️:              at[centerpoint].valence != at[centerpoint].chem_bonds_valence) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1; /* not a carbon with all single bonds */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (at[at_no].num_H == 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *s_subtype |= SALT_p_DONOR;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at[at_no].charge == -1)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *s_subtype |= SALT_p_ACCEPTOR;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return -1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return type;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetOtherSaltType
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetOtherSaltType
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_S  ((U_CHAR)16)
    // INCHI✔️✔️: #define EL_NUMBER_SE ((U_CHAR)34)
    // INCHI✔️✔️: #define EL_NUMBER_TE ((U_CHAR)52)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define SALT_p_DONOR 8
    // INCHI✔️✔️: #define SALT_p_ACCEPTOR 16
    // END INCHI ACTIVE MACRO CONFIGURATION: GetOtherSaltType

    let at_index = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = at
        .get(at_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.valence != 1
        || atom.chem_bonds_valence != 1
        || 1 != i32::from(atom.num_H == 1) + i32::from(atom.charge == -1)
    {
        return Ok(-1);
    }

    *s_subtype = 0;
    if !matches!(atom.el_number, 16 | 34 | 52) {
        return Ok(-1);
    }

    let type_ = 2;
    let mut eif = ENDPOINT_INFO::default();
    if nGetEndpointInfo(at, at_no, &mut eif) == 0
        || (eif.cMoveableCharge != 0 && atom.c_point == 0)
        || eif.cDonor == 0
        || eif.cAcceptor != 0
    {
        return Ok(-1);
    }

    let centerpoint_index = usize::from(atom.neighbor[0]);
    let centerpoint = at
        .get(centerpoint_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if centerpoint.el_number != 6
        || centerpoint.charge != 0
        || (centerpoint.radical != 0 && centerpoint.radical != RADICAL_SINGLET as S_CHAR)
        || centerpoint.valence != centerpoint.chem_bonds_valence
    {
        return Ok(-1);
    }
    if atom.num_H == 1 {
        *s_subtype |= SALT_p_DONOR as i32;
    } else if atom.charge == -1 {
        *s_subtype |= SALT_p_ACCEPTOR as i32;
    } else {
        return Ok(-1);
    }

    Ok(type_)
}

const REGISTER_END_POINTS_STACK_LEN: usize = MAX_STACK_ARRAY_LEN as usize + 1;

struct RegisterEndPointsScratch<T> {
    stack: [MaybeUninit<T>; REGISTER_END_POINTS_STACK_LEN],
    heap: Option<SourceMutPointer<T>>,
}

impl<T: Copy + Default + 'static> RegisterEndPointsScratch<T> {
    fn new() -> Self {
        Self {
            // The three source stack arrays are deliberately uninitialized.
            // RegisterEndPoints writes every element before reading it.
            stack: [MaybeUninit::uninit(); REGISTER_END_POINTS_STACK_LEN],
            heap: None,
        }
    }

    fn ensure_heap(
        &mut self,
        heap: &mut SourceHeap,
        len: usize,
        initialized: usize,
    ) -> Result<(), SourceHeapError> {
        if self.heap.is_none() {
            let pointer = inchi_calloc::<T>(heap, len as u64, std::mem::size_of::<T>() as u64)?;
            {
                let target = heap.slice_mut(pointer)?;
                let target = target
                    .get_mut(..initialized)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                // SAFETY: callers pass the initialized prefix length maintained
                // by the source loops. MaybeUninit<T> has the same layout as T,
                // and the newly allocated target cannot alias the stack array.
                unsafe {
                    std::ptr::copy_nonoverlapping(
                        self.stack.as_ptr().cast::<T>(),
                        target.as_mut_ptr(),
                        initialized,
                    );
                }
            }
            self.heap = Some(pointer);
        }
        Ok(())
    }

    fn write(
        &mut self,
        heap: &mut SourceHeap,
        index: usize,
        value: T,
    ) -> Result<(), SourceHeapError> {
        if let Some(pointer) = self.heap {
            heap.slice_mut(pointer)?
                .get_mut(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from(&value);
        } else {
            self.stack
                .get_mut(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .write(value);
        }
        Ok(())
    }

    fn read(&self, heap: &SourceHeap, index: usize) -> Result<T, SourceHeapError> {
        if let Some(pointer) = self.heap {
            heap.slice(pointer.as_const())?
                .get(index)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)
        } else {
            self.stack
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)
                // SAFETY: every source read is bounded by the number of values
                // already written, or follows the explicit zero initialization.
                .map(|value| unsafe { value.assume_init_read() })
        }
    }

    fn clear_prefix(&mut self, heap: &mut SourceHeap, len: usize) -> Result<(), SourceHeapError> {
        if let Some(pointer) = self.heap {
            heap.slice_mut(pointer)?
                .get_mut(..len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .fill(T::default());
        } else {
            let stack = self
                .stack
                .get_mut(..len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for value in stack {
                value.write(T::default());
            }
        }
        Ok(())
    }

    fn free(&mut self, heap: &mut SourceHeap) -> Result<(), SourceHeapError> {
        if let Some(pointer) = self.heap.take() {
            inchi_free(heap, pointer)?;
        }
        Ok(())
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RegisterEndPoints(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    t_group_info: &mut T_GROUP_INFO,
    EndPoint: &mut [T_ENDPOINT],
    nNumEndPoints: i32,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    cgi: Option<&mut C_GROUP_INFO>,
    mut pBNS: Option<&mut BN_STRUCT>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1021 RegisterEndPoints
    // INCHI✔️❌: int RegisterEndPoints( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                        T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                        T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                        int nNumEndPoints,
    // INCHI✔️❌:                        inp_ATOM *at,
    // INCHI✔️❌:                        int num_atoms,
    // INCHI✔️❌:                        C_GROUP_INFO *cgi,
    // INCHI✔️❌:                        struct BalancedNetworkStructure *pBNS )
    // INCHI✔️❌: {
    // INCHI✔️❌:     T_GROUP  *t_group = t_group_info->t_group;
    // INCHI✔️❌:     int      *pnum_t = &t_group_info->num_t_groups;
    // INCHI✔️❌:     int       max_num_t = t_group_info->max_num_t_groups;
    // INCHI✔️❌:     int       nNumZeroEqu, nNumNewTGroups;
    // INCHI✔️❌:     AT_NUMB   group, prev_group, prev_eqnum, nNextGroupNumber, nLeastGroupNumber;
    // INCHI✔️❌:     int       nNumGroups, num_t, difference;
    // INCHI✔️❌:     int       i, j, k, ret;
    // INCHI✔️❌:     AT_NUMB   nNewTgNumberStackArray[MAX_STACK_ARRAY_LEN + 1];
    // INCHI✔️❌:     AT_NUMB   nGroupNumberStackArray[MAX_STACK_ARRAY_LEN + 1];
    // INCHI✔️❌:     AT_NUMB   nGroupNewNumberStackArray[MAX_STACK_ARRAY_LEN + 1];
    // INCHI✔️❌:     AT_NUMB  *nNewTgNumber = nNewTgNumberStackArray;
    // INCHI✔️❌:     AT_NUMB  *nGroupNumber = nGroupNumberStackArray;
    // INCHI✔️❌:     AT_NUMB  *nGroupNewNumber = nGroupNewNumberStackArray;
    // INCHI✔️❌: ...
    // INCHI✔️❌:     return difference;
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌: ...
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: RegisterEndPoints
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RegisterEndPoints
    // INCHI✔️❌: #define MAX_STACK_ARRAY_LEN 127
    // INCHI✔️❌: #define ALWAYS_ADD_TG_ON_THE_FLY 1
    // INCHI✔️❌: #define KETO_ENOL_TAUT 1
    // INCHI✔️❌: #define TAUT_PT_22_00 1
    // INCHI✔️❌: #define TAUT_PT_16_00 1
    // INCHI✔️❌: #define TAUT_PT_06_00 1
    // INCHI✔️❌: #define TAUT_PT_39_00 1
    // INCHI✔️❌: #define TAUT_PT_13_00 1
    // INCHI✔️❌: #define TAUT_PT_18_00 1
    // INCHI✔️❌: #define TG_FLAG_MOVE_POS_CHARGES 8
    // END INCHI ACTIVE MACRO CONFIGURATION: RegisterEndPoints

    fn is_bns_error(value: i32) -> bool {
        BNS_ERR <= value && value <= BNS_MAX_ERR_VALUE
    }

    macro_rules! c_alloc {
        ($expr:expr) => {
            match $expr {
                Ok(value) => value,
                Err(
                    SourceHeapError::AllocationFailed
                    | SourceHeapError::AllocationSizeOverflow
                    | SourceHeapError::AllocationElementCountOutOfRange,
                ) => return Ok(-1),
                Err(error) => return Err(error),
            }
        };
    }

    let mut n_new_tg_number = RegisterEndPointsScratch::<AT_NUMB>::new();
    let mut n_group_number = RegisterEndPointsScratch::<AT_NUMB>::new();
    let mut n_group_new_number = RegisterEndPointsScratch::<AT_NUMB>::new();
    let result = (|| -> Result<i32, SourceHeapError> {
        let t_group = t_group_info.t_group;
        let pnum_t = &mut t_group_info.num_t_groups;
        let max_num_t = t_group_info.max_num_t_groups;
        let mut n_num_zero_equ = 0_i32;
        let mut n_num_new_t_groups = 0_i32;
        let mut n_next_group_number: AT_NUMB = 0;
        let mut n_least_group_number: AT_NUMB;
        let mut n_num_groups = 0_i32;
        let mut num_t = *pnum_t;
        let mut difference = 0_i32;
        let mut i = 0_i32;
        let mut j = 0_i32;
        let mut k = 0_i32;
        let mut ret = 0_i32;

        if nNumEndPoints <= 0 {
            return Ok(0);
        }

        let num_endpoints =
            usize::try_from(nNumEndPoints).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atom_count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let num_t_usize =
            usize::try_from(num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let max_num_t_usize =
            usize::try_from(max_num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let t_group_slice = heap.slice(t_group.as_const())?;
        if t_group_slice.len() < max_num_t_usize {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let mut t_group_snapshot = t_group_slice[..num_t_usize].to_vec();

        for group in t_group_slice.iter().take(num_t_usize) {
            if n_next_group_number < group.nGroupNumber {
                n_next_group_number = group.nGroupNumber;
            }
        }
        n_next_group_number = n_next_group_number.wrapping_add(1);
        n_least_group_number = n_next_group_number;
        let first = &EndPoint[0];
        let prev_group = first.nGroupNumber;
        let prev_eqnum = first.nEquNumber;
        for idx in 0..num_endpoints {
            let endpoint = &EndPoint[idx];
            if endpoint.nGroupNumber != 0 && endpoint.nGroupNumber < n_least_group_number {
                n_least_group_number = endpoint.nGroupNumber;
            }
            j += i32::from(prev_group == endpoint.nGroupNumber);
            k += i32::from(prev_eqnum == endpoint.nEquNumber);
            n_num_zero_equ += i32::from(endpoint.nEquNumber == 0);
        }
        if j == nNumEndPoints && prev_group != 0 && k == nNumEndPoints {
            return Ok(0);
        }

        if n_num_zero_equ == 0 {
            for idx in 0..num_endpoints {
                let group = EndPoint[idx].nEquNumber;
                if group >= n_next_group_number {
                    let mut found_pos = None;
                    for pos in 0..n_num_new_t_groups as usize {
                        if group == n_group_new_number.read(heap, pos)? {
                            found_pos = Some(pos);
                            break;
                        }
                    }
                    let mapped_offset = if let Some(pos) = found_pos {
                        pos
                    } else {
                        if n_num_new_t_groups as usize == MAX_STACK_ARRAY_LEN as usize {
                            c_alloc!(n_group_new_number.ensure_heap(
                                heap,
                                num_endpoints,
                                n_num_new_t_groups as usize
                            ));
                        }
                        n_group_new_number.write(heap, n_num_new_t_groups as usize, group)?;
                        let pos = n_num_new_t_groups as usize;
                        n_num_new_t_groups += 1;
                        pos
                    };
                    let mapped = n_next_group_number.wrapping_add(mapped_offset as AT_NUMB);
                    EndPoint[idx].nEquNumber = mapped;
                }
            }
        } else if n_num_zero_equ == nNumEndPoints {
            if n_least_group_number == n_next_group_number {
                n_num_new_t_groups = 1;
            }
            for endpoint in EndPoint.iter_mut().take(num_endpoints) {
                endpoint.nEquNumber = n_least_group_number;
            }
        } else {
            ret = -1;
            return Ok(ret);
        }

        if n_num_new_t_groups != 0 {
            if num_t.wrapping_add(n_num_new_t_groups) > max_num_t {
                ret = -1;
                return Ok(ret);
            }
            let t_groups = heap.slice_mut(t_group)?;
            let new_count = usize::try_from(n_num_new_t_groups)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let start = num_t_usize;
            for offset in 0..new_count {
                t_groups[start + offset] = T_GROUP::default();
                t_groups[start + offset].nGroupNumber =
                    n_next_group_number.wrapping_add(offset as AT_NUMB);
            }
            t_group_snapshot = t_groups[..usize::try_from(num_t)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                + new_count]
                .to_vec();
        }

        n_num_groups = 0;
        i = 0;
        j = 0;
        while i < nNumEndPoints {
            let endpoint = &EndPoint[i as usize];
            let group = endpoint.nGroupNumber;
            if group != 0 {
                if group == endpoint.nEquNumber {
                    i += 1;
                    continue;
                }
                let mut found = false;
                for pos in 0..n_num_groups as usize {
                    if group == n_group_number.read(heap, pos)? {
                        if endpoint.nEquNumber != n_group_new_number.read(heap, pos)? {
                            ret = -1;
                            return Ok(ret);
                        }
                        found = true;
                        break;
                    }
                }
                if !found {
                    if n_num_groups as usize == MAX_STACK_ARRAY_LEN as usize {
                        c_alloc!(n_group_new_number.ensure_heap(
                            heap,
                            num_endpoints,
                            n_num_groups as usize
                        ));
                        c_alloc!(n_group_number.ensure_heap(
                            heap,
                            num_endpoints,
                            n_num_groups as usize
                        ));
                    }
                    n_group_number.write(heap, n_num_groups as usize, group)?;
                    n_group_new_number.write(heap, n_num_groups as usize, endpoint.nEquNumber)?;
                    n_num_groups += 1;
                }
            } else {
                let group = endpoint.nEquNumber;
                if group >= n_next_group_number {
                    j = num_t + i32::from(group - n_next_group_number);
                } else if j >= num_t
                    || group
                        != t_group_snapshot[usize::try_from(j)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?]
                        .nGroupNumber
                {
                    j = 0;
                    while j < num_t {
                        if group
                            == t_group_snapshot[usize::try_from(j)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?]
                            .nGroupNumber
                        {
                            break;
                        }
                        j += 1;
                    }
                    if j == num_t {
                        ret = -1;
                        return Ok(ret);
                    }
                }
                {
                    let t_groups = heap.slice_mut(t_group)?;
                    let group_index =
                        usize::try_from(j).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    t_groups[group_index].nNumEndpoints =
                        t_groups[group_index].nNumEndpoints.wrapping_add(1);
                    for idx in 0..t_groups[group_index].num.len() {
                        t_groups[group_index].num[idx] =
                            t_groups[group_index].num[idx].wrapping_add(endpoint.num[idx]);
                    }
                    for idx in 0..t_groups[group_index].num_DA.len() {
                        t_groups[group_index].num_DA[idx] =
                            t_groups[group_index].num_DA[idx].wrapping_add(endpoint.num_DA[idx]);
                    }
                }
                let atom_index = usize::try_from(endpoint.nAtomNumber)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                heap.slice_mut(at)?
                    .get_mut(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .endpoint = group;
                difference += 1;
            }
            i += 1;
        }

        difference += n_num_groups;
        num_t = num_t.wrapping_add(n_num_new_t_groups);
        if difference == 0 {
            return Ok(0);
        }

        if n_num_groups != 0 {
            t_group_snapshot = heap.slice(t_group.as_const())?
                [..usize::try_from(num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?]
                .to_vec();
            for group in t_group_snapshot.iter() {
                if n_next_group_number < group.nGroupNumber {
                    n_next_group_number = group.nGroupNumber;
                }
            }
        }

        for idx in
            0..usize::try_from(n_num_groups).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        {
            let group1 = n_group_number.read(heap, idx)?;
            let group2 = n_group_new_number.read(heap, idx)?;
            let mut i1 = -1_i32;
            let mut i2 = -1_i32;
            for jdx in
                0..usize::try_from(num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            {
                let group = t_group_snapshot[jdx].nGroupNumber;
                if i1 < 0 && group1 == group {
                    i1 = jdx as i32;
                }
                if i2 < 0 && group2 == group {
                    i2 = jdx as i32;
                }
            }
            if i1 < 0 || i2 < 0 {
                ret = -1;
                return Ok(ret);
            }
            let i1_index =
                usize::try_from(i1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let i2_index =
                usize::try_from(i2).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            {
                let t_groups = heap.slice_mut(t_group)?;
                for idx2 in 0..t_groups[i1_index].num.len() {
                    t_groups[i2_index].num[idx2] =
                        t_groups[i2_index].num[idx2].wrapping_add(t_groups[i1_index].num[idx2]);
                }
                for idx2 in 0..t_groups[i1_index].num_DA.len() {
                    t_groups[i2_index].num_DA[idx2] = t_groups[i2_index].num_DA[idx2]
                        .wrapping_add(t_groups[i1_index].num_DA[idx2]);
                }
                t_groups[i2_index].nNumEndpoints = t_groups[i2_index]
                    .nNumEndpoints
                    .wrapping_add(t_groups[i1_index].nNumEndpoints);
                num_t -= 1;
                if num_t > i1 {
                    let end = usize::try_from(num_t)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    for shift in i1_index..end {
                        t_groups[shift] = t_groups[shift + 1].clone();
                    }
                }
            }
            t_group_snapshot = heap.slice(t_group.as_const())?
                [..usize::try_from(num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?]
                .to_vec();
        }

        if n_num_groups != 0 {
            if n_next_group_number >= MAX_STACK_ARRAY_LEN as AT_NUMB {
                c_alloc!(
                    n_new_tg_number.ensure_heap(
                        heap,
                        usize::try_from(n_next_group_number)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                            + 1,
                        0
                    )
                );
            }
            let clear_len = usize::try_from(n_next_group_number)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                + 1;
            // INCHI✔️✔️:         memset( nNewTgNumber, 0, (nNextGroupNumber + 1) * sizeof(nNewTgNumber[0]) );
            n_new_tg_number.clear_prefix(heap, clear_len)?;
            for idx in
                0..usize::try_from(num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            {
                n_new_tg_number.write(
                    heap,
                    t_group_snapshot[idx].nGroupNumber as usize,
                    (idx + 1) as AT_NUMB,
                )?;
            }
            for idx in 0..usize::try_from(n_num_groups)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            {
                let old_group = n_group_number.read(heap, idx)?;
                let new_group = n_group_new_number.read(heap, idx)?;
                let old_slot = n_new_tg_number.read(heap, old_group as usize)?;
                let new_slot = n_new_tg_number.read(heap, new_group as usize)?;
                if old_slot == 0 && new_slot != 0 {
                    n_new_tg_number.write(heap, old_group as usize, new_slot)?;
                } else {
                    ret = -1;
                    return Ok(ret);
                }
            }
            let group_count =
                usize::try_from(num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let mut mapped_group_numbers = Vec::with_capacity(group_count);
            for group in t_group_snapshot.iter().take(group_count) {
                mapped_group_numbers.push(n_new_tg_number.read(heap, group.nGroupNumber as usize)?);
            }
            {
                let t_groups = heap.slice_mut(t_group)?;
                for idx in 0..group_count {
                    t_groups[idx].nGroupNumber = mapped_group_numbers[idx];
                }
            }
            for idx in 0..atom_count {
                let endpoint_value = {
                    let atom = heap
                        .slice(at.as_const())?
                        .get(idx)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    atom.endpoint
                };
                if endpoint_value != 0 {
                    let mapped = n_new_tg_number.read(heap, endpoint_value as usize)?;
                    if mapped == 0 || n_next_group_number <= mapped {
                        ret = -1;
                        return Ok(ret);
                    }
                    heap.slice_mut(at)?
                        .get_mut(idx)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint = mapped;
                }
            }
        }

        n_new_tg_number.free(heap)?;
        n_group_number.free(heap)?;
        n_group_new_number.free(heap)?;
        if t_group_info.tGroupNumber.is_null() {
            let alloc_len = usize::try_from(max_num_t)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                .checked_mul(2)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?;
            t_group_info.tGroupNumber = c_alloc!(inchi_calloc::<AT_NUMB>(
                heap,
                alloc_len as u64,
                std::mem::size_of::<AT_NUMB>() as u64,
            ));
        }
        heap.slice_mut(t_group_info.tGroupNumber)?
            .fill(AT_NUMB::default());
        {
            let group_count =
                usize::try_from(num_t).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let t_groups = heap.slice(t_group.as_const())?[..group_count].to_vec();
            let t_group_number = heap.slice_mut(t_group_info.tGroupNumber)?;
            for (idx, group) in t_groups.iter().enumerate() {
                if group.nNumEndpoints != 0 && group.nGroupNumber != 0 {
                    t_group_number
                        .get_mut(group.nGroupNumber as usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone_from(&((idx + 1) as AT_NUMB));
                }
            }
        }

        if let Some(mut pBNS) = pBNS.as_deref_mut() {
            if pBNS.tot_st_cap == pBNS.tot_st_flow || ALWAYS_ADD_TG_ON_THE_FLY == 1 {
                let mut tgi = T_GROUP_INFO {
                    t_group,
                    num_t_groups: num_t,
                    ..T_GROUP_INFO::default()
                };
                if KETO_ENOL_TAUT == 1 {
                    tgi.bTautFlags |= t_group_info.bTautFlags & TG_FLAG_KETO_ENOL_TAUT as u64;
                }
                if TAUT_PT_22_00 == 1 {
                    tgi.bTautFlags |= t_group_info.bTautFlags & TG_FLAG_PT_22_00 as u64;
                }
                if TAUT_PT_16_00 == 1 {
                    tgi.bTautFlags |= t_group_info.bTautFlags & TG_FLAG_PT_16_00 as u64;
                }
                if TAUT_PT_06_00 == 1 {
                    tgi.bTautFlags |= t_group_info.bTautFlags & TG_FLAG_PT_06_00 as u64;
                }
                if TAUT_PT_39_00 == 1 {
                    tgi.bTautFlags |= t_group_info.bTautFlags & TG_FLAG_PT_39_00 as u64;
                }
                if TAUT_PT_13_00 == 1 {
                    tgi.bTautFlags |= t_group_info.bTautFlags & TG_FLAG_PT_13_00 as u64;
                }
                if TAUT_PT_18_00 == 1 {
                    tgi.bTautFlags |= t_group_info.bTautFlags & TG_FLAG_PT_18_00 as u64;
                }
                let ret_bns = ReInitBnStruct(heap, Some(pBNS), at, num_atoms, 0)?;
                if is_bns_error(ret_bns) {
                    return Ok(ret_bns);
                }
                // INCHI✔️🔝: AddCGroups2BnStruct and AddTGroups2BnStruct receive
                // the same atom allocation that the C source passes directly.
                // ReInitBnStruct is called first with bRemoveGroupsFromAtoms=0,
                // so this stable read-only view is valid for both helpers and
                // avoids copying the complete inp_ATOM array on every call.
                let atoms_view = unsafe { heap.stable_slice(at.as_const())? };
                let atoms = atoms_view.prefix(atoms_view.len())?;
                if let Some(pb) = heap.slice(pBNS.pbTautFlags.as_const())?.first() {
                    if (*pb & TG_FLAG_MOVE_POS_CHARGES as u64) != 0 {
                        let ret_bns = AddCGroups2BnStruct(heap, pCG, pBNS, atoms, num_atoms, cgi)?;
                        if is_bns_error(ret_bns) {
                            return Ok(ret_bns);
                        }
                    }
                }
                let ret_bns =
                    AddTGroups2BnStruct(heap, pCG, pBNS, atoms, num_atoms, Some(&mut tgi))?;
                if is_bns_error(ret_bns) {
                    return Ok(ret_bns);
                }
            }
        }

        *pnum_t = num_t;
        Ok(difference)
    })();

    let cleanup = (|| -> Result<(), SourceHeapError> {
        n_new_tg_number.free(heap)?;
        n_group_number.free(heap)?;
        n_group_new_number.free(heap)?;
        Ok(())
    })();

    match (result, cleanup) {
        (Err(error), _) => Err(error),
        (Ok(value), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MergeSaltTautGroups(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    s_group_info: Option<&mut S_GROUP_INFO>,
    t_group_info: Option<&mut T_GROUP_INFO>,
    c_group_info: Option<&mut C_GROUP_INFO>,
    pBNS: Option<&mut BN_STRUCT>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:3953 MergeSaltTautGroups
    // INCHI✔️❌: int MergeSaltTautGroups( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                          inp_ATOM *at,
    // INCHI✔️❌:                          int num_atoms,
    // INCHI✔️❌:                          S_GROUP_INFO *s_group_info,
    // INCHI✔️❌:                          T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                          C_GROUP_INFO *c_group_info,
    // INCHI✔️❌:                          struct BalancedNetworkStructure *pBNS )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* Count candidates to be connected: exclude pure donors that do not belong to any t-group */
    // INCHI✔️❌:     AT_NUMB    nCurTGroupNumber;
    // INCHI✔️❌:     int        i, j, /*k,*/ ret, iat, /*nMobile,*/ nMinNumEndpoints;
    // INCHI✔️❌:     int        s_subtype_all, s_subtype_taut;
    // INCHI✔️❌:     int        nMaxNumCandidates, nNumCandidates; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     T_ENDPOINT EndPointStackArray[MAX_STACK_ARRAY_LEN]; /* will be reallocated if too short */
    // INCHI✔️❌:     T_ENDPOINT  *EndPoint = EndPointStackArray;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!s_group_info || !s_group_info->s_candidate || /*s_group_info->num_candidates <= 0 ||*/
    // INCHI✔️❌:          !t_group_info || !t_group_info->t_group || !c_group_info)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     nMinNumEndpoints = 0;
    // INCHI✔️❌:     nMaxNumCandidates = s_group_info->max_num_candidates;
    // INCHI✔️❌:     nCurTGroupNumber = MAX_ATOMS;  /* impossible t-group number */
    // INCHI✔️❌:     s_subtype_all = s_subtype_taut = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Collect tautomeric acidic O and previously non-tautomeric C-OH, C-SH, C-O(-), C-S(-)  */
    // INCHI✔️❌:     /* find whether previously found tautomeric atoms have both mobile H and (-) */
    // INCHI✔️❌:     if (1 || ( s_group_info->num_candidates < 0 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Can be only -O(-)  and -OH */
    // INCHI✔️❌:         int          s_type, s_subtype;
    // INCHI✔️❌:         S_CANDIDATE *s_candidate = s_group_info->s_candidate;
    // INCHI✔️❌:         for (i = 0, nNumCandidates = 0; i < num_atoms; i++) /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             s_subtype = 0;
    // INCHI✔️❌:             if (0 == ( s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype ) ) ||
    // INCHI✔️❌:                  /* -C=O or =C-OH, O = S, Se, Te */
    // INCHI✔️❌:
    // INCHI✔️❌:                  /*(t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT) &&*/
    // INCHI✔️❌:                  1 == ( s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1/* bAccept_O*/ ) ) ||
    // INCHI✔️❌:                  /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above. M may be N */
    // INCHI✔️❌:
    // INCHI✔️❌:                  2 == ( s_type = GetOtherSaltType( at, i, &s_subtype ) ) ||
    // INCHI✔️❌:                  /* >C-SH, >C-S(-); S=S,Se,Te */
    // INCHI✔️❌:
    // INCHI✔️❌:                  /* other proton donor or acceptor */
    // INCHI✔️❌:                 (bHasAcidicHydrogen(at, i) && ((s_type = 3), (s_subtype = SALT_p_DONOR))) ||
    // INCHI✔️❌:                 (bHasAcidicMinus(at, i) && ((s_type = 3), (s_subtype = SALT_p_ACCEPTOR)))
    // INCHI✔️❌:                 )
    // INCHI✔️❌:             {
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (nNumCandidates >= nMaxNumCandidates)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[i].endpoint)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     s_subtype_taut |= s_subtype;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bDoNotMergeNonTautAtom( at, i ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue; /* ignore non-tautomeric N */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (!( s_subtype & SALT_DONOR_ALL ) ||
    // INCHI✔️❌:                     (( s_subtype & SALT_ACCEPTOR ) && !at[i].endpoint)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;  /* do not include non-taut acceptors like -C=O */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 s_candidate[nNumCandidates].atnumber = i;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].type = s_type;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].subtype = s_subtype;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].endpoint = at[i].endpoint;
    // INCHI✔️❌:                 nNumCandidates++;
    // INCHI✔️❌:                 s_subtype_all |= s_subtype;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         Forced merging occurs upon:
    // INCHI✔️❌:         ===========================
    // INCHI✔️❌:         (t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O) or
    // INCHI✔️❌:         (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT)
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         Allow forced merging in cases:
    // INCHI✔️❌:         {t-groups}  (H, (-)}  {H, (-), t-groups}
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         Normal salt merging in cases:
    // INCHI✔️❌:         (H, (-)} {H, (-), t-groups},
    // INCHI✔️❌:
    // INCHI✔️❌:         Cannot merge H into t-groups if no (-) is present
    // INCHI✔️❌:         */
    // INCHI✔️❌:         if (( t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O ) ||
    // INCHI✔️❌:             ( t_group_info->bTautFlagsDone & TG_FLAG_FOUND_SALT_CHARGES_DONE ) ||
    // INCHI✔️❌:              ( t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Force merge even though no negative charges are present */
    // INCHI✔️❌:             if (nNumCandidates <= 1 ||
    // INCHI✔️❌:                 (( !( s_subtype_all & SALT_DONOR_Neg2 ) || !( s_subtype_all & SALT_DONOR_H2 ) ) &&
    // INCHI✔️❌:                  !t_group_info->num_t_groups)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 s_group_info->num_candidates = -1; /* no candidate exists */
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* normal salt mode: merge if both -XH and -X(-) are present */
    // INCHI✔️❌:             if (nNumCandidates <= 1 ||
    // INCHI✔️❌:                 ( !( s_subtype_all & SALT_DONOR_Neg2 ) || !( s_subtype_all & SALT_DONOR_H2 ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 s_group_info->num_candidates = -1; /* no candidate exists */
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* -- old code --
    // INCHI✔️❌:         if ( nNumCandidates <= 1 ||
    // INCHI✔️❌:         (((t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O) ||
    // INCHI✔️❌:         (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT)) ?
    // INCHI✔️❌:         !(s_subtype_all & SALT_DONOR_ALL):
    // INCHI✔️❌:         !(s_subtype_all & SALT_DONOR_Neg2)
    // INCHI✔️❌:         )
    // INCHI✔️❌:         ) {
    // INCHI✔️❌:         s_group_info->num_candidates = -1;
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         */
    // INCHI✔️❌:         if (!( s_subtype_all & ( SALT_DONOR_Neg2 ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t_group_info->bTautFlagsDone |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         s_group_info->num_candidates = nNumCandidates;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < s_group_info->num_candidates; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         iat = s_group_info->s_candidate[i].atnumber;
    // INCHI✔️❌:         if (( s_group_info->s_candidate[i].subtype & SALT_ACCEPTOR ) && !at[iat].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* should not happen */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         s_subtype_all |= s_group_info->s_candidate[i].subtype;
    // INCHI✔️❌:         if (at[iat].endpoint != nCurTGroupNumber || !at[iat].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nMinNumEndpoints++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nCurTGroupNumber = (int) at[iat].endpoint;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nMinNumEndpoints <= 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* too few endpoints */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Make sure we have enough memory */
    // INCHI✔️❌:     if (nMinNumEndpoints > MAX_STACK_ARRAY_LEN)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!( EndPoint = (T_ENDPOINT *) inchi_calloc( nMinNumEndpoints, sizeof( EndPoint[0] ) ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*printf("BNS_OUT_OF_RAM-8\n");*/
    // INCHI✔️❌:             return BNS_OUT_OF_RAM;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nCurTGroupNumber = MAX_ATOMS;  /* impossible t-group number */
    // INCHI✔️❌:     for (i = j = 0; i < s_group_info->num_candidates; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         iat = s_group_info->s_candidate[i].atnumber;
    // INCHI✔️❌:         if (s_group_info->s_candidate[i].subtype == SALT_ACCEPTOR && !at[iat].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[iat].endpoint != nCurTGroupNumber || !at[iat].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             AddEndPoint( EndPoint + j, at, iat );
    // INCHI✔️❌:             j++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nCurTGroupNumber = (int) at[iat].endpoint;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = RegisterEndPoints( pCG, t_group_info, EndPoint, j, at,
    // INCHI✔️❌:                              num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ret == -1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = BNS_PROGRAM_ERR;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (EndPoint != EndPointStackArray)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( EndPoint );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MergeSaltTautGroups
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MergeSaltTautGroups
    // INCHI✔️❌: #define MAX_STACK_ARRAY_LEN 127
    // INCHI✔️❌: #define MAX_ATOMS 32766
    // END INCHI ACTIVE MACRO CONFIGURATION: MergeSaltTautGroups

    let Some(s_group_info) = s_group_info else {
        return Ok(0);
    };
    let Some(t_group_info) = t_group_info else {
        return Ok(0);
    };
    let Some(c_group_info) = c_group_info else {
        return Ok(0);
    };
    if s_group_info.s_candidate.is_null() || t_group_info.t_group.is_null() {
        return Ok(0);
    }

    let atoms = heap.slice(at.as_const())?.to_vec();
    let n_max_num_candidates = s_group_info.max_num_candidates;
    let mut n_cur_t_group_number = MAX_ATOMS as AT_NUMB;
    let mut s_subtype_all = 0_i32;
    let mut _s_subtype_taut = 0_i32;
    let mut n_num_candidates = 0_i32;

    for i in 0..num_atoms {
        let mut s_subtype = 0_i32;
        let mut s_type = GetSaltChargeType(heap, &atoms, i, Some(&*t_group_info), &mut s_subtype)?;
        let mut found = s_type == 0;
        if !found {
            s_type =
                GetOtherSaltChargeType(heap, &atoms, i, Some(&*t_group_info), &mut s_subtype, 1)?;
            found = s_type == 1;
        }
        if !found {
            s_type = GetOtherSaltType(&atoms, i, &mut s_subtype)?;
            found = s_type == 2;
        }
        if !found && bHasAcidicHydrogen(&atoms, i)? != 0 {
            s_type = 3;
            s_subtype = SALT_p_DONOR as i32;
            found = true;
        }
        if !found && bHasAcidicMinus(&atoms, i)? != 0 {
            s_type = 3;
            s_subtype = SALT_p_ACCEPTOR as i32;
            found = true;
        }
        if !found {
            continue;
        }
        if n_num_candidates >= n_max_num_candidates {
            return Ok(BNS_VERT_EDGE_OVFL);
        }
        let atom = atoms
            .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.endpoint != 0 {
            _s_subtype_taut |= s_subtype;
        } else if bDoNotMergeNonTautAtom(&atoms, i)? != 0 {
            continue;
        }
        if (s_subtype & SALT_DONOR_ALL as i32) == 0
            || ((s_subtype & SALT_ACCEPTOR as i32) != 0 && atom.endpoint == 0)
        {
            continue;
        }
        let candidate = S_CANDIDATE {
            atnumber: i as AT_NUMB,
            type_: s_type as S_CHAR,
            subtype: s_subtype as S_CHAR,
            endpoint: atom.endpoint,
        };
        let candidate_index =
            usize::try_from(n_num_candidates).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(s_group_info.s_candidate)?
            .get_mut(candidate_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = candidate;
        n_num_candidates += 1;
        s_subtype_all |= s_subtype;
    }

    let force_merge = (t_group_info.bTautFlags & u64::from(TG_FLAG_ALLOW_NO_NEGTV_O)) != 0
        || (t_group_info.bTautFlagsDone & u64::from(TG_FLAG_FOUND_SALT_CHARGES_DONE)) != 0
        || (t_group_info.tni.bNormalizationFlags & u64::from(FLAG_FORCE_SALT_TAUT)) != 0;
    if force_merge {
        if n_num_candidates <= 1
            || (((s_subtype_all & SALT_DONOR_Neg2 as i32) == 0
                || (s_subtype_all & SALT_DONOR_H2 as i32) == 0)
                && t_group_info.num_t_groups == 0)
        {
            s_group_info.num_candidates = -1;
            return Ok(0);
        }
    } else if n_num_candidates <= 1
        || (s_subtype_all & SALT_DONOR_Neg2 as i32) == 0
        || (s_subtype_all & SALT_DONOR_H2 as i32) == 0
    {
        s_group_info.num_candidates = -1;
        return Ok(0);
    }
    if (s_subtype_all & SALT_DONOR_Neg2 as i32) == 0 {
        t_group_info.bTautFlagsDone |= u64::from(TG_FLAG_ALLOW_NO_NEGTV_O_DONE);
    }
    s_group_info.num_candidates = n_num_candidates;

    let mut n_min_num_endpoints = 0_i32;
    for i in 0..s_group_info.num_candidates {
        let candidate = heap
            .slice(s_group_info.s_candidate.as_const())?
            .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let atom = atoms
            .get(usize::from(candidate.atnumber))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if (i32::from(candidate.subtype) & SALT_ACCEPTOR as i32) != 0 && atom.endpoint == 0 {
            continue;
        }
        s_subtype_all |= i32::from(candidate.subtype);
        if atom.endpoint != n_cur_t_group_number || atom.endpoint == 0 {
            n_min_num_endpoints += 1;
        }
        n_cur_t_group_number = atom.endpoint;
    }
    if n_min_num_endpoints <= 1 {
        return Ok(0);
    }

    let dynamic_endpoint = if n_min_num_endpoints > MAX_STACK_ARRAY_LEN as i32 {
        match inchi_calloc::<T_ENDPOINT>(
            heap,
            n_min_num_endpoints as u64,
            std::mem::size_of::<T_ENDPOINT>() as u64,
        ) {
            Ok(pointer) => Some(pointer),
            Err(
                SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange,
            ) => return Ok(BNS_OUT_OF_RAM),
            Err(error) => return Err(error),
        }
    } else {
        None
    };
    let mut endpoint_stack: [T_ENDPOINT; MAX_STACK_ARRAY_LEN as usize] =
        std::array::from_fn(|_| T_ENDPOINT::default());
    let mut dynamic_endpoint_values = if let Some(pointer) = dynamic_endpoint {
        heap.slice(pointer.as_const())?.to_vec()
    } else {
        Vec::new()
    };
    let endpoints = if dynamic_endpoint.is_some() {
        dynamic_endpoint_values.as_mut_slice()
    } else {
        &mut endpoint_stack[..usize::try_from(n_min_num_endpoints)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
    };

    n_cur_t_group_number = MAX_ATOMS as AT_NUMB;
    let mut j = 0_i32;
    for i in 0..s_group_info.num_candidates {
        let candidate = heap
            .slice(s_group_info.s_candidate.as_const())?
            .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let atom = atoms
            .get(usize::from(candidate.atnumber))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if i32::from(candidate.subtype) == SALT_ACCEPTOR as i32 && atom.endpoint == 0 {
            continue;
        }
        if atom.endpoint != n_cur_t_group_number || atom.endpoint == 0 {
            AddEndPoint(
                endpoints
                    .get_mut(usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                &atoms,
                i32::from(candidate.atnumber),
            )?;
            j += 1;
        }
        n_cur_t_group_number = atom.endpoint;
    }

    let register_result = RegisterEndPoints(
        heap,
        pCG,
        t_group_info,
        endpoints,
        j,
        at,
        num_atoms,
        Some(c_group_info),
        pBNS,
    );
    let cleanup_result = if let Some(pointer) = dynamic_endpoint {
        inchi_free(heap, pointer)
    } else {
        Ok(())
    };
    let mut ret = match (register_result, cleanup_result) {
        (Err(error), _) => return Err(error),
        (Ok(_), Err(error)) => return Err(error),
        (Ok(value), Ok(())) => value,
    };
    if ret == -1 {
        ret = BNS_PROGRAM_ERR;
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn MakeIsotopicHGroup(
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    s_group_info: Option<&mut S_GROUP_INFO>,
    t_group_info: Option<&mut T_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:4156 MakeIsotopicHGroup
    // INCHI✔️❌: int MakeIsotopicHGroup( inp_ATOM *at,
    // INCHI✔️❌:                         int num_atoms,
    // INCHI✔️❌:                         S_GROUP_INFO *s_group_info,
    // INCHI✔️❌:                         T_GROUP_INFO *t_group_info )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* All tautomeric atoms and all possible H+ donors and acceptors that have H */
    // INCHI✔️❌:     int        i, j, k, n, bHasH, tg, nError = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int        nMaxNumCandidates, nNumCandidates, nNumNonTautCandidates;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!s_group_info || !s_group_info->s_candidate || /*s_group_info->num_candidates <= 0 ||*/
    // INCHI✔️❌:          !t_group_info || !t_group_info->t_group)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     nMaxNumCandidates = s_group_info->max_num_candidates;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     memset( t_group_info->num_iso_H, 0, sizeof( t_group_info->num_iso_H ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     if (1 || ( s_group_info->num_candidates < 0 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int s_type, s_subtype;
    // INCHI✔️❌:         S_CANDIDATE *s_candidate = s_group_info->s_candidate;
    // INCHI✔️❌:         for (i = 0, nNumCandidates = nNumNonTautCandidates = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             s_subtype = 0;
    // INCHI✔️❌:             s_type = 0;
    // INCHI✔️❌:             if (at[i].endpoint)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (( tg = t_group_info->tGroupNumber[at[i].endpoint] ) &&
    // INCHI✔️❌:                      at[i].endpoint == t_group_info->t_group[tg -= 1].nGroupNumber)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bHasH = (int) t_group_info->t_group[tg].num[0] - (int) t_group_info->t_group[tg].num[1];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nError = BNS_PROGRAM_ERR;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bHasH = (int) at[i].num_H;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((bHasH && at[i].endpoint) || /* tautomeric atoms */
    // INCHI✔️❌:                                            /* non-tautomeric heteroatoms that
    // INCHI✔️❌:                                            (a) have H and
    // INCHI✔️❌:                                            (b) may be donors of H
    // INCHI✔️❌:                                            therefore may exchange isotopic-non-isotopic H */
    // INCHI✔️❌:                  (bHasH &&
    // INCHI✔️❌:                  ( 0 == ( s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype ) ) ||
    // INCHI✔️❌:                    /* -C=O or =C-OH, O = S, Se, Te */
    // INCHI✔️❌:
    // INCHI✔️❌:                    /*(t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT) &&*/
    // INCHI✔️❌:                    1 == ( s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1/* bAccept_O*/ ) ) ||
    // INCHI✔️❌:                    /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above. M may be N */
    // INCHI✔️❌:
    // INCHI✔️❌:                    2 == ( s_type = GetOtherSaltType( at, i, &s_subtype ) ) ||
    // INCHI✔️❌:                    /* >C-SH, >C-S(-); S=S,Se,Te */
    // INCHI✔️❌:
    // INCHI✔️❌:                    /* other proton donor or acceptor */
    // INCHI✔️❌:                      (bHasAcidicHydrogen(at, i) && ((s_type = 3), (s_subtype = SALT_p_DONOR))) ||
    // INCHI✔️❌:                      (bHasAcidicMinus(at, i) && ((s_type = 3), (s_subtype = SALT_p_ACCEPTOR))) ||
    // INCHI✔️❌:                      (bHasOtherExchangableH(at, i) && ((s_type = 3), (s_subtype = SALT_DONOR_H)))))
    // INCHI✔️❌:
    // INCHI✔️❌:                      )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nNumCandidates >= nMaxNumCandidates)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 s_candidate[nNumCandidates].atnumber = i;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].type = s_type;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].subtype = s_subtype;
    // INCHI✔️❌:                 s_candidate[nNumCandidates].endpoint = at[i].endpoint;
    // INCHI✔️❌:                 nNumCandidates++;
    // INCHI✔️❌:                 nNumNonTautCandidates += !at[i].endpoint;
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nError)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return nError;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumCandidates > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t_group_info->nIsotopicEndpointAtomNumber = (AT_NUMB *) inchi_calloc( (long long)nNumNonTautCandidates + 1, sizeof( t_group_info->nIsotopicEndpointAtomNumber[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:             if (t_group_info->nIsotopicEndpointAtomNumber) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 t_group_info->nIsotopicEndpointAtomNumber[0] = nNumNonTautCandidates;
    // INCHI✔️❌:                 for (i = 0, n = 1; i < nNumCandidates; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     k = s_candidate[i].atnumber;
    // INCHI✔️❌:                     if (!at[k].endpoint)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         t_group_info->nIsotopicEndpointAtomNumber[n++] = k;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     for (j = 0; j < NUM_H_ISOTOPES; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         t_group_info->num_iso_H[j] += at[k].num_iso_H[j];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     at[k].cFlags |= AT_FLAG_ISO_H_POINT;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 t_group_info->nNumIsotopicEndpoints = nNumNonTautCandidates + 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumCandidates;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeIsotopicHGroup
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MakeIsotopicHGroup
    // INCHI✔️❌: typedef unsigned short AT_NUMB;
    // INCHI✔️❌: typedef signed short NUM_H;
    // INCHI✔️❌: typedef signed char S_CHAR;
    // INCHI✔️❌: #define NUM_H_ISOTOPES        3  /* number of hydrogen isotopes: protium, deuterium, tritium */
    // INCHI✔️❌: #define AT_FLAG_ISO_H_POINT 0x01 /**< may have isotopic H */
    // INCHI✔️❌: #define BNS_ERR            -9999
    // INCHI✔️❌: #define BNS_PROGRAM_ERR    (BNS_ERR +  2) /*(-9997)*/
    // INCHI✔️❌: #define BNS_VERT_EDGE_OVFL (BNS_ERR +  6) /*(-9993)*/
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MakeIsotopicHGroup

    let Some(s_group_info) = s_group_info else {
        return Ok(0);
    };
    let Some(t_group_info) = t_group_info else {
        return Ok(0);
    };
    if s_group_info.s_candidate.is_null() || t_group_info.t_group.is_null() {
        return Ok(0);
    }

    let n_max_num_candidates = s_group_info.max_num_candidates;
    t_group_info.num_iso_H = [0; NUM_H_ISOTOPES as usize];
    let atoms = heap.slice(at.as_const())?.to_vec();
    let mut n_num_candidates = 0_i32;
    let mut n_num_non_taut_candidates = 0_i32;
    let mut n_error = 0_i32;

    for i in 0..num_atoms {
        let atom_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = atoms
            .get(atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut s_subtype = 0_i32;
        let mut s_type = 0_i32;
        let b_has_h = if atom.endpoint != 0 {
            let tg_number = *heap
                .slice(t_group_info.tGroupNumber.as_const())?
                .get(usize::from(atom.endpoint))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if tg_number == 0 {
                n_error = BNS_PROGRAM_ERR;
                break;
            }
            let tg = usize::from(tg_number.wrapping_sub(1));
            let group = heap
                .slice(t_group_info.t_group.as_const())?
                .get(tg)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom.endpoint != group.nGroupNumber {
                n_error = BNS_PROGRAM_ERR;
                break;
            }
            i32::from(group.num[0]) - i32::from(group.num[1])
        } else {
            i32::from(atom.num_H)
        };

        let mut accepted = b_has_h != 0 && atom.endpoint != 0;
        if !accepted && b_has_h != 0 {
            s_type = GetSaltChargeType(heap, &atoms, i, Some(&*t_group_info), &mut s_subtype)?;
            accepted = s_type == 0;
            if !accepted {
                s_type = GetOtherSaltChargeType(
                    heap,
                    &atoms,
                    i,
                    Some(&*t_group_info),
                    &mut s_subtype,
                    1,
                )?;
                accepted = s_type == 1;
            }
            if !accepted {
                s_type = GetOtherSaltType(&atoms, i, &mut s_subtype)?;
                accepted = s_type == 2;
            }
            if !accepted && bHasAcidicHydrogen(&atoms, i)? != 0 {
                s_type = 3;
                s_subtype = SALT_p_DONOR as i32;
                accepted = true;
            }
            if !accepted && bHasAcidicMinus(&atoms, i)? != 0 {
                s_type = 3;
                s_subtype = SALT_p_ACCEPTOR as i32;
                accepted = true;
            }
            if !accepted && bHasOtherExchangableH(&atoms, i)? != 0 {
                s_type = 3;
                s_subtype = SALT_DONOR_H as i32;
                accepted = true;
            }
        }
        if !accepted {
            continue;
        }
        if n_num_candidates >= n_max_num_candidates {
            return Ok(BNS_VERT_EDGE_OVFL);
        }
        let candidate_index =
            usize::try_from(n_num_candidates).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(s_group_info.s_candidate)?
            .get_mut(candidate_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = S_CANDIDATE {
            atnumber: i as AT_NUMB,
            type_: s_type as S_CHAR,
            subtype: s_subtype as S_CHAR,
            endpoint: atom.endpoint,
        };
        n_num_candidates = n_num_candidates.wrapping_add(1);
        n_num_non_taut_candidates =
            n_num_non_taut_candidates.wrapping_add(i32::from(atom.endpoint == 0));
    }

    if n_error != 0 {
        return Ok(n_error);
    }

    if n_num_candidates > 0 {
        let allocation_count = u64::try_from(i64::from(n_num_non_taut_candidates) + 1)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let isotopic_endpoints = match inchi_calloc::<AT_NUMB>(
            heap,
            allocation_count,
            std::mem::size_of::<AT_NUMB>() as u64,
        ) {
            Ok(pointer) => Some(pointer),
            Err(
                SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange,
            ) => None,
            Err(error) => return Err(error),
        };
        t_group_info.nIsotopicEndpointAtomNumber = isotopic_endpoints.unwrap_or_default();

        if let Some(isotopic_endpoints) = isotopic_endpoints {
            heap.slice_mut(isotopic_endpoints)?[0] = n_num_non_taut_candidates as AT_NUMB;
            let candidates = heap
                .slice(s_group_info.s_candidate.as_const())?
                .get(
                    ..usize::try_from(n_num_candidates)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let mut n = 1_i32;
            for candidate in candidates {
                let k = usize::from(candidate.atnumber);
                let atom = heap
                    .slice(at.as_const())?
                    .get(k)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                if atom.endpoint == 0 {
                    *heap
                        .slice_mut(isotopic_endpoints)?
                        .get_mut(
                            usize::try_from(n).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = candidate.atnumber;
                    n = n.wrapping_add(1);
                }
                for j in 0..NUM_H_ISOTOPES as usize {
                    t_group_info.num_iso_H[j] =
                        t_group_info.num_iso_H[j].wrapping_add(i16::from(atom.num_iso_H[j]));
                }
                let atom = heap
                    .slice_mut(at)?
                    .get_mut(k)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                atom.cFlags |= AT_FLAG_ISO_H_POINT as S_CHAR;
            }
            t_group_info.nNumIsotopicEndpoints = n_num_non_taut_candidates.wrapping_add(1);
        }
    }

    Ok(n_num_candidates)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RegisterCPoints(
    c_group: &mut [C_GROUP],
    pnum_c: &mut i32,
    max_num_c: i32,
    t_group_info: Option<&T_GROUP_INFO>,
    mut point1: i32,
    mut point2: i32,
    ctype: i32,
    at: &mut [inp_ATOM],
    num_atoms: i32,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2249 RegisterCPoints
    // INCHI✔️❌: int RegisterCPoints( C_GROUP *c_group,
    // INCHI✔️❌:                      int *pnum_c,
    // INCHI✔️❌:                      int max_num_c,
    // INCHI✔️❌:                      T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                      int point1,
    // INCHI✔️❌:                      int point2,
    // INCHI✔️❌:                      int ctype,
    // INCHI✔️❌:                      inp_ATOM *at,
    // INCHI✔️❌:                      int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int num_c = *pnum_c, i, i1, i2;
    // INCHI✔️❌:     AT_NUMB nGroupNumber = 0, nNewGroupNumber;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[point1].c_point == at[point2].c_point)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[point1].c_point)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         memset( c_group + num_c, 0, sizeof( c_group[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (num_c < max_num_c)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             c_group[num_c].num[0] = CHARGED_CPOINT( at, point1 ) + CHARGED_CPOINT( at, point2 );
    // INCHI✔️❌:             c_group[num_c].num_CPoints += 2;
    // INCHI✔️❌:             c_group[num_c].cGroupType = ctype;
    // INCHI✔️❌:             /* get next available c-group number */
    // INCHI✔️❌:             for (i = 0; i < num_c; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nGroupNumber < c_group[i].nGroupNumber)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nGroupNumber = c_group[i].nGroupNumber;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nGroupNumber++;
    // INCHI✔️❌:             c_group[num_c].nGroupNumber =
    // INCHI✔️❌:                 at[point1].c_point =
    // INCHI✔️❌:                 at[point2].c_point = nGroupNumber;
    // INCHI✔️❌:             *pnum_c = num_c + 1;
    // INCHI✔️❌:             /* count protons */
    // INCHI✔️❌:             if (at[point1].num_H)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 c_group[num_c].num[1] ++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[point2].num_H)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     c_group[num_c].num[1] ++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (( at[point1].endpoint || at[point2].endpoint ) && t_group_info && t_group_info->t_group && t_group_info->num_t_groups)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* !!! add later !!! */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return BNS_CPOINT_ERR; /* overflow */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[point1].c_point > at[point2].c_point)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* make sure at[point1].c_point < at[point2].c_point */
    // INCHI✔️❌:         i = point1;
    // INCHI✔️❌:         point1 = point2;
    // INCHI✔️❌:         point2 = i;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!at[point1].c_point)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* add a new c-endpoint to an existing c-group */
    // INCHI✔️❌:         nGroupNumber = at[point2].c_point;
    // INCHI✔️❌:         for (i = 0; i < num_c; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nGroupNumber == c_group[i].nGroupNumber)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[point1].c_point = at[point2].c_point;
    // INCHI✔️❌:                 c_group[i].num_CPoints++;
    // INCHI✔️❌:                 c_group[i].num[0] += CHARGED_CPOINT( at, point1 );
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return BNS_CPOINT_ERR; /* program error: c-group not found */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* merge two c-groups */
    // INCHI✔️❌:         nNewGroupNumber = at[point1].c_point;
    // INCHI✔️❌:         nGroupNumber = at[point2].c_point;
    // INCHI✔️❌:         for (i = 0, i1 = i2 = -1; i < num_c && ( i1 < 0 || i2 < 0 ); i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nNewGroupNumber == c_group[i].nGroupNumber)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 i1 = i;
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nGroupNumber == c_group[i].nGroupNumber)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 i2 = i;
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i1 < 0 || i2 < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return BNS_CPOINT_ERR; /* at least one not found */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         c_group[i1].num[0] += c_group[i2].num[0];
    // INCHI✔️❌:         c_group[i1].num_CPoints += c_group[i2].num_CPoints;
    // INCHI✔️❌:         num_c--;
    // INCHI✔️❌:         if (num_c > i2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memmove(c_group + i2, c_group + i2 + 1, ((long long)num_c - i2) * sizeof(c_group[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         *pnum_c = num_c;
    // INCHI✔️❌:         /* renumber c-groups */
    // INCHI✔️❌:         for (i = 0; i < num_c; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (c_group[i].nGroupNumber > nGroupNumber)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 c_group[i].nGroupNumber--;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* renumber c-points */
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[i].c_point > nGroupNumber)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[i].c_point--;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[i].c_point == nGroupNumber)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at[i].c_point = nNewGroupNumber;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: RegisterCPoints
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: RegisterCPoints
    // INCHI✔️❌: #define CHARGED_CPOINT(X,i) ((X)[i].charge==1)
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: RegisterCPoints

    let mut num_c = *pnum_c;
    let point1_index = point1 as usize;
    let point2_index = point2 as usize;

    if at[point1_index].c_point == at[point2_index].c_point {
        if at[point1_index].c_point != 0 {
            return 0;
        }
        c_group[num_c as usize] = C_GROUP::default();
        if num_c < max_num_c {
            c_group[num_c as usize].num[0] =
                charged_cpoint(at, point1).wrapping_add(charged_cpoint(at, point2));
            c_group[num_c as usize].num_CPoints =
                c_group[num_c as usize].num_CPoints.wrapping_add(2);
            c_group[num_c as usize].cGroupType = ctype as u8;

            let mut nGroupNumber: AT_NUMB = 0;
            for i in 0..num_c as usize {
                if nGroupNumber < c_group[i].nGroupNumber {
                    nGroupNumber = c_group[i].nGroupNumber;
                }
            }
            nGroupNumber = nGroupNumber.wrapping_add(1);
            c_group[num_c as usize].nGroupNumber = nGroupNumber;
            at[point1_index].c_point = nGroupNumber;
            at[point2_index].c_point = nGroupNumber;
            *pnum_c = num_c.wrapping_add(1);

            if at[point1_index].num_H != 0 {
                c_group[num_c as usize].num[1] = c_group[num_c as usize].num[1].wrapping_add(1);
            } else if at[point2_index].num_H != 0 {
                c_group[num_c as usize].num[1] = c_group[num_c as usize].num[1].wrapping_add(1);
            } else if (at[point1_index].endpoint != 0 || at[point2_index].endpoint != 0)
                && t_group_info
                    .is_some_and(|info| !info.t_group.is_null() && info.num_t_groups != 0)
            {
            }
            return 1;
        }

        return BNS_CPOINT_ERR;
    }

    if at[point1 as usize].c_point > at[point2 as usize].c_point {
        std::mem::swap(&mut point1, &mut point2);
    }

    if at[point1 as usize].c_point == 0 {
        let nGroupNumber = at[point2 as usize].c_point;
        for i in 0..num_c as usize {
            if nGroupNumber == c_group[i].nGroupNumber {
                at[point1 as usize].c_point = at[point2 as usize].c_point;
                c_group[i].num_CPoints = c_group[i].num_CPoints.wrapping_add(1);
                c_group[i].num[0] = c_group[i].num[0].wrapping_add(charged_cpoint(at, point1));
                return 1;
            }
        }
        return BNS_CPOINT_ERR;
    }

    let nNewGroupNumber = at[point1 as usize].c_point;
    let nGroupNumber = at[point2 as usize].c_point;
    let mut i1 = -1_i32;
    let mut i2 = -1_i32;
    let mut i = 0_i32;
    while i < num_c && (i1 < 0 || i2 < 0) {
        if nNewGroupNumber == c_group[i as usize].nGroupNumber {
            i1 = i;
            i += 1;
            continue;
        }
        if nGroupNumber == c_group[i as usize].nGroupNumber {
            i2 = i;
            i += 1;
            continue;
        }
        i += 1;
    }
    if i1 < 0 || i2 < 0 {
        return BNS_CPOINT_ERR;
    }

    let i1_index = i1 as usize;
    let i2_index = i2 as usize;
    c_group[i1_index].num[0] = c_group[i1_index].num[0].wrapping_add(c_group[i2_index].num[0]);
    c_group[i1_index].num_CPoints = c_group[i1_index]
        .num_CPoints
        .wrapping_add(c_group[i2_index].num_CPoints);
    num_c = num_c.wrapping_sub(1);
    if num_c > i2 {
        for index in i2_index..num_c as usize {
            c_group[index] = c_group[index + 1].clone();
        }
    }
    *pnum_c = num_c;

    for group in c_group.iter_mut().take(num_c as usize) {
        if group.nGroupNumber > nGroupNumber {
            group.nGroupNumber = group.nGroupNumber.wrapping_sub(1);
        }
    }
    for atom in at.iter_mut().take(num_atoms as usize) {
        if atom.c_point > nGroupNumber {
            atom.c_point = atom.c_point.wrapping_sub(1);
        } else if atom.c_point == nGroupNumber {
            atom.c_point = nNewGroupNumber;
        }
    }
    1
}

#[allow(non_snake_case)]
pub(crate) fn CmpCCandidates(a1: &C_CANDIDATE, a2: &C_CANDIDATE) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2229 CmpCCandidates
    // INCHI✔️❌: int CmpCCandidates( const void *a1, const void *a2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     const C_CANDIDATE *c1 = (const C_CANDIDATE *) a1;
    // INCHI✔️❌:     const C_CANDIDATE *c2 = (const C_CANDIDATE *) a2;
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     if ((ret = (int) c1->type - (int) c2->type)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((ret = (int) c1->subtype - (int) c2->subtype)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     ret = (int) c1->atnumber - (int) c2->atnumber;
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CmpCCandidates

    let ret = i32::from(a1.type_) - i32::from(a2.type_);
    if ret != 0 {
        return ret;
    }
    let ret = i32::from(a1.subtype) - i32::from(a2.subtype);
    if ret != 0 {
        return ret;
    }
    i32::from(a1.atnumber) - i32::from(a2.atnumber)
}

fn ichitaut_is_bns_error(value: i32) -> bool {
    (BNS_ERR..=BNS_MAX_ERR_VALUE).contains(&value)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MarkChargeGroups(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    c_group_info: Option<&mut C_GROUP_INFO>,
    t_group_info: Option<&T_GROUP_INFO>,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2397 MarkChargeGroups
    // INCHI✔️❌: int MarkChargeGroups( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                       inp_ATOM *at,
    // INCHI✔️❌:                       int num_atoms,
    // INCHI✔️❌:                       C_GROUP_INFO *c_group_info,
    // INCHI✔️❌:                       T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                       struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                       struct BalancedNetworkData *pBD )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nNumChanges = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (c_group_info && c_group_info->c_candidate && c_group_info->max_num_candidates > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i, i1, i2, i3, j, num_tested;
    // INCHI✔️❌:         C_CANDIDATE *c_candidate = c_group_info->c_candidate;
    // INCHI✔️❌:         int          nMaxNumCandidates = c_group_info->max_num_candidates;
    // INCHI✔️❌:         int          nNumCandidates = c_group_info->num_candidates;
    // INCHI✔️❌:         S_CHAR       c_type, c_subtype;
    // INCHI✔️❌:         int          iat1, iat2, ret, nDelta;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumCandidates == -1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nNumCandidates = 0; /* 2004-02-26 they could appear after t-group discovery */
    // INCHI✔️❌:                                 /*return 0;*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nNumCandidates == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0, nNumCandidates = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (0 <= ( c_type = GetChargeType( at, i, &c_subtype ) ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (nNumCandidates >= nMaxNumCandidates)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return BNS_VERT_EDGE_OVFL;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     c_candidate[nNumCandidates].atnumber = i;
    // INCHI✔️❌:                     c_candidate[nNumCandidates].type = c_type;
    // INCHI✔️❌:                     c_candidate[nNumCandidates].subtype = c_subtype;
    // INCHI✔️❌:                     nNumCandidates++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nNumCandidates <= 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 c_group_info->num_candidates = -1; /* no candidate exists */
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         qsort( c_candidate, nNumCandidates, sizeof( c_candidate[0] ), CmpCCandidates );
    // INCHI✔️❌:
    // INCHI✔️❌:         i1 = 0;
    // INCHI✔️❌:         num_tested = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         while (i1 < nNumCandidates)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (; i1 < nNumCandidates && ( c_candidate[i1].subtype & C_SUBTYPE_NEUTRAL ); i1++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i1 == nNumCandidates)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break; /* not found */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             for (i2 = i1 + 1; i2 < nNumCandidates &&
    // INCHI✔️❌:                   c_candidate[i2].type == c_candidate[i1].type &&
    // INCHI✔️❌:                   !( c_candidate[i2].subtype & C_SUBTYPE_NEUTRAL );
    // INCHI✔️❌:                   i2++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i2 == nNumCandidates)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break; /* no neutral candidates */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             for (i3 = i2;  i3 < nNumCandidates &&
    // INCHI✔️❌:                   c_candidate[i3].type == c_candidate[i1].type;
    // INCHI✔️❌:                   i3++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (i3 == i2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2 < nNumCandidates)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     i1 = i3;
    // INCHI✔️❌:                     continue;  /* move to the next atom type */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 break; /* nothing more to do */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             for (i = i1; i < i2; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 iat1 = c_candidate[i].atnumber;
    // INCHI✔️❌:                 for (j = i2; j < i3; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_tested++;
    // INCHI✔️❌:                     iat2 = c_candidate[j].atnumber;
    // INCHI✔️❌:                     if (at[iat1].c_point && at[iat1].c_point == at[iat2].c_point)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     ret = bExistsAltPath( pCG, pBNS, pBD, NULL, at, num_atoms, iat1, iat2, ALT_PATH_MODE_CHARGE );
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return ret;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (ret & 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nDelta = ( ret & ~3 ) >> 2;
    // INCHI✔️❌:                         nNumChanges += ( ret & 2 );
    // INCHI✔️❌:
    // INCHI✔️❌:                         ret = RegisterCPoints( c_group_info->c_group, &c_group_info->num_c_groups,
    // INCHI✔️❌:                                                c_group_info->max_num_c_groups, t_group_info,
    // INCHI✔️❌:                                                iat1, iat2, c_candidate[i1].type, at, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return ret;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nDelta)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto quick_exit;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             i1 = i3;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     quick_exit:
    // INCHI✔️❌:         if (c_group_info->num_candidates == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             c_group_info->num_candidates = num_tested ? nNumCandidates : -1; /* no candidate exists */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumChanges;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkChargeGroups
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MarkChargeGroups
    // INCHI✔️❌: #define IS_BNS_ERROR(X) (BNS_ERR <= (X) && (X) <= BNS_MAX_ERR_VALUE)
    // INCHI✔️❌: #define ALT_PATH_MODE_CHARGE 2
    // END INCHI ACTIVE MACRO CONFIGURATION: MarkChargeGroups

    let mut n_num_changes = 0_i32;
    let Some(c_group_info) = c_group_info else {
        return Ok(n_num_changes);
    };
    if c_group_info.c_candidate.is_null() || c_group_info.max_num_candidates <= 0 {
        return Ok(n_num_changes);
    }

    let c_candidate = c_group_info.c_candidate;
    let n_max_num_candidates = c_group_info.max_num_candidates;
    let mut n_num_candidates = c_group_info.num_candidates;
    let mut num_tested = 0_i32;

    if n_num_candidates == -1 {
        n_num_candidates = 0;
    }
    if n_num_candidates == 0 {
        n_num_candidates = 0;
        for i in 0..num_atoms {
            let mut c_subtype = 0 as S_CHAR;
            let c_type = {
                let atoms = heap.slice(at.as_const())?;
                GetChargeType(atoms, i, &mut c_subtype)
            };
            if c_type >= 0 {
                if n_num_candidates >= n_max_num_candidates {
                    return Ok(BNS_VERT_EDGE_OVFL);
                }
                let candidate = heap
                    .slice_mut(c_candidate)?
                    .get_mut(n_num_candidates as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                candidate.atnumber = i as AT_NUMB;
                candidate.type_ = c_type as S_CHAR;
                candidate.subtype = c_subtype;
                n_num_candidates += 1;
            }
        }
        if n_num_candidates <= 1 {
            c_group_info.num_candidates = -1;
            return Ok(0);
        }
    }

    let candidate_count =
        usize::try_from(n_num_candidates).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut candidates = heap
        .slice(c_candidate.as_const())?
        .get(..candidate_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    candidates.sort_by(|left, right| CmpCCandidates(left, right).cmp(&0));
    heap.slice_mut(c_candidate)?
        .get_mut(..candidate_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone_from_slice(&candidates);

    let mut i1 = 0_i32;
    'quick_exit: while i1 < n_num_candidates {
        while i1 < n_num_candidates
            && (candidates[i1 as usize].subtype & C_SUBTYPE_NEUTRAL as S_CHAR) != 0
        {
            i1 += 1;
        }
        if i1 == n_num_candidates {
            break;
        }

        let mut i2 = i1 + 1;
        while i2 < n_num_candidates
            && candidates[i2 as usize].type_ == candidates[i1 as usize].type_
            && (candidates[i2 as usize].subtype & C_SUBTYPE_NEUTRAL as S_CHAR) == 0
        {
            i2 += 1;
        }
        if i2 == n_num_candidates {
            break;
        }

        let mut i3 = i2;
        while i3 < n_num_candidates
            && candidates[i3 as usize].type_ == candidates[i1 as usize].type_
        {
            i3 += 1;
        }

        if i3 == i2 {
            if i2 < n_num_candidates {
                i1 = i3;
                continue;
            }
            break;
        }

        for i in i1..i2 {
            let iat1 = i32::from(candidates[i as usize].atnumber);
            for j in i2..i3 {
                num_tested += 1;
                let iat2 = i32::from(candidates[j as usize].atnumber);
                {
                    let atoms = heap.slice(at.as_const())?;
                    if atoms[iat1 as usize].c_point != 0
                        && atoms[iat1 as usize].c_point == atoms[iat2 as usize].c_point
                    {
                        continue;
                    }
                }

                let ret = bExistsAltPath(
                    heap,
                    pCG,
                    pBNS,
                    pBD,
                    None,
                    at,
                    num_atoms,
                    iat1,
                    iat2,
                    ALT_PATH_MODE_CHARGE as i32,
                    clock_result,
                )?;
                if ichitaut_is_bns_error(ret) {
                    return Ok(ret);
                }
                if (ret & 1) != 0 {
                    let n_delta = (ret & !3) >> 2;
                    n_num_changes += ret & 2;
                    let mut c_groups = heap.slice(c_group_info.c_group.as_const())?.to_vec();
                    let mut atoms = heap.slice(at.as_const())?.to_vec();
                    let ret = RegisterCPoints(
                        &mut c_groups,
                        &mut c_group_info.num_c_groups,
                        c_group_info.max_num_c_groups,
                        t_group_info,
                        iat1,
                        iat2,
                        i32::from(candidates[i1 as usize].type_),
                        &mut atoms,
                        num_atoms,
                    );
                    heap.slice_mut(c_group_info.c_group)?
                        .clone_from_slice(&c_groups);
                    heap.slice_mut(at)?.clone_from_slice(&atoms);
                    if ichitaut_is_bns_error(ret) {
                        return Ok(ret);
                    }
                    if n_delta != 0 {
                        break 'quick_exit;
                    }
                }
            }
        }
        i1 = i3;
    }

    if c_group_info.num_candidates == 0 {
        c_group_info.num_candidates = if num_tested != 0 {
            n_num_candidates
        } else {
            -1
        };
    }

    Ok(n_num_changes)
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo(atom: &[inp_ATOM], iat: i32, eif: &mut ENDPOINT_INFO) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:359 nGetEndpointInfo
    // INCHI✔️✔️: int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (!( nEndpointValence = get_endpoint_valence( atom[iat].el_number ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* not an endpoint */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nEndpointValence <= atom[iat].valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         if (nEndpointValence < atom[iat].chem_bonds_valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0; /* abnormal valence > standard endpoint valence */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         nMobile = atom[iat].num_H + ( atom[iat].charge == -1 );
    // INCHI✔️✔️:         if (nMobile + atom[iat].chem_bonds_valence != nEndpointValence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0; /* non-standard endpoint valence */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         switch (atom[iat].chem_bonds_valence - atom[iat].valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             case 0:
    // INCHI✔️✔️:                 eif->cDonor = 1;
    // INCHI✔️✔️:                 eif->cAcceptor = 0;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             case 1:
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             default:
    // INCHI✔️✔️:                 return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:              0 <= GetChargeType( atom, iat, &cChargeSubtype ) &&
    // INCHI✔️✔️:              ( (int) cChargeSubtype & ( C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_ACCEPT)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     eif->cDonor = 0;
    // INCHI✔️✔️:                     eif->cAcceptor = 1;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (cChargeSubtype & C_SUBTYPE_H_DONOR)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         eif->cDonor = 1;
    // INCHI✔️✔️:                         eif->cAcceptor = 0;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         return 0;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:             eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:             eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:             eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:             return nEndpointValence;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = get_endpoint_valence(at.el_number);
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(at.valence) {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        if n_endpoint_valence < i32::from(at.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(at.num_H) + i32::from(at.charge == -1);
        if n_mobile + i32::from(at.chem_bonds_valence) != n_endpoint_valence {
            return 0;
        }
        match i32::from(at.chem_bonds_valence) - i32::from(at.valence) {
            0 => {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            }
            1 => {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            }
            _ => return 0,
        }
        eif.cMobile = n_mobile as S_CHAR;
        eif.cNeutralBondsValence = (n_endpoint_valence - n_mobile) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            } else {
                return 0;
            }
            eif.cMobile = at.num_H;
            eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
            eif.cMoveableCharge = at.charge;
            eif.cKetoEnolCode = 0;
            return n_endpoint_valence;
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo_PT_22_00(
    atom: &[inp_ATOM],
    iat: i32,
    eif: &mut ENDPOINT_INFO,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:452 nGetEndpointInfo_PT_22_00
    // INCHI✔️✔️: int nGetEndpointInfo_PT_22_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     nEndpointValence = atom[iat].el_number == EL_NUMBER_C ? 4 : 0;
    // INCHI✔️✔️:     if (!nEndpointValence)
    // INCHI✔️✔️:         return 0; /* not an endpoint */
    // INCHI✔️✔️:     if (nEndpointValence <= atom[iat].valence)
    // INCHI✔️✔️:         return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0) {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         if (nEndpointValence < atom[iat].chem_bonds_valence)
    // INCHI✔️✔️:             return 0; /* abnormal valence > standard endpoint valence */
    // INCHI✔️✔️:         nMobile = atom[iat].num_H + (atom[iat].charge == -1);
    // INCHI✔️✔️:         if (nMobile + atom[iat].chem_bonds_valence != nEndpointValence)
    // INCHI✔️✔️:             return 0; /* non-standard endpoint valence */
    // INCHI✔️✔️:         switch (atom[iat].chem_bonds_valence - atom[iat].valence) {
    // INCHI✔️✔️:         case 0:
    // INCHI✔️✔️:             eif->cDonor = 1;
    // INCHI✔️✔️:             eif->cAcceptor = 0;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case 1:
    // INCHI✔️✔️:             eif->cDonor = 0;
    // INCHI✔️✔️:             eif->cAcceptor = 1;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:             0 <= GetChargeType(atom, iat, &cChargeSubtype) &&
    // INCHI✔️✔️:             ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR))
    // INCHI✔️✔️:             ) {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             if (cChargeSubtype & C_SUBTYPE_H_ACCEPT) {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_DONOR) {
    // INCHI✔️✔️:                     eif->cDonor = 1;
    // INCHI✔️✔️:                     eif->cAcceptor = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else {
    // INCHI✔️✔️:                     return 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:                 eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:                 return nEndpointValence;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo_PT_22_00
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_22_00
    // INCHI✔️✔️: #define TAUT_PT_22_00              1 /* tautomerism rule PT_22_00 */
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_22_00

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = if at.el_number == 6 { 4 } else { 0 };
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(at.valence) {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        if n_endpoint_valence < i32::from(at.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(at.num_H) + i32::from(at.charge == -1);
        if n_mobile + i32::from(at.chem_bonds_valence) != n_endpoint_valence {
            return 0;
        }
        match i32::from(at.chem_bonds_valence) - i32::from(at.valence) {
            0 => {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            }
            1 => {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            }
            _ => return 0,
        }
        eif.cMobile = n_mobile as S_CHAR;
        eif.cNeutralBondsValence = (n_endpoint_valence - n_mobile) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            } else {
                return 0;
            }
            eif.cMobile = at.num_H;
            eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
            eif.cMoveableCharge = at.charge;
            eif.cKetoEnolCode = 0;
            return n_endpoint_valence;
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo_PT_16_00(
    atom: &[inp_ATOM],
    iat: i32,
    eif: &mut ENDPOINT_INFO,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:524 nGetEndpointInfo_PT_16_00
    // INCHI✔️✔️: int nGetEndpointInfo_PT_16_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     nEndpointValence = get_endpoint_valence_KET( atom[iat].el_number );
    // INCHI✔️✔️:     if (!nEndpointValence)
    // INCHI✔️✔️:         return 0; /* not an endpoint */
    // INCHI✔️✔️:     if (nEndpointValence <= atom[iat].valence)
    // INCHI✔️✔️:         return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    // INCHI✔️✔️:     if (nEndpointValence == 4 && atom[iat].valence < 2)
    // INCHI✔️✔️:         return 0; /* exclude O==N--CH3  <=> HO--N==CH2 */
    // INCHI✔️✔️:     if (nEndpointValence == 2 && atom[iat].valence > 1)
    // INCHI✔️✔️:         return 0; /* exclude --O--N==CH-- */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0) {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         if (nEndpointValence < atom[iat].chem_bonds_valence)
    // INCHI✔️✔️:             return 0; /* abnormal valence > standard endpoint valence */
    // INCHI✔️✔️:         nMobile = atom[iat].num_H + (atom[iat].charge == -1);
    // INCHI✔️✔️:         if (nMobile + atom[iat].chem_bonds_valence != nEndpointValence)
    // INCHI✔️✔️:             return 0; /* non-standard endpoint valence */
    // INCHI✔️✔️:         switch (atom[iat].chem_bonds_valence - atom[iat].valence) {
    // INCHI✔️✔️:         case 0:
    // INCHI✔️✔️:             eif->cDonor = 1;
    // INCHI✔️✔️:             eif->cAcceptor = 0;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case 1:
    // INCHI✔️✔️:             eif->cDonor = 0;
    // INCHI✔️✔️:             eif->cAcceptor = 1;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:             0 <= GetChargeType(atom, iat, &cChargeSubtype) &&
    // INCHI✔️✔️:             ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR))
    // INCHI✔️✔️:             ) {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             if (cChargeSubtype & C_SUBTYPE_H_ACCEPT) {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_DONOR) {
    // INCHI✔️✔️:                     eif->cDonor = 1;
    // INCHI✔️✔️:                     eif->cAcceptor = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else {
    // INCHI✔️✔️:                     return 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:                 eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:                 return nEndpointValence;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo_PT_16_00
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_16_00
    // INCHI✔️✔️: #define TAUT_PT_16_00              1 /* tautomerism rule PT_16_00 */
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_16_00

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = get_endpoint_valence_KET(at.el_number);
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(at.valence) {
        return 0;
    }
    if n_endpoint_valence == 4 && i32::from(at.valence) < 2 {
        return 0;
    }
    if n_endpoint_valence == 2 && i32::from(at.valence) > 1 {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        if n_endpoint_valence < i32::from(at.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(at.num_H) + i32::from(at.charge == -1);
        if n_mobile + i32::from(at.chem_bonds_valence) != n_endpoint_valence {
            return 0;
        }
        match i32::from(at.chem_bonds_valence) - i32::from(at.valence) {
            0 => {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            }
            1 => {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            }
            _ => return 0,
        }
        eif.cMobile = n_mobile as S_CHAR;
        eif.cNeutralBondsValence = (n_endpoint_valence - n_mobile) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            } else {
                return 0;
            }
            eif.cMobile = at.num_H;
            eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
            eif.cMoveableCharge = at.charge;
            eif.cKetoEnolCode = 0;
            return n_endpoint_valence;
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo_PT_06_00(
    atom: &[inp_ATOM],
    iat: i32,
    eif: &mut ENDPOINT_INFO,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:600 nGetEndpointInfo_PT_06_00
    // INCHI✔️✔️: int nGetEndpointInfo_PT_06_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     nEndpointValence = atom[iat].el_number == EL_NUMBER_C ? 4 :
    // INCHI✔️✔️:         get_endpoint_valence( atom[iat].el_number );
    // INCHI✔️✔️:     /*printf("Connectivity: %d\n", atom[iat].valence);
    // INCHI✔️✔️:     printf("Charge: %d\n", atom[iat].charge);
    // INCHI✔️✔️:     printf("Actual valence: %d\n", atom[iat].chem_bonds_valence);
    // INCHI✔️✔️:     printf("Num H: %d\n", atom[iat].num_H);*/
    // INCHI✔️✔️:     if (!nEndpointValence)
    // INCHI✔️✔️:         return 0; /* not an endpoint */
    // INCHI✔️✔️:     if (nEndpointValence <= atom[iat].valence)
    // INCHI✔️✔️:         return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0) {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         if (nEndpointValence < atom[iat].chem_bonds_valence)
    // INCHI✔️✔️:             return 0; /* abnormal valence > standard endpoint valence */
    // INCHI✔️✔️:         nMobile = atom[iat].num_H + (atom[iat].charge == -1);
    // INCHI✔️✔️:         if (nMobile + atom[iat].chem_bonds_valence != nEndpointValence)
    // INCHI✔️✔️:             return 0; /* non-standard endpoint valence */
    // INCHI✔️✔️:         switch (atom[iat].chem_bonds_valence - atom[iat].valence) {
    // INCHI✔️✔️:         case 0:
    // INCHI✔️✔️:             eif->cDonor = 1;
    // INCHI✔️✔️:             eif->cAcceptor = 0;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case 1:
    // INCHI✔️✔️:             eif->cDonor = 0;
    // INCHI✔️✔️:             eif->cAcceptor = 1;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:             0 <= GetChargeType(atom, iat, &cChargeSubtype) &&
    // INCHI✔️✔️:             ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR))
    // INCHI✔️✔️:             ) {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             if (cChargeSubtype & C_SUBTYPE_H_ACCEPT) {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_DONOR) {
    // INCHI✔️✔️:                     eif->cDonor = 1;
    // INCHI✔️✔️:                     eif->cAcceptor = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else {
    // INCHI✔️✔️:                     return 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:                 eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:                 return nEndpointValence;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo_PT_06_00
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_06_00
    // INCHI✔️✔️: #define TAUT_PT_06_00              1 /* tautomerism rule PT_06_00 */
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_06_00

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = if at.el_number == 6 {
        4
    } else {
        get_endpoint_valence(at.el_number)
    };
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(at.valence) {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        if n_endpoint_valence < i32::from(at.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(at.num_H) + i32::from(at.charge == -1);
        if n_mobile + i32::from(at.chem_bonds_valence) != n_endpoint_valence {
            return 0;
        }
        match i32::from(at.chem_bonds_valence) - i32::from(at.valence) {
            0 => {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            }
            1 => {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            }
            _ => return 0,
        }
        eif.cMobile = n_mobile as S_CHAR;
        eif.cNeutralBondsValence = (n_endpoint_valence - n_mobile) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            } else {
                return 0;
            }
            eif.cMobile = at.num_H;
            eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
            eif.cMoveableCharge = at.charge;
            eif.cKetoEnolCode = 0;
            return n_endpoint_valence;
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo_PT_39_00(
    atom: &[inp_ATOM],
    iat: i32,
    eif: &mut ENDPOINT_INFO,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:677 nGetEndpointInfo_PT_39_00
    // INCHI✔️✔️: int nGetEndpointInfo_PT_39_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     nEndpointValence = atom[iat].el_number == EL_NUMBER_C ? 4 :
    // INCHI✔️✔️:         atom[iat].el_number == EL_NUMBER_N ? 3 : 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!nEndpointValence)
    // INCHI✔️✔️:         return 0; /* not an endpoint */
    // INCHI✔️✔️:     if (nEndpointValence <= atom[iat].valence)
    // INCHI✔️✔️:         return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    // INCHI✔️✔️:                   /*
    // INCHI✔️✔️:                   if ( nEndpointValence == 4 && !(atom[iat].valence == 3 || atom[iat].valence == 4))
    // INCHI✔️✔️:                   return 0;
    // INCHI✔️✔️:                   if ( nEndpointValence == 3 && !(atom[iat].valence == 2 || atom[iat].valence == 3))
    // INCHI✔️✔️:                   return 0;
    // INCHI✔️✔️:                   */
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0) {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         if (nEndpointValence < atom[iat].chem_bonds_valence)
    // INCHI✔️✔️:             return 0; /* abnormal valence > standard endpoint valence */
    // INCHI✔️✔️:         nMobile = atom[iat].num_H + (atom[iat].charge == -1);
    // INCHI✔️✔️:         if (nMobile + atom[iat].chem_bonds_valence != nEndpointValence)
    // INCHI✔️✔️:             return 0; /* non-standard endpoint valence */
    // INCHI✔️✔️:         switch (atom[iat].chem_bonds_valence - atom[iat].valence) {
    // INCHI✔️✔️:         case 0:
    // INCHI✔️✔️:             eif->cDonor = 1;
    // INCHI✔️✔️:             eif->cAcceptor = 0;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case 1:
    // INCHI✔️✔️:             eif->cDonor = 0;
    // INCHI✔️✔️:             eif->cAcceptor = 1;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:             0 <= GetChargeType(atom, iat, &cChargeSubtype) &&
    // INCHI✔️✔️:             ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR))
    // INCHI✔️✔️:             ) {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             if (cChargeSubtype & C_SUBTYPE_H_ACCEPT) {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_DONOR) {
    // INCHI✔️✔️:                     eif->cDonor = 1;
    // INCHI✔️✔️:                     eif->cAcceptor = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else {
    // INCHI✔️✔️:                     return 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:                 eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:                 return nEndpointValence;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo_PT_39_00
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_39_00
    // INCHI✔️✔️: #define TAUT_PT_39_00              1 /* tautomerism rule PT_39_00 */
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_N  ((U_CHAR)7)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_39_00

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = if at.el_number == 6 {
        4
    } else if at.el_number == 7 {
        3
    } else {
        0
    };
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(at.valence) {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        if n_endpoint_valence < i32::from(at.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(at.num_H) + i32::from(at.charge == -1);
        if n_mobile + i32::from(at.chem_bonds_valence) != n_endpoint_valence {
            return 0;
        }
        match i32::from(at.chem_bonds_valence) - i32::from(at.valence) {
            0 => {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            }
            1 => {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            }
            _ => return 0,
        }
        eif.cMobile = n_mobile as S_CHAR;
        eif.cNeutralBondsValence = (n_endpoint_valence - n_mobile) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            } else {
                return 0;
            }
            eif.cMobile = at.num_H;
            eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
            eif.cMoveableCharge = at.charge;
            eif.cKetoEnolCode = 0;
            return n_endpoint_valence;
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo_PT_13_00(
    atom: &[inp_ATOM],
    iat: i32,
    eif: &mut ENDPOINT_INFO,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:756 nGetEndpointInfo_PT_13_00
    // INCHI✔️✔️: int nGetEndpointInfo_PT_13_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     nEndpointValence = atom[iat].el_number == EL_NUMBER_C ? 4 :
    // INCHI✔️✔️:         atom[iat].el_number == EL_NUMBER_S ? 2 :
    // INCHI✔️✔️:         atom[iat].el_number == EL_NUMBER_O ? 2 :
    // INCHI✔️✔️:         atom[iat].el_number == EL_NUMBER_SE ? 2 :
    // INCHI✔️✔️:         atom[iat].el_number == EL_NUMBER_TE ? 2 : 0;
    // INCHI✔️✔️:     /*printf("Connectivity: %d\n", atom[iat].valence);
    // INCHI✔️✔️:     printf("Charge: %d\n", atom[iat].charge);
    // INCHI✔️✔️:     printf("Actual valence: %d\n", atom[iat].chem_bonds_valence);
    // INCHI✔️✔️:     printf("Num H: %d\n", atom[iat].num_H);*/
    // INCHI✔️✔️:     if (!nEndpointValence)
    // INCHI✔️✔️:         return 0; /* not an endpoint */
    // INCHI✔️✔️:     if (nEndpointValence <= atom[iat].valence)
    // INCHI✔️✔️:         return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0) {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         if (nEndpointValence < atom[iat].chem_bonds_valence)
    // INCHI✔️✔️:             return 0; /* abnormal valence > standard endpoint valence */
    // INCHI✔️✔️:         nMobile = atom[iat].num_H + (atom[iat].charge == -1);
    // INCHI✔️✔️:         if (nMobile + atom[iat].chem_bonds_valence != nEndpointValence)
    // INCHI✔️✔️:             return 0; /* non-standard endpoint valence */
    // INCHI✔️✔️:         if (nMobile > 0) {
    // INCHI✔️✔️:             eif->cDonor = 1;
    // INCHI✔️✔️:             eif->cAcceptor = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else {
    // INCHI✔️✔️:             eif->cDonor = 0;
    // INCHI✔️✔️:             eif->cAcceptor = 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:             0 <= GetChargeType(atom, iat, &cChargeSubtype) &&
    // INCHI✔️✔️:             ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR))
    // INCHI✔️✔️:             ) {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             if (cChargeSubtype & C_SUBTYPE_H_ACCEPT) {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_DONOR) {
    // INCHI✔️✔️:                     eif->cDonor = 1;
    // INCHI✔️✔️:                     eif->cAcceptor = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else {
    // INCHI✔️✔️:                     return 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:                 eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:                 eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:                 return nEndpointValence;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo_PT_13_00
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_13_00
    // INCHI✔️✔️: #define TAUT_PT_13_00              1 /* tautomerism rule PT_13_00 */
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define EL_NUMBER_S  ((U_CHAR)16)
    // INCHI✔️✔️: #define EL_NUMBER_SE ((U_CHAR)34)
    // INCHI✔️✔️: #define EL_NUMBER_TE ((U_CHAR)52)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_13_00

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = match at.el_number {
        6 => 4,
        16 | 8 | 34 | 52 => 2,
        _ => 0,
    };
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(at.valence) {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        if n_endpoint_valence < i32::from(at.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(at.num_H) + i32::from(at.charge == -1);
        if n_mobile + i32::from(at.chem_bonds_valence) != n_endpoint_valence {
            return 0;
        }
        if n_mobile > 0 {
            eif.cDonor = 1;
            eif.cAcceptor = 0;
        } else {
            eif.cDonor = 0;
            eif.cAcceptor = 1;
        }
        eif.cMobile = n_mobile as S_CHAR;
        eif.cNeutralBondsValence = (n_endpoint_valence - n_mobile) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            } else {
                return 0;
            }
            eif.cMobile = at.num_H;
            eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
            eif.cMoveableCharge = at.charge;
            eif.cKetoEnolCode = 0;
            return n_endpoint_valence;
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo_PT_18_00(
    atom: &[inp_ATOM],
    iat: i32,
    eif: &mut ENDPOINT_INFO,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:832 nGetEndpointInfo_PT_18_00
    // INCHI✔️✔️: int nGetEndpointInfo_PT_18_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:     /* int res; removed */
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     nEndpointValence = atom[iat].el_number == EL_NUMBER_O ? 2 :
    // INCHI✔️✔️:         atom[iat].el_number == EL_NUMBER_N ? 3 : 0;
    // INCHI✔️✔️:     /*printf("Connectivity: %d\n", atom[iat].valence);
    // INCHI✔️✔️:     printf("Charge: %d\n", atom[iat].charge);
    // INCHI✔️✔️:     printf("Actual valence: %d\n", atom[iat].chem_bonds_valence);
    // INCHI✔️✔️:     printf("Num H: %d\n", atom[iat].num_H);
    // INCHI✔️✔️:     printf("c-point: %d\n", atom[iat].c_point);
    // INCHI✔️✔️:     res = GetChargeType( atom, iat, &cChargeSubtype );
    // INCHI✔️✔️:     printf("Charge subtype: %d %d\n", res, cChargeSubtype);*/
    // INCHI✔️✔️:     if (!nEndpointValence)
    // INCHI✔️✔️:         return 0; /* not an endpoint */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0) {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         nMobile = atom[iat].num_H;
    // INCHI✔️✔️:         if (nMobile > 0) {
    // INCHI✔️✔️:             eif->cDonor = 1;
    // INCHI✔️✔️:             eif->cAcceptor = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else {
    // INCHI✔️✔️:             eif->cDonor = 0;
    // INCHI✔️✔️:             eif->cAcceptor = 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else {
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:             0 <= GetChargeType(atom, iat, &cChargeSubtype) &&
    // INCHI✔️✔️:             ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR))
    // INCHI✔️✔️:             ) {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             if (cChargeSubtype & C_SUBTYPE_H_ACCEPT) {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_DONOR) {
    // INCHI✔️✔️:                     eif->cDonor = 1;
    // INCHI✔️✔️:                     eif->cAcceptor = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else {
    // INCHI✔️✔️:             /* charge-point */
    // INCHI✔️✔️:             if (atom[iat].num_H) {
    // INCHI✔️✔️:                 eif->cDonor = 1;
    // INCHI✔️✔️:                 eif->cAcceptor = 0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:         eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔️✔️:         eif->cKetoEnolCode = 0;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo_PT_18_00
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_18_00
    // INCHI✔️✔️: #define TAUT_PT_18_00              1 /* tautomerism rule PT_18_00 */
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_N  ((U_CHAR)7)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_PT_18_00

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = if at.el_number == 8 {
        2
    } else if at.el_number == 7 {
        3
    } else {
        0
    };
    if n_endpoint_valence == 0 {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        let n_mobile = at.num_H;
        if n_mobile > 0 {
            eif.cDonor = 1;
            eif.cAcceptor = 0;
        } else {
            eif.cDonor = 0;
            eif.cAcceptor = 1;
        }
        eif.cMobile = n_mobile;
        eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(n_mobile)) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            }
        } else if at.num_H != 0 {
            eif.cDonor = 1;
            eif.cAcceptor = 0;
        } else {
            eif.cDonor = 0;
            eif.cAcceptor = 1;
        }
        eif.cMobile = at.num_H;
        eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
        eif.cMoveableCharge = at.charge;
        eif.cKetoEnolCode = 0;
        return n_endpoint_valence;
    }
}

#[allow(non_snake_case)]
pub(crate) fn nGetEndpointInfo_KET(atom: &[inp_ATOM], iat: i32, eif: &mut ENDPOINT_INFO) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:916 nGetEndpointInfo_KET
    // INCHI✔️✔️: int nGetEndpointInfo_KET( inp_ATOM *atom,
    // INCHI✔️✔️:                           int iat,
    // INCHI✔️✔️:                           ENDPOINT_INFO *eif )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nEndpointValence;
    // INCHI✔️✔️:     int  nMobile;
    // INCHI✔️✔️:     S_CHAR cChargeSubtype;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].radical && atom[iat].radical != RADICAL_SINGLET)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* a radical */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (!( nEndpointValence = get_endpoint_valence_KET( atom[iat].el_number ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* not an endpoint; only O and C can be an endpoint for keto-enol tautomerism */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nEndpointValence <= atom[iat].valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nEndpointValence == 4 && atom[iat].valence < 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* exclude O==C--CH3  <=> HO--C==CH2 */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nEndpointValence == 2 && atom[iat].valence > 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* exclude --O--C==CH-- */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[iat].charge == -1 || atom[iat].charge == 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* not a positive charge-point */
    // INCHI✔️✔️:         if (nEndpointValence < atom[iat].chem_bonds_valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0; /* abnormal valence > standard endpoint valence */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         nMobile = atom[iat].num_H + ( atom[iat].charge == -1 );
    // INCHI✔️✔️:         if (nMobile + atom[iat].chem_bonds_valence != nEndpointValence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0; /* non-standard endpoint valence */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         switch (atom[iat].chem_bonds_valence - atom[iat].valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             case 0:
    // INCHI✔️✔️:                 eif->cDonor = 1;
    // INCHI✔️✔️:                 eif->cAcceptor = 0;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             case 1:
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             default:
    // INCHI✔️✔️:                 return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         eif->cMobile = nMobile;
    // INCHI✔️✔️:         eif->cNeutralBondsValence = nEndpointValence - nMobile;
    // INCHI✔️✔️:         eif->cMoveableCharge = 0;
    // INCHI✔️✔️:         eif->cKetoEnolCode = ( nEndpointValence == 2 ) ? 1 : ( nEndpointValence == 4 ) ? 2 : 0;
    // INCHI✔️✔️:         return nEndpointValence;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (atom[iat].c_point &&
    // INCHI✔️✔️:              0 <= GetChargeType( atom, iat, &cChargeSubtype ) &&
    // INCHI✔️✔️:              ( (int) cChargeSubtype & ( C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* charge-point; currently only O for keto-enol tautomerism */
    // INCHI✔️✔️:             if (cChargeSubtype & C_SUBTYPE_H_ACCEPT)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 eif->cDonor = 0;
    // INCHI✔️✔️:                 eif->cAcceptor = 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (cChargeSubtype & C_SUBTYPE_H_DONOR)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     eif->cDonor = 1;
    // INCHI✔️✔️:                     eif->cAcceptor = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     return 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             eif->cMobile = atom[iat].num_H;
    // INCHI✔️✔️:             eif->cNeutralBondsValence = nEndpointValence - atom[iat].num_H;
    // INCHI✔️✔️:             eif->cMoveableCharge = atom[iat].charge;
    // INCHI✔️✔️:             eif->cKetoEnolCode = ( nEndpointValence == 2 ) ? 1 : ( nEndpointValence == 4 ) ? 2 : 0;
    // INCHI✔️✔️:             return nEndpointValence;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nGetEndpointInfo_KET
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_KET
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nGetEndpointInfo_KET

    let at = &atom[iat as usize];
    if at.radical != 0 && at.radical != RADICAL_SINGLET as S_CHAR {
        return 0;
    }
    let n_endpoint_valence = get_endpoint_valence_KET(at.el_number);
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(at.valence) {
        return 0;
    }
    if n_endpoint_valence == 4 && i32::from(at.valence) < 2 {
        return 0;
    }
    if n_endpoint_valence == 2 && i32::from(at.valence) > 1 {
        return 0;
    }

    if at.charge == -1 || at.charge == 0 {
        if n_endpoint_valence < i32::from(at.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(at.num_H) + i32::from(at.charge == -1);
        if n_mobile + i32::from(at.chem_bonds_valence) != n_endpoint_valence {
            return 0;
        }
        match i32::from(at.chem_bonds_valence) - i32::from(at.valence) {
            0 => {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            }
            1 => {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            }
            _ => return 0,
        }
        eif.cMobile = n_mobile as S_CHAR;
        eif.cNeutralBondsValence = (n_endpoint_valence - n_mobile) as S_CHAR;
        eif.cMoveableCharge = 0;
        eif.cKetoEnolCode = if n_endpoint_valence == 2 {
            1
        } else if n_endpoint_valence == 4 {
            2
        } else {
            0
        };
        return n_endpoint_valence;
    } else {
        let mut c_charge_subtype = 0;
        if at.c_point != 0
            && GetChargeType(atom, iat, &mut c_charge_subtype) >= 0
            && (i32::from(c_charge_subtype)
                & (C_SUBTYPE_H_ACCEPT as i32 | C_SUBTYPE_H_DONOR as i32))
                != 0
        {
            if (i32::from(c_charge_subtype) & C_SUBTYPE_H_ACCEPT as i32) != 0 {
                eif.cDonor = 0;
                eif.cAcceptor = 1;
            } else if (i32::from(c_charge_subtype) & C_SUBTYPE_H_DONOR as i32) != 0 {
                eif.cDonor = 1;
                eif.cAcceptor = 0;
            } else {
                return 0;
            }
            eif.cMobile = at.num_H;
            eif.cNeutralBondsValence = (n_endpoint_valence - i32::from(at.num_H)) as S_CHAR;
            eif.cMoveableCharge = at.charge;
            eif.cKetoEnolCode = if n_endpoint_valence == 2 {
                1
            } else if n_endpoint_valence == 4 {
                2
            } else {
                0
            };
            return n_endpoint_valence;
        }
    }

    0
}

#[allow(non_snake_case)]
pub(crate) fn SetTautomericBonds(
    at: &mut [inp_ATOM],
    nNumBondPos: i32,
    BondPos: &[T_BONDPOS],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1523 SetTautomericBonds
    // INCHI✔️✔️: int SetTautomericBonds( inp_ATOM *at,
    // INCHI✔️✔️:                         int nNumBondPos,
    // INCHI✔️✔️:                         T_BONDPOS *BondPos )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k, n;
    // INCHI✔️✔️:     for (k = n = 0; k < nNumBondPos; k++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int neighbor_index = BondPos[k].neighbor_index;
    // INCHI✔️✔️:         int center = BondPos[k].nAtomNumber;
    // INCHI✔️✔️:         int bond_mark = at[center].bond_type[neighbor_index];
    // INCHI✔️✔️:         int bond_type = bond_mark & ~BOND_MARK_ALL;
    // INCHI✔️✔️:         int neighbor;
    // INCHI✔️✔️: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️✔️:         if (bond_type != BOND_TAUTOM)
    // INCHI❌❌: #else
    // INCHI❌❌:         if (bond_type != BOND_ALTERN && bond_type != BOND_TAUTOM)
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             int ii;
    // INCHI✔️✔️:             /*  change bond type to BOND_TAUTOM presering higher bits marks */
    // INCHI✔️✔️:             bond_type = ( bond_mark & BOND_MARK_ALL ) | BOND_TAUTOM;
    // INCHI✔️✔️:             /*  change center-neighbor bond */
    // INCHI✔️✔️:             at[center].bond_type[neighbor_index] = bond_type;
    // INCHI✔️✔️:             neighbor = at[center].neighbor[neighbor_index];
    // INCHI✔️✔️:             for (ii = 0; ii < at[neighbor].valence; ii++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (at[neighbor].neighbor[ii] == center)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /*  neighbor-center bond found */
    // INCHI✔️✔️:                     at[neighbor].bond_type[ii] = bond_type;
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             n++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return n;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SetTautomericBonds
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: SetTautomericBonds
    // INCHI✔️✔️: #define REPLACE_ALT_WITH_TAUT 1
    // INCHI✔️✔️: #define BOND_MARK_ALL 0xf0
    // END INCHI ACTIVE MACRO CONFIGURATION: SetTautomericBonds

    let mut k = 0_i32;
    let mut n = 0_i32;
    while k < nNumBondPos {
        let position = BondPos
            .get(usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let neighbor_index = usize::from(position.neighbor_index);
        let center = usize::from(position.nAtomNumber);
        let center_atom = at.get(center).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bond_mark = *center_atom
            .bond_type
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let source_bond_type = bond_mark & !(BOND_MARK_ALL as u8);
        if source_bond_type != BOND_TAUTOM as u8 {
            let bond_type = (bond_mark & BOND_MARK_ALL as u8) | BOND_TAUTOM as u8;
            let neighbor = usize::from(
                *center_atom
                    .neighbor
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            at.get_mut(center)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bond_type[neighbor_index] = bond_type;

            let neighbor_atom = at
                .get_mut(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut ii = 0_i32;
            while ii < i32::from(neighbor_atom.valence) {
                let reverse_index =
                    usize::try_from(ii).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if usize::from(neighbor_atom.neighbor[reverse_index]) == center {
                    neighbor_atom.bond_type[reverse_index] = bond_type;
                    break;
                }
                ii = ii.wrapping_add(1);
            }
            n = n.wrapping_add(1);
        }
        k = k.wrapping_add(1);
    }

    Ok(n)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetNeutralRepsIfNeeded(
    heap: &SourceHeap,
    pri: &mut AT_NUMB,
    prj: &mut AT_NUMB,
    at: &[inp_ATOM],
    num_atoms: i32,
    EndPoint: &[T_ENDPOINT],
    nNumEndPoints: i32,
    cgi: Option<&C_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1565 GetNeutralRepsIfNeeded
    // INCHI✔️✔️: int GetNeutralRepsIfNeeded( AT_NUMB      *pri,
    // INCHI✔️✔️:                             AT_NUMB      *prj,
    // INCHI✔️✔️:                             inp_ATOM     *at,
    // INCHI✔️✔️:                             int          num_atoms,
    // INCHI✔️✔️:                             T_ENDPOINT   *EndPoint,
    // INCHI✔️✔️:                             int          nNumEndPoints,
    // INCHI✔️✔️:                             C_GROUP_INFO *cgi )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_NUMB ri = *pri;
    // INCHI✔️✔️:     AT_NUMB rj = *prj;
    // INCHI✔️✔️:     int     i, k;
    // INCHI✔️✔️:     AT_NUMB c_point, endpoint, r;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (( c_point = at[ri].c_point ) &&
    // INCHI✔️✔️:         ( at[rj].c_point == c_point ) &&
    // INCHI✔️✔️:          ( at[ri].charge == 1 || at[rj].charge == 1 ) &&
    // INCHI✔️✔️:          cgi                                        &&
    // INCHI✔️✔️:          cgi->num_c_groups > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* at[ri] and at[rj] belong to the same charge group, at least one is charged   */
    // INCHI✔️✔️:         for (k = 0; k < cgi->num_c_groups; k++) /* MS VC++ 2008 reports unreachable code here ??? */ /* djb-rwth: addressing coverity ID #499559 -- read the previous comment; can cgi->num_c_groups only have values 0 and 1? */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (cgi->c_group[k].nGroupNumber == c_point)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* cgi->c_group[k] is found to be this charge group */
    // INCHI✔️✔️:                 if (cgi->c_group[k].num_CPoints - cgi->c_group[k].num[0] < 2)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* Only one neutral in the c-group: we will not be able to neutralize both */
    // INCHI✔️✔️:                     /* when looking for the alt path to discover the tautomerism.              */
    // INCHI✔️✔️:                     /* Therefore we need to find a neutral t-group representative at[rj]       */
    // INCHI✔️✔️:                     if ((endpoint = at[rj].endpoint)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         for (i = 0; i < nNumEndPoints; i++)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (( r = EndPoint[i].nAtomNumber ) == *prj)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 continue; /* ignore at[*prj] */
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             if (at[r].endpoint != endpoint)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 continue; /* at[r] does not belong to the same t-group as at[*prj]; ignore the atom */
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             if (!at[r].c_point)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 rj = r; /* found a neutral t-group representative */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             if (at[r].c_point != c_point && c_point == at[rj].c_point)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 /* replace only once because of (c_point == at[rj].c_point) condition  */
    // INCHI✔️✔️:                                 rj = r;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         if (rj == *prj /*&& at[ri].endpoint*/)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             /* !!! "&& at[ri].endpoint": only between 2 t-groups 2004-02-27;
    // INCHI✔️✔️:                             the change disabled due to undiscovered yet possibility of ambiguity*/
    // INCHI✔️✔️:                             /* no replacement has been found in EndPoint[]; try all atoms in the t-group */
    // INCHI✔️✔️:                             for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 if (at[i].endpoint != endpoint)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     continue;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                                 if (i == (int) *prj)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     continue;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                                 if (!at[i].c_point)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     rj = (AT_NUMB) i; /* found neutral t-group representative */
    // INCHI✔️✔️:                                     break;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                                 if (at[i].c_point != c_point && c_point == at[rj].c_point)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     /* replace only once */
    // INCHI✔️✔️:                                     rj = (AT_NUMB) i;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     /* at[ri] */
    // INCHI✔️✔️:                     if ((endpoint = at[ri].endpoint)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         for (i = 0; i < nNumEndPoints; i++)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (( r = EndPoint[i].nAtomNumber ) == *pri)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 continue;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             if (at[r].endpoint != endpoint)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 continue;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             if (!at[r].c_point)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 ri = r; /* found neutral t-group representative */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             if (at[r].c_point != c_point && c_point == at[ri].c_point &&
    // INCHI✔️✔️:                                  at[r].c_point != at[rj].c_point)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 /* replace only once */
    // INCHI✔️✔️:                                 ri = r;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         if (ri == *pri && at[rj].endpoint)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             /* !!! "&& at[rj].endpoint": only between 2 t-groups 2004-02-27;
    // INCHI✔️✔️:                             the change disabled due to undiscovered yet possibility of ambiguity */
    // INCHI✔️✔️:                             for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 if (at[i].endpoint != endpoint)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     continue;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                                 if (i == (int) *pri)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     continue;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                                 if (!at[i].c_point)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     ri = (AT_NUMB) i; /* found neutral t-group representative */
    // INCHI✔️✔️:                                     break;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                                 if (at[i].c_point != c_point && c_point == at[ri].c_point &&
    // INCHI✔️✔️:                                      at[i].c_point != at[rj].c_point)
    // INCHI✔️✔️:                                 {
    // INCHI✔️✔️:                                     /* replace only once */
    // INCHI✔️✔️:                                     ri = (AT_NUMB) i;
    // INCHI✔️✔️:                                 }
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         *prj = rj;
    // INCHI✔️✔️:         *pri = ri;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetNeutralRepsIfNeeded

    let mut ri = *pri;
    let mut rj = *prj;
    let ri_index = usize::from(ri);
    let ri_atom = at
        .get(ri_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let c_point = ri_atom.c_point;

    if c_point != 0 {
        let rj_atom = at
            .get(usize::from(rj))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if rj_atom.c_point == c_point
            && (ri_atom.charge == 1 || rj_atom.charge == 1)
            && cgi.is_some_and(|value| value.num_c_groups > 0)
        {
            let cgi = cgi.expect("checked above");
            let groups = heap.slice(cgi.c_group.as_const())?;
            let mut k = 0_i32;
            while k < cgi.num_c_groups {
                let group = groups
                    .get(usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if group.nGroupNumber == c_point
                    && i32::from(group.num_CPoints) - i32::from(group.num[0]) < 2
                {
                    let endpoint = at
                        .get(usize::from(rj))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint;
                    if endpoint != 0 {
                        let mut i = 0_i32;
                        while i < nNumEndPoints {
                            let r = EndPoint
                                .get(
                                    usize::try_from(i)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nAtomNumber;
                            if r == *prj {
                                i = i.wrapping_add(1);
                                continue;
                            }
                            let candidate = at
                                .get(usize::from(r))
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if candidate.endpoint != endpoint {
                                i = i.wrapping_add(1);
                                continue;
                            }
                            if candidate.c_point == 0 {
                                rj = r;
                                break;
                            }
                            if candidate.c_point != c_point
                                && c_point
                                    == at
                                        .get(usize::from(rj))
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                        .c_point
                            {
                                rj = r;
                            }
                            i = i.wrapping_add(1);
                        }
                        if rj == *prj {
                            let mut i = 0_i32;
                            while i < num_atoms {
                                let index = usize::try_from(i)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let candidate =
                                    at.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if candidate.endpoint != endpoint || i == i32::from(*prj) {
                                    i = i.wrapping_add(1);
                                    continue;
                                }
                                if candidate.c_point == 0 {
                                    rj = i as AT_NUMB;
                                    break;
                                }
                                if candidate.c_point != c_point
                                    && c_point
                                        == at
                                            .get(usize::from(rj))
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                                            .c_point
                                {
                                    rj = i as AT_NUMB;
                                }
                                i = i.wrapping_add(1);
                            }
                        }
                    }

                    let endpoint = at
                        .get(usize::from(ri))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint;
                    if endpoint != 0 {
                        let mut i = 0_i32;
                        while i < nNumEndPoints {
                            let r = EndPoint
                                .get(
                                    usize::try_from(i)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nAtomNumber;
                            if r == *pri {
                                i = i.wrapping_add(1);
                                continue;
                            }
                            let candidate = at
                                .get(usize::from(r))
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if candidate.endpoint != endpoint {
                                i = i.wrapping_add(1);
                                continue;
                            }
                            if candidate.c_point == 0 {
                                ri = r;
                                break;
                            }
                            if candidate.c_point != c_point
                                && c_point
                                    == at
                                        .get(usize::from(ri))
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                        .c_point
                                && candidate.c_point
                                    != at
                                        .get(usize::from(rj))
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                        .c_point
                            {
                                ri = r;
                            }
                            i = i.wrapping_add(1);
                        }
                        if ri == *pri
                            && at
                                .get(usize::from(rj))
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .endpoint
                                != 0
                        {
                            let mut i = 0_i32;
                            while i < num_atoms {
                                let index = usize::try_from(i)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let candidate =
                                    at.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if candidate.endpoint != endpoint || i == i32::from(*pri) {
                                    i = i.wrapping_add(1);
                                    continue;
                                }
                                if candidate.c_point == 0 {
                                    ri = i as AT_NUMB;
                                    break;
                                }
                                if candidate.c_point != c_point
                                    && c_point
                                        == at
                                            .get(usize::from(ri))
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                                            .c_point
                                    && candidate.c_point
                                        != at
                                            .get(usize::from(rj))
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                                            .c_point
                                {
                                    ri = i as AT_NUMB;
                                }
                                i = i.wrapping_add(1);
                            }
                        }
                    }
                }
                k = k.wrapping_add(1);
                break;
            }
            *prj = rj;
            *pri = ri;
        }
    }

    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FindAccessibleEndPoints(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    EndPoint: &mut [T_ENDPOINT],
    nNumEndPoints: &mut i32,
    BondPos: &mut [T_BONDPOS],
    nNumBondPos: &mut i32,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    cgi: Option<&C_GROUP_INFO>,
    taut_mode: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1712 FindAccessibleEndPoints
    // INCHI✔️❌: int FindAccessibleEndPoints( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                              T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                              int *nNumEndPoints,
    // INCHI✔️❌:                              T_BONDPOS *BondPos,
    // INCHI✔️❌:                              int *nNumBondPos,
    // INCHI✔️❌:                              struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                              struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                              inp_ATOM *at,
    // INCHI✔️❌:                              int num_atoms,
    // INCHI✔️❌:                              C_GROUP_INFO *cgi,
    // INCHI✔️❌:                              int taut_mode )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB nTGroupRepresenative[MAXVAL], nTGroupEqu[MAXVAL], nTGEndPointNo[MAXVAL], ri, rj;
    // INCHI✔️❌:     AT_NUMB nCurTGroupNumber, nMaxTGroupNumber, nNumTgroupNumbers, nMaxEquNumber;
    // INCHI✔️❌:     int   i, j, k, nNumDiffTGroupNumbers = 0, nNumFoundEqu, nErr;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (*nNumEndPoints != *nNumBondPos)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     /* collect all group numbers. Fill EndPoint[i].nEquNumber */
    // INCHI✔️❌:     for (i = 0; i < *nNumEndPoints; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nCurTGroupNumber = EndPoint[i].nEquNumber = EndPoint[i].nGroupNumber; /* initial equivalence */
    // INCHI✔️❌:         if (nCurTGroupNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* found endpoint that already belongs to a t-group */
    // INCHI✔️❌:             for (j = 0; j < nNumDiffTGroupNumbers; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nTGroupEqu[j] == nCurTGroupNumber)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (j == nNumDiffTGroupNumbers)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nTGroupRepresenative[nNumDiffTGroupNumbers] = EndPoint[i].nAtomNumber;
    // INCHI✔️❌:                 nTGroupEqu[nNumDiffTGroupNumbers] = EndPoint[i].nGroupNumber;
    // INCHI✔️❌:                 nTGEndPointNo[nNumDiffTGroupNumbers] = i;
    // INCHI✔️❌:                 nNumDiffTGroupNumbers++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check whether each pair belongs to the same t-group and establish the equivalence(s) */
    // INCHI✔️❌:     for (i = 0, nNumFoundEqu = 0; i < nNumDiffTGroupNumbers; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = i + 1; j < nNumDiffTGroupNumbers; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ri = nTGroupRepresenative[i];
    // INCHI✔️❌:             rj = nTGroupRepresenative[j];
    // INCHI✔️❌:             /* both at[ri] and at[rj] are known to belong to tautomeric groups */
    // INCHI✔️❌:             GetNeutralRepsIfNeeded( &ri, &rj, at, num_atoms, EndPoint, *nNumEndPoints, cgi );
    // INCHI✔️❌:             nErr = bExistsAnyAltPath( pCG, pBNS, pBD, at, num_atoms, ri, rj, taut_mode );
    // INCHI✔️❌:             if (IS_BNS_ERROR( nErr ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return nErr;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (0 == nErr)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue; /* alt path between at[ri] and at[rj] not found */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nCurTGroupNumber = inchi_min( nTGroupEqu[i], nTGroupEqu[j] );
    // INCHI✔️❌:             nMaxTGroupNumber = inchi_max( nTGroupEqu[i], nTGroupEqu[j] );
    // INCHI✔️❌:             for (k = 0; k < nNumDiffTGroupNumbers; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nTGroupEqu[k] == nMaxTGroupNumber)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nTGroupEqu[k] = nCurTGroupNumber;
    // INCHI✔️❌:                     nNumFoundEqu++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (k = 0; k < *nNumEndPoints; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (EndPoint[k].nEquNumber == nMaxTGroupNumber)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     EndPoint[k].nEquNumber = nCurTGroupNumber;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumFoundEqu)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* leave in only non-equivalent representatives */
    // INCHI✔️❌:         for (i = 1; i < nNumDiffTGroupNumbers; i++) /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = 0; j < i; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nTGroupEqu[j] == nTGroupEqu[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nTGroupEqu[i] = 0;  /* i > j; mark equivalent for removal*/
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (i = j = 0; i < nNumDiffTGroupNumbers; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nTGroupEqu[i])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i != j)
    // INCHI✔️❌:                 { /* remove the marked */
    // INCHI✔️❌:                     nTGroupEqu[j] = nTGroupEqu[i];
    // INCHI✔️❌:                     nTGroupRepresenative[j] = nTGroupRepresenative[i];
    // INCHI✔️❌:                     nTGEndPointNo[j] = nTGEndPointNo[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 j++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nNumDiffTGroupNumbers = j; /* number of known t-group representatives */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* collect endpoints that have not been assigned to t-groups */
    // INCHI✔️❌:     for (i = 0, j = nNumDiffTGroupNumbers; i < *nNumEndPoints; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (EndPoint[i].nEquNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nTGroupEqu[j] = 0;
    // INCHI✔️❌:         nTGroupRepresenative[j] = EndPoint[i].nAtomNumber;
    // INCHI✔️❌:         nTGEndPointNo[j] = i;
    // INCHI✔️❌:         j++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     nNumTgroupNumbers = j;
    // INCHI✔️❌:     nMaxEquNumber = num_atoms + 1; /* impossible atom or t-group number */
    // INCHI✔️❌:
    // INCHI✔️❌:                                    /* check whether each pair belongs to the same group and establish the equivalence(s) */
    // INCHI✔️❌:     for (i = 0, nNumFoundEqu = 0; i < nNumTgroupNumbers; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = i + 1; j < nNumTgroupNumbers; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((nTGroupEqu[i] != nTGroupEqu[j] && ( i >= nNumDiffTGroupNumbers || j >= nNumDiffTGroupNumbers )) ||
    // INCHI✔️❌:                  /* equivalence of a t-group and a non-t-group atom */
    // INCHI✔️❌:                  (!nTGroupEqu[i] && !nTGroupEqu[j])
    // INCHI✔️❌:                  /* equivalence of two non-t-group atoms */
    // INCHI✔️❌:                  ) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ri = nTGroupRepresenative[i];
    // INCHI✔️❌:                 rj = nTGroupRepresenative[j];
    // INCHI✔️❌:                 /*------------------------------!!!---------------------------------------------
    // INCHI✔️❌:                 Explanation why GetNeutralRepsIfNeeded() may need to be changed 2004-02-27
    // INCHI✔️❌:                 The change has been disabled due to undiscovered yet possibility of ambiguity
    // INCHI✔️❌:                 to search for neutral only among EndPoint[] in case taut-not_taut pairs
    // INCHI✔️❌:
    // INCHI✔️❌:                 Counterexample:   O=C-NH(+)=C-NH2
    // INCHI✔️❌:                 1   2       3
    // INCHI✔️❌:
    // INCHI✔️❌:                 Has already been found: 2-3 (+)-charge exchange
    // INCHI✔️❌:                 1-2 tautomerism (charge removed to 3)
    // INCHI✔️❌:                 Now testing: 2-3 tautomerism. If not commented out,
    // INCHI✔️❌:                 GetNeutralRepsIfNeeded() would replace 2-3 test with 1-3 test because:
    // INCHI✔️❌:                 o Charge group has only one neutral and both 2 and 3 belong to it,
    // INCHI✔️❌:                 therefore we cannot neutralize both; search for neutral representative;
    // INCHI✔️❌:                 o Since 1 and 2 belong to the same t-group and 1 is neutral,
    // INCHI✔️❌:                 test 1-3 instead of 2-3.
    // INCHI✔️❌:                 This breaks our condition:
    // INCHI✔️❌:                 Test tautomeric H movement only between neutral atoms.
    // INCHI✔️❌:                 -----------------------------------------------------------------------------*/
    // INCHI✔️❌:                 GetNeutralRepsIfNeeded( &ri, &rj, at, num_atoms, EndPoint, *nNumEndPoints, cgi );
    // INCHI✔️❌:
    // INCHI✔️❌:                 nErr = bExistsAnyAltPath( pCG, pBNS, pBD, at, num_atoms, ri, rj, taut_mode );
    // INCHI✔️❌:                 if (IS_BNS_ERROR( nErr ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return nErr;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (nErr <= 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (nTGroupEqu[i] && nTGroupEqu[j])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* found equivalence of two t-groups; at least one of them must be a new one */
    // INCHI✔️❌:                     nCurTGroupNumber = inchi_min( nTGroupEqu[i], nTGroupEqu[j] );
    // INCHI✔️❌:                     nMaxTGroupNumber = inchi_max( nTGroupEqu[i], nTGroupEqu[j] );
    // INCHI✔️❌:                     for (k = 0; k < nNumTgroupNumbers; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nTGroupEqu[k] == nMaxTGroupNumber)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nTGroupEqu[k] = nCurTGroupNumber;
    // INCHI✔️❌:                             nNumFoundEqu++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     for (k = 0; k < *nNumEndPoints; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (EndPoint[k].nEquNumber == nMaxTGroupNumber)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             EndPoint[k].nEquNumber = nCurTGroupNumber;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (nTGroupEqu[i])
    // INCHI✔️❌:                     { /* extend existing t-group */
    // INCHI✔️❌:                         nTGroupEqu[j] = nTGroupEqu[i];
    // INCHI✔️❌:                         EndPoint[nTGEndPointNo[j]].nEquNumber = nTGroupEqu[i];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nTGroupEqu[j])
    // INCHI✔️❌:                         { /* extend existing t-group */
    // INCHI✔️❌:                             nTGroupEqu[i] = nTGroupEqu[j];
    // INCHI✔️❌:                             EndPoint[nTGEndPointNo[i]].nEquNumber = nTGroupEqu[j];
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         { /* establis a new t-group */
    // INCHI✔️❌:                             nTGroupEqu[i] =
    // INCHI✔️❌:                                 nTGroupEqu[j] = nMaxEquNumber; /* assign a fict. ID to establish equivalence */
    // INCHI✔️❌:                             EndPoint[nTGEndPointNo[i]].nEquNumber =
    // INCHI✔️❌:                                 EndPoint[nTGEndPointNo[j]].nEquNumber = nMaxEquNumber;
    // INCHI✔️❌:                             nMaxEquNumber++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* eliminate endpoints and bonds that do not belong to t-group(s)
    // INCHI✔️❌:     (they have not been found connected by an alt path to any other endpoint)
    // INCHI✔️❌:     */
    // INCHI✔️❌:     for (i = 0, j = 0; i < *nNumEndPoints; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (EndPoint[i].nEquNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌: #if ( IGNORE_SINGLE_ENDPOINTS == 1 )  /* 1-28-2003 */
    // INCHI✔️❌:             for (k = 0, nNumFoundEqu = 0; k < *nNumEndPoints; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nNumFoundEqu += ( EndPoint[i].nEquNumber == EndPoint[k].nEquNumber );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nNumFoundEqu <= 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* one time it is equal to itself when i == k above */
    // INCHI✔️❌:                 /* if EndPoint[i] is not "equivalent" to any other EndPoint then ignore it */
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             if (i != j)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* save endpoints that are found to be connected to other endpoints by alt paths */
    // INCHI✔️❌:                 EndPoint[j] = EndPoint[i];
    // INCHI✔️❌:                 BondPos[j] = BondPos[i];
    // INCHI✔️❌:             }
    // INCHI✔️❌:             j++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( IGNORE_SINGLE_ENDPOINTS != 1 )  /* 1-28-2003 */
    // INCHI❌❌:     /* Do not allow a centerpoint to have only one tautomeric bond */
    // INCHI❌❌:     /* Hack: we may have only one centerpoint */
    // INCHI❌❌:     /* BondPos[*].nAtomNumber are centerpoints */
    // INCHI❌❌:     if (j == 1)
    // INCHI❌❌:     {
    // INCHI❌❌:         /* check if there exist other centerpoint neighbors
    // INCHI❌❌:         * connected to it by another tautomeric-bond
    // INCHI❌❌:         */
    // INCHI❌❌:         for (i = 0, k = 0; i < at[BondPos[0].nAtomNumber].valence; i++)
    // INCHI❌❌:         {
    // INCHI❌❌:             k += ( i != BondPos[0].neighbor_index &&
    // INCHI❌❌:                    BOND_TAUTOM == ( at[BondPos[0].nAtomNumber].bond_type[i] & ~BOND_MARK_ALL ) );
    // INCHI❌❌:         }
    // INCHI❌❌:         if (!k)
    // INCHI❌❌:         {
    // INCHI❌❌:             j = 0;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     *nNumEndPoints = *nNumBondPos = j;
    // INCHI✔️❌:
    // INCHI✔️❌:     return j;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FindAccessibleEndPoints
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FindAccessibleEndPoints
    // INCHI✔️❌: #define MAXVAL 20
    // INCHI✔️❌: #define IGNORE_SINGLE_ENDPOINTS 1
    // END INCHI ACTIVE MACRO CONFIGURATION: FindAccessibleEndPoints

    if *nNumEndPoints != *nNumBondPos {
        return Ok(0);
    }

    let mut representatives = [0 as AT_NUMB; MAXVAL as usize];
    let mut equivalences = [0 as AT_NUMB; MAXVAL as usize];
    let mut endpoint_numbers = [0 as AT_NUMB; MAXVAL as usize];
    let mut different_group_count = 0_i32;
    let mut i = 0_i32;
    while i < *nNumEndPoints {
        let endpoint_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let endpoint = EndPoint
            .get_mut(endpoint_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        endpoint.nEquNumber = endpoint.nGroupNumber;
        let current_group = endpoint.nEquNumber;
        if current_group != 0 {
            let mut j = 0_i32;
            while j < different_group_count {
                let index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if equivalences[index] == current_group {
                    break;
                }
                j = j.wrapping_add(1);
            }
            if j == different_group_count {
                let index = usize::try_from(different_group_count)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                *representatives
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = endpoint.nAtomNumber;
                *equivalences
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = endpoint.nGroupNumber;
                *endpoint_numbers
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = i as AT_NUMB;
                different_group_count = different_group_count.wrapping_add(1);
            }
        }
        i = i.wrapping_add(1);
    }

    let is_bns_error = |value: i32| BNS_ERR <= value && value <= BNS_MAX_ERR_VALUE;
    let mut found_equivalences = 0_i32;
    i = 0;
    while i < different_group_count {
        let mut j = i.wrapping_add(1);
        while j < different_group_count {
            let i_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let j_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut ri = representatives[i_index];
            let mut rj = representatives[j_index];
            {
                let atoms = heap.slice(at.as_const())?;
                GetNeutralRepsIfNeeded(
                    heap,
                    &mut ri,
                    &mut rj,
                    atoms,
                    num_atoms,
                    EndPoint,
                    *nNumEndPoints,
                    cgi,
                )?;
            }
            let error = bExistsAnyAltPath(
                heap,
                pCG,
                pBNS,
                pBD,
                at,
                num_atoms,
                i32::from(ri),
                i32::from(rj),
                taut_mode,
                clock_result,
            )?;
            if is_bns_error(error) {
                return Ok(error);
            }
            if error != 0 {
                let current_group = equivalences[i_index].min(equivalences[j_index]);
                let maximum_group = equivalences[i_index].max(equivalences[j_index]);
                let mut k = 0_i32;
                while k < different_group_count {
                    let index =
                        usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    if equivalences[index] == maximum_group {
                        equivalences[index] = current_group;
                        found_equivalences = found_equivalences.wrapping_add(1);
                    }
                    k = k.wrapping_add(1);
                }
                k = 0;
                while k < *nNumEndPoints {
                    let index =
                        usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let endpoint = EndPoint
                        .get_mut(index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if endpoint.nEquNumber == maximum_group {
                        endpoint.nEquNumber = current_group;
                    }
                    k = k.wrapping_add(1);
                }
            }
            j = j.wrapping_add(1);
        }
        i = i.wrapping_add(1);
    }

    if found_equivalences != 0 {
        i = 1;
        while i < different_group_count {
            let i_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut j = 0_i32;
            while j < i {
                let j_index =
                    usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if equivalences[j_index] == equivalences[i_index] {
                    equivalences[i_index] = 0;
                    break;
                }
                j = j.wrapping_add(1);
            }
            i = i.wrapping_add(1);
        }
        i = 0;
        let mut j = 0_i32;
        while i < different_group_count {
            let i_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if equivalences[i_index] != 0 {
                let j_index =
                    usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if i != j {
                    equivalences[j_index] = equivalences[i_index];
                    representatives[j_index] = representatives[i_index];
                    endpoint_numbers[j_index] = endpoint_numbers[i_index];
                }
                j = j.wrapping_add(1);
            }
            i = i.wrapping_add(1);
        }
        different_group_count = j;
    }

    i = 0;
    let mut j = different_group_count;
    while i < *nNumEndPoints {
        let endpoint_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let endpoint = EndPoint
            .get(endpoint_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if endpoint.nEquNumber == 0 {
            let index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *equivalences
                .get_mut(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            *representatives
                .get_mut(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = endpoint.nAtomNumber;
            *endpoint_numbers
                .get_mut(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = i as AT_NUMB;
            j = j.wrapping_add(1);
        }
        i = i.wrapping_add(1);
    }
    let group_number_count = j;
    let mut maximum_equivalence = num_atoms.wrapping_add(1) as AT_NUMB;

    i = 0;
    while i < group_number_count {
        j = i.wrapping_add(1);
        while j < group_number_count {
            let i_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let j_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if (equivalences[i_index] != equivalences[j_index]
                && (i >= different_group_count || j >= different_group_count))
                || (equivalences[i_index] == 0 && equivalences[j_index] == 0)
            {
                let mut ri = representatives[i_index];
                let mut rj = representatives[j_index];
                {
                    let atoms = heap.slice(at.as_const())?;
                    GetNeutralRepsIfNeeded(
                        heap,
                        &mut ri,
                        &mut rj,
                        atoms,
                        num_atoms,
                        EndPoint,
                        *nNumEndPoints,
                        cgi,
                    )?;
                }
                let error = bExistsAnyAltPath(
                    heap,
                    pCG,
                    pBNS,
                    pBD,
                    at,
                    num_atoms,
                    i32::from(ri),
                    i32::from(rj),
                    taut_mode,
                    clock_result,
                )?;
                if is_bns_error(error) {
                    return Ok(error);
                }
                if error > 0 {
                    if equivalences[i_index] != 0 && equivalences[j_index] != 0 {
                        let current_group = equivalences[i_index].min(equivalences[j_index]);
                        let maximum_group = equivalences[i_index].max(equivalences[j_index]);
                        let mut k = 0_i32;
                        while k < group_number_count {
                            let index = usize::try_from(k)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            if equivalences[index] == maximum_group {
                                equivalences[index] = current_group;
                                found_equivalences = found_equivalences.wrapping_add(1);
                            }
                            k = k.wrapping_add(1);
                        }
                        k = 0;
                        while k < *nNumEndPoints {
                            let index = usize::try_from(k)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let endpoint = EndPoint
                                .get_mut(index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if endpoint.nEquNumber == maximum_group {
                                endpoint.nEquNumber = current_group;
                            }
                            k = k.wrapping_add(1);
                        }
                    } else if equivalences[i_index] != 0 {
                        equivalences[j_index] = equivalences[i_index];
                        let endpoint_index = usize::from(endpoint_numbers[j_index]);
                        EndPoint
                            .get_mut(endpoint_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nEquNumber = equivalences[i_index];
                    } else if equivalences[j_index] != 0 {
                        equivalences[i_index] = equivalences[j_index];
                        let endpoint_index = usize::from(endpoint_numbers[i_index]);
                        EndPoint
                            .get_mut(endpoint_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nEquNumber = equivalences[j_index];
                    } else {
                        equivalences[i_index] = maximum_equivalence;
                        equivalences[j_index] = maximum_equivalence;
                        let i_endpoint = usize::from(endpoint_numbers[i_index]);
                        let j_endpoint = usize::from(endpoint_numbers[j_index]);
                        EndPoint
                            .get_mut(i_endpoint)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nEquNumber = maximum_equivalence;
                        EndPoint
                            .get_mut(j_endpoint)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nEquNumber = maximum_equivalence;
                        maximum_equivalence = maximum_equivalence.wrapping_add(1);
                    }
                }
            }
            j = j.wrapping_add(1);
        }
        i = i.wrapping_add(1);
    }

    i = 0;
    j = 0;
    while i < *nNumEndPoints {
        let i_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let equivalence = EndPoint
            .get(i_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nEquNumber;
        if equivalence != 0 {
            let mut k = 0_i32;
            found_equivalences = 0;
            while k < *nNumEndPoints {
                let k_index =
                    usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                found_equivalences = found_equivalences.wrapping_add(i32::from(
                    equivalence
                        == EndPoint
                            .get(k_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nEquNumber,
                ));
                k = k.wrapping_add(1);
            }
            if found_equivalences > 1 {
                let j_index =
                    usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if i != j {
                    let endpoint = EndPoint
                        .get(i_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let bond = BondPos
                        .get(i_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    *EndPoint
                        .get_mut(j_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = endpoint;
                    *BondPos
                        .get_mut(j_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = bond;
                }
                j = j.wrapping_add(1);
            }
        }
        i = i.wrapping_add(1);
    }

    *nNumEndPoints = j;
    *nNumBondPos = j;
    Ok(j)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MarkTautomerGroups(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    t_group_info: Option<&mut T_GROUP_INFO>,
    mut c_group_info: Option<&mut C_GROUP_INFO>,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:4336 MarkTautomerGroups
    // INCHI✔️❌: int MarkTautomerGroups( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                         inp_ATOM *at,
    // INCHI✔️❌:                         int num_atoms,
    // INCHI✔️❌:                         T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                         C_GROUP_INFO *c_group_info,
    // INCHI✔️❌:                         struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                         struct BalancedNetworkData *pBD )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, m, endpoint_valence, centerpoint, endpoint, bond_type, nMobile, num_changes = 0, tot_changes = 0;
    // INCHI✔️❌:     T_ENDPOINT EndPoint[MAXVAL];
    // INCHI✔️❌:     T_BONDPOS  BondPos[MAXVAL];
    // INCHI✔️❌:     AT_NUMB    nGroupNumber;
    // INCHI✔️❌:     int        bDiffGroups;
    // INCHI✔️❌:     int  nNumEndPoints, nNumBondPos, nNumPossibleMobile;
    // INCHI✔️❌:     int  bNonTautBond, bAltBond; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int  nNumDonor, nNumAcceptor, bPossiblyEndpoint;
    // INCHI✔️❌:     int bIgnoreIsotopic; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     ENDPOINT_INFO eif1, eif2;
    // INCHI✔️❌:     int nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #define ALLOWED_EDGE(PBNS, IAT,IBOND)  ( !(PBNS) || !(PBNS)->edge || !(PBNS)->vert || !(PBNS)->edge[(PBNS)->vert[IAT].iedge[IBOND]].forbidden)
    // INCHI✔️❌: #define ACTUAL_ORDER(PBNS, IAT,IBOND, BTYPE)  ( ((PBNS) && (PBNS)->edge && (PBNS)->vert &&\
    // INCHI✔️❌:     ((BTYPE)==BOND_ALT_123 || (BTYPE)==BOND_ALT_13 || (BTYPE)==BOND_ALT_23))? (PBNS)->edge[(PBNS)->vert[IAT].iedge[IBOND]].flow+BOND_TYPE_SINGLE:(BTYPE))
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!t_group_info || !( t_group_info->bTautFlags & TG_FLAG_TEST_TAUT__ATOMS ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Initial t_group allocation */
    // INCHI✔️❌:     if (!t_group_info->t_group && !t_group_info->max_num_t_groups)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         INCHI_MODE bTautFlags = t_group_info->bTautFlags;       /*  save initial setting */
    // INCHI✔️❌:         INCHI_MODE bTautFlagsDone = t_group_info->bTautFlagsDone;   /*  save previous findings, if any */
    // INCHI✔️❌:         TNI       tni = t_group_info->tni;
    // INCHI✔️❌:         AT_NUMB   *tGroupNumber = t_group_info->tGroupNumber;
    // INCHI✔️❌:         T_GROUP* tgi_tgr = NULL;  /* copied from below 2024-09-01 DT */
    // INCHI✔️❌:
    // INCHI✔️❌:         bIgnoreIsotopic = t_group_info->bIgnoreIsotopic;
    // INCHI✔️❌:         memset( t_group_info, 0, sizeof( *t_group_info ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         t_group_info->bIgnoreIsotopic = bIgnoreIsotopic; /*  restore initial setting */
    // INCHI✔️❌:         t_group_info->bTautFlags = bTautFlags;
    // INCHI✔️❌:         t_group_info->bTautFlagsDone = bTautFlagsDone;
    // INCHI✔️❌:         t_group_info->tni = tni;
    // INCHI✔️❌:         t_group_info->tGroupNumber = tGroupNumber;
    // INCHI✔️❌:         t_group_info->max_num_t_groups = num_atoms / 2 + 1; /*  upper limit */
    // INCHI✔️❌:         /* djb-rwth: fixing oss-fuzz issue #52978 */
    // INCHI✔️❌:         tgi_tgr = (T_GROUP*)inchi_calloc((long long)t_group_info->max_num_t_groups + 1, sizeof(t_group_info->t_group[0]));
    // INCHI✔️❌:         if (!tgi_tgr)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t_group_info->max_num_t_groups = -1;
    // INCHI✔️❌:             t_group_info->t_group = NULL;
    // INCHI✔️❌:             return (-1); /*  failed, out of RAM */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t_group_info->t_group = tgi_tgr;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Check if t_group_info exists */
    // INCHI✔️❌:     if (!t_group_info->t_group || !t_group_info->max_num_t_groups)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 > t_group_info->max_num_t_groups)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return t_group_info->max_num_t_groups;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  1-3 tautomers */
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  Find possible endpoint Z = at[i] */
    // INCHI✔️❌:         if ((endpoint_valence = nGetEndpointInfo( at, i, &eif1 ))) /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:             for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bond_type = (int) at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                 bond_type = ACTUAL_ORDER( pBNS, i, j, bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:                 centerpoint = (int) at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                 if (( bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                       bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                       bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                       bond_type == BOND_TAUTOM ) && is_centerpoint_elem( at[centerpoint].el_number )
    // INCHI✔️❌:                      && ALLOWED_EDGE( pBNS, i, j )
    // INCHI✔️❌:                      )
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  Test a centerpoint candidate. */
    // INCHI✔️❌:                     /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                     nNumPossibleMobile = 0;
    // INCHI✔️❌:                     nGroupNumber = (AT_NUMB) num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                     bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                     nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                     for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                         bond_type = (int) at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                         bond_type = ACTUAL_ORDER( pBNS, centerpoint, k, bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             bNonTautBond =
    // INCHI✔️❌:                             bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                             bPossiblyEndpoint = 0;
    // INCHI✔️❌:                         if (!ALLOWED_EDGE( pBNS, centerpoint, k ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                 bAltBond = ( bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS ); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     bNonTautBond = 1;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     continue;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (!( endpoint_valence = nGetEndpointInfo( at, endpoint, &eif1 ) )) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                             continue; /*  not an endpoint element or can't have mobile groups */
    // INCHI✔️❌:                                       /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                                       /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                         nMobile =
    // INCHI✔️❌:                             AddAtom2num( EndPoint[nNumEndPoints].num, at, endpoint, 2 ); /* fill out */
    // INCHI✔️❌:                         AddAtom2DA( EndPoint[nNumEndPoints].num_DA, at, endpoint, 2 );
    // INCHI✔️❌:                         /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                         nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                         nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                         */
    // INCHI✔️❌:                         if (bNonTautBond)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             m = ( bond_type == BOND_SINGLE && ( nMobile || at[endpoint].endpoint ) );
    // INCHI✔️❌:                             nNumDonor += m;
    // INCHI✔️❌:                             bPossiblyEndpoint += m;
    // INCHI✔️❌:                             m = ( bond_type == BOND_DOUBLE );
    // INCHI✔️❌:                             nNumAcceptor += m;
    // INCHI✔️❌:                             bPossiblyEndpoint += m;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  Tautomeric or alternating bond */
    // INCHI✔️❌:                             m = ( 0 != at[endpoint].endpoint || eif1.cDonor );
    // INCHI✔️❌:                             nNumDonor += m;
    // INCHI✔️❌:                             bPossiblyEndpoint += m;
    // INCHI✔️❌:                             m = ( at[endpoint].endpoint ||
    // INCHI✔️❌:                                   eif1.cNeutralBondsValence > at[endpoint].valence );
    // INCHI✔️❌:                             nNumAcceptor += m;
    // INCHI✔️❌:                             bPossiblyEndpoint += m;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (!bPossiblyEndpoint)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                         EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                         EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB) endpoint;
    // INCHI✔️❌:                         if (nGroupNumber != at[endpoint].endpoint)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bDiffGroups++;
    // INCHI✔️❌:                             nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         /*  save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                         if (bNonTautBond || bAltBond)
    // INCHI❌❌:                         {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             BondPos[nNumBondPos].nAtomNumber = (AT_NUMB) centerpoint;
    // INCHI✔️❌:                             BondPos[nNumBondPos].neighbor_index = (AT_NUMB) k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                             nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                         }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                         /*  mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                         /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                         nNumPossibleMobile += ( nMobile > 0 || at[endpoint].endpoint );
    // INCHI✔️❌:                         nNumEndPoints++;
    // INCHI✔️❌:                         /*printf("Found %d %d %d %d\n", centerpoint+1, at[centerpoint].el_number, endpoint+1, at[endpoint].el_number);*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*
    // INCHI✔️❌:                         * a tautomeric group has been found
    // INCHI✔️❌:                         *
    // INCHI✔️❌:                         * at this point:
    // INCHI✔️❌:                         * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                         * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                         * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                         * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                         */
    // INCHI✔️❌:
    // INCHI✔️❌:                         nErr = FindAccessibleEndPoints( pCG, EndPoint,
    // INCHI✔️❌:                                                         &nNumEndPoints,
    // INCHI✔️❌:                                                         BondPos, &nNumBondPos,
    // INCHI✔️❌:                                                         pBNS, pBD, at,
    // INCHI✔️❌:                                                         num_atoms, c_group_info,
    // INCHI✔️❌:                                                         ALT_PATH_MODE_TAUTOM );
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (IS_BNS_ERROR( nErr ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return nErr;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         nErr = 0;
    // INCHI✔️❌:                         if (nNumEndPoints > 0) {
    // INCHI✔️❌:                             if (!nGroupNumber || bDiffGroups > 0) {
    // INCHI✔️❌:
    // INCHI✔️❌:                                 num_changes = RegisterEndPoints(pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
    // INCHI✔️❌:                                 if (num_changes == -1) {
    // INCHI✔️❌:                                     nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (num_changes < 0) {
    // INCHI✔️❌:                                     nErr = num_changes;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nErr)
    // INCHI✔️❌:                                     goto exit_function;
    // INCHI✔️❌:                                 tot_changes += (num_changes>0);
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (nNumBondPos > 0) {
    // INCHI✔️❌:                                 /*  some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                 num_changes = SetTautomericBonds(at, nNumBondPos, BondPos);
    // INCHI✔️❌:                                 tot_changes += (num_changes>0);
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_PT_22_00 == 1 ) /******* BEGIN PT_22_00 ********/
    // INCHI✔️❌:     if (t_group_info->bTautFlags & TG_FLAG_PT_22_00) {
    // INCHI✔️❌:         /*** [#1:1][CX4:2][NX2:3]=[CX3:4]>>[CX3:2]=[NX2:3][CX4:4][#1:1] ***/
    // INCHI✔️❌:         /*** Similar to the previous case of M=Q-ZH >> MH-Q=Z, with M,Z = "C" and Q = "N" ***/
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++) {
    // INCHI✔️❌:             /*  find possible endpoint Z = at[i] */
    // INCHI✔️❌:             if ((endpoint_valence = nGetEndpointInfo_PT_22_00(at, i, &eif1))) { /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                 /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:                 for (j = 0; j < at[i].valence; j++) {
    // INCHI✔️❌:                     bond_type = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                     bond_type = ACTUAL_ORDER(pBNS, i, j, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                     if ((bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                         bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                         bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                         bond_type == BOND_TAUTOM) && at[centerpoint].el_number == EL_NUMBER_N
    // INCHI✔️❌:                         && ALLOWED_EDGE(pBNS, i, j)
    // INCHI✔️❌:                         ) {
    // INCHI✔️❌:                         /*  test a centerpoint candidate. */
    // INCHI✔️❌:                         /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                         nNumPossibleMobile = 0;
    // INCHI✔️❌:                         nGroupNumber = (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                         bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                         nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                         for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++) {
    // INCHI✔️❌:                             endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                             bond_type = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                             bond_type = ACTUAL_ORDER(pBNS, centerpoint, k, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 bNonTautBond =
    // INCHI✔️❌:                                 bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                                 bPossiblyEndpoint = 0;
    // INCHI✔️❌:                             if (!ALLOWED_EDGE(pBNS, centerpoint, k)) {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM) {
    // INCHI✔️❌:                                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                     bAltBond = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE)
    // INCHI✔️❌:                                         bNonTautBond = 1;
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (!(endpoint_valence = nGetEndpointInfo_PT_22_00(at, endpoint, &eif1))) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                                 continue; /*  not an endpoint element or can't have mobile groups */
    // INCHI✔️❌:                                           /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                                           /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                             nMobile =
    // INCHI✔️❌:                                 AddAtom2num(EndPoint[nNumEndPoints].num, at, endpoint, 2); /* fill out */
    // INCHI✔️❌:                             AddAtom2DA(EndPoint[nNumEndPoints].num_DA, at, endpoint, 2);
    // INCHI✔️❌:                             /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             if (bNonTautBond) {
    // INCHI✔️❌:                                 m = (bond_type == BOND_SINGLE && (nMobile || at[endpoint].endpoint));
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (bond_type == BOND_DOUBLE);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else {
    // INCHI✔️❌:                                 /*  tautomeric or alternating bond */
    // INCHI✔️❌:                                 m = (0 != at[endpoint].endpoint || eif1.cDonor);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (at[endpoint].endpoint ||
    // INCHI✔️❌:                                     eif1.cNeutralBondsValence > at[endpoint].valence);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!bPossiblyEndpoint)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB)endpoint;
    // INCHI✔️❌:                             if (nGroupNumber != at[endpoint].endpoint) {
    // INCHI✔️❌:                                 bDiffGroups++;
    // INCHI✔️❌:                                 nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             if (bNonTautBond || bAltBond) {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                                 BondPos[nNumBondPos].nAtomNumber = (AT_NUMB)centerpoint;
    // INCHI✔️❌:                                 BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                                 nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /*  mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                             /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                             nNumPossibleMobile += (nMobile>0 || at[endpoint].endpoint);
    // INCHI✔️❌:                             nNumEndPoints++;
    // INCHI✔️❌:                             /*printf("Found %d %d %d %d\n", centerpoint+1, at[centerpoint].el_number, endpoint+1, at[endpoint].el_number);*/
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor) {
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             * a tautomeric group has been found
    // INCHI✔️❌:                             *
    // INCHI✔️❌:                             * at this point:
    // INCHI✔️❌:                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                             */
    // INCHI✔️❌:
    // INCHI✔️❌:                             nErr = FindAccessibleEndPoints(pCG,
    // INCHI✔️❌:                                 EndPoint, &nNumEndPoints,
    // INCHI✔️❌:                                 BondPos, &nNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, at,
    // INCHI✔️❌:                                 num_atoms, c_group_info,
    // INCHI✔️❌:                                 ALT_PATH_MODE_TAUTOM_PT_22_00);
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (IS_BNS_ERROR(nErr)) {
    // INCHI✔️❌:                                 return nErr;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nNumEndPoints > 0) {
    // INCHI✔️❌:                                 if (!nGroupNumber || bDiffGroups > 0) {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     num_changes = RegisterEndPoints(pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
    // INCHI✔️❌:                                     if (num_changes == -1) {
    // INCHI✔️❌:                                         nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (num_changes < 0) {
    // INCHI✔️❌:                                         nErr = num_changes;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                                     if (nErr)
    // INCHI✔️❌:                                         goto exit_function;
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nNumBondPos > 0) {
    // INCHI✔️❌:                                     /*  some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                     num_changes = SetTautomericBonds(at, nNumBondPos, BondPos);
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /********** END PT_22_00 ************/
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_PT_16_00 == 1 )    /******* BEGIN PT_16_00 ********/
    // INCHI✔️❌:     if (t_group_info->bTautFlags & TG_FLAG_PT_16_00) {
    // INCHI✔️❌:         /*** [#1:1][O;!R:2][N+0z1:3]=[CX3:4]>>[O;!R:2]=[N+0z1:3][CX4:4][#1:1] ***/
    // INCHI✔️❌:         /*** Similar to the previous case of M=Q-ZH >> MH-Q=Z, with M,Z = "C, O" and Q = "N" ***/
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++) {
    // INCHI✔️❌:             /*  find possible endpoint Z = at[i] */
    // INCHI✔️❌:             if ((endpoint_valence = nGetEndpointInfo_PT_16_00(at, i, &eif1))) { /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                 /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:                 for (j = 0; j < at[i].valence; j++) {
    // INCHI✔️❌:                     bond_type = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                     bond_type = ACTUAL_ORDER(pBNS, i, j, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                     if ((bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                         bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                         bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                         bond_type == BOND_TAUTOM)
    // INCHI✔️❌:                         && at[centerpoint].el_number == EL_NUMBER_N
    // INCHI✔️❌:                         && at[centerpoint].valence == 2
    // INCHI✔️❌:                         && at[centerpoint].charge == 0
    // INCHI✔️❌:                         && ALLOWED_EDGE(pBNS, i, j)
    // INCHI✔️❌:                         ) {
    // INCHI✔️❌:                         /*  test a centerpoint candidate. */
    // INCHI✔️❌:                         /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                         nNumPossibleMobile = 0;
    // INCHI✔️❌:                         nGroupNumber = (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                         bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                         nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                         for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++) {
    // INCHI✔️❌:                             endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                             bond_type = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                             bond_type = ACTUAL_ORDER(pBNS, centerpoint, k, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 bNonTautBond =
    // INCHI✔️❌:                                 bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                                 bPossiblyEndpoint = 0;
    // INCHI✔️❌:                             if (!ALLOWED_EDGE(pBNS, centerpoint, k)) {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM) {
    // INCHI✔️❌:                                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                     bAltBond = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE)
    // INCHI✔️❌:                                         bNonTautBond = 1;
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (!(endpoint_valence = nGetEndpointInfo_PT_16_00(at, endpoint, &eif1))) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                                 continue; /*  not an endpoint element or can't have mobile groups */
    // INCHI✔️❌:                             if (at[endpoint].el_number == EL_NUMBER_O &&
    // INCHI✔️❌:                                 at[endpoint].nNumAtInRingSystem != 1)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             if (endpoint != i && at[endpoint].el_number == at[i].el_number)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                             /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                             nMobile =
    // INCHI✔️❌:                                 AddAtom2num(EndPoint[nNumEndPoints].num, at, endpoint, 2); /* fill out */
    // INCHI✔️❌:                             AddAtom2DA(EndPoint[nNumEndPoints].num_DA, at, endpoint, 2);
    // INCHI✔️❌:                             /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             if (bNonTautBond) {
    // INCHI✔️❌:                                 m = (bond_type == BOND_SINGLE && (nMobile || at[endpoint].endpoint));
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (bond_type == BOND_DOUBLE);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else {
    // INCHI✔️❌:                                 /*  tautomeric or alternating bond */
    // INCHI✔️❌:                                 m = (0 != at[endpoint].endpoint || eif1.cDonor);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (at[endpoint].endpoint ||
    // INCHI✔️❌:                                     eif1.cNeutralBondsValence > at[endpoint].valence);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!bPossiblyEndpoint)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB)endpoint;
    // INCHI✔️❌:                             if (nGroupNumber != at[endpoint].endpoint) {
    // INCHI✔️❌:                                 bDiffGroups++;
    // INCHI✔️❌:                                 nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             if (bNonTautBond || bAltBond) {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                                 BondPos[nNumBondPos].nAtomNumber = (AT_NUMB)centerpoint;
    // INCHI✔️❌:                                 BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                                 nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /*  mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                             /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                             nNumPossibleMobile += (nMobile>0 || at[endpoint].endpoint);
    // INCHI✔️❌:                             nNumEndPoints++;
    // INCHI✔️❌:                             /*printf("Found %d %d %d %d\n", centerpoint+1, at[centerpoint].el_number, endpoint+1, at[endpoint].el_number);*/
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor) {
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             * a tautomeric group has been found
    // INCHI✔️❌:                             *
    // INCHI✔️❌:                             * at this point:
    // INCHI✔️❌:                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             nErr = FindAccessibleEndPoints(pCG,
    // INCHI✔️❌:                                 EndPoint, &nNumEndPoints,
    // INCHI✔️❌:                                 BondPos, &nNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, at,
    // INCHI✔️❌:                                 num_atoms, c_group_info,
    // INCHI✔️❌:                                 ALT_PATH_MODE_TAUTOM_PT_16_00);
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (IS_BNS_ERROR(nErr)) {
    // INCHI✔️❌:                                 return nErr;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nNumEndPoints > 0) {
    // INCHI✔️❌:                                 if (!nGroupNumber || bDiffGroups > 0) {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     num_changes = RegisterEndPoints(pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
    // INCHI✔️❌:                                     if (num_changes == -1) {
    // INCHI✔️❌:                                         nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (num_changes < 0) {
    // INCHI✔️❌:                                         nErr = num_changes;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                                     if (nErr)
    // INCHI✔️❌:                                         goto exit_function;
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nNumBondPos > 0) {
    // INCHI✔️❌:                                     /*  some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                     num_changes = SetTautomericBonds(at, nNumBondPos, BondPos);
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /********** END PT_16_00 ************/
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_PT_06_00 == 1 ) /******* BEGIN PT_06_00 ********/
    // INCHI✔️❌:     if (t_group_info->bTautFlags & TG_FLAG_PT_06_00) {
    // INCHI✔️❌:         /*** [CX{2-3}z{0-1},N,n,S,s,O,o,Se,Te:1]=[NX2,nX2,CX3,c,P,p:2][N,n,S,O,Se,Te:3][#1:4] >> [#1:4][CX4z{0-1},N,n,S,O,Se,Te:1][NX2,nX2,CX3z{0-1},c,P,p:2]=[N,n,S,s,O,o,Se,Te:3] ***/
    // INCHI✔️❌:         /*** Similar to the previous case of M=Q-ZH >> MH-Q=Z, with M = "C,N,S,O,Se,Te", Q = "N,C,P", and Z = "N,S,O,Se,Te" ***/
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++) {
    // INCHI✔️❌:             /*  find possible endpoint Z = at[i] */
    // INCHI✔️❌:             if ((endpoint_valence = nGetEndpointInfo_PT_06_00(at, i, &eif1))) { /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                 /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:                 for (j = 0; j < at[i].valence; j++) {
    // INCHI✔️❌:                     bond_type = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                     bond_type = ACTUAL_ORDER(pBNS, i, j, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                     if ((bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                         bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                         bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                         bond_type == BOND_TAUTOM) &&
    // INCHI✔️❌:                         (at[centerpoint].el_number == EL_NUMBER_N ||
    // INCHI✔️❌:                             at[centerpoint].el_number == EL_NUMBER_C ||
    // INCHI✔️❌:                             at[centerpoint].el_number == EL_NUMBER_P)
    // INCHI✔️❌:                         && ALLOWED_EDGE(pBNS, i, j)
    // INCHI✔️❌:                         ) {
    // INCHI✔️❌:                         /*  test a centerpoint candidate. */
    // INCHI✔️❌:                         /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                         nNumPossibleMobile = 0;
    // INCHI✔️❌:                         nGroupNumber = (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                         bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                         nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                         for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++) {
    // INCHI✔️❌:                             endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                             bond_type = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                             bond_type = ACTUAL_ORDER(pBNS, centerpoint, k, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 bNonTautBond =
    // INCHI✔️❌:                                 bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                                 bPossiblyEndpoint = 0;
    // INCHI✔️❌:                             if (!ALLOWED_EDGE(pBNS, centerpoint, k)) {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM) {
    // INCHI✔️❌:                                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                     bAltBond = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE)
    // INCHI✔️❌:                                         bNonTautBond = 1;
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:                             if (!(endpoint_valence = nGetEndpointInfo_PT_06_00(at, endpoint, &eif1))) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                                 continue; /*  not an endpoint element or can't have mobile groups */
    // INCHI✔️❌:                             if (i != endpoint &&
    // INCHI✔️❌:                                 at[endpoint].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:                                 at[i].el_number == EL_NUMBER_C)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                             /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                             nMobile =
    // INCHI✔️❌:                                 AddAtom2num(EndPoint[nNumEndPoints].num, at, endpoint, 2); /* fill out */
    // INCHI✔️❌:                             AddAtom2DA(EndPoint[nNumEndPoints].num_DA, at, endpoint, 2);
    // INCHI✔️❌:                             /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             if (bNonTautBond) {
    // INCHI✔️❌:                                 m = (bond_type == BOND_SINGLE && (nMobile || at[endpoint].endpoint));
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (bond_type == BOND_DOUBLE);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else {
    // INCHI✔️❌:                                 /*  tautomeric or alternating bond */
    // INCHI✔️❌:                                 m = (0 != at[endpoint].endpoint || eif1.cDonor);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (at[endpoint].endpoint ||
    // INCHI✔️❌:                                     eif1.cNeutralBondsValence > at[endpoint].valence);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!bPossiblyEndpoint)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB)endpoint;
    // INCHI✔️❌:                             if (nGroupNumber != at[endpoint].endpoint) {
    // INCHI✔️❌:                                 bDiffGroups++;
    // INCHI✔️❌:                                 nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             if (bNonTautBond || bAltBond) {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                                 BondPos[nNumBondPos].nAtomNumber = (AT_NUMB)centerpoint;
    // INCHI✔️❌:                                 BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                                 nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /*  mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                             /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                             nNumPossibleMobile += (nMobile>0 || at[endpoint].endpoint);
    // INCHI✔️❌:                             nNumEndPoints++;
    // INCHI✔️❌:                             /*printf("Found %d %d %d %d\n", centerpoint+1, at[centerpoint].el_number, endpoint+1, at[endpoint].el_number);*/
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor) {
    // INCHI✔️❌:                             /*printf("Real %d\n", nNumEndPoints);*/
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             * a tautomeric group has been found
    // INCHI✔️❌:                             *
    // INCHI✔️❌:                             * at this point:
    // INCHI✔️❌:                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                             */
    // INCHI✔️❌:
    // INCHI✔️❌:                             nErr = FindAccessibleEndPoints(pCG,
    // INCHI✔️❌:                                 EndPoint, &nNumEndPoints,
    // INCHI✔️❌:                                 BondPos, &nNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, at,
    // INCHI✔️❌:                                 num_atoms, c_group_info,
    // INCHI✔️❌:                                 ALT_PATH_MODE_TAUTOM_PT_06_00);
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (IS_BNS_ERROR(nErr)) {
    // INCHI✔️❌:                                 return nErr;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nNumEndPoints > 0) {
    // INCHI✔️❌:                                 if (!nGroupNumber || bDiffGroups > 0) {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     num_changes = RegisterEndPoints(pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
    // INCHI✔️❌:                                     if (num_changes == -1) {
    // INCHI✔️❌:                                         nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (num_changes < 0) {
    // INCHI✔️❌:                                         nErr = num_changes;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                                     if (nErr)
    // INCHI✔️❌:                                         goto exit_function;
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nNumBondPos > 0) {
    // INCHI✔️❌:                                     /*  some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                     num_changes = SetTautomericBonds(at, nNumBondPos, BondPos);
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /********** END PT_06_00 ************/
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_PT_39_00 == 1 )    /******* BEGIN PT_39_00 ********/
    // INCHI✔️❌:     if (t_group_info->bTautFlags & TG_FLAG_PT_39_00) {
    // INCHI✔️❌:         /*** [CX3,NX2:1]=[NX3+:2]([O-:3])[CX4:4][#1:5]>>[#1:5][CX4,NX3:1][NX3+:2]([O-:3])=[CX3:4] ***/
    // INCHI✔️❌:         /*** Similar to the previous case of M=Q-ZH >> MH-Q=Z, with M,Z = "C, N" and Q = "N+" ***/
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++) {
    // INCHI✔️❌:             /*  find possible endpoint Z = at[i] */
    // INCHI✔️❌:             if ((endpoint_valence = nGetEndpointInfo_PT_39_00(at, i, &eif1))) { /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                 /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:                 for (j = 0; j < at[i].valence; j++) {
    // INCHI✔️❌:                     bond_type = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                     bond_type = ACTUAL_ORDER(pBNS, i, j, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                     if ((bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                         bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                         bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                         bond_type == BOND_TAUTOM)
    // INCHI✔️❌:                         && at[centerpoint].el_number == EL_NUMBER_N
    // INCHI✔️❌:                         && at[centerpoint].valence == 3
    // INCHI✔️❌:                         && ALLOWED_EDGE(pBNS, i, j)
    // INCHI✔️❌:                         ) {
    // INCHI✔️❌:                         int num_O = 0;
    // INCHI✔️❌:                         int num_N = 0;
    // INCHI✔️❌:                         /*  test a centerpoint candidate. */
    // INCHI✔️❌:                         /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                         nNumPossibleMobile = 0;
    // INCHI✔️❌:                         nGroupNumber = (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                         bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                         nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                         for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++) {
    // INCHI✔️❌:                             endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                             num_O += (at[endpoint].el_number == EL_NUMBER_O);
    // INCHI✔️❌:                             num_N += (at[endpoint].el_number == EL_NUMBER_N);
    // INCHI✔️❌:                             bond_type = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                             bond_type = ACTUAL_ORDER(pBNS, centerpoint, k, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 bNonTautBond =
    // INCHI✔️❌:                                 bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                                 bPossiblyEndpoint = 0;
    // INCHI✔️❌:                             if (!ALLOWED_EDGE(pBNS, centerpoint, k)) {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM) {
    // INCHI✔️❌:                                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                     bAltBond = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE)
    // INCHI✔️❌:                                         bNonTautBond = 1;
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (!(endpoint_valence = nGetEndpointInfo_PT_39_00(at, endpoint, &eif1))) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                                 continue; /*  not an endpoint element or can't have mobile groups */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                                           /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                                           /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                             nMobile =
    // INCHI✔️❌:                                 AddAtom2num(EndPoint[nNumEndPoints].num, at, endpoint, 2); /* fill out */
    // INCHI✔️❌:                             AddAtom2DA(EndPoint[nNumEndPoints].num_DA, at, endpoint, 2);
    // INCHI✔️❌:                             /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             if (bNonTautBond) {
    // INCHI✔️❌:                                 m = (bond_type == BOND_SINGLE && (nMobile || at[endpoint].endpoint));
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (bond_type == BOND_DOUBLE);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else {
    // INCHI✔️❌:                                 /*  tautomeric or alternating bond */
    // INCHI✔️❌:                                 m = (0 != at[endpoint].endpoint || eif1.cDonor);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (at[endpoint].endpoint ||
    // INCHI✔️❌:                                     eif1.cNeutralBondsValence > at[endpoint].valence);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!bPossiblyEndpoint)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB)endpoint;
    // INCHI✔️❌:                             if (nGroupNumber != at[endpoint].endpoint) {
    // INCHI✔️❌:                                 bDiffGroups++;
    // INCHI✔️❌:                                 nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             if (bNonTautBond || bAltBond) {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                                 BondPos[nNumBondPos].nAtomNumber = (AT_NUMB)centerpoint;
    // INCHI✔️❌:                                 BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                                 nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /*  mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                             /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                             nNumPossibleMobile += (nMobile>0 || at[endpoint].endpoint);
    // INCHI✔️❌:                             nNumEndPoints++;
    // INCHI✔️❌:                             /*printf("Found %d %d %d %d\n", centerpoint+1, at[centerpoint].el_number, endpoint+1, at[endpoint].el_number);*/
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor && num_O == 1 && num_N < 2) {
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             * a tautomeric group has been found
    // INCHI✔️❌:                             *
    // INCHI✔️❌:                             * at this point:
    // INCHI✔️❌:                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             nErr = FindAccessibleEndPoints(pCG,
    // INCHI✔️❌:                                 EndPoint, &nNumEndPoints,
    // INCHI✔️❌:                                 BondPos, &nNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, at,
    // INCHI✔️❌:                                 num_atoms, c_group_info,
    // INCHI✔️❌:                                 ALT_PATH_MODE_TAUTOM_PT_39_00);
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (IS_BNS_ERROR(nErr)) {
    // INCHI✔️❌:                                 return nErr;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nNumEndPoints > 0) {
    // INCHI✔️❌:                                 if (!nGroupNumber || bDiffGroups > 0) {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     num_changes = RegisterEndPoints(pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
    // INCHI✔️❌:                                     if (num_changes == -1) {
    // INCHI✔️❌:                                         nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (num_changes < 0) {
    // INCHI✔️❌:                                         nErr = num_changes;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                                     if (nErr)
    // INCHI✔️❌:                                         goto exit_function;
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nNumBondPos > 0) {
    // INCHI✔️❌:                                     /*  some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                     num_changes = SetTautomericBonds(at, nNumBondPos, BondPos);
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /********** END PT_39_00 ************/
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_PT_13_00 == 1 ) /******* BEGIN PT_13_00 ********/
    // INCHI✔️❌:     if (t_group_info->bTautFlags & TG_FLAG_PT_13_00) {
    // INCHI✔️❌:         /*** [O,S,Se,Te;X1:1]=[C:2]=[C:3][#1:4]>>[#1:4][O,S,Se,Te;X2:1][C:2]#[C:3] ***/
    // INCHI✔️❌:         /*** Similar to the previous case of M=Q-ZH >> MH-Q=Z, with M = "S,O,Se,Te", Q = "C", and Z = "C" ***/
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++) {
    // INCHI✔️❌:             /*  find possible endpoint Z = at[i] */
    // INCHI✔️❌:             if ((endpoint_valence = nGetEndpointInfo_PT_13_00(at, i, &eif1))) { /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                 /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:                 for (j = 0; j < at[i].valence; j++) {
    // INCHI✔️❌:                     bond_type = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                     bond_type = ACTUAL_ORDER(pBNS, i, j, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                     if ((bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                         bond_type == BOND_SINGLE ||
    // INCHI✔️❌:                         bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                         bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                         bond_type == BOND_ALT_13 ||
    // INCHI✔️❌:                         bond_type == BOND_TAUTOM) &&
    // INCHI✔️❌:                         at[centerpoint].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:                         at[centerpoint].valence == 2 &&
    // INCHI✔️❌:                         ALLOWED_EDGE(pBNS, i, j)
    // INCHI✔️❌:                         ) {
    // INCHI✔️❌:                         /*  test a centerpoint candidate. */
    // INCHI✔️❌:                         /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                         int num_O = 0;
    // INCHI✔️❌:                         int num_C = 0;
    // INCHI✔️❌:                         nNumPossibleMobile = 0;
    // INCHI✔️❌:                         nGroupNumber = (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                         bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                         nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                         for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++) {
    // INCHI✔️❌:                             endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                             bond_type = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                             bond_type = ACTUAL_ORDER(pBNS, centerpoint, k, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 bNonTautBond =
    // INCHI✔️❌:                                 bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                                 bPossiblyEndpoint = 0;
    // INCHI✔️❌:                             if (!ALLOWED_EDGE(pBNS, centerpoint, k)) {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_ALT_13 || bond_type == BOND_TAUTOM) {
    // INCHI✔️❌:                                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                     bAltBond = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_ALT_13); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE || bond_type == BOND_TRIPLE)
    // INCHI✔️❌:                                         bNonTautBond = 1;
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:                             if (!(endpoint_valence = nGetEndpointInfo_PT_13_00(at, endpoint, &eif1)))
    // INCHI✔️❌:                                 continue; /*  not an endpoint element or can't have mobile groups */
    // INCHI✔️❌:
    // INCHI✔️❌:                                           /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                                           /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                             nMobile =
    // INCHI✔️❌:                                 AddAtom2num(EndPoint[nNumEndPoints].num, at, endpoint, 2); /* fill out */
    // INCHI✔️❌:                             AddAtom2DA(EndPoint[nNumEndPoints].num_DA, at, endpoint, 2);
    // INCHI✔️❌:                             /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             if (bNonTautBond) {
    // INCHI✔️❌:                                 m = (nMobile || at[endpoint].endpoint);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (nMobile == 0);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else {
    // INCHI✔️❌:                                 /*  tautomeric or alternating bond */
    // INCHI✔️❌:                                 m = (eif1.cDonor != 0);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (eif1.cNeutralBondsValence > at[endpoint].valence);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (!bPossiblyEndpoint)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             num_O += (endpoint_valence == 2);
    // INCHI✔️❌:                             num_C += (endpoint_valence == 4);
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB)endpoint;
    // INCHI✔️❌:                             if (nGroupNumber != at[endpoint].endpoint) {
    // INCHI✔️❌:                                 bDiffGroups++;
    // INCHI✔️❌:                                 nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             if (bNonTautBond || bAltBond) {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                                 BondPos[nNumBondPos].nAtomNumber = (AT_NUMB)centerpoint;
    // INCHI✔️❌:                                 BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                                 nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /*  mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                             /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                             nNumPossibleMobile += (nMobile>0);
    // INCHI✔️❌:                             nNumEndPoints++;
    // INCHI✔️❌:                             /*printf("Found %d %d %d %d\n", centerpoint+1, at[centerpoint].el_number, endpoint+1, at[endpoint].el_number);*/
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor && num_O == 1 && num_C == 1) {
    // INCHI✔️❌:                             /*printf("Real %d\n", nNumEndPoints);*/
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             * a tautomeric group has been found
    // INCHI✔️❌:                             *
    // INCHI✔️❌:                             * at this point:
    // INCHI✔️❌:                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                             */
    // INCHI✔️❌:
    // INCHI✔️❌:                             nErr = FindAccessibleEndPoints(pCG,
    // INCHI✔️❌:                                 EndPoint, &nNumEndPoints,
    // INCHI✔️❌:                                 BondPos, &nNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, at,
    // INCHI✔️❌:                                 num_atoms, c_group_info,
    // INCHI✔️❌:                                 ALT_PATH_MODE_TAUTOM_PT_13_00);
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (IS_BNS_ERROR(nErr)) {
    // INCHI✔️❌:                                 return nErr;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nNumEndPoints > 0) {
    // INCHI✔️❌:                                 if (!nGroupNumber || bDiffGroups > 0) {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     num_changes = RegisterEndPoints(pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
    // INCHI✔️❌:                                     if (num_changes == -1) {
    // INCHI✔️❌:                                         nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (num_changes < 0) {
    // INCHI✔️❌:                                         nErr = num_changes;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                                     if (nErr)
    // INCHI✔️❌:                                         goto exit_function;
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nNumBondPos > 0) {
    // INCHI✔️❌:                                     /*  some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                     num_changes = SetTautomericBonds(at, nNumBondPos, BondPos);
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /********** END PT_13_00 ************/
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_PT_18_00 == 1 ) /******* BEGIN PT_18_00 ********/
    // INCHI✔️❌:     if (t_group_info->bTautFlags & TG_FLAG_PT_18_00) {
    // INCHI✔️❌:         /*** [#1:1][O:2][C:3]#[N:4]>>[O:2]=[C:3]=[N:4][#1:1] ***/
    // INCHI✔️❌:         /*** Similar to the previous case of M=Q-ZH >> MH-Q=Z, with M = "O", Q = "C", and Z = "N" ***/
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++) {
    // INCHI✔️❌:             /*  find possible endpoint Z = at[i] */
    // INCHI✔️❌:             if ((endpoint_valence = nGetEndpointInfo_PT_18_00(at, i, &eif1))) { /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                 /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:                 for (j = 0; j < at[i].valence; j++) {
    // INCHI✔️❌:                     bond_type = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                     bond_type = ACTUAL_ORDER(pBNS, i, j, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                                                            /*printf("Centerpoint: %d\n", centerpoint+1);*/
    // INCHI✔️❌:                     if ((bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                         bond_type == BOND_SINGLE ||
    // INCHI✔️❌:                         bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                         bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                         bond_type == BOND_ALT_13 ||
    // INCHI✔️❌:                         bond_type == BOND_TAUTOM) &&
    // INCHI✔️❌:                         at[centerpoint].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:                         at[centerpoint].valence == 2 &&
    // INCHI✔️❌:                         ALLOWED_EDGE(pBNS, i, j)
    // INCHI✔️❌:                         ) {
    // INCHI✔️❌:                         /*  test a centerpoint candidate. */
    // INCHI✔️❌:                         /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                         int num_O = 0;
    // INCHI✔️❌:                         int num_N = 0;
    // INCHI✔️❌:                         nNumPossibleMobile = 0;
    // INCHI✔️❌:                         nGroupNumber = (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                         bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                         nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                         for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++) {
    // INCHI✔️❌:                             endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                                                                     /*printf("Endpoint: %d\n", endpoint+1);*/
    // INCHI✔️❌:                             bond_type = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                             bond_type = ACTUAL_ORDER(pBNS, centerpoint, k, bond_type);
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 bNonTautBond =
    // INCHI✔️❌:                                 bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                                 bPossiblyEndpoint = 0;
    // INCHI✔️❌:                             if (!ALLOWED_EDGE(pBNS, centerpoint, k)) {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_ALT_13 || bond_type == BOND_TAUTOM) {
    // INCHI✔️❌:                                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                     bAltBond = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_ALT_13); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE || bond_type == BOND_TRIPLE)
    // INCHI✔️❌:                                         bNonTautBond = 1;
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:                             if (!(endpoint_valence = nGetEndpointInfo_PT_18_00(at, endpoint, &eif1)))
    // INCHI✔️❌:                                 continue; /*  not an endpoint element or can't have mobile groups */
    // INCHI✔️❌:
    // INCHI✔️❌:                                           /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                                           /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                             nMobile = eif1.cMobile;
    // INCHI✔️❌:                             AddAtom2num(EndPoint[nNumEndPoints].num, at, endpoint, 2); /* fill out */
    // INCHI✔️❌:                             AddAtom2DA(EndPoint[nNumEndPoints].num_DA, at, endpoint, 2);
    // INCHI✔️❌:                             /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (bond_type == BOND_DOUBLE && endpoint_valence == 3 && nMobile == 0)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (bNonTautBond) {
    // INCHI✔️❌:                                 /*printf("AAA %d\n", nMobile);*/
    // INCHI✔️❌:                                 m = (nMobile || at[endpoint].endpoint);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (nMobile == 0);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else {
    // INCHI✔️❌:                                 /*  tautomeric or alternating bond */
    // INCHI✔️❌:                                 /*printf("BBB %d %d %d\n", eif1.cDonor, eif1.cNeutralBondsValence, at[endpoint].valence);*/
    // INCHI✔️❌:                                 m = (eif1.cDonor != 0);
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = (eif1.cAcceptor != 0);
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (!bPossiblyEndpoint)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                             num_O += (endpoint_valence == 2);
    // INCHI✔️❌:                             num_N += (endpoint_valence == 3);
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB)endpoint;
    // INCHI✔️❌:                             if (nGroupNumber != at[endpoint].endpoint) {
    // INCHI✔️❌:                                 bDiffGroups++;
    // INCHI✔️❌:                                 nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             if (bNonTautBond || bAltBond) {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                                 BondPos[nNumBondPos].nAtomNumber = (AT_NUMB)centerpoint;
    // INCHI✔️❌:                                 BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                                 nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /*  mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                             /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                             nNumPossibleMobile += (nMobile>0);
    // INCHI✔️❌:                             nNumEndPoints++;
    // INCHI✔️❌:                             /*printf("Found %d %d %d %d\n", centerpoint+1, at[centerpoint].el_number, endpoint+1, at[endpoint].el_number);*/
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor && num_O == 1 && num_N == 1) {
    // INCHI✔️❌:                             /*printf("Real %d\n", nNumEndPoints);*/
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             * a tautomeric group has been found
    // INCHI✔️❌:                             *
    // INCHI✔️❌:                             * at this point:
    // INCHI✔️❌:                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                             */
    // INCHI✔️❌:
    // INCHI✔️❌:                             nErr = FindAccessibleEndPoints(pCG,
    // INCHI✔️❌:                                 EndPoint, &nNumEndPoints,
    // INCHI✔️❌:                                 BondPos, &nNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, at,
    // INCHI✔️❌:                                 num_atoms, c_group_info,
    // INCHI✔️❌:                                 ALT_PATH_MODE_TAUTOM_PT_18_00);
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (IS_BNS_ERROR(nErr)) {
    // INCHI✔️❌:                                 return nErr;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nNumEndPoints > 0) {
    // INCHI✔️❌:                                 if (!nGroupNumber || bDiffGroups > 0) {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     num_changes = RegisterEndPoints(pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
    // INCHI✔️❌:                                     if (num_changes == -1) {
    // INCHI✔️❌:                                         nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (num_changes < 0) {
    // INCHI✔️❌:                                         nErr = num_changes;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                                     if (nErr)
    // INCHI✔️❌:                                         goto exit_function;
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nNumBondPos > 0) {
    // INCHI✔️❌:                                     /*  some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                     num_changes = SetTautomericBonds(at, nNumBondPos, BondPos);
    // INCHI✔️❌:                                     tot_changes += (num_changes>0);
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /********** END PT_18_00 ************/
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( KETO_ENOL_TAUT == 1 )  /***** post v.1 feature *****/
    // INCHI✔️❌:     if (t_group_info->bTautFlags & TG_FLAG_KETO_ENOL_TAUT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* 1,3 keto-enol tautomerism */
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  Find possible endpoint Z = at[i] */
    // INCHI✔️❌:             if ((endpoint_valence = nGetEndpointInfo_KET( at, i, &eif1 ))) /* djb-rwth: addressing LLVM warning; ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*  1st endpoint candidate found. Find centerpoint candidate */
    // INCHI✔️❌:                 for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bond_type = (int) at[i].bond_type[j] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                     bond_type = ACTUAL_ORDER( pBNS, i, j, bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     centerpoint = (int) at[i].neighbor[j];  /*  a centerpoint candidate */
    // INCHI✔️❌:                     if (( bond_type == BOND_DOUBLE ||
    // INCHI✔️❌:                           bond_type == BOND_ALTERN ||
    // INCHI✔️❌:                           bond_type == BOND_ALT12NS ||
    // INCHI✔️❌:                           bond_type == BOND_TAUTOM ) &&
    // INCHI✔️❌:                          is_centerpoint_elem_KET( at[centerpoint].el_number ) &&
    // INCHI✔️❌:                          !at[centerpoint].charge && !at[centerpoint].radical &&
    // INCHI✔️❌:                          /* only normal carbon is allowed */
    // INCHI✔️❌:                          4 == at[centerpoint].chem_bonds_valence + at[centerpoint].num_H
    // INCHI✔️❌:                          && ALLOWED_EDGE( pBNS, i, j ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_O = 0;
    // INCHI✔️❌:                         int num_C = 0;
    // INCHI✔️❌:                         /*  Test a centerpoint candidate. */
    // INCHI✔️❌:                         /*  find all endpoints including at[i] and store them into EndPoint[] */
    // INCHI✔️❌:                         nNumPossibleMobile = 0;
    // INCHI✔️❌:                         nGroupNumber = (AT_NUMB) num_atoms; /*  greater than any tautomeric group number */
    // INCHI✔️❌:                         bDiffGroups = -1;         /*  ignore the first difference */
    // INCHI✔️❌:                         nNumDonor = nNumAcceptor = 0;
    // INCHI✔️❌:                         for (k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
    // INCHI✔️❌:                             bond_type = (int) at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:                             bond_type = ACTUAL_ORDER( pBNS, centerpoint, k, bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                 bNonTautBond =
    // INCHI✔️❌:                                 bAltBond = /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:                                 bPossiblyEndpoint = 0;
    // INCHI✔️❌:                             if (!ALLOWED_EDGE( pBNS, centerpoint, k ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( REPLACE_ALT_WITH_TAUT == 1 )
    // INCHI✔️❌:                                     bAltBond = ( bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS ); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         bNonTautBond = 1;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (!( endpoint_valence = nGetEndpointInfo_KET( at, endpoint, &eif2 ) ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             if ( 3 != eif1.cKetoEnolCode + eif2.cKetoEnolCode && endpoint != i )
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             /*  save information about the found possible tautomeric endpoint */
    // INCHI✔️❌:                             /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
    // INCHI✔️❌:                             nMobile = AddAtom2num( EndPoint[nNumEndPoints].num, at, endpoint, 2 ); /* fill out */
    // INCHI✔️❌:                             AddAtom2DA( EndPoint[nNumEndPoints].num_DA, at, endpoint, 2 );
    // INCHI✔️❌:                             /* --- why is isitopic info missing ? -- see below
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
    // INCHI✔️❌:                             nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             if (bNonTautBond)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 m = ( bond_type == BOND_SINGLE && ( nMobile || at[endpoint].endpoint ) );
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = ( bond_type == BOND_DOUBLE );
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /*  Tautomeric or alternating bond */
    // INCHI✔️❌:                                 m = ( 0 != at[endpoint].endpoint || eif1.cDonor );
    // INCHI✔️❌:                                 nNumDonor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                                 m = ( at[endpoint].endpoint ||
    // INCHI✔️❌:                                       eif1.cNeutralBondsValence > at[endpoint].valence );
    // INCHI✔️❌:                                 nNumAcceptor += m;
    // INCHI✔️❌:                                 bPossiblyEndpoint += m;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!bPossiblyEndpoint)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             num_O += ( endpoint_valence == 2 );
    // INCHI✔️❌:                             num_C += ( endpoint_valence == 4 );
    // INCHI✔️❌:
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nGroupNumber = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nEquNumber = 0;
    // INCHI✔️❌:                             EndPoint[nNumEndPoints].nAtomNumber = (AT_NUMB) endpoint;
    // INCHI✔️❌:                             if (nGroupNumber != at[endpoint].endpoint)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 bDiffGroups++;
    // INCHI✔️❌:                                 nGroupNumber = at[endpoint].endpoint;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  Save positions of all, not only possibly tautomeric bonds */
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             if (bNonTautBond || bAltBond)
    // INCHI❌❌:                             {
    // INCHI❌❌: #endif
    // INCHI✔️❌:                                 BondPos[nNumBondPos].nAtomNumber = (AT_NUMB) centerpoint;
    // INCHI✔️❌:                                 BondPos[nNumBondPos].neighbor_index = (AT_NUMB) k; /* bond ordering number; used to change bonds to tautomeric only  */
    // INCHI✔️❌:                                 nNumBondPos++;
    // INCHI❌❌: #if ( REPLACE_ALT_WITH_TAUT != 1 )
    // INCHI❌❌:                             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                             /*  Mobile group is possible if (a) the endpoint has a mobile group or */
    // INCHI✔️❌:                             /*                              (b) the centerpoint is adjacent to another endpoint */
    // INCHI✔️❌:                             nNumPossibleMobile += ( nMobile > 0 || at[endpoint].endpoint );
    // INCHI✔️❌:                             nNumEndPoints++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumEndPoints > 1 && nNumPossibleMobile &&
    // INCHI✔️❌:                              nNumDonor && nNumAcceptor &&
    // INCHI✔️❌:                              num_O == 1 && num_C)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             * A tautomeric group has been found
    // INCHI✔️❌:                             *
    // INCHI✔️❌:                             * at this point:
    // INCHI✔️❌:                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
    // INCHI✔️❌:                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
    // INCHI✔️❌:                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
    // INCHI✔️❌:                             */
    // INCHI✔️❌:
    // INCHI✔️❌:                             nErr = FindAccessibleEndPoints( pCG, EndPoint, &nNumEndPoints, BondPos, &nNumBondPos,
    // INCHI✔️❌:                                                             pBNS, pBD, at, num_atoms, c_group_info, ALT_PATH_MODE_TAUTOM_KET );
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (IS_BNS_ERROR( nErr ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 return nErr;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nNumEndPoints > 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (!nGroupNumber || bDiffGroups > 0)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     num_changes = RegisterEndPoints( pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:                                     if (num_changes == -1)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (num_changes < 0)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         nErr = num_changes;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (nErr)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         goto exit_function;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     tot_changes += ( num_changes > 0 );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (nNumBondPos > 0)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     /*  Some of the bonds have not been marked as tautomeric yet */
    // INCHI✔️❌:                                     num_changes = SetTautomericBonds( at, nNumBondPos, BondPos );
    // INCHI✔️❌:                                     tot_changes += ( num_changes > 0 );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif  /* KETO_ENOL_TAUT */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_OTHER == 1 ) /* { */
    // INCHI✔️❌:     if (!tot_changes)
    // INCHI✔️❌:     {
    // INCHI✔️❌: #define MAX_ALT_PATH_LEN 8
    // INCHI✔️❌:         int nMaxLenDfsPath = MAX_ALT_PATH_LEN;
    // INCHI✔️❌:         int i1, i2;
    // INCHI✔️❌:         AT_RANK *nDfsPathPos = (AT_RANK  *) inchi_calloc( num_atoms, sizeof( nDfsPathPos[0] ) );
    // INCHI✔️❌:         DFS_PATH DfsPath[MAX_ALT_PATH_LEN];
    // INCHI✔️❌:         int      ret;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!nDfsPathPos) /* djb-rwth: removing redundant code as address of DfsPath will always evaluate to true */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             tot_changes = CT_OUT_OF_RAM;  /*   <BRKPT> */
    // INCHI✔️❌:             goto free_memory;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_15_NON_RING      == 1 ) /***** post v.1 feature *****/
    // INCHI✔️❌:         if (t_group_info->bTautFlags & TG_FLAG_1_5_TAUT)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  1,5 tautomerism; one of the endpoints should no be on a ring  */
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             O                OH                 O
    // INCHI✔️❌:             ||               |                  ||
    // INCHI✔️❌:             A--pos-          A--pos-            A--pos-
    // INCHI✔️❌:             /   sib-        //   sib-     ?     /   sib-
    // INCHI✔️❌:             C    ly          C    ly            CH   ly
    // INCHI✔️❌:             \\   a     <-->   \   a       <-->   \   a
    // INCHI✔️❌:             B--ring          B--ring            B--ring
    // INCHI✔️❌:             |                ||                 ||
    // INCHI✔️❌:             NH               N                  N
    // INCHI✔️❌:
    // INCHI✔️❌:             Note: few recent modifications now allow the terminal N be in a ring, too
    // INCHI✔️❌:             */
    // INCHI✔️❌:             for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*  Find possible endpoint Z = at[i1] */
    // INCHI✔️❌:                 if (!( endpoint_valence = nGetEndpointInfo( at, i1, &eif1 ) ) /*||
    // INCHI✔️❌:                                                                               at[i1].nNumAtInRingSystem > 1*/) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue; /*  not a possibly endpoint */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nNumEndPoints = 0;
    // INCHI✔️❌:                     nNumBondPos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                     ret = nGet15TautInAltPath( pCG, at, i1, nDfsPathPos,
    // INCHI✔️❌:                                                DfsPath, nMaxLenDfsPath,
    // INCHI✔️❌:                                                EndPoint, sizeof( EndPoint ) / sizeof( EndPoint[0] ),
    // INCHI✔️❌:                                                BondPos, sizeof( BondPos ) / sizeof( BondPos[0] ),
    // INCHI✔️❌:                                                &nNumEndPoints, &nNumBondPos,
    // INCHI✔️❌:                                                pBNS, pBD, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (ret > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nNumEndPoints)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             num_changes = RegisterEndPoints( pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:                             if (num_changes == -1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (num_changes < 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 nErr = num_changes;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (nErr)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 goto free_memory;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             tot_changes += ( num_changes > 0 );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumBondPos)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nErr = ret;
    // INCHI✔️❌:                             goto free_memory;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_4PYRIDINOL_RINGS == 1 )
    // INCHI✔️❌:         /*  6-member rings */
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         O              OH             OH
    // INCHI✔️❌:         ||             |              |
    // INCHI✔️❌:         /  \          //  \           /  \\
    // INCHI✔️❌:         ||   ||  <-->  |    ||  <-->  ||   |
    // INCHI✔️❌:         \  /          \\  /           \  //
    // INCHI✔️❌:         NH             N              N
    // INCHI✔️❌:         */
    // INCHI✔️❌:         for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  Find possible endpoint Z = at[i1] */
    // INCHI✔️❌:             if (3 != ( endpoint_valence = nGetEndpointInfo( at, i1, &eif1 ) ) ||
    // INCHI✔️❌:                  2 != at[i1].valence) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue; /*  not a nitrogen atom or a wrong valence */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (at[i1].nNumAtInRingSystem >= 6)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nNumEndPoints = 0;
    // INCHI✔️❌:                 nNumBondPos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                 ret = nGet15TautIn6MembAltRing( pCG, at, i1, nDfsPathPos,
    // INCHI✔️❌:                                                 DfsPath, nMaxLenDfsPath, EndPoint,
    // INCHI✔️❌:                                                 sizeof( EndPoint ) / sizeof( EndPoint[0] ),
    // INCHI✔️❌:                                                 BondPos,
    // INCHI✔️❌:                                                 sizeof( BondPos ) / sizeof( BondPos[0] ),
    // INCHI✔️❌:                                                 &nNumEndPoints, &nNumBondPos,
    // INCHI✔️❌:                                                 pBNS, pBD, num_atoms );
    // INCHI✔️❌:                 if (ret > 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (nNumEndPoints)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_changes = RegisterEndPoints( pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:                         if (num_changes == -1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (num_changes < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nErr = num_changes;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nErr)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto free_memory;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         tot_changes += ( num_changes > 0 );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (nNumBondPos)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nErr = ret;
    // INCHI✔️❌:                         goto free_memory;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif /* TAUT_4PYRIDINOL_RINGS */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_PYRAZOLE_RINGS == 1 )
    // INCHI✔️❌:         /* 5-member rings:
    // INCHI✔️❌:
    // INCHI✔️❌:         Z               Z
    // INCHI✔️❌:         /  \\           //  \
    // INCHI✔️❌:         X     Y  <-->   X     Y
    // INCHI✔️❌:         \\   /          \    //
    // INCHI✔️❌:         N--NH           HN--N
    // INCHI✔️❌:
    // INCHI✔️❌:         ^             ^
    // INCHI✔️❌:         search for these NH
    // INCHI✔️❌:         */
    // INCHI✔️❌:         /*  5-member rings (pyrazole derivatives): look for the neighboring N */
    // INCHI✔️❌:         for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (2 == at[i1].valence &&
    // INCHI✔️❌:                  at[i1].nNumAtInRingSystem >= 5 &&
    // INCHI✔️❌:                  3 == ( endpoint_valence = nGetEndpointInfo( at, i1, &eif1 ) )
    // INCHI✔️❌:                  ) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nMobile = at[i1].num_H + ( at[i1].charge == -1 );
    // INCHI✔️❌:                 for (j = 0; j < at[i1].valence; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int nMobile2, endpoint_valence2;
    // INCHI✔️❌:                     i2 = at[i1].neighbor[j];
    // INCHI✔️❌:
    // INCHI✔️❌:                     /*  may be important */
    // INCHI✔️❌:                     if (i2 >= i1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue; /*  do not try same pair 2 times */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (at[i2].nRingSystem != at[i1].nRingSystem)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     bond_type = ( at[i1].bond_type[j] & ~BOND_MARK_ALL );
    // INCHI✔️❌:                     if ((bond_type != BOND_SINGLE &&
    // INCHI✔️❌:                          bond_type != BOND_TAUTOM &&
    // INCHI✔️❌:                          bond_type != BOND_ALT12NS &&
    // INCHI✔️❌:                          bond_type != BOND_ALTERN) ||  /* added 1-15-2002 */
    // INCHI✔️❌:                          2 != at[i2].valence ||
    // INCHI✔️❌:                          3 != ( endpoint_valence2 = nGetEndpointInfo( at, i2, &eif2 ) )) /* djb-rwth: addressing LLVM warning; variable used to store function return value */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue; /*  not a nitrogen atom or a wrong valence or not a single bond */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     nMobile2 = at[i2].num_H + ( at[i2].charge == -1 );  /*  number of mobile groups */
    // INCHI❌❌: #if ( TAUT_IGNORE_EQL_ENDPOINTS == 1 )
    // INCHI❌❌:                     if (at[i1].endpoint && at[i1].endpoint == at[i2].endpoint)
    // INCHI❌❌:                     {
    // INCHI❌❌:                         continue; /* atoms already belong to the same t-group */
    // INCHI❌❌:                     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     if (!at[i1].endpoint && !at[i2].endpoint && 1 != nMobile + nMobile2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     ret = nGet12TautIn5MembAltRing( pCG, at, i1, j, nDfsPathPos,
    // INCHI✔️❌:                                                     DfsPath, nMaxLenDfsPath,
    // INCHI✔️❌:                                                     EndPoint,
    // INCHI✔️❌:                                                     sizeof( EndPoint ) / sizeof( EndPoint[0] ),
    // INCHI✔️❌:                                                     BondPos,
    // INCHI✔️❌:                                                     sizeof( BondPos ) / sizeof( BondPos[0] ),
    // INCHI✔️❌:                                                     &nNumEndPoints, &nNumBondPos,
    // INCHI✔️❌:                                                     pBNS, pBD, num_atoms );
    // INCHI✔️❌:                     if (ret > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nNumEndPoints)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             num_changes = RegisterEndPoints( pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:                             if (num_changes == -1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (num_changes < 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 nErr = num_changes;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (nErr)
    // INCHI✔️❌:                                 goto free_memory;
    // INCHI✔️❌:                             tot_changes += ( num_changes > 0 );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nNumBondPos)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nErr = ret;
    // INCHI✔️❌:                             goto free_memory;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif /* TAUT_PYRAZOLE_RINGS */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_TROPOLONE_7 == 1 || TAUT_TROPOLONE_5 == 1 ) /* { */
    // INCHI✔️❌:         /********************************************************
    // INCHI✔️❌:         *                                         A  B
    // INCHI✔️❌:         *                                         | ||
    // INCHI✔️❌:         * 7-member rings (tropolones): look for M=Q--R--ZH,
    // INCHI✔️❌:         *                                       ^ ^  ^  ^
    // INCHI✔️❌:         *                               endpoint1 i1 i2 endpoint2
    // INCHI✔️❌:         * where A-Q-R=B belong to a 7-member alt. (except Q-R bond) ring: ..=A-(Q-R)=B-..
    // INCHI✔️❌:         * Bond Q-R should be single or tautomeric or alternating
    // INCHI✔️❌:         * M=Q and R-ZH should be chain (non-ring) bonds
    // INCHI✔️❌:         * Same for 5-member rings
    // INCHI✔️❌:         */
    // INCHI✔️❌:         for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[i1].nNumAtInRingSystem >=
    // INCHI✔️❌: #if ( TAUT_TROPOLONE_5 == 1 )
    // INCHI✔️❌:                  5
    // INCHI❌❌: #else
    // INCHI❌❌:                  7
    // INCHI❌❌: #endif
    // INCHI✔️❌:                  &&
    // INCHI✔️❌:                  bIsCenterPointStrict( at, i1 ) &&
    // INCHI✔️❌: #if ( TAUT_RINGS_ATTACH_CHAIN == 1 )
    // INCHI✔️❌:                  at[i1].bCutVertex &&
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                  at[i1].valence == 3 && !at[i1].endpoint)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int nMobile1, endpoint1, endpoint1_valence, bond_type1;
    // INCHI✔️❌:                 int nMobile2, endpoint2, endpoint2_valence, bond_type2;
    // INCHI✔️❌:                 for (j = 0; j < at[i1].valence; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     i2 = at[i1].neighbor[j];
    // INCHI✔️❌:                     /* may be important
    // INCHI✔️❌:                     if ( i2 > i1 )
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                     do not try same pair 2 times
    // INCHI✔️❌:                     */
    // INCHI✔️❌:                     if (at[i2].nRingSystem != at[i1].nRingSystem ||
    // INCHI✔️❌:                          !bIsCenterPointStrict( at, i2 ) ||
    // INCHI✔️❌: #if ( TAUT_RINGS_ATTACH_CHAIN == 1 )
    // INCHI✔️❌:                          !at[i2].bCutVertex ||
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                          at[i2].valence != 3 || at[i2].endpoint)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     bond_type = ( at[i1].bond_type[j] & ~BOND_MARK_ALL );
    // INCHI✔️❌:                     if (bond_type != BOND_SINGLE &&
    // INCHI✔️❌:                          bond_type != BOND_TAUTOM &&
    // INCHI✔️❌:                          bond_type != BOND_ALT12NS &&
    // INCHI✔️❌:                          bond_type != BOND_ALTERN)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue; /*  not a single bond between Q-R */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /*  Find endpoints */
    // INCHI✔️❌:                     for (k = 0; k < at[i1].valence; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         endpoint1 = at[i1].neighbor[k];
    // INCHI✔️❌:                         if (endpoint1 == i2)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue; /*  j == k */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (!( endpoint1_valence = nGetEndpointInfo( at, endpoint1, &eif1 ) ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue; /*  not an endpoint1 element or can't have mobile groups */
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #if ( TAUT_RINGS_ATTACH_CHAIN == 1 )
    // INCHI✔️❌:                         if (at[endpoint1].nRingSystem == at[i1].nRingSystem)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         nMobile1 = at[endpoint1].num_H + ( at[endpoint1].charge == -1 );  /*  number of mobile groups */
    // INCHI✔️❌:                         if (nMobile1 + at[endpoint1].chem_bonds_valence != endpoint1_valence)
    // INCHI✔️❌:                             continue; /*  abnormal endpoint1 valence; ignore. */
    // INCHI✔️❌:                         bond_type1 = ( at[i1].bond_type[k] & ~BOND_MARK_ALL );
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (bond_type1 != BOND_SINGLE &&
    // INCHI✔️❌:                              bond_type1 != BOND_DOUBLE &&
    // INCHI✔️❌:                              bond_type1 != BOND_TAUTOM &&
    // INCHI✔️❌:                              bond_type1 != BOND_ALT12NS &&
    // INCHI✔️❌:                              bond_type1 != BOND_ALTERN)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         for (m = 0; m < at[i2].valence; m++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             endpoint2 = at[i2].neighbor[m];
    // INCHI✔️❌:                             if (endpoint2 == i1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!( endpoint2_valence = nGetEndpointInfo( at, endpoint2, &eif2 ) )) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue; /*  not an endpoint2 element or can't have mobile groups */
    // INCHI✔️❌:                             }
    // INCHI✔️❌: #if ( TAUT_RINGS_ATTACH_CHAIN == 1 )
    // INCHI✔️❌:                             if (at[endpoint2].nRingSystem == at[i2].nRingSystem)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                             nMobile2 = at[endpoint2].num_H + ( at[endpoint2].charge == -1 );  /*  number of mobile groups */
    // INCHI✔️❌:                             bond_type2 = ( at[i2].bond_type[m] & ~BOND_MARK_ALL );
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (bond_type2 != BOND_SINGLE &&
    // INCHI✔️❌:                                  bond_type2 != BOND_DOUBLE &&
    // INCHI✔️❌:                                  bond_type2 != BOND_TAUTOM &&
    // INCHI✔️❌:                                  bond_type2 != BOND_ALT12NS &&
    // INCHI✔️❌:                                  bond_type2 != BOND_ALTERN)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  Final test for possible tautomerism */
    // INCHI✔️❌:                             nMobile = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (ALLOWED_EDGE( pBNS, i1, k ) && ALLOWED_EDGE( pBNS, i2, m ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /*  Can mobile group move from 1 to 2? */
    // INCHI✔️❌:                                 nMobile += ( at[endpoint1].endpoint || nMobile1 ) &&  /*  from endpoint1 */
    // INCHI✔️❌:                                     ( bond_type1 != BOND_DOUBLE ) &&
    // INCHI✔️❌:
    // INCHI✔️❌:                                     ( at[endpoint2].endpoint ||          /*  to endpoint2 */
    // INCHI✔️❌:                                       eif2.cNeutralBondsValence > at[endpoint2].valence ) &&
    // INCHI✔️❌:                                       ( bond_type2 != BOND_SINGLE );
    // INCHI✔️❌:
    // INCHI✔️❌:                                 /*  Can mobile group move from 2 to 1? */
    // INCHI✔️❌:                                 nMobile += ( at[endpoint2].endpoint || nMobile2 ) &&  /*  from endpoint2 */
    // INCHI✔️❌:                                     ( bond_type2 != BOND_DOUBLE ) && /*changed from BOND_SINGLE 2004-02-26 */
    // INCHI✔️❌:
    // INCHI✔️❌:                                     ( at[endpoint1].endpoint ||          /*  to endpoint1 */
    // INCHI✔️❌:                                       eif1.cNeutralBondsValence > at[endpoint1].valence ) &&
    // INCHI✔️❌:                                       ( bond_type1 != BOND_SINGLE );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!nMobile)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (bond_type1 == bond_type2 &&
    // INCHI✔️❌:                                 ( bond_type1 == BOND_SINGLE || bond_type1 == BOND_DOUBLE ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             /* -- old --
    // INCHI✔️❌:                             if ( !at[endpoint1].endpoint && !at[endpoint2].endpoint && 1 != nMobile1 + nMobile2 )
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             /* -- new --
    // INCHI✔️❌:
    // INCHI✔️❌:                             if ( !at[endpoint1].endpoint && !at[endpoint2].endpoint ) {
    // INCHI✔️❌:                             if ( !(bond_type1 == BOND_SINGLE || bond_type1 == BOND_DOUBLE) ||
    // INCHI✔️❌:                             !(bond_type2 == BOND_SINGLE || bond_type2 == BOND_DOUBLE) ) {
    // INCHI✔️❌:                             // at this point bond_type1 != bond_type2
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if ( bond_type1 == BOND_SINGLE && !nMobile1 ||
    // INCHI✔️❌:                             bond_type2 == BOND_SINGLE && !nMobile2 ||
    // INCHI✔️❌:                             0 == nMobile1 + nMobile2 ) {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_TROPOLONE_7 == 1 )
    // INCHI✔️❌:                             if (at[i1].nNumAtInRingSystem >= 7)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:
    // INCHI✔️❌:                                 ret = nGet14TautIn7MembAltRing( pCG, at, i1,
    // INCHI✔️❌:                                                                 j, k, m, nDfsPathPos,
    // INCHI✔️❌:                                                                 DfsPath, nMaxLenDfsPath,
    // INCHI✔️❌:                                                                 EndPoint,
    // INCHI✔️❌:                                                                 sizeof( EndPoint ) / sizeof( EndPoint[0] ),
    // INCHI✔️❌:                                                                 BondPos,
    // INCHI✔️❌:                                                                 sizeof( BondPos ) / sizeof( BondPos[0] ),
    // INCHI✔️❌:                                                                 &nNumEndPoints, &nNumBondPos,
    // INCHI✔️❌:                                                                 pBNS, pBD, num_atoms );
    // INCHI✔️❌:                                 if (ret > 0)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (nNumEndPoints)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         num_changes = RegisterEndPoints( pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:                                         if (num_changes == -1)
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         if (num_changes < 0)
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             nErr = num_changes;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         if (nErr)
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             goto free_memory;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         tot_changes += ( num_changes > 0 );
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (nNumBondPos)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         nErr = ret;
    // INCHI✔️❌:                                         goto free_memory;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( TAUT_TROPOLONE_5 == 1 )
    // INCHI✔️❌:                             if (at[i1].nNumAtInRingSystem >= 5)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:
    // INCHI✔️❌:                                 ret = nGet14TautIn5MembAltRing( pCG, at, i1, j, k, m, nDfsPathPos,
    // INCHI✔️❌:                                                                 DfsPath, nMaxLenDfsPath,
    // INCHI✔️❌:                                                                 EndPoint, sizeof( EndPoint ) / sizeof( EndPoint[0] ),
    // INCHI✔️❌:                                                                 BondPos, sizeof( BondPos ) / sizeof( BondPos[0] ),
    // INCHI✔️❌:                                                                 &nNumEndPoints, &nNumBondPos,
    // INCHI✔️❌:                                                                 pBNS, pBD, num_atoms );
    // INCHI✔️❌:                                 if (ret > 0)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (nNumEndPoints)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         num_changes = RegisterEndPoints( pCG, t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
    // INCHI✔️❌:                                         if (num_changes == -1)
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             nErr = CT_TAUCOUNT_ERR;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         if (num_changes < 0)
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             nErr = num_changes;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         if (nErr)
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             goto free_memory;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         tot_changes += ( num_changes > 0 );
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     if (nNumBondPos)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (IS_BNS_ERROR( ret ))
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         nErr = ret;
    // INCHI✔️❌:                                         goto free_memory;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif /* } TAUT_TROPOLONE */
    // INCHI✔️❌:
    // INCHI✔️❌:     free_memory:
    // INCHI✔️❌:         if (nDfsPathPos)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( nDfsPathPos );
    // INCHI✔️❌:         }
    // INCHI✔️❌: #undef MAX_ALT_PATH_LEN
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif  /* } FIND_RING_SYSTEMS */
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return nErr < 0 ? nErr : tot_changes;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkTautomerGroups
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MarkTautomerGroups
    // INCHI✔️❌: #define FIX_BOND23_IN_TAUT 0
    // INCHI✔️❌: #define REPLACE_ALT_WITH_TAUT 1
    // INCHI✔️❌: #define TAUT_PT_22_00 1
    // INCHI✔️❌: #define TAUT_PT_16_00 1
    // INCHI✔️❌: #define TAUT_PT_06_00 1
    // INCHI✔️❌: #define TAUT_PT_39_00 1
    // INCHI✔️❌: #define TAUT_PT_13_00 1
    // INCHI✔️❌: #define TAUT_PT_18_00 1
    // INCHI✔️❌: #define KETO_ENOL_TAUT 1
    // INCHI✔️❌: #define TAUT_OTHER 1
    // INCHI✔️❌: #define TAUT_15_NON_RING 1
    // INCHI✔️❌: #define TAUT_4PYRIDINOL_RINGS 1
    // INCHI✔️❌: #define TAUT_PYRAZOLE_RINGS 1
    // INCHI✔️❌: #define TAUT_IGNORE_EQL_ENDPOINTS 0
    // INCHI✔️❌: #define TAUT_TROPOLONE_7 1
    // INCHI✔️❌: #define TAUT_TROPOLONE_5 1
    // INCHI✔️❌: #define TAUT_RINGS_ATTACH_CHAIN 1
    // INCHI✔️❌: #define MAX_ALT_PATH_LEN 8
    // INCHI✔️❌: #define MAXVAL 20
    // END INCHI ACTIVE MACRO CONFIGURATION: MarkTautomerGroups

    let Some(t_group_info) = t_group_info else {
        return Ok(0);
    };
    if t_group_info.bTautFlags & TG_FLAG_TEST_TAUT__ATOMS as u64 == 0 {
        return Ok(0);
    }

    if t_group_info.t_group.is_null() && t_group_info.max_num_t_groups == 0 {
        let b_taut_flags = t_group_info.bTautFlags;
        let b_taut_flags_done = t_group_info.bTautFlagsDone;
        let tni = t_group_info.tni.clone();
        let t_group_number = t_group_info.tGroupNumber;
        let b_ignore_isotopic = t_group_info.bIgnoreIsotopic;
        *t_group_info = T_GROUP_INFO::default();
        t_group_info.bIgnoreIsotopic = b_ignore_isotopic;
        t_group_info.bTautFlags = b_taut_flags;
        t_group_info.bTautFlagsDone = b_taut_flags_done;
        t_group_info.tni = tni;
        t_group_info.tGroupNumber = t_group_number;
        t_group_info.max_num_t_groups = num_atoms.wrapping_div(2).wrapping_add(1);
        let count = i64::from(t_group_info.max_num_t_groups).wrapping_add(1) as u64;
        match inchi_calloc::<T_GROUP>(heap, count, 40) {
            Ok(groups) => t_group_info.t_group = groups,
            Err(
                SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange,
            ) => {
                t_group_info.max_num_t_groups = -1;
                t_group_info.t_group = SourceMutPointer::null();
                return Ok(-1);
            }
            Err(error) => return Err(error),
        }
    }

    if t_group_info.t_group.is_null() || t_group_info.max_num_t_groups == 0 {
        return Ok(0);
    }
    if t_group_info.max_num_t_groups < 0 {
        return Ok(t_group_info.max_num_t_groups);
    }

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Rule {
        Standard,
        Pt22,
        Pt16,
        Pt06,
        Pt39,
        Pt13,
        Pt18,
    }

    macro_rules! endpoint_info {
        ($rule:expr, $atoms:expr, $index:expr, $eif:expr) => {
            match $rule {
                Rule::Standard => nGetEndpointInfo($atoms, $index, $eif),
                Rule::Pt22 => nGetEndpointInfo_PT_22_00($atoms, $index, $eif),
                Rule::Pt16 => nGetEndpointInfo_PT_16_00($atoms, $index, $eif),
                Rule::Pt06 => nGetEndpointInfo_PT_06_00($atoms, $index, $eif),
                Rule::Pt39 => nGetEndpointInfo_PT_39_00($atoms, $index, $eif),
                Rule::Pt13 => nGetEndpointInfo_PT_13_00($atoms, $index, $eif),
                Rule::Pt18 => nGetEndpointInfo_PT_18_00($atoms, $index, $eif),
            }
        };
    }

    macro_rules! allowed_edge {
        ($atom_number:expr, $bond_number:expr) => {{
            if pBNS.edge.is_null() || pBNS.vert.is_null() {
                true
            } else {
                let vertex = heap
                    .slice(pBNS.vert.as_const())?
                    .get(
                        usize::try_from($atom_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let edge_number = *heap
                    .slice(vertex.iedge.as_const())?
                    .get(
                        usize::try_from($bond_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                heap.slice(pBNS.edge.as_const())?
                    .get(
                        usize::try_from(edge_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .forbidden
                    == 0
            }
        }};
    }

    let is_bns_error = |value: i32| BNS_ERR <= value && value <= BNS_MAX_ERR_VALUE;
    let mut tot_changes = 0_i32;
    let mut n_err = 0_i32;

    macro_rules! register_and_mark {
        ($endpoints:expr, $endpoint_count:expr, $bonds:expr, $bond_count:expr,
         $group_number:expr, $different_groups:expr, $mode:expr) => {{
            let result = FindAccessibleEndPoints(
                heap,
                pCG,
                $endpoints,
                $endpoint_count,
                $bonds,
                $bond_count,
                pBNS,
                pBD,
                at,
                num_atoms,
                c_group_info.as_deref(),
                $mode,
                clock_result,
            )?;
            if is_bns_error(result) {
                return Ok(result);
            }
            n_err = 0;
            if *$endpoint_count > 0 {
                if $group_number == 0 || $different_groups > 0 {
                    let num_changes = RegisterEndPoints(
                        heap,
                        pCG,
                        t_group_info,
                        $endpoints,
                        *$endpoint_count,
                        at,
                        num_atoms,
                        c_group_info.as_deref_mut(),
                        Some(&mut *pBNS),
                    )?;
                    if num_changes == -1 {
                        n_err = CT_TAUCOUNT_ERR;
                    }
                    if num_changes < 0 {
                        n_err = num_changes;
                    }
                    if n_err != 0 {
                        return Ok(n_err);
                    }
                    tot_changes = tot_changes.wrapping_add(i32::from(num_changes > 0));
                }
                if *$bond_count > 0 {
                    let num_changes =
                        SetTautomericBonds(heap.slice_mut(at)?, *$bond_count, $bonds)?;
                    tot_changes = tot_changes.wrapping_add(i32::from(num_changes > 0));
                }
            }
        }};
    }

    macro_rules! scan_rule {
        ($rule:expr, $mode:expr) => {{
            let rule = $rule;
            let mut i = 0_i32;
            while i < num_atoms {
                let mut eif1 = ENDPOINT_INFO::default();
                let endpoint_valence = {
                    let atoms = heap.slice(at.as_const())?;
                    endpoint_info!(rule, atoms, i, &mut eif1)
                };
                if endpoint_valence != 0 {
                    let initial_atom = heap
                        .slice(at.as_const())?
                        .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let mut j = 0_i32;
                    while j < i32::from(initial_atom.valence) {
                        let j_index =
                            usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let bond_type =
                            i32::from(initial_atom.bond_type[j_index]) & !(BOND_MARK_ALL as i32);
                        let centerpoint = i32::from(initial_atom.neighbor[j_index]);
                        let center = heap
                            .slice(at.as_const())?
                            .get(
                                usize::try_from(centerpoint)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone();
                        let broad_center_bonds = rule == Rule::Pt13 || rule == Rule::Pt18;
                        let center_bond_ok = bond_type == BOND_DOUBLE as i32
                            || (broad_center_bonds && bond_type == BOND_SINGLE as i32)
                            || bond_type == BOND_ALTERN as i32
                            || bond_type == BOND_ALT12NS as i32
                            || (broad_center_bonds && bond_type == BOND_ALT_13 as i32)
                            || bond_type == BOND_TAUTOM as i32;
                        let center_ok = match rule {
                            Rule::Standard => is_centerpoint_elem(center.el_number) != 0,
                            Rule::Pt22 => center.el_number == 7,
                            Rule::Pt16 => {
                                center.el_number == 7 && center.valence == 2 && center.charge == 0
                            }
                            Rule::Pt06 => matches!(center.el_number, 6 | 7 | 15),
                            Rule::Pt39 => center.el_number == 7 && center.valence == 3,
                            Rule::Pt13 | Rule::Pt18 => center.el_number == 6 && center.valence == 2,
                        };
                        if center_bond_ok && center_ok && allowed_edge!(i, j) {
                            let mut endpoints =
                                std::array::from_fn::<_, { MAXVAL as usize }, _>(|_| {
                                    T_ENDPOINT::default()
                                });
                            let mut bonds =
                                std::array::from_fn::<_, { MAXVAL as usize }, _>(|_| {
                                    T_BONDPOS::default()
                                });
                            let mut endpoint_count = 0_i32;
                            let mut bond_count = 0_i32;
                            let mut possible_mobile = 0_i32;
                            let mut group_number = num_atoms as AT_NUMB;
                            let mut different_groups = -1_i32;
                            let mut donors = 0_i32;
                            let mut acceptors = 0_i32;
                            let mut num_o = 0_i32;
                            let mut num_n = 0_i32;
                            let mut num_c = 0_i32;
                            let mut k = 0_i32;
                            while k < i32::from(center.valence) {
                                let k_index = usize::try_from(k)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let endpoint = i32::from(center.neighbor[k_index]);
                                let endpoint_atom = heap
                                    .slice(at.as_const())?
                                    .get(
                                        usize::try_from(endpoint)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .clone();
                                if rule == Rule::Pt39 {
                                    num_o =
                                        num_o.wrapping_add(i32::from(endpoint_atom.el_number == 8));
                                    num_n =
                                        num_n.wrapping_add(i32::from(endpoint_atom.el_number == 7));
                                }
                                let endpoint_bond_type =
                                    i32::from(center.bond_type[k_index]) & !(BOND_MARK_ALL as i32);
                                let broad_endpoint_bonds = rule == Rule::Pt13 || rule == Rule::Pt18;
                                let alt_bond = endpoint_bond_type == BOND_ALTERN as i32
                                    || endpoint_bond_type == BOND_ALT12NS as i32
                                    || (broad_endpoint_bonds
                                        && endpoint_bond_type == BOND_ALT_13 as i32)
                                    || endpoint_bond_type == BOND_TAUTOM as i32;
                                let non_taut_bond = endpoint_bond_type == BOND_SINGLE as i32
                                    || endpoint_bond_type == BOND_DOUBLE as i32
                                    || (broad_endpoint_bonds
                                        && endpoint_bond_type == BOND_TRIPLE as i32);
                                if !allowed_edge!(centerpoint, k) || (!alt_bond && !non_taut_bond) {
                                    k = k.wrapping_add(1);
                                    continue;
                                }
                                let endpoint_valence = {
                                    let atoms = heap.slice(at.as_const())?;
                                    endpoint_info!(rule, atoms, endpoint, &mut eif1)
                                };
                                if endpoint_valence == 0 {
                                    k = k.wrapping_add(1);
                                    continue;
                                }
                                if rule == Rule::Pt16
                                    && ((endpoint_atom.el_number == 8
                                        && endpoint_atom.nNumAtInRingSystem != 1)
                                        || (endpoint != i
                                            && endpoint_atom.el_number == initial_atom.el_number))
                                {
                                    k = k.wrapping_add(1);
                                    continue;
                                }
                                if rule == Rule::Pt06
                                    && endpoint != i
                                    && endpoint_atom.el_number == 6
                                    && initial_atom.el_number == 6
                                {
                                    k = k.wrapping_add(1);
                                    continue;
                                }

                                let endpoint_slot = usize::try_from(endpoint_count)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let mobile = if rule == Rule::Pt18 {
                                    let value = i32::from(eif1.cMobile);
                                    let _ = AddAtom2num(
                                        &mut endpoints
                                            .get_mut(endpoint_slot)
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                                            .num,
                                        heap.slice(at.as_const())?,
                                        endpoint,
                                        2,
                                    )?;
                                    value
                                } else {
                                    AddAtom2num(
                                        &mut endpoints
                                            .get_mut(endpoint_slot)
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                                            .num,
                                        heap.slice(at.as_const())?,
                                        endpoint,
                                        2,
                                    )?
                                };
                                AddAtom2DA(
                                    &mut endpoints
                                        .get_mut(endpoint_slot)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                        .num_DA,
                                    heap.slice(at.as_const())?,
                                    endpoint,
                                    2,
                                )?;
                                if rule == Rule::Pt18
                                    && endpoint_bond_type == BOND_DOUBLE as i32
                                    && endpoint_valence == 3
                                    && mobile == 0
                                {
                                    k = k.wrapping_add(1);
                                    continue;
                                }

                                let (donor, acceptor) = match rule {
                                    Rule::Pt13 => {
                                        if non_taut_bond {
                                            (
                                                i32::from(
                                                    mobile != 0 || endpoint_atom.endpoint != 0,
                                                ),
                                                i32::from(mobile == 0),
                                            )
                                        } else {
                                            (
                                                i32::from(eif1.cDonor != 0),
                                                i32::from(
                                                    i32::from(eif1.cNeutralBondsValence)
                                                        > i32::from(endpoint_atom.valence),
                                                ),
                                            )
                                        }
                                    }
                                    Rule::Pt18 => {
                                        if non_taut_bond {
                                            (
                                                i32::from(
                                                    mobile != 0 || endpoint_atom.endpoint != 0,
                                                ),
                                                i32::from(mobile == 0),
                                            )
                                        } else {
                                            (
                                                i32::from(eif1.cDonor != 0),
                                                i32::from(eif1.cAcceptor != 0),
                                            )
                                        }
                                    }
                                    _ => {
                                        if non_taut_bond {
                                            (
                                                i32::from(
                                                    endpoint_bond_type == BOND_SINGLE as i32
                                                        && (mobile != 0
                                                            || endpoint_atom.endpoint != 0),
                                                ),
                                                i32::from(endpoint_bond_type == BOND_DOUBLE as i32),
                                            )
                                        } else {
                                            (
                                                i32::from(
                                                    endpoint_atom.endpoint != 0 || eif1.cDonor != 0,
                                                ),
                                                i32::from(
                                                    endpoint_atom.endpoint != 0
                                                        || i32::from(eif1.cNeutralBondsValence)
                                                            > i32::from(endpoint_atom.valence),
                                                ),
                                            )
                                        }
                                    }
                                };
                                let possibly_endpoint = donor.wrapping_add(acceptor);
                                donors = donors.wrapping_add(donor);
                                acceptors = acceptors.wrapping_add(acceptor);
                                if possibly_endpoint == 0 {
                                    k = k.wrapping_add(1);
                                    continue;
                                }
                                if rule == Rule::Pt13 {
                                    num_o = num_o.wrapping_add(i32::from(endpoint_valence == 2));
                                    num_c = num_c.wrapping_add(i32::from(endpoint_valence == 4));
                                } else if rule == Rule::Pt18 {
                                    num_o = num_o.wrapping_add(i32::from(endpoint_valence == 2));
                                    num_n = num_n.wrapping_add(i32::from(endpoint_valence == 3));
                                }
                                let endpoint_record = endpoints
                                    .get_mut(endpoint_slot)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                endpoint_record.nGroupNumber = endpoint_atom.endpoint;
                                endpoint_record.nEquNumber = 0;
                                endpoint_record.nAtomNumber = endpoint as AT_NUMB;
                                if group_number != endpoint_atom.endpoint {
                                    different_groups = different_groups.wrapping_add(1);
                                    group_number = endpoint_atom.endpoint;
                                }
                                let bond_slot = usize::try_from(bond_count)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                *bonds
                                    .get_mut(bond_slot)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)? = T_BONDPOS {
                                    nAtomNumber: centerpoint as AT_NUMB,
                                    neighbor_index: k as AT_NUMB,
                                };
                                bond_count = bond_count.wrapping_add(1);
                                possible_mobile = possible_mobile.wrapping_add(i32::from(
                                    mobile > 0
                                        || ((rule != Rule::Pt13 && rule != Rule::Pt18)
                                            && endpoint_atom.endpoint != 0),
                                ));
                                endpoint_count = endpoint_count.wrapping_add(1);
                                k = k.wrapping_add(1);
                            }
                            let special_counts = match rule {
                                Rule::Pt39 => num_o == 1 && num_n < 2,
                                Rule::Pt13 => num_o == 1 && num_c == 1,
                                Rule::Pt18 => num_o == 1 && num_n == 1,
                                _ => true,
                            };
                            if endpoint_count > 1
                                && possible_mobile != 0
                                && donors != 0
                                && acceptors != 0
                                && special_counts
                            {
                                register_and_mark!(
                                    &mut endpoints,
                                    &mut endpoint_count,
                                    &mut bonds,
                                    &mut bond_count,
                                    group_number,
                                    different_groups,
                                    $mode
                                );
                            }
                        }
                        j = j.wrapping_add(1);
                    }
                }
                i = i.wrapping_add(1);
            }
        }};
    }

    scan_rule!(Rule::Standard, ALT_PATH_MODE_TAUTOM as i32);
    if t_group_info.bTautFlags & TG_FLAG_PT_22_00 as u64 != 0 {
        scan_rule!(Rule::Pt22, ALT_PATH_MODE_TAUTOM_PT_22_00 as i32);
    }
    if t_group_info.bTautFlags & TG_FLAG_PT_16_00 as u64 != 0 {
        scan_rule!(Rule::Pt16, ALT_PATH_MODE_TAUTOM_PT_16_00 as i32);
    }
    if t_group_info.bTautFlags & TG_FLAG_PT_06_00 as u64 != 0 {
        scan_rule!(Rule::Pt06, ALT_PATH_MODE_TAUTOM_PT_06_00 as i32);
    }
    if t_group_info.bTautFlags & TG_FLAG_PT_39_00 as u64 != 0 {
        scan_rule!(Rule::Pt39, ALT_PATH_MODE_TAUTOM_PT_39_00 as i32);
    }
    if t_group_info.bTautFlags & TG_FLAG_PT_13_00 as u64 != 0 {
        scan_rule!(Rule::Pt13, ALT_PATH_MODE_TAUTOM_PT_13_00 as i32);
    }
    if t_group_info.bTautFlags & TG_FLAG_PT_18_00 as u64 != 0 {
        scan_rule!(Rule::Pt18, ALT_PATH_MODE_TAUTOM_PT_18_00 as i32);
    }

    if t_group_info.bTautFlags & TG_FLAG_KETO_ENOL_TAUT as u64 != 0 {
        let mut i = 0_i32;
        while i < num_atoms {
            let mut eif1 = ENDPOINT_INFO::default();
            let endpoint_valence = nGetEndpointInfo_KET(heap.slice(at.as_const())?, i, &mut eif1);
            if endpoint_valence != 0 {
                let initial_atom = heap
                    .slice(at.as_const())?
                    .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let mut j = 0_i32;
                while j < i32::from(initial_atom.valence) {
                    let j_index =
                        usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let bond_type =
                        i32::from(initial_atom.bond_type[j_index]) & !(BOND_MARK_ALL as i32);
                    let centerpoint = i32::from(initial_atom.neighbor[j_index]);
                    let center = heap
                        .slice(at.as_const())?
                        .get(
                            usize::try_from(centerpoint)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    if matches!(
                        bond_type,
                        value if value == BOND_DOUBLE as i32
                            || value == BOND_ALTERN as i32
                            || value == BOND_ALT12NS as i32
                            || value == BOND_TAUTOM as i32
                    ) && is_centerpoint_elem_KET(center.el_number) != 0
                        && center.charge == 0
                        && center.radical == 0
                        && i32::from(center.chem_bonds_valence)
                            .wrapping_add(i32::from(center.num_H))
                            == 4
                        && allowed_edge!(i, j)
                    {
                        let mut endpoints =
                            std::array::from_fn::<_, { MAXVAL as usize }, _>(|_| {
                                T_ENDPOINT::default()
                            });
                        let mut bonds = std::array::from_fn::<_, { MAXVAL as usize }, _>(|_| {
                            T_BONDPOS::default()
                        });
                        let mut endpoint_count = 0_i32;
                        let mut bond_count = 0_i32;
                        let mut possible_mobile = 0_i32;
                        let mut group_number = num_atoms as AT_NUMB;
                        let mut different_groups = -1_i32;
                        let mut donors = 0_i32;
                        let mut acceptors = 0_i32;
                        let mut num_o = 0_i32;
                        let mut num_c = 0_i32;
                        let mut k = 0_i32;
                        while k < i32::from(center.valence) {
                            let k_index = usize::try_from(k)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let endpoint = i32::from(center.neighbor[k_index]);
                            let endpoint_atom = heap
                                .slice(at.as_const())?
                                .get(
                                    usize::try_from(endpoint)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .clone();
                            let endpoint_bond_type =
                                i32::from(center.bond_type[k_index]) & !(BOND_MARK_ALL as i32);
                            let alt_bond = endpoint_bond_type == BOND_ALTERN as i32
                                || endpoint_bond_type == BOND_ALT12NS as i32
                                || endpoint_bond_type == BOND_TAUTOM as i32;
                            let non_taut_bond = endpoint_bond_type == BOND_SINGLE as i32
                                || endpoint_bond_type == BOND_DOUBLE as i32;
                            if !allowed_edge!(centerpoint, k) || (!alt_bond && !non_taut_bond) {
                                k = k.wrapping_add(1);
                                continue;
                            }
                            let mut eif2 = ENDPOINT_INFO::default();
                            let endpoint_valence = nGetEndpointInfo_KET(
                                heap.slice(at.as_const())?,
                                endpoint,
                                &mut eif2,
                            );
                            if endpoint_valence == 0 {
                                k = k.wrapping_add(1);
                                continue;
                            }
                            let endpoint_slot = usize::try_from(endpoint_count)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let mobile = AddAtom2num(
                                &mut endpoints
                                    .get_mut(endpoint_slot)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .num,
                                heap.slice(at.as_const())?,
                                endpoint,
                                2,
                            )?;
                            AddAtom2DA(
                                &mut endpoints
                                    .get_mut(endpoint_slot)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .num_DA,
                                heap.slice(at.as_const())?,
                                endpoint,
                                2,
                            )?;
                            let (donor, acceptor) = if non_taut_bond {
                                (
                                    i32::from(
                                        endpoint_bond_type == BOND_SINGLE as i32
                                            && (mobile != 0 || endpoint_atom.endpoint != 0),
                                    ),
                                    i32::from(endpoint_bond_type == BOND_DOUBLE as i32),
                                )
                            } else {
                                (
                                    i32::from(endpoint_atom.endpoint != 0 || eif1.cDonor != 0),
                                    i32::from(
                                        endpoint_atom.endpoint != 0
                                            || i32::from(eif1.cNeutralBondsValence)
                                                > i32::from(endpoint_atom.valence),
                                    ),
                                )
                            };
                            donors = donors.wrapping_add(donor);
                            acceptors = acceptors.wrapping_add(acceptor);
                            if donor.wrapping_add(acceptor) == 0 {
                                k = k.wrapping_add(1);
                                continue;
                            }
                            num_o = num_o.wrapping_add(i32::from(endpoint_valence == 2));
                            num_c = num_c.wrapping_add(i32::from(endpoint_valence == 4));
                            let record = endpoints
                                .get_mut(endpoint_slot)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            record.nGroupNumber = endpoint_atom.endpoint;
                            record.nEquNumber = 0;
                            record.nAtomNumber = endpoint as AT_NUMB;
                            if group_number != endpoint_atom.endpoint {
                                different_groups = different_groups.wrapping_add(1);
                                group_number = endpoint_atom.endpoint;
                            }
                            let bond_slot = usize::try_from(bond_count)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            *bonds
                                .get_mut(bond_slot)
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = T_BONDPOS {
                                nAtomNumber: centerpoint as AT_NUMB,
                                neighbor_index: k as AT_NUMB,
                            };
                            bond_count = bond_count.wrapping_add(1);
                            possible_mobile = possible_mobile
                                .wrapping_add(i32::from(mobile > 0 || endpoint_atom.endpoint != 0));
                            endpoint_count = endpoint_count.wrapping_add(1);
                            k = k.wrapping_add(1);
                        }
                        if endpoint_count > 1
                            && possible_mobile != 0
                            && donors != 0
                            && acceptors != 0
                            && num_o == 1
                            && num_c != 0
                        {
                            register_and_mark!(
                                &mut endpoints,
                                &mut endpoint_count,
                                &mut bonds,
                                &mut bond_count,
                                group_number,
                                different_groups,
                                ALT_PATH_MODE_TAUTOM_KET as i32
                            );
                        }
                    }
                    j = j.wrapping_add(1);
                }
            }
            i = i.wrapping_add(1);
        }
    }

    if tot_changes == 0 {
        let path_positions = match inchi_calloc::<AT_RANK>(heap, num_atoms as u64, 2) {
            Ok(pointer) => pointer,
            Err(
                SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange,
            ) => {
                tot_changes = CT_OUT_OF_RAM;
                SourceMutPointer::null()
            }
            Err(error) => return Err(error),
        };
        if !path_positions.is_null() {
            let ring_result = (|| -> Result<(), SourceHeapError> {
                let mut dfs_path = std::array::from_fn::<_, 8, _>(|_| DFS_PATH::default());
                let mut endpoints =
                    std::array::from_fn::<_, { MAXVAL as usize }, _>(|_| T_ENDPOINT::default());
                let mut bonds =
                    std::array::from_fn::<_, { MAXVAL as usize }, _>(|_| T_BONDPOS::default());

                macro_rules! accept_ring_result {
                    ($ret:expr, $endpoint_count:expr, $bond_count:expr) => {{
                        let ret = $ret;
                        if ret > 0 {
                            if $endpoint_count != 0 {
                                let num_changes = RegisterEndPoints(
                                    heap,
                                    pCG,
                                    t_group_info,
                                    &mut endpoints,
                                    $endpoint_count,
                                    at,
                                    num_atoms,
                                    c_group_info.as_deref_mut(),
                                    Some(&mut *pBNS),
                                )?;
                                if num_changes == -1 {
                                    n_err = CT_TAUCOUNT_ERR;
                                }
                                if num_changes < 0 {
                                    n_err = num_changes;
                                }
                                if n_err != 0 {
                                    return Ok(());
                                }
                                tot_changes = tot_changes.wrapping_add(i32::from(num_changes > 0));
                            }
                            if $bond_count != 0 {
                                let changed =
                                    SetTautomericBonds(heap.slice_mut(at)?, $bond_count, &bonds)?;
                                tot_changes = tot_changes.wrapping_add(i32::from(changed > 0));
                            }
                        } else if is_bns_error(ret) {
                            n_err = ret;
                            return Ok(());
                        }
                    }};
                }

                if t_group_info.bTautFlags & TG_FLAG_1_5_TAUT as u64 != 0 {
                    let mut i1 = 0_i32;
                    while i1 < num_atoms {
                        let mut eif1 = ENDPOINT_INFO::default();
                        if nGetEndpointInfo(heap.slice(at.as_const())?, i1, &mut eif1) != 0 {
                            let mut endpoint_count = 0_i32;
                            let mut bond_count = 0_i32;
                            let ret = nGet15TautInAltPath(
                                heap,
                                pCG,
                                at,
                                i1,
                                path_positions,
                                &mut dfs_path,
                                8,
                                &mut endpoints,
                                MAXVAL as i32,
                                &mut bonds,
                                MAXVAL as i32,
                                &mut endpoint_count,
                                &mut bond_count,
                                pBNS,
                                pBD,
                                num_atoms,
                                clock_result,
                            )?;
                            accept_ring_result!(ret, endpoint_count, bond_count);
                            if n_err != 0 {
                                return Ok(());
                            }
                        }
                        i1 = i1.wrapping_add(1);
                    }
                }

                let mut i1 = 0_i32;
                while i1 < num_atoms {
                    let mut eif1 = ENDPOINT_INFO::default();
                    let atom_i1 = heap
                        .slice(at.as_const())?
                        .get(usize::try_from(i1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    if nGetEndpointInfo(heap.slice(at.as_const())?, i1, &mut eif1) == 3
                        && atom_i1.valence == 2
                        && atom_i1.nNumAtInRingSystem >= 6
                    {
                        let mut endpoint_count = 0_i32;
                        let mut bond_count = 0_i32;
                        let ret = nGet15TautIn6MembAltRing(
                            heap,
                            pCG,
                            at,
                            i1,
                            path_positions,
                            &mut dfs_path,
                            8,
                            &mut endpoints,
                            MAXVAL as i32,
                            &mut bonds,
                            MAXVAL as i32,
                            &mut endpoint_count,
                            &mut bond_count,
                            pBNS,
                            pBD,
                            num_atoms,
                            clock_result,
                        )?;
                        accept_ring_result!(ret, endpoint_count, bond_count);
                        if n_err != 0 {
                            return Ok(());
                        }
                    }
                    i1 = i1.wrapping_add(1);
                }

                i1 = 0;
                while i1 < num_atoms {
                    let atom_i1 = heap
                        .slice(at.as_const())?
                        .get(usize::try_from(i1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let mut eif1 = ENDPOINT_INFO::default();
                    if atom_i1.valence == 2
                        && atom_i1.nNumAtInRingSystem >= 5
                        && nGetEndpointInfo(heap.slice(at.as_const())?, i1, &mut eif1) == 3
                    {
                        let mobile1 =
                            i32::from(atom_i1.num_H).wrapping_add(i32::from(atom_i1.charge == -1));
                        let mut j = 0_i32;
                        while j < i32::from(atom_i1.valence) {
                            let j_index = usize::try_from(j)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let i2 = i32::from(atom_i1.neighbor[j_index]);
                            if i2 < i1 {
                                let atom_i2 = heap
                                    .slice(at.as_const())?
                                    .get(
                                        usize::try_from(i2)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .clone();
                                let bond_type =
                                    i32::from(atom_i1.bond_type[j_index]) & !(BOND_MARK_ALL as i32);
                                let mut eif2 = ENDPOINT_INFO::default();
                                if atom_i2.nRingSystem == atom_i1.nRingSystem
                                    && matches!(
                                        bond_type,
                                        value if value == BOND_SINGLE as i32
                                            || value == BOND_TAUTOM as i32
                                            || value == BOND_ALT12NS as i32
                                            || value == BOND_ALTERN as i32
                                    )
                                    && atom_i2.valence == 2
                                    && nGetEndpointInfo(heap.slice(at.as_const())?, i2, &mut eif2)
                                        == 3
                                {
                                    let mobile2 = i32::from(atom_i2.num_H)
                                        .wrapping_add(i32::from(atom_i2.charge == -1));
                                    if (atom_i1.endpoint != 0 || atom_i2.endpoint != 0)
                                        || mobile1.wrapping_add(mobile2) == 1
                                    {
                                        let mut endpoint_count = 0_i32;
                                        let mut bond_count = 0_i32;
                                        let ret = nGet12TautIn5MembAltRing(
                                            heap,
                                            pCG,
                                            at,
                                            i1,
                                            j,
                                            path_positions,
                                            &mut dfs_path,
                                            8,
                                            &mut endpoints,
                                            MAXVAL as i32,
                                            &mut bonds,
                                            MAXVAL as i32,
                                            &mut endpoint_count,
                                            &mut bond_count,
                                            pBNS,
                                            pBD,
                                            num_atoms,
                                            clock_result,
                                        )?;
                                        accept_ring_result!(ret, endpoint_count, bond_count);
                                        if n_err != 0 {
                                            return Ok(());
                                        }
                                    }
                                }
                            }
                            j = j.wrapping_add(1);
                        }
                    }
                    i1 = i1.wrapping_add(1);
                }

                i1 = 0;
                while i1 < num_atoms {
                    let atom_i1 = heap
                        .slice(at.as_const())?
                        .get(usize::try_from(i1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    if atom_i1.nNumAtInRingSystem >= 5
                        && bIsCenterPointStrict(heap, at, i1)? != 0
                        && atom_i1.bCutVertex != 0
                        && atom_i1.valence == 3
                        && atom_i1.endpoint == 0
                    {
                        let mut j = 0_i32;
                        while j < i32::from(atom_i1.valence) {
                            let j_index = usize::try_from(j)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let i2 = i32::from(atom_i1.neighbor[j_index]);
                            let atom_i2 = heap
                                .slice(at.as_const())?
                                .get(
                                    usize::try_from(i2)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .clone();
                            let qr_bond =
                                i32::from(atom_i1.bond_type[j_index]) & !(BOND_MARK_ALL as i32);
                            if atom_i2.nRingSystem == atom_i1.nRingSystem
                                && bIsCenterPointStrict(heap, at, i2)? != 0
                                && atom_i2.bCutVertex != 0
                                && atom_i2.valence == 3
                                && atom_i2.endpoint == 0
                                && matches!(
                                    qr_bond,
                                    value if value == BOND_SINGLE as i32
                                        || value == BOND_TAUTOM as i32
                                        || value == BOND_ALT12NS as i32
                                        || value == BOND_ALTERN as i32
                                )
                            {
                                let mut k = 0_i32;
                                while k < i32::from(atom_i1.valence) {
                                    let k_index = usize::try_from(k)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                    let endpoint1 = i32::from(atom_i1.neighbor[k_index]);
                                    if endpoint1 != i2 {
                                        let endpoint1_atom =
                                            heap.slice(at.as_const())?
                                                .get(usize::try_from(endpoint1).map_err(|_| {
                                                    SourceHeapError::PointerOutOfBounds
                                                })?)
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                                .clone();
                                        let mut eif1 = ENDPOINT_INFO::default();
                                        let endpoint1_valence = nGetEndpointInfo(
                                            heap.slice(at.as_const())?,
                                            endpoint1,
                                            &mut eif1,
                                        );
                                        let mobile1 = i32::from(endpoint1_atom.num_H)
                                            .wrapping_add(i32::from(endpoint1_atom.charge == -1));
                                        let bond1 = i32::from(atom_i1.bond_type[k_index])
                                            & !(BOND_MARK_ALL as i32);
                                        if endpoint1_valence != 0
                                            && endpoint1_atom.nRingSystem != atom_i1.nRingSystem
                                            && mobile1.wrapping_add(i32::from(
                                                endpoint1_atom.chem_bonds_valence,
                                            )) == endpoint1_valence
                                            && matches!(
                                                bond1,
                                                value if value == BOND_SINGLE as i32
                                                    || value == BOND_DOUBLE as i32
                                                    || value == BOND_TAUTOM as i32
                                                    || value == BOND_ALT12NS as i32
                                                    || value == BOND_ALTERN as i32
                                            )
                                        {
                                            let mut m = 0_i32;
                                            while m < i32::from(atom_i2.valence) {
                                                let m_index = usize::try_from(m).map_err(|_| {
                                                    SourceHeapError::PointerOutOfBounds
                                                })?;
                                                let endpoint2 =
                                                    i32::from(atom_i2.neighbor[m_index]);
                                                if endpoint2 != i1 {
                                                    let endpoint2_atom = heap
                                                        .slice(at.as_const())?
                                                        .get(usize::try_from(endpoint2).map_err(
                                                            |_| SourceHeapError::PointerOutOfBounds,
                                                        )?)
                                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                                        .clone();
                                                    let mut eif2 = ENDPOINT_INFO::default();
                                                    let endpoint2_valence = nGetEndpointInfo(
                                                        heap.slice(at.as_const())?,
                                                        endpoint2,
                                                        &mut eif2,
                                                    );
                                                    let mobile2 = i32::from(endpoint2_atom.num_H)
                                                        .wrapping_add(i32::from(
                                                            endpoint2_atom.charge == -1,
                                                        ));
                                                    let bond2 =
                                                        i32::from(atom_i2.bond_type[m_index])
                                                            & !(BOND_MARK_ALL as i32);
                                                    if endpoint2_valence != 0
                                                        && endpoint2_atom.nRingSystem
                                                            != atom_i2.nRingSystem
                                                        && matches!(
                                                            bond2,
                                                            value if value == BOND_SINGLE as i32
                                                                || value == BOND_DOUBLE as i32
                                                                || value == BOND_TAUTOM as i32
                                                                || value == BOND_ALT12NS as i32
                                                                || value == BOND_ALTERN as i32
                                                        )
                                                    {
                                                        let can_move_1_to_2 = allowed_edge!(i1, k)
                                                            && allowed_edge!(i2, m)
                                                            && (endpoint1_atom.endpoint != 0
                                                                || mobile1 != 0)
                                                            && bond1 != BOND_DOUBLE as i32
                                                            && (endpoint2_atom.endpoint != 0
                                                                || i32::from(
                                                                    eif2.cNeutralBondsValence,
                                                                ) > i32::from(
                                                                    endpoint2_atom.valence,
                                                                ))
                                                            && bond2 != BOND_SINGLE as i32;
                                                        let can_move_2_to_1 = allowed_edge!(i1, k)
                                                            && allowed_edge!(i2, m)
                                                            && (endpoint2_atom.endpoint != 0
                                                                || mobile2 != 0)
                                                            && bond2 != BOND_DOUBLE as i32
                                                            && (endpoint1_atom.endpoint != 0
                                                                || i32::from(
                                                                    eif1.cNeutralBondsValence,
                                                                ) > i32::from(
                                                                    endpoint1_atom.valence,
                                                                ))
                                                            && bond1 != BOND_SINGLE as i32;
                                                        if (can_move_1_to_2 || can_move_2_to_1)
                                                            && !(bond1 == bond2
                                                                && matches!(
                                                                    bond1,
                                                                    value if value
                                                                        == BOND_SINGLE as i32
                                                                        || value
                                                                            == BOND_DOUBLE as i32
                                                                ))
                                                        {
                                                            if atom_i1.nNumAtInRingSystem >= 7 {
                                                                let mut endpoint_count = 0_i32;
                                                                let mut bond_count = 0_i32;
                                                                let ret = nGet14TautIn7MembAltRing(
                                                                    heap,
                                                                    pCG,
                                                                    at,
                                                                    i1,
                                                                    j,
                                                                    k,
                                                                    m,
                                                                    path_positions,
                                                                    &mut dfs_path,
                                                                    8,
                                                                    &mut endpoints,
                                                                    MAXVAL as i32,
                                                                    &mut bonds,
                                                                    MAXVAL as i32,
                                                                    &mut endpoint_count,
                                                                    &mut bond_count,
                                                                    pBNS,
                                                                    pBD,
                                                                    num_atoms,
                                                                    clock_result,
                                                                )?;
                                                                accept_ring_result!(
                                                                    ret,
                                                                    endpoint_count,
                                                                    bond_count
                                                                );
                                                                if n_err != 0 {
                                                                    return Ok(());
                                                                }
                                                            }
                                                            if atom_i1.nNumAtInRingSystem >= 5 {
                                                                let mut endpoint_count = 0_i32;
                                                                let mut bond_count = 0_i32;
                                                                let ret = nGet14TautIn5MembAltRing(
                                                                    heap,
                                                                    pCG,
                                                                    at,
                                                                    i1,
                                                                    j,
                                                                    k,
                                                                    m,
                                                                    path_positions,
                                                                    &mut dfs_path,
                                                                    8,
                                                                    &mut endpoints,
                                                                    MAXVAL as i32,
                                                                    &mut bonds,
                                                                    MAXVAL as i32,
                                                                    &mut endpoint_count,
                                                                    &mut bond_count,
                                                                    pBNS,
                                                                    pBD,
                                                                    num_atoms,
                                                                    clock_result,
                                                                )?;
                                                                accept_ring_result!(
                                                                    ret,
                                                                    endpoint_count,
                                                                    bond_count
                                                                );
                                                                if n_err != 0 {
                                                                    return Ok(());
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                m = m.wrapping_add(1);
                                            }
                                        }
                                    }
                                    k = k.wrapping_add(1);
                                }
                            }
                            j = j.wrapping_add(1);
                        }
                    }
                    i1 = i1.wrapping_add(1);
                }
                Ok(())
            })();
            let free_result = inchi_free(heap, path_positions);
            ring_result?;
            free_result?;
        }
    }

    Ok(if n_err < 0 { n_err } else { tot_changes })
}

#[allow(non_snake_case)]
pub(crate) fn CompRankTautomer(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pn_tRankForSort: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:7008 CompRankTautomer
    // INCHI✔️✔️: complete active source frame follows verbatim.
    /*
    int CompRankTautomer( const void* a1, const void* a2, void *p )
    {
        AT_RANK *pn_tRankForSort = (AT_RANK *) p;
        int ret = (int) pn_tRankForSort[(int) ( *(const AT_RANK*) a1 )] -
            (int) pn_tRankForSort[(int) ( *(const AT_RANK*) a2 )];

        return ret;
    }
    */
    // END INCHI C FUNCTION: CompRankTautomer

    let ranks = heap.slice(pn_tRankForSort.as_const())?;
    let first = *ranks
        .get(usize::from(a1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = *ranks
        .get(usize::from(a2))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(first) - i32::from(second))
}

#[allow(non_snake_case)]
pub(crate) fn SortTautomerGroupsAndEndpoints(
    heap: &mut SourceHeap,
    _pCG: &mut CANON_GLOBALS,
    t_group_info: &T_GROUP_INFO,
    num_atoms: i32,
    num_at_tg: i32,
    nRank: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:7019 SortTautomerGroupsAndEndpoints
    // INCHI✔️❌: complete active source frame follows verbatim.
    /*
    int SortTautomerGroupsAndEndpoints( CANON_GLOBALS *pCG,
                                        T_GROUP_INFO *t_group_info,
                                        int num_atoms,
                                        int num_at_tg,
                                        AT_RANK *nRank )
    {
        int i, nFirstEndpointAtNoPos, nNumEndpoints;
        AT_RANK *pn_tRankForSort;
        AT_NUMB  *nEndpointAtomNumber;
        int       num_t_groups = num_at_tg - num_atoms;
        T_GROUP   *t_group = NULL;

        /*  Check if sorting is required */
        if (num_t_groups <= 0 || t_group_info->nNumEndpoints < 2)
        {
            return 0; /*  no tautomer data */
        }

        t_group = t_group_info->t_group;

        /*  Sort endpoints within the groups */
        for (i = 0; i < num_t_groups; i++)
        {
            if (t_group[i].nNumEndpoints < 2)
            {
                continue;  /*  program error; should not happen */ /*   <BRKPT> */
            }

            /*  Set globals for sorting */
            nFirstEndpointAtNoPos = t_group[i].nFirstEndpointAtNoPos;
            nNumEndpoints = t_group[i].nNumEndpoints;
            if (nNumEndpoints + nFirstEndpointAtNoPos > t_group_info->nNumEndpoints)
            {
                /*  for debug only */
                return CT_TAUCOUNT_ERR; /*  program error */ /*   <BRKPT> */
            }
            nEndpointAtomNumber = t_group_info->nEndpointAtomNumber + (int) nFirstEndpointAtNoPos;
            pn_tRankForSort = nRank;

            insertions_sort( pn_tRankForSort, nEndpointAtomNumber, nNumEndpoints,
                             sizeof( nEndpointAtomNumber[0] ), CompRankTautomer );
        }

        /*  Sort the tautomeric groups according to their ranks only
        (that is, ignoring the isotopic composition of the mobile groups and ranks of the endpoints) */
        if (t_group_info->num_t_groups > 1)
        {
            /*  Set globals for sorting */
            /*  a hack: the ranks of all tautomeric groups are */
            /*  located at nRank[num_atoms..num_at_tg-1] */
            pn_tRankForSort = nRank + num_atoms;
            /*  Sort */
            /*  ordering numbers to sort : t_group_info->tGroupNumber; */

            insertions_sort( pn_tRankForSort, t_group_info->tGroupNumber,
                             num_t_groups, sizeof( t_group_info->tGroupNumber[0] ),
                             CompRankTautomer );
        }

        return t_group_info->num_t_groups;
    }
    */
    // END INCHI C FUNCTION: SortTautomerGroupsAndEndpoints

    let num_t_groups = num_at_tg.wrapping_sub(num_atoms);
    if num_t_groups <= 0 || t_group_info.nNumEndpoints < 2 {
        return Ok(0);
    }

    let group_count = usize::try_from(num_t_groups)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let groups = heap
        .slice(t_group_info.t_group.as_const())?
        .get(..group_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();

    for group in groups {
        if group.nNumEndpoints < 2 {
            continue;
        }

        let first_endpoint = i32::from(group.nFirstEndpointAtNoPos);
        let endpoint_count = i32::from(group.nNumEndpoints);
        if endpoint_count + first_endpoint > t_group_info.nNumEndpoints {
            return Ok(CT_TAUCOUNT_ERR);
        }
        let endpoint_pointer = t_group_info
            .nEndpointAtomNumber
            .offset(i64::from(first_endpoint))?;
        heap.with_slice_mut_and_heap(endpoint_pointer, |endpoints, heap| {
            insertions_sort(
                bytemuck::cast_slice_mut::<AT_NUMB, u8>(endpoints),
                usize::try_from(endpoint_count)
                    .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
                std::mem::size_of::<AT_NUMB>(),
                &mut |left, right| {
                    CompRankTautomer(
                        heap,
                        AT_NUMB::from_ne_bytes([left[0], left[1]]),
                        AT_NUMB::from_ne_bytes([right[0], right[1]]),
                        nRank,
                    )
                },
            )
            .map(|_| ())
        })?;
    }

    if t_group_info.num_t_groups > 1 {
        let rank_pointer = nRank.offset(i64::from(num_atoms))?;
        heap.with_slice_mut_and_heap(t_group_info.tGroupNumber, |group_numbers, heap| {
            insertions_sort(
                bytemuck::cast_slice_mut::<AT_NUMB, u8>(group_numbers),
                group_count,
                std::mem::size_of::<AT_NUMB>(),
                &mut |left, right| {
                    CompRankTautomer(
                        heap,
                        AT_NUMB::from_ne_bytes([left[0], left[1]]),
                        AT_NUMB::from_ne_bytes([right[0], right[1]]),
                        rank_pointer,
                    )
                },
            )
            .map(|_| ())
        })?;
    }

    Ok(t_group_info.num_t_groups)
}

pub(crate) fn free_t_group_info(
    heap: &mut SourceHeap,
    t_group_info: Option<&mut T_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6336 free_t_group_info
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap frees add overhead.
    /*
    int free_t_group_info( T_GROUP_INFO *t_group_info )
    {
        if (t_group_info)
        {
            if (t_group_info->t_group)
            {
                inchi_free( t_group_info->t_group );
            }
            if (t_group_info->nEndpointAtomNumber)
            {
                inchi_free( t_group_info->nEndpointAtomNumber );
            }
            if (t_group_info->tGroupNumber)
            {
                inchi_free( t_group_info->tGroupNumber );
            }
            if (t_group_info->nIsotopicEndpointAtomNumber)
            {
                inchi_free( t_group_info->nIsotopicEndpointAtomNumber );
            }
            memset( t_group_info, 0, sizeof( *t_group_info ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }

        return 0;
    }
        */
    // END INCHI C FUNCTION: free_t_group_info

    if let Some(t_group_info) = t_group_info {
        if !t_group_info.t_group.is_null() {
            inchi_free(heap, t_group_info.t_group)?;
        }
        if !t_group_info.nEndpointAtomNumber.is_null() {
            inchi_free(heap, t_group_info.nEndpointAtomNumber)?;
        }
        if !t_group_info.tGroupNumber.is_null() {
            inchi_free(heap, t_group_info.tGroupNumber)?;
        }
        if !t_group_info.nIsotopicEndpointAtomNumber.is_null() {
            inchi_free(heap, t_group_info.nIsotopicEndpointAtomNumber)?;
        }
        *t_group_info = T_GROUP_INFO::default();
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CountTautomerGroups(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    t_group_info: Option<&mut T_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6519 CountTautomerGroups
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int CountTautomerGroups( sp_ATOM *at,
                             int num_atoms,
                             T_GROUP_INFO *t_group_info )
    {
        int i, j, ret = 0, nNumEndpoints, max_t_group, num_groups_noH;

        AT_NUMB    nGroupNumber, nNewGroupNumber, *nCurrEndpointAtNoPos = NULL;

        T_GROUP   *t_group;
        int        num_t;
        /* int bIgnoreIsotopic, max_num_t; */
        AT_NUMB   *nTautomerGroupNumber = NULL;
        AT_NUMB   *nEndpointAtomNumber = NULL;
        AT_NUMB   *tGroupNumber = NULL;

        if (!t_group_info || !t_group_info->t_group || 0 >= t_group_info->max_num_t_groups)
        {
            return 0; /* empty t-groups */
        }
        num_t = t_group_info->num_t_groups;
        t_group = t_group_info->t_group;
        /*
        max_num_t       =  t_group_info->max_num_t_groups;
        bIgnoreIsotopic =  t_group_info->bIgnoreIsotopic;
        */
        num_groups_noH = 0;

        /* the following 2 arrays are to be rebuilt here */
        /* djb-rwth: fixing oss-fuzz issue #70006 */
        if (t_group_info->nEndpointAtomNumber)
        {
            /* free( t_group_info->nEndpointAtomNumber ); */
            t_group_info->nEndpointAtomNumber = NULL;
        }
        if (t_group_info->tGroupNumber)
        {
            /* inchi_free( t_group_info->tGroupNumber ); */
            t_group_info->tGroupNumber = NULL;
        }

        /*  Find max_t_group */
        for (i = 0, max_t_group = 0; i < t_group_info->num_t_groups; i++)
        {
            if (max_t_group < t_group[i].nGroupNumber)
                max_t_group = t_group[i].nGroupNumber;
        }
        /*  Allocate memory for temp storage of numbers of endpoints  */
        if (max_t_group &&
             !( nTautomerGroupNumber = (AT_NUMB*) inchi_calloc( (long long)max_t_group + 1, sizeof( nTautomerGroupNumber[0] ) ) /*temp*/ )) /* djb-rwth: cast operator added */
        {
            goto err_exit_function; /*  program error: out of RAM */ /*   <BRKPT> */
        }

        /*  Count endpoints for each tautomer group */
        for (i = 0, nNumEndpoints = 0; i < num_atoms; i++)
        {
            if (( j = at[i].endpoint ) == 0)
            {
                continue;
            }
            if (j > max_t_group) /*  debug only */
            {
                goto err_exit_function; /*  program error */ /*   <BRKPT> */
            }
            if (nTautomerGroupNumber) /* djb-rwth: fixing a NULL pointer dereference */
                nTautomerGroupNumber[j] ++;
            nNumEndpoints++;
        }

        if (!nNumEndpoints)
        {
            goto exit_function; /*  not a tautomer */
        }

        /*  allocate temporary array */
        if (!( nEndpointAtomNumber = (AT_NUMB*) inchi_calloc( nNumEndpoints, sizeof( nEndpointAtomNumber[0] ) ) ) ||
             !( nCurrEndpointAtNoPos = (AT_NUMB*) inchi_calloc( num_t, sizeof( nCurrEndpointAtNoPos[0] ) ) /*temp*/ ))
        {
            goto err_exit_function; /*   program error: out of RAM */ /*   <BRKPT> */
        }

        /*
        * Remove missing endpoints from t_group. Since only one
        * disconnected part is processed, some endpoints groups may have disappeared.
        * Mark t_groups containing charges only for subsequent removal
        */
        for (i = 0, nNewGroupNumber = 0; i < num_t; /*i ++*/)
        {
            int bNoH = 0, nNumH;
            nGroupNumber = t_group[i].nGroupNumber;
            for (j = 1, nNumH = t_group[i].num[0]; j < T_NUM_NO_ISOTOPIC; j++)
            {
                nNumH -= (int) t_group[i].num[j];
            }
            if (nTautomerGroupNumber) /* djb-rwth: fixing a NULL pointer dereference */
            {
                if (t_group[i].nNumEndpoints != nTautomerGroupNumber[(int)nGroupNumber]
    #if ( IGNORE_TGROUP_WITHOUT_H == 1 )
                    || (bNoH = (t_group[i].num[0] == t_group[i].num[1]))  /* only for (H,-) t-groups; (+) t-groups are not removed */
    #endif
                    ) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    if (!nTautomerGroupNumber[(int)nGroupNumber] || bNoH)
                    {
                        /*  The group belongs to another disconnected part of the structure or has only charges */
                        /*  Remove the group */
                        num_t--;
                        if (i < num_t)
                            memmove(t_group + i, t_group + i + 1, ((long long)num_t - (long long)i) * sizeof(t_group[0])); /* djb-rwth: cast operators added */
                        if (bNoH)
                        {
                            /*  Group contains no mobile hydrogen atoms, only charges. Prepare to remove it. */
                            nTautomerGroupNumber[(int)nGroupNumber] = 0;
                            num_groups_noH++;
                        }
                        /*i --;*/
                    }
                    else
                    {
                        /*  Different number of endpoints */
                        goto err_exit_function; /*  program error */ /*   <BRKPT> */
                    }
                }
                else
                {
                    /*  Renumber t_group and prepare to renumber at[i].endpoint */
                    nTautomerGroupNumber[(int)nGroupNumber] =
                        t_group[i].nGroupNumber = ++nNewGroupNumber; /*  = i+1 */
                    /*  get first group atom orig. number position in the nEndpointAtomNumber[] */
                    /*  and in the tautomer endpoint canon numbers part of the connection table */
                    t_group[i].nFirstEndpointAtNoPos = nCurrEndpointAtNoPos[i] =
                        i ? (t_group[i - 1].nFirstEndpointAtNoPos + t_group[i - 1].nNumEndpoints) : 0;
                    t_group[i].num[0] = nNumH;
    #if ( REMOVE_TGROUP_CHARGE == 1 )
                    t_group[i].num[1] = 0;  /* remove only (-) charges */
    #endif
                    /* -- wrong condition. Disabled.
                    if ( t_group[i].nGroupNumber != i + 1 ) { // for debug only
                    goto err_exit_function; // program error
                    }
                    */
                    i++;
                }
            }
        }
        if (num_t != nNewGroupNumber)
        {
            /*  for debug only */
            goto err_exit_function; /*  program error */ /*   <BRKPT> */
        }

        /*  Check if any tautomer group was left */
        if (!nNewGroupNumber)
        {
            if (!num_groups_noH)
            {
                goto err_exit_function; /*  program error: not a tautomer */ /*   <BRKPT> */
            }
            else
            {
                goto exit_function;
            }
        }

        /*
        * An array for tautomer group sorting later, at the time of storing Connection Table
        * Later the sorting consists out of 2 steps:
        * 1) Sort t_group[i].nNumEndpoints endpoint atom ranks within each endpoint group
        *    starting from t_group[i].nFirstEndpointAtNoPos; i = 0..t_group_info->num_t_groups-1
        * 2) Sort the groups indexes t_group_info->tGroupNumber[]
        */
        if (!( tGroupNumber =
            (AT_NUMB*) inchi_calloc( (long long)nNewGroupNumber * (long long)TGSO_TOTAL_LEN, sizeof( tGroupNumber[0] ) ) )) /* djb-rwth: cast operator added */
        {
            goto err_exit_function; /*  out of RAM */
        }
        for (i = 0; i < nNewGroupNumber; i++)
        {
            tGroupNumber[i] = (AT_NUMB) i; /*  initialization: original t_group number = (at[i]->endpoint-1) */
        }
        /*
        * Renumber endpoint atoms and save their orig. atom
        * numbers for filling out the tautomer part of the LinearCT.
        * nCurrEndpointAtNoPos[j] is an index of the atom number in the nEndpointAtomNumber[]
        */
        for (i = 0; i < num_atoms; i++)
        {
            if ((j = (int) at[i].endpoint)) /* djb-rwth: addressing LLVM warning */
            {
                j = (int) ( at[i].endpoint = nTautomerGroupNumber[j] ) - 1; /*  new t_group number */
                if (j >= 0 && j < num_t) /* djb-rwth: fixing buffer overflow */
                {
                    /*  j=-1 in case of no mobile hydrogen atoms (charges only), group being removed */
                    if (nCurrEndpointAtNoPos[j] >=   /*  debug only */
                         t_group[j].nFirstEndpointAtNoPos + t_group[j].nNumEndpoints)
                    {
                        goto err_exit_function; /*  program error */ /*   <BRKPT> */
                    }
                    nEndpointAtomNumber[(int) nCurrEndpointAtNoPos[j] ++] = (AT_NUMB) i;
                }
                else
                {
                    nNumEndpoints--; /*  endpoint has been removed */
                }
            }
        }
        t_group_info->num_t_groups = nNewGroupNumber;
        t_group_info->nNumEndpoints = nNumEndpoints;
        t_group_info->nEndpointAtomNumber = nEndpointAtomNumber;
        t_group_info->tGroupNumber = tGroupNumber; /* only the 1st segment filled */
        inchi_free( nTautomerGroupNumber );
        inchi_free( nCurrEndpointAtNoPos );
        return nNumEndpoints + T_GROUP_HDR_LEN * nNewGroupNumber + 1; /*  nLenLinearCTTautomer */

    err_exit_function:
        ret = CT_TAUCOUNT_ERR;

    exit_function:

        /*  release allocated memory; set "no tautomeric group" */
        if (nEndpointAtomNumber)
        {
            inchi_free( nEndpointAtomNumber );
        }
        if (nTautomerGroupNumber)
        {
            inchi_free( nTautomerGroupNumber );
        }
        if (tGroupNumber)
        {
            inchi_free( tGroupNumber );
        }
        if (nCurrEndpointAtNoPos)
        {
            inchi_free( nCurrEndpointAtNoPos );
        }
        t_group_info->nNumEndpoints = 0;
        t_group_info->num_t_groups = 0;
        if (!ret && ( ( t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT ) ||
                      (t_group_info->nNumIsotopicEndpoints > 1 && ( t_group_info->bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ) )) )) /* djb-rwth: addressing LLVM warning */
        {
            ret = 1; /* only protons have been (re)moved or neitralization happened */
        }

        return ret;
    }
        */
    // END INCHI C FUNCTION: CountTautomerGroups
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CountTautomerGroups
    // INCHI✔️❌: IGNORE_TGROUP_WITHOUT_H == 1 and REMOVE_TGROUP_CHARGE == 0.
    // INCHI✔️❌: TGSO_TOTAL_LEN == 4 and T_GROUP_HDR_LEN == 3.
    // END INCHI ACTIVE MACRO CONFIGURATION: CountTautomerGroups

    let Some(info) = t_group_info else {
        return Ok(0);
    };
    if info.t_group.is_null() || info.max_num_t_groups <= 0 {
        return Ok(0);
    }

    let mut number_of_groups = info.num_t_groups;
    let group_count =
        usize::try_from(number_of_groups).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if heap.slice(info.t_group.as_const())?.len() < group_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    // The active 2024 fix intentionally orphans these old allocations.
    info.nEndpointAtomNumber = SourceMutPointer::null();
    info.tGroupNumber = SourceMutPointer::null();

    let mut maximum_group = 0_i32;
    for group in &heap.slice(info.t_group.as_const())?[..group_count] {
        maximum_group = maximum_group.max(i32::from(group.nGroupNumber));
    }

    let mut group_number_map = SourceMutPointer::<AT_NUMB>::null();
    let mut endpoint_atom_numbers = SourceMutPointer::<AT_NUMB>::null();
    let mut current_endpoint_positions = SourceMutPointer::<AT_NUMB>::null();
    let mut sorted_group_numbers = SourceMutPointer::<AT_NUMB>::null();
    let mut return_code = 0_i32;

    let process = (|| -> Result<Option<i32>, SourceHeapError> {
        if maximum_group != 0 {
            group_number_map = match inchi_calloc::<AT_NUMB>(
                heap,
                u64::try_from(maximum_group)
                    .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
                    + 1,
                std::mem::size_of::<AT_NUMB>() as u64,
            ) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
            if group_number_map.is_null() {
                return_code = CT_TAUCOUNT_ERR;
                return Ok(None);
            }
        }

        let atom_count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if heap.slice(at.as_const())?.len() < atom_count {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let mut number_of_endpoints = 0_i32;
        for index in 0..atom_count {
            let endpoint = i32::from(heap.slice(at.as_const())?[index].endpoint);
            if endpoint == 0 {
                continue;
            }
            if endpoint > maximum_group {
                return_code = CT_TAUCOUNT_ERR;
                return Ok(None);
            }
            if !group_number_map.is_null() {
                let endpoint =
                    usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let value = heap.slice(group_number_map.as_const())?[endpoint];
                heap.slice_mut(group_number_map)?[endpoint] = value.wrapping_add(1);
            }
            number_of_endpoints = number_of_endpoints.wrapping_add(1);
        }

        if number_of_endpoints == 0 {
            return Ok(None);
        }

        endpoint_atom_numbers = match inchi_calloc::<AT_NUMB>(
            heap,
            u64::try_from(number_of_endpoints)
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
            std::mem::size_of::<AT_NUMB>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if endpoint_atom_numbers.is_null() {
            return_code = CT_TAUCOUNT_ERR;
            return Ok(None);
        }
        current_endpoint_positions = match inchi_calloc::<AT_NUMB>(
            heap,
            number_of_groups as u64,
            std::mem::size_of::<AT_NUMB>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if current_endpoint_positions.is_null() {
            return_code = CT_TAUCOUNT_ERR;
            return Ok(None);
        }

        let mut groups_without_hydrogen = 0_i32;
        let mut new_group_number = 0_u16;
        let mut index = 0_i32;
        while index < number_of_groups {
            let group_index =
                usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let group = heap.slice(info.t_group.as_const())?[group_index].clone();
            let old_group_number = usize::from(group.nGroupNumber);
            let number_of_hydrogens = i32::from(group.num[0]).wrapping_sub(i32::from(group.num[1]));
            let no_hydrogen = group.num[0] == group.num[1];
            let actual_endpoints = heap.slice(group_number_map.as_const())?[old_group_number];

            if group.nNumEndpoints != actual_endpoints || no_hydrogen {
                if actual_endpoints == 0 || no_hydrogen {
                    number_of_groups -= 1;
                    if index < number_of_groups {
                        let groups = heap.slice_mut(info.t_group)?;
                        for moved in group_index
                            ..usize::try_from(number_of_groups)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                        {
                            groups[moved] = groups[moved + 1].clone();
                        }
                    }
                    if no_hydrogen {
                        heap.slice_mut(group_number_map)?[old_group_number] = 0;
                        groups_without_hydrogen = groups_without_hydrogen.wrapping_add(1);
                    }
                } else {
                    return_code = CT_TAUCOUNT_ERR;
                    return Ok(None);
                }
            } else {
                new_group_number = new_group_number.wrapping_add(1);
                heap.slice_mut(group_number_map)?[old_group_number] = new_group_number;
                let first_position = if index != 0 {
                    let previous = &heap.slice(info.t_group.as_const())?[group_index - 1];
                    previous
                        .nFirstEndpointAtNoPos
                        .wrapping_add(previous.nNumEndpoints)
                } else {
                    0
                };
                {
                    let groups = heap.slice_mut(info.t_group)?;
                    groups[group_index].nGroupNumber = new_group_number;
                    groups[group_index].nFirstEndpointAtNoPos = first_position;
                    groups[group_index].num[0] = number_of_hydrogens as AT_RANK;
                }
                heap.slice_mut(current_endpoint_positions)?[group_index] = first_position;
                index += 1;
            }
        }

        if number_of_groups != i32::from(new_group_number) {
            return_code = CT_TAUCOUNT_ERR;
            return Ok(None);
        }
        if new_group_number == 0 {
            if groups_without_hydrogen == 0 {
                return_code = CT_TAUCOUNT_ERR;
            }
            return Ok(None);
        }

        let sorted_length = u64::from(new_group_number)
            .checked_mul(u64::from(crate::source_types::TGSO_TOTAL_LEN))
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        sorted_group_numbers = match inchi_calloc::<AT_NUMB>(
            heap,
            sorted_length,
            std::mem::size_of::<AT_NUMB>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if sorted_group_numbers.is_null() {
            return_code = CT_TAUCOUNT_ERR;
            return Ok(None);
        }
        for index in 0..usize::from(new_group_number) {
            heap.slice_mut(sorted_group_numbers)?[index] = index as AT_NUMB;
        }

        for atom_index in 0..atom_count {
            let old_endpoint = usize::from(heap.slice(at.as_const())?[atom_index].endpoint);
            if old_endpoint == 0 {
                continue;
            }
            let mapped = heap.slice(group_number_map.as_const())?[old_endpoint];
            heap.slice_mut(at)?[atom_index].endpoint = mapped;
            let group_index = i32::from(mapped) - 1;
            if group_index >= 0 && group_index < number_of_groups {
                let group_index = usize::try_from(group_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let position =
                    usize::from(heap.slice(current_endpoint_positions.as_const())?[group_index]);
                let group = &heap.slice(info.t_group.as_const())?[group_index];
                if position
                    >= usize::from(
                        group
                            .nFirstEndpointAtNoPos
                            .wrapping_add(group.nNumEndpoints),
                    )
                {
                    return_code = CT_TAUCOUNT_ERR;
                    return Ok(None);
                }
                heap.slice_mut(endpoint_atom_numbers)?[position] = atom_index as AT_NUMB;
                heap.slice_mut(current_endpoint_positions)?[group_index] =
                    (position as AT_NUMB).wrapping_add(1);
            } else {
                number_of_endpoints = number_of_endpoints.wrapping_sub(1);
            }
        }

        info.num_t_groups = i32::from(new_group_number);
        info.nNumEndpoints = number_of_endpoints;
        info.nEndpointAtomNumber = endpoint_atom_numbers;
        info.tGroupNumber = sorted_group_numbers;
        endpoint_atom_numbers = SourceMutPointer::null();
        sorted_group_numbers = SourceMutPointer::null();
        Ok(Some(
            number_of_endpoints
                .wrapping_add((T_GROUP_HDR_LEN as i32).wrapping_mul(i32::from(new_group_number)))
                .wrapping_add(1),
        ))
    })();

    let process = match process {
        Ok(value) => value,
        Err(error) => {
            inchi_free(heap, endpoint_atom_numbers)?;
            inchi_free(heap, group_number_map)?;
            inchi_free(heap, sorted_group_numbers)?;
            inchi_free(heap, current_endpoint_positions)?;
            return Err(error);
        }
    };
    inchi_free(heap, endpoint_atom_numbers)?;
    inchi_free(heap, group_number_map)?;
    inchi_free(heap, sorted_group_numbers)?;
    inchi_free(heap, current_endpoint_positions)?;

    if let Some(success) = process {
        return Ok(success);
    }

    info.nNumEndpoints = 0;
    info.num_t_groups = 0;
    if return_code == 0
        && ((info.tni.bNormalizationFlags & u64::from(FLAG_NORM_CONSIDER_TAUT)) != 0
            || (info.nNumIsotopicEndpoints > 1
                && (info.bTautFlagsDone
                    & u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))
                    != 0))
    {
        return_code = 1;
    }
    Ok(return_code)
}

pub(crate) fn set_tautomer_iso_sort_keys(
    heap: &mut SourceHeap,
    t_group_info: Option<&mut T_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6477 set_tautomer_iso_sort_keys
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int set_tautomer_iso_sort_keys( T_GROUP_INFO *t_group_info )
    {
        T_GROUP       *t_group;
        T_GROUP_ISOWT Mult = 1;
        int     i, j, num_t_groups, num_iso_t_groups = 0;

        if (!t_group_info || !( t_group = t_group_info->t_group ) ||
             0 >= ( num_t_groups = t_group_info->num_t_groups ) || t_group_info->nNumIsotopicEndpoints)
        {
            return 0;
        }

        for (i = 0; i < num_t_groups; i++)
        {
            t_group[i].iWeight = 0;
            j = T_NUM_ISOTOPIC - 1;
            Mult = 1;
            do
            {
                t_group[i].iWeight += Mult * (T_GROUP_ISOWT) t_group[i].num[T_NUM_NO_ISOTOPIC + j];
            } while (--j >= 0 && ( Mult *= T_GROUP_ISOWT_MULT ));

            num_iso_t_groups += ( t_group[i].iWeight != 0 );
        }

        return num_iso_t_groups;
    }
        */
    // END INCHI C FUNCTION: set_tautomer_iso_sort_keys

    let Some(t_group_info) = t_group_info else {
        return Ok(0);
    };
    if t_group_info.t_group.is_null()
        || t_group_info.num_t_groups <= 0
        || t_group_info.nNumIsotopicEndpoints != 0
    {
        return Ok(0);
    }

    let count = usize::try_from(t_group_info.num_t_groups)
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let groups = heap.slice_mut(t_group_info.t_group)?;
    if groups.len() < count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut number_of_isotopic_groups = 0_i32;
    for group in &mut groups[..count] {
        group.iWeight = 0;
        let mut isotope = T_NUM_ISOTOPIC as i32 - 1;
        let mut multiplier = 1_i64;
        loop {
            let index = usize::try_from(T_NUM_NO_ISOTOPIC as i32 + isotope)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            group.iWeight = group
                .iWeight
                .wrapping_add(multiplier.wrapping_mul(i64::from(group.num[index])));
            isotope -= 1;
            if isotope < 0 {
                break;
            }
            multiplier =
                multiplier.wrapping_mul(i64::from(crate::source_types::T_GROUP_ISOWT_MULT));
            if multiplier == 0 {
                break;
            }
        }
        number_of_isotopic_groups =
            number_of_isotopic_groups.wrapping_add(i32::from(group.iWeight != 0));
    }
    Ok(number_of_isotopic_groups)
}

pub(crate) fn make_a_copy_of_t_group_info(
    heap: &mut SourceHeap,
    t_group_info: Option<&mut T_GROUP_INFO>,
    t_group_info_orig: Option<&mut T_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6364 make_a_copy_of_t_group_info
    // INCHI✔️❌: complete active source frame follows verbatim; checked typed allocations and copies add overhead.
    /*
    int make_a_copy_of_t_group_info( T_GROUP_INFO *t_group_info,
                                     T_GROUP_INFO *t_group_info_orig )
    {
        int err = 0, len;
        free_t_group_info( t_group_info );
        if (t_group_info_orig && t_group_info)
        {
            if (( len = t_group_info_orig->max_num_t_groups ) > 0)
            {
                T_GROUP* tgi_tg = NULL;  /* Copied from below 2024-09-01 DT */
                /* djb-rwth: fixing oss-fuzz issue #52978 */
                tgi_tg = (T_GROUP*)inchi_malloc(((long long)len+1) * sizeof(t_group_info->t_group[0])); /* djb-rwth: cast operator added */
                if (tgi_tg && t_group_info_orig->t_group)
                {
                    t_group_info->t_group = tgi_tg;
                    memcpy(tgi_tg,
                        t_group_info_orig->t_group,
                        len * sizeof(t_group_info->t_group[0])); /* djb-rwth: cast operator added */
                }
                else
                {
                    if (tgi_tg) /* djb-rwth: avoiding memory leak */
                    {
                        inchi_free(tgi_tg);
                    }
                    err++;
                }
            }
            if (( len = t_group_info_orig->nNumEndpoints ) > 0)
            {
                if ((t_group_info->nEndpointAtomNumber =
                    (AT_NUMB*) inchi_malloc( len * sizeof( t_group_info->nEndpointAtomNumber[0] ) ))) /* djb-rwth: addressing LLVM warning */
                {
                    memcpy(t_group_info->nEndpointAtomNumber,
                        t_group_info_orig->nEndpointAtomNumber,
                        len * sizeof(t_group_info->nEndpointAtomNumber[0]));
                }
                else
                {
                    err++;
                }
            }
            if (( len = t_group_info_orig->num_t_groups ) > 0)
            {
                if ((t_group_info->tGroupNumber =
                    (AT_NUMB*) inchi_malloc( (long long)len * TGSO_TOTAL_LEN * sizeof( t_group_info->tGroupNumber[0] ) ))) /* djb-rwth: cast operator added; djb-rwth: addressing LLVM warning */
                {
                    memcpy(t_group_info->tGroupNumber,
                        t_group_info_orig->tGroupNumber,
                        (long long)len * TGSO_TOTAL_LEN * sizeof(t_group_info->tGroupNumber[0])); /* djb-rwth: cast operator added */
                }
                else
                {
                    err++;
                }
            }
            if (( len = t_group_info_orig->nNumIsotopicEndpoints ) > 0)
            {
                /* djb-rwth: fixing oss-fuzz issue #53519 */
                AT_NUMB* tgi_niean = (AT_NUMB*)inchi_malloc(len * sizeof(t_group_info->nIsotopicEndpointAtomNumber[0]));
                AT_NUMB* tgior_niean = (AT_NUMB*)inchi_realloc(t_group_info_orig->nIsotopicEndpointAtomNumber, len * sizeof(t_group_info_orig->nIsotopicEndpointAtomNumber[0]));
                if (tgi_niean && tgior_niean) /* djb-rwth: addressing LLVM warning */
                {
                    t_group_info->nIsotopicEndpointAtomNumber = tgi_niean;
                    t_group_info_orig->nIsotopicEndpointAtomNumber = tgior_niean;
                    memcpy(tgi_niean,
                        tgior_niean,
                        len * sizeof(t_group_info->nIsotopicEndpointAtomNumber[0]));
                }
                else
                {
                    /* djb-rwth: avoiding memory leaks */
                    if (tgi_niean)
                    {
                        inchi_free(tgi_niean);
                    }
                    if (tgior_niean)
                    {
                        inchi_free(tgior_niean);
                    }
                    err++;
                }
            }
            if (!err)
            {
                t_group_info->nNumEndpoints = t_group_info_orig->nNumEndpoints;
                t_group_info->num_t_groups = t_group_info_orig->num_t_groups;
                t_group_info->max_num_t_groups = t_group_info_orig->max_num_t_groups;
                t_group_info->bIgnoreIsotopic = t_group_info_orig->bIgnoreIsotopic;
                t_group_info->nNumIsotopicEndpoints = t_group_info_orig->nNumIsotopicEndpoints;
                t_group_info->tni = t_group_info_orig->tni;
                /*
                t_group_info->nNumRemovedExplicitH       = t_group_info_orig->nNumRemovedExplicitH;
                t_group_info->nNumRemovedProtons         = t_group_info_orig->nNumRemovedProtons;
                t_group_info->bNormalizationFlags        = t_group_info_orig->bNormalizationFlags;
                */
                /*
                t_group_info->bHardAddedRemovedProtons   = t_group_info_orig->bHardAddedRemovedProtons;
                t_group_info->bSimpleAddedRemovedProtons = t_group_info_orig->bSimpleAddedRemovedProtons;
                t_group_info->nNumCanceledCharges        = t_group_info_orig->nNumCanceledCharges;
                */
            }
            t_group_info->bTautFlags = t_group_info_orig->bTautFlags;
            t_group_info->bTautFlagsDone = t_group_info_orig->bTautFlagsDone;
        }

        return err;
    }
        */
    // END INCHI C FUNCTION: make_a_copy_of_t_group_info
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: make_a_copy_of_t_group_info
    // INCHI✔️❌: inchi_malloc/inchi_realloc/inchi_free select GCC/Linux libc behavior.
    // INCHI✔️❌: TGSO_TOTAL_LEN == 4.
    // END INCHI ACTIVE MACRO CONFIGURATION: make_a_copy_of_t_group_info

    let mut destination = t_group_info;
    free_t_group_info(heap, destination.as_deref_mut())?;
    let (Some(source), Some(destination)) = (t_group_info_orig, destination) else {
        return Ok(0);
    };

    let mut error_count = 0_i32;
    if source.max_num_t_groups > 0 {
        let length = u64::try_from(source.max_num_t_groups)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let allocated = match inchi_calloc::<T_GROUP>(
            heap,
            length + 1,
            std::mem::size_of::<T_GROUP>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if !allocated.is_null() && !source.t_group.is_null() {
            let copied = heap.slice(source.t_group.as_const())?
                [..usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .to_vec();
            destination.t_group = allocated;
            heap.slice_mut(allocated)?[..copied.len()].clone_from_slice(&copied);
        } else {
            inchi_free(heap, allocated)?;
            error_count = error_count.wrapping_add(1);
        }
    }

    if source.nNumEndpoints > 0 {
        let length = u64::try_from(source.nNumEndpoints)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let allocated =
            match inchi_calloc::<AT_NUMB>(heap, length, std::mem::size_of::<AT_NUMB>() as u64) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
        if !allocated.is_null() {
            destination.nEndpointAtomNumber = allocated;
            let copied = heap.slice(source.nEndpointAtomNumber.as_const())?
                [..usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .to_vec();
            heap.slice_mut(allocated)?.clone_from_slice(&copied);
        } else {
            error_count = error_count.wrapping_add(1);
        }
    }

    if source.num_t_groups > 0 {
        let length = u64::try_from(source.num_t_groups)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
            .checked_mul(u64::from(crate::source_types::TGSO_TOTAL_LEN))
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        let allocated =
            match inchi_calloc::<AT_NUMB>(heap, length, std::mem::size_of::<AT_NUMB>() as u64) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
        if !allocated.is_null() {
            destination.tGroupNumber = allocated;
            let copied = heap.slice(source.tGroupNumber.as_const())?
                [..usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .to_vec();
            heap.slice_mut(allocated)?.clone_from_slice(&copied);
        } else {
            error_count = error_count.wrapping_add(1);
        }
    }

    if source.nNumIsotopicEndpoints > 0 {
        let length = u64::try_from(source.nNumIsotopicEndpoints)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let allocated =
            match inchi_calloc::<AT_NUMB>(heap, length, std::mem::size_of::<AT_NUMB>() as u64) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
        let reallocated =
            match inchi_realloc::<AT_NUMB>(heap, source.nIsotopicEndpointAtomNumber, length) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
        if !allocated.is_null() && !reallocated.is_null() {
            destination.nIsotopicEndpointAtomNumber = allocated;
            source.nIsotopicEndpointAtomNumber = reallocated;
            let copied = heap.slice(reallocated.as_const())?.to_vec();
            heap.slice_mut(allocated)?.clone_from_slice(&copied);
        } else {
            inchi_free(heap, allocated)?;
            inchi_free(heap, reallocated)?;
            error_count = error_count.wrapping_add(1);
        }
    }

    if error_count == 0 {
        destination.nNumEndpoints = source.nNumEndpoints;
        destination.num_t_groups = source.num_t_groups;
        destination.max_num_t_groups = source.max_num_t_groups;
        destination.bIgnoreIsotopic = source.bIgnoreIsotopic;
        destination.nNumIsotopicEndpoints = source.nNumIsotopicEndpoints;
        destination.tni = source.tni.clone();
    }
    destination.bTautFlags = source.bTautFlags;
    destination.bTautFlagsDone = source.bTautFlagsDone;

    if error_count != 0 {
        free_t_group_info(heap, Some(destination))?;
    }
    Ok(error_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::ichi_bns::{AllocateAndInitBnData, AllocateAndInitBnStruct};
    use crate::source_types::{ALT_PATH_MODE_TAUTOM, INCHI_CLOCK};

    #[test]
    fn source_port__ichitaut__compranktautomer__line_7008() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![0, AT_RANK::MAX, 7, 7, 1])
            .unwrap();

        assert_eq!(CompRankTautomer(&heap, 2, 3, ranks), Ok(0));
        assert_eq!(CompRankTautomer(&heap, 1, 4, ranks), Ok(65_534));
        assert_eq!(CompRankTautomer(&heap, 4, 1, ranks), Ok(-65_534));
        assert_eq!(CompRankTautomer(&heap, 0, 2, ranks), Ok(-7));
        assert_eq!(
            CompRankTautomer(&heap, 5, 0, ranks),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichitaut__sorttautomergroupsandendpoints__line_7019() {
        let mut heap = SourceHeap::default();
        let mut globals = CANON_GLOBALS::default();
        let empty = T_GROUP_INFO::default();
        assert_eq!(
            SortTautomerGroupsAndEndpoints(
                &mut heap,
                &mut globals,
                &empty,
                4,
                4,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let endpoints = heap.allocate_model_storage(vec![2, 0, 1, 4, 3, 1]).unwrap();
        let groups = heap
            .allocate_model_storage(vec![
                T_GROUP {
                    nNumEndpoints: 3,
                    nFirstEndpointAtNoPos: 0,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    nNumEndpoints: 1,
                    nFirstEndpointAtNoPos: 3,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    nNumEndpoints: 2,
                    nFirstEndpointAtNoPos: 4,
                    ..T_GROUP::default()
                },
            ])
            .unwrap();
        let group_numbers = heap.allocate_model_storage(vec![0, 1, 2]).unwrap();
        let ranks = heap
            .allocate_model_storage(vec![10, 20, 30, 40, 50, 20, 20, 10])
            .unwrap();
        let info = T_GROUP_INFO {
            t_group: groups,
            nEndpointAtomNumber: endpoints,
            tGroupNumber: group_numbers,
            nNumEndpoints: 6,
            num_t_groups: 3,
            ..T_GROUP_INFO::default()
        };

        assert_eq!(
            SortTautomerGroupsAndEndpoints(&mut heap, &mut globals, &info, 5, 8, ranks),
            Ok(3)
        );
        assert_eq!(
            heap.slice(endpoints.as_const()).unwrap(),
            &[0, 1, 2, 4, 1, 3]
        );
        assert_eq!(
            heap.slice(group_numbers.as_const()).unwrap(),
            &[2, 0, 1],
            "equal-rank group numbers retain insertion order"
        );

        let short_info = T_GROUP_INFO {
            nNumEndpoints: 1,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            SortTautomerGroupsAndEndpoints(
                &mut heap,
                &mut globals,
                &short_info,
                0,
                1,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let overflow_groups = heap
            .allocate_model_storage(vec![T_GROUP {
                nNumEndpoints: 2,
                nFirstEndpointAtNoPos: 1,
                ..T_GROUP::default()
            }])
            .unwrap();
        let overflow_endpoints = heap.allocate_model_storage(vec![1, 0]).unwrap();
        let overflow_ranks = heap.allocate_model_storage(vec![1, 2]).unwrap();
        let overflow_info = T_GROUP_INFO {
            t_group: overflow_groups,
            nEndpointAtomNumber: overflow_endpoints,
            nNumEndpoints: 2,
            num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            SortTautomerGroupsAndEndpoints(
                &mut heap,
                &mut globals,
                &overflow_info,
                1,
                2,
                overflow_ranks,
            ),
            Ok(CT_TAUCOUNT_ERR)
        );
        assert_eq!(heap.slice(overflow_endpoints.as_const()).unwrap(), &[1, 0]);
    }

    fn bonded_pair(bond_type: u8) -> Vec<inp_ATOM> {
        let mut atoms = vec![inp_ATOM::default(); 2];
        atoms[0].valence = 1;
        atoms[0].neighbor[0] = 1;
        atoms[0].bond_type[0] = bond_type;
        atoms[1].valence = 1;
        atoms[1].neighbor[0] = 0;
        atoms[1].bond_type[0] = bond_type;
        atoms
    }

    #[test]
    fn source_port__ichitaut__marktautomergroups__line_4336() {
        fn add_bond(atoms: &mut [inp_ATOM], left: usize, right: usize, bond_type: u8) {
            let left_pos = usize::try_from(atoms[left].valence).unwrap();
            let right_pos = usize::try_from(atoms[right].valence).unwrap();
            atoms[left].neighbor[left_pos] = right as AT_NUMB;
            atoms[left].bond_type[left_pos] = bond_type;
            atoms[left].valence += 1;
            atoms[right].neighbor[right_pos] = left as AT_NUMB;
            atoms[right].bond_type[right_pos] = bond_type;
            atoms[right].valence += 1;
        }

        fn standard_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 3];
            atoms[0].el_number = 8;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].num_H = 1;
            atoms[0].nNumAtInRingSystem = 1;
            atoms[1].el_number = 6;
            atoms[1].chem_bonds_valence = 3;
            atoms[2].el_number = 8;
            atoms[2].chem_bonds_valence = 2;
            atoms[2].nNumAtInRingSystem = 1;
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_DOUBLE as u8);
            atoms
        }

        fn network(
            heap: &mut SourceHeap,
            atoms: SourceMutPointer<inp_ATOM>,
            atom_count: i32,
        ) -> (BN_STRUCT, BN_DATA) {
            let mut changed_bonds = 0;
            let bns_pointer =
                AllocateAndInitBnStruct(heap, atoms, atom_count, 0, 0, 1, &mut changed_bonds)
                    .unwrap();
            assert_eq!(changed_bonds, 0);
            let mut bns = heap.slice(bns_pointer.as_const()).unwrap()[0].clone();
            bns.pbTautFlags = heap.allocate_model_storage(vec![0_u64]).unwrap();
            bns.ic = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let bd_pointer = AllocateAndInitBnData(heap, bns.max_vertices).unwrap();
            let bd = heap.slice(bd_pointer.as_const()).unwrap()[0].clone();
            (bns, bd)
        }

        fn assert_rule_case(atoms: Vec<inp_ATOM>, rule_flag: u64, expected_endpoints: &[usize]) {
            let atom_count = atoms.len() as i32;
            let mut heap = SourceHeap::default();
            let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
            let (mut bns, mut bd) = network(&mut heap, atom_pointer, atom_count);
            let mut tgi = T_GROUP_INFO {
                bTautFlags: TG_FLAG_TEST_TAUT__ATOMS as u64 | rule_flag,
                ..T_GROUP_INFO::default()
            };
            let result = MarkTautomerGroups(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atom_pointer,
                atom_count,
                Some(&mut tgi),
                None,
                &mut bns,
                &mut bd,
                0,
            )
            .unwrap();
            assert!(result > 0, "rule flag {rule_flag:#x} did not mutate state");
            let atoms_after = heap.slice(atom_pointer.as_const()).unwrap();
            let endpoint_group = atoms_after[expected_endpoints[0]].endpoint;
            assert_ne!(endpoint_group, 0);
            for &index in expected_endpoints {
                assert_eq!(atoms_after[index].endpoint, endpoint_group);
            }
            let tautomeric_directions = atoms_after
                .iter()
                .flat_map(|atom| atom.bond_type[..usize::try_from(atom.valence).unwrap()].iter())
                .filter(|&&bond| bond & !(BOND_MARK_ALL as u8) == BOND_TAUTOM as u8)
                .count();
            assert!(tautomeric_directions >= 4);
        }

        fn pt22_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 7];
            atoms[0].el_number = 6;
            atoms[1].el_number = 7;
            atoms[2].el_number = 6;
            for atom in &mut atoms[3..] {
                atom.el_number = 1;
                atom.chem_bonds_valence = 1;
            }
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_DOUBLE as u8);
            add_bond(&mut atoms, 0, 3, BOND_SINGLE as u8);
            add_bond(&mut atoms, 0, 4, BOND_SINGLE as u8);
            add_bond(&mut atoms, 2, 5, BOND_SINGLE as u8);
            add_bond(&mut atoms, 2, 6, BOND_SINGLE as u8);
            atoms[0].chem_bonds_valence = 3;
            atoms[0].num_H = 1;
            atoms[1].chem_bonds_valence = 3;
            atoms[2].chem_bonds_valence = 4;
            atoms
        }

        fn pt16_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 5];
            atoms[0].el_number = 8;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].num_H = 1;
            atoms[0].nNumAtInRingSystem = 1;
            atoms[1].el_number = 7;
            atoms[2].el_number = 6;
            for atom in &mut atoms[3..] {
                atom.el_number = 1;
                atom.chem_bonds_valence = 1;
            }
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_DOUBLE as u8);
            add_bond(&mut atoms, 2, 3, BOND_SINGLE as u8);
            add_bond(&mut atoms, 2, 4, BOND_SINGLE as u8);
            atoms[1].chem_bonds_valence = 3;
            atoms[2].chem_bonds_valence = 4;
            atoms
        }

        fn pt06_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 6];
            atoms[0].el_number = 7;
            atoms[1].el_number = 6;
            atoms[2].el_number = 6;
            for atom in &mut atoms[3..] {
                atom.el_number = 1;
                atom.chem_bonds_valence = 1;
            }
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 0, 3, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_DOUBLE as u8);
            add_bond(&mut atoms, 2, 4, BOND_SINGLE as u8);
            add_bond(&mut atoms, 2, 5, BOND_SINGLE as u8);
            atoms[0].chem_bonds_valence = 2;
            atoms[0].num_H = 1;
            atoms[1].chem_bonds_valence = 3;
            atoms[2].chem_bonds_valence = 4;
            atoms
        }

        fn pt39_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 8];
            atoms[0].el_number = 6;
            atoms[1].el_number = 7;
            atoms[2].el_number = 6;
            atoms[3].el_number = 8;
            atoms[3].chem_bonds_valence = 1;
            atoms[3].num_H = 1;
            for atom in &mut atoms[4..] {
                atom.el_number = 1;
                atom.chem_bonds_valence = 1;
            }
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_DOUBLE as u8);
            add_bond(&mut atoms, 1, 3, BOND_SINGLE as u8);
            add_bond(&mut atoms, 0, 4, BOND_SINGLE as u8);
            add_bond(&mut atoms, 0, 5, BOND_SINGLE as u8);
            add_bond(&mut atoms, 2, 6, BOND_SINGLE as u8);
            add_bond(&mut atoms, 2, 7, BOND_SINGLE as u8);
            atoms[0].chem_bonds_valence = 3;
            atoms[0].num_H = 1;
            atoms[1].chem_bonds_valence = 4;
            atoms[2].chem_bonds_valence = 4;
            atoms
        }

        fn pt13_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 3];
            atoms[0].el_number = 8;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].num_H = 1;
            atoms[1].el_number = 6;
            atoms[2].el_number = 6;
            atoms[2].chem_bonds_valence = 4;
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_TRIPLE as u8);
            atoms[1].chem_bonds_valence = 4;
            atoms
        }

        fn pt18_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 3];
            atoms[0].el_number = 8;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].num_H = 1;
            atoms[1].el_number = 6;
            atoms[2].el_number = 7;
            atoms[2].chem_bonds_valence = 3;
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_TRIPLE as u8);
            atoms[1].chem_bonds_valence = 4;
            atoms
        }

        fn keto_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 5];
            atoms[0].el_number = 8;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].num_H = 1;
            atoms[1].el_number = 6;
            atoms[2].el_number = 6;
            for atom in &mut atoms[3..] {
                atom.el_number = 1;
                atom.chem_bonds_valence = 1;
            }
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, 2, BOND_DOUBLE as u8);
            add_bond(&mut atoms, 2, 3, BOND_SINGLE as u8);
            add_bond(&mut atoms, 2, 4, BOND_SINGLE as u8);
            atoms[1].chem_bonds_valence = 3;
            atoms[1].num_H = 1;
            atoms[2].chem_bonds_valence = 4;
            atoms
        }

        fn assert_ring_case(
            case_name: &str,
            mut atoms: Vec<inp_ATOM>,
            rule_flag: u64,
            expected_endpoints: &[usize],
        ) {
            let one_five = rule_flag == TG_FLAG_1_5_TAUT as u64;
            for (position, &index) in expected_endpoints.iter().enumerate() {
                atoms[index].endpoint = if one_five {
                    (position + 1) as AT_NUMB
                } else {
                    9
                };
            }
            let atom_count = atoms.len() as i32;
            let mut heap = SourceHeap::default();
            let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
            let (mut bns, mut bd) = network(&mut heap, atom_pointer, atom_count);
            let initial_groups = if one_five {
                vec![
                    T_GROUP {
                        nGroupNumber: 1,
                        nNumEndpoints: 1,
                        ..T_GROUP::default()
                    },
                    T_GROUP {
                        nGroupNumber: 2,
                        nNumEndpoints: 1,
                        ..T_GROUP::default()
                    },
                ]
            } else {
                vec![T_GROUP {
                    nGroupNumber: 9,
                    nNumEndpoints: expected_endpoints.len() as AT_NUMB,
                    ..T_GROUP::default()
                }]
            };
            let group_count = initial_groups.len() as i32;
            let groups = heap.allocate_model_storage(initial_groups).unwrap();
            let mut tgi = T_GROUP_INFO {
                t_group: groups,
                num_t_groups: group_count,
                max_num_t_groups: group_count,
                bTautFlags: TG_FLAG_TEST_TAUT__ATOMS as u64 | rule_flag,
                ..T_GROUP_INFO::default()
            };
            let result = MarkTautomerGroups(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atom_pointer,
                atom_count,
                Some(&mut tgi),
                None,
                &mut bns,
                &mut bd,
                0,
            )
            .unwrap();
            assert!(
                result > 0,
                "ring case {case_name} with flag {rule_flag:#x} did not mutate state"
            );
            let atoms_after = heap.slice(atom_pointer.as_const()).unwrap();
            let endpoint_group = atoms_after[expected_endpoints[0]].endpoint;
            assert_ne!(endpoint_group, 0);
            for &index in expected_endpoints {
                assert_eq!(atoms_after[index].endpoint, endpoint_group);
            }
        }

        fn one_five_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 5];
            for (index, bond_type) in [
                BOND_SINGLE as u8,
                BOND_DOUBLE as u8,
                BOND_SINGLE as u8,
                BOND_DOUBLE as u8,
            ]
            .into_iter()
            .enumerate()
            {
                add_bond(&mut atoms, index, index + 1, bond_type);
            }
            atoms[0].el_number = 8;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].num_H = 1;
            atoms[0].nNumAtInRingSystem = 1;
            atoms[4].el_number = 8;
            atoms[4].chem_bonds_valence = 2;
            atoms[4].nNumAtInRingSystem = 1;
            for atom in &mut atoms[1..4] {
                atom.el_number = 6;
                atom.chem_bonds_valence = 3;
            }
            atoms
        }

        fn six_member_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 7];
            let bond_types = [
                BOND_SINGLE as u8,
                BOND_DOUBLE as u8,
                BOND_SINGLE as u8,
                BOND_SINGLE as u8,
                BOND_DOUBLE as u8,
                BOND_SINGLE as u8,
            ];
            for (index, bond_type) in bond_types.into_iter().enumerate() {
                add_bond(&mut atoms, index, (index + 1) % 6, bond_type);
            }
            add_bond(&mut atoms, 3, 6, BOND_DOUBLE as u8);
            for atom in &mut atoms[..6] {
                atom.el_number = 6;
                atom.chem_bonds_valence = 3;
                atom.nRingSystem = 1;
                atom.nNumAtInRingSystem = 6;
            }
            atoms[0].el_number = 7;
            atoms[0].chem_bonds_valence = 2;
            atoms[0].num_H = 1;
            atoms[3].chem_bonds_valence = 4;
            atoms[3].bCutVertex = 1;
            atoms[6].el_number = 8;
            atoms[6].chem_bonds_valence = 2;
            atoms
        }

        fn pyrazole_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 5];
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            for index in 1..=4 {
                let next = if index == 4 { 0 } else { index + 1 };
                let bond_type = if index % 2 == 1 {
                    BOND_SINGLE as u8
                } else {
                    BOND_DOUBLE as u8
                };
                add_bond(&mut atoms, index, next, bond_type);
            }
            atoms[0].el_number = 7;
            atoms[0].chem_bonds_valence = 3;
            atoms[1].el_number = 7;
            atoms[1].chem_bonds_valence = 2;
            atoms[1].num_H = 1;
            for atom in &mut atoms[2..] {
                atom.el_number = 6;
                atom.chem_bonds_valence = 3;
            }
            for atom in &mut atoms {
                atom.nRingSystem = 1;
                atom.nNumAtInRingSystem = 5;
            }
            atoms
        }

        fn tropolone_atoms(path_len: usize) -> Vec<inp_ATOM> {
            let first_endpoint = path_len + 1;
            let second_endpoint = path_len + 2;
            let mut atoms = vec![inp_ATOM::default(); path_len + 3];
            add_bond(&mut atoms, 0, 1, BOND_SINGLE as u8);
            add_bond(&mut atoms, 1, first_endpoint, BOND_DOUBLE as u8);
            add_bond(&mut atoms, 0, second_endpoint, BOND_SINGLE as u8);
            for index in 1..=path_len {
                let next = if index == path_len { 0 } else { index + 1 };
                let bond_type = if index % 2 == 1 {
                    BOND_SINGLE as u8
                } else {
                    BOND_DOUBLE as u8
                };
                add_bond(&mut atoms, index, next, bond_type);
            }
            for atom in &mut atoms[..=path_len] {
                atom.el_number = 6;
                atom.chem_bonds_valence = atom.valence.wrapping_add(1);
                atom.nRingSystem = 1;
                atom.nNumAtInRingSystem = (path_len + 1) as AT_NUMB;
            }
            atoms[0].bCutVertex = 1;
            atoms[1].bCutVertex = 1;
            atoms[first_endpoint].el_number = 8;
            atoms[first_endpoint].chem_bonds_valence = 2;
            atoms[second_endpoint].el_number = 8;
            atoms[second_endpoint].chem_bonds_valence = 1;
            atoms[second_endpoint].num_H = 1;
            atoms
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty_bns = BN_STRUCT::default();
        let mut empty_bd = BN_DATA::default();
        assert_eq!(
            MarkTautomerGroups(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                0,
                None,
                None,
                &mut empty_bns,
                &mut empty_bd,
                0,
            ),
            Ok(0)
        );
        let mut disabled = T_GROUP_INFO::default();
        assert_eq!(
            MarkTautomerGroups(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                0,
                Some(&mut disabled),
                None,
                &mut empty_bns,
                &mut empty_bd,
                0,
            ),
            Ok(0)
        );

        let existing_groups = empty_heap
            .allocate_model_storage(vec![T_GROUP::default()])
            .unwrap();
        let mut zero_capacity = T_GROUP_INFO {
            t_group: existing_groups,
            bTautFlags: TG_FLAG_TEST_TAUT__ATOMS as u64,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            MarkTautomerGroups(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                0,
                Some(&mut zero_capacity),
                None,
                &mut empty_bns,
                &mut empty_bd,
                0,
            ),
            Ok(0)
        );
        zero_capacity.max_num_t_groups = -7;
        assert_eq!(
            MarkTautomerGroups(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                0,
                Some(&mut zero_capacity),
                None,
                &mut empty_bns,
                &mut empty_bd,
                0,
            ),
            Ok(-7)
        );

        let mut allocation_heap = SourceHeap::default();
        allocation_heap.fail_after_allocations(0);
        let mut allocation_tgi = T_GROUP_INFO {
            bIgnoreIsotopic: 17,
            bTautFlags: TG_FLAG_TEST_TAUT__ATOMS as u64,
            bTautFlagsDone: 23,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            MarkTautomerGroups(
                &mut allocation_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                0,
                Some(&mut allocation_tgi),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(-1)
        );
        assert_eq!(allocation_tgi.max_num_t_groups, -1);
        assert!(allocation_tgi.t_group.is_null());
        assert_eq!(allocation_tgi.bIgnoreIsotopic, 17);
        assert_eq!(allocation_tgi.bTautFlagsDone, 23);

        let mut path_allocation_heap = SourceHeap::default();
        let path_groups = path_allocation_heap
            .allocate_model_storage(vec![T_GROUP::default()])
            .unwrap();
        let path_atoms = path_allocation_heap
            .allocate_model_storage(Vec::<inp_ATOM>::new())
            .unwrap();
        let allocations_before = path_allocation_heap.live_allocation_count();
        path_allocation_heap.fail_after_allocations(0);
        let mut path_tgi = T_GROUP_INFO {
            t_group: path_groups,
            max_num_t_groups: 1,
            bTautFlags: TG_FLAG_TEST_TAUT__ATOMS as u64,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            MarkTautomerGroups(
                &mut path_allocation_heap,
                &mut CANON_GLOBALS::default(),
                path_atoms,
                0,
                Some(&mut path_tgi),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert_eq!(
            path_allocation_heap.live_allocation_count(),
            allocations_before
        );

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(standard_atoms()).unwrap();
        let (mut bns, mut bd) = network(&mut heap, atoms, 3);
        let mut tgi = T_GROUP_INFO {
            bTautFlags: TG_FLAG_TEST_TAUT__ATOMS as u64,
            ..T_GROUP_INFO::default()
        };
        let result = MarkTautomerGroups(
            &mut heap,
            &mut CANON_GLOBALS::default(),
            atoms,
            3,
            Some(&mut tgi),
            None,
            &mut bns,
            &mut bd,
            0,
        )
        .unwrap();
        assert!(result > 0, "standard 1,3 rule must mutate the source state");
        assert_eq!(tgi.num_t_groups, 1);
        assert_eq!(tgi.nNumEndpoints, 0);
        let atoms_after = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(atoms_after[0].endpoint, 1);
        assert_eq!(atoms_after[2].endpoint, 1);
        assert_eq!(
            atoms_after[1].bond_type[0] & !(BOND_MARK_ALL as u8),
            BOND_TAUTOM as u8
        );
        assert_eq!(
            atoms_after[1].bond_type[1] & !(BOND_MARK_ALL as u8),
            BOND_TAUTOM as u8
        );

        assert_rule_case(pt22_atoms(), TG_FLAG_PT_22_00 as u64, &[0, 2]);
        assert_rule_case(pt16_atoms(), TG_FLAG_PT_16_00 as u64, &[0, 2]);
        assert_rule_case(pt06_atoms(), TG_FLAG_PT_06_00 as u64, &[0, 2]);
        assert_rule_case(pt39_atoms(), TG_FLAG_PT_39_00 as u64, &[0, 2]);
        assert_rule_case(pt13_atoms(), TG_FLAG_PT_13_00 as u64, &[0, 2]);
        assert_rule_case(pt18_atoms(), TG_FLAG_PT_18_00 as u64, &[0, 2]);
        assert_rule_case(keto_atoms(), TG_FLAG_KETO_ENOL_TAUT as u64, &[0, 2]);
        assert_ring_case(
            "1,5 non-ring",
            one_five_atoms(),
            TG_FLAG_1_5_TAUT as u64,
            &[0, 4],
        );
        assert_ring_case("six-member", six_member_atoms(), 0, &[0, 6]);
        assert_ring_case("pyrazole", pyrazole_atoms(), 0, &[0, 1]);
        assert_ring_case("tropolone-7", tropolone_atoms(6), 0, &[7, 8]);
        assert_ring_case("tropolone-5", tropolone_atoms(4), 0, &[5, 6]);

        let mut forbidden_heap = SourceHeap::default();
        let forbidden_atoms = forbidden_heap
            .allocate_model_storage(standard_atoms())
            .unwrap();
        let (mut forbidden_bns, mut forbidden_bd) =
            network(&mut forbidden_heap, forbidden_atoms, 3);
        let center_vertex = forbidden_heap.slice(forbidden_bns.vert.as_const()).unwrap()[1].clone();
        let edge_number = forbidden_heap
            .slice(center_vertex.iedge.as_const())
            .unwrap()[0];
        forbidden_heap.slice_mut(forbidden_bns.edge).unwrap()
            [usize::try_from(edge_number).unwrap()]
        .forbidden = 1;
        let mut forbidden_tgi = T_GROUP_INFO {
            bTautFlags: TG_FLAG_TEST_TAUT__ATOMS as u64,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            MarkTautomerGroups(
                &mut forbidden_heap,
                &mut CANON_GLOBALS::default(),
                forbidden_atoms,
                3,
                Some(&mut forbidden_tgi),
                None,
                &mut forbidden_bns,
                &mut forbidden_bd,
                0,
            ),
            Ok(0)
        );
        assert!(
            forbidden_heap
                .slice(forbidden_atoms.as_const())
                .unwrap()
                .iter()
                .all(|atom| atom.endpoint == 0)
        );
    }

    #[test]
    fn source_port__ichitaut__settautomericbonds__line_1523() {
        let first_direction = T_BONDPOS {
            nAtomNumber: 0,
            neighbor_index: 0,
        };
        let reverse_direction = T_BONDPOS {
            nAtomNumber: 1,
            neighbor_index: 0,
        };

        let marked_single = 0x20 | BOND_SINGLE as u8;
        let mut atoms = bonded_pair(marked_single);
        assert_eq!(
            SetTautomericBonds(&mut atoms, 1, &[first_direction.clone()]),
            Ok(1)
        );
        assert_eq!(atoms[0].bond_type[0], 0x20 | BOND_TAUTOM as u8);
        assert_eq!(atoms[1].bond_type[0], 0x20 | BOND_TAUTOM as u8);

        let mut duplicate = bonded_pair(BOND_DOUBLE as u8);
        assert_eq!(
            SetTautomericBonds(
                &mut duplicate,
                2,
                &[first_direction.clone(), reverse_direction]
            ),
            Ok(1)
        );
        assert_eq!(duplicate[0].bond_type[0], BOND_TAUTOM as u8);
        assert_eq!(duplicate[1].bond_type[0], BOND_TAUTOM as u8);

        let mut alternating = bonded_pair(BOND_ALTERN as u8);
        assert_eq!(
            SetTautomericBonds(&mut alternating, 1, &[first_direction.clone()]),
            Ok(1)
        );
        assert_eq!(alternating[0].bond_type[0], BOND_TAUTOM as u8);
        assert_eq!(alternating[1].bond_type[0], BOND_TAUTOM as u8);

        let marked_tautomeric = 0xa0 | BOND_TAUTOM as u8;
        let mut unchanged = bonded_pair(marked_tautomeric);
        assert_eq!(
            SetTautomericBonds(&mut unchanged, 1, &[first_direction.clone()]),
            Ok(0)
        );
        assert_eq!(unchanged[0].bond_type[0], marked_tautomeric);
        assert_eq!(unchanged[1].bond_type[0], marked_tautomeric);

        let mut missing_reverse = bonded_pair(BOND_SINGLE as u8);
        missing_reverse[1].valence = 0;
        missing_reverse[1].bond_type[0] = 0x55;
        assert_eq!(
            SetTautomericBonds(&mut missing_reverse, 1, &[first_direction.clone()]),
            Ok(1)
        );
        assert_eq!(missing_reverse[0].bond_type[0], BOND_TAUTOM as u8);
        assert_eq!(missing_reverse[1].bond_type[0], 0x55);

        let mut chain = vec![inp_ATOM::default(); 3];
        chain[0].valence = 1;
        chain[0].neighbor[0] = 1;
        chain[0].bond_type[0] = BOND_SINGLE as u8;
        chain[1].valence = 2;
        chain[1].neighbor[0] = 0;
        chain[1].bond_type[0] = BOND_SINGLE as u8;
        chain[1].neighbor[1] = 2;
        chain[1].bond_type[1] = BOND_DOUBLE as u8;
        chain[2].valence = 1;
        chain[2].neighbor[0] = 1;
        chain[2].bond_type[0] = BOND_DOUBLE as u8;
        assert_eq!(
            SetTautomericBonds(
                &mut chain,
                2,
                &[
                    first_direction.clone(),
                    T_BONDPOS {
                        nAtomNumber: 1,
                        neighbor_index: 1,
                    },
                ],
            ),
            Ok(2)
        );
        assert_eq!(chain[0].bond_type[0], BOND_TAUTOM as u8);
        assert_eq!(chain[1].bond_type[0], BOND_TAUTOM as u8);
        assert_eq!(chain[1].bond_type[1], BOND_TAUTOM as u8);
        assert_eq!(chain[2].bond_type[0], BOND_TAUTOM as u8);

        let mut untouched = bonded_pair(BOND_SINGLE as u8);
        let original = untouched.clone();
        assert_eq!(SetTautomericBonds(&mut untouched, -1, &[]), Ok(0));
        assert_eq!(untouched, original);
        assert_eq!(SetTautomericBonds(&mut untouched, 0, &[]), Ok(0));
        assert_eq!(untouched, original);

        assert_eq!(
            SetTautomericBonds(&mut untouched, 1, &[]),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            SetTautomericBonds(
                &mut untouched,
                1,
                &[T_BONDPOS {
                    nAtomNumber: 9,
                    neighbor_index: 0,
                }],
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichitaut__getneutralrepsifneeded__line_1565() {
        fn atom(c_point: AT_NUMB, charge: S_CHAR, endpoint: AT_NUMB) -> inp_ATOM {
            inp_ATOM {
                c_point,
                charge,
                endpoint,
                ..inp_ATOM::default()
            }
        }

        fn group_info(heap: &mut SourceHeap, groups: Vec<C_GROUP>) -> C_GROUP_INFO {
            let num_c_groups = i32::try_from(groups.len()).unwrap();
            C_GROUP_INFO {
                c_group: heap.allocate_model_storage(groups).unwrap(),
                num_c_groups,
                max_num_c_groups: num_c_groups,
                ..C_GROUP_INFO::default()
            }
        }

        let mut heap = SourceHeap::default();
        let cgi = group_info(&mut heap, vec![c_group(1, 0, 1, 0)]);
        let atoms = vec![
            atom(1, 1, 10),
            atom(1, 0, 20),
            atom(0, 0, 20),
            atom(0, 0, 10),
            atom(0, 0, 99),
        ];
        let endpoints = vec![
            T_ENDPOINT {
                nAtomNumber: 1,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 4,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 2,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 0,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 3,
                ..T_ENDPOINT::default()
            },
        ];
        let mut ri = 0;
        let mut rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(
                &heap,
                &mut ri,
                &mut rj,
                &atoms,
                atoms.len() as i32,
                &endpoints,
                endpoints.len() as i32,
                Some(&cgi),
            ),
            Ok(0)
        );
        assert_eq!((ri, rj), (3, 2));

        let charged_alternatives = vec![
            atom(1, 1, 10),
            atom(1, 0, 20),
            atom(2, 0, 20),
            atom(3, 0, 20),
            atom(2, 0, 10),
            atom(3, 0, 10),
        ];
        let alternative_endpoints = [
            T_ENDPOINT {
                nAtomNumber: 2,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 3,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 4,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 5,
                ..T_ENDPOINT::default()
            },
        ];
        ri = 0;
        rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(
                &heap,
                &mut ri,
                &mut rj,
                &charged_alternatives,
                charged_alternatives.len() as i32,
                &alternative_endpoints,
                alternative_endpoints.len() as i32,
                Some(&cgi),
            ),
            Ok(0)
        );
        assert_eq!(rj, 2, "the first alternate charge group is retained");
        assert_eq!(ri, 5, "the rj charge group is excluded for ri");

        ri = 0;
        rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(
                &heap,
                &mut ri,
                &mut rj,
                &atoms,
                atoms.len() as i32,
                &[],
                0,
                Some(&cgi),
            ),
            Ok(0)
        );
        assert_eq!((ri, rj), (3, 2), "the full atom scan is the fallback");

        let mut second_group_heap = SourceHeap::default();
        let second_group_cgi = group_info(
            &mut second_group_heap,
            vec![c_group(2, 0, 1, 0), c_group(1, 0, 1, 0)],
        );
        ri = 0;
        rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(
                &second_group_heap,
                &mut ri,
                &mut rj,
                &atoms,
                atoms.len() as i32,
                &endpoints,
                endpoints.len() as i32,
                Some(&second_group_cgi),
            ),
            Ok(0)
        );
        assert_eq!(
            (ri, rj),
            (0, 1),
            "the source's unconditional loop break checks only c_group[0]"
        );

        let mut enough_neutral_heap = SourceHeap::default();
        let enough_neutral = group_info(&mut enough_neutral_heap, vec![c_group(1, 0, 2, 0)]);
        ri = 0;
        rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(
                &enough_neutral_heap,
                &mut ri,
                &mut rj,
                &atoms,
                atoms.len() as i32,
                &endpoints,
                endpoints.len() as i32,
                Some(&enough_neutral),
            ),
            Ok(0)
        );
        assert_eq!((ri, rj), (0, 1));

        for (case_atoms, case_cgi) in [
            (vec![atom(0, 1, 10)], Some(&cgi)),
            (vec![atom(1, 1, 10), atom(2, 0, 20)], Some(&cgi)),
            (vec![atom(1, 0, 10), atom(1, 0, 20)], Some(&cgi)),
            (vec![atom(1, 1, 10), atom(1, 0, 20)], None),
        ] {
            ri = 0;
            rj = if case_atoms.len() > 1 { 1 } else { u16::MAX };
            assert_eq!(
                GetNeutralRepsIfNeeded(
                    &heap,
                    &mut ri,
                    &mut rj,
                    &case_atoms,
                    case_atoms.len() as i32,
                    &[],
                    0,
                    case_cgi,
                ),
                Ok(0)
            );
            assert_eq!(ri, 0);
        }

        let empty_cgi = C_GROUP_INFO::default();
        ri = 0;
        rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(
                &heap,
                &mut ri,
                &mut rj,
                &atoms,
                atoms.len() as i32,
                &[],
                0,
                Some(&empty_cgi),
            ),
            Ok(0)
        );
        assert_eq!((ri, rj), (0, 1));

        let zero_endpoint_atoms = vec![atom(1, 1, 0), atom(1, 0, 0)];
        ri = 0;
        rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(
                &heap,
                &mut ri,
                &mut rj,
                &zero_endpoint_atoms,
                2,
                &[],
                -1,
                Some(&cgi),
            ),
            Ok(0)
        );
        assert_eq!((ri, rj), (0, 1));

        ri = 0;
        rj = 9;
        assert_eq!(
            GetNeutralRepsIfNeeded(&heap, &mut ri, &mut rj, &atoms, 5, &[], 0, Some(&cgi)),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        ri = 0;
        rj = 1;
        assert_eq!(
            GetNeutralRepsIfNeeded(&heap, &mut ri, &mut rj, &atoms, 5, &[], 1, Some(&cgi)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((ri, rj), (0, 1));
    }

    #[test]
    fn source_port__ichitaut__findaccessibleendpoints__line_1712() {
        fn fixture() -> (SourceHeap, SourceMutPointer<inp_ATOM>, BN_STRUCT, BN_DATA) {
            fn add_bond(atoms: &mut [inp_ATOM], left: usize, right: usize, bond_type: u8) {
                let left_pos = usize::try_from(atoms[left].valence).unwrap();
                let right_pos = usize::try_from(atoms[right].valence).unwrap();
                atoms[left].neighbor[left_pos] = right as AT_NUMB;
                atoms[left].bond_type[left_pos] = bond_type;
                atoms[left].valence += 1;
                atoms[right].neighbor[right_pos] = left as AT_NUMB;
                atoms[right].bond_type[right_pos] = bond_type;
                atoms[right].valence += 1;
            }

            let mut atoms = vec![inp_ATOM::default(); 5];
            for (index, bond_type) in [
                BOND_SINGLE as u8,
                BOND_DOUBLE as u8,
                BOND_SINGLE as u8,
                BOND_DOUBLE as u8,
            ]
            .into_iter()
            .enumerate()
            {
                add_bond(&mut atoms, index, index + 1, bond_type);
            }
            atoms[0].el_number = 8;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].num_H = 1;
            atoms[0].endpoint = 1;
            atoms[0].nNumAtInRingSystem = 1;
            atoms[4].el_number = 8;
            atoms[4].chem_bonds_valence = 2;
            atoms[4].endpoint = 2;
            atoms[4].nNumAtInRingSystem = 1;
            for atom in &mut atoms[1..4] {
                atom.el_number = 6;
                atom.chem_bonds_valence = 3;
            }

            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let mut changed_bonds = 0;
            let bns_pointer =
                AllocateAndInitBnStruct(&mut heap, atoms, 5, 0, 0, 1, &mut changed_bonds).unwrap();
            assert_eq!(changed_bonds, 0);
            let mut bns = heap.slice(bns_pointer.as_const()).unwrap()[0].clone();
            bns.pbTautFlags = heap.allocate_model_storage(vec![0_u64]).unwrap();
            bns.ic = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let bd_pointer = AllocateAndInitBnData(&mut heap, bns.max_vertices).unwrap();
            let bd = heap.slice(bd_pointer.as_const()).unwrap()[0].clone();
            (heap, atoms, bns, bd)
        }

        fn disconnected_fixture() -> (SourceHeap, SourceMutPointer<inp_ATOM>, BN_STRUCT, BN_DATA) {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM::default(); 2])
                .unwrap();
            let mut changed_bonds = 0;
            let bns_pointer =
                AllocateAndInitBnStruct(&mut heap, atoms, 2, 0, 0, 1, &mut changed_bonds).unwrap();
            assert_eq!(changed_bonds, 0);
            let mut bns = heap.slice(bns_pointer.as_const()).unwrap()[0].clone();
            bns.pbTautFlags = heap.allocate_model_storage(vec![0_u64]).unwrap();
            bns.ic = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let bd_pointer = AllocateAndInitBnData(&mut heap, bns.max_vertices).unwrap();
            let bd = heap.slice(bd_pointer.as_const()).unwrap()[0].clone();
            (heap, atoms, bns, bd)
        }

        let (mut heap, atoms, mut bns, mut bd) = fixture();
        let mut endpoints = vec![
            T_ENDPOINT {
                nAtomNumber: 0,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 4,
                ..T_ENDPOINT::default()
            },
        ];
        let mut bonds = vec![
            T_BONDPOS {
                nAtomNumber: 0,
                neighbor_index: 0,
            },
            T_BONDPOS {
                nAtomNumber: 3,
                neighbor_index: 1,
            },
        ];
        let mut endpoint_count = 2;
        let mut bond_count = 2;
        assert_eq!(
            FindAccessibleEndPoints(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                &mut endpoints,
                &mut endpoint_count,
                &mut bonds,
                &mut bond_count,
                &mut bns,
                &mut bd,
                atoms,
                5,
                None,
                ALT_PATH_MODE_TAUTOM as i32,
                0,
            ),
            Ok(2)
        );
        assert_eq!((endpoint_count, bond_count), (2, 2));
        assert_eq!(endpoints[0].nEquNumber, 6);
        assert_eq!(endpoints[1].nEquNumber, 6);

        let (mut heap, atoms, mut bns, mut bd) = disconnected_fixture();
        let mut grouped_endpoints = vec![
            T_ENDPOINT {
                nAtomNumber: 0,
                nGroupNumber: 1,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 1,
                nGroupNumber: 2,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 1,
                nGroupNumber: 2,
                ..T_ENDPOINT::default()
            },
        ];
        let mut grouped_bonds = vec![
            T_BONDPOS {
                nAtomNumber: 10,
                neighbor_index: 0,
            },
            T_BONDPOS {
                nAtomNumber: 20,
                neighbor_index: 1,
            },
            T_BONDPOS {
                nAtomNumber: 30,
                neighbor_index: 2,
            },
        ];
        endpoint_count = 3;
        bond_count = 3;
        assert_eq!(
            FindAccessibleEndPoints(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                &mut grouped_endpoints,
                &mut endpoint_count,
                &mut grouped_bonds,
                &mut bond_count,
                &mut bns,
                &mut bd,
                atoms,
                2,
                None,
                ALT_PATH_MODE_TAUTOM as i32,
                0,
            ),
            Ok(2)
        );
        assert_eq!((endpoint_count, bond_count), (2, 2));
        assert_eq!(
            grouped_endpoints[..2]
                .iter()
                .map(|endpoint| (endpoint.nAtomNumber, endpoint.nEquNumber))
                .collect::<Vec<_>>(),
            vec![(1, 2), (1, 2)]
        );
        assert_eq!(
            grouped_bonds[..2]
                .iter()
                .map(|bond| (bond.nAtomNumber, bond.neighbor_index))
                .collect::<Vec<_>>(),
            vec![(20, 1), (30, 2)]
        );

        let (mut heap, atoms, mut bns, mut bd) = disconnected_fixture();
        let mut no_path_endpoints = vec![
            T_ENDPOINT {
                nAtomNumber: 0,
                ..T_ENDPOINT::default()
            },
            T_ENDPOINT {
                nAtomNumber: 1,
                ..T_ENDPOINT::default()
            },
        ];
        let mut no_path_bonds = bonds.clone();
        endpoint_count = 2;
        bond_count = 2;
        assert_eq!(
            FindAccessibleEndPoints(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                &mut no_path_endpoints,
                &mut endpoint_count,
                &mut no_path_bonds,
                &mut bond_count,
                &mut bns,
                &mut bd,
                atoms,
                2,
                None,
                ALT_PATH_MODE_TAUTOM as i32,
                0,
            ),
            Ok(0)
        );
        assert_eq!((endpoint_count, bond_count), (0, 0));

        let (mut heap, atoms, mut bns, mut bd) = fixture();
        let mut mismatched_endpoints = vec![T_ENDPOINT {
            nAtomNumber: 0,
            nGroupNumber: 7,
            nEquNumber: 99,
            ..T_ENDPOINT::default()
        }];
        let mut mismatched_bonds = Vec::new();
        endpoint_count = 1;
        bond_count = 0;
        assert_eq!(
            FindAccessibleEndPoints(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                &mut mismatched_endpoints,
                &mut endpoint_count,
                &mut mismatched_bonds,
                &mut bond_count,
                &mut bns,
                &mut bd,
                atoms,
                5,
                None,
                ALT_PATH_MODE_TAUTOM as i32,
                0,
            ),
            Ok(0)
        );
        assert_eq!(mismatched_endpoints[0].nEquNumber, 99);
        assert_eq!((endpoint_count, bond_count), (1, 0));

        let (mut heap, atoms, mut bns, mut bd) = fixture();
        let mut overflow_endpoints = (0..=MAXVAL)
            .map(|index| T_ENDPOINT {
                nAtomNumber: (index % 5) as AT_NUMB,
                nGroupNumber: (index + 1) as AT_NUMB,
                ..T_ENDPOINT::default()
            })
            .collect::<Vec<_>>();
        let mut overflow_bonds = vec![T_BONDPOS::default(); overflow_endpoints.len()];
        endpoint_count = overflow_endpoints.len() as i32;
        bond_count = endpoint_count;
        assert_eq!(
            FindAccessibleEndPoints(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                &mut overflow_endpoints,
                &mut endpoint_count,
                &mut overflow_bonds,
                &mut bond_count,
                &mut bns,
                &mut bd,
                atoms,
                5,
                None,
                i32::MAX,
                0,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichitaut__is_centerpoint_elem__line_157() {
        let accepted = [6_u8, 7, 15, 16, 53, 33, 51, 34, 52, 17, 35];
        for element in 0_u8..=118 {
            let expected = i32::from(accepted.contains(&element));
            assert_eq!(is_centerpoint_elem(element), expected, "element {element}");
        }
        assert_eq!(is_centerpoint_elem(u8::MAX), 0);
    }

    #[test]
    fn source_port__ichitaut__is_centerpoint_elem_ket__line_182() {
        for element in u8::MIN..=u8::MAX {
            assert_eq!(
                is_centerpoint_elem_KET(element),
                i32::from(element == 6),
                "element {element}",
            );
        }
    }

    #[test]
    fn source_port__ichitaut__is_centerpoint_elem_strict__line_190() {
        let accepted = [6_u8, 7, 15, 33, 51];
        for element in u8::MIN..=u8::MAX {
            assert_eq!(
                is_centerpoint_elem_strict(element),
                i32::from(accepted.contains(&element)),
                "element {element}",
            );
        }
    }

    #[test]
    fn source_port__ichitaut__addatom2num__line_211() {
        let mut atom = named_atom(b"O", 8, -1, 1, 1, 2, 0);
        atom.num_iso_H = [3, 5, 7];
        let atoms = vec![atom];

        let mut add = [10_u16, 20, 30, 40, 50];
        assert_eq!(AddAtom2num(&mut add, &atoms, 0, 0).unwrap(), 3);
        assert_eq!(add, [13, 21, 37, 45, 53]);

        let mut subtract = [10_u16, 20, 30, 40, 50];
        assert_eq!(AddAtom2num(&mut subtract, &atoms, 0, 1).unwrap(), 3);
        assert_eq!(subtract, [7, 19, 23, 35, 47]);

        let mut fill = [999_u16, 888, 777, 666, 555];
        assert_eq!(AddAtom2num(&mut fill, &atoms, 0, 2).unwrap(), 3);
        assert_eq!(fill, [3, 1, 7, 5, 3]);

        let mut underflow = [0_u16, 0, 0, 0, 0];
        assert_eq!(AddAtom2num(&mut underflow, &atoms, 0, 1).unwrap(), 3);
        assert_eq!(
            underflow,
            [
                u16::MAX - 2,
                u16::MAX,
                u16::MAX - 6,
                u16::MAX - 4,
                u16::MAX - 2
            ]
        );

        let mut short = [0_u16; 4];
        assert_eq!(
            AddAtom2num(&mut short, &atoms, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichitaut__addatom2da__line_250() {
        let mut acidic_donor = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        acidic_donor.at_type = ATT_ACIDIC_CO as u16;
        let atoms = vec![acidic_donor];

        let mut add = [10_u16; 6];
        AddAtom2DA(&mut add, &atoms, 0, 0).unwrap();
        assert_eq!(add, [11, 10, 10, 10, 11, 10]);

        let mut subtract = [10_u16; 6];
        AddAtom2DA(&mut subtract, &atoms, 0, 1).unwrap();
        assert_eq!(subtract, [9, 10, 10, 10, 9, 10]);

        let mut fill = [1_u16, 2, 3, 4, 5, 6];
        AddAtom2DA(&mut fill, &atoms, 0, 2).unwrap();
        assert_eq!(fill, [1, 0, 0, 0, 1, 0]);

        let mut acceptor = named_atom(b"O", 8, -1, 1, 2, 0, 0);
        acceptor.at_type = ATT_ACIDIC_CO as u16;
        let atoms = vec![acceptor];
        let mut acceptor_slots = [0_u16; 6];
        AddAtom2DA(&mut acceptor_slots, &atoms, 0, 0).unwrap();
        assert_eq!(acceptor_slots, [0, 0, 0, 1, 0, 0]);

        let mut invalid = named_atom(b"O", 8, 2, 1, 1, 0, 0);
        invalid.at_type = ATT_ACIDIC_CO as u16;
        let atoms = vec![invalid];
        let mut unchanged = [7_u16; 6];
        AddAtom2DA(&mut unchanged, &atoms, 0, 0).unwrap();
        assert_eq!(unchanged, [7; 6]);
    }

    #[test]
    fn source_port__ichitaut__addendpoint__line_330() {
        let mut existing_atom = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        existing_atom.endpoint = 42;
        let mut existing = taut_endpoint(99, 88, 77, 6, 9);
        existing.num = [1, 2, 3, 4, 5];
        existing.num_DA = [6, 7, 8, 9, 10, 11];
        assert_eq!(AddEndPoint(&mut existing, &[existing_atom], 0), Ok(0));
        assert_eq!(existing.nAtomNumber, 0);
        assert_eq!(existing.nEquNumber, 0);
        assert_eq!(existing.nGroupNumber, 42);
        assert_eq!(existing.num, [0; 5]);
        assert_eq!(existing.num_DA, [6, 7, 8, 9, 10, 11]);

        let mut new_atom = named_atom(b"O", 8, -1, 1, 1, 2, 0);
        new_atom.num_iso_H = [3, 5, 7];
        new_atom.at_type = ATT_ACIDIC_CO as u16;
        let mut new_endpoint = taut_endpoint(99, 88, 77, 6, 9);
        new_endpoint.num = [99; 5];
        new_endpoint.num_DA = [99; 6];
        assert_eq!(AddEndPoint(&mut new_endpoint, &[new_atom], 0), Ok(0));
        assert_eq!(new_endpoint.nAtomNumber, 0);
        assert_eq!(new_endpoint.nEquNumber, 0);
        assert_eq!(new_endpoint.nGroupNumber, 0);
        assert_eq!(new_endpoint.num, [3, 1, 7, 5, 3]);
        assert_eq!(new_endpoint.num_DA, [0, 1, 0, 0, 1, 0]);

        let before = new_endpoint.clone();
        assert_eq!(
            AddEndPoint(&mut new_endpoint, &[], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(new_endpoint, before);
    }

    #[test]
    fn source_port__ichitaut__comp_candidates__line_2930() {
        let candidate = |atnumber, type_, endpoint| S_CANDIDATE {
            atnumber,
            type_,
            endpoint,
            ..S_CANDIDATE::default()
        };

        assert_eq!(
            comp_candidates(&candidate(9, 0, 0), &candidate(1, -1, 7)),
            -1
        );
        assert_eq!(
            comp_candidates(&candidate(9, -1, 7), &candidate(1, 0, 0)),
            1
        );
        assert_eq!(
            comp_candidates(&candidate(9, 0, 7), &candidate(1, 0, 0)),
            -1
        );
        assert_eq!(comp_candidates(&candidate(9, 0, 0), &candidate(1, 0, 7)), 1);
        assert_eq!(
            comp_candidates(&candidate(9, 0, 3), &candidate(1, 0, 7)),
            -4
        );
        assert_eq!(comp_candidates(&candidate(1, 0, 7), &candidate(9, 0, 3)), 4);
        assert_eq!(
            comp_candidates(&candidate(3, 0, 7), &candidate(8, 0, 7)),
            -5
        );
        assert_eq!(comp_candidates(&candidate(8, 0, 0), &candidate(3, 0, 0)), 5);
    }

    #[test]
    fn source_port__ichitaut__marksaltchargegroups2__line_2961() {
        fn salt_atoms(donor_h: bool) -> Vec<inp_ATOM> {
            let mut donor = if donor_h {
                named_atom(b"O", 8, 0, 1, 1, 1, 0)
            } else {
                named_atom(b"O", 8, -1, 1, 1, 0, 0)
            };
            donor.neighbor[0] = 1;
            let carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);
            let mut acceptor = named_atom(b"O", 8, 0, 1, 2, 0, 0);
            acceptor.neighbor[0] = 3;
            vec![donor, carbon.clone(), acceptor, carbon]
        }

        fn allocate_case(
            heap: &mut SourceHeap,
            donor_h: bool,
            capacity: usize,
        ) -> (SourceMutPointer<inp_ATOM>, S_GROUP_INFO, T_GROUP_INFO) {
            let atoms = heap.allocate_model_storage(salt_atoms(donor_h)).unwrap();
            let candidates = heap
                .allocate_model_storage(vec![S_CANDIDATE::default(); capacity])
                .unwrap();
            let groups = heap
                .allocate_model_storage(vec![T_GROUP::default(); capacity.max(1)])
                .unwrap();
            (
                atoms,
                S_GROUP_INFO {
                    s_candidate: candidates,
                    max_num_candidates: capacity as i32,
                    ..S_GROUP_INFO::default()
                },
                T_GROUP_INFO {
                    t_group: groups,
                    max_num_t_groups: capacity.max(1) as i32,
                    ..T_GROUP_INFO::default()
                },
            )
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                1,
                None,
                None,
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );

        let mut null_candidates = S_GROUP_INFO {
            max_num_candidates: 4,
            ..S_GROUP_INFO::default()
        };
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                1,
                Some(&mut null_candidates),
                None,
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );

        let mut invalid_heap = SourceHeap::default();
        let (invalid_atoms, mut invalid_salt, mut invalid_taut) =
            allocate_case(&mut invalid_heap, false, 4);
        invalid_salt.num_candidates = -2;
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut invalid_heap,
                &mut CANON_GLOBALS::default(),
                invalid_atoms,
                4,
                Some(&mut invalid_salt),
                Some(&mut invalid_taut),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(invalid_salt.num_candidates, -2);

        let mut overflow_heap = SourceHeap::default();
        let (overflow_atoms, mut overflow_salt, mut overflow_taut) =
            allocate_case(&mut overflow_heap, false, 1);
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut overflow_heap,
                &mut CANON_GLOBALS::default(),
                overflow_atoms,
                4,
                Some(&mut overflow_salt),
                Some(&mut overflow_taut),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!(overflow_salt.num_candidates, 0);
        assert_eq!(
            overflow_heap
                .slice(overflow_salt.s_candidate.as_const())
                .unwrap()[0]
                .atnumber,
            0
        );

        let mut strict_heap = SourceHeap::default();
        let (strict_atoms, mut strict_salt, mut strict_taut) =
            allocate_case(&mut strict_heap, true, 4);
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut strict_heap,
                &mut CANON_GLOBALS::default(),
                strict_atoms,
                4,
                Some(&mut strict_salt),
                Some(&mut strict_taut),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(strict_salt.num_candidates, 0);
        assert_eq!(strict_taut.bTautFlagsDone, 0);

        let mut permissive_heap = SourceHeap::default();
        let (permissive_atoms, mut permissive_salt, mut permissive_taut) =
            allocate_case(&mut permissive_heap, true, 4);
        permissive_taut.bTautFlags = u64::from(TG_FLAG_ALLOW_NO_NEGTV_O);
        let live_before = permissive_heap.live_allocation_count();
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut permissive_heap,
                &mut CANON_GLOBALS::default(),
                permissive_atoms,
                4,
                Some(&mut permissive_salt),
                Some(&mut permissive_taut),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );
        assert_ne!(
            permissive_taut.bTautFlagsDone & u64::from(TG_FLAG_ALLOW_NO_NEGTV_O_DONE),
            0
        );
        assert_eq!(permissive_heap.live_allocation_count(), live_before);
        assert_eq!(permissive_heap.source_allocation_calls(), 1);

        let mut allocation_heap = SourceHeap::default();
        let (allocation_atoms, mut allocation_salt, mut allocation_taut) =
            allocate_case(&mut allocation_heap, false, 4);
        allocation_heap.fail_after_allocations(0);
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut allocation_heap,
                &mut CANON_GLOBALS::default(),
                allocation_atoms,
                4,
                Some(&mut allocation_salt),
                Some(&mut allocation_taut),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(BNS_OUT_OF_RAM)
        );
        assert_eq!(allocation_heap.source_allocation_calls(), 1);
        assert_eq!(allocation_salt.num_candidates, 0);

        let mut stale_heap = SourceHeap::default();
        let (stale_atoms, mut stale_salt, mut stale_taut) =
            allocate_case(&mut stale_heap, false, 4);
        stale_heap.slice_mut(stale_salt.s_candidate).unwrap()[2].type_ = 2;
        stale_heap.trace_source_allocations();
        assert_eq!(
            MarkSaltChargeGroups2(
                &mut stale_heap,
                &mut CANON_GLOBALS::default(),
                stale_atoms,
                4,
                Some(&mut stale_salt),
                Some(&mut stale_taut),
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(BNS_OUT_OF_RAM)
        );
        assert_eq!(stale_heap.source_allocation_calls(), 0);
        let stale_candidates = stale_heap.slice(stale_salt.s_candidate.as_const()).unwrap();
        assert_eq!(stale_candidates[0].type_, -10);
        assert_eq!(stale_candidates[1].type_, -10);
        assert_eq!(stale_candidates[2].type_, 2);
        assert_eq!(stale_salt.num_candidates, 0);
    }

    #[test]
    fn source_port__ichitaut__marksaltchargegroups__line_3483() {
        use crate::source_types::{BNS_ALT_PATH, BNS_EDGE, BNS_VERTEX};

        fn salt_atoms(donor_h: bool) -> Vec<inp_ATOM> {
            let mut donor = if donor_h {
                named_atom(b"O", 8, 0, 1, 1, 1, 0)
            } else {
                named_atom(b"O", 8, -1, 1, 1, 0, 0)
            };
            donor.neighbor[0] = 1;
            let carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);
            let mut acceptor = named_atom(b"O", 8, 0, 1, 2, 0, 0);
            acceptor.neighbor[0] = 3;
            vec![donor, carbon.clone(), acceptor, carbon]
        }

        fn empty_bns(heap: &mut SourceHeap, atom_count: usize) -> BN_STRUCT {
            let mut vertices = Vec::with_capacity(atom_count);
            for _ in 0..atom_count {
                vertices.push(BNS_VERTEX {
                    iedge: heap.allocate_model_storage(Vec::<i32>::new()).unwrap(),
                    ..BNS_VERTEX::default()
                });
            }
            let mut bns = BN_STRUCT {
                num_atoms: atom_count as i32,
                num_vertices: atom_count as i32,
                max_vertices: atom_count as i32,
                vert: heap.allocate_model_storage(vertices).unwrap(),
                edge: heap.allocate_model_storage(Vec::<BNS_EDGE>::new()).unwrap(),
                max_altp: 1,
                pbTautFlags: heap.allocate_model_storage(vec![0_u64]).unwrap(),
                ..BN_STRUCT::default()
            };
            bns.altp[0] = heap
                .allocate_model_storage(vec![BNS_ALT_PATH::default(); 8])
                .unwrap();
            bns
        }

        fn allocate_case(
            heap: &mut SourceHeap,
            donor_h: bool,
            capacity: usize,
        ) -> (
            SourceMutPointer<inp_ATOM>,
            S_GROUP_INFO,
            T_GROUP_INFO,
            BN_STRUCT,
        ) {
            let atoms = heap.allocate_model_storage(salt_atoms(donor_h)).unwrap();
            let candidates = heap
                .allocate_model_storage(vec![S_CANDIDATE::default(); capacity])
                .unwrap();
            let groups = heap
                .allocate_model_storage(vec![T_GROUP::default(); capacity.max(1)])
                .unwrap();
            (
                atoms,
                S_GROUP_INFO {
                    s_candidate: candidates,
                    max_num_candidates: capacity as i32,
                    ..S_GROUP_INFO::default()
                },
                T_GROUP_INFO {
                    t_group: groups,
                    max_num_t_groups: capacity.max(1) as i32,
                    ..T_GROUP_INFO::default()
                },
                empty_bns(heap, 4),
            )
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        assert_eq!(
            MarkSaltChargeGroups(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                1,
                None,
                None,
                None,
                &mut BN_STRUCT::default(),
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );

        let mut invalid_heap = SourceHeap::default();
        let (invalid_atoms, mut invalid_salt, mut invalid_taut, mut invalid_bns) =
            allocate_case(&mut invalid_heap, false, 4);
        invalid_salt.num_candidates = -1;
        assert_eq!(
            MarkSaltChargeGroups(
                &mut invalid_heap,
                &mut CANON_GLOBALS::default(),
                invalid_atoms,
                4,
                Some(&mut invalid_salt),
                Some(&mut invalid_taut),
                None,
                &mut invalid_bns,
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(invalid_salt.num_candidates, -1);

        let mut overflow_heap = SourceHeap::default();
        let (overflow_atoms, mut overflow_salt, mut overflow_taut, mut overflow_bns) =
            allocate_case(&mut overflow_heap, false, 1);
        assert_eq!(
            MarkSaltChargeGroups(
                &mut overflow_heap,
                &mut CANON_GLOBALS::default(),
                overflow_atoms,
                4,
                Some(&mut overflow_salt),
                Some(&mut overflow_taut),
                None,
                &mut overflow_bns,
                &mut BN_DATA::default(),
                0,
            ),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!(overflow_salt.num_candidates, 0);

        let mut strict_heap = SourceHeap::default();
        let (strict_atoms, mut strict_salt, mut strict_taut, mut strict_bns) =
            allocate_case(&mut strict_heap, true, 4);
        assert_eq!(
            MarkSaltChargeGroups(
                &mut strict_heap,
                &mut CANON_GLOBALS::default(),
                strict_atoms,
                4,
                Some(&mut strict_salt),
                Some(&mut strict_taut),
                None,
                &mut strict_bns,
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(strict_salt.num_candidates, -1);

        let mut tested_heap = SourceHeap::default();
        let (tested_atoms, mut tested_salt, mut tested_taut, mut tested_bns) =
            allocate_case(&mut tested_heap, false, 4);
        assert_eq!(
            MarkSaltChargeGroups(
                &mut tested_heap,
                &mut CANON_GLOBALS::default(),
                tested_atoms,
                4,
                Some(&mut tested_salt),
                Some(&mut tested_taut),
                None,
                &mut tested_bns,
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(tested_salt.num_candidates, 2);
        let candidates = tested_heap
            .slice(tested_salt.s_candidate.as_const())
            .unwrap();
        assert_eq!(candidates[0].atnumber, 0);
        assert_eq!(candidates[0].subtype, SALT_DONOR_Neg as S_CHAR);
        assert_eq!(candidates[1].atnumber, 2);
        assert_eq!(candidates[1].subtype, SALT_ACCEPTOR as S_CHAR);

        let mut permissive_heap = SourceHeap::default();
        let (permissive_atoms, mut permissive_salt, mut permissive_taut, mut permissive_bns) =
            allocate_case(&mut permissive_heap, true, 4);
        permissive_taut.bTautFlags = u64::from(TG_FLAG_ALLOW_NO_NEGTV_O);
        assert_eq!(
            MarkSaltChargeGroups(
                &mut permissive_heap,
                &mut CANON_GLOBALS::default(),
                permissive_atoms,
                4,
                Some(&mut permissive_salt),
                Some(&mut permissive_taut),
                None,
                &mut permissive_bns,
                &mut BN_DATA::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(permissive_salt.num_candidates, 2);
        assert_ne!(
            permissive_taut.bTautFlagsDone & u64::from(TG_FLAG_ALLOW_NO_NEGTV_O_DONE),
            0
        );
    }

    fn atom(
        charge: S_CHAR,
        valence: S_CHAR,
        chem_bonds_valence: S_CHAR,
        num_h: S_CHAR,
    ) -> inp_ATOM {
        inp_ATOM {
            charge,
            valence,
            chem_bonds_valence,
            num_H: num_h,
            ..inp_ATOM::default()
        }
    }

    fn named_atom(
        symbol: &[u8],
        el_number: u8,
        charge: S_CHAR,
        valence: S_CHAR,
        chem_bonds_valence: S_CHAR,
        num_h: S_CHAR,
        n_num_at_in_ring_system: u16,
    ) -> inp_ATOM {
        let mut atom = atom(charge, valence, chem_bonds_valence, num_h);
        for (target, source) in atom.elname.iter_mut().zip(symbol.iter().copied()) {
            *target = source as i8;
        }
        atom.el_number = el_number;
        atom.nNumAtInRingSystem = n_num_at_in_ring_system;
        atom
    }

    fn atoms_with_neighbors(mut center: inp_ATOM, neighbors: &[inp_ATOM]) -> Vec<inp_ATOM> {
        center.valence = neighbors.len() as S_CHAR;
        for index in 0..neighbors.len() {
            center.neighbor[index] = (index + 1) as _;
        }
        let mut atoms = vec![center];
        atoms.extend_from_slice(neighbors);
        atoms
    }

    fn endpoint_info_sentinel() -> ENDPOINT_INFO {
        ENDPOINT_INFO {
            cMoveableCharge: 11,
            cNeutralBondsValence: 12,
            cMobile: 13,
            cDonor: 14,
            cAcceptor: 15,
            cKetoEnolCode: 16,
        }
    }

    #[test]
    fn source_port__ichitaut__bcanbeacpoint__line_2047() {
        let mut subtype = 99;
        let donor = atom(1, 2, 2, 1);
        assert_eq!(bCanBeACPoint(&donor, 1, 1, 2, 3, 3, &mut subtype), 0);
        assert_eq!(subtype, C_SUBTYPE_CHARGED_p_DONOR as S_CHAR);

        subtype = 77;
        assert_eq!(bCanBeACPoint(&donor, 1, 1, 2, 3, 0, &mut subtype), 0);
        assert_eq!(subtype, 77);

        subtype = 0;
        let charged_h_donor = atom(1, 2, 3, 1);
        assert_eq!(
            bCanBeACPoint(&charged_h_donor, 1, 1, 3, 3, 3, &mut subtype),
            1
        );
        assert_eq!(subtype, C_SUBTYPE_CHARGED_H_DONOR as S_CHAR);

        subtype = 0;
        let charged_non_taut = atom(1, 3, 4, 0);
        assert_eq!(
            bCanBeACPoint(&charged_non_taut, 1, 1, 3, 3, 3, &mut subtype),
            1
        );
        assert_eq!(subtype, C_SUBTYPE_CHARGED_NON_TAUT as S_CHAR);

        subtype = 0;
        let charged_acceptor_h = atom(1, 1, 3, 1);
        assert_eq!(
            bCanBeACPoint(&charged_acceptor_h, 1, 1, 3, 3, 3, &mut subtype),
            1
        );
        assert_eq!(subtype, C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR as S_CHAR);

        subtype = 0;
        let charged_acceptor = atom(1, 2, 4, 0);
        assert_eq!(
            bCanBeACPoint(&charged_acceptor, 1, 1, 3, 3, 3, &mut subtype),
            1
        );
        assert_eq!(subtype, C_SUBTYPE_CHARGED_H_ACCEPT as S_CHAR);

        subtype = 0;
        let neutral_non_taut = atom(0, 3, 3, 0);
        assert_eq!(
            bCanBeACPoint(&neutral_non_taut, 1, 1, 3, 3, 3, &mut subtype),
            1
        );
        assert_eq!(subtype, C_SUBTYPE_NEUTRAL_NON_TAUT as S_CHAR);

        subtype = 0;
        let neutral_h_donor = atom(0, 2, 2, 1);
        assert_eq!(
            bCanBeACPoint(&neutral_h_donor, 1, 1, 3, 3, 3, &mut subtype),
            1
        );
        assert_eq!(subtype, C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR);

        subtype = 0;
        let neutral_h_accept = atom(0, 2, 3, 0);
        assert_eq!(
            bCanBeACPoint(&neutral_h_accept, 1, 1, 3, 3, 3, &mut subtype),
            1
        );
        assert_eq!(subtype, C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT as S_CHAR);

        subtype = 12;
        let no_match = atom(2, 1, 1, 0);
        assert_eq!(bCanBeACPoint(&no_match, 1, 1, 3, 3, 3, &mut subtype), 0);
        assert_eq!(subtype, 12);
    }

    #[test]
    fn source_port__ichitaut__getchargetype__line_2181() {
        let neutral_neighbor = named_atom(b"C", 6, 0, 0, 0, 0, 0);

        let nitrogen = named_atom(b"N", 7, 1, 2, 3, 1, 0);
        let atoms = atoms_with_neighbors(
            nitrogen,
            &[neutral_neighbor.clone(), neutral_neighbor.clone()],
        );
        let mut subtype = 99;
        assert_eq!(GetChargeType(&atoms, 0, &mut subtype), 0);
        assert_eq!(subtype, C_SUBTYPE_CHARGED_H_DONOR as S_CHAR);

        let phosphorus = named_atom(b"P", 15, 1, 3, 4, 0, 0);
        let atoms = atoms_with_neighbors(
            phosphorus,
            &[
                neutral_neighbor.clone(),
                neutral_neighbor.clone(),
                neutral_neighbor.clone(),
            ],
        );
        subtype = 88;
        assert_eq!(GetChargeType(&atoms, 0, &mut subtype), 1);
        assert_eq!(subtype, C_SUBTYPE_CHARGED_NON_TAUT as S_CHAR);

        for (symbol, el_number, charge_type) in [
            (b"O".as_slice(), 8, 2),
            (b"S".as_slice(), 16, 3),
            (b"Se".as_slice(), 34, 4),
            (b"Te".as_slice(), 52, 5),
        ] {
            let center = named_atom(symbol, el_number, 1, 2, 3, 0, 5);
            let atoms = atoms_with_neighbors(
                center,
                &[neutral_neighbor.clone(), neutral_neighbor.clone()],
            );
            subtype = 77;
            assert_eq!(GetChargeType(&atoms, 0, &mut subtype), charge_type);
            assert_eq!(subtype, C_SUBTYPE_CHARGED_NON_TAUT as S_CHAR);
        }

        let oxygen_outside_ring = named_atom(b"O", 8, 1, 2, 3, 0, 4);
        let atoms = atoms_with_neighbors(
            oxygen_outside_ring,
            &[neutral_neighbor.clone(), neutral_neighbor.clone()],
        );
        subtype = 66;
        assert_eq!(GetChargeType(&atoms, 0, &mut subtype), -1);
        assert_eq!(subtype, 0);

        let opposite_charge_neighbor = named_atom(b"C", 6, -1, 0, 0, 0, 0);
        let nitrogen = named_atom(b"N", 7, 1, 1, 2, 0, 0);
        let atoms = atoms_with_neighbors(nitrogen, &[opposite_charge_neighbor]);
        subtype = 55;
        assert_eq!(GetChargeType(&atoms, 0, &mut subtype), -1);
        assert_eq!(subtype, 0);

        let mut endpoint_opposite_charge_neighbor = named_atom(b"C", 6, -1, 0, 0, 0, 0);
        endpoint_opposite_charge_neighbor.endpoint = 1;
        let nitrogen = named_atom(b"N", 7, 1, 1, 2, 0, 0);
        let atoms = atoms_with_neighbors(nitrogen, &[endpoint_opposite_charge_neighbor]);
        subtype = 44;
        assert_eq!(GetChargeType(&atoms, 0, &mut subtype), -1);
        assert_eq!(subtype, 0);

        let doubly_charged = named_atom(b"N", 7, 2, 0, 0, 0, 0);
        subtype = 33;
        assert_eq!(GetChargeType(&[doubly_charged], 0, &mut subtype), -1);
        assert_eq!(subtype, 0);

        let neutral = named_atom(b"N", 7, 0, 2, 2, 1, 0);
        subtype = 22;
        assert_eq!(GetChargeType(&[neutral], 0, &mut subtype), 0);
        assert_eq!(subtype, C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR);

        let unknown = named_atom(b"Xx", 0, 0, 0, 0, 0, 0);
        subtype = 11;
        assert_eq!(GetChargeType(&[unknown], 0, &mut subtype), -1);
        assert_eq!(subtype, 0);
    }

    #[test]
    fn source_port__ichitaut__getsaltchargetype__line_2565() {
        let mut heap = SourceHeap::default();
        let mut donor_neg_group = taut_group(7, 1, 5, 0);
        donor_neg_group.num[1] = 2;
        let t_groups = heap
            .allocate_model_storage(vec![donor_neg_group, taut_group(8, 1, 1, 0)])
            .unwrap();
        let t_group_info = T_GROUP_INFO {
            t_group: t_groups,
            num_t_groups: 2,
            max_num_t_groups: 2,
            ..T_GROUP_INFO::default()
        };
        let carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);

        let mut subtype = 99;
        let invalid_valence = vec![named_atom(b"O", 8, 0, 2, 2, 0, 0)];
        assert_eq!(
            GetSaltChargeType(&heap, &invalid_valence, 0, None, &mut subtype).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 88;
        let mut radical = named_atom(b"O", 8, 0, 1, 1, 0, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        assert_eq!(
            GetSaltChargeType(&heap, &[radical], 0, None, &mut subtype).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 77;
        let charge_too_low = vec![named_atom(b"O", 8, -2, 1, 1, 0, 0)];
        assert_eq!(
            GetSaltChargeType(&heap, &charge_too_low, 0, None, &mut subtype).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 66;
        let positive_without_c_point = vec![named_atom(b"O", 8, 1, 1, 1, 0, 0)];
        assert_eq!(
            GetSaltChargeType(&heap, &positive_without_c_point, 0, None, &mut subtype,).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        for (symbol, el_number) in [(b"N".as_slice(), 7), (b"P".as_slice(), 15)] {
            subtype = 55;
            let atom = vec![named_atom(symbol, el_number, 0, 1, 1, 0, 0)];
            assert_eq!(
                GetSaltChargeType(&heap, &atom, 0, None, &mut subtype).unwrap(),
                -1
            );
            assert_eq!(subtype, 0);
        }

        subtype = 44;
        let wrong_valence = vec![named_atom(b"O", 8, 0, 1, 1, 0, 0)];
        assert_eq!(
            GetSaltChargeType(&heap, &wrong_valence, 0, None, &mut subtype).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        for atom in [
            named_atom(b"O", 8, 0, 1, 2, 0, 0),
            named_atom(b"S", 16, 0, 1, 2, 0, 0),
            named_atom(b"Se", 34, 0, 1, 2, 0, 0),
            named_atom(b"Te", 52, 0, 1, 2, 0, 0),
        ] {
            subtype = 33;
            let atoms = {
                let mut atoms = vec![atom];
                atoms[0].neighbor[0] = 1;
                atoms.push(carbon.clone());
                atoms
            };
            assert_eq!(
                GetSaltChargeType(&heap, &atoms, 0, None, &mut subtype).unwrap(),
                0
            );
            assert_eq!(subtype, 4);
        }

        let salt_atom = named_atom(b"O", 8, 0, 1, 2, 0, 0);
        for wrong_carbon in [
            named_atom(b"N", 7, 0, 2, 3, 1, 0),
            named_atom(b"C", 6, 0, 2, 2, 1, 0),
            named_atom(b"C", 6, 1, 2, 3, 1, 0),
            {
                let mut atom = named_atom(b"C", 6, 0, 2, 3, 1, 0);
                atom.radical = RADICAL_SINGLET as S_CHAR + 1;
                atom
            },
            named_atom(b"C", 6, 0, 3, 3, 1, 0),
        ] {
            subtype = 32;
            let atoms = {
                let mut atoms = vec![salt_atom.clone(), wrong_carbon];
                atoms[0].neighbor[0] = 1;
                atoms
            };
            assert_eq!(
                GetSaltChargeType(&heap, &atoms, 0, None, &mut subtype).unwrap(),
                -1
            );
            assert_eq!(subtype, 0);
        }

        subtype = 22;
        let mut taut_atom = named_atom(b"O", 8, 0, 1, 2, 0, 0);
        taut_atom.endpoint = 7;
        let taut_atoms = {
            let mut atoms = vec![taut_atom, carbon.clone()];
            atoms[0].neighbor[0] = 1;
            atoms
        };
        assert_eq!(
            GetSaltChargeType(&heap, &taut_atoms, 0, Some(&t_group_info), &mut subtype,).unwrap(),
            0
        );
        assert_eq!(
            subtype,
            (SALT_DONOR_H | SALT_DONOR_Neg | SALT_ACCEPTOR) as i32
        );

        subtype = 11;
        let mut missing_group_atom = named_atom(b"O", 8, 0, 1, 2, 0, 0);
        missing_group_atom.endpoint = 9;
        let missing_group_atoms = {
            let mut atoms = vec![missing_group_atom, carbon.clone()];
            atoms[0].neighbor[0] = 1;
            atoms
        };
        assert_eq!(
            GetSaltChargeType(
                &heap,
                &missing_group_atoms,
                0,
                Some(&t_group_info),
                &mut subtype,
            )
            .unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 10;
        let mut neutral_h = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        neutral_h.endpoint = 0;
        let neutral_h_atoms = {
            let mut atoms = vec![neutral_h, carbon.clone()];
            atoms[0].neighbor[0] = 1;
            atoms
        };
        assert_eq!(
            GetSaltChargeType(&heap, &neutral_h_atoms, 0, None, &mut subtype).unwrap(),
            0
        );
        assert_eq!(subtype, SALT_DONOR_H as i32);

        subtype = 9;
        let mut neutral_acceptor = named_atom(b"O", 8, 0, 1, 2, 0, 0);
        neutral_acceptor.endpoint = 0;
        let neutral_acceptor_atoms = {
            let mut atoms = vec![neutral_acceptor, carbon.clone()];
            atoms[0].neighbor[0] = 1;
            atoms
        };
        assert_eq!(
            GetSaltChargeType(&heap, &neutral_acceptor_atoms, 0, None, &mut subtype,).unwrap(),
            0
        );
        assert_eq!(subtype, SALT_ACCEPTOR as i32);

        subtype = 8;
        let mut anion = named_atom(b"O", 8, -1, 1, 1, 0, 0);
        anion.endpoint = 0;
        let anion_atoms = {
            let mut atoms = vec![anion, carbon];
            atoms[0].neighbor[0] = 1;
            atoms
        };
        assert_eq!(
            GetSaltChargeType(&heap, &anion_atoms, 0, None, &mut subtype).unwrap(),
            0
        );
        assert_eq!(subtype, SALT_DONOR_Neg as i32);

        subtype = 7;
        let mut positive = named_atom(b"O", 8, 1, 1, 2, 1, 0);
        positive.c_point = 1;
        let positive_atoms = {
            let mut atoms = vec![positive, named_atom(b"C", 6, 0, 2, 3, 1, 0)];
            atoms[0].neighbor[0] = 1;
            atoms
        };
        assert_eq!(
            GetSaltChargeType(&heap, &positive_atoms, 0, None, &mut subtype).unwrap(),
            0
        );
        assert_eq!(subtype, SALT_DONOR_H as i32);
    }

    #[test]
    fn source_port__ichitaut__bdonotmergenontautatom__line_2691() {
        let nitrogen = inp_ATOM {
            el_number: 7,
            ..inp_ATOM::default()
        };
        let nitrogen_atoms = [nitrogen];
        assert_eq!(bDoNotMergeNonTautAtom(&nitrogen_atoms, 0), Ok(1));

        for el_number in [0, 1, 6, 8, u8::MAX] {
            let atom = inp_ATOM {
                el_number,
                ..inp_ATOM::default()
            };
            assert_eq!(bDoNotMergeNonTautAtom(&[atom], 0), Ok(0));
        }

        assert_eq!(
            bDoNotMergeNonTautAtom(&nitrogen_atoms, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            bDoNotMergeNonTautAtom(&nitrogen_atoms, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichitaut__getothersaltchargetype__line_2702() {
        let mut heap = SourceHeap::default();
        let mut donor_group = taut_group(7, 1, 2, 1);
        donor_group.num[1] = 1;
        let t_groups = heap
            .allocate_model_storage(vec![donor_group, taut_group(8, 1, 0, 0)])
            .unwrap();
        let t_group_info = T_GROUP_INFO {
            t_group: t_groups,
            num_t_groups: 2,
            max_num_t_groups: 2,
            ..T_GROUP_INFO::default()
        };

        let center_double = named_atom(b"C", 6, 0, 1, 2, 0, 0);
        let center_single = named_atom(b"C", 6, 0, 1, 2, 0, 0);

        let mut subtype = 99;
        let oxygen_rejected = {
            let mut atom = named_atom(b"O", 8, 0, 1, 1, 0, 0);
            atom.endpoint = 7;
            atom.bond_type[0] = BOND_DOUBLE as u8;
            atom.neighbor[0] = 1;
            vec![atom, center_double.clone()]
        };
        assert_eq!(
            GetOtherSaltChargeType(&heap, &oxygen_rejected, 0, None, &mut subtype, 0).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 88;
        let not_endpoint = {
            let mut atom = named_atom(b"Xx", 0, 0, 1, 1, 0, 0);
            atom.bond_type[0] = BOND_SINGLE as u8;
            atom.neighbor[0] = 1;
            vec![atom, center_single.clone()]
        };
        assert_eq!(
            GetOtherSaltChargeType(&heap, &not_endpoint, 0, None, &mut subtype, 1).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 77;
        let donor = {
            let mut atom = named_atom(b"O", 8, 0, 1, 1, 1, 0);
            atom.bond_type[0] = BOND_SINGLE as u8;
            atom.neighbor[0] = 1;
            vec![atom, center_single.clone()]
        };
        assert_eq!(
            GetOtherSaltChargeType(&heap, &donor, 0, None, &mut subtype, 1).unwrap(),
            1
        );
        assert_eq!(subtype, SALT_DONOR_H as i32);

        subtype = 66;
        let acceptor = {
            let mut atom = named_atom(b"O", 8, 0, 1, 2, 0, 0);
            atom.bond_type[0] = BOND_DOUBLE as u8;
            atom.neighbor[0] = 1;
            vec![atom, center_double.clone()]
        };
        assert_eq!(
            GetOtherSaltChargeType(&heap, &acceptor, 0, None, &mut subtype, 1).unwrap(),
            1
        );
        assert_eq!(subtype, SALT_ACCEPTOR as i32);

        subtype = 55;
        let tgroup = {
            let mut atom = named_atom(b"O", 8, 0, 1, 2, 0, 0);
            atom.endpoint = 7;
            atom.bond_type[0] = BOND_DOUBLE as u8;
            atom.neighbor[0] = 1;
            vec![atom, center_double.clone()]
        };
        assert_eq!(
            GetOtherSaltChargeType(&heap, &tgroup, 0, Some(&t_group_info), &mut subtype, 1)
                .unwrap(),
            1
        );
        assert_eq!(
            subtype,
            (SALT_DONOR_H | SALT_DONOR_Neg | SALT_ACCEPTOR) as i32
        );

        subtype = 44;
        let missing_group = {
            let mut atom = named_atom(b"O", 8, 0, 1, 2, 0, 0);
            atom.endpoint = 9;
            atom.bond_type[0] = BOND_DOUBLE as u8;
            atom.neighbor[0] = 1;
            vec![atom, center_double]
        };
        assert_eq!(
            GetOtherSaltChargeType(
                &heap,
                &missing_group,
                0,
                Some(&t_group_info),
                &mut subtype,
                1,
            )
            .unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 33;
        let wrong_bond_for_acceptor = {
            let mut atom = named_atom(b"O", 8, 0, 1, 2, 0, 0);
            atom.bond_type[0] = BOND_SINGLE as u8;
            atom.neighbor[0] = 1;
            vec![atom, center_single]
        };
        assert_eq!(
            GetOtherSaltChargeType(&heap, &wrong_bond_for_acceptor, 0, None, &mut subtype, 1,)
                .unwrap(),
            -1
        );
        assert_eq!(subtype, 0);
    }

    #[test]
    fn source_port__ichitaut__getothersalttype__line_2828() {
        let mut subtype = 99;
        let wrong_valence = vec![named_atom(b"S", 16, 0, 1, 2, 1, 0)];
        assert_eq!(
            GetOtherSaltType(&wrong_valence, 0, &mut subtype).unwrap(),
            -1
        );
        assert_eq!(subtype, 99);

        subtype = 88;
        let wrong_element = vec![named_atom(b"O", 8, 0, 1, 1, 1, 0)];
        assert_eq!(
            GetOtherSaltType(&wrong_element, 0, &mut subtype).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);

        subtype = 77;
        let mut donor = named_atom(b"S", 16, 0, 1, 1, 1, 0);
        donor.bond_type[0] = BOND_SINGLE as u8;
        donor.neighbor[0] = 1;
        let donor_atoms = vec![donor, named_atom(b"C", 6, 0, 1, 1, 0, 0)];
        assert_eq!(GetOtherSaltType(&donor_atoms, 0, &mut subtype).unwrap(), 2);
        assert_eq!(subtype, SALT_p_DONOR as i32);

        subtype = 66;
        let mut acceptor = named_atom(b"S", 16, -1, 1, 1, 0, 0);
        acceptor.bond_type[0] = BOND_SINGLE as u8;
        acceptor.neighbor[0] = 1;
        let acceptor_atoms = vec![acceptor, named_atom(b"C", 6, 0, 1, 1, 0, 0)];
        assert_eq!(
            GetOtherSaltType(&acceptor_atoms, 0, &mut subtype).unwrap(),
            2
        );
        assert_eq!(subtype, SALT_p_ACCEPTOR as i32);

        subtype = 55;
        let mut wrong_center = named_atom(b"S", 16, 0, 1, 1, 1, 0);
        wrong_center.bond_type[0] = BOND_SINGLE as u8;
        wrong_center.neighbor[0] = 1;
        let mut center = named_atom(b"C", 6, 0, 1, 1, 0, 0);
        center.charge = 1;
        let wrong_center_atoms = vec![wrong_center, center];
        assert_eq!(
            GetOtherSaltType(&wrong_center_atoms, 0, &mut subtype).unwrap(),
            -1
        );
        assert_eq!(subtype, 0);
    }

    fn c_group(
        group_number: AT_NUMB,
        charged_points: u16,
        cpoints: u16,
        proton_points: u16,
    ) -> C_GROUP {
        C_GROUP {
            num: [charged_points, proton_points],
            num_CPoints: cpoints,
            nGroupNumber: group_number,
            cGroupType: 9,
        }
    }

    fn cpoint_atom(c_point: AT_NUMB, charge: S_CHAR, num_h: S_CHAR) -> inp_ATOM {
        inp_ATOM {
            c_point,
            charge,
            num_H: num_h,
            ..inp_ATOM::default()
        }
    }

    fn taut_group(group_number: AT_NUMB, endpoints: AT_NUMB, num0: u16, num_da0: u16) -> T_GROUP {
        let mut group = T_GROUP {
            nGroupNumber: group_number,
            nNumEndpoints: endpoints,
            ..T_GROUP::default()
        };
        group.num[0] = num0;
        group.num_DA[0] = num_da0;
        group
    }

    fn taut_endpoint(
        atom_number: AT_NUMB,
        group_number: AT_NUMB,
        equ_number: AT_NUMB,
        num0: u16,
        num_da0: u16,
    ) -> T_ENDPOINT {
        let mut endpoint = T_ENDPOINT {
            nAtomNumber: atom_number,
            nGroupNumber: group_number,
            nEquNumber: equ_number,
            ..T_ENDPOINT::default()
        };
        endpoint.num[0] = num0;
        endpoint.num_DA[0] = num_da0;
        endpoint
    }

    fn endpoint_atom(endpoint: AT_NUMB) -> inp_ATOM {
        inp_ATOM {
            endpoint,
            ..inp_ATOM::default()
        }
    }

    fn c_candidate(atnumber: AT_NUMB, type_: S_CHAR, subtype: S_CHAR) -> C_CANDIDATE {
        C_CANDIDATE {
            atnumber,
            type_,
            subtype,
        }
    }

    #[test]
    fn source_port__ichitaut__cmpccandidates__line_2229() {
        let base = c_candidate(9, 3, 7);
        assert_eq!(CmpCCandidates(&base, &c_candidate(9, 5, 1)), -2);
        assert_eq!(CmpCCandidates(&base, &c_candidate(9, 1, 20)), 2);
        assert_eq!(CmpCCandidates(&base, &c_candidate(9, 3, 10)), -3);
        assert_eq!(CmpCCandidates(&base, &c_candidate(9, 3, 4)), 3);
        assert_eq!(CmpCCandidates(&base, &c_candidate(12, 3, 7)), -3);
        assert_eq!(CmpCCandidates(&base, &c_candidate(2, 3, 7)), 7);
        assert_eq!(CmpCCandidates(&base, &c_candidate(9, 3, 7)), 0);

        let mut candidates = vec![
            c_candidate(3, 1, 8),
            c_candidate(2, 0, 8),
            c_candidate(4, 1, 4),
            c_candidate(1, 0, 8),
        ];
        candidates.sort_by(|left, right| CmpCCandidates(left, right).cmp(&0));
        assert_eq!(
            candidates,
            vec![
                c_candidate(1, 0, 8),
                c_candidate(2, 0, 8),
                c_candidate(4, 1, 4),
                c_candidate(3, 1, 8),
            ]
        );
    }

    #[test]
    fn source_port__ichitaut__markchargegroups__line_2397() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                named_atom(b"N", 7, 0, 2, 2, 1, 0),
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
            ])
            .unwrap();
        let c_groups = heap
            .allocate_model_storage(vec![C_GROUP::default(); 2])
            .unwrap();
        let c_candidates = heap
            .allocate_model_storage(vec![C_CANDIDATE::default(); 2])
            .unwrap();
        let mut cgi = C_GROUP_INFO {
            c_group: c_groups,
            max_num_c_groups: 2,
            c_candidate: c_candidates,
            max_num_candidates: 2,
            num_candidates: 0,
            ..C_GROUP_INFO::default()
        };
        let mut cg = CANON_GLOBALS::default();
        let mut bns = BN_STRUCT::default();
        let mut bd = BN_DATA::default();

        assert_eq!(
            MarkChargeGroups(
                &mut heap,
                &mut cg,
                atoms,
                2,
                Some(&mut cgi),
                None,
                &mut bns,
                &mut bd,
                0,
            ),
            Ok(0)
        );
        assert_eq!(cgi.num_candidates, -1);
        assert_eq!(
            heap.slice(c_candidates.as_const()).unwrap()[0],
            c_candidate(0, 0, C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR)
        );

        {
            let candidates = heap.slice_mut(c_candidates).unwrap();
            candidates[0] = c_candidate(1, 1, C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR);
            candidates[1] = c_candidate(0, 1, C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR);
        }
        cgi.num_candidates = 2;
        assert_eq!(
            MarkChargeGroups(
                &mut heap,
                &mut cg,
                atoms,
                2,
                Some(&mut cgi),
                None,
                &mut bns,
                &mut bd,
                0,
            ),
            Ok(0)
        );
        assert_eq!(cgi.num_candidates, 2);
        assert_eq!(
            heap.slice(c_candidates.as_const()).unwrap()[..2],
            [
                c_candidate(0, 1, C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR),
                c_candidate(1, 1, C_SUBTYPE_NEUTRAL_H_DONOR as S_CHAR),
            ]
        );

        let mut null_candidate_info = C_GROUP_INFO {
            max_num_candidates: 2,
            ..C_GROUP_INFO::default()
        };
        assert_eq!(
            MarkChargeGroups(
                &mut heap,
                &mut cg,
                atoms,
                2,
                Some(&mut null_candidate_info),
                None,
                &mut bns,
                &mut bd,
                0,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichitaut__registerendpoints__line_1021() {
        let mut cg = CANON_GLOBALS::default();

        let mut empty_heap = SourceHeap::default();
        let empty_groups = empty_heap
            .allocate_model_storage(Vec::<T_GROUP>::new())
            .unwrap();
        let empty_atoms = empty_heap
            .allocate_model_storage(Vec::<inp_ATOM>::new())
            .unwrap();
        let mut empty_tgi = T_GROUP_INFO {
            t_group: empty_groups,
            max_num_t_groups: 0,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            RegisterEndPoints(
                &mut empty_heap,
                &mut cg,
                &mut empty_tgi,
                &mut [],
                0,
                empty_atoms,
                0,
                None,
                None,
            )
            .unwrap(),
            0
        );

        let mut same_heap = SourceHeap::default();
        let same_groups = same_heap
            .allocate_model_storage(vec![taut_group(7, 2, 4, 5)])
            .unwrap();
        let same_atoms = same_heap
            .allocate_model_storage(vec![endpoint_atom(7), endpoint_atom(7)])
            .unwrap();
        let mut same_tgi = T_GROUP_INFO {
            t_group: same_groups,
            num_t_groups: 1,
            max_num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        let mut same_endpoints = vec![taut_endpoint(0, 7, 7, 1, 1), taut_endpoint(1, 7, 7, 2, 2)];
        assert_eq!(
            RegisterEndPoints(
                &mut same_heap,
                &mut cg,
                &mut same_tgi,
                &mut same_endpoints,
                2,
                same_atoms,
                2,
                None,
                None,
            )
            .unwrap(),
            0
        );
        assert_eq!(same_tgi.num_t_groups, 1);
        assert!(same_tgi.tGroupNumber.is_null());

        let mut create_heap = SourceHeap::default();
        let create_groups = create_heap
            .allocate_model_storage(vec![T_GROUP::default(); 3])
            .unwrap();
        let create_atoms = create_heap
            .allocate_model_storage(vec![endpoint_atom(0), endpoint_atom(0)])
            .unwrap();
        let mut create_tgi = T_GROUP_INFO {
            t_group: create_groups,
            max_num_t_groups: 3,
            ..T_GROUP_INFO::default()
        };
        let mut create_endpoints = vec![taut_endpoint(0, 0, 0, 2, 3), taut_endpoint(1, 0, 0, 5, 7)];
        assert_eq!(
            RegisterEndPoints(
                &mut create_heap,
                &mut cg,
                &mut create_tgi,
                &mut create_endpoints,
                2,
                create_atoms,
                2,
                None,
                None,
            )
            .unwrap(),
            2
        );
        assert_eq!(create_tgi.num_t_groups, 1);
        assert_eq!(create_endpoints[0].nEquNumber, 1);
        assert_eq!(create_endpoints[1].nEquNumber, 1);
        let create_group = &create_heap.slice(create_groups.as_const()).unwrap()[0];
        assert_eq!(create_group.nGroupNumber, 1);
        assert_eq!(create_group.nNumEndpoints, 2);
        assert_eq!(create_group.num[0], 7);
        assert_eq!(create_group.num_DA[0], 10);
        assert_eq!(
            create_heap
                .slice(create_atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.endpoint)
                .collect::<Vec<_>>(),
            vec![1, 1]
        );
        assert_eq!(
            create_heap
                .slice(create_tgi.tGroupNumber.as_const())
                .unwrap()[1],
            1
        );

        let mut fict_heap = SourceHeap::default();
        let fict_groups = fict_heap
            .allocate_model_storage(vec![
                taut_group(3, 1, 11, 13),
                T_GROUP::default(),
                T_GROUP::default(),
            ])
            .unwrap();
        let fict_atoms = fict_heap
            .allocate_model_storage(vec![endpoint_atom(0), endpoint_atom(0), endpoint_atom(0)])
            .unwrap();
        let mut fict_tgi = T_GROUP_INFO {
            t_group: fict_groups,
            num_t_groups: 1,
            max_num_t_groups: 3,
            ..T_GROUP_INFO::default()
        };
        let mut fict_endpoints = vec![
            taut_endpoint(0, 0, 20, 1, 1),
            taut_endpoint(1, 0, 22, 2, 3),
            taut_endpoint(2, 0, 20, 4, 5),
        ];
        assert_eq!(
            RegisterEndPoints(
                &mut fict_heap,
                &mut cg,
                &mut fict_tgi,
                &mut fict_endpoints,
                3,
                fict_atoms,
                3,
                None,
                None,
            )
            .unwrap(),
            3
        );
        assert_eq!(
            fict_endpoints
                .iter()
                .map(|endpoint| endpoint.nEquNumber)
                .collect::<Vec<_>>(),
            vec![4, 5, 4]
        );
        assert_eq!(
            fict_heap
                .slice(fict_atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.endpoint)
                .collect::<Vec<_>>(),
            vec![4, 5, 4]
        );

        let mut merge_heap = SourceHeap::default();
        let merge_groups = merge_heap
            .allocate_model_storage(vec![
                taut_group(1, 1, 10, 1),
                taut_group(3, 2, 20, 2),
                taut_group(4, 3, 30, 3),
                T_GROUP::default(),
                T_GROUP::default(),
                T_GROUP::default(),
            ])
            .unwrap();
        let merge_atoms = merge_heap
            .allocate_model_storage(vec![
                endpoint_atom(1),
                endpoint_atom(3),
                endpoint_atom(4),
                endpoint_atom(3),
            ])
            .unwrap();
        let mut merge_tgi = T_GROUP_INFO {
            t_group: merge_groups,
            num_t_groups: 3,
            max_num_t_groups: 6,
            ..T_GROUP_INFO::default()
        };
        let mut merge_endpoints = vec![
            taut_endpoint(0, 1, 1, 0, 0),
            taut_endpoint(1, 3, 1, 0, 0),
            taut_endpoint(2, 4, 1, 0, 0),
            taut_endpoint(3, 3, 1, 0, 0),
        ];
        assert_eq!(
            RegisterEndPoints(
                &mut merge_heap,
                &mut cg,
                &mut merge_tgi,
                &mut merge_endpoints,
                4,
                merge_atoms,
                4,
                None,
                None,
            )
            .unwrap(),
            2
        );
        assert_eq!(merge_tgi.num_t_groups, 1);
        let merged = &merge_heap.slice(merge_groups.as_const()).unwrap()[0];
        assert_eq!(merged.nGroupNumber, 1);
        assert_eq!(merged.nNumEndpoints, 6);
        assert_eq!(merged.num[0], 60);
        assert_eq!(merged.num_DA[0], 6);
        assert_eq!(
            merge_heap
                .slice(merge_atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.endpoint)
                .collect::<Vec<_>>(),
            vec![1, 1, 1, 1]
        );

        let mut mixed_heap = SourceHeap::default();
        let mixed_groups = mixed_heap
            .allocate_model_storage(vec![T_GROUP::default(); 2])
            .unwrap();
        let mixed_atoms = mixed_heap
            .allocate_model_storage(vec![endpoint_atom(0), endpoint_atom(0)])
            .unwrap();
        let mut mixed_tgi = T_GROUP_INFO {
            t_group: mixed_groups,
            max_num_t_groups: 2,
            ..T_GROUP_INFO::default()
        };
        let mut mixed_endpoints = vec![taut_endpoint(0, 0, 0, 0, 0), taut_endpoint(1, 0, 2, 0, 0)];
        assert_eq!(
            RegisterEndPoints(
                &mut mixed_heap,
                &mut cg,
                &mut mixed_tgi,
                &mut mixed_endpoints,
                2,
                mixed_atoms,
                2,
                None,
                None,
            )
            .unwrap(),
            -1
        );
        assert_eq!(mixed_tgi.num_t_groups, 0);

        let mut overflow_heap = SourceHeap::default();
        let overflow_groups = overflow_heap
            .allocate_model_storage(vec![T_GROUP::default()])
            .unwrap();
        let overflow_atoms = overflow_heap
            .allocate_model_storage(vec![endpoint_atom(0), endpoint_atom(0)])
            .unwrap();
        let mut overflow_tgi = T_GROUP_INFO {
            t_group: overflow_groups,
            max_num_t_groups: 0,
            ..T_GROUP_INFO::default()
        };
        let mut overflow_endpoints =
            vec![taut_endpoint(0, 0, 0, 0, 0), taut_endpoint(1, 0, 0, 0, 0)];
        assert_eq!(
            RegisterEndPoints(
                &mut overflow_heap,
                &mut cg,
                &mut overflow_tgi,
                &mut overflow_endpoints,
                2,
                overflow_atoms,
                2,
                None,
                None,
            )
            .unwrap(),
            -1
        );

        let mut fail_index_heap = SourceHeap::default();
        let fail_index_groups = fail_index_heap
            .allocate_model_storage(vec![T_GROUP::default(); 2])
            .unwrap();
        let fail_index_atoms = fail_index_heap
            .allocate_model_storage(vec![endpoint_atom(0)])
            .unwrap();
        let mut fail_index_tgi = T_GROUP_INFO {
            t_group: fail_index_groups,
            max_num_t_groups: 2,
            ..T_GROUP_INFO::default()
        };
        let mut fail_index_endpoints = vec![taut_endpoint(0, 0, 0, 1, 1)];
        fail_index_heap.fail_after_allocations(0);
        assert_eq!(
            RegisterEndPoints(
                &mut fail_index_heap,
                &mut cg,
                &mut fail_index_tgi,
                &mut fail_index_endpoints,
                1,
                fail_index_atoms,
                1,
                None,
                None,
            )
            .unwrap(),
            -1
        );
        assert!(fail_index_tgi.tGroupNumber.is_null());

        let mut fail_scratch_heap = SourceHeap::default();
        let fail_scratch_groups = fail_scratch_heap
            .allocate_model_storage(vec![T_GROUP::default(); 200])
            .unwrap();
        let fail_scratch_atoms = fail_scratch_heap
            .allocate_model_storage(vec![endpoint_atom(0); REGISTER_END_POINTS_STACK_LEN])
            .unwrap();
        let mut fail_scratch_tgi = T_GROUP_INFO {
            t_group: fail_scratch_groups,
            max_num_t_groups: 200,
            ..T_GROUP_INFO::default()
        };
        let mut fail_scratch_endpoints = (0..REGISTER_END_POINTS_STACK_LEN)
            .map(|idx| taut_endpoint(idx as AT_NUMB, 0, (1000 + idx) as AT_NUMB, 1, 0))
            .collect::<Vec<_>>();
        fail_scratch_heap.fail_after_allocations(0);
        assert_eq!(
            RegisterEndPoints(
                &mut fail_scratch_heap,
                &mut cg,
                &mut fail_scratch_tgi,
                &mut fail_scratch_endpoints,
                REGISTER_END_POINTS_STACK_LEN as i32,
                fail_scratch_atoms,
                REGISTER_END_POINTS_STACK_LEN as i32,
                None,
                None,
            )
            .unwrap(),
            -1
        );
    }

    #[test]
    fn source_port__ichitaut__mergesalttautgroups__line_3953() {
        fn salt_pair_atoms() -> Vec<inp_ATOM> {
            let mut donor = named_atom(b"O", 8, 0, 1, 1, 1, 0);
            donor.neighbor[0] = 1;
            donor.bond_type[0] = BOND_SINGLE as u8;
            let mut donor_center = named_atom(b"C", 6, 0, 2, 3, 1, 0);
            donor_center.neighbor[0] = 0;
            donor_center.bond_type[0] = BOND_SINGLE as u8;

            let mut acceptor = named_atom(b"O", 8, -1, 1, 1, 0, 0);
            acceptor.neighbor[0] = 3;
            acceptor.bond_type[0] = BOND_SINGLE as u8;
            let mut acceptor_center = named_atom(b"C", 6, 0, 2, 3, 1, 0);
            acceptor_center.neighbor[0] = 2;
            acceptor_center.bond_type[0] = BOND_SINGLE as u8;

            vec![donor, donor_center, acceptor, acceptor_center]
        }

        fn storage(
            heap: &mut SourceHeap,
            atoms: Vec<inp_ATOM>,
            candidate_capacity: usize,
            group_capacity: usize,
        ) -> (SourceMutPointer<inp_ATOM>, S_GROUP_INFO, T_GROUP_INFO) {
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let candidates = heap
                .allocate_model_storage(vec![S_CANDIDATE::default(); candidate_capacity])
                .unwrap();
            let groups = heap
                .allocate_model_storage(vec![T_GROUP::default(); group_capacity])
                .unwrap();
            (
                atoms,
                S_GROUP_INFO {
                    s_candidate: candidates,
                    max_num_candidates: candidate_capacity as i32,
                    ..S_GROUP_INFO::default()
                },
                T_GROUP_INFO {
                    t_group: groups,
                    max_num_t_groups: group_capacity as i32,
                    ..T_GROUP_INFO::default()
                },
            )
        }

        let mut cg = CANON_GLOBALS::default();
        assert_eq!(
            MergeSaltTautGroups(
                &mut SourceHeap::default(),
                &mut cg,
                SourceMutPointer::default(),
                0,
                None,
                None,
                None,
                None,
            ),
            Ok(0)
        );

        let mut empty_heap = SourceHeap::default();
        let (empty_atoms, mut empty_sgi, mut empty_tgi) = storage(
            &mut empty_heap,
            vec![named_atom(b"C", 6, 0, 0, 0, 0, 0)],
            1,
            1,
        );
        let mut empty_cgi = C_GROUP_INFO::default();
        assert_eq!(
            MergeSaltTautGroups(
                &mut empty_heap,
                &mut cg,
                empty_atoms,
                1,
                Some(&mut empty_sgi),
                Some(&mut empty_tgi),
                Some(&mut empty_cgi),
                None,
            ),
            Ok(0)
        );
        assert_eq!(empty_sgi.num_candidates, -1);

        let mut overflow_heap = SourceHeap::default();
        let overflow_atoms = overflow_heap
            .allocate_model_storage(salt_pair_atoms())
            .unwrap();
        let overflow_candidates = overflow_heap
            .allocate_model_storage(vec![S_CANDIDATE::default()])
            .unwrap();
        let overflow_groups = overflow_heap
            .allocate_model_storage(vec![T_GROUP::default(); 2])
            .unwrap();
        let mut overflow_sgi = S_GROUP_INFO {
            s_candidate: overflow_candidates,
            max_num_candidates: 0,
            ..S_GROUP_INFO::default()
        };
        let mut overflow_tgi = T_GROUP_INFO {
            t_group: overflow_groups,
            max_num_t_groups: 2,
            ..T_GROUP_INFO::default()
        };
        let mut overflow_cgi = C_GROUP_INFO::default();
        assert_eq!(
            MergeSaltTautGroups(
                &mut overflow_heap,
                &mut cg,
                overflow_atoms,
                4,
                Some(&mut overflow_sgi),
                Some(&mut overflow_tgi),
                Some(&mut overflow_cgi),
                None,
            ),
            Ok(BNS_VERT_EDGE_OVFL)
        );

        let mut merge_heap = SourceHeap::default();
        let (merge_atoms, mut merge_sgi, mut merge_tgi) =
            storage(&mut merge_heap, salt_pair_atoms(), 4, 4);
        let merge_candidates = merge_sgi.s_candidate;
        let merge_groups = merge_tgi.t_group;
        let mut merge_cgi = C_GROUP_INFO::default();
        assert_eq!(
            MergeSaltTautGroups(
                &mut merge_heap,
                &mut cg,
                merge_atoms,
                4,
                Some(&mut merge_sgi),
                Some(&mut merge_tgi),
                Some(&mut merge_cgi),
                None,
            ),
            Ok(2)
        );
        assert_eq!(merge_sgi.num_candidates, 2);
        assert_eq!(
            merge_heap.slice(merge_candidates.as_const()).unwrap()[..2]
                .iter()
                .map(|candidate| (candidate.atnumber, candidate.type_, candidate.subtype))
                .collect::<Vec<_>>(),
            vec![
                (0, 0, SALT_DONOR_H as S_CHAR),
                (2, 0, SALT_DONOR_Neg as S_CHAR)
            ]
        );
        assert_eq!(merge_tgi.num_t_groups, 1);
        assert_eq!(
            merge_heap
                .slice(merge_atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.endpoint)
                .collect::<Vec<_>>(),
            vec![1, 0, 1, 0]
        );
        assert_eq!(
            merge_heap.slice(merge_groups.as_const()).unwrap()[0].nNumEndpoints,
            2
        );

        let mut many_atoms = Vec::new();
        for index in 0..128_u16 {
            let base = index * 2;
            let mut oxygen = if index % 2 == 0 {
                named_atom(b"O", 8, 0, 1, 1, 1, 0)
            } else {
                named_atom(b"O", 8, -1, 1, 1, 0, 0)
            };
            oxygen.neighbor[0] = base + 1;
            oxygen.bond_type[0] = BOND_SINGLE as u8;
            let mut carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);
            carbon.neighbor[0] = base;
            carbon.bond_type[0] = BOND_SINGLE as u8;
            many_atoms.push(oxygen);
            many_atoms.push(carbon);
        }
        let mut allocation_heap = SourceHeap::default();
        let (allocation_atoms, mut allocation_sgi, mut allocation_tgi) =
            storage(&mut allocation_heap, many_atoms, 128, 130);
        let mut allocation_cgi = C_GROUP_INFO::default();
        allocation_heap.fail_after_allocations(0);
        assert_eq!(
            MergeSaltTautGroups(
                &mut allocation_heap,
                &mut cg,
                allocation_atoms,
                256,
                Some(&mut allocation_sgi),
                Some(&mut allocation_tgi),
                Some(&mut allocation_cgi),
                None,
            ),
            Ok(BNS_OUT_OF_RAM)
        );
        assert_eq!(allocation_sgi.num_candidates, 128);
    }

    #[test]
    fn source_port__ichitaut__makeisotopichgroup__line_4156() {
        fn storage(
            heap: &mut SourceHeap,
            atoms: Vec<inp_ATOM>,
            candidate_capacity: usize,
        ) -> (SourceMutPointer<inp_ATOM>, S_GROUP_INFO, T_GROUP_INFO) {
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let candidates = heap
                .allocate_model_storage(vec![S_CANDIDATE::default(); candidate_capacity.max(1)])
                .unwrap();
            let groups = heap
                .allocate_model_storage(vec![T_GROUP::default()])
                .unwrap();
            let group_numbers = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
            (
                atoms,
                S_GROUP_INFO {
                    s_candidate: candidates,
                    max_num_candidates: candidate_capacity as i32,
                    ..S_GROUP_INFO::default()
                },
                T_GROUP_INFO {
                    t_group: groups,
                    tGroupNumber: group_numbers,
                    num_iso_H: [91, 92, 93],
                    ..T_GROUP_INFO::default()
                },
            )
        }

        let mut null_heap = SourceHeap::default();
        assert_eq!(
            MakeIsotopicHGroup(&mut null_heap, SourceMutPointer::default(), 0, None, None,),
            Ok(0)
        );
        let mut null_candidates = S_GROUP_INFO::default();
        let mut null_groups = T_GROUP_INFO {
            num_iso_H: [7, 8, 9],
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            MakeIsotopicHGroup(
                &mut null_heap,
                SourceMutPointer::default(),
                0,
                Some(&mut null_candidates),
                Some(&mut null_groups),
            ),
            Ok(0)
        );
        assert_eq!(null_groups.num_iso_H, [7, 8, 9]);

        let mut empty_heap = SourceHeap::default();
        let (empty_atoms, mut empty_sgi, mut empty_tgi) =
            storage(&mut empty_heap, vec![named_atom(b"C", 6, 0, 0, 0, 0, 0)], 1);
        assert_eq!(
            MakeIsotopicHGroup(
                &mut empty_heap,
                empty_atoms,
                1,
                Some(&mut empty_sgi),
                Some(&mut empty_tgi),
            ),
            Ok(0)
        );
        assert_eq!(empty_tgi.num_iso_H, [0, 0, 0]);
        assert!(empty_tgi.nIsotopicEndpointAtomNumber.is_null());

        let mut invalid_heap = SourceHeap::default();
        let mut invalid_atom = named_atom(b"O", 8, 0, 0, 0, 0, 0);
        invalid_atom.endpoint = 1;
        let (invalid_atoms, mut invalid_sgi, mut invalid_tgi) =
            storage(&mut invalid_heap, vec![invalid_atom], 1);
        invalid_heap.slice_mut(invalid_tgi.tGroupNumber).unwrap()[1] = 0;
        assert_eq!(
            MakeIsotopicHGroup(
                &mut invalid_heap,
                invalid_atoms,
                1,
                Some(&mut invalid_sgi),
                Some(&mut invalid_tgi),
            ),
            Ok(BNS_PROGRAM_ERR)
        );
        assert_eq!(invalid_tgi.num_iso_H, [0, 0, 0]);

        let mut overflow_heap = SourceHeap::default();
        let mut overflow_atom = named_atom(b"O", 8, 0, 0, 0, 0, 0);
        overflow_atom.endpoint = 1;
        let (overflow_atoms, mut overflow_sgi, mut overflow_tgi) =
            storage(&mut overflow_heap, vec![overflow_atom], 0);
        let mut group = T_GROUP {
            nGroupNumber: 1,
            ..T_GROUP::default()
        };
        group.num[0] = 1;
        overflow_heap.slice_mut(overflow_tgi.t_group).unwrap()[0] = group;
        assert_eq!(
            MakeIsotopicHGroup(
                &mut overflow_heap,
                overflow_atoms,
                1,
                Some(&mut overflow_sgi),
                Some(&mut overflow_tgi),
            ),
            Ok(BNS_VERT_EDGE_OVFL)
        );
        assert_eq!(overflow_tgi.num_iso_H, [0, 0, 0]);

        fn candidate_atoms() -> Vec<inp_ATOM> {
            let mut tautomer = named_atom(b"O", 8, 0, 0, 0, 0, 0);
            tautomer.endpoint = 1;
            tautomer.num_iso_H = [1, 2, 3];
            tautomer.cFlags = 2;

            let mut oxygen = named_atom(b"O", 8, 0, 1, 1, 1, 0);
            oxygen.neighbor[0] = 2;
            oxygen.bond_type[0] = BOND_SINGLE as u8;
            oxygen.num_iso_H = [4, 5, 6];
            oxygen.cFlags = 2;
            let mut oxygen_center = named_atom(b"C", 6, 0, 2, 3, 1, 0);
            oxygen_center.neighbor[0] = 1;
            oxygen_center.bond_type[0] = BOND_SINGLE as u8;

            let mut sulfur = named_atom(b"S", 16, 0, 1, 1, 1, 0);
            sulfur.neighbor[0] = 4;
            sulfur.bond_type[0] = BOND_SINGLE as u8;
            sulfur.num_iso_H = [7, 8, 9];
            sulfur.cFlags = 2;
            let mut sulfur_center = named_atom(b"C", 6, 0, 1, 1, 0, 0);
            sulfur_center.neighbor[0] = 3;
            sulfur_center.bond_type[0] = BOND_SINGLE as u8;

            let mut nitrogen = named_atom(b"N", 7, 0, 0, 2, 1, 0);
            nitrogen.num_iso_H = [10, 11, 12];
            nitrogen.cFlags = 2;

            vec![
                tautomer,
                oxygen,
                oxygen_center,
                sulfur,
                sulfur_center,
                nitrogen,
            ]
        }

        let mut success_heap = SourceHeap::default();
        let (success_atoms, mut success_sgi, mut success_tgi) =
            storage(&mut success_heap, candidate_atoms(), 6);
        let mut success_group = T_GROUP {
            nGroupNumber: 1,
            ..T_GROUP::default()
        };
        success_group.num[0] = 3;
        success_group.num[1] = 1;
        success_heap.slice_mut(success_tgi.t_group).unwrap()[0] = success_group;
        assert_eq!(
            MakeIsotopicHGroup(
                &mut success_heap,
                success_atoms,
                6,
                Some(&mut success_sgi),
                Some(&mut success_tgi),
            ),
            Ok(4)
        );
        let candidates = success_heap
            .slice(success_sgi.s_candidate.as_const())
            .unwrap();
        assert_eq!(
            candidates[..4]
                .iter()
                .map(|candidate| (
                    candidate.atnumber,
                    candidate.type_,
                    candidate.subtype,
                    candidate.endpoint,
                ))
                .collect::<Vec<_>>(),
            vec![
                (0, 0, 0, 1),
                (1, 0, SALT_DONOR_H as S_CHAR, 0),
                (3, 2, SALT_p_DONOR as S_CHAR, 0),
                (5, 3, SALT_DONOR_H as S_CHAR, 0),
            ]
        );
        assert_eq!(success_tgi.num_iso_H, [22, 26, 30]);
        assert_eq!(success_tgi.nNumIsotopicEndpoints, 4);
        assert_eq!(
            success_heap
                .slice(success_tgi.nIsotopicEndpointAtomNumber.as_const())
                .unwrap(),
            &[3, 1, 3, 5]
        );
        let success_atoms = success_heap.slice(success_atoms.as_const()).unwrap();
        assert_eq!(success_atoms[0].cFlags, 3);
        assert_eq!(success_atoms[1].cFlags, 3);
        assert_eq!(success_atoms[3].cFlags, 3);
        assert_eq!(success_atoms[5].cFlags, 3);
        assert_eq!(success_atoms[2].cFlags, 0);
        assert_eq!(success_atoms[4].cFlags, 0);

        let mut failure_heap = SourceHeap::default();
        let (failure_atoms, mut failure_sgi, mut failure_tgi) =
            storage(&mut failure_heap, candidate_atoms(), 6);
        let mut failure_group = T_GROUP {
            nGroupNumber: 1,
            ..T_GROUP::default()
        };
        failure_group.num[0] = 3;
        failure_group.num[1] = 1;
        failure_heap.slice_mut(failure_tgi.t_group).unwrap()[0] = failure_group;
        failure_tgi.nIsotopicEndpointAtomNumber = failure_heap
            .allocate_model_storage(vec![99_u16, 98])
            .unwrap();
        failure_tgi.nNumIsotopicEndpoints = 77;
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            MakeIsotopicHGroup(
                &mut failure_heap,
                failure_atoms,
                6,
                Some(&mut failure_sgi),
                Some(&mut failure_tgi),
            ),
            Ok(4)
        );
        assert!(failure_tgi.nIsotopicEndpointAtomNumber.is_null());
        assert_eq!(failure_tgi.nNumIsotopicEndpoints, 77);
        assert_eq!(failure_tgi.num_iso_H, [0, 0, 0]);
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert!(
            failure_heap
                .slice(failure_atoms.as_const())
                .unwrap()
                .iter()
                .all(|atom| atom.cFlags != 3)
        );

        let mut wrapping_heap = SourceHeap::default();
        let mut wrapping_atoms = Vec::new();
        for _ in 0..259 {
            let mut atom = named_atom(b"O", 8, 0, 0, 0, 0, 0);
            atom.endpoint = 1;
            atom.num_iso_H = [127, -128, 1];
            wrapping_atoms.push(atom);
        }
        let (wrapping_atoms, mut wrapping_sgi, mut wrapping_tgi) =
            storage(&mut wrapping_heap, wrapping_atoms, 259);
        let mut wrapping_group = T_GROUP {
            nGroupNumber: 1,
            ..T_GROUP::default()
        };
        wrapping_group.num[0] = 1;
        wrapping_heap.slice_mut(wrapping_tgi.t_group).unwrap()[0] = wrapping_group;
        assert_eq!(
            MakeIsotopicHGroup(
                &mut wrapping_heap,
                wrapping_atoms,
                259,
                Some(&mut wrapping_sgi),
                Some(&mut wrapping_tgi),
            ),
            Ok(259)
        );
        assert_eq!(wrapping_tgi.num_iso_H, [-32643, 32384, 259]);
        assert_eq!(wrapping_tgi.nNumIsotopicEndpoints, 1);
        assert_eq!(
            wrapping_heap
                .slice(wrapping_tgi.nIsotopicEndpointAtomNumber.as_const())
                .unwrap(),
            &[0]
        );
    }

    #[test]
    fn source_port__ichitaut__registercpoints__line_2249() {
        let mut same_group_atoms = vec![cpoint_atom(7, 1, 0), cpoint_atom(7, 1, 0)];
        let mut same_groups = vec![c_group(7, 2, 2, 0)];
        let mut same_num_c = 1;
        assert_eq!(
            RegisterCPoints(
                &mut same_groups,
                &mut same_num_c,
                1,
                None,
                0,
                1,
                3,
                &mut same_group_atoms,
                2,
            ),
            0
        );
        assert_eq!(same_num_c, 1);
        assert_eq!(same_groups[0], c_group(7, 2, 2, 0));

        let mut new_groups = vec![
            c_group(4, 2, 2, 0),
            C_GROUP {
                num: [99, 88],
                num_CPoints: 77,
                nGroupNumber: 66,
                cGroupType: 55,
            },
        ];
        let mut new_atoms = vec![cpoint_atom(0, 1, 1), cpoint_atom(0, 0, 0)];
        let mut new_num_c = 1;
        assert_eq!(
            RegisterCPoints(
                &mut new_groups,
                &mut new_num_c,
                2,
                None,
                0,
                1,
                513,
                &mut new_atoms,
                2,
            ),
            1
        );
        assert_eq!(new_num_c, 2);
        assert_eq!(new_atoms[0].c_point, 5);
        assert_eq!(new_atoms[1].c_point, 5);
        assert_eq!(new_groups[1].num, [1, 1]);
        assert_eq!(new_groups[1].num_CPoints, 2);
        assert_eq!(new_groups[1].nGroupNumber, 5);
        assert_eq!(new_groups[1].cGroupType, 1);

        let mut overflow_groups = vec![
            c_group(4, 1, 2, 0),
            C_GROUP {
                num: [8, 8],
                num_CPoints: 8,
                nGroupNumber: 8,
                cGroupType: 8,
            },
        ];
        let mut overflow_atoms = vec![cpoint_atom(0, 1, 0), cpoint_atom(0, 1, 0)];
        let mut overflow_num_c = 1;
        assert_eq!(
            RegisterCPoints(
                &mut overflow_groups,
                &mut overflow_num_c,
                1,
                None,
                0,
                1,
                3,
                &mut overflow_atoms,
                2,
            ),
            BNS_CPOINT_ERR
        );
        assert_eq!(overflow_num_c, 1);
        assert_eq!(overflow_groups[1], C_GROUP::default());
        assert_eq!(overflow_atoms[0].c_point, 0);
        assert_eq!(overflow_atoms[1].c_point, 0);

        let mut add_groups = vec![c_group(4, 1, 2, 0)];
        let mut add_atoms = vec![cpoint_atom(4, 0, 0), cpoint_atom(0, 1, 0)];
        let mut add_num_c = 1;
        assert_eq!(
            RegisterCPoints(
                &mut add_groups,
                &mut add_num_c,
                1,
                None,
                0,
                1,
                3,
                &mut add_atoms,
                2,
            ),
            1
        );
        assert_eq!(add_atoms[1].c_point, 4);
        assert_eq!(add_groups[0].num, [2, 0]);
        assert_eq!(add_groups[0].num_CPoints, 3);

        let mut missing_groups = vec![c_group(4, 1, 2, 0)];
        let mut missing_atoms = vec![cpoint_atom(0, 1, 0), cpoint_atom(5, 0, 0)];
        let mut missing_num_c = 1;
        assert_eq!(
            RegisterCPoints(
                &mut missing_groups,
                &mut missing_num_c,
                1,
                None,
                0,
                1,
                3,
                &mut missing_atoms,
                2,
            ),
            BNS_CPOINT_ERR
        );
        assert_eq!(missing_num_c, 1);
        assert_eq!(missing_atoms[0].c_point, 0);

        let mut merge_groups = vec![
            c_group(1, 1, 2, 7),
            c_group(3, 2, 4, 9),
            c_group(4, 0, 1, 5),
        ];
        let mut merge_atoms = vec![
            cpoint_atom(1, 0, 0),
            cpoint_atom(3, 0, 0),
            cpoint_atom(4, 0, 0),
            cpoint_atom(3, 0, 0),
        ];
        let mut merge_num_c = 3;
        assert_eq!(
            RegisterCPoints(
                &mut merge_groups,
                &mut merge_num_c,
                3,
                None,
                0,
                1,
                3,
                &mut merge_atoms,
                4,
            ),
            1
        );
        assert_eq!(merge_num_c, 2);
        assert_eq!(merge_groups[0].num, [3, 7]);
        assert_eq!(merge_groups[0].num_CPoints, 6);
        assert_eq!(merge_groups[0].nGroupNumber, 1);
        assert_eq!(merge_groups[1].nGroupNumber, 3);
        assert_eq!(
            merge_atoms
                .iter()
                .map(|atom| atom.c_point)
                .collect::<Vec<_>>(),
            vec![1, 1, 3, 1]
        );
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo__line_359() {
        let mut radical = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let carbon = named_atom(b"C", 6, 0, 0, 0, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[carbon], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let saturated_nitrogen = named_atom(b"N", 7, 0, 3, 3, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[saturated_nitrogen], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let abnormal_oxygen = named_atom(b"O", 8, 0, 1, 3, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[abnormal_oxygen], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let nonstandard_oxygen = named_atom(b"O", 8, 0, 1, 1, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[nonstandard_oxygen], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let donor_oxygen = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[donor_oxygen], 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 1,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let acceptor_oxygen = named_atom(b"O", 8, 0, 1, 2, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[acceptor_oxygen], 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 2,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let anion_donor = named_atom(b"O", 8, -1, 1, 1, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&[anion_donor], 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 1,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let unreachable_switch_delta = named_atom(b"O", 8, 0, 0, 2, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo(&[unreachable_switch_delta], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let neutral_neighbor = named_atom(b"C", 6, 0, 0, 0, 0, 0);
        let mut charged_h_donor = named_atom(b"N", 7, 1, 2, 3, 1, 0);
        charged_h_donor.c_point = 1;
        let atoms = atoms_with_neighbors(
            charged_h_donor,
            &[neutral_neighbor.clone(), neutral_neighbor.clone()],
        );
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 2,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let mut charged_h_acceptor = named_atom(b"N", 7, 1, 2, 4, 0, 0);
        charged_h_acceptor.c_point = 1;
        let atoms = atoms_with_neighbors(
            charged_h_acceptor,
            &[neutral_neighbor.clone(), neutral_neighbor],
        );
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 3,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let positive_without_c_point = named_atom(b"N", 7, 1, 2, 3, 1, 0);
        let atoms = atoms_with_neighbors(
            positive_without_c_point,
            &[
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
            ],
        );
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo(&atoms, 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo_pt_22_00__line_452() {
        let mut radical = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_22_00(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let oxygen = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_22_00(&[oxygen], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let saturated_carbon = named_atom(b"C", 6, 0, 4, 4, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_22_00(&[saturated_carbon], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let donor_carbon = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_22_00(&[donor_carbon], 0, &mut eif), 4);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 1,
                cMobile: 3,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let acceptor_carbon = named_atom(b"C", 6, 0, 1, 2, 2, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_22_00(&[acceptor_carbon], 0, &mut eif),
            4
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 2,
                cMobile: 2,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let anion_carbon = named_atom(b"C", 6, -1, 1, 1, 2, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_22_00(&[anion_carbon], 0, &mut eif), 4);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 1,
                cMobile: 3,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let default_delta = named_atom(b"C", 6, 0, 0, 2, 2, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_22_00(&[default_delta], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let mut positive_carbon = named_atom(b"C", 6, 1, 1, 3, 1, 0);
        positive_carbon.c_point = 1;
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_22_00(&[positive_carbon], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo_pt_16_00__line_524() {
        let mut radical = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_16_00(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let excluded_carbon = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_16_00(&[excluded_carbon], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let excluded_oxygen = named_atom(b"O", 8, 0, 2, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_16_00(&[excluded_oxygen], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let donor_carbon = named_atom(b"C", 6, 0, 2, 2, 2, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_16_00(&[donor_carbon], 0, &mut eif), 4);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 2,
                cMobile: 2,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let acceptor_carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_16_00(&[acceptor_carbon], 0, &mut eif),
            4
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 1,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_point = named_atom(b"N", 8, 1, 1, 2, 2, 0);
        charge_point.c_point = 1;
        let atoms = atoms_with_neighbors(charge_point, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_16_00(&atoms, 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 0,
                cMobile: 2,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo_pt_06_00__line_600() {
        let mut radical = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_06_00(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let non_endpoint = named_atom(b"P", 15, 0, 1, 1, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_06_00(&[non_endpoint], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let saturated_carbon = named_atom(b"C", 6, 0, 4, 4, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_06_00(&[saturated_carbon], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let abnormal_valence = named_atom(b"O", 8, 0, 1, 3, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_06_00(&[abnormal_valence], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let non_standard = named_atom(b"C", 6, 0, 3, 3, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_06_00(&[non_standard], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let donor_carbon = named_atom(b"C", 6, 0, 3, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_06_00(&[donor_carbon], 0, &mut eif), 4);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let acceptor_carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_06_00(&[acceptor_carbon], 0, &mut eif),
            4
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 1,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_donor = named_atom(b"N", 7, 1, 1, 2, 2, 0);
        charge_donor.c_point = 1;
        let atoms = atoms_with_neighbors(charge_donor, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_06_00(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 1,
                cMobile: 2,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_acceptor = named_atom(b"N", 7, 1, 2, 4, 0, 0);
        charge_acceptor.c_point = 1;
        let atoms = atoms_with_neighbors(
            charge_acceptor,
            &[
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
            ],
        );
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_06_00(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 3,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo_pt_39_00__line_677() {
        let mut radical = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_39_00(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let oxygen_not_endpoint = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_39_00(&[oxygen_not_endpoint], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let saturated_nitrogen = named_atom(b"N", 7, 0, 3, 3, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_39_00(&[saturated_nitrogen], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let abnormal_valence = named_atom(b"N", 7, 0, 2, 4, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_39_00(&[abnormal_valence], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let non_standard = named_atom(b"N", 7, 0, 2, 2, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_39_00(&[non_standard], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let donor_nitrogen = named_atom(b"N", 7, 0, 2, 2, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_39_00(&[donor_nitrogen], 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 2,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let acceptor_carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_39_00(&[acceptor_carbon], 0, &mut eif),
            4
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 1,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_donor = named_atom(b"N", 7, 1, 1, 2, 2, 0);
        charge_donor.c_point = 1;
        let atoms = atoms_with_neighbors(charge_donor, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_39_00(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 1,
                cMobile: 2,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_acceptor = named_atom(b"N", 7, 1, 2, 4, 0, 0);
        charge_acceptor.c_point = 1;
        let atoms = atoms_with_neighbors(
            charge_acceptor,
            &[
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
            ],
        );
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_39_00(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 3,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo_pt_13_00__line_756() {
        let mut radical = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_13_00(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let nitrogen_not_endpoint = named_atom(b"N", 7, 0, 1, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_13_00(&[nitrogen_not_endpoint], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        for (symbol, el_number) in [
            (b"O".as_slice(), 8),
            (b"S".as_slice(), 16),
            (b"Se".as_slice(), 34),
            (b"Te".as_slice(), 52),
        ] {
            let donor = named_atom(symbol, el_number, 0, 1, 1, 1, 0);
            eif = endpoint_info_sentinel();
            assert_eq!(nGetEndpointInfo_PT_13_00(&[donor], 0, &mut eif), 2);
            assert_eq!(
                eif,
                ENDPOINT_INFO {
                    cMoveableCharge: 0,
                    cNeutralBondsValence: 1,
                    cMobile: 1,
                    cDonor: 1,
                    cAcceptor: 0,
                    cKetoEnolCode: 0,
                }
            );
        }

        let saturated_oxygen = named_atom(b"O", 8, 0, 2, 2, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_13_00(&[saturated_oxygen], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let abnormal_valence = named_atom(b"O", 8, 0, 1, 3, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_13_00(&[abnormal_valence], 0, &mut eif),
            0
        );
        assert_eq!(eif, endpoint_info_sentinel());

        let non_standard = named_atom(b"C", 6, 0, 3, 3, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_13_00(&[non_standard], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let acceptor_oxygen = named_atom(b"O", 8, 0, 1, 2, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_13_00(&[acceptor_oxygen], 0, &mut eif),
            2
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 2,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let donor_carbon = named_atom(b"C", 6, 0, 3, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_13_00(&[donor_carbon], 0, &mut eif), 4);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_donor = named_atom(b"N", 8, 1, 1, 2, 2, 0);
        charge_donor.c_point = 1;
        let atoms = atoms_with_neighbors(charge_donor, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_13_00(&atoms, 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 0,
                cMobile: 2,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_acceptor = named_atom(b"N", 8, 1, 1, 3, 1, 0);
        charge_acceptor.c_point = 1;
        let atoms = atoms_with_neighbors(charge_acceptor, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_13_00(&atoms, 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 1,
                cMobile: 1,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo_pt_18_00__line_832() {
        let mut radical = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_18_00(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let carbon = named_atom(b"C", 6, 0, 3, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_18_00(&[carbon], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let neutral_oxygen_donor = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_18_00(&[neutral_oxygen_donor], 0, &mut eif),
            2
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 1,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let neutral_nitrogen_acceptor = named_atom(b"N", 7, -1, 1, 1, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_18_00(&[neutral_nitrogen_acceptor], 0, &mut eif),
            3
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_donor = named_atom(b"N", 7, 1, 1, 2, 1, 0);
        charge_donor.c_point = 1;
        let atoms = atoms_with_neighbors(charge_donor, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_18_00(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 2,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let mut charge_acceptor = named_atom(b"N", 7, 1, 2, 3, 0, 0);
        charge_acceptor.c_point = 1;
        let atoms = atoms_with_neighbors(
            charge_acceptor,
            &[
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
                named_atom(b"C", 6, 0, 0, 0, 0, 0),
            ],
        );
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_18_00(&atoms, 0, &mut eif), 3);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 3,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );

        let mut fallback_donor = named_atom(b"O", 8, 1, 1, 1, 1, 0);
        fallback_donor.c_point = 0;
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_PT_18_00(&[fallback_donor], 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 1,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 0,
            }
        );

        let mut fallback_acceptor = named_atom(b"O", 8, 1, 1, 1, 0, 0);
        fallback_acceptor.c_point = 0;
        eif = endpoint_info_sentinel();
        assert_eq!(
            nGetEndpointInfo_PT_18_00(&[fallback_acceptor], 0, &mut eif),
            2
        );
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 2,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 0,
            }
        );
    }

    #[test]
    fn source_port__ichitaut__ngetendpointinfo_ket__line_916() {
        let mut radical = named_atom(b"C", 6, 0, 1, 1, 3, 0);
        radical.radical = RADICAL_SINGLET as S_CHAR + 1;
        let mut eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[radical], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let non_endpoint = named_atom(b"N", 7, 0, 1, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[non_endpoint], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let excluded_o = named_atom(b"O", 8, 0, 2, 2, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[excluded_o], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let excluded_c = named_atom(b"C", 6, 0, 1, 1, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[excluded_c], 0, &mut eif), 0);
        assert_eq!(eif, endpoint_info_sentinel());

        let donor_oxygen = named_atom(b"O", 8, 0, 1, 1, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[donor_oxygen], 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 1,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 1,
            }
        );

        let acceptor_oxygen = named_atom(b"O", 8, 0, 1, 2, 0, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[acceptor_oxygen], 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 2,
                cMobile: 0,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 1,
            }
        );

        let donor_carbon = named_atom(b"C", 6, 0, 3, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[donor_carbon], 0, &mut eif), 4);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 1,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 2,
            }
        );

        let acceptor_carbon = named_atom(b"C", 6, 0, 2, 3, 1, 0);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&[acceptor_carbon], 0, &mut eif), 4);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 0,
                cNeutralBondsValence: 3,
                cMobile: 1,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 2,
            }
        );

        let mut charge_donor = named_atom(b"N", 8, 1, 1, 2, 2, 0);
        charge_donor.c_point = 1;
        let atoms = atoms_with_neighbors(charge_donor, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&atoms, 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 0,
                cMobile: 2,
                cDonor: 1,
                cAcceptor: 0,
                cKetoEnolCode: 1,
            }
        );

        let mut charge_acceptor = named_atom(b"N", 8, 1, 1, 3, 1, 0);
        charge_acceptor.c_point = 1;
        let atoms = atoms_with_neighbors(charge_acceptor, &[named_atom(b"C", 6, 0, 0, 0, 0, 0)]);
        eif = endpoint_info_sentinel();
        assert_eq!(nGetEndpointInfo_KET(&atoms, 0, &mut eif), 2);
        assert_eq!(
            eif,
            ENDPOINT_INFO {
                cMoveableCharge: 1,
                cNeutralBondsValence: 1,
                cMobile: 1,
                cDonor: 0,
                cAcceptor: 1,
                cKetoEnolCode: 1,
            }
        );
    }

    #[test]
    fn source_port__ichitaut__free_t_group_info__line_6336() {
        let mut empty_heap = SourceHeap::default();
        assert_eq!(free_t_group_info(&mut empty_heap, None), Ok(0));
        assert_eq!(empty_heap.live_allocation_count(), 0);
        let mut empty_info = T_GROUP_INFO::default();
        assert_eq!(
            free_t_group_info(&mut empty_heap, Some(&mut empty_info)),
            Ok(0)
        );
        assert_eq!(empty_info, T_GROUP_INFO::default());
        assert_eq!(empty_heap.live_allocation_count(), 0);

        let mut heap = SourceHeap::default();
        let groups = heap
            .allocate_model_storage(vec![T_GROUP::default(), T_GROUP::default()])
            .unwrap();
        let endpoints = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let group_numbers = heap.allocate_model_storage(vec![4_u16, 5, 6, 7]).unwrap();
        let isotopic_endpoints = heap.allocate_model_storage(vec![8_u16, 9]).unwrap();
        let mut info = T_GROUP_INFO {
            t_group: groups,
            nEndpointAtomNumber: endpoints,
            tGroupNumber: group_numbers,
            nNumEndpoints: 3,
            num_t_groups: 2,
            max_num_t_groups: 17,
            bIgnoreIsotopic: -1,
            nIsotopicEndpointAtomNumber: isotopic_endpoints,
            nNumIsotopicEndpoints: 2,
            num_iso_H: [1, 2, 3],
            tni: crate::source_types::TNI {
                nNumRemovedExplicitH: 4,
                nNumRemovedProtons: 5,
                nNumRemovedProtonsIsotopic: [6, 7, 8],
                bNormalizationFlags: u64::MAX,
            },
            bTautFlags: 0x1234,
            bTautFlagsDone: 0x5678,
        };
        assert_eq!(heap.live_allocation_count(), 4);
        assert_eq!(free_t_group_info(&mut heap, Some(&mut info)), Ok(0));
        assert_eq!(info, T_GROUP_INFO::default());
        assert_eq!(heap.live_allocation_count(), 0);

        let mut partial_heap = SourceHeap::default();
        let endpoint = partial_heap.allocate_model_storage(vec![10_u16]).unwrap();
        let mut partial = T_GROUP_INFO {
            nEndpointAtomNumber: endpoint,
            nNumEndpoints: 1,
            bTautFlags: u64::MAX,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            free_t_group_info(&mut partial_heap, Some(&mut partial)),
            Ok(0)
        );
        assert_eq!(partial, T_GROUP_INFO::default());
        assert_eq!(partial_heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__ichitaut__make_a_copy_of_t_group_info__line_6364() {
        fn source_info(heap: &mut SourceHeap) -> T_GROUP_INFO {
            let mut first = T_GROUP::default();
            first.nGroupNumber = 11;
            first.num = [1, 2, 3, 4, 5];
            let mut second = T_GROUP::default();
            second.nGroupNumber = 22;
            second.num = [6, 7, 8, 9, 10];
            T_GROUP_INFO {
                t_group: heap.allocate_model_storage(vec![first, second]).unwrap(),
                nEndpointAtomNumber: heap.allocate_model_storage(vec![3_u16, 5, 7]).unwrap(),
                tGroupNumber: heap
                    .allocate_model_storage(vec![10_u16, 11, 12, 13, 20, 21, 22, 23])
                    .unwrap(),
                nNumEndpoints: 3,
                num_t_groups: 2,
                max_num_t_groups: 2,
                bIgnoreIsotopic: 1,
                nIsotopicEndpointAtomNumber: heap.allocate_model_storage(vec![31_u16, 32]).unwrap(),
                nNumIsotopicEndpoints: 2,
                num_iso_H: [4, 5, 6],
                tni: crate::source_types::TNI {
                    nNumRemovedExplicitH: 7,
                    nNumRemovedProtons: 8,
                    nNumRemovedProtonsIsotopic: [9, 10, 11],
                    bNormalizationFlags: 0x1234,
                },
                bTautFlags: 0x5678,
                bTautFlagsDone: 0x9abc,
            }
        }

        let mut null_heap = SourceHeap::default();
        let mut source = source_info(&mut null_heap);
        let original_isotopic = source.nIsotopicEndpointAtomNumber;
        assert_eq!(
            make_a_copy_of_t_group_info(&mut null_heap, None, Some(&mut source)),
            Ok(0)
        );
        assert_eq!(source.nIsotopicEndpointAtomNumber, original_isotopic);
        assert_eq!(null_heap.live_allocation_count(), 4);

        let old_destination_group = null_heap
            .allocate_model_storage(vec![T_GROUP::default()])
            .unwrap();
        let mut destination = T_GROUP_INFO {
            t_group: old_destination_group,
            max_num_t_groups: 99,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            make_a_copy_of_t_group_info(&mut null_heap, Some(&mut destination), None),
            Ok(0)
        );
        assert_eq!(destination, T_GROUP_INFO::default());
        assert_eq!(null_heap.live_allocation_count(), 4);

        let mut heap = SourceHeap::default();
        let mut source = source_info(&mut heap);
        let old_source_isotopic = source.nIsotopicEndpointAtomNumber;
        let old_destination_endpoint = heap.allocate_model_storage(vec![99_u16, 98]).unwrap();
        let mut destination = T_GROUP_INFO {
            nEndpointAtomNumber: old_destination_endpoint,
            nNumEndpoints: 2,
            bTautFlags: u64::MAX,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            make_a_copy_of_t_group_info(&mut heap, Some(&mut destination), Some(&mut source)),
            Ok(0)
        );
        assert_ne!(source.nIsotopicEndpointAtomNumber, old_source_isotopic);
        assert_eq!(
            heap.slice(source.nIsotopicEndpointAtomNumber.as_const())
                .unwrap(),
            [31, 32]
        );
        assert_eq!(destination.nNumEndpoints, 3);
        assert_eq!(destination.num_t_groups, 2);
        assert_eq!(destination.max_num_t_groups, 2);
        assert_eq!(destination.bIgnoreIsotopic, 1);
        assert_eq!(destination.nNumIsotopicEndpoints, 2);
        assert_eq!(destination.tni, source.tni);
        assert_eq!(destination.num_iso_H, [0, 0, 0]);
        assert_eq!(destination.bTautFlags, 0x5678);
        assert_eq!(destination.bTautFlagsDone, 0x9abc);
        assert_eq!(
            heap.slice(destination.nEndpointAtomNumber.as_const())
                .unwrap(),
            [3, 5, 7]
        );
        assert_eq!(
            heap.slice(destination.tGroupNumber.as_const()).unwrap(),
            [10, 11, 12, 13, 20, 21, 22, 23]
        );
        let copied_groups = heap.slice(destination.t_group.as_const()).unwrap();
        assert_eq!(copied_groups.len(), 3);
        assert_eq!(copied_groups[0].nGroupNumber, 11);
        assert_eq!(copied_groups[1].nGroupNumber, 22);
        assert_eq!(
            heap.slice(destination.nIsotopicEndpointAtomNumber.as_const())
                .unwrap(),
            [31, 32]
        );
        assert_eq!(heap.live_allocation_count(), 8);

        for successful_allocations in 0..5 {
            let mut failure_heap = SourceHeap::default();
            let mut source = source_info(&mut failure_heap);
            let mut destination = T_GROUP_INFO {
                bTautFlags: u64::MAX,
                ..T_GROUP_INFO::default()
            };
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                make_a_copy_of_t_group_info(
                    &mut failure_heap,
                    Some(&mut destination),
                    Some(&mut source),
                ),
                Ok(1),
                "allocation failure after {successful_allocations} successes"
            );
            assert_eq!(
                destination,
                T_GROUP_INFO::default(),
                "destination cleanup after failure {successful_allocations}"
            );
            if successful_allocations == 3 {
                assert_eq!(failure_heap.live_allocation_count(), 3);
                assert_eq!(
                    failure_heap.slice(source.nIsotopicEndpointAtomNumber.as_const()),
                    Err(SourceHeapError::MissingAllocation)
                );
            } else {
                assert_eq!(
                    failure_heap.live_allocation_count(),
                    4,
                    "source ownership after failure {successful_allocations}"
                );
                assert_eq!(
                    failure_heap
                        .slice(source.nIsotopicEndpointAtomNumber.as_const())
                        .unwrap(),
                    [31, 32]
                );
            }
        }
    }

    #[test]
    fn source_port__ichitaut__set_tautomer_iso_sort_keys__line_6477() {
        let mut empty_heap = SourceHeap::default();
        assert_eq!(set_tautomer_iso_sort_keys(&mut empty_heap, None), Ok(0));
        let mut empty_info = T_GROUP_INFO::default();
        assert_eq!(
            set_tautomer_iso_sort_keys(&mut empty_heap, Some(&mut empty_info)),
            Ok(0)
        );

        let mut heap = SourceHeap::default();
        let groups = heap
            .allocate_model_storage(vec![
                T_GROUP {
                    num: [0, 0, 1, 2, 3],
                    iWeight: -1,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    num: [9, 8, 0, 0, 0],
                    iWeight: -2,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    num: [0, 0, u16::MAX, u16::MAX, u16::MAX],
                    iWeight: -3,
                    ..T_GROUP::default()
                },
            ])
            .unwrap();
        let mut info = T_GROUP_INFO {
            t_group: groups,
            num_t_groups: 3,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            set_tautomer_iso_sort_keys(&mut heap, Some(&mut info)),
            Ok(2)
        );
        let groups_after = heap.slice(groups.as_const()).unwrap();
        assert_eq!(groups_after[0].iWeight, 1_050_627);
        assert_eq!(groups_after[1].iWeight, 0);
        assert_eq!(
            groups_after[2].iWeight,
            i64::from(u16::MAX) * (1 + 1024 + 1_048_576)
        );

        let unchanged = groups_after.to_vec();
        info.nNumIsotopicEndpoints = 1;
        assert_eq!(
            set_tautomer_iso_sort_keys(&mut heap, Some(&mut info)),
            Ok(0)
        );
        assert_eq!(heap.slice(groups.as_const()).unwrap(), unchanged);

        info.nNumIsotopicEndpoints = 0;
        info.num_t_groups = 0;
        assert_eq!(
            set_tautomer_iso_sort_keys(&mut heap, Some(&mut info)),
            Ok(0)
        );
        assert_eq!(heap.slice(groups.as_const()).unwrap(), unchanged);

        info.num_t_groups = 4;
        assert_eq!(
            set_tautomer_iso_sort_keys(&mut heap, Some(&mut info)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichitaut__counttautomergroups__line_6519() {
        fn group(number: u16, endpoints: u16, hydrogens: u16, charges: u16) -> T_GROUP {
            let mut group = T_GROUP {
                nGroupNumber: number,
                nNumEndpoints: endpoints,
                ..T_GROUP::default()
            };
            group.num[0] = hydrogens;
            group.num[1] = charges;
            group
        }

        fn normal_fixture(heap: &mut SourceHeap) -> (SourceMutPointer<sp_ATOM>, T_GROUP_INFO) {
            let atoms = heap
                .allocate_model_storage(vec![
                    sp_ATOM {
                        endpoint: 2,
                        ..sp_ATOM::default()
                    },
                    sp_ATOM {
                        endpoint: 5,
                        ..sp_ATOM::default()
                    },
                    sp_ATOM {
                        endpoint: 2,
                        ..sp_ATOM::default()
                    },
                    sp_ATOM::default(),
                ])
                .unwrap();
            let groups = heap
                .allocate_model_storage(vec![group(2, 2, 4, 1), group(5, 1, 2, 0)])
                .unwrap();
            (
                atoms,
                T_GROUP_INFO {
                    t_group: groups,
                    num_t_groups: 2,
                    max_num_t_groups: 2,
                    ..T_GROUP_INFO::default()
                },
            )
        }

        let mut null_heap = SourceHeap::default();
        assert_eq!(
            CountTautomerGroups(&mut null_heap, SourceMutPointer::null(), 0, None),
            Ok(0)
        );
        let mut empty = T_GROUP_INFO::default();
        assert_eq!(
            CountTautomerGroups(
                &mut null_heap,
                SourceMutPointer::null(),
                0,
                Some(&mut empty)
            ),
            Ok(0)
        );

        let mut heap = SourceHeap::default();
        let (atoms, mut info) = normal_fixture(&mut heap);
        let orphaned_endpoints = heap.allocate_model_storage(vec![91_u16, 92]).unwrap();
        let orphaned_groups = heap
            .allocate_model_storage(vec![81_u16, 82, 83, 84])
            .unwrap();
        info.nEndpointAtomNumber = orphaned_endpoints;
        info.tGroupNumber = orphaned_groups;
        assert_eq!(
            CountTautomerGroups(&mut heap, atoms, 4, Some(&mut info)),
            Ok(10)
        );
        assert_eq!((info.num_t_groups, info.nNumEndpoints), (2, 3));
        assert_eq!(
            heap.slice(atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.endpoint)
                .collect::<Vec<_>>(),
            [1, 2, 1, 0]
        );
        let groups = heap.slice(info.t_group.as_const()).unwrap();
        assert_eq!(
            (
                groups[0].nGroupNumber,
                groups[0].nFirstEndpointAtNoPos,
                groups[0].num[0],
                groups[1].nGroupNumber,
                groups[1].nFirstEndpointAtNoPos,
                groups[1].num[0],
            ),
            (1, 0, 3, 2, 2, 2)
        );
        assert_eq!(
            heap.slice(info.nEndpointAtomNumber.as_const()).unwrap(),
            [0, 2, 1]
        );
        assert_eq!(
            heap.slice(info.tGroupNumber.as_const()).unwrap(),
            [0, 1, 0, 0, 0, 0, 0, 0]
        );
        assert_eq!(heap.slice(orphaned_endpoints.as_const()).unwrap(), [91, 92]);
        assert_eq!(
            heap.slice(orphaned_groups.as_const()).unwrap(),
            [81, 82, 83, 84]
        );
        assert_eq!(heap.live_allocation_count(), 6);

        let mut no_h_heap = SourceHeap::default();
        let no_h_atoms = no_h_heap
            .allocate_model_storage(vec![sp_ATOM {
                endpoint: 1,
                ..sp_ATOM::default()
            }])
            .unwrap();
        let no_h_groups = no_h_heap
            .allocate_model_storage(vec![group(1, 1, 1, 1)])
            .unwrap();
        let mut no_h_info = T_GROUP_INFO {
            t_group: no_h_groups,
            num_t_groups: 1,
            max_num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            CountTautomerGroups(&mut no_h_heap, no_h_atoms, 1, Some(&mut no_h_info)),
            Ok(0)
        );
        assert_eq!((no_h_info.num_t_groups, no_h_info.nNumEndpoints), (0, 0));
        assert_eq!(
            no_h_heap.slice(no_h_atoms.as_const()).unwrap()[0].endpoint,
            1
        );
        assert_eq!(no_h_heap.live_allocation_count(), 2);

        let mut no_endpoint_heap = SourceHeap::default();
        let no_endpoint_atoms = no_endpoint_heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        let no_endpoint_groups = no_endpoint_heap
            .allocate_model_storage(vec![group(1, 1, 2, 0)])
            .unwrap();
        let mut no_endpoint_info = T_GROUP_INFO {
            t_group: no_endpoint_groups,
            num_t_groups: 1,
            max_num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        no_endpoint_info.tni.bNormalizationFlags = u64::from(FLAG_NORM_CONSIDER_TAUT);
        assert_eq!(
            CountTautomerGroups(
                &mut no_endpoint_heap,
                no_endpoint_atoms,
                1,
                Some(&mut no_endpoint_info),
            ),
            Ok(1)
        );
        assert_eq!(
            (
                no_endpoint_info.num_t_groups,
                no_endpoint_info.nNumEndpoints
            ),
            (0, 0)
        );

        for invalid_endpoint in [3_u16, u16::MAX] {
            let mut invalid_heap = SourceHeap::default();
            let invalid_atoms = invalid_heap
                .allocate_model_storage(vec![sp_ATOM {
                    endpoint: invalid_endpoint,
                    ..sp_ATOM::default()
                }])
                .unwrap();
            let invalid_groups = invalid_heap
                .allocate_model_storage(vec![group(2, 1, 2, 0)])
                .unwrap();
            let mut invalid_info = T_GROUP_INFO {
                t_group: invalid_groups,
                num_t_groups: 1,
                max_num_t_groups: 1,
                ..T_GROUP_INFO::default()
            };
            assert_eq!(
                CountTautomerGroups(&mut invalid_heap, invalid_atoms, 1, Some(&mut invalid_info),),
                Ok(CT_TAUCOUNT_ERR)
            );
            assert_eq!(invalid_heap.live_allocation_count(), 2);
        }

        let mut mismatch_heap = SourceHeap::default();
        let mismatch_atoms = mismatch_heap
            .allocate_model_storage(vec![sp_ATOM {
                endpoint: 1,
                ..sp_ATOM::default()
            }])
            .unwrap();
        let mismatch_groups = mismatch_heap
            .allocate_model_storage(vec![group(1, 2, 3, 0)])
            .unwrap();
        let mut mismatch_info = T_GROUP_INFO {
            t_group: mismatch_groups,
            num_t_groups: 1,
            max_num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        assert_eq!(
            CountTautomerGroups(
                &mut mismatch_heap,
                mismatch_atoms,
                1,
                Some(&mut mismatch_info),
            ),
            Ok(CT_TAUCOUNT_ERR)
        );

        for successful_allocations in 0..4 {
            let mut failure_heap = SourceHeap::default();
            let (atoms, mut info) = normal_fixture(&mut failure_heap);
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                CountTautomerGroups(&mut failure_heap, atoms, 4, Some(&mut info)),
                Ok(CT_TAUCOUNT_ERR),
                "allocation failure after {successful_allocations} successes"
            );
            assert_eq!(
                failure_heap.source_allocation_calls(),
                successful_allocations + 1
            );
            assert_eq!(failure_heap.live_allocation_count(), 2);
            assert_eq!((info.num_t_groups, info.nNumEndpoints), (0, 0));
        }
    }
}
