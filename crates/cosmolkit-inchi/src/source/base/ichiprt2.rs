use crate::source::base::ichi_io::inchi_strbuf_printf;
use crate::source::base::ichimake::{GetDfsOrder4CT, mystrrev};
use crate::source::base::ichiparm::source_strtod_with_end;
use crate::source::base::util::inchi_free;
use crate::source_types::{
    CT_MODE_ABC_NUM_CLOSURES, CT_MODE_ABC_NUMBERS, CT_MODE_PREDECESSORS, EQL_EQU_ISO, EQL_EQU_TG,
    EQL_EXISTS, EQL_NUM, EQL_NUM_INV, EQL_NUM_ISO, EQL_SP2, EQL_SP3, EQL_SP3_INV, INCHI_IOS_STRING,
    INChI, INChI_Stereo, MAX_ATOMS, ORIG_INFO, SourceConstPointer, SourceFormatArgument,
    SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList,
};

#[allow(non_snake_case)]
pub(crate) fn Eql_INChI_Stereo(
    heap: &SourceHeap,
    stereo1: Option<&INChI_Stereo>,
    equal1: i32,
    stereo2: Option<&INChI_Stereo>,
    equal2: i32,
    relative_or_racemic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:54 Eql_INChI_Stereo
    // INCHI✔️❌: int Eql_INChI_Stereo( INChI_Stereo  *s1,
    // INCHI✔️❌:                       int           eql1,
    // INCHI✔️❌:                       INChI_Stereo  *s2,
    // INCHI✔️❌:                       int           eql2,
    // INCHI✔️❌:                       int           bRelRac )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int inv1 = 0, inv2 = 0, len;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!s1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌: #else
    // INCHI✔️❌:     bRelRac = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (EQL_SP2 == eql1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( len = s1->nNumberOfStereoBonds ) > 0 && s1->b_parity && s1->nBondAtom1 && s1->nBondAtom2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!s2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (EQL_EXISTS == eql2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* find whether double bond stereo exists*/
    // INCHI✔️❌:                     return 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (EQL_SP2 == eql2 &&
    // INCHI✔️❌:                  len == s2->nNumberOfStereoBonds && s2->b_parity && s2->nBondAtom1 && s2->nBondAtom2 &&
    // INCHI✔️❌:                  !memcmp( s1->nBondAtom1, s2->nBondAtom1, len * sizeof( s1->nBondAtom1[0] ) ) &&
    // INCHI✔️❌:                  !memcmp( s1->nBondAtom2, s2->nBondAtom2, len * sizeof( s1->nBondAtom2[0] ) ) &&
    // INCHI✔️❌:                  !memcmp( s1->b_parity, s2->b_parity, len * sizeof( s1->b_parity[0] ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( eql1 == EQL_SP3 || ( inv1 = ( eql1 == EQL_SP3_INV ) ) ) &&
    // INCHI✔️❌:             ( len = s1->nNumberOfStereoCenters ) > ( bRelRac ? 1 : 0 )) /* djb-rwth: addressing coverity ID #499484 -- bRelRac does not have to be 0 */
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌:             S_CHAR  *t_parity1, *t_parity2;
    // INCHI✔️❌:             AT_NUMB *nNumber1, *nNumber2;
    // INCHI✔️❌:             if (inv1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (s1->nCompInv2Abs)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     t_parity1 = s1->t_parityInv;
    // INCHI✔️❌:                     nNumber1 = s1->nNumberInv;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 t_parity1 = s1->t_parity;
    // INCHI✔️❌:                 nNumber1 = s1->nNumber;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (t_parity1 && nNumber1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!s2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (EQL_EXISTS == eql2 && ( !inv1 || s1->nCompInv2Abs ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* the 1st sp3 (inverted if requested) stereo exists*/
    // INCHI✔️❌:                         return 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     return 0;  /* both sp3 do not exist */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (( eql2 == EQL_SP3 || ( inv2 = ( eql2 == EQL_SP3_INV ) ) ) &&
    // INCHI✔️❌:                     len == s2->nNumberOfStereoCenters)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (inv2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (s2->nCompInv2Abs && s1->nCompInv2Abs)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             t_parity2 = s2->t_parityInv;
    // INCHI✔️❌:                             nNumber2 = s2->nNumberInv;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* if one sp3 is inverted then another should have non-trivial inverted stereo */
    // INCHI✔️❌:                             return 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (inv1 && !s2->nCompInv2Abs)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* if one sp3 is inverted then another should have non-trivial inverted stereo */
    // INCHI✔️❌:                             return 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         t_parity2 = s2->t_parity;
    // INCHI✔️❌:                         nNumber2 = s2->nNumber;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (t_parity2 && nNumber2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (inv1 ^ inv2)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int i, num_inv;
    // INCHI✔️❌:                             for (i = 0, num_inv = 0; i < len; i++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (nNumber1[i] != nNumber2[i])
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 if (ATOM_PARITY_WELL_DEF( t_parity1[i] ) &&
    // INCHI✔️❌:                                      ATOM_PARITY_WELL_DEF( t_parity2[i] ))
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (3 == t_parity1[i] + t_parity2[i])
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         num_inv++;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         break;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (t_parity1[i] != t_parity2[i])
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         break;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             return ( len == i && num_inv > 0 );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return !memcmp( t_parity1, t_parity2, len * sizeof( t_parity1[0] ) ) &&
    // INCHI✔️❌:                                 !memcmp( nNumber1, nNumber2, len * sizeof( nNumber1[0] ) );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Eql_INChI_Stereo

    let Some(stereo1) = stereo1 else {
        return Ok(0);
    };
    if equal1 == EQL_SP2 as i32 {
        let length = stereo1.nNumberOfStereoBonds;
        if length > 0
            && !stereo1.b_parity.is_null()
            && !stereo1.nBondAtom1.is_null()
            && !stereo1.nBondAtom2.is_null()
        {
            let Some(stereo2) = stereo2 else {
                return Ok(i32::from(equal2 == EQL_EXISTS as i32));
            };
            if equal2 == EQL_SP2 as i32
                && length == stereo2.nNumberOfStereoBonds
                && !stereo2.b_parity.is_null()
                && !stereo2.nBondAtom1.is_null()
                && !stereo2.nBondAtom2.is_null()
            {
                let length =
                    usize::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let atoms11 = heap
                    .slice(stereo1.nBondAtom1.as_const())?
                    .get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atoms12 = heap
                    .slice(stereo2.nBondAtom1.as_const())?
                    .get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atoms21 = heap
                    .slice(stereo1.nBondAtom2.as_const())?
                    .get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atoms22 = heap
                    .slice(stereo2.nBondAtom2.as_const())?
                    .get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let parity1 = heap
                    .slice(stereo1.b_parity.as_const())?
                    .get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let parity2 = heap
                    .slice(stereo2.b_parity.as_const())?
                    .get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                return Ok(i32::from(
                    atoms11 == atoms12 && atoms21 == atoms22 && parity1 == parity2,
                ));
            }
        }
        return Ok(0);
    }

    let inverted1 = equal1 == EQL_SP3_INV as i32;
    if !(equal1 == EQL_SP3 as i32 || inverted1)
        || stereo1.nNumberOfStereoCenters <= i32::from(relative_or_racemic != 0)
    {
        return Ok(0);
    }
    let length = usize::try_from(stereo1.nNumberOfStereoCenters)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let (parity1_pointer, number1_pointer) = if inverted1 {
        if stereo1.nCompInv2Abs == 0 {
            return Ok(0);
        }
        (stereo1.t_parityInv, stereo1.nNumberInv)
    } else {
        (stereo1.t_parity, stereo1.nNumber)
    };
    if parity1_pointer.is_null() || number1_pointer.is_null() {
        return Ok(0);
    }
    let Some(stereo2) = stereo2 else {
        return Ok(i32::from(
            equal2 == EQL_EXISTS as i32 && (!inverted1 || stereo1.nCompInv2Abs != 0),
        ));
    };
    let inverted2 = equal2 == EQL_SP3_INV as i32;
    if !(equal2 == EQL_SP3 as i32 || inverted2)
        || stereo1.nNumberOfStereoCenters != stereo2.nNumberOfStereoCenters
    {
        return Ok(0);
    }
    let (parity2_pointer, number2_pointer) = if inverted2 {
        if stereo2.nCompInv2Abs == 0 || stereo1.nCompInv2Abs == 0 {
            return Ok(0);
        }
        (stereo2.t_parityInv, stereo2.nNumberInv)
    } else {
        if inverted1 && stereo2.nCompInv2Abs == 0 {
            return Ok(0);
        }
        (stereo2.t_parity, stereo2.nNumber)
    };
    if parity2_pointer.is_null() || number2_pointer.is_null() {
        return Ok(0);
    }
    let parity1 = heap
        .slice(parity1_pointer.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let parity2 = heap
        .slice(parity2_pointer.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let numbers1 = heap
        .slice(number1_pointer.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let numbers2 = heap
        .slice(number2_pointer.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if inverted1 ^ inverted2 {
        let mut inversions = 0;
        for index in 0..length {
            if numbers1[index] != numbers2[index] {
                return Ok(0);
            }
            let parity1_well_defined = (1..=2).contains(&parity1[index]);
            let parity2_well_defined = (1..=2).contains(&parity2[index]);
            if parity1_well_defined && parity2_well_defined {
                if parity1[index] + parity2[index] != 3 {
                    return Ok(0);
                }
                inversions += 1;
            } else if parity1[index] != parity2[index] {
                return Ok(0);
            }
        }
        return Ok(i32::from(inversions > 0));
    }
    Ok(i32::from(parity1 == parity2 && numbers1 == numbers2))
}

#[allow(non_snake_case)]
pub(crate) fn Eql_INChI_Isotopic(
    heap: &SourceHeap,
    inchi1: Option<&INChI>,
    inchi2: Option<&INChI>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:203 Eql_INChI_Isotopic
    // INCHI✔️✔️: int Eql_INChI_Isotopic( INChI *i1, INChI *i2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int eq = i1
    // INCHI✔️✔️:             && i2
    // INCHI✔️✔️:             && !i1->bDeleted
    // INCHI✔️✔️:             && !i2->bDeleted
    // INCHI✔️✔️:             && ( i1->nNumberOfIsotopicAtoms > 0 || i1->nNumberOfIsotopicTGroups > 0 )
    // INCHI✔️✔️:             && i1->nNumberOfIsotopicAtoms == i2->nNumberOfIsotopicAtoms
    // INCHI✔️✔️:             && i1->nNumberOfIsotopicTGroups == i2->nNumberOfIsotopicTGroups
    // INCHI✔️✔️:             && ( !i1->nNumberOfIsotopicAtoms ||
    // INCHI✔️✔️:                 (i1->IsotopicAtom && i2->IsotopicAtom &&
    // INCHI✔️✔️:                 !memcmp( i1->IsotopicAtom, i2->IsotopicAtom,
    // INCHI✔️✔️:                         i1->nNumberOfIsotopicAtoms * sizeof( i1->IsotopicAtom[0] ) )) )
    // INCHI✔️✔️:             && ( !i1->nNumberOfIsotopicTGroups ||
    // INCHI✔️✔️:                 (i1->IsotopicTGroup && i2->IsotopicTGroup &&
    // INCHI✔️✔️:                 !memcmp( i1->IsotopicTGroup, i2->IsotopicTGroup,
    // INCHI✔️✔️:                     i1->nNumberOfIsotopicTGroups * sizeof( i1->IsotopicAtom[0] )) )
    // INCHI✔️✔️:             ); /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return eq;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: Eql_INChI_Isotopic
    let (Some(inchi1), Some(inchi2)) = (inchi1, inchi2) else {
        return Ok(0);
    };
    if inchi1.bDeleted != 0
        || inchi2.bDeleted != 0
        || (inchi1.nNumberOfIsotopicAtoms <= 0 && inchi1.nNumberOfIsotopicTGroups <= 0)
        || inchi1.nNumberOfIsotopicAtoms != inchi2.nNumberOfIsotopicAtoms
        || inchi1.nNumberOfIsotopicTGroups != inchi2.nNumberOfIsotopicTGroups
    {
        return Ok(0);
    }
    if inchi1.nNumberOfIsotopicAtoms != 0 {
        if inchi1.IsotopicAtom.is_null() || inchi2.IsotopicAtom.is_null() {
            return Ok(0);
        }
        let count = usize::try_from(inchi1.nNumberOfIsotopicAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atoms1 = heap
            .slice(inchi1.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms2 = heap
            .slice(inchi2.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atoms1 != atoms2 {
            return Ok(0);
        }
    }
    if inchi1.nNumberOfIsotopicTGroups != 0 {
        if inchi1.IsotopicTGroup.is_null() || inchi2.IsotopicTGroup.is_null() {
            return Ok(0);
        }
        let count = usize::try_from(inchi1.nNumberOfIsotopicTGroups)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let groups1 = heap
            .slice(inchi1.IsotopicTGroup.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let groups2 = heap
            .slice(inchi2.IsotopicTGroup.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if groups1 != groups2 {
            return Ok(0);
        }
    }
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn Eql_INChI_Aux_Equ(
    heap: &SourceHeap,
    aux1: Option<&crate::source_types::INChI_Aux>,
    equal1: i32,
    aux2: Option<&crate::source_types::INChI_Aux>,
    equal2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:227 Eql_INChI_Aux_Equ
    // INCHI✔️✔️: int Eql_INChI_Aux_Equ( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int t1 = 0, t2 = 0, len = 0;
    // INCHI✔️✔️:     AT_NUMB *n1 = NULL, *n2 = NULL;
    // INCHI✔️✔️:     if (!a1 || !a2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     t1 = ( eql1 & EQL_EQU_TG );
    // INCHI✔️✔️:     t2 = ( eql2 & EQL_EQU_TG );
    // INCHI✔️✔️:     if (t1 && t2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (( len = a1->nNumberOfTGroups ) > 0 && len == a2->nNumberOfTGroups && !a1->bDeleted && !a2->bDeleted)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (eql1 & EQL_EQU_ISO)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (a1->bIsIsotopic)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     n1 = a1->nConstitEquIsotopicTGroupNumbers;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 n1 = a1->nConstitEquTGroupNumbers;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (eql2 & EQL_EQU_ISO)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (a2->bIsIsotopic)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     n2 = a2->nConstitEquIsotopicTGroupNumbers;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 n2 = a2->nConstitEquTGroupNumbers;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:         if (!t1 && !t2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (( len = a1->nNumberOfAtoms ) > 0 && len == a2->nNumberOfAtoms && !a1->bDeleted && !a2->bDeleted)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (eql1 & EQL_EQU_ISO)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (a1->bIsIsotopic)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         n1 = a1->nConstitEquIsotopicNumbers;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     n1 = a1->nConstitEquNumbers;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (eql2 & EQL_EQU_ISO)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (a2->bIsIsotopic)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         n2 = a2->nConstitEquIsotopicNumbers;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     n2 = a2->nConstitEquNumbers;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     if (n1 && n2 && !memcmp( n1, n2, len * sizeof( n1[0] ) ) && bHasEquString( n1, len ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: Eql_INChI_Aux_Equ

    let (Some(aux1), Some(aux2)) = (aux1, aux2) else {
        return Ok(0);
    };
    let tgroup1 = equal1 & EQL_EQU_TG as i32;
    let tgroup2 = equal2 & EQL_EQU_TG as i32;
    let (length, numbers1, numbers2) = if tgroup1 != 0 && tgroup2 != 0 {
        if aux1.nNumberOfTGroups <= 0
            || aux1.nNumberOfTGroups != aux2.nNumberOfTGroups
            || aux1.bDeleted != 0
            || aux2.bDeleted != 0
        {
            return Ok(0);
        }
        let first = if equal1 & EQL_EQU_ISO as i32 != 0 {
            if aux1.bIsIsotopic == 0 {
                return Ok(0);
            }
            aux1.nConstitEquIsotopicTGroupNumbers
        } else {
            aux1.nConstitEquTGroupNumbers
        };
        let second = if equal2 & EQL_EQU_ISO as i32 != 0 {
            if aux2.bIsIsotopic == 0 {
                return Ok(0);
            }
            aux2.nConstitEquIsotopicTGroupNumbers
        } else {
            aux2.nConstitEquTGroupNumbers
        };
        (aux1.nNumberOfTGroups, first, second)
    } else if tgroup1 == 0 && tgroup2 == 0 {
        if aux1.nNumberOfAtoms <= 0
            || aux1.nNumberOfAtoms != aux2.nNumberOfAtoms
            || aux1.bDeleted != 0
            || aux2.bDeleted != 0
        {
            return Ok(0);
        }
        let first = if equal1 & EQL_EQU_ISO as i32 != 0 {
            if aux1.bIsIsotopic == 0 {
                return Ok(0);
            }
            aux1.nConstitEquIsotopicNumbers
        } else {
            aux1.nConstitEquNumbers
        };
        let second = if equal2 & EQL_EQU_ISO as i32 != 0 {
            if aux2.bIsIsotopic == 0 {
                return Ok(0);
            }
            aux2.nConstitEquIsotopicNumbers
        } else {
            aux2.nConstitEquNumbers
        };
        (aux1.nNumberOfAtoms, first, second)
    } else {
        return Ok(0);
    };
    if numbers1.is_null() || numbers2.is_null() {
        return Ok(0);
    }
    let count = usize::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let first = heap
        .slice(numbers1.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = heap
        .slice(numbers2.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if first != second {
        return Ok(0);
    }
    bHasEquString(heap, numbers1.as_const(), length)
}

#[allow(non_snake_case)]
pub(crate) fn Eql_INChI_Aux_Num(
    heap: &SourceHeap,
    aux1: Option<&crate::source_types::INChI_Aux>,
    equal1: i32,
    aux2: Option<&crate::source_types::INChI_Aux>,
    equal2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:305 Eql_INChI_Aux_Num
    // INCHI✔️✔️: int Eql_INChI_Aux_Num( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int len;
    // INCHI✔️✔️:     AT_NUMB *n1 = NULL, *n2 = NULL;
    // INCHI✔️✔️:     if (!a1 || !a2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (( len = a1->nNumberOfAtoms ) <= 0 || len != a2->nNumberOfAtoms || a1->bDeleted || a2->bDeleted)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if ((( eql1 & EQL_NUM_ISO ) && !a1->bIsIsotopic) ||
    // INCHI✔️✔️:         (( eql2 & EQL_NUM_ISO ) && !a2->bIsIsotopic)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     switch (eql1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         case EQL_NUM:
    // INCHI✔️✔️:             n1 = a1->nOrigAtNosInCanonOrd;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case EQL_NUM_ISO:
    // INCHI✔️✔️:             n1 = a1->nIsotopicOrigAtNosInCanonOrd;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case EQL_NUM_INV:
    // INCHI✔️✔️:             n1 = a1->nOrigAtNosInCanonOrdInv;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case ( EQL_NUM_INV | EQL_NUM_ISO ):
    // INCHI✔️✔️:             n1 = a1->nIsotopicOrigAtNosInCanonOrdInv;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     switch (eql2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         case EQL_NUM:
    // INCHI✔️✔️:             n2 = a2->nOrigAtNosInCanonOrd;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case EQL_NUM_ISO:
    // INCHI✔️✔️:             n2 = a2->nIsotopicOrigAtNosInCanonOrd;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case EQL_NUM_INV:
    // INCHI✔️✔️:             n2 = a2->nOrigAtNosInCanonOrdInv;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case ( EQL_NUM_INV | EQL_NUM_ISO ):
    // INCHI✔️✔️:             n2 = a2->nIsotopicOrigAtNosInCanonOrdInv;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (n1 && n2 && !memcmp( n1, n2, len * sizeof( n1[0] ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: Eql_INChI_Aux_Num

    let (Some(aux1), Some(aux2)) = (aux1, aux2) else {
        return Ok(0);
    };
    if aux1.nNumberOfAtoms <= 0
        || aux1.nNumberOfAtoms != aux2.nNumberOfAtoms
        || aux1.bDeleted != 0
        || aux2.bDeleted != 0
    {
        return Ok(0);
    }
    if ((equal1 & EQL_NUM_ISO as i32) != 0 && aux1.bIsIsotopic == 0)
        || ((equal2 & EQL_NUM_ISO as i32) != 0 && aux2.bIsIsotopic == 0)
    {
        return Ok(0);
    }
    let number_pointer = |aux: &crate::source_types::INChI_Aux,
                          equal: i32|
     -> Option<SourceMutPointer<crate::source_types::AT_NUMB>> {
        match equal {
            value if value == EQL_NUM as i32 => Some(aux.nOrigAtNosInCanonOrd),
            value if value == EQL_NUM_ISO as i32 => Some(aux.nIsotopicOrigAtNosInCanonOrd),
            value if value == EQL_NUM_INV as i32 => Some(aux.nOrigAtNosInCanonOrdInv),
            value if value == (EQL_NUM_INV | EQL_NUM_ISO) as i32 => {
                Some(aux.nIsotopicOrigAtNosInCanonOrdInv)
            }
            _ => None,
        }
    };
    let (Some(numbers1), Some(numbers2)) =
        (number_pointer(aux1, equal1), number_pointer(aux2, equal2))
    else {
        return Ok(0);
    };
    if numbers1.is_null() || numbers2.is_null() {
        return Ok(0);
    }
    let length =
        usize::try_from(aux1.nNumberOfAtoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let numbers1 = heap
        .slice(numbers1.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let numbers2 = heap
        .slice(numbers2.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(numbers1 == numbers2))
}

#[allow(non_snake_case)]
pub(crate) fn MakeCtStringNew(
    heap: &mut SourceHeap,
    linear_ct: SourceMutPointer<u16>,
    length_ct: i32,
    add_delimiter: i32,
    number_hydrogens: SourceConstPointer<i8>,
    number_of_atoms: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:551 MakeCtStringNew
    // INCHI✔️❌: int MakeCtStringNew( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                      AT_NUMB          *LinearCT,
    // INCHI✔️❌:                      int              nLenCT,
    // INCHI✔️❌:                      int              bAddDelim,
    // INCHI✔️❌:                      S_CHAR           *nNum_H,
    // INCHI✔️❌:                      int              num_atoms,
    // INCHI✔️❌:                      INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                      int              nCtMode,
    // INCHI✔️❌:                      int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, i, bOvfl = *bOverflow;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     int   nValue, nDelim, num_H;
    // INCHI✔️❌:     AT_NUMB *nDfsOrderCT = NULL;
    // INCHI✔️❌:     int      bNoNum_H = ( NULL == nNum_H );
    // INCHI✔️❌:     int      nNumRingClosures;
    // INCHI✔️❌:     int      bAbcNumbers = ( 0 != ( nCtMode & CT_MODE_ABC_NUMBERS ) );
    // INCHI✔️❌:     int      bPredecessors = ( 0 != ( nCtMode & CT_MODE_PREDECESSORS ) );
    // INCHI✔️❌:     int      bCountRingClosures = bAbcNumbers && bPredecessors && ( nCtMode & CT_MODE_ABC_NUM_CLOSURES );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nLenCT <= 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;  /*  no atoms or a single atom: no connection table */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  make array containing connection string data */
    // INCHI✔️❌:     if (!( nDfsOrderCT = GetDfsOrder4CT( pCG, LinearCT, nLenCT,
    // INCHI✔️❌:         nNum_H, num_atoms, nCtMode ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ( *bOverflow )++;
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add connection table string */
    // INCHI✔️❌:     if (!bOvfl && bAddDelim)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "," );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!bOvfl)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nNumRingClosures = 0;
    // INCHI✔️❌:         for (i = 0; nDfsOrderCT[i]/* && nLen < nLen_szLinearCT*/; i += 3)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nValue = ( nDfsOrderCT[i] > MAX_ATOMS ) ? 0 : nDfsOrderCT[i];
    // INCHI✔️❌:             num_H = nDfsOrderCT[i + 1] ? nDfsOrderCT[i + 1] - 16 : 0;
    // INCHI✔️❌:             nDelim = nDfsOrderCT[i + 2];
    // INCHI✔️❌:             len = 0;
    // INCHI✔️❌:             /*  delimiter */
    // INCHI✔️❌:             if (bPredecessors)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bCountRingClosures)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (nDelim == '-' && i > 3 && bNoNum_H)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (!nNumRingClosures)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int j;
    // INCHI✔️❌:                             for (j = i; nDfsOrderCT[j] && '-' == nDfsOrderCT[j + 2]; j += 3)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 nNumRingClosures++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (nNumRingClosures)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, nNumRingClosures );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             nNumRingClosures--;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nNumRingClosures--;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nNumRingClosures = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (nDelim && !( bAbcNumbers && nDelim == ',' ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (nNum_H || i > 3)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         szValue[len++] = nDelim;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nDelim && !( bAbcNumbers && nDelim == '-' ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     szValue[len++] = nDelim;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (bAbcNumbers)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nValue || i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* the 1st value may be zero in case of presdecessor list */
    // INCHI✔️❌:                     len += MakeAbcNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, nValue );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (num_H)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, num_H );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nValue || i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* the 1st value may be zero in case of presdecessor list */
    // INCHI✔️❌:                     len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, nValue );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (num_H)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     szValue[len] = 'H';
    // INCHI✔️❌:                     len++;
    // INCHI✔️❌:                     if (num_H > 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, num_H );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (len > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:     if (nDfsOrderCT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nDfsOrderCT );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeCtStringNew

    if length_ct <= 1 {
        return Ok(0);
    }
    let initial_used_length = string_buffer.nUsedLength;
    let initial_overflow = *overflow;
    let dfs = GetDfsOrder4CT(
        heap,
        linear_ct,
        length_ct,
        number_hydrogens,
        number_of_atoms,
        ct_mode,
    )?;
    if dfs.is_null() {
        *overflow = overflow.wrapping_add(1);
        return Ok(0);
    }

    let computation = (|| -> Result<i32, SourceHeapError> {
        let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
        let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
        let string_format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
        if initial_overflow == 0 && add_delimiter != 0 {
            match inchi_strbuf_printf(
                heap,
                Some(string_buffer),
                comma.as_const(),
                &SourceVaList::default(),
            ) {
                Ok(_) | Err(SourceHeapError::AllocationFailed) => {}
                Err(error) => return Err(error),
            }
        }
        if initial_overflow == 0 {
            let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0;
            let predecessors = ct_mode & CT_MODE_PREDECESSORS as i32 != 0;
            let count_ring_closures =
                abc_numbers && predecessors && ct_mode & CT_MODE_ABC_NUM_CLOSURES as i32 != 0;
            let no_hydrogens = number_hydrogens.is_null();
            let mut remaining_ring_closures = 0_i32;
            let mut index = 0_usize;
            while heap.slice(dfs.as_const())?[index] != 0 {
                let value = {
                    let raw = heap.slice(dfs.as_const())?[index];
                    if u32::from(raw) > MAX_ATOMS {
                        0
                    } else {
                        i32::from(raw)
                    }
                };
                let hydrogens = {
                    let raw = heap.slice(dfs.as_const())?[index + 1];
                    if raw != 0 { i32::from(raw) - 16 } else { 0 }
                };
                let delimiter = heap.slice(dfs.as_const())?[index + 2];
                let mut length = 0_i32;

                if predecessors {
                    if count_ring_closures {
                        if delimiter == b'-' as u16 && index > 3 && no_hydrogens {
                            if remaining_ring_closures == 0 {
                                let mut scan = index;
                                while heap.slice(dfs.as_const())?[scan] != 0
                                    && heap.slice(dfs.as_const())?[scan + 2] == b'-' as u16
                                {
                                    remaining_ring_closures =
                                        remaining_ring_closures.wrapping_add(1);
                                    scan += 3;
                                }
                                if remaining_ring_closures != 0 {
                                    length += MakeDecNumber(
                                        heap,
                                        value_buffer.offset(i64::from(length))?,
                                        2048 - length,
                                        SourceConstPointer::null(),
                                        remaining_ring_closures,
                                    )?;
                                }
                            }
                            remaining_ring_closures = remaining_ring_closures.wrapping_sub(1);
                        } else {
                            remaining_ring_closures = 0;
                        }
                    } else if delimiter != 0 && !(abc_numbers && delimiter == b',' as u16) {
                        if !number_hydrogens.is_null() || index > 3 {
                            heap.slice_mut(value_buffer)?[length as usize] = delimiter as i8;
                            length += 1;
                        }
                    }
                } else if delimiter != 0 && !(abc_numbers && delimiter == b'-' as u16) {
                    heap.slice_mut(value_buffer)?[length as usize] = delimiter as i8;
                    length += 1;
                }

                if abc_numbers {
                    if value != 0 || index != 0 {
                        length += MakeAbcNumber(
                            heap,
                            value_buffer.offset(i64::from(length))?,
                            2048 - length,
                            SourceConstPointer::null(),
                            value,
                        )?;
                    }
                    if hydrogens != 0 {
                        length += MakeDecNumber(
                            heap,
                            value_buffer.offset(i64::from(length))?,
                            2048 - length,
                            SourceConstPointer::null(),
                            hydrogens,
                        )?;
                    }
                } else {
                    if value != 0 || index != 0 {
                        length += MakeDecNumber(
                            heap,
                            value_buffer.offset(i64::from(length))?,
                            2048 - length,
                            SourceConstPointer::null(),
                            value,
                        )?;
                    }
                    if hydrogens != 0 {
                        heap.slice_mut(value_buffer)?[length as usize] = b'H' as i8;
                        length += 1;
                        if hydrogens > 1 {
                            length += MakeDecNumber(
                                heap,
                                value_buffer.offset(i64::from(length))?,
                                2048 - length,
                                SourceConstPointer::null(),
                                hydrogens,
                            )?;
                        } else {
                            heap.slice_mut(value_buffer)?[length as usize] = 0;
                        }
                    }
                }

                if length > 0 {
                    match inchi_strbuf_printf(
                        heap,
                        Some(string_buffer),
                        string_format.as_const(),
                        &SourceVaList {
                            arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                            ..SourceVaList::default()
                        },
                    ) {
                        Ok(_) | Err(SourceHeapError::AllocationFailed) => {}
                        Err(error) => return Err(error),
                    }
                }
                index += 3;
            }
        }
        *overflow |= initial_overflow;
        Ok(string_buffer.nUsedLength - initial_used_length)
    })();
    inchi_free(heap, dfs)?;
    computation
}

#[allow(non_snake_case)]
pub(crate) fn MakeCtStringOld(
    heap: &mut SourceHeap,
    linear_ct: SourceConstPointer<u16>,
    length_ct: i32,
    add_delimiter: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:703 MakeCtStringOld
    // INCHI✔️❌: int MakeCtStringOld( AT_NUMB          *LinearCT,
    // INCHI✔️❌:                      int              nLenCT,
    // INCHI✔️❌:                      int              bAddDelim,
    // INCHI✔️❌:                      INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                      int              nCtMode,
    // INCHI✔️❌:                      int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, i, bLessThanPrev, bOvfl = *bOverflow;
    // INCHI✔️❌:     AT_NUMB nMax = 0;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     int   nValue, bNext = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add connection table string */
    // INCHI✔️❌:     if (!( nCtMode & CT_MODE_ABC_NUMBERS ) && !bOvfl && bAddDelim)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "," );
    // INCHI✔️❌:                                         /*
    // INCHI✔️❌:                                         if ( nLen_szLinearCT > 1 ) {
    // INCHI✔️❌:                                             strcpy( szLinearCT, "," );
    // INCHI✔️❌:                                             nLen ++;
    // INCHI✔️❌:                                         } else {
    // INCHI✔️❌:                                             bOvfl = 1;
    // INCHI✔️❌:                                         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!bOvfl)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < nLenCT/* && nLen < nLen_szLinearCT*/; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             bLessThanPrev = 0;
    // INCHI✔️❌:             if (!( nCtMode & CT_MODE_NO_ORPHANS ) || ( ( bLessThanPrev = LinearCT[i] < nMax ) ||
    // INCHI✔️❌:                 (i + 1 < nLenCT && LinearCT[i + 1] < ( nMax = LinearCT[i] )) )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nValue = LinearCT[i];
    // INCHI✔️❌:                 if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     len = MakeAbcNumber( szValue, ( int )sizeof( szValue ), ( !bNext && bAddDelim ) ? ITEM_DELIMETER : NULL, nValue );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (nCtMode & CT_MODE_NO_ORPHANS)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*  censored CT */
    // INCHI✔️❌:                         /*  output '-' as a delimiter to show a bonding for decimal output of the connection table */
    // INCHI✔️❌:                         len = MakeDecNumber( szValue, ( int )sizeof( szValue ), bLessThanPrev ? "-" : ITEM_DELIMETER, nValue );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = MakeDecNumber( szValue, ( int )sizeof( szValue ), i ? ITEM_DELIMETER : NULL, nValue );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (len > 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:                     bNext++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 /*if ( 0 <= len && nLen+len < nLen_szLinearCT) {
    // INCHI✔️❌:                     if ( len ) {
    // INCHI✔️❌:                         strcpy( szLinearCT+nLen, szValue );
    // INCHI✔️❌:                         nLen += len;
    // INCHI✔️❌:                         bNext ++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     } else {
    // INCHI✔️❌:                         bOvfl = 1;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }*/
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeCtStringOld

    let initial_length = string_buffer.nUsedLength;
    let local_overflow = *overflow;
    let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
    let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
    let hyphen = heap.allocate_model_storage(vec![b'-' as i8, 0])?;
    let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
    let append = |heap: &mut SourceHeap,
                  string_buffer: &mut INCHI_IOS_STRING,
                  value: SourceConstPointer<i8>|
     -> Result<(), SourceHeapError> {
        match inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(value)],
                ..SourceVaList::default()
            },
        ) {
            Ok(_) | Err(SourceHeapError::AllocationFailed) => Ok(()),
            Err(error) => Err(error),
        }
    };
    let abc = ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0;
    let no_orphans = ct_mode & crate::source_types::CT_MODE_NO_ORPHANS as i32 != 0;
    if !abc && local_overflow == 0 && add_delimiter != 0 {
        append(heap, string_buffer, comma.as_const())?;
    }
    if local_overflow == 0 {
        let length = if length_ct > 0 {
            usize::try_from(length_ct).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        } else {
            0
        };
        let values = if length == 0 {
            Vec::new()
        } else {
            heap.slice(linear_ct)?
                .get(..length)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        };
        let mut maximum = 0_u16;
        let mut next = 0_i32;
        for index in 0..length {
            let less_than_previous = values[index] < maximum;
            let include = if !no_orphans || less_than_previous {
                true
            } else if index + 1 < length {
                maximum = values[index];
                values[index + 1] < maximum
            } else {
                false
            };
            if !include {
                continue;
            }
            let value = i32::from(values[index]);
            let written = if abc {
                MakeAbcNumber(
                    heap,
                    value_buffer,
                    2048,
                    if next == 0 && add_delimiter != 0 {
                        comma.as_const()
                    } else {
                        SourceConstPointer::null()
                    },
                    value,
                )?
            } else if no_orphans {
                MakeDecNumber(
                    heap,
                    value_buffer,
                    2048,
                    if less_than_previous {
                        hyphen.as_const()
                    } else {
                        comma.as_const()
                    },
                    value,
                )?
            } else {
                MakeDecNumber(
                    heap,
                    value_buffer,
                    2048,
                    if index != 0 {
                        comma.as_const()
                    } else {
                        SourceConstPointer::null()
                    },
                    value,
                )?
            };
            if written > 0 {
                append(heap, string_buffer, value_buffer.as_const())?;
                next = next.wrapping_add(1);
            }
        }
    }
    *overflow |= local_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeCtString(
    heap: &mut SourceHeap,
    linear_ct: SourceMutPointer<u16>,
    length_ct: i32,
    add_delimiter: i32,
    number_hydrogens: SourceConstPointer<i8>,
    number_of_atoms: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1069 MakeCtString
    // INCHI✔️✔️: int MakeCtString( CANON_GLOBALS    *pCG,
    // INCHI✔️✔️:                   AT_NUMB          *LinearCT,
    // INCHI✔️✔️:                   int              nLenCT,
    // INCHI✔️✔️:                   int              bAddDelim,
    // INCHI✔️✔️:                   S_CHAR           *nNum_H,
    // INCHI✔️✔️:                   int              num_atoms,
    // INCHI✔️✔️:                   INCHI_IOS_STRING *strbuf,
    // INCHI✔️✔️:                   int              nCtMode,
    // INCHI✔️✔️:                   int              *bOverflow )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!nNum_H || !( nCtMode & CT_MODE_NO_ORPHANS ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return MakeCtStringOld( LinearCT, nLenCT, bAddDelim, strbuf, nCtMode, bOverflow );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return MakeCtStringNew( pCG, LinearCT, nLenCT, bAddDelim, nNum_H,
    // INCHI✔️✔️:                                 num_atoms, strbuf, nCtMode, bOverflow );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: MakeCtString

    if number_hydrogens.is_null() || ct_mode & crate::source_types::CT_MODE_NO_ORPHANS as i32 == 0 {
        MakeCtStringOld(
            heap,
            linear_ct.as_const(),
            length_ct,
            add_delimiter,
            string_buffer,
            ct_mode,
            overflow,
        )
    } else {
        MakeCtStringNew(
            heap,
            linear_ct,
            length_ct,
            add_delimiter,
            number_hydrogens,
            number_of_atoms,
            string_buffer,
            ct_mode,
            overflow,
        )
    }
}

#[allow(non_snake_case)]
pub(crate) fn MakeCRVString(
    heap: &mut SourceHeap,
    original_info: SourceConstPointer<ORIG_INFO>,
    length_ct: i32,
    add_delimiter: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1308 MakeCRVString
    // INCHI✔️❌: int MakeCRVString( ORIG_INFO        *OrigInfo,
    // INCHI✔️❌:                    int              nLenCT,
    // INCHI✔️❌:                    int              bAddDelim,
    // INCHI✔️❌:                    INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                    int              nCtMode,
    // INCHI✔️❌:                    int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, k, bAbcNumbers;
    // INCHI✔️❌:     int bOvfl = *bOverflow;
    // INCHI✔️❌:     char szValue[2048] = { 0 };
    // INCHI✔️❌:     int  bNext = 0;
    // INCHI✔️❌:     bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add connection table string */
    // INCHI✔️❌:     if (!bOvfl && bAddDelim)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, ", " );
    // INCHI✔️❌:         /*if ( nLen_szLinearCT > 2 ) {
    // INCHI✔️❌:             strcpy( szLinearCT, ", " );
    // INCHI✔️❌:             nLen += 2;
    // INCHI✔️❌:         } else {
    // INCHI✔️❌:             bOvfl = 1;
    // INCHI✔️❌:         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k = 0; !bOvfl && k < nLenCT/* && nLen < nLen_szLinearCT*/; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  find the next non-empty entry */
    // INCHI✔️❌:         if (OrigInfo[k].cCharge || OrigInfo[k].cRadical || OrigInfo[k].cUnusualValence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (bAbcNumbers)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                 3 items: Ad+3d4 (canon. numb=Ad, charge=+3, doublet, valence = 4
    // INCHI✔️❌:                 2 items: Ad.d4  Ad+3.4  Ad+3d
    // INCHI✔️❌:                 1 item:  Ad+3   Ad.d     Ad4
    // INCHI✔️❌:
    // INCHI✔️❌:                 dot output before radical: no charge, radical is present
    // INCHI✔️❌:                 dot before valence:        charge is present, no radical, valence is present
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 len = MakeAbcNumber( szValue, ( int )sizeof( szValue ), NULL, k + 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* charge */
    // INCHI✔️❌:                 if (OrigInfo[k].cCharge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (OrigInfo[k].cCharge > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, "+", OrigInfo[k].cCharge );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, OrigInfo[k].cCharge );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* radical */
    // INCHI✔️❌:                 if (OrigInfo[k].cRadical)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!OrigInfo[k].cCharge)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (len >= 2047) /* djb-rwth: fixing coverity ID #499515 */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = 2047;
    // INCHI✔️❌:                             goto early_break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else if (len < 0) /* djb-rwth: fixing coverity ID #500400 */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = 0;
    // INCHI✔️❌:                             goto early_break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             szValue[len] = '.';
    // INCHI✔️❌:                             len++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (len >= 2047) /* djb-rwth: fixing coverity ID #499515 */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = 2047;
    // INCHI✔️❌:                         goto early_break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (len < 0) /* djb-rwth: fixing coverity ID #500382 */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = 0;
    // INCHI✔️❌:                         goto early_break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         switch (OrigInfo[k].cRadical)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case 1:
    // INCHI✔️❌:                                 /* djb-rwth: fixing coverity ID #499515 -- false positive, len tested for overflow */
    // INCHI✔️❌:                                 szValue[len] = 'd';
    // INCHI✔️❌:                                 len++;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case 2:
    // INCHI✔️❌:                                 szValue[len] = 't';
    // INCHI✔️❌:                                 len++;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             default:
    // INCHI✔️❌:                                 szValue[len] = 'u';
    // INCHI✔️❌:                                 len++;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* valence */
    // INCHI✔️❌:                 if (OrigInfo[k].cUnusualValence)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (OrigInfo[k].cCharge && !OrigInfo[k].cRadical)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (len >= 2047) /* djb-rwth: fixing coverity ID #499515 */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = 2047;
    // INCHI✔️❌:                             goto early_break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else if (len < 0) /* djb-rwth: fixing coverity ID #500382 */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = 0;
    // INCHI✔️❌:                             goto early_break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             szValue[len] = '.';
    // INCHI✔️❌:                             len++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, OrigInfo[k].cUnusualValence );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                     3 items: 22+3d4 (canon. numb=22, charge=+3, doublet, valence = 4
    // INCHI✔️❌:                     2 items: 22d4  22+3.4  22+3d
    // INCHI✔️❌:                     1 item:  22+3  22d     22.4
    // INCHI✔️❌:
    // INCHI✔️❌:                     dot output before valence:
    // INCHI✔️❌:                                 (a) charge,    no radical, valence
    // INCHI✔️❌:                                 (b) no charge, no radical, valence
    // INCHI✔️❌:                     that is, whenever valence is present and no radical
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 len = MakeDecNumber( szValue, ( int )sizeof( szValue ), bNext ? ITEM_DELIMETER : NULL, k + 1 );
    // INCHI✔️❌:                 /* charge */
    // INCHI✔️❌:                 if (OrigInfo[k].cCharge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (OrigInfo[k].cCharge > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, "+", OrigInfo[k].cCharge );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, OrigInfo[k].cCharge );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* radical */
    // INCHI✔️❌:                 if (OrigInfo[k].cRadical)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (len >= 2047) /* djb-rwth: fixing coverity ID #499515 */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = 2047;
    // INCHI✔️❌:                         goto early_break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (len < 0) /* djb-rwth: fixing coverity ID #500382 */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = 0;
    // INCHI✔️❌:                         goto early_break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         switch (OrigInfo[k].cRadical)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case 1:
    // INCHI✔️❌:                                 szValue[len] = 'd'; /* djb-rwth: GCC 14 false positive */
    // INCHI✔️❌:                                 len++;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case 2:
    // INCHI✔️❌:                                 szValue[len] = 't';
    // INCHI✔️❌:                                 len++;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             default:
    // INCHI✔️❌:                                 szValue[len] = 'u';
    // INCHI✔️❌:                                 len++;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* valence */
    // INCHI✔️❌:                 if (OrigInfo[k].cUnusualValence)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!OrigInfo[k].cRadical)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (len >= 2047) /* djb-rwth: fixing coverity ID #499515 */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = 2047;
    // INCHI✔️❌:                             goto early_break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else if (len < 0) /* djb-rwth: fixing coverity ID #500382 */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = 0;
    // INCHI✔️❌:                             goto early_break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             szValue[len] = '.';
    // INCHI✔️❌:                             len++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     len += MakeDecNumber(szValue + len, (int)sizeof(szValue) - len, NULL, OrigInfo[k].cUnusualValence);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: early_break:
    // INCHI✔️❌:         if (len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szValue[len] = '\0';
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:             bNext++;
    // INCHI✔️❌:             szValue[0] = '\0';
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*if ( len && nLen+len < nLen_szLinearCT ) {
    // INCHI✔️❌:             strcpy( szLinearCT+nLen, szValue );
    // INCHI✔️❌:             nLen += len;
    // INCHI✔️❌:             bNext ++;
    // INCHI✔️❌:         } else
    // INCHI✔️❌:         if ( len ) {
    // INCHI✔️❌:             bOvfl = 1;
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeCRVString

    let initial_length = string_buffer.nUsedLength;
    let local_overflow = *overflow;
    let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS as i32;
    let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
    let comma_space = heap.allocate_model_storage(vec![b',' as i8, b' ' as i8, 0])?;
    let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
    let plus = heap.allocate_model_storage(vec![b'+' as i8, 0])?;
    let string_format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;

    if local_overflow == 0 && add_delimiter != 0 {
        inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            comma_space.as_const(),
            &SourceVaList::default(),
        )?;
    }

    let mut next = 0_i32;
    let mut index = 0_i32;
    while local_overflow == 0 && index < length_ct {
        let info = heap
            .slice(original_info.offset(i64::from(index))?)?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut length = 0_i32;
        let mut early_break = false;
        if info.cCharge != 0 || info.cRadical != 0 || info.cUnusualValence != 0 {
            length = if abc_numbers != 0 {
                MakeAbcNumber(
                    heap,
                    value_buffer,
                    2048,
                    SourceConstPointer::null(),
                    index.wrapping_add(1),
                )?
            } else {
                MakeDecNumber(
                    heap,
                    value_buffer,
                    2048,
                    if next != 0 {
                        comma.as_const()
                    } else {
                        SourceConstPointer::null()
                    },
                    index.wrapping_add(1),
                )?
            };

            if info.cCharge != 0 {
                let destination = value_buffer.offset(i64::from(length))?;
                let written = MakeDecNumber(
                    heap,
                    destination,
                    2048_i32.wrapping_sub(length),
                    if info.cCharge > 0 {
                        plus.as_const()
                    } else {
                        SourceConstPointer::null()
                    },
                    i32::from(info.cCharge),
                )?;
                length = length.wrapping_add(written);
            }

            if info.cRadical != 0 {
                if abc_numbers != 0 && info.cCharge == 0 {
                    if length >= 2047 {
                        length = 2047;
                        early_break = true;
                    } else if length < 0 {
                        length = 0;
                        early_break = true;
                    } else {
                        heap.slice_mut(value_buffer)?[length as usize] = b'.' as i8;
                        length = length.wrapping_add(1);
                    }
                }
                if !early_break {
                    if length >= 2047 {
                        length = 2047;
                        early_break = true;
                    } else if length < 0 {
                        length = 0;
                        early_break = true;
                    } else {
                        heap.slice_mut(value_buffer)?[length as usize] = match info.cRadical {
                            1 => b'd' as i8,
                            2 => b't' as i8,
                            _ => b'u' as i8,
                        };
                        length = length.wrapping_add(1);
                    }
                }
            }

            if !early_break && info.cUnusualValence != 0 {
                let needs_dot = if abc_numbers != 0 {
                    info.cCharge != 0 && info.cRadical == 0
                } else {
                    info.cRadical == 0
                };
                if needs_dot {
                    if length >= 2047 {
                        length = 2047;
                        early_break = true;
                    } else if length < 0 {
                        length = 0;
                        early_break = true;
                    } else {
                        heap.slice_mut(value_buffer)?[length as usize] = b'.' as i8;
                        length = length.wrapping_add(1);
                    }
                }
                if !early_break {
                    let destination = value_buffer.offset(i64::from(length))?;
                    let written = MakeDecNumber(
                        heap,
                        destination,
                        2048_i32.wrapping_sub(length),
                        SourceConstPointer::null(),
                        i32::from(info.cUnusualValence),
                    )?;
                    length = length.wrapping_add(written);
                }
            }
        }

        if length != 0 {
            let length_index =
                usize::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            *heap
                .slice_mut(value_buffer)?
                .get_mut(length_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            inchi_strbuf_printf(
                heap,
                Some(string_buffer),
                string_format.as_const(),
                &SourceVaList {
                    arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                    ..SourceVaList::default()
                },
            )?;
            next = next.wrapping_add(1);
            heap.slice_mut(value_buffer)?[0] = 0;
        }
        index = index.wrapping_add(1);
    }
    *overflow |= local_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn bHasOrigInfo(
    heap: &SourceHeap,
    original_info: SourceConstPointer<ORIG_INFO>,
    number_of_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:369 bHasOrigInfo
    // INCHI✔️✔️: int bHasOrigInfo( ORIG_INFO *OrigInfo, int num_atoms )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, bFound = 0;
    // INCHI✔️✔️:     if (OrigInfo && num_atoms > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0; !bFound && i < num_atoms; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bFound |= ( 0 != OrigInfo[i].cCharge ) ||
    // INCHI✔️✔️:                 ( 0 != OrigInfo[i].cRadical ) ||
    // INCHI✔️✔️:                 ( 0 != OrigInfo[i].cUnusualValence );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return bFound;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: bHasOrigInfo

    let mut found = 0_i32;
    if !original_info.is_null() && number_of_atoms > 0 {
        let mut index = 0_i32;
        while found == 0 && index < number_of_atoms {
            let info = heap
                .slice(original_info.offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            found |=
                i32::from(info.cCharge != 0 || info.cRadical != 0 || info.cUnusualValence != 0);
            index = index.wrapping_add(1);
        }
    }
    Ok(found)
}

#[allow(non_snake_case)]
pub(crate) fn EqlOrigInfo(
    heap: &SourceHeap,
    aux1: Option<&crate::source_types::INChI_Aux>,
    aux2: Option<&crate::source_types::INChI_Aux>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:387 EqlOrigInfo
    // INCHI✔️✔️: int EqlOrigInfo( INChI_Aux *a1, INChI_Aux *a2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int ret = a1 && a2 && a1->nNumberOfAtoms == a2->nNumberOfAtoms &&
    // INCHI✔️✔️:         bHasOrigInfo( a1->OrigInfo, a1->nNumberOfAtoms ) && a2->OrigInfo &&
    // INCHI✔️✔️:         !memcmp( a1->OrigInfo, a2->OrigInfo, a1->nNumberOfAtoms * sizeof( a1->OrigInfo[0] ) );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: EqlOrigInfo

    let (Some(aux1), Some(aux2)) = (aux1, aux2) else {
        return Ok(0);
    };
    if aux1.nNumberOfAtoms != aux2.nNumberOfAtoms
        || bHasOrigInfo(heap, aux1.OrigInfo.as_const(), aux1.nNumberOfAtoms)? == 0
        || aux2.OrigInfo.is_null()
    {
        return Ok(0);
    }
    let length =
        usize::try_from(aux1.nNumberOfAtoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let first = heap
        .slice(aux1.OrigInfo.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = heap
        .slice(aux2.OrigInfo.as_const())?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(first == second))
}

#[allow(non_snake_case)]
pub(crate) fn bHasEquString(
    heap: &SourceHeap,
    linear_ct: SourceConstPointer<u16>,
    linear_ct_length: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:398 bHasEquString
    // INCHI✔️✔️: int bHasEquString( AT_NUMB *LinearCT, int nLenCT )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /*  produce output string; */
    // INCHI✔️✔️:     int i, k;
    // INCHI✔️✔️:     if (!LinearCT)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (k = 0; k < nLenCT; k++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  find the first equivalence number */
    // INCHI✔️✔️:         if (k != (int) LinearCT[k] - 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         for (i = k; i < nLenCT; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (k != (int) LinearCT[i] - 1)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (k < i)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: bHasEquString

    if linear_ct.is_null() {
        return Ok(0);
    }
    let mut group = 0_i32;
    while group < linear_ct_length {
        let representative = i32::from(
            *heap
                .slice(linear_ct.offset(i64::from(group))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if group == representative.wrapping_sub(1) {
            let mut index = group;
            while index < linear_ct_length {
                let member_representative = i32::from(
                    *heap
                        .slice(linear_ct.offset(i64::from(index))?)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                if group == member_representative.wrapping_sub(1) && group < index {
                    return Ok(1);
                }
                index = index.wrapping_add(1);
            }
        }
        group = group.wrapping_add(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn MakeMult(
    heap: &mut SourceHeap,
    multiplier: i32,
    trailing_delimiter: SourceConstPointer<i8>,
    buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:431 MakeMult
    // INCHI✔️❌: int MakeMult( int mult,
    // INCHI✔️❌:               const char       *szTailingDelim,
    // INCHI✔️❌:               INCHI_IOS_STRING *buf,
    // INCHI✔️❌:               int              nCtMode,
    // INCHI✔️❌:               int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     int  len = 0, len_delim, n;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (mult == 1 || *bOverflow)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len += MakeAbcNumber( szValue, ( int )sizeof( szValue ), NULL, mult );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len += MakeDecNumber( szValue, ( int )sizeof( szValue ), NULL, mult );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     len_delim = (int) strlen( szTailingDelim );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (len + len_delim < ( int )sizeof( szValue ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         strcpy(szValue + len, szTailingDelim);
    // INCHI✔️❌:         n = inchi_strbuf_printf( buf, "%s", szValue );
    // INCHI✔️❌:         if (-1 == n) *bOverflow |= 1;
    // INCHI✔️❌:         return n;
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         len += len_delim;
    // INCHI✔️❌:         if ( len < nLen_szLinearCT )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             strcpy( szLinearCT, szValue );
    // INCHI✔️❌:             return len;
    // INCHI✔️❌:         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *bOverflow |= 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeMult

    if multiplier == 1 || *overflow != 0 {
        return Ok(0);
    }
    let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
    let length = if ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0 {
        MakeAbcNumber(
            heap,
            value_buffer,
            2048,
            SourceConstPointer::null(),
            multiplier,
        )?
    } else {
        MakeDecNumber(
            heap,
            value_buffer,
            2048,
            SourceConstPointer::null(),
            multiplier,
        )?
    };
    let delimiter = heap.slice(trailing_delimiter)?;
    let delimiter_length = delimiter
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let delimiter = delimiter[..delimiter_length].to_vec();
    if length.wrapping_add(
        i32::try_from(delimiter_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
    ) < 2048
    {
        let length = usize::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let output = heap.slice_mut(value_buffer)?;
        output[length..length + delimiter.len()].copy_from_slice(&delimiter);
        output[length + delimiter.len()] = 0;
        let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
        return match inchi_strbuf_printf(
            heap,
            Some(buffer),
            format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                ..SourceVaList::default()
            },
        ) {
            Ok(value) => {
                if value == -1 {
                    *overflow |= 1;
                }
                Ok(value)
            }
            Err(SourceHeapError::AllocationFailed) => {
                *overflow |= 1;
                Ok(-1)
            }
            Err(error) => Err(error),
        };
    }
    *overflow |= 1;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn MakeDelim(
    heap: &mut SourceHeap,
    trailing_delimiter: SourceConstPointer<i8>,
    buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:476 MakeDelim
    // INCHI✔️❌: int MakeDelim( const char       *szTailingDelim,
    // INCHI✔️❌:                INCHI_IOS_STRING *buf,
    // INCHI✔️❌:                int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int n;
    // INCHI✔️❌:     if (!szTailingDelim || !*szTailingDelim || *bOverflow)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     n = inchi_strbuf_printf( buf, szTailingDelim );
    // INCHI✔️❌:     if (-1 == n)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *bOverflow |= 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return n;
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     len_delim = (int) strlen(szTailingDelim);
    // INCHI✔️❌:     if ( len_delim < nLen_szLinearCT ) {
    // INCHI✔️❌:     strcpy( szLinearCT, szTailingDelim );
    // INCHI✔️❌:     return len_delim;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= 1;
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌:     */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeDelim

    if trailing_delimiter.is_null() || *overflow != 0 {
        return Ok(0);
    }
    let delimiter = heap.slice(trailing_delimiter)?;
    if *delimiter
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        == 0
    {
        return Ok(0);
    }
    match inchi_strbuf_printf(
        heap,
        Some(buffer),
        trailing_delimiter,
        &SourceVaList::default(),
    ) {
        Ok(value) => {
            if value == -1 {
                *overflow |= 1;
            }
            Ok(value)
        }
        Err(SourceHeapError::AllocationFailed) => {
            *overflow |= 1;
            Ok(-1)
        }
        Err(error) => Err(error),
    }
}

#[allow(non_snake_case)]
pub(crate) fn MakeEqStr(
    heap: &mut SourceHeap,
    trailing_delimiter: SourceConstPointer<i8>,
    multiplier: i32,
    buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:506 MakeEqStr
    // INCHI✔️❌: int MakeEqStr( const char       *szTailingDelim,
    // INCHI✔️❌:                int              mult,
    // INCHI✔️❌:                INCHI_IOS_STRING *buf,
    // INCHI✔️❌:                int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int  n = 0, n0;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!szTailingDelim || !*szTailingDelim || *bOverflow)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     n0 = buf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (mult != 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         n = MakeDecNumber( szValue, ( int )sizeof( szValue ), NULL, mult );
    // INCHI✔️❌:         if (-1 == n)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *bOverflow |= 1;
    // INCHI✔️❌:             return -1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (n > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         n = inchi_strbuf_printf( buf, "%-s", szValue );
    // INCHI✔️❌:         if (-1 == n) *bOverflow |= 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     n = inchi_strbuf_printf( buf, "%-s", szTailingDelim );
    // INCHI✔️❌:     if (-1 == n)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *bOverflow |= 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( buf->nUsedLength - n0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeEqStr

    if trailing_delimiter.is_null() || *overflow != 0 {
        return Ok(0);
    }
    if *heap
        .slice(trailing_delimiter)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        == 0
    {
        return Ok(0);
    }

    let initial_length = buffer.nUsedLength;
    let format = heap.allocate_model_storage(vec![b'%' as i8, b'-' as i8, b's' as i8, 0])?;
    if multiplier != 1 {
        let value = heap.allocate_model_storage(vec![0_i8; 2048])?;
        let length = MakeDecNumber(heap, value, 2048, SourceConstPointer::null(), multiplier)?;
        if length == -1 {
            *overflow |= 1;
            return Ok(-1);
        }
        if length > 0 {
            match inchi_strbuf_printf(
                heap,
                Some(buffer),
                format.as_const(),
                &SourceVaList {
                    arguments: vec![SourceFormatArgument::Bytes(value.as_const())],
                    ..SourceVaList::default()
                },
            ) {
                Ok(-1) | Err(SourceHeapError::AllocationFailed) => *overflow |= 1,
                Ok(_) => {}
                Err(error) => return Err(error),
            }
        }
    }
    match inchi_strbuf_printf(
        heap,
        Some(buffer),
        format.as_const(),
        &SourceVaList {
            arguments: vec![SourceFormatArgument::Bytes(trailing_delimiter)],
            ..SourceVaList::default()
        },
    ) {
        Ok(-1) | Err(SourceHeapError::AllocationFailed) => *overflow |= 1,
        Ok(_) => {}
        Err(error) => return Err(error),
    }
    Ok(buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn MakeTautString(
    heap: &mut SourceHeap,
    mut linear_ct: SourceMutPointer<u16>,
    mut length_ct: i32,
    add_delimiter: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1111 MakeTautString
    // INCHI✔️❌: int MakeTautString( AT_NUMB          *LinearCT,
    // INCHI✔️❌:                     int              nLenCT,
    // INCHI✔️❌:                     int              bAddDelim,
    // INCHI✔️❌:                     INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                     int              nCtMode,
    // INCHI✔️❌:                     int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, i, bOvfl = *bOverflow;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     const char *p;
    // INCHI✔️❌:     int   nValue, nGroupLen, iGroupOutputCount, bCompressed;
    // INCHI✔️❌:     /*  make tautomer string */
    // INCHI✔️❌:     if (!nLenCT || !LinearCT || !*LinearCT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return nLen;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     bCompressed = ( nCtMode & CT_MODE_ABC_NUMBERS );
    // INCHI✔️❌:     if (!bCompressed && !bOvfl && bAddDelim)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "%s", COMMA_EXTRA_SPACE );
    // INCHI✔️❌:         /*if ( nLen_szLinearCT > 1+LEN_EXTRA_SPACE ) {
    // INCHI✔️❌:             strcpy( szLinearCT, COMMA_EXTRA_SPACE);
    // INCHI✔️❌:             nLen += 1+LEN_EXTRA_SPACE;
    // INCHI✔️❌:         } else {
    // INCHI✔️❌:             bOvfl = 1;
    // INCHI✔️❌:         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     LinearCT++; /*  bypass number of tautomeric groups */
    // INCHI✔️❌:     nLenCT--;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!bOvfl)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = nGroupLen = iGroupOutputCount = 0; i < nLenCT/* && nLen < nLen_szLinearCT*/; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nValue = (int) LinearCT[i];
    // INCHI✔️❌:             if (nGroupLen == iGroupOutputCount)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nGroupLen = nValue;
    // INCHI✔️❌:                 iGroupOutputCount = 0;
    // INCHI✔️❌:                 /* group delimiter (uncompressed) */
    // INCHI✔️❌:                 if (!bCompressed)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!i)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         strcpy( szValue, "(" );
    // INCHI✔️❌:                         len = 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         strcpy( szValue, ")(" );
    // INCHI✔️❌:                         len = 2;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     len = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bCompressed && iGroupOutputCount >= INCHI_T_NUM_MOVABLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  compressed canon number in Abc */
    // INCHI✔️❌:                     len = MakeAbcNumber( szValue, ( int )sizeof( szValue ), NULL, nValue );
    // INCHI✔️❌:                     iGroupOutputCount++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  always output number of hydrogen atoms as a decimal */
    // INCHI✔️❌:                     /*  output leading space if: */
    // INCHI✔️❌:                     /*  (a) this is the first output value in compressed mode (i==1 && bCompressed) */
    // INCHI✔️❌:                     /*  (b) this is not the first output value in non-compressed mode ( iGroupOutputCount && !bCompressed) */
    // INCHI✔️❌:                     if (bCompressed)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         p = NULL;
    // INCHI✔️❌:                         len = 0;
    // INCHI✔️❌:                         switch (iGroupOutputCount)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case 0:
    // INCHI✔️❌:                                 len = MakeDecNumber( szValue, ( int )sizeof( szValue ), ( i == 1 ) ? ITEM_DELIMETER : NULL, nValue );
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case 1:
    // INCHI✔️❌:                                 p = "-";
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case 2:
    // INCHI✔️❌:                                 p = "+";
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (p)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             switch (nValue)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 case 0:
    // INCHI✔️❌:                                     len = 0;
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 case 1:
    // INCHI✔️❌:                                     strcpy(szValue, p);
    // INCHI✔️❌:                                     len = (int) strlen( szValue );
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 default:
    // INCHI✔️❌:                                     len = MakeDecNumber( szValue, ( int )sizeof( szValue ), p, nValue );
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (iGroupOutputCount >= INCHI_T_NUM_MOVABLE)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  canonical number of the atom in the tautomeric group */
    // INCHI✔️❌:                             len = MakeDecNumber( szValue, ( int )sizeof( szValue ), ITEM_DELIMETER, nValue );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             p = NULL;
    // INCHI✔️❌:                             len = 0;
    // INCHI✔️❌:                             if (nValue)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 switch (iGroupOutputCount)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     case 0:
    // INCHI✔️❌:                                         p = "H";
    // INCHI✔️❌:                                         break;
    // INCHI✔️❌:                                     case 1:
    // INCHI✔️❌:                                         p = "-";
    // INCHI✔️❌:                                         break;
    // INCHI✔️❌:                                     case 2:
    // INCHI✔️❌:                                         p = "+";
    // INCHI✔️❌:                                         break;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (p)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     /*  number of hydrogens */
    // INCHI✔️❌:                                     if (nValue == 1)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         strcpy(szValue, p);
    // INCHI✔️❌:                                         len = (int) strlen( szValue );
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         len = MakeDecNumber( szValue, ( int )sizeof( szValue ), p, nValue );
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     iGroupOutputCount++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (len > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /*if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
    // INCHI✔️❌:                 if ( len ) {
    // INCHI✔️❌:                     strcpy( szLinearCT+nLen, szValue );
    // INCHI✔️❌:                     nLen += len;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             } else {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!bOvfl && !bCompressed && i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, ")" );
    // INCHI✔️❌:             /*if ( nLen + 1 < nLen_szLinearCT ) {
    // INCHI✔️❌:                 strcpy( szLinearCT+nLen, ")" );
    // INCHI✔️❌:                 nLen ++;
    // INCHI✔️❌:             } else {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:             }*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeTautString

    if length_ct == 0 || linear_ct.is_null() {
        return Ok(0);
    }
    let first = *heap
        .slice(linear_ct.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if first == 0 {
        return Ok(0);
    }
    let initial_used_length = string_buffer.nUsedLength;
    let initial_overflow = *overflow;
    let compressed = ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0;
    let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
    let string_format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
    let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
    let minus = heap.allocate_model_storage(vec![b'-' as i8, 0])?;
    let plus = heap.allocate_model_storage(vec![b'+' as i8, 0])?;
    let hydrogen = heap.allocate_model_storage(vec![b'H' as i8, 0])?;
    let close = heap.allocate_model_storage(vec![b')' as i8, 0])?;

    let append = |heap: &mut SourceHeap,
                  buffer: &mut INCHI_IOS_STRING,
                  value: SourceConstPointer<i8>|
     -> Result<(), SourceHeapError> {
        match inchi_strbuf_printf(
            heap,
            Some(buffer),
            string_format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(value)],
                ..SourceVaList::default()
            },
        ) {
            Ok(_) | Err(SourceHeapError::AllocationFailed) => Ok(()),
            Err(error) => Err(error),
        }
    };

    if !compressed && initial_overflow == 0 && add_delimiter != 0 {
        append(heap, string_buffer, comma.as_const())?;
    }
    linear_ct = linear_ct.offset(1)?;
    length_ct = length_ct.wrapping_sub(1);
    let mut index = 0_i32;
    if initial_overflow == 0 {
        let mut group_length = 0_i32;
        let mut group_output_count = 0_i32;
        while index < length_ct {
            let value = i32::from(
                *heap
                    .slice(linear_ct.as_const())?
                    .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let length;
            if group_length == group_output_count {
                group_length = value;
                group_output_count = 0;
                if !compressed {
                    let bytes: &[i8] = if index == 0 {
                        &[b'(' as i8, 0]
                    } else {
                        &[b')' as i8, b'(' as i8, 0]
                    };
                    heap.slice_mut(value_buffer)?[..bytes.len()].copy_from_slice(bytes);
                    length = i32::try_from(bytes.len() - 1)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                } else {
                    length = 0;
                }
            } else if compressed
                && group_output_count >= crate::source_types::INCHI_T_NUM_MOVABLE as i32
            {
                length =
                    MakeAbcNumber(heap, value_buffer, 2048, SourceConstPointer::null(), value)?;
                group_output_count = group_output_count.wrapping_add(1);
            } else {
                if compressed {
                    length = match group_output_count {
                        0 => MakeDecNumber(
                            heap,
                            value_buffer,
                            2048,
                            if index == 1 {
                                comma.as_const()
                            } else {
                                SourceConstPointer::null()
                            },
                            value,
                        )?,
                        1 | 2 => {
                            let prefix = if group_output_count == 1 { minus } else { plus };
                            match value {
                                0 => 0,
                                1 => {
                                    let prefix_bytes = heap.slice(prefix.as_const())?.to_vec();
                                    heap.slice_mut(value_buffer)?[..prefix_bytes.len()]
                                        .copy_from_slice(&prefix_bytes);
                                    1
                                }
                                _ => MakeDecNumber(
                                    heap,
                                    value_buffer,
                                    2048,
                                    prefix.as_const(),
                                    value,
                                )?,
                            }
                        }
                        _ => 0,
                    };
                } else if group_output_count >= crate::source_types::INCHI_T_NUM_MOVABLE as i32 {
                    length = MakeDecNumber(heap, value_buffer, 2048, comma.as_const(), value)?;
                } else if value != 0 {
                    let prefix = match group_output_count {
                        0 => hydrogen,
                        1 => minus,
                        2 => plus,
                        _ => SourceMutPointer::null(),
                    };
                    if prefix.is_null() {
                        length = 0;
                    } else if value == 1 {
                        let prefix_bytes = heap.slice(prefix.as_const())?.to_vec();
                        heap.slice_mut(value_buffer)?[..prefix_bytes.len()]
                            .copy_from_slice(&prefix_bytes);
                        length = i32::try_from(prefix_bytes.len() - 1)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    } else {
                        length = MakeDecNumber(heap, value_buffer, 2048, prefix.as_const(), value)?;
                    }
                } else {
                    length = 0;
                }
                group_output_count = group_output_count.wrapping_add(1);
            }
            if length > 0 {
                append(heap, string_buffer, value_buffer.as_const())?;
            }
            index = index.wrapping_add(1);
        }
        if !compressed && index != 0 {
            append(heap, string_buffer, close.as_const())?;
        }
    }
    *overflow |= initial_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_used_length))
}

#[allow(non_snake_case)]
pub(crate) fn MakeHString(
    heap: &mut SourceHeap,
    add_delimiter: i32,
    linear_ct: SourceConstPointer<i8>,
    length_ct: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:789 MakeHString
    // INCHI✔️❌: int MakeHString( int              bAddDelim,
    // INCHI✔️❌:                  S_CHAR           *LinearCT,
    // INCHI✔️❌:                  int              nLenCT,
    // INCHI✔️❌:                  INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                  int              nCtMode,
    // INCHI✔️❌:                  int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌: #define INIT_MIN_NUM_H (-4)
    // INCHI✔️❌: #define INIT_MAX_NUM_H 16
    // INCHI✔️❌: #define INIT_LEN_NUM_H (INIT_MAX_NUM_H - INIT_MIN_NUM_H + 1)
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, i, iFirst, nVal, bOvfl = *bOverflow;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     const char *pH;
    // INCHI✔️❌:     int  bNext = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add connection table string */
    // INCHI✔️❌:     if (!( nCtMode & CT_MODE_ABC_NUMBERS ) && !bOvfl && bAddDelim)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "," );
    // INCHI✔️❌:                             /*if ( nLen_szLinearCT > 1 ) {
    // INCHI✔️❌:                                 strcpy( szLinearCT, "," );
    // INCHI✔️❌:                                 nLen ++;
    // INCHI✔️❌:                                 } else {
    // INCHI✔️❌:                                     bOvfl = 1;
    // INCHI✔️❌:                             }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!bOvfl && 0 < nLenCT && LinearCT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (nCtMode & CT_MODE_EQL_H_TOGETHER)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int  curMinH = INIT_MIN_NUM_H;
    // INCHI✔️❌:             int  curMaxH = INIT_MAX_NUM_H;
    // INCHI✔️❌:             int  curLenH = INIT_LEN_NUM_H;
    // INCHI✔️❌:             int  nInitNumH[INIT_LEN_NUM_H];
    // INCHI✔️❌:             int *nNumH = nInitNumH;
    // INCHI✔️❌:             int  numAt, curNumH;
    // INCHI✔️❌:             int      j, bOutOfRange, tot_num_no_H;
    // INCHI✔️❌:             /* count atoms H */
    // INCHI✔️❌:             do
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bOutOfRange = 0;
    // INCHI✔️❌:                 tot_num_no_H = 0; /* number of atoms that have no H */
    // INCHI✔️❌:                 memset( nNumH, 0, curLenH * sizeof( nNumH[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                 for (i = 0; i < nLenCT; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     curNumH = LinearCT[i];
    // INCHI✔️❌:                     if (curNumH < curMinH)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         curMinH = curNumH;
    // INCHI✔️❌:                         bOutOfRange++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (curNumH > curMaxH)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             curMaxH = curNumH;
    // INCHI✔️❌:                             bOutOfRange++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                             if (!bOutOfRange)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 nNumH[curNumH - curMinH] ++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     tot_num_no_H += !curNumH;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (tot_num_no_H == nLenCT)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return nLen; /* empty string */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (bOutOfRange)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* for debug only */
    // INCHI✔️❌:                     if (nNumH != nInitNumH)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *bOverflow |= 1;
    // INCHI✔️❌:                         inchi_free( nNumH );
    // INCHI✔️❌:                         return nLen;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* end debug */
    // INCHI✔️❌:                     curLenH = curMaxH - curMinH + 1;
    // INCHI✔️❌:                     nNumH = (int*) inchi_malloc( curLenH * sizeof( nNumH[0] ) );
    // INCHI✔️❌:                     if (!nNumH)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *bOverflow |= 1;
    // INCHI✔️❌:                         return nLen;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             while (bOutOfRange); /* the loop may be executed 1 or 2 times only */
    // INCHI✔️❌:
    // INCHI✔️❌:             for (curNumH = curMinH; curNumH <= curMaxH; curNumH++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 numAt = nNumH[curNumH - curMinH]; /* number of atoms that have curNumH atoms H */
    // INCHI✔️❌:                 if (!numAt || !curNumH)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue; /* no atom has this number of H or number of H = 0 */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 j = 0;
    // INCHI✔️❌:                 while (j < nLenCT && numAt)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (curNumH == LinearCT[j])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iFirst = ++j;
    // INCHI✔️❌:                         numAt--;
    // INCHI✔️❌:                         for (; j < nLenCT && curNumH == LinearCT[j] && numAt; j++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             numAt--;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = MakeAbcNumber( szValue, ( int )sizeof( szValue ), NULL, iFirst );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = MakeDecNumber( szValue, ( int )sizeof( szValue ), bNext ? ITEM_DELIMETER : NULL, iFirst );
    // INCHI✔️❌:                             bNext++; /* add a delimiter (comma) before all except the first */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (iFirst < j)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* output last canonical number */
    // INCHI✔️❌:                             if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 len += MakeAbcNumber( szValue + len, ( int )sizeof( szValue ), NULL, j );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, "-", j );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (!numAt || ( nCtMode & CT_MODE_ABC_NUMBERS ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* add number of H */
    // INCHI✔️❌:                             /* output number of H */
    // INCHI✔️❌:                             nVal = curNumH;
    // INCHI✔️❌:                             if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, nVal );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 pH = nVal > 0 ? "H" : "h";
    // INCHI✔️❌:                                 nVal = abs( nVal );
    // INCHI✔️❌:                                 if (nVal > 1)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, pH, nVal );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     strcpy(szValue + len, pH);
    // INCHI✔️❌:                                     len++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /* add to the output */
    // INCHI✔️❌:                         if (len > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             inchi_strbuf_printf( strbuf, "%-s", szValue );
    // INCHI✔️❌:                             bNext++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /*
    // INCHI✔️❌:                         if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
    // INCHI✔️❌:                             if ( len ) {
    // INCHI✔️❌:                                 strcpy( szLinearCT+nLen, szValue );
    // INCHI✔️❌:                                 nLen += len;
    // INCHI✔️❌:                                 bNext ++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         } else {
    // INCHI✔️❌:                             bOvfl = 1;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         j++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nNumH != nInitNumH)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( nNumH );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iFirst = 0;
    // INCHI✔️❌:             for (i = iFirst + 1; i <= nLenCT/* && nLen < nLen_szLinearCT*/; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i < nLenCT && LinearCT[i] == LinearCT[iFirst])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* output identical values located at i = iFirst..i-1 */
    // INCHI✔️❌:                 if (LinearCT[iFirst])
    // INCHI✔️❌:                 { /* output only non-zero values */
    // INCHI✔️❌:                   /* first canonical number */
    // INCHI✔️❌:                     nVal = LinearCT[iFirst];
    // INCHI✔️❌:                     iFirst++;
    // INCHI✔️❌:                     if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = MakeAbcNumber( szValue, ( int )sizeof( szValue ), NULL, iFirst );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = MakeDecNumber( szValue, ( int )sizeof( szValue ), bNext ? ITEM_DELIMETER : NULL, iFirst );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (iFirst < i)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* output last canonical number */
    // INCHI✔️❌:                         if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len += MakeAbcNumber( szValue + len, ( int )sizeof( szValue ), NULL, i );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, "-", i );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* output number of H */
    // INCHI✔️❌:                     if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, NULL, nVal );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         pH = nVal > 0 ? "H" : "h";
    // INCHI✔️❌:                         nVal = abs( nVal );
    // INCHI✔️❌:                         if (nVal > 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len += MakeDecNumber( szValue + len, ( int )sizeof( szValue ) - len, pH, nVal );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             strcpy(szValue + len, pH); /* djb-rwth: GCC 14 false positive */
    // INCHI✔️❌:                             len++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (0 <= len/* && nLen+len < nLen_szLinearCT*/)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (len)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             inchi_strbuf_printf( strbuf, "%-s", szValue );
    // INCHI✔️❌:                                         /*strcpy( szLinearCT+nLen, szValue );
    // INCHI✔️❌:                                         nLen += len;*/
    // INCHI✔️❌:                             bNext++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bOvfl = 1;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 iFirst = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌:
    // INCHI✔️❌: #undef INIT_MIN_NUM_H
    // INCHI✔️❌: #undef INIT_MAX_NUM_H
    // INCHI✔️❌: #undef INIT_LEN_NUM_H
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeHString

    const INITIAL_MIN_H: i32 = -4;
    const INITIAL_MAX_H: i32 = 16;
    let initial_length = string_buffer.nUsedLength;
    let initial_overflow = *overflow;
    let abc = ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0;
    let together = ct_mode & crate::source_types::CT_MODE_EQL_H_TOGETHER as i32 != 0;
    let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
    let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
    let minus = heap.allocate_model_storage(vec![b'-' as i8, 0])?;
    let upper_h = heap.allocate_model_storage(vec![b'H' as i8, 0])?;
    let lower_h = heap.allocate_model_storage(vec![b'h' as i8, 0])?;
    let string_format = heap.allocate_model_storage(vec![b'%' as i8, b'-' as i8, b's' as i8, 0])?;
    let append =
        |heap: &mut SourceHeap, buffer: &mut INCHI_IOS_STRING, value| match inchi_strbuf_printf(
            heap,
            Some(buffer),
            string_format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(value)],
                ..SourceVaList::default()
            },
        ) {
            Ok(_) | Err(SourceHeapError::AllocationFailed) => Ok(()),
            Err(error) => Err(error),
        };
    if !abc && initial_overflow == 0 && add_delimiter != 0 {
        append(heap, string_buffer, comma.as_const())?;
    }
    if initial_overflow != 0 || length_ct <= 0 || linear_ct.is_null() {
        *overflow |= initial_overflow;
        return Ok(string_buffer.nUsedLength.wrapping_sub(initial_length));
    }
    let length = usize::try_from(length_ct).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let values = heap
        .slice(linear_ct)?
        .get(..length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let mut next = 0_i32;

    let write_run = |heap: &mut SourceHeap,
                     buffer: &mut INCHI_IOS_STRING,
                     first: i32,
                     last: i32,
                     h_value: i32,
                     include_h: bool,
                     next: &mut i32|
     -> Result<(), SourceHeapError> {
        let mut written = if abc {
            MakeAbcNumber(heap, value_buffer, 2048, SourceConstPointer::null(), first)?
        } else {
            let result = MakeDecNumber(
                heap,
                value_buffer,
                2048,
                if *next != 0 {
                    comma.as_const()
                } else {
                    SourceConstPointer::null()
                },
                first,
            )?;
            if together {
                *next = next.wrapping_add(1);
            }
            result
        };
        if first < last {
            written += if abc {
                MakeAbcNumber(
                    heap,
                    value_buffer.offset(i64::from(written))?,
                    2048,
                    SourceConstPointer::null(),
                    last,
                )?
            } else {
                MakeDecNumber(
                    heap,
                    value_buffer.offset(i64::from(written))?,
                    2048 - written,
                    minus.as_const(),
                    last,
                )?
            };
        }
        if include_h {
            if abc {
                written += MakeDecNumber(
                    heap,
                    value_buffer.offset(i64::from(written))?,
                    2048 - written,
                    SourceConstPointer::null(),
                    h_value,
                )?;
            } else {
                let h = if h_value > 0 { upper_h } else { lower_h };
                let magnitude = h_value.abs();
                if magnitude > 1 {
                    written += MakeDecNumber(
                        heap,
                        value_buffer.offset(i64::from(written))?,
                        2048 - written,
                        h.as_const(),
                        magnitude,
                    )?;
                } else {
                    let h_bytes = heap.slice(h.as_const())?.to_vec();
                    heap.slice_mut(value_buffer.offset(i64::from(written))?)?[..h_bytes.len()]
                        .copy_from_slice(&h_bytes);
                    written += 1;
                }
            }
        }
        if written > 0 {
            append(heap, buffer, value_buffer.as_const())?;
            *next = next.wrapping_add(1);
        }
        Ok(())
    };

    if together {
        let zero_count = values.iter().filter(|value| **value == 0).count();
        if zero_count == length {
            return Ok(0);
        }
        let mut minimum = INITIAL_MIN_H;
        let mut maximum = INITIAL_MAX_H;
        for value in &values {
            minimum = minimum.min(i32::from(*value));
            maximum = maximum.max(i32::from(*value));
        }
        let count_length = usize::try_from(maximum - minimum + 1)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let dynamic = minimum < INITIAL_MIN_H || maximum > INITIAL_MAX_H;
        let counts_pointer = if dynamic {
            match heap.allocate(vec![0_i32; count_length]) {
                Ok(pointer) => Some(pointer),
                Err(SourceHeapError::AllocationFailed) => {
                    *overflow |= 1;
                    return Ok(0);
                }
                Err(error) => return Err(error),
            }
        } else {
            None
        };
        let mut counts = vec![0_i32; count_length];
        for value in &values {
            let index = usize::try_from(i32::from(*value) - minimum)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            counts[index] = counts[index].wrapping_add(1);
        }
        if let Some(pointer) = counts_pointer {
            heap.slice_mut(pointer)?.copy_from_slice(&counts);
        }
        for current_h in minimum..=maximum {
            let mut remaining = counts[usize::try_from(current_h - minimum)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
            if remaining == 0 || current_h == 0 {
                continue;
            }
            let mut position = 0_usize;
            while position < length && remaining != 0 {
                if i32::from(values[position]) == current_h {
                    let first = i32::try_from(position + 1)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    position += 1;
                    remaining -= 1;
                    while position < length
                        && i32::from(values[position]) == current_h
                        && remaining != 0
                    {
                        position += 1;
                        remaining -= 1;
                    }
                    write_run(
                        heap,
                        string_buffer,
                        first,
                        i32::try_from(position)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                        current_h,
                        remaining == 0 || abc,
                        &mut next,
                    )?;
                } else {
                    position += 1;
                }
            }
        }
        if let Some(pointer) = counts_pointer {
            inchi_free(heap, pointer)?;
        }
    } else {
        let mut first = 0_usize;
        for index in 1..=length {
            if index < length && values[index] == values[first] {
                continue;
            }
            let value = i32::from(values[first]);
            if value != 0 {
                write_run(
                    heap,
                    string_buffer,
                    i32::try_from(first + 1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                    i32::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                    value,
                    true,
                    &mut next,
                )?;
            }
            first = index;
        }
    }
    *overflow |= initial_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn MakeAbcNumber(
    heap: &mut SourceHeap,
    string: SourceMutPointer<i8>,
    mut string_length: i32,
    leading_delimiter: SourceConstPointer<i8>,
    mut value: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2143 MakeAbcNumber
    // INCHI✔️✔️: int MakeAbcNumber( char      *szString,
    // INCHI✔️✔️:                   int        nStringLen,
    // INCHI✔️✔️:                   const char *szLeadingDelim,
    // INCHI✔️✔️:                   int        nValue )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     char *p = szString;
    // INCHI✔️✔️:     char *q;
    // INCHI✔️✔️:     int  nChar;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (nStringLen < 2)
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     while (szLeadingDelim && *szLeadingDelim && --nStringLen)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *p++ = *szLeadingDelim++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nStringLen < 2)
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     if (!nValue)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *p++ = ALPHA_ZERO_VAL;  /*  zero value (cannot use 0) */
    // INCHI✔️✔️:         *p = '\0';
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nValue < 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *p++ = ALPHA_MINUS;
    // INCHI✔️✔️:         nStringLen--;
    // INCHI✔️✔️:         nValue = -nValue;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (q = p; nValue && --nStringLen; nValue /= ALPHA_BASE)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if ((nChar = nValue % ALPHA_BASE)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nChar = ALPHA_ONE + nChar - 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nChar = ALPHA_ZERO;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         *q++ = nChar;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nStringLen <= 0)
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     *q = '\0';
    // INCHI✔️✔️:     mystrrev( p );
    // INCHI✔️✔️:     p[0] = toupper( p[0] );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return (int) ( q - szString );
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: MakeAbcNumber

    if string_length < 2 {
        return Ok(-1);
    }
    let mut position = 0_usize;
    if !leading_delimiter.is_null() {
        let delimiter = heap.slice(leading_delimiter)?;
        let delimiter_nul = delimiter
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let delimiter = delimiter[..delimiter_nul].to_vec();
        for byte in delimiter {
            string_length = string_length.wrapping_sub(1);
            if string_length == 0 {
                break;
            }
            *heap
                .slice_mut(string)?
                .get_mut(position)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = byte;
            position += 1;
        }
    }
    if string_length < 2 {
        return Ok(-1);
    }
    if value == 0 {
        let output = heap.slice_mut(string)?;
        *output
            .get_mut(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = b'.' as i8;
        position += 1;
        *output
            .get_mut(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        return Ok(1);
    }
    if value < 0 {
        *heap
            .slice_mut(string)?
            .get_mut(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = b'-' as i8;
        position += 1;
        string_length = string_length.wrapping_sub(1);
        value = value
            .checked_neg()
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    }
    let digit_start = position;
    while value != 0 {
        string_length = string_length.wrapping_sub(1);
        if string_length == 0 {
            break;
        }
        let remainder = value % 27;
        let character = if remainder != 0 {
            b'a' as i32 + remainder - 1
        } else {
            b'@' as i32
        };
        *heap
            .slice_mut(string)?
            .get_mut(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = character as i8;
        position += 1;
        value /= 27;
    }
    if string_length <= 0 {
        return Ok(-1);
    }
    *heap
        .slice_mut(string)?
        .get_mut(position)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    let digits = string
        .offset(i64::try_from(digit_start).map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;
    mystrrev(heap, digits)?;
    let first = heap
        .slice_mut(digits)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *first = (*first as u8).to_ascii_uppercase() as i8;
    i32::try_from(position).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

#[allow(non_snake_case)]
pub(crate) fn MakeDecNumber(
    heap: &mut SourceHeap,
    string: SourceMutPointer<i8>,
    mut string_length: i32,
    leading_delimiter: SourceConstPointer<i8>,
    mut value: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2302 MakeDecNumber
    // INCHI✔️✔️: int MakeDecNumber( char        *szString,
    // INCHI✔️✔️:                    int         nStringLen,
    // INCHI✔️✔️:                    const char  *szLeadingDelim,
    // INCHI✔️✔️:                    int         nValue )
    // INCHI✔️✔️: {
    // INCHI✔️✔️: #define DECIMAL_BASE     10
    // INCHI✔️✔️: #define DECIMAL_MINUS    '-'
    // INCHI✔️✔️: #define DECIMAL_ZERO_VAL '0'
    // INCHI✔️✔️: #define DECIMAL_ONE      '1'
    // INCHI✔️✔️: #define DECIMAL_ZERO     '0'
    // INCHI✔️✔️:
    // INCHI✔️✔️:     char *p = szString;
    // INCHI✔️✔️:     char *q;
    // INCHI✔️✔️:     int  nChar;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (nStringLen < 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     while (szLeadingDelim && *szLeadingDelim && --nStringLen)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *p++ = *szLeadingDelim++; /* djb-rwth: GCC 14 false positive */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nStringLen < 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (!nValue)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *p++ = DECIMAL_ZERO_VAL;  /*  zero value (cannot use 0) */
    // INCHI✔️✔️:         *p = '\0';
    // INCHI✔️✔️:         return (int) ( p - szString );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nValue < 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *p++ = DECIMAL_MINUS;
    // INCHI✔️✔️:         nStringLen--;
    // INCHI✔️✔️:         nValue = -nValue;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (q = p; nValue && --nStringLen; nValue /= DECIMAL_BASE)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if ((nChar = nValue % DECIMAL_BASE)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nChar = DECIMAL_ONE + nChar - 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nChar = DECIMAL_ZERO;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         *q++ = nChar;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nStringLen <= 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     *q = '\0';
    // INCHI✔️✔️:     mystrrev( p );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return (int) ( q - szString );
    // INCHI✔️✔️:
    // INCHI✔️✔️: #undef DECIMAL_BASE
    // INCHI✔️✔️: #undef DECIMAL_MINUS
    // INCHI✔️✔️: #undef DECIMAL_ZERO_VAL
    // INCHI✔️✔️: #undef DECIMAL_ONE
    // INCHI✔️✔️: #undef DECIMAL_ZERO
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: MakeDecNumber

    if string_length < 2 {
        return Ok(-1);
    }
    let mut position = 0_usize;
    if !leading_delimiter.is_null() {
        let delimiter = heap.slice(leading_delimiter)?;
        let delimiter_nul = delimiter
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let delimiter = delimiter[..delimiter_nul].to_vec();
        for byte in delimiter {
            string_length = string_length.wrapping_sub(1);
            if string_length == 0 {
                break;
            }
            *heap
                .slice_mut(string)?
                .get_mut(position)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = byte;
            position += 1;
        }
    }
    if string_length < 2 {
        return Ok(-1);
    }
    if value == 0 {
        let output = heap.slice_mut(string)?;
        output[position] = b'0' as i8;
        position += 1;
        output[position] = 0;
        return i32::try_from(position).map_err(|_| SourceHeapError::SourceIntegerOverflow);
    }
    if value < 0 {
        *heap
            .slice_mut(string)?
            .get_mut(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = b'-' as i8;
        position += 1;
        string_length = string_length.wrapping_sub(1);
        value = value
            .checked_neg()
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    }
    let digit_start = position;
    while value != 0 {
        string_length = string_length.wrapping_sub(1);
        if string_length == 0 {
            break;
        }
        let remainder = value % 10;
        let character = if remainder != 0 {
            b'1' as i32 + remainder - 1
        } else {
            b'0' as i32
        };
        *heap
            .slice_mut(string)?
            .get_mut(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = character as i8;
        position += 1;
        value /= 10;
    }
    if string_length <= 0 {
        return Ok(-1);
    }
    *heap
        .slice_mut(string)?
        .get_mut(position)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    mystrrev(
        heap,
        string.offset(
            i64::try_from(digit_start).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
        )?,
    )?;
    i32::try_from(position).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

#[allow(non_snake_case)]
pub(crate) fn MakeEquString(
    heap: &mut SourceHeap,
    linear_ct: SourceConstPointer<u16>,
    linear_ct_length: i32,
    add_delimiter: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1559 MakeEquString
    // INCHI✔️❌: int MakeEquString( AT_NUMB          *LinearCT,
    // INCHI✔️❌:                    int              nLenCT,
    // INCHI✔️❌:                    int              bAddDelim,
    // INCHI✔️❌:                    INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                    int              nCtMode,
    // INCHI✔️❌:                    int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, i, k, bAbcNumbers; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     int bOvfl = *bOverflow;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     int  bNext = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    // INCHI✔️❌:     /*  add connection table string */
    // INCHI✔️❌:     if (!bOvfl && bAddDelim)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, ", " );
    // INCHI✔️❌:         /*if ( nLen_szLinearCT > 2 ) {
    // INCHI✔️❌:             strcpy( szLinearCT, ", " );
    // INCHI✔️❌:             nLen += 2;
    // INCHI✔️❌:         } else {
    // INCHI✔️❌:             bOvfl = 1;
    // INCHI✔️❌:         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (k = 0; !bOvfl && k < nLenCT/* && nLen < nLen_szLinearCT*/; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  find the first equivalence number */
    // INCHI✔️❌:         if (k != (int) LinearCT[k] - 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (i = k; i < nLenCT/* && nLen < nLen_szLinearCT*/; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (k != (int) LinearCT[i] - 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /*  equivalence number: a minimal canon_number out of a group of equivalent atoms */
    // INCHI✔️❌:             /*  is at canon_number-1 position of each equivalent atom.  */
    // INCHI✔️❌:             if (bAbcNumbers)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len = MakeAbcNumber( szValue, ( int )sizeof( szValue ), ( i == k && bNext ) ? ITEM_DELIMETER : NULL, i + 1 ); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len = MakeDecNumber( szValue, ( int )sizeof( szValue ), ( i == k ) ? "(" : ITEM_DELIMETER, i + 1 ); /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:             bNext++;
    // INCHI✔️❌:             /*if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
    // INCHI✔️❌:                 strcpy( szLinearCT+nLen, szValue );
    // INCHI✔️❌:                 nLen += len;
    // INCHI✔️❌:                 bNext ++;
    // INCHI✔️❌:             } else
    // INCHI✔️❌:             if ( 0 > len ) {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, ")" );
    // INCHI✔️❌:         /*if ( !bOvfl && !bAbcNumbers ) {
    // INCHI✔️❌:             if ( nLen + 2 < nLen_szLinearCT ) {
    // INCHI✔️❌:                 strcpy( szLinearCT+nLen, ")" );
    // INCHI✔️❌:                 nLen ++;
    // INCHI✔️❌:             } else {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeEquString

    let initial_length = string_buffer.nUsedLength;
    let mut source_overflow = *overflow;
    let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
    let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
    let space_comma = heap.allocate_model_storage(vec![b',' as i8, b' ' as i8, 0])?;
    let open = heap.allocate_model_storage(vec![b'(' as i8, 0])?;
    let close = heap.allocate_model_storage(vec![b')' as i8, 0])?;
    let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
    if source_overflow == 0 && add_delimiter != 0 {
        inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(space_comma.as_const())],
                ..SourceVaList::default()
            },
        )?;
    }
    let mut next = 0_i32;
    let mut group = 0_i32;
    while source_overflow == 0 && group < linear_ct_length {
        let first = i32::from(
            *heap
                .slice(linear_ct.offset(i64::from(group))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if group != first.wrapping_sub(1) {
            group = group.wrapping_add(1);
            continue;
        }
        let mut index = group;
        while index < linear_ct_length {
            let value = i32::from(
                *heap
                    .slice(linear_ct.offset(i64::from(index))?)?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            if group == value.wrapping_sub(1) {
                let leading = if ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0 {
                    if index == group && next != 0 {
                        comma.as_const()
                    } else {
                        SourceConstPointer::null()
                    }
                } else if index == group {
                    open.as_const()
                } else {
                    comma.as_const()
                };
                if ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0 {
                    let _ =
                        MakeAbcNumber(heap, value_buffer, 2048, leading, index.wrapping_add(1))?;
                } else {
                    let _ =
                        MakeDecNumber(heap, value_buffer, 2048, leading, index.wrapping_add(1))?;
                }
                inchi_strbuf_printf(
                    heap,
                    Some(string_buffer),
                    format.as_const(),
                    &SourceVaList {
                        arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                        ..SourceVaList::default()
                    },
                )?;
                next = next.wrapping_add(1);
            }
            index = index.wrapping_add(1);
        }
        inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(close.as_const())],
                ..SourceVaList::default()
            },
        )?;
        group = group.wrapping_add(1);
    }
    *overflow |= source_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn MakeIsoAtomString(
    heap: &mut SourceHeap,
    isotopic_atoms: SourceConstPointer<crate::source_types::INChI_IsotopicAtom>,
    number_of_isotopic_atoms: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1645 MakeIsoAtomString
    // INCHI✔️❌: int MakeIsoAtomString( INChI_IsotopicAtom *IsotopicAtom,
    // INCHI✔️❌:                        int                nNumberOfIsotopicAtoms,
    // INCHI✔️❌:                        INCHI_IOS_STRING   *strbuf,
    // INCHI✔️❌:                        int                nCtMode,
    // INCHI✔️❌:                        int                *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, tot_len, ret, i, j, bOvfl = *bOverflow;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     char *p;
    // INCHI✔️❌:     int   nValue;
    // INCHI✔️❌:     int   bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    // INCHI✔️❌:     static const char letter[] = "itdh";
    // INCHI✔️❌:     static const char *h[] = { "T", "D", "H" };
    // INCHI✔️❌:     static const char *sign[] = { "-", "+" };
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!bOvfl)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < nNumberOfIsotopicAtoms/* && nLen < nLen_szLinearCT*/; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             p = szValue;
    // INCHI✔️❌:             tot_len = 0;
    // INCHI✔️❌:             for (j = 0; j < 5; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len = 0;
    // INCHI✔️❌:                 switch (j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     case 0:
    // INCHI✔️❌:                         nValue = (int) IsotopicAtom[i].nAtomNumber;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case 1:
    // INCHI✔️❌:                         nValue = (int) IsotopicAtom[i].nIsoDifference;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case 2:
    // INCHI✔️❌:                         nValue = (int) IsotopicAtom[i].nNum_T;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case 3:
    // INCHI✔️❌:                         nValue = (int) IsotopicAtom[i].nNum_D;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case 4:
    // INCHI✔️❌:                         nValue = (int) IsotopicAtom[i].nNum_H;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (!j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  atom canonical number */
    // INCHI✔️❌:                     len = ( bAbcNumbers ? MakeAbcNumber : MakeDecNumber )
    // INCHI✔️❌:                         ( p, ( int )sizeof( szValue ) - tot_len,
    // INCHI✔️❌:                           bAbcNumbers ? NULL : ( i ? ITEM_DELIMETER : EXTRA_SPACE ), nValue
    // INCHI✔️❌:                         );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bAbcNumbers)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*  Abc output */
    // INCHI✔️❌:                         switch (j)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case 1: /* nIsoDifference */
    // INCHI✔️❌:                                 len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, NULL, nValue );
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case 2: /* nNum_T */
    // INCHI✔️❌:                             case 3: /* nNum_D */
    // INCHI✔️❌:                             case 4: /* nNum_H */
    // INCHI✔️❌:                                 if (nValue)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (( int )sizeof( szValue ) - tot_len > 1)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         p[len++] = letter[j - 1];
    // INCHI✔️❌:                                         if (1 == nValue)
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             p[len] = '\0';
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         else
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             ret = MakeDecNumber( p + len, ( int )sizeof( szValue ) - tot_len - len, NULL, nValue );
    // INCHI✔️❌:                                             len = ( ret >= 0 ) ? len + ret : ret;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         len = -1; /* overflow */
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nValue)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (j == 1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /*  Decimal output */
    // INCHI✔️❌:                                 /*  signed isotopic mass difference */
    // INCHI✔️❌:                                 int subtract = ( nValue > 0 );
    // INCHI✔️❌:                                 /*  (n = mass difference) > 0 corresponds to nValue = n+1 */
    // INCHI✔️❌:                                 /*  subtract 1 from it so that mass difference for 35Cl or 12C is zero */
    // INCHI✔️❌:                                 len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, sign[nValue >= 0], abs( nValue - subtract ) );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /*  hydrogen isotope */
    // INCHI✔️❌:                                 if (nValue != 1)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, h[j - 2], nValue );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (( int )sizeof( szValue ) - tot_len > 1)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         strcpy(p, h[j - 2]);
    // INCHI✔️❌:                                         len = 1;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         len = -1; /*  overflow */
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue; /*  do not write zeroes */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (len < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bOvfl = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 tot_len += len;
    // INCHI✔️❌:                 p += len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:             /*if ( nLen+tot_len < nLen_szLinearCT )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memcpy( szLinearCT+nLen, szValue, tot_len+1 );
    // INCHI✔️❌:                 nLen += tot_len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeIsoAtomString
    let initial_length = string_buffer.nUsedLength;
    let mut local_overflow = *overflow;
    if local_overflow == 0 {
        let count = if number_of_isotopic_atoms > 0 {
            usize::try_from(number_of_isotopic_atoms)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        } else {
            0
        };
        let atoms = if count == 0 {
            Vec::new()
        } else {
            heap.slice(isotopic_atoms)?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        };
        let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
        let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
        let empty = heap.allocate_model_storage(vec![0_i8])?;
        let minus = heap.allocate_model_storage(vec![b'-' as i8, 0])?;
        let plus = heap.allocate_model_storage(vec![b'+' as i8, 0])?;
        let isotope_letters = [b't' as i8, b'd' as i8, b'h' as i8];
        let hydrogen_letters = [
            heap.allocate_model_storage(vec![b'T' as i8, 0])?,
            heap.allocate_model_storage(vec![b'D' as i8, 0])?,
            heap.allocate_model_storage(vec![b'H' as i8, 0])?,
        ];
        let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
        let abc = ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0;

        for (atom_index, atom) in atoms.iter().enumerate() {
            let values = [
                i32::from(atom.nAtomNumber),
                i32::from(atom.nIsoDifference),
                i32::from(atom.nNum_T),
                i32::from(atom.nNum_D),
                i32::from(atom.nNum_H),
            ];
            let mut total_length = 0_i32;
            for field in 0..5_usize {
                let output = value_buffer.offset(i64::from(total_length))?;
                let remaining = 2048_i32.wrapping_sub(total_length);
                let value = values[field];
                let mut length = if field == 0 {
                    let delimiter = if abc {
                        SourceConstPointer::null()
                    } else if atom_index == 0 {
                        empty.as_const()
                    } else {
                        comma.as_const()
                    };
                    if abc {
                        MakeAbcNumber(heap, output, remaining, delimiter, value)?
                    } else {
                        MakeDecNumber(heap, output, remaining, delimiter, value)?
                    }
                } else if abc {
                    if field == 1 {
                        MakeDecNumber(heap, output, remaining, SourceConstPointer::null(), value)?
                    } else if value != 0 {
                        if remaining > 1 {
                            let target = heap.slice_mut(output)?;
                            target[0] = isotope_letters[field - 2];
                            let mut written = 1_i32;
                            if value == 1 {
                                target[1] = 0;
                            } else {
                                let digits = output.offset(1)?;
                                let result = MakeDecNumber(
                                    heap,
                                    digits,
                                    remaining.wrapping_sub(1),
                                    SourceConstPointer::null(),
                                    value,
                                )?;
                                written = if result >= 0 {
                                    written.wrapping_add(result)
                                } else {
                                    result
                                };
                            }
                            written
                        } else {
                            -1
                        }
                    } else {
                        0
                    }
                } else if value == 0 {
                    0
                } else if field == 1 {
                    let subtract = i32::from(value > 0);
                    let adjusted = value.wrapping_sub(subtract);
                    let magnitude = adjusted
                        .checked_abs()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    MakeDecNumber(
                        heap,
                        output,
                        remaining,
                        if value >= 0 {
                            plus.as_const()
                        } else {
                            minus.as_const()
                        },
                        magnitude,
                    )?
                } else if value != 1 {
                    MakeDecNumber(
                        heap,
                        output,
                        remaining,
                        hydrogen_letters[field - 2].as_const(),
                        value,
                    )?
                } else if remaining > 1 {
                    let target = heap.slice_mut(output)?;
                    target[0] = [b'T' as i8, b'D' as i8, b'H' as i8][field - 2];
                    target[1] = 0;
                    1
                } else {
                    -1
                };
                if length < 0 {
                    local_overflow = 1;
                    break;
                }
                total_length = total_length.wrapping_add(length);
            }
            match inchi_strbuf_printf(
                heap,
                Some(string_buffer),
                format.as_const(),
                &SourceVaList {
                    arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                    ..SourceVaList::default()
                },
            ) {
                Ok(_) | Err(SourceHeapError::AllocationFailed) => {}
                Err(error) => return Err(error),
            }
        }
    }
    *overflow |= local_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn MakeIsoTautString(
    heap: &mut SourceHeap,
    isotopic_groups: SourceConstPointer<crate::source_types::INChI_IsotopicTGroup>,
    number_of_isotopic_groups: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1803 MakeIsoTautString
    // INCHI✔️❌: int MakeIsoTautString( INChI_IsotopicTGroup *IsotopicTGroup,
    // INCHI✔️❌:                        int                  nNumberOfIsotopicTGroups,
    // INCHI✔️❌:                        INCHI_IOS_STRING     *strbuf,
    // INCHI✔️❌:                        int                  nCtMode,
    // INCHI✔️❌:                        int                  *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, tot_len, i, j, bOvfl = *bOverflow;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     char *p;
    // INCHI✔️❌:     int   nValue;
    // INCHI✔️❌:     int   bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    // INCHI✔️❌:     static const char letter[] = "tdh";
    // INCHI✔️❌:     static const char *h[] = { "T", "D", "H" };
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add connection table string */
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     if (!bOvfl)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < nNumberOfIsotopicTGroups/* && nLen < nLen_szLinearCT*/; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             p = szValue;
    // INCHI✔️❌:             tot_len = 0;
    // INCHI✔️❌:             for (j = 0; j < 4; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 switch (j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     case 0:
    // INCHI✔️❌:                         nValue = (int) IsotopicTGroup[i].nTGroupNumber;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case 1:
    // INCHI✔️❌:                         nValue = (int) IsotopicTGroup[i].nNum_T;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case 2:
    // INCHI✔️❌:                         nValue = (int) IsotopicTGroup[i].nNum_D;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case 3:
    // INCHI✔️❌:                         nValue = (int) IsotopicTGroup[i].nNum_H;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (!j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  atom canonical number */
    // INCHI✔️❌:                     len = ( bAbcNumbers ? MakeAbcNumber : MakeDecNumber )
    // INCHI✔️❌:                         ( p, ( int )sizeof( szValue ) - tot_len,
    // INCHI✔️❌:                           bAbcNumbers ? NULL : ( i ? ITEM_DELIMETER : EXTRA_SPACE ),
    // INCHI✔️❌:                           nValue
    // INCHI✔️❌:                         );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                     if (nValue)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (bAbcNumbers)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, NULL, nValue );
    // INCHI✔️❌:                             if (len > 0)
    // INCHI✔️❌:                             { /*  make sure overflow has not happened */
    // INCHI✔️❌:                                 if (( int )sizeof( szValue ) - tot_len - len > 1)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     p[len++] = letter[j - 1];
    // INCHI✔️❌:                                     p[len] = '\0';
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     len = -1; /*  overflow */
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  hydrogen isotope */
    // INCHI✔️❌:                             if (nValue != 1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, h[j - 1], nValue );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (( int )sizeof( szValue ) - tot_len > 1)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     strcpy(p, h[j - 1]);
    // INCHI✔️❌:                                     len = 1;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     len = -1; /*  overflow */
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue; /*  do not write zeroes */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 if (len < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bOvfl = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 p += len;
    // INCHI✔️❌:                 tot_len += len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:             /*if ( nLen+tot_len < nLen_szLinearCT ) {
    // INCHI✔️❌:                 memcpy( szLinearCT+nLen, szValue, tot_len+1 );
    // INCHI✔️❌:                 nLen += tot_len;
    // INCHI✔️❌:             } else {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeIsoTautString
    let initial_length = string_buffer.nUsedLength;
    let mut local_overflow = *overflow;
    if local_overflow == 0 {
        let count = if number_of_isotopic_groups > 0 {
            usize::try_from(number_of_isotopic_groups)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        } else {
            0
        };
        let groups = if count == 0 {
            Vec::new()
        } else {
            heap.slice(isotopic_groups)?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        };
        let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
        let comma = heap.allocate_model_storage(vec![b',' as i8, 0])?;
        let empty = heap.allocate_model_storage(vec![0_i8])?;
        let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
        let isotope_letters = [b't' as i8, b'd' as i8, b'h' as i8];
        let isotope_names = [b'T' as i8, b'D' as i8, b'H' as i8];
        let abc = ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0;

        for (group_index, group) in groups.iter().enumerate() {
            let values = [
                i32::from(group.nTGroupNumber),
                i32::from(group.nNum_T),
                i32::from(group.nNum_D),
                i32::from(group.nNum_H),
            ];
            let mut total_length = 0_i32;
            for field in 0..4_usize {
                let output = value_buffer.offset(i64::from(total_length))?;
                let remaining = 2048_i32.wrapping_sub(total_length);
                let value = values[field];
                let mut length = if field == 0 {
                    let delimiter = if abc {
                        SourceConstPointer::null()
                    } else if group_index == 0 {
                        empty.as_const()
                    } else {
                        comma.as_const()
                    };
                    if abc {
                        MakeAbcNumber(heap, output, remaining, delimiter, value)?
                    } else {
                        MakeDecNumber(heap, output, remaining, delimiter, value)?
                    }
                } else if value == 0 {
                    continue;
                } else if abc {
                    let mut written =
                        MakeDecNumber(heap, output, remaining, SourceConstPointer::null(), value)?;
                    if written > 0 {
                        if remaining.wrapping_sub(written) > 1 {
                            let target = heap.slice_mut(output)?;
                            let index = usize::try_from(written)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                            target[index] = isotope_letters[field - 1];
                            target[index + 1] = 0;
                            written = written.wrapping_add(1);
                        } else {
                            written = -1;
                        }
                    }
                    written
                } else if value != 1 {
                    let delimiter =
                        heap.allocate_model_storage(vec![isotope_names[field - 1], 0])?;
                    MakeDecNumber(heap, output, remaining, delimiter.as_const(), value)?
                } else if remaining > 1 {
                    let target = heap.slice_mut(output)?;
                    target[0] = isotope_names[field - 1];
                    target[1] = 0;
                    1
                } else {
                    -1
                };
                if length < 0 {
                    local_overflow = 1;
                    break;
                }
                total_length = total_length.wrapping_add(length);
            }
            match inchi_strbuf_printf(
                heap,
                Some(string_buffer),
                format.as_const(),
                &SourceVaList {
                    arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                    ..SourceVaList::default()
                },
            ) {
                Ok(_) | Err(SourceHeapError::AllocationFailed) => {}
                Err(error) => return Err(error),
            }
        }
    }
    *overflow |= local_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn MakeIsoHString(
    heap: &mut SourceHeap,
    isotope_hydrogens: &[i32; 3],
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1924 MakeIsoHString
    // INCHI✔️❌: int MakeIsoHString( int              num_iso_H[],
    // INCHI✔️❌:                     INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                     int              nCtMode,
    // INCHI✔️❌:                     int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, tot_len, j, bOvfl = *bOverflow;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     char *p;
    // INCHI✔️❌:     int   nValue;
    // INCHI✔️❌:     int   bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    // INCHI✔️❌:     static const char letter[] = "tdh";
    // INCHI✔️❌:     static const char *h[] = { "T", "D", "H" };
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add connection table string */
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     if (!bOvfl)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         p = szValue;
    // INCHI✔️❌:         tot_len = 0;
    // INCHI✔️❌:         for (j = 1; j < 4; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nValue = num_iso_H[NUM_H_ISOTOPES - j];/* j: 1=>T, 2=>D, 3=>1H */
    // INCHI✔️❌:             if (nValue)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bAbcNumbers)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, NULL, nValue );
    // INCHI✔️❌:                     if (len > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*  make sure overflow has not happened */
    // INCHI✔️❌:                         if (( int )sizeof( szValue ) - tot_len - len > 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             p[len++] = letter[j - 1];
    // INCHI✔️❌:                             p[len] = '\0';
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = -1; /*  overflow */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  hydrogen isotope */
    // INCHI✔️❌:                     if (nValue != 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, h[j - 1], nValue );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (( int )sizeof( szValue ) - tot_len > 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             strcpy(p, h[j - 1]);
    // INCHI✔️❌:                             len = 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = -1; /*  overflow */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue; /*  do not write zeroes */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (len < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             p += len;
    // INCHI✔️❌:             tot_len += len;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:         /*if ( nLen+tot_len < nLen_szLinearCT ) {
    // INCHI✔️❌:             memcpy( szLinearCT+nLen, szValue, tot_len+1 );
    // INCHI✔️❌:             nLen += tot_len;
    // INCHI✔️❌:         } else {
    // INCHI✔️❌:             bOvfl = 1;
    // INCHI✔️❌:         }*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeIsoHString
    let initial_length = string_buffer.nUsedLength;
    let mut local_overflow = *overflow;
    if local_overflow == 0 {
        let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
        let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
        let delimiters = [
            heap.allocate_model_storage(vec![b'T' as i8, 0])?,
            heap.allocate_model_storage(vec![b'D' as i8, 0])?,
            heap.allocate_model_storage(vec![b'H' as i8, 0])?,
        ];
        let lowercase = [b't' as i8, b'd' as i8, b'h' as i8];
        let abc = ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0;
        let mut total_length = 0_i32;
        for isotope_index in 0..3_usize {
            let value = isotope_hydrogens[2 - isotope_index];
            if value == 0 {
                continue;
            }
            let output = value_buffer.offset(i64::from(total_length))?;
            let remaining = 2048_i32.wrapping_sub(total_length);
            let mut length = if abc {
                MakeDecNumber(heap, output, remaining, SourceConstPointer::null(), value)?
            } else if value != 1 {
                MakeDecNumber(
                    heap,
                    output,
                    remaining,
                    delimiters[isotope_index].as_const(),
                    value,
                )?
            } else if remaining > 1 {
                let target = heap.slice_mut(output)?;
                target[0] = [b'T' as i8, b'D' as i8, b'H' as i8][isotope_index];
                target[1] = 0;
                1
            } else {
                -1
            };
            if abc && length > 0 {
                if remaining.wrapping_sub(length) > 1 {
                    let position = usize::try_from(length)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    let target = heap.slice_mut(output)?;
                    *target
                        .get_mut(position)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = lowercase[isotope_index];
                    *target
                        .get_mut(position + 1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                    length = length.wrapping_add(1);
                } else {
                    length = -1;
                }
            }
            if length < 0 {
                local_overflow = 1;
                break;
            }
            total_length = total_length.wrapping_add(length);
        }
        match inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                ..SourceVaList::default()
            },
        ) {
            Ok(_) | Err(SourceHeapError::AllocationFailed) => {}
            Err(error) => return Err(error),
        }
    }
    *overflow |= local_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeStereoString(
    heap: &mut SourceHeap,
    atom1: SourceConstPointer<u16>,
    atom2: SourceConstPointer<u16>,
    parity: SourceConstPointer<i8>,
    mut add_delimiter: i32,
    linear_ct_length: i32,
    string_buffer: &mut INCHI_IOS_STRING,
    ct_mode: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2019 MakeStereoString
    // INCHI✔️❌: int MakeStereoString( AT_NUMB          *at1,
    // INCHI✔️❌:                       AT_NUMB          *at2,
    // INCHI✔️❌:                       S_CHAR           *parity,
    // INCHI✔️❌:                       int              bAddDelim,
    // INCHI✔️❌:                       int              nLenCT,
    // INCHI✔️❌:                       INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                       int              nCtMode,
    // INCHI✔️❌:                       int              *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  produce output string; */
    // INCHI✔️❌:     int nUsedLength0 = 0, nLen = 0, len, tot_len, i, j, bOvfl = *bOverflow;
    // INCHI✔️❌:     char szValue[2048];
    // INCHI✔️❌:     char *p;
    // INCHI✔️❌:     int   nValue;
    // INCHI✔️❌:     static const char parity_char[] = "!-+u?";
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:     bAddDelim = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!bOvfl)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < nLenCT/* && nLen < nLen_szLinearCT*/; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             p = szValue;
    // INCHI✔️❌:             tot_len = 0;
    // INCHI✔️❌:             for (j = 0; j < 3; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (j == 0 && at1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nValue = (int) at1[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (j == 1 && at2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nValue = (int) at2[i];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (j == 2 && parity)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nValue = (int) parity[i];
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (nCtMode & CT_MODE_ABC_NUMBERS)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     len = ( j == 2 ? MakeDecNumber : MakeAbcNumber )( p, ( int )sizeof( szValue ) - tot_len, NULL, nValue );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (j < 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len = MakeDecNumber( p, ( int )sizeof( szValue ) - tot_len, tot_len ? "-" : ( i || bAddDelim ) ? ITEM_DELIMETER : NULL, nValue );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (tot_len + 1 < ( int )sizeof( szValue ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             *p++ = ( 0 <= nValue && nValue <= 4 ) ? parity_char[nValue] : parity_char[0];
    // INCHI✔️❌:                             *p = '\0';
    // INCHI✔️❌:                             len = 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len = -1; /*  Overflow */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (len < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bOvfl = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 p += len;
    // INCHI✔️❌:                 tot_len += len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             if ( nLen+tot_len < nLen_szLinearCT ) {
    // INCHI✔️❌:                 memcpy(  szLinearCT+nLen, szValue, tot_len+1 );
    // INCHI✔️❌:                 nLen += tot_len;
    // INCHI✔️❌:             } else {
    // INCHI✔️❌:                 bOvfl = 1;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= bOvfl;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = strbuf->nUsedLength - nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeStereoString

    let initial_length = string_buffer.nUsedLength;
    let mut local_overflow = *overflow;
    add_delimiter = 0;
    let value_buffer = heap.allocate_model_storage(vec![0_i8; 2048])?;
    let item_delimiter = heap.allocate_model_storage(vec![b',' as i8, 0])?;
    let hyphen = heap.allocate_model_storage(vec![b'-' as i8, 0])?;
    let string_format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
    let parity_char = [b'!' as i8, b'-' as i8, b'+' as i8, b'u' as i8, b'?' as i8];

    if local_overflow == 0 {
        let mut index = 0_i32;
        while index < linear_ct_length {
            let mut total_length = 0_i32;
            let mut field = 0_i32;
            while field < 3 {
                let source_value = if field == 0 && !atom1.is_null() {
                    Some(i32::from(
                        *heap
                            .slice(atom1.offset(i64::from(index))?)?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ))
                } else if field == 1 && !atom2.is_null() {
                    Some(i32::from(
                        *heap
                            .slice(atom2.offset(i64::from(index))?)?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ))
                } else if field == 2 && !parity.is_null() {
                    Some(i32::from(
                        *heap
                            .slice(parity.offset(i64::from(index))?)?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ))
                } else {
                    None
                };
                let Some(source_value) = source_value else {
                    field = field.wrapping_add(1);
                    continue;
                };
                let position = value_buffer.offset(i64::from(total_length))?;
                let remaining = 2048_i32.wrapping_sub(total_length);
                let length = if ct_mode & CT_MODE_ABC_NUMBERS as i32 != 0 {
                    if field == 2 {
                        MakeDecNumber(
                            heap,
                            position,
                            remaining,
                            SourceConstPointer::null(),
                            source_value,
                        )?
                    } else {
                        MakeAbcNumber(
                            heap,
                            position,
                            remaining,
                            SourceConstPointer::null(),
                            source_value,
                        )?
                    }
                } else if field < 2 {
                    let delimiter = if total_length != 0 {
                        hyphen.as_const()
                    } else if index != 0 || add_delimiter != 0 {
                        item_delimiter.as_const()
                    } else {
                        SourceConstPointer::null()
                    };
                    MakeDecNumber(heap, position, remaining, delimiter, source_value)?
                } else if total_length.wrapping_add(1) < 2048 {
                    let parity_index = usize::try_from(source_value).ok();
                    let character = parity_index
                        .and_then(|parity_index| parity_char.get(parity_index))
                        .copied()
                        .unwrap_or(parity_char[0]);
                    let destination = heap.slice_mut(position)?;
                    destination[0] = character;
                    destination[1] = 0;
                    1
                } else {
                    -1
                };
                if length < 0 {
                    local_overflow = 1;
                    break;
                }
                total_length = total_length.wrapping_add(length);
                field = field.wrapping_add(1);
            }
            match inchi_strbuf_printf(
                heap,
                Some(string_buffer),
                string_format.as_const(),
                &SourceVaList {
                    arguments: vec![SourceFormatArgument::Bytes(value_buffer.as_const())],
                    ..SourceVaList::default()
                },
            ) {
                Ok(_) | Err(SourceHeapError::AllocationFailed) => {}
                Err(error) => return Err(error),
            }
            index = index.wrapping_add(1);
        }
    }
    *overflow |= local_overflow;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

pub(crate) fn abctol(
    heap: &SourceHeap,
    sz_string: SourceConstPointer<i8>,
    q: Option<&mut SourceConstPointer<i8>>,
) -> Result<i64, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2203 abctol
    // INCHI✔❌: long abctol( const char *szString, char **q )
    // INCHI✔❌: {
    // INCHI✔❌: #define __MYTOLOWER(c) ( ((c) >= 'A') && ((c) <= 'Z') ? ((c) - 'A' + 'a') : (c) )
    // INCHI✔❌:
    // INCHI✔❌:     long        val = 0;
    // INCHI✔❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔❌:     const char *p = szString, *pp = szString;
    // INCHI✔❌:     long treshold = ( LONG_MAX - 1 ) / ALPHA_BASE;
    // INCHI✔❌:
    // INCHI✔❌:     if (*p == ALPHA_MINUS)
    // INCHI✔❌:     {
    // INCHI✔❌:         p++;
    // INCHI✔❌:         /* djb-rwth: removing redundant code */
    // INCHI✔❌:     }
    // INCHI✔❌:     if (*p == ALPHA_ZERO)
    // INCHI✔❌:     {
    // INCHI✔❌:         p++;
    // INCHI✔❌:         goto exit_function;
    // INCHI✔❌:     }
    // INCHI✔❌:     if (!isupper( UCINT *p ))
    // INCHI✔❌:     {
    // INCHI✔❌:         p = szString;
    // INCHI✔❌:         goto exit_function; /* not an abc-number */
    // INCHI✔❌:     }
    // INCHI✔❌:     val = __MYTOLOWER( *p ) - ALPHA_ONE + 1;
    // INCHI✔❌:     p++;
    // INCHI✔❌:     while (*p)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (islower( UCINT *p ))
    // INCHI✔❌:         {
    // INCHI✔❌:             val *= ALPHA_BASE;
    // INCHI✔❌:             val += *p - ALPHA_ONE + 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:             if (*p == ALPHA_ZERO)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /*^^^ 1.06 */
    // INCHI✔❌:                 if (val > treshold)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     p = pp;
    // INCHI✔❌:                     val = 0;
    // INCHI✔❌:                     goto exit_function;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 /* Software version 1.07 ^^^*/
    // INCHI✔❌:                     val *= ALPHA_BASE;  /* @1.06 fuzz testing caused overflow here! */
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 break;
    // INCHI✔❌:             }
    // INCHI✔❌:         p++;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌: exit_function:
    // INCHI✔❌:     if (q)
    // INCHI✔❌:     {
    // INCHI✔❌:         *q = (char *) p;  /* cast deliberately discards const qualifier */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return val;
    // INCHI✔❌: #undef __MYTOLOWER
    // INCHI✔❌: }
    // END INCHI C FUNCTION: abctol

    const ALPHA_BASE: i64 = 27;
    const THRESHOLD: i64 = (i64::MAX - 1) / ALPHA_BASE;

    let bytes = heap.slice(sz_string)?;
    let nul = bytes
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut val = 0_i64;
    let mut position = 0_usize;
    if position < nul && bytes[position] as u8 == b'-' {
        position += 1;
    }
    if position < nul && bytes[position] as u8 == b'@' {
        position += 1;
    } else if position >= nul || !(bytes[position] as u8).is_ascii_uppercase() {
        position = 0;
    } else {
        val = i64::from(bytes[position] as u8 - b'A' + 1);
        position += 1;
        while position < nul {
            let byte = bytes[position] as u8;
            if byte.is_ascii_lowercase() {
                val = val.wrapping_mul(ALPHA_BASE);
                val = val.wrapping_add(i64::from(byte - b'a' + 1));
            } else if byte == b'@' {
                if val > THRESHOLD {
                    position = 0;
                    val = 0;
                    break;
                }
                val *= ALPHA_BASE;
            } else {
                break;
            }
            position += 1;
        }
    }
    if let Some(q) = q {
        *q = sz_string
            .offset(i64::try_from(position).map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;
    }
    Ok(val)
}

pub(crate) fn inchi_strtol(
    heap: &mut SourceHeap,
    str_: SourceConstPointer<i8>,
    p: Option<&mut SourceConstPointer<i8>>,
    base: i32,
) -> Result<i64, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2269 inchi_strtol
    // INCHI✔❌: long inchi_strtol(const char *str, const char **p, int base)
    // INCHI✔❌: {
    // INCHI✔❌:     if (base == ALPHA_BASE) {
    // INCHI✔❌:         return abctol(str, (char **)p); /* cast deliberately discards const qualifier */
    // INCHI✔❌:     }
    // INCHI✔❌:     else {
    // INCHI✔❌:         return strtol(str, (char **)p, base); /* cast deliberately discards const qualifier */
    // INCHI✔❌:     }
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_strtol

    if base == 27 {
        return abctol(heap, str_, p);
    }
    if base != 0 && !(2..=36).contains(&base) {
        heap.set_source_errno(22);
        return Ok(0);
    }

    let bytes = heap.slice(str_)?;
    let nul = bytes
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut position = 0_usize;
    while position < nul
        && matches!(
            bytes[position] as u8,
            b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'
        )
    {
        position += 1;
    }
    let negative = if position < nul && matches!(bytes[position] as u8, b'+' | b'-') {
        let negative = bytes[position] as u8 == b'-';
        position += 1;
        negative
    } else {
        false
    };

    let digit_value = |byte: u8| -> Option<u32> {
        match byte {
            b'0'..=b'9' => Some(u32::from(byte - b'0')),
            b'a'..=b'z' => Some(u32::from(byte - b'a') + 10),
            b'A'..=b'Z' => Some(u32::from(byte - b'A') + 10),
            _ => None,
        }
    };
    let mut effective_base = base;
    if effective_base == 0 {
        if position < nul && bytes[position] as u8 == b'0' {
            if position + 2 < nul
                && matches!(bytes[position + 1] as u8, b'x' | b'X')
                && digit_value(bytes[position + 2] as u8).is_some_and(|digit| digit < 16)
            {
                effective_base = 16;
                position += 2;
            } else {
                effective_base = 8;
            }
        } else {
            effective_base = 10;
        }
    } else if effective_base == 16
        && position + 2 < nul
        && bytes[position] as u8 == b'0'
        && matches!(bytes[position + 1] as u8, b'x' | b'X')
        && digit_value(bytes[position + 2] as u8).is_some_and(|digit| digit < 16)
    {
        position += 2;
    }

    let digit_start = position;
    let limit = if negative {
        (i64::MAX as u128) + 1
    } else {
        i64::MAX as u128
    };
    let radix = effective_base as u128;
    let mut magnitude = 0_u128;
    let mut overflowed = false;
    while position < nul {
        let Some(digit) = digit_value(bytes[position] as u8) else {
            break;
        };
        if digit >= effective_base as u32 {
            break;
        }
        let digit = u128::from(digit);
        if magnitude > (limit - digit) / radix {
            overflowed = true;
            magnitude = limit;
        } else if !overflowed {
            magnitude = magnitude * radix + digit;
        }
        position += 1;
    }

    if position == digit_start {
        if let Some(p) = p {
            *p = str_;
        }
        return Ok(0);
    }
    if let Some(p) = p {
        *p = str_
            .offset(i64::try_from(position).map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;
    }
    if overflowed {
        heap.set_source_errno(34);
        return Ok(if negative { i64::MIN } else { i64::MAX });
    }
    if negative {
        if magnitude == (i64::MAX as u128) + 1 {
            Ok(i64::MIN)
        } else {
            Ok(-(magnitude as i64))
        }
    } else {
        Ok(magnitude as i64)
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn inchi_strtod(
    heap: &mut SourceHeap,
    str_: SourceConstPointer<i8>,
    p: Option<&mut SourceConstPointer<i8>>,
) -> Result<f64, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2288 inchi_strtod
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
double inchi_strtod(const char *str, const char **p)
{
    return strtod( str, (char **) p );
}
    */
    // END INCHI C FUNCTION: inchi_strtod
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: inchi_strtod
    // INCHI✔️❌: GCC/Linux TARGET_API_LIB resolves strtod to the C runtime implementation.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and READ_INCHI_STRING do not alter this wrapper.
    // END INCHI ACTIVE MACRO CONFIGURATION: inchi_strtod
    // The shared source-runtime parser clones the remaining allocation and its
    // numeric token, unlike libc strtod's direct scan.

    let (value, end_offset) = source_strtod_with_end(heap, str_)?;
    if let Some(p) = p {
        *p = str_.offset(
            i64::try_from(end_offset).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
        )?;
    }
    Ok(value)
}

pub(crate) fn print_sequence_of_nums_compressing_ranges(
    heap: &mut SourceHeap,
    numbers: &[i32],
    string_buffer: &mut INCHI_IOS_STRING,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2376 print_sequence_of_nums_compressing_ranges
    // INCHI✔️❌: void print_sequence_of_nums_compressing_ranges( int n,
    // INCHI✔️❌:                                                 int *num,
    // INCHI✔️❌:                                                 INCHI_IOS_STRING *strbuf )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int k, range;
    // INCHI✔️❌:
    // INCHI✔️❌:     range = 0;
    // INCHI✔️❌:     for (k = 0; k < n - 1; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (num[k + 1] == num[k] + 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (range)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 range++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_strbuf_printf( strbuf, "%d-", num[k] );
    // INCHI✔️❌:                 range = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (range)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 range = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%d,", num[k] );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_strbuf_printf( strbuf, "%d", num[n - 1] );
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: print_sequence_of_nums_compressing_ranges

    let last = *numbers.last().ok_or(SourceHeapError::PointerOutOfBounds)?;
    let range_format =
        heap.allocate_model_storage(b"%d-\0".iter().map(|byte| *byte as i8).collect())?;
    let item_format =
        heap.allocate_model_storage(b"%d,\0".iter().map(|byte| *byte as i8).collect())?;
    let last_format =
        heap.allocate_model_storage(b"%d\0".iter().map(|byte| *byte as i8).collect())?;
    let mut range = 0_i32;
    for index in 0..numbers.len() - 1 {
        let successor = numbers[index]
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        if numbers[index + 1] == successor {
            if range != 0 {
                range = range.wrapping_add(1);
            } else {
                inchi_strbuf_printf(
                    heap,
                    Some(string_buffer),
                    range_format.as_const(),
                    &SourceVaList {
                        arguments: vec![SourceFormatArgument::Signed(i64::from(numbers[index]))],
                        ..SourceVaList::default()
                    },
                )?;
                range = 1;
            }
        } else {
            if range != 0 {
                range = 0;
            }
            inchi_strbuf_printf(
                heap,
                Some(string_buffer),
                item_format.as_const(),
                &SourceVaList {
                    arguments: vec![SourceFormatArgument::Signed(i64::from(numbers[index]))],
                    ..SourceVaList::default()
                },
            )?;
        }
    }
    inchi_strbuf_printf(
        heap,
        Some(string_buffer),
        last_format.as_const(),
        &SourceVaList {
            arguments: vec![SourceFormatArgument::Signed(i64::from(last))],
            ..SourceVaList::default()
        },
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stereo_fixture(
        heap: &mut SourceHeap,
        centers: &[u16],
        parity: &[i8],
        inverse_centers: &[u16],
        inverse_parity: &[i8],
        bond1: &[u16],
        bond2: &[u16],
        bond_parity: &[i8],
    ) -> INChI_Stereo {
        INChI_Stereo {
            nNumberOfStereoCenters: centers.len() as i32,
            nNumber: heap.allocate(centers.to_vec()).unwrap(),
            t_parity: heap.allocate(parity.to_vec()).unwrap(),
            nNumberInv: heap.allocate(inverse_centers.to_vec()).unwrap(),
            t_parityInv: heap.allocate(inverse_parity.to_vec()).unwrap(),
            nCompInv2Abs: i32::from(!inverse_centers.is_empty()),
            nNumberOfStereoBonds: bond1.len() as i32,
            nBondAtom1: heap.allocate(bond1.to_vec()).unwrap(),
            nBondAtom2: heap.allocate(bond2.to_vec()).unwrap(),
            b_parity: heap.allocate(bond_parity.to_vec()).unwrap(),
            ..INChI_Stereo::default()
        }
    }

    fn allocate_c_string(heap: &mut SourceHeap, bytes: &[u8]) -> SourceConstPointer<i8> {
        let mut allocation: Vec<i8> = bytes.iter().map(|&byte| byte as i8).collect();
        allocation.push(0);
        heap.allocate(allocation).unwrap().as_const()
    }

    fn abc_case(value: i32, delimiter: Option<&str>, length: i32) -> (i32, Vec<i8>) {
        let mut heap = SourceHeap::default();
        let output = heap.allocate_model_storage(vec![83_i8; 64]).unwrap();
        let delimiter = delimiter.map_or(SourceConstPointer::null(), |text| {
            allocate_c_string(&mut heap, text.as_bytes())
        });
        let result = MakeAbcNumber(&mut heap, output, length, delimiter, value).unwrap();
        (result, heap.slice(output.as_const()).unwrap().to_vec())
    }

    #[test]
    fn source_port__ichiprt2__makeabcnumber__line_2143() {
        for (value, expected) in [
            (0, "."),
            (1, "A"),
            (26, "Z"),
            (27, "A@"),
            (28, "Aa"),
            (729, "A@@"),
            (-28, "-Aa"),
            (i32::MAX, "Enqwltj"),
        ] {
            let (result, bytes) = abc_case(value, None, 64);
            assert_eq!(result, expected.len() as i32, "value={value}");
            assert_eq!(
                &bytes[..=expected.len()],
                &expected
                    .bytes()
                    .map(|byte| byte as i8)
                    .chain(std::iter::once(0))
                    .collect::<Vec<_>>(),
                "value={value}"
            );
            assert_eq!(bytes[expected.len() + 1], 83);
        }

        let (result, bytes) = abc_case(28, Some("/x"), 64);
        assert_eq!(result, 4);
        assert_eq!(
            &bytes[..5],
            &[b'/' as i8, b'x' as i8, b'A' as i8, b'a' as i8, 0]
        );
        let (result, bytes) = abc_case(0, Some("/"), 64);
        assert_eq!(result, 1);
        assert_eq!(&bytes[..3], &[b'/' as i8, b'.' as i8, 0]);

        for length in [i32::MIN, -1, 0, 1] {
            let (result, bytes) = abc_case(1, None, length);
            assert_eq!(result, -1);
            assert_eq!(bytes[0], 83);
        }
        let (result, bytes) = abc_case(1, Some("abc"), 3);
        assert_eq!(result, -1);
        assert_eq!(&bytes[..2], &[b'a' as i8, b'b' as i8]);
        assert_eq!(bytes[2], 83);
        let (result, bytes) = abc_case(28, None, 2);
        assert_eq!(result, -1);
        assert_eq!(bytes[0], b'a' as i8);
        assert_eq!(bytes[1], 83);

        let mut heap = SourceHeap::default();
        let output = heap.allocate_model_storage(vec![83_i8; 64]).unwrap();
        assert_eq!(
            MakeAbcNumber(&mut heap, output, 64, SourceConstPointer::null(), i32::MIN,),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], b'-' as i8);
    }

    fn dec_case(value: i32, delimiter: Option<&str>, length: i32) -> (i32, Vec<i8>) {
        let mut heap = SourceHeap::default();
        let output = heap.allocate_model_storage(vec![83_i8; 64]).unwrap();
        let delimiter = delimiter.map_or(SourceConstPointer::null(), |text| {
            allocate_c_string(&mut heap, text.as_bytes())
        });
        let result = MakeDecNumber(&mut heap, output, length, delimiter, value).unwrap();
        (result, heap.slice(output.as_const()).unwrap().to_vec())
    }

    #[test]
    fn source_port__ichiprt2__makedecnumber__line_2302() {
        for (value, expected) in [
            (0, "0"),
            (1, "1"),
            (10, "10"),
            (123456789, "123456789"),
            (-42, "-42"),
            (i32::MAX, "2147483647"),
        ] {
            let (result, bytes) = dec_case(value, None, 64);
            assert_eq!(result, expected.len() as i32, "value={value}");
            assert_eq!(
                &bytes[..=expected.len()],
                &expected
                    .bytes()
                    .map(|byte| byte as i8)
                    .chain(std::iter::once(0))
                    .collect::<Vec<_>>()
            );
            assert_eq!(bytes[expected.len() + 1], 83);
        }

        let (result, bytes) = dec_case(42, Some("/x"), 64);
        assert_eq!(result, 4);
        assert_eq!(
            &bytes[..5],
            &[b'/' as i8, b'x' as i8, b'4' as i8, b'2' as i8, 0]
        );
        let (result, bytes) = dec_case(0, Some("/"), 64);
        assert_eq!(result, 2);
        assert_eq!(&bytes[..3], &[b'/' as i8, b'0' as i8, 0]);

        for length in [i32::MIN, -1, 0, 1] {
            let (result, bytes) = dec_case(1, None, length);
            assert_eq!(result, -1);
            assert_eq!(bytes[0], 83);
        }
        let (result, bytes) = dec_case(1, Some("abc"), 3);
        assert_eq!(result, -1);
        assert_eq!(&bytes[..2], &[b'a' as i8, b'b' as i8]);
        assert_eq!(bytes[2], 83);
        let (result, bytes) = dec_case(12, None, 2);
        assert_eq!(result, -1);
        assert_eq!(bytes[0], b'2' as i8);
        assert_eq!(bytes[1], 83);

        let mut heap = SourceHeap::default();
        let output = heap.allocate_model_storage(vec![83_i8; 64]).unwrap();
        assert_eq!(
            MakeDecNumber(&mut heap, output, 64, SourceConstPointer::null(), i32::MIN,),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], b'-' as i8);
    }

    #[test]
    fn source_port__ichiprt2__makeisoatomstring__line_1645() {
        fn run(
            atoms: &[crate::source_types::INChI_IsotopicAtom],
            count: i32,
            mode: i32,
            initial_overflow: i32,
            prefix: &str,
        ) -> (Result<i32, SourceHeapError>, String, i32) {
            let mut heap = SourceHeap::default();
            let pointer = if atoms.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(atoms.to_vec()).unwrap()
            };
            let mut buffer = mult_buffer(&mut heap, prefix, 256);
            let mut overflow = initial_overflow;
            let result = MakeIsoAtomString(
                &mut heap,
                pointer.as_const(),
                count,
                &mut buffer,
                mode,
                &mut overflow,
            );
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            let output =
                String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap();
            (result, output, overflow)
        }

        let atoms = [
            crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 1,
                nIsoDifference: 1,
                nNum_T: 1,
                nNum_D: 2,
                nNum_H: 0,
            },
            crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 27,
                nIsoDifference: -2,
                nNum_T: 0,
                nNum_D: 0,
                nNum_H: 3,
            },
        ];
        assert_eq!(
            run(&atoms, 2, 0, 0, ""),
            (Ok(13), "1+0TD2,27-2H3".to_owned(), 0)
        );
        assert_eq!(
            run(&atoms, 2, CT_MODE_ABC_NUMBERS as i32, 0, ""),
            (Ok(11), "A1td2A@-2h3".to_owned(), 0)
        );
        let no_isotopes = [crate::source_types::INChI_IsotopicAtom {
            nAtomNumber: u16::MAX,
            ..crate::source_types::INChI_IsotopicAtom::default()
        }];
        assert_eq!(
            run(&no_isotopes, 1, 0, 0, "p"),
            (Ok(5), "p65535".to_owned(), 0)
        );
        assert_eq!(
            run(&no_isotopes, 1, CT_MODE_ABC_NUMBERS as i32, 0, ""),
            (Ok(5), "Chxf0".to_owned(), 0)
        );
        assert_eq!(run(&[], 0, 0, 0, "held"), (Ok(0), "held".to_owned(), 0));
        assert_eq!(run(&atoms, -1, 0, 0, "held"), (Ok(0), "held".to_owned(), 0));
        assert_eq!(run(&atoms, 2, 0, 9, "held"), (Ok(0), "held".to_owned(), 9));
        assert!(matches!(
            run(&atoms[..1], 2, 0, 0, "").0,
            Err(SourceHeapError::PointerOutOfBounds)
        ));
    }

    #[test]
    fn source_port__ichiprt2__makeisotautstring__line_1803() {
        fn run(
            groups: &[crate::source_types::INChI_IsotopicTGroup],
            count: i32,
            mode: i32,
            initial_overflow: i32,
            prefix: &str,
        ) -> (Result<i32, SourceHeapError>, String, i32) {
            let mut heap = SourceHeap::default();
            let pointer = if groups.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(groups.to_vec()).unwrap()
            };
            let mut buffer = mult_buffer(&mut heap, prefix, 256);
            let mut overflow = initial_overflow;
            let result = MakeIsoTautString(
                &mut heap,
                pointer.as_const(),
                count,
                &mut buffer,
                mode,
                &mut overflow,
            );
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            let output =
                String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap();
            (result, output, overflow)
        }

        let groups = [
            crate::source_types::INChI_IsotopicTGroup {
                nTGroupNumber: 1,
                nNum_T: 1,
                nNum_D: 2,
                nNum_H: 0,
            },
            crate::source_types::INChI_IsotopicTGroup {
                nTGroupNumber: 27,
                nNum_T: 0,
                nNum_D: 0,
                nNum_H: 3,
            },
        ];
        assert_eq!(
            run(&groups, 2, 0, 0, ""),
            (Ok(9), "1TD2,27H3".to_owned(), 0)
        );
        assert_eq!(
            run(&groups, 2, CT_MODE_ABC_NUMBERS as i32, 0, ""),
            (Ok(9), "A1t2dA@3h".to_owned(), 0)
        );
        let boundary = [crate::source_types::INChI_IsotopicTGroup {
            nTGroupNumber: u16::MAX,
            nNum_H: u16::MAX,
            ..crate::source_types::INChI_IsotopicTGroup::default()
        }];
        assert_eq!(
            run(&boundary, 1, 0, 0, "p"),
            (Ok(11), "p65535H65535".to_owned(), 0)
        );
        assert_eq!(run(&[], 0, 0, 0, "held"), (Ok(0), "held".to_owned(), 0));
        assert_eq!(
            run(&groups, -1, 0, 0, "held"),
            (Ok(0), "held".to_owned(), 0)
        );
        assert_eq!(run(&groups, 2, 0, 7, "held"), (Ok(0), "held".to_owned(), 7));
        assert!(matches!(
            run(&groups[..1], 2, 0, 0, "").0,
            Err(SourceHeapError::PointerOutOfBounds)
        ));
    }

    #[test]
    fn source_port__ichiprt2__makeisohstring__line_1924() {
        fn run(
            values: [i32; 3],
            mode: i32,
            initial_overflow: i32,
            prefix: &str,
        ) -> (i32, String, i32) {
            let mut heap = SourceHeap::default();
            let mut buffer = mult_buffer(&mut heap, prefix, 256);
            let mut overflow = initial_overflow;
            let length =
                MakeIsoHString(&mut heap, &values, &mut buffer, mode, &mut overflow).unwrap();
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            let output =
                String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap();
            (length, output, overflow)
        }

        assert_eq!(run([3, 2, 1], 0, 0, ""), (5, "TD2H3".to_owned(), 0));
        assert_eq!(
            run([3, 2, 1], CT_MODE_ABC_NUMBERS as i32, 0, "p"),
            (6, "p1t2d3h".to_owned(), 0)
        );
        assert_eq!(run([-3, -2, -1], 0, 0, ""), (9, "T-1D-2H-3".to_owned(), 0));
        assert_eq!(
            run([-3, -2, -1], CT_MODE_ABC_NUMBERS as i32, 0, ""),
            (9, "-1t-2d-3h".to_owned(), 0)
        );
        assert_eq!(run([0, 0, 0], 0, 0, "held"), (0, "held".to_owned(), 0));
        assert_eq!(
            run([i32::MAX, i32::MAX, i32::MAX], 0, 0, ""),
            (33, "T2147483647D2147483647H2147483647".to_owned(), 0)
        );
        assert_eq!(run([3, 2, 1], 0, 7, "held"), (0, "held".to_owned(), 7));
    }

    #[test]
    fn source_port__ichiprt2__makestereostring__line_2019() {
        fn run(
            atom1: Option<&[u16]>,
            atom2: Option<&[u16]>,
            parity: Option<&[i8]>,
            add_delimiter: i32,
            length: i32,
            mode: i32,
            initial_overflow: i32,
            fail_allocation: bool,
            prefix: &str,
        ) -> (Result<i32, SourceHeapError>, String, i32) {
            let mut heap = SourceHeap::default();
            let atom1 = atom1.map_or(SourceConstPointer::null(), |values| {
                heap.allocate_model_storage(values.to_vec())
                    .unwrap()
                    .as_const()
            });
            let atom2 = atom2.map_or(SourceConstPointer::null(), |values| {
                heap.allocate_model_storage(values.to_vec())
                    .unwrap()
                    .as_const()
            });
            let parity = parity.map_or(SourceConstPointer::null(), |values| {
                heap.allocate_model_storage(values.to_vec())
                    .unwrap()
                    .as_const()
            });
            let mut buffer = if fail_allocation {
                let bytes = prefix
                    .bytes()
                    .map(|byte| byte as i8)
                    .chain(std::iter::once(0))
                    .collect::<Vec<_>>();
                INCHI_IOS_STRING {
                    pStr: heap.allocate_model_storage(bytes).unwrap(),
                    nAllocatedLength: i32::try_from(prefix.len() + 1).unwrap(),
                    nUsedLength: i32::try_from(prefix.len()).unwrap(),
                    nPtr: 0,
                }
            } else {
                mult_buffer(&mut heap, prefix, 128)
            };
            let mut overflow = initial_overflow;
            if fail_allocation {
                heap.fail_after_allocations(0);
            }
            let result = MakeStereoString(
                &mut heap,
                atom1,
                atom2,
                parity,
                add_delimiter,
                length,
                &mut buffer,
                mode,
                &mut overflow,
            );
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let visible = bytes.iter().position(|byte| *byte == 0).unwrap();
            (
                result,
                String::from_utf8(bytes[..visible].iter().map(|byte| *byte as u8).collect())
                    .unwrap(),
                overflow,
            )
        }

        for add_delimiter in [0, 1, -1, i32::MAX] {
            let expected = "pre:1-2!,10-20-";
            assert_eq!(
                run(
                    Some(&[1, 10]),
                    Some(&[2, 20]),
                    Some(&[0, 1]),
                    add_delimiter,
                    2,
                    0,
                    0,
                    false,
                    "pre:",
                ),
                (Ok((expected.len() - 4) as i32), expected.to_owned(), 0)
            );
        }
        assert_eq!(
            run(None, Some(&[2, 20]), None, 0, 2, 0, 0, false, ""),
            (Ok(4), "2,20".to_owned(), 0)
        );
        assert_eq!(
            run(
                None,
                None,
                Some(&[0, 1, 2, 3, 4, 5, -1]),
                0,
                7,
                0,
                0,
                false,
                "",
            ),
            (Ok(7), "!-+u?!!".to_owned(), 0)
        );
        assert_eq!(
            run(
                Some(&[1, 27]),
                Some(&[2, 28]),
                Some(&[1, -2]),
                0,
                2,
                CT_MODE_ABC_NUMBERS as i32,
                0,
                false,
                "",
            ),
            (Ok(9), "AB1A@Aa-2".to_owned(), 0)
        );
        assert_eq!(
            run(
                Some(&[1]),
                None,
                None,
                0,
                i32::MAX,
                0,
                i32::MIN,
                false,
                "held",
            ),
            (Ok(0), "held".to_owned(), i32::MIN)
        );
        assert_eq!(
            run(Some(&[1]), None, None, 0, -1, 0, 0, false, "held"),
            (Ok(0), "held".to_owned(), 0)
        );
        assert_eq!(
            run(Some(&[1]), None, None, 0, 1, 0, 0, true, ""),
            (Ok(0), "".to_owned(), 0)
        );
        assert!(matches!(
            run(Some(&[]), None, None, 0, 1, 0, 0, false, "").0,
            Err(SourceHeapError::PointerOutOfBounds)
        ));
    }

    fn mult_buffer(heap: &mut SourceHeap, text: &str, allocated: usize) -> INCHI_IOS_STRING {
        let mut bytes = vec![87_i8; allocated];
        for (target, source) in bytes.iter_mut().zip(text.bytes()) {
            *target = source as i8;
        }
        bytes[text.len()] = 0;
        INCHI_IOS_STRING {
            pStr: heap.allocate_model_storage(bytes).unwrap(),
            nAllocatedLength: allocated as i32,
            nUsedLength: text.len() as i32,
            nPtr: 8,
        }
    }

    #[test]
    fn source_port__ichiprt2__makecrvstring__line_1308() {
        fn run(
            values: Option<&[ORIG_INFO]>,
            length: i32,
            add_delimiter: i32,
            mode: i32,
            initial_overflow: i32,
            prefix: &str,
        ) -> (Result<i32, SourceHeapError>, String, i32) {
            let mut heap = SourceHeap::default();
            let pointer = values.map_or(SourceConstPointer::null(), |items| {
                heap.allocate_model_storage(items.to_vec())
                    .unwrap()
                    .as_const()
            });
            let mut buffer = mult_buffer(&mut heap, prefix, 512);
            let mut overflow = initial_overflow;
            let result = MakeCRVString(
                &mut heap,
                pointer,
                length,
                add_delimiter,
                &mut buffer,
                mode,
                &mut overflow,
            );
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let used = usize::try_from(buffer.nUsedLength).unwrap();
            let output =
                String::from_utf8(bytes[..used].iter().map(|byte| *byte as u8).collect()).unwrap();
            (result, output, overflow)
        }

        let info = |charge, radical, valence| ORIG_INFO {
            cCharge: charge,
            cRadical: radical,
            cUnusualValence: valence,
        };
        let values = [
            info(0, 0, 0),
            info(3, 0, 0),
            info(-2, 0, 0),
            info(0, 1, 0),
            info(0, 2, 0),
            info(0, 3, 0),
            info(0, 0, 4),
            info(1, 0, 3),
            info(0, 1, 4),
            info(2, 2, 0),
            info(2, 1, 3),
            info(-128, -128, -128),
            info(127, 127, 127),
        ];
        assert_eq!(
            run(Some(&values), values.len() as i32, 1, 0, 0, "P"),
            (
                Ok(68),
                "P, 2+3,3-2,4d,5t,6u,7.4,8+1.3,9d4,10+2t,11+2d3,12-128u-128,13+127u127".to_owned(),
                0,
            )
        );
        assert_eq!(
            run(
                Some(&values),
                values.len() as i32,
                0,
                CT_MODE_ABC_NUMBERS as i32,
                0,
                "",
            ),
            (
                Ok(54),
                "B+3C-2D.dE.tF.uG4H+1.3I.d4J+2tK+2d3L-128u-128M+127u127".to_owned(),
                0,
            )
        );

        let empty = std::array::from_fn::<_, 4, _>(|_| ORIG_INFO::default());
        assert_eq!(
            run(Some(&empty), 4, 0, 0, 0, "held"),
            (Ok(0), "held".to_owned(), 0)
        );
        assert_eq!(
            run(None, 0, 1, 0, 0, "held"),
            (Ok(2), "held, ".to_owned(), 0)
        );
        assert_eq!(
            run(None, -1, 0, 0, 0, "held"),
            (Ok(0), "held".to_owned(), 0)
        );
        assert_eq!(
            run(None, i32::MAX, 1, 0, 7, "held"),
            (Ok(0), "held".to_owned(), 7)
        );
        assert!(matches!(
            run(Some(&[]), 1, 0, 0, 0, "").0,
            Err(SourceHeapError::PointerOutOfBounds)
        ));
    }

    #[test]
    fn source_port__ichiprt2__bhasoriginfo__line_369() {
        fn run(values: Option<&[ORIG_INFO]>, count: i32) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let pointer = values.map_or(SourceConstPointer::null(), |items| {
                heap.allocate_model_storage(items.to_vec())
                    .unwrap()
                    .as_const()
            });
            bHasOrigInfo(&heap, pointer, count)
        }

        let info = |charge, radical, valence| ORIG_INFO {
            cCharge: charge,
            cRadical: radical,
            cUnusualValence: valence,
        };
        assert_eq!(run(None, i32::MAX), Ok(0));
        assert_eq!(run(Some(&[]), 0), Ok(0));
        assert_eq!(run(Some(&[]), -1), Ok(0));
        assert_eq!(run(Some(&[info(0, 0, 0), info(0, 0, 0)]), 2), Ok(0));
        assert_eq!(run(Some(&[info(-128, 0, 0)]), 1), Ok(1));
        assert_eq!(run(Some(&[info(0, 127, 0)]), 1), Ok(1));
        assert_eq!(run(Some(&[info(0, 0, -128)]), 1), Ok(1));
        assert_eq!(run(Some(&[info(1, 0, 0)]), i32::MAX), Ok(1));
        assert!(matches!(
            run(Some(&[]), 1),
            Err(SourceHeapError::PointerOutOfBounds)
        ));
    }

    #[test]
    fn source_port__ichiprt2__eqloriginfo__line_387() {
        use crate::source_types::INChI_Aux;

        let info = |charge, radical, valence| ORIG_INFO {
            cCharge: charge,
            cRadical: radical,
            cUnusualValence: valence,
        };
        let mut heap = SourceHeap::default();
        let left = heap
            .allocate_model_storage(vec![info(1, 0, 0), info(0, 2, -3), info(9, 9, 9)])
            .unwrap();
        let equal = heap
            .allocate_model_storage(vec![info(1, 0, 0), info(0, 2, -3), info(-9, -9, -9)])
            .unwrap();
        let unequal_charge = heap
            .allocate_model_storage(vec![info(2, 0, 0), info(0, 2, -3)])
            .unwrap();
        let unequal_radical = heap
            .allocate_model_storage(vec![info(1, 0, 0), info(0, 1, -3)])
            .unwrap();
        let unequal_valence = heap
            .allocate_model_storage(vec![info(1, 0, 0), info(0, 2, 3)])
            .unwrap();
        let aux = |pointer, count| INChI_Aux {
            OrigInfo: pointer,
            nNumberOfAtoms: count,
            ..INChI_Aux::default()
        };
        let aux_left = aux(left, 2);
        let aux_equal = aux(equal, 2);
        assert_eq!(EqlOrigInfo(&heap, None, Some(&aux_equal)), Ok(0));
        assert_eq!(EqlOrigInfo(&heap, Some(&aux_left), None), Ok(0));
        assert_eq!(EqlOrigInfo(&heap, Some(&aux_left), Some(&aux_equal)), Ok(1));
        assert_eq!(
            EqlOrigInfo(&heap, Some(&aux_left), Some(&aux(unequal_charge, 2))),
            Ok(0)
        );
        assert_eq!(
            EqlOrigInfo(&heap, Some(&aux_left), Some(&aux(unequal_radical, 2))),
            Ok(0)
        );
        assert_eq!(
            EqlOrigInfo(&heap, Some(&aux_left), Some(&aux(unequal_valence, 2))),
            Ok(0)
        );
        assert_eq!(
            EqlOrigInfo(&heap, Some(&aux_left), Some(&aux(equal, 1))),
            Ok(0)
        );

        let empty = heap
            .allocate_model_storage(vec![ORIG_INFO::default()])
            .unwrap();
        let empty_aux = aux(empty, 1);
        assert_eq!(
            EqlOrigInfo(&heap, Some(&empty_aux), Some(&empty_aux)),
            Ok(0)
        );
        let zero_aux = aux(SourceMutPointer::null(), 0);
        assert_eq!(EqlOrigInfo(&heap, Some(&zero_aux), Some(&zero_aux)), Ok(0));
        assert_eq!(
            EqlOrigInfo(
                &heap,
                Some(&aux(left, 2)),
                Some(&aux(SourceMutPointer::null(), 2)),
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichiprt2__print_sequence_of_nums_compressing_ranges__line_2376() {
        fn run(numbers: &[i32], prefix: &str) -> Result<String, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let mut buffer = mult_buffer(&mut heap, prefix, 128);
            print_sequence_of_nums_compressing_ranges(&mut heap, numbers, &mut buffer)?;
            let used = usize::try_from(buffer.nUsedLength)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let bytes = heap
                .slice(buffer.pStr.as_const())?
                .get(..used)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            String::from_utf8(bytes.iter().map(|byte| *byte as u8).collect())
                .map_err(|_| SourceHeapError::InvalidSourceTextEncoding)
        }

        for (numbers, expected) in [
            (vec![1], "1"),
            (vec![1, 2], "1-2"),
            (vec![1, 2, 3], "1-3"),
            (vec![1, 2, 4], "1-2,4"),
            (vec![1, 3, 4, 5, 7], "1,3-5,7"),
            (vec![1, 1], "1,1"),
            (vec![3, 2, 1], "3,2,1"),
            (vec![-3, -2, -1, 1], "-3--1,1"),
            (vec![i32::MAX], "2147483647"),
            (vec![i32::MAX - 1, i32::MAX], "2147483646-2147483647"),
        ] {
            assert_eq!(run(&numbers, "").unwrap(), expected, "numbers={numbers:?}");
        }
        assert_eq!(run(&[5, 6, 9], "prefix:").unwrap(), "prefix:5-6,9");
        assert!(matches!(
            run(&[], ""),
            Err(SourceHeapError::PointerOutOfBounds)
        ));
        assert!(matches!(
            run(&[i32::MAX, i32::MAX], ""),
            Err(SourceHeapError::SourceIntegerOverflow)
        ));
    }

    fn taut_case(
        ct: Option<&[u16]>,
        length: i32,
        add_delimiter: i32,
        mode: i32,
        initial_overflow: i32,
        prefix: &str,
    ) -> (i32, String, i32) {
        let mut heap = SourceHeap::default();
        let pointer = ct.map_or(SourceMutPointer::null(), |values| {
            heap.allocate_model_storage(values.to_vec()).unwrap()
        });
        let mut buffer = mult_buffer(&mut heap, prefix, 128);
        let mut overflow = initial_overflow;
        let result = MakeTautString(
            &mut heap,
            pointer,
            length,
            add_delimiter,
            &mut buffer,
            mode,
            &mut overflow,
        )
        .unwrap();
        let bytes = &heap.slice(buffer.pStr.as_const()).unwrap()
            [..usize::try_from(buffer.nUsedLength).unwrap()];
        (
            result,
            String::from_utf8(bytes.iter().map(|byte| *byte as u8).collect()).unwrap(),
            overflow,
        )
    }

    fn h_case(
        values: Option<&[i8]>,
        length: i32,
        add_delimiter: i32,
        mode: i32,
        initial_overflow: i32,
        fail_after_allocations: Option<u64>,
        prefix: &str,
    ) -> (i32, String, i32) {
        let mut heap = SourceHeap::default();
        let pointer = values.map_or(SourceConstPointer::null(), |items| {
            heap.allocate_model_storage(items.to_vec())
                .unwrap()
                .as_const()
        });
        let mut buffer = mult_buffer(&mut heap, prefix, 128);
        let mut overflow = initial_overflow;
        if let Some(count) = fail_after_allocations {
            heap.fail_after_allocations(count);
        }
        let result = MakeHString(
            &mut heap,
            add_delimiter,
            pointer,
            length,
            &mut buffer,
            mode,
            &mut overflow,
        )
        .unwrap();
        let bytes = &heap.slice(buffer.pStr.as_const()).unwrap()
            [..usize::try_from(buffer.nUsedLength).unwrap()];
        (
            result,
            String::from_utf8(bytes.iter().map(|byte| *byte as u8).collect()).unwrap(),
            overflow,
        )
    }

    #[test]
    fn source_port__ichiprt2__makectstring__line_1069() {
        fn run(
            hydrogens: Option<&[i8]>,
            mode: i32,
            add_delimiter: i32,
            initial_overflow: i32,
        ) -> (Result<i32, SourceHeapError>, String, i32) {
            let mut heap = SourceHeap::default();
            let chain = heap.allocate_model_storage(vec![2_u16, 1, 3, 2]).unwrap();
            let hydrogen_pointer = hydrogens.map_or(SourceConstPointer::null(), |values| {
                heap.allocate_model_storage(values.to_vec())
                    .unwrap()
                    .as_const()
            });
            let mut buffer = mult_buffer(&mut heap, "", 128);
            let mut overflow = initial_overflow;
            let result = MakeCtString(
                &mut heap,
                chain,
                4,
                add_delimiter,
                hydrogen_pointer,
                3,
                &mut buffer,
                mode,
                &mut overflow,
            );
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            let output =
                String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap();
            (result, output, overflow)
        }

        let no_orphans = crate::source_types::CT_MODE_NO_ORPHANS as i32;
        assert_eq!(
            run(None, no_orphans, 0, 0),
            (Ok(8), ",2-1,3-2".to_owned(), 0)
        );
        assert_eq!(
            run(Some(&[1, 2, 0]), 0, 0, 0),
            (Ok(7), "2,1,3,2".to_owned(), 0)
        );
        assert_eq!(
            run(Some(&[1, 2, 0]), no_orphans, 0, 0),
            (Ok(8), "1H-2H2-3".to_owned(), 0)
        );
        assert_eq!(
            run(Some(&[1, 2, 0]), no_orphans, 1, 0),
            (Ok(9), ",1H-2H2-3".to_owned(), 0)
        );
        assert_eq!(
            run(Some(&[1, 2, 0]), no_orphans, 1, 9),
            (Ok(0), "".to_owned(), 9)
        );
    }

    #[test]
    fn source_port__ichiprt2__makectstringold__line_703() {
        fn run(
            values: &[u16],
            count: i32,
            add_delimiter: i32,
            mode: i32,
            initial_overflow: i32,
            prefix: &str,
        ) -> (Result<i32, SourceHeapError>, String, i32) {
            let mut heap = SourceHeap::default();
            let pointer = if values.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(values.to_vec()).unwrap()
            };
            let mut buffer = mult_buffer(&mut heap, prefix, 256);
            let mut overflow = initial_overflow;
            let result = MakeCtStringOld(
                &mut heap,
                pointer.as_const(),
                count,
                add_delimiter,
                &mut buffer,
                mode,
                &mut overflow,
            );
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            let output =
                String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap();
            (result, output, overflow)
        }

        assert_eq!(
            run(&[1, 2, 3], 3, 1, 0, 0, "pre"),
            (Ok(6), "pre,1,2,3".to_owned(), 0)
        );
        assert_eq!(
            run(&[1, 28, 0], 3, 1, CT_MODE_ABC_NUMBERS as i32, 0, "",),
            (Ok(5), ",AAa.".to_owned(), 0)
        );
        let no_orphans = crate::source_types::CT_MODE_NO_ORPHANS as i32;
        assert_eq!(
            run(&[1, 2, 1, 3, 1], 5, 0, no_orphans, 0, ""),
            (Ok(8), ",2-1,3-1".to_owned(), 0)
        );
        assert_eq!(
            run(&[1, 2, 1, 3, 1], 5, 1, no_orphans, 0, ""),
            (Ok(9), ",,2-1,3-1".to_owned(), 0)
        );
        assert_eq!(
            run(
                &[1, 2, 1, 3, 1],
                5,
                0,
                no_orphans | CT_MODE_ABC_NUMBERS as i32,
                0,
                "",
            ),
            (Ok(4), "BACA".to_owned(), 0)
        );
        assert_eq!(
            run(&[1, 2, 3], 3, 0, no_orphans, 0, "held"),
            (Ok(0), "held".to_owned(), 0)
        );
        assert_eq!(
            run(&[u16::MAX], 1, 0, 0, 0, ""),
            (Ok(5), "65535".to_owned(), 0)
        );
        assert_eq!(run(&[], 0, 1, 0, 0, ""), (Ok(1), ",".to_owned(), 0));
        assert_eq!(run(&[], -1, 0, 0, 0, "held"), (Ok(0), "held".to_owned(), 0));
        assert_eq!(run(&[1], 1, 1, 0, 7, "held"), (Ok(0), "held".to_owned(), 7));
        assert!(matches!(
            run(&[1], 2, 0, 0, 0, "").0,
            Err(SourceHeapError::PointerOutOfBounds)
        ));
    }

    #[test]
    fn source_port__ichiprt2__makemult__line_431() {
        let mut heap = SourceHeap::default();
        let delimiter = allocate_c_string(&mut heap, b"*");
        let mut buffer = mult_buffer(&mut heap, "pre", 32);
        let before = buffer.clone();
        let mut overflow = 0;
        assert_eq!(
            MakeMult(&mut heap, 1, delimiter, &mut buffer, 0, &mut overflow),
            Ok(0)
        );
        assert_eq!(buffer, before);
        overflow = 2;
        assert_eq!(
            MakeMult(&mut heap, 2, delimiter, &mut buffer, 0, &mut overflow),
            Ok(0)
        );
        assert_eq!(buffer, before);

        overflow = 0;
        assert_eq!(
            MakeMult(&mut heap, 28, delimiter, &mut buffer, 0, &mut overflow),
            Ok(3)
        );
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..7],
            &[
                b'p' as i8, b'r' as i8, b'e' as i8, b'2' as i8, b'8' as i8, b'*' as i8, 0
            ]
        );
        assert_eq!(overflow, 0);

        assert_eq!(
            MakeMult(
                &mut heap,
                28,
                delimiter,
                &mut buffer,
                CT_MODE_ABC_NUMBERS as i32,
                &mut overflow,
            ),
            Ok(3)
        );
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..10],
            &[
                b'p' as i8, b'r' as i8, b'e' as i8, b'2' as i8, b'8' as i8, b'*' as i8, b'A' as i8,
                b'a' as i8, b'*' as i8, 0
            ]
        );

        let long_delimiter = allocate_c_string(&mut heap, &vec![b'x'; 2047]);
        let used_before = buffer.nUsedLength;
        assert_eq!(
            MakeMult(&mut heap, 2, long_delimiter, &mut buffer, 0, &mut overflow),
            Ok(0)
        );
        assert_eq!(overflow, 1);
        assert_eq!(buffer.nUsedLength, used_before);

        let mut failing_heap = SourceHeap::default();
        let delimiter = allocate_c_string(&mut failing_heap, b"*");
        let mut failing = mult_buffer(&mut failing_heap, "abc", 4);
        let mut overflow = 0;
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            MakeMult(
                &mut failing_heap,
                2,
                delimiter,
                &mut failing,
                0,
                &mut overflow
            ),
            Ok(-1)
        );
        assert_eq!(overflow, 1);
        assert_eq!(failing.nUsedLength, 3);
    }

    #[test]
    fn source_port__ichiprt2__makedelim__line_476() {
        let mut heap = SourceHeap::default();
        let mut buffer = mult_buffer(&mut heap, "pre", 32);
        let before = buffer.clone();
        let mut overflow = 0;
        assert_eq!(
            MakeDelim(
                &mut heap,
                SourceConstPointer::null(),
                &mut buffer,
                &mut overflow,
            ),
            Ok(0)
        );
        assert_eq!(buffer, before);
        let empty = allocate_c_string(&mut heap, b"");
        assert_eq!(
            MakeDelim(&mut heap, empty, &mut buffer, &mut overflow),
            Ok(0)
        );
        assert_eq!(buffer, before);
        let slash = allocate_c_string(&mut heap, b"/c");
        overflow = 2;
        assert_eq!(
            MakeDelim(&mut heap, slash, &mut buffer, &mut overflow),
            Ok(0)
        );
        assert_eq!(buffer, before);

        overflow = 0;
        assert_eq!(
            MakeDelim(&mut heap, slash, &mut buffer, &mut overflow),
            Ok(2)
        );
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..6],
            &[
                b'p' as i8, b'r' as i8, b'e' as i8, b'/' as i8, b'c' as i8, 0
            ]
        );
        let percent = allocate_c_string(&mut heap, b"%%");
        assert_eq!(
            MakeDelim(&mut heap, percent, &mut buffer, &mut overflow),
            Ok(1)
        );
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..7],
            &[
                b'p' as i8, b'r' as i8, b'e' as i8, b'/' as i8, b'c' as i8, b'%' as i8, 0
            ]
        );

        let mut failing_heap = SourceHeap::default();
        let slash = allocate_c_string(&mut failing_heap, b"/c");
        let mut failing = mult_buffer(&mut failing_heap, "abc", 4);
        let mut overflow = 0;
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            MakeDelim(&mut failing_heap, slash, &mut failing, &mut overflow),
            Ok(-1)
        );
        assert_eq!(overflow, 1);
        assert_eq!(failing.nUsedLength, 3);
    }

    #[test]
    fn source_port__ichiprt2__makeeqstr__line_506() {
        for (multiplier, delimiter_text, expected_suffix) in [
            (1, "=", "="),
            (0, "=", "0="),
            (2, "=", "2="),
            (-2, "=", "-2="),
            (i32::MAX, "%%", "2147483647%%"),
        ] {
            let mut heap = SourceHeap::default();
            let delimiter = allocate_c_string(&mut heap, delimiter_text.as_bytes());
            let mut buffer = mult_buffer(&mut heap, "pre", 64);
            let mut overflow = 0;
            assert_eq!(
                MakeEqStr(&mut heap, delimiter, multiplier, &mut buffer, &mut overflow,),
                Ok(expected_suffix.len() as i32),
                "multiplier={multiplier}"
            );
            let bytes = &heap.slice(buffer.pStr.as_const()).unwrap()
                [..usize::try_from(buffer.nUsedLength).unwrap()];
            assert_eq!(
                bytes,
                format!("pre{expected_suffix}")
                    .bytes()
                    .map(|byte| byte as i8)
                    .collect::<Vec<_>>(),
                "multiplier={multiplier}"
            );
            assert_eq!(overflow, 0);
        }

        let mut heap = SourceHeap::default();
        let empty = allocate_c_string(&mut heap, b"");
        let delimiter = allocate_c_string(&mut heap, b"=");
        let mut buffer = mult_buffer(&mut heap, "pre", 32);
        let before = buffer.clone();
        let mut overflow = 0;
        assert_eq!(
            MakeEqStr(
                &mut heap,
                SourceConstPointer::null(),
                2,
                &mut buffer,
                &mut overflow,
            ),
            Ok(0)
        );
        assert_eq!(buffer, before);
        assert_eq!(
            MakeEqStr(&mut heap, empty, 2, &mut buffer, &mut overflow),
            Ok(0)
        );
        assert_eq!(buffer, before);
        overflow = 7;
        assert_eq!(
            MakeEqStr(&mut heap, delimiter, 2, &mut buffer, &mut overflow),
            Ok(0)
        );
        assert_eq!(buffer, before);
        assert_eq!(overflow, 7);

        let mut partial_heap = SourceHeap::default();
        let delimiter = allocate_c_string(&mut partial_heap, b"=");
        let mut partial = mult_buffer(&mut partial_heap, "abc", 5);
        let mut overflow = 0;
        partial_heap.fail_after_allocations(0);
        assert_eq!(
            MakeEqStr(&mut partial_heap, delimiter, 2, &mut partial, &mut overflow,),
            Ok(1)
        );
        assert_eq!(overflow, 1);
        assert_eq!(partial.nUsedLength, 4);
        assert_eq!(
            &partial_heap.slice(partial.pStr.as_const()).unwrap()[..5],
            &[b'a' as i8, b'b' as i8, b'c' as i8, b'2' as i8, 0]
        );

        let mut failing_heap = SourceHeap::default();
        let delimiter = allocate_c_string(&mut failing_heap, b"=");
        let mut failing = mult_buffer(&mut failing_heap, "abc", 4);
        let mut overflow = 0;
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            MakeEqStr(&mut failing_heap, delimiter, 2, &mut failing, &mut overflow,),
            Ok(1)
        );
        assert_eq!(overflow, 1);
        assert_eq!(failing.nUsedLength, 4);
        assert_eq!(
            &failing_heap.slice(failing.pStr.as_const()).unwrap()[..5],
            &[b'a' as i8, b'b' as i8, b'c' as i8, b'=' as i8, 0]
        );

        let mut boundary_heap = SourceHeap::default();
        let delimiter = allocate_c_string(&mut boundary_heap, b"=");
        let mut boundary = mult_buffer(&mut boundary_heap, "pre", 32);
        let mut overflow = 0;
        assert_eq!(
            MakeEqStr(
                &mut boundary_heap,
                delimiter,
                i32::MIN,
                &mut boundary,
                &mut overflow,
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(boundary.nUsedLength, 3);
        assert_eq!(overflow, 0);
    }

    #[test]
    fn source_port__ichiprt2__maketautstring__line_1111() {
        for (ct, length) in [(None, 1), (Some(&[][..]), 0), (Some(&[0][..]), 1)] {
            assert_eq!(taut_case(ct, length, 1, 0, 0, "pre"), (0, "pre".into(), 0));
        }
        assert_eq!(
            taut_case(Some(&[1]), -1, 1, 0, 0, "pre"),
            (1, "pre,".into(), 0)
        );
        assert_eq!(
            taut_case(Some(&[1, 4, 2, 0, 1, 5]), 6, 1, 0, 0, "pre"),
            (9, "pre,(H2,1,5)".into(), 0)
        );
        assert_eq!(
            taut_case(Some(&[1, 4, 2, 1, 28, 5]), 6, 0, 0, 0, ""),
            (10, "(H2-,28,5)".into(), 0)
        );
        assert_eq!(
            taut_case(
                Some(&[1, 4, 2, 1, 28, 5]),
                6,
                1,
                CT_MODE_ABC_NUMBERS as i32,
                0,
                "pre",
            ),
            (6, "pre,2-AaE".into(), 0)
        );
        assert_eq!(
            taut_case(Some(&[2, 2, 0, 1, 2, 1, 2]), 7, 0, 0, 0, ""),
            (8, "(-)(H-2)".into(), 0)
        );
        assert_eq!(
            taut_case(Some(&[1, 3, 0, 0, 0]), 5, 0, 0, 0, ""),
            (4, "(,0)".into(), 0)
        );
        assert_eq!(
            taut_case(
                Some(&[1, 3, 0, 0, 0]),
                5,
                1,
                CT_MODE_ABC_NUMBERS as i32,
                0,
                "",
            ),
            (3, ",0.".into(), 0)
        );
        assert_eq!(
            taut_case(Some(&[1, 0, 2]), 3, 0, 0, 0, ""),
            (4, "()()".into(), 0)
        );
        assert_eq!(
            taut_case(Some(&[1, 4, 2, 1, 28, 5]), 6, 1, 0, 3, "pre"),
            (0, "pre".into(), 3)
        );
    }

    #[test]
    fn source_port__ichiprt2__makehstring__line_789() {
        assert_eq!(h_case(None, 3, 1, 0, 0, None, "pre"), (1, "pre,".into(), 0));
        assert_eq!(
            h_case(Some(&[1]), 0, 1, 0, 0, None, "pre"),
            (1, "pre,".into(), 0)
        );
        assert_eq!(
            h_case(Some(&[1]), -1, 1, 0, 0, None, "pre"),
            (1, "pre,".into(), 0)
        );
        assert_eq!(
            h_case(Some(&[1, 1, 0, 2, -1, -1]), 6, 1, 0, 0, None, "pre"),
            (14, "pre,1-2H,4H2,5-6h".into(), 0)
        );
        assert_eq!(
            h_case(
                Some(&[1, 1, 0, 2, -1, -1]),
                6,
                1,
                CT_MODE_ABC_NUMBERS as i32,
                0,
                None,
                "pre",
            ),
            (9, "preAB1D2EF-1".into(), 0)
        );
        let together = crate::source_types::CT_MODE_EQL_H_TOGETHER as i32;
        assert_eq!(
            h_case(Some(&[1, 1, 0, 2, -1, -1]), 6, 0, together, 0, None, ""),
            (13, "5-6h,1-2H,4H2".into(), 0)
        );
        assert_eq!(
            h_case(Some(&[1, 0, 1]), 3, 0, together, 0, None, ""),
            (4, "1,3H".into(), 0)
        );
        assert_eq!(
            h_case(
                Some(&[1, 0, 1]),
                3,
                1,
                together | CT_MODE_ABC_NUMBERS as i32,
                0,
                None,
                "",
            ),
            (4, "A1C1".into(), 0)
        );
        assert_eq!(
            h_case(Some(&[-128, 127]), 2, 0, together, 0, None, ""),
            (11, "1h128,2H127".into(), 0)
        );
        assert_eq!(
            h_case(Some(&[0, 0]), 2, 1, together, 0, None, "pre"),
            (0, "pre,".into(), 0)
        );
        assert_eq!(
            h_case(Some(&[0, 0]), 2, 1, 0, 0, None, "pre"),
            (1, "pre,".into(), 0)
        );
        assert_eq!(
            h_case(Some(&[-128, 127]), 2, 1, together, 0, Some(0), "pre"),
            (0, "pre,".into(), 1)
        );
        assert_eq!(
            h_case(Some(&[1, 2]), 2, 1, 0, 7, None, "pre"),
            (0, "pre".into(), 7)
        );
    }

    #[test]
    fn source_port__ichiprt2__eql_inchi_stereo__line_54() {
        let mut heap = SourceHeap::default();
        assert_eq!(Eql_INChI_Stereo(&heap, None, 0, None, 0, 0), Ok(0));

        let sp2 = stereo_fixture(&mut heap, &[], &[], &[], &[], &[1, 3], &[2, 4], &[1, 2]);
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&sp2),
                EQL_SP2 as i32,
                None,
                EQL_EXISTS as i32,
                0
            ),
            Ok(1)
        );
        assert_eq!(
            Eql_INChI_Stereo(&heap, Some(&sp2), EQL_SP2 as i32, None, EQL_SP2 as i32, 0),
            Ok(0)
        );
        let sp2_equal = stereo_fixture(&mut heap, &[], &[], &[], &[], &[1, 3], &[2, 4], &[1, 2]);
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&sp2),
                EQL_SP2 as i32,
                Some(&sp2_equal),
                EQL_SP2 as i32,
                0,
            ),
            Ok(1)
        );
        for changed in [
            stereo_fixture(&mut heap, &[], &[], &[], &[], &[1, 5], &[2, 4], &[1, 2]),
            stereo_fixture(&mut heap, &[], &[], &[], &[], &[1, 3], &[2, 5], &[1, 2]),
            stereo_fixture(&mut heap, &[], &[], &[], &[], &[1, 3], &[2, 4], &[1, 1]),
            stereo_fixture(&mut heap, &[], &[], &[], &[], &[1], &[2], &[1]),
        ] {
            assert_eq!(
                Eql_INChI_Stereo(
                    &heap,
                    Some(&sp2),
                    EQL_SP2 as i32,
                    Some(&changed),
                    EQL_SP2 as i32,
                    0,
                ),
                Ok(0)
            );
        }
        let mut missing = sp2_equal.clone();
        missing.b_parity = Default::default();
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&sp2),
                EQL_SP2 as i32,
                Some(&missing),
                EQL_SP2 as i32,
                0,
            ),
            Ok(0)
        );

        let normal = stereo_fixture(&mut heap, &[1, 2], &[1, 4], &[1, 2], &[2, 4], &[], &[], &[]);
        let equal = stereo_fixture(&mut heap, &[1, 2], &[1, 4], &[1, 2], &[2, 4], &[], &[], &[]);
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&normal),
                EQL_SP3 as i32,
                Some(&equal),
                EQL_SP3 as i32,
                0,
            ),
            Ok(1)
        );
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&normal),
                EQL_SP3 as i32,
                None,
                EQL_EXISTS as i32,
                0,
            ),
            Ok(1)
        );
        let one = stereo_fixture(&mut heap, &[1], &[1], &[1], &[2], &[], &[], &[]);
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&one),
                EQL_SP3 as i32,
                None,
                EQL_EXISTS as i32,
                1,
            ),
            Ok(0)
        );

        let inverse = stereo_fixture(&mut heap, &[1, 2], &[1, 4], &[1, 2], &[2, 4], &[], &[], &[]);
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&normal),
                EQL_SP3 as i32,
                Some(&inverse),
                EQL_SP3_INV as i32,
                0,
            ),
            Ok(1)
        );
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&normal),
                EQL_SP3_INV as i32,
                Some(&inverse),
                EQL_SP3_INV as i32,
                0,
            ),
            Ok(1)
        );
        let no_actual_inversion =
            stereo_fixture(&mut heap, &[1, 2], &[4, 4], &[1, 2], &[4, 4], &[], &[], &[]);
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&no_actual_inversion),
                EQL_SP3 as i32,
                Some(&no_actual_inversion),
                EQL_SP3_INV as i32,
                0,
            ),
            Ok(0)
        );
        let bad_number =
            stereo_fixture(&mut heap, &[1, 2], &[1, 4], &[1, 3], &[2, 4], &[], &[], &[]);
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&normal),
                EQL_SP3 as i32,
                Some(&bad_number),
                EQL_SP3_INV as i32,
                0,
            ),
            Ok(0)
        );
        let mut no_inverse = normal.clone();
        no_inverse.nCompInv2Abs = 0;
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&no_inverse),
                EQL_SP3_INV as i32,
                None,
                EQL_EXISTS as i32,
                0,
            ),
            Ok(0)
        );
        let mut missing_center = normal.clone();
        missing_center.t_parity = Default::default();
        assert_eq!(
            Eql_INChI_Stereo(
                &heap,
                Some(&missing_center),
                EQL_SP3 as i32,
                Some(&equal),
                EQL_SP3 as i32,
                0,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichiprt2__eql_inchi_isotopic__line_203() {
        let mut heap = SourceHeap::default();
        assert_eq!(Eql_INChI_Isotopic(&heap, None, None), Ok(0));

        let mut empty1 = INChI::default();
        let empty2 = INChI::default();
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&empty1), Some(&empty2)),
            Ok(0)
        );
        empty1.bDeleted = 1;
        empty1.nNumberOfIsotopicAtoms = 1;
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&empty1), Some(&empty2)),
            Ok(0)
        );

        let atoms = vec![
            crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 1,
                nIsoDifference: -2,
                nNum_H: 3,
                nNum_D: 4,
                nNum_T: 5,
            },
            crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: u16::MAX,
                nIsoDifference: i16::MAX,
                nNum_H: i16::MIN,
                nNum_D: 0,
                nNum_T: 1,
            },
        ];
        let atoms1 = heap.allocate_model_storage(atoms.clone()).unwrap();
        let atoms2 = heap.allocate_model_storage(atoms.clone()).unwrap();
        let mut atom_inchi1 = INChI {
            nNumberOfIsotopicAtoms: 2,
            IsotopicAtom: atoms1,
            ..INChI::default()
        };
        let mut atom_inchi2 = INChI {
            nNumberOfIsotopicAtoms: 2,
            IsotopicAtom: atoms2,
            nTotalCharge: 99,
            ..INChI::default()
        };
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&atom_inchi1), Some(&atom_inchi2)),
            Ok(1)
        );
        atom_inchi2.nNumberOfIsotopicAtoms = 1;
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&atom_inchi1), Some(&atom_inchi2)),
            Ok(0)
        );
        atom_inchi2.nNumberOfIsotopicAtoms = 2;
        atom_inchi2.IsotopicAtom = SourceMutPointer::null();
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&atom_inchi1), Some(&atom_inchi2)),
            Ok(0)
        );
        let mut changed_atoms = atoms;
        changed_atoms[1].nNum_T = 2;
        atom_inchi2.IsotopicAtom = heap.allocate_model_storage(changed_atoms).unwrap();
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&atom_inchi1), Some(&atom_inchi2)),
            Ok(0)
        );

        let groups = vec![crate::source_types::INChI_IsotopicTGroup {
            nTGroupNumber: u16::MAX,
            nNum_H: 1,
            nNum_D: 2,
            nNum_T: 3,
        }];
        let groups1 = heap.allocate_model_storage(groups.clone()).unwrap();
        let groups2 = heap.allocate_model_storage(groups).unwrap();
        let group_inchi1 = INChI {
            nNumberOfIsotopicTGroups: 1,
            IsotopicTGroup: groups1,
            ..INChI::default()
        };
        let mut group_inchi2 = INChI {
            nNumberOfIsotopicTGroups: 1,
            IsotopicTGroup: groups2,
            ..INChI::default()
        };
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&group_inchi1), Some(&group_inchi2)),
            Ok(1)
        );
        heap.slice_mut(group_inchi2.IsotopicTGroup).unwrap()[0].nNum_D = 7;
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&group_inchi1), Some(&group_inchi2)),
            Ok(0)
        );
        group_inchi2.bDeleted = -1;
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&group_inchi1), Some(&group_inchi2)),
            Ok(0)
        );

        atom_inchi1.IsotopicAtom = heap
            .allocate_model_storage(vec![crate::source_types::INChI_IsotopicAtom::default()])
            .unwrap();
        assert_eq!(
            Eql_INChI_Isotopic(&heap, Some(&atom_inchi1), Some(&atom_inchi1)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiprt2__abctol__line_2203() {
        let mut heap = SourceHeap::default();
        let cases: &[(&[u8], i64, i64)] = &[
            (b"", 0, 0),
            (b"@tail", 0, 1),
            (b"-@tail", 0, 2),
            (b"A", 1, 1),
            (b"Z", 26, 1),
            (b"Aa", 28, 2),
            (b"Az", 53, 2),
            (b"Aa@", 756, 3),
            (b"-Aa!", 28, 3),
            (b"AA", 1, 1),
            (b"a", 0, 0),
            (b"-", 0, 0),
            (&[0x80], 0, 0),
        ];
        for &(bytes, expected_value, expected_offset) in cases {
            let input = allocate_c_string(&mut heap, bytes);
            let mut end = SourceConstPointer::null();
            assert_eq!(abctol(&heap, input, Some(&mut end)), Ok(expected_value));
            assert_eq!(end, input.offset(expected_offset).unwrap());
        }

        let overflow = allocate_c_string(&mut heap, b"Zzzzzzzzzzzzz@");
        let mut end = SourceConstPointer::null();
        assert_eq!(abctol(&heap, overflow, Some(&mut end)), Ok(0));
        assert_eq!(end, overflow);

        let input = allocate_c_string(&mut heap, b"Bc");
        assert_eq!(abctol(&heap, input, None), Ok(57));

        let unterminated = heap.allocate(vec![b'A' as i8]).unwrap().as_const();
        let mut end = SourceConstPointer::null();
        assert_eq!(
            abctol(&heap, unterminated, Some(&mut end)),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert!(end.is_null());
    }

    #[test]
    fn source_port__ichiprt2__inchi_strtol__line_2269() {
        let mut heap = SourceHeap::default();
        let cases: &[(&[u8], i32, i64, i64)] = &[
            (b"123tail", 10, 123, 3),
            (b"  -42!", 10, -42, 5),
            (b"0x10z", 0, 16, 4),
            (b"0778", 0, 63, 3),
            (b"0x", 0, 0, 1),
            (b"0xg", 16, 0, 1),
            (b"10102", 2, 10, 4),
            (b"zZ!", 36, 1_295, 2),
            (b"+", 10, 0, 0),
            (b"  nope", 10, 0, 0),
            (b"Aa@!tail", 27, 756, 3),
        ];
        for &(bytes, base, expected_value, expected_offset) in cases {
            let input = allocate_c_string(&mut heap, bytes);
            let mut end = SourceConstPointer::null();
            heap.set_source_errno(77);
            assert_eq!(
                inchi_strtol(&mut heap, input, Some(&mut end), base),
                Ok(expected_value)
            );
            assert_eq!(end, input.offset(expected_offset).unwrap());
            assert_eq!(heap.source_errno(), 77);
        }

        let invalid = allocate_c_string(&mut heap, b"123");
        let sentinel = invalid.offset(2).unwrap();
        let mut end = sentinel;
        heap.set_source_errno(77);
        assert_eq!(inchi_strtol(&mut heap, invalid, Some(&mut end), 1), Ok(0));
        assert_eq!(end, sentinel);
        assert_eq!(heap.source_errno(), 22);

        let positive_overflow = allocate_c_string(&mut heap, b"9223372036854775808x");
        let mut end = SourceConstPointer::null();
        heap.set_source_errno(0);
        assert_eq!(
            inchi_strtol(&mut heap, positive_overflow, Some(&mut end), 10),
            Ok(i64::MAX)
        );
        assert_eq!(end, positive_overflow.offset(19).unwrap());
        assert_eq!(heap.source_errno(), 34);

        let negative_overflow = allocate_c_string(&mut heap, b"-9223372036854775809x");
        let mut end = SourceConstPointer::null();
        heap.set_source_errno(0);
        assert_eq!(
            inchi_strtol(&mut heap, negative_overflow, Some(&mut end), 10),
            Ok(i64::MIN)
        );
        assert_eq!(end, negative_overflow.offset(20).unwrap());
        assert_eq!(heap.source_errno(), 34);

        let exact_min = allocate_c_string(&mut heap, b"-9223372036854775808");
        heap.set_source_errno(9);
        assert_eq!(inchi_strtol(&mut heap, exact_min, None, 10), Ok(i64::MIN));
        assert_eq!(heap.source_errno(), 9);

        let unterminated = heap.allocate(vec![b'1' as i8]).unwrap().as_const();
        assert_eq!(
            inchi_strtol(&mut heap, unterminated, None, 10),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichiprt2__inchi_strtod__line_2288() {
        let mut heap = SourceHeap::default();
        let finite_cases: &[(&[u8], f64, i64)] = &[
            (b"123.5tail", 123.5, 5),
            (b" \t-0.25e+2!", -25.0, 10),
            (b".5x", 0.5, 2),
            (b"1.e2x", 100.0, 4),
            (b"1e+x", 1.0, 1),
            (b"nope", 0.0, 0),
            (b"0x1.8p+1z", 3.0, 8),
            (b"0x1pz", 1.0, 3),
        ];
        for &(bytes, expected, expected_offset) in finite_cases {
            let input = allocate_c_string(&mut heap, bytes);
            let mut end = SourceConstPointer::null();
            heap.set_source_errno(77);
            let actual = inchi_strtod(&mut heap, input, Some(&mut end)).unwrap();
            assert_eq!(actual.to_bits(), expected.to_bits(), "{bytes:?}");
            assert_eq!(end, input.offset(expected_offset).unwrap(), "{bytes:?}");
            assert_eq!(heap.source_errno(), 77, "{bytes:?}");
        }

        let positive_infinity = allocate_c_string(&mut heap, b"INFz");
        let mut end = SourceConstPointer::null();
        assert_eq!(
            inchi_strtod(&mut heap, positive_infinity, Some(&mut end)),
            Ok(f64::INFINITY)
        );
        assert_eq!(end, positive_infinity.offset(3).unwrap());

        let negative_infinity = allocate_c_string(&mut heap, b"-infinity!");
        assert_eq!(
            inchi_strtod(&mut heap, negative_infinity, Some(&mut end)),
            Ok(f64::NEG_INFINITY)
        );
        assert_eq!(end, negative_infinity.offset(9).unwrap());

        let nan = allocate_c_string(&mut heap, b"NAN(payload)!");
        let value = inchi_strtod(&mut heap, nan, Some(&mut end)).unwrap();
        assert!(value.is_nan());
        assert!(!value.is_sign_negative());
        assert_eq!(end, nan.offset(12).unwrap());

        let numeric_nan = allocate_c_string(&mut heap, b"nan(123)!");
        let value = inchi_strtod(&mut heap, numeric_nan, Some(&mut end)).unwrap();
        assert_eq!(value.to_bits(), 0x7ff8_0000_0000_007b);
        assert_eq!(end, numeric_nan.offset(8).unwrap());

        let hexadecimal_nan = allocate_c_string(&mut heap, b"nan(0x123)!");
        let value = inchi_strtod(&mut heap, hexadecimal_nan, Some(&mut end)).unwrap();
        assert_eq!(value.to_bits(), 0x7ff8_0000_0000_0123);
        assert_eq!(end, hexadecimal_nan.offset(10).unwrap());

        let octal_nan = allocate_c_string(&mut heap, b"nan(0123)!");
        let value = inchi_strtod(&mut heap, octal_nan, Some(&mut end)).unwrap();
        assert_eq!(value.to_bits(), 0x7ff8_0000_0000_0053);
        assert_eq!(end, octal_nan.offset(9).unwrap());

        let invalid_nan_payload = allocate_c_string(&mut heap, b"nan(x-y)!");
        let value = inchi_strtod(&mut heap, invalid_nan_payload, Some(&mut end)).unwrap();
        assert_eq!(value.to_bits(), f64::NAN.to_bits());
        assert_eq!(end, invalid_nan_payload.offset(3).unwrap());

        let negative_nan = allocate_c_string(&mut heap, b"-nan!");
        let value = inchi_strtod(&mut heap, negative_nan, Some(&mut end)).unwrap();
        assert!(value.is_nan());
        assert!(value.is_sign_negative());
        assert_eq!(end, negative_nan.offset(4).unwrap());

        let negative_zero = allocate_c_string(&mut heap, b"-0tail");
        let value = inchi_strtod(&mut heap, negative_zero, Some(&mut end)).unwrap();
        assert_eq!(value.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(end, negative_zero.offset(2).unwrap());

        let overflow = allocate_c_string(&mut heap, b"1e309x");
        heap.set_source_errno(0);
        assert_eq!(
            inchi_strtod(&mut heap, overflow, Some(&mut end)),
            Ok(f64::INFINITY)
        );
        assert_eq!(end, overflow.offset(5).unwrap());
        assert_eq!(heap.source_errno(), 34);

        let underflow = allocate_c_string(&mut heap, b"1e-4000x");
        heap.set_source_errno(0);
        let value = inchi_strtod(&mut heap, underflow, Some(&mut end)).unwrap();
        assert_eq!(value.to_bits(), 0.0_f64.to_bits());
        assert_eq!(end, underflow.offset(7).unwrap());
        assert_eq!(heap.source_errno(), 34);

        let subnormal = allocate_c_string(&mut heap, b"5e-324x");
        heap.set_source_errno(0);
        let value = inchi_strtod(&mut heap, subnormal, Some(&mut end)).unwrap();
        assert_eq!(value.to_bits(), 1);
        assert_eq!(end, subnormal.offset(6).unwrap());
        assert_eq!(heap.source_errno(), 34);

        let no_end = allocate_c_string(&mut heap, b"2.5");
        assert_eq!(inchi_strtod(&mut heap, no_end, None), Ok(2.5));

        let unterminated = heap.allocate(vec![b'1' as i8]).unwrap().as_const();
        assert_eq!(
            inchi_strtod(&mut heap, unterminated, Some(&mut end)),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            inchi_strtod(&mut heap, SourceConstPointer::null(), Some(&mut end)),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichiprt2__makectstringnew__line_551() {
        fn text(heap: &SourceHeap, buffer: &INCHI_IOS_STRING) -> Vec<u8> {
            heap.slice(buffer.pStr.as_const()).unwrap()[..buffer.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        let mut heap = SourceHeap::default();
        let chain = heap.allocate_model_storage(vec![2_u16, 1, 3, 2]).unwrap();
        let hydrogens = heap.allocate_model_storage(vec![1_i8, 2, 0]).unwrap();

        let mut buffer = mult_buffer(&mut heap, "", 128);
        let mut overflow = 0;
        assert_eq!(
            MakeCtStringNew(
                &mut heap,
                chain,
                4,
                0,
                SourceConstPointer::null(),
                3,
                &mut buffer,
                0,
                &mut overflow,
            ),
            Ok(5)
        );
        assert_eq!(text(&heap, &buffer), b"1-2-3");
        assert_eq!(overflow, 0);

        buffer = mult_buffer(&mut heap, "p", 128);
        assert_eq!(
            MakeCtStringNew(
                &mut heap,
                chain,
                4,
                1,
                hydrogens.as_const(),
                3,
                &mut buffer,
                0,
                &mut overflow,
            ),
            Ok(9)
        );
        assert_eq!(text(&heap, &buffer), b"p,1H-2H2-3");

        buffer = mult_buffer(&mut heap, "", 128);
        assert_eq!(
            MakeCtStringNew(
                &mut heap,
                chain,
                4,
                0,
                SourceConstPointer::null(),
                3,
                &mut buffer,
                CT_MODE_ABC_NUMBERS as i32,
                &mut overflow,
            ),
            Ok(3)
        );
        assert_eq!(text(&heap, &buffer), b"ABC");

        buffer = mult_buffer(&mut heap, "", 128);
        assert_eq!(
            MakeCtStringNew(
                &mut heap,
                chain,
                4,
                0,
                SourceConstPointer::null(),
                3,
                &mut buffer,
                CT_MODE_PREDECESSORS as i32,
                &mut overflow,
            ),
            Ok(3)
        );
        assert_eq!(text(&heap, &buffer), b"1,2");

        buffer = mult_buffer(&mut heap, "x", 128);
        overflow = 7;
        let before = buffer.clone();
        assert_eq!(
            MakeCtStringNew(
                &mut heap,
                chain,
                4,
                1,
                SourceConstPointer::null(),
                3,
                &mut buffer,
                0,
                &mut overflow,
            ),
            Ok(0)
        );
        assert_eq!(buffer, before);
        assert_eq!(overflow, 7);

        overflow = 0;
        assert_eq!(
            MakeCtStringNew(
                &mut heap,
                chain,
                1,
                1,
                SourceConstPointer::null(),
                3,
                &mut buffer,
                0,
                &mut overflow,
            ),
            Ok(0)
        );
        let mut failing_heap = SourceHeap::default();
        let chain = failing_heap
            .allocate_model_storage(vec![2_u16, 1, 3, 2])
            .unwrap();
        let mut buffer = mult_buffer(&mut failing_heap, "", 128);
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            MakeCtStringNew(
                &mut failing_heap,
                chain,
                4,
                0,
                SourceConstPointer::null(),
                3,
                &mut buffer,
                0,
                &mut overflow,
            ),
            Ok(0)
        );
        assert_eq!(overflow, 1);
        assert_eq!(buffer.nUsedLength, 0);
    }

    #[test]
    fn source_port__ichiprt2__eql_inchi_aux_num__line_305() {
        use crate::source_types::INChI_Aux;

        let mut heap = SourceHeap::default();
        let ordinary = heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let isotopic = heap.allocate(vec![3_u16, 2, 1]).unwrap();
        let inverted = heap.allocate(vec![2_u16, 1, 3]).unwrap();
        let isotopic_inverted = heap.allocate(vec![2_u16, 3, 1]).unwrap();
        let same_ordinary = heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let same_isotopic = heap.allocate(vec![3_u16, 2, 1]).unwrap();
        let same_inverted = heap.allocate(vec![2_u16, 1, 3]).unwrap();
        let same_isotopic_inverted = heap.allocate(vec![2_u16, 3, 1]).unwrap();
        let aux1 = INChI_Aux {
            nNumberOfAtoms: 3,
            bIsIsotopic: 1,
            nOrigAtNosInCanonOrd: ordinary,
            nIsotopicOrigAtNosInCanonOrd: isotopic,
            nOrigAtNosInCanonOrdInv: inverted,
            nIsotopicOrigAtNosInCanonOrdInv: isotopic_inverted,
            ..INChI_Aux::default()
        };
        let aux2 = INChI_Aux {
            nNumberOfAtoms: 3,
            bIsIsotopic: 1,
            nOrigAtNosInCanonOrd: same_ordinary,
            nIsotopicOrigAtNosInCanonOrd: same_isotopic,
            nOrigAtNosInCanonOrdInv: same_inverted,
            nIsotopicOrigAtNosInCanonOrdInv: same_isotopic_inverted,
            ..INChI_Aux::default()
        };

        for selector in [
            EQL_NUM as i32,
            EQL_NUM_ISO as i32,
            EQL_NUM_INV as i32,
            (EQL_NUM_INV | EQL_NUM_ISO) as i32,
        ] {
            assert_eq!(
                Eql_INChI_Aux_Num(&heap, Some(&aux1), selector, Some(&aux2), selector),
                Ok(1),
                "selector={selector}"
            );
        }
        assert_eq!(
            Eql_INChI_Aux_Num(
                &heap,
                Some(&aux1),
                EQL_NUM as i32,
                Some(&aux2),
                EQL_NUM_ISO as i32,
            ),
            Ok(0)
        );
        assert_eq!(Eql_INChI_Aux_Num(&heap, None, 0, Some(&aux2), 0), Ok(0));
        assert_eq!(
            Eql_INChI_Aux_Num(&heap, Some(&aux1), 99, Some(&aux2), 0),
            Ok(0)
        );

        let non_isotopic = INChI_Aux {
            bIsIsotopic: 0,
            ..aux1.clone()
        };
        assert_eq!(
            Eql_INChI_Aux_Num(
                &heap,
                Some(&non_isotopic),
                EQL_NUM_ISO as i32,
                Some(&aux2),
                EQL_NUM_ISO as i32,
            ),
            Ok(0)
        );
        let deleted = INChI_Aux {
            bDeleted: 1,
            ..aux1.clone()
        };
        assert_eq!(
            Eql_INChI_Aux_Num(&heap, Some(&deleted), 0, Some(&aux2), 0),
            Ok(0)
        );
        let wrong_length = INChI_Aux {
            nNumberOfAtoms: 2,
            ..aux2.clone()
        };
        assert_eq!(
            Eql_INChI_Aux_Num(&heap, Some(&aux1), 0, Some(&wrong_length), 0),
            Ok(0)
        );
        let null_numbers = INChI_Aux {
            nOrigAtNosInCanonOrd: SourceMutPointer::null(),
            ..aux2.clone()
        };
        assert_eq!(
            Eql_INChI_Aux_Num(&heap, Some(&aux1), 0, Some(&null_numbers), 0),
            Ok(0)
        );
        let short = heap.allocate(vec![1_u16, 2]).unwrap();
        let short_aux = INChI_Aux {
            nOrigAtNosInCanonOrd: short,
            ..aux2
        };
        assert_eq!(
            Eql_INChI_Aux_Num(&heap, Some(&aux1), 0, Some(&short_aux), 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiprt2__makeequstring__line_1559() {
        fn text(heap: &SourceHeap, buffer: &INCHI_IOS_STRING) -> Vec<u8> {
            let used = usize::try_from(buffer.nUsedLength).unwrap();
            heap.slice(buffer.pStr.as_const())
                .unwrap()
                .get(..used)
                .unwrap()
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        let mut heap = SourceHeap::default();
        let linear = heap
            .allocate_model_storage(vec![1_u16, 1_u16, 3_u16])
            .unwrap();
        let mut buffer = mult_buffer(&mut heap, "", 8);
        let mut overflow = 0;
        assert_eq!(
            MakeEquString(
                &mut heap,
                linear.as_const(),
                3,
                0,
                &mut buffer,
                0,
                &mut overflow
            ),
            Ok(8)
        );
        assert_eq!(text(&heap, &buffer), b"(1,2)(3)");
        assert_eq!(overflow, 0);

        buffer = mult_buffer(&mut heap, "seed", 8);
        overflow = 0;
        assert_eq!(
            MakeEquString(
                &mut heap,
                linear.as_const(),
                3,
                1,
                &mut buffer,
                0,
                &mut overflow
            ),
            Ok(10)
        );
        assert_eq!(text(&heap, &buffer), b"seed, (1,2)(3)");
        assert_eq!(overflow, 0);

        buffer = mult_buffer(&mut heap, "", 8);
        overflow = 0;
        assert_eq!(
            MakeEquString(
                &mut heap,
                linear.as_const(),
                3,
                0,
                &mut buffer,
                CT_MODE_ABC_NUMBERS as i32,
                &mut overflow,
            ),
            Ok(6)
        );
        assert_eq!(text(&heap, &buffer), b"AB),C)");
        assert_eq!(overflow, 0);

        buffer = mult_buffer(&mut heap, "held", 8);
        overflow = 1;
        let before = text(&heap, &buffer);
        assert_eq!(
            MakeEquString(
                &mut heap,
                linear.as_const(),
                3,
                1,
                &mut buffer,
                0,
                &mut overflow
            ),
            Ok(0)
        );
        assert_eq!(text(&heap, &buffer), before);
        assert_eq!(overflow, 1);

        buffer = mult_buffer(&mut heap, "", 8);
        overflow = 0;
        assert_eq!(
            MakeEquString(
                &mut heap,
                linear.as_const(),
                0,
                1,
                &mut buffer,
                0,
                &mut overflow
            ),
            Ok(2)
        );
        assert_eq!(text(&heap, &buffer), b", ");
        assert_eq!(overflow, 0);

        buffer = mult_buffer(&mut heap, "", 8);
        overflow = 0;
        assert_eq!(
            MakeEquString(
                &mut heap,
                linear.as_const(),
                -1,
                1,
                &mut buffer,
                0,
                &mut overflow
            ),
            Ok(2)
        );
        assert_eq!(text(&heap, &buffer), b", ");
        assert_eq!(overflow, 0);
    }

    #[test]
    fn source_port__ichiprt2__bhasequstring__line_398() {
        let mut heap = SourceHeap::default();
        let duplicate = heap.allocate_model_storage(vec![1_u16, 1, 3, 4]).unwrap();
        let identity = heap.allocate_model_storage(vec![1_u16, 2, 3, 4]).unwrap();
        let no_representative = heap.allocate_model_storage(vec![2_u16, 2, 4, 4]).unwrap();
        let late_duplicate = heap.allocate_model_storage(vec![1_u16, 2, 3, 3]).unwrap();

        assert_eq!(bHasEquString(&heap, SourceConstPointer::null(), 4), Ok(0));
        assert_eq!(bHasEquString(&heap, duplicate.as_const(), 4), Ok(1));
        assert_eq!(bHasEquString(&heap, identity.as_const(), 4), Ok(0));
        assert_eq!(bHasEquString(&heap, no_representative.as_const(), 4), Ok(0));
        assert_eq!(bHasEquString(&heap, late_duplicate.as_const(), 4), Ok(1));
        assert_eq!(bHasEquString(&heap, duplicate.as_const(), 0), Ok(0));
        assert_eq!(bHasEquString(&heap, duplicate.as_const(), -1), Ok(0));

        let short = heap.allocate_model_storage(vec![1_u16]).unwrap();
        assert_eq!(
            bHasEquString(&heap, short.as_const(), 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiprt2__eql_inchi_aux_equ__line_227() {
        use crate::source_types::INChI_Aux;

        let mut heap = SourceHeap::default();
        let ordinary1 = heap.allocate_model_storage(vec![1_u16, 1, 3]).unwrap();
        let ordinary2 = heap.allocate_model_storage(vec![1_u16, 1, 3]).unwrap();
        let isotopic1 = heap.allocate_model_storage(vec![1_u16, 2, 2]).unwrap();
        let isotopic2 = heap.allocate_model_storage(vec![1_u16, 2, 2]).unwrap();
        let tgroups1 = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let tgroups2 = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let mut first = INChI_Aux {
            nNumberOfAtoms: 3,
            nNumberOfTGroups: 2,
            bIsIsotopic: 1,
            nConstitEquNumbers: ordinary1,
            nConstitEquIsotopicNumbers: isotopic1,
            nConstitEquTGroupNumbers: tgroups1,
            ..INChI_Aux::default()
        };
        let mut second = INChI_Aux {
            nNumberOfAtoms: 3,
            nNumberOfTGroups: 2,
            bIsIsotopic: 1,
            nConstitEquNumbers: ordinary2,
            nConstitEquIsotopicNumbers: isotopic2,
            nConstitEquTGroupNumbers: tgroups2,
            ..INChI_Aux::default()
        };

        assert_eq!(Eql_INChI_Aux_Equ(&heap, None, 0, Some(&second), 0), Ok(0));
        assert_eq!(
            Eql_INChI_Aux_Equ(&heap, Some(&first), 0, Some(&second), 0),
            Ok(1)
        );
        assert_eq!(
            Eql_INChI_Aux_Equ(
                &heap,
                Some(&first),
                EQL_EQU_ISO as i32,
                Some(&second),
                EQL_EQU_ISO as i32,
            ),
            Ok(1)
        );
        assert_eq!(
            Eql_INChI_Aux_Equ(
                &heap,
                Some(&first),
                EQL_EQU_TG as i32,
                Some(&second),
                EQL_EQU_TG as i32,
            ),
            Ok(1)
        );
        assert_eq!(
            Eql_INChI_Aux_Equ(&heap, Some(&first), 0, Some(&second), EQL_EQU_TG as i32),
            Ok(0)
        );

        let identity = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        second.nConstitEquNumbers = identity;
        assert_eq!(
            Eql_INChI_Aux_Equ(&heap, Some(&first), 0, Some(&second), 0),
            Ok(0)
        );
        second.nConstitEquNumbers = ordinary2;
        second.bIsIsotopic = 0;
        assert_eq!(
            Eql_INChI_Aux_Equ(
                &heap,
                Some(&first),
                EQL_EQU_ISO as i32,
                Some(&second),
                EQL_EQU_ISO as i32,
            ),
            Ok(0)
        );
        second.bIsIsotopic = 1;
        first.bDeleted = 1;
        assert_eq!(
            Eql_INChI_Aux_Equ(&heap, Some(&first), 0, Some(&second), 0),
            Ok(0)
        );
        first.bDeleted = 0;
        second.nNumberOfAtoms = 2;
        assert_eq!(
            Eql_INChI_Aux_Equ(&heap, Some(&first), 0, Some(&second), 0),
            Ok(0)
        );
        second.nNumberOfAtoms = 3;
        second.nConstitEquNumbers = SourceMutPointer::null();
        assert_eq!(
            Eql_INChI_Aux_Equ(&heap, Some(&first), 0, Some(&second), 0),
            Ok(0)
        );
    }
}
