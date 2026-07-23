use crate::source::base::ichi_io::inchi_strbuf_printf;
use crate::source::base::util::{get_element_or_pseudoelement_symbol, inchi_calloc, inchi_free};
use crate::source_types::{
    AT_NUMB, INCHI_IOS_STRING, INChI, S_CHAR, SourceConstPointer, SourceFormatArgument, SourceHeap,
    SourceHeapError, SourceMutPointer, SourceVaList, U_CHAR,
};

#[allow(non_snake_case)]
pub(crate) fn MakeHillFormulaString(
    heap: &mut SourceHeap,
    hill_formula: SourceConstPointer<i8>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:121 MakeHillFormulaString
    // INCHI✔️❌: int MakeHillFormulaString( char *szHillFormula,
    // INCHI✔️❌:                            INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                            int *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:     if (szHillFormula && !*bOverflow)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (-1 == inchi_strbuf_printf( strbuf, "%s", szHillFormula ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *bOverflow |= 1;
    // INCHI✔️❌:             return  nUsedLength0 + 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeHillFormulaString

    let initial_used_length = string_buffer.nUsedLength;
    if !hill_formula.is_null() && *overflow == 0 {
        let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
        let result = inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(hill_formula)],
                ..SourceVaList::default()
            },
        );
        match result {
            Ok(-1) | Err(SourceHeapError::AllocationFailed) => {
                *overflow |= 1;
                return Ok(initial_used_length.wrapping_add(1));
            }
            Ok(_) => {}
            Err(error) => return Err(error),
        }
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_used_length))
}

#[allow(non_snake_case)]
pub(crate) fn GetHillFormulaIndexLength(count: i32) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:144 GetHillFormulaIndexLength
    // INCHI✔️❌: int GetHillFormulaIndexLength( int count )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szCount[16];
    // INCHI✔️❌:     if (count > 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return sprintf(szCount, "%d", count);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetHillFormulaIndexLength

    if count > 1 {
        count.to_string().len() as i32
    } else {
        0
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetHillFormulaCounts(
    atom: &[U_CHAR],
    num_hydrogens: &[S_CHAR],
    num_atoms: i32,
    tautomer: Option<&[AT_NUMB]>,
    tautomer_length: i32,
    num_carbons: &mut i32,
    num_hydrogens_total: &mut i32,
    formula_length: &mut i32,
    num_non_hydrogen_atoms: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:159 GetHillFormulaCounts
    // INCHI✔️❌: int GetHillFormulaCounts( U_CHAR *nAtom,
    // INCHI✔️❌:                           S_CHAR *nNum_H,
    // INCHI✔️❌:                           int num_atoms,
    // INCHI✔️❌:                           AT_NUMB *nTautomer,
    // INCHI✔️❌:                           int lenTautomer,
    // INCHI✔️❌:                           int *pnum_C,
    // INCHI✔️❌:                           int *pnum_H,
    // INCHI✔️❌:                           int *pnLen,
    // INCHI✔️❌:                           int *pnNumNonHAtoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szElement[4];
    // INCHI✔️❌:     U_CHAR nPrevAtom = (U_CHAR) -2;
    // INCHI✔️❌:     int  bCarbon, bHydrogen, nElemLen, nFormLen, nNumNonHAtoms;
    // INCHI✔️❌:     int  mult, i, num_H, num_C;
    // INCHI✔️❌:
    // INCHI✔️❌:     num_H = 0;
    // INCHI✔️❌:     num_C = 0;
    // INCHI✔️❌:     bCarbon = 0;
    // INCHI✔️❌:     bHydrogen = 0;
    // INCHI✔️❌:     nElemLen = 0;
    // INCHI✔️❌:     nFormLen = 0;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     nNumNonHAtoms = num_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (nPrevAtom != nAtom[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (mult)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bHydrogen)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_H += mult;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (bCarbon)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_C += mult;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nFormLen += nElemLen;
    // INCHI✔️❌:                     nFormLen += GetHillFormulaIndexLength( mult );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /*if (-1 == get_element_chemical_symbol( (int) nAtom[i], szElement ))*/
    // INCHI✔️❌:             if (-1 == get_element_or_pseudoelement_symbol( (int) nAtom[i], szElement ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             mult = 1;
    // INCHI✔️❌:             nElemLen = (int) strlen( szElement );
    // INCHI✔️❌:             nPrevAtom = nAtom[i];
    // INCHI✔️❌:             bCarbon = !strcmp( szElement, "C" );
    // INCHI✔️❌:             bHydrogen = !strcmp( szElement, "H" );
    // INCHI✔️❌:             if (bHydrogen)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nNumNonHAtoms = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         num_H += nNum_H[i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* NumGroups; ((NumAt+1, NumH, At1..AtNumAt),...) */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nTautomer && lenTautomer > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int num_groups = nTautomer[0];
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 1; i < lenTautomer && num_groups > 0; i += nTautomer[i] + 1, num_groups--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_H += nTautomer[i + 1];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (mult)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bHydrogen)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_H += mult;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (bCarbon)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_C += mult;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nFormLen += nElemLen;
    // INCHI✔️❌:             nFormLen += GetHillFormulaIndexLength( mult );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_C)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nFormLen += (int) strlen( "C" );
    // INCHI✔️❌:         nFormLen += GetHillFormulaIndexLength( num_C );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_H)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nFormLen += (int) strlen( "H" );
    // INCHI✔️❌:         nFormLen += GetHillFormulaIndexLength( num_H );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *pnum_C = num_C;
    // INCHI✔️❌:     *pnum_H = num_H;
    // INCHI✔️❌:     *pnLen = nFormLen;
    // INCHI✔️❌:     *pnNumNonHAtoms = nNumNonHAtoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetHillFormulaCounts

    let mut local_num_h = 0_i32;
    let mut local_num_c = 0_i32;
    let mut is_carbon = false;
    let mut is_hydrogen = false;
    let mut element_length = 0_i32;
    let mut local_formula_length = 0_i32;
    let mut multiplicity = 0_i32;
    let mut local_num_non_hydrogen_atoms = num_atoms;
    let mut previous_atom = (-2_i32) as U_CHAR;

    let atom_count = if num_atoms > 0 {
        usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    if atom.len() < atom_count || num_hydrogens.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    for index in 0..atom_count {
        if previous_atom != atom[index] {
            if multiplicity != 0 {
                if is_hydrogen {
                    local_num_h = local_num_h.wrapping_add(multiplicity);
                } else if is_carbon {
                    local_num_c = local_num_c.wrapping_add(multiplicity);
                } else {
                    local_formula_length = local_formula_length.wrapping_add(element_length);
                    local_formula_length =
                        local_formula_length.wrapping_add(GetHillFormulaIndexLength(multiplicity));
                }
            }
            let mut symbol = [0_i8; 4];
            if get_element_or_pseudoelement_symbol(i32::from(atom[index]), &mut symbol)? == -1 {
                return Ok(-1);
            }
            multiplicity = 1;
            element_length = symbol.iter().position(|byte| *byte == 0).unwrap_or(4) as i32;
            previous_atom = atom[index];
            is_carbon = symbol[..element_length as usize] == [b'C' as i8];
            is_hydrogen = symbol[..element_length as usize] == [b'H' as i8];
            if is_hydrogen {
                local_num_non_hydrogen_atoms = index as i32;
            }
        } else {
            multiplicity = multiplicity.wrapping_add(1);
        }
        local_num_h = local_num_h.wrapping_add(i32::from(num_hydrogens[index]));
    }

    if let Some(tautomer) = tautomer.filter(|_| tautomer_length > 0) {
        let active_length =
            usize::try_from(tautomer_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if tautomer.len() < active_length || active_length == 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let mut groups = i32::from(tautomer[0]);
        let mut index = 1_usize;
        while index < active_length && groups > 0 {
            let hydrogen_index = index
                .checked_add(1)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            let hydrogen_count = *tautomer
                .get(hydrogen_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            local_num_h = local_num_h.wrapping_add(i32::from(hydrogen_count));
            index = index
                .checked_add(usize::from(tautomer[index]) + 1)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            groups -= 1;
        }
    }

    if multiplicity != 0 {
        if is_hydrogen {
            local_num_h = local_num_h.wrapping_add(multiplicity);
        } else if is_carbon {
            local_num_c = local_num_c.wrapping_add(multiplicity);
        } else {
            local_formula_length = local_formula_length.wrapping_add(element_length);
            local_formula_length =
                local_formula_length.wrapping_add(GetHillFormulaIndexLength(multiplicity));
        }
    }
    if local_num_c != 0 {
        local_formula_length = local_formula_length.wrapping_add(1);
        local_formula_length =
            local_formula_length.wrapping_add(GetHillFormulaIndexLength(local_num_c));
    }
    if local_num_h != 0 {
        local_formula_length = local_formula_length.wrapping_add(1);
        local_formula_length =
            local_formula_length.wrapping_add(GetHillFormulaIndexLength(local_num_h));
    }

    *num_carbons = local_num_c;
    *num_hydrogens_total = local_num_h;
    *formula_length = local_formula_length;
    *num_non_hydrogen_atoms = local_num_non_hydrogen_atoms;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn AddElementAndCount(
    element: &[i8],
    multiplicity: i32,
    linear_ct: &mut [i8],
    linear_ct_length: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:279 AddElementAndCount
    // INCHI✔️❌: int AddElementAndCount( const char *szElement, int mult, char *szLinearCT, int nLenLinearCT, int *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szMult[16];
    // INCHI✔️❌:     int len1, len2;
    // INCHI✔️❌:     if (mult > 0 && !*bOverflow && 0 < ( len1 = strlen( szElement ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (mult > 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len2 = sprintf(szMult, "%d", mult);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len2 = 0;
    // INCHI✔️❌:             szMult[0] = '\0';
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (len1 + len2 < nLenLinearCT)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memcpy(szLinearCT, szElement, len1);
    // INCHI✔️❌:             memcpy(szLinearCT + len1, szMult, (long long)len2 + 1); /*  adding zero termination */ /* djb-rwth: added cast operator */
    // INCHI✔️❌:             return len1 + len2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ( *bOverflow )++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AddElementAndCount

    if multiplicity <= 0 || *overflow != 0 {
        return Ok(0);
    }
    let element_length = element
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if element_length == 0 {
        return Ok(0);
    }
    let multiplicity_text = if multiplicity > 1 {
        multiplicity.to_string()
    } else {
        String::new()
    };
    let visible_length = element_length
        .checked_add(multiplicity_text.len())
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let visible_length_i32 =
        i32::try_from(visible_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if visible_length_i32 < linear_ct_length {
        let required = visible_length
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        let destination = linear_ct
            .get_mut(..required)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        destination[..element_length].copy_from_slice(&element[..element_length]);
        for (target, source) in destination[element_length..visible_length]
            .iter_mut()
            .zip(multiplicity_text.bytes())
        {
            *target = source as i8;
        }
        destination[visible_length] = 0;
        return Ok(visible_length_i32);
    }
    *overflow = overflow.wrapping_add(1);
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeHillFormula(
    atom: &[U_CHAR],
    num_atoms: i32,
    linear_ct: &mut [i8],
    linear_ct_length: i32,
    num_carbons: i32,
    mut num_hydrogens: i32,
    overflow: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:316 MakeHillFormula
    // INCHI✔️❌: int MakeHillFormula( U_CHAR *nAtom,
    // INCHI✔️❌:                      int num_atoms,
    // INCHI✔️❌:                      char *szLinearCT,
    // INCHI✔️❌:                      int nLen_szLinearCT,
    // INCHI✔️❌:                      int num_C,
    // INCHI✔️❌:                      int num_H,
    // INCHI✔️❌:                      int *bOverflow )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szElement[4];
    // INCHI✔️❌:     int  mult, compare2H;
    // INCHI✔️❌:     int  i, nLen, bOvfl;
    // INCHI✔️❌:     U_CHAR nPrevAtom;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLen = 0;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bOvfl = 0;
    // INCHI✔️❌:     nPrevAtom = (U_CHAR) -2; /*  non-existent number */
    // INCHI✔️❌:     memset(szElement, '\0', sizeof(szElement)); /* djb-rwth: fixing coverity ID #499542 */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_C)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nLen += AddElementAndCount( "C", num_C, szLinearCT + nLen, nLen_szLinearCT - nLen, &bOvfl );
    // INCHI✔️❌:         if (num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nLen += AddElementAndCount( "H", num_H, szLinearCT + nLen, nLen_szLinearCT - nLen, &bOvfl );
    // INCHI✔️❌:             num_H = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nPrevAtom != nAtom[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (mult)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nLen += AddElementAndCount( szElement, mult, szLinearCT + nLen, nLen_szLinearCT - nLen, &bOvfl );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             mult = 1;
    // INCHI✔️❌:             /*if (-1 == get_element_chemical_symbol( (int) nAtom[i], szElement ))*/
    // INCHI✔️❌:             if (-1 == get_element_or_pseudoelement_symbol( (int) nAtom[i], szElement ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  wrong element */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             nPrevAtom = nAtom[i];
    // INCHI✔️❌:             if (!strcmp( "C", szElement ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             compare2H = strcmp( "H", szElement );
    // INCHI✔️❌:             if (!compare2H)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (compare2H < 0 && num_H)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*   H-atom should be located in front of szElement */
    // INCHI✔️❌:                 nLen += AddElementAndCount( "H", num_H, szLinearCT + nLen, nLen_szLinearCT - nLen, &bOvfl );
    // INCHI✔️❌:                 num_H = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (mult)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  the last element if any */
    // INCHI✔️❌:         nLen += AddElementAndCount( szElement, mult, szLinearCT + nLen, nLen_szLinearCT - nLen, &bOvfl );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (num_H)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  if H has not been output... */
    // INCHI✔️❌:         nLen += AddElementAndCount( "H", num_H, szLinearCT + nLen, nLen_szLinearCT - nLen, &bOvfl );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *bOverflow |= ( 0 != bOvfl );
    // INCHI✔️❌:
    // INCHI✔️❌:     return bOvfl ? nLen_szLinearCT + 1 : nLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MakeHillFormula

    let mut length = 0_i32;
    let mut multiplicity = 0_i32;
    let mut local_overflow = 0_i32;
    let mut previous_atom = (-2_i32) as U_CHAR;
    let mut symbol = [0_i8; 4];

    if num_carbons != 0 {
        let offset = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let destination = linear_ct
            .get_mut(offset..)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        length = length.wrapping_add(AddElementAndCount(
            &[b'C' as i8, 0],
            num_carbons,
            destination,
            linear_ct_length.wrapping_sub(length),
            &mut local_overflow,
        )?);
        if num_hydrogens != 0 {
            let offset =
                usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let destination = linear_ct
                .get_mut(offset..)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            length = length.wrapping_add(AddElementAndCount(
                &[b'H' as i8, 0],
                num_hydrogens,
                destination,
                linear_ct_length.wrapping_sub(length),
                &mut local_overflow,
            )?);
            num_hydrogens = 0;
        }
    }

    let atom_count = if num_atoms > 0 {
        usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    if atom.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    for &current_atom in &atom[..atom_count] {
        if previous_atom != current_atom {
            if multiplicity != 0 {
                let offset =
                    usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let destination = linear_ct
                    .get_mut(offset..)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                length = length.wrapping_add(AddElementAndCount(
                    &symbol,
                    multiplicity,
                    destination,
                    linear_ct_length.wrapping_sub(length),
                    &mut local_overflow,
                )?);
            }
            multiplicity = 1;
            if get_element_or_pseudoelement_symbol(i32::from(current_atom), &mut symbol)? == -1 {
                return Ok(-1);
            }
            previous_atom = current_atom;
            let symbol_length = symbol.iter().position(|byte| *byte == 0).unwrap_or(4);
            let symbol_bytes = symbol.map(|byte| byte as u8);
            let symbol_text = &symbol_bytes[..symbol_length];
            if symbol_text == b"C" || symbol_text == b"H" {
                return Ok(-1);
            }
            if b"H".as_slice() < symbol_text && num_hydrogens != 0 {
                let offset =
                    usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let destination = linear_ct
                    .get_mut(offset..)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                length = length.wrapping_add(AddElementAndCount(
                    &[b'H' as i8, 0],
                    num_hydrogens,
                    destination,
                    linear_ct_length.wrapping_sub(length),
                    &mut local_overflow,
                )?);
                num_hydrogens = 0;
            }
        } else {
            multiplicity = multiplicity.wrapping_add(1);
        }
    }
    if multiplicity != 0 {
        let offset = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let destination = linear_ct
            .get_mut(offset..)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        length = length.wrapping_add(AddElementAndCount(
            &symbol,
            multiplicity,
            destination,
            linear_ct_length.wrapping_sub(length),
            &mut local_overflow,
        )?);
    }
    if num_hydrogens != 0 {
        let offset = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let destination = linear_ct
            .get_mut(offset..)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        length = length.wrapping_add(AddElementAndCount(
            &[b'H' as i8, 0],
            num_hydrogens,
            destination,
            linear_ct_length.wrapping_sub(length),
            &mut local_overflow,
        )?);
    }
    *overflow |= i32::from(local_overflow != 0);
    Ok(if local_overflow != 0 {
        linear_ct_length.wrapping_add(1)
    } else {
        length
    })
}

#[allow(non_snake_case)]
pub(crate) fn AllocateAndFillHillFormula(
    heap: &mut SourceHeap,
    inchi: &INChI,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // Active libinchi configuration: FIX_GAF_2019_2 == 1 selects inchi_calloc.
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:402 AllocateAndFillHillFormula
    // INCHI✔️❌: char *AllocateAndFillHillFormula( INChI *pINChI )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int num_C, num_H, nLen, nNumNonHAtoms, ret, bOverflow;
    // INCHI✔️❌:     char *pHillFormula = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     bOverflow = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!GetHillFormulaCounts( pINChI->nAtom,
    // INCHI✔️❌:         pINChI->nNum_H,
    // INCHI✔️❌:         pINChI->nNumberOfAtoms,
    // INCHI✔️❌:         pINChI->nTautomer,
    // INCHI✔️❌:         pINChI->lenTautomer,
    // INCHI✔️❌:         &num_C, &num_H,
    // INCHI✔️❌:         &nLen,
    // INCHI✔️❌:         &nNumNonHAtoms ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( FIX_GAF_2019_2==1 )
    // INCHI✔️❌:         pHillFormula = (char*)inchi_calloc((long long )nLen + 1, sizeof(char)); /* djb-rwth: cast operator added */
    // INCHI✔️❌: #else
    // INCHI✔️❌:         pHillFormula = (char*)inchi_malloc(nLen + 1);
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         if (pHillFormula)
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌:             ret = MakeHillFormula( pINChI->nAtom + num_C,
    // INCHI✔️❌:                                    nNumNonHAtoms - num_C,
    // INCHI✔️❌:                                    pHillFormula,
    // INCHI✔️❌:                                    nLen + 1,
    // INCHI✔️❌:                                    num_C, num_H,
    // INCHI✔️❌:                                    &bOverflow );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (ret != nLen || bOverflow)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( pHillFormula );
    // INCHI✔️❌:                 pHillFormula = NULL;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return pHillFormula;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AllocateAndFillHillFormula

    let mut num_carbons = 0_i32;
    let mut num_hydrogens = 0_i32;
    let mut formula_length = 0_i32;
    let mut num_non_hydrogen_atoms = 0_i32;
    let count_result = {
        let atom_count = if inchi.nNumberOfAtoms > 0 {
            usize::try_from(inchi.nNumberOfAtoms)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        } else {
            0
        };
        let atoms = if atom_count == 0 && inchi.nAtom.is_null() {
            &[][..]
        } else {
            heap.slice(inchi.nAtom.as_const())?
        };
        let atom_hydrogens = if atom_count == 0 && inchi.nNum_H.is_null() {
            &[][..]
        } else {
            heap.slice(inchi.nNum_H.as_const())?
        };
        let tautomer = if inchi.nTautomer.is_null() {
            None
        } else {
            Some(heap.slice(inchi.nTautomer.as_const())?)
        };
        GetHillFormulaCounts(
            atoms,
            atom_hydrogens,
            inchi.nNumberOfAtoms,
            tautomer,
            inchi.lenTautomer,
            &mut num_carbons,
            &mut num_hydrogens,
            &mut formula_length,
            &mut num_non_hydrogen_atoms,
        )?
    };
    if count_result != 0 {
        return Ok(SourceMutPointer::null());
    }

    let allocation_count = i64::from(formula_length) + 1;
    let allocation_count =
        u64::try_from(allocation_count).map_err(|_| SourceHeapError::AllocationSizeOverflow)?;
    let hill_formula = match inchi_calloc::<i8>(heap, allocation_count, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };

    let make_result = heap.with_slice_mut_and_heap(hill_formula, |output, heap| {
        let atom_offset =
            usize::try_from(num_carbons).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atoms = if inchi.nAtom.is_null() {
            if atom_offset == 0 {
                &[][..]
            } else {
                return Err(SourceHeapError::NullPointer);
            }
        } else {
            heap.slice(inchi.nAtom.as_const())?
                .get(atom_offset..)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let mut overflow = 0_i32;
        let result = MakeHillFormula(
            atoms,
            num_non_hydrogen_atoms.wrapping_sub(num_carbons),
            output,
            formula_length.wrapping_add(1),
            num_carbons,
            num_hydrogens,
            &mut overflow,
        )?;
        Ok((result, overflow))
    });
    let (result, overflow) = match make_result {
        Ok(values) => values,
        Err(error) => {
            inchi_free(heap, hill_formula)?;
            return Err(error);
        }
    };
    if result != formula_length || overflow != 0 {
        inchi_free(heap, hill_formula)?;
        return Ok(SourceMutPointer::null());
    }
    Ok(hill_formula)
}

fn c_fixed(value: f64, precision: usize) -> String {
    let value = if value.is_nan() {
        if value.is_sign_negative() {
            "-nan".to_owned()
        } else {
            "nan".to_owned()
        }
    } else if value == f64::INFINITY {
        "inf".to_owned()
    } else if value == f64::NEG_INFINITY {
        "-inf".to_owned()
    } else {
        format!("{value:.precision$}")
    };
    format!("{value:>10}")
}

fn c_scientific(value: f64, precision: usize) -> String {
    let value = if value.is_nan() {
        if value.is_sign_negative() {
            "-nan".to_owned()
        } else {
            "nan".to_owned()
        }
    } else if value == f64::INFINITY {
        "inf".to_owned()
    } else if value == f64::NEG_INFINITY {
        "-inf".to_owned()
    } else {
        let rust = format!("{value:.precision$e}");
        let (mantissa, exponent) = rust
            .split_once('e')
            .expect("Rust scientific formatting always contains an exponent");
        let exponent = exponent
            .parse::<i32>()
            .expect("Rust scientific formatting emits a decimal exponent");
        format!("{mantissa}e{exponent:+03}")
    };
    format!("{value:>10}")
}

pub(crate) fn WriteCoord(output: &mut [i8], value: f64) -> Result<(), SourceHeapError> {
    // Active libinchi configuration: TARGET_API_LIB includes this function;
    // GHI100_FIX is undefined, so every active leaf is the trailing sprintf.
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:890 WriteCoord
    // INCHI✔❌: void WriteCoord( char *str, double x )
    // INCHI✔❌: {
    // INCHI✔❌:     if (x < -9999999.9)
    // INCHI✔❌:     {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:         dbl2int(str, 10, 2, 'e', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:         stbsp_sprintf(str, "%10.2e", x);
    // INCHI✔❌: #else
    // INCHI✔❌:         sprintf(str, "%10.2e", x);
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:         sprintf(str, "%10.2e", x);
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         if (x < -999999.99)
    // INCHI✔❌:         {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:             dbl2int(str, 10, 2, 'f', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:             stbsp_sprintf(str, "%10.2f", x);
    // INCHI✔❌: #else
    // INCHI✔❌:             sprintf( str, "%10.2f", x );
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:             sprintf( str, "%10.2f", x );
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             if (x < -99999.999)
    // INCHI✔❌:             {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:                 dbl2int(str, 10, 3, 'f', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:                 stbsp_sprintf(str, "%10.3f", x);
    // INCHI✔❌: #else
    // INCHI✔❌:                 sprintf( str, "%10.3f", x );
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:                 sprintf( str, "%10.3f", x );
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 if (x < 99999.9999)
    // INCHI✔❌:                 {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:                     dbl2int(str, 10, 4, 'f', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:                     stbsp_sprintf(str, "%10.4f", x);
    // INCHI✔❌: #else
    // INCHI✔❌:                     sprintf( str, "%10.4f", x );
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:                     sprintf( str, "%10.4f", x );
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (x < 999999.999)
    // INCHI✔❌:                     {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:                         dbl2int(str, 10, 3, 'f', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:                         stbsp_sprintf(str, "%10.3f", x);
    // INCHI✔❌: #else
    // INCHI✔❌:                         sprintf( str, "%10.3f", x );
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:                         sprintf( str, "%10.3f", x );
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (x < 9999999.99)
    // INCHI✔❌:                         {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:                             dbl2int(str, 10, 2, 'f', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:                             stbsp_sprintf(str, "%10.2f", x);
    // INCHI✔❌: #else
    // INCHI✔❌:                             sprintf( str, "%10.2f", x );
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:                             sprintf( str, "%10.2f", x );
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (x < 99999999.9)
    // INCHI✔❌:                             {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:                                 dbl2int(str, 10, 1, 'f', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:                                 stbsp_sprintf(str, "%10.1f", x);
    // INCHI✔❌: #else
    // INCHI✔❌:                                 sprintf( str, "%10.1f", x );
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:                                 sprintf( str, "%10.1f", x );
    // INCHI✔❌:                             }
    // INCHI✔❌:                             else
    // INCHI✔❌:                             {
    // INCHI✔❌: #ifdef GHI100_FIX
    // INCHI✔❌: #if (SPRINTF_FLAG == 2)
    // INCHI✔❌:                                 dbl2int(str, 10, 3, 'e', x);
    // INCHI✔❌: #elif (SPRINTF_FLAG == 1)
    // INCHI✔❌:                                 stbsp_sprintf(str, "%10.3e", x);
    // INCHI✔❌: #else
    // INCHI✔❌:                                 sprintf( str, "%10.3e", x );
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:                                 sprintf( str, "%10.3e", x );
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌: }
    // END INCHI C FUNCTION: WriteCoord

    let formatted = if value < -9_999_999.9 {
        c_scientific(value, 2)
    } else if value < -999_999.99 {
        c_fixed(value, 2)
    } else if value < -99_999.999 {
        c_fixed(value, 3)
    } else if value < 99_999.9999 {
        c_fixed(value, 4)
    } else if value < 999_999.999 {
        c_fixed(value, 3)
    } else if value < 9_999_999.99 {
        c_fixed(value, 2)
    } else if value < 99_999_999.9 {
        c_fixed(value, 1)
    } else {
        c_scientific(value, 3)
    };
    let required = formatted
        .len()
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    if output.len() < required {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    for (destination, source) in output.iter_mut().zip(formatted.bytes()) {
        *destination = source as i8;
    }
    output[formatted.len()] = 0;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn format_coordinate(value: f64) -> String {
        let mut output = [b'X' as i8; 32];
        WriteCoord(&mut output, value).unwrap();
        let length = output.iter().position(|byte| *byte == 0).unwrap();
        let text =
            String::from_utf8(output[..length].iter().map(|byte| *byte as u8).collect()).unwrap();
        assert_eq!(output[length + 1], b'X' as i8);
        text
    }

    #[test]
    fn source_port__ichimak2__writecoord__line_890() {
        for (value, expected) in [
            (-10_000_000.0, " -1.00e+07"),
            (-9_999_999.9, "-9999999.90"),
            (-999_999.99, "-999999.990"),
            (-99_999.999, "-99999.9990"),
            (-0.0, "   -0.0000"),
            (0.0, "    0.0000"),
            (99_999.9999, "100000.000"),
            (999_999.999, "1000000.00"),
            (9_999_999.99, "10000000.0"),
            (99_999_999.9, " 1.000e+08"),
            (f64::INFINITY, "       inf"),
            (f64::NEG_INFINITY, "      -inf"),
        ] {
            assert_eq!(format_coordinate(value), expected, "{value:?}");
        }
        assert_eq!(format_coordinate(f64::NAN), "       nan");

        // These exact strings were emitted by the independent official C
        // WriteCoord oracle under LC_NUMERIC=C.
        for (bits, expected) in [
            (13935002356523257038, " -1.00e+07"),
            (13935002356523257037, "-9999999.90"),
            (13935002356523257036, "-9999999.90"),
            (13920209183453562799, "-999999.99"),
            (13920209183453562798, "-999999.990"),
            (13920209183453562797, "-999999.990"),
            (13904980397670231180, "-99999.999"),
            (13904980397670231179, "-99999.9990"),
            (13904980397670231178, "-99999.9990"),
            (4681608360877302899, "99999.9999"),
            (4681608360877302900, "100000.000"),
            (4681608360877302901, "100000.000"),
            (4696837146676096400, "999999.999"),
            (4696837146676096401, "1000000.00"),
            (4696837146676096402, "1000000.00"),
            (4711630319716799610, "9999999.99"),
            (4711630319716799611, "10000000.0"),
            (4711630319716799612, "10000000.0"),
            (4726483295877568921, "99999999.9"),
            (4726483295877568922, " 1.000e+08"),
            (4726483295877568923, " 1.000e+08"),
            (0, "    0.0000"),
            (9223372036854775808, "   -0.0000"),
            (1, "    0.0000"),
            (9223372036854775809, "   -0.0000"),
            (4503599627370495, "    0.0000"),
            (9227875636482146303, "   -0.0000"),
            (4503599627370496, "    0.0000"),
            (9227875636482146304, "   -0.0000"),
            (9218868437227405311, "1.798e+308"),
            (18442240474082181119, "-1.80e+308"),
            (9218868437227405312, "       inf"),
            (18442240474082181120, "      -inf"),
            (9221120237041090560, "       nan"),
            (18444492273895866368, "      -nan"),
            (9218868437227405313, "       nan"),
            (18442240474082181121, "      -nan"),
            (4608238287732654420, "    1.2344"),
            (4608238287732654421, "    1.2345"),
            (4608238287732654422, "    1.2345"),
            (4608238738092617157, "    1.2345"),
            (4608238738092617158, "    1.2346"),
            (4608238738092617159, "    1.2346"),
            (13831610324587430230, "   -1.2345"),
            (13831610324587430229, "   -1.2345"),
            (13831610324587430228, "   -1.2344"),
            (13831610774947392967, "   -1.2346"),
            (13831610774947392966, "   -1.2346"),
            (13831610774947392965, "   -1.2345"),
            (13904980397773310395, "-100000.001"),
            (13904980397773310394, "-100000.000"),
            (13904980397773310393, "-100000.000"),
            (4681608360918534585, "100000.000"),
            (4681608360918534586, "100000.000"),
            (4681608360918534587, "100000.001"),
            (13920209183582411818, "-1000000.01"),
            (13920209183582411817, "-1000000.01"),
            (13920209183582411816, "-1000000.00"),
            (4696837146727636008, "1000000.00"),
            (4696837146727636009, "1000000.01"),
            (4696837146727636010, "1000000.01"),
            (4711630319749011865, "10000000.0"),
            (4711630319749011866, "10000000.1"),
            (4711630319749011867, "10000000.1"),
            (13935029200122544129, " -1.01e+07"),
            (13935029200122544128, " -1.00e+07"),
            (13935029200122544127, " -1.00e+07"),
            (13935082887213744129, " -1.02e+07"),
            (13935082887213744128, " -1.02e+07"),
            (13935082887213744127, " -1.01e+07"),
            (4726486651327479807, " 1.000e+08"),
            (4726486651327479808, " 1.000e+08"),
            (4726486651327479809, " 1.001e+08"),
            (4726493362213879807, " 1.001e+08"),
            (4726493362213879808, " 1.002e+08"),
            (4726493362213879809, " 1.002e+08"),
        ] {
            assert_eq!(
                format_coordinate(f64::from_bits(bits)),
                expected,
                "{bits:#018x}"
            );
        }

        let mut short = [0_i8; 4];
        assert_eq!(
            WriteCoord(&mut short, 0.0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn official_c_oracle__ichimak2__writecoord__line_890() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--write-coord-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );

        let mut record_count = 0;
        for line in String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
        {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["operation"], "write_coord");
            assert_eq!(official["input"]["locale"], "C");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let bits = official["input"]["bits"]
                .as_u64()
                .expect("input bits must be u64");
            let official_bytes = official["output"]["bytes"]
                .as_array()
                .expect("output bytes must be an array")
                .iter()
                .map(|value| value.as_u64().expect("output byte must be u8") as u8)
                .collect::<Vec<_>>();
            assert_eq!(official_bytes.len(), 64);
            assert!(
                official_bytes.contains(&0),
                "{case_id} has no NUL terminator"
            );

            let mut rust_output = [0xa5_u8 as i8; 64];
            WriteCoord(&mut rust_output, f64::from_bits(bits)).unwrap();
            let rust_bytes = rust_output.map(|byte| byte as u8);
            assert_eq!(official_bytes, rust_bytes, "byte mismatch for {case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 4172);
    }

    fn formula_text(heap: &mut SourceHeap, text: &str) -> SourceConstPointer<i8> {
        heap.allocate_model_storage(
            text.bytes()
                .map(|byte| byte as i8)
                .chain(std::iter::once(0))
                .collect(),
        )
        .unwrap()
        .as_const()
    }

    fn formula_buffer(heap: &mut SourceHeap, text: &str, allocated: usize) -> INCHI_IOS_STRING {
        let mut bytes = vec![71_i8; allocated];
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
    fn source_port__ichimak2__makehillformulastring__line_121() {
        let mut heap = SourceHeap::default();
        let mut buffer = formula_buffer(&mut heap, "pre", 32);
        let mut overflow = 0;
        let formula = formula_text(&mut heap, "C2H6");
        assert_eq!(
            MakeHillFormulaString(&mut heap, formula, &mut buffer, &mut overflow),
            Ok(4)
        );
        assert_eq!(overflow, 0);
        assert_eq!(buffer.nUsedLength, 7);
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..8],
            &[
                b'p' as i8, b'r' as i8, b'e' as i8, b'C' as i8, b'2' as i8, b'H' as i8, b'6' as i8,
                0
            ]
        );

        let before = buffer.clone();
        assert_eq!(
            MakeHillFormulaString(
                &mut heap,
                SourceConstPointer::null(),
                &mut buffer,
                &mut overflow,
            ),
            Ok(0)
        );
        assert_eq!(buffer, before);

        overflow = 2;
        assert_eq!(
            MakeHillFormulaString(&mut heap, formula, &mut buffer, &mut overflow),
            Ok(0)
        );
        assert_eq!(overflow, 2);
        assert_eq!(buffer, before);

        let mut failing_heap = SourceHeap::default();
        let formula = formula_text(&mut failing_heap, "C2H6");
        let mut failing = formula_buffer(&mut failing_heap, "abc", 4);
        let pointer = failing.pStr;
        let mut overflow = 0;
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            MakeHillFormulaString(&mut failing_heap, formula, &mut failing, &mut overflow,),
            Ok(4)
        );
        assert_eq!(overflow, 1);
        assert_eq!(failing.pStr, pointer);
        assert_eq!(failing.nUsedLength, 3);
    }

    #[test]
    fn source_port__ichimak2__gethillformulaindexlength__line_144() {
        for count in [i32::MIN, -100, -1, 0, 1] {
            assert_eq!(GetHillFormulaIndexLength(count), 0, "{count}");
        }
        for (count, expected) in [
            (2, 1),
            (9, 1),
            (10, 2),
            (99, 2),
            (100, 3),
            (999, 3),
            (1_000, 4),
            (9_999, 4),
            (10_000, 5),
            (99_999, 5),
            (100_000, 6),
            (999_999, 6),
            (1_000_000, 7),
            (9_999_999, 7),
            (10_000_000, 8),
            (99_999_999, 8),
            (100_000_000, 9),
            (999_999_999, 9),
            (1_000_000_000, 10),
            (i32::MAX, 10),
        ] {
            assert_eq!(GetHillFormulaIndexLength(count), expected, "{count}");
        }
    }

    #[test]
    fn source_port__ichimak2__gethillformulacounts__line_159() {
        fn counts(
            atom: &[U_CHAR],
            atom_h: &[S_CHAR],
            num_atoms: i32,
            tautomer: Option<&[AT_NUMB]>,
            tautomer_length: i32,
        ) -> Result<(i32, [i32; 4]), SourceHeapError> {
            let mut output = [91, 92, 93, 94];
            let (carbon, rest) = output.split_at_mut(1);
            let (hydrogen, rest) = rest.split_at_mut(1);
            let (length, non_hydrogen) = rest.split_at_mut(1);
            let result = GetHillFormulaCounts(
                atom,
                atom_h,
                num_atoms,
                tautomer,
                tautomer_length,
                &mut carbon[0],
                &mut hydrogen[0],
                &mut length[0],
                &mut non_hydrogen[0],
            )?;
            Ok((result, output))
        }

        assert_eq!(counts(&[], &[], 0, None, 0), Ok((0, [0, 0, 0, 0])));
        assert_eq!(counts(&[], &[], -7, None, 0), Ok((0, [0, 0, 0, -7])));

        let atoms = [6, 6, 8, 8, 8, 1, 1];
        let attached_h = [3, 2, 0, 0, 0, 0, 0];
        assert_eq!(
            counts(&atoms, &attached_h, 7, None, 0),
            Ok((0, [2, 7, 6, 5]))
        );

        let tautomer = [2, 3, 2, 10, 11, 2, 1, 12];
        assert_eq!(
            counts(&atoms, &attached_h, 7, Some(&tautomer), 8),
            Ok((0, [2, 10, 7, 5]))
        );
        assert_eq!(
            counts(&atoms, &attached_h, 7, None, 8),
            Ok((0, [2, 7, 6, 5]))
        );
        assert_eq!(
            counts(&atoms, &attached_h, 7, Some(&tautomer), 0),
            Ok((0, [2, 7, 6, 5]))
        );

        assert_eq!(
            counts(&[17, 17, 17], &[-1, 0, 0], 3, None, 0),
            Ok((0, [0, -1, 4, 3]))
        );
        assert_eq!(counts(&[120], &[0], 1, None, 0), Ok((0, [0, 0, 2, 1])));

        assert_eq!(counts(&[0], &[0], 1, None, 0), Ok((-1, [91, 92, 93, 94])));
        assert_eq!(
            counts(&[6], &[], 1, None, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            counts(&atoms, &attached_h, 7, Some(&tautomer[..3]), 8),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            counts(&[6], &[0], 1, Some(&[1, 1]), 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        // The source initializes nPrevAtom to `(U_CHAR)-2`, so an otherwise
        // invalid first atomic number 254 bypasses symbol lookup.
        assert_eq!(counts(&[254], &[0], 1, None, 0), Ok((0, [0, 0, 0, 1])));
    }

    #[test]
    fn source_port__ichimak2__addelementandcount__line_279() {
        let carbon = [b'C' as i8, 0];
        let chlorine = [b'C' as i8, b'l' as i8, 0];
        let empty = [0_i8];

        let mut output = [b'X' as i8; 8];
        let mut overflow = 0;
        assert_eq!(
            AddElementAndCount(&carbon, 1, &mut output, 8, &mut overflow),
            Ok(1)
        );
        assert_eq!(output, [b'C' as i8, 0, 88, 88, 88, 88, 88, 88]);
        assert_eq!(overflow, 0);

        output.fill(b'X' as i8);
        assert_eq!(
            AddElementAndCount(&chlorine, 12, &mut output, 8, &mut overflow),
            Ok(4)
        );
        assert_eq!(
            output,
            [
                b'C' as i8, b'l' as i8, b'1' as i8, b'2' as i8, 0, 88, 88, 88
            ]
        );

        for multiplicity in [i32::MIN, -1, 0] {
            let before = output;
            assert_eq!(
                AddElementAndCount(&[b'Q' as i8], multiplicity, &mut output, 8, &mut overflow),
                Ok(0)
            );
            assert_eq!(output, before);
        }
        let before = output;
        assert_eq!(
            AddElementAndCount(&empty, 2, &mut output, 8, &mut overflow),
            Ok(0)
        );
        assert_eq!(output, before);

        overflow = 9;
        assert_eq!(
            AddElementAndCount(&[b'Q' as i8], 2, &mut output, 8, &mut overflow),
            Ok(0)
        );
        assert_eq!(overflow, 9);
        assert_eq!(output, before);

        overflow = 0;
        output.fill(b'X' as i8);
        assert_eq!(
            AddElementAndCount(&chlorine, 12, &mut output, 4, &mut overflow),
            Ok(0)
        );
        assert_eq!(overflow, 1);
        assert_eq!(output, [b'X' as i8; 8]);

        overflow = i32::MAX;
        assert_eq!(
            AddElementAndCount(&carbon, 1, &mut output, 1, &mut overflow),
            Ok(0)
        );
        assert_eq!(overflow, i32::MAX);

        overflow = 0;
        assert_eq!(
            AddElementAndCount(&carbon, 1, &mut [], 2, &mut overflow),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(overflow, 0);
        assert_eq!(
            AddElementAndCount(&[b'C' as i8], 1, &mut output, 8, &mut overflow),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichimak2__makehillformula__line_316() {
        fn run(
            atoms: &[U_CHAR],
            num_atoms: i32,
            capacity: i32,
            num_carbons: i32,
            num_hydrogens: i32,
            initial_overflow: i32,
        ) -> (Result<i32, SourceHeapError>, i32, [i8; 16]) {
            let mut output = [b'X' as i8; 16];
            let mut overflow = initial_overflow;
            let result = MakeHillFormula(
                atoms,
                num_atoms,
                &mut output,
                capacity,
                num_carbons,
                num_hydrogens,
                &mut overflow,
            );
            (result, overflow, output)
        }

        let (result, overflow, output) = run(&[8, 8], 2, 16, 2, 6, 0);
        assert_eq!(result, Ok(6));
        assert_eq!(overflow, 0);
        assert_eq!(
            &output[..8],
            &[
                b'C' as i8, b'2' as i8, b'H' as i8, b'6' as i8, b'O' as i8, b'2' as i8, 0, 88
            ]
        );

        let (result, overflow, output) = run(&[9, 9, 8], 3, 16, 0, 3, 0);
        assert_eq!(result, Ok(5));
        assert_eq!(overflow, 0);
        assert_eq!(
            &output[..7],
            &[
                b'F' as i8, b'2' as i8, b'H' as i8, b'3' as i8, b'O' as i8, 0, 88
            ]
        );

        let (result, _, output) = run(&[9, 9], 2, 16, 0, 2, 0);
        assert_eq!(result, Ok(4));
        assert_eq!(
            &output[..6],
            &[b'F' as i8, b'2' as i8, b'H' as i8, b'2' as i8, 0, 88]
        );
        let (result, _, output) = run(&[], -1, 16, 0, 2, 0);
        assert_eq!(result, Ok(2));
        assert_eq!(&output[..4], &[b'H' as i8, b'2' as i8, 0, 88]);
        let (result, _, output) = run(&[], 0, 16, -1, 2, 0);
        assert_eq!(result, Ok(2));
        assert_eq!(&output[..4], &[b'H' as i8, b'2' as i8, 0, 88]);

        for invalid in [0, 1, 6] {
            let (result, overflow, _) = run(&[invalid], 1, 16, 0, 0, 7);
            assert_eq!(result, Ok(-1), "{invalid}");
            assert_eq!(overflow, 7, "{invalid}");
        }
        let (result, overflow, output) = run(&[6], 1, 16, 2, 0, 0);
        assert_eq!(result, Ok(-1));
        assert_eq!(overflow, 0);
        assert_eq!(&output[..4], &[b'C' as i8, b'2' as i8, 0, 88]);

        let (result, overflow, output) = run(&[8], 1, 3, 2, 0, 0);
        assert_eq!(result, Ok(4));
        assert_eq!(overflow, 1);
        assert_eq!(&output[..4], &[b'C' as i8, b'2' as i8, 0, 88]);
        let (result, overflow, output) = run(&[8], 1, 3, 2, 0, 2);
        assert_eq!(result, Ok(4));
        assert_eq!(overflow, 3);
        assert_eq!(&output[..4], &[b'C' as i8, b'2' as i8, 0, 88]);

        let (result, overflow, output) = run(&[0], 1, 2, 2, 0, 4);
        assert_eq!(result, Ok(-1));
        assert_eq!(overflow, 4);
        assert_eq!(output, [b'X' as i8; 16]);

        let (result, overflow, output) = run(&[254], 1, 16, 0, 0, 0);
        assert_eq!(result, Ok(0));
        assert_eq!(overflow, 0);
        assert_eq!(output, [b'X' as i8; 16]);
        assert_eq!(
            MakeHillFormula(&[8], 2, &mut [0_i8; 8], 8, 0, 0, &mut 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimak2__allocateandfillhillformula__line_402() {
        fn fixture(atoms: &[U_CHAR], hydrogens: &[S_CHAR]) -> (SourceHeap, INChI) {
            let mut heap = SourceHeap::default();
            let atom_pointer = if atoms.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(atoms.to_vec()).unwrap()
            };
            let hydrogen_pointer = if hydrogens.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(hydrogens.to_vec()).unwrap()
            };
            (
                heap,
                INChI {
                    nNumberOfAtoms: atoms.len() as i32,
                    nAtom: atom_pointer,
                    nNum_H: hydrogen_pointer,
                    ..INChI::default()
                },
            )
        }

        let (mut heap, inchi) = fixture(&[6, 6, 8, 8], &[3, 2, 0, 0]);
        let formula = AllocateAndFillHillFormula(&mut heap, &inchi).unwrap();
        assert!(!formula.is_null());
        assert_eq!(
            &heap.slice(formula.as_const()).unwrap()[..7],
            &[
                b'C' as i8, b'2' as i8, b'H' as i8, b'5' as i8, b'O' as i8, b'2' as i8, 0
            ]
        );

        let (mut heap, mut inchi) = fixture(&[6, 8], &[0, 0]);
        inchi.nTautomer = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        inchi.lenTautomer = 3;
        let formula = AllocateAndFillHillFormula(&mut heap, &inchi).unwrap();
        assert_eq!(
            &heap.slice(formula.as_const()).unwrap()[..5],
            &[b'C' as i8, b'H' as i8, b'3' as i8, b'O' as i8, 0]
        );

        let (mut heap, inchi) = fixture(&[], &[]);
        let formula = AllocateAndFillHillFormula(&mut heap, &inchi).unwrap();
        assert!(!formula.is_null());
        assert_eq!(heap.slice(formula.as_const()).unwrap(), &[0]);

        let (mut heap, inchi) = fixture(&[0], &[0]);
        heap.trace_source_allocations();
        let formula = AllocateAndFillHillFormula(&mut heap, &inchi).unwrap();
        assert!(formula.is_null());
        assert_eq!(heap.source_allocation_calls(), 0);

        let (mut heap, inchi) = fixture(&[6], &[4]);
        heap.trace_source_allocations();
        heap.fail_after_allocations(0);
        let formula = AllocateAndFillHillFormula(&mut heap, &inchi).unwrap();
        assert!(formula.is_null());
        assert_eq!(heap.source_allocation_calls(), 1);

        // This source-invalid ordering makes the counting pass, then makes
        // MakeHillFormula reject C after the source's `nAtom + num_C` offset.
        let (mut heap, inchi) = fixture(&[8, 6], &[0, 0]);
        heap.trace_source_allocations();
        let formula = AllocateAndFillHillFormula(&mut heap, &inchi).unwrap();
        assert!(formula.is_null());
        assert_eq!(heap.source_allocation_calls(), 1);
    }
}
