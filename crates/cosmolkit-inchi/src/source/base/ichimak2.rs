use crate::source::base::ichi_io::inchi_strbuf_printf;
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichimap1::CompareLinCtStereoDble;
use crate::source::base::ichimap2::switch_ptrs;
use crate::source::base::ichisort::{CompRanksOrd, inchi_qsort};
use crate::source::base::ichitaut::SortTautomerGroupsAndEndpoints;
use crate::source::base::util::{
    get_element_or_pseudoelement_symbol, get_unusual_el_valence, inchi_calloc, inchi_free,
};
use crate::source_types::{
    AB_MAX_KNOWN_PARITY, AB_MAX_WELL_DEFINED_PARITY, AB_MIN_KNOWN_PARITY,
    AB_MIN_WELL_DEFINED_PARITY, AB_PARITY_UNKN, AMBIGUOUS_STEREO_ATOM, AMBIGUOUS_STEREO_ATOM_ISO,
    AMBIGUOUS_STEREO_BOND, AMBIGUOUS_STEREO_BOND_ISO, AT_NUMB, AT_RANK, AT_STEREO_CARB,
    AT_STEREO_DBLE, CANON_GLOBALS, CANON_STAT, CT_OUT_OF_RAM, CT_WRONG_FORMULA,
    ERR_NO_CANON_RESULTS, FLAG_FORCE_SALT_TAUT, FLAG_NORM_CONSIDER_TAUT, FLAG_PROTON_CHARGE_CANCEL,
    INCHI_FLAG_ACID_TAUT, INCHI_FLAG_HARD_ADD_REM_PROTON, INCHI_FLAG_RAC_STEREO,
    INCHI_FLAG_REL_STEREO, INCHI_FLAG_SB_IGN_ALL_UU, INCHI_FLAG_SC_IGN_ALL_ISO_UU,
    INCHI_FLAG_SC_IGN_ALL_UU, INCHI_IOS_STRING, INCHI_MODE, INCHI_T_NUM_MOVABLE, INChI, INChI_Aux,
    INChI_Stereo, MASK_CUMULENE_LEN, MAX_ATOMS, MAX_NUM_STEREO_BONDS, MULT_STEREOBOND,
    NUM_H_ISOTOPES, RADICAL_DOUBLET, RADICAL_SINGLET, RADICAL_TRIPLET, REQ_MODE_RACEMIC_STEREO,
    REQ_MODE_RELATIVE_STEREO, REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU, S_CHAR,
    SourceConstPointer, SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer,
    SourceVaList, TG_FLAG_ALL_SALT_DONE, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE,
    TG_FLAG_FOUND_ISOTOPIC_H_DONE, U_CHAR, WARN_FAILED_ISOTOPIC, WARN_FAILED_ISOTOPIC_STEREO,
    WARN_FAILED_STEREO, inp_ATOM, sp_ATOM,
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

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Copy2StereoBondOrAllene(
    heap: &mut SourceHeap,
    Stereo: &INChI_Stereo,
    nNumberOfStereoCenters: &mut i32,
    mut nNumberOfStereoBonds: Option<&mut i32>,
    LinearCTStereoDble: &AT_STEREO_DBLE,
    pCanonOrd: SourceConstPointer<AT_NUMB>,
    pCanonRank: SourceConstPointer<AT_RANK>,
    at: SourceConstPointer<sp_ATOM>,
    bIsotopic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:451 Copy2StereoBondOrAllene
    // INCHI✔️❌: int Copy2StereoBondOrAllene( INChI_Stereo *Stereo,
    // INCHI✔️❌:                              int *nNumberOfStereoCenters,
    // INCHI✔️❌:                              int *nNumberOfStereoBonds,
    // INCHI✔️❌:                              AT_STEREO_DBLE *LinearCTStereoDble,
    // INCHI✔️❌:                              AT_NUMB *pCanonOrd,
    // INCHI✔️❌:                              AT_RANK *pCanonRank,
    // INCHI✔️❌:                              sp_ATOM *at,
    // INCHI✔️❌:                              int bIsotopic )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int cumulene_len,
    // INCHI✔️❌:         j,
    // INCHI✔️❌:         next_j = 0 /* ordering number of the central allene atom */,
    // INCHI✔️❌:         next_neigh;
    // INCHI✔️❌:     AT_RANK
    // INCHI✔️❌:         at_num;
    // INCHI✔️❌:     int parity;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pCanonOrd && pCanonRank)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = pCanonOrd[(int) LinearCTStereoDble->at_num1 - 1];
    // INCHI✔️❌:
    // INCHI✔️❌:         /* if allene then find the central atom, at[next_j] */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (bIsotopic)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cumulene_len = BOND_CHAIN_LEN( at[j].stereo_bond_parity2[0] );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (cumulene_len % 2 && ( 1 >= MAX_NUM_STEREO_BONDS || !at[j].stereo_bond_neighbor2[1] ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 next_j = at[j].neighbor[(int) at[j].stereo_bond_ord2[0]];
    // INCHI✔️❌:
    // INCHI✔️❌:                 for (cumulene_len = ( cumulene_len - 1 ) / 2; cumulene_len && 2 == at[next_j].valence; cumulene_len--)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     next_neigh = ( j == at[next_j].neighbor[0] );
    // INCHI✔️❌:                     j = next_j;
    // INCHI✔️❌:                     next_j = at[next_j].neighbor[next_neigh];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* next_j is the central atom */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 cumulene_len = -1; /* not an allene */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cumulene_len = BOND_CHAIN_LEN( at[j].stereo_bond_parity[0] );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (cumulene_len % 2 && ( 1 >= MAX_NUM_STEREO_BONDS || !at[j].stereo_bond_neighbor[1] ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 next_j = at[j].neighbor[(int) at[j].stereo_bond_ord[0]];
    // INCHI✔️❌:
    // INCHI✔️❌:                 for (cumulene_len = ( cumulene_len - 1 ) / 2; cumulene_len && 2 == at[next_j].valence; cumulene_len--)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     next_neigh = ( j == at[next_j].neighbor[0] );
    // INCHI✔️❌:                     j = next_j;
    // INCHI✔️❌:                     next_j = at[next_j].neighbor[next_neigh];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 cumulene_len = -1; /* not an allene */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!cumulene_len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* allene has been found; insert new stereocenter and parity */
    // INCHI✔️❌:
    // INCHI✔️❌:             AT_NUMB *nNumber;
    // INCHI✔️❌:             S_CHAR  *t_parity;
    // INCHI✔️❌:
    // INCHI✔️❌:             nNumber = nNumberOfStereoBonds ? Stereo->nNumber : Stereo->nNumberInv;
    // INCHI✔️❌:             t_parity = nNumberOfStereoBonds ? Stereo->t_parity : Stereo->t_parityInv;
    // INCHI✔️❌:
    // INCHI✔️❌:             at_num = pCanonRank[next_j];
    // INCHI✔️❌:             parity = LinearCTStereoDble->parity;
    // INCHI✔️❌:
    // INCHI✔️❌:             /* free room for the new stereocenter */
    // INCHI✔️❌:             for (j = 0; j < *nNumberOfStereoCenters && Stereo->nNumber[j] < at_num; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (j < *nNumberOfStereoCenters)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memmove(nNumber + j + 1, nNumber + j, (*nNumberOfStereoCenters - (long long)j) * sizeof(nNumber[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                 memmove(t_parity + j + 1, t_parity + j, (*nNumberOfStereoCenters - (long long)j) * sizeof(t_parity[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* fill the new stereo center info */
    // INCHI✔️❌:
    // INCHI✔️❌:             nNumber[j] = at_num;
    // INCHI✔️❌:             t_parity[j] = parity;
    // INCHI✔️❌:             ( *nNumberOfStereoCenters )++;
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Save the stereo bond info */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = *nNumberOfStereoBonds;
    // INCHI✔️❌:         Stereo->b_parity[j] = LinearCTStereoDble->parity;
    // INCHI✔️❌:         Stereo->nBondAtom1[j] = LinearCTStereoDble->at_num1;
    // INCHI✔️❌:         Stereo->nBondAtom2[j] = LinearCTStereoDble->at_num2;
    // INCHI✔️❌:         ( *nNumberOfStereoBonds )++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Copy2StereoBondOrAllene
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: Copy2StereoBondOrAllene
    // INCHI✔️❌: #define MASK_CUMULENE_LEN 0x38
    // INCHI✔️❌: #define MULT_STEREOBOND 0x08
    // INCHI✔️❌: #define BOND_CHAIN_LEN(X) (((X)&MASK_CUMULENE_LEN)/MULT_STEREOBOND)
    // INCHI✔️❌: #define MAX_NUM_STEREO_BONDS 3
    // END INCHI ACTIVE MACRO CONFIGURATION: Copy2StereoBondOrAllene

    let use_primary_arrays = nNumberOfStereoBonds.is_some();
    if !pCanonOrd.is_null() && !pCanonRank.is_null() {
        let allene = {
            let canonical_index = usize::from(LinearCTStereoDble.at_num1)
                .checked_sub(1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut j = usize::from(
                *heap
                    .slice(pCanonOrd)?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let atoms = heap.slice(at)?;
            let first = atoms.get(j).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let (encoded_parity, stereo_neighbor, stereo_order) = if bIsotopic != 0 {
                (
                    first.stereo_bond_parity2[0],
                    first.stereo_bond_neighbor2,
                    first.stereo_bond_ord2[0],
                )
            } else {
                (
                    first.stereo_bond_parity[0],
                    first.stereo_bond_neighbor,
                    first.stereo_bond_ord[0],
                )
            };
            let mut cumulene_len =
                (i32::from(encoded_parity) & MASK_CUMULENE_LEN as i32) / MULT_STEREOBOND as i32;
            if cumulene_len % 2 != 0 && (1 >= MAX_NUM_STEREO_BONDS || stereo_neighbor[1] == 0) {
                let order = usize::try_from(stereo_order)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut next_j = usize::from(
                    *first
                        .neighbor
                        .get(order)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                cumulene_len = (cumulene_len - 1) / 2;
                while cumulene_len != 0 {
                    let next_atom = atoms
                        .get(next_j)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if next_atom.valence != 2 {
                        break;
                    }
                    let next_neigh = usize::from(j == usize::from(next_atom.neighbor[0]));
                    j = next_j;
                    next_j = usize::from(next_atom.neighbor[next_neigh]);
                    cumulene_len -= 1;
                }
                if cumulene_len == 0 {
                    Some(
                        *heap
                            .slice(pCanonRank)?
                            .get(next_j)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    )
                } else {
                    None
                }
            } else {
                None
            }
        };

        if let Some(at_num) = allene {
            let count = usize::try_from(*nNumberOfStereoCenters)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let insertion_index = {
                let primary_numbers = heap.slice(Stereo.nNumber.as_const())?;
                if primary_numbers.len() < count {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                let mut index = 0;
                while index < count && primary_numbers[index] < at_num {
                    index += 1;
                }
                index
            };
            let number_pointer = if use_primary_arrays {
                Stereo.nNumber
            } else {
                Stereo.nNumberInv
            };
            let parity_pointer = if use_primary_arrays {
                Stereo.t_parity
            } else {
                Stereo.t_parityInv
            };
            {
                let numbers = heap.slice_mut(number_pointer)?;
                if numbers.len() <= count {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                if insertion_index < count {
                    numbers.copy_within(insertion_index..count, insertion_index + 1);
                }
                numbers[insertion_index] = at_num;
            }
            {
                let parities = heap.slice_mut(parity_pointer)?;
                if parities.len() <= count {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                if insertion_index < count {
                    parities.copy_within(insertion_index..count, insertion_index + 1);
                }
                parities[insertion_index] = LinearCTStereoDble.parity as S_CHAR;
            }
            *nNumberOfStereoCenters = (*nNumberOfStereoCenters).wrapping_add(1);
            return Ok(1);
        }
    }

    if let Some(number_of_bonds) = nNumberOfStereoBonds.as_deref_mut() {
        let index =
            usize::try_from(*number_of_bonds).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(Stereo.b_parity)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = LinearCTStereoDble.parity as S_CHAR;
        *heap
            .slice_mut(Stereo.nBondAtom1)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = LinearCTStereoDble.at_num1;
        *heap
            .slice_mut(Stereo.nBondAtom2)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = LinearCTStereoDble.at_num2;
        *number_of_bonds = (*number_of_bonds).wrapping_add(1);
    }

    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CopyLinearCTStereoToINChIStereo(
    heap: &mut SourceHeap,
    Stereo: Option<&mut INChI_Stereo>,
    LinearCTStereoCarb: SourceConstPointer<AT_STEREO_CARB>,
    nLenLinearCTStereoCarb: i32,
    LinearCTStereoDble: SourceConstPointer<AT_STEREO_DBLE>,
    nLenLinearCTStereoDble: i32,
    pCanonOrd: SourceConstPointer<AT_NUMB>,
    pCanonRank: SourceConstPointer<AT_RANK>,
    at: SourceConstPointer<sp_ATOM>,
    bIsotopic: i32,
    LinearCTStereoCarbInv: SourceConstPointer<AT_STEREO_CARB>,
    LinearCTStereoDbleInv: SourceConstPointer<AT_STEREO_DBLE>,
    pCanonOrdInv: SourceConstPointer<AT_NUMB>,
    pCanonRankInv: SourceConstPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:566 CopyLinearCTStereoToINChIStereo
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int CopyLinearCTStereoToINChIStereo( INChI_Stereo *Stereo,
                                         AT_STEREO_CARB *LinearCTStereoCarb,
                                         int nLenLinearCTStereoCarb,
                                         AT_STEREO_DBLE *LinearCTStereoDble,
                                         int nLenLinearCTStereoDble,
                                         AT_NUMB *pCanonOrd, AT_RANK *pCanonRank,
                                         sp_ATOM *at,
                                         int bIsotopic,
                                         AT_STEREO_CARB *LinearCTStereoCarbInv,
                                         AT_STEREO_DBLE *LinearCTStereoDbleInv,
                                         AT_NUMB *pCanonOrdInv, AT_RANK *pCanonRankInv )
    {
        int n, i, nErrorCode = 0, len;
        int bAllene;
        int diff;
        int lenInv, bAlleneInv;

        /* Stereo centers */

        /* djb-rwth: fixing oss-fuzz issue #68271 */
        if (Stereo)
        {
            n = Stereo->nNumberOfStereoCenters = nLenLinearCTStereoCarb;

            for (i = 0; i < n; i++)
            {
                if (LinearCTStereoCarb) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    Stereo->nNumber[i] = LinearCTStereoCarb[i].at_num;
                    Stereo->t_parity[i] = LinearCTStereoCarb[i].parity;
                }
                if (LinearCTStereoCarbInv) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    Stereo->nNumberInv[i] = LinearCTStereoCarbInv[i].at_num;
                    Stereo->t_parityInv[i] = LinearCTStereoCarbInv[i].parity;
                }
            }

            /* Stereo bonds */

            n = nLenLinearCTStereoDble;
            lenInv = Stereo->nNumberOfStereoCenters;

            for (i = len = 0; i < n; i++)
            {
                bAllene =
                    Copy2StereoBondOrAllene(Stereo,
                        &Stereo->nNumberOfStereoCenters,
                        &len,
                        LinearCTStereoDble + i,
                        pCanonOrd, pCanonRank,
                        at,
                        bIsotopic);

                bAlleneInv =
                    Copy2StereoBondOrAllene(Stereo,
                        &lenInv,
                        NULL,
                        LinearCTStereoDbleInv + i,
                        pCanonOrdInv, pCanonRankInv,
                        at,
                        bIsotopic);

                /* make sure double bond stereo is identical in original and inverted geometry */
                /* Note: all allenes are AFTER double bonds in LinearCTStereoDble... */
                if (bAllene != bAlleneInv || (!bAllene &&
                    CompareLinCtStereoDble(LinearCTStereoDble + i, 1,
                        LinearCTStereoDbleInv + i, 1))) /* djb-rwth: addressing LLVM warning */
                {
                    /* double bond stereo Inv is NOT identical to Abs */
                    nErrorCode = -4;
                    goto exit_function;
                }
            }

            Stereo->nNumberOfStereoBonds = len;

            if (lenInv != Stereo->nNumberOfStereoCenters)
            {
                nErrorCode = -5; /* different number of stereo centers in Abs and Inv */
                goto exit_function;
            }


            /* compare inverted stereocenters to absolute */

            n = Stereo->nNumberOfStereoCenters;
            /* djb-rwth: removing redundant code */

            for (i = 0, diff = 0; i < n; i++)
            {
                if (Stereo->nNumberInv[i] != Stereo->nNumber[i])
                {
                    diff = (Stereo->nNumberInv[i] > Stereo->nNumber[i]) ? 2 : -2;
                    break; /* Abs != Inv */
                }
                if (Stereo->t_parityInv[i] != Stereo->t_parity[i])
                {
                    diff = (Stereo->t_parityInv[i] > Stereo->t_parity[i]) ? 1 : -1;
                    break; /* Abs != Inv */
                }
            }

            Stereo->nCompInv2Abs =
                (diff > 0) ? 1 : (diff < 0) ? -1 : 0;

            if (diff == -1 || diff == 1)
            {
                /* The first found difference was in parities */
                for (i = 0, diff = 0; i < n; i++)
                {
                    if (Stereo->nNumberInv[i] != Stereo->nNumber[i])
                    {
                        diff = 2; /* difference in stereo center numbering */
                        break;
                    }

                    /*  Parities can be only 1, 2, 3, 4. Therefore only mutually inverted pairs
                     *  (t_parityInv, t_parity) = (1,2) or (2,1) statisfy conditions
                     *  (t_parityInv != t_parity) && (t_parityInv + t_parity == 3)
                     */

                    if (Stereo->t_parityInv[i] == Stereo->t_parity[i] ||
                        Stereo->t_parityInv[i] + Stereo->t_parity[i] != 3)
                    {
                        diff = 1; /* parities are same or different and cannot be obtained by simple inversion */
                        break;
                    }
                }
                Stereo->bTrivialInv = !diff;
            }
            else
            {
                Stereo->bTrivialInv = 0;
            }
        }
        else
        {
            nErrorCode = 1;
        }

    exit_function:
        return nErrorCode;
    }
    */
    // END INCHI C FUNCTION: CopyLinearCTStereoToINChIStereo

    let Some(stereo) = Stereo else {
        return Ok(1);
    };
    stereo.nNumberOfStereoCenters = nLenLinearCTStereoCarb;
    let center_count = nLenLinearCTStereoCarb.max(0) as usize;
    for i in 0..center_count {
        if !LinearCTStereoCarb.is_null() {
            let value = heap
                .slice(LinearCTStereoCarb)?
                .get(i)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(stereo.nNumber)?
                .get_mut(i)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = value.at_num;
            *heap
                .slice_mut(stereo.t_parity)?
                .get_mut(i)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = value.parity as S_CHAR;
        }
        if !LinearCTStereoCarbInv.is_null() {
            let value = heap
                .slice(LinearCTStereoCarbInv)?
                .get(i)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(stereo.nNumberInv)?
                .get_mut(i)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = value.at_num;
            *heap
                .slice_mut(stereo.t_parityInv)?
                .get_mut(i)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = value.parity as S_CHAR;
        }
    }

    let mut len = 0_i32;
    let mut len_inv = stereo.nNumberOfStereoCenters;
    let bond_count = nLenLinearCTStereoDble.max(0) as usize;
    for i in 0..bond_count {
        let direct = heap
            .slice(LinearCTStereoDble)?
            .get(i)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let inverted = heap
            .slice(LinearCTStereoDbleInv)?
            .get(i)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut updated_centers = stereo.nNumberOfStereoCenters;
        let allene = Copy2StereoBondOrAllene(
            heap,
            stereo,
            &mut updated_centers,
            Some(&mut len),
            &direct,
            pCanonOrd,
            pCanonRank,
            at,
            bIsotopic,
        )?;
        stereo.nNumberOfStereoCenters = updated_centers;
        let allene_inv = Copy2StereoBondOrAllene(
            heap,
            stereo,
            &mut len_inv,
            None,
            &inverted,
            pCanonOrdInv,
            pCanonRankInv,
            at,
            bIsotopic,
        )?;
        if allene != allene_inv
            || (allene == 0
                && CompareLinCtStereoDble(
                    heap,
                    LinearCTStereoDble.offset(i as i64)?,
                    1,
                    LinearCTStereoDbleInv.offset(i as i64)?,
                    1,
                )? != 0)
        {
            return Ok(-4);
        }
    }

    stereo.nNumberOfStereoBonds = len;
    if len_inv != stereo.nNumberOfStereoCenters {
        return Ok(-5);
    }

    let n = stereo.nNumberOfStereoCenters.max(0) as usize;
    let mut diff = 0_i32;
    for i in 0..n {
        let direct_number = *heap
            .slice(stereo.nNumber.as_const())?
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let inverted_number = *heap
            .slice(stereo.nNumberInv.as_const())?
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if inverted_number != direct_number {
            diff = if inverted_number > direct_number {
                2
            } else {
                -2
            };
            break;
        }
        let direct_parity = *heap
            .slice(stereo.t_parity.as_const())?
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let inverted_parity = *heap
            .slice(stereo.t_parityInv.as_const())?
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if inverted_parity != direct_parity {
            diff = if inverted_parity > direct_parity {
                1
            } else {
                -1
            };
            break;
        }
    }
    stereo.nCompInv2Abs = diff.signum();
    if diff == -1 || diff == 1 {
        diff = 0;
        for i in 0..n {
            let direct_number = heap.slice(stereo.nNumber.as_const())?[i];
            let inverted_number = heap.slice(stereo.nNumberInv.as_const())?[i];
            if inverted_number != direct_number {
                diff = 2;
                break;
            }
            let direct_parity = i32::from(heap.slice(stereo.t_parity.as_const())?[i]);
            let inverted_parity = i32::from(heap.slice(stereo.t_parityInv.as_const())?[i]);
            if inverted_parity == direct_parity || inverted_parity + direct_parity != 3 {
                diff = 1;
                break;
            }
        }
        stereo.bTrivialInv = i32::from(diff == 0);
    } else {
        stereo.bTrivialInv = 0;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MarkAmbiguousStereo(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    norm_at: SourceMutPointer<inp_ATOM>,
    bIsotopic: i32,
    pCanonOrd: SourceConstPointer<AT_NUMB>,
    LinearCTStereoCarb: SourceConstPointer<AT_STEREO_CARB>,
    nLenLinearCTStereoCarb: i32,
    LinearCTStereoDble: SourceConstPointer<AT_STEREO_DBLE>,
    nLenLinearCTStereoDble: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:713 MarkAmbiguousStereo
    // INCHI✔️❌: complete behavior; SourceHeap requires target-array clones.
    /*
    int MarkAmbiguousStereo( sp_ATOM *at,
                             inp_ATOM *norm_at,
                             int bIsotopic,
                             AT_NUMB *pCanonOrd,
                             AT_STEREO_CARB *LinearCTStereoCarb,
                             int nLenLinearCTStereoCarb,
                             AT_STEREO_DBLE *LinearCTStereoDble,
                             int nLenLinearCTStereoDble )
    {
        int n, i, j1, j2, num, mark_atom, mark_bond;

        if (!pCanonOrd)
        {
            return -1;
        }

        num = 0;

        n = nLenLinearCTStereoCarb;

        mark_atom = bIsotopic ? AMBIGUOUS_STEREO_ATOM_ISO : AMBIGUOUS_STEREO_ATOM;

        for (i = 0; i < n; i++)
        {
            /*  Mark ambiguous stereo centers (for displaying and "Ambiguous stereo" message) */
            if (ATOM_PARITY_NOT_UNKN( LinearCTStereoCarb[i].parity ) &&
                 at[j1 = pCanonOrd[(int) LinearCTStereoCarb[i].at_num - 1]].bAmbiguousStereo)
            {
                at[j1].bAmbiguousStereo |= mark_atom;
                norm_at[j1].bAmbiguousStereo |= mark_atom;
                num++;
            }
        }

        n = nLenLinearCTStereoDble;

        mark_bond = bIsotopic ? AMBIGUOUS_STEREO_BOND_ISO : AMBIGUOUS_STEREO_BOND;

        for (i = 0; i < n; i++)
        {
            /*  Mark ambiguous stereo bonds or allenes (for displaying and "Ambiguous stereo" message) */

            if (ATOM_PARITY_WELL_DEF( LinearCTStereoDble[i].parity ))
            {
                j1 = pCanonOrd[(int) LinearCTStereoDble[i].at_num1 - 1];
                j2 = pCanonOrd[(int) LinearCTStereoDble[i].at_num2 - 1];

                if (at[j1].bAmbiguousStereo || at[j2].bAmbiguousStereo)
                {
                    /* If it is an allene then mark the central atom only
                       because the bonds should not be marked to avoid misleading
                       message "Ambiguous stereo: bond(s)": Allene makes a stereocenter
                    */

                    int j1_parity = bIsotopic ? at[j1].stereo_bond_parity2[0] : at[j1].stereo_bond_parity[0];

                    int cumulene_len = BOND_CHAIN_LEN( j1_parity ); /* 0 => double bond, 1 => allene, 2 => cumulene,..*/

                    if (cumulene_len % 2 && ( 1 >= MAX_NUM_STEREO_BONDS ||
                        !( bIsotopic ? at[j1].stereo_bond_neighbor2[1] : at[j1].stereo_bond_neighbor[1] ) )
                        )
                    {
                        /*  found an allene; locate its central atom */

                        int next_j, next_neigh;
                        int j = j1;

                        next_j = at[j].neighbor[bIsotopic ? at[j].stereo_bond_ord2[0] : at[j].stereo_bond_ord[0]];

                        for (cumulene_len = ( cumulene_len - 1 ) / 2;
                             cumulene_len && 2 == at[next_j].valence;
                             cumulene_len--)
                        {
                            next_neigh = ( j == at[next_j].neighbor[0] );
                            j = next_j;
                            next_j = at[next_j].neighbor[next_neigh];
                        }
                        /* next_j is the central atom */
                        if (2 == at[next_j].valence)
                        {
                            at[next_j].bAmbiguousStereo |= mark_atom;
                            norm_at[next_j].bAmbiguousStereo |= mark_atom;
                            num++;
                            continue; /* do not mark the cumulene "bond" endpoints */
                        }
                    }

                    /* Not an allene, mark double bond or cumulene end atoms */
                    if (at[j1].bAmbiguousStereo)
                    {
                        at[j1].bAmbiguousStereo |= mark_bond; /*  ??? */
                        norm_at[j1].bAmbiguousStereo |= mark_bond;
                        num++;
                    }

                    if (at[j2].bAmbiguousStereo)
                    {
                        at[j2].bAmbiguousStereo |= mark_bond; /*  ??? */
                        norm_at[j2].bAmbiguousStereo |= mark_bond;
                        num++;
                    }
                }
            }
        }

        return num;
    }
    */
    // END INCHI C FUNCTION: MarkAmbiguousStereo
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MarkAmbiguousStereo
    // INCHI✔️❌: #define ATOM_PARITY_KNOWN(X)        (AB_MIN_KNOWN_PARITY <= (X) && (X) <= AB_MAX_KNOWN_PARITY)
    // INCHI✔️❌: #define ATOM_PARITY_WELL_DEF(X)     (AB_MIN_WELL_DEFINED_PARITY <= (X) && (X) <= AB_MAX_WELL_DEFINED_PARITY)
    // INCHI✔️❌: #define ATOM_PARITY_NOT_UNKN(X)     (ATOM_PARITY_KNOWN(X) && (X) != AB_PARITY_UNKN)
    // INCHI✔️❌: #define BOND_CHAIN_LEN(X)           (GET_BITS_CUMULENE_LEN(X)/MULT_STEREOBOND) /* 0 => double bond, 1 => allene, 2 => cumulene,..*/
    // INCHI✔️❌: #define GET_BITS_CUMULENE_LEN(X)    ((X)&MASK_CUMULENE_LEN)
    // END INCHI ACTIVE MACRO CONFIGURATION: MarkAmbiguousStereo

    if pCanonOrd.is_null() {
        return Ok(-1);
    }

    let canonical_order = heap.slice(pCanonOrd)?.to_vec();
    let center_count = nLenLinearCTStereoCarb.max(0) as usize;
    let centers = if center_count == 0 {
        Vec::new()
    } else {
        heap.slice(LinearCTStereoCarb)?
            .get(..center_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let bond_count = nLenLinearCTStereoDble.max(0) as usize;
    let bonds = if bond_count == 0 {
        Vec::new()
    } else {
        heap.slice(LinearCTStereoDble)?
            .get(..bond_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let mut atoms = heap.slice(at.as_const())?.to_vec();
    let mut normalized_atoms = heap.slice(norm_at.as_const())?.to_vec();

    let canonical_atom = |number: AT_NUMB| -> Result<usize, SourceHeapError> {
        let index = usize::from(number)
            .checked_sub(1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        canonical_order
            .get(index)
            .copied()
            .map(usize::from)
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };
    let mark_atom = if bIsotopic != 0 {
        AMBIGUOUS_STEREO_ATOM_ISO as i32
    } else {
        AMBIGUOUS_STEREO_ATOM as i32
    };
    let mark_bond = if bIsotopic != 0 {
        AMBIGUOUS_STEREO_BOND_ISO as i32
    } else {
        AMBIGUOUS_STEREO_BOND as i32
    };
    let mut num = 0_i32;

    for center in centers {
        let parity = i32::from(center.parity);
        let parity_not_unknown = (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32)
            .contains(&parity)
            && parity != AB_PARITY_UNKN as i32;
        if parity_not_unknown {
            let j1 = canonical_atom(center.at_num)?;
            let atom = atoms
                .get_mut(j1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom.bAmbiguousStereo != 0 {
                atom.bAmbiguousStereo = (i32::from(atom.bAmbiguousStereo) | mark_atom) as S_CHAR;
                let normalized = normalized_atoms
                    .get_mut(j1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                normalized.bAmbiguousStereo =
                    (i32::from(normalized.bAmbiguousStereo) | mark_atom) as S_CHAR;
                num = num.wrapping_add(1);
            }
        }
    }

    for bond in bonds {
        let parity = i32::from(bond.parity);
        if !(AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
            .contains(&parity)
        {
            continue;
        }
        let j1 = canonical_atom(bond.at_num1)?;
        let j2 = canonical_atom(bond.at_num2)?;
        let endpoint_is_ambiguous = atoms
            .get(j1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .bAmbiguousStereo
            != 0
            || atoms
                .get(j2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bAmbiguousStereo
                != 0;
        if !endpoint_is_ambiguous {
            continue;
        }

        let first = atoms.get(j1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let (encoded_parity, stereo_neighbor, stereo_order) = if bIsotopic != 0 {
            (
                first.stereo_bond_parity2[0],
                first.stereo_bond_neighbor2,
                first.stereo_bond_ord2[0],
            )
        } else {
            (
                first.stereo_bond_parity[0],
                first.stereo_bond_neighbor,
                first.stereo_bond_ord[0],
            )
        };
        let mut cumulene_len =
            (i32::from(encoded_parity) & MASK_CUMULENE_LEN as i32) / MULT_STEREOBOND as i32;
        if cumulene_len % 2 != 0 && (1 >= MAX_NUM_STEREO_BONDS || stereo_neighbor[1] == 0) {
            let order =
                usize::try_from(stereo_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut j = j1;
            let mut next_j = usize::from(
                *atoms
                    .get(j)
                    .and_then(|atom| atom.neighbor.get(order))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            cumulene_len = (cumulene_len - 1) / 2;
            while cumulene_len != 0
                && atoms
                    .get(next_j)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .valence
                    == 2
            {
                let next_atom = atoms
                    .get(next_j)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let next_neigh = usize::from(j == usize::from(next_atom.neighbor[0]));
                j = next_j;
                next_j = usize::from(next_atom.neighbor[next_neigh]);
                cumulene_len -= 1;
            }
            if atoms
                .get(next_j)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .valence
                == 2
            {
                let atom = atoms
                    .get_mut(next_j)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                atom.bAmbiguousStereo = (i32::from(atom.bAmbiguousStereo) | mark_atom) as S_CHAR;
                let normalized = normalized_atoms
                    .get_mut(next_j)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                normalized.bAmbiguousStereo =
                    (i32::from(normalized.bAmbiguousStereo) | mark_atom) as S_CHAR;
                num = num.wrapping_add(1);
                continue;
            }
        }

        if atoms[j1].bAmbiguousStereo != 0 {
            atoms[j1].bAmbiguousStereo =
                (i32::from(atoms[j1].bAmbiguousStereo) | mark_bond) as S_CHAR;
            normalized_atoms[j1].bAmbiguousStereo =
                (i32::from(normalized_atoms[j1].bAmbiguousStereo) | mark_bond) as S_CHAR;
            num = num.wrapping_add(1);
        }
        if atoms[j2].bAmbiguousStereo != 0 {
            atoms[j2].bAmbiguousStereo =
                (i32::from(atoms[j2].bAmbiguousStereo) | mark_bond) as S_CHAR;
            normalized_atoms[j2].bAmbiguousStereo =
                (i32::from(normalized_atoms[j2].bAmbiguousStereo) | mark_bond) as S_CHAR;
            num = num.wrapping_add(1);
        }
    }

    heap.slice_mut(at)?.clone_from_slice(&atoms);
    heap.slice_mut(norm_at)?.clone_from_slice(&normalized_atoms);
    Ok(num)
}

#[allow(non_snake_case)]
pub(crate) fn UnmarkAllUndefinedUnknownStereo(
    heap: &mut SourceHeap,
    Stereo: Option<&mut INChI_Stereo>,
    nUserMode: INCHI_MODE,
) -> Result<INCHI_MODE, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:823 UnmarkAllUndefinedUnknownStereo
    // INCHI✔️❌: complete behavior; SourceHeap adds repeated checked array lookups.
    /*
    INCHI_MODE UnmarkAllUndefinedUnknownStereo( INChI_Stereo *Stereo,
                                                INCHI_MODE nUserMode )
    {
        INCHI_MODE nRet = 0;
        int   i, n;

        if (!Stereo || (Stereo && !Stereo->nNumberOfStereoCenters && !Stereo->nNumberOfStereoBonds)) /* djb-rwth: addressing LLVM warning */
        {
            return nRet;
        }

        /* Stereocenters */
        if (!Stereo->nCompInv2Abs &&
            ( n = Stereo->nNumberOfStereoCenters ) > 0 && ( nUserMode & REQ_MODE_SC_IGN_ALL_UU ))
        {

            for (i = 0; i < n && !ATOM_PARITY_WELL_DEF( Stereo->t_parity[i] ); i++)
            {
                ;
            }

            if (i == n)
            {
                Stereo->nNumberOfStereoCenters = 0;

                for (i = 0; i < n; i++)
                {
                    Stereo->t_parity[i] = 0;
                    Stereo->nNumber[i] = 0;
                    Stereo->t_parityInv[i] = 0;
                    Stereo->nNumberInv[i] = 0;
                }

                nRet |= REQ_MODE_SC_IGN_ALL_UU;
            }
        }

        /* Stereobonds */

        if (( n = Stereo->nNumberOfStereoBonds ) > 0 && ( nUserMode & REQ_MODE_SB_IGN_ALL_UU ))
        {
            for (i = 0; i < n && !ATOM_PARITY_WELL_DEF( Stereo->b_parity[i] ); i++)
            {
                ;
            }

            if (i == n)
            {
                Stereo->nNumberOfStereoBonds = 0;
                for (i = 0; i < n; i++)
                {
                    Stereo->b_parity[i] = 0;
                    Stereo->nBondAtom1[i] = 0;
                    Stereo->nBondAtom2[i] = 0;
                }
                nRet |= REQ_MODE_SB_IGN_ALL_UU;
            }
        }

        return nRet;
    }
    */
    // END INCHI C FUNCTION: UnmarkAllUndefinedUnknownStereo
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: UnmarkAllUndefinedUnknownStereo
    // INCHI✔️❌: #define ATOM_PARITY_WELL_DEF(X)     (AB_MIN_WELL_DEFINED_PARITY <= (X) && (X) <= AB_MAX_WELL_DEFINED_PARITY)
    // INCHI✔️❌: #define REQ_MODE_SC_IGN_ALL_UU      0x000800    /* IAUSC Ignore stereocenters if All Undef/Unknown */
    // INCHI✔️❌: #define REQ_MODE_SB_IGN_ALL_UU      0x001000    /* IAUSC Ignore stereobonds if All Undef/Unknown */
    // END INCHI ACTIVE MACRO CONFIGURATION: UnmarkAllUndefinedUnknownStereo

    let Some(stereo) = Stereo else {
        return Ok(0);
    };
    if stereo.nNumberOfStereoCenters == 0 && stereo.nNumberOfStereoBonds == 0 {
        return Ok(0);
    }

    let mut result = 0 as INCHI_MODE;
    let center_count = stereo.nNumberOfStereoCenters;
    if stereo.nCompInv2Abs == 0
        && center_count > 0
        && nUserMode & REQ_MODE_SC_IGN_ALL_UU as INCHI_MODE != 0
    {
        let center_count = center_count as usize;
        let all_undefined_or_unknown = heap
            .slice(stereo.t_parity.as_const())?
            .get(..center_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .iter()
            .all(|parity| {
                !(AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
                    .contains(&i32::from(*parity))
            });
        if all_undefined_or_unknown {
            stereo.nNumberOfStereoCenters = 0;
            for index in 0..center_count {
                *heap
                    .slice_mut(stereo.t_parity)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                *heap
                    .slice_mut(stereo.nNumber)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                *heap
                    .slice_mut(stereo.t_parityInv)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                *heap
                    .slice_mut(stereo.nNumberInv)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            }
            result |= REQ_MODE_SC_IGN_ALL_UU as INCHI_MODE;
        }
    }

    let bond_count = stereo.nNumberOfStereoBonds;
    if bond_count > 0 && nUserMode & REQ_MODE_SB_IGN_ALL_UU as INCHI_MODE != 0 {
        let bond_count = bond_count as usize;
        let all_undefined_or_unknown = heap
            .slice(stereo.b_parity.as_const())?
            .get(..bond_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .iter()
            .all(|parity| {
                !(AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
                    .contains(&i32::from(*parity))
            });
        if all_undefined_or_unknown {
            stereo.nNumberOfStereoBonds = 0;
            for index in 0..bond_count {
                *heap
                    .slice_mut(stereo.b_parity)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                *heap
                    .slice_mut(stereo.nBondAtom1)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                *heap
                    .slice_mut(stereo.nBondAtom2)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            }
            result |= REQ_MODE_SB_IGN_ALL_UU as INCHI_MODE;
        }
    }

    Ok(result)
}

fn fill_out_source_get<T: Clone + 'static>(
    heap: &SourceHeap,
    pointer: SourceConstPointer<T>,
    index: usize,
) -> Result<T, SourceHeapError> {
    heap.slice(pointer)?
        .get(index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

fn fill_out_source_set<T: Clone + 'static>(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<T>,
    index: usize,
    value: T,
) -> Result<(), SourceHeapError> {
    *heap
        .slice_mut(pointer)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
    Ok(())
}

fn fill_out_copy_prefix<T: Clone + 'static>(
    heap: &mut SourceHeap,
    destination: SourceMutPointer<T>,
    source: SourceConstPointer<T>,
    count: usize,
) -> Result<(), SourceHeapError> {
    let values = heap
        .slice(source)?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(destination)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone_from_slice(&values);
    Ok(())
}

fn fill_out_sort_equivalence(
    heap: &mut SourceHeap,
    globals: &mut CANON_GLOBALS,
    equivalence: SourceMutPointer<AT_NUMB>,
    sort_order: SourceMutPointer<AT_NUMB>,
    count: usize,
) -> Result<(), SourceHeapError> {
    globals.m_pn_RankForSort = equivalence.as_const();
    if count > 1 {
        heap.with_slice_mut_and_heap(sort_order, |order, heap| {
            let order = order
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let bytes = bytemuck::cast_slice_mut::<AT_NUMB, u8>(order);
            inchi_qsort(
                bytes,
                count,
                std::mem::size_of::<AT_NUMB>(),
                &mut |first, second| {
                    let first = AT_NUMB::from_ne_bytes(
                        first
                            .try_into()
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    );
                    let second = AT_NUMB::from_ne_bytes(
                        second
                            .try_into()
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    );
                    CompRanksOrd(heap, first, second, globals)
                },
            )
        })?;
    }

    let mut i = 0_usize;
    let mut minimum_order = fill_out_source_get(heap, sort_order.as_const(), 0)?;
    for j in 1..=count {
        let boundary = if j == count {
            true
        } else {
            let left_order = fill_out_source_get(heap, sort_order.as_const(), i)?;
            let right_order = fill_out_source_get(heap, sort_order.as_const(), j)?;
            fill_out_source_get(heap, equivalence.as_const(), usize::from(left_order))?
                != fill_out_source_get(heap, equivalence.as_const(), usize::from(right_order))?
        };
        if boundary {
            minimum_order = minimum_order.wrapping_add(1);
            if j - i > 1 {
                while i < j {
                    let order = fill_out_source_get(heap, sort_order.as_const(), i)?;
                    fill_out_source_set(heap, equivalence, usize::from(order), minimum_order)?;
                    i += 1;
                }
            } else {
                let order = fill_out_source_get(heap, sort_order.as_const(), i)?;
                fill_out_source_set(heap, equivalence, usize::from(order), 0)?;
                i += 1;
            }
            if j < count {
                minimum_order = fill_out_source_get(heap, sort_order.as_const(), j)?;
            }
        }
    }
    Ok(())
}

fn fill_out_warning(
    error_buffer: Option<&mut [i8]>,
    message: &[u8],
) -> Result<(), SourceHeapError> {
    let message: Vec<i8> = message.iter().map(|byte| *byte as i8).collect();
    AddErrorMessage(error_buffer, Some(&message))?;
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillOutINChI(
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
    mut pStrErrStruct: Option<&mut [i8]>,
    bNoWarnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:1072 FillOutINChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int FillOutINChI( INChI *pINChI,
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
                      char *pStrErrStruct,
                      int bNoWarnings )
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
    
        /*  Total charge */
        /* djb-rwth: fixing oss-fuzz issue #69656 */
        for (i = 0, n = 0; i < num_atoms + num_removed_H; i++)
        {
            n += at[i].charge;
        }
        pINChI->nTotalCharge = n;
    
        /*  Number of atoms */
        pINChI->nNumberOfAtoms = num_atoms;
        pINChI_Aux->nNumberOfAtoms = num_atoms;
    
        /* Removed protons and detachable isotopic H */
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
            if (pINChI_Aux->bNormalizationFlags & ( FLAG_NORM_CONSIDER_TAUT &~FLAG_PROTON_CHARGE_CANCEL ))
            {
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErrStruct, "Proton(s) added/removed" );
                }
            }
            if (pINChI_Aux->bNormalizationFlags & FLAG_PROTON_CHARGE_CANCEL)
            {
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErrStruct, "Charges neutralized" );
                }
            }
        }
    
        /* Abs or rel stereo may establish one of two canonical numberings */
        if (( pCS->nLenLinearCTStereoCarb > 0 || pCS->nLenLinearCTStereoDble > 0 ) &&
              pCS->nLenCanonOrdStereo > 0 &&
              ( (pCS->LinearCTStereoCarb && pCS->LinearCTStereoCarbInv) ||
                  (pCS->LinearCTStereoDble && pCS->LinearCTStereoDbleInv) ) &&
              pCS->nCanonOrdStereo && pCS->nCanonOrdStereoInv
           ) /* djb-rwth: addressing LLVM warning */
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
            /* Copy stereo bonds and stereo centers; compare Inv and Abs stereo */
            /********************************************************************/
    
            nErrorCode = CopyLinearCTStereoToINChIStereo( Stereo,
                                                          pCS->LinearCTStereoCarb,
                                                          pCS->nLenLinearCTStereoCarb,
                                                          pCS->LinearCTStereoDble,
                                                          pCS->nLenLinearCTStereoDble,
                                                          pCanonOrd, pCanonRank,
                                                          at, 0 /* non-isotopic */,
                                                          pCS->LinearCTStereoCarbInv,
                                                          pCS->LinearCTStereoDbleInv,
                                                          pCanonOrdInv, pCanonRankInv );
    
            /* djb-rwth: fixing oss-fuzz issue #68271 */
            if (nErrorCode == 1)
            {
                nErrorCode = 0;
                ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
                pINChI->nErrorCode = pINChI_Aux->nErrorCode = CT_OUT_OF_RAM;
                goto exit_function;
            }
    
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
    
            LOG_NO_ARGS("************************** Canonical Ordering with Stereo (L1183:ichimak2.c) ***************************\n");
            for (i = 0; i < num_atoms; i++)
            {
                pINChI_Aux->nOrigAtNosInCanonOrdInv[i] = at[pCanonOrdInv[i]].orig_at_number;
                pINChI_Aux->nOrigAtNosInCanonOrd[i] = at[pCanonOrd[i]].orig_at_number;
                LOG_MULT_ARGS("Atom Nr: %d, Canonical Numbering Normal: %d, Element Name: %s\n", i + 1, at[pCanonOrd[i]].orig_at_number, at[i].elname);
            }
    
            LOG_NO_ARGS("\n******************************************************************************************************\n");
    
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
        {
            /*------------------------------ no stereo */
            pCanonOrd = pCS->nLenCanonOrdStereo > 0 ? pCS->nCanonOrdStereo
                                                    : pCS->nLenCanonOrd > 0  ? pCS->nCanonOrd
                                                                             : NULL;
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
    
        /* Non-isotopic canonical numbers and equivalence of atoms (Aux) */
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
    
    
        /*  Atomic numbers from the Periodic Table */
        for (i = 0; i < num_atoms; i++)
        {
            pINChI->nAtom[i] = (int) at[pCanonOrd[i]].el_number;
        }
    
    
        /*  Connection table: atoms only (before 7-29-2003 pCS->LinearCT2 contained non-isotopic CT) */
        if (pCS->nLenLinearCTAtOnly <= 0 || !pCS->LinearCT || !pINChI->nConnTable)
        {
            nErrorCode |= ERR_NO_CANON_RESULTS;
            ret = -2;
            goto exit_function;
        }
    
        memcpy(pINChI->nConnTable, pCS->LinearCT, sizeof(pINChI->nConnTable[0]) * pCS->nLenLinearCTAtOnly);
    
        pINChI->lenConnTable = pCS->nLenLinearCTAtOnly;
    
        /*  Tautomeric group(s) canonical representation */
        len = 0;
        if (bTautomeric && 0 < ( n = SortTautomerGroupsAndEndpoints(
            pCG, t_group_info,
            num_atoms, num_at_tg, pCanonRank ) ))
        {
            /* SortTautomerGroupsAndEndpoints() produces canonically ordered t-groups */
            pINChI->nFlags |=
                ( t_group_info->bTautFlagsDone & TG_FLAG_ALL_SALT_DONE ) ? INCHI_FLAG_ACID_TAUT
                : 0;
    
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
    
                for (j = 0; j < INCHI_T_NUM_MOVABLE; j++) /* djb-rwth: redundant condition; && j < T_NUM_NO_ISOTOPIC part should be checked; can INCHI_T_NUM_MOVABLE and T_NUM_NO_ISOTOPIC change values from 2? */
                {
                    pINChI->nTautomer[len++] = t_group->num[j];
                }
                
                
                /* djb-rwth: erroneous loop execution condition
                
                for (j = T_NUM_NO_ISOTOPIC; j < INCHI_T_NUM_MOVABLE; j++)
                {
                    pINChI->nTautomer[len++] = 0; // should not happen 
                }
                
                */
    
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
    
        /*  Number of H (excluding tautomeric) */
        if (pCS->nNum_H)
        {
            for (i = 0; i < num_atoms; i++)
            {
                pINChI->nNum_H[i] = pCS->nNum_H[i];
            }
        }
    
    
        /*  Number of fixed H (tautomeric H in non-tautomeric representation) */
        if (pCS->nNum_H_fixed && !pINChI->lenTautomer)
        {
            for (i = 0; i < num_atoms; i++)
            {
                pINChI->nNum_H_fixed[i] = pCS->nNum_H_fixed[i];
                pINChI->nNum_H[i] += pCS->nNum_H_fixed[i];
            }
        }
    
    
        /***********************************************************
         *  Tautomeric group(s) numbering and symmetry;
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
    
        pINChI->szHillFormula = AllocateAndFillHillFormula( pINChI );
    
        if (!pINChI->szHillFormula)
        {
            nErrorCode = 0;
            ret = CT_WRONG_FORMULA; /* CT_OUT_OF_RAM;*/  /*   <BRKPT> */
            pINChI->nErrorCode = pINChI_Aux->nErrorCode = ret;
            goto exit_function;
        }
    
        nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo( pINChI->Stereo, nUserMode );
    
        if (nStereoUnmarkMode)
        {
            pINChI->nFlags |=
                ( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ? INCHI_FLAG_SC_IGN_ALL_UU
                : 0;
            pINChI->nFlags |=
                ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ) ? INCHI_FLAG_SB_IGN_ALL_UU
                : 0;
    
            if (( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ||
                ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ))
            {
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErrStruct, "Omitted undefined stereo" );
                }
            }
        }
    
        /*************************/
        /* Mark ambiguous stereo */
        /*************************/
    
        MarkAmbiguousStereo( at, norm_at, 0 /* non-isotopic */, pCanonOrd,
                             pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb,
                             pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble );
    
    
        /************************************************************************
         *
         *  Isotopic part
         */
    
        /* abs or rel stereo may establish one of two canonical numberings */
    
        if (( pCS->nLenLinearCTIsotopicStereoCarb > 0 || pCS->nLenLinearCTIsotopicStereoDble > 0 ) &&
              pCS->nLenCanonOrdIsotopicStereo > 0 &&
              ( (pCS->LinearCTIsotopicStereoCarb && pCS->LinearCTIsotopicStereoCarbInv) ||
                  (pCS->LinearCTIsotopicStereoDble && pCS->LinearCTIsotopicStereoDbleInv) ) &&
              pCS->nCanonOrdIsotopicStereo && pCS->nCanonOrdIsotopicStereoInv
              ) /* djb-rwth: addressing LLVM warning */
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
                                                          pCS->LinearCTIsotopicStereoCarb,
                                                          pCS->nLenLinearCTIsotopicStereoCarb,
                                                          pCS->LinearCTIsotopicStereoDble,
                                                          pCS->nLenLinearCTIsotopicStereoDble,
                                                          pCanonOrd, pCanonRank, at, 1 /* isotopic */,
                                                          pCS->LinearCTIsotopicStereoCarbInv,
                                                          pCS->LinearCTIsotopicStereoDbleInv,
                                                          pCanonOrdInv, pCanonRankInv );
    
            /* djb-rwth: fixing oss-fuzz issue #68271 */
            if (nErrorCode == 1)
            {
                nErrorCode = 0;
                ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
                pINChI->nErrorCode = pINChI_Aux->nErrorCode = CT_OUT_OF_RAM;
                goto exit_function;
            }
            
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
    
            pCanonOrd = pCS->nLenCanonOrdIsotopicStereo > 0 ? pCS->nCanonOrdIsotopicStereo
                : pCS->nLenCanonOrdIsotopic > 0 ? pCS->nCanonOrdIsotopic
                : NULL;
            pCanonRank = pCanonRankAtoms;
            pOrigNosInCanonOrd = pINChI_Aux->nIsotopicOrigAtNosInCanonOrd;
    
            /* djb-rwth: fixing oss-fuzz issue #30496 */
            if (pCanonOrd && pCanonRank && pOrigNosInCanonOrd)
            {
                for (i = 0; i < num_atoms; i++)
                {
                    /* Fix13 -- out of bounds */
                    pCanonRank[pCanonOrd[i]] = (AT_NUMB) ( i + 1 );
                    pOrigNosInCanonOrd[i] = at[pCanonOrd[i]].orig_at_number;
                }
                for (; i < num_at_tg; i++)
                {
                    /* Fix13 -- out of bounds */
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
    
    
        /*  Isotopic atoms */
    
        n = pINChI->nNumberOfIsotopicAtoms = pCS->nLenLinearCTIsotopic;
    
        for (i = 0; i < n; i++)
        {
            pINChI->IsotopicAtom[i].nAtomNumber = pCS->LinearCTIsotopic[i].at_num;
            pINChI->IsotopicAtom[i].nIsoDifference = pCS->LinearCTIsotopic[i].iso_atw_diff;
            pINChI->IsotopicAtom[i].nNum_H = pCS->LinearCTIsotopic[i].num_1H;
            pINChI->IsotopicAtom[i].nNum_D = pCS->LinearCTIsotopic[i].num_D;
            pINChI->IsotopicAtom[i].nNum_T = pCS->LinearCTIsotopic[i].num_T;
        }
    
        /*  Isotopic tautomeric groups */
    
        n = pINChI->nNumberOfIsotopicTGroups = pCS->nLenLinearCTIsotopicTautomer;
    
        for (i = 0; i < n; i++)
        {
            pINChI->IsotopicTGroup[i].nTGroupNumber = pCS->LinearCTIsotopicTautomer[i].tgroup_num;
            pINChI->IsotopicTGroup[i].nNum_H = pCS->LinearCTIsotopicTautomer[i].num[2];
            pINChI->IsotopicTGroup[i].nNum_D = pCS->LinearCTIsotopicTautomer[i].num[1];
            pINChI->IsotopicTGroup[i].nNum_T = pCS->LinearCTIsotopicTautomer[i].num[0];
        }
    
        /* Atoms that may exchange isotopic H-atoms */
    
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
            pINChI->nFlags |=
                ( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ? INCHI_FLAG_SC_IGN_ALL_ISO_UU
                : 0;
            pINChI->nFlags |=
                ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ) ? INCHI_FLAG_SC_IGN_ALL_ISO_UU
                : 0;
            if (( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ||
                ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ))
            {
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErrStruct, "Omitted undefined stereo" );
                }
            }
        }
    
        /* Mark ambiguous stereo */
    
        MarkAmbiguousStereo( at, norm_at, 1 /* isotopic */, pCanonOrd,
                             pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb,
                             pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble );
    
    
        /***********************************************************
         *  Isotopic tautomeric group(s) numbering and symmetry;
         *  should not depend on switching to rel. stereo numbering
         */
    
        if (pINChI->lenTautomer &&
             pINChI_Aux->nConstitEquIsotopicTGroupNumbers &&
             pCS->nSymmRankIsotopicTaut &&
             ( pCS->nLenLinearCTIsotopic || pCS->nLenLinearCTIsotopicTautomer ) &&
             t_group_info && t_group_info->num_t_groups > 0)
        {
    
            /* djb-rwth: removing redundant code */
    
            pCanonOrdTaut =
                pCS->nLenCanonOrdIsotopicStereoTaut > 0 ? ( n = pCS->nLenCanonOrdIsotopicStereoTaut, pCS->nCanonOrdIsotopicStereoTaut )
                : pCS->nLenCanonOrdIsotopicTaut > 0 ? ( n = pCS->nLenCanonOrdIsotopicTaut, pCS->nCanonOrdIsotopicTaut )
                : ( n = 0, (AT_RANK*) NULL );
    
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
    // END INCHI C FUNCTION: FillOutINChI

    let mut pCanonRankAtoms = SourceMutPointer::<AT_NUMB>::null();
    let mut pSortOrd = SourceMutPointer::<AT_NUMB>::null();
    let mut nErrorCode = 0_i32;

    let operation = (|| -> Result<i32, SourceHeapError> {
        let mut ret = 0_i32;
        if pCS.nLenLinearCTStereoCarb < 0
            || pCS.nLenLinearCTStereoDble < 0
            || pCS.nLenCanonOrdStereo < 0
            || pCS.nLenCanonOrdStereoTaut < 0
        {
            nErrorCode |= WARN_FAILED_STEREO as i32;
        }
        if pCS.nLenLinearCTIsotopic < 0
            || pCS.nLenLinearCTIsotopicTautomer < 0
            || pCS.nLenCanonOrdIsotopic < 0
            || pCS.nLenCanonOrdIsotopicTaut < 0
        {
            nErrorCode |= WARN_FAILED_ISOTOPIC as i32;
        }
        if pCS.nLenLinearCTIsotopicStereoCarb < 0
            || pCS.nLenLinearCTIsotopicStereoDble < 0
            || pCS.nLenCanonOrdIsotopicStereo < 0
            || pCS.nLenCanonOrdIsotopicStereoTaut < 0
        {
            nErrorCode |= WARN_FAILED_ISOTOPIC_STEREO as i32;
        }

        let allocation_count = i64::from(num_at_tg)
            .checked_add(1)
            .and_then(|value| u64::try_from(value).ok())
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        pCanonRankAtoms = match inchi_calloc::<AT_NUMB>(heap, allocation_count, 2) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        pSortOrd = match inchi_calloc::<AT_NUMB>(heap, allocation_count, 2) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if pCanonRankAtoms.is_null() || pSortOrd.is_null() {
            nErrorCode = 0;
            pINChI.nErrorCode = CT_OUT_OF_RAM;
            pINChI_Aux.nErrorCode = CT_OUT_OF_RAM;
            return Ok(CT_OUT_OF_RAM);
        }

        let atom_count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atom_and_removed_count = num_atoms
            .checked_add(num_removed_H)
            .and_then(|value| usize::try_from(value).ok())
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let atoms = heap
            .slice(at.as_const())?
            .get(..atom_and_removed_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut total_charge = 0_i32;
        for atom in atoms {
            total_charge = total_charge
                .checked_add(i32::from(atom.charge))
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
        pINChI.nTotalCharge = total_charge;
        pINChI.nNumberOfAtoms = num_atoms;
        pINChI_Aux.nNumberOfAtoms = num_atoms;

        let t_group_info = if pCS.t_group_info.is_null() {
            None
        } else {
            Some(
                heap.slice(pCS.t_group_info.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        if bTautomeric != 0
            && let Some(t_group_info) = t_group_info.as_ref()
        {
            pINChI_Aux.nNumRemovedProtons = t_group_info.tni.nNumRemovedProtons;
            for i in 0..NUM_H_ISOTOPES as usize {
                pINChI_Aux.nNumRemovedIsotopicH[i] = t_group_info.num_iso_H[i]
                    .wrapping_add(t_group_info.tni.nNumRemovedProtonsIsotopic[i]);
            }
            if pINChI_Aux.bNormalizationFlags & u64::from(FLAG_FORCE_SALT_TAUT) != 0 {
                pINChI.nFlags |= u64::from(INCHI_FLAG_HARD_ADD_REM_PROTON);
            }
            if pINChI_Aux.bNormalizationFlags
                & u64::from(FLAG_NORM_CONSIDER_TAUT & !FLAG_PROTON_CHARGE_CANCEL)
                != 0
                && bNoWarnings == 0
            {
                fill_out_warning(pStrErrStruct.as_deref_mut(), b"Proton(s) added/removed\0")?;
            }
            if pINChI_Aux.bNormalizationFlags & u64::from(FLAG_PROTON_CHARGE_CANCEL) != 0
                && bNoWarnings == 0
            {
                fill_out_warning(pStrErrStruct.as_deref_mut(), b"Charges neutralized\0")?;
            }
        }

        let num_at_tg_usize =
            usize::try_from(num_at_tg).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let mut pCanonRank = pCanonRankAtoms;
        let mut pCanonOrd = SourceMutPointer::<AT_NUMB>::null();
        let mut pCanonRankInv = SourceMutPointer::<AT_NUMB>::null();
        let mut pCanonOrdInv = SourceMutPointer::<AT_NUMB>::null();
        let mut bUseNumberingInv = 0_i32;

        let has_nonisotopic_stereo = (pCS.nLenLinearCTStereoCarb > 0
            || pCS.nLenLinearCTStereoDble > 0)
            && pCS.nLenCanonOrdStereo > 0
            && ((!pCS.LinearCTStereoCarb.is_null() && !pCS.LinearCTStereoCarbInv.is_null())
                || (!pCS.LinearCTStereoDble.is_null() && !pCS.LinearCTStereoDbleInv.is_null()))
            && !pCS.nCanonOrdStereo.is_null()
            && !pCS.nCanonOrdStereoInv.is_null();
        if has_nonisotopic_stereo {
            pCanonOrd = pCS.nCanonOrdStereo;
            pCanonRankInv = pSortOrd;
            pCanonOrdInv = pCS.nCanonOrdStereoInv;
            for i in 0..num_at_tg_usize {
                let rank = (i as AT_NUMB).wrapping_add(1);
                let direct = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                let inverted = fill_out_source_get(heap, pCanonOrdInv.as_const(), i)?;
                fill_out_source_set(heap, pCanonRank, usize::from(direct), rank)?;
                fill_out_source_set(heap, pCanonRankInv, usize::from(inverted), rank)?;
            }
            let mut stereo = heap
                .slice(pINChI.Stereo.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            nErrorCode = CopyLinearCTStereoToINChIStereo(
                heap,
                Some(&mut stereo),
                pCS.LinearCTStereoCarb.as_const(),
                pCS.nLenLinearCTStereoCarb,
                pCS.LinearCTStereoDble.as_const(),
                pCS.nLenLinearCTStereoDble,
                pCanonOrd.as_const(),
                pCanonRank.as_const(),
                at.as_const(),
                0,
                pCS.LinearCTStereoCarbInv.as_const(),
                pCS.LinearCTStereoDbleInv.as_const(),
                pCanonOrdInv.as_const(),
                pCanonRankInv.as_const(),
            )?;
            if nErrorCode == 1 {
                nErrorCode = 0;
                pINChI.nErrorCode = CT_OUT_OF_RAM;
                pINChI_Aux.nErrorCode = CT_OUT_OF_RAM;
                return Ok(CT_OUT_OF_RAM);
            }
            if !stereo.t_parityInv.is_null() && !stereo.nNumberInv.is_null() {
                if nUserMode & u64::from(REQ_MODE_RELATIVE_STEREO) != 0 {
                    pINChI.nFlags |= u64::from(INCHI_FLAG_REL_STEREO);
                }
                if nUserMode & u64::from(REQ_MODE_RACEMIC_STEREO) != 0 {
                    pINChI.nFlags |= u64::from(INCHI_FLAG_RAC_STEREO);
                }
                if stereo.nCompInv2Abs == -1 {
                    std::mem::swap(&mut stereo.nNumber, &mut stereo.nNumberInv);
                    std::mem::swap(&mut stereo.t_parity, &mut stereo.t_parityInv);
                    switch_ptrs(&mut pCanonRank, &mut pCanonRankInv);
                    switch_ptrs(&mut pCanonOrd, &mut pCanonOrdInv);
                    bUseNumberingInv = 1;
                }
            }
            heap.slice_mut(pINChI.Stereo)?[0] = stereo;
            for i in 0..atom_count {
                let direct = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                let inverted = fill_out_source_get(heap, pCanonOrdInv.as_const(), i)?;
                let direct_original =
                    fill_out_source_get(heap, at.as_const(), usize::from(direct))?.orig_at_number;
                let inverted_original =
                    fill_out_source_get(heap, at.as_const(), usize::from(inverted))?.orig_at_number;
                fill_out_source_set(heap, pINChI_Aux.nOrigAtNosInCanonOrd, i, direct_original)?;
                fill_out_source_set(
                    heap,
                    pINChI_Aux.nOrigAtNosInCanonOrdInv,
                    i,
                    inverted_original,
                )?;
            }
            if bUseNumberingInv != 0 {
                switch_ptrs(&mut pCanonRank, &mut pCanonRankInv);
                switch_ptrs(&mut pCanonOrd, &mut pCanonOrdInv);
                fill_out_copy_prefix(heap, pCanonRank, pCanonRankInv.as_const(), num_at_tg_usize)?;
                fill_out_copy_prefix(heap, pCanonOrd, pCanonOrdInv.as_const(), num_at_tg_usize)?;
            }
        } else {
            pCanonOrd = if pCS.nLenCanonOrdStereo > 0 {
                pCS.nCanonOrdStereo
            } else if pCS.nLenCanonOrd > 0 {
                pCS.nCanonOrd
            } else {
                SourceMutPointer::null()
            };
            if !pCanonOrd.is_null() {
                for i in 0..atom_count {
                    let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                    fill_out_source_set(
                        heap,
                        pCanonRank,
                        usize::from(order),
                        (i as AT_NUMB).wrapping_add(1),
                    )?;
                    let original = fill_out_source_get(heap, at.as_const(), usize::from(order))?
                        .orig_at_number;
                    fill_out_source_set(heap, pINChI_Aux.nOrigAtNosInCanonOrd, i, original)?;
                }
                for i in atom_count..num_at_tg_usize {
                    let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                    fill_out_source_set(
                        heap,
                        pCanonRank,
                        usize::from(order),
                        (i as AT_NUMB).wrapping_add(1),
                    )?;
                }
            }
        }

        if !pINChI_Aux.OrigInfo.is_null() {
            for i in 0..atom_count {
                if pCanonOrd.is_null() {
                    continue;
                }
                let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                let atom = fill_out_source_get(heap, norm_at.as_const(), usize::from(order))?;
                let radical = if atom.radical == RADICAL_SINGLET as i8 {
                    0
                } else if atom.radical == RADICAL_DOUBLET as i8 {
                    1
                } else if atom.radical == RADICAL_TRIPLET as i8 {
                    2
                } else if atom.radical != 0 {
                    3
                } else {
                    0
                };
                let mut info = fill_out_source_get(heap, pINChI_Aux.OrigInfo.as_const(), i)?;
                info.cRadical = radical;
                if atom.valence != 0 || atom.num_H != 0 {
                    info.cCharge = atom.charge;
                    info.cUnusualValence = get_unusual_el_valence(
                        i32::from(atom.el_number),
                        i32::from(atom.charge),
                        i32::from(atom.radical),
                        i32::from(atom.chem_bonds_valence),
                        i32::from(atom.num_H),
                        i32::from(atom.valence),
                    )? as i8;
                }
                fill_out_source_set(heap, pINChI_Aux.OrigInfo, i, info)?;
            }
        }

        let pConstitEquNumb = pINChI_Aux.nConstitEquNumbers;
        if pCanonOrd.is_null()
            || pCanonRank.is_null()
            || pCS.nSymmRank.is_null()
            || pConstitEquNumb.is_null()
        {
            nErrorCode |= ERR_NO_CANON_RESULTS as i32;
            return Ok(-1);
        }
        for i in 0..atom_count {
            let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
            let rank = fill_out_source_get(heap, pCS.nSymmRank.as_const(), usize::from(order))?;
            fill_out_source_set(heap, pConstitEquNumb, i, rank)?;
            fill_out_source_set(heap, pSortOrd, i, i as AT_NUMB)?;
        }
        for i in atom_count..num_at_tg_usize {
            fill_out_source_set(heap, pSortOrd, i, MAX_ATOMS as AT_NUMB)?;
        }
        fill_out_sort_equivalence(heap, pCG, pConstitEquNumb, pSortOrd, atom_count)?;

        for i in 0..atom_count {
            let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
            let atomic_number =
                fill_out_source_get(heap, at.as_const(), usize::from(order))?.el_number;
            fill_out_source_set(heap, pINChI.nAtom, i, atomic_number)?;
        }
        if pCS.nLenLinearCTAtOnly <= 0 || pCS.LinearCT.is_null() || pINChI.nConnTable.is_null() {
            nErrorCode |= ERR_NO_CANON_RESULTS as i32;
            return Ok(-2);
        }
        let connection_count = usize::try_from(pCS.nLenLinearCTAtOnly)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        fill_out_copy_prefix(
            heap,
            pINChI.nConnTable,
            pCS.LinearCT.as_const(),
            connection_count,
        )?;
        pINChI.lenConnTable = pCS.nLenLinearCTAtOnly;

        let mut tautomer_length = 0_usize;
        let tautomer_groups = if bTautomeric != 0 {
            if let Some(info) = t_group_info.as_ref() {
                SortTautomerGroupsAndEndpoints(heap, pCG, info, num_atoms, num_at_tg, pCanonRank)?
            } else {
                0
            }
        } else {
            0
        };
        if tautomer_groups > 0 {
            let info = t_group_info
                .as_ref()
                .expect("positive group count requires info");
            if info.bTautFlagsDone & u64::from(TG_FLAG_ALL_SALT_DONE) != 0 {
                pINChI.nFlags |= u64::from(INCHI_FLAG_ACID_TAUT);
            }
            fill_out_source_set(
                heap,
                pINChI.nTautomer,
                tautomer_length,
                tautomer_groups as AT_NUMB,
            )?;
            tautomer_length += 1;
            for i in 0..usize::try_from(tautomer_groups)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            {
                let group_index =
                    usize::from(fill_out_source_get(heap, info.tGroupNumber.as_const(), i)?);
                let group = fill_out_source_get(heap, info.t_group.as_const(), group_index)?;
                fill_out_source_set(
                    heap,
                    pINChI.nTautomer,
                    tautomer_length,
                    group
                        .nNumEndpoints
                        .wrapping_add(INCHI_T_NUM_MOVABLE as AT_NUMB),
                )?;
                tautomer_length += 1;
                for j in 0..INCHI_T_NUM_MOVABLE as usize {
                    fill_out_source_set(heap, pINChI.nTautomer, tautomer_length, group.num[j])?;
                    tautomer_length += 1;
                }
                let first = usize::from(group.nFirstEndpointAtNoPos);
                let end = first + usize::from(group.nNumEndpoints);
                for j in first..end {
                    let endpoint =
                        fill_out_source_get(heap, info.nEndpointAtomNumber.as_const(), j)?;
                    let rank =
                        fill_out_source_get(heap, pCanonRank.as_const(), usize::from(endpoint))?;
                    fill_out_source_set(heap, pINChI.nTautomer, tautomer_length, rank)?;
                    tautomer_length += 1;
                }
            }
            pINChI.lenTautomer = tautomer_length as i32;
            pINChI_Aux.nNumberOfTGroups = tautomer_groups;
        } else {
            pINChI.lenTautomer = 0;
            pINChI_Aux.nNumberOfTGroups = 0;
            if let Some(info) = t_group_info.as_ref()
                && (info.tni.bNormalizationFlags & u64::from(FLAG_NORM_CONSIDER_TAUT) != 0
                    || (info.nNumIsotopicEndpoints > 1
                        && info.bTautFlagsDone
                            & u64::from(
                                TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE,
                            )
                            != 0))
            {
                pINChI.lenTautomer = 1;
                fill_out_source_set(heap, pINChI.nTautomer, 0, 0)?;
            }
        }

        if !pCS.nNum_H.is_null() {
            fill_out_copy_prefix(heap, pINChI.nNum_H, pCS.nNum_H.as_const(), atom_count)?;
        }
        if !pCS.nNum_H_fixed.is_null() && pINChI.lenTautomer == 0 {
            for i in 0..atom_count {
                let fixed = fill_out_source_get(heap, pCS.nNum_H_fixed.as_const(), i)?;
                fill_out_source_set(heap, pINChI.nNum_H_fixed, i, fixed)?;
                let normal = fill_out_source_get(heap, pINChI.nNum_H.as_const(), i)?;
                fill_out_source_set(heap, pINChI.nNum_H, i, normal.wrapping_add(fixed))?;
            }
        }

        if pINChI.lenTautomer != 0 && pINChI_Aux.nNumberOfTGroups != 0 {
            let group_count = usize::try_from(pINChI_Aux.nNumberOfTGroups)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let pCanonOrdTaut = if pCS.nLenCanonOrdStereoTaut > 0 {
                pCS.nCanonOrdStereoTaut
            } else if pCS.nLenCanonOrdTaut > 0 {
                pCS.nCanonOrdTaut
            } else {
                SourceMutPointer::null()
            };
            if !pCanonOrdTaut.is_null()
                && !pCS.nSymmRankTaut.is_null()
                && !pINChI_Aux.nConstitEquTGroupNumbers.is_null()
            {
                for i in 0..group_count {
                    let order = fill_out_source_get(heap, pCanonOrdTaut.as_const(), i)?;
                    let rank = fill_out_source_get(
                        heap,
                        pCS.nSymmRankTaut.as_const(),
                        usize::from(order),
                    )?;
                    fill_out_source_set(heap, pINChI_Aux.nConstitEquTGroupNumbers, i, rank)?;
                    fill_out_source_set(heap, pSortOrd, i, i as AT_NUMB)?;
                }
                fill_out_sort_equivalence(
                    heap,
                    pCG,
                    pINChI_Aux.nConstitEquTGroupNumbers,
                    pSortOrd,
                    group_count,
                )?;
            }
        }

        pINChI.szHillFormula = AllocateAndFillHillFormula(heap, pINChI)?;
        if pINChI.szHillFormula.is_null() {
            nErrorCode = 0;
            ret = CT_WRONG_FORMULA;
            pINChI.nErrorCode = ret;
            pINChI_Aux.nErrorCode = ret;
            return Ok(ret);
        }

        let mut stereo = if pINChI.Stereo.is_null() {
            None
        } else {
            Some(
                heap.slice(pINChI.Stereo.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo(heap, stereo.as_mut(), nUserMode)?;
        if let Some(stereo) = stereo {
            heap.slice_mut(pINChI.Stereo)?[0] = stereo;
        }
        if nStereoUnmarkMode != 0 {
            if nStereoUnmarkMode & u64::from(REQ_MODE_SC_IGN_ALL_UU) != 0 {
                pINChI.nFlags |= u64::from(INCHI_FLAG_SC_IGN_ALL_UU);
            }
            if nStereoUnmarkMode & u64::from(REQ_MODE_SB_IGN_ALL_UU) != 0 {
                pINChI.nFlags |= u64::from(INCHI_FLAG_SB_IGN_ALL_UU);
            }
            if bNoWarnings == 0 {
                fill_out_warning(pStrErrStruct.as_deref_mut(), b"Omitted undefined stereo\0")?;
            }
        }
        MarkAmbiguousStereo(
            heap,
            at,
            norm_at,
            0,
            pCanonOrd.as_const(),
            pCS.LinearCTStereoCarb.as_const(),
            pCS.nLenLinearCTStereoCarb,
            pCS.LinearCTStereoDble.as_const(),
            pCS.nLenLinearCTStereoDble,
        )?;

        let mut bUseIsotopicNumberingInv = 0_i32;
        let has_isotopic_stereo = (pCS.nLenLinearCTIsotopicStereoCarb > 0
            || pCS.nLenLinearCTIsotopicStereoDble > 0)
            && pCS.nLenCanonOrdIsotopicStereo > 0
            && ((!pCS.LinearCTIsotopicStereoCarb.is_null()
                && !pCS.LinearCTIsotopicStereoCarbInv.is_null())
                || (!pCS.LinearCTIsotopicStereoDble.is_null()
                    && !pCS.LinearCTIsotopicStereoDbleInv.is_null()))
            && !pCS.nCanonOrdIsotopicStereo.is_null()
            && !pCS.nCanonOrdIsotopicStereoInv.is_null();
        if has_isotopic_stereo {
            pCanonRank = pCanonRankAtoms;
            pCanonOrd = pCS.nCanonOrdIsotopicStereo;
            pCanonRankInv = pSortOrd;
            pCanonOrdInv = pCS.nCanonOrdIsotopicStereoInv;
            for i in 0..num_at_tg_usize {
                let rank = (i as AT_NUMB).wrapping_add(1);
                let direct = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                let inverted = fill_out_source_get(heap, pCanonOrdInv.as_const(), i)?;
                fill_out_source_set(heap, pCanonRank, usize::from(direct), rank)?;
                fill_out_source_set(heap, pCanonRankInv, usize::from(inverted), rank)?;
            }
            let mut stereo = heap
                .slice(pINChI.StereoIsotopic.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            nErrorCode = CopyLinearCTStereoToINChIStereo(
                heap,
                Some(&mut stereo),
                pCS.LinearCTIsotopicStereoCarb.as_const(),
                pCS.nLenLinearCTIsotopicStereoCarb,
                pCS.LinearCTIsotopicStereoDble.as_const(),
                pCS.nLenLinearCTIsotopicStereoDble,
                pCanonOrd.as_const(),
                pCanonRank.as_const(),
                at.as_const(),
                1,
                pCS.LinearCTIsotopicStereoCarbInv.as_const(),
                pCS.LinearCTIsotopicStereoDbleInv.as_const(),
                pCanonOrdInv.as_const(),
                pCanonRankInv.as_const(),
            )?;
            if nErrorCode == 1 {
                nErrorCode = 0;
                pINChI.nErrorCode = CT_OUT_OF_RAM;
                pINChI_Aux.nErrorCode = CT_OUT_OF_RAM;
                return Ok(CT_OUT_OF_RAM);
            }
            if !stereo.t_parityInv.is_null() && !stereo.nNumberInv.is_null() {
                if nUserMode & u64::from(REQ_MODE_RELATIVE_STEREO) != 0 {
                    pINChI.nFlags |= u64::from(INCHI_FLAG_REL_STEREO);
                }
                if nUserMode & u64::from(REQ_MODE_RACEMIC_STEREO) != 0 {
                    pINChI.nFlags |= u64::from(INCHI_FLAG_RAC_STEREO);
                }
                if stereo.nCompInv2Abs == -1 {
                    std::mem::swap(&mut stereo.nNumber, &mut stereo.nNumberInv);
                    std::mem::swap(&mut stereo.t_parity, &mut stereo.t_parityInv);
                    switch_ptrs(&mut pCanonRank, &mut pCanonRankInv);
                    switch_ptrs(&mut pCanonOrd, &mut pCanonOrdInv);
                    bUseIsotopicNumberingInv = 1;
                }
            }
            heap.slice_mut(pINChI.StereoIsotopic)?[0] = stereo;
            for i in 0..atom_count {
                let direct = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                let inverted = fill_out_source_get(heap, pCanonOrdInv.as_const(), i)?;
                let direct_original =
                    fill_out_source_get(heap, at.as_const(), usize::from(direct))?.orig_at_number;
                let inverted_original =
                    fill_out_source_get(heap, at.as_const(), usize::from(inverted))?.orig_at_number;
                fill_out_source_set(
                    heap,
                    pINChI_Aux.nIsotopicOrigAtNosInCanonOrd,
                    i,
                    direct_original,
                )?;
                fill_out_source_set(
                    heap,
                    pINChI_Aux.nIsotopicOrigAtNosInCanonOrdInv,
                    i,
                    inverted_original,
                )?;
            }
            if bUseIsotopicNumberingInv != 0 {
                switch_ptrs(&mut pCanonRank, &mut pCanonRankInv);
                switch_ptrs(&mut pCanonOrd, &mut pCanonOrdInv);
                fill_out_copy_prefix(heap, pCanonRank, pCanonRankInv.as_const(), num_at_tg_usize)?;
                fill_out_copy_prefix(heap, pCanonOrd, pCanonOrdInv.as_const(), num_at_tg_usize)?;
            }
        } else {
            pCanonOrd = if pCS.nLenCanonOrdIsotopicStereo > 0 {
                pCS.nCanonOrdIsotopicStereo
            } else if pCS.nLenCanonOrdIsotopic > 0 {
                pCS.nCanonOrdIsotopic
            } else {
                SourceMutPointer::null()
            };
            pCanonRank = pCanonRankAtoms;
            if !pCanonOrd.is_null() && !pINChI_Aux.nIsotopicOrigAtNosInCanonOrd.is_null() {
                for i in 0..atom_count {
                    let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                    fill_out_source_set(
                        heap,
                        pCanonRank,
                        usize::from(order),
                        (i as AT_NUMB).wrapping_add(1),
                    )?;
                    let original = fill_out_source_get(heap, at.as_const(), usize::from(order))?
                        .orig_at_number;
                    fill_out_source_set(
                        heap,
                        pINChI_Aux.nIsotopicOrigAtNosInCanonOrd,
                        i,
                        original,
                    )?;
                }
                for i in atom_count..num_at_tg_usize {
                    let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
                    fill_out_source_set(
                        heap,
                        pCanonRank,
                        usize::from(order),
                        (i as AT_NUMB).wrapping_add(1),
                    )?;
                }
            }
        }

        if pCanonOrd.is_null()
            || pCanonRank.is_null()
            || pINChI_Aux.nConstitEquIsotopicNumbers.is_null()
            || pCS.nSymmRankIsotopic.is_null()
        {
            return Ok(ret);
        }
        for i in 0..atom_count {
            let order = fill_out_source_get(heap, pCanonOrd.as_const(), i)?;
            let rank =
                fill_out_source_get(heap, pCS.nSymmRankIsotopic.as_const(), usize::from(order))?;
            fill_out_source_set(heap, pINChI_Aux.nConstitEquIsotopicNumbers, i, rank)?;
            fill_out_source_set(heap, pSortOrd, i, i as AT_NUMB)?;
        }
        for i in atom_count..num_at_tg_usize {
            fill_out_source_set(heap, pSortOrd, i, i as AT_NUMB)?;
        }
        fill_out_sort_equivalence(
            heap,
            pCG,
            pINChI_Aux.nConstitEquIsotopicNumbers,
            pSortOrd,
            atom_count,
        )?;

        let isotope_atom_count = usize::try_from(pCS.nLenLinearCTIsotopic.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        pINChI.nNumberOfIsotopicAtoms = pCS.nLenLinearCTIsotopic;
        for i in 0..isotope_atom_count {
            let source = fill_out_source_get(heap, pCS.LinearCTIsotopic.as_const(), i)?;
            let mut destination = fill_out_source_get(heap, pINChI.IsotopicAtom.as_const(), i)?;
            destination.nAtomNumber = source.at_num;
            destination.nIsoDifference = source.iso_atw_diff;
            destination.nNum_H = source.num_1H;
            destination.nNum_D = source.num_D;
            destination.nNum_T = source.num_T;
            fill_out_source_set(heap, pINChI.IsotopicAtom, i, destination)?;
        }
        let isotope_group_count = usize::try_from(pCS.nLenLinearCTIsotopicTautomer.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        pINChI.nNumberOfIsotopicTGroups = pCS.nLenLinearCTIsotopicTautomer;
        for i in 0..isotope_group_count {
            let source = fill_out_source_get(heap, pCS.LinearCTIsotopicTautomer.as_const(), i)?;
            let mut destination = fill_out_source_get(heap, pINChI.IsotopicTGroup.as_const(), i)?;
            destination.nTGroupNumber = source.tgroup_num;
            destination.nNum_H = source.num[2];
            destination.nNum_D = source.num[1];
            destination.nNum_T = source.num[0];
            fill_out_source_set(heap, pINChI.IsotopicTGroup, i, destination)?;
        }
        if !pCS.nExchgIsoH.is_null() && !pINChI.nPossibleLocationsOfIsotopicH.is_null() {
            let mut j = 1_usize;
            for i in 0..atom_count {
                if fill_out_source_get(heap, pCS.nExchgIsoH.as_const(), i)? != 0 {
                    fill_out_source_set(
                        heap,
                        pINChI.nPossibleLocationsOfIsotopicH,
                        j,
                        (i as AT_NUMB).wrapping_add(1),
                    )?;
                    j += 1;
                }
            }
            fill_out_source_set(heap, pINChI.nPossibleLocationsOfIsotopicH, 0, j as AT_NUMB)?;
        }

        let mut isotope_stereo = if pINChI.StereoIsotopic.is_null() {
            None
        } else {
            Some(
                heap.slice(pINChI.StereoIsotopic.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let isotope_unmark =
            UnmarkAllUndefinedUnknownStereo(heap, isotope_stereo.as_mut(), nUserMode)?;
        if let Some(stereo) = isotope_stereo {
            heap.slice_mut(pINChI.StereoIsotopic)?[0] = stereo;
        }
        if isotope_unmark != 0 {
            if isotope_unmark & u64::from(REQ_MODE_SC_IGN_ALL_UU) != 0 {
                pINChI.nFlags |= u64::from(INCHI_FLAG_SC_IGN_ALL_ISO_UU);
            }
            if isotope_unmark & u64::from(REQ_MODE_SB_IGN_ALL_UU) != 0 {
                pINChI.nFlags |= u64::from(INCHI_FLAG_SC_IGN_ALL_ISO_UU);
            }
            if bNoWarnings == 0 {
                fill_out_warning(pStrErrStruct.as_deref_mut(), b"Omitted undefined stereo\0")?;
            }
        }
        MarkAmbiguousStereo(
            heap,
            at,
            norm_at,
            1,
            pCanonOrd.as_const(),
            pCS.LinearCTIsotopicStereoCarb.as_const(),
            pCS.nLenLinearCTIsotopicStereoCarb,
            pCS.LinearCTIsotopicStereoDble.as_const(),
            pCS.nLenLinearCTIsotopicStereoDble,
        )?;

        if pINChI.lenTautomer != 0
            && !pINChI_Aux.nConstitEquIsotopicTGroupNumbers.is_null()
            && !pCS.nSymmRankIsotopicTaut.is_null()
            && (pCS.nLenLinearCTIsotopic != 0 || pCS.nLenLinearCTIsotopicTautomer != 0)
            && let Some(info) = t_group_info.as_ref()
            && info.num_t_groups > 0
        {
            let (count, order) = if pCS.nLenCanonOrdIsotopicStereoTaut > 0 {
                (
                    pCS.nLenCanonOrdIsotopicStereoTaut,
                    pCS.nCanonOrdIsotopicStereoTaut,
                )
            } else if pCS.nLenCanonOrdIsotopicTaut > 0 {
                (pCS.nLenCanonOrdIsotopicTaut, pCS.nCanonOrdIsotopicTaut)
            } else {
                (0, SourceMutPointer::null())
            };
            if !order.is_null() && count > 0 {
                let count =
                    usize::try_from(count).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                for i in 0..count {
                    let order = fill_out_source_get(heap, order.as_const(), i)?;
                    let rank = fill_out_source_get(
                        heap,
                        pCS.nSymmRankIsotopicTaut.as_const(),
                        usize::from(order),
                    )?;
                    fill_out_source_set(
                        heap,
                        pINChI_Aux.nConstitEquIsotopicTGroupNumbers,
                        i,
                        rank,
                    )?;
                    fill_out_source_set(heap, pSortOrd, i, i as AT_NUMB)?;
                }
                fill_out_sort_equivalence(
                    heap,
                    pCG,
                    pINChI_Aux.nConstitEquIsotopicTGroupNumbers,
                    pSortOrd,
                    count,
                )?;
            }
        }
        Ok(ret)
    })();

    let rank_cleanup = inchi_free(heap, pCanonRankAtoms);
    let order_cleanup = inchi_free(heap, pSortOrd);
    pINChI.nErrorCode |= nErrorCode;
    pINChI_Aux.nErrorCode |= nErrorCode;
    let ret = operation?;
    rank_cleanup?;
    order_cleanup?;
    Ok(ret)
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

    struct FillOutFixture {
        heap: SourceHeap,
        inchi: INChI,
        aux: INChI_Aux,
        canon: CANON_STAT,
        globals: CANON_GLOBALS,
        atoms: SourceMutPointer<sp_ATOM>,
        normalized: SourceMutPointer<inp_ATOM>,
        errors: [i8; 256],
    }

    impl FillOutFixture {
        fn new() -> Self {
            let mut heap = SourceHeap::default();
            let mut atoms = vec![sp_ATOM::default(); 2];
            atoms[0].el_number = 6;
            atoms[0].orig_at_number = 11;
            atoms[0].charge = 1;
            atoms[1].el_number = 8;
            atoms[1].orig_at_number = 22;
            atoms[1].charge = -1;
            let atoms = heap.allocate_model_storage(atoms).unwrap();

            let mut normalized = vec![inp_ATOM::default(); 2];
            normalized[0].el_number = 6;
            normalized[0].valence = 4;
            normalized[0].chem_bonds_valence = 4;
            normalized[0].radical = RADICAL_DOUBLET as i8;
            normalized[1].el_number = 8;
            normalized[1].valence = 1;
            normalized[1].chem_bonds_valence = 1;
            normalized[1].charge = -1;
            let normalized = heap.allocate_model_storage(normalized).unwrap();

            let canonical_order = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
            let isotope_order = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
            let symmetry = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
            let isotope_symmetry = heap.allocate_model_storage(vec![7_u16, 7]).unwrap();
            let connection = heap
                .allocate_model_storage(vec![2_u16, 1_u16, 2_u16])
                .unwrap();
            let hydrogens = heap.allocate_model_storage(vec![4_i8, 0_i8]).unwrap();
            let fixed_hydrogens = heap.allocate_model_storage(vec![0_i8, 1_i8]).unwrap();
            let exchange = heap.allocate_model_storage(vec![0_i8, 1_i8]).unwrap();

            let inchi = INChI {
                nAtom: heap.allocate_model_storage(vec![0_u8; 2]).unwrap(),
                nConnTable: heap.allocate_model_storage(vec![0_u16; 3]).unwrap(),
                nTautomer: heap.allocate_model_storage(vec![0_u16; 8]).unwrap(),
                nNum_H: heap.allocate_model_storage(vec![0_i8; 2]).unwrap(),
                nNum_H_fixed: heap.allocate_model_storage(vec![0_i8; 2]).unwrap(),
                IsotopicAtom: heap.allocate_model_storage(Vec::new()).unwrap(),
                IsotopicTGroup: heap.allocate_model_storage(Vec::new()).unwrap(),
                nPossibleLocationsOfIsotopicH: heap.allocate_model_storage(vec![0_u16; 3]).unwrap(),
                ..INChI::default()
            };
            let aux = INChI_Aux {
                nOrigAtNosInCanonOrd: heap.allocate_model_storage(vec![0_u16; 2]).unwrap(),
                nIsotopicOrigAtNosInCanonOrd: heap.allocate_model_storage(vec![0_u16; 2]).unwrap(),
                nOrigAtNosInCanonOrdInv: heap.allocate_model_storage(vec![0_u16; 2]).unwrap(),
                nIsotopicOrigAtNosInCanonOrdInv: heap
                    .allocate_model_storage(vec![0_u16; 2])
                    .unwrap(),
                nConstitEquNumbers: heap.allocate_model_storage(vec![0_u16; 2]).unwrap(),
                nConstitEquTGroupNumbers: heap.allocate_model_storage(vec![0_u16; 2]).unwrap(),
                nConstitEquIsotopicNumbers: heap.allocate_model_storage(vec![0_u16; 2]).unwrap(),
                nConstitEquIsotopicTGroupNumbers: heap
                    .allocate_model_storage(vec![0_u16; 2])
                    .unwrap(),
                OrigInfo: heap
                    .allocate_model_storage(vec![crate::source_types::ORIG_INFO::default(); 2])
                    .unwrap(),
                ..INChI_Aux::default()
            };
            let canon = CANON_STAT {
                LinearCT: connection,
                nLenLinearCTAtOnly: 3,
                nCanonOrd: canonical_order,
                nLenCanonOrd: 2,
                nSymmRank: symmetry,
                nCanonOrdIsotopic: isotope_order,
                nLenCanonOrdIsotopic: 2,
                nSymmRankIsotopic: isotope_symmetry,
                nNum_H: hydrogens,
                nNum_H_fixed: fixed_hydrogens,
                nExchgIsoH: exchange,
                ..CANON_STAT::default()
            };
            Self {
                heap,
                inchi,
                aux,
                canon,
                globals: CANON_GLOBALS::default(),
                atoms,
                normalized,
                errors: [0; 256],
            }
        }

        fn call(&mut self) -> Result<i32, SourceHeapError> {
            FillOutINChI(
                &mut self.heap,
                &mut self.inchi,
                &mut self.aux,
                2,
                2,
                0,
                self.atoms,
                self.normalized,
                &mut self.canon,
                &mut self.globals,
                0,
                0,
                Some(&mut self.errors),
                0,
            )
        }
    }

    #[test]
    fn source_port__ichimak2__filloutinchi__line_1072() {
        let mut fixture = FillOutFixture::new();
        let allocations_before = fixture.heap.live_allocation_count();
        assert_eq!(fixture.call(), Ok(0));
        assert_eq!(fixture.heap.live_allocation_count(), allocations_before + 1);
        assert_eq!(fixture.inchi.nErrorCode, 0);
        assert_eq!(fixture.aux.nErrorCode, 0);
        assert_eq!(fixture.inchi.nTotalCharge, 0);
        assert_eq!(fixture.inchi.nNumberOfAtoms, 2);
        assert_eq!(fixture.aux.nNumberOfAtoms, 2);
        assert_eq!(fixture.inchi.lenConnTable, 3);
        assert_eq!(
            fixture.heap.slice(fixture.inchi.nAtom.as_const()).unwrap(),
            &[6, 8]
        );
        assert_eq!(
            fixture
                .heap
                .slice(fixture.inchi.nConnTable.as_const())
                .unwrap(),
            &[2, 1, 2]
        );
        assert_eq!(
            fixture.heap.slice(fixture.inchi.nNum_H.as_const()).unwrap(),
            &[4, 1]
        );
        assert_eq!(
            fixture
                .heap
                .slice(fixture.aux.nOrigAtNosInCanonOrd.as_const())
                .unwrap(),
            &[11, 22]
        );
        assert_eq!(
            fixture
                .heap
                .slice(fixture.aux.nIsotopicOrigAtNosInCanonOrd.as_const())
                .unwrap(),
            &[22, 11]
        );
        assert_eq!(
            fixture
                .heap
                .slice(fixture.aux.nConstitEquNumbers.as_const())
                .unwrap(),
            &[0, 0]
        );
        assert_eq!(
            fixture
                .heap
                .slice(fixture.aux.nConstitEquIsotopicNumbers.as_const())
                .unwrap(),
            &[1, 1]
        );
        assert_eq!(
            fixture
                .heap
                .slice(fixture.inchi.nPossibleLocationsOfIsotopicH.as_const())
                .unwrap(),
            &[2, 2, 0]
        );
        assert_eq!(
            &fixture
                .heap
                .slice(fixture.inchi.szHillFormula.as_const())
                .unwrap()[..5],
            &[b'C' as i8, b'H' as i8, b'5' as i8, b'O' as i8, 0]
        );
        let original_info = fixture.heap.slice(fixture.aux.OrigInfo.as_const()).unwrap();
        assert_eq!(original_info[0].cRadical, 1);
        assert_eq!(original_info[1].cCharge, -1);

        let mut warnings = FillOutFixture::new();
        warnings.canon.nLenLinearCTStereoCarb = -1;
        warnings.canon.nLenLinearCTIsotopic = -1;
        warnings.canon.nLenLinearCTIsotopicStereoCarb = -1;
        assert_eq!(warnings.call(), Ok(0));
        assert_eq!(
            warnings.inchi.nErrorCode,
            (WARN_FAILED_STEREO | WARN_FAILED_ISOTOPIC | WARN_FAILED_ISOTOPIC_STEREO) as i32
        );
        assert_eq!(warnings.aux.nErrorCode, warnings.inchi.nErrorCode);

        let mut no_canonical_results = FillOutFixture::new();
        no_canonical_results.canon.nSymmRank = SourceMutPointer::null();
        assert_eq!(no_canonical_results.call(), Ok(-1));
        assert_eq!(
            no_canonical_results.inchi.nErrorCode,
            ERR_NO_CANON_RESULTS as i32
        );

        let mut no_connection_table = FillOutFixture::new();
        no_connection_table.canon.nLenLinearCTAtOnly = 0;
        assert_eq!(no_connection_table.call(), Ok(-2));
        assert_eq!(
            no_connection_table.inchi.nErrorCode,
            ERR_NO_CANON_RESULTS as i32
        );

        let mut wrong_formula = FillOutFixture::new();
        wrong_formula.heap.slice_mut(wrong_formula.atoms).unwrap()[0].el_number = 0;
        assert_eq!(wrong_formula.call(), Ok(CT_WRONG_FORMULA));
        assert_eq!(wrong_formula.inchi.nErrorCode, CT_WRONG_FORMULA);
        assert_eq!(wrong_formula.aux.nErrorCode, CT_WRONG_FORMULA);

        for successful_allocations in [0, 1] {
            let mut out_of_memory = FillOutFixture::new();
            let allocations_before = out_of_memory.heap.live_allocation_count();
            out_of_memory
                .heap
                .fail_after_allocations(successful_allocations);
            assert_eq!(out_of_memory.call(), Ok(CT_OUT_OF_RAM));
            assert_eq!(out_of_memory.inchi.nErrorCode, CT_OUT_OF_RAM);
            assert_eq!(out_of_memory.aux.nErrorCode, CT_OUT_OF_RAM);
            assert_eq!(
                out_of_memory.heap.live_allocation_count(),
                allocations_before
            );
        }
    }

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
    fn source_port__ichimak2__unmarkallundefinedunknownstereo__line_823() {
        fn stereo_fixture(
            heap: &mut SourceHeap,
            center_parities: Vec<i8>,
            bond_parities: Vec<i8>,
        ) -> INChI_Stereo {
            let center_count = center_parities.len() as i32;
            let bond_count = bond_parities.len() as i32;
            INChI_Stereo {
                nNumberOfStereoCenters: center_count,
                nNumber: heap
                    .allocate_model_storage((0..center_count).map(|i| i as u16 + 11).collect())
                    .unwrap(),
                t_parity: heap.allocate_model_storage(center_parities).unwrap(),
                nNumberInv: heap
                    .allocate_model_storage((0..center_count).map(|i| i as u16 + 21).collect())
                    .unwrap(),
                t_parityInv: heap
                    .allocate_model_storage((0..center_count).map(|i| i as i8 + 31).collect())
                    .unwrap(),
                nNumberOfStereoBonds: bond_count,
                nBondAtom1: heap
                    .allocate_model_storage((0..bond_count).map(|i| i as u16 + 41).collect())
                    .unwrap(),
                nBondAtom2: heap
                    .allocate_model_storage((0..bond_count).map(|i| i as u16 + 51).collect())
                    .unwrap(),
                b_parity: heap.allocate_model_storage(bond_parities).unwrap(),
                ..INChI_Stereo::default()
            }
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            UnmarkAllUndefinedUnknownStereo(&mut heap, None, INCHI_MODE::MAX),
            Ok(0)
        );
        let mut empty = INChI_Stereo::default();
        assert_eq!(
            UnmarkAllUndefinedUnknownStereo(&mut heap, Some(&mut empty), INCHI_MODE::MAX),
            Ok(0)
        );
        empty.nNumberOfStereoCenters = -1;
        empty.nNumberOfStereoBonds = -2;
        assert_eq!(
            UnmarkAllUndefinedUnknownStereo(&mut heap, Some(&mut empty), INCHI_MODE::MAX),
            Ok(0)
        );

        let mut no_mode = stereo_fixture(&mut heap, vec![3, 4], vec![3, 4]);
        assert_eq!(
            UnmarkAllUndefinedUnknownStereo(&mut heap, Some(&mut no_mode), 0),
            Ok(0)
        );
        assert_eq!(no_mode.nNumberOfStereoCenters, 2);
        assert_eq!(no_mode.nNumberOfStereoBonds, 2);

        let mut mixed = stereo_fixture(&mut heap, vec![3, 1, 4], vec![4, 2, 3]);
        assert_eq!(
            UnmarkAllUndefinedUnknownStereo(
                &mut heap,
                Some(&mut mixed),
                (REQ_MODE_SC_IGN_ALL_UU | REQ_MODE_SB_IGN_ALL_UU) as INCHI_MODE,
            ),
            Ok(0)
        );
        assert_eq!(mixed.nNumberOfStereoCenters, 3);
        assert_eq!(mixed.nNumberOfStereoBonds, 3);
        assert_eq!(heap.slice(mixed.t_parity.as_const()).unwrap(), &[3, 1, 4]);
        assert_eq!(heap.slice(mixed.b_parity.as_const()).unwrap(), &[4, 2, 3]);

        let mut inversion_comparison = stereo_fixture(&mut heap, vec![0, 3, 4], vec![3, 4]);
        inversion_comparison.nCompInv2Abs = -1;
        assert_eq!(
            UnmarkAllUndefinedUnknownStereo(
                &mut heap,
                Some(&mut inversion_comparison),
                INCHI_MODE::MAX,
            ),
            Ok(REQ_MODE_SB_IGN_ALL_UU as INCHI_MODE)
        );
        assert_eq!(inversion_comparison.nNumberOfStereoCenters, 3);
        assert_eq!(inversion_comparison.nNumberOfStereoBonds, 0);

        let mut all_undefined =
            stereo_fixture(&mut heap, vec![0, 3, 4, 5, -1], vec![0, 3, 4, 5, -1]);
        let result =
            UnmarkAllUndefinedUnknownStereo(&mut heap, Some(&mut all_undefined), INCHI_MODE::MAX)
                .unwrap();
        assert_eq!(
            result,
            (REQ_MODE_SC_IGN_ALL_UU | REQ_MODE_SB_IGN_ALL_UU) as INCHI_MODE
        );
        assert_eq!(all_undefined.nNumberOfStereoCenters, 0);
        assert_eq!(all_undefined.nNumberOfStereoBonds, 0);
        assert_eq!(
            heap.slice(all_undefined.t_parity.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(
            heap.slice(all_undefined.nNumber.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(
            heap.slice(all_undefined.t_parityInv.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(
            heap.slice(all_undefined.nNumberInv.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(
            heap.slice(all_undefined.b_parity.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(
            heap.slice(all_undefined.nBondAtom1.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(
            heap.slice(all_undefined.nBondAtom2.as_const()).unwrap(),
            &[0; 5]
        );
    }

    #[test]
    fn source_port__ichimak2__markambiguousstereo__line_713() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            MarkAmbiguousStereo(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                i32::MAX,
                SourceConstPointer::null(),
                i32::MAX,
            ),
            Ok(-1)
        );

        let canonical_order = heap
            .allocate_model_storage(vec![0_u16, 1, 2, 3, 4, 5])
            .unwrap();
        let mut atoms = vec![sp_ATOM::default(); 6];
        atoms[0].bAmbiguousStereo = 1;
        atoms[1].bAmbiguousStereo = 1;
        atoms[2].bAmbiguousStereo = 1;
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let mut normalized = vec![inp_ATOM::default(); 6];
        normalized[0].bAmbiguousStereo = 32;
        let normalized = heap.allocate_model_storage(normalized).unwrap();
        let centers = heap
            .allocate_model_storage(vec![
                AT_STEREO_CARB {
                    at_num: 1,
                    parity: 1,
                },
                AT_STEREO_CARB {
                    at_num: 2,
                    parity: AB_PARITY_UNKN as u8,
                },
                AT_STEREO_CARB {
                    at_num: 3,
                    parity: 4,
                },
                AT_STEREO_CARB {
                    at_num: 4,
                    parity: 0,
                },
                AT_STEREO_CARB {
                    at_num: 5,
                    parity: 5,
                },
            ])
            .unwrap();
        let bonds = heap
            .allocate_model_storage(vec![
                AT_STEREO_DBLE {
                    at_num1: 1,
                    at_num2: 2,
                    parity: 1,
                },
                AT_STEREO_DBLE {
                    at_num1: 3,
                    at_num2: 4,
                    parity: 2,
                },
                AT_STEREO_DBLE {
                    at_num1: 1,
                    at_num2: 2,
                    parity: 3,
                },
                AT_STEREO_DBLE {
                    at_num1: 5,
                    at_num2: 6,
                    parity: 1,
                },
            ])
            .unwrap();
        assert_eq!(
            MarkAmbiguousStereo(
                &mut heap,
                atoms,
                normalized,
                0,
                canonical_order.as_const(),
                centers.as_const(),
                5,
                bonds.as_const(),
                4,
            ),
            Ok(5)
        );
        assert_eq!(
            heap.slice(atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.bAmbiguousStereo)
                .collect::<Vec<_>>(),
            vec![7, 5, 7, 0, 0, 0]
        );
        assert_eq!(
            heap.slice(normalized.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.bAmbiguousStereo)
                .collect::<Vec<_>>(),
            vec![38, 4, 6, 0, 0, 0]
        );
        assert_eq!(
            MarkAmbiguousStereo(
                &mut heap,
                atoms,
                normalized,
                0,
                canonical_order.as_const(),
                SourceConstPointer::null(),
                -1,
                SourceConstPointer::null(),
                i32::MIN,
            ),
            Ok(0)
        );

        fn allene_fixture(
            heap: &mut SourceHeap,
            encoded_parity: i8,
            second_stereo_neighbor: AT_NUMB,
            first_internal_valence: i8,
        ) -> (
            SourceMutPointer<sp_ATOM>,
            SourceMutPointer<inp_ATOM>,
            SourceMutPointer<AT_NUMB>,
            SourceMutPointer<AT_STEREO_DBLE>,
        ) {
            let mut atoms = vec![sp_ATOM::default(); 4];
            atoms[0].bAmbiguousStereo = 1;
            atoms[3].bAmbiguousStereo = 1;
            atoms[0].stereo_bond_parity2[0] = encoded_parity;
            atoms[0].stereo_bond_neighbor2[1] = second_stereo_neighbor;
            atoms[0].stereo_bond_ord2[0] = 0;
            atoms[0].neighbor[0] = 1;
            atoms[1].valence = first_internal_valence;
            atoms[1].neighbor[0] = 0;
            atoms[1].neighbor[1] = 2;
            atoms[2].valence = 2;
            (
                heap.allocate_model_storage(atoms).unwrap(),
                heap.allocate_model_storage(vec![inp_ATOM::default(); 4])
                    .unwrap(),
                heap.allocate_model_storage(vec![0_u16, 1, 2, 3]).unwrap(),
                heap.allocate_model_storage(vec![AT_STEREO_DBLE {
                    at_num1: 1,
                    at_num2: 4,
                    parity: 2,
                }])
                .unwrap(),
            )
        }

        let (atoms, normalized, order, bond) = allene_fixture(&mut heap, 24, 0, 2);
        heap.slice_mut(normalized).unwrap()[2].bAmbiguousStereo = 16;
        assert_eq!(
            MarkAmbiguousStereo(
                &mut heap,
                atoms,
                normalized,
                -1,
                order.as_const(),
                SourceConstPointer::null(),
                0,
                bond.as_const(),
                1,
            ),
            Ok(1)
        );
        assert_eq!(
            heap.slice(atoms.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.bAmbiguousStereo)
                .collect::<Vec<_>>(),
            vec![1, 0, 8, 1]
        );
        assert_eq!(
            heap.slice(normalized.as_const()).unwrap()[2].bAmbiguousStereo,
            24
        );

        let (atoms, normalized, order, bond) = allene_fixture(&mut heap, 8, 0, 2);
        assert_eq!(
            MarkAmbiguousStereo(
                &mut heap,
                atoms,
                normalized,
                1,
                order.as_const(),
                SourceConstPointer::null(),
                0,
                bond.as_const(),
                1,
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap()[1].bAmbiguousStereo, 8);

        for (second_neighbor, internal_valence) in [(0_u16, 1_i8), (9, 2)] {
            let (atoms, normalized, order, bond) =
                allene_fixture(&mut heap, 8, second_neighbor, internal_valence);
            assert_eq!(
                MarkAmbiguousStereo(
                    &mut heap,
                    atoms,
                    normalized,
                    1,
                    order.as_const(),
                    SourceConstPointer::null(),
                    0,
                    bond.as_const(),
                    1,
                ),
                Ok(2)
            );
            assert_eq!(
                heap.slice(atoms.as_const()).unwrap()[0].bAmbiguousStereo,
                17
            );
            assert_eq!(
                heap.slice(atoms.as_const()).unwrap()[3].bAmbiguousStereo,
                17
            );
            assert_eq!(
                heap.slice(normalized.as_const()).unwrap()[0].bAmbiguousStereo,
                16
            );
            assert_eq!(
                heap.slice(normalized.as_const()).unwrap()[3].bAmbiguousStereo,
                16
            );
        }
    }

    #[test]
    fn source_port__ichimak2__copylinearctstereotoinchistereo__line_566() {
        fn stereo_fixture(heap: &mut SourceHeap) -> INChI_Stereo {
            INChI_Stereo {
                nNumber: heap.allocate_model_storage(vec![0_u16; 8]).unwrap(),
                t_parity: heap.allocate_model_storage(vec![0_i8; 8]).unwrap(),
                nNumberInv: heap.allocate_model_storage(vec![0_u16; 8]).unwrap(),
                t_parityInv: heap.allocate_model_storage(vec![0_i8; 8]).unwrap(),
                nBondAtom1: heap.allocate_model_storage(vec![0_u16; 8]).unwrap(),
                nBondAtom2: heap.allocate_model_storage(vec![0_u16; 8]).unwrap(),
                b_parity: heap.allocate_model_storage(vec![0_i8; 8]).unwrap(),
                ..INChI_Stereo::default()
            }
        }

        let null_carb: SourceConstPointer<AT_STEREO_CARB> = SourceConstPointer::null();
        let null_dble: SourceConstPointer<AT_STEREO_DBLE> = SourceConstPointer::null();
        let null_rank: SourceConstPointer<AT_RANK> = SourceConstPointer::null();
        let null_atom: SourceConstPointer<sp_ATOM> = SourceConstPointer::null();
        let mut heap = SourceHeap::default();
        assert_eq!(
            CopyLinearCTStereoToINChIStereo(
                &mut heap, None, null_carb, 0, null_dble, 0, null_rank, null_rank, null_atom, 0,
                null_carb, null_dble, null_rank, null_rank,
            ),
            Ok(1)
        );

        let mut stereo = stereo_fixture(&mut heap);
        let direct_centers = heap
            .allocate_model_storage(vec![
                AT_STEREO_CARB {
                    at_num: 2,
                    parity: 1,
                },
                AT_STEREO_CARB {
                    at_num: 7,
                    parity: 2,
                },
            ])
            .unwrap();
        let inverted_centers = heap
            .allocate_model_storage(vec![
                AT_STEREO_CARB {
                    at_num: 2,
                    parity: 2,
                },
                AT_STEREO_CARB {
                    at_num: 7,
                    parity: 1,
                },
            ])
            .unwrap();
        assert_eq!(
            CopyLinearCTStereoToINChIStereo(
                &mut heap,
                Some(&mut stereo),
                direct_centers.as_const(),
                2,
                null_dble,
                0,
                null_rank,
                null_rank,
                null_atom,
                0,
                inverted_centers.as_const(),
                null_dble,
                null_rank,
                null_rank,
            ),
            Ok(0)
        );
        assert_eq!(stereo.nNumberOfStereoCenters, 2);
        assert_eq!(stereo.nCompInv2Abs, 1);
        assert_eq!(stereo.bTrivialInv, 1);

        let number_difference = heap
            .allocate_model_storage(vec![
                AT_STEREO_CARB {
                    at_num: 3,
                    parity: 1,
                },
                AT_STEREO_CARB {
                    at_num: 7,
                    parity: 2,
                },
            ])
            .unwrap();
        assert_eq!(
            CopyLinearCTStereoToINChIStereo(
                &mut heap,
                Some(&mut stereo),
                direct_centers.as_const(),
                2,
                null_dble,
                0,
                null_rank,
                null_rank,
                null_atom,
                0,
                number_difference.as_const(),
                null_dble,
                null_rank,
                null_rank,
            ),
            Ok(0)
        );
        assert_eq!(stereo.nCompInv2Abs, 1);
        assert_eq!(stereo.bTrivialInv, 0);

        let bond = heap
            .allocate_model_storage(vec![AT_STEREO_DBLE {
                at_num1: 4,
                at_num2: 9,
                parity: 3,
            }])
            .unwrap();
        assert_eq!(
            CopyLinearCTStereoToINChIStereo(
                &mut heap,
                Some(&mut stereo),
                null_carb,
                0,
                bond.as_const(),
                1,
                null_rank,
                null_rank,
                null_atom,
                0,
                null_carb,
                bond.as_const(),
                null_rank,
                null_rank,
            ),
            Ok(0)
        );
        assert_eq!(stereo.nNumberOfStereoBonds, 1);
        assert_eq!(heap.slice(stereo.nBondAtom1.as_const()).unwrap()[0], 4);
        assert_eq!(heap.slice(stereo.nBondAtom2.as_const()).unwrap()[0], 9);
        assert_eq!(heap.slice(stereo.b_parity.as_const()).unwrap()[0], 3);

        let different_bond = heap
            .allocate_model_storage(vec![AT_STEREO_DBLE {
                at_num1: 4,
                at_num2: 10,
                parity: 3,
            }])
            .unwrap();
        stereo.nNumberOfStereoBonds = 77;
        assert_eq!(
            CopyLinearCTStereoToINChIStereo(
                &mut heap,
                Some(&mut stereo),
                null_carb,
                0,
                bond.as_const(),
                1,
                null_rank,
                null_rank,
                null_atom,
                0,
                null_carb,
                different_bond.as_const(),
                null_rank,
                null_rank,
            ),
            Ok(-4)
        );
        assert_eq!(stereo.nNumberOfStereoBonds, 77);

        assert_eq!(
            CopyLinearCTStereoToINChIStereo(
                &mut heap,
                Some(&mut stereo),
                null_carb,
                -2,
                null_dble,
                -1,
                null_rank,
                null_rank,
                null_atom,
                0,
                null_carb,
                null_dble,
                null_rank,
                null_rank,
            ),
            Ok(0)
        );
        assert_eq!(stereo.nNumberOfStereoCenters, -2);
        assert_eq!(stereo.nNumberOfStereoBonds, 0);
        assert_eq!(stereo.nCompInv2Abs, 0);
        assert_eq!(stereo.bTrivialInv, 0);
    }

    #[test]
    fn source_port__ichimak2__copy2stereobondorallene__line_451() {
        fn stereo_fixture(heap: &mut SourceHeap) -> INChI_Stereo {
            INChI_Stereo {
                nNumberOfStereoCenters: 2,
                nNumber: heap
                    .allocate_model_storage(vec![2_u16, 8, 0, 0, 0, 0])
                    .unwrap(),
                t_parity: heap
                    .allocate_model_storage(vec![11_i8, 22, 0, 0, 0, 0])
                    .unwrap(),
                nNumberInv: heap
                    .allocate_model_storage(vec![20_u16, 50, 80, 0, 0, 0])
                    .unwrap(),
                t_parityInv: heap
                    .allocate_model_storage(vec![-11_i8, -22, -33, 0, 0, 0])
                    .unwrap(),
                nNumberOfStereoBonds: 0,
                nBondAtom1: heap
                    .allocate_model_storage(vec![91_u16, 92, 93, 94, 95, 96])
                    .unwrap(),
                nBondAtom2: heap
                    .allocate_model_storage(vec![81_u16, 82, 83, 84, 85, 86])
                    .unwrap(),
                b_parity: heap
                    .allocate_model_storage(vec![71_i8, 72, 73, 74, 75, 76])
                    .unwrap(),
                ..INChI_Stereo::default()
            }
        }

        let mut heap = SourceHeap::default();
        let stereo = stereo_fixture(&mut heap);
        let direct_bond = AT_STEREO_DBLE {
            at_num1: 3,
            at_num2: 9,
            parity: u8::MAX,
        };
        let mut centers = 2;
        let mut bonds = 1;
        assert_eq!(
            Copy2StereoBondOrAllene(
                &mut heap,
                &stereo,
                &mut centers,
                Some(&mut bonds),
                &direct_bond,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(centers, 2);
        assert_eq!(bonds, 2);
        assert_eq!(heap.slice(stereo.b_parity.as_const()).unwrap()[1], -1);
        assert_eq!(heap.slice(stereo.nBondAtom1.as_const()).unwrap()[1], 3);
        assert_eq!(heap.slice(stereo.nBondAtom2.as_const()).unwrap()[1], 9);

        let canonical_order_only = heap.allocate_model_storage(vec![0_u16]).unwrap();
        assert_eq!(
            Copy2StereoBondOrAllene(
                &mut heap,
                &stereo,
                &mut centers,
                Some(&mut bonds),
                &AT_STEREO_DBLE {
                    at_num1: 7,
                    at_num2: 6,
                    parity: 4,
                },
                canonical_order_only.as_const(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(bonds, 3);
        assert_eq!(heap.slice(stereo.b_parity.as_const()).unwrap()[2], 4);
        let bonds_before = heap.slice(stereo.nBondAtom1.as_const()).unwrap().to_vec();
        assert_eq!(
            Copy2StereoBondOrAllene(
                &mut heap,
                &stereo,
                &mut centers,
                None,
                &direct_bond,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(
            heap.slice(stereo.nBondAtom1.as_const()).unwrap(),
            bonds_before
        );

        let atoms = heap
            .allocate_model_storage(vec![
                sp_ATOM {
                    neighbor: {
                        let mut values = [0; 20];
                        values[0] = 1;
                        values
                    },
                    stereo_bond_parity: [(MULT_STEREOBOND + 1) as i8, 0, 0],
                    stereo_bond_ord: [0, 0, 0],
                    ..sp_ATOM::default()
                },
                sp_ATOM::default(),
                sp_ATOM::default(),
            ])
            .unwrap();
        let canonical_order = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let canonical_rank = heap.allocate_model_storage(vec![10_u16, 5, 4]).unwrap();
        let mut allene_bonds = 0;
        assert_eq!(
            Copy2StereoBondOrAllene(
                &mut heap,
                &stereo,
                &mut centers,
                Some(&mut allene_bonds),
                &AT_STEREO_DBLE {
                    at_num1: 1,
                    at_num2: 3,
                    parity: u8::MAX,
                },
                canonical_order.as_const(),
                canonical_rank.as_const(),
                atoms.as_const(),
                0,
            ),
            Ok(1)
        );
        assert_eq!(centers, 3);
        assert_eq!(allene_bonds, 0);
        assert_eq!(
            &heap.slice(stereo.nNumber.as_const()).unwrap()[..4],
            &[2, 5, 8, 0]
        );
        assert_eq!(
            &heap.slice(stereo.t_parity.as_const()).unwrap()[..4],
            &[11, -1, 22, 0]
        );

        {
            let atoms = heap.slice_mut(atoms).unwrap();
            atoms[0].stereo_bond_parity2[0] = (3 * MULT_STEREOBOND + 2) as i8;
            atoms[0].stereo_bond_ord2[0] = 0;
            atoms[1].valence = 2;
            atoms[1].neighbor[0] = 0;
            atoms[1].neighbor[1] = 2;
        }
        let primary_before = heap.slice(stereo.nNumber.as_const()).unwrap().to_vec();
        let mut inverted_centers = 3;
        assert_eq!(
            Copy2StereoBondOrAllene(
                &mut heap,
                &stereo,
                &mut inverted_centers,
                None,
                &AT_STEREO_DBLE {
                    at_num1: 1,
                    at_num2: 3,
                    parity: 7,
                },
                canonical_order.as_const(),
                canonical_rank.as_const(),
                atoms.as_const(),
                -1,
            ),
            Ok(1)
        );
        assert_eq!(inverted_centers, 4);
        assert_eq!(
            heap.slice(stereo.nNumber.as_const()).unwrap(),
            primary_before
        );
        assert_eq!(
            &heap.slice(stereo.nNumberInv.as_const()).unwrap()[..5],
            &[20, 4, 50, 80, 0]
        );
        assert_eq!(
            &heap.slice(stereo.t_parityInv.as_const()).unwrap()[..5],
            &[-11, 7, -22, -33, 0]
        );

        let mut rejected_bonds = 3;
        for (encoded_parity, second_neighbor, middle_valence) in [
            ((2 * MULT_STEREOBOND + 1) as i8, 0_u16, 2_i8),
            ((MULT_STEREOBOND + 1) as i8, 9_u16, 2_i8),
            ((3 * MULT_STEREOBOND + 1) as i8, 0_u16, 1_i8),
        ] {
            {
                let atoms = heap.slice_mut(atoms).unwrap();
                atoms[0].stereo_bond_parity[0] = encoded_parity;
                atoms[0].stereo_bond_neighbor[1] = second_neighbor;
                atoms[1].valence = middle_valence;
            }
            assert_eq!(
                Copy2StereoBondOrAllene(
                    &mut heap,
                    &stereo,
                    &mut centers,
                    Some(&mut rejected_bonds),
                    &AT_STEREO_DBLE {
                        at_num1: 1,
                        at_num2: 2,
                        parity: 6,
                    },
                    canonical_order.as_const(),
                    canonical_rank.as_const(),
                    atoms.as_const(),
                    0,
                ),
                Ok(0)
            );
        }
        assert_eq!(rejected_bonds, 6);
        assert_eq!(
            &heap.slice(stereo.b_parity.as_const()).unwrap()[3..6],
            &[6, 6, 6]
        );
        assert_eq!(
            &heap.slice(stereo.nBondAtom1.as_const()).unwrap()[3..6],
            &[1, 1, 1]
        );
        assert_eq!(
            &heap.slice(stereo.nBondAtom2.as_const()).unwrap()[3..6],
            &[2, 2, 2]
        );
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
