use crate::source::base::ichiring::is_bond_in_Nmax_memb_ring;
use crate::source::base::ichirvr1::{
    AddToEdgeList, AllocEdgeList, FindInEdgeList, GetChargeFlowerUpperEdge,
    MakeOneInChIOutOfStrFromINChI2, RemoveForbiddenEdgeMask, RemoveFromEdgeListByValue,
    RunBnsRestoreOnce, RunBnsTestOnce, SetForbiddenEdgeMask,
};
use crate::source::base::ichirvr4::{FillOutCMP2FHINCHI, FillOutExtraFixedHDataRestr};
use crate::source_types::{
    ALL_TC_GROUPS, AT_NUMB, BN_DATA, BN_STRUCT, BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, CANON_GLOBALS,
    CMP2FHINCHI, EDGE_LIST, EDGE_LIST_CLEAR, EDGE_LIST_FREE, INCHI_CLOCK, INChI, INPUT_PARMS,
    MAX_DIFF_FIXH, NO_VERTEX, STRUCT_DATA, SourceHeap, SourceHeapError, SourceMutPointer,
    SourceTGroupInfoPointer, StrFromINChI, T_GROUP_INFO, TAUT_YES, VAL_AT, clock_t, inp_ATOM,
    local_ichirvr3::{fNumRNeutrlH, fNumRPosChgH},
};

const INC_ADD_EDGE: i32 = 64;
const EL_NUMBER_N: u8 = 7;

const F_NUM_R_POS_CHG_H: usize = 0;
const F_NUM_R_POS_CHG_U: usize = 1;
const F_NUM_R_NEG_CHG_O: usize = 2;
const F_NUM_R_NEG_CHG_N: usize = 3;
const F_NUM_R_NEUTRL_H: usize = 4;
const F_NUM_N_POS_CHG_H: usize = 5;
const F_NUM_N_POS_CHG_U: usize = 6;
const F_NUM_N_NEG_CHG_O: usize = 7;
const F_NUM_N_NEG_CHG_N: usize = 8;
const F_NUM_N_NEUTRL_H: usize = 9;
const F_NUM_ALL_CHG_T: usize = 10;

// BEGIN INCHI C TYPE: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:57 TgDiffHChgFH
// INCHI✔️✔️: #define fNumRPosChgH 0 /* number of positive charges on endpoints that have H in at2[] */
// INCHI✔️✔️: #define fNumRPosChgU 1 /* number of positive charges on endpoints that have no H in at2[] */
// INCHI✔️✔️: #define fNumRNegChgO 2 /* number of negative charges on O endpoints */
// INCHI✔️✔️: #define fNumRNegChgN 3 /* number of negative charges on N endpoints */
// INCHI✔️✔️: #define fNumRNeutrlH 4 /* number of neutral endp that have H in at2[] */
// INCHI✔️✔️:
// INCHI✔️✔️: #define fNumNPosChgH 5 /* number of positive charges on endpoints that have H in atf[] */
// INCHI✔️✔️: #define fNumNPosChgU 6 /* number of positive charges on endpoints that have no H in atf[] */
// INCHI✔️✔️: #define fNumNNegChgO 7 /* number of negative charges on O endpoints */
// INCHI✔️✔️: #define fNumNNegChgN 8 /* number of negative charges on N endpoints */
// INCHI✔️✔️: #define fNumNNeutrlH 9 /* number of neutral endp that have H in atf[] */
// INCHI✔️✔️:
// INCHI✔️✔️: #define fNumAllChgT 10 /* total  number of fNum... */
// INCHI✔️✔️:
// INCHI✔️✔️: typedef struct tagTgDiffHChgFH {
// INCHI✔️✔️:     short  itg; /* t-group index; endpoint = itg+1 */
// INCHI✔️✔️:     short  nNumHInchi;  /* number of H in t-group from orig. InChI */
// INCHI✔️✔️:     short  nNumHRevrs;  /* number of H in at2[] */
// INCHI✔️✔️:     short  nNumHNorml;  /* number of H in Normalized atfMobile_H_Revrs[] */
// INCHI✔️✔️:     short  nNumMInchi;  /* number of (-) in InChI */
// INCHI✔️✔️:     short  nNumMRevrs;  /* number of (-) in at2[] */
// INCHI✔️✔️:     short  nNumMNorml;  /* number of (-) in atf[] */
// INCHI✔️✔️:     short  nNumPRevrs;  /* number of (+) in at2[] */
// INCHI✔️✔️:     short  nNumPNorml;  /* number of (+) in Normalized atfMobile_H_Revrs[] */
// INCHI✔️✔️:     short n[fNumAllChgT]; /* all numbers */
// INCHI✔️✔️:     short i[fNumAllChgT]; /* all indices */
// INCHI✔️✔️: } TgDiffHChgFH;
// END INCHI C TYPE: TgDiffHChgFH
#[repr(C)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct TgDiffHChgFH {
    pub(crate) itg: i16,
    pub(crate) nNumHInchi: i16,
    pub(crate) nNumHRevrs: i16,
    pub(crate) nNumHNorml: i16,
    pub(crate) nNumMInchi: i16,
    pub(crate) nNumMRevrs: i16,
    pub(crate) nNumMNorml: i16,
    pub(crate) nNumPRevrs: i16,
    pub(crate) nNumPNorml: i16,
    pub(crate) n: [i16; F_NUM_ALL_CHG_T],
    pub(crate) i: [i16; F_NUM_ALL_CHG_T],
}

impl TgDiffHChgFH {
    fn word(&self, index: usize) -> i16 {
        match index {
            0 => self.itg,
            1 => self.nNumHInchi,
            2 => self.nNumHRevrs,
            3 => self.nNumHNorml,
            4 => self.nNumMInchi,
            5 => self.nNumMRevrs,
            6 => self.nNumMNorml,
            7 => self.nNumPRevrs,
            8 => self.nNumPNorml,
            9..=18 => self.n[index - 9],
            19..=28 => self.i[index - 19],
            _ => unreachable!(),
        }
    }

    fn set_word(&mut self, index: usize, value: i16) {
        match index {
            0 => self.itg = value,
            1 => self.nNumHInchi = value,
            2 => self.nNumHRevrs = value,
            3 => self.nNumHNorml = value,
            4 => self.nNumMInchi = value,
            5 => self.nNumMRevrs = value,
            6 => self.nNumMNorml = value,
            7 => self.nNumPRevrs = value,
            8 => self.nNumPNorml = value,
            9..=18 => self.n[index - 9] = value,
            19..=28 => self.i[index - 19] = value,
            _ => unreachable!(),
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn bHas_N_V(at2: &[inp_ATOM], num_atoms: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:92 bHas_N_V
    // INCHI✔️✔️: int bHas_N_V( inp_ATOM *at2, int num_atoms )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, num_found = 0;
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (at2[i].el_number == EL_NUMBER_N &&
    // INCHI✔️✔️:              !at2[i].charge &&
    // INCHI✔️✔️:              !at2[i].num_H &&
    // INCHI✔️✔️:              !at2[i].radical &&
    // INCHI✔️✔️:              at2[i].chem_bonds_valence == 5 &&
    // INCHI✔️✔️:              ( at2[i].valence == 3 ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             num_found++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_found;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: bHas_N_V
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: bHas_N_V
    // INCHI✔️✔️: READ_INCHI_STRING=1; COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // INCHI✔️✔️: EL_NUMBER_N=((U_CHAR)7); el_number is unsigned 8-bit and the tested
    // INCHI✔️✔️: charge, num_H, radical, valence fields are signed 8-bit values.
    // END INCHI ACTIVE MACRO CONFIGURATION: bHas_N_V

    if num_atoms <= 0 {
        return Ok(0);
    }
    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = at2
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut number_found = 0_i32;
    for atom in atoms {
        if atom.el_number == EL_NUMBER_N
            && atom.charge == 0
            && atom.num_H == 0
            && atom.radical == 0
            && atom.chem_bonds_valence == 5
            && atom.valence == 3
        {
            number_found = number_found.wrapping_add(1);
        }
    }
    Ok(number_found)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillTgDiffHChgFH(
    heap: &mut SourceHeap,
    tdhc: &mut [TgDiffHChgFH],
    max_tdhc: i32,
    at2: &[inp_ATOM],
    atf: &[inp_ATOM],
    nCanon2AtnoRevrs: &[AT_NUMB],
    pVA: &[VAL_AT],
    ti: &T_GROUP_INFO,
    pAtomIndList: &mut EDGE_LIST,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:113 FillTgDiffHChgFH
    // INCHI✔️❌: complete active source frame follows verbatim; SourceHeap pointer checks and
    // INCHI✔️❌: temporary value copies add overhead compared with direct C pointer traversal.
    /*
    int FillTgDiffHChgFH( TgDiffHChgFH tdhc[],
                          int max_tdhc,
                          inp_ATOM at2[],
                          inp_ATOM atf[],
                          AT_NUMB  *nCanon2AtnoRevrs,
                          VAL_AT *pVA,
                          T_GROUP_INFO *ti,
                          EDGE_LIST *pAtomIndList )
    {

        int i, j, iat, itg, num, itg_out; /* djb-rwth: removing redundant variables */
        EDGE_LIST IndList;   /* type, itg */
        TgDiffHChgFH cur_tdhc;
        AT_NUMB    *pEndp0;
        inp_ATOM *at2i, *atfi;
        int         typeR, typeN, type, ret = 0, nCurIndListLen;

        AllocEdgeList( &IndList, EDGE_LIST_CLEAR );
        pAtomIndList->num_edges = 0;
        itg_out = 0;
        /* djb-rwth: removing redundant code */
        memset( tdhc, 0, max_tdhc ); /* djb-rwth: memset_s C11/Annex K variant? */

        for (itg = 0; itg < ti->num_t_groups; itg++)
        {
            memset( &cur_tdhc, 0, sizeof( cur_tdhc ) ); /* djb-rwth: memset_s C11/Annex K variant? */

            cur_tdhc.itg = itg;
            cur_tdhc.nNumHInchi = ti->t_group[itg].num[0] - ti->t_group[itg].num[1];
            cur_tdhc.nNumMInchi = ti->t_group[itg].num[1];

            pEndp0 = ti->nEndpointAtomNumber + ti->t_group[itg].nFirstEndpointAtNoPos;
            nCurIndListLen = IndList.num_edges;
            for (j = 0; j < ti->t_group[itg].nNumEndpoints; j++)
            {
                i = pEndp0[j];
                iat = nCanon2AtnoRevrs[i];

                at2i = at2 + iat;
                atfi = atf + iat;

                typeR = typeN = -1;
                if (at2i->charge == 1)
                {
                    if (at2i->num_H)
                    {
                        typeR = fNumRPosChgH;
                    }
                    else
                    {
                        typeR = fNumRPosChgU;
                    }
                    cur_tdhc.nNumPRevrs++;
                }
                else
                {
                    if (at2i->charge == -1)
                    {
                        if (pVA[iat].cNumValenceElectrons == 6)
                        {
                            typeR = fNumRNegChgO;
                        }
                        else
                        {
                            if (pVA[iat].cNumValenceElectrons == 5)
                            {
                                typeR = fNumRNegChgN;
                            }
                        }
                        cur_tdhc.nNumMRevrs++;
                    }
                    else
                    {
                        if (at2i->num_H && at2i->valence == at2i->chem_bonds_valence)
                        {
                            typeR = fNumRNeutrlH;
                        }
                    }
                }
                cur_tdhc.nNumHRevrs += at2i->num_H;

                if (atfi->charge == 1)
                {
                    if (atfi->num_H)
                    {
                        typeN = fNumNPosChgH;
                    }
                    else
                    {
                        typeN = fNumNPosChgU;
                    }
                    cur_tdhc.nNumPNorml++;
                }
                else
                {
                    if (atfi->charge == -1)
                    {
                        if (pVA[iat].cNumValenceElectrons == 6)
                        {
                            typeN = fNumNNegChgO;
                        }
                        else
                            if (pVA[iat].cNumValenceElectrons == 5)
                            {
                                typeN = fNumNNegChgN;
                            }
                        cur_tdhc.nNumMNorml++;
                    }
                    else
                    {
                        if (atfi->num_H && atfi->valence == atfi->chem_bonds_valence)
                        {
                            typeN = fNumNNeutrlH;
                        }
                    }
                }
                cur_tdhc.nNumHNorml += atfi->num_H;
                if (at2[iat].charge < 0 || 0 < pVA[iat].nCPlusGroupEdge)
                {
                    if (typeR >= 0 && (
                        ( ret = AddToEdgeList( &IndList, typeR, INC_ADD_EDGE ) ) ||
                        ( ret = AddToEdgeList( &IndList, itg, INC_ADD_EDGE ) ) ||
                        ( ret = AddToEdgeList( &IndList, iat, INC_ADD_EDGE ) ) ))
                    {
                        goto exit_function;
                    }
                    if (typeN >= 0 && (
                        ( ret = AddToEdgeList( &IndList, typeN, INC_ADD_EDGE ) ) ||
                        ( ret = AddToEdgeList( &IndList, itg, INC_ADD_EDGE ) ) ||
                        ( ret = AddToEdgeList( &IndList, iat, INC_ADD_EDGE ) ) ))
                    {
                        goto exit_function;
                    }
                }
            }
            if (cur_tdhc.nNumHNorml == cur_tdhc.nNumHInchi &&
                 cur_tdhc.nNumMNorml == cur_tdhc.nNumMInchi)
            {
                IndList.num_edges = nCurIndListLen; /* t-group seems to be correct */
                continue;
            }
            if (itg_out < max_tdhc)
            {
                tdhc[itg_out++] = cur_tdhc;
            }
            else
            {
                /* djb-rwth: removing redundant code */
                IndList.num_edges = nCurIndListLen;
                break;
            }
        }

        /* fill out atom index list */
        if (itg_out)
        {
            itg_out--; /* djb-rwth: fixing undefined index value / buffer overflow */
            /* djb-rwth: removing redundant code */
            for (type = 0; type < fNumAllChgT; type++)
            {
                j = 0;
                for (i = 0; i < itg_out; i++)
                {
                    num = 0;
                    itg = tdhc[i].itg;
                    tdhc[i].i[type] = -999; /* empty */
                    while (IndList.pnEdges[j + 1] == itg)
                    {
                        if (IndList.pnEdges[j] == type)
                        {
                            if (!num++)
                            {
                                tdhc[i].i[type] = pAtomIndList->num_edges;
                            }
                            if ((ret = AddToEdgeList( pAtomIndList, IndList.pnEdges[j + 2], INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                        j += 3;
                    }
                    tdhc[i].n[type] = num;
                }
            }
        }
        ret = itg_out;

    exit_function:

        AllocEdgeList( &IndList, EDGE_LIST_FREE );

        return ret;
        /*
        #undef fNumRPosChgH
        #undef fNumRPosChgU
        #undef fNumRNegChgO
        #undef fNumRNegChgN

        #undef fNumNPosChgH
        #undef fNumNPosChgU
        #undef fNumNNegChgO
        #undef fNumNNegChgN

        #undef fNumAllChgT
        */
    }
    */
    // END INCHI C FUNCTION: FillTgDiffHChgFH
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillTgDiffHChgFH
    // INCHI✔️❌: READ_INCHI_STRING=1; COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // INCHI✔️❌: short is signed 16-bit, AT_NUMB is unsigned 16-bit, and EdgeIndex is int.
    // INCHI✔️❌: The active memset clears exactly max_tdhc bytes, not max_tdhc structures.
    // END INCHI ACTIVE MACRO CONFIGURATION: FillTgDiffHChgFH

    let max_tdhc_usize =
        usize::try_from(max_tdhc).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if max_tdhc_usize > tdhc.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let mut ind_list = EDGE_LIST::default();
    AllocEdgeList(heap, &mut ind_list, EDGE_LIST_CLEAR)?;
    pAtomIndList.num_edges = 0;

    const WORDS_PER_RECORD: usize = 29;
    const BYTES_PER_WORD: usize = std::mem::size_of::<i16>();
    for byte_index in 0..max_tdhc_usize {
        let record_index = byte_index / (WORDS_PER_RECORD * BYTES_PER_WORD);
        let byte_in_record = byte_index % (WORDS_PER_RECORD * BYTES_PER_WORD);
        let word_index = byte_in_record / BYTES_PER_WORD;
        let byte_in_word = byte_in_record % BYTES_PER_WORD;
        let record = tdhc
            .get_mut(record_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut bytes = record.word(word_index).to_ne_bytes();
        bytes[byte_in_word] = 0;
        record.set_word(word_index, i16::from_ne_bytes(bytes));
    }

    let result = (|| -> Result<i32, SourceHeapError> {
        let mut itg_out = 0_i32;
        let group_count =
            usize::try_from(ti.num_t_groups).map_err(|_| SourceHeapError::PointerOutOfBounds)?;

        for group_index in 0..group_count {
            let group = heap
                .slice(ti.t_group.as_const())?
                .get(group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let mut current = TgDiffHChgFH {
                itg: group_index as i16,
                nNumHInchi: i32::from(group.num[0]).wrapping_sub(i32::from(group.num[1])) as i16,
                nNumMInchi: group.num[1] as i16,
                ..TgDiffHChgFH::default()
            };
            let first_endpoint = usize::from(group.nFirstEndpointAtNoPos);
            let endpoint_count = usize::from(group.nNumEndpoints);
            let current_ind_list_len = ind_list.num_edges;

            for endpoint_offset in 0..endpoint_count {
                let endpoint_position = first_endpoint
                    .checked_add(endpoint_offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let canonical_index = usize::from(
                    *heap
                        .slice(ti.nEndpointAtomNumber.as_const())?
                        .get(endpoint_position)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::from(
                    *nCanon2AtnoRevrs
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let reverse_atom = at2
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let normalized_atom = atf
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence_atom = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;

                let mut reverse_type = -1_i32;
                let mut normalized_type = -1_i32;
                if reverse_atom.charge == 1 {
                    reverse_type = if reverse_atom.num_H != 0 {
                        F_NUM_R_POS_CHG_H as i32
                    } else {
                        F_NUM_R_POS_CHG_U as i32
                    };
                    current.nNumPRevrs = current.nNumPRevrs.wrapping_add(1);
                } else if reverse_atom.charge == -1 {
                    if valence_atom.cNumValenceElectrons == 6 {
                        reverse_type = F_NUM_R_NEG_CHG_O as i32;
                    } else if valence_atom.cNumValenceElectrons == 5 {
                        reverse_type = F_NUM_R_NEG_CHG_N as i32;
                    }
                    current.nNumMRevrs = current.nNumMRevrs.wrapping_add(1);
                } else if reverse_atom.num_H != 0
                    && reverse_atom.valence == reverse_atom.chem_bonds_valence
                {
                    reverse_type = F_NUM_R_NEUTRL_H as i32;
                }
                current.nNumHRevrs = current
                    .nNumHRevrs
                    .wrapping_add(i16::from(reverse_atom.num_H));

                if normalized_atom.charge == 1 {
                    normalized_type = if normalized_atom.num_H != 0 {
                        F_NUM_N_POS_CHG_H as i32
                    } else {
                        F_NUM_N_POS_CHG_U as i32
                    };
                    current.nNumPNorml = current.nNumPNorml.wrapping_add(1);
                } else if normalized_atom.charge == -1 {
                    if valence_atom.cNumValenceElectrons == 6 {
                        normalized_type = F_NUM_N_NEG_CHG_O as i32;
                    } else if valence_atom.cNumValenceElectrons == 5 {
                        normalized_type = F_NUM_N_NEG_CHG_N as i32;
                    }
                    current.nNumMNorml = current.nNumMNorml.wrapping_add(1);
                } else if normalized_atom.num_H != 0
                    && normalized_atom.valence == normalized_atom.chem_bonds_valence
                {
                    normalized_type = F_NUM_N_NEUTRL_H as i32;
                }
                current.nNumHNorml = current
                    .nNumHNorml
                    .wrapping_add(i16::from(normalized_atom.num_H));

                if reverse_atom.charge < 0 || valence_atom.nCPlusGroupEdge > 0 {
                    if reverse_type >= 0 {
                        for value in [reverse_type, group_index as i32, atom_index as i32] {
                            let ret = AddToEdgeList(heap, &mut ind_list, value, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(ret);
                            }
                        }
                    }
                    if normalized_type >= 0 {
                        for value in [normalized_type, group_index as i32, atom_index as i32] {
                            let ret = AddToEdgeList(heap, &mut ind_list, value, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(ret);
                            }
                        }
                    }
                }
            }

            if current.nNumHNorml == current.nNumHInchi && current.nNumMNorml == current.nNumMInchi
            {
                ind_list.num_edges = current_ind_list_len;
                continue;
            }
            if itg_out < max_tdhc {
                *tdhc
                    .get_mut(
                        usize::try_from(itg_out)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = current;
                itg_out = itg_out.wrapping_add(1);
            } else {
                ind_list.num_edges = current_ind_list_len;
                break;
            }
        }

        if itg_out != 0 {
            itg_out = itg_out.wrapping_sub(1);
            for type_index in 0..F_NUM_ALL_CHG_T {
                let mut list_index = 0_i32;
                for output_index in 0..itg_out {
                    let output_index_usize = usize::try_from(output_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let group_index = i32::from(
                        tdhc.get(output_index_usize)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .itg,
                    );
                    tdhc[output_index_usize].i[type_index] = -999;
                    let mut number = 0_i32;
                    loop {
                        let list_index_usize = usize::try_from(list_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let values = heap.slice(ind_list.pnEdges.as_const())?;
                        let record_group = *values
                            .get(list_index_usize + 1)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if record_group != group_index {
                            break;
                        }
                        let record_type = *values
                            .get(list_index_usize)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if record_type == type_index as i32 {
                            if number == 0 {
                                tdhc[output_index_usize].i[type_index] =
                                    pAtomIndList.num_edges as i16;
                            }
                            number = number.wrapping_add(1);
                            let atom_index = *values
                                .get(list_index_usize + 2)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let ret = AddToEdgeList(heap, pAtomIndList, atom_index, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(ret);
                            }
                        }
                        list_index = list_index.wrapping_add(3);
                    }
                    tdhc[output_index_usize].n[type_index] = number as i16;
                }
            }
        }
        Ok(itg_out)
    })();

    let cleanup = AllocEdgeList(heap, &mut ind_list, EDGE_LIST_FREE);
    match (result, cleanup) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(_)) => Ok(value),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixFixedHRestoredStructure(
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
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:333 FixFixedHRestoredStructure
    // INCHI✔️❌: complete active source frame follows verbatim; source-ordered behavior is
    // INCHI✔️❌: reproduced, while SourceHeap pointer and ownership checks add known overhead.
    /*
    int FixFixedHRestoredStructure( CANON_GLOBALS *pCG,
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
        int i, j, k, delta, num_try, tot_succes, cur_success, ret = 0, bAllowedNFlowerEdges = 0, num_zero_ret;
        CMP2FHINCHI c2i;
        CMP2FHINCHI *pc2i = &c2i;

        EDGE_LIST AllChargeEdges, CurrEdges, SFlowerEdges, NFlowerEdges, OtherNFlowerEdges, FixedLargeRingStereoEdges;
        EDGE_LIST AllBondEdges;

        EdgeIndex e;
        BNS_EDGE  *pe;
        Vertex v1, v2;
        BNS_VERTEX *pv1, *pv2;

        Vertex     vPathStart, vPathEnd;
        int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

        int        nNumRunBNS = 0, forbidden_edge_mask_inv = ~forbidden_edge_mask; /* djb-rwth: ignoring LLVM warning: variable used */

        INCHI_HEAPCHK

        AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &CurrEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &NFlowerEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &SFlowerEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &AllBondEdges, EDGE_LIST_CLEAR );

        tot_succes = 0;

        if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
        {
            goto exit_function;  /* no fixed-H found */
        }

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

                    if (pBNS->edge[j].forbidden)
                    {
                        continue;
                    }

                    if (pBNS->edge[j].flow)
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
                else
                {
                    /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                    if (pVA[i].cNumValenceElectrons == 6 && !pVA[i].cMetal && /* N, P, As */
                         NO_VERTEX != ( j = GetChargeFlowerUpperEdge( pBNS, pVA, e ) ))
                    {

                        if (pBNS->edge[j].forbidden)
                        {
                            continue;
                        }

                        if (pBNS->edge[j].flow)
                        {
                            if ((ret = AddToEdgeList( &SFlowerEdges, j, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                    }
                }
            }
            for (j = 0; j < at2[i].valence; j++)
            {
                k = at2[i].neighbor[j];
                if (k < i && !pBNS->edge[e = pBNS->vert[i].iedge[j]].forbidden)
                {
                    if ((ret = AddToEdgeList( &AllBondEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
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
            if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }

        INCHI_HEAPCHK

            if (!pc2i->bHasDifference ||
                 (!pc2i->len_c2at && pc2i->nNumTgRevrs == pc2i->nNumTgInChI &&
                 pc2i->nNumEndpRevrs == pc2i->nNumRemHInChI &&
                 pc2i->nNumEndpRevrs == pc2i->nNumEndpInChI &&
                 !pc2i->nNumTgDiffMinus && !pc2i->nNumTgDiffH)) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function; /* nothing to do */
            }

        /*goto exit_function;*/ /* debug only*/

        if (pc2i->len_c2at >= 2)
        {
            /*----------------------------------------------------*/
            /* case 01: restored: O=AB-O(-)  original:  (-)O-AB=O */
            /* FixH:              0    -1                 -1    0 */
            /* MobH:              0     1                  1    0 */
            /*                         non-taut      non-taut     */
            /* O = O, S, Se; charged atoms O are not tautomeric   */
            /* Solution: move (-) from B-O(-) to O=A              */
            /*----------------------------------------------------*/
            int num_DB_O = 0, num_SB_O_Minus = 0, iat;
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if (pc2i->c2at[i].nValElectr == 6 /* && !pc2i->c2at[i].endptInChI -- mod#1*/ &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                {
                    if ( /* orig. InChI info: */
                        num_SB_O_Minus < MAX_DIFF_FIXH &&
                        pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                        /* reversed structure info: */
                        pc2i->c2at[i].nFixHRevrs == -1 && pc2i->c2at[i].nMobHRevrs == 1 &&
                        pc2i->c2at[i].nAtChargeRevrs == -1 && !at2[iat].num_H && /* at2 is Fixed-H */
                        at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 1)
                    {
                        iat_SB_O_Minus[num_SB_O_Minus++] = iat;
                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                    else
                    {
                        if ( /* orig. InChI info: */
                            num_DB_O < MAX_DIFF_FIXH &&
                            pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI == 1 &&
                            /* reversed structure info: */
                            pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                            pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                            at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2)
                        {
                            iat_DB_O[num_DB_O++] = iat;
                            if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_O_Minus, num_DB_O))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_O_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_O_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta; /* remove (-) from AB-O(-) */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Added (-)charge to O=AB => nDeltaCharge == -1 */
                        /* Flow change on pe (-)charge edge (atom B-O(-)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 01 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                CurrEdges.num_edges = 0; /* clear current edge list */
            }
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic,
                    ip, sd, pBNS, pStruct,
                    at, at2, at3, pVA,
                    pTCGroups,
                    ppt_group_info, ppat_norm,
                    ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 1)
        {
            /*--------------------------------------------------------------*/
            /* case 02: restored: -O(+)=AB-NH2  original:  -O-AB=NH2(+)     */
            /* FixH:               0        0               0      1        */
            /* MobH:               0        2               0      1        */
            /* O = P, As, Sb, O, S, Se, F, Cl, Br, I; not taut. in InChI    */
            /* N = N, O, S, Se, Te; has H; tautomeric or not tautomeric     */
            /* Solution: move (+) from O(+) to NH2                          */
            /*--------------------------------------------------------------*/
            int num_DB_O_Plus = 0, num_SB_NH = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : NULL;
            cur_success = 0;
            num_zero_ret = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: =NH2(+), =OH(+) */
                    num_SB_NH < MAX_DIFF_FIXH &&
                    ((pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) ||
                        pc2i->c2at[i].nValElectr == 6) /* N, O, S, Se, Te */ &&
                    /*!pc2i->c2at[i].endptInChI &&*/ /* <=== relaxation */
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI>0 /*== 1 --modification#2*/ && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                    /* reversed structure info: */
                    pc2i->c2at[i].nFixHRevrs == 0 && /* pc2i->c2at[i].nMobHRevrs == 0 &&*/
                    pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H &&
                    at2[iat].valence == at2[iat].chem_bonds_valence) /* djb-rwth: addressing LLVM warning */
                {
                    iat_SB_NH[num_SB_NH++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=+1, no H, has double bond, P, As, O, S, Se, Te, F, Cl, Br, I */
                    num_DB_O_Plus < MAX_DIFF_FIXH &&
                    at2[iat].charge == 1 && !at2[iat].num_H &&
                    at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    (pVA[iat].cNumValenceElectrons == 6 || pVA[iat].cNumValenceElectrons == 7 ||
                        (pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber > 1)) &&
                    /* in orig.InChI: not an endpoint, has no H */
                    !pStruct->endpoint[i] &&
                    !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    !(nMobHInChI && nMobHInChI[i]) &&
                    /* has (+) edge */
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                {
                    iat_DB_O_Plus[num_DB_O_Plus++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if ((num_try = inchi_min(num_DB_O_Plus, num_SB_NH))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
            repeat_02_allow_NV:
                for (i = 0; i < num_SB_NH && cur_success < num_try; i++)
                {
                    iat = iat_SB_NH[i];
                    pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Removed charge from O(+) => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 02 */
                        }
                    }
                    else
                    {
                        num_zero_ret += !ret;
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                if (num_zero_ret == num_try && !bAllowedNFlowerEdges && NFlowerEdges.num_edges)
                {
                    RemoveForbiddenEdgeMask(pBNS, &NFlowerEdges, forbidden_edge_mask);
                    bAllowedNFlowerEdges = 1;
                    goto repeat_02_allow_NV;
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                /* djb-rwth: removing redundant code */
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 1 && pc2i->nNumTgRevrs == 1 &&
            ( pc2i->nNumEndpRevrs > pc2i->nNumEndpInChI || pc2i->nNumTgInChI > 1 ) /* ADP in Revrs */)
        {
            /*--------------------------------------------------------------*/
            /* case 03: restored: -N(-)-AB=O    original:  -N=AB-O(-)       */
            /* FixH:               0       0                0     -1        */
            /* MobH:               0       0                0      1        */
            /* O = O, S, Se; N = N;                                         */
            /* restored atoms are tautomeric; original atoms are not taut.  */
            /* restored struct has 1 t-group; original has less endpoints   */
            /*                                and possibly >1 t-groups      */
            /* Solution: move (-) from N(-) to =O                           */
            /*           these atoms are tautomeric in restored structure   */
            /*--------------------------------------------------------------*/
            int num_SB_N_Minus = 0, num_DB_O = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: -O(-) */
                    num_DB_O < MAX_DIFF_FIXH &&
                    pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                    !pc2i->c2at[i].endptInChI &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI == 1 &&
                    /* reversed structure info: */
                    pc2i->c2at[i].endptRevrs &&
                    pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                    pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                    at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2)
                {
                    iat_DB_O[num_DB_O++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom N: charge=-1, no H, has no double bond, endpoint */
                    num_SB_N_Minus < MAX_DIFF_FIXH &&
                    at2[iat].charge == -1 && /*!at2[iat].num_H &&*/
                    at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                    at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint &&
                    /* in orig.InChI: not an endpoint, has no H */
                    /* !pStruct->endpoint[i] && */
                    /*
                    !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    !(nMobHInChI && nMobHInChI[i] ) &&
                    */
                    /* has (-) edge */
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                {
                    iat_SB_N_Minus[num_SB_N_Minus++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_N_Minus, num_DB_O))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_N_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_N_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1; /* 2006-03-03: changed from CPlusGroupEdge */
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Added (-) charge to =O => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 03 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->nNumTgRevrs == 1 && /* pc2i->nNumRemHInChI < 0 &&*/
            ( pc2i->nNumEndpRevrs > pc2i->nNumEndpInChI || pc2i->nNumTgInChI > 1 ) /* ADP in Revrs */)
        {
            /*--------------------------------------------------------------*/
            /* case 03a:restored: -N(-)-AB=O    original:  -N=AB-O(-)       */
            /* FixH:               0       0                0      0        */
            /* MobH:               0       0                0      0        */
            /* O = O, S, Se; N = N;                              taut       */
            /* restored atoms are tautomeric; original atom is; N may be.   */
            /* restored struct has 1 t-group; original has less endpoints   */
            /*                                and possibly >1 t-groups      */
            /* Solution: move (-) from N(-) to =O                           */
            /*           these atoms are tautomeric in restored structure   */
            /*--------------------------------------------------------------*/
            int num_SB_N_Minus = 0, num_DB_O = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            S_CHAR   *pnMobHInChI = ( pInChI[1] && pInChI[1]->nNum_H ) ? pInChI[1]->nNum_H :
                ( pInChI[0] && pInChI[0]->nNum_H ) ? pInChI[0]->nNum_H : NULL;
            S_CHAR   *pnFixHInChI = pStruct->fixed_H;

            cur_success = 0;
            CurrEdges.num_edges = 0;

            for (i = 0; i < pStruct->num_atoms; i++)
            { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom N: charge=-1, no H, has no double bond, endpoint */
                    num_SB_N_Minus < MAX_DIFF_FIXH &&
                    at2[iat].charge == -1 && /*!at2[iat].num_H &&*/
                    at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                    at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint &&
                    /* in orig.InChI: may be an endpoint, has no H */
                    /* has (-) edge */
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                {
                    iat_SB_N_Minus[num_SB_N_Minus++] = iat;
                }
                else
                    if (num_DB_O < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && /*!at2[iat].num_H &&*/
                        at2[iat].valence + 1 == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        pVA[iat].cNumValenceElectrons == 6 &&
                        at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint && /* endpoint in Reconstructed */
                        (pStruct->endpoint[i] || /* endpoint or H(+) acceptor in original */
                            (pnMobHInChI && pnMobHInChI[i] == 1 && pnFixHInChI && pnFixHInChI[i] == -1)) &&
                        /* has (-) edge */
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                    {
                        iat_DB_O[num_DB_O++] = iat;
                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
            }
            if ((num_try = inchi_min(num_SB_N_Minus, num_DB_O))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                /* allow charge transfer to all found =O */
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_N_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_N_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Added (-) charge to =O => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 03a */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 1 && pc2i->nNumTgInChI == 1 && /* ADP in InChI */
            ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ))
        {
            /*--------------------------------------------------------------*/
            /* case 04: restored: OH(+)=AB-O- OH- orig.  HO-AB=O(+)- OH-    */
            /* FixH:               1       0   0          1      0   1      */
            /* MobH:               0       0   1          0      0   0      */
            /*                 non-taut.                taut        taut    */
            /*                                    ADP: one t-group or more endpoints */
            /* O(+) = N, P, As, As, O, S, Se; OH = N, O, S, Se, Te          */
            /* Solution: move (+) from O(+) to NH2                          */
            /*--------------------------------------------------------------*/
            int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0;

            for (i = 0; i < pStruct->num_atoms; i++)
            { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=+1, has H, has double bond, N, O, S, Se, Te */
                    num_DB_Charged < MAX_DIFF_FIXH &&
                    at2[iat].charge == 1 && at2[iat].num_H &&
                    at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    (pVA[iat].cNumValenceElectrons == 6 ||
                        (pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1)) &&
                    /* in orig.InChI: an endpoint, has fixed-H */
                    pStruct->endpoint[i] &&
                    (pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                    /* has (+) edge */
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                {

                    iat_DB_Charged[num_DB_Charged++] = iat;

                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
                else
                {
                    if ( /* in restored atom: charge=0, has no H, has no double bond, N, P, O, S, Se, Te */
                        num_SB_Neutr < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H &&
                        at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        (pVA[iat].cNumValenceElectrons == 6 ||
                            pVA[iat].cNumValenceElectrons == 5) &&
                        /* in orig.InChI: an endpoint, has fixed-H */
                        /* pStruct->endpoint[i] && */
                        !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                        !(nMobHInChI && nMobHInChI[i]) &&
                        /* has (+) edge */
                        (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 &&
                        0 == pBNS->edge[e].forbidden)
                    {

                        iat_SB_Neutr[num_SB_Neutr++] = iat;

                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_Neutr, num_DB_Charged))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_Neutr && cur_success < num_try; i++)
                {
                    iat = iat_SB_Neutr[i];
                    pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Removed charge from O(+) => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 04 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at > 1)
        {
            /*--------------------------------------------------------------*/
            /* case 05: restored:  O=AB-NH      original:(-)O-AB=NH(+)      */
            /* FixH:               0     0                 -1     1         */
            /* MobH:               0     1                  1     0         */
            /* O = O, S, Se; N = N, O, S, Se, Te; all atoms not tautomeric  */
            /* Solution: Separate charges                                   */
            /*--------------------------------------------------------------*/
            int num_DB_O = 0, num_SB_NH = 0, iat;
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: =NH2(+), =OH(+) */
                    num_SB_NH < MAX_DIFF_FIXH &&
                    ((pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) ||
                        pc2i->c2at[i].nValElectr == 6) /* N, O, S, Se, Te */ &&
                    !pc2i->c2at[i].endptInChI &&
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI == 1 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                    /* reversed structure info: */
                    pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs &&
                    pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H &&
                    !pc2i->c2at[i].endptRevrs &&
                    at2[iat].valence == at2[iat].chem_bonds_valence)
                {
                    iat_SB_NH[num_SB_NH++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
                else
                {
                    if ( /* orig. InChI info: -O(-) */
                        num_DB_O < MAX_DIFF_FIXH &&
                        (pc2i->c2at[i].nValElectr == 6) /* O, S, Se, Te */ &&
                        !pc2i->c2at[i].endptInChI &&
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                        pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI == 1 &&
                        /* reversed structure info: */
                        pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                        pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                        !pc2i->c2at[i].endptRevrs &&
                        at2[iat].valence + 1 == at2[iat].chem_bonds_valence)
                    {
                        iat_DB_O[num_DB_O++] = iat;
                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_DB_O, num_SB_NH))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_NH && cur_success < num_try; i++)
                {
                    iat = iat_SB_NH[i];
                    pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Added charge to =O => nDeltaCharge == 1 */
                        /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 05 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pStruct->fixed_H && pStruct->endpoint && pc2i->nChargeFixHInChI > 0 && pc2i->nChargeFixHInChI > pc2i->nChargeMobHInChI)
        {
            /*----------------------------------------------------------*/
            /* case 06c: restored -NH- or -NH(+)  orig: -NH-            */
            /*  Fixed-H            1       1             0              */
            /*  Mobile-H           0       0             1              */
            /*                     not tautomeric    not tautomeric     */
            /*           has adjacent (+)                               */
            /*           charges                                        */
            /*  Solution: move (+) charges to the -NH- unless it already*/
            /*            N = N, O, S, Se, Te                           */
            /*            has (+) charge blocked by adjacent (+)        */
            /*----------------------------------------------------------*/
            int iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            /*
            inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
            inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
            pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds?
            pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds : NULL;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : NULL;
            */
            EDGE_LIST CurChargeEdges;
            EdgeIndex e2;
            cur_success = 0;
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
            CurrEdges.num_edges = 0;
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                /* atoms -NH- from which H(+) were removed by the Normalization in orig. InChI */
                iat = pc2i->c2at[i].atomNumber;
                if (( pc2i->c2at[i].nValElectr == 6 ||
                      (pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) ) &&
                     !pc2i->c2at[i].endptInChI &&
                     ( e = pVA[iat].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                {
                    if ( /* orig. InChI info: -NH- */
                         pc2i->c2at[i].nFixHInChI == 1 && pc2i->c2at[i].nMobHInChI == 0 &&
                         /* reversed structure info: */
                         pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 1 && /* was not removed */
                         /*pc2i->c2at[i].nAtChargeRevrs == 0 &&*/ at2[iat].num_H && /* at2 is Fixed-H */
                         at2[iat].valence == at2[iat].chem_bonds_valence)
                    {
                        if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                /* find adjacent charged atoms */
                iat = nCanon2AtnoRevrs[i];
                if (pStruct->endpoint[i] || at2[iat].charge != 1 || at2[iat].radical || pVA[iat].cMetal)
                {
                    continue;
                }
                if (0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow && pVA[iat].cNumValenceElectrons >= 5)
                {
                    /* positively charged atom */
                    for (j = 0; j < at2[iat].valence; j++)
                    {
                        if (at2[k = (int) at2[iat].neighbor[j]].charge == 1 && !pVA[k].cMetal &&
                             0 <= ( e2 = pVA[k].nCPlusGroupEdge - 1 ) && !pBNS->edge[e2].forbidden && !pBNS->edge[e2].flow)
                        {
                            if (0 > FindInEdgeList( &CurrEdges, e ) &&
                                 0 > FindInEdgeList( &CurChargeEdges, e ) &&
                                 ( ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ) ))
                            {
                                goto exit_case_06c;
                            }
                            if (0 > FindInEdgeList( &CurrEdges, e2 ) &&
                                 0 > FindInEdgeList( &CurChargeEdges, e2 ) &&
                                 ( ret = AddToEdgeList( &CurChargeEdges, e2, INC_ADD_EDGE ) ))
                            {
                                goto exit_case_06c;
                            }
                        }
                    }
                }
            }
            if ((num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask );
                delta = 1;
                for (i = 0; i < CurrEdges.num_edges && cur_success < num_try; i++)
                {
                    e = CurrEdges.pnEdges[i];
                    pe = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                    pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                    pe->flow -= delta; /* add (+) to -NHm */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                      (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 06c */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            }
        exit_case_06c:
            CurrEdges.num_edges = 0; /* clear current edge list */
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
            if (ret < 0)
            {
                goto exit_function;
            }
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 >( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 2)
        {
            /*------------------------------------------------------------*/
            /* case 06d: restored: XH(+)=-AB-NH    orig.: XH-=AB=NH(+)    */
            /* FixH:                1       1 0          0 1      1       */
            /* MobH:                0    taut 1          1 taut   0       */
            /*                                                            */
            /*                                                            */
            /* N  = N, O, S, Se; atoms N are not tautomeric in orig InChI */
            /* X  = N, O, S, Se, Te, F, Cl, Br, I; atom X is non-taut     */
            /* Solution: move (+) from X  to NH                           */
            /*------------------------------------------------------------*/
            int iat;
            /*
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
            inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
            pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds?
            pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds :
            pStruct->pOne_norm_data[TAUT_NON]->at;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            EDGE_LIST CurChargeEdges;
            cur_success = 0;
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
            CurrEdges.num_edges = 0;
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                /* XH(+) */
                if ( /* reconstructed: non-taut and (+) */
                    ( pc2i->c2at[i].nMobHRevrs + 1 == pc2i->c2at[i].nFixHRevrs &&
                      pc2i->c2at[i].nFixHRevrs > 0 && !pc2i->c2at[i].endptRevrs &&
                      pc2i->c2at[i].nAtChargeRevrs == 1 &&
                      /* original InChI: non-taut & has H or an endpoint, has Fixed H */
                      ( (!pc2i->c2at[i].nFixHInChI && pc2i->c2at[i].nMobHInChI == pc2i->c2at[i].nFixHRevrs) ||
                        (pc2i->c2at[i].nFixHInChI == pc2i->c2at[i].nFixHRevrs && pc2i->c2at[i].endptInChI) ) ) &&
                     0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow) /* djb-rwth: addressing LLVM warnings */
                {

                    if ((ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_case_06d;
                    }
                }
                else
                {
                    /* -NH- */
                    if ( /* original InChI: has H and is not an endpoint */
                        ( pc2i->c2at[i].nMobHInChI + 1 == pc2i->c2at[i].nFixHInChI &&
                          pc2i->c2at[i].nFixHInChI > 0 && !pc2i->c2at[i].endptInChI &&
                          pc2i->c2at[i].nAtChargeRevrs == 0 &&
                          /* reconstructed InChI: non-taut & has H or an endpoint, has Fixed H */
                          ( (!pc2i->c2at[i].nFixHRevrs && pc2i->c2at[i].nMobHRevrs == pc2i->c2at[i].nFixHInChI) ||
                            (pc2i->c2at[i].nFixHRevrs == pc2i->c2at[i].nFixHInChI && pc2i->c2at[i].endptRevrs) ) ) &&
                         0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden &&
                         pBNS->edge[e].flow) /* djb-rwth: addressing LLVM warnings */
                    {

                        if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_case_06d;
                        }
                    }
                }
            }
            if ((num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                int bSFlowerEdgesMayBeForbidden = ( SFlowerEdges.num_edges > 0 );
                int bSFlowerEdgesIsForbidden;
                for (bSFlowerEdgesIsForbidden = bSFlowerEdgesMayBeForbidden;
                      0 <= bSFlowerEdgesIsForbidden; bSFlowerEdgesIsForbidden--)
                {
                    if (bSFlowerEdgesIsForbidden)
                    {
                        /* on the 1st pass disallow -S(+)= => =S=, allow only -S(+)= => -S- */
                        SetForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
                    }
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask );
                    delta = 1;
                    for (i = 0; i < CurrEdges.num_edges && cur_success < num_try; i++)
                    {
                        e = CurrEdges.pnEdges[i];
                        pe = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                        if (!pe->flow)
                            continue;
                        pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                        pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                        pe->flow -= delta; /* add (+) to -NHm */
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;

                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                        if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                          (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                            /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if (ret > 0)
                            {
                                nNumRunBNS++;
                                cur_success++; /* 06d */
                            }
                        }
                        else
                        {
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                        INCHI_HEAPCHK
                    }
                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
                }
            }
        exit_case_06d:
            CurrEdges.num_edges = 0; /* clear current edge list */
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
            if (cur_success)
            {
                tot_succes += cur_success;
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 2)
        {
            /*--------------------------------------------------------*/
            /* case 06: restored: NHn(+)=AB-NHm  orig.: NHn-AB=NHm(+) */
            /* FixH:               1        0            0     1      */
            /* MobH:              n-1       m            n    m-1     */
            /* N = N, O, S, Se; atoms N are not tautomeric            */
            /* Solution: move (+) from NHn(+) to NHn                  */
            /*--------------------------------------------------------*/
            int num_DB_NHn_Plus = 0, num_SB_NHm_Neutr = 0, iat;
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ((pc2i->c2at[i].nValElectr == 6 ||
                    (pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1)) &&
                    !pc2i->c2at[i].endptInChI &&
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                {
                    if ( /* orig. InChI info: NHm */
                        num_SB_NHm_Neutr < MAX_DIFF_FIXH &&
                        pc2i->c2at[i].nFixHInChI == 1 && /*pc2i->c2at[i].nMobHInChI == 0 &&*/
                        /* reversed structure info: */
                        pc2i->c2at[i].nFixHRevrs == 0 && /*pc2i->c2at[i].nMobHRevrs == 1 &&*/
                        pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H && /* at2 is Fixed-H */
                        at2[iat].valence == at2[iat].chem_bonds_valence)
                    {
                        iat_SB_NHm_Neutr[num_SB_NHm_Neutr++] = iat;
                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                    else
                    {
                        if ( /* orig. InChI info: */
                            num_DB_NHn_Plus < MAX_DIFF_FIXH &&
                            pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI &&*/
                            /* reversed structure info: */
                            pc2i->c2at[i].nFixHRevrs == 1 && /*pc2i->c2at[i].nMobHRevrs ==  0 &&*/
                            pc2i->c2at[i].nAtChargeRevrs == 1 && at2[iat].num_H &&
                            at2[iat].valence < at2[iat].chem_bonds_valence)
                        {
                            iat_DB_NHn_Plus[num_DB_NHn_Plus++] = iat;
                            if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_NHm_Neutr, num_DB_NHn_Plus))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_NHm_Neutr && cur_success < num_try; i++)
                {
                    iat = iat_SB_NHm_Neutr[i];
                    pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta; /* add (+) to -NHm */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Removed (+)charge from -NHn => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 06 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (( (pc2i->nNumTgInChI > pc2i->nNumTgRevrs && pc2i->nNumTgRevrs == 1) ||
              pc2i->nNumEndpInChI < pc2i->nNumEndpRevrs ) &&
             pStruct->nNumRemovedProtonsMobHInChI == pStruct->One_ti.tni.nNumRemovedProtons &&
             pStruct->fixed_H && pStruct->endpoint && pStruct->pOne_norm_data[TAUT_YES] && pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds) /* djb-rwth: addressing LLVM warning; fixing a NULL pointer dereference */
        {
            /*----------------------------------------------------------*/
            /* case 06a: restored: N'(+)=-AB-NH    orig.: N'-=AB=NH(+)  */
            /* FixH:               0         1            0      1      */
            /* MobH:               0         0            0      0      */
            /*                    single t-group      multiple t-groups */
            /* N  = N, O, S, Se; atoms N are not tautomeric             */
            /* N' = N            atom N' is not tautomeric              */
            /* Solution: move (+) from N' to NH                         */
            /*----------------------------------------------------------*/
            int iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            /*
            inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
            */
            inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds ?
                pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds : NULL;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            EDGE_LIST CurChargeEdges;
            cur_success = 0;
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
            CurrEdges.num_edges = 0;
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                iat = nCanon2AtnoRevrs[i];
                if (pStruct->endpoint[i])
                {
                    continue;
                }
                /* -NH-, -OH */
                if (pStruct->fixed_H[i] && nMobHInChI && !nMobHInChI[i] &&
                     at2[iat].charge == 0 && at2[iat].radical == 0 &&
                     0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                     ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) )) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    goto exit_case_06a;
                }
                else
                {
                    /* >N(+)= */
                    if (at2[iat].charge == 1 && !at2[iat].num_H &&
                         pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                         atfMobile_H_Revrs && atfMobile_H_Revrs[iat].charge == 0 &&
                         0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow &&
                         ( ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ) ))
                    {
                        goto exit_case_06a;
                    }
                }
            }
            if ((num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask );
                delta = 1;
                for (i = 0; i < CurrEdges.num_edges && cur_success < num_try; i++)
                {
                    e = CurrEdges.pnEdges[i];
                    pe = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                    pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                    pe->flow -= delta; /* add (+) to -NHm */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                      (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 06a */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            }
        exit_case_06a:
            CurrEdges.num_edges = 0; /* clear current edge list */
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
            if (cur_success)
            {
                tot_succes += cur_success;
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }
        if (( (pc2i->nNumTgInChI > pc2i->nNumTgRevrs && pc2i->nNumTgRevrs == 1) ||
              pc2i->nNumEndpInChI < pc2i->nNumEndpRevrs ) &&
              ( pStruct->nNumRemovedProtonsMobHInChI == pStruct->One_ti.tni.nNumRemovedProtons ||
                pStruct->nNumRemovedProtonsMobHInChI > pStruct->One_ti.tni.nNumRemovedProtons ) &&
             pStruct->fixed_H && pStruct->endpoint && pStruct->pOne_norm_data[TAUT_YES] && pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds) /* djb-rwth: addressing LLVM warning; fixing a NULL pointer dereference */
        {
            /*----------------------------------------------------------*/
            /* case 06b: restored: X(+)=-AB-NH    orig.: X-=AB=NH(+)    */
            /* FixH:               0        1 1          0      1       */
            /* MobH:               0        0 t          0      0       */
            /*                    single t-group      multiple t-groups */
            /*                                        or no t-groupd    */
            /* N  = N, O, S, Se; atoms N are not tautomeric             */
            /* X  = O, S, Se, Te, F, Cl, Br, I; atom X is not tautomeric*/
            /* Solution: move (+) from X  to NH                         */
            /*----------------------------------------------------------*/
            int iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            /*
            inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
            */
            inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds ?
                pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds : NULL;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            EDGE_LIST CurChargeEdges;
            cur_success = 0;
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
            CurrEdges.num_edges = 0;
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                iat = nCanon2AtnoRevrs[i];
                if (pStruct->endpoint[i])
                {
                    continue;
                }
                /* -NH-, -OH */
                if (pStruct->fixed_H[i] && nMobHInChI && !nMobHInChI[i] &&
                     at2[iat].charge == 0 && at2[iat].radical == 0 &&
                     0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                     ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) )) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    goto exit_case_06b;
                }
                else
                {
                    /* X(+)= */
                    if (at2[iat].charge == 1 && !at2[iat].num_H &&
                        ( pVA[iat].cNumValenceElectrons == 6 || pVA[iat].cPeriodicRowNumber == 7 ) &&
                         atfMobile_H_Revrs && atfMobile_H_Revrs[iat].charge == 1 &&
                         0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow &&
                         ( ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ) ))
                    {
                        goto exit_case_06b;
                    }
                }
            }
            if ((num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                int bSFlowerEdgesMayBeForbidden = ( SFlowerEdges.num_edges > 0 );
                int bSFlowerEdgesIsForbidden;
                for (bSFlowerEdgesIsForbidden = bSFlowerEdgesMayBeForbidden;
                      0 <= bSFlowerEdgesIsForbidden; bSFlowerEdgesIsForbidden--)
                {
                    if (bSFlowerEdgesIsForbidden)
                    {
                        /* on the 1st pass disallow -S(+)= => =S=, allow only -S(+)= => -S- */
                        SetForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
                    }
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask );
                    delta = 1;
                    for (i = 0; i < CurrEdges.num_edges && cur_success < num_try; i++)
                    {
                        e = CurrEdges.pnEdges[i];
                        pe = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                        if (!pe->flow)
                            continue;
                        pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                        pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                        pe->flow -= delta; /* add (+) to -NHm */
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;

                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                        if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                          (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                            /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if (ret > 0)
                            {
                                nNumRunBNS++;
                                cur_success++; /* 06b */
                            }
                        }
                        else
                        {
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                        INCHI_HEAPCHK
                    }
                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
                }
            }
        exit_case_06b:
            CurrEdges.num_edges = 0; /* clear current edge list */
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
            if (cur_success)
            {
                tot_succes += cur_success;
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->nNumTgInChI > 1 &&
            ( pStruct->nNumRemovedProtonsMobHInChI > 0 || pStruct->ti.tni.nNumRemovedProtons > 0 ) &&
             pStruct->fixed_H && pStruct->endpoint &&
             pStruct->pOne_norm_data[TAUT_YES] && pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds)
        {
            /*----------------------------------------------------------*/
            /* case 06e:restored: XHn(+)=-AB-YHm  orig.: XHn-=AB=YHm(+) */
            /* FixH:               1          0           1      1      */
            /* MobH:               0          1           t      t      */
            /*                   non-taut atoms       multiple t-groups */
            /*                                                          */
            /* 1. orig. t-group has more H on its endpoints counted     */
            /*          in atf and has no (+) on endpoint that has H    */
            /* 2. orig. t-group has less H on its endpoints counted     */
            /*          in atf and has (+) on endpoint that has H       */
            /*          in reconstructed struct and less H in atf       */
            /* Solution: move (+) from (2) to atom in (1) that has H    */
            /*                                                          */
            /*   tg1  reconstr:   XHn and more H than in orig t-group   */
            /*             atf:   XHn                                   */
            /*   tg2  reconstr:   XHm(+) and less H than in             */
            /*             atf:   XH(m-1)           orig in t-group     */
            /*                                                          */
            /* N  = N, O, S, Se; atoms N are not tautomeric             */
            /* X  = O, S, Se, Te, F, Cl, Br, I; atom X is not tautomeric*/
            /* Solution: move (+) from X  to NH                         */
            /*----------------------------------------------------------*/

            int iat, nNumWrongTg, jjoffs, jj, nNum2RemovePlus, nNum2AddPlus, nNum2MovePlus;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            /*
            inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
            */
            inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds ?
                pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds :
                pStruct->pOne_norm_data[TAUT_YES] &&
                pStruct->pOne_norm_data[TAUT_YES]->at ?
                pStruct->pOne_norm_data[TAUT_YES]->at : NULL;
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            EDGE_LIST CurChargeEdges /* source of (+)*/, EndpList;
            BNS_VERTEX *pv1n, *pv2n;
            BNS_EDGE   *pe1n, *pe2n;
            Vertex      v1n, v2n;

            cur_success = 0;
            AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
            AllocEdgeList( &EndpList, EDGE_LIST_CLEAR );
            CurrEdges.num_edges = 0; /* receptors of (+) */

            if (!atfMobile_H_Revrs)
            {
                goto exit_case_06e;
            }
            nNumWrongTg = FillTgDiffHChgFH(tdhc, MAX_DIFF_FIXH, at2, atfMobile_H_Revrs,
                nCanon2AtnoRevrs, pVA, &pStruct->ti, &EndpList);
            if (nNumWrongTg < 1)
            {
                goto exit_case_06e; /* for now only transfer (+) from one Mobile-H group to another */
            }
            nNum2RemovePlus = nNum2AddPlus = 0; /* djb-rwth: removing redundant code */
            for (i = 0; i < nNumWrongTg; i++)
            {
                /* detect t-group that has extra (+) on H */
                if (tdhc[i].nNumHInchi > tdhc[i].nNumHNorml &&
                    tdhc[i].nNumPRevrs > tdhc[i].nNumPNorml && tdhc[i].n[fNumRPosChgH])
                {
                    /* count how many (+) to remove */
                    /* store XH(+) atom numbers */
                    int nNumNeeded = inchi_min(tdhc[i].nNumHInchi - tdhc[i].nNumHNorml, tdhc[i].n[fNumRPosChgH]);
                    nNum2RemovePlus += nNumNeeded;
                    jjoffs = tdhc[i].i[fNumRPosChgH];
                    for (jj = 0; jj < tdhc[i].n[fNumRPosChgH]; jj++)
                    {
                        iat = EndpList.pnEdges[jjoffs + jj];
                        e = pVA[iat].nCPlusGroupEdge - 1;
                        if ((ret = AddToEdgeList(&CurChargeEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_case_06e;
                        }
                    }
                }
                else
                {
                    /* detect t-group that needs (+) on XH to reduce number of H */
                    if (tdhc[i].nNumHInchi < tdhc[i].nNumHNorml && tdhc[i].n[fNumRNeutrlH])
                    {
                        /* store XH atom numbers */
                        int nNumNeeded = inchi_min(tdhc[i].nNumHNorml - tdhc[i].nNumHInchi, tdhc[i].n[fNumRNeutrlH]);
                        nNum2AddPlus += nNumNeeded;
                        jjoffs = tdhc[i].i[fNumRNeutrlH];
                        for (jj = 0; jj < tdhc[i].n[fNumRNeutrlH]; jj++)
                        {
                            iat = EndpList.pnEdges[jjoffs + jj];
                            e = pVA[iat].nCPlusGroupEdge - 1;
                            if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_case_06e;
                            }
                        }
                    }
                }
            }
            nNum2MovePlus = inchi_min(nNum2RemovePlus, nNum2AddPlus);
            if (CurrEdges.num_edges > 0 && CurChargeEdges.num_edges > 0)
            {
                for (i = 0; 0 < nNum2MovePlus && i < nNumWrongTg; i++)
                {
                    /* detect t-group that has extra (+) on H */
                    if (tdhc[i].nNumHInchi > tdhc[i].nNumHNorml &&
                        tdhc[i].nNumPRevrs > tdhc[i].nNumPNorml && tdhc[i].n[fNumRPosChgH])
                    {
                        int nNum2Remove = tdhc[i].nNumHInchi - tdhc[i].nNumHNorml;
                        if (nNum2Remove < tdhc[i].n[fNumRPosChgH])
                        {
                            nNum2Remove = tdhc[i].n[fNumRPosChgH];
                        }
                        /* store XH(+) atom numbers */
                        jjoffs = tdhc[i].i[fNumRPosChgH];
                        SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                        RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                        for (jj = 0; 0 < nNum2MovePlus && 0 < nNum2Remove && jj < tdhc[i].n[fNumRPosChgH]; jj++)
                        {
                            iat = EndpList.pnEdges[jjoffs + jj];
                            e = pVA[iat].nCPlusGroupEdge - 1;
                            pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                            if (pe->flow)
                                continue;
                            pv1 = pBNS->vert + (v1 = pe->neighbor1);
                            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                            for (j = pv1->num_adj_edges - 1; 0 <= j; j--)
                            {
                                pe1n = pBNS->edge + pv1->iedge[j];
                                if (pe1n->flow && !pe1n->forbidden)
                                {
                                    pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                                    break;
                                }
                            }
                            if (j < 0)
                            {
                                continue; /* not found */
                            }
                            for (j = pv2->num_adj_edges - 2; 0 <= j; j--)
                            {
                                pe2n = pBNS->edge + pv2->iedge[j];
                                if (pe2n->flow && !pe2n->forbidden)
                                {
                                    pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                                    break;
                                }
                            }
                            if (j < 0)
                            {
                                continue; /* not found */
                            }
                            delta = 1;
                            pe->flow += delta;
                            pe1n->flow -= delta;
                            pe2n->flow -= delta;
                            pv1n->st_edge.flow -= delta;
                            pv2n->st_edge.flow -= delta;
                            pBNS->tot_st_flow -= 2 * delta;

                            ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                            if (ret == 1 && ((vPathEnd == v1n && vPathStart == v2n) ||
                                (vPathEnd == v2n && vPathStart == v1n)) &&
                                (nDeltaCharge == 0 || nDeltaCharge == 1)) /* djb-rwth: addressing LLVM warnings */
                            {
                                ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                                if (ret > 0)
                                {
                                    nNumRunBNS++;
                                    nNum2Remove--;
                                    nNum2MovePlus--;
                                    cur_success++; /* 06e */
                                }
                            }
                            else
                            {
                                pe->flow -= delta;
                                pe1n->flow += delta;
                                pe2n->flow += delta;
                                pv1n->st_edge.flow += delta;
                                pv2n->st_edge.flow += delta;
                                pBNS->tot_st_flow += 2 * delta;
                            }
                            if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE)))
                            {
                                goto exit_case_06e;
                            }
                        }
                        RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                    }
                }
            }
        exit_case_06e:
            CurrEdges.num_edges = 0; /* clear current edge list */
            AllocEdgeList(&CurChargeEdges, EDGE_LIST_FREE);
            AllocEdgeList(&EndpList, EDGE_LIST_FREE);
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 1)
        {
            /*--------------------------------------------------------------*/
            /* case 07: restored:  O(-)-AB=O  original:  O=AB-O(-)          */
            /* FixH:               0       0             0     -1           */
            /* MobH:               0       0             0      1           */
            /*                    taut  (non-taut)     (taut) non-taut      */
            /*                    taut  (taut)     (non-taut) non-taut      */
            /* O = O, S, Se, Te                                             */
            /* Solution: move (-) from O(-)-AB to AB=O                      */
            /*--------------------------------------------------------------*/
            int num_SB_O_Minus = 0, num_DB_O_Neutr = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: -O(-), non-taut */
                    num_DB_O_Neutr < MAX_DIFF_FIXH &&
                    pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                    !pc2i->c2at[i].endptInChI &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI == 1 &&
                    /* reversed structure info: */
                    pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                    pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                    at2[iat].valence < at2[iat].chem_bonds_valence)
                {
                    iat_DB_O_Neutr[num_DB_O_Neutr++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=-1, no H, has single bond, O, S, Se, Te */
                    num_SB_O_Minus < MAX_DIFF_FIXH &&
                    at2[iat].charge == -1 && !at2[iat].num_H &&
                    at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    pVA[iat].cNumValenceElectrons == 6 &&
                    at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint &&
                    /* in orig.InChI: not an endpoint, has no H */
                    /*pStruct->endpoint[i] && -- modificatuion#1 */
                    !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    !(nMobHInChI && nMobHInChI[i]) &&
                    /* has (-) edge */
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                {
                    iat_SB_O_Minus[num_SB_O_Minus++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_O_Minus, num_DB_O_Neutr))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_O_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_O_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Moved (-) charge to AB=O => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (O(-)-AB) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 07 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 1)
        {
            /*--------------------------------------------------------------*/
            /* case 07a: restored: O(-)-N(V)B=O  original:  O=N(V)B-O(-)    */
            /* FixH:               0          0             0      -1       */
            /* MobH:               0          0             0       1       */
            /*                non-taut  (non-taut)  non-taut  non-taut      */
            /*                non-taut  (taut)      non-taut  non-taut      */
            /* O = O, S, Se, Te                                             */
            /* Solution: move (-) from O(-)-AB to AB=O                      */
            /*--------------------------------------------------------------*/
            int num_SB_O_Minus = 0, num_DB_O_Neutr = 0, iat, iN;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: -O(-), non-taut */
                    num_DB_O_Neutr < MAX_DIFF_FIXH &&
                    pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                    !pc2i->c2at[i].endptInChI &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI == 1 &&
                    /* reversed structure info: */
                    pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                    pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                    at2[iat].valence < at2[iat].chem_bonds_valence)
                {
                    iat_DB_O_Neutr[num_DB_O_Neutr++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=-1, no H, has single bond, O, S, Se, Te */
                    num_SB_O_Minus < MAX_DIFF_FIXH &&
                    at2[iat].charge == -1 && !at2[iat].num_H &&
                    at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    pVA[iat].cNumValenceElectrons == 6 &&
                    /*at_Mobile_H_Revrs && !at_Mobile_H_Revrs[iat].endpoint &&*/
                    /* in orig.InChI: not an endpoint, has no H */
                    !pStruct->endpoint[i] &&
                    !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    !(nMobHInChI && nMobHInChI[i]) &&
                    /* has N(V) neighbor */
                    1 == at2[iat].valence && at2[iN = at2[iat].neighbor[0]].chem_bonds_valence == 5 &&
                    !at2[iN].charge && pVA[iN].cNumValenceElectrons == 5 &&
                    /* has (-) edge */
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                {
                    iat_SB_O_Minus[num_SB_O_Minus++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_O_Minus, num_DB_O_Neutr))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_O_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_O_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Moved (-) charge to AB=O => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (O(-)-AB) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 07 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if ( /*(pc2i->len_c2at >= 1 || pc2i->nNumRemHRevrs) &&*/ pc2i->nNumTgInChI == 1 && /* ADP in InChI */
            (pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1))
        {
            /*----------------------------------------------------------------*/
            /* case 08: restored: O(-)-AB=N- OH- orig.   O=AB-N(-)- OH-       */
            /* FixH:               1      0   0          0      0   1         */
            /* MobH:               0      0   1          0      0   0         */
            /*           may be taut or not  non-taut   taut  taut taut       */
            /*                                    ADP: one t-group or more endpoints */
            /* O(-) = S, Se, Te; N = N;                                       */
            /* Solution: move (-) from O(-) to =N-; avoid stereogenic DB on N */
            /*----------------------------------------------------------------*/
            int num_DB_N_Neutr = 0, num_SB_O_Minus = 0, iat;

            AT_NUMB* nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR* nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0;

            for (i = 0; i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=-1, has no H, has single bond, O, S, Se, Te */
                    num_SB_O_Minus < MAX_DIFF_FIXH &&
                    at2[iat].charge == -1 && !at2[iat].num_H &&
                    at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    pVA[iat].cNumValenceElectrons == 6 &&
                    /* in orig.InChI: an endpoint, may have fixed-H */
                    pStruct->endpoint[i] &&
                    /*!(pStruct->fixed_H && pStruct->fixed_H[i]) &&*/
                    !(nMobHInChI && nMobHInChI[i]) &&
                    /* has (-) edge */
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                {

                    iat_SB_O_Minus[num_SB_O_Minus++] = iat;

                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
                else
                {
                    if ( /* in restored atom: charge=0, has no H, has double non-stereogenic bond, N */
                        num_DB_N_Neutr < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H && !at2[iat].sb_parity[0] &&
                        at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                        /* in orig.InChI: an endpoint, has no fixed-H */
                        pStruct->endpoint[i] &&
                        !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                        !(nMobHInChI && nMobHInChI[i]) &&
                        /* has (-) edge */
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 &&
                        0 == pBNS->edge[e].forbidden)
                    {

                        iat_DB_N_Neutr[num_DB_N_Neutr++] = iat;

                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_DB_N_Neutr, num_SB_O_Minus))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                /* allow stereobonds in rings change */
                if (forbidden_stereo_edge_mask)
                    RemoveForbiddenEdgeMask(pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask);

                delta = 1;
                for (i = 0; i < num_SB_O_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_O_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Moved (-) charge to =N- => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom (-)O-) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 08 */
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                if (forbidden_stereo_edge_mask)
                    SetForbiddenEdgeMask(pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 2)
        {
            /*--------------------------------------------------------*/
            /* case 09: restored: NH2(+)=C--NH2 orig.: NH2-C(+)-NH2   */
            /* FixH:               2     |  2            0  |   0     */
            /* MobH:               0        0            2      2     */
            /* N = N,            taut      taut     non-taut  non-taut*/
            /* Solution: move (+) from NH2(+) to C                    */
            /*--------------------------------------------------------*/
            int iat;
            cur_success = 0;
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if (( pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1 ) &&
                     /* orig. InChI info: */
                     !pc2i->c2at[i].endptInChI &&
                     pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI &&
                     /* reversed structure info: */
                     pc2i->c2at[i].endptRevrs &&
                     pc2i->c2at[i].nFixHRevrs && !pc2i->c2at[i].nMobHRevrs &&
                     pc2i->c2at[i].nAtChargeRevrs == 1 &&
                     at2[iat].valence + 1 == at2[iat].chem_bonds_valence &&
                     ( e = pVA[iat].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden)
                {
                    EdgeIndex eNC = NO_VERTEX, eCPlusC;
                    int       iNH2, iatC, iatNH2; /* djb-rwth: removing redundant variables */
                    /* found NH2(+)=; locate =C< and find whether it has -NH2 neighbor */
                    for (j = 0; j < at2[iat].valence; j++)
                    {
                        if (at2[iat].bond_type[j] == BOND_TYPE_DOUBLE)
                            break;
                    }
                    if (j == at2[iat].valence)
                    {
                        continue;
                    }
                    eNC = pBNS->vert[iat].iedge[j]; /* edge NH2(+)=C */
                    iatC = at2[iat].neighbor[j];
                    if (pVA[iatC].cNumValenceElectrons != 4 || pVA[iatC].cMetal || at2[iatC].charge ||
                         at2[iatC].valence != 3 || at2[iatC].valence + 1 != at2[iatC].chem_bonds_valence ||
                         ( eCPlusC = pVA[iatC].nCPlusGroupEdge - 1 ) < 0 || pBNS->edge[eCPlusC].forbidden)
                    {
                        continue;
                    }
                    for (j = 0; j < at2[iatC].valence; j++)
                    {
                        iatNH2 = at2[iatC].neighbor[j];
                        if (iatNH2 == iat || pVA[iatNH2].cNumValenceElectrons != 5 ||
                             pVA[iatNH2].cPeriodicRowNumber != 1 || !at2[iatNH2].num_H || at2[iatNH2].charge)
                            continue;
                        /* djb-rwth: removing redundant code */
                        for (iNH2 = 0; iNH2 < pc2i->len_c2at; iNH2++)
                        {
                            if (iatNH2 == pc2i->c2at[iNH2].atomNumber)
                                break;
                        }
                        if (iNH2 == pc2i->len_c2at)
                        {
                            continue;
                        }

                        if (( pc2i->c2at[iNH2].nValElectr == 5 && pc2i->c2at[iNH2].nPeriodNum == 1 ) &&
                             /* orig. InChI info: */
                             !pc2i->c2at[iNH2].endptInChI &&
                             pc2i->c2at[iNH2].nFixHInChI == 0 && pc2i->c2at[iNH2].nMobHInChI &&
                             /* reversed structure info: */
                             pc2i->c2at[iNH2].endptRevrs &&
                             pc2i->c2at[iNH2].nFixHRevrs && !pc2i->c2at[iNH2].nMobHRevrs &&
                             pc2i->c2at[iNH2].nAtChargeRevrs == 0 &&
                             at2[iatNH2].valence == at2[iatNH2].chem_bonds_valence)
                        {
                            /* we have found NH2(+)=, =C<, and bond between them */

                            if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                            if ((ret = AddToEdgeList( &CurrEdges, eCPlusC, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
                            delta = 1;

                            pe = pBNS->edge + eNC;
                            if (!pe->flow)
                                continue;
                            pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                            pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                            pe->forbidden |= forbidden_edge_mask;
                            pe->flow -= delta; /* add (+) to -NHm */
                            pv1->st_edge.flow -= delta;
                            pv2->st_edge.flow -= delta;
                            pBNS->tot_st_flow -= 2 * delta;

                            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                            if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                              (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 0) /* djb-rwth: addressing LLVM warnings */
                            {
                                /* Removed (+)charge from -NHn => nDeltaCharge == -1 */
                                /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                                if (ret > 0)
                                {
                                    nNumRunBNS++;
                                    cur_success++; /* 09 */
                                }
                            }
                            else
                            {
                                pe->flow += delta;
                                pv1->st_edge.flow += delta;
                                pv2->st_edge.flow += delta;
                                pBNS->tot_st_flow += 2 * delta;
                            }
                            INCHI_HEAPCHK
                                pe->forbidden &= forbidden_edge_mask_inv;
                            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                            CurrEdges.num_edges = 0; /* clear current edge list */
                            break;
                        }
                    }
                }
            }
            if (cur_success)
            {
                tot_succes += cur_success;
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 2)
        {
            /*--------------------------------------------------------*/
            /* case 10: restored: NH2-X(+)-NH-  orig.: NH2(+)=X-NH-   */
            /* FixH:               0        0            2      1     */
            /* MobH:               2        1            0      0     */
            /* N = N,O,S,Se,Te non-taut  non-taut       taut   taut   */
            /* Solution: move (+) from X(+) to NH2 or NH              */
            /*--------------------------------------------------------*/
            int iat;
            cur_success = 0;
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                if (pc2i->c2at[i].nValue)
                    continue;
                iat = pc2i->c2at[i].atomNumber;
                if (( pc2i->c2at[i].nValElectr == 6 ||
                      (pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) ) &&
                     /* orig. InChI info: */
                     pc2i->c2at[i].endptInChI &&
                     pc2i->c2at[i].nFixHInChI && !pc2i->c2at[i].nMobHInChI &&
                     /* reversed structure info: */
                     !pc2i->c2at[i].endptRevrs &&
                     !pc2i->c2at[i].nFixHRevrs && pc2i->c2at[i].nMobHRevrs &&
                     pc2i->c2at[i].nAtChargeRevrs == 0 &&
                     at2[iat].valence == at2[iat].chem_bonds_valence &&
                     ( e = pVA[iat].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                {

                    EdgeIndex eCPlusC, eCPlusNH2, bContinue = 1;
                    int       iNH2, iatC, iatNH2, j1, j2; /* djb-rwth: removing redundant variables */
                    BNS_EDGE *pe_iat, *pe_iNH2;
                    /* found NH2- locate -X(+) and find whether it has another -NH2 neighbor */
                    for (j1 = 0; j1 < at2[iat].valence && bContinue; j1++)
                    {
                        if (at2[iat].bond_type[j1] == BOND_TYPE_SINGLE &&
                             at2[iatC = at2[iat].neighbor[j1]].charge == 1 &&
                             ( 4 <= pVA[iatC].cNumValenceElectrons && pVA[iatC].cNumValenceElectrons <= 6 ) &&
                             at2[iatC].valence == at2[iatC].chem_bonds_valence &&
                             ( eCPlusC = pVA[iatC].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[eCPlusC].forbidden)
                        {
                            /* found a candidate for X; find another NH2 */
                            for (j2 = 0; j2 < at2[iatC].valence && bContinue; j2++)
                            {
                                if (at2[iatC].bond_type[j2] == BOND_TYPE_SINGLE &&
                                     iat != ( iatNH2 = at2[iatC].neighbor[j2] ) &&
                                     at2[iatNH2].charge == 0 && at2[iatNH2].num_H &&
                                     ( pVA[iatNH2].cNumValenceElectrons == 5 || pVA[iatNH2].cNumValenceElectrons == 6 ) &&
                                     at2[iatNH2].valence == at2[iatNH2].chem_bonds_valence &&
                                     ( eCPlusNH2 = pVA[iatNH2].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[eCPlusNH2].forbidden)
                                {
                                    for (iNH2 = 0; iNH2 < pc2i->len_c2at; iNH2++)
                                    {
                                        if (iatNH2 != pc2i->c2at[iNH2].atomNumber || pc2i->c2at[iNH2].nValue)
                                            continue;
                                        /* check the second -NH */
                                        /* djb-rwth: removing redundant code */
                                        if ( /* orig. InChI info: */
                                             pc2i->c2at[iNH2].endptInChI &&
                                             pc2i->c2at[iNH2].nFixHInChI && !pc2i->c2at[iNH2].nMobHInChI &&
                                             /* reversed structure info: */
                                             !pc2i->c2at[iNH2].endptRevrs &&
                                             !pc2i->c2at[iNH2].nFixHRevrs && pc2i->c2at[iNH2].nMobHRevrs &&
                                             pc2i->c2at[iNH2].nAtChargeRevrs == 0)
                                        {
                                            /* we have found NH-X(+)-NH; remove charge from X(+) */
                                            pe_iat = pBNS->edge + pBNS->vert[iat].iedge[j1];
                                            pe_iNH2 = pBNS->edge + pBNS->vert[iatC].iedge[j2];
                                            /* pick up one of -NH to move (+) to it */
                                            if (!pe_iat->forbidden && pBNS->edge[e].flow)
                                            {
                                                pe = pBNS->edge + e;
                                            }
                                            else
                                                if (!pe_iNH2->forbidden && pBNS->edge[eCPlusNH2].flow)
                                                {
                                                    pe = pBNS->edge + eCPlusNH2;
                                                }
                                                else
                                                {
                                                    continue; /* none of the two -X(+)- bonds may be changed */
                                                }
                                            if ((ret = AddToEdgeList( &CurrEdges, eCPlusC, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                            {
                                                goto exit_function;
                                            }
                                            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                                            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
                                            delta = 1;

                                            pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                                            pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                                            /*pe->forbidden |= forbidden_edge_mask;*/
                                            pe->flow -= delta; /* add (+) to -NHm */
                                            pv1->st_edge.flow -= delta;
                                            pv2->st_edge.flow -= delta;
                                            pBNS->tot_st_flow -= 2 * delta;

                                            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                                            if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                                              (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == -1) /* djb-rwth: addressing LLVM warnings */
                                            {
                                                /* Removed (+)charge from -NHn => nDeltaCharge == -1 */
                                                /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                                                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                                                if (ret > 0)
                                                {
                                                    nNumRunBNS++;
                                                    cur_success++; /* 10 */
                                                    bContinue = 0;
                                                    pc2i->c2at[i].nValue = 1;    /* mark as used */
                                                    pc2i->c2at[iNH2].nValue = 1; /* mark as used */
                                                }
                                            }
                                            else
                                            {
                                                pe->flow += delta;
                                                pv1->st_edge.flow += delta;
                                                pv2->st_edge.flow += delta;
                                                pBNS->tot_st_flow += 2 * delta;
                                            }
                                            INCHI_HEAPCHK

                                                /*pe->forbidden &= forbidden_edge_mask_inv;*/
                                                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                                            CurrEdges.num_edges = 0; /* clear current edge list */
                                            break;
                                        }
                                    } /* iNH2: pc2i->c2at[iNH2] cycle */
                                }
                            } /* j2: iatC neighbors cycle */
                        }
                    } /* j1: iat neighbors cycle */
                }
            } /* i: pc2i->c2at[i] cycle */

            if (cur_success)
            {
                /*
                for ( i = 0; i < pc2i->len_c2at; i ++ ) {
                pc2i->c2at[i].nValue = 0;
                }
                */
                tot_succes += cur_success;
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if ( /*pc2i->len_c2at >= 1 &&*/ pc2i->nNumTgInChI == 1 && /* ADP in InChI */
            ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ))
        {
            /*--------------------------------------------------------------*/
            /* case 11: restored: NH(+)=AB-N< OH- orig.  NH-AB=N(+)< OH-    */
            /* FixH:               0       0   0          1      0   1      */
            /* MobH:               1       0   1          0      0   0      */
            /*                 non-taut.                taut        taut    */
            /*                                    ADP: one t-group or more endpoints */
            /* NH(+)= => N, O, S, Se; -N< => N                              */
            /* Solution: move (+) from NH(+) to -N<                         */
            /*--------------------------------------------------------------*/
            int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0;

            /* search for NH(+)= */
    /* search for -N< */
            for (i = 0; i < pStruct->num_atoms; i++)
            { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=0, has no H, has no double bond, N only */
                    num_DB_Charged < MAX_DIFF_FIXH &&
                    at2[iat].charge == 1 && at2[iat].num_H &&
                    at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    ((pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1) ||
                        pVA[iat].cNumValenceElectrons == 6) &&
                    /* in orig.InChI: an endpoint, has fixed-H */
                    /*pStruct->endpoint[i] &&*/
                    (pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                    /* has (+) edge */
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && 0 == pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                {

                    iat_DB_Charged[num_DB_Charged++] = iat;
                    /*
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                    }
                    */
                }
                else
                {
                    if ( /* in restored atom: charge=0, has no H, has no double bond, N only */
                        num_SB_Neutr < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H &&
                        at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        (pVA[iat].cNumValenceElectrons == 5 &&
                            pVA[iat].cPeriodicRowNumber == 1) &&
                        /* in orig.InChI: an endpoint, has fixed-H */
                        /*pStruct->endpoint[i] &&*/
                        !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                        !(nMobHInChI && nMobHInChI[i]) &&
                        /* has (+) edge */
                        (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && 0 == pBNS->edge[e].forbidden)
                    {

                        iat_SB_Neutr[num_SB_Neutr++] = iat;
                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_Neutr, num_DB_Charged))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                BNS_VERTEX* pv1n, * pv2n;
                BNS_EDGE* pe1n, * pe2n;
                Vertex      v1n, v2n;
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_DB_Charged && cur_success < num_try; i++)
                {
                    iat = iat_DB_Charged[i];
                    pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                    if (pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    for (j = pv1->num_adj_edges - 1; 0 <= j; j--)
                    {
                        pe1n = pBNS->edge + pv1->iedge[j];
                        if (pe1n->flow && !pe1n->forbidden)
                        {
                            pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                            break;
                        }
                    }
                    if (j < 0)
                    {
                        continue; /* not found */
                    }

                    for (j = pv2->num_adj_edges - 2; 0 <= j; j--)
                    {
                        pe2n = pBNS->edge + pv2->iedge[j];
                        if (pe2n->flow && !pe2n->forbidden)
                        {
                            pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                            break;
                        }
                    }
                    if (j < 0)
                    {
                        continue; /* not found */
                    }

                    pe->flow += delta;
                    pe1n->flow -= delta;
                    pe2n->flow -= delta;
                    pv1n->st_edge.flow -= delta;
                    pv2n->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1n && vPathStart == v2n) ||
                        (vPathEnd == v2n && vPathStart == v1n)) &&
                        (nDeltaCharge == 0 || nDeltaCharge == 1)) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* before setting flows the structure could be:
                        [NH+ neigh, v1n]=e1n=[NH+,v1]-pe-[+,v2]=e2n=[another at or its chargeStruct]
                        or

                        [NH+ or ChStr, v1n]=pe1n=[NH+ or ChStr, v1]-pe-[+,v2]=pe2n=[at2 or ChStr, v2n]
                        ^    ^    ^
                        NH+(+)edge |  N (+) edge: only
                        |  these are not forbidden
                        |
                        hetero (+) vertex

                        After setting flows (* mark radicals, =pe= is forbidden):

                        *[NH+ or ChStr, v1n]-pe1n-[NH+ or ChStr, v1]=pe=[+,v2]-pe2n-[at2 or ChStr, v2n]*
                        ^    ^    ^
                        NH+(+)edge |  N (+) edge: only
                        |  these are not forbidden
                        |
                        hetero (+) vertex

                        Flow in
                        pe1n and pe2n will or will not change, depending on the structure.

                        Consider what happens if pe2n changes. It may only increment.
                        If pe2n flow increments then another (+)edge flow dectrements. If
                        [at2 or ChStr, v2n] is at2 then at2 charge would change from (+) to 0,
                        and another N charge would change from 0 to (+), giving tot. change of
                        number of charges  (-1)+(+1)=0. However, if [at2 or ChStr, v2n] is
                        ChargeStruct then at2 will not be on the alt path and only the creation
                        of another (+) will be detected.
                        */
                        /* Removed charge from O(+) => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 11 */
                        }
                    }
                    else
                    {
                        pe->flow -= delta;
                        pe1n->flow += delta;
                        pe2n->flow += delta;
                        pv1n->st_edge.flow += delta;
                        pv2n->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 1 && pc2i->nNumTgInChI == 1 &&
             pc2i->nNumRemHInChI >= -1 && /* 2006-03-03 */
             ( pc2i->nNumEndpInChI > pc2i->nNumEndpRevrs || pc2i->nNumTgRevrs > 1 ) /* ADP in InChI */)
        {
            /*--------------------------------------------------------------*/
            /* case 12: restored: O=AB-N<         original: (-)O-AB=N(+)<   */
            /* FixH:              0    0                     0        0     */
            /* MobH:              0    0                     0        0     */
            /*                   non-taut                   taut            */
            /* O = O, S, Se, N; N = N;                                         */
            /* restored atom O is not tautomeric; original atom O is taut.  */
            /* original struct has 1 t-group; restored has less endpoints   */
            /*                             and/or possibly >1 t-groups      */
            /* Solution: separate charges between O= and -N<                */
            /*           allow moving charge to N(V) to make it N(IV)(+)    */
            /*--------------------------------------------------------------*/
            int bOnly_N_V = 1;
            /* djb-rwth: removing redundant code */
            while (1)
            {
                int num_SB_N_Neutr = 0, num_DB_O = 0, iat, num_N_V = 0, bN_V;
                AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
                inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                                pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
                S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                    pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
                cur_success = 0;

                for (i = 0; i < pc2i->len_c2at; i++)
                {
                    iat = pc2i->c2at[i].atomNumber;
                    if ( /* orig. InChI info: -O(-) */
                        num_DB_O < MAX_DIFF_FIXH &&
                        (pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ ||
                            (pc2i->c2at[i].nValElectr == 5 &&
                            pc2i->c2at[i].nPeriodNum == 1) /* N */) && /* djb-rwth: addressing LLVM warning */
                        pc2i->c2at[i].endptInChI &&
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                        pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                        /* reversed structure info: */
                        !pc2i->c2at[i].endptRevrs &&
                        pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                        pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                        ((pc2i->c2at[i].nValElectr == 6) ?
                            (at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2) :
                            (pc2i->c2at[i].nValElectr == 5) ?
                            (at2[iat].valence == 2 && at2[iat].chem_bonds_valence == 3) :
                            0))
                    {

                        iat_DB_O[num_DB_O++] = iat;
                        /*
                        if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                        }
                        */
                    }
                }
                for (i = 0; i < pStruct->num_atoms; i++)
                {
                    /* i = canonical number - 1 */
                    iat = nCanon2AtnoRevrs[i];
                    bN_V = 0;
                    if ( /* in restored atom N: charge=0, no H, has no double bond, not an endpoint */
                        num_SB_N_Neutr < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H &&
                        (at2[iat].valence == at2[iat].chem_bonds_valence ||
                            (bN_V = at2[iat].valence + 2 == at2[iat].chem_bonds_valence)) &&
                        !pVA[iat].cMetal &&
                        pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                        !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                        /* in orig.InChI: not an endpoint, has no H */
                        !pStruct->endpoint[i] &&
                        !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                        !(nMobHInChI && nMobHInChI[i]) &&
                        /* has (+) edge */
                        (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                    {

                        if (bOnly_N_V && bN_V &&
                            NO_VERTEX != (j = GetChargeFlowerUpperEdge(pBNS, pVA, e)) &&
                            !pBNS->edge[j].forbidden && !pBNS->edge[j].flow)
                        {
                            if (!num_N_V)
                            {
                                /* switch to N(V) only mode */
                                CurrEdges.num_edges = 0;
                                num_SB_N_Neutr = 0;
                            }
                            iat_SB_N_Neutr[num_SB_N_Neutr++] = iat;
                            num_N_V++;
                            if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                            if ((ret = AddToEdgeList(&CurrEdges, j, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                        else
                        {
                            if (!num_N_V)
                            {
                                iat_SB_N_Neutr[num_SB_N_Neutr++] = iat;
                                if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                                /* in addition, permit N(V)=>N(IV)(+) change by allowing charge flower edge change flow */
                                if (bN_V && NO_VERTEX != (j = GetChargeFlowerUpperEdge(pBNS, pVA, e)) &&
                                    !pBNS->edge[j].forbidden && !pBNS->edge[j].flow)
                                {
                                    if ((ret = AddToEdgeList(&CurrEdges, j, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                            }
                        }
                    }
                }
                if ((num_try = inchi_min(num_SB_N_Neutr, num_DB_O))) /* djb-rwth: addressing LLVM warning */
                {
                    /* detected; attempt to fix */
                    BNS_EDGE* pe_CMinus;
                    SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                    RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                    delta = 1;
                    for (i = 0; i < num_DB_O && cur_success < num_try; i++)
                    {
                        iat = iat_DB_O[i];
                        pe_CMinus = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                        pe_CMinus->forbidden &= forbidden_edge_mask_inv;

                        pe = pBNS->edge + pBNS->vert[iat].iedge[0]; /* double bond O=...*/
                        if (!pe->flow)
                            continue;
                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                        pe->forbidden |= forbidden_edge_mask; /* change bond O=X to O(rad)-X(rad) */
                        pe->flow -= delta;
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;

                        ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                            &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                        if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 2) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* Added (-) charge to =O and (+) charge to N => nDeltaCharge == 2 */
                            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                            if (ret > 0)
                            {
                                nNumRunBNS++;
                                cur_success++; /* 12 */
                            }
                        }
                        else
                        {
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                        pe->forbidden &= forbidden_edge_mask_inv; /* allow changes to O=X bond */
                        INCHI_HEAPCHK
                    }
                    RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                }
                CurrEdges.num_edges = 0; /* clear current edge list */
                if (cur_success)
                {
                    tot_succes += cur_success;
                    /* recalculate InChI from the structure */
                    if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                        ppt_group_info, ppat_norm, ppat_prep)))
                    {
                        goto exit_function;
                    }
                    if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                    if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                    {
                        goto exit_function;  /* no fixed-H found */
                    }
                    if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                    if (!pc2i->bHasDifference)
                    {
                        goto exit_function; /* nothing to do */
                    }
                    break;
                }
                else
                {
                    if (bOnly_N_V)
                    {
                        bOnly_N_V = 0;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }


        if (pc2i->nNumTgDiffMinus /*|| pc2i->nNumTgDiffH */ /* no ADP in InChI needed */)
        {
            /*--------------------------------------------------------------*/
            /*                         |                            |       */
            /* case 13: restored: O=AB=N=         original: (-)O-AB-N(+)=   */
            /* FixH:              0    0                     0        0     */
            /* MobH:              0    0                     0        0     */
            /*                        non-taut              taut   non-taut */
            /* O = O, S, Se, N; N = N, P, ...                               */
            /* t-group in original has same num. endpoints                  */
            /*       same num_H and less (-) than in the restored structure */
            /* original atom O is tautomeric, N is not taut in both         */
            /* original struct has 1 t-group; restored has less endpoints   */
            /*                             and/or possibly >1 t-groups      */
            /* Solution: separate charges between O= and -N<                */
            /*           allow moving charge to N(V) to make it N(IV)(+)    */
            /*--------------------------------------------------------------*/
            int itg;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;

            S_CHAR   *num_Fixed_H_Revrs = pStruct->pOneINChI[0]->nNum_H_fixed ? pStruct->pOneINChI[0]->nNum_H_fixed : NULL;
            S_CHAR   *pnMobHRevrs = ( pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H ) ?
                pStruct->pOneINChI[1]->nNum_H :
                ( pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H ) ?
                pStruct->pOneINChI[0]->nNum_H : NULL;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            /* djb-rwth: removing redundant code */
            /* find whether this may help */
            for (itg = 0; itg < pStruct->ti.num_t_groups && itg < pStruct->One_ti.num_t_groups; itg++)
            {
                if (pStruct->ti.t_group[itg].nNumEndpoints == pStruct->One_ti.t_group[itg].nNumEndpoints &&
                     pStruct->ti.t_group[itg].num[0] - pStruct->ti.t_group[itg].num[1] ==
                     pStruct->One_ti.t_group[itg].num[0] - pStruct->One_ti.t_group[itg].num[1] &&
                     pStruct->ti.t_group[itg].num[1] > pStruct->One_ti.t_group[itg].num[1])
                {
                    /* restored InChI t-group has more (-) and same number of H */
                    int num_SB_N_Neutr = 0, num_DB_O = 0, iat;
                    cur_success = 0;

                    for (i = 0; i < pStruct->num_atoms; i++)
                    { /* i = canonical number - 1 */
                        iat = nCanon2AtnoRevrs[i];
                        if ( /* orig. InChI info: -O(-) */
                            num_DB_O < MAX_DIFF_FIXH &&
                            (pVA[i].cNumValenceElectrons == 6 /* O, S, Se, Te */) &&
                            pStruct->endpoint[i] == itg + 1 &&
                            (e = pVA[i].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                            !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                            !(nMobHInChI && nMobHInChI[i]) &&
                            /* reversed structure info: */
                            /*!pc2i->c2at[i].endptRevrs &&*/
                            !(num_Fixed_H_Revrs && num_Fixed_H_Revrs[iat]) &&
                            !(pnMobHRevrs && pnMobHRevrs[iat]) &&
                            at2[iat].charge == 0 && at2[iat].num_H == 0 &&
                            at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2)
                        {

                            iat_DB_O[num_DB_O++] = iat;
                            /*
                            if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_function;
                            }
                            */
                        }
                        else
                        {
                            if ( /* in restored atom N: charge=0, no H, has no double bond, not an endpoint */
                                num_SB_N_Neutr < MAX_DIFF_FIXH &&
                                at2[iat].charge == 0 && !at2[iat].num_H &&
                                /*at2[iat].valence == at2[iat].chem_bonds_valence ||*/
                                (at2[iat].valence == 4 && at2[iat].chem_bonds_valence == 5) &&
                                !pVA[iat].cMetal &&
                                pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber >= 1 &&
                                !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                                /* in orig.InChI: not an endpoint, has no H */
                                !pStruct->endpoint[i] &&
                                !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                                !(nMobHInChI && nMobHInChI[i]) &&
                                /* has (+) edge */
                                (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                            {

                                iat_SB_N_Neutr[num_SB_N_Neutr++] = iat;
                                if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                            }
                        }
                    }
                    if ((num_try = inchi_min(num_SB_N_Neutr, num_DB_O))) /* djb-rwth: addressing LLVM warning */
                    {
                        /* detected; attempt to fix */
                        BNS_EDGE* pe_CMinus;
                        SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                        RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                        delta = 1;
                        for (i = 0; i < num_DB_O && cur_success < num_try; i++)
                        {
                            iat = iat_DB_O[i];
                            pe_CMinus = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                            pe_CMinus->forbidden &= forbidden_edge_mask_inv;

                            pe = pBNS->edge + pBNS->vert[iat].iedge[0]; /* double bond O=...*/
                            if (!pe->flow)
                                continue;
                            pv1 = pBNS->vert + (v1 = pe->neighbor1);
                            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                            pe->forbidden |= forbidden_edge_mask; /* change bond O=X to O(rad)-X(rad) */
                            pe->flow -= delta;
                            pv1->st_edge.flow -= delta;
                            pv2->st_edge.flow -= delta;
                            pBNS->tot_st_flow -= 2 * delta;

                            ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                            if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                                (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 2) /* djb-rwth: addressing LLVM warning */
                            {
                                /* Added (-) charge to =O and (+) charge to N => nDeltaCharge == 2 */
                                ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                                if (ret > 0)
                                {
                                    nNumRunBNS++;
                                    cur_success++; /* 13 */
                                }
                            }
                            else
                            {
                                pe->flow += delta;
                                pv1->st_edge.flow += delta;
                                pv2->st_edge.flow += delta;
                                pBNS->tot_st_flow += 2 * delta;
                            }
                            pe->forbidden &= forbidden_edge_mask_inv; /* allow changes to O=X bond */
                            INCHI_HEAPCHK
                        }
                        RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                    }
                    CurrEdges.num_edges = 0; /* clear current edge list */
                    if (cur_success)
                    {
                        tot_succes += cur_success;
                        /* recalculate InChI from the structure */
                        if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                            ppt_group_info, ppat_norm, ppat_prep)))
                        {
                            goto exit_function;
                        }
                        if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                        if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                        {
                            goto exit_function;  /* no fixed-H found */
                        }
                        if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                        if (!pc2i->bHasDifference)
                        {
                            goto exit_function; /* nothing to do */
                        }
                        break;
                    }/* else
                        if ( bOnly_N_V ) {
                        bOnly_N_V = 0;
                        }
                        */
                }
            }
        }

        if (( (pc2i->nNumTgInChI <= 1 &&
              pc2i->nNumRemHInChI > pc2i->nNumRemHRevrs) || pc2i->len_c2at ) &&
             bHas_N_V( at2, pStruct->num_atoms )) /* djb-rwth: addressing LLVM warnings */
        {
            /*-----------------------------------------------------------------*/
            /*                         |                         |             */
            /* case 14: restored:-N=AB=N=CD-XH original: (-)N-AB-N(+)=CD-XH    */
            /* FixH:              0    0   0/1              0            1     */
            /* MobH:              0    0   1/0              0            0     */
            /*                   non-taut  n/t             any  non     any    */
            /*                                                  taut           */
            /* X = O, S, Se, N; N = N                                          */
            /* t-group in original may have more (-) than in restored          */
            /*       same num_H and less (-) than in the restored structure    */
            /* atom N(V)/N(IV)(+) is not taut in both                          */
            /* The following transformation should be possible:                */
            /*        |                         |                              */
            /*   N=AB=N=CD-XH    ->     (-)N-AB-N-CD=XH(+)                     */
            /* This allows ADP to remove H(+) from -XH                         */
            /* As the result, the original structure has 0 or 1 t-group        */
            /* Solution: separate charges between -N(III)= and  N(V)           */
            /*-----------------------------------------------------------------*/
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;

            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            int num_N_V = 0, iat, i1, i2, i3, e1Flower, e1Plus, e2Plus, e2Minus, e3Plus;
            int max_success = pc2i->nNumRemHInChI - pc2i->nNumRemHRevrs;
            EDGE_LIST iat_X_List, iat_N_III_List;
            AllocEdgeList( &iat_X_List, EDGE_LIST_CLEAR );
            AllocEdgeList( &iat_N_III_List, EDGE_LIST_CLEAR );
            cur_success = 0;
            ret = 0;

            for (i = 0; i < pStruct->num_atoms; i++)
            {
                iat = nCanon2AtnoRevrs[i];
                /* search for N(V), 3 bonds */
                if ( /* restored structure */
                    num_N_V < MAX_DIFF_FIXH &&
                    at2[iat].chem_bonds_valence == 5 && at2[iat].valence == 3 &&
                    !at2[iat].charge && !at2[iat].radical &&
                    pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                    !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint) &&
                    !at2[iat].num_H &&
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pBNS->edge[e].flow /* no charge */ &&
                    NO_VERTEX != (j = GetChargeFlowerUpperEdge(pBNS, pVA, e)) && !pBNS->edge[j].forbidden &&
                    !pBNS->edge[j].flow /* neutral, valence=5 */ &&
                    /* orig. InChI */
                    !pStruct->endpoint[i] &&
                    !(nMobHInChI && nMobHInChI[i]) && pStruct->fixed_H && !pStruct->fixed_H[i]) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    iat_N_V_Array[num_N_V++] = iat;
                }
                else
                {
                    /* search for -N= */
                    if ( /* restored structure */
                        at2[iat].chem_bonds_valence == 3 && at2[iat].valence == 2 &&
                        !at2[iat].charge && !at2[iat].radical &&
                        pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                        !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint) &&
                        !at2[iat].num_H &&
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                        !pBNS->edge[e].flow /* no charge */ &&
                        /* orig. InChI */
                        /*!pStruct->endpoint[i] &&*/
                        !(nMobHInChI && nMobHInChI[i]) && pStruct->fixed_H && !pStruct->fixed_H[i]) /* djb-rwth: fixing a NULL pointer dereference */
                    {

                        if ((ret = AddToEdgeList(&iat_N_III_List, iat, 32))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_case_14;
                        }
                    }
                    else
                    {
                        /* search for -OH -NH-, -NH2 */
                        if ( /* restored structure */
                            at2[iat].chem_bonds_valence == at2[iat].valence &&
                            !at2[iat].charge && !at2[iat].radical &&
                            ((pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1) ||
                                pVA[iat].cNumValenceElectrons == 6) &&
                            at2[iat].num_H &&
                            (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                            pBNS->edge[e].flow /* no charge */ &&
                            /* orig. InChI */
                            !(nMobHInChI && nMobHInChI[i]) && pStruct->fixed_H && pStruct->fixed_H[i]) /* djb-rwth: addressing LLVM warning; fixing a NULL pointer dereference */
                        {

                            if ((ret = AddToEdgeList(&iat_X_List, iat, 32))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_case_14;
                            }
                        }
                    }
                }
            }
            if (!max_success)
            {
                max_success = inchi_min(num_N_V, iat_N_III_List.num_edges);
                max_success = inchi_min(max_success, iat_X_List.num_edges);
            }
            if (num_N_V && iat_N_III_List.num_edges && iat_X_List.num_edges)
            {
                for (i1 = 0; i1 < num_N_V && cur_success < max_success; i1++)
                {
                    int iat_N_V = iat_N_V_Array[i1];
                    if (NO_VERTEX == iat_N_V ||
                        0 >= (e1Plus = pVA[iat_N_V].nCPlusGroupEdge - 1) ||
                        NO_VERTEX == (e1Flower = GetChargeFlowerUpperEdge(pBNS, pVA, e1Plus)) ||
                        1 != pBNS->edge[e1Plus].flow ||
                        0 != pBNS->edge[e1Flower].flow)
                    {
                        continue;
                    }
                    for (i2 = iat_N_III_List.num_edges - 1; 0 <= i2 && cur_success < max_success; i2--)
                    {
                        int iat_N_III = iat_N_III_List.pnEdges[i2];
                        if (NO_VERTEX == iat_N_III ||
                            0 >= (e2Minus = pVA[iat_N_III].nCMinusGroupEdge - 1) ||
                            0 >= (e2Plus = pVA[iat_N_III].nCPlusGroupEdge - 1) ||
                            0 != pBNS->edge[e2Minus].flow ||
                            1 != pBNS->edge[e2Plus].flow)
                        {
                            /* do not consider this atom anymore */
                            iat_N_III_List.pnEdges[i2] = NO_VERTEX;
                            continue;
                        }
                        for (i3 = iat_X_List.num_edges - 1; 0 <= i3 && cur_success < max_success; i3--)
                        {
                            int iat_X = iat_X_List.pnEdges[i3];
                            BNS_VERTEX* pv1n, * pv2n;
                            BNS_EDGE* pe1n, * pe2n, * pe1Plus, * pe2Minus, * pe3Plus;
                            Vertex      v1n, v2n;
                            ret = 0;
                            if (NO_VERTEX == iat_X ||
                                0 >= (e3Plus = pVA[iat_X].nCPlusGroupEdge - 1) ||
                                1 != pBNS->edge[e3Plus].flow)
                            {
                                /* do not consider this atom anymore */
                                iat_X_List.pnEdges[i3] = NO_VERTEX;
                                continue;
                            }
                            /* all is ready to check whether the following applies:
                            forbid changes of all charges and N,P,... flowers
                            allow to change edges: e2Minus, e3Plus
                            Increment flow in e1Flower
                            The result should be: increase in number of charges by 2
                            */
                            pe1Plus = pBNS->edge + e1Plus;  /* N(V) positive charge edge */
                            pe2Minus = pBNS->edge + e2Minus; /* =N- negative charge edge */
                            pe3Plus = pBNS->edge + e3Plus;  /* -XH positive charge edge */
                            pe = pBNS->edge + e1Flower; /* N(V) flower edge */
                            pv1 = pBNS->vert + (v1 = pe->neighbor1);
                            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);
                            for (j = pv1->num_adj_edges - 1; 0 <= j; j--)
                            {
                                pe1n = pBNS->edge + pv1->iedge[j];
                                if (pe1n->flow && !pe1n->forbidden && pe1n != pe1Plus)
                                {
                                    pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                                    break;
                                }
                            }
                            if (j < 0)
                            {
                                continue; /* not found -- should not happen */
                            }
                            for (j = pv2->num_adj_edges - 1; 0 <= j; j--)
                            { /* was -2; changed 2006-2-28 12:35pm*/
                                pe2n = pBNS->edge + pv2->iedge[j];
                                if (pe2n->flow && !pe2n->forbidden && pe2n != pe1Plus)
                                {
                                    pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                                    break;
                                }
                            }
                            if (j < 0)
                            {
                                continue; /* not found -- should not happen */
                            }
                            delta = 1;
                            pe->flow += delta;
                            pe1n->flow -= delta;
                            pe2n->flow -= delta;
                            pv1n->st_edge.flow -= delta;
                            pv2n->st_edge.flow -= delta;
                            pBNS->tot_st_flow -= 2 * delta;

                            SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                            SetForbiddenEdgeMask(pBNS, &OtherNFlowerEdges, forbidden_edge_mask);

                            /* allow two charges to change */
                            pe2Minus->forbidden &= forbidden_edge_mask_inv;
                            pe3Plus->forbidden &= forbidden_edge_mask_inv;
                            /* test #1 */
                            ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);
                            INCHI_HEAPCHK
                                if (ret < 0)
                                {
                                    goto exit_case_14;
                                }
                                else
                                {
                                    if (ret == 1 && ((vPathEnd == v1n && vPathStart == v2n) ||
                                        (vPathEnd == v2n && vPathStart == v1n)) &&
                                        nDeltaCharge == 2) /* djb-rwth: addressing LLVM warning */
                                    {
                                        ; /* success */
                                    }
                                    else
                                    {
                                        ret = 0;
                                    }
                                }
                            /* restore BNS */
                            pe2Minus->forbidden |= forbidden_edge_mask;
                            pe3Plus->forbidden |= forbidden_edge_mask;
                            pe->flow -= delta;
                            pe1n->flow += delta;
                            pe2n->flow += delta;
                            pv1n->st_edge.flow += delta;
                            pv2n->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                            if (ret == 1)
                            {
                                /* test #2: check if charge separation is possible */
                                pe->flow += delta;
                                pe1n->flow -= delta;
                                pe2n->flow -= delta;
                                pv1n->st_edge.flow -= delta;
                                pv2n->st_edge.flow -= delta;
                                pBNS->tot_st_flow -= 2 * delta;

                                /* allow two charges (N(V) and N(III)) to change */
                                pe2Minus->forbidden &= forbidden_edge_mask_inv;
                                pe1Plus->forbidden &= forbidden_edge_mask_inv;
                                /* test #2 */
                                ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                    &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);
                                if (ret == 1 && ((vPathEnd == v1n && vPathStart == v2n) ||
                                    (vPathEnd == v2n && vPathStart == v1n)) &&
                                    nDeltaCharge == 2) /* djb-rwth: addressing LLVM warnings */
                                {
                                    /* success; actually change charges */
                                    ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                                    if (ret > 0)
                                    {
                                        nNumRunBNS++;
                                        cur_success++; /* 14 */
                                    }
                                }
                                if (ret <= 0)
                                {
                                    /* failed: restore BNS flow */
                                    pe->flow -= delta;
                                    pe1n->flow += delta;
                                    pe2n->flow += delta;
                                    pv1n->st_edge.flow += delta;
                                    pv2n->st_edge.flow += delta;
                                    pBNS->tot_st_flow += 2 * delta;
                                }
                                INCHI_HEAPCHK
                            }
                            RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                            RemoveForbiddenEdgeMask(pBNS, &OtherNFlowerEdges, forbidden_edge_mask);
                            if (ret > 0)
                            {
                                /* do not repeat for the same atoms */
                                iat_N_V_Array[i1] = NO_VERTEX;
                                iat_N_III_List.pnEdges[i2] = NO_VERTEX;
                                iat_X_List.pnEdges[i3] = NO_VERTEX;
                            }
                            if (ret < 0)
                            {
                                goto exit_case_14;
                            }
                            if (ret > 0)
                            {
                                break;
                            }
                        } /* i3 cycle */
                        if (ret > 0)
                        {
                            break;
                        }
                    } /* i2 cycle */
                }
            }
        exit_case_14:
            AllocEdgeList(&iat_X_List, EDGE_LIST_FREE);
            AllocEdgeList(&iat_N_III_List, EDGE_LIST_FREE);
            RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            RemoveForbiddenEdgeMask(pBNS, &OtherNFlowerEdges, forbidden_edge_mask);
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (ret < 0)
            {
                goto exit_function;
            }
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->nNumTgMRevrs > pc2i->nNumTgMInChI ||
             pc2i->nNumRemHRevrs < pc2i->nNumRemHInChI ||
             pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI ||
             (pc2i->nNumTgInChI <= 1 && pc2i->nNumTgRevrs > pc2i->nNumTgInChI)) /* djb-rwth: addressing LLVM warning */
        {
            /*--------------------------------------------------------------*/
            /* case 15: restored: -(+)O=AB-N<  orig: -O-AB=N(+)<            */
            /* (a) restored t-groups have more (-) than in original InChI   */
            /* (b) Mobile-H    charge: restored > original InChI *and*      */
            /*              removed H: restored < original InChI            */
            /* (c) restored t-groups have less endpnoits than in orig InChI */
            /* O = O, S, Se, Te; N = N                                      */
            /* Solution: move (+) from -O(+)= to -N<                        */
            /*--------------------------------------------------------------*/
            int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0;

            /* search for -O(+)= */
            /* search for -N< */
            for (i = 0; i < pStruct->num_atoms; i++)
            { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* -O(+)= in restored atom: charge=1, has no H, a double bond */
                    num_DB_Charged < MAX_DIFF_FIXH &&
                    at2[iat].charge == 1 && !at2[iat].num_H &&
                    at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    (pVA[iat].cNumValenceElectrons == 6) &&
                    /* in orig.InChI: an endpoint, has fixed-H */
                    /*pStruct->endpoint[i] &&*/
                    !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    !(nMobHInChI && nMobHInChI[i]) &&
                    /* has (+) edge */
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && 0 == pBNS->edge[e].forbidden)
                {

                    iat_DB_Charged[num_DB_Charged++] = iat;
                    /*
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                    }
                    */
                }
                else
                {
                    if ( /* -N< in restored atom: charge=0, has no H, has no double bond, N only */
                        num_SB_Neutr < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H &&
                        at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        (pVA[iat].cNumValenceElectrons == 5 &&
                            pVA[iat].cPeriodicRowNumber == 1) &&
                        /* in orig.InChI: an endpoint, has fixed-H */
                        /*pStruct->endpoint[i] &&*/
                        !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                        !(nMobHInChI && nMobHInChI[i]) &&
                        /* has (+) edge */
                        (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && 0 == pBNS->edge[e].forbidden)
                    {

                        iat_SB_Neutr[num_SB_Neutr++] = iat;
                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_Neutr, num_DB_Charged))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                BNS_VERTEX* pv1n, * pv2n;
                BNS_EDGE* pe1n, * pe2n;
                Vertex      v1n, v2n;
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_DB_Charged && cur_success < num_try; i++)
                {
                    iat = iat_DB_Charged[i];
                    pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                    if (pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    for (j = pv1->num_adj_edges - 1; 0 <= j; j--)
                    {
                        pe1n = pBNS->edge + pv1->iedge[j];
                        if (pe1n->flow && !pe1n->forbidden)
                        {
                            pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                            break;
                        }
                    }
                    if (j < 0)
                    {
                        continue; /* not found */
                    }

                    for (j = pv2->num_adj_edges - 1; 0 <= j; j--)
                    { /* was -2; changed 2006-2-28 12:35pm*/
                        pe2n = pBNS->edge + pv2->iedge[j];
                        if (pe2n->flow && !pe2n->forbidden)
                        {
                            pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                            break;
                        }
                    }
                    if (j < 0)
                    {
                        continue; /* not found */
                    }

                    pe->flow += delta;
                    pe1n->flow -= delta;
                    pe2n->flow -= delta;
                    pv1n->st_edge.flow -= delta;
                    pv2n->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1n && vPathStart == v2n) ||
                        (vPathEnd == v2n && vPathStart == v1n)) &&
                        (nDeltaCharge == 0 || nDeltaCharge == 1)) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Moved charge from O(+) to -N< => nDeltaCharge == 1 or 0 if pe2n = -N< charge edge */
                        /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 15 */
                        }
                    }
                    else
                    {
                        pe->flow -= delta;
                        pe1n->flow += delta;
                        pe2n->flow += delta;
                        pv1n->st_edge.flow += delta;
                        pv2n->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->nNumTgDiffMinus)
        {
            /*----------------------------------------------------------------*/
            /* case 16: restored: O=X-NH(-)      orig.:  O(-)-X=NH            */
            /*            t-group: (H,-)                  (2H)                */
            /* O(-) = S, Se, Te; N = N;                                       */
            /* Solution: move (-) from O(-) to -NH(-)                         */
            /*----------------------------------------------------------------*/
            int num_SB_N_Minus = 0, num_DB_O_Neutr = 0, iat, itg;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0; /* djb-rwth: ignoring LLVM warning: variable used */

            for (itg = 0; itg < pStruct->ti.num_t_groups && itg < pStruct->One_ti.num_t_groups; itg++)
            {
                if (pStruct->ti.t_group[itg].nNumEndpoints != pStruct->One_ti.t_group[itg].nNumEndpoints ||
                    pStruct->ti.t_group[itg].num[1] >= pStruct->One_ti.t_group[itg].num[1])
                {
                    continue;
                }
                CurrEdges.num_edges = num_SB_N_Minus = num_DB_O_Neutr = 0;
                cur_success = 0;
                for (j = 0, k = pStruct->One_ti.t_group[itg].nFirstEndpointAtNoPos;
                    j < pStruct->One_ti.t_group[itg].nNumEndpoints; j++)
                {
                    i = pStruct->One_ti.nEndpointAtomNumber[k + j]; /* canonical number in restored struct. */
                    iat = nCanon2AtnoRevrs[i];
                    if ( /* in restored atom: charge=0, has no H, has double bond, O, S, Se, Te */
                        num_DB_O_Neutr < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H &&
                        at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        pVA[iat].cNumValenceElectrons == 6 &&
                        /* in orig.InChI: an endpoint, may have fixed-H */
                        pStruct->endpoint[i] &&
                        /*!(pStruct->fixed_H && pStruct->fixed_H[i]) &&*/
                        !(nMobHInChI && nMobHInChI[i]) &&
                        /* has (-) edge */
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden)
                    {

                        iat_DB_O_Neutr[num_DB_O_Neutr++] = iat;

                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                    else
                    {
                        if ( /* in restored atom: charge=-1, has H, has double bond, N */
                            num_SB_N_Minus < MAX_DIFF_FIXH &&
                            at2[iat].charge == -1 && at2[iat].num_H &&
                            at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                            pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                            /* in orig.InChI: an endpoint, has no fixed-H */
                            pStruct->endpoint[i] &&
                            (pStruct->fixed_H && pStruct->fixed_H[i]) &&
                            !(nMobHInChI && nMobHInChI[i]) &&
                            /* has (-) edge */
                            (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 &&
                            0 == pBNS->edge[e].forbidden)
                        {

                            iat_SB_N_Minus[num_SB_N_Minus++] = iat;
                            /*
                            if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_function;
                            }
                            */
                        }
                    }
                }
                if ((num_try = inchi_min(num_SB_N_Minus, num_DB_O_Neutr))) /* djb-rwth: addressing LLVM warning */
                {
                    /* detected; attempt to fix */
                    SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                    RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                    /* allow stereobonds in rings change */
                    /*
                    if ( forbidden_stereo_edge_mask )
                    RemoveForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );
                    */
                    delta = 1;
                    for (i = 0; i < num_SB_N_Minus && cur_success < num_try; i++)
                    {
                        iat = iat_SB_N_Minus[i];
                        pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                        if (!pe->flow)
                            continue;
                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                        /*pe->forbidden |= forbidden_edge_mask;*/
                        pe->flow -= delta;
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;

                        ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                            &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                        if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                        {
                            /* Moved (-) charge to =O => nDeltaCharge == 1 */
                            /* Flow change on pe (-)charge edge (atom -NH(-)) is not known to RunBnsTestOnce()) */
                            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                            if (ret > 0)
                            {
                                nNumRunBNS++;
                                cur_success++; /* 16 */
                            }
                        }
                        else
                        {
                            pe->forbidden &= forbidden_edge_mask_inv;
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                        INCHI_HEAPCHK
                    }
                    RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                    /*
                    if ( forbidden_stereo_edge_mask )
                    SetForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );
                    */
                }
                CurrEdges.num_edges = 0; /* clear current edge list */
                if (cur_success)
                {
                    tot_succes += cur_success;
                    /* recalculate InChI from the structure */
                    if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                        ppt_group_info, ppat_norm, ppat_prep)))
                    {
                        goto exit_function;
                    }
                    if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                    if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                    {
                        goto exit_function;  /* no fixed-H found */
                    }
                    if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
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

        if (pc2i->nNumRemHInChI < pc2i->nNumRemHRevrs)
        {
            /*--------------------------------------------------------------*/
            /* case 17: restored: OH(+)=AB-O-     orig.  HO-AB=O(+)-        */
            /* number of removed H:  n+m                     n              */
            /* OH(+) = N, O, S, Se; -O- = P,As,O,S,Se,Te,F,Cl,Br,I          */
            /* Solution: move (+) from OH(+) to -O-                         */
            /*--------------------------------------------------------------*/
            int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H ? pInChI[1]->nNum_H :
                pInChI[0] && pInChI[0]->nNum_H ? pInChI[0]->nNum_H : 0;
            cur_success = 0;

            for (i = 0; i < pStruct->num_atoms; i++)
            { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=+1, has H, has double bond, N, O, S, Se, Te */
                    num_DB_Charged < MAX_DIFF_FIXH &&
                    at2[iat].charge == 1 && at2[iat].num_H &&
                    at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    (pVA[iat].cNumValenceElectrons == 6 ||
                        (pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1)) &&
                    /* has (+) edge */
                    (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warnings */
                {

                    iat_DB_Charged[num_DB_Charged++] = iat;
                    /*
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                    }
                    */
                }
                else
                {
                    if ( /* in restored atom: charge=0, has no H, has no double bond, N, P, O, S, Se, Te */
                        num_SB_Neutr < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H &&
                        at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        (pVA[iat].cNumValenceElectrons == 6 || pVA[iat].cNumValenceElectrons == 7 ||
                            (pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber > 1)) &&
                        /* in orig.InChI: not an endpoint */
                        !pStruct->endpoint[i] &&
                        !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                        !(nMobHInChI && nMobHInChI[i]) &&
                        /* has (+) edge */
                        (e = pVA[iat].nCPlusGroupEdge - 1) >= 0 &&
                        0 == pBNS->edge[e].forbidden) /* djb-rwth: addressing LLVM warning */
                    {

                        iat_SB_Neutr[num_SB_Neutr++] = iat;

                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_Neutr, num_DB_Charged))) /* djb-rwth: addressing LLVM warning */
            {
                BNS_VERTEX* pv1n, * pv2n;
                BNS_EDGE* pe1n, * pe2n;
                Vertex      v1n, v2n;

                num_try = inchi_min(num_try, pc2i->nNumRemHRevrs - pc2i->nNumRemHInChI);
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_DB_Charged && cur_success < num_try; i++)
                {
                    iat = iat_DB_Charged[i];
                    pe = pBNS->edge + pVA[iat].nCPlusGroupEdge - 1;
                    if (pe->flow)
                    {
                        continue;
                    }
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    for (j = pv1->num_adj_edges - 1; 0 <= j; j--)
                    {
                        pe1n = pBNS->edge + pv1->iedge[j];
                        if (pe1n->flow && !pe1n->forbidden)
                        {
                            pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                            break;
                        }
                    }
                    if (j < 0)
                    {
                        continue; /* not found */
                    }

                    for (j = pv2->num_adj_edges - 1; 0 <= j; j--)
                    { /* was -2; changed 2006-2-28 12:35pm*/
                        pe2n = pBNS->edge + pv2->iedge[j];
                        if (pe2n->flow && !pe2n->forbidden)
                        {
                            pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                            break;
                        }
                    }
                    if (j < 0)
                        continue; /* not found */

                    pe->flow += delta;
                    pe1n->flow -= delta;
                    pe2n->flow -= delta;
                    pv1n->st_edge.flow -= delta;
                    pv2n->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1n && vPathStart == v2n) ||
                        (vPathEnd == v2n && vPathStart == v1n)) &&
                        (nDeltaCharge == 0 || nDeltaCharge == 1)) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Moved charge from OH(+) to -O- => nDeltaCharge == 1 or 0 if pe2n = -O- charge edge */
                        /* Flow change on pe (+)charge edge (atom OH(+)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 17 */
                        }
                    }
                    else
                    {
                        pe->flow -= delta;
                        pe1n->flow += delta;
                        pe2n->flow += delta;
                        pv1n->st_edge.flow += delta;
                        pv2n->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (( pc2i->nNumTgInChI && pStruct->endpoint &&
              pc2i->nNumTgMInChI > pc2i->nNumTgMRevrs && pc2i->nNumEndpInChI > pc2i->nNumEndpRevrs ))
        {
            /*-----------------------------------------------------------------*/
            /*                                                                 */
            /* case 18: restored:-N=AB-X                -(-)N-AB-X(+)          */
            /* FixH:              0    0                    0    0             */
            /* MobH:              0    0                    0    0             */
            /*                   non  non                 taut  non            */
            /*                  taut  taut                      taut           */
            /* X = any heteroatom   N=N                                        */
            /* t-group in original has (Hn,-m) in the restored: (Hn,-m+1)      */
            /*       same num_H and more (-) than in the restored structure    */
            /* atom X is not taut in both                                      */
            /* Solution: separate charges between -N(III)= and  X              */
            /*-----------------------------------------------------------------*/
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            int iat, e1, itg, max_success;
            CurrEdges.num_edges = 0;
            cur_success = 0;
            ret = 0;
            /* search for -N= */
            for (itg = 0; itg < pStruct->ti.num_t_groups && itg < pStruct->One_ti.num_t_groups; itg++)
            {
                if (pStruct->ti.t_group[itg].nNumEndpoints <= pStruct->One_ti.t_group[itg].nNumEndpoints ||
                     pStruct->ti.t_group[itg].num[1] <= pStruct->One_ti.t_group[itg].num[1])
                {
                    continue;
                }
                CurrEdges.num_edges = 0;
                cur_success = 0;
                for (j = 0, k = pStruct->ti.t_group[itg].nFirstEndpointAtNoPos;
                      j < pStruct->ti.t_group[itg].nNumEndpoints; j++)
                {
                    i = pStruct->ti.nEndpointAtomNumber[k + j]; /* canonical number in restored struct. */
                    iat = nCanon2AtnoRevrs[i];
                    if (!pStruct->endpoint[i] || !at_Mobile_H_Revrs || at_Mobile_H_Revrs[iat].endpoint ||
                         pVA[i].cNumValenceElectrons != 5 || pVA[i].cPeriodicRowNumber != 1 ||
                         2 != at2[iat].valence || at2[iat].num_H || at2[iat].radical ||
                         (0 <= ( e1 = pVA[iat].nCPlusGroupEdge - 1 ) && !pBNS->edge[e1].flow) ||
                         0 >( e = pVA[iat].nCMinusGroupEdge - 1 ) || pBNS->edge[e].forbidden || pBNS->edge[e].flow) /* djb-rwth: addressing LLVM warning */
                    {
                        continue;
                    }
                    /* found -N= */
                    if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if (!( max_success = CurrEdges.num_edges ))
            {
                goto exit_case_18;
            }
            /* search for X */
            for (i = 0; i < pStruct->num_atoms && cur_success < max_success; i++)
            {
                iat = nCanon2AtnoRevrs[i];
                if (pStruct->endpoint[i] || !pVA[i].cNumValenceElectrons || pVA[i].cNumValenceElectrons == 4 ||
                     at2[iat].num_H || at2[iat].radical ||
                     (0 <= ( e1 = pVA[iat].nCMinusGroupEdge - 1 ) && !pBNS->edge[e1].flow) ||
                     0 >( e = pVA[iat].nCPlusGroupEdge - 1 ) || pBNS->edge[e].forbidden || pBNS->edge[e].flow != 1) /* djb-rwth: addressing LLVM warnings */
                {
                    continue;
                }
                /* try to move the charge */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                SetForbiddenEdgeMask( pBNS, &OtherNFlowerEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );

                pe = pBNS->edge + e;
                if (!pe->flow)
                    continue;
                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                delta = 1;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2 * delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                  (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                {
                    /* Created (-) charge on -N= => nDeltaCharge == 1 */
                    /* Flow change on pe (+)charge edge (atom X) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if (ret > 0)
                    {
                        nNumRunBNS++;
                        cur_success++; /* 18 */
                    }
                }
                else
                {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2 * delta;
                }
                INCHI_HEAPCHK
            }
        exit_case_18:
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &OtherNFlowerEdges, forbidden_edge_mask );
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (ret < 0)
            {
                goto exit_function;
            }
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 >( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at >= 1)
        {
            /*--------------------------------------------------------------*/
            /* case 19 restored:       M--OH   original:  M(-)==OH(+)       */
            /* FixH:               metal  0                      1          */
            /* MobH:                      1                      0          */
            /* O =  O, S, Se, Te; not taut. in InChI                        */
            /* In restored structure has H; tautomeric or not tautomeric    */
            /* Solution: move (+) from -OH to M; charhe on M may vary       */
            /*--------------------------------------------------------------*/
            int iat;
            EdgeIndex eMPlus, eMMinus; /* djb-rwth: removing redundant variables/code */
            BNS_EDGE  *peOHPlus, *peMPlus, *peMMinus, *peOMBond;
            int       iatMetal, ChargeOnMetal, DeltaChargeExpected;
            cur_success = 0;
            /* djb-rwth: removing redundant code */
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: =NH2(+), =OH(+) */
                    ( pc2i->c2at[i].nValElectr == 6 ) /* N, O, S, Se, Te */ &&
                     /*!pc2i->c2at[i].endptInChI &&*/ /* <=== relaxation */
                     ( e = pVA[iat].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                     pc2i->c2at[i].nFixHInChI == 1 && pc2i->c2at[i].nMobHInChI == 0 &&
                     /* reversed structure info: */
                     pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 1 &&
                     pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H &&
                     at2[iat].valence == 1 &&
                     at2[iat].valence == at2[iat].chem_bonds_valence &&
                     /* metal atom */
                     pVA[iatMetal = at2[iat].neighbor[0]].cMetal &&
                     ( eMPlus = pVA[iatMetal].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[eMPlus].forbidden &&
                     ( eMMinus = pVA[iatMetal].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[eMMinus].forbidden &&
                     !pBNS->edge[pBNS->vert[iat].iedge[0]].forbidden
                     ) /* djb-rwth: removing redundant code */
                {
                    /* -OH charge edges */
                    if ((ret = AddToEdgeList( &CurrEdges, iat, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if (CurrEdges.num_edges)
            {
                /* detected; fix */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                SetForbiddenEdgeMask( pBNS, &NFlowerEdges, forbidden_edge_mask );
                SetForbiddenEdgeMask( pBNS, &AllBondEdges, forbidden_edge_mask );
                for (i = 0; i < CurrEdges.num_edges; i++)
                {
                    /* v1 is -OH, v2 is adjacent to it Metal */
                    iat = CurrEdges.pnEdges[i];
                    iatMetal = at2[iat].neighbor[0];
                    peOHPlus = pBNS->edge + ( (long long)pVA[iat].nCPlusGroupEdge - 1 ); /* djb-rwth: removing redundant variables/code */
                    peMPlus = pBNS->edge + ( (long long)pVA[iatMetal].nCPlusGroupEdge - 1 ); /* djb-rwth: removing redundant variables/code */
                    peMMinus = pBNS->edge + ( (long long)pVA[iatMetal].nCMinusGroupEdge - 1 ); /* djb-rwth: removing redundant variables/code */
                    peOMBond = pBNS->edge + ( pBNS->vert[iat].iedge[0] ); /* djb-rwth: removing redundant variables/code */
                    /* remove forbidden edge masks */
                    peMPlus->forbidden &= forbidden_edge_mask_inv;
                    peMMinus->forbidden &= forbidden_edge_mask_inv;
                    peOMBond->forbidden &= forbidden_edge_mask_inv;

                    ChargeOnMetal = ( peMPlus->cap - peMPlus->flow ) - peMMinus->flow;
                    if (1 == ChargeOnMetal)
                    {
                        /* We are going to subtract 1 from the charge on Metal */
                        /* Added (+)charge to -OH is not known to RunBnsTestOnce() */
                        DeltaChargeExpected = -1; /* charge will become = 0 */
                    }
                    else
                    {
                        if (0 == ChargeOnMetal)
                        {
                            DeltaChargeExpected = 1; /* charge on Metal will be created */
                        }
                        else
                        {
                            DeltaChargeExpected = 0;
                        }
                    }

                    delta = 1;
                    pe = peOHPlus;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                    pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                    pe->flow -= delta; /* remove (-) from AB-O(-) */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                      (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == DeltaChargeExpected) /* djb-rwth: addressing LLVM warnings */
                    {
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 19 */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                        /* set forbidden edge masks back */
                        peMPlus->forbidden |= forbidden_edge_mask;
                    peMMinus->forbidden |= forbidden_edge_mask;
                    peOMBond->forbidden |= forbidden_edge_mask;
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &NFlowerEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &AllBondEdges, forbidden_edge_mask );

                CurrEdges.num_edges = 0; /* clear current edge list */
                if (cur_success)
                {
                    tot_succes += cur_success;
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
                    if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
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

        if (pc2i->len_c2at > 1 && pc2i->nNumTgRevrs && pc2i->nNumTgInChI)
        {
            /*--------------------------------------------------------------*/
            /* case 20: restored:  O(-)-AB=N-   original:   O=AB-N(-)-      */
            /* FixH:               0       0                0     -1        */
            /* MobH:               0       0                0      1        */
            /*                   taut    non-taut       non-taut taut       */
            /*                           or taut                  no H      */
            /*                           no H                               */
            /* O = O, S, Se; N = N, O, S, Se, Te;                           */
            /* restored atoms are taut/non-taut; original are opposite.     */
            /* Solution: move (-) from O(-) to =N-                          */
            /*--------------------------------------------------------------*/
            int num_SB_O_Minus = 0, num_DB_N = 0, iat;

            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            /*
            inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            cur_success = 0;

            CurrEdges.num_edges = 0; /* clear current edge list */
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: =O or -N= */
                    num_DB_N < MAX_DIFF_FIXH &&
                    pc2i->c2at[i].endptInChI &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pBNS->edge[e].flow == 0 &&
                    pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                    /* if  more than 1 t-group are in orig. InChI then do not move (-) to N */
                    (pc2i->nNumTgInChI == 1 || pc2i->c2at[i].nValElectr == 6) &&
                    /* reversed structure info: */
                    !pc2i->c2at[i].endptRevrs &&
                    pc2i->c2at[i].nFixHRevrs == 0 && /*pc2i->c2at[i].nMobHRevrs == 0 &&*/
                    pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                    at2[iat].valence + 1 == at2[iat].chem_bonds_valence)
                {
                    iat_DB_N[num_DB_N++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
                else
                {
                    if ( /* orig. InChI info: -O(-) */
                        num_SB_O_Minus < MAX_DIFF_FIXH &&
                        !pc2i->c2at[i].endptInChI &&
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                        pBNS->edge[e].flow == 1 &&
                        pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                        pc2i->c2at[i].nValElectr == 6 &&
                        /* reversed structure info: */
                        pc2i->c2at[i].endptRevrs &&
                        pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                        pc2i->c2at[i].nAtChargeRevrs == -1 && !at2[iat].num_H &&
                        at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 1)
                    {
                        iat_SB_O_Minus[num_SB_O_Minus++] = iat;
                    }
                }
            }
            if (!num_DB_N)
            {
                /* search among N that are tautomeric in both cases */
                for (i = 0; i < pStruct->num_atoms; i++)
                { /* i = canonical number - 1 */
                    if (!pStruct->endpoint[i])
                    {
                        continue;
                    }
                    iat = nCanon2AtnoRevrs[i];
                    if ( /* in restored atom O: charge=-1, no H, has no double bond, endpoint */
                        num_DB_N < MAX_DIFF_FIXH &&
                        at2[iat].charge == 0 && !at2[iat].num_H &&
                        at2[iat].valence + 1 == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                        /* in orig.InChI: an endpoint, has no H */
                        !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                        /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                        /* has (-) edge */
                        (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                        !pBNS->edge[e].flow)
                    {

                        iat_DB_N[num_DB_N++] = iat;
                        if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
            if ((num_try = inchi_min(num_SB_O_Minus, num_DB_N))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_O_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_O_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Added (-) charge to =N- => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom -O(-)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 20 */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at && pc2i->nNumTgRevrs && pc2i->nNumTgHInChI && pStruct->endpoint)
        {
            /*--------------------------------------------------------------*/
            /*                      O(-)                      O             */
            /*                      |                         ||            */
            /* case 21: restored: R=S=O         original:   R-S=O           */
            /*                      |                         |             */
            /*                      O(-)                      O(-)          */
            /*                           All O are taut     R is not taut   */
            /*                                                              */
            /* In addition, another atom O that should have been tautomeric */
            /* or has H(+) added in Mobile-H layer is not like that         */
            /* O = O, S, Se;  S=S, Se, Te                                  */
            /* Solution: move (-) from O(-) to =O                           */
            /*           these atoms are tautomeric in restored structure   */
            /*--------------------------------------------------------------*/
            int num_SB_O_Minus = 0, num_DB_O = 0, iat, iS;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            CurrEdges.num_edges = 0; /* clear current edge list */
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: =O    */
                    num_DB_O < MAX_DIFF_FIXH &&
                    pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                    (pc2i->c2at[i].endptInChI || pc2i->c2at[i].nMobHInChI) &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                    /* reversed structure info: */
                    !(pc2i->c2at[i].endptRevrs || pc2i->c2at[i].nMobHRevrs) &&
                    pc2i->c2at[i].nFixHRevrs == 0 &&
                    pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                    at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2)
                {
                    iat_DB_O[num_DB_O++] = iat;
                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            for (i = 0; num_DB_O && i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                if (!pStruct->endpoint[i])
                {
                    continue;
                }
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom O: charge=-1, no H, has no double bond, endpoint */
                    num_SB_O_Minus < MAX_DIFF_FIXH &&
                    at2[iat].charge == -1 && !at2[iat].num_H &&
                    at2[iat].valence == 1 && at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                    pVA[iat].cNumValenceElectrons == 6 &&
                    (at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                    /* in orig.InChI: an endpoint, has no H */
                    !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                    /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                    /* has (-) edge */
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pBNS->edge[e].flow)
                {
                    int nNumTautSB = 0, nNumTautDB = 0, nNumOtherDB = 0, nNumOtherSB = 0, nNumOthers = 0, nNumNegEndp = 0; /* djb-rwth: ignoring LLVM warning: variables used */
                    /* traverse neighbors of the centerpoint iS */
                    iS = at2[i].neighbor[0];
                    for (j = 0; j < num_SB_O_Minus; j++)
                    {
                        if (iat_Central[j] == iS)
                            break;
                    }
                    if (j < num_SB_O_Minus)
                    {
                        continue;  /* have already been there */
                    }
                    for (j = 0; j < at[iS].valence; j++)
                    {
                        int bond_type = at2[iS].bond_type[j];
                        k = at2[iS].neighbor[j];
                        if (k == i)
                        {
                            continue;
                        }
                        if (pStruct->endpoint[k] == pStruct->endpoint[i])
                        {
                            nNumTautSB += (bond_type == BOND_TYPE_SINGLE);
                            nNumTautDB += (bond_type == BOND_TYPE_DOUBLE);
                        }
                        else
                        {
                            if (bond_type == BOND_TYPE_DOUBLE)
                            {
                                nNumOtherDB++;
                            }
                            else
                            {
                                if (bond_type == BOND_TYPE_SINGLE)
                                {
                                    nNumOtherSB++;
                                }
                                else
                                {
                                    nNumOthers++;
                                }
                            }
                        }
                        if (at2[k].endpoint == at2[i].endpoint && at2[k].valence == 1 &&
                            at2[k].charge == -1 && pVA[k].cNumValenceElectrons == 6)
                        {
                            nNumNegEndp++;
                        }
                    }
                    if (!nNumTautSB)
                    {
                        continue;
                    }
                    if (!(nNumOtherDB && nNumTautDB))
                    {
                        continue; /* ignore */
                    }

                    iat_SB_O_Minus[num_SB_O_Minus] = iat;
                    iat_Central[num_SB_O_Minus++] = iS;
                }
            }
            if ((num_try = inchi_min(num_SB_O_Minus, num_DB_O))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_SB_O_Minus && cur_success < num_try; i++)
                {
                    iat = iat_SB_O_Minus[i];
                    pe = pBNS->edge + pVA[iat].nCMinusGroupEdge - 1;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Added (-) charge to =O => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 21 */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at && pc2i->nNumTgRevrs && pc2i->nNumEndpInChI < pc2i->nNumEndpRevrs)
        {
            /*--------------------------------------------------------------*/
            /*                      O                         O             */
            /*                      ||                        ||            */
            /* case 21a:restored: R=S-R' =X     original:   R-S-R' -X(-)    */
            /*                      |                         ||            */
            /*                      O(-)                      O(-)          */
            /*             All O and X are taut      O and X are not taut   */
            /*             it is possible that X is R                       */
            /*                                                              */
            /* O = O, S, Se;  S=S, Se, Te; X = N, O, S, Se, Te              */
            /* Solution: move (-) from O(-) to =X                           */
            /*           these atoms are tautomeric in restored structure   */
            /*--------------------------------------------------------------*/
            int iat, iS;
            /*
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            */
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            EDGE_LIST  OtherSO, CentralS, SOMinus, MinusAcceptord;
            CurrEdges.num_edges = 0; /* clear current edge list */
            AllocEdgeList( &OtherSO, EDGE_LIST_CLEAR );
            AllocEdgeList( &CentralS, EDGE_LIST_CLEAR );
            AllocEdgeList( &SOMinus, EDGE_LIST_CLEAR );
            AllocEdgeList( &MinusAcceptord, EDGE_LIST_CLEAR );
            cur_success = 0;
            if (!at_Mobile_H_Revrs)
            {
                goto exit_case_21a;
            }
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: -X(-)    */
                     /*num_DB_O < MAX_DIFF_FIXH &&*/
                     /*pc2i->c2at[i].nValElectr == 6 */ /* O, S, Se, Te */
                     !pc2i->c2at[i].endptInChI &&
                     ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                     pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                                                      /* reversed structure info: */
                     ( pc2i->c2at[i].endptRevrs || pc2i->c2at[i].nMobHRevrs ) &&
                     pc2i->c2at[i].nFixHRevrs == 0 &&
                     /*pc2i->c2at[i].nAtChargeRevrs == 0 &&*/ !at2[iat].num_H)
                {
                    if (pVA[iat].cNumValenceElectrons == 6 && at2[iat].charge == -1 &&
                         pBNS->edge[e].flow &&
                         at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 1 &&
                         pVA[iS = (int) at2[iat].neighbor[0]].cNumValenceElectrons == 6 && pVA[iS].cPeriodicRowNumber > 1 &&
                         at2[iS].valence >= 4)
                    {
                        /* a candidate for S in -SO2- */
                        int nNumTautSB = 0, nNumTautDB = 0, nNumOtherDB = 0, nNumOtherSB = 0; /* djb-rwth: ignoring LLVM warning: variables used */
                        int nNumOthers = 0, nNumNegEndp = 0, nNumEndpO = 0;
                        /* check whether we have already found it */
                        if (0 <= FindInEdgeList( &CentralS, iS ))
                        {
                            continue;
                        }
                        for (j = 0; j < at[iS].valence; j++)
                        {
                            int bond_type = at2[iS].bond_type[j];
                            k = at2[iS].neighbor[j];
                            if (k == iat)
                            {
                                continue;
                            }
                            if (pc2i->c2at[i].endptRevrs == at_Mobile_H_Revrs[k].endpoint && !at2[k].endpoint)
                            {
                                nNumTautSB += ( bond_type == BOND_TYPE_SINGLE );
                                nNumTautDB += ( bond_type == BOND_TYPE_DOUBLE );
                                nNumEndpO += ( pVA[k].cNumValenceElectrons == 6 && at2[k].valence == 1 );
                            }
                            else
                            {
                                if (bond_type == BOND_TYPE_DOUBLE)
                                {
                                    nNumOtherDB++;
                                }
                                else
                                {
                                    if (bond_type == BOND_TYPE_SINGLE)
                                    {
                                        nNumOtherSB++;
                                    }
                                    else
                                    {
                                        nNumOthers++;
                                    }
                                }
                            }
                            if (at2[k].endpoint == at2[i].endpoint && at2[k].valence == 1 &&
                                 at2[k].charge == -1 && pVA[k].cNumValenceElectrons == 6)
                            {
                                nNumNegEndp++;
                            }
                        }
                        if (!nNumEndpO)
                        {
                            continue;
                        }
                        if (nNumTautSB + nNumTautDB + nNumOtherDB <= nNumEndpO)
                        {
                            continue; /* ignore */
                        }
                        /* collect double bond taut =O */
                        for (j = 0; j < at[iS].valence; j++)
                        {
                            int bond_type = at2[iS].bond_type[j];
                            k = at2[iS].neighbor[j];
                            if (pc2i->c2at[i].endptRevrs == at_Mobile_H_Revrs[k].endpoint &&
                                 !at2[k].endpoint && pVA[k].cNumValenceElectrons == 6 && at2[k].valence == 1 &&
                                 0 <= ( e = pVA[k].nCMinusGroupEdge - 1 ) && !pBNS->edge[e].forbidden)
                            {
                                if (bond_type == BOND_TYPE_DOUBLE && !at2[k].charge && !pBNS->edge[e].flow)
                                {
                                    /* charges to be unchanged */
                                    if ((ret = AddToEdgeList( &OtherSO, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_case_21a;
                                    }
                                }
                                else
                                {
                                    if (bond_type == BOND_TYPE_SINGLE && at2[k].charge == -1 && pBNS->edge[e].flow)
                                    {
                                        /* charges to be removed */
                                        if ((ret = AddToEdgeList( &SOMinus, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                        {
                                            goto exit_case_21a;
                                        }
                                    }
                                }
                            }
                        }
                        if ((ret = AddToEdgeList( &CentralS, iS, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_case_21a;
                        }
                    }
                    else
                        if (at2[iat].charge == 0 && !pBNS->edge[e].flow &&
                             at2[iat].valence + 1 == at2[iat].chem_bonds_valence)
                        {
                            /* changeable charges */
                            if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                }
            }
            /* remove unchangeable from changeable */
            for (i = 0; i < OtherSO.num_edges; i++)
            {
                RemoveFromEdgeListByValue( &CurrEdges, OtherSO.pnEdges[i] );
            }

            if ((num_try = inchi_min( SOMinus.num_edges, CurrEdges.num_edges ))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
                delta = 1;
                for (i = 0; i < SOMinus.num_edges && cur_success < num_try; i++)
                {
                    pe = pBNS->edge + SOMinus.pnEdges[i];
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                    pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                    /*pe->forbidden |= forbidden_edge_mask;*/
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                      (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Added (-) charge to =O => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 21a */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            }
        exit_case_21a:
            CurrEdges.num_edges = 0; /* clear current edge list */
            AllocEdgeList( &OtherSO, EDGE_LIST_FREE );
            AllocEdgeList( &CentralS, EDGE_LIST_FREE );
            AllocEdgeList( &SOMinus, EDGE_LIST_FREE );
            AllocEdgeList( &MinusAcceptord, EDGE_LIST_FREE );
            if (cur_success)
            {
                tot_succes += cur_success;
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at)
        {
            /*------------------------------------------------------------------*/
            /* case 22: restored: N(-)=N(+)=C...=O orig: N#N-N=...-O(-)         */
            /*     im InChI        -O(-) may have H(+) added by Normalization   */
            /*                           or may be tautomeric                   */
            /* Solution: move (-) from N(-) to =O                               */
            /*                                                                  */
            /*------------------------------------------------------------------*/
            int num_DB_O = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            int iN2, iC;
            BNS_EDGE *peDB_O_Minus;
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            CurrEdges.num_edges = 0; /* clear current edge list */
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: =O    */
                    num_DB_O < MAX_DIFF_FIXH &&
                    pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                    (pc2i->c2at[i].endptInChI || pc2i->c2at[i].nMobHInChI) &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                    /* reversed structure info: */
                    !(pc2i->c2at[i].endptRevrs || pc2i->c2at[i].nMobHRevrs) &&
                    pc2i->c2at[i].nFixHRevrs == 0 &&
                    pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                    at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2)
                {
                    iat_DB_O[num_DB_O++] = iat;
                }
            }
            for (i = 0; num_DB_O && i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom O: charge=-1, no H, has no double bond, endpoint */
                    at2[iat].charge == -1 && !at2[iat].num_H &&
                    at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 && !pVA[iat].cMetal &&
                    pVA[iat].cNumValenceElectrons == 5 &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pBNS->edge[e].flow &&
                    !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                    pVA[iN2 = at2[iat].neighbor[0]].cNumValenceElectrons == 5 &&
                    at2[iat].bond_type[0] == BOND_TYPE_DOUBLE &&
                    at2[iN2].charge == 1 && at2[iN2].valence == 2 && at2[iN2].chem_bonds_valence == 4 &&
                    pVA[iC = at2[iN2].neighbor[at2[iN2].neighbor[0] == iN2]].cNumValenceElectrons == 4) /* djb-rwth: ignoring LLVM warning: variable used */
                {

                    if ((ret = AddToEdgeList(&CurrEdges, e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if ((num_try = inchi_min(CurrEdges.num_edges, num_DB_O))) /* djb-rwth: addressing LLVM warning */
            {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                RemoveForbiddenEdgeMask(pBNS, &CurrEdges, forbidden_edge_mask);
                delta = 1;
                for (i = 0; i < num_DB_O && cur_success < num_try; i++)
                {
                    iat = iat_DB_O[i];

                    peDB_O_Minus = pBNS->edge + ((long long)pVA[iat].nCMinusGroupEdge - 1); /* djb-rwth: cast operator added */
                    pe = pBNS->edge + pBNS->vert[iat].iedge[0];
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask;
                    peDB_O_Minus->forbidden &= forbidden_edge_mask_inv;

                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                        &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                    if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 0) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Added (-) charge to =O and removed from =N(-) => nDeltaCharge == 0 */
                        ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                        if (ret > 0)
                        {
                            nNumRunBNS++;
                            cur_success++; /* 22 */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                    INCHI_HEAPCHK
                        pe->forbidden &= forbidden_edge_mask_inv;
                    peDB_O_Minus->forbidden |= forbidden_edge_mask;
                }
                RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->len_c2at && pc2i->nNumTgInChI == 1)
        {
            /*------------------------------------------------------------------*/
            /* case 23: -NO2 are to be tautomeric but they are not AND          */
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
            int num_DB_O = 0, iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            /*
            inp_ATOM *atfMobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
            pStruct->pOne_norm_data[1]->at_fixed_bonds)?
            pStruct->pOne_norm_data[1]->at_fixed_bonds : NULL;
            */
            S_CHAR   *num_Fixed_H_Revrs = pStruct->pOneINChI[0]->nNum_H_fixed ? pStruct->pOneINChI[0]->nNum_H_fixed : NULL;
            S_CHAR   *pnMobHRevrs = ( pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H ) ?
                pStruct->pOneINChI[1]->nNum_H :
                ( pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H ) ?
                pStruct->pOneINChI[0]->nNum_H : NULL;
            int iN, one_success;
            BNS_EDGE *peDB_O_Minus;
            int neigh, nNumO, nNumOthers;
    #define CHG_SET_NOOH         0
    #define CHG_SET_WRONG_TAUT   1
    #define CHG_SET_TAUT         2
    #define CHG_LAST_SET         2 /* the last index in trying */
    #define CHG_SET_O_FIXED      3
    #define CHG_SET_NUM          4
            EDGE_LIST ChangeableEdges[CHG_SET_NUM];
            memset( ChangeableEdges, 0, sizeof( ChangeableEdges ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            CurrEdges.num_edges = 0; /* clear current edge list */
            cur_success = 0;

            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: taut in orig. InChI =O located in -NO2 that is not taut in Reconstructed InChI */
                    num_DB_O < MAX_DIFF_FIXH &&
                    pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                    (pc2i->c2at[i].endptInChI /*|| pc2i->c2at[i].nMobHInChI*/) &&
                    (e = pVA[iat].nCMinusGroupEdge - 1) >= 0 && !pBNS->edge[e].forbidden &&
                    pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                    /* reversed structure info: */
                    !(pc2i->c2at[i].endptRevrs /*|| pc2i->c2at[i].nMobHRevrs*/) &&
                    pc2i->c2at[i].nFixHRevrs == 0 &&
                    pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                    at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 &&
                    /* find whether it belongs to NO2 */
                    pVA[iN = at2[iat].neighbor[0]].cNumValenceElectrons == 5 &&
                    at2[iN].valence == 3 && (at2[iN].charge == 0 || at2[iN].charge == 1) &&
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
                            pStruct->endpoint[neigh] &&
                            !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[neigh].endpoint) &&
                            at2[neigh].valence == 1 && at2[neigh].num_H == 0 &&
                            at2[neigh].radical == 0 && (at2[neigh].charge == 0 || at2[neigh].charge == -1) &&
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
                    /* save the edge to avoid interference */
                    if ((ret = AddToEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_case_23;
                    }
                }
            }
            if (num_DB_O)
            {
                /* 1. search for =N(=O)-OH; assume =N(+)(-O(-))(-OH) does not happen */
                for (i = 0; i < pStruct->num_atoms; i++)
                {
                    /* i = canonical number - 1 */
                    /* find O=N(V) */
                    iat = nCanon2AtnoRevrs[i];
                    if (!pStruct->endpoint[i] || pVA[i].cNumValenceElectrons != 6 ||
                        at2[iat].valence != 1 || at2[iat].charge ||
                        0 > (e = pVA[iat].nCMinusGroupEdge - 1) ||
                        at2[iat].num_H + at2[iat].chem_bonds_valence != 2 ||
                        pVA[iN = at2[iat].neighbor[0]].cNumValenceElectrons != 5 ||
                        0 > (e = pVA[iN].nCPlusGroupEdge - 1) ||
                        pBNS->edge[e].forbidden || !pBNS->edge[e].flow ||
                        at2[iN].charge || at2[iN].valence != 3 || at2[iN].chem_bonds_valence != 5) /* djb-rwth: ignoring LLVM warning: variable used */
                    {
                        continue;
                    }
                    /* find the second O, -OH */
                    nNumO = nNumOthers = 0;
                    for (k = 0; k < at2[iN].valence; k++)
                    {
                        neigh = at2[iN].neighbor[k];
                        if (neigh == iat)
                        {
                            continue;
                        }
                        if (pVA[neigh].cNumValenceElectrons == 6 &&
                            pStruct->endpoint[neigh] &&
                            at2[neigh].valence == 1 && at2[neigh].num_H == 1 &&
                            at2[neigh].radical == 0 && (at2[neigh].charge == 0))
                        {
                            nNumO++;
                        }
                        else
                            if (at2[iN].bond_type[k] == BOND_TYPE_DOUBLE &&
                                at2[neigh].valence >= 2 &&
                                at2[neigh].valence < at2[neigh].chem_bonds_valence)
                            {
                                nNumOthers++;
                            }
                    }
                    if (nNumO != 1 || nNumOthers != 1)
                    {
                        continue;
                    }
                    /* save edges to be changed */
                    if ((ret = AddToEdgeList(&ChangeableEdges[CHG_SET_NOOH], e, INC_ADD_EDGE)) ||
                        (ret = AddToEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE)))
                    {
                        goto exit_case_23;
                    }
                    if (NO_VERTEX != (j = GetChargeFlowerUpperEdge(pBNS, pVA, e)) &&
                        ((ret = AddToEdgeList(&ChangeableEdges[CHG_SET_NOOH], j, INC_ADD_EDGE)) ||
                            (ret = AddToEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE))))
                    {
                        goto exit_case_23;
                    }
                }
                /* 2. search for (-) atoms that are tautomeric but should not be  */
                /*           or that got H from Normalization but they shouldn't  */
                for (i = 0; i < pStruct->num_atoms; i++)
                { /* i = canonical number - 1 */
                    iat = nCanon2AtnoRevrs[i];
                    if (at2[iat].charge == -1 &&
                        !pStruct->endpoint[i] &&
                        (at_Mobile_H_Revrs &&
                            (at_Mobile_H_Revrs[i].endpoint || at2[iat].num_H < at_Mobile_H_Revrs[i].num_H)))
                    {

                        if (0 <= (e = pVA[iat].nCMinusGroupEdge - 1) &&
                            0 > FindInEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e) &&
                            !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                            (
                                (ret = AddToEdgeList(&ChangeableEdges[CHG_SET_WRONG_TAUT], e, INC_ADD_EDGE)) ||
                                (ret = AddToEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE))
                                ))
                        {
                            goto exit_case_23;
                        }
                    }
                    else
                    {
                        /* negatively charged atom in Reconstructed structure got H(+) from Normalization */
                        /* and is not tautomeric; in the original structure it is tautomeric */
                        if (at2[iat].charge == -1 &&
                            pStruct->endpoint[i] &&
                            !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint) &&
                            (num_Fixed_H_Revrs && num_Fixed_H_Revrs[i] == -1) &&
                            (pnMobHRevrs && pnMobHRevrs[i] == 1) &&
                            pStruct->fixed_H[i] == 0)
                        {

                            if (0 <= (e = pVA[iat].nCMinusGroupEdge - 1) &&
                                0 > FindInEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e) &&
                                !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                                (
                                    (ret = AddToEdgeList(&ChangeableEdges[CHG_SET_WRONG_TAUT], e, INC_ADD_EDGE)) ||
                                    (ret = AddToEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE))
                                    ))
                            {
                                goto exit_case_23;
                            }
                        }
                    }
                }
                /* 3. Search for (-) atoms that are tautomeric */
                for (i = 0; i < pStruct->num_atoms; i++)
                {
                    /* i = canonical number - 1 */
                    iat = nCanon2AtnoRevrs[i];
                    if (pStruct->endpoint[i] &&
                        (at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint) &&
                        at2[iat].charge == -1
                        /*&& pVA[i].cNumValenceElectrons == 6*/)
                    {
                        if (0 <= (e = pVA[iat].nCMinusGroupEdge - 1) &&
                            !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                            0 > FindInEdgeList(&ChangeableEdges[CHG_SET_O_FIXED], e) &&
                            (ret = AddToEdgeList(&ChangeableEdges[CHG_SET_TAUT], e, INC_ADD_EDGE)))
                        {
                            goto exit_case_23;
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
                    peDB_O_Minus = pBNS->edge + ((long long)pVA[iat].nCMinusGroupEdge - 1);/* djb-rwth: cast operator added */
                    pe = pBNS->edge + pBNS->vert[iat].iedge[0];

                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

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
                        nDeltaChargeExpected = (k == CHG_SET_NOOH) ? 2 : 0;

                        SetForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
                        RemoveForbiddenEdgeMask(pBNS, &ChangeableEdges[k], forbidden_edge_mask);
                        /* allow (-) charge to move to N=O */
                        peDB_O_Minus->forbidden &= forbidden_edge_mask_inv;

                        ret = RunBnsTestOnce(pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                            &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms);

                        if (ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1)) &&
                            nDeltaCharge == nDeltaChargeExpected) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* Move (-) charge to =O and remove it an endpoint => nDeltaCharge == 0 */
                            ret = RunBnsRestoreOnce(pBNS, pBD, pVA, pTCGroups);
                            if (ret > 0)
                            {
                                nNumRunBNS++;
                                one_success++; /* 23 */
                            }
                        }
                        INCHI_HEAPCHK
                    }
                    cur_success += one_success;

                    RemoveForbiddenEdgeMask(pBNS, &AllChargeEdges, forbidden_edge_mask);
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
        exit_case_23:
            for (i = 0; i < CHG_SET_NUM; i++)
            {
                AllocEdgeList(&ChangeableEdges[i], EDGE_LIST_FREE);
            }

            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 > (ret = MakeOneInChIOutOfStrFromINChI2(pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep)))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr(pStruct)))/* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
                {
                    goto exit_function;  /* no fixed-H found */
                }
                if ((ret = FillOutCMP2FHINCHI(pStruct, at2, pVA, pInChI, pc2i))) /* djb-rwth: addressing LLVM warning */
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

        if (pc2i->len_c2at && pc2i->nNumTgInChI == 1)
        {
            /*------------------------------------------------------------------*/
            /* case 24: InChI norm. -N(-)-N(+)(IV) => -N=N(V) prevents tauto-   */
            /*          merism on -N(-)- in case of ADP                         */
            /*                                                                  */
            /* Solution: convert       N(V)=N-   ...=X    -> N(IV)(+)-N=...-X(-)*/
            /*                     N(IV)(+)-N(-)-...=X                          */
            /*                                                                  */
            /*      Orig InChI            taut      taut, 1 t-group only(ADP?)  */
            /*   Reconstructed struct   non-taut    possibly not taut           */
            /*                                                                  */
            /*   Details: 1a. store next to N(V) (+)edge its flower edge        */
            /*            1b. store next to N(-) edge NO_VERTEX                 */
            /*            2.  Release (-) edges of other missing endpoints or   */
            /*                all endpoints if no other is missing              */
            /*            3.  Decrement flow on (+) edge                        */
            /*                if flower edge is stored then expect DeltaCharge=2*/
            /*                otherwise DeltaCharge = 0                         */
            /*------------------------------------------------------------------*/
            int iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            inp_ATOM *atf = ( pStruct->pOne_norm_data[1] && pStruct->pOne_norm_data[1]->at_fixed_bonds ) ?
                pStruct->pOne_norm_data[1]->at_fixed_bonds : NULL;
            int iN, one_success;
            EdgeIndex  ef, e1;
            BNS_EDGE *pef;
    #define CHG_SET_MISSED_TAUT   0
    #define CHG_SET_OTHER_TAUT_O  1
    #define CHG_SET_OTHER_TAUT_N  2
    #define CHG_LAST_SET          2 /* the last index in trying */
    #define CHG_SET_NN            3
    #define CHG_SET_AVOID         4
    #define CHG_SET_NUM           5
            EDGE_LIST ChangeableEdges[CHG_SET_NUM];
            memset( ChangeableEdges, 0, sizeof( ChangeableEdges ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            CurrEdges.num_edges = 0; /* clear current edge list */
            cur_success = 0;
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: -N=N(V)    */
                     pc2i->c2at[i].nValElectr == 5 /* N or P */ &&
                     ( pc2i->c2at[i].endptInChI /* only N */ ) &&
                     ( e1 = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e1].forbidden &&
                     pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                     /* reversed structure info: */
                     !pc2i->c2at[i].endptRevrs &&
                     pc2i->c2at[i].nFixHRevrs == 0 &&
                     pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                     at2[iat].valence == 2 && at2[iat].chem_bonds_valence == 3 &&
                     /* find whether -N= has =N(V) neighbor; Note: operator comma: (A,B) returns B */
                     ( iN = at2[iat].neighbor[at2[iat].bond_type[0] != BOND_TYPE_DOUBLE],
                       pVA[iN].cNumValenceElectrons == 5 ) &&
                     at2[iN].chem_bonds_valence == 5 &&
                     at2[iN].charge == 0 && !at2[iN].num_H && !at2[iN].radical &&
                     0 <= ( e = pVA[iN].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                     0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ))
                {

                    ef = GetChargeFlowerUpperEdge( pBNS, pVA, e ); /* == NO_VERTEX if N(V) has 4 bonds */
                    if (( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], e, INC_ADD_EDGE ) ) ||
                        ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], ef, INC_ADD_EDGE ) ) ||
                         ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], 1, INC_ADD_EDGE ) ) || /* expected nDeltaCharge */
                         ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e1, INC_ADD_EDGE ) ) ||
                         ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ) ) ||
                         ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], ef, INC_ADD_EDGE ) ))
                    {
                        goto exit_case_24;
                    }
                    /* mark -N= so that (-) will not be moved to it */
                    if (0 <= ( e = pVA[iat].nCMinusGroupEdge ) && !pBNS->edge[e].forbidden &&
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) &&
                         ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ) ))
                    {
                        goto exit_case_24;
                    }
                }
                else
                {
                    if ( /* orig. InChI info: -N(-)N(IV)(+)    */
                         atf &&
                         pc2i->c2at[i].nValElectr == 5 /* N or P */ &&
                         pc2i->c2at[i].endptInChI /* only N */ &&
                         ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                         pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                         /* reversed structure info: */
                         !pc2i->c2at[i].endptRevrs &&
                         pc2i->c2at[i].nFixHRevrs == 0 &&
                         pc2i->c2at[i].nAtChargeRevrs == -1 && !at2[iat].num_H &&
                         at2[iat].valence == 2 && at2[iat].chem_bonds_valence == 2 &&
                         atf[iat].valence == 2 && atf[iat].chem_bonds_valence == 3 &&
                         /* find whether -N= has =N(V) neighbor; Note: operator comma: (A,B) returns B */
                         ( iN = atf[iat].neighbor[atf[iat].bond_type[0] != BOND_TYPE_DOUBLE],
                           pVA[iN].cNumValenceElectrons == 5 ) &&
                         at2[iN].charge == 1 && /* double bond neighbor */
                         at2[iN].chem_bonds_valence == 4 &&
                         atf[iN].charge == 0 &&
                         atf[iN].chem_bonds_valence == 5 &&  /* InChI normalization created N(V)=N- out of N(IV)(+)-N(-)- */
                         !at2[iN].num_H && !at2[iN].radical &&
                         0 <= ( e = pVA[iat].nCMinusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ))
                    {
                        /* save (-) edge */
                        if (( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], e, INC_ADD_EDGE ) ) ||
                            ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], NO_VERTEX, INC_ADD_EDGE ) ) ||
                             ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], 1, INC_ADD_EDGE ) ) || /* expected nDeltaCharge */
                             ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ) ))
                        {
                            goto exit_case_24;
                        }
                    }
                }
            }
            if (!ChangeableEdges[CHG_SET_NN].num_edges)
            {
                goto  exit_case_24;
            }
            /* Collect all relevant tautomeric atoms */
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                /* i = canonical number - 1 */
                if (!pStruct->endpoint[i])
                {
                    continue;
                }
                iat = nCanon2AtnoRevrs[i];
                if (at2[iat].charge || at2[iat].radical || at2[iat].valence == at2[iat].chem_bonds_valence)
                {
                    continue; /* cannot be an acceptor of (-) */
                }
                if (0 > ( e = pVA[iat].nCMinusGroupEdge - 1 ) || pBNS->edge[e].forbidden || pBNS->edge[e].flow)
                {
                    continue;
                }
                if (0 <= FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ))
                {
                    continue; /* has already been used */
                }
                /* missing endpoint */
                if (!( at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint ))
                {
                    if (0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) && (
                        ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_MISSED_TAUT], e, INC_ADD_EDGE ) ) ||
                        ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ) ) ))
                    {
                        goto exit_case_24;
                    }
                }
                else
                {
                    /* endpoint O */
                    if (pVA[iat].cNumValenceElectrons == 6)
                    {
                        if (0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) && (
                            ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_O], e, INC_ADD_EDGE ) ) ||
                            ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ) ) ))
                        {
                            goto exit_case_24;
                        }
                    }
                    else
                    {
                        /* endpoint N */
                        if (pVA[iat].cNumValenceElectrons == 5)
                        {
                            if (0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) && (
                                ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_N], e, INC_ADD_EDGE ) ) ||
                                ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ) ) ))
                            {
                                goto exit_case_24;
                            }
                        }
                    }
                }
            }
            /* ------- finally, try to move charges from -N(-)-N(+) or to N(V) --------------*/
            for (i = 0; i < ChangeableEdges[CHG_SET_NN].num_edges; i += 3)
            {
                int nDeltaChargeExpected;
                one_success = 0;
                delta = 1;
                pe = pBNS->edge + ChangeableEdges[CHG_SET_NN].pnEdges[i];
                pef = ( NO_VERTEX != ChangeableEdges[CHG_SET_NN].pnEdges[i + 1] ) ?
                    pBNS->edge + ChangeableEdges[CHG_SET_NN].pnEdges[i + 1] : NULL;
                nDeltaChargeExpected = ChangeableEdges[CHG_SET_NN].pnEdges[i + 2];

                if (!pe->flow)
                    continue;
                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

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
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &ChangeableEdges[k], forbidden_edge_mask );
                    /* allow change of N(V) flower edge */
                    if (pef)
                    {
                        pef->forbidden &= forbidden_edge_mask_inv;
                    }

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
                            nNumRunBNS++;
                            one_success++; /* 24 */
                        }
                    }
                    INCHI_HEAPCHK
                }
                cur_success += one_success;

                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );

                if (!one_success)
                {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2 * delta;
                }
            }
        exit_case_24:
            for (i = 0; i < CHG_SET_NUM; i++)
            {
                AllocEdgeList( &ChangeableEdges[i], EDGE_LIST_FREE );
            }

            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 >( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
    #undef CHG_SET_NN
    #undef CHG_SET_MISSED_TAUT
    #undef CHG_SET_OTHER_TAUT_O
    #undef CHG_SET_OTHER_TAUT_N
    #undef CHG_LAST_SET
    #undef CHG_SET_AVOID
    #undef CHG_SET_NUM
        }

        /* pStruct->nNumRemovedProtonsMobHInChI == pc2i->nNumRemHInChI */

        if (pc2i->len_c2at && pc2i->nNumTgInChI == 1 &&
             pc2i->nNumRemHRevrs > pc2i->nNumRemHInChI && 0 > pc2i->nNumRemHInChI &&
             ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI ||
               pc2i->nNumTgRevrs > pc2i->nNumTgInChI ))
        {
            /*------------------------------------------------------------------*/
            /* case 25: Restored InChI does not have 2 or more added protons    */
            /*                         possibly taut. endpoints are missing     */
            /*                         has -N(-O(-))-O(-) group(s)              */
            /*          Original InChI has only one t-group                     */
            /*                                                                  */
            /* Solution: convert       -N(-O(-))-O(-) -> -N(+)(=O)-O(-)         */
            /*                         and direct 2(-) to the missing taut atoms*/
            /*           at first attempt try to move (-) to N only             */
            /*                                                                  */
            /*------------------------------------------------------------------*/
            int iat;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            AT_NUMB  *nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[1] &&
                                            pStruct->pOne_norm_data[1]->at ) ? pStruct->pOne_norm_data[1]->at : NULL;
            /*
            inp_ATOM *atf  = (pStruct->pOne_norm_data[1] && pStruct->pOne_norm_data[1]->at_fixed_bonds)?
            pStruct->pOne_norm_data[1]->at_fixed_bonds : NULL;
            */
            int iN, neigh, one_success;
            EdgeIndex  e1, bFirst;
            BNS_EDGE *pef;
    #define CHG_SET_MISSED_TAUT_1   0
    #define CHG_SET_MISSED_TAUT_ALL 1
    #define CHG_SET_OTHER_TAUT_1    2
    #define CHG_SET_OTHER_TAUT_ALL  3
    #define CHG_LAST_SET            3 /* the last index in trying */
    #define CHG_SET_NO_IN_NO2M2     4
    #define CHG_SET_AVOID           5
    #define CHG_SET_NUM             6
            EDGE_LIST ChangeableEdges[CHG_SET_NUM];
            memset( ChangeableEdges, 0, sizeof( ChangeableEdges ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
            pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            CurrEdges.num_edges = 0; /* clear current edge list */
            cur_success = 0;
            /* find all -N(-O(-))-O(-) */
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                iat = nCanon2AtnoRevrs[i];
                if (pStruct->endpoint[i])
                {
                    if (0 >( e = pVA[iat].nCMinusGroupEdge - 1 ) || pBNS->edge[e].forbidden ||
                         0 <= FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ))
                    {
                        continue;
                    }
                    bFirst = ( (pVA[iat].cNumValenceElectrons == 5 && pc2i->nNumTgInChI == 1) ||
                               (pVA[iat].cNumValenceElectrons == 6 && pc2i->nNumTgInChI != 1) ); /* djb-rwth: addressing LLVM warnings */
                    /* many or no t-groups -> try O only first */
                    /* single t-group -> try only N first */
                    if (!( at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint ))
                    {
                        /* missed tautomeric endpoint */
                        if (bFirst &&
                            ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_MISSED_TAUT_1], e, INC_ADD_EDGE ) ))
                        {
                            goto exit_case_25;
                        }
                        if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_MISSED_TAUT_ALL], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_case_25;
                        }
                    }
                    if (bFirst &&
                        ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_1], e, INC_ADD_EDGE ) ))
                    {
                        goto exit_case_25;
                    }
                    if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_ALL], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_case_25;
                    }
                    if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_case_25;
                    }
                }
                else
                {
                    if (at2[iat].valence == 1 && at2[iat].charge == -1 &&
                         pVA[iat].cNumValenceElectrons == 6 &&
                         pVA[iN = at2[iat].neighbor[0]].cNumValenceElectrons == 5 && /* -O(-) */
                         !pStruct->endpoint[nAtno2CanonRevrs[iN]] &&
                         at2[iN].valence == 3 && at2[iN].chem_bonds_valence == 3 &&
                         !at2[iN].charge && !at2[iN].radical &&
                         0 <= ( e = pVA[iN].nCPlusGroupEdge - 1 ) && !pBNS->edge[e].forbidden &&
                         pBNS->edge[e].flow && /* NPlus edge */
                         0 <= ( e1 = pVA[iat].nCMinusGroupEdge - 1 ) && !pBNS->edge[e1].forbidden &&
                         pBNS->edge[e1].flow &&  /* OMinus edge */
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) &&
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e1 ))
                    {
                        /* found >N-O(-) */
                        int nNumO = 0, nNumOthers = 0;
                        for (k = 0; k < at2[iN].valence; k++)
                        {
                            neigh = at2[iN].neighbor[k];
                            if (neigh == iat)
                            {
                                continue;
                            }
                            if (pVA[neigh].cNumValenceElectrons == 6 &&
                                 !pStruct->endpoint[neigh] &&
                                 at2[neigh].valence == 1 && at2[neigh].num_H == 0 &&
                                 at2[neigh].radical == 0 && at2[neigh].charge == -1 &&
                                 at2[neigh].chem_bonds_valence == 1)
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
                        if (nNumO != 1 && nNumOthers != 1)
                        {
                            continue;
                        }
                        /* save charge edges: NPlus first, OMinus second */
                        if (( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NO_IN_NO2M2], e, INC_ADD_EDGE ) ) ||
                            ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NO_IN_NO2M2], e1, INC_ADD_EDGE ) ) ||
                             ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ) ) ||
                             ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e1, INC_ADD_EDGE ) ))
                        {
                            goto exit_case_25;
                        }
                    }
                }
            }
            if (!ChangeableEdges[CHG_SET_NO_IN_NO2M2].num_edges ||
                 !ChangeableEdges[CHG_SET_OTHER_TAUT_ALL].num_edges)
            {
                goto exit_case_25;
            }
            /* ------- finally, try to move charges from -NO2(2-) or to tautomeric endpoints ----*/
            for (i = 0; i < ChangeableEdges[CHG_SET_NO_IN_NO2M2].num_edges; i += 2)
            {
                int nDeltaChargeExpected = 3;
                /* change flow on O(-) to make it neutral; 3 new charges will be created:
                N(+), and two (-) on InChI endpoints
                alternatively, if we change flow on N to make N(+) then O(-) will
                be nutralized (-1 charge) and two (-) charges on taut. endpoints will be
                created (+2); the total change in this case would be (-1)+(+2) = +1
                */
                one_success = 0;
                delta = 1;
                pe = pBNS->edge + ChangeableEdges[CHG_SET_NO_IN_NO2M2].pnEdges[i + 1]; /* O(-) edge */
                pef = pBNS->edge + ChangeableEdges[CHG_SET_NO_IN_NO2M2].pnEdges[i]; /* >N- (+) edge */

                if (!pe->flow)
                    continue;
                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

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
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &ChangeableEdges[k], forbidden_edge_mask );
                    /* allow change of N(V) flower edge */
                    pef->forbidden &= forbidden_edge_mask_inv;

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
                            nNumRunBNS++;
                            one_success++; /* 24 */
                        }
                    }
                    INCHI_HEAPCHK
                }
                cur_success += one_success;

                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );

                if (!one_success)
                {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2 * delta;
                }
            }
        exit_case_25:
            for (i = 0; i < CHG_SET_NUM; i++)
            {
                AllocEdgeList( &ChangeableEdges[i], EDGE_LIST_FREE );
            }

            CurrEdges.num_edges = 0; /* clear current edge list */
            if (cur_success)
            {
                tot_succes += cur_success;
                /* recalculate InChI from the structure */
                if (0 >( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
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
                if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
    #undef CHG_SET_NN
    #undef CHG_SET_MISSED_TAUT
    #undef CHG_SET_OTHER_TAUT_O
    #undef CHG_SET_OTHER_TAUT_N
    #undef CHG_LAST_SET
    #undef CHG_SET_AVOID
    #undef CHG_SET_NUM
        }


    exit_function:

        AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
        AllocEdgeList( &CurrEdges, EDGE_LIST_FREE );
        AllocEdgeList( &NFlowerEdges, EDGE_LIST_FREE );
        AllocEdgeList( &SFlowerEdges, EDGE_LIST_FREE );
        AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_FREE );
        AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_FREE );
        AllocEdgeList( &AllBondEdges, EDGE_LIST_FREE );

        return ret < 0 ? ret : ( tot_succes && pc2i->bHasDifference );
    }
        */
    // END INCHI C FUNCTION: FixFixedHRestoredStructure
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FixFixedHRestoredStructure
    // INCHI✔️❌: READ_INCHI_STRING=1; COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // INCHI✔️❌: The complete 25-case body and function-level cleanup are active.
    // END INCHI ACTIVE MACRO CONFIGURATION: FixFixedHRestoredStructure

    let mut ret = 0_i32;
    let mut tot_success = 0_i32;
    let mut n_num_run_bns = 0_i32;
    let mut allowed_nitrogen_flower_edges = false;
    let forbidden_edge_mask_inv = !forbidden_edge_mask;
    let mut comparison = CMP2FHINCHI::default();

    let mut all_charge_edges = EDGE_LIST::default();
    let mut current_edges = EDGE_LIST::default();
    let mut sulfur_flower_edges = EDGE_LIST::default();
    let mut nitrogen_flower_edges = EDGE_LIST::default();
    let mut other_nitrogen_flower_edges = EDGE_LIST::default();
    let mut fixed_large_ring_stereo_edges = EDGE_LIST::default();
    let mut all_bond_edges = EDGE_LIST::default();

    let _ = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut current_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut nitrogen_flower_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut sulfur_flower_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut other_nitrogen_flower_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut fixed_large_ring_stereo_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut all_bond_edges, EDGE_LIST_CLEAR)?;

    let execution = (|| -> Result<(), SourceHeapError> {
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

        let mut atom_number = 0_i32;
        while atom_number < pStruct.num_atoms {
            let atom_index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
            if minus_edge >= 0 {
                let edge = heap
                    .slice(pBNS.edge.as_const())?
                    .get(
                        usize::try_from(minus_edge)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if edge.forbidden == 0 {
                    ret = AddToEdgeList(heap, &mut all_charge_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
            }

            let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
            let mut continue_outer = false;
            if plus_edge >= 0 {
                let edge_forbidden = heap
                    .slice(pBNS.edge.as_const())?
                    .get(
                        usize::try_from(plus_edge)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .forbidden;
                if edge_forbidden == 0 {
                    ret = AddToEdgeList(heap, &mut all_charge_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                    let valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if valence.cNumValenceElectrons == 5 && valence.cMetal == 0 {
                        let upper = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                        if upper != NO_VERTEX {
                            let flower = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(upper)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if flower.forbidden != 0 {
                                continue_outer = true;
                            } else if flower.flow != 0 {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut all_charge_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                                ret = AddToEdgeList(
                                    heap,
                                    &mut nitrogen_flower_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            } else {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut other_nitrogen_flower_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    } else if valence.cNumValenceElectrons == 6 && valence.cMetal == 0 {
                        let upper = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                        if upper != NO_VERTEX {
                            let flower = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(upper)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if flower.forbidden != 0 {
                                continue_outer = true;
                            } else if flower.flow != 0 {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut sulfur_flower_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    }
                }
            }
            if continue_outer {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }

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
            let mut neighbor_order = 0_i32;
            while neighbor_order < i32::from(atom.valence) {
                let order = usize::try_from(neighbor_order)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let neighbor = i32::from(atom.neighbor[order]);
                if neighbor < atom_number {
                    let edge_number = *heap
                        .slice(vertex.iedge.as_const())?
                        .get(order)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let forbidden = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden;
                    if forbidden == 0 {
                        ret = AddToEdgeList(heap, &mut all_bond_edges, edge_number, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                }
                neighbor_order = neighbor_order.wrapping_add(1);
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
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
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
                let mut neighbor_order = 0_i32;
                while neighbor_order < i32::from(atom.valence) {
                    let order = usize::try_from(neighbor_order)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = *heap
                        .slice(vertex.iedge.as_const())?
                        .get(order)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let forbidden = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden;
                    if i32::from(forbidden) == forbidden_stereo_edge_mask {
                        let ring_size = is_bond_in_Nmax_memb_ring(
                            heap,
                            at2,
                            atom_number,
                            neighbor_order,
                            bfs.q,
                            bfs.nAtomLevel,
                            bfs.cSource,
                            99,
                        )?;
                        if ring_size > 0 {
                            ret = AddToEdgeList(
                                heap,
                                &mut fixed_large_ring_stereo_edges,
                                edge_number,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                return Ok(());
                            }
                        }
                    }
                    neighbor_order = neighbor_order.wrapping_add(1);
                }
                atom_number = atom_number.wrapping_add(1);
            }
        }

        ret = FillOutExtraFixedHDataRestr(heap, pStruct)?;
        if ret != 0 {
            return Ok(());
        }
        let at2_snapshot = heap.slice(at2.as_const())?.to_vec();
        ret = FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
        if ret != 0 {
            return Ok(());
        }
        if comparison.bHasDifference == 0
            || (comparison.len_c2at == 0
                && comparison.nNumTgRevrs == comparison.nNumTgInChI
                && comparison.nNumEndpRevrs == comparison.nNumRemHInChI
                && comparison.nNumEndpRevrs == comparison.nNumEndpInChI
                && comparison.nNumTgDiffMinus == 0
                && comparison.nNumTgDiffH == 0)
        {
            return Ok(());
        }

        // Complete translation of source case 01, ichirvr3.c:514-646.
        if comparison.len_c2at >= 2 {
            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut single_bond_oxygen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_oxygen = 0_i32;
            let mut number_single_bond_oxygen_minus = 0_i32;
            let mut current_success = 0_i32;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let edge_number = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                let eligible_edge = edge_number >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if difference.nValElectr == 6 && eligible_edge {
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if number_single_bond_oxygen_minus < MAX_DIFF_FIXH as i32
                        && difference.nFixHInChI == 0
                        && difference.nMobHInChI == 0
                        && difference.nFixHRevrs == -1
                        && difference.nMobHRevrs == 1
                        && difference.nAtChargeRevrs == -1
                        && atom.num_H == 0
                        && atom.valence == 1
                        && atom.chem_bonds_valence == 1
                    {
                        single_bond_oxygen_minus[number_single_bond_oxygen_minus as usize] =
                            atom_number as i16;
                        number_single_bond_oxygen_minus =
                            number_single_bond_oxygen_minus.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, edge_number, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    } else if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                        && difference.nFixHInChI == -1
                        && difference.nMobHInChI == 1
                        && difference.nFixHRevrs == 0
                        && difference.nMobHRevrs == 0
                        && difference.nAtChargeRevrs == 0
                        && atom.num_H == 0
                        && atom.valence == 1
                        && atom.chem_bonds_valence == 2
                    {
                        double_bond_oxygen[number_double_bond_oxygen as usize] = atom_number as i16;
                        number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, edge_number, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let number_to_try = number_single_bond_oxygen_minus.min(number_double_bond_oxygen);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let delta = 1_i32;
                let mut index = 0_i32;
                while index < number_single_bond_oxygen_minus && current_success < number_to_try {
                    let atom_number = i32::from(single_bond_oxygen_minus[index as usize]);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(delta);
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
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));

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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            edge.flow = edge.flow.wrapping_add(delta);
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
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                    }
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                current_edges.num_edges = 0;
            }
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 02, ichirvr3.c:648-795.
        if comparison.len_c2at >= 1 {
            let mut double_bond_positive_hetero = [0_i16; MAX_DIFF_FIXH as usize];
            let mut single_bond_nh = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_positive_hetero = 0_i32;
            let mut number_single_bond_nh = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let input_mobile_h = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut number_zero_returns = 0_i32;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCPlusGroupEdge
                    .wrapping_sub(1);
                let edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if number_single_bond_nh < MAX_DIFF_FIXH as i32
                    && ((difference.nValElectr == 5 && difference.nPeriodNum == 1)
                        || difference.nValElectr == 6)
                    && edge_available
                    && difference.nFixHInChI > 0
                    && difference.nFixHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H != 0
                    && atom.valence == atom.chem_bonds_valence
                {
                    single_bond_nh[number_single_bond_nh as usize] = atom_number as i16;
                    number_single_bond_nh = number_single_bond_nh.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let fixed_h = if pStruct.fixed_H.is_null() {
                    0
                } else {
                    *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                };
                let mobile_h = if input_mobile_h.is_null() {
                    0
                } else {
                    *heap
                        .slice(input_mobile_h.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                };
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_double_bond_positive_hetero < MAX_DIFF_FIXH as i32
                    && atom.charge == 1
                    && atom.num_H == 0
                    && atom.valence < atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && (valence.cNumValenceElectrons == 6
                        || valence.cNumValenceElectrons == 7
                        || (valence.cNumValenceElectrons == 5 && valence.cPeriodicRowNumber > 1))
                    && endpoint == 0
                    && fixed_h == 0
                    && mobile_h == 0
                    && edge_available
                {
                    double_bond_positive_hetero[number_double_bond_positive_hetero as usize] =
                        atom_number as i16;
                    number_double_bond_positive_hetero =
                        number_double_bond_positive_hetero.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            let number_to_try = number_double_bond_positive_hetero.min(number_single_bond_nh);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let delta = 1_i32;
                loop {
                    let mut index = 0_i32;
                    while index < number_single_bond_nh && current_success < number_to_try {
                        let atom_number = i32::from(single_bond_nh[index as usize]);
                        let atom_index = usize::try_from(atom_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let edge_number = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCPlusGroupEdge
                            .wrapping_sub(1);
                        let edge_index = usize::try_from(edge_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let edge_before = heap
                            .slice(pBNS.edge.as_const())?
                            .get(edge_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if edge_before.flow == 0 {
                            index = index.wrapping_add(1);
                            continue;
                        }
                        let first_vertex = i32::from(edge_before.neighbor1);
                        let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                            edge.flow = edge.flow.wrapping_sub(delta);
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
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));

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
                            && delta_charge == -1
                        {
                            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                            if ret > 0 {
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                current_success = current_success.wrapping_add(1);
                            }
                        } else {
                            number_zero_returns =
                                number_zero_returns.wrapping_add(i32::from(ret == 0));
                            {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden =
                                    (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                                edge.flow = edge.flow.wrapping_add(delta);
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
                                    vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                                }
                            }
                            pBNS.tot_st_flow =
                                pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                        }
                        index = index.wrapping_add(1);
                    }
                    if number_zero_returns == number_to_try
                        && !allowed_nitrogen_flower_edges
                        && nitrogen_flower_edges.num_edges != 0
                    {
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &nitrogen_flower_edges,
                            forbidden_edge_mask,
                        )?;
                        allowed_nitrogen_flower_edges = true;
                        continue;
                    }
                    break;
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;

            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 03, ichirvr3.c:796-944.
        if comparison.len_c2at >= 1
            && comparison.nNumTgRevrs == 1
            && (comparison.nNumEndpRevrs > comparison.nNumEndpInChI || comparison.nNumTgInChI > 1)
        {
            let mut single_bond_nitrogen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_nitrogen_minus = 0_i32;
            let mut number_double_bond_oxygen = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let mut current_success = 0_i32;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                    && difference.nValElectr == 6
                    && difference.endptInChI == 0
                    && edge_available
                    && difference.nFixHInChI == -1
                    && difference.nMobHInChI == 1
                    && difference.endptRevrs != 0
                    && difference.nFixHRevrs == 0
                    && difference.nMobHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
                    && atom.valence == 1
                    && atom.chem_bonds_valence == 2
                {
                    double_bond_oxygen[number_double_bond_oxygen as usize] = atom_number as i16;
                    number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_endpoint = !mobile_h_reversed.is_null()
                    && heap
                        .slice(mobile_h_reversed.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint
                        != 0;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_single_bond_nitrogen_minus < MAX_DIFF_FIXH as i32
                    && atom.charge == -1
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 5
                    && valence.cPeriodicRowNumber == 1
                    && reversed_endpoint
                    && edge_available
                {
                    single_bond_nitrogen_minus[number_single_bond_nitrogen_minus as usize] =
                        atom_number as i16;
                    number_single_bond_nitrogen_minus =
                        number_single_bond_nitrogen_minus.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            let number_to_try = number_single_bond_nitrogen_minus.min(number_double_bond_oxygen);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let delta = 1_i32;
                let mut index = 0_i32;
                while index < number_single_bond_nitrogen_minus && current_success < number_to_try {
                    let atom_number = i32::from(single_bond_nitrogen_minus[index as usize]);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(delta);
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
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));

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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            edge.flow = edge.flow.wrapping_add(delta);
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
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                    }
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;

            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 03a, ichirvr3.c:945-1083.
        if comparison.nNumTgRevrs == 1
            && (comparison.nNumEndpRevrs > comparison.nNumEndpInChI || comparison.nNumTgInChI > 1)
        {
            let mut single_bond_nitrogen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_nitrogen_minus = 0_i32;
            let mut number_double_bond_oxygen = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let fixed_h_input = pStruct.fixed_H;
            let mut current_success = 0_i32;
            current_edges.num_edges = 0;

            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_endpoint = !mobile_h_reversed.is_null()
                    && heap
                        .slice(mobile_h_reversed.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint
                        != 0;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_single_bond_nitrogen_minus < MAX_DIFF_FIXH as i32
                    && atom.charge == -1
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 5
                    && valence.cPeriodicRowNumber == 1
                    && reversed_endpoint
                    && edge_available
                {
                    single_bond_nitrogen_minus[number_single_bond_nitrogen_minus as usize] =
                        atom_number as i16;
                    number_single_bond_nitrogen_minus =
                        number_single_bond_nitrogen_minus.wrapping_add(1);
                } else {
                    let original_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                    let original_h_acceptor = !mobile_h_input.is_null()
                        && *heap
                            .slice(mobile_h_input.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            == 1
                        && !fixed_h_input.is_null()
                        && *heap
                            .slice(fixed_h_input.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            == -1;
                    if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                        && atom.charge == 0
                        && i32::from(atom.valence) + 1 == i32::from(atom.chem_bonds_valence)
                        && valence.cMetal == 0
                        && valence.cNumValenceElectrons == 6
                        && reversed_endpoint
                        && (original_endpoint || original_h_acceptor)
                        && edge_available
                    {
                        double_bond_oxygen[number_double_bond_oxygen as usize] = atom_number as i16;
                        number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            let number_to_try = number_single_bond_nitrogen_minus.min(number_double_bond_oxygen);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let delta = 1_i32;
                let mut index = 0_i32;
                while index < number_single_bond_nitrogen_minus && current_success < number_to_try {
                    let atom_number = i32::from(single_bond_nitrogen_minus[index as usize]);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(delta);
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
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));

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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            edge.flow = edge.flow.wrapping_add(delta);
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
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                    }
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;

            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 04, ichirvr3.c:1080-1227.
        if comparison.len_c2at >= 1
            && comparison.nNumTgInChI == 1
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI || comparison.nNumTgRevrs > 1)
        {
            let mut single_bond_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_charged = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_neutral = 0_i32;
            let mut number_double_bond_charged = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;

            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                let fixed_h = !pStruct.fixed_H.is_null()
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let mobile_h = !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0;
                if number_double_bond_charged < MAX_DIFF_FIXH as i32
                    && atom.charge == 1
                    && atom.num_H != 0
                    && atom.valence < atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && (valence.cNumValenceElectrons == 6
                        || (valence.cNumValenceElectrons == 5 && valence.cPeriodicRowNumber == 1))
                    && endpoint
                    && fixed_h
                    && edge_available
                {
                    double_bond_charged[number_double_bond_charged as usize] = atom_number as i16;
                    number_double_bond_charged = number_double_bond_charged.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                } else if number_single_bond_neutral < MAX_DIFF_FIXH as i32
                    && atom.charge == 0
                    && atom.num_H == 0
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && (valence.cNumValenceElectrons == 6 || valence.cNumValenceElectrons == 5)
                    && !fixed_h
                    && !mobile_h
                    && edge_available
                {
                    single_bond_neutral[number_single_bond_neutral as usize] = atom_number as i16;
                    number_single_bond_neutral = number_single_bond_neutral.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            let number_to_try = number_single_bond_neutral.min(number_double_bond_charged);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let delta = 1_i32;
                let mut index = 0_i32;
                while index < number_single_bond_neutral && current_success < number_to_try {
                    let atom_number = i32::from(single_bond_neutral[index as usize]);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(delta);
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
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));

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
                        && delta_charge == -1
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            edge.flow = edge.flow.wrapping_add(delta);
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
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                    }
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;

            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 05, ichirvr3.c:1225-1360.
        if comparison.len_c2at > 1 {
            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut single_bond_nh = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_oxygen = 0_i32;
            let mut number_single_bond_nh = 0_i32;
            let mut current_success = 0_i32;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let plus_edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_single_bond_nh < MAX_DIFF_FIXH as i32
                    && ((difference.nValElectr == 5 && difference.nPeriodNum == 1)
                        || difference.nValElectr == 6)
                    && difference.endptInChI == 0
                    && plus_edge_available
                    && difference.nFixHInChI == 1
                    && difference.nFixHRevrs == 0
                    && difference.nMobHRevrs != 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H != 0
                    && difference.endptRevrs == 0
                    && atom.valence == atom.chem_bonds_valence
                {
                    single_bond_nh[number_single_bond_nh as usize] = atom_number as i16;
                    number_single_bond_nh = number_single_bond_nh.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                } else {
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    let minus_edge_available = minus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .forbidden
                            == 0;
                    if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                        && difference.nValElectr == 6
                        && difference.endptInChI == 0
                        && minus_edge_available
                        && difference.nFixHInChI == -1
                        && difference.nMobHInChI == 1
                        && difference.nFixHRevrs == 0
                        && difference.nMobHRevrs == 0
                        && difference.nAtChargeRevrs == 0
                        && atom.num_H == 0
                        && difference.endptRevrs == 0
                        && i32::from(atom.valence) + 1 == i32::from(atom.chem_bonds_valence)
                    {
                        double_bond_oxygen[number_double_bond_oxygen as usize] = atom_number as i16;
                        number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let number_to_try = number_double_bond_oxygen.min(number_single_bond_nh);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let delta = 1_i32;
                let mut index = 0_i32;
                while index < number_single_bond_nh && current_success < number_to_try {
                    let atom_number = i32::from(single_bond_nh[index as usize]);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(delta);
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
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));

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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            edge.flow = edge.flow.wrapping_add(delta);
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
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                    }
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;

            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 06c, ichirvr3.c:1358-1525.
        if !pStruct.fixed_H.is_null()
            && !pStruct.endpoint.is_null()
            && comparison.nChargeFixHInChI > 0
            && comparison.nChargeFixHInChI > comparison.nChargeMobHInChI
        {
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mut current_charge_edges = EDGE_LIST::default();
            let mut current_success = 0_i32;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_CLEAR)?;
            current_edges.num_edges = 0;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCPlusGroupEdge
                    .wrapping_sub(1);
                let edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if (difference.nValElectr == 6
                    || (difference.nValElectr == 5 && difference.nPeriodNum == 1))
                    && difference.endptInChI == 0
                    && edge_available
                    && difference.nFixHInChI == 1
                    && difference.nMobHInChI == 0
                    && difference.nFixHRevrs == 0
                    && difference.nMobHRevrs == 1
                    && atom.num_H != 0
                    && atom.valence == atom.chem_bonds_valence
                {
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        // The source jumps directly to exit_function here and does not
                        // release CurChargeEdges on this allocation-failure path.
                        return Ok(());
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let mut leave_case = false;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms && !leave_case {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if endpoint != 0 || atom.charge != 1 || atom.radical != 0 || valence.cMetal != 0 {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let positively_charged = plus_edge >= 0
                    && {
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden == 0 && edge.flow == 0
                    }
                    && valence.cNumValenceElectrons >= 5;
                if positively_charged {
                    let mut neighbor_number = 0_i32;
                    while neighbor_number < i32::from(atom.valence) {
                        let adjacent_atom_number = i32::from(
                            *atom
                                .neighbor
                                .get(
                                    usize::try_from(neighbor_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        let adjacent_index = usize::try_from(adjacent_atom_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_atom = heap
                            .slice(at2.as_const())?
                            .get(adjacent_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_valence = pVA
                            .get(adjacent_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_plus_edge = adjacent_valence.nCPlusGroupEdge.wrapping_sub(1);
                        let adjacent_positive = adjacent_plus_edge >= 0 && {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(adjacent_plus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden == 0 && edge.flow == 0
                        };
                        if adjacent_atom.charge == 1
                            && adjacent_valence.cMetal == 0
                            && adjacent_positive
                        {
                            if FindInEdgeList(heap, &current_edges, plus_edge)? < 0
                                && FindInEdgeList(heap, &current_charge_edges, plus_edge)? < 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut current_charge_edges,
                                    plus_edge,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    leave_case = true;
                                    break;
                                }
                            }
                            if FindInEdgeList(heap, &current_edges, adjacent_plus_edge)? < 0
                                && FindInEdgeList(heap, &current_charge_edges, adjacent_plus_edge)?
                                    < 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut current_charge_edges,
                                    adjacent_plus_edge,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    leave_case = true;
                                    break;
                                }
                            }
                        }
                        neighbor_number = neighbor_number.wrapping_add(1);
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            if !leave_case {
                let number_to_try = current_edges.num_edges.min(current_charge_edges.num_edges);
                if number_to_try != 0 {
                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    RemoveForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &current_charge_edges,
                        forbidden_edge_mask,
                    )?;
                    let delta = 1_i32;
                    let mut index = 0_i32;
                    while index < current_edges.num_edges && current_success < number_to_try {
                        let edge_number = *heap
                            .slice(current_edges.pnEdges.as_const())?
                            .get(
                                usize::try_from(index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let edge_index = usize::try_from(edge_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let edge_before = heap
                            .slice(pBNS.edge.as_const())?
                            .get(edge_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if edge_before.flow == 0 {
                            index = index.wrapping_add(1);
                            continue;
                        }
                        let first_vertex = i32::from(edge_before.neighbor1);
                        let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                        heap.slice_mut(pBNS.edge)?[edge_index].flow =
                            edge_before.flow.wrapping_sub(delta);
                        {
                            let vertices = heap.slice_mut(pBNS.vert)?;
                            for vertex_number in [first_vertex, second_vertex] {
                                let vertex = vertices
                                    .get_mut(
                                        usize::try_from(vertex_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));

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
                            && delta_charge == -1
                        {
                            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                            if ret > 0 {
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                current_success = current_success.wrapping_add(1);
                            }
                        } else {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.flow = edge.flow.wrapping_add(delta);
                            {
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [first_vertex, second_vertex] {
                                    let vertex = vertices
                                        .get_mut(
                                            usize::try_from(vertex_number)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                                }
                            }
                            pBNS.tot_st_flow =
                                pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                        }
                        index = index.wrapping_add(1);
                    }
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                }
            }

            current_edges.num_edges = 0;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_FREE)?;
            if ret < 0 {
                return Ok(());
            }
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 06d, ichirvr3.c:1526-1682.
        if comparison.len_c2at >= 2 {
            let mut current_charge_edges = EDGE_LIST::default();
            let mut current_success = 0_i32;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_CLEAR)?;
            current_edges.num_edges = 0;
            let mut leave_case = false;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) && !leave_case {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_index = usize::try_from(i32::from(difference.atomNumber))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCPlusGroupEdge
                    .wrapping_sub(1);
                let edge = if plus_edge >= 0 {
                    Some(
                        heap.slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    )
                } else {
                    None
                };
                let reconstructed_charged = i32::from(difference.nMobHRevrs) + 1
                    == i32::from(difference.nFixHRevrs)
                    && difference.nFixHRevrs > 0
                    && difference.endptRevrs == 0
                    && difference.nAtChargeRevrs == 1
                    && ((difference.nFixHInChI == 0
                        && difference.nMobHInChI == difference.nFixHRevrs)
                        || (difference.nFixHInChI == difference.nFixHRevrs
                            && difference.endptInChI != 0));
                if reconstructed_charged
                    && edge.is_some_and(|edge| edge.forbidden == 0 && edge.flow == 0)
                {
                    ret = AddToEdgeList(heap, &mut current_charge_edges, plus_edge, INC_ADD_EDGE)?;
                    leave_case = ret != 0;
                } else {
                    let original_has_h = i32::from(difference.nMobHInChI) + 1
                        == i32::from(difference.nFixHInChI)
                        && difference.nFixHInChI > 0
                        && difference.endptInChI == 0
                        && difference.nAtChargeRevrs == 0
                        && ((difference.nFixHRevrs == 0
                            && difference.nMobHRevrs == difference.nFixHInChI)
                            || (difference.nFixHRevrs == difference.nFixHInChI
                                && difference.endptRevrs != 0));
                    if original_has_h
                        && edge.is_some_and(|edge| edge.forbidden == 0 && edge.flow != 0)
                    {
                        ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                        leave_case = ret != 0;
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            if !leave_case {
                let number_to_try = current_edges.num_edges.min(current_charge_edges.num_edges);
                if number_to_try != 0 {
                    let sulfur_may_be_forbidden = i32::from(sulfur_flower_edges.num_edges > 0);
                    let mut sulfur_is_forbidden = sulfur_may_be_forbidden;
                    while sulfur_is_forbidden >= 0 {
                        if sulfur_is_forbidden != 0 {
                            SetForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &sulfur_flower_edges,
                                forbidden_edge_mask,
                            )?;
                        }
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &current_charge_edges,
                            forbidden_edge_mask,
                        )?;
                        let delta = 1_i32;
                        let mut index = 0_i32;
                        while index < current_edges.num_edges && current_success < number_to_try {
                            let edge_number = *heap
                                .slice(current_edges.pnEdges.as_const())?
                                .get(
                                    usize::try_from(index)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let edge_index = usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let edge_before = heap
                                .slice(pBNS.edge.as_const())?
                                .get(edge_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge_before.flow == 0 {
                                index = index.wrapping_add(1);
                                continue;
                            }
                            let first_vertex = i32::from(edge_before.neighbor1);
                            let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                            heap.slice_mut(pBNS.edge)?[edge_index].flow =
                                edge_before.flow.wrapping_sub(delta);
                            {
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [first_vertex, second_vertex] {
                                    let vertex = vertices
                                        .get_mut(
                                            usize::try_from(vertex_number)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                                }
                            }
                            pBNS.tot_st_flow =
                                pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));
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
                                && delta_charge == -1
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
                                    n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                    current_success = current_success.wrapping_add(1);
                                }
                            } else {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.flow = edge.flow.wrapping_add(delta);
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [first_vertex, second_vertex] {
                                    let vertex = vertices
                                        .get_mut(
                                            usize::try_from(vertex_number)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                                }
                                pBNS.tot_st_flow =
                                    pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                            }
                            index = index.wrapping_add(1);
                        }
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &all_charge_edges,
                            forbidden_edge_mask,
                        )?;
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &sulfur_flower_edges,
                            forbidden_edge_mask,
                        )?;
                        sulfur_is_forbidden = sulfur_is_forbidden.wrapping_sub(1);
                    }
                }
            }
            current_edges.num_edges = 0;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_FREE)?;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 06, ichirvr3.c:1683-1815.
        if comparison.len_c2at >= 2 {
            let mut double_bond_nhn_plus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut single_bond_nhm_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_nhn_plus = 0_i32;
            let mut number_single_bond_nhm_neutral = 0_i32;
            let mut current_success = 0_i32;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCPlusGroupEdge
                    .wrapping_sub(1);
                let edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if (difference.nValElectr == 6
                    || (difference.nValElectr == 5 && difference.nPeriodNum == 1))
                    && difference.endptInChI == 0
                    && edge_available
                {
                    if number_single_bond_nhm_neutral < MAX_DIFF_FIXH as i32
                        && difference.nFixHInChI == 1
                        && difference.nFixHRevrs == 0
                        && difference.nAtChargeRevrs == 0
                        && atom.num_H != 0
                        && atom.valence == atom.chem_bonds_valence
                    {
                        single_bond_nhm_neutral[number_single_bond_nhm_neutral as usize] =
                            atom_number as i16;
                        number_single_bond_nhm_neutral =
                            number_single_bond_nhm_neutral.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    } else if number_double_bond_nhn_plus < MAX_DIFF_FIXH as i32
                        && difference.nFixHInChI == 0
                        && difference.nFixHRevrs == 1
                        && difference.nAtChargeRevrs == 1
                        && atom.num_H != 0
                        && atom.valence < atom.chem_bonds_valence
                    {
                        double_bond_nhn_plus[number_double_bond_nhn_plus as usize] =
                            atom_number as i16;
                        number_double_bond_nhn_plus = number_double_bond_nhn_plus.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let number_to_try = number_single_bond_nhm_neutral.min(number_double_bond_nhn_plus);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let delta = 1_i32;
                let mut index = 0_i32;
                while index < number_single_bond_nhm_neutral && current_success < number_to_try {
                    let atom_number = i32::from(single_bond_nhm_neutral[index as usize]);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(delta);
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
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));
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
                        && delta_charge == -1
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        edge.flow = edge.flow.wrapping_add(delta);
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_vertex, second_vertex] {
                            let vertex = vertices
                                .get_mut(
                                    usize::try_from(vertex_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                    }
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 06a, ichirvr3.c:1810-1949.
        let normalized_tautomeric = pStruct.pOne_norm_data[TAUT_YES as usize];
        let fixed_bond_atoms = if normalized_tautomeric.is_null() {
            SourceMutPointer::null()
        } else {
            heap.slice(normalized_tautomeric.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .at_fixed_bonds
        };
        if ((comparison.nNumTgInChI > comparison.nNumTgRevrs && comparison.nNumTgRevrs == 1)
            || comparison.nNumEndpInChI < comparison.nNumEndpRevrs)
            && pStruct.nNumRemovedProtonsMobHInChI
                == i32::from(pStruct.One_ti.tni.nNumRemovedProtons)
            && !pStruct.fixed_H.is_null()
            && !pStruct.endpoint.is_null()
            && !normalized_tautomeric.is_null()
            && !fixed_bond_atoms.is_null()
        {
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_charge_edges = EDGE_LIST::default();
            let mut current_success = 0_i32;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_CLEAR)?;
            current_edges.num_edges = 0;
            let mut leave_case = false;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms && !leave_case {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
                {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let edge = if plus_edge >= 0 {
                    Some(
                        heap.slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone(),
                    )
                } else {
                    None
                };
                let neutral = *heap
                    .slice(pStruct.fixed_H.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
                    && !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        == 0
                    && atom.charge == 0
                    && atom.radical == 0
                    && edge
                        .as_ref()
                        .is_some_and(|edge| edge.forbidden == 0 && edge.flow != 0);
                if neutral {
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        leave_case = true;
                    }
                }
                if !leave_case {
                    let fixed_bond_charge = heap
                        .slice(fixed_bond_atoms.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .charge;
                    let charged = atom.charge == 1
                        && atom.num_H == 0
                        && valence.cNumValenceElectrons == 5
                        && valence.cPeriodicRowNumber == 1
                        && fixed_bond_charge == 0
                        && edge
                            .as_ref()
                            .is_some_and(|edge| edge.forbidden == 0 && edge.flow == 0);
                    if charged {
                        ret = AddToEdgeList(
                            heap,
                            &mut current_charge_edges,
                            plus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            leave_case = true;
                        }
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            if !leave_case {
                let number_to_try = current_edges.num_edges.min(current_charge_edges.num_edges);
                if number_to_try != 0 {
                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    RemoveForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &current_charge_edges,
                        forbidden_edge_mask,
                    )?;
                    let delta = 1_i32;
                    let mut index = 0_i32;
                    while index < current_edges.num_edges && current_success < number_to_try {
                        let edge_number = *heap
                            .slice(current_edges.pnEdges.as_const())?
                            .get(
                                usize::try_from(index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let edge_index = usize::try_from(edge_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let edge_before = heap
                            .slice(pBNS.edge.as_const())?
                            .get(edge_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if edge_before.flow == 0 {
                            index = index.wrapping_add(1);
                            continue;
                        }
                        let first_vertex = i32::from(edge_before.neighbor1);
                        let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                        heap.slice_mut(pBNS.edge)?[edge_index].flow =
                            edge_before.flow.wrapping_sub(delta);
                        {
                            let vertices = heap.slice_mut(pBNS.vert)?;
                            for vertex_number in [first_vertex, second_vertex] {
                                let vertex = vertices
                                    .get_mut(
                                        usize::try_from(vertex_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(delta);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2_i32.wrapping_mul(delta));
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
                            && delta_charge == -1
                        {
                            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                            if ret > 0 {
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                current_success = current_success.wrapping_add(1);
                            }
                        } else {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.flow = edge.flow.wrapping_add(delta);
                            let vertices = heap.slice_mut(pBNS.vert)?;
                            for vertex_number in [first_vertex, second_vertex] {
                                let vertex = vertices
                                    .get_mut(
                                        usize::try_from(vertex_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(delta);
                            }
                            pBNS.tot_st_flow =
                                pBNS.tot_st_flow.wrapping_add(2_i32.wrapping_mul(delta));
                        }
                        index = index.wrapping_add(1);
                    }
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                }
            }
            current_edges.num_edges = 0;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_FREE)?;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 06b, ichirvr3.c:1950-2096.
        let normalized_tautomeric = pStruct.pOne_norm_data[TAUT_YES as usize];
        let fixed_bond_atoms = if normalized_tautomeric.is_null() {
            SourceMutPointer::null()
        } else {
            heap.slice(normalized_tautomeric.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .at_fixed_bonds
        };
        if ((comparison.nNumTgInChI > comparison.nNumTgRevrs && comparison.nNumTgRevrs == 1)
            || comparison.nNumEndpInChI < comparison.nNumEndpRevrs)
            && (pStruct.nNumRemovedProtonsMobHInChI
                == i32::from(pStruct.One_ti.tni.nNumRemovedProtons)
                || pStruct.nNumRemovedProtonsMobHInChI
                    > i32::from(pStruct.One_ti.tni.nNumRemovedProtons))
            && !pStruct.fixed_H.is_null()
            && !pStruct.endpoint.is_null()
            && !normalized_tautomeric.is_null()
            && !fixed_bond_atoms.is_null()
        {
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_charge_edges = EDGE_LIST::default();
            let mut current_success = 0_i32;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_CLEAR)?;
            current_edges.num_edges = 0;
            let mut leave_case = false;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms && !leave_case {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
                {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let edge = if plus_edge >= 0 {
                    Some(
                        heap.slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone(),
                    )
                } else {
                    None
                };
                let neutral = *heap
                    .slice(pStruct.fixed_H.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
                    && !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        == 0
                    && atom.charge == 0
                    && atom.radical == 0
                    && edge
                        .as_ref()
                        .is_some_and(|edge| edge.forbidden == 0 && edge.flow != 0);
                if neutral {
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        leave_case = true;
                    }
                }
                if !leave_case {
                    let fixed_bond_charge = heap
                        .slice(fixed_bond_atoms.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .charge;
                    let charged = atom.charge == 1
                        && atom.num_H == 0
                        && (valence.cNumValenceElectrons == 6 || valence.cPeriodicRowNumber == 7)
                        && fixed_bond_charge == 1
                        && edge
                            .as_ref()
                            .is_some_and(|edge| edge.forbidden == 0 && edge.flow == 0);
                    if charged {
                        ret = AddToEdgeList(
                            heap,
                            &mut current_charge_edges,
                            plus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            leave_case = true;
                        }
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            if !leave_case {
                let number_to_try = current_edges.num_edges.min(current_charge_edges.num_edges);
                if number_to_try != 0 {
                    let sulfur_may_be_forbidden = i32::from(sulfur_flower_edges.num_edges > 0);
                    let mut sulfur_is_forbidden = sulfur_may_be_forbidden;
                    while sulfur_is_forbidden >= 0 {
                        if sulfur_is_forbidden != 0 {
                            SetForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &sulfur_flower_edges,
                                forbidden_edge_mask,
                            )?;
                        }
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &current_charge_edges,
                            forbidden_edge_mask,
                        )?;
                        let delta = 1_i32;
                        let mut index = 0_i32;
                        while index < current_edges.num_edges && current_success < number_to_try {
                            let edge_number = *heap
                                .slice(current_edges.pnEdges.as_const())?
                                .get(
                                    usize::try_from(index)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let edge_index = usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let edge_before = heap
                                .slice(pBNS.edge.as_const())?
                                .get(edge_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge_before.flow == 0 {
                                index = index.wrapping_add(1);
                                continue;
                            }
                            let first_vertex = i32::from(edge_before.neighbor1);
                            let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                            heap.slice_mut(pBNS.edge)?[edge_index].flow =
                                edge_before.flow.wrapping_sub(1);
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
                                && delta_charge == -1
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
                                    n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                    current_success = current_success.wrapping_add(1);
                                }
                            } else {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.flow = edge.flow.wrapping_add(1);
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
                            index = index.wrapping_add(1);
                        }
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &all_charge_edges,
                            forbidden_edge_mask,
                        )?;
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &sulfur_flower_edges,
                            forbidden_edge_mask,
                        )?;
                        sulfur_is_forbidden = sulfur_is_forbidden.wrapping_sub(1);
                    }
                }
            }
            current_edges.num_edges = 0;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_FREE)?;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 06e, ichirvr3.c:2091-2328.
        let normalized_tautomeric = pStruct.pOne_norm_data[TAUT_YES as usize];
        let normalized_at_fixed_bonds = if normalized_tautomeric.is_null() {
            SourceMutPointer::null()
        } else {
            heap.slice(normalized_tautomeric.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .at_fixed_bonds
        };
        let fixed_bond_atoms = if normalized_tautomeric.is_null() {
            SourceMutPointer::null()
        } else {
            let normalized = heap
                .slice(normalized_tautomeric.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if !normalized.at_fixed_bonds.is_null() {
                normalized.at_fixed_bonds
            } else {
                normalized.at
            }
        };
        if comparison.nNumTgInChI > 1
            && (pStruct.nNumRemovedProtonsMobHInChI > 0 || pStruct.ti.tni.nNumRemovedProtons > 0)
            && !pStruct.fixed_H.is_null()
            && !pStruct.endpoint.is_null()
            && !normalized_tautomeric.is_null()
            && !normalized_at_fixed_bonds.is_null()
        {
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mut tgroup_differences: [TgDiffHChgFH; MAX_DIFF_FIXH as usize] =
                std::array::from_fn(|_| TgDiffHChgFH::default());
            let mut current_charge_edges = EDGE_LIST::default();
            let mut endpoint_list = EDGE_LIST::default();
            let mut current_success = 0_i32;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_CLEAR)?;
            let _ = AllocEdgeList(heap, &mut endpoint_list, EDGE_LIST_CLEAR)?;
            current_edges.num_edges = 0;
            let mut leave_case = false;
            let mut number_wrong_tgroups = 0_i32;
            if fixed_bond_atoms.is_null() {
                leave_case = true;
            } else {
                let at2_snapshot = heap.slice(at2.as_const())?.to_vec();
                let fixed_bond_snapshot = heap.slice(fixed_bond_atoms.as_const())?.to_vec();
                let canonical_snapshot = heap.slice(canonical_to_atom.as_const())?.to_vec();
                number_wrong_tgroups = FillTgDiffHChgFH(
                    heap,
                    &mut tgroup_differences,
                    MAX_DIFF_FIXH as i32,
                    &at2_snapshot,
                    &fixed_bond_snapshot,
                    &canonical_snapshot,
                    pVA,
                    &pStruct.ti,
                    &mut endpoint_list,
                )?;
                if number_wrong_tgroups < 1 {
                    leave_case = true;
                }
            }
            let mut number_remove_plus = 0_i32;
            let mut number_add_plus = 0_i32;
            if !leave_case {
                let mut group_index = 0_i32;
                while group_index < number_wrong_tgroups && !leave_case {
                    let difference = &tgroup_differences[usize::try_from(group_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                    if difference.nNumHInchi > difference.nNumHNorml
                        && difference.nNumPRevrs > difference.nNumPNorml
                        && difference.n[fNumRPosChgH as usize] != 0
                    {
                        let needed = (i32::from(difference.nNumHInchi)
                            - i32::from(difference.nNumHNorml))
                        .min(i32::from(difference.n[fNumRPosChgH as usize]));
                        number_remove_plus = number_remove_plus.wrapping_add(needed);
                        let offset = i32::from(difference.i[fNumRPosChgH as usize]);
                        let mut endpoint_number = 0_i32;
                        while endpoint_number < i32::from(difference.n[fNumRPosChgH as usize]) {
                            let endpoint_index = offset.wrapping_add(endpoint_number);
                            let atom_number = *heap
                                .slice(endpoint_list.pnEdges.as_const())?
                                .get(
                                    usize::try_from(endpoint_index)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let atom_index = usize::try_from(atom_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let edge_number = pVA
                                .get(atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nCPlusGroupEdge
                                .wrapping_sub(1);
                            ret = AddToEdgeList(
                                heap,
                                &mut current_charge_edges,
                                edge_number,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                leave_case = true;
                                break;
                            }
                            endpoint_number = endpoint_number.wrapping_add(1);
                        }
                    } else if difference.nNumHInchi < difference.nNumHNorml
                        && difference.n[fNumRNeutrlH as usize] != 0
                    {
                        let needed = (i32::from(difference.nNumHNorml)
                            - i32::from(difference.nNumHInchi))
                        .min(i32::from(difference.n[fNumRNeutrlH as usize]));
                        number_add_plus = number_add_plus.wrapping_add(needed);
                        let offset = i32::from(difference.i[fNumRNeutrlH as usize]);
                        let mut endpoint_number = 0_i32;
                        while endpoint_number < i32::from(difference.n[fNumRNeutrlH as usize]) {
                            let endpoint_index = offset.wrapping_add(endpoint_number);
                            let atom_number = *heap
                                .slice(endpoint_list.pnEdges.as_const())?
                                .get(
                                    usize::try_from(endpoint_index)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let atom_index = usize::try_from(atom_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let edge_number = pVA
                                .get(atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nCPlusGroupEdge
                                .wrapping_sub(1);
                            ret =
                                AddToEdgeList(heap, &mut current_edges, edge_number, INC_ADD_EDGE)?;
                            if ret != 0 {
                                leave_case = true;
                                break;
                            }
                            endpoint_number = endpoint_number.wrapping_add(1);
                        }
                    }
                    group_index = group_index.wrapping_add(1);
                }
            }
            let mut number_move_plus = number_remove_plus.min(number_add_plus);
            if !leave_case && current_edges.num_edges > 0 && current_charge_edges.num_edges > 0 {
                let mut group_index = 0_i32;
                while number_move_plus > 0 && group_index < number_wrong_tgroups && !leave_case {
                    let difference = tgroup_differences[usize::try_from(group_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                    .clone();
                    if difference.nNumHInchi > difference.nNumHNorml
                        && difference.nNumPRevrs > difference.nNumPNorml
                        && difference.n[fNumRPosChgH as usize] != 0
                    {
                        let mut number_remove =
                            i32::from(difference.nNumHInchi) - i32::from(difference.nNumHNorml);
                        if number_remove < i32::from(difference.n[fNumRPosChgH as usize]) {
                            number_remove = i32::from(difference.n[fNumRPosChgH as usize]);
                        }
                        let offset = i32::from(difference.i[fNumRPosChgH as usize]);
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                        let mut endpoint_number = 0_i32;
                        while number_move_plus > 0
                            && number_remove > 0
                            && endpoint_number < i32::from(difference.n[fNumRPosChgH as usize])
                        {
                            let atom_number = *heap
                                .slice(endpoint_list.pnEdges.as_const())?
                                .get(
                                    usize::try_from(offset.wrapping_add(endpoint_number))
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let atom_index = usize::try_from(atom_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let edge_number = pVA
                                .get(atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nCPlusGroupEdge
                                .wrapping_sub(1);
                            let edge_index = usize::try_from(edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let edge_before = heap
                                .slice(pBNS.edge.as_const())?
                                .get(edge_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge_before.flow != 0 {
                                endpoint_number = endpoint_number.wrapping_add(1);
                                continue;
                            }
                            let first_vertex = i32::from(edge_before.neighbor1);
                            let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;

                            let first_vertex_data = heap
                                .slice(pBNS.vert.as_const())?
                                .get(
                                    usize::try_from(first_vertex)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut adjacent_index = i32::from(first_vertex_data.num_adj_edges) - 1;
                            let mut first_neighbor = None;
                            while adjacent_index >= 0 {
                                let adjacent_edge_number = *heap
                                    .slice(first_vertex_data.iedge.as_const())?
                                    .get(
                                        usize::try_from(adjacent_index)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let adjacent_edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(adjacent_edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                                    first_neighbor = Some((
                                        adjacent_edge_number,
                                        i32::from(adjacent_edge.neighbor12) ^ first_vertex,
                                    ));
                                    break;
                                }
                                adjacent_index -= 1;
                            }
                            let Some((first_neighbor_edge, first_neighbor_vertex)) = first_neighbor
                            else {
                                endpoint_number = endpoint_number.wrapping_add(1);
                                continue;
                            };

                            let second_vertex_data = heap
                                .slice(pBNS.vert.as_const())?
                                .get(
                                    usize::try_from(second_vertex)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut adjacent_index =
                                i32::from(second_vertex_data.num_adj_edges) - 2;
                            let mut second_neighbor = None;
                            while adjacent_index >= 0 {
                                let adjacent_edge_number = *heap
                                    .slice(second_vertex_data.iedge.as_const())?
                                    .get(
                                        usize::try_from(adjacent_index)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let adjacent_edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(adjacent_edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                                    second_neighbor = Some((
                                        adjacent_edge_number,
                                        i32::from(adjacent_edge.neighbor12) ^ second_vertex,
                                    ));
                                    break;
                                }
                                adjacent_index -= 1;
                            }
                            let Some((second_neighbor_edge, second_neighbor_vertex)) =
                                second_neighbor
                            else {
                                endpoint_number = endpoint_number.wrapping_add(1);
                                continue;
                            };

                            heap.slice_mut(pBNS.edge)?[edge_index].flow =
                                edge_before.flow.wrapping_add(1);
                            for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                                let adjacent_edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(
                                        usize::try_from(neighbor_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                adjacent_edge.flow = adjacent_edge.flow.wrapping_sub(1);
                            }
                            {
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [first_neighbor_vertex, second_neighbor_vertex]
                                {
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
                                && ((path_end == first_neighbor_vertex
                                    && path_start == second_neighbor_vertex)
                                    || (path_end == second_neighbor_vertex
                                        && path_start == first_neighbor_vertex))
                                && (delta_charge == 0 || delta_charge == 1)
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
                                    n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                    number_remove = number_remove.wrapping_sub(1);
                                    number_move_plus = number_move_plus.wrapping_sub(1);
                                    current_success = current_success.wrapping_add(1);
                                }
                            } else {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.flow = edge.flow.wrapping_sub(1);
                                for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                                    let adjacent_edge = heap
                                        .slice_mut(pBNS.edge)?
                                        .get_mut(
                                            usize::try_from(neighbor_edge)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    adjacent_edge.flow = adjacent_edge.flow.wrapping_add(1);
                                }
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [first_neighbor_vertex, second_neighbor_vertex]
                                {
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
                            ret =
                                AddToEdgeList(heap, &mut current_edges, edge_number, INC_ADD_EDGE)?;
                            if ret != 0 {
                                leave_case = true;
                                break;
                            }
                            endpoint_number = endpoint_number.wrapping_add(1);
                        }
                        if !leave_case {
                            RemoveForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &all_charge_edges,
                                forbidden_edge_mask,
                            )?;
                        }
                    }
                    group_index = group_index.wrapping_add(1);
                }
            }
            current_edges.num_edges = 0;
            let _ = AllocEdgeList(heap, &mut current_charge_edges, EDGE_LIST_FREE)?;
            let _ = AllocEdgeList(heap, &mut endpoint_list, EDGE_LIST_FREE)?;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 07, ichirvr3.c:2329-2469.
        if comparison.len_c2at >= 1 {
            let mut single_bond_oxygen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_oxygen_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_oxygen_minus = 0_i32;
            let mut number_double_bond_oxygen_neutral = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_double_bond_oxygen_neutral < MAX_DIFF_FIXH as i32
                    && difference.nValElectr == 6
                    && difference.endptInChI == 0
                    && edge_available
                    && difference.nFixHInChI == -1
                    && difference.nMobHInChI == 1
                    && difference.nFixHRevrs == 0
                    && difference.nMobHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
                    && atom.valence < atom.chem_bonds_valence
                {
                    double_bond_oxygen_neutral[number_double_bond_oxygen_neutral as usize] =
                        atom_number as i16;
                    number_double_bond_oxygen_neutral =
                        number_double_bond_oxygen_neutral.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_endpoint = !mobile_h_reversed.is_null()
                    && heap
                        .slice(mobile_h_reversed.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint
                        != 0;
                let fixed_h = !pStruct.fixed_H.is_null()
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let mobile_h = !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_single_bond_oxygen_minus < MAX_DIFF_FIXH as i32
                    && atom.charge == -1
                    && atom.num_H == 0
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 6
                    && reversed_endpoint
                    && !fixed_h
                    && !mobile_h
                    && edge_available
                {
                    single_bond_oxygen_minus[number_single_bond_oxygen_minus as usize] =
                        atom_number as i16;
                    number_single_bond_oxygen_minus =
                        number_single_bond_oxygen_minus.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            let number_to_try =
                number_single_bond_oxygen_minus.min(number_double_bond_oxygen_neutral);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut index = 0_i32;
                while index < number_single_bond_oxygen_minus && current_success < number_to_try {
                    let atom_number = i32::from(single_bond_oxygen_minus[index as usize]);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        edge.flow = edge.flow.wrapping_add(1);
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
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 07a, ichirvr3.c:2470-2612.
        if comparison.len_c2at >= 1 {
            let mut single_bond_oxygen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_oxygen_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_oxygen_minus = 0_i32;
            let mut number_double_bond_oxygen_neutral = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_double_bond_oxygen_neutral < MAX_DIFF_FIXH as i32
                    && difference.nValElectr == 6
                    && difference.endptInChI == 0
                    && edge_available
                    && difference.nFixHInChI == -1
                    && difference.nMobHInChI == 1
                    && difference.nFixHRevrs == 0
                    && difference.nMobHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
                    && atom.valence < atom.chem_bonds_valence
                {
                    double_bond_oxygen_neutral[number_double_bond_oxygen_neutral as usize] =
                        atom_number as i16;
                    number_double_bond_oxygen_neutral =
                        number_double_bond_oxygen_neutral.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let fixed_h = !pStruct.fixed_H.is_null()
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let mobile_h = !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                let mut nitrogen_five_neighbor = false;
                if number_single_bond_oxygen_minus < MAX_DIFF_FIXH as i32
                    && atom.charge == -1
                    && atom.num_H == 0
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 6
                    && endpoint == 0
                    && !fixed_h
                    && !mobile_h
                    && atom.valence == 1
                {
                    let nitrogen_number = i32::from(atom.neighbor[0]);
                    let nitrogen_index = usize::try_from(nitrogen_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let nitrogen = heap
                        .slice(at2.as_const())?
                        .get(nitrogen_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    nitrogen_five_neighbor = nitrogen.chem_bonds_valence == 5
                        && nitrogen.charge == 0
                        && pVA
                            .get(nitrogen_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .cNumValenceElectrons
                            == 5;
                }
                if number_single_bond_oxygen_minus < MAX_DIFF_FIXH as i32
                    && atom.charge == -1
                    && atom.num_H == 0
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 6
                    && endpoint == 0
                    && !fixed_h
                    && !mobile_h
                    && atom.valence == 1
                    && nitrogen_five_neighbor
                    && edge_available
                {
                    single_bond_oxygen_minus[number_single_bond_oxygen_minus as usize] =
                        atom_number as i16;
                    number_single_bond_oxygen_minus =
                        number_single_bond_oxygen_minus.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            let number_to_try =
                number_single_bond_oxygen_minus.min(number_double_bond_oxygen_neutral);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut index = 0_i32;
                while index < number_single_bond_oxygen_minus && current_success < number_to_try {
                    let atom_index =
                        usize::try_from(i32::from(single_bond_oxygen_minus[index as usize]))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        edge.flow = edge.flow.wrapping_add(1);
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
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 08, ichirvr3.c:2609-2762.
        if comparison.nNumTgInChI == 1
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI || comparison.nNumTgRevrs > 1)
        {
            let mut double_bond_nitrogen_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut single_bond_oxygen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_nitrogen_neutral = 0_i32;
            let mut number_single_bond_oxygen_minus = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let fixed_h = !pStruct.fixed_H.is_null()
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let mobile_h = !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let edge_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_single_bond_oxygen_minus < MAX_DIFF_FIXH as i32
                    && atom.charge == -1
                    && atom.num_H == 0
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 6
                    && endpoint != 0
                    && !mobile_h
                    && edge_available
                {
                    single_bond_oxygen_minus[number_single_bond_oxygen_minus as usize] =
                        atom_number as i16;
                    number_single_bond_oxygen_minus =
                        number_single_bond_oxygen_minus.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                } else if number_double_bond_nitrogen_neutral < MAX_DIFF_FIXH as i32
                    && atom.charge == 0
                    && atom.num_H == 0
                    && atom.sb_parity[0] == 0
                    && atom.valence < atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 5
                    && valence.cPeriodicRowNumber == 1
                    && endpoint != 0
                    && !fixed_h
                    && !mobile_h
                    && edge_available
                {
                    double_bond_nitrogen_neutral[number_double_bond_nitrogen_neutral as usize] =
                        atom_number as i16;
                    number_double_bond_nitrogen_neutral =
                        number_double_bond_nitrogen_neutral.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            let number_to_try =
                number_double_bond_nitrogen_neutral.min(number_single_bond_oxygen_minus);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                if forbidden_stereo_edge_mask != 0 {
                    RemoveForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &fixed_large_ring_stereo_edges,
                        forbidden_stereo_edge_mask,
                    )?;
                }
                let mut index = 0_i32;
                while index < number_single_bond_oxygen_minus && current_success < number_to_try {
                    let atom_index =
                        usize::try_from(i32::from(single_bond_oxygen_minus[index as usize]))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        index = index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        edge.flow = edge.flow.wrapping_add(1);
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
                    index = index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                if forbidden_stereo_edge_mask != 0 {
                    SetForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &fixed_large_ring_stereo_edges,
                        forbidden_stereo_edge_mask,
                    )?;
                }
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 09, ichirvr3.c:2763-2919.
        if comparison.len_c2at >= 2 {
            let mut current_success = 0_i32;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCPlusGroupEdge
                    .wrapping_sub(1);
                let plus_edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if difference.nValElectr == 5
                    && difference.nPeriodNum == 1
                    && difference.endptInChI == 0
                    && difference.nFixHInChI == 0
                    && difference.nMobHInChI != 0
                    && difference.endptRevrs != 0
                    && difference.nFixHRevrs != 0
                    && difference.nMobHRevrs == 0
                    && difference.nAtChargeRevrs == 1
                    && i32::from(atom.valence) + 1 == i32::from(atom.chem_bonds_valence)
                    && plus_edge_available
                {
                    let mut bond_number = 0_i32;
                    while bond_number < i32::from(atom.valence) {
                        let bond_index = usize::try_from(bond_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if atom.bond_type[bond_index] == BOND_TYPE_DOUBLE as u8 {
                            break;
                        }
                        bond_number = bond_number.wrapping_add(1);
                    }
                    if bond_number == i32::from(atom.valence) {
                        difference_index = difference_index.wrapping_add(1);
                        continue;
                    }
                    let atom_vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let nitrogen_carbon_edge = *heap
                        .slice(atom_vertex.iedge.as_const())?
                        .get(
                            usize::try_from(bond_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let carbon_number = i32::from(
                        *atom
                            .neighbor
                            .get(
                                usize::try_from(bond_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let carbon_index = usize::try_from(carbon_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let carbon = heap
                        .slice(at2.as_const())?
                        .get(carbon_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let carbon_valence = pVA
                        .get(carbon_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let carbon_plus_edge = carbon_valence.nCPlusGroupEdge.wrapping_sub(1);
                    let carbon_plus_available = carbon_plus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(carbon_plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .forbidden
                            == 0;
                    if carbon_valence.cNumValenceElectrons != 4
                        || carbon_valence.cMetal != 0
                        || carbon.charge != 0
                        || carbon.valence != 3
                        || i32::from(carbon.valence) + 1 != i32::from(carbon.chem_bonds_valence)
                        || !carbon_plus_available
                    {
                        difference_index = difference_index.wrapping_add(1);
                        continue;
                    }
                    let mut carbon_neighbor_number = 0_i32;
                    while carbon_neighbor_number < i32::from(carbon.valence) {
                        let carbon_neighbor_index = usize::try_from(carbon_neighbor_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let second_nitrogen_number =
                            i32::from(carbon.neighbor[carbon_neighbor_index]);
                        let second_nitrogen_index = usize::try_from(second_nitrogen_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let second_nitrogen = heap
                            .slice(at2.as_const())?
                            .get(second_nitrogen_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let second_valence = pVA
                            .get(second_nitrogen_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if second_nitrogen_number == atom_number
                            || second_valence.cNumValenceElectrons != 5
                            || second_valence.cPeriodicRowNumber != 1
                            || second_nitrogen.num_H == 0
                            || second_nitrogen.charge != 0
                        {
                            carbon_neighbor_number = carbon_neighbor_number.wrapping_add(1);
                            continue;
                        }
                        let mut second_difference_index = 0_i32;
                        while second_difference_index < i32::from(comparison.len_c2at)
                            && second_nitrogen_number
                                != i32::from(
                                    comparison.c2at[usize::try_from(second_difference_index)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                                    .atomNumber,
                                )
                        {
                            second_difference_index = second_difference_index.wrapping_add(1);
                        }
                        if second_difference_index == i32::from(comparison.len_c2at) {
                            carbon_neighbor_number = carbon_neighbor_number.wrapping_add(1);
                            continue;
                        }
                        let second_difference =
                            comparison.c2at[usize::try_from(second_difference_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                            .clone();
                        if second_difference.nValElectr == 5
                            && second_difference.nPeriodNum == 1
                            && second_difference.endptInChI == 0
                            && second_difference.nFixHInChI == 0
                            && second_difference.nMobHInChI != 0
                            && second_difference.endptRevrs != 0
                            && second_difference.nFixHRevrs != 0
                            && second_difference.nMobHRevrs == 0
                            && second_difference.nAtChargeRevrs == 0
                            && second_nitrogen.valence == second_nitrogen.chem_bonds_valence
                        {
                            ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                            ret = AddToEdgeList(
                                heap,
                                &mut current_edges,
                                carbon_plus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                return Ok(());
                            }
                            SetForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &all_charge_edges,
                                forbidden_edge_mask,
                            )?;
                            RemoveForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &current_edges,
                                forbidden_edge_mask,
                            )?;
                            let edge_index = usize::try_from(nitrogen_carbon_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let edge_before = heap
                                .slice(pBNS.edge.as_const())?
                                .get(edge_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge_before.flow == 0 {
                                // Source `continue` targets the carbon-neighbor loop and
                                // deliberately skips mask removal and CurrEdges clearing.
                                carbon_neighbor_number = carbon_neighbor_number.wrapping_add(1);
                                continue;
                            }
                            let first_vertex = i32::from(edge_before.neighbor1);
                            let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                            {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden =
                                    (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
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
                                    n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                    current_success = current_success.wrapping_add(1);
                                }
                            } else {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.flow = edge.flow.wrapping_add(1);
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
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            RemoveForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &all_charge_edges,
                                forbidden_edge_mask,
                            )?;
                            current_edges.num_edges = 0;
                            break;
                        }
                        carbon_neighbor_number = carbon_neighbor_number.wrapping_add(1);
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 10, ichirvr3.c:2920-3092.
        if comparison.len_c2at >= 2 {
            let mut current_success = 0_i32;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let comparison_index = usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let difference = comparison.c2at[comparison_index].clone();
                if difference.nValue != 0 {
                    difference_index = difference_index.wrapping_add(1);
                    continue;
                }
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atom_valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atom_plus_edge = atom_valence.nCPlusGroupEdge.wrapping_sub(1);
                let atom_plus_available = atom_plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(atom_plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if (difference.nValElectr == 6
                    || (difference.nValElectr == 5 && difference.nPeriodNum == 1))
                    && difference.endptInChI != 0
                    && difference.nFixHInChI != 0
                    && difference.nMobHInChI == 0
                    && difference.endptRevrs == 0
                    && difference.nFixHRevrs == 0
                    && difference.nMobHRevrs != 0
                    && difference.nAtChargeRevrs == 0
                    && atom.valence == atom.chem_bonds_valence
                    && atom_plus_available
                {
                    let atom_vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut keep_searching = true;
                    let mut first_bond_number = 0_i32;
                    while first_bond_number < i32::from(atom.valence) && keep_searching {
                        let first_bond_index = usize::try_from(first_bond_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if atom.bond_type[first_bond_index] == BOND_TYPE_SINGLE as u8 {
                            let center_number = i32::from(atom.neighbor[first_bond_index]);
                            let center_index = usize::try_from(center_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let center = heap
                                .slice(at2.as_const())?
                                .get(center_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let center_valence = pVA
                                .get(center_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let center_plus_edge = center_valence.nCPlusGroupEdge.wrapping_sub(1);
                            let center_plus_available = center_plus_edge >= 0
                                && heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(center_plus_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .forbidden
                                    == 0;
                            if center.charge == 1
                                && (4..=6).contains(&center_valence.cNumValenceElectrons)
                                && center.valence == center.chem_bonds_valence
                                && center_plus_available
                            {
                                let center_vertex = heap
                                    .slice(pBNS.vert.as_const())?
                                    .get(center_index)
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let mut second_bond_number = 0_i32;
                                while second_bond_number < i32::from(center.valence)
                                    && keep_searching
                                {
                                    let second_bond_index = usize::try_from(second_bond_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                    if center.bond_type[second_bond_index] == BOND_TYPE_SINGLE as u8
                                    {
                                        let second_atom_number =
                                            i32::from(center.neighbor[second_bond_index]);
                                        let second_atom_index = usize::try_from(second_atom_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                        let second_atom = heap
                                            .slice(at2.as_const())?
                                            .get(second_atom_index)
                                            .cloned()
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        let second_valence = pVA
                                            .get(second_atom_index)
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        let second_plus_edge =
                                            second_valence.nCPlusGroupEdge.wrapping_sub(1);
                                        let second_plus_available = second_plus_edge >= 0
                                            && heap
                                                .slice(pBNS.edge.as_const())?
                                                .get(usize::try_from(second_plus_edge).map_err(
                                                    |_| SourceHeapError::PointerOutOfBounds,
                                                )?)
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                                .forbidden
                                                == 0;
                                        if second_atom_number != atom_number
                                            && second_atom.charge == 0
                                            && second_atom.num_H != 0
                                            && (second_valence.cNumValenceElectrons == 5
                                                || second_valence.cNumValenceElectrons == 6)
                                            && second_atom.valence == second_atom.chem_bonds_valence
                                            && second_plus_available
                                        {
                                            let mut second_difference_index = 0_i32;
                                            while second_difference_index
                                                < i32::from(comparison.len_c2at)
                                            {
                                                let second_comparison_index =
                                                    usize::try_from(second_difference_index)
                                                        .map_err(|_| {
                                                            SourceHeapError::PointerOutOfBounds
                                                        })?;
                                                let second_difference = comparison.c2at
                                                    [second_comparison_index]
                                                    .clone();
                                                if second_atom_number
                                                    != i32::from(second_difference.atomNumber)
                                                    || second_difference.nValue != 0
                                                {
                                                    second_difference_index =
                                                        second_difference_index.wrapping_add(1);
                                                    continue;
                                                }
                                                if second_difference.endptInChI != 0
                                                    && second_difference.nFixHInChI != 0
                                                    && second_difference.nMobHInChI == 0
                                                    && second_difference.endptRevrs == 0
                                                    && second_difference.nFixHRevrs == 0
                                                    && second_difference.nMobHRevrs != 0
                                                    && second_difference.nAtChargeRevrs == 0
                                                {
                                                    let first_bond_edge = *heap
                                                        .slice(atom_vertex.iedge.as_const())?
                                                        .get(first_bond_index)
                                                        .ok_or(
                                                            SourceHeapError::PointerOutOfBounds,
                                                        )?;
                                                    let second_bond_edge = *heap
                                                        .slice(center_vertex.iedge.as_const())?
                                                        .get(second_bond_index)
                                                        .ok_or(
                                                            SourceHeapError::PointerOutOfBounds,
                                                        )?;
                                                    let first_bond_forbidden = heap
                                                        .slice(pBNS.edge.as_const())?
                                                        .get(
                                                            usize::try_from(first_bond_edge)
                                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                        )
                                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                                        .forbidden
                                                        != 0;
                                                    let first_charge_flow = heap
                                                        .slice(pBNS.edge.as_const())?
                                                        .get(
                                                            usize::try_from(atom_plus_edge)
                                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                        )
                                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                                        .flow;
                                                    let second_bond_forbidden = heap
                                                        .slice(pBNS.edge.as_const())?
                                                        .get(
                                                            usize::try_from(second_bond_edge)
                                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                        )
                                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                                        .forbidden
                                                        != 0;
                                                    let second_charge_flow = heap
                                                        .slice(pBNS.edge.as_const())?
                                                        .get(
                                                            usize::try_from(second_plus_edge)
                                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                        )
                                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                                        .flow;
                                                    let selected_edge = if !first_bond_forbidden
                                                        && first_charge_flow != 0
                                                    {
                                                        Some(atom_plus_edge)
                                                    } else if !second_bond_forbidden
                                                        && second_charge_flow != 0
                                                    {
                                                        Some(second_plus_edge)
                                                    } else {
                                                        None
                                                    };
                                                    let Some(selected_edge) = selected_edge else {
                                                        second_difference_index =
                                                            second_difference_index.wrapping_add(1);
                                                        continue;
                                                    };
                                                    ret = AddToEdgeList(
                                                        heap,
                                                        &mut current_edges,
                                                        center_plus_edge,
                                                        INC_ADD_EDGE,
                                                    )?;
                                                    if ret != 0 {
                                                        return Ok(());
                                                    }
                                                    SetForbiddenEdgeMask(
                                                        heap,
                                                        pBNS,
                                                        &all_charge_edges,
                                                        forbidden_edge_mask,
                                                    )?;
                                                    RemoveForbiddenEdgeMask(
                                                        heap,
                                                        pBNS,
                                                        &current_edges,
                                                        forbidden_edge_mask,
                                                    )?;
                                                    let selected_index =
                                                        usize::try_from(selected_edge).map_err(
                                                            |_| SourceHeapError::PointerOutOfBounds,
                                                        )?;
                                                    let selected_before = heap
                                                        .slice(pBNS.edge.as_const())?
                                                        .get(selected_index)
                                                        .cloned()
                                                        .ok_or(
                                                            SourceHeapError::PointerOutOfBounds,
                                                        )?;
                                                    let first_vertex =
                                                        i32::from(selected_before.neighbor1);
                                                    let second_vertex =
                                                        i32::from(selected_before.neighbor12)
                                                            ^ first_vertex;
                                                    heap.slice_mut(pBNS.edge)?[selected_index]
                                                        .flow =
                                                        selected_before.flow.wrapping_sub(1);
                                                    {
                                                        let vertices = heap.slice_mut(pBNS.vert)?;
                                                        for vertex_number in
                                                            [first_vertex, second_vertex]
                                                        {
                                                            let vertex =
                                                                vertices
                                                                    .get_mut(usize::try_from(vertex_number).map_err(
                                                                        |_| SourceHeapError::PointerOutOfBounds,
                                                                    )?)
                                                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                                            vertex.st_edge.flow =
                                                                vertex.st_edge.flow.wrapping_sub(1);
                                                        }
                                                    }
                                                    pBNS.tot_st_flow =
                                                        pBNS.tot_st_flow.wrapping_sub(2);
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
                                                        && ((path_end == first_vertex
                                                            && path_start == second_vertex)
                                                            || (path_end == second_vertex
                                                                && path_start == first_vertex))
                                                        && delta_charge == -1
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
                                                            n_num_run_bns =
                                                                n_num_run_bns.wrapping_add(1);
                                                            current_success =
                                                                current_success.wrapping_add(1);
                                                            keep_searching = false;
                                                            comparison.c2at[comparison_index]
                                                                .nValue = 1;
                                                            comparison.c2at
                                                                [second_comparison_index]
                                                                .nValue = 1;
                                                        }
                                                    } else {
                                                        let selected = heap
                                                            .slice_mut(pBNS.edge)?
                                                            .get_mut(selected_index)
                                                            .ok_or(
                                                                SourceHeapError::PointerOutOfBounds,
                                                            )?;
                                                        selected.flow =
                                                            selected.flow.wrapping_add(1);
                                                        let vertices = heap.slice_mut(pBNS.vert)?;
                                                        for vertex_number in
                                                            [first_vertex, second_vertex]
                                                        {
                                                            let vertex =
                                                                vertices
                                                                    .get_mut(usize::try_from(vertex_number).map_err(
                                                                        |_| SourceHeapError::PointerOutOfBounds,
                                                                    )?)
                                                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                                            vertex.st_edge.flow =
                                                                vertex.st_edge.flow.wrapping_add(1);
                                                        }
                                                        pBNS.tot_st_flow =
                                                            pBNS.tot_st_flow.wrapping_add(2);
                                                    }
                                                    RemoveForbiddenEdgeMask(
                                                        heap,
                                                        pBNS,
                                                        &all_charge_edges,
                                                        forbidden_edge_mask,
                                                    )?;
                                                    current_edges.num_edges = 0;
                                                    break;
                                                }
                                                second_difference_index =
                                                    second_difference_index.wrapping_add(1);
                                            }
                                        }
                                    }
                                    second_bond_number = second_bond_number.wrapping_add(1);
                                }
                            }
                        }
                        first_bond_number = first_bond_number.wrapping_add(1);
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 11, ichirvr3.c:3089-3303.
        if comparison.nNumTgInChI == 1
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI || comparison.nNumTgRevrs > 1)
        {
            let mut single_bond_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_charged = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_neutral = 0_i32;
            let mut number_double_bond_charged = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                let fixed_h = !pStruct.fixed_H.is_null()
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let mobile_h = !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                if number_double_bond_charged < MAX_DIFF_FIXH as i32
                    && atom.charge == 1
                    && atom.num_H != 0
                    && atom.valence < atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && ((valence.cNumValenceElectrons == 5 && valence.cPeriodicRowNumber == 1)
                        || valence.cNumValenceElectrons == 6)
                    && fixed_h
                    && edge_available
                {
                    double_bond_charged[number_double_bond_charged as usize] = atom_number as i16;
                    number_double_bond_charged = number_double_bond_charged.wrapping_add(1);
                } else if number_single_bond_neutral < MAX_DIFF_FIXH as i32
                    && atom.charge == 0
                    && atom.num_H == 0
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 5
                    && valence.cPeriodicRowNumber == 1
                    && !fixed_h
                    && !mobile_h
                    && edge_available
                {
                    single_bond_neutral[number_single_bond_neutral as usize] = atom_number as i16;
                    number_single_bond_neutral = number_single_bond_neutral.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            let number_to_try = number_single_bond_neutral.min(number_double_bond_charged);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < number_double_bond_charged
                    && current_success < number_to_try
                {
                    let atom_index =
                        usize::try_from(i32::from(double_bond_charged[candidate_index as usize]))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow != 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    let first_vertex_data = heap
                        .slice(pBNS.vert.as_const())?
                        .get(
                            usize::try_from(first_vertex)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut adjacency_index = i32::from(first_vertex_data.num_adj_edges) - 1;
                    let mut first_neighbor = None;
                    while adjacency_index >= 0 {
                        let adjacent_edge_number = *heap
                            .slice(first_vertex_data.iedge.as_const())?
                            .get(
                                usize::try_from(adjacency_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(adjacent_edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                            first_neighbor = Some((
                                adjacent_edge_number,
                                i32::from(adjacent_edge.neighbor12) ^ first_vertex,
                            ));
                            break;
                        }
                        adjacency_index -= 1;
                    }
                    let Some((first_neighbor_edge, first_neighbor_vertex)) = first_neighbor else {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    };
                    let second_vertex_data = heap
                        .slice(pBNS.vert.as_const())?
                        .get(
                            usize::try_from(second_vertex)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut adjacency_index = i32::from(second_vertex_data.num_adj_edges) - 2;
                    let mut second_neighbor = None;
                    while adjacency_index >= 0 {
                        let adjacent_edge_number = *heap
                            .slice(second_vertex_data.iedge.as_const())?
                            .get(
                                usize::try_from(adjacency_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(adjacent_edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                            second_neighbor = Some((
                                adjacent_edge_number,
                                i32::from(adjacent_edge.neighbor12) ^ second_vertex,
                            ));
                            break;
                        }
                        adjacency_index -= 1;
                    }
                    let Some((second_neighbor_edge, second_neighbor_vertex)) = second_neighbor
                    else {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    };
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow.wrapping_add(1);
                    for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                        let adjacent_edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(
                                usize::try_from(neighbor_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        adjacent_edge.flow = adjacent_edge.flow.wrapping_sub(1);
                    }
                    {
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_neighbor_vertex, second_neighbor_vertex] {
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
                        && ((path_end == first_neighbor_vertex
                            && path_start == second_neighbor_vertex)
                            || (path_end == second_neighbor_vertex
                                && path_start == first_neighbor_vertex))
                        && (delta_charge == 0 || delta_charge == 1)
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.flow = edge.flow.wrapping_sub(1);
                        for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                            let adjacent_edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(
                                    usize::try_from(neighbor_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            adjacent_edge.flow = adjacent_edge.flow.wrapping_add(1);
                        }
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_neighbor_vertex, second_neighbor_vertex] {
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
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 12, ichirvr3.c:3299-3518.
        if comparison.len_c2at >= 1
            && comparison.nNumTgInChI == 1
            && comparison.nNumRemHInChI >= -1
            && (comparison.nNumEndpInChI > comparison.nNumEndpRevrs || comparison.nNumTgRevrs > 1)
        {
            let mut only_nitrogen_five = true;
            loop {
                let mut single_bond_nitrogen_neutral = [0_i16; MAX_DIFF_FIXH as usize];
                let mut double_bond_hetero = [0_i16; MAX_DIFF_FIXH as usize];
                let mut number_single_bond_nitrogen_neutral = 0_i32;
                let mut number_double_bond_hetero = 0_i32;
                let mut number_nitrogen_five = 0_i32;
                let canonical_to_atom = pStruct.nCanon2Atno[0];
                let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                    SourceMutPointer::null()
                } else {
                    heap.slice(pStruct.pOne_norm_data[1].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .at
                };
                let mobile_h_input = if !pInChI[1].is_null() {
                    let inchi = heap
                        .slice(pInChI[1].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if !inchi.nNum_H.is_null() {
                        inchi.nNum_H
                    } else if !pInChI[0].is_null() {
                        heap.slice(pInChI[0].as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nNum_H
                    } else {
                        SourceMutPointer::null()
                    }
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                };
                let mut current_success = 0_i32;

                let mut difference_index = 0_i32;
                while difference_index < i32::from(comparison.len_c2at) {
                    let difference = comparison.c2at[usize::try_from(difference_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                    .clone();
                    let atom_number = i32::from(difference.atomNumber);
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let minus_edge = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_available = minus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .forbidden
                            == 0;
                    let source_valence_matches = if difference.nValElectr == 6 {
                        atom.valence == 1 && atom.chem_bonds_valence == 2
                    } else if difference.nValElectr == 5 {
                        atom.valence == 2 && atom.chem_bonds_valence == 3
                    } else {
                        false
                    };
                    if number_double_bond_hetero < MAX_DIFF_FIXH as i32
                        && (difference.nValElectr == 6
                            || (difference.nValElectr == 5 && difference.nPeriodNum == 1))
                        && difference.endptInChI != 0
                        && edge_available
                        && difference.nFixHInChI == 0
                        && difference.nMobHInChI == 0
                        && difference.endptRevrs == 0
                        && difference.nFixHRevrs == 0
                        && difference.nMobHRevrs == 0
                        && difference.nAtChargeRevrs == 0
                        && atom.num_H == 0
                        && source_valence_matches
                    {
                        double_bond_hetero[number_double_bond_hetero as usize] = atom_number as i16;
                        number_double_bond_hetero = number_double_bond_hetero.wrapping_add(1);
                    }
                    difference_index = difference_index.wrapping_add(1);
                }

                let mut canonical_number = 0_i32;
                while canonical_number < pStruct.num_atoms {
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut nitrogen_five = false;
                    let valence_matches = atom.valence == atom.chem_bonds_valence || {
                        nitrogen_five =
                            i32::from(atom.valence) + 2 == i32::from(atom.chem_bonds_valence);
                        nitrogen_five
                    };
                    let reversed_endpoint = !mobile_h_reversed.is_null()
                        && heap
                            .slice(mobile_h_reversed.as_const())?
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .endpoint
                            != 0;
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                    let fixed_h = !pStruct.fixed_H.is_null()
                        && *heap
                            .slice(pStruct.fixed_H.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            != 0;
                    let mobile_h = !mobile_h_input.is_null()
                        && *heap
                            .slice(mobile_h_input.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            != 0;
                    let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                    let edge_available = plus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .forbidden
                            == 0;
                    if number_single_bond_nitrogen_neutral < MAX_DIFF_FIXH as i32
                        && atom.charge == 0
                        && atom.num_H == 0
                        && valence_matches
                        && valence.cMetal == 0
                        && valence.cNumValenceElectrons == 5
                        && valence.cPeriodicRowNumber == 1
                        && !reversed_endpoint
                        && !input_endpoint
                        && !fixed_h
                        && !mobile_h
                        && edge_available
                    {
                        let upper_edge = if nitrogen_five {
                            GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?
                        } else {
                            NO_VERTEX
                        };
                        let upper_available = upper_edge != NO_VERTEX && {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(upper_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden == 0 && edge.flow == 0
                        };
                        if only_nitrogen_five && nitrogen_five && upper_available {
                            if number_nitrogen_five == 0 {
                                current_edges.num_edges = 0;
                                number_single_bond_nitrogen_neutral = 0;
                            }
                            single_bond_nitrogen_neutral
                                [number_single_bond_nitrogen_neutral as usize] = atom_number as i16;
                            number_single_bond_nitrogen_neutral =
                                number_single_bond_nitrogen_neutral.wrapping_add(1);
                            number_nitrogen_five = number_nitrogen_five.wrapping_add(1);
                            ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                            ret =
                                AddToEdgeList(heap, &mut current_edges, upper_edge, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                        } else if number_nitrogen_five == 0 {
                            single_bond_nitrogen_neutral
                                [number_single_bond_nitrogen_neutral as usize] = atom_number as i16;
                            number_single_bond_nitrogen_neutral =
                                number_single_bond_nitrogen_neutral.wrapping_add(1);
                            ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                            if nitrogen_five && upper_available {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut current_edges,
                                    upper_edge,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    }
                    canonical_number = canonical_number.wrapping_add(1);
                }

                let number_to_try =
                    number_single_bond_nitrogen_neutral.min(number_double_bond_hetero);
                if number_to_try != 0 {
                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                    let mut donor_index = 0_i32;
                    while donor_index < number_double_bond_hetero && current_success < number_to_try
                    {
                        let atom_index =
                            usize::try_from(i32::from(double_bond_hetero[donor_index as usize]))
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let minus_edge = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCMinusGroupEdge
                            .wrapping_sub(1);
                        let minus_edge_index = usize::try_from(minus_edge)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(minus_edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        }
                        let atom_vertex = heap
                            .slice(pBNS.vert.as_const())?
                            .get(atom_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let bond_edge_number = *heap
                            .slice(atom_vertex.iedge.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let bond_edge_index = usize::try_from(bond_edge_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let bond_before = heap
                            .slice(pBNS.edge.as_const())?
                            .get(bond_edge_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if bond_before.flow == 0 {
                            donor_index = donor_index.wrapping_add(1);
                            continue;
                        }
                        let first_vertex = i32::from(bond_before.neighbor1);
                        let second_vertex = i32::from(bond_before.neighbor12) ^ first_vertex;
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(bond_edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
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
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                current_success = current_success.wrapping_add(1);
                            }
                        } else {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(bond_edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.flow = edge.flow.wrapping_add(1);
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
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(bond_edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        donor_index = donor_index.wrapping_add(1);
                    }
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                }
                current_edges.num_edges = 0;
                if current_success != 0 {
                    tot_success = tot_success.wrapping_add(current_success);
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
                    ret = FillOutCMP2FHINCHI(
                        heap,
                        pStruct,
                        &at2_snapshot,
                        pVA,
                        pInChI,
                        &mut comparison,
                    )?;
                    if ret != 0 || comparison.bHasDifference == 0 {
                        return Ok(());
                    }
                    break;
                } else if only_nitrogen_five {
                    only_nitrogen_five = false;
                } else {
                    break;
                }
            }
        }

        // Complete translation of source case 13, ichirvr3.c:3515-3700.
        if comparison.nNumTgDiffMinus != 0 {
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed_atoms = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let fixed_h_reversed = heap
                .slice(pStruct.pOneINChI[0].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .nNum_H_fixed;
            let mobile_h_reversed = if !pStruct.pOneINChI[1].is_null() {
                let inchi = heap
                    .slice(pStruct.pOneINChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pStruct.pOneINChI[0].is_null() {
                    heap.slice(pStruct.pOneINChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pStruct.pOneINChI[0].is_null() {
                heap.slice(pStruct.pOneINChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };

            let mut tgroup_index = 0_i32;
            while tgroup_index < pStruct.ti.num_t_groups
                && tgroup_index < pStruct.One_ti.num_t_groups
            {
                let group_index = usize::try_from(tgroup_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let input_group = heap
                    .slice(pStruct.ti.t_group.as_const())?
                    .get(group_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_group = heap
                    .slice(pStruct.One_ti.t_group.as_const())?
                    .get(group_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if input_group.nNumEndpoints == reversed_group.nNumEndpoints
                    && i32::from(input_group.num[0]) - i32::from(input_group.num[1])
                        == i32::from(reversed_group.num[0]) - i32::from(reversed_group.num[1])
                    && input_group.num[1] > reversed_group.num[1]
                {
                    let mut single_bond_nitrogen_neutral = [0_i16; MAX_DIFF_FIXH as usize];
                    let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
                    let mut number_single_bond_nitrogen_neutral = 0_i32;
                    let mut number_double_bond_oxygen = 0_i32;
                    let mut current_success = 0_i32;
                    let mut canonical_number = 0_i32;
                    while canonical_number < pStruct.num_atoms {
                        let canonical_index = usize::try_from(canonical_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let atom_number = i32::from(
                            *heap
                                .slice(canonical_to_atom.as_const())?
                                .get(canonical_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        let atom_index = usize::try_from(atom_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let atom = heap
                            .slice(at2.as_const())?
                            .get(atom_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let canonical_valence = pVA
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let atom_valence = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let input_endpoint = *heap
                            .slice(pStruct.endpoint.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let input_fixed_h = !pStruct.fixed_H.is_null()
                            && *heap
                                .slice(pStruct.fixed_H.as_const())?
                                .get(canonical_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                != 0;
                        let input_mobile_h = !mobile_h_input.is_null()
                            && *heap
                                .slice(mobile_h_input.as_const())?
                                .get(canonical_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                != 0;
                        let reversed_fixed = !fixed_h_reversed.is_null()
                            && *heap
                                .slice(fixed_h_reversed.as_const())?
                                .get(atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                != 0;
                        let reversed_mobile = !mobile_h_reversed.is_null()
                            && *heap
                                .slice(mobile_h_reversed.as_const())?
                                .get(atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                != 0;
                        let canonical_minus_edge =
                            canonical_valence.nCMinusGroupEdge.wrapping_sub(1);
                        let canonical_minus_available = canonical_minus_edge >= 0
                            && heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(canonical_minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .forbidden
                                == 0;
                        if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                            && canonical_valence.cNumValenceElectrons == 6
                            && i32::from(input_endpoint) == tgroup_index.wrapping_add(1)
                            && canonical_minus_available
                            && !input_fixed_h
                            && !input_mobile_h
                            && !reversed_fixed
                            && !reversed_mobile
                            && atom.charge == 0
                            && atom.num_H == 0
                            && atom.valence == 1
                            && atom.chem_bonds_valence == 2
                        {
                            double_bond_oxygen[number_double_bond_oxygen as usize] =
                                atom_number as i16;
                            number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                        } else {
                            let reversed_endpoint = !mobile_h_reversed_atoms.is_null()
                                && heap
                                    .slice(mobile_h_reversed_atoms.as_const())?
                                    .get(atom_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .endpoint
                                    != 0;
                            let plus_edge = atom_valence.nCPlusGroupEdge.wrapping_sub(1);
                            let plus_available = plus_edge >= 0
                                && heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(plus_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .forbidden
                                    == 0;
                            if number_single_bond_nitrogen_neutral < MAX_DIFF_FIXH as i32
                                && atom.charge == 0
                                && atom.num_H == 0
                                && atom.valence == 4
                                && atom.chem_bonds_valence == 5
                                && atom_valence.cMetal == 0
                                && atom_valence.cNumValenceElectrons == 5
                                && atom_valence.cPeriodicRowNumber >= 1
                                && !reversed_endpoint
                                && input_endpoint == 0
                                && !input_fixed_h
                                && !input_mobile_h
                                && plus_available
                            {
                                single_bond_nitrogen_neutral
                                    [number_single_bond_nitrogen_neutral as usize] =
                                    atom_number as i16;
                                number_single_bond_nitrogen_neutral =
                                    number_single_bond_nitrogen_neutral.wrapping_add(1);
                                ret = AddToEdgeList(
                                    heap,
                                    &mut current_edges,
                                    plus_edge,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                        canonical_number = canonical_number.wrapping_add(1);
                    }
                    let number_to_try =
                        number_single_bond_nitrogen_neutral.min(number_double_bond_oxygen);
                    if number_to_try != 0 {
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                        let mut donor_index = 0_i32;
                        while donor_index < number_double_bond_oxygen
                            && current_success < number_to_try
                        {
                            let atom_index = usize::try_from(i32::from(
                                double_bond_oxygen[donor_index as usize],
                            ))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let minus_edge = pVA
                                .get(atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nCMinusGroupEdge
                                .wrapping_sub(1);
                            let minus_edge_index = usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(minus_edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden =
                                    (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            }
                            let atom_vertex = heap
                                .slice(pBNS.vert.as_const())?
                                .get(atom_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let bond_edge_number = *heap
                                .slice(atom_vertex.iedge.as_const())?
                                .first()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let bond_edge_index = usize::try_from(bond_edge_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let bond_before = heap
                                .slice(pBNS.edge.as_const())?
                                .get(bond_edge_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if bond_before.flow == 0 {
                                donor_index = donor_index.wrapping_add(1);
                                continue;
                            }
                            let first_vertex = i32::from(bond_before.neighbor1);
                            let second_vertex = i32::from(bond_before.neighbor12) ^ first_vertex;
                            {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(bond_edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden =
                                    (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
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
                                ret = RunBnsRestoreOnce(
                                    heap,
                                    pBNS,
                                    pBD,
                                    pVA,
                                    pTCGroups,
                                    clock_result,
                                )?;
                                if ret > 0 {
                                    n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                    current_success = current_success.wrapping_add(1);
                                }
                            } else {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(bond_edge_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.flow = edge.flow.wrapping_add(1);
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
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(bond_edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            donor_index = donor_index.wrapping_add(1);
                        }
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &all_charge_edges,
                            forbidden_edge_mask,
                        )?;
                    }
                    current_edges.num_edges = 0;
                    if current_success != 0 {
                        tot_success = tot_success.wrapping_add(current_success);
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
                        ret = FillOutCMP2FHINCHI(
                            heap,
                            pStruct,
                            &at2_snapshot,
                            pVA,
                            pInChI,
                            &mut comparison,
                        )?;
                        if ret != 0 || comparison.bHasDifference == 0 {
                            return Ok(());
                        }
                        break;
                    }
                }
                tgroup_index = tgroup_index.wrapping_add(1);
            }
        }

        // Complete translation of source case 14, ichirvr3.c:3696-4032.
        let case_14_guard = (comparison.nNumTgInChI <= 1
            && comparison.nNumRemHInChI > comparison.nNumRemHRevrs)
            || comparison.len_c2at != 0;
        let has_nitrogen_five = if case_14_guard {
            let at2_snapshot = heap.slice(at2.as_const())?.to_vec();
            bHas_N_V(&at2_snapshot, pStruct.num_atoms)? != 0
        } else {
            false
        };
        if case_14_guard && has_nitrogen_five {
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut nitrogen_five_atoms = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_nitrogen_five = 0_i32;
            let mut max_success =
                i32::from(comparison.nNumRemHInChI) - i32::from(comparison.nNumRemHRevrs);
            let mut x_atoms = EDGE_LIST::default();
            let mut nitrogen_three_atoms = EDGE_LIST::default();
            let _ = AllocEdgeList(heap, &mut x_atoms, EDGE_LIST_CLEAR)?;
            let _ = AllocEdgeList(heap, &mut nitrogen_three_atoms, EDGE_LIST_CLEAR)?;
            let mut current_success = 0_i32;
            ret = 0;
            let mut leave_case = false;

            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms && !leave_case {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_endpoint = !mobile_h_reversed.is_null()
                    && heap
                        .slice(mobile_h_reversed.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint
                        != 0;
                let input_mobile_h = !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let input_endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let fixed_h_pointer_present = !pStruct.fixed_H.is_null();
                let input_fixed_h = fixed_h_pointer_present
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let plus_state = if plus_edge >= 0 {
                    Some(
                        heap.slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone(),
                    )
                } else {
                    None
                };
                let mut is_nitrogen_five = false;
                if number_nitrogen_five < MAX_DIFF_FIXH as i32
                    && atom.chem_bonds_valence == 5
                    && atom.valence == 3
                    && atom.charge == 0
                    && atom.radical == 0
                    && valence.cNumValenceElectrons == 5
                    && valence.cPeriodicRowNumber == 1
                    && !reversed_endpoint
                    && atom.num_H == 0
                    && plus_state
                        .as_ref()
                        .is_some_and(|edge| edge.forbidden == 0 && edge.flow != 0)
                {
                    let upper_edge = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                    let upper_available = upper_edge != NO_VERTEX && {
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(upper_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden == 0 && edge.flow == 0
                    };
                    is_nitrogen_five = upper_available
                        && input_endpoint == 0
                        && !input_mobile_h
                        && fixed_h_pointer_present
                        && !input_fixed_h;
                }
                if is_nitrogen_five {
                    nitrogen_five_atoms[number_nitrogen_five as usize] = atom_number as i16;
                    number_nitrogen_five = number_nitrogen_five.wrapping_add(1);
                } else {
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    let minus_state = if minus_edge >= 0 {
                        Some(
                            heap.slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .clone(),
                        )
                    } else {
                        None
                    };
                    let is_nitrogen_three = atom.chem_bonds_valence == 3
                        && atom.valence == 2
                        && atom.charge == 0
                        && atom.radical == 0
                        && valence.cNumValenceElectrons == 5
                        && valence.cPeriodicRowNumber == 1
                        && !reversed_endpoint
                        && atom.num_H == 0
                        && minus_state
                            .as_ref()
                            .is_some_and(|edge| edge.forbidden == 0 && edge.flow == 0)
                        && !input_mobile_h
                        && fixed_h_pointer_present
                        && !input_fixed_h;
                    if is_nitrogen_three {
                        ret = AddToEdgeList(heap, &mut nitrogen_three_atoms, atom_number, 32)?;
                        if ret != 0 {
                            leave_case = true;
                        }
                    } else {
                        let is_x = atom.chem_bonds_valence == atom.valence
                            && atom.charge == 0
                            && atom.radical == 0
                            && ((valence.cNumValenceElectrons == 5
                                && valence.cPeriodicRowNumber == 1)
                                || valence.cNumValenceElectrons == 6)
                            && atom.num_H != 0
                            && plus_state
                                .as_ref()
                                .is_some_and(|edge| edge.forbidden == 0 && edge.flow != 0)
                            && !input_mobile_h
                            && fixed_h_pointer_present
                            && input_fixed_h;
                        if is_x {
                            ret = AddToEdgeList(heap, &mut x_atoms, atom_number, 32)?;
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            if max_success == 0 {
                max_success = number_nitrogen_five.min(nitrogen_three_atoms.num_edges);
                max_success = max_success.min(x_atoms.num_edges);
            }
            if !leave_case
                && number_nitrogen_five != 0
                && nitrogen_three_atoms.num_edges != 0
                && x_atoms.num_edges != 0
            {
                let mut nitrogen_five_index = 0_i32;
                while nitrogen_five_index < number_nitrogen_five
                    && current_success < max_success
                    && !leave_case
                {
                    let nitrogen_five_atom =
                        i32::from(nitrogen_five_atoms[nitrogen_five_index as usize]);
                    let nitrogen_five_atom_index = usize::try_from(nitrogen_five_atom)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let nitrogen_five_plus = pVA
                        .get(nitrogen_five_atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let nitrogen_five_flower = if nitrogen_five_plus > 0 {
                        GetChargeFlowerUpperEdge(heap, pBNS, pVA, nitrogen_five_plus)?
                    } else {
                        NO_VERTEX
                    };
                    let valid_nitrogen_five = nitrogen_five_atom != NO_VERTEX
                        && nitrogen_five_plus > 0
                        && nitrogen_five_flower != NO_VERTEX
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(nitrogen_five_plus)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            == 1
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(nitrogen_five_flower)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            == 0;
                    if !valid_nitrogen_five {
                        nitrogen_five_index = nitrogen_five_index.wrapping_add(1);
                        continue;
                    }
                    let mut nitrogen_three_list_index = nitrogen_three_atoms.num_edges - 1;
                    while nitrogen_three_list_index >= 0
                        && current_success < max_success
                        && !leave_case
                    {
                        let list_index = usize::try_from(nitrogen_three_list_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let nitrogen_three_atom = *heap
                            .slice(nitrogen_three_atoms.pnEdges.as_const())?
                            .get(list_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let nitrogen_three_atom_index = if nitrogen_three_atom != NO_VERTEX {
                            Some(
                                usize::try_from(nitrogen_three_atom)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                        } else {
                            None
                        };
                        let (nitrogen_three_minus, nitrogen_three_plus) =
                            if let Some(atom_index) = nitrogen_three_atom_index {
                                let valence = pVA
                                    .get(atom_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                (
                                    valence.nCMinusGroupEdge.wrapping_sub(1),
                                    valence.nCPlusGroupEdge.wrapping_sub(1),
                                )
                            } else {
                                (NO_VERTEX, NO_VERTEX)
                            };
                        let valid_nitrogen_three = nitrogen_three_atom != NO_VERTEX
                            && nitrogen_three_minus > 0
                            && nitrogen_three_plus > 0
                            && heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(nitrogen_three_minus)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .flow
                                == 0
                            && heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(nitrogen_three_plus)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .flow
                                == 1;
                        if !valid_nitrogen_three {
                            *heap
                                .slice_mut(nitrogen_three_atoms.pnEdges)?
                                .get_mut(list_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = NO_VERTEX;
                            nitrogen_three_list_index -= 1;
                            continue;
                        }
                        let mut x_list_index = x_atoms.num_edges - 1;
                        while x_list_index >= 0 && current_success < max_success && !leave_case {
                            let x_index = usize::try_from(x_list_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let x_atom = *heap
                                .slice(x_atoms.pnEdges.as_const())?
                                .get(x_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            ret = 0;
                            let x_atom_index = if x_atom != NO_VERTEX {
                                Some(
                                    usize::try_from(x_atom)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                            } else {
                                None
                            };
                            let x_plus = if let Some(atom_index) = x_atom_index {
                                pVA.get(atom_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .nCPlusGroupEdge
                                    .wrapping_sub(1)
                            } else {
                                NO_VERTEX
                            };
                            let valid_x = x_atom != NO_VERTEX
                                && x_plus > 0
                                && heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(x_plus)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .flow
                                    == 1;
                            if !valid_x {
                                *heap
                                    .slice_mut(x_atoms.pnEdges)?
                                    .get_mut(x_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)? = NO_VERTEX;
                                x_list_index -= 1;
                                continue;
                            }

                            let flower_index = usize::try_from(nitrogen_five_flower)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let flower_before = heap
                                .slice(pBNS.edge.as_const())?
                                .get(flower_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let first_vertex = i32::from(flower_before.neighbor1);
                            let second_vertex = i32::from(flower_before.neighbor12) ^ first_vertex;
                            let first_vertex_data = heap
                                .slice(pBNS.vert.as_const())?
                                .get(
                                    usize::try_from(first_vertex)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut adjacency_index =
                                i32::from(first_vertex_data.num_adj_edges) - 1;
                            let mut first_neighbor = None;
                            while adjacency_index >= 0 {
                                let adjacent_edge_number = *heap
                                    .slice(first_vertex_data.iedge.as_const())?
                                    .get(
                                        usize::try_from(adjacency_index)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let adjacent_edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(adjacent_edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if adjacent_edge.flow != 0
                                    && adjacent_edge.forbidden == 0
                                    && adjacent_edge_number != nitrogen_five_plus
                                {
                                    first_neighbor = Some((
                                        adjacent_edge_number,
                                        i32::from(adjacent_edge.neighbor12) ^ first_vertex,
                                    ));
                                    break;
                                }
                                adjacency_index -= 1;
                            }
                            let Some((first_neighbor_edge, first_neighbor_vertex)) = first_neighbor
                            else {
                                x_list_index -= 1;
                                continue;
                            };
                            let second_vertex_data = heap
                                .slice(pBNS.vert.as_const())?
                                .get(
                                    usize::try_from(second_vertex)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut adjacency_index =
                                i32::from(second_vertex_data.num_adj_edges) - 1;
                            let mut second_neighbor = None;
                            while adjacency_index >= 0 {
                                let adjacent_edge_number = *heap
                                    .slice(second_vertex_data.iedge.as_const())?
                                    .get(
                                        usize::try_from(adjacency_index)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let adjacent_edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(adjacent_edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if adjacent_edge.flow != 0
                                    && adjacent_edge.forbidden == 0
                                    && adjacent_edge_number != nitrogen_five_plus
                                {
                                    second_neighbor = Some((
                                        adjacent_edge_number,
                                        i32::from(adjacent_edge.neighbor12) ^ second_vertex,
                                    ));
                                    break;
                                }
                                adjacency_index -= 1;
                            }
                            let Some((second_neighbor_edge, second_neighbor_vertex)) =
                                second_neighbor
                            else {
                                x_list_index -= 1;
                                continue;
                            };

                            heap.slice_mut(pBNS.edge)?[flower_index].flow =
                                flower_before.flow.wrapping_add(1);
                            for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(
                                        usize::try_from(neighbor_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.flow = edge.flow.wrapping_sub(1);
                            }
                            {
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [first_neighbor_vertex, second_neighbor_vertex]
                                {
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
                            SetForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &all_charge_edges,
                                forbidden_edge_mask,
                            )?;
                            SetForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &other_nitrogen_flower_edges,
                                forbidden_edge_mask,
                            )?;
                            for edge_number in [nitrogen_three_minus, x_plus] {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(
                                        usize::try_from(edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden =
                                    (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
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
                            if ret < 0 {
                                leave_case = true;
                                break;
                            }
                            if !(ret == 1
                                && ((path_end == first_neighbor_vertex
                                    && path_start == second_neighbor_vertex)
                                    || (path_end == second_neighbor_vertex
                                        && path_start == first_neighbor_vertex))
                                && delta_charge == 2)
                            {
                                ret = 0;
                            }
                            for edge_number in [nitrogen_three_minus, x_plus] {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(
                                        usize::try_from(edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden =
                                    (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                            }
                            {
                                let flower = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(flower_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                flower.flow = flower.flow.wrapping_sub(1);
                            }
                            for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(
                                        usize::try_from(neighbor_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.flow = edge.flow.wrapping_add(1);
                            }
                            {
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [first_neighbor_vertex, second_neighbor_vertex]
                                {
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

                            if ret == 1 {
                                {
                                    let flower = heap
                                        .slice_mut(pBNS.edge)?
                                        .get_mut(flower_index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    flower.flow = flower.flow.wrapping_add(1);
                                }
                                for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                                    let edge = heap
                                        .slice_mut(pBNS.edge)?
                                        .get_mut(
                                            usize::try_from(neighbor_edge)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    edge.flow = edge.flow.wrapping_sub(1);
                                }
                                {
                                    let vertices = heap.slice_mut(pBNS.vert)?;
                                    for vertex_number in
                                        [first_neighbor_vertex, second_neighbor_vertex]
                                    {
                                        let vertex =
                                            vertices
                                                .get_mut(usize::try_from(vertex_number).map_err(
                                                    |_| SourceHeapError::PointerOutOfBounds,
                                                )?)
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                                    }
                                }
                                pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
                                for edge_number in [nitrogen_three_minus, nitrogen_five_plus] {
                                    let edge = heap
                                        .slice_mut(pBNS.edge)?
                                        .get_mut(
                                            usize::try_from(edge_number)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    edge.forbidden =
                                        (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                                }
                                path_start = 0;
                                path_end = 0;
                                path_length = 0;
                                delta_h = 0;
                                delta_charge = 0;
                                number_visited_atoms = 0;
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
                                    && ((path_end == first_neighbor_vertex
                                        && path_start == second_neighbor_vertex)
                                        || (path_end == second_neighbor_vertex
                                            && path_start == first_neighbor_vertex))
                                    && delta_charge == 2
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
                                        n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                        current_success = current_success.wrapping_add(1);
                                    }
                                }
                                if ret <= 0 {
                                    let flower = heap
                                        .slice_mut(pBNS.edge)?
                                        .get_mut(flower_index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    flower.flow = flower.flow.wrapping_sub(1);
                                    for neighbor_edge in [first_neighbor_edge, second_neighbor_edge]
                                    {
                                        let edge =
                                            heap.slice_mut(pBNS.edge)?
                                                .get_mut(usize::try_from(neighbor_edge).map_err(
                                                    |_| SourceHeapError::PointerOutOfBounds,
                                                )?)
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        edge.flow = edge.flow.wrapping_add(1);
                                    }
                                    let vertices = heap.slice_mut(pBNS.vert)?;
                                    for vertex_number in
                                        [first_neighbor_vertex, second_neighbor_vertex]
                                    {
                                        let vertex =
                                            vertices
                                                .get_mut(usize::try_from(vertex_number).map_err(
                                                    |_| SourceHeapError::PointerOutOfBounds,
                                                )?)
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(1);
                                    }
                                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                                }
                            }
                            RemoveForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &all_charge_edges,
                                forbidden_edge_mask,
                            )?;
                            RemoveForbiddenEdgeMask(
                                heap,
                                pBNS,
                                &other_nitrogen_flower_edges,
                                forbidden_edge_mask,
                            )?;
                            if ret > 0 {
                                nitrogen_five_atoms[nitrogen_five_index as usize] =
                                    NO_VERTEX as i16;
                                *heap
                                    .slice_mut(nitrogen_three_atoms.pnEdges)?
                                    .get_mut(list_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)? = NO_VERTEX;
                                *heap
                                    .slice_mut(x_atoms.pnEdges)?
                                    .get_mut(x_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)? = NO_VERTEX;
                            }
                            if ret < 0 {
                                leave_case = true;
                                break;
                            }
                            if ret > 0 {
                                break;
                            }
                            x_list_index -= 1;
                        }
                        if ret > 0 {
                            break;
                        }
                        nitrogen_three_list_index -= 1;
                    }
                    nitrogen_five_index = nitrogen_five_index.wrapping_add(1);
                }
            }

            let _ = AllocEdgeList(heap, &mut x_atoms, EDGE_LIST_FREE)?;
            let _ = AllocEdgeList(heap, &mut nitrogen_three_atoms, EDGE_LIST_FREE)?;
            RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            RemoveForbiddenEdgeMask(
                heap,
                pBNS,
                &other_nitrogen_flower_edges,
                forbidden_edge_mask,
            )?;
            current_edges.num_edges = 0;
            if ret < 0 {
                return Ok(());
            }
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 15, ichirvr3.c:4029-4209.
        if comparison.nNumTgMRevrs > comparison.nNumTgMInChI
            || comparison.nNumRemHRevrs < comparison.nNumRemHInChI
            || comparison.nNumEndpRevrs < comparison.nNumEndpInChI
            || (comparison.nNumTgInChI <= 1 && comparison.nNumTgRevrs > comparison.nNumTgInChI)
        {
            let mut single_bond_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_charged = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_neutral = 0_i32;
            let mut number_double_bond_charged = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let fixed_h = !pStruct.fixed_H.is_null()
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let mobile_h = !mobile_h_input.is_null()
                    && *heap
                        .slice(mobile_h_input.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let edge_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_double_bond_charged < MAX_DIFF_FIXH as i32
                    && atom.charge == 1
                    && atom.num_H == 0
                    && atom.valence < atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 6
                    && !fixed_h
                    && !mobile_h
                    && edge_available
                {
                    double_bond_charged[number_double_bond_charged as usize] = atom_number as i16;
                    number_double_bond_charged = number_double_bond_charged.wrapping_add(1);
                } else if number_single_bond_neutral < MAX_DIFF_FIXH as i32
                    && atom.charge == 0
                    && atom.num_H == 0
                    && atom.valence == atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 5
                    && valence.cPeriodicRowNumber == 1
                    && !fixed_h
                    && !mobile_h
                    && edge_available
                {
                    single_bond_neutral[number_single_bond_neutral as usize] = atom_number as i16;
                    number_single_bond_neutral = number_single_bond_neutral.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            let number_to_try = number_single_bond_neutral.min(number_double_bond_charged);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < number_double_bond_charged
                    && current_success < number_to_try
                {
                    let atom_index =
                        usize::try_from(i32::from(double_bond_charged[candidate_index as usize]))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow != 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    let first_vertex_data = heap
                        .slice(pBNS.vert.as_const())?
                        .get(
                            usize::try_from(first_vertex)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut adjacency_index = i32::from(first_vertex_data.num_adj_edges) - 1;
                    let mut first_neighbor = None;
                    while adjacency_index >= 0 {
                        let adjacent_edge_number = *heap
                            .slice(first_vertex_data.iedge.as_const())?
                            .get(
                                usize::try_from(adjacency_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(adjacent_edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                            first_neighbor = Some((
                                adjacent_edge_number,
                                i32::from(adjacent_edge.neighbor12) ^ first_vertex,
                            ));
                            break;
                        }
                        adjacency_index -= 1;
                    }
                    let Some((first_neighbor_edge, first_neighbor_vertex)) = first_neighbor else {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    };
                    let second_vertex_data = heap
                        .slice(pBNS.vert.as_const())?
                        .get(
                            usize::try_from(second_vertex)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut adjacency_index = i32::from(second_vertex_data.num_adj_edges) - 1;
                    let mut second_neighbor = None;
                    while adjacency_index >= 0 {
                        let adjacent_edge_number = *heap
                            .slice(second_vertex_data.iedge.as_const())?
                            .get(
                                usize::try_from(adjacency_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(adjacent_edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                            second_neighbor = Some((
                                adjacent_edge_number,
                                i32::from(adjacent_edge.neighbor12) ^ second_vertex,
                            ));
                            break;
                        }
                        adjacency_index -= 1;
                    }
                    let Some((second_neighbor_edge, second_neighbor_vertex)) = second_neighbor
                    else {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    };
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow.wrapping_add(1);
                    for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                        let adjacent_edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(
                                usize::try_from(neighbor_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        adjacent_edge.flow = adjacent_edge.flow.wrapping_sub(1);
                    }
                    {
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_neighbor_vertex, second_neighbor_vertex] {
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
                        && ((path_end == first_neighbor_vertex
                            && path_start == second_neighbor_vertex)
                            || (path_end == second_neighbor_vertex
                                && path_start == first_neighbor_vertex))
                        && (delta_charge == 0 || delta_charge == 1)
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.flow = edge.flow.wrapping_sub(1);
                        for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                            let adjacent_edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(
                                    usize::try_from(neighbor_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            adjacent_edge.flow = adjacent_edge.flow.wrapping_add(1);
                        }
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_neighbor_vertex, second_neighbor_vertex] {
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
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 16, ichirvr3.c:4207-4365.
        if comparison.nNumTgDiffMinus != 0 {
            let mut single_bond_nitrogen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_oxygen_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut tgroup_index = 0_i32;
            while tgroup_index < pStruct.ti.num_t_groups
                && tgroup_index < pStruct.One_ti.num_t_groups
            {
                let group_index = usize::try_from(tgroup_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let input_group = heap
                    .slice(pStruct.ti.t_group.as_const())?
                    .get(group_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_group = heap
                    .slice(pStruct.One_ti.t_group.as_const())?
                    .get(group_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if input_group.nNumEndpoints != reversed_group.nNumEndpoints
                    || input_group.num[1] >= reversed_group.num[1]
                {
                    tgroup_index = tgroup_index.wrapping_add(1);
                    continue;
                }

                current_edges.num_edges = 0;
                let mut number_single_bond_nitrogen_minus = 0_i32;
                let mut number_double_bond_oxygen_neutral = 0_i32;
                current_success = 0;
                let first_endpoint = i32::from(reversed_group.nFirstEndpointAtNoPos);
                let mut endpoint_index = 0_i32;
                while endpoint_index < i32::from(reversed_group.nNumEndpoints) {
                    let endpoint_position =
                        usize::try_from(first_endpoint.wrapping_add(endpoint_index))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let canonical_number = i32::from(
                        *heap
                            .slice(pStruct.One_ti.nEndpointAtomNumber.as_const())?
                            .get(endpoint_position)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let input_mobile_h = !mobile_h_input.is_null()
                        && *heap
                            .slice(mobile_h_input.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            != 0;
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    let minus_available = minus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .forbidden
                            == 0;
                    if number_double_bond_oxygen_neutral < MAX_DIFF_FIXH as i32
                        && atom.charge == 0
                        && atom.num_H == 0
                        && atom.valence < atom.chem_bonds_valence
                        && valence.cMetal == 0
                        && valence.cNumValenceElectrons == 6
                        && input_endpoint != 0
                        && !input_mobile_h
                        && minus_available
                    {
                        double_bond_oxygen_neutral[number_double_bond_oxygen_neutral as usize] =
                            atom_number as i16;
                        number_double_bond_oxygen_neutral =
                            number_double_bond_oxygen_neutral.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    } else {
                        let input_fixed_h = !pStruct.fixed_H.is_null()
                            && *heap
                                .slice(pStruct.fixed_H.as_const())?
                                .get(canonical_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                != 0;
                        if number_single_bond_nitrogen_minus < MAX_DIFF_FIXH as i32
                            && atom.charge == -1
                            && atom.num_H != 0
                            && atom.valence == atom.chem_bonds_valence
                            && valence.cMetal == 0
                            && valence.cNumValenceElectrons == 5
                            && valence.cPeriodicRowNumber == 1
                            && input_endpoint != 0
                            && input_fixed_h
                            && !input_mobile_h
                            && minus_available
                        {
                            single_bond_nitrogen_minus
                                [number_single_bond_nitrogen_minus as usize] = atom_number as i16;
                            number_single_bond_nitrogen_minus =
                                number_single_bond_nitrogen_minus.wrapping_add(1);
                        }
                    }
                    endpoint_index = endpoint_index.wrapping_add(1);
                }

                let number_to_try =
                    number_single_bond_nitrogen_minus.min(number_double_bond_oxygen_neutral);
                if number_to_try != 0 {
                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                    let mut candidate_index = 0_i32;
                    while candidate_index < number_single_bond_nitrogen_minus
                        && current_success < number_to_try
                    {
                        let atom_index = usize::try_from(i32::from(
                            single_bond_nitrogen_minus[candidate_index as usize],
                        ))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let minus_edge = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCMinusGroupEdge
                            .wrapping_sub(1);
                        let edge_index = usize::try_from(minus_edge)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let edge_before = heap
                            .slice(pBNS.edge.as_const())?
                            .get(edge_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if edge_before.flow == 0 {
                            candidate_index = candidate_index.wrapping_add(1);
                            continue;
                        }
                        let first_vertex = i32::from(edge_before.neighbor1);
                        let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
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
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                current_success = current_success.wrapping_add(1);
                            }
                        } else {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(edge_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                            edge.flow = edge.flow.wrapping_add(1);
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
                        candidate_index = candidate_index.wrapping_add(1);
                    }
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                }
                current_edges.num_edges = 0;
                if current_success != 0 {
                    tot_success = tot_success.wrapping_add(current_success);
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
                    ret = FillOutCMP2FHINCHI(
                        heap,
                        pStruct,
                        &at2_snapshot,
                        pVA,
                        pInChI,
                        &mut comparison,
                    )?;
                    if ret != 0 || comparison.bHasDifference == 0 {
                        return Ok(());
                    }
                }
                tgroup_index = tgroup_index.wrapping_add(1);
            }
        }

        // Complete translation of source case 17, ichirvr3.c:4367-4538.
        if comparison.nNumRemHInChI < comparison.nNumRemHRevrs {
            let mut single_bond_neutral = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_charged = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_neutral = 0_i32;
            let mut number_double_bond_charged = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_input = if !pInChI[1].is_null() {
                let inchi = heap
                    .slice(pInChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pInChI[0].is_null() {
                    heap.slice(pInChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pInChI[0].is_null() {
                heap.slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut current_success = 0_i32;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                let plus_available = plus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_double_bond_charged < MAX_DIFF_FIXH as i32
                    && atom.charge == 1
                    && atom.num_H != 0
                    && atom.valence < atom.chem_bonds_valence
                    && valence.cMetal == 0
                    && (valence.cNumValenceElectrons == 6
                        || (valence.cNumValenceElectrons == 5 && valence.cPeriodicRowNumber == 1))
                    && plus_available
                {
                    double_bond_charged[number_double_bond_charged as usize] = atom_number as i16;
                    number_double_bond_charged = number_double_bond_charged.wrapping_add(1);
                } else {
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let input_fixed_h = !pStruct.fixed_H.is_null()
                        && *heap
                            .slice(pStruct.fixed_H.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            != 0;
                    let input_mobile_h = !mobile_h_input.is_null()
                        && *heap
                            .slice(mobile_h_input.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            != 0;
                    if number_single_bond_neutral < MAX_DIFF_FIXH as i32
                        && atom.charge == 0
                        && atom.num_H == 0
                        && atom.valence == atom.chem_bonds_valence
                        && valence.cMetal == 0
                        && (valence.cNumValenceElectrons == 6
                            || valence.cNumValenceElectrons == 7
                            || (valence.cNumValenceElectrons == 5
                                && valence.cPeriodicRowNumber > 1))
                        && input_endpoint == 0
                        && !input_fixed_h
                        && !input_mobile_h
                        && plus_available
                    {
                        single_bond_neutral[number_single_bond_neutral as usize] =
                            atom_number as i16;
                        number_single_bond_neutral = number_single_bond_neutral.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }
            let mut number_to_try = number_single_bond_neutral.min(number_double_bond_charged);
            if number_to_try != 0 {
                number_to_try = number_to_try.min(
                    i32::from(comparison.nNumRemHRevrs)
                        .wrapping_sub(i32::from(comparison.nNumRemHInChI)),
                );
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < number_double_bond_charged
                    && current_success < number_to_try
                {
                    let atom_index =
                        usize::try_from(i32::from(double_bond_charged[candidate_index as usize]))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_number = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow != 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    let first_vertex_data = heap
                        .slice(pBNS.vert.as_const())?
                        .get(
                            usize::try_from(first_vertex)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut adjacency_index = i32::from(first_vertex_data.num_adj_edges) - 1;
                    let mut first_neighbor = None;
                    while adjacency_index >= 0 {
                        let adjacent_edge_number = *heap
                            .slice(first_vertex_data.iedge.as_const())?
                            .get(
                                usize::try_from(adjacency_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(adjacent_edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                            first_neighbor = Some((
                                adjacent_edge_number,
                                i32::from(adjacent_edge.neighbor12) ^ first_vertex,
                            ));
                            break;
                        }
                        adjacency_index -= 1;
                    }
                    let Some((first_neighbor_edge, first_neighbor_vertex)) = first_neighbor else {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    };
                    let second_vertex_data = heap
                        .slice(pBNS.vert.as_const())?
                        .get(
                            usize::try_from(second_vertex)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut adjacency_index = i32::from(second_vertex_data.num_adj_edges) - 1;
                    let mut second_neighbor = None;
                    while adjacency_index >= 0 {
                        let adjacent_edge_number = *heap
                            .slice(second_vertex_data.iedge.as_const())?
                            .get(
                                usize::try_from(adjacency_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let adjacent_edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(adjacent_edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if adjacent_edge.flow != 0 && adjacent_edge.forbidden == 0 {
                            second_neighbor = Some((
                                adjacent_edge_number,
                                i32::from(adjacent_edge.neighbor12) ^ second_vertex,
                            ));
                            break;
                        }
                        adjacency_index -= 1;
                    }
                    let Some((second_neighbor_edge, second_neighbor_vertex)) = second_neighbor
                    else {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    };
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow.wrapping_add(1);
                    for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                        let adjacent_edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(
                                usize::try_from(neighbor_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        adjacent_edge.flow = adjacent_edge.flow.wrapping_sub(1);
                    }
                    {
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_neighbor_vertex, second_neighbor_vertex] {
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
                        && ((path_end == first_neighbor_vertex
                            && path_start == second_neighbor_vertex)
                            || (path_end == second_neighbor_vertex
                                && path_start == first_neighbor_vertex))
                        && (delta_charge == 0 || delta_charge == 1)
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.flow = edge.flow.wrapping_sub(1);
                        for neighbor_edge in [first_neighbor_edge, second_neighbor_edge] {
                            let adjacent_edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(
                                    usize::try_from(neighbor_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            adjacent_edge.flow = adjacent_edge.flow.wrapping_add(1);
                        }
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_neighbor_vertex, second_neighbor_vertex] {
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
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 18, ichirvr3.c:4540-4683.
        if comparison.nNumTgInChI != 0
            && !pStruct.endpoint.is_null()
            && comparison.nNumTgMInChI > comparison.nNumTgMRevrs
            && comparison.nNumEndpInChI > comparison.nNumEndpRevrs
        {
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;
            ret = 0;
            let mut tgroup_index = 0_i32;
            while tgroup_index < pStruct.ti.num_t_groups
                && tgroup_index < pStruct.One_ti.num_t_groups
            {
                let group_index = usize::try_from(tgroup_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let input_group = heap
                    .slice(pStruct.ti.t_group.as_const())?
                    .get(group_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_group = heap
                    .slice(pStruct.One_ti.t_group.as_const())?
                    .get(group_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if input_group.nNumEndpoints <= reversed_group.nNumEndpoints
                    || input_group.num[1] <= reversed_group.num[1]
                {
                    tgroup_index = tgroup_index.wrapping_add(1);
                    continue;
                }
                current_edges.num_edges = 0;
                current_success = 0;
                let first_endpoint = i32::from(input_group.nFirstEndpointAtNoPos);
                let mut endpoint_index = 0_i32;
                while endpoint_index < i32::from(input_group.nNumEndpoints) {
                    let endpoint_position =
                        usize::try_from(first_endpoint.wrapping_add(endpoint_index))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let canonical_number = i32::from(
                        *heap
                            .slice(pStruct.ti.nEndpointAtomNumber.as_const())?
                            .get(endpoint_position)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let reversed_endpoint = !mobile_h_reversed.is_null()
                        && heap
                            .slice(mobile_h_reversed.as_const())?
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .endpoint
                            != 0;
                    let canonical_valence = pVA
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let atom_valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let plus_edge = atom_valence.nCPlusGroupEdge.wrapping_sub(1);
                    let plus_rejects = plus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            == 0;
                    let minus_edge = atom_valence.nCMinusGroupEdge.wrapping_sub(1);
                    let minus_available = minus_edge >= 0 && {
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden == 0 && edge.flow == 0
                    };
                    if input_endpoint != 0
                        && !mobile_h_reversed.is_null()
                        && !reversed_endpoint
                        && canonical_valence.cNumValenceElectrons == 5
                        && canonical_valence.cPeriodicRowNumber == 1
                        && atom.valence == 2
                        && atom.num_H == 0
                        && atom.radical == 0
                        && !plus_rejects
                        && minus_available
                    {
                        ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    endpoint_index = endpoint_index.wrapping_add(1);
                }
                tgroup_index = tgroup_index.wrapping_add(1);
            }

            let max_success = current_edges.num_edges;
            if max_success != 0 {
                let mut canonical_number = 0_i32;
                while canonical_number < pStruct.num_atoms && current_success < max_success {
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let canonical_valence = pVA
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let atom_valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let minus_edge = atom_valence.nCMinusGroupEdge.wrapping_sub(1);
                    let minus_rejects = minus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            == 0;
                    let plus_edge = atom_valence.nCPlusGroupEdge.wrapping_sub(1);
                    let plus_available = plus_edge >= 0 && {
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden == 0 && edge.flow == 1
                    };
                    if input_endpoint != 0
                        || canonical_valence.cNumValenceElectrons == 0
                        || canonical_valence.cNumValenceElectrons == 4
                        || atom.num_H != 0
                        || atom.radical != 0
                        || minus_rejects
                        || !plus_available
                    {
                        canonical_number = canonical_number.wrapping_add(1);
                        continue;
                    }

                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    SetForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &other_nitrogen_flower_edges,
                        forbidden_edge_mask,
                    )?;
                    RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                    let edge_index = usize::try_from(plus_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        canonical_number = canonical_number.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow.wrapping_sub(1);
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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[edge_index].flow =
                            edge_before.flow.wrapping_add(0);
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
                    canonical_number = canonical_number.wrapping_add(1);
                }
            }

            RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            RemoveForbiddenEdgeMask(
                heap,
                pBNS,
                &other_nitrogen_flower_edges,
                forbidden_edge_mask,
            )?;
            current_edges.num_edges = 0;
            if ret < 0 {
                return Ok(());
            }
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 19, ichirvr3.c:4685-4837.
        if comparison.len_c2at >= 1 {
            let mut current_success = 0_i32;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let mut is_candidate = false;
                if difference.nValElectr == 6 {
                    let valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let oxygen_plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                    let oxygen_plus_available = oxygen_plus_edge >= 0 && {
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(oxygen_plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden == 0 && edge.flow != 0
                    };
                    if oxygen_plus_available
                        && difference.nFixHInChI == 1
                        && difference.nMobHInChI == 0
                        && difference.nFixHRevrs == 0
                        && difference.nMobHRevrs == 1
                        && difference.nAtChargeRevrs == 0
                        && atom.num_H != 0
                        && atom.valence == 1
                        && atom.valence == atom.chem_bonds_valence
                    {
                        let metal_index = usize::from(atom.neighbor[0]);
                        let metal_valence = pVA
                            .get(metal_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if metal_valence.cMetal != 0 {
                            let metal_plus_edge = metal_valence.nCPlusGroupEdge.wrapping_sub(1);
                            let metal_plus_available = metal_plus_edge >= 0
                                && heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(metal_plus_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .forbidden
                                    == 0;
                            if metal_plus_available {
                                let metal_minus_edge =
                                    metal_valence.nCMinusGroupEdge.wrapping_sub(1);
                                let metal_minus_available = metal_minus_edge >= 0
                                    && heap
                                        .slice(pBNS.edge.as_const())?
                                        .get(
                                            usize::try_from(metal_minus_edge)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                        .forbidden
                                        == 0;
                                if metal_minus_available {
                                    let atom_vertex = heap
                                        .slice(pBNS.vert.as_const())?
                                        .get(atom_index)
                                        .cloned()
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    let metal_bond_edge = *heap
                                        .slice(atom_vertex.iedge.as_const())?
                                        .first()
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    is_candidate = heap
                                        .slice(pBNS.edge.as_const())?
                                        .get(
                                            usize::try_from(metal_bond_edge)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                        .forbidden
                                        == 0;
                                }
                            }
                        }
                    }
                }
                if is_candidate {
                    ret = AddToEdgeList(heap, &mut current_edges, atom_number, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            if current_edges.num_edges != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                SetForbiddenEdgeMask(heap, pBNS, &nitrogen_flower_edges, forbidden_edge_mask)?;
                SetForbiddenEdgeMask(heap, pBNS, &all_bond_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < current_edges.num_edges {
                    let atom_number = *heap
                        .slice(current_edges.pnEdges.as_const())?
                        .get(
                            usize::try_from(candidate_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let metal_number = i32::from(atom.neighbor[0]);
                    let metal_index = usize::try_from(metal_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let oxygen_plus_edge = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let metal_plus_edge = pVA
                        .get(metal_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCPlusGroupEdge
                        .wrapping_sub(1);
                    let metal_minus_edge = pVA
                        .get(metal_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let atom_vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let metal_bond_edge = *heap
                        .slice(atom_vertex.iedge.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    for edge_number in [metal_plus_edge, metal_minus_edge, metal_bond_edge] {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(
                                usize::try_from(edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                    }
                    let metal_plus = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(metal_plus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let metal_minus = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(metal_minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let charge_on_metal = (i32::from(metal_plus.cap)
                        .wrapping_sub(i32::from(metal_plus.flow)))
                    .wrapping_sub(i32::from(metal_minus.flow));
                    let expected_delta_charge = if charge_on_metal == 1 {
                        -1
                    } else if charge_on_metal == 0 {
                        1
                    } else {
                        0
                    };
                    let oxygen_plus_index = usize::try_from(oxygen_plus_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(oxygen_plus_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[oxygen_plus_index].flow =
                        edge_before.flow.wrapping_sub(1);
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
                        && delta_charge == expected_delta_charge
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[oxygen_plus_index].flow =
                            edge_before.flow.wrapping_add(0);
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
                    for edge_number in [metal_plus_edge, metal_minus_edge, metal_bond_edge] {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(
                                usize::try_from(edge_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                    }
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &nitrogen_flower_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &all_bond_edges, forbidden_edge_mask)?;
                current_edges.num_edges = 0;
                if current_success != 0 {
                    tot_success = tot_success.wrapping_add(current_success);
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
                    ret = FillOutCMP2FHINCHI(
                        heap,
                        pStruct,
                        &at2_snapshot,
                        pVA,
                        pInChI,
                        &mut comparison,
                    )?;
                    if ret != 0 || comparison.bHasDifference == 0 {
                        return Ok(());
                    }
                }
            }
        }

        // Complete translation of source case 20, ichirvr3.c:4839-5009.
        if comparison.len_c2at > 1 && comparison.nNumTgRevrs != 0 && comparison.nNumTgInChI != 0 {
            let mut single_bond_oxygen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_nitrogen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_oxygen_minus = 0_i32;
            let mut number_double_bond_nitrogen = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mut current_success = 0_i32;
            current_edges.num_edges = 0;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                let minus_state = if minus_edge >= 0 {
                    Some(
                        heap.slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    )
                } else {
                    None
                };
                let double_bond_candidate = number_double_bond_nitrogen < MAX_DIFF_FIXH as i32
                    && difference.endptInChI != 0
                    && minus_state
                        .as_ref()
                        .is_some_and(|edge| edge.forbidden == 0 && edge.flow == 0)
                    && difference.nFixHInChI == 0
                    && difference.nMobHInChI == 0
                    && (comparison.nNumTgInChI == 1 || difference.nValElectr == 6)
                    && difference.endptRevrs == 0
                    && difference.nFixHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
                    && i32::from(atom.valence).wrapping_add(1)
                        == i32::from(atom.chem_bonds_valence);
                if double_bond_candidate {
                    double_bond_nitrogen[number_double_bond_nitrogen as usize] = atom_number as i16;
                    number_double_bond_nitrogen = number_double_bond_nitrogen.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                } else if number_single_bond_oxygen_minus < MAX_DIFF_FIXH as i32
                    && difference.endptInChI == 0
                    && minus_state
                        .as_ref()
                        .is_some_and(|edge| edge.forbidden == 0 && edge.flow == 1)
                    && difference.nFixHInChI == 0
                    && difference.nMobHInChI == 0
                    && difference.nValElectr == 6
                    && difference.endptRevrs != 0
                    && difference.nFixHRevrs == 0
                    && difference.nMobHRevrs == 0
                    && difference.nAtChargeRevrs == -1
                    && atom.num_H == 0
                    && atom.valence == 1
                    && atom.chem_bonds_valence == 1
                {
                    single_bond_oxygen_minus[number_single_bond_oxygen_minus as usize] =
                        atom_number as i16;
                    number_single_bond_oxygen_minus =
                        number_single_bond_oxygen_minus.wrapping_add(1);
                }
                difference_index = difference_index.wrapping_add(1);
            }

            if number_double_bond_nitrogen == 0 {
                let mut canonical_number = 0_i32;
                while canonical_number < pStruct.num_atoms {
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if input_endpoint == 0 {
                        canonical_number = canonical_number.wrapping_add(1);
                        continue;
                    }
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let input_fixed_h = !pStruct.fixed_H.is_null()
                        && *heap
                            .slice(pStruct.fixed_H.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            != 0;
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    let minus_available = minus_edge >= 0 && {
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden == 0 && edge.flow == 0
                    };
                    if number_double_bond_nitrogen < MAX_DIFF_FIXH as i32
                        && atom.charge == 0
                        && atom.num_H == 0
                        && i32::from(atom.valence).wrapping_add(1)
                            == i32::from(atom.chem_bonds_valence)
                        && valence.cMetal == 0
                        && !input_fixed_h
                        && minus_available
                    {
                        double_bond_nitrogen[number_double_bond_nitrogen as usize] =
                            atom_number as i16;
                        number_double_bond_nitrogen = number_double_bond_nitrogen.wrapping_add(1);
                        ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    canonical_number = canonical_number.wrapping_add(1);
                }
            }

            let number_to_try = number_single_bond_oxygen_minus.min(number_double_bond_nitrogen);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < number_single_bond_oxygen_minus
                    && current_success < number_to_try
                {
                    let atom_index = usize::try_from(i32::from(
                        single_bond_oxygen_minus[candidate_index as usize],
                    ))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let minus_edge = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(minus_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow.wrapping_sub(1);
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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow;
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
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 21, ichirvr3.c:5011-5217.
        if comparison.len_c2at != 0
            && comparison.nNumTgRevrs != 0
            && comparison.nNumTgHInChI != 0
            && !pStruct.endpoint.is_null()
        {
            let mut single_bond_oxygen_minus = [0_i16; MAX_DIFF_FIXH as usize];
            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut central_atoms = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_single_bond_oxygen_minus = 0_i32;
            let mut number_double_bond_oxygen = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                let minus_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                    && difference.nValElectr == 6
                    && (difference.endptInChI != 0 || difference.nMobHInChI != 0)
                    && minus_available
                    && difference.nFixHInChI == 0
                    && difference.endptRevrs == 0
                    && difference.nMobHRevrs == 0
                    && difference.nFixHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
                    && atom.valence == 1
                    && atom.chem_bonds_valence == 2
                {
                    double_bond_oxygen[number_double_bond_oxygen as usize] = atom_number as i16;
                    number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let mut canonical_number = 0_i32;
            while number_double_bond_oxygen != 0 && canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let input_endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if input_endpoint == 0 {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reversed_endpoint = !mobile_h_reversed.is_null()
                    && heap
                        .slice(mobile_h_reversed.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint
                        != 0;
                let input_fixed_h = !pStruct.fixed_H.is_null()
                    && *heap
                        .slice(pStruct.fixed_H.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let minus_available = minus_edge >= 0 && {
                    let edge = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.forbidden == 0 && edge.flow != 0
                };
                if number_single_bond_oxygen_minus < MAX_DIFF_FIXH as i32
                    && atom.charge == -1
                    && atom.num_H == 0
                    && atom.valence == 1
                    && atom.chem_bonds_valence != 0
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 6
                    && reversed_endpoint
                    && !input_fixed_h
                    && minus_available
                {
                    let canonical_atom = heap
                        .slice(at2.as_const())?
                        .get(canonical_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let central_number = i32::from(canonical_atom.neighbor[0]);
                    let mut duplicate_index = 0_i32;
                    while duplicate_index < number_single_bond_oxygen_minus
                        && i32::from(central_atoms[duplicate_index as usize]) != central_number
                    {
                        duplicate_index = duplicate_index.wrapping_add(1);
                    }
                    if duplicate_index < number_single_bond_oxygen_minus {
                        canonical_number = canonical_number.wrapping_add(1);
                        continue;
                    }
                    let central_index = usize::try_from(central_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let central_at = heap
                        .slice(at.as_const())?
                        .get(central_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let central_at2 = heap
                        .slice(at2.as_const())?
                        .get(central_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut number_tautomer_single_bonds = 0_i32;
                    let mut number_tautomer_double_bonds = 0_i32;
                    let mut number_other_double_bonds = 0_i32;
                    let mut number_other_single_bonds = 0_i32;
                    let mut number_other_bonds = 0_i32;
                    let mut number_negative_endpoints = 0_i32;
                    let mut neighbor_index = 0_i32;
                    while neighbor_index < i32::from(central_at.valence) {
                        let neighbor_position = usize::try_from(neighbor_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let bond_type = i32::from(central_at2.bond_type[neighbor_position]);
                        let neighbor_number = i32::from(central_at2.neighbor[neighbor_position]);
                        if neighbor_number == canonical_number {
                            neighbor_index = neighbor_index.wrapping_add(1);
                            continue;
                        }
                        let neighbor_index_usize = usize::try_from(neighbor_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let neighbor_input_endpoint = *heap
                            .slice(pStruct.endpoint.as_const())?
                            .get(neighbor_index_usize)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if neighbor_input_endpoint == input_endpoint {
                            number_tautomer_single_bonds = number_tautomer_single_bonds
                                .wrapping_add(i32::from(bond_type == BOND_TYPE_SINGLE as i32));
                            number_tautomer_double_bonds = number_tautomer_double_bonds
                                .wrapping_add(i32::from(bond_type == BOND_TYPE_DOUBLE as i32));
                        } else if bond_type == BOND_TYPE_DOUBLE as i32 {
                            number_other_double_bonds = number_other_double_bonds.wrapping_add(1);
                        } else if bond_type == BOND_TYPE_SINGLE as i32 {
                            number_other_single_bonds = number_other_single_bonds.wrapping_add(1);
                        } else {
                            number_other_bonds = number_other_bonds.wrapping_add(1);
                        }
                        let neighbor_atom = heap
                            .slice(at2.as_const())?
                            .get(neighbor_index_usize)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if neighbor_atom.endpoint == canonical_atom.endpoint
                            && neighbor_atom.valence == 1
                            && neighbor_atom.charge == -1
                            && pVA
                                .get(neighbor_index_usize)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .cNumValenceElectrons
                                == 6
                        {
                            number_negative_endpoints = number_negative_endpoints.wrapping_add(1);
                        }
                        neighbor_index = neighbor_index.wrapping_add(1);
                    }
                    let _ = (
                        number_other_single_bonds,
                        number_other_bonds,
                        number_negative_endpoints,
                    );
                    if number_tautomer_single_bonds == 0
                        || !(number_other_double_bonds != 0 && number_tautomer_double_bonds != 0)
                    {
                        canonical_number = canonical_number.wrapping_add(1);
                        continue;
                    }
                    single_bond_oxygen_minus[number_single_bond_oxygen_minus as usize] =
                        atom_number as i16;
                    central_atoms[number_single_bond_oxygen_minus as usize] = central_number as i16;
                    number_single_bond_oxygen_minus =
                        number_single_bond_oxygen_minus.wrapping_add(1);
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            let number_to_try = number_single_bond_oxygen_minus.min(number_double_bond_oxygen);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < number_single_bond_oxygen_minus
                    && current_success < number_to_try
                {
                    let atom_index = usize::try_from(i32::from(
                        single_bond_oxygen_minus[candidate_index as usize],
                    ))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let minus_edge = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let edge_index = usize::try_from(minus_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
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
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow;
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
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 21a, ichirvr3.c:5219-5460.
        if comparison.len_c2at != 0
            && comparison.nNumTgRevrs != 0
            && comparison.nNumEndpInChI < comparison.nNumEndpRevrs
        {
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let mut other_sulfur_oxygen = EDGE_LIST::default();
            let mut central_sulfur = EDGE_LIST::default();
            let mut sulfur_oxygen_minus = EDGE_LIST::default();
            let mut minus_acceptor = EDGE_LIST::default();
            current_edges.num_edges = 0;
            let _ = AllocEdgeList(heap, &mut other_sulfur_oxygen, EDGE_LIST_CLEAR)?;
            let _ = AllocEdgeList(heap, &mut central_sulfur, EDGE_LIST_CLEAR)?;
            let _ = AllocEdgeList(heap, &mut sulfur_oxygen_minus, EDGE_LIST_CLEAR)?;
            let _ = AllocEdgeList(heap, &mut minus_acceptor, EDGE_LIST_CLEAR)?;
            let mut current_success = 0_i32;
            let mut leave_case = mobile_h_reversed.is_null();

            let mut difference_index = 0_i32;
            while !leave_case && difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let minus_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if difference.endptInChI == 0
                    && minus_available
                    && difference.nFixHInChI == 0
                    && (difference.endptRevrs != 0 || difference.nMobHRevrs != 0)
                    && difference.nFixHRevrs == 0
                    && atom.num_H == 0
                {
                    let edge_flow = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .flow;
                    let central_number = i32::from(atom.neighbor[0]);
                    let central_index = usize::try_from(central_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let central_valence = pVA
                        .get(central_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let central_atom = heap
                        .slice(at2.as_const())?
                        .get(central_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if valence.cNumValenceElectrons == 6
                        && atom.charge == -1
                        && edge_flow != 0
                        && atom.valence == 1
                        && atom.chem_bonds_valence == 1
                        && central_valence.cNumValenceElectrons == 6
                        && central_valence.cPeriodicRowNumber > 1
                        && central_atom.valence >= 4
                    {
                        if FindInEdgeList(heap, &central_sulfur, central_number)? >= 0 {
                            difference_index = difference_index.wrapping_add(1);
                            continue;
                        }
                        let central_original = heap
                            .slice(at.as_const())?
                            .get(central_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let difference_atom = heap
                            .slice(at2.as_const())?
                            .get(
                                usize::try_from(difference_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let mut number_tautomer_single_bonds = 0_i32;
                        let mut number_tautomer_double_bonds = 0_i32;
                        let mut number_other_double_bonds = 0_i32;
                        let mut number_other_single_bonds = 0_i32;
                        let mut number_other_bonds = 0_i32;
                        let mut number_negative_endpoints = 0_i32;
                        let mut number_endpoint_oxygen = 0_i32;
                        let mut neighbor_index = 0_i32;
                        while neighbor_index < i32::from(central_original.valence) {
                            let neighbor_position = usize::try_from(neighbor_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let bond_type = i32::from(central_atom.bond_type[neighbor_position]);
                            let neighbor_number =
                                i32::from(central_atom.neighbor[neighbor_position]);
                            if neighbor_number == atom_number {
                                neighbor_index = neighbor_index.wrapping_add(1);
                                continue;
                            }
                            let neighbor_atom_index = usize::try_from(neighbor_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let neighbor_mobile_endpoint = heap
                                .slice(mobile_h_reversed.as_const())?
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .endpoint;
                            let neighbor_atom = heap
                                .slice(at2.as_const())?
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let neighbor_valence = pVA
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if difference.endptRevrs == neighbor_mobile_endpoint
                                && neighbor_atom.endpoint == 0
                            {
                                number_tautomer_single_bonds = number_tautomer_single_bonds
                                    .wrapping_add(i32::from(bond_type == BOND_TYPE_SINGLE as i32));
                                number_tautomer_double_bonds = number_tautomer_double_bonds
                                    .wrapping_add(i32::from(bond_type == BOND_TYPE_DOUBLE as i32));
                                number_endpoint_oxygen =
                                    number_endpoint_oxygen.wrapping_add(i32::from(
                                        neighbor_valence.cNumValenceElectrons == 6
                                            && neighbor_atom.valence == 1,
                                    ));
                            } else if bond_type == BOND_TYPE_DOUBLE as i32 {
                                number_other_double_bonds =
                                    number_other_double_bonds.wrapping_add(1);
                            } else if bond_type == BOND_TYPE_SINGLE as i32 {
                                number_other_single_bonds =
                                    number_other_single_bonds.wrapping_add(1);
                            } else {
                                number_other_bonds = number_other_bonds.wrapping_add(1);
                            }
                            if neighbor_atom.endpoint == difference_atom.endpoint
                                && neighbor_atom.valence == 1
                                && neighbor_atom.charge == -1
                                && neighbor_valence.cNumValenceElectrons == 6
                            {
                                number_negative_endpoints =
                                    number_negative_endpoints.wrapping_add(1);
                            }
                            neighbor_index = neighbor_index.wrapping_add(1);
                        }
                        let _ = (
                            number_other_single_bonds,
                            number_other_bonds,
                            number_negative_endpoints,
                        );
                        if number_endpoint_oxygen == 0
                            || number_tautomer_single_bonds
                                .wrapping_add(number_tautomer_double_bonds)
                                .wrapping_add(number_other_double_bonds)
                                <= number_endpoint_oxygen
                        {
                            difference_index = difference_index.wrapping_add(1);
                            continue;
                        }

                        neighbor_index = 0;
                        while neighbor_index < i32::from(central_original.valence) && !leave_case {
                            let neighbor_position = usize::try_from(neighbor_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let bond_type = i32::from(central_atom.bond_type[neighbor_position]);
                            let neighbor_number =
                                i32::from(central_atom.neighbor[neighbor_position]);
                            let neighbor_atom_index = usize::try_from(neighbor_number)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let neighbor_mobile_endpoint = heap
                                .slice(mobile_h_reversed.as_const())?
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .endpoint;
                            let neighbor_atom = heap
                                .slice(at2.as_const())?
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let neighbor_valence = pVA
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let neighbor_minus_edge =
                                neighbor_valence.nCMinusGroupEdge.wrapping_sub(1);
                            let neighbor_minus_available = neighbor_minus_edge >= 0
                                && heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(neighbor_minus_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .forbidden
                                    == 0;
                            if difference.endptRevrs == neighbor_mobile_endpoint
                                && neighbor_atom.endpoint == 0
                                && neighbor_valence.cNumValenceElectrons == 6
                                && neighbor_atom.valence == 1
                                && neighbor_minus_available
                            {
                                let neighbor_minus_flow = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(neighbor_minus_edge)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .flow;
                                if bond_type == BOND_TYPE_DOUBLE as i32
                                    && neighbor_atom.charge == 0
                                    && neighbor_minus_flow == 0
                                {
                                    ret = AddToEdgeList(
                                        heap,
                                        &mut other_sulfur_oxygen,
                                        neighbor_minus_edge,
                                        INC_ADD_EDGE,
                                    )?;
                                    if ret != 0 {
                                        leave_case = true;
                                    }
                                } else if bond_type == BOND_TYPE_SINGLE as i32
                                    && neighbor_atom.charge == -1
                                    && neighbor_minus_flow != 0
                                {
                                    ret = AddToEdgeList(
                                        heap,
                                        &mut sulfur_oxygen_minus,
                                        neighbor_minus_edge,
                                        INC_ADD_EDGE,
                                    )?;
                                    if ret != 0 {
                                        leave_case = true;
                                    }
                                }
                            }
                            neighbor_index = neighbor_index.wrapping_add(1);
                        }
                        if !leave_case {
                            ret = AddToEdgeList(
                                heap,
                                &mut central_sulfur,
                                central_number,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                    } else if atom.charge == 0
                        && edge_flow == 0
                        && i32::from(atom.valence).wrapping_add(1)
                            == i32::from(atom.chem_bonds_valence)
                    {
                        ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            if !leave_case {
                let mut index = 0_i32;
                while index < other_sulfur_oxygen.num_edges {
                    let edge_number = *heap
                        .slice(other_sulfur_oxygen.pnEdges.as_const())?
                        .get(
                            usize::try_from(index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let _ = RemoveFromEdgeListByValue(heap, &mut current_edges, edge_number)?;
                    index = index.wrapping_add(1);
                }
                let number_to_try = sulfur_oxygen_minus.num_edges.min(current_edges.num_edges);
                if number_to_try != 0 {
                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                    let mut candidate_index = 0_i32;
                    while candidate_index < sulfur_oxygen_minus.num_edges
                        && current_success < number_to_try
                    {
                        let edge_number = *heap
                            .slice(sulfur_oxygen_minus.pnEdges.as_const())?
                            .get(
                                usize::try_from(candidate_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let edge_index = usize::try_from(edge_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let edge_before = heap
                            .slice(pBNS.edge.as_const())?
                            .get(edge_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if edge_before.flow == 0 {
                            candidate_index = candidate_index.wrapping_add(1);
                            continue;
                        }
                        let first_vertex = i32::from(edge_before.neighbor1);
                        let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                        heap.slice_mut(pBNS.edge)?[edge_index].flow =
                            edge_before.flow.wrapping_sub(1);
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
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                current_success = current_success.wrapping_add(1);
                            }
                        } else {
                            heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow;
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
                        candidate_index = candidate_index.wrapping_add(1);
                    }
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                }
            }

            current_edges.num_edges = 0;
            let _ = AllocEdgeList(heap, &mut other_sulfur_oxygen, EDGE_LIST_FREE)?;
            let _ = AllocEdgeList(heap, &mut central_sulfur, EDGE_LIST_FREE)?;
            let _ = AllocEdgeList(heap, &mut sulfur_oxygen_minus, EDGE_LIST_FREE)?;
            let _ = AllocEdgeList(heap, &mut minus_acceptor, EDGE_LIST_FREE)?;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 22, ichirvr3.c:5462-5604.
        if comparison.len_c2at != 0 {
            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_oxygen = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;
            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nCMinusGroupEdge
                    .wrapping_sub(1);
                let minus_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                    && difference.nValElectr == 6
                    && (difference.endptInChI != 0 || difference.nMobHInChI != 0)
                    && minus_available
                    && difference.nFixHInChI == 0
                    && difference.endptRevrs == 0
                    && difference.nMobHRevrs == 0
                    && difference.nFixHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
                    && atom.valence == 1
                    && atom.chem_bonds_valence == 2
                {
                    double_bond_oxygen[number_double_bond_oxygen as usize] = atom_number as i16;
                    number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                }
                difference_index = difference_index.wrapping_add(1);
            }

            let mut canonical_number = 0_i32;
            while number_double_bond_oxygen != 0 && canonical_number < pStruct.num_atoms {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let minus_available = minus_edge >= 0 && {
                    let edge = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    edge.forbidden == 0 && edge.flow != 0
                };
                let reversed_endpoint = !mobile_h_reversed.is_null()
                    && heap
                        .slice(mobile_h_reversed.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint
                        != 0;
                let mut chain_matches = false;
                if atom.charge == -1
                    && atom.num_H == 0
                    && atom.valence == 1
                    && atom.chem_bonds_valence == 2
                    && valence.cMetal == 0
                    && valence.cNumValenceElectrons == 5
                    && minus_available
                    && !reversed_endpoint
                {
                    let second_nitrogen_index = usize::from(atom.neighbor[0]);
                    let second_nitrogen_valence = pVA
                        .get(second_nitrogen_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let second_nitrogen = heap
                        .slice(at2.as_const())?
                        .get(second_nitrogen_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if second_nitrogen_valence.cNumValenceElectrons == 5
                        && atom.bond_type[0] == BOND_TYPE_DOUBLE as u8
                        && second_nitrogen.charge == 1
                        && second_nitrogen.valence == 2
                        && second_nitrogen.chem_bonds_valence == 4
                    {
                        let carbon_neighbor_index = usize::from(
                            second_nitrogen.neighbor[usize::from(
                                second_nitrogen.neighbor[0] as usize == second_nitrogen_index,
                            )],
                        );
                        chain_matches = pVA
                            .get(carbon_neighbor_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .cNumValenceElectrons
                            == 4;
                    }
                }
                if chain_matches {
                    ret = AddToEdgeList(heap, &mut current_edges, minus_edge, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            let number_to_try = current_edges.num_edges.min(number_double_bond_oxygen);
            if number_to_try != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < number_double_bond_oxygen && current_success < number_to_try
                {
                    let atom_index =
                        usize::try_from(i32::from(double_bond_oxygen[candidate_index as usize]))
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
                    let bond_edge = *heap
                        .slice(atom_vertex.iedge.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let bond_edge_index = usize::try_from(bond_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(bond_edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    {
                        let edges = heap.slice_mut(pBNS.edge)?;
                        let edge = edges
                            .get_mut(bond_edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(1);
                        let target_edge = edges
                            .get_mut(
                                usize::try_from(target_minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        target_edge.forbidden =
                            (i32::from(target_edge.forbidden) & forbidden_edge_mask_inv) as i8;
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
                        && delta_charge == 0
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            n_num_run_bns = n_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[bond_edge_index].flow = edge_before.flow;
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
                    {
                        let edges = heap.slice_mut(pBNS.edge)?;
                        let edge = edges
                            .get_mut(bond_edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        let target_edge = edges
                            .get_mut(
                                usize::try_from(target_minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        target_edge.forbidden =
                            (i32::from(target_edge.forbidden) | forbidden_edge_mask) as i8;
                    }
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 23, ichirvr3.c:5606-5973.
        if comparison.len_c2at != 0 && comparison.nNumTgInChI == 1 {
            const CHG_SET_NOOH: usize = 0;
            const CHG_SET_WRONG_TAUT: usize = 1;
            const CHG_SET_TAUT: usize = 2;
            const CHG_LAST_SET: usize = 2;
            const CHG_SET_O_FIXED: usize = 3;
            const CHG_SET_NUM: usize = 4;

            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut nitrogen_dioxide = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_oxygen = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let reversed_fixed_h = heap
                .slice(pStruct.pOneINChI[0].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .nNum_H_fixed;
            let reversed_mobile_h = if !pStruct.pOneINChI[1].is_null() {
                let inchi = heap
                    .slice(pStruct.pOneINChI[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !inchi.nNum_H.is_null() {
                    inchi.nNum_H
                } else if !pStruct.pOneINChI[0].is_null() {
                    heap.slice(pStruct.pOneINChI[0].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNum_H
                } else {
                    SourceMutPointer::null()
                }
            } else if !pStruct.pOneINChI[0].is_null() {
                heap.slice(pStruct.pOneINChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H
            } else {
                SourceMutPointer::null()
            };
            let mut changeable_edges: [EDGE_LIST; CHG_SET_NUM] =
                std::array::from_fn(|_| EDGE_LIST::default());
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;
            let mut leave_case = false;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) && !leave_case {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let minus_available = minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                let mut dioxide_center = None;
                if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                    && difference.nValElectr == 6
                    && difference.endptInChI != 0
                    && minus_available
                    && difference.nFixHInChI == 0
                    && difference.endptRevrs == 0
                    && difference.nFixHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
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
                        && i32::from(nitrogen.chem_bonds_valence)
                            == 5_i32.wrapping_sub(i32::from(nitrogen.charge))
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
                    let mut neighbor_index = 0_i32;
                    while neighbor_index < i32::from(nitrogen.valence) {
                        let position = usize::try_from(neighbor_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let neighbor_number = i32::from(nitrogen.neighbor[position]);
                        if neighbor_number == atom_number {
                            neighbor_index = neighbor_index.wrapping_add(1);
                            continue;
                        }
                        let neighbor_atom_index = usize::try_from(neighbor_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let neighbor_atom = heap
                            .slice(at2.as_const())?
                            .get(neighbor_atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let neighbor_valence = pVA
                            .get(neighbor_atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let input_endpoint = *heap
                            .slice(pStruct.endpoint.as_const())?
                            .get(neighbor_atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let reversed_endpoint = !mobile_h_reversed.is_null()
                            && heap
                                .slice(mobile_h_reversed.as_const())?
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .endpoint
                                != 0;
                        if neighbor_valence.cNumValenceElectrons == 6
                            && input_endpoint != 0
                            && !reversed_endpoint
                            && neighbor_atom.valence == 1
                            && neighbor_atom.num_H == 0
                            && neighbor_atom.radical == 0
                            && (neighbor_atom.charge == 0 || neighbor_atom.charge == -1)
                            && i32::from(neighbor_atom.chem_bonds_valence)
                                .wrapping_sub(i32::from(neighbor_atom.charge))
                                == 2
                        {
                            number_oxygen = number_oxygen.wrapping_add(1);
                        } else if nitrogen.bond_type[position] == BOND_TYPE_SINGLE as u8
                            && neighbor_atom.valence > 1
                            && neighbor_atom.valence < neighbor_atom.chem_bonds_valence
                        {
                            number_others = number_others.wrapping_add(1);
                        }
                        neighbor_index = neighbor_index.wrapping_add(1);
                    }
                    if number_oxygen == 1 && number_others == 1 {
                        let nitrogen_number = nitrogen_index as i32;
                        let mut duplicate_index = 0_i32;
                        while duplicate_index < number_double_bond_oxygen
                            && i32::from(nitrogen_dioxide[duplicate_index as usize])
                                != nitrogen_number
                        {
                            duplicate_index = duplicate_index.wrapping_add(1);
                        }
                        if duplicate_index == number_double_bond_oxygen {
                            nitrogen_dioxide[number_double_bond_oxygen as usize] =
                                nitrogen_number as i16;
                            double_bond_oxygen[number_double_bond_oxygen as usize] =
                                atom_number as i16;
                            number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                        }
                        ret = AddToEdgeList(
                            heap,
                            &mut changeable_edges[CHG_SET_O_FIXED],
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

            if number_double_bond_oxygen != 0 && !leave_case {
                let mut canonical_number = 0_i32;
                while canonical_number < pStruct.num_atoms && !leave_case {
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let canonical_valence = pVA
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut nitrogen_match = None;
                    if *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0
                        && canonical_valence.cNumValenceElectrons == 6
                        && atom.valence == 1
                        && atom.charge == 0
                    {
                        let atom_minus_edge = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCMinusGroupEdge
                            .wrapping_sub(1);
                        if atom_minus_edge >= 0
                            && i32::from(atom.num_H)
                                .wrapping_add(i32::from(atom.chem_bonds_valence))
                                == 2
                        {
                            let nitrogen_index = usize::from(atom.neighbor[0]);
                            if pVA
                                .get(nitrogen_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .cNumValenceElectrons
                                == 5
                            {
                                let nitrogen_plus_edge = pVA
                                    .get(nitrogen_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .nCPlusGroupEdge
                                    .wrapping_sub(1);
                                if nitrogen_plus_edge >= 0 {
                                    let plus_edge = heap
                                        .slice(pBNS.edge.as_const())?
                                        .get(
                                            usize::try_from(nitrogen_plus_edge)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    let nitrogen = heap
                                        .slice(at2.as_const())?
                                        .get(nitrogen_index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    if plus_edge.forbidden == 0
                                        && plus_edge.flow != 0
                                        && nitrogen.charge == 0
                                        && nitrogen.valence == 3
                                        && nitrogen.chem_bonds_valence == 5
                                    {
                                        nitrogen_match = Some(nitrogen_plus_edge);
                                    }
                                }
                            }
                        }
                    }
                    if let Some(nitrogen_plus_edge) = nitrogen_match {
                        let nitrogen_index = usize::from(atom.neighbor[0]);
                        let nitrogen = heap
                            .slice(at2.as_const())?
                            .get(nitrogen_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let mut number_oxygen = 0_i32;
                        let mut number_others = 0_i32;
                        let mut neighbor_index = 0_i32;
                        while neighbor_index < i32::from(nitrogen.valence) {
                            let position = usize::try_from(neighbor_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let neighbor_atom_index = usize::from(nitrogen.neighbor[position]);
                            if neighbor_atom_index == atom_index {
                                neighbor_index = neighbor_index.wrapping_add(1);
                                continue;
                            }
                            let neighbor = heap
                                .slice(at2.as_const())?
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if pVA
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .cNumValenceElectrons
                                == 6
                                && *heap
                                    .slice(pStruct.endpoint.as_const())?
                                    .get(neighbor_atom_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    != 0
                                && neighbor.valence == 1
                                && neighbor.num_H == 1
                                && neighbor.radical == 0
                                && neighbor.charge == 0
                            {
                                number_oxygen = number_oxygen.wrapping_add(1);
                            } else if nitrogen.bond_type[position] == BOND_TYPE_DOUBLE as u8
                                && neighbor.valence >= 2
                                && neighbor.valence < neighbor.chem_bonds_valence
                            {
                                number_others = number_others.wrapping_add(1);
                            }
                            neighbor_index = neighbor_index.wrapping_add(1);
                        }
                        if number_oxygen == 1 && number_others == 1 {
                            ret = AddToEdgeList(
                                heap,
                                &mut changeable_edges[CHG_SET_NOOH],
                                nitrogen_plus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret == 0 {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut changeable_edges[CHG_SET_O_FIXED],
                                    nitrogen_plus_edge,
                                    INC_ADD_EDGE,
                                )?;
                            }
                            if ret != 0 {
                                leave_case = true;
                            } else {
                                let flower_edge =
                                    GetChargeFlowerUpperEdge(heap, pBNS, pVA, nitrogen_plus_edge)?;
                                if flower_edge != NO_VERTEX {
                                    ret = AddToEdgeList(
                                        heap,
                                        &mut changeable_edges[CHG_SET_NOOH],
                                        flower_edge,
                                        INC_ADD_EDGE,
                                    )?;
                                    if ret == 0 {
                                        ret = AddToEdgeList(
                                            heap,
                                            &mut changeable_edges[CHG_SET_O_FIXED],
                                            nitrogen_plus_edge,
                                            INC_ADD_EDGE,
                                        )?;
                                    }
                                    if ret != 0 {
                                        leave_case = true;
                                    }
                                }
                            }
                        }
                    }
                    canonical_number = canonical_number.wrapping_add(1);
                }

                canonical_number = 0;
                while canonical_number < pStruct.num_atoms && !leave_case {
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let reversed_endpoint = !mobile_h_reversed.is_null()
                        && heap
                            .slice(mobile_h_reversed.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .endpoint
                            != 0;
                    let reversed_mobile_atom = if mobile_h_reversed.is_null() {
                        None
                    } else {
                        Some(
                            heap.slice(mobile_h_reversed.as_const())?
                                .get(canonical_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                    };
                    let wrong_tautomer = atom.charge == -1
                        && input_endpoint == 0
                        && reversed_mobile_atom.is_some_and(|mobile| {
                            mobile.endpoint != 0 || atom.num_H < mobile.num_H
                        });
                    let normalized_fixed = atom.charge == -1
                        && input_endpoint != 0
                        && !reversed_endpoint
                        && !reversed_fixed_h.is_null()
                        && *heap
                            .slice(reversed_fixed_h.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            == -1
                        && !reversed_mobile_h.is_null()
                        && *heap
                            .slice(reversed_mobile_h.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            == 1
                        && *heap
                            .slice(pStruct.fixed_H.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            == 0;
                    if wrong_tautomer || normalized_fixed {
                        let minus_edge = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCMinusGroupEdge
                            .wrapping_sub(1);
                        if minus_edge >= 0
                            && FindInEdgeList(heap, &changeable_edges[CHG_SET_O_FIXED], minus_edge)?
                                < 0
                        {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge.forbidden == 0 && edge.flow != 0 {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut changeable_edges[CHG_SET_WRONG_TAUT],
                                    minus_edge,
                                    INC_ADD_EDGE,
                                )?;
                                if ret == 0 {
                                    ret = AddToEdgeList(
                                        heap,
                                        &mut changeable_edges[CHG_SET_O_FIXED],
                                        minus_edge,
                                        INC_ADD_EDGE,
                                    )?;
                                }
                                if ret != 0 {
                                    leave_case = true;
                                }
                            }
                        }
                    }
                    canonical_number = canonical_number.wrapping_add(1);
                }

                canonical_number = 0;
                while canonical_number < pStruct.num_atoms && !leave_case {
                    let canonical_index = usize::try_from(canonical_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let reversed_endpoint = !mobile_h_reversed.is_null()
                        && heap
                            .slice(mobile_h_reversed.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .endpoint
                            != 0;
                    let atom_number = i32::from(
                        *heap
                            .slice(canonical_to_atom.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap
                        .slice(at2.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if input_endpoint != 0 && reversed_endpoint && atom.charge == -1 {
                        let minus_edge = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCMinusGroupEdge
                            .wrapping_sub(1);
                        if minus_edge >= 0 {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge.forbidden == 0
                                && edge.flow != 0
                                && FindInEdgeList(
                                    heap,
                                    &changeable_edges[CHG_SET_O_FIXED],
                                    minus_edge,
                                )? < 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut changeable_edges[CHG_SET_TAUT],
                                    minus_edge,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    leave_case = true;
                                }
                            }
                        }
                    }
                    canonical_number = canonical_number.wrapping_add(1);
                }

                let mut target_index = 0_i32;
                while target_index < number_double_bond_oxygen && !leave_case {
                    let atom_index =
                        usize::try_from(i32::from(double_bond_oxygen[target_index as usize]))
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
                    let bond_edge = *heap
                        .slice(atom_vertex.iedge.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let bond_edge_index = usize::try_from(bond_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(bond_edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        target_index = target_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
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
                                    usize::try_from(vertex_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
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
                        let expected_delta_charge = if set_index == CHG_SET_NOOH { 2 } else { 0 };
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &changeable_edges[set_index],
                            forbidden_edge_mask,
                        )?;
                        let target_edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(
                                usize::try_from(target_minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        target_edge.forbidden =
                            (i32::from(target_edge.forbidden) & forbidden_edge_mask_inv) as i8;
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
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
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
                        edge.forbidden =
                            (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                    }
                    if one_success == 0 {
                        heap.slice_mut(pBNS.edge)?[bond_edge_index].flow = edge_before.flow;
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
                    target_index = target_index.wrapping_add(1);
                }
            }

            for edges in &mut changeable_edges {
                let _ = AllocEdgeList(heap, edges, EDGE_LIST_FREE)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 24, ichirvr3.c:5975-6262.
        if comparison.len_c2at != 0 && comparison.nNumTgInChI == 1 {
            const CHG_SET_MISSED_TAUT: usize = 0;
            const CHG_SET_OTHER_TAUT_O: usize = 1;
            const CHG_SET_OTHER_TAUT_N: usize = 2;
            const CHG_LAST_SET: usize = 2;
            const CHG_SET_NN: usize = 3;
            const CHG_SET_AVOID: usize = 4;
            const CHG_SET_NUM: usize = 5;

            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let fixed_bond_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at_fixed_bonds
            };
            let mut changeable_edges: [EDGE_LIST; CHG_SET_NUM] =
                std::array::from_fn(|_| EDGE_LIST::default());
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;
            let mut leave_case = false;

            let mut difference_index = 0_i32;
            while difference_index < i32::from(comparison.len_c2at) && !leave_case {
                let difference = comparison.c2at[usize::try_from(difference_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                let atom_number = i32::from(difference.atomNumber);
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let input_minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                let input_minus_available = input_minus_edge >= 0
                    && heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(input_minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .forbidden
                        == 0;
                let mut first_kind = None;
                if difference.nValElectr == 5
                    && difference.endptInChI != 0
                    && input_minus_available
                    && difference.nFixHInChI == 0
                    && difference.nMobHInChI == 0
                    && difference.endptRevrs == 0
                    && difference.nFixHRevrs == 0
                    && difference.nAtChargeRevrs == 0
                    && atom.num_H == 0
                    && atom.valence == 2
                    && atom.chem_bonds_valence == 3
                {
                    let nitrogen_index = usize::from(
                        atom.neighbor[usize::from(atom.bond_type[0] != BOND_TYPE_DOUBLE as u8)],
                    );
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
                        && nitrogen.chem_bonds_valence == 5
                        && nitrogen.charge == 0
                        && nitrogen.num_H == 0
                        && nitrogen.radical == 0
                    {
                        let plus_edge = pVA
                            .get(nitrogen_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCPlusGroupEdge
                            .wrapping_sub(1);
                        if plus_edge >= 0 {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(plus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge.forbidden == 0
                                && edge.flow != 0
                                && FindInEdgeList(
                                    heap,
                                    &changeable_edges[CHG_SET_AVOID],
                                    plus_edge,
                                )? < 0
                            {
                                first_kind = Some(plus_edge);
                            }
                        }
                    }
                }
                if let Some(plus_edge) = first_kind {
                    let flower_edge = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                    for (list_index, edge_number) in [
                        (CHG_SET_NN, plus_edge),
                        (CHG_SET_NN, flower_edge),
                        (CHG_SET_NN, 1),
                        (CHG_SET_AVOID, input_minus_edge),
                        (CHG_SET_AVOID, plus_edge),
                        (CHG_SET_AVOID, flower_edge),
                    ] {
                        if ret == 0 {
                            ret = AddToEdgeList(
                                heap,
                                &mut changeable_edges[list_index],
                                edge_number,
                                INC_ADD_EDGE,
                            )?;
                        }
                    }
                    if ret != 0 {
                        leave_case = true;
                    } else {
                        let raw_minus_edge = valence.nCMinusGroupEdge;
                        if raw_minus_edge >= 0 {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(
                                    usize::try_from(raw_minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if edge.forbidden == 0
                                && FindInEdgeList(
                                    heap,
                                    &changeable_edges[CHG_SET_AVOID],
                                    raw_minus_edge,
                                )? < 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut changeable_edges[CHG_SET_AVOID],
                                    raw_minus_edge,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    leave_case = true;
                                }
                            }
                        }
                    }
                } else if !fixed_bond_reversed.is_null()
                    && difference.nValElectr == 5
                    && difference.endptInChI != 0
                    && input_minus_available
                    && difference.nFixHInChI == 0
                    && difference.nMobHInChI == 0
                    && difference.endptRevrs == 0
                    && difference.nFixHRevrs == 0
                    && difference.nAtChargeRevrs == -1
                    && atom.num_H == 0
                    && atom.valence == 2
                    && atom.chem_bonds_valence == 2
                {
                    let fixed_atom = heap
                        .slice(fixed_bond_reversed.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if fixed_atom.valence == 2 && fixed_atom.chem_bonds_valence == 3 {
                        let nitrogen_index = usize::from(
                            fixed_atom.neighbor
                                [usize::from(fixed_atom.bond_type[0] != BOND_TYPE_DOUBLE as u8)],
                        );
                        let nitrogen = heap
                            .slice(at2.as_const())?
                            .get(nitrogen_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let fixed_nitrogen = heap
                            .slice(fixed_bond_reversed.as_const())?
                            .get(nitrogen_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(
                                usize::try_from(input_minus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if pVA
                            .get(nitrogen_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .cNumValenceElectrons
                            == 5
                            && nitrogen.charge == 1
                            && nitrogen.chem_bonds_valence == 4
                            && fixed_nitrogen.charge == 0
                            && fixed_nitrogen.chem_bonds_valence == 5
                            && nitrogen.num_H == 0
                            && nitrogen.radical == 0
                            && edge.forbidden == 0
                            && edge.flow != 0
                            && FindInEdgeList(
                                heap,
                                &changeable_edges[CHG_SET_AVOID],
                                input_minus_edge,
                            )? < 0
                        {
                            for (list_index, edge_number) in [
                                (CHG_SET_NN, input_minus_edge),
                                (CHG_SET_NN, NO_VERTEX),
                                (CHG_SET_NN, 1),
                                (CHG_SET_AVOID, input_minus_edge),
                            ] {
                                if ret == 0 {
                                    ret = AddToEdgeList(
                                        heap,
                                        &mut changeable_edges[list_index],
                                        edge_number,
                                        INC_ADD_EDGE,
                                    )?;
                                }
                            }
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }

            if changeable_edges[CHG_SET_NN].num_edges == 0 {
                leave_case = true;
            }
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms && !leave_case {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    == 0
                {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if atom.charge != 0 || atom.radical != 0 || atom.valence == atom.chem_bonds_valence
                {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                if minus_edge < 0 {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let edge = heap
                    .slice(pBNS.edge.as_const())?
                    .get(
                        usize::try_from(minus_edge)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if edge.forbidden != 0
                    || edge.flow != 0
                    || FindInEdgeList(heap, &changeable_edges[CHG_SET_AVOID], minus_edge)? >= 0
                {
                    canonical_number = canonical_number.wrapping_add(1);
                    continue;
                }
                let reversed_endpoint = !mobile_h_reversed.is_null()
                    && heap
                        .slice(mobile_h_reversed.as_const())?
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .endpoint
                        != 0;
                let destination_set = if !reversed_endpoint {
                    Some(CHG_SET_MISSED_TAUT)
                } else if valence.cNumValenceElectrons == 6 {
                    Some(CHG_SET_OTHER_TAUT_O)
                } else if valence.cNumValenceElectrons == 5 {
                    Some(CHG_SET_OTHER_TAUT_N)
                } else {
                    None
                };
                if let Some(destination_set) = destination_set
                    && FindInEdgeList(heap, &changeable_edges[CHG_SET_AVOID], minus_edge)? < 0
                {
                    ret = AddToEdgeList(
                        heap,
                        &mut changeable_edges[destination_set],
                        minus_edge,
                        INC_ADD_EDGE,
                    )?;
                    if ret == 0 {
                        ret = AddToEdgeList(
                            heap,
                            &mut changeable_edges[CHG_SET_AVOID],
                            minus_edge,
                            INC_ADD_EDGE,
                        )?;
                    }
                    if ret != 0 {
                        leave_case = true;
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            if !leave_case {
                let mut triple_index = 0_i32;
                while triple_index < changeable_edges[CHG_SET_NN].num_edges {
                    let edge_number = *heap
                        .slice(changeable_edges[CHG_SET_NN].pnEdges.as_const())?
                        .get(
                            usize::try_from(triple_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let flower_edge = *heap
                        .slice(changeable_edges[CHG_SET_NN].pnEdges.as_const())?
                        .get(
                            usize::try_from(triple_index.wrapping_add(1))
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let expected_delta_charge = *heap
                        .slice(changeable_edges[CHG_SET_NN].pnEdges.as_const())?
                        .get(
                            usize::try_from(triple_index.wrapping_add(2))
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        triple_index = triple_index.wrapping_add(3);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow.wrapping_sub(1);
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
                    let mut one_success = 0_i32;
                    let mut set_index = 0_usize;
                    while one_success == 0 && set_index <= CHG_LAST_SET {
                        if changeable_edges[set_index].num_edges == 0 {
                            set_index += 1;
                            continue;
                        }
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &changeable_edges[set_index],
                            forbidden_edge_mask,
                        )?;
                        if flower_edge != NO_VERTEX {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(
                                    usize::try_from(flower_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden =
                                (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
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
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                one_success = one_success.wrapping_add(1);
                            }
                        }
                        set_index += 1;
                    }
                    current_success = current_success.wrapping_add(one_success);
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    if one_success == 0 {
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow;
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
                    triple_index = triple_index.wrapping_add(3);
                }
            }

            for edges in &mut changeable_edges {
                let _ = AllocEdgeList(heap, edges, EDGE_LIST_FREE)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        // Complete translation of source case 25, ichirvr3.c:6264-6520.
        if comparison.len_c2at != 0
            && comparison.nNumTgInChI == 1
            && comparison.nNumRemHRevrs > comparison.nNumRemHInChI
            && comparison.nNumRemHInChI < 0
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI
                || comparison.nNumTgRevrs > comparison.nNumTgInChI)
        {
            const CHG_SET_MISSED_TAUT_1: usize = 0;
            const CHG_SET_MISSED_TAUT_ALL: usize = 1;
            const CHG_SET_OTHER_TAUT_1: usize = 2;
            const CHG_SET_OTHER_TAUT_ALL: usize = 3;
            const CHG_LAST_SET: usize = 3;
            const CHG_SET_NO_IN_NO2M2: usize = 4;
            const CHG_SET_AVOID: usize = 5;
            const CHG_SET_NUM: usize = 6;

            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let atom_to_canonical = pStruct.nAtno2Canon[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[1].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[1].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let mut changeable_edges: [EDGE_LIST; CHG_SET_NUM] =
                std::array::from_fn(|_| EDGE_LIST::default());
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;
            let mut leave_case = false;
            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms && !leave_case {
                let canonical_index = usize::try_from(canonical_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let input_endpoint = *heap
                    .slice(pStruct.endpoint.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if input_endpoint != 0 {
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    if minus_edge < 0 {
                        canonical_number = canonical_number.wrapping_add(1);
                        continue;
                    }
                    let edge = heap
                        .slice(pBNS.edge.as_const())?
                        .get(
                            usize::try_from(minus_edge)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge.forbidden != 0
                        || FindInEdgeList(heap, &changeable_edges[CHG_SET_AVOID], minus_edge)? >= 0
                    {
                        canonical_number = canonical_number.wrapping_add(1);
                        continue;
                    }
                    let first = (valence.cNumValenceElectrons == 5 && comparison.nNumTgInChI == 1)
                        || (valence.cNumValenceElectrons == 6 && comparison.nNumTgInChI != 1);
                    let reversed_endpoint = !mobile_h_reversed.is_null()
                        && heap
                            .slice(mobile_h_reversed.as_const())?
                            .get(canonical_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .endpoint
                            != 0;
                    if !reversed_endpoint {
                        if first {
                            ret = AddToEdgeList(
                                heap,
                                &mut changeable_edges[CHG_SET_MISSED_TAUT_1],
                                minus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                        if !leave_case {
                            ret = AddToEdgeList(
                                heap,
                                &mut changeable_edges[CHG_SET_MISSED_TAUT_ALL],
                                minus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                    }
                    if !leave_case && first {
                        ret = AddToEdgeList(
                            heap,
                            &mut changeable_edges[CHG_SET_OTHER_TAUT_1],
                            minus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            leave_case = true;
                        }
                    }
                    if !leave_case {
                        ret = AddToEdgeList(
                            heap,
                            &mut changeable_edges[CHG_SET_OTHER_TAUT_ALL],
                            minus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            leave_case = true;
                        }
                    }
                    if !leave_case {
                        ret = AddToEdgeList(
                            heap,
                            &mut changeable_edges[CHG_SET_AVOID],
                            minus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            leave_case = true;
                        }
                    }
                } else if atom.valence == 1
                    && atom.charge == -1
                    && valence.cNumValenceElectrons == 6
                {
                    let nitrogen_index = usize::from(atom.neighbor[0]);
                    let nitrogen_valence = pVA
                        .get(nitrogen_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let nitrogen = heap
                        .slice(at2.as_const())?
                        .get(nitrogen_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let nitrogen_canonical = usize::from(
                        *heap
                            .slice(atom_to_canonical.as_const())?
                            .get(nitrogen_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let nitrogen_input_endpoint = *heap
                        .slice(pStruct.endpoint.as_const())?
                        .get(nitrogen_canonical)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let nitrogen_plus_edge = nitrogen_valence.nCPlusGroupEdge.wrapping_sub(1);
                    let oxygen_minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    let charge_edges_match = nitrogen_plus_edge >= 0
                        && oxygen_minus_edge >= 0
                        && {
                            let edges = heap.slice(pBNS.edge.as_const())?;
                            let plus = edges
                                .get(
                                    usize::try_from(nitrogen_plus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let minus = edges
                                .get(
                                    usize::try_from(oxygen_minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            plus.forbidden == 0
                                && plus.flow != 0
                                && minus.forbidden == 0
                                && minus.flow != 0
                        }
                        && FindInEdgeList(
                            heap,
                            &changeable_edges[CHG_SET_AVOID],
                            nitrogen_plus_edge,
                        )? < 0
                        && FindInEdgeList(
                            heap,
                            &changeable_edges[CHG_SET_AVOID],
                            oxygen_minus_edge,
                        )? < 0;
                    if nitrogen_valence.cNumValenceElectrons == 5
                        && nitrogen_input_endpoint == 0
                        && nitrogen.valence == 3
                        && nitrogen.chem_bonds_valence == 3
                        && nitrogen.charge == 0
                        && nitrogen.radical == 0
                        && charge_edges_match
                    {
                        let mut number_oxygen = 0_i32;
                        let mut number_others = 0_i32;
                        let mut neighbor_index = 0_i32;
                        while neighbor_index < i32::from(nitrogen.valence) {
                            let position = usize::try_from(neighbor_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let neighbor_atom_index = usize::from(nitrogen.neighbor[position]);
                            if neighbor_atom_index == atom_index {
                                neighbor_index = neighbor_index.wrapping_add(1);
                                continue;
                            }
                            let neighbor = heap
                                .slice(at2.as_const())?
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if pVA
                                .get(neighbor_atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .cNumValenceElectrons
                                == 6
                                && *heap
                                    .slice(pStruct.endpoint.as_const())?
                                    .get(neighbor_atom_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    == 0
                                && neighbor.valence == 1
                                && neighbor.num_H == 0
                                && neighbor.radical == 0
                                && neighbor.charge == -1
                                && neighbor.chem_bonds_valence == 1
                            {
                                number_oxygen = number_oxygen.wrapping_add(1);
                            } else if nitrogen.bond_type[position] == BOND_TYPE_SINGLE as u8
                                && neighbor.valence > 1
                                && neighbor.valence < neighbor.chem_bonds_valence
                            {
                                number_others = number_others.wrapping_add(1);
                            }
                            neighbor_index = neighbor_index.wrapping_add(1);
                        }
                        if !(number_oxygen != 1 && number_others != 1) {
                            for (list_index, edge_number) in [
                                (CHG_SET_NO_IN_NO2M2, nitrogen_plus_edge),
                                (CHG_SET_NO_IN_NO2M2, oxygen_minus_edge),
                                (CHG_SET_AVOID, nitrogen_plus_edge),
                                (CHG_SET_AVOID, oxygen_minus_edge),
                            ] {
                                if ret == 0 {
                                    ret = AddToEdgeList(
                                        heap,
                                        &mut changeable_edges[list_index],
                                        edge_number,
                                        INC_ADD_EDGE,
                                    )?;
                                }
                            }
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            if changeable_edges[CHG_SET_NO_IN_NO2M2].num_edges == 0
                || changeable_edges[CHG_SET_OTHER_TAUT_ALL].num_edges == 0
            {
                leave_case = true;
            }
            if !leave_case {
                let mut pair_index = 0_i32;
                while pair_index < changeable_edges[CHG_SET_NO_IN_NO2M2].num_edges {
                    let nitrogen_plus_edge = *heap
                        .slice(changeable_edges[CHG_SET_NO_IN_NO2M2].pnEdges.as_const())?
                        .get(
                            usize::try_from(pair_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let oxygen_minus_edge = *heap
                        .slice(changeable_edges[CHG_SET_NO_IN_NO2M2].pnEdges.as_const())?
                        .get(
                            usize::try_from(pair_index.wrapping_add(1))
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let oxygen_minus_index = usize::try_from(oxygen_minus_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(oxygen_minus_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        pair_index = pair_index.wrapping_add(2);
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    heap.slice_mut(pBNS.edge)?[oxygen_minus_index].flow =
                        edge_before.flow.wrapping_sub(1);
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
                    let mut one_success = 0_i32;
                    let mut set_index = 0_usize;
                    while one_success == 0 && set_index <= CHG_LAST_SET {
                        if changeable_edges[set_index].num_edges == 0 {
                            set_index += 1;
                            continue;
                        }
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(
                            heap,
                            pBNS,
                            &changeable_edges[set_index],
                            forbidden_edge_mask,
                        )?;
                        let nitrogen_edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(
                                usize::try_from(nitrogen_plus_edge)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        nitrogen_edge.forbidden =
                            (i32::from(nitrogen_edge.forbidden) & forbidden_edge_mask_inv) as i8;
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
                            && delta_charge == 3
                        {
                            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                            if ret > 0 {
                                n_num_run_bns = n_num_run_bns.wrapping_add(1);
                                one_success = one_success.wrapping_add(1);
                            }
                        }
                        set_index += 1;
                    }
                    current_success = current_success.wrapping_add(one_success);
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    if one_success == 0 {
                        heap.slice_mut(pBNS.edge)?[oxygen_minus_index].flow = edge_before.flow;
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
                    pair_index = pair_index.wrapping_add(2);
                }
            }

            for edges in &mut changeable_edges {
                let _ = AllocEdgeList(heap, edges, EDGE_LIST_FREE)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                tot_success = tot_success.wrapping_add(current_success);
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
                ret =
                    FillOutCMP2FHINCHI(heap, pStruct, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }
        Ok(())
    })();

    let cleanup_result = (|| -> Result<(), SourceHeapError> {
        let _ = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut current_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut nitrogen_flower_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut sulfur_flower_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut other_nitrogen_flower_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut fixed_large_ring_stereo_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut all_bond_edges, EDGE_LIST_FREE)?;
        Ok(())
    })();

    execution?;
    cleanup_result?;
    Ok(if ret < 0 {
        ret
    } else {
        i32::from(tot_success != 0 && comparison.bHasDifference != 0)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        BNS_EDGE, BNS_VERTEX, INChI_Aux, RI_ERR_ALLOC, RI_ERR_PROGR, T_GROUP,
    };

    fn call_fix_fixed_h_restored_structure(
        heap: &mut SourceHeap,
        structure: &mut StrFromINChI,
        atoms: SourceMutPointer<inp_ATOM>,
        bns: &mut BN_STRUCT,
        valence: &mut [VAL_AT],
        input: [SourceMutPointer<INChI>; 2],
        runs: &mut i32,
        total_delta: &mut i32,
    ) -> Result<i32, SourceHeapError> {
        FixFixedHRestoredStructure(
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
            &mut ALL_TC_GROUPS::default(),
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

    fn atom(charge: i8, num_h: i8, valence: i8) -> inp_ATOM {
        inp_ATOM {
            charge,
            num_H: num_h,
            valence,
            chem_bonds_valence: valence,
            ..inp_ATOM::default()
        }
    }

    fn make_info(
        heap: &mut SourceHeap,
        groups: Vec<T_GROUP>,
        endpoints: Vec<AT_NUMB>,
    ) -> T_GROUP_INFO {
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
    fn source_port__ichirvr3__filltgdiffhchgfh__line_113() {
        let mut empty_heap = SourceHeap::default();
        let empty_info = T_GROUP_INFO::default();
        let mut byte_exact = vec![TgDiffHChgFH {
            itg: i16::from_ne_bytes([0x7f, 0x7f]),
            nNumHInchi: i16::from_ne_bytes([0x7f, 0x7f]),
            ..TgDiffHChgFH::default()
        }];
        let mut empty_indices = EDGE_LIST::default();
        assert_eq!(
            FillTgDiffHChgFH(
                &mut empty_heap,
                &mut byte_exact,
                1,
                &[],
                &[],
                &[],
                &[],
                &empty_info,
                &mut empty_indices,
            ),
            Ok(0)
        );
        assert_eq!(byte_exact[0].itg.to_ne_bytes(), [0, 0x7f]);
        assert_eq!(byte_exact[0].nNumHInchi.to_ne_bytes(), [0x7f, 0x7f]);
        assert_eq!(empty_indices.num_edges, 0);

        let mut heap = SourceHeap::default();
        let groups = vec![
            T_GROUP {
                num: [99, 2, 0, 0, 0],
                nNumEndpoints: 5,
                nFirstEndpointAtNoPos: 0,
                ..T_GROUP::default()
            },
            T_GROUP {
                num: [9, 1, 0, 0, 0],
                nNumEndpoints: 1,
                nFirstEndpointAtNoPos: 5,
                ..T_GROUP::default()
            },
        ];
        let info = make_info(&mut heap, groups, vec![0, 1, 2, 3, 4, 5]);
        let at2 = vec![
            atom(1, 1, 1),
            atom(1, 0, 1),
            atom(-1, 0, 1),
            atom(-1, 0, 1),
            atom(0, 1, 1),
            atom(1, 0, 1),
        ];
        let atf = at2.clone();
        let mut valence = vec![VAL_AT::default(); 6];
        valence[0].nCPlusGroupEdge = 1;
        valence[1].nCPlusGroupEdge = 1;
        valence[2].cNumValenceElectrons = 6;
        valence[3].cNumValenceElectrons = 5;
        valence[4].nCPlusGroupEdge = 1;
        valence[5].nCPlusGroupEdge = 1;
        let canonical = vec![0, 1, 2, 3, 4, 5];
        let mut differences = vec![TgDiffHChgFH::default(); 2];
        let mut atom_indices = EDGE_LIST::default();
        assert_eq!(
            FillTgDiffHChgFH(
                &mut heap,
                &mut differences,
                2,
                &at2,
                &atf,
                &canonical,
                &valence,
                &info,
                &mut atom_indices,
            ),
            Ok(1)
        );
        assert_eq!(differences[0].itg, 0);
        assert_eq!(differences[0].nNumHInchi, 97);
        assert_eq!(differences[0].nNumMInchi, 2);
        assert_eq!(differences[0].nNumHRevrs, 2);
        assert_eq!(differences[0].nNumHNorml, 2);
        assert_eq!(differences[0].nNumMRevrs, 2);
        assert_eq!(differences[0].nNumMNorml, 2);
        assert_eq!(differences[0].nNumPRevrs, 2);
        assert_eq!(differences[0].nNumPNorml, 2);
        assert_eq!(differences[0].n, [1; F_NUM_ALL_CHG_T]);
        assert_eq!(differences[0].i, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(atom_indices.num_edges, 10);
        assert_eq!(
            &heap.slice(atom_indices.pnEdges.as_const()).unwrap()[..10],
            &[0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
        );

        let mut correct_heap = SourceHeap::default();
        let correct_info = make_info(
            &mut correct_heap,
            vec![T_GROUP {
                num: [2, 1, 0, 0, 0],
                nNumEndpoints: 1,
                ..T_GROUP::default()
            }],
            vec![0],
        );
        let reverse = vec![atom(-1, 7, 1)];
        let normal = vec![atom(-1, 1, 1)];
        let mut oxygen = vec![VAL_AT::default()];
        oxygen[0].cNumValenceElectrons = 6;
        let mut correct_differences = vec![TgDiffHChgFH::default()];
        let mut correct_indices = EDGE_LIST::default();
        assert_eq!(
            FillTgDiffHChgFH(
                &mut correct_heap,
                &mut correct_differences,
                1,
                &reverse,
                &normal,
                &[0],
                &oxygen,
                &correct_info,
                &mut correct_indices,
            ),
            Ok(0)
        );
        assert_eq!(correct_indices.num_edges, 0);

        let mut failure_heap = SourceHeap::default();
        let failure_info = make_info(
            &mut failure_heap,
            vec![T_GROUP {
                num: [8, 0, 0, 0, 0],
                nNumEndpoints: 1,
                ..T_GROUP::default()
            }],
            vec![0],
        );
        failure_heap.fail_after_allocations(0);
        let mut failure_differences = vec![TgDiffHChgFH::default()];
        let mut failure_indices = EDGE_LIST::default();
        assert_eq!(
            FillTgDiffHChgFH(
                &mut failure_heap,
                &mut failure_differences,
                1,
                &[atom(1, 1, 1)],
                &[atom(1, 1, 1)],
                &[0],
                &[VAL_AT {
                    nCPlusGroupEdge: 1,
                    ..VAL_AT::default()
                }],
                &failure_info,
                &mut failure_indices,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(failure_indices.num_edges, 0);
    }

    #[test]
    fn source_port__ichirvr3__bhas_n_v__line_92() {
        let matching = inp_ATOM {
            el_number: EL_NUMBER_N,
            chem_bonds_valence: 5,
            valence: 3,
            ..inp_ATOM::default()
        };
        let mut atoms = vec![matching.clone(), matching.clone()];

        let mut wrong_element = matching.clone();
        wrong_element.el_number = EL_NUMBER_N - 1;
        atoms.push(wrong_element);
        let mut charged = matching.clone();
        charged.charge = -1;
        atoms.push(charged);
        let mut hydrogenated = matching.clone();
        hydrogenated.num_H = -1;
        atoms.push(hydrogenated);
        let mut radical = matching.clone();
        radical.radical = 1;
        atoms.push(radical);
        let mut wrong_chemical_valence = matching.clone();
        wrong_chemical_valence.chem_bonds_valence = 4;
        atoms.push(wrong_chemical_valence);
        let mut wrong_valence = matching;
        wrong_valence.valence = 4;
        atoms.push(wrong_valence);

        assert_eq!(bHas_N_V(&atoms, -1), Ok(0));
        assert_eq!(bHas_N_V(&atoms, 0), Ok(0));
        assert_eq!(bHas_N_V(&atoms, 1), Ok(1));
        assert_eq!(bHas_N_V(&atoms, 2), Ok(2));
        assert_eq!(bHas_N_V(&atoms, atoms.len() as i32), Ok(2));
        assert_eq!(
            bHas_N_V(&atoms, atoms.len() as i32 + 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr3__fixfixedhrestoredstructure__line_333() {
        fn fixed_h_inchi(fixed_h: SourceMutPointer<i8>) -> INChI {
            INChI {
                nNum_H_fixed: fixed_h,
                ..INChI::default()
            }
        }

        fn zero_atom_fixture(
            heap: &mut SourceHeap,
            with_mobile_layer: bool,
        ) -> (
            StrFromINChI,
            SourceMutPointer<inp_ATOM>,
            [SourceMutPointer<INChI>; 2],
        ) {
            let fixed_h = heap.allocate_model_storage(vec![0_i8]).unwrap();
            let input_fixed = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h)])
                .unwrap();
            let reversed_fixed = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h)])
                .unwrap();
            let input_mobile = if with_mobile_layer {
                heap.allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: 1,
                    ..INChI::default()
                }])
                .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let canonical_numbers = heap.allocate_model_storage(vec![1_u16]).unwrap();
            let auxiliary = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nOrigAtNosInCanonOrd: canonical_numbers,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let atoms = heap.allocate_model_storage(Vec::<inp_ATOM>::new()).unwrap();
            (
                StrFromINChI {
                    num_atoms: 0,
                    pOneINChI: [reversed_fixed, SourceMutPointer::null()],
                    pOneINChI_Aux: [auxiliary, SourceMutPointer::null()],
                    ..StrFromINChI::default()
                },
                atoms,
                [input_fixed, input_mobile],
            )
        }

        let mut no_fixed_heap = SourceHeap::default();
        let input = no_fixed_heap
            .allocate_model_storage(vec![INChI::default()])
            .unwrap();
        let reversed = no_fixed_heap
            .allocate_model_storage(vec![INChI::default()])
            .unwrap();
        let mut no_fixed_structure = StrFromINChI {
            pOneINChI: [reversed, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let no_fixed_live = no_fixed_heap.live_allocation_count();
        let mut runs = i32::MIN;
        let mut total_delta = i32::MAX;
        assert_eq!(
            call_fix_fixed_h_restored_structure(
                &mut no_fixed_heap,
                &mut no_fixed_structure,
                SourceMutPointer::null(),
                &mut BN_STRUCT::default(),
                &mut [],
                [input, SourceMutPointer::null()],
                &mut runs,
                &mut total_delta,
            ),
            Ok(0)
        );
        assert_eq!((runs, total_delta), (i32::MIN, i32::MAX));
        assert_eq!(no_fixed_heap.live_allocation_count(), no_fixed_live);

        for with_mobile_layer in [false, true] {
            let mut heap = SourceHeap::default();
            let (mut structure, atoms, input) = zero_atom_fixture(&mut heap, with_mobile_layer);
            let fixture_allocations = heap.live_allocation_count();
            heap.trace_source_allocations();
            let mut runs = 17;
            let mut total_delta = -19;
            assert_eq!(
                call_fix_fixed_h_restored_structure(
                    &mut heap,
                    &mut structure,
                    atoms,
                    &mut BN_STRUCT::default(),
                    &mut [],
                    input,
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(0),
                "mobile layer={with_mobile_layer}"
            );
            assert_eq!((runs, total_delta), (17, -19));
            assert_eq!(heap.source_allocation_calls(), 2);
            assert_eq!(heap.live_allocation_count(), fixture_allocations + 2);
            assert!(!structure.nCanon2Atno[0].is_null());
            assert!(!structure.nAtno2Canon[0].is_null());
        }

        for successful_allocations in [0_u64, 1] {
            let mut heap = SourceHeap::default();
            let (mut structure, atoms, input) = zero_atom_fixture(&mut heap, false);
            let fixture_allocations = heap.live_allocation_count();
            heap.fail_after_allocations(successful_allocations);
            let mut runs = 23;
            let mut total_delta = 29;
            assert_eq!(
                call_fix_fixed_h_restored_structure(
                    &mut heap,
                    &mut structure,
                    atoms,
                    &mut BN_STRUCT::default(),
                    &mut [],
                    input,
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(RI_ERR_ALLOC),
                "allocation ordinal={successful_allocations}"
            );
            assert_eq!((runs, total_delta), (23, 29));
            assert_eq!(
                heap.live_allocation_count(),
                fixture_allocations + successful_allocations as usize
            );
            assert_eq!(
                structure.nCanon2Atno[0].is_null(),
                successful_allocations == 0
            );
            assert!(structure.nAtno2Canon[0].is_null());
        }

        for fail_charge_list_allocation in [true, false] {
            let mut heap = SourceHeap::default();
            let fixed_h = heap.allocate_model_storage(vec![0_i8]).unwrap();
            let input = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h)])
                .unwrap();
            let reversed = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h)])
                .unwrap();
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let vertices = heap
                .allocate_model_storage(vec![BNS_VERTEX::default()])
                .unwrap();
            let edges = heap
                .allocate_model_storage(vec![BNS_EDGE::default()])
                .unwrap();
            let mut structure = StrFromINChI {
                num_atoms: 1,
                pOneINChI: [reversed, SourceMutPointer::null()],
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
            if fail_charge_list_allocation {
                heap.fail_after_allocations(0);
            } else {
                heap.trace_source_allocations();
            }
            let mut runs = 31;
            let mut total_delta = 37;
            assert_eq!(
                call_fix_fixed_h_restored_structure(
                    &mut heap,
                    &mut structure,
                    atoms,
                    &mut bns,
                    &mut valence,
                    [input, SourceMutPointer::null()],
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(if fail_charge_list_allocation {
                    RI_ERR_ALLOC
                } else {
                    RI_ERR_PROGR
                })
            );
            assert_eq!((runs, total_delta), (31, 37));
            assert_eq!(heap.source_allocation_calls(), 1);
            assert_eq!(heap.live_allocation_count(), fixture_allocations);
            assert_eq!(
                heap.slice(edges.as_const()).unwrap(),
                &[BNS_EDGE::default()]
            );
        }
    }
}
