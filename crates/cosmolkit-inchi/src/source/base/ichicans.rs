use crate::source::base::ichimap2::HalfStereoBondParity;
use crate::source::base::ichisort::{
    CompNeighborsAT_NUMBER, CompNeighborsATNumberContext, comp_AT_RANK, insertions_sort,
};
use crate::source_types::{
    AB_INV_PARITY_BITS, AB_PARITY_CALC, AB_PARITY_EVEN, AB_PARITY_NONE, AB_PARITY_UNDF, AT_NUMB,
    AT_RANK, AT_STEREO_CARB, AT_STEREO_DBLE, CANON_GLOBALS, CANON_STAT, CT_ERR_MAX, CT_ERR_MIN,
    CT_OUT_OF_RAM, CT_OVERFLOW, CT_STEREOBOND_ERROR, CT_STEREOCOUNT_ERR, MAX_ATOMS,
    MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BOND_NEIGH, MAX_NUM_STEREO_BONDS, MAXVAL,
    MIN_DOT_PROD, NUM_H_ISOTOPES, S_CHAR, SourceConstPointer, SourceHeap, SourceHeapError,
    SourceMutPointer, sp_ATOM,
};

const BITS_PARITY: i32 = 0x07;
const MASK_CUMULENE_LEN: i32 = 0x38;
const MULT_STEREOBOND: i32 = 0x08;
const AB_MIN_WELL_DEFINED_PARITY: i32 = 1;
const AB_MAX_WELL_DEFINED_PARITY: i32 = 2;

fn atom_parity_well_defined(parity: i32) -> bool {
    (AB_MIN_WELL_DEFINED_PARITY..=AB_MAX_WELL_DEFINED_PARITY).contains(&parity)
}

fn parity_well_defined(parity: i32) -> bool {
    let parity = parity & BITS_PARITY;
    (AB_MIN_WELL_DEFINED_PARITY..=AB_MAX_WELL_DEFINED_PARITY).contains(&parity)
}

fn parity_calculate(parity: i32) -> bool {
    parity & BITS_PARITY == AB_PARITY_CALC as i32
}

fn parity_known(parity: i32) -> bool {
    (1..=4).contains(&(parity & BITS_PARITY))
}

fn bond_chain_len(parity: i32) -> i32 {
    (parity & MASK_CUMULENE_LEN) / MULT_STEREOBOND
}

pub(crate) fn find_atoms_with_parity(
    heap: &mut SourceHeap,
    at: SourceConstPointer<sp_ATOM>,
    visited: SourceMutPointer<S_CHAR>,
    from_atom: i32,
    cur_atom: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:155 find_atoms_with_parity
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int find_atoms_with_parity( sp_ATOM *at,
                                S_CHAR  *visited,
                                int     from_atom,
                                int     cur_atom )
    {
        int i, next_atom;

        if (visited[cur_atom])
        {
            return 0;
        }
        if (at[cur_atom].parity)
        {
            return 1;
        }

        visited[cur_atom] = 1;

        for (i = 0; i < at[cur_atom].valence; i++)
        {
            next_atom = at[cur_atom].neighbor[i];

            if (next_atom != from_atom &&
                 find_atoms_with_parity( at, visited, cur_atom, next_atom ))
            {
                return 1;
            }
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: find_atoms_with_parity

    let current = usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if heap
        .slice(visited.as_const())?
        .get(current)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        != 0
    {
        return Ok(0);
    }
    let atom = heap
        .slice(at)?
        .get(current)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.parity != 0 {
        return Ok(1);
    }
    heap.slice_mut(visited)?[current] = 1;
    for i in 0..i32::from(atom.valence) {
        let next_atom = i32::from(atom.neighbor[i as usize]);
        if next_atom != from_atom
            && find_atoms_with_parity(heap, at, visited, cur_atom, next_atom)? != 0
        {
            return Ok(1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn SetHalfStereoBondIllDefPariy(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    jn: i32,
    k1: i32,
    new_parity: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:189 SetHalfStereoBondIllDefPariy
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int SetHalfStereoBondIllDefPariy( sp_ATOM *at,
                                      int     jn, /* atom number*/
                                      int     k1 /* stereo bond number*/,
                                      int     new_parity )
    {
        int parity;
        if (k1 < MAX_NUM_STEREO_BOND_NEIGH && at[jn].stereo_bond_neighbor[k1])
        {
            parity = at[jn].stereo_bond_parity[k1] ^ PARITY_VAL( at[jn].stereo_bond_parity[k1] );
            at[jn].stereo_bond_parity[k1] = parity | PARITY_VAL( new_parity );
            at[jn].parity = PARITY_VAL( new_parity );
            return 1;  /*  success */
        }

        return 0; /*  failed             */
    }
    */
    // END INCHI C FUNCTION: SetHalfStereoBondIllDefPariy
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/extr_ct.h:297 PARITY_VAL
    // INCHI✔️❌: #define PARITY_VAL(X)               ((X) & BITS_PARITY)
    // END INCHI ACTIVE HEADER MACRO: PARITY_VAL

    if k1 >= MAX_NUM_STEREO_BOND_NEIGH as i32 {
        return Ok(0);
    }
    let atom_index = usize::try_from(jn).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let stereo_index = usize::try_from(k1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = heap
        .slice_mut(at)?
        .get_mut(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.stereo_bond_neighbor[stereo_index] == 0 {
        return Ok(0);
    }
    let encoded = i32::from(atom.stereo_bond_parity[stereo_index]);
    let parity = encoded ^ (encoded & BITS_PARITY);
    atom.stereo_bond_parity[stereo_index] = (parity | (new_parity & BITS_PARITY)) as i8;
    atom.parity = (new_parity & BITS_PARITY) as i8;
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveHalfStereoBond(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    jn: i32,
    k1: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:208 RemoveHalfStereoBond
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int RemoveHalfStereoBond( sp_ATOM *at,
                              int     jn, /* atom number*/
                              int     k1 /* stereo bond number*/ )
    {
        int k2;
        if (k1 < MAX_NUM_STEREO_BOND_NEIGH && at[jn].stereo_bond_neighbor[k1])
        {
            for (k2 = k1; k2 < MAX_NUM_STEREO_BOND_NEIGH - 1; k2++) /* djb-rwth: loop condition corrected (buffer error) */
            {
                at[jn].stereo_bond_neighbor[k2] = at[jn].stereo_bond_neighbor[k2 + 1];
                at[jn].stereo_bond_ord[k2] = at[jn].stereo_bond_ord[k2 + 1];
                at[jn].stereo_bond_z_prod[k2] = at[jn].stereo_bond_z_prod[k2 + 1];
                at[jn].stereo_bond_parity[k2] = at[jn].stereo_bond_parity[k2 + 1];
            }
            at[jn].stereo_bond_neighbor[k2] = 0;
            at[jn].stereo_bond_ord[k2] = 0;
            at[jn].stereo_bond_z_prod[k2] = 0;
            at[jn].stereo_bond_parity[k2] = 0;

            if (!at[jn].stereo_bond_neighbor[0])
            {   /*  curled braces added 6-6-2002 */
                at[jn].parity = 0;
                at[jn].stereo_atom_parity = 0;
                at[jn].final_parity = 0;
                /* at[jn].bHasStereoOrEquToStereo = 0; */
            }
            return 1; /*  success            */
        }

        return 0; /*  failed             */
    }
    */
    // END INCHI C FUNCTION: RemoveHalfStereoBond

    let atom_index = usize::try_from(jn).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let stereo_index = usize::try_from(k1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if stereo_index >= MAX_NUM_STEREO_BOND_NEIGH as usize {
        return Ok(0);
    }
    let atom = heap
        .slice_mut(at)?
        .get_mut(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.stereo_bond_neighbor[stereo_index] == 0 {
        return Ok(0);
    }
    let last = MAX_NUM_STEREO_BOND_NEIGH as usize - 1;
    for k2 in stereo_index..last {
        atom.stereo_bond_neighbor[k2] = atom.stereo_bond_neighbor[k2 + 1];
        atom.stereo_bond_ord[k2] = atom.stereo_bond_ord[k2 + 1];
        atom.stereo_bond_z_prod[k2] = atom.stereo_bond_z_prod[k2 + 1];
        atom.stereo_bond_parity[k2] = atom.stereo_bond_parity[k2 + 1];
    }
    atom.stereo_bond_neighbor[last] = 0;
    atom.stereo_bond_ord[last] = 0;
    atom.stereo_bond_z_prod[last] = 0;
    atom.stereo_bond_parity[last] = 0;
    if atom.stereo_bond_neighbor[0] == 0 {
        atom.parity = 0;
        atom.stereo_atom_parity = 0;
        atom.final_parity = 0;
    }
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn SetOneStereoBondIllDefParity(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    jc: i32,
    k: i32,
    new_parity: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:242 SetOneStereoBondIllDefParity
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int SetOneStereoBondIllDefParity( sp_ATOM *at,
                                      int     jc,       /* atom number              */
                                      int     k,        /* stereo bond ord. number  */
                                      int     new_parity )
    {
        int k1, ret = 0, kn, jn = (int) at[jc].stereo_bond_neighbor[k] - 1;

        /*  opposite end */
        for (k1 = ret = 0;
                k1 < MAX_NUM_STEREO_BOND_NEIGH && ( kn = at[jn].stereo_bond_neighbor[k1] );
                  k1++) /* djb-rwth: removing redundant code */
        {
            if (kn - 1 == jc)
            {
                ret = SetHalfStereoBondIllDefPariy( at, jn, /* atom number*/ k1 /* stereo bond number*/, new_parity );
                break;
            }
        }

        if (ret)
        {
            ret = SetHalfStereoBondIllDefPariy( at, jc, k, new_parity );
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: SetOneStereoBondIllDefParity

    let current_index = usize::try_from(jc).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let stereo_index = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let jn = i32::from(
        *heap
            .slice(at.as_const())?
            .get(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .stereo_bond_neighbor
            .get(stereo_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    )
    .wrapping_sub(1);

    let mut ret = 0_i32;
    let mut k1 = 0_i32;
    while k1 < MAX_NUM_STEREO_BOND_NEIGH as i32 {
        let opposite_index =
            usize::try_from(jn).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let kn = heap
            .slice(at.as_const())?
            .get(opposite_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .stereo_bond_neighbor[k1 as usize];
        if kn == 0 {
            break;
        }
        if i32::from(kn).wrapping_sub(1) == jc {
            ret = SetHalfStereoBondIllDefPariy(heap, at, jn, k1, new_parity)?;
            break;
        }
        k1 = k1.wrapping_add(1);
    }
    if ret != 0 {
        ret = SetHalfStereoBondIllDefPariy(heap, at, jc, k, new_parity)?;
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveOneStereoBond(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    jc: i32,
    k: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:271 RemoveOneStereoBond
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int RemoveOneStereoBond( sp_ATOM *at,
                             int     jc,    /* atom number          */
                             int     k      /* stereo bond number   */
    )
    {
        int k1, ret = 0, kn, jn = (int) at[jc].stereo_bond_neighbor[k] - 1;

        /*  opposite end */
        for (k1 = ret = 0;
                k1 < MAX_NUM_STEREO_BOND_NEIGH && ( kn = at[jn].stereo_bond_neighbor[k1] );
                  k1++) /* djb-rwth: removing redundant code */
        {
            if (kn - 1 == jc)
            {
                ret = RemoveHalfStereoBond( at, jn, k1 );
                break;
            }
        }

        if (ret)
        {
            ret = RemoveHalfStereoBond( at, jc, k );
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: RemoveOneStereoBond

    let current_index = usize::try_from(jc).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let stereo_index = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let opposite_index = {
        let atoms = heap.slice(at.as_const())?;
        let atom = atoms
            .get(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        usize::from(
            *atom
                .stereo_bond_neighbor
                .get(stereo_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
        .checked_sub(1)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
    };

    let mut ret = 0_i32;
    for opposite_stereo_index in 0..MAX_NUM_STEREO_BOND_NEIGH as usize {
        let opposite_neighbor = {
            let atoms = heap.slice(at.as_const())?;
            atoms
                .get(opposite_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .stereo_bond_neighbor[opposite_stereo_index]
        };
        if opposite_neighbor == 0 {
            break;
        }
        if i32::from(opposite_neighbor).wrapping_sub(1) == jc {
            ret = RemoveHalfStereoBond(
                heap,
                at,
                opposite_index as i32,
                opposite_stereo_index as i32,
            )?;
            break;
        }
    }
    if ret != 0 {
        ret = RemoveHalfStereoBond(heap, at, jc, k)?;
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveOneStereoCenter(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    jc: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:300 RemoveOneStereoCenter
    // INCHI✔️✔️: int RemoveOneStereoCenter( sp_ATOM *at,
    // INCHI✔️✔️:                            int     jc /* atom number*/ )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (at[jc].parity)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         at[jc].parity = 0; /*  remove parity */
    // INCHI✔️✔️:         at[jc].stereo_atom_parity = 0;
    // INCHI✔️✔️:         at[jc].final_parity = 0;
    // INCHI✔️✔️:         /*  at[jc].bHasStereoOrEquToStereo = 0; */
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0; /*  failed: not a stereo center */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: RemoveOneStereoCenter

    let index = usize::try_from(jc).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = heap
        .slice_mut(at)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.parity != 0 {
        atom.parity = 0;
        atom.stereo_atom_parity = 0;
        atom.final_parity = 0;
        Ok(1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn InvertStereo(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_at_tg: i32,
    nCanonRank: SourceConstPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
    bInvertLinearCTStereo: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2016 InvertStereo
    // INCHI✔️✔️: int InvertStereo( sp_ATOM    *at,
    // INCHI✔️✔️:                   int        num_at_tg,
    // INCHI✔️✔️:                   AT_RANK    *nCanonRank,
    // INCHI✔️✔️:                   AT_RANK    *nAtomNumberCanon,
    // INCHI✔️✔️:                   CANON_STAT *pCS,
    // INCHI✔️✔️:                   int        bInvertLinearCTStereo )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j, j1, j2, num_changes, parity, cumulene_len;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     num_changes = 0;
    // INCHI✔️✔️:     for (i = 0; i < num_at_tg; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nAtomNumberCanon[(int) nCanonRank[i] - 1] = i;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i = 0; i < pCS->nLenLinearCTStereoCarb; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         parity = pCS->LinearCTStereoCarb[i].parity;
    // INCHI✔️✔️:         if (ATOM_PARITY_WELL_DEF( parity ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             j = nAtomNumberCanon[(int) pCS->LinearCTStereoCarb[i].at_num - 1];
    // INCHI✔️✔️:             if (PARITY_WELL_DEF( at[j].parity ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[j].parity ^= AB_INV_PARITY_BITS;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 goto exit_error; /* inconsistency */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (bInvertLinearCTStereo)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️: #ifdef FIX_STEREOCOUNT_ERR
    // INCHI✔️✔️:                 pCS->LinearCTStereoCarb[i].parity = AB_PARITY_EVEN; /* deliberately worse */
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:                 pCS->LinearCTStereoCarb[i].parity ^= AB_INV_PARITY_BITS;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             num_changes++;
    // INCHI✔️✔️:             if (PARITY_WELL_DEF( at[j].stereo_atom_parity ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[j].stereo_atom_parity ^= AB_INV_PARITY_BITS;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (PARITY_WELL_DEF( at[j].final_parity ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[j].final_parity ^= AB_INV_PARITY_BITS;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i = 0; i < pCS->nLenLinearCTStereoDble; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         parity = pCS->LinearCTStereoDble[i].parity;
    // INCHI✔️✔️:         if (ATOM_PARITY_WELL_DEF( parity ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             j1 = nAtomNumberCanon[(int) pCS->LinearCTStereoDble[i].at_num1 - 1];
    // INCHI✔️✔️:             cumulene_len = BOND_CHAIN_LEN( at[j1].stereo_bond_parity[0] );
    // INCHI✔️✔️:             if (cumulene_len % 2)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* invert only in case of allene */
    // INCHI✔️✔️:                 j2 = nAtomNumberCanon[(int) pCS->LinearCTStereoDble[i].at_num2 - 1];
    // INCHI✔️✔️:                 /* checks for debug only */
    // INCHI✔️✔️:                 if (1 < MAX_NUM_STEREO_BONDS)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (at[j1].stereo_bond_neighbor[1] ||
    // INCHI✔️✔️:                          at[j2].stereo_bond_neighbor[1])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         goto exit_error; /* inconsitency: atom has more than one cumulene bond */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (cumulene_len != BOND_CHAIN_LEN( at[j2].stereo_bond_parity[0] ) ||
    // INCHI✔️✔️:                      j1 + 1 != at[j2].stereo_bond_neighbor[0] ||
    // INCHI✔️✔️:                      j2 + 1 != at[j1].stereo_bond_neighbor[0])
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     goto exit_error; /* inconsitency: atoms should refer to each other */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /* invert parities */
    // INCHI✔️✔️:                 if (PARITY_WELL_DEF( at[j1].parity ) && PARITY_WELL_DEF( at[j2].parity ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     j = inchi_min( j1, j2 );
    // INCHI✔️✔️:                     at[j].parity ^= AB_INV_PARITY_BITS; /* for reversability always invert only atom with the smaller number */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     goto exit_error; /* inconsistency */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (bInvertLinearCTStereo)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     pCS->LinearCTStereoDble[i].parity ^= AB_INV_PARITY_BITS;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 num_changes++;
    // INCHI✔️✔️:                 if (PARITY_WELL_DEF( at[j1].stereo_bond_parity[0] ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     at[j1].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (PARITY_WELL_DEF( at[j2].stereo_bond_parity[0] ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     at[j2].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_changes;
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_error:
    // INCHI✔️✔️:     return CT_STEREOCOUNT_ERR;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: InvertStereo

    for i in 0..num_at_tg {
        let rank = *heap
            .slice(nCanonRank)?
            .get(i as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let output_index = usize::from(rank)
            .checked_sub(1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(nAtomNumberCanon)?
            .get_mut(output_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = i as AT_RANK;
    }

    let mut num_changes = 0_i32;
    for i in 0..pCS.nLenLinearCTStereoCarb {
        let index = i as usize;
        let record = heap
            .slice(pCS.LinearCTStereoCarb.as_const())?
            .get(index)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom_parity_well_defined(i32::from(record.parity)) {
            let canonical_index = usize::from(record.at_num)
                .checked_sub(1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom_index = usize::from(
                *heap
                    .slice(nAtomNumberCanon.as_const())?
                    .get(canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let atoms = heap.slice_mut(at)?;
            let atom = atoms
                .get_mut(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if !parity_well_defined(i32::from(atom.parity)) {
                return Ok(CT_STEREOCOUNT_ERR);
            }
            atom.parity ^= AB_INV_PARITY_BITS as i8;
            if bInvertLinearCTStereo != 0 {
                heap.slice_mut(pCS.LinearCTStereoCarb)?[index].parity = AB_PARITY_EVEN as u8;
            }
            num_changes = num_changes.wrapping_add(1);
            let atom = heap
                .slice_mut(at)?
                .get_mut(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if parity_well_defined(i32::from(atom.stereo_atom_parity)) {
                atom.stereo_atom_parity ^= AB_INV_PARITY_BITS as i8;
            }
            if parity_well_defined(i32::from(atom.final_parity)) {
                atom.final_parity ^= AB_INV_PARITY_BITS as i8;
            }
        }
    }

    for i in 0..pCS.nLenLinearCTStereoDble {
        let index = i as usize;
        let record = heap
            .slice(pCS.LinearCTStereoDble.as_const())?
            .get(index)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom_parity_well_defined(i32::from(record.parity)) {
            let first_canonical_index = usize::from(record.at_num1)
                .checked_sub(1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let j1 = usize::from(
                *heap
                    .slice(nAtomNumberCanon.as_const())?
                    .get(first_canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let first = heap
                .slice(at.as_const())?
                .get(j1)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let cumulene_len = bond_chain_len(i32::from(first.stereo_bond_parity[0]));
            if cumulene_len % 2 != 0 {
                let second_canonical_index = usize::from(record.at_num2)
                    .checked_sub(1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let j2 = usize::from(
                    *heap
                        .slice(nAtomNumberCanon.as_const())?
                        .get(second_canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let second = heap
                    .slice(at.as_const())?
                    .get(j2)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if 1 < MAX_NUM_STEREO_BONDS
                    && (first.stereo_bond_neighbor[1] != 0 || second.stereo_bond_neighbor[1] != 0)
                {
                    return Ok(CT_STEREOCOUNT_ERR);
                }
                if cumulene_len != bond_chain_len(i32::from(second.stereo_bond_parity[0]))
                    || j1 + 1 != usize::from(second.stereo_bond_neighbor[0])
                    || j2 + 1 != usize::from(first.stereo_bond_neighbor[0])
                {
                    return Ok(CT_STEREOCOUNT_ERR);
                }
                if !parity_well_defined(i32::from(first.parity))
                    || !parity_well_defined(i32::from(second.parity))
                {
                    return Ok(CT_STEREOCOUNT_ERR);
                }
                heap.slice_mut(at)?[j1.min(j2)].parity ^= AB_INV_PARITY_BITS as i8;
                if bInvertLinearCTStereo != 0 {
                    heap.slice_mut(pCS.LinearCTStereoDble)?[index].parity ^=
                        AB_INV_PARITY_BITS as u8;
                }
                num_changes = num_changes.wrapping_add(1);
                let atoms = heap.slice_mut(at)?;
                if parity_well_defined(i32::from(atoms[j1].stereo_bond_parity[0])) {
                    atoms[j1].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS as i8;
                }
                if parity_well_defined(i32::from(atoms[j2].stereo_bond_parity[0])) {
                    atoms[j2].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS as i8;
                }
            }
        }
    }

    Ok(num_changes)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillSingleStereoDescriptors(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    i: i32,
    mut num_trans: i32,
    nRank: SourceMutPointer<AT_RANK>,
    LinearCTStereoCarb: SourceMutPointer<AT_STEREO_CARB>,
    nStereoCarbLen: &mut i32,
    nMaxStereoCarbLen: i32,
    LinearCTStereoDble: SourceMutPointer<AT_STEREO_DBLE>,
    nStereoDbleLen: &mut i32,
    nMaxStereoDbleLen: i32,
    bAllene: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:525 FillSingleStereoDescriptors
    // INCHI✔️✔️: complete active source frame follows verbatim; the direct slice comparator view preserves the C stack-array behavior without allocation.
    /*
    int FillSingleStereoDescriptors( CANON_GLOBALS  *pCG,
                                     sp_ATOM        *at,
                                     int            i,
                                     int            num_trans,
                                     const AT_RANK  *nRank,
                                     AT_STEREO_CARB *LinearCTStereoCarb,
                                     int            *nStereoCarbLen,
                                     int            nMaxStereoCarbLen,
                                     AT_STEREO_DBLE *LinearCTStereoDble,
                                     int            *nStereoDbleLen,
                                     int            nMaxStereoDbleLen,
                                     int            bAllene )
    {

        if (!LinearCTStereoDble && !LinearCTStereoCarb)
        {
            return 0; /*  return immediately if no stereo have been requested */
        }

        /***************************************************
        add stereo centers and stereo bonds to the CT
        ***************************************************/
        if (at[i].parity || at[i].stereo_bond_neighbor[0])
        {
            AT_RANK r_neigh, rank = nRank[i];
            AT_NUMB nNeighborNumber2[MAXVAL];
            unsigned parity;
            int      k;
            int num_allene = 0;

            if (ATOM_PARITY_WELL_DEF( at[i].parity ) && num_trans < 0)
            {
                /*  number of neighbors transpositions to the sorted order is unknown. Find it. */
                /*  If parity is not well-defined then doing this is a waste of time */
                int num_neigh = at[i].valence;
                for (k = 0; k < num_neigh; k++)
                {
                    nNeighborNumber2[k] = k;
                }

                pCG->m_pNeighborsForSort = at[i].neighbor;
                pCG->m_pn_RankForSort = nRank;
                num_trans = insertions_sort( pCG, nNeighborNumber2, num_neigh, sizeof( nNeighborNumber2[0] ), CompNeighborsAT_NUMBER );

    #ifndef CT_NEIGH_INCREASE
                num_trans += ( ( num_neigh*( num_neigh - 1 ) ) / 2 ) % 2;  /*  get correct parity for ascending order */
    #endif
            }

            /*  stereo bonds */
            if (LinearCTStereoDble && at[i].stereo_bond_neighbor[0])
            {

                /* HalfStereoBondParity( sp_ATOM *at, int at_no1, int i_sb_neigh, AT_RANK *nRank ) */
                AT_NUMB nStereoNeighNumber[MAX_NUM_STEREO_BONDS], nStereoNeigh[MAX_NUM_STEREO_BONDS], n;
                int       num_stereo, stereo_neigh, stereo_neigh_ord, stereo_bond_parity;
                for (num_stereo = 0;
                          num_stereo < MAX_NUM_STEREO_BONDS &&
                          ( n = at[i].stereo_bond_neighbor[num_stereo] ); num_stereo++)
                {
                    nStereoNeighNumber[num_stereo] = num_stereo;
                    nStereoNeigh[num_stereo] = n - 1;
                    num_allene += IS_ALLENE_CHAIN( at[i].stereo_bond_parity[num_stereo] );
                }
                if ((bAllene > 0 && !num_allene) || (bAllene == 0 && num_allene)) /* djb-rwth: addressing LLVM warning */
                {
                    return 0;
                }

                /*  sort stereo bonds according to the ranks of the neighbors */
                pCG->m_pNeighborsForSort = nStereoNeigh;
                pCG->m_pn_RankForSort = nRank;
                insertions_sort( pCG, nStereoNeighNumber, num_stereo, sizeof( nStereoNeighNumber[0] ), CompNeighborsAT_NUMBER );

                /*  process stereo bonds one by one */
                for (k = 0; k < num_stereo; k++)
                {
                    stereo_neigh = nStereoNeigh[stereo_neigh_ord = (int) nStereoNeighNumber[k]];

                    if (( r_neigh = (AT_NUMB) nRank[stereo_neigh] ) CT_NEIGH_SMALLER_THAN rank)
                    {
                        /* accept only neighbors that have smaller ranks */
                        stereo_bond_parity = PARITY_VAL( at[i].stereo_bond_parity[stereo_neigh_ord] );
                        if (stereo_bond_parity == AB_PARITY_NONE)
                        {
                            continue;
                        }

                        /* stereo_neigh      = at[i].stereo_bond_neighbor[nStereoNeighNumber[k]]-1; */
                        if (ATOM_PARITY_KNOWN( stereo_bond_parity ))
                        {
                            parity = stereo_bond_parity;
                        }
                        else if (ATOM_PARITY_WELL_DEF( at[i].parity ) &&
                                  ATOM_PARITY_WELL_DEF( at[stereo_neigh].parity ) &&
                                  MIN_DOT_PROD <= abs( at[i].stereo_bond_z_prod[stereo_neigh_ord] ))
                        {
                            /*  bond parity can be calculated */
                            int half_parity1, half_parity2, j, nn, stereo_neigh_ord2;
                            stereo_neigh_ord2 = -1;
                            for (j = 0; j < MAX_NUM_STEREO_BONDS &&
                                ( nn = (int) at[stereo_neigh].stereo_bond_neighbor[j] );
                                           j++)
                            {
                                if (i + 1 == nn)
                                {
                                    /* found the opposite end of the stereo bond */
                                    stereo_neigh_ord2 = j;
                                    break;
                                }
                            }
                            if (stereo_neigh_ord2 >= 0)
                            {
                                half_parity1 = HalfStereoBondParity( at, i, stereo_neigh_ord, nRank );
                                half_parity2 = HalfStereoBondParity( at, stereo_neigh, stereo_neigh_ord2, nRank );
                                if (ATOM_PARITY_WELL_DEF( half_parity1 ) &&
                                     ATOM_PARITY_WELL_DEF( half_parity2 ))
                                {
                                    parity = 2 - ( half_parity1 + half_parity2
                                             + ( at[i].stereo_bond_z_prod[stereo_neigh_ord] < 0 ) ) % 2;
                                }
                                else
                                {
                                    return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                                }
                            }
                            else
                            {
                                return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                            }
                        }
                        else
                        {
                            /*  parity cannot be calculated: not enough info or 'unknown' */
                            if (AB_PARITY_NONE == ( parity = inchi_max( at[i].parity, at[stereo_neigh].parity ) ))
                            {
                                continue;
                            }
                            if (ATOM_PARITY_WELL_DEF( parity ))
                            {
                                parity = AB_PARITY_UNDF; /*  should not happen */
                            }
                        }
                        if (CHECK_OVERFLOW( *nStereoDbleLen, nMaxStereoDbleLen ))
                            return CT_OVERFLOW;  /*   <BRKPT> */
                        /*  first stereo bond atom */
                        LinearCTStereoDble[*nStereoDbleLen].at_num1 = rank;
                        /*  second stereo bond atom (opposite end) */
                        LinearCTStereoDble[*nStereoDbleLen].at_num2 = r_neigh;
                        /*  bond parity */
                        LinearCTStereoDble[*nStereoDbleLen].parity = parity;
                        ( *nStereoDbleLen )++;
                    }
                }
            }

            /*  stereo carbon */
            if (bAllene > 0)
            {
                return 0;
            }

            if (LinearCTStereoCarb && !at[i].stereo_bond_neighbor[0])
            {
                if (CHECK_OVERFLOW( *nStereoCarbLen, nMaxStereoCarbLen ))
                    return CT_OVERFLOW;  /*   <BRKPT> */
                /*  stereo atom rank */
                LinearCTStereoCarb[*nStereoCarbLen].at_num = rank;
                /*  stereo atom parity */
                parity = ATOM_PARITY_WELL_DEF( at[i].parity ) ? ( 2 - ( at[i].parity + num_trans ) % 2 ) : at[i].parity;
                LinearCTStereoCarb[*nStereoCarbLen].parity = parity;
                ( *nStereoCarbLen )++;
            }
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: FillSingleStereoDescriptors

    fn sort_neighbor_orders(
        values: &mut [AT_NUMB],
        count: usize,
        neighbors: &[AT_NUMB],
        ranks: &[AT_RANK],
    ) -> Result<i32, SourceHeapError> {
        let bytes = bytemuck::cast_slice_mut::<AT_NUMB, u8>(values);
        insertions_sort(
            bytes,
            count,
            std::mem::size_of::<AT_NUMB>(),
            &mut |left, right| {
                let left = AT_NUMB::from_ne_bytes([left[0], left[1]]);
                let right = AT_NUMB::from_ne_bytes([right[0], right[1]]);
                CompNeighborsAT_NUMBER(
                    left,
                    right,
                    CompNeighborsATNumberContext::Slices { neighbors, ranks },
                )
            },
        )
    }

    if LinearCTStereoDble.is_null() && LinearCTStereoCarb.is_null() {
        return Ok(0);
    }

    let atom_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom_snapshot = heap
        .slice(at.as_const())?
        .get(atom_index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom_snapshot.parity == 0 && atom_snapshot.stereo_bond_neighbor[0] == 0 {
        return Ok(0);
    }

    let rank = *heap
        .slice(nRank.as_const())?
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut neighbor_numbers = [0_u16; MAXVAL as usize];
    if atom_parity_well_defined(i32::from(atom_snapshot.parity)) && num_trans < 0 {
        let neighbor_count = usize::try_from(i32::from(atom_snapshot.valence))
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        if neighbor_count > neighbor_numbers.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        for (index, value) in neighbor_numbers[..neighbor_count].iter_mut().enumerate() {
            *value = index as AT_NUMB;
        }
        pCG.m_pn_RankForSort = nRank.as_const();
        num_trans = sort_neighbor_orders(
            &mut neighbor_numbers,
            neighbor_count,
            &atom_snapshot.neighbor,
            heap.slice(nRank.as_const())?,
        )?;
    }

    if !LinearCTStereoDble.is_null() && atom_snapshot.stereo_bond_neighbor[0] != 0 {
        let mut stereo_neighbor_numbers = [0_u16; MAX_NUM_STEREO_BONDS as usize];
        let mut stereo_neighbors = [0_u16; MAX_NUM_STEREO_BONDS as usize];
        let mut stereo_count = 0_usize;
        let mut num_allene = 0_i32;
        while stereo_count < MAX_NUM_STEREO_BONDS as usize {
            let neighbor = atom_snapshot.stereo_bond_neighbor[stereo_count];
            if neighbor == 0 {
                break;
            }
            stereo_neighbor_numbers[stereo_count] = stereo_count as AT_NUMB;
            stereo_neighbors[stereo_count] = neighbor.wrapping_sub(1);
            let encoded_parity = i32::from(atom_snapshot.stereo_bond_parity[stereo_count]);
            num_allene = num_allene.wrapping_add(
                ((encoded_parity & MASK_CUMULENE_LEN as i32) / MULT_STEREOBOND as i32) % 2,
            );
            stereo_count += 1;
        }
        if (bAllene > 0 && num_allene == 0) || (bAllene == 0 && num_allene != 0) {
            return Ok(0);
        }

        pCG.m_pn_RankForSort = nRank.as_const();
        sort_neighbor_orders(
            &mut stereo_neighbor_numbers,
            stereo_count,
            &stereo_neighbors,
            heap.slice(nRank.as_const())?,
        )?;

        for &stereo_order_value in &stereo_neighbor_numbers[..stereo_count] {
            let stereo_order = usize::from(stereo_order_value);
            let stereo_neighbor = usize::from(stereo_neighbors[stereo_order]);
            let neighbor_rank = *heap
                .slice(nRank.as_const())?
                .get(stereo_neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor_rank >= rank {
                continue;
            }

            let stereo_bond_parity =
                i32::from(atom_snapshot.stereo_bond_parity[stereo_order]) & BITS_PARITY;
            if stereo_bond_parity == AB_PARITY_NONE as i32 {
                continue;
            }
            let parity = if (1..=4).contains(&stereo_bond_parity) {
                stereo_bond_parity
            } else {
                let opposite_atom = heap
                    .slice(at.as_const())?
                    .get(stereo_neighbor)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if atom_parity_well_defined(i32::from(atom_snapshot.parity))
                    && atom_parity_well_defined(i32::from(opposite_atom.parity))
                    && i32::from(atom_snapshot.stereo_bond_z_prod[stereo_order]).abs()
                        >= MIN_DOT_PROD as i32
                {
                    let mut opposite_order = None;
                    for reverse_order in 0..MAX_NUM_STEREO_BONDS as usize {
                        let reverse_neighbor = opposite_atom.stereo_bond_neighbor[reverse_order];
                        if reverse_neighbor == 0 {
                            break;
                        }
                        if i.wrapping_add(1) == i32::from(reverse_neighbor) {
                            opposite_order = Some(reverse_order);
                            break;
                        }
                    }
                    let opposite_order = match opposite_order {
                        Some(value) => value,
                        None => return Ok(CT_STEREOBOND_ERROR),
                    };
                    let half_parity1 =
                        HalfStereoBondParity(heap, at, i, stereo_order as i32, nRank)?;
                    let half_parity2 = HalfStereoBondParity(
                        heap,
                        at,
                        stereo_neighbor as i32,
                        opposite_order as i32,
                        nRank,
                    )?;
                    if !atom_parity_well_defined(half_parity1)
                        || !atom_parity_well_defined(half_parity2)
                    {
                        return Ok(CT_STEREOBOND_ERROR);
                    }
                    2_i32.wrapping_sub(
                        half_parity1
                            .wrapping_add(half_parity2)
                            .wrapping_add(i32::from(
                                atom_snapshot.stereo_bond_z_prod[stereo_order] < 0,
                            ))
                            % 2,
                    )
                } else {
                    let mut fallback =
                        i32::from(atom_snapshot.parity).max(i32::from(opposite_atom.parity));
                    if fallback == AB_PARITY_NONE as i32 {
                        continue;
                    }
                    if atom_parity_well_defined(fallback) {
                        fallback = AB_PARITY_UNDF as i32;
                    }
                    fallback
                }
            };

            if *nStereoDbleLen >= nMaxStereoDbleLen {
                return Ok(CT_OVERFLOW);
            }
            let descriptor_index = usize::try_from(*nStereoDbleLen)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let descriptor = heap
                .slice_mut(LinearCTStereoDble)?
                .get_mut(descriptor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            descriptor.at_num1 = rank;
            descriptor.at_num2 = neighbor_rank;
            descriptor.parity = parity as u8;
            *nStereoDbleLen = nStereoDbleLen.wrapping_add(1);
        }
    }

    if bAllene > 0 {
        return Ok(0);
    }
    if !LinearCTStereoCarb.is_null() && atom_snapshot.stereo_bond_neighbor[0] == 0 {
        if *nStereoCarbLen >= nMaxStereoCarbLen {
            return Ok(CT_OVERFLOW);
        }
        let descriptor_index =
            usize::try_from(*nStereoCarbLen).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let descriptor = heap
            .slice_mut(LinearCTStereoCarb)?
            .get_mut(descriptor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        descriptor.at_num = rank;
        descriptor.parity = if atom_parity_well_defined(i32::from(atom_snapshot.parity)) {
            2_i32.wrapping_sub(i32::from(atom_snapshot.parity).wrapping_add(num_trans) % 2) as u8
        } else {
            atom_snapshot.parity as u8
        };
        *nStereoCarbLen = nStereoCarbLen.wrapping_add(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn SwitchAtomStereoAndIsotopicStereo(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    bSwitched: &mut i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:705 SwitchAtomStereoAndIsotopicStereo
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    void SwitchAtomStereoAndIsotopicStereo( sp_ATOM *at,
                                            int     num_atoms,
                                            int     *bSwitched )
    {
        int i;
        /*  switch atom stereo data */
        for (i = 0; i < num_atoms; i++)
        {
            inchi_swap( (char*) &at[i].parity, (char*) &at[i].parity2, sizeof( at[i].parity ) );
            inchi_swap( (char*) &at[i].final_parity, (char*) &at[i].final_parity2, sizeof( at[i].final_parity ) );
            inchi_swap( (char*) &at[i].stereo_atom_parity, (char*) &at[i].stereo_atom_parity2, sizeof( at[i].stereo_atom_parity ) );
            inchi_swap( (char*) &at[i].bHasStereoOrEquToStereo, (char*) &at[i].bHasStereoOrEquToStereo2, sizeof( at[i].bHasStereoOrEquToStereo ) );

            inchi_swap( (char*) at[i].stereo_bond_neighbor, (char*) at[i].stereo_bond_neighbor2, sizeof( at[i].stereo_bond_neighbor ) );
            inchi_swap( (char*) at[i].stereo_bond_ord, (char*) at[i].stereo_bond_ord2, sizeof( at[i].stereo_bond_ord ) );
            inchi_swap( (char*) at[i].stereo_bond_z_prod, (char*) at[i].stereo_bond_z_prod2, sizeof( at[i].stereo_bond_z_prod ) );
            inchi_swap( (char*) at[i].stereo_bond_parity, (char*) at[i].stereo_bond_parity2, sizeof( at[i].stereo_bond_parity ) );
        }

        *bSwitched = !*bSwitched;
    }
    */
    // END INCHI C FUNCTION: SwitchAtomStereoAndIsotopicStereo

    for i in 0..num_atoms {
        let atom = heap
            .slice_mut(at)?
            .get_mut(i as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        std::mem::swap(&mut atom.parity, &mut atom.parity2);
        std::mem::swap(&mut atom.final_parity, &mut atom.final_parity2);
        std::mem::swap(&mut atom.stereo_atom_parity, &mut atom.stereo_atom_parity2);
        std::mem::swap(
            &mut atom.bHasStereoOrEquToStereo,
            &mut atom.bHasStereoOrEquToStereo2,
        );
        std::mem::swap(
            &mut atom.stereo_bond_neighbor,
            &mut atom.stereo_bond_neighbor2,
        );
        std::mem::swap(&mut atom.stereo_bond_ord, &mut atom.stereo_bond_ord2);
        std::mem::swap(&mut atom.stereo_bond_z_prod, &mut atom.stereo_bond_z_prod2);
        std::mem::swap(&mut atom.stereo_bond_parity, &mut atom.stereo_bond_parity2);
    }
    *bSwitched = i32::from(*bSwitched == 0);
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn SetCtToIsotopicStereo(pCS: &mut CANON_STAT, pCS2: &CANON_STAT) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:729 SetCtToIsotopicStereo
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    void SetCtToIsotopicStereo( CANON_STAT *pCS,
                                CANON_STAT *pCS2 )
    {
        pCS->LinearCTStereoDble = pCS2->LinearCTIsotopicStereoDble; /*  enable stereo */
        pCS->LinearCTStereoCarb = pCS2->LinearCTIsotopicStereoCarb;

        pCS->LinearCTStereoDbleInv = pCS2->LinearCTIsotopicStereoDbleInv; /*  enable inv. stereo */
        pCS->LinearCTStereoCarbInv = pCS2->LinearCTIsotopicStereoCarbInv;
        pCS->nMaxLenLinearCTStereoDble = pCS2->nMaxLenLinearCTIsotopicStereoDble;
        pCS->nMaxLenLinearCTStereoCarb = pCS2->nMaxLenLinearCTIsotopicStereoCarb;

        pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
        pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
    }
    */
    // END INCHI C FUNCTION: SetCtToIsotopicStereo

    pCS.LinearCTStereoDble = pCS2.LinearCTIsotopicStereoDble;
    pCS.LinearCTStereoCarb = pCS2.LinearCTIsotopicStereoCarb;
    pCS.LinearCTStereoDbleInv = pCS2.LinearCTIsotopicStereoDbleInv;
    pCS.LinearCTStereoCarbInv = pCS2.LinearCTIsotopicStereoCarbInv;
    pCS.nMaxLenLinearCTStereoDble = pCS2.nMaxLenLinearCTIsotopicStereoDble;
    pCS.nMaxLenLinearCTStereoCarb = pCS2.nMaxLenLinearCTIsotopicStereoCarb;
    pCS.nLenLinearCTStereoDble = pCS2.nLenLinearCTIsotopicStereoDble;
    pCS.nLenLinearCTStereoCarb = pCS2.nLenLinearCTIsotopicStereoCarb;
}

#[allow(non_snake_case)]
pub(crate) fn SetCtToNonIsotopicStereo(pCS: &mut CANON_STAT, pCS2: &CANON_STAT) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:746 SetCtToNonIsotopicStereo
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    void SetCtToNonIsotopicStereo( CANON_STAT *pCS,
                                   CANON_STAT *pCS2 )
    {
        pCS->LinearCTStereoDble = pCS2->LinearCTStereoDble; /*  enable stereo */
        pCS->LinearCTStereoCarb = pCS2->LinearCTStereoCarb;

        pCS->LinearCTStereoDbleInv = pCS2->LinearCTStereoDbleInv; /*  enable inv. stereo */
        pCS->LinearCTStereoCarbInv = pCS2->LinearCTStereoCarbInv;
        pCS->nMaxLenLinearCTStereoDble = pCS2->nMaxLenLinearCTStereoDble;
        pCS->nMaxLenLinearCTStereoCarb = pCS2->nMaxLenLinearCTStereoCarb;

        pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTStereoDble;
        pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTStereoCarb;

        pCS->nLenLinearCTIsotopicStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
        pCS->nLenLinearCTIsotopicStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
    }
    */
    // END INCHI C FUNCTION: SetCtToNonIsotopicStereo

    pCS.LinearCTStereoDble = pCS2.LinearCTStereoDble;
    pCS.LinearCTStereoCarb = pCS2.LinearCTStereoCarb;
    pCS.LinearCTStereoDbleInv = pCS2.LinearCTStereoDbleInv;
    pCS.LinearCTStereoCarbInv = pCS2.LinearCTStereoCarbInv;
    pCS.nMaxLenLinearCTStereoDble = pCS2.nMaxLenLinearCTStereoDble;
    pCS.nMaxLenLinearCTStereoCarb = pCS2.nMaxLenLinearCTStereoCarb;
    pCS.nLenLinearCTStereoDble = pCS2.nLenLinearCTStereoDble;
    pCS.nLenLinearCTStereoCarb = pCS2.nLenLinearCTStereoCarb;
    pCS.nLenLinearCTIsotopicStereoDble = pCS2.nLenLinearCTIsotopicStereoDble;
    pCS.nLenLinearCTIsotopicStereoCarb = pCS2.nLenLinearCTIsotopicStereoCarb;
}

#[allow(non_snake_case)]
pub(crate) fn FillAllStereoDescriptors(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:766 FillAllStereoDescriptors
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int FillAllStereoDescriptors( CANON_GLOBALS *pCG,
                                  sp_ATOM       *at,
                                  int           num_atoms,
                                  const         AT_RANK *nCanonRank,
                                  const         AT_RANK *nAtomNumberCanon,
                                  CANON_STAT    *pCS )
    {
        int ret = 0, i;
        /*  initialize zero lengths */
        pCS->nLenLinearCTStereoCarb = 0;
        pCS->nLenLinearCTStereoDble = 0;

        /*  fill atom by atom */
        for (i = 0; !ret && i < num_atoms; i++)
        {
            ret = FillSingleStereoDescriptors( pCG, at, (int) nAtomNumberCanon[i], -1, nCanonRank
                              , pCS->LinearCTStereoCarb, &pCS->nLenLinearCTStereoCarb, pCS->nMaxLenLinearCTStereoCarb
                              , pCS->LinearCTStereoDble, &pCS->nLenLinearCTStereoDble, pCS->nMaxLenLinearCTStereoDble
                              , 0 /* bAllene */ );
        }
        for (i = 0; !ret && i < num_atoms; i++)
        {
            ret = FillSingleStereoDescriptors( pCG, at, (int) nAtomNumberCanon[i], -1, nCanonRank
                              , pCS->LinearCTStereoCarb, &pCS->nLenLinearCTStereoCarb, pCS->nMaxLenLinearCTStereoCarb
                              , pCS->LinearCTStereoDble, &pCS->nLenLinearCTStereoDble, pCS->nMaxLenLinearCTStereoDble
                              , 1 /* bAllene */ );
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: FillAllStereoDescriptors

    pCS.nLenLinearCTStereoCarb = 0;
    pCS.nLenLinearCTStereoDble = 0;
    let mut carb_len = 0_i32;
    let mut double_len = 0_i32;
    let mut ret = 0_i32;
    for allene in [0_i32, 1_i32] {
        let mut index = 0_i32;
        while ret == 0 && index < num_atoms {
            let atom_number = *heap
                .slice(nAtomNumberCanon.as_const())?
                .get(index as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let call_result = FillSingleStereoDescriptors(
                heap,
                pCG,
                at,
                i32::from(atom_number),
                -1,
                nCanonRank,
                pCS.LinearCTStereoCarb,
                &mut carb_len,
                pCS.nMaxLenLinearCTStereoCarb,
                pCS.LinearCTStereoDble,
                &mut double_len,
                pCS.nMaxLenLinearCTStereoDble,
                allene,
            );
            pCS.nLenLinearCTStereoCarb = carb_len;
            pCS.nLenLinearCTStereoDble = double_len;
            ret = call_result?;
            index = index.wrapping_add(1);
        }
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn SetKnownStereoBondParities(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nCanonRank: SourceConstPointer<AT_RANK>,
    nRank: SourceConstPointer<AT_RANK>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:801 SetKnownStereoBondParities
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int SetKnownStereoBondParities( CANON_GLOBALS *pCG,
                                    sp_ATOM       *at,
                                    int           num_atoms,
                                    const AT_RANK *nCanonRank,
                                    const AT_RANK *nRank,
                                    const AT_RANK *nAtomNumber )
    {
        int i, j, n, m, j1, k, num_neigh1, num_neigh2, iMax1, parity;
        int trans_i1, trans_i2, trans_k1, trans_k2, prev_trans, trans_k, num_set;
        int i1, i2, k1, k2, n1, n2, m1, m2, /*stereo_bond_parity,*/ cumulene_len;

        AT_RANK nAtomRank1, nAtomRank2, nAtom1NeighRank;
        AT_RANK nNeighRank1[MAX_NUM_STEREO_BONDS],
            nNeighRank2[MAX_NUM_STEREO_BONDS];
        AT_RANK nNeighCanonRank1[MAX_NUM_STEREO_BONDS],
            nNeighCanonRank2[MAX_NUM_STEREO_BONDS];

        for (i1 = 0, num_set = 0; i1 < num_atoms; i1++)
        {
            if (!at[i1].parity || !at[i1].stereo_bond_neighbor[0])
            {
                continue;
            }

            if (!PARITY_WELL_DEF( at[i1].parity ))
            {
                continue;
            }

            nAtomRank1 = nRank[i1];
            iMax1 = (int) nAtomRank1 - 1;
            num_neigh1 = at[i1].valence;

            for (n1 = 0;    n1 < MAX_NUM_STEREO_BONDS &&
                ( i2 = (int) at[i1].stereo_bond_neighbor[n1] );
                               n1++)
            {
                i2--;

                /*  found a stereo bond at[i1]-at[i2] adjacent to at[i1] */
                for (n2 = 0, m = 0;
                        n2 < MAX_NUM_STEREO_BONDS &&
                        ( m = (int) at[i2].stereo_bond_neighbor[n2] ) && m - 1 != i1;
                            n2++)
                    ; /* locate stereo bond (#n2) at the opposite atom at[i2] */

                if (m - 1 != i1 || at[i1].stereo_bond_parity[n1] != at[i2].stereo_bond_parity[n2])
                {
                    return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                }
                if (i1 < i2)
                {
                    continue; /* do not process same bond 2 times */
                }
                if (PARITY_KNOWN( at[i1].stereo_bond_parity[n1] ) ||
                    !PARITY_VAL( at[i1].stereo_bond_parity[n1] ))
                {
                    continue;
                }
                if (!PARITY_WELL_DEF( at[i1].parity ) ||
                     !PARITY_WELL_DEF( at[i2].parity ))
                {
                    continue;
                }
                if (PARITY_VAL( at[i1].stereo_bond_parity[n1] ) != AB_PARITY_CALC)
                {
                    continue;  /*  ?? program error ?? should not happen */ /*   <BRKPT> */
                }

                /*stereo_bond_parity = PARITY_VAL(at[i1].stereo_bond_parity[n1]);*/
                cumulene_len = BOND_CHAIN_LEN( at[i1].stereo_bond_parity[n1] );
                nAtomRank2 = nRank[i2];
                nAtom1NeighRank = nRank[(int) at[i1].neighbor[(int) at[i1].stereo_bond_ord[n1]]];
                num_neigh2 = at[i2].valence;

                /*  store ranks of at[i1] stereo bond neighbors except one connected by a stereo bond */
                k = (int) at[i1].stereo_bond_ord[n1];
                trans_i1 = 0;
                for (i = j = 0; i < num_neigh1; i++)
                {
                    if (i != k)
                    {
                        nNeighRank1[j] = nRank[(int) at[i1].neighbor[i]];
                        j++;
                    }
                }
                if (j == 2)
                {
                    if (nNeighRank1[0] == nNeighRank1[1])
                    {
                        /*  neighbors are constitutionally identical, can't find bond parity */
                        continue;
                    }
                    trans_i1 = insertions_sort( pCG, nNeighRank1, j, sizeof( nNeighRank1[0] ), comp_AT_RANK );
                }

                /*  store ranks of at[i2] stereo bond neighbors except one connected by a stereo bond */
                k = (int) at[i2].stereo_bond_ord[n2];
                trans_i2 = 0;
                for (i = j = 0; i < num_neigh2; i++)
                {
                    if (i != k)
                    {
                        nNeighRank2[j] = nRank[(int) at[i2].neighbor[i]];
                        j++;
                    }
                }

                if (j == 2)
                {
                    if (nNeighRank2[0] == nNeighRank2[1])
                    {
                        /*  neighbors are constitutionally identical, can't find bond parity */
                        continue;
                    }
                    trans_i2 = insertions_sort( pCG, nNeighRank2, j, sizeof( nNeighRank2[0] ), comp_AT_RANK );
                }

                prev_trans = -1;
                trans_k1 = -2; /* djb-rwth: ignoring LLVM warning: value used */
                trans_k = -4; /* 2004-04-28 */

                /*  find all pairs of atoms that can be mapped on at[i1], at[i2] pair */
                for (j1 = 0;
                        j1 <= iMax1 && nAtomRank1 == nRank[k1 = (int) nAtomNumber[iMax1 - j1]];
                            j1++)
                {
                    /*  at[k1] is constitutionally equivalent to at[i1] */
                    /*  find all at[k1] neighbors that have rank nAtomRank2; */
                    /*  then find at[k2] constitutionally equivalent at at[i2] */
                    if (at[k1].valence != num_neigh1)
                    {
                        return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                    }
                    for (m1 = 0; m1 < num_neigh1; m1++)
                    {
                        int prev, next, len;
                        if (nAtom1NeighRank != nRank[k2 = (int) at[k1].neighbor[m1]])
                        {
                            continue;
                        }
                        m2 = -1; /*  undefined yet */
                        prev = k1;
                        /* djb-rwth: removing redundant code */
                        if (cumulene_len)
                        {
                            for (len = 0, next = (int) at[k1].neighbor[m1]; len < cumulene_len; len++)
                            {
                                if (at[next].valence == 2 && !at[next].num_H)
                                {
                                    j = ( (int) at[next].neighbor[0] == prev );
                                    prev = next;
                                    next = at[next].neighbor[j];
                                }
                                else
                                {
                                    break; /*  cannot continue */
                                }
                            }
                            if (len != cumulene_len || nAtomRank2 != nRank[next])
                            {
                                continue;  /*  not found */
                            }
                            k2 = next;
                        }
                        if (at[k2].valence != num_neigh2)
                        {
                            return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                        }

                        /*  store canon. ranks of at[k1] neighbors */ /*  use i,j,k,m,n */
                        for (n = j = 0; n < num_neigh1; n++)
                        {
                            if (n != m1)
                            {
                                i = (int) at[k1].neighbor[n];
                                for (m = 0; m < num_neigh1 - 1; m++)
                                {
                                    if (nRank[i] == nNeighRank1[m])
                                    {
                                        nNeighCanonRank1[m] = nCanonRank[i];
                                        j++;
                                        break;
                                    }
                                }
                            }
                        }
                        if (j != num_neigh1 - 1)
                        {
                            return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                        }
                        if (j == 2)
                        {
                            trans_k1 = insertions_sort( pCG, nNeighCanonRank1, j, sizeof( nNeighCanonRank1[0] ), comp_AT_RANK );
                        }
                        else
                        {
                            trans_k1 = 0;
                        }

                        /*  store canon. ranks of at[k2] neighbors */ /*  use i,j,k,m,n */
                        for (n = j = 0; n < num_neigh2; n++)
                        {
                            i = (int) at[k2].neighbor[n];
                            if (i == prev)
                            {
                                /* neighbor belongs to the stereobond */
                                m2 = n;
                            }
                            else
                            {
                                for (m = 0; m < num_neigh2 - 1; m++)
                                {
                                    if (nRank[i] == nNeighRank2[m])
                                    {
                                        nNeighCanonRank2[m] = nCanonRank[i];
                                        j++;
                                        break;
                                    }
                                }
                            }
                        }
                        if (j != num_neigh2 - 1 || m2 < 0)
                        {
                            return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                        }
                        if (j == 2)
                        {
                            trans_k2 = insertions_sort( pCG, nNeighCanonRank2, j, sizeof( nNeighCanonRank2[0] ), comp_AT_RANK );
                        }
                        else
                        {
                            trans_k2 = 0;
                        }
                        trans_k = ( trans_k1 + trans_k2 ) % 2;
                        if (prev_trans < 0)
                        {
                            prev_trans = trans_k;
                        }
                        else if (prev_trans != trans_k)
                        {
                            /* was != trans_k1, changed 9-23-2003 */
                            break; /*  different number of transpositions */
                        }
                    } /* end of the second atom mapping cycle */
                    if (prev_trans >= 0 && prev_trans != trans_k)
                    { /* was != trans_k1, changed 9-23-2003 */
                        break;
                    }
                } /* end of the first atom mapping cycle */

                if (prev_trans == trans_k)
                {
                    /* was == trans_k1, changed 9-23-2003 */
                    int z_prod;

                    /*  all mappings of canonical numbers on the */
                    /*  stereo bond at[i1]-at[i2] produce equivalent numberings. */
                    /*  Therefore the stereo bond parity is known at this time. */
                    /*  parity_1 = at[i1].parity + (trans_i1 + trans_k1 + num_neigh1 - 1) + (int)at[i1].stereo_bond_ord[n1] */
                    /*  expression in parentheses is equivalent to rank[first neigh] > rank[second neigh] */
                    /*  same for parity_2. */
                    /*  parity_2 = at[i2].parity + (trans_i2 + trans_k2 + num_neigh2 - 1) + (int)at[i2].stereo_bond_ord[n2] */
                    /*  Sum of the two parities (without stereo_bond_z_prod) is: */

                    parity = ( at[i1].parity + at[i2].parity + prev_trans + trans_i1 + trans_i2
                                  + num_neigh1 + num_neigh2
                                  + (int) at[i1].stereo_bond_ord[n1] + (int) at[i2].stereo_bond_ord[n2] ) % 2;

                    z_prod = at[i1].stereo_bond_z_prod[n1];
                    if (MIN_DOT_PROD > abs( z_prod ))
                    {
                        parity = AB_PARITY_UNDF; /*  undefined because of geometry */
                    }
                    else
                    {
                        parity = ( z_prod > 0 ) ? 2 - parity : 1 + parity;
                    }
                    at[i1].stereo_bond_parity[n1] = ALL_BUT_PARITY( at[i1].stereo_bond_parity[n1] ) | parity;
                    at[i2].stereo_bond_parity[n2] = ALL_BUT_PARITY( at[i2].stereo_bond_parity[n2] ) | parity;
                    num_set++;
                }
            }
        }

        return num_set;
    }
    */
    // END INCHI C FUNCTION: SetKnownStereoBondParities

    fn sort_ranks(values: &mut [AT_RANK], count: usize) -> Result<i32, SourceHeapError> {
        let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(values);
        insertions_sort(
            bytes,
            count,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                Ok(comp_AT_RANK(
                    AT_RANK::from_ne_bytes([left[0], left[1]]),
                    AT_RANK::from_ne_bytes([right[0], right[1]]),
                ))
            },
        )
    }

    let _ = pCG;
    let mut atoms = heap.slice(at.as_const())?.to_vec();
    let canonical_ranks = heap.slice(nCanonRank)?.to_vec();
    let ranks = heap.slice(nRank)?.to_vec();
    let atom_numbers = heap.slice(nAtomNumber)?.to_vec();
    macro_rules! finish {
        ($value:expr) => {{
            heap.slice_mut(at)?
                .get_mut(..atoms.len())
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&atoms);
            return Ok($value);
        }};
    }

    let mut num_set = 0_i32;
    for i1_raw in 0..num_atoms {
        let i1 = i1_raw as usize;
        if atoms[i1].parity == 0
            || atoms[i1].stereo_bond_neighbor[0] == 0
            || !parity_well_defined(i32::from(atoms[i1].parity))
        {
            continue;
        }
        let atom_rank1 = ranks[i1];
        let maximum1 = i32::from(atom_rank1) - 1;
        let neighbor_count1 = i32::from(atoms[i1].valence);
        let mut n1 = 0_usize;
        while n1 < MAX_NUM_STEREO_BONDS as usize {
            let opposite_number = atoms[i1].stereo_bond_neighbor[n1];
            if opposite_number == 0 {
                break;
            }
            let i2 = usize::from(opposite_number - 1);
            let mut n2 = 0_usize;
            let mut reverse_number = 0_i32;
            while n2 < MAX_NUM_STEREO_BONDS as usize {
                reverse_number = i32::from(atoms[i2].stereo_bond_neighbor[n2]);
                if reverse_number == 0 || reverse_number - 1 == i1_raw {
                    break;
                }
                n2 += 1;
            }
            if reverse_number - 1 != i1_raw
                || n2 >= MAX_NUM_STEREO_BONDS as usize
                || atoms[i1].stereo_bond_parity[n1] != atoms[i2].stereo_bond_parity[n2]
            {
                finish!(CT_STEREOCOUNT_ERR);
            }
            if i1 < i2 {
                n1 += 1;
                continue;
            }
            let encoded_parity = i32::from(atoms[i1].stereo_bond_parity[n1]);
            let parity_value = encoded_parity & BITS_PARITY;
            if parity_known(encoded_parity) || parity_value == 0 {
                n1 += 1;
                continue;
            }
            if !parity_well_defined(i32::from(atoms[i1].parity))
                || !parity_well_defined(i32::from(atoms[i2].parity))
                || parity_value != AB_PARITY_CALC as i32
            {
                n1 += 1;
                continue;
            }

            let cumulene_len = bond_chain_len(encoded_parity);
            let atom_rank2 = ranks[i2];
            let order1 = usize::try_from(i32::from(atoms[i1].stereo_bond_ord[n1]))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom1_neighbor_rank = ranks[usize::from(atoms[i1].neighbor[order1])];
            let neighbor_count2 = i32::from(atoms[i2].valence);

            let mut neighbor_ranks1 = [0_u16; MAX_NUM_STEREO_BONDS as usize];
            let mut count1 = 0_usize;
            for neighbor_order in 0..neighbor_count1 {
                if neighbor_order as usize != order1 {
                    if count1 >= neighbor_ranks1.len() {
                        return Err(SourceHeapError::PointerOutOfBounds);
                    }
                    neighbor_ranks1[count1] =
                        ranks[usize::from(atoms[i1].neighbor[neighbor_order as usize])];
                    count1 += 1;
                }
            }
            let mut trans_i1 = 0_i32;
            if count1 == 2 {
                if neighbor_ranks1[0] == neighbor_ranks1[1] {
                    n1 += 1;
                    continue;
                }
                trans_i1 = sort_ranks(&mut neighbor_ranks1, count1)?;
            }

            let order2 = usize::try_from(i32::from(atoms[i2].stereo_bond_ord[n2]))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut neighbor_ranks2 = [0_u16; MAX_NUM_STEREO_BONDS as usize];
            let mut count2 = 0_usize;
            for neighbor_order in 0..neighbor_count2 {
                if neighbor_order as usize != order2 {
                    if count2 >= neighbor_ranks2.len() {
                        return Err(SourceHeapError::PointerOutOfBounds);
                    }
                    neighbor_ranks2[count2] =
                        ranks[usize::from(atoms[i2].neighbor[neighbor_order as usize])];
                    count2 += 1;
                }
            }
            let mut trans_i2 = 0_i32;
            if count2 == 2 {
                if neighbor_ranks2[0] == neighbor_ranks2[1] {
                    n1 += 1;
                    continue;
                }
                trans_i2 = sort_ranks(&mut neighbor_ranks2, count2)?;
            }

            let mut previous_trans = -1_i32;
            let mut trans_k = -4_i32;
            let mut j1 = 0_i32;
            while j1 <= maximum1 {
                let k1 = usize::from(atom_numbers[(maximum1 - j1) as usize]);
                if atom_rank1 != ranks[k1] {
                    break;
                }
                if i32::from(atoms[k1].valence) != neighbor_count1 {
                    finish!(CT_STEREOCOUNT_ERR);
                }
                for m1 in 0..neighbor_count1 {
                    let mut k2 = usize::from(atoms[k1].neighbor[m1 as usize]);
                    if atom1_neighbor_rank != ranks[k2] {
                        continue;
                    }
                    let mut previous = k1;
                    if cumulene_len != 0 {
                        let mut length = 0_i32;
                        let mut next = k2;
                        while length < cumulene_len {
                            if atoms[next].valence == 2 && atoms[next].num_H == 0 {
                                let branch =
                                    usize::from(usize::from(atoms[next].neighbor[0]) == previous);
                                previous = next;
                                next = usize::from(atoms[next].neighbor[branch]);
                                length += 1;
                            } else {
                                break;
                            }
                        }
                        if length != cumulene_len || atom_rank2 != ranks[next] {
                            continue;
                        }
                        k2 = next;
                    }
                    if i32::from(atoms[k2].valence) != neighbor_count2 {
                        finish!(CT_STEREOCOUNT_ERR);
                    }

                    let mut canonical_neighbor_ranks1 = [0_u16; MAX_NUM_STEREO_BONDS as usize];
                    let mut matched1 = 0_i32;
                    for neighbor_order in 0..neighbor_count1 {
                        if neighbor_order == m1 {
                            continue;
                        }
                        let neighbor = usize::from(atoms[k1].neighbor[neighbor_order as usize]);
                        for rank_slot in 0..neighbor_count1 - 1 {
                            if ranks[neighbor] == neighbor_ranks1[rank_slot as usize] {
                                canonical_neighbor_ranks1[rank_slot as usize] =
                                    canonical_ranks[neighbor];
                                matched1 += 1;
                                break;
                            }
                        }
                    }
                    if matched1 != neighbor_count1 - 1 {
                        finish!(CT_STEREOCOUNT_ERR);
                    }
                    let trans_k1 = if matched1 == 2 {
                        sort_ranks(&mut canonical_neighbor_ranks1, matched1 as usize)?
                    } else {
                        0
                    };

                    let mut canonical_neighbor_ranks2 = [0_u16; MAX_NUM_STEREO_BONDS as usize];
                    let mut matched2 = 0_i32;
                    let mut m2 = -1_i32;
                    for neighbor_order in 0..neighbor_count2 {
                        let neighbor = usize::from(atoms[k2].neighbor[neighbor_order as usize]);
                        if neighbor == previous {
                            m2 = neighbor_order;
                        } else {
                            for rank_slot in 0..neighbor_count2 - 1 {
                                if ranks[neighbor] == neighbor_ranks2[rank_slot as usize] {
                                    canonical_neighbor_ranks2[rank_slot as usize] =
                                        canonical_ranks[neighbor];
                                    matched2 += 1;
                                    break;
                                }
                            }
                        }
                    }
                    if matched2 != neighbor_count2 - 1 || m2 < 0 {
                        finish!(CT_STEREOCOUNT_ERR);
                    }
                    let trans_k2 = if matched2 == 2 {
                        sort_ranks(&mut canonical_neighbor_ranks2, matched2 as usize)?
                    } else {
                        0
                    };
                    trans_k = trans_k1.wrapping_add(trans_k2) % 2;
                    if previous_trans < 0 {
                        previous_trans = trans_k;
                    } else if previous_trans != trans_k {
                        break;
                    }
                }
                if previous_trans >= 0 && previous_trans != trans_k {
                    break;
                }
                j1 += 1;
            }

            if previous_trans == trans_k {
                let mut parity = (i32::from(atoms[i1].parity)
                    .wrapping_add(i32::from(atoms[i2].parity))
                    .wrapping_add(previous_trans)
                    .wrapping_add(trans_i1)
                    .wrapping_add(trans_i2)
                    .wrapping_add(neighbor_count1)
                    .wrapping_add(neighbor_count2)
                    .wrapping_add(i32::from(atoms[i1].stereo_bond_ord[n1]))
                    .wrapping_add(i32::from(atoms[i2].stereo_bond_ord[n2])))
                    % 2;
                let z_product = i32::from(atoms[i1].stereo_bond_z_prod[n1]);
                if z_product.abs() < MIN_DOT_PROD as i32 {
                    parity = AB_PARITY_UNDF as i32;
                } else {
                    parity = if z_product > 0 {
                        2_i32.wrapping_sub(parity)
                    } else {
                        1_i32.wrapping_add(parity)
                    };
                }
                atoms[i1].stereo_bond_parity[n1] = ((encoded_parity & !BITS_PARITY) | parity) as i8;
                let opposite_encoded = i32::from(atoms[i2].stereo_bond_parity[n2]);
                atoms[i2].stereo_bond_parity[n2] =
                    ((opposite_encoded & !BITS_PARITY) | parity) as i8;
                num_set = num_set.wrapping_add(1);
            }
            n1 += 1;
        }
    }
    finish!(num_set)
}

#[allow(non_snake_case)]
pub(crate) fn MarkKnownEqualStereoBondParities(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nRank: SourceConstPointer<AT_RANK>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1093 MarkKnownEqualStereoBondParities
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int MarkKnownEqualStereoBondParities( sp_ATOM       *at,
                                          int           num_atoms,
                                          const AT_RANK *nRank,
                                          const AT_RANK *nAtomNumber )
    {
        int j, n, m, j1, num_neigh1, num_neigh2, iMax1;
        int num_set, /*num_sb1, num_sb2,*/ bDifferentParities;
        int i1, i2, k1, k2, n1, n2, m1, m2, s1, s2, stereo_bond_parity, stereo_bond_parity2, cumulene_len;
        AT_RANK nAtomRank1, nAtomRank2, nAtom1NeighRank, nAtom2NeighRank;

        /* djb-rwth: removing redundant code */

        for (i1 = 0, num_set = 0; i1 < num_atoms; i1++)
        {
            if (!at[i1].parity || !at[i1].stereo_bond_neighbor[0])
            {
                continue;
            }

            nAtomRank1 = nRank[i1];
            iMax1 = (int) nAtomRank1 - 1;
            num_neigh1 = at[i1].valence;

            /*  count stereogenic bonds adjacent to at[i1] */
            for (n1 = 0;
                    n1 < MAX_NUM_STEREO_BONDS && at[i1].stereo_bond_neighbor[n1];
                        n1++);


            /*num_sb1 = n1;*/
            /*  search for bonds possibly constitutionally equivalent to each of the adjacent bonds */
            /*  and find if all of them have same already known parity */

            for (n1 = 0;
                    n1 < MAX_NUM_STEREO_BONDS && ( i2 = (int) at[i1].stereo_bond_neighbor[n1] );
                        n1++) /* djb-rwth: removing redundant code */
            {
                i2--;

                nAtomRank2 = nRank[i2];
                if (nAtomRank2 < nAtomRank1 || (nAtomRank2 == nAtomRank1 && i1 < i2)) /* djb-rwth: addressing LLVM warning */
                {
                    /*  An attempt to reduce unnecessary repetitions. */
                    /*  We still have repetitions because we do not accumulate a list of */
                    /*  processed (nAtomRank2, nAtomRank1) pairs. */
                    continue;
                }

                bDifferentParities = -1;   /*  parities have not been compared yet */

                /*  found a stereo bond at[i1]-at[i2] (adjacent to at[i1]) */
                /*
                if ( !PARITY_KNOWN(at[i1].stereo_bond_parity[n1]) || (at[i1].stereo_bond_parity[n1] & KNOWN_PARITIES_EQL) )
                {
                    continue;
                }
                */
                if (at[i1].stereo_bond_parity[n1] & KNOWN_PARITIES_EQL)
                {
                    continue;
                }


                /*  stereo bond has known or unknown parity; we have not checked it yet */

                for (n2 = 0;
                     n2 < MAX_NUM_STEREO_BONDS && at[i2].stereo_bond_neighbor[n2];
                     n2++)
                {
                    ;
                }

                /*num_sb2 = n2;*/
                for (n2 = 0, m = 0;
                     n2 < MAX_NUM_STEREO_BONDS &&
                     ( m = (int) at[i2].stereo_bond_neighbor[n2] ) &&
                        m - 1 != i1;
                     n2++)
                {
                    ;
                }

                if (m - 1 != i1 || at[i1].stereo_bond_parity[n1] != at[i2].stereo_bond_parity[n2])
                {
                    return CT_STEREOCOUNT_ERR; /*  program error: stereo bonds data in two directions are different */ /*   <BRKPT> */
                }

                stereo_bond_parity = PARITY_VAL( at[i1].stereo_bond_parity[n1] );
                cumulene_len = BOND_CHAIN_LEN( at[i1].stereo_bond_parity[n1] );
                nAtom1NeighRank = nRank[(int) at[i1].neighbor[(int) at[i1].stereo_bond_ord[n1]]];
                nAtom2NeighRank = nRank[(int) at[i2].neighbor[(int) at[i2].stereo_bond_ord[n2]]];
                num_neigh2 = at[i2].valence;

                /*  find all pairs of atoms that possibly can be mapped on at[i1], at[i2] pair */
                /*  (we may also find pairs that cannot be mapped, but we cannot miss any pair */
                /*  that can be mapped) */

                for (j1 = 0; j1 <= iMax1 && nAtomRank1 == nRank[k1 = (int) nAtomNumber[iMax1 - j1]]; j1++)
                {
                    /*  at[k1] is constitutionally equivalent to at[i1] */
                    /*  find all at[k1] stereo bond neighbors at[k2] that have rank nAtomRank2; */
                    /*  then find at[k2] constitutionally equivalent at at[i2] */

                    if (at[k1].valence != num_neigh1)
                    {
                        return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                    }

                    if (!at[k1].bHasStereoOrEquToStereo)
                    {
                        at[k1].bHasStereoOrEquToStereo = 1;
                    }

                    /* -- do not check number of stereo bonds, check bonds themselves --
                    for ( s1 = 0; s1 < MAX_NUM_STEREO_BONDS && at[k1].stereo_bond_neighbor[s1]; s1++ )
                    {
                        ;
                    }
                    if ( num_sb1 != s1 )
                    {
                        bDifferentParities = 1;
                    }
                    */

                    for (m1 = 0; m1 < num_neigh1; m1++)
                    {
                        /*  Looking for at[k1] neighbor with nRank=nAtom1NeighRank. */
                        /*  This neighbor may be on the bond constit. equivalent to at[i1]-at[i2] stereo bond */
                        /*  (or may be constit. equivalent an adjacent to at[i1] atom in a stereogenic cumulene chain) */
                        int prev, next, len;
                        if (nAtom1NeighRank != nRank[k2 = (int) at[k1].neighbor[m1]])
                            continue;

                        /*  found at[k1] neighbor with nRank=nAtom1NeighRank */

                        m2 = -1; /*  undefined yet */
                        prev = k1;
                        /* djb-rwth: removing redundant code */

                        /*  if cumulene then bypass the cumulene chain */

                        if (cumulene_len)
                        {

                            for (len = 0, next = (int) at[k1].neighbor[m1]; len < cumulene_len; len++)
                            {
                                if (at[next].valence == 2 && !at[next].num_H)
                                {
                                    j = ( (int) at[next].neighbor[0] == prev );
                                    prev = next;
                                    next = at[next].neighbor[j];
                                }
                                else
                                {
                                    break; /*  cannot continue: end of cumulene chain */
                                }
                            }

                            if (len != cumulene_len || nAtomRank2 != nRank[next])
                            {
                                continue;  /*  cumulene chain not found at this neighbor */
                            }

                            if (nAtom2NeighRank != nRank[prev])
                            {
                                /* continue; */ /*  ??? program error ??? If not, must be a very rare event */
                                return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                            }

                            k2 = next;
                        }

                        /*  a connected pair of constit. equivalent atoms found */

                        if (at[k2].valence != num_neigh2)
                        {
                            return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                        }

                        for (n = 0; n < num_neigh2; n++)
                        {
                            if (prev == (int) at[k2].neighbor[n])
                            {
                                m2 = n; /*  found bond from the opposite end of a possibly stereogenic bond */
                                break;
                            }
                        }

                        if (m2 < 0)
                        {
                            return CT_STEREOCOUNT_ERR; /*  program error: opposite direction bond not found */ /*   <BRKPT> */
                        }

                        if (!at[k2].bHasStereoOrEquToStereo)
                        {
                            at[k2].bHasStereoOrEquToStereo = 1;
                        }


                        /*  check if atoms at[k1] and at[k2] are connected by a stereo bond */
                        for (s1 = 0, m = 0;
                             s1 < MAX_NUM_STEREO_BONDS &&
                             ( m = (int) at[k1].stereo_bond_neighbor[s1] ) &&
                                        m - 1 != k2;
                             s1++)
                        {
                            ;
                        }
                        if (m - 1 != k2)
                        {
                            bDifferentParities = 1; /*  cannot find the stereo bond */
                            at[k1].bHasStereoOrEquToStereo =
                                at[k2].bHasStereoOrEquToStereo = 2;
                            continue;
                        }

                        /*  -- do not check number of stereo bonds, check bonds themselves --
                        for ( s2 = 0; s2 < MAX_NUM_STEREO_BONDS && at[k2].stereo_bond_neighbor[s2]; s2++ )
                         {
                            ;
                         }
                        if ( num_sb2 != s2 )
                        {
                            bDifferentParities = 1;
                            continue;
                        }
                        */

                        for (s2 = 0, m = 0;
                             s2 < MAX_NUM_STEREO_BONDS &&
                             ( m = (int) at[k2].stereo_bond_neighbor[s2] ) &&
                                        m - 1 != k1;
                             s2++)
                        {
                            ;
                        }

                        if (m - 1 != k1)
                        {
                            /*
                            bDifferentParities = 1; // cannot find the stereo bond
                            continue;
                            */
                            return CT_STEREOCOUNT_ERR; /*  program error: opposite direction bond not found */ /*   <BRKPT> */
                        }

                        if (at[k1].stereo_bond_parity[s1] != at[k2].stereo_bond_parity[s2])
                        {
                            bDifferentParities = 1;
                            continue;
                        }
                        stereo_bond_parity2 = PARITY_VAL( at[k1].stereo_bond_parity[s1] );
                        if (stereo_bond_parity2 != stereo_bond_parity)
                        {
                            bDifferentParities = 1;
                            continue;
                        }
                        if (stereo_bond_parity2 == stereo_bond_parity && bDifferentParities < 0)
                        {
                            bDifferentParities = 0;
                        }
                    }
                }

                /*  mark equal parities */
                if (0 == bDifferentParities && PARITY_KNOWN( stereo_bond_parity ))
                {
                    for (j1 = 0;
                            j1 <= iMax1 && nAtomRank1 == nRank[k1 = (int) nAtomNumber[iMax1 - j1]];
                                j1++)
                    {
                        /*  at[k1] is constitutionally equivalent to at[i1] */
                        for (s1 = 0; s1 < MAX_NUM_STEREO_BONDS && ( k2 = (int) at[k1].stereo_bond_neighbor[s1] ); s1++) /* djb-rwth: removing redundant code */
                        {
                            k2--;
                            if (nRank[k2] == nAtomRank2)
                            {
                                int b1, b2;
                                for (s2 = 0, m = 0;
                                     s2 < MAX_NUM_STEREO_BONDS &&
                                        ( m = (int) at[k2].stereo_bond_neighbor[s2] ) &&
                                        m - 1 != k1;
                                     s2++)
                                {
                                    ;
                                }

                                if (m - 1 != k1)
                                {
                                    return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                                }
                                /*  mark the stereo bonds */
                                b1 = !( at[k1].stereo_bond_parity[s1] & KNOWN_PARITIES_EQL );
                                b2 = !( at[k2].stereo_bond_parity[s2] & KNOWN_PARITIES_EQL );
                                if (2 == b1 + b2)
                                {
                                    at[k1].stereo_bond_parity[s1] |= KNOWN_PARITIES_EQL;
                                    at[k2].stereo_bond_parity[s2] |= KNOWN_PARITIES_EQL;
                                    num_set++;
                                }
                                else if (b1 || b2)
                                {
                                    return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                                }
                            }
                        }
                    }
                }
            }
        }

        return num_set;
    }
    */
    // END INCHI C FUNCTION: MarkKnownEqualStereoBondParities

    if num_atoms <= 0 {
        return Ok(0);
    }
    let mut atoms = heap.slice(at.as_const())?.to_vec();
    let ranks = heap.slice(nRank)?.to_vec();
    let atom_numbers = heap.slice(nAtomNumber)?.to_vec();
    macro_rules! finish {
        ($value:expr) => {{
            let destination = heap.slice_mut(at)?;
            destination
                .get_mut(..atoms.len())
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&atoms);
            return Ok($value);
        }};
    }

    let mut num_set = 0_i32;
    for i1_raw in 0..num_atoms {
        let i1 = i1_raw as usize;
        if atoms[i1].parity == 0 || atoms[i1].stereo_bond_neighbor[0] == 0 {
            continue;
        }
        let atom_rank1 = ranks[i1];
        let i_max1 = i32::from(atom_rank1) - 1;
        let num_neigh1 = i32::from(atoms[i1].valence);

        let mut n1 = 0_usize;
        while n1 < MAX_NUM_STEREO_BONDS as usize && atoms[i1].stereo_bond_neighbor[n1] != 0 {
            n1 += 1;
        }
        n1 = 0;
        while n1 < MAX_NUM_STEREO_BONDS as usize {
            let neighbor_number = atoms[i1].stereo_bond_neighbor[n1];
            if neighbor_number == 0 {
                break;
            }
            let i2 = usize::from(neighbor_number - 1);
            let atom_rank2 = ranks[i2];
            if atom_rank2 < atom_rank1 || (atom_rank2 == atom_rank1 && i1 < i2) {
                n1 += 1;
                continue;
            }
            let mut different_parities = -1_i32;
            if i32::from(atoms[i1].stereo_bond_parity[n1]) & 0x40 != 0 {
                n1 += 1;
                continue;
            }

            let mut n2 = 0_usize;
            while n2 < MAX_NUM_STEREO_BONDS as usize && atoms[i2].stereo_bond_neighbor[n2] != 0 {
                n2 += 1;
            }
            n2 = 0;
            let mut m = 0_i32;
            while n2 < MAX_NUM_STEREO_BONDS as usize {
                m = i32::from(atoms[i2].stereo_bond_neighbor[n2]);
                if m == 0 || m - 1 == i1_raw {
                    break;
                }
                n2 += 1;
            }
            if m - 1 != i1_raw
                || n2 >= MAX_NUM_STEREO_BONDS as usize
                || atoms[i1].stereo_bond_parity[n1] != atoms[i2].stereo_bond_parity[n2]
            {
                finish!(CT_STEREOCOUNT_ERR);
            }

            let stereo_bond_parity = i32::from(atoms[i1].stereo_bond_parity[n1]) & BITS_PARITY;
            let cumulene_len = bond_chain_len(i32::from(atoms[i1].stereo_bond_parity[n1]));
            let order1 = atoms[i1].stereo_bond_ord[n1] as usize;
            let order2 = atoms[i2].stereo_bond_ord[n2] as usize;
            let atom1_neighbor_rank = ranks[usize::from(atoms[i1].neighbor[order1])];
            let atom2_neighbor_rank = ranks[usize::from(atoms[i2].neighbor[order2])];
            let num_neigh2 = i32::from(atoms[i2].valence);

            let mut j1 = 0_i32;
            while j1 <= i_max1 {
                let k1 = usize::from(atom_numbers[(i_max1 - j1) as usize]);
                if atom_rank1 != ranks[k1] {
                    break;
                }
                if i32::from(atoms[k1].valence) != num_neigh1 {
                    finish!(CT_STEREOCOUNT_ERR);
                }
                if atoms[k1].bHasStereoOrEquToStereo == 0 {
                    atoms[k1].bHasStereoOrEquToStereo = 1;
                }
                for m1_raw in 0..num_neigh1 {
                    let m1 = m1_raw as usize;
                    let mut k2 = usize::from(atoms[k1].neighbor[m1]);
                    if atom1_neighbor_rank != ranks[k2] {
                        continue;
                    }
                    let mut m2 = -1_i32;
                    let mut prev = k1;
                    if cumulene_len != 0 {
                        let mut len = 0_i32;
                        let mut next = k2;
                        while len < cumulene_len {
                            if atoms[next].valence == 2 && atoms[next].num_H == 0 {
                                let branch =
                                    usize::from(usize::from(atoms[next].neighbor[0]) == prev);
                                prev = next;
                                next = usize::from(atoms[next].neighbor[branch]);
                                len += 1;
                            } else {
                                break;
                            }
                        }
                        if len != cumulene_len || atom_rank2 != ranks[next] {
                            continue;
                        }
                        if atom2_neighbor_rank != ranks[prev] {
                            finish!(CT_STEREOCOUNT_ERR);
                        }
                        k2 = next;
                    }
                    if i32::from(atoms[k2].valence) != num_neigh2 {
                        finish!(CT_STEREOCOUNT_ERR);
                    }
                    for n_raw in 0..num_neigh2 {
                        if prev == usize::from(atoms[k2].neighbor[n_raw as usize]) {
                            m2 = n_raw;
                            break;
                        }
                    }
                    if m2 < 0 {
                        finish!(CT_STEREOCOUNT_ERR);
                    }
                    if atoms[k2].bHasStereoOrEquToStereo == 0 {
                        atoms[k2].bHasStereoOrEquToStereo = 1;
                    }

                    let mut s1 = 0_usize;
                    m = 0;
                    while s1 < MAX_NUM_STEREO_BONDS as usize {
                        m = i32::from(atoms[k1].stereo_bond_neighbor[s1]);
                        if m == 0 || m - 1 == k2 as i32 {
                            break;
                        }
                        s1 += 1;
                    }
                    if m - 1 != k2 as i32 {
                        different_parities = 1;
                        atoms[k1].bHasStereoOrEquToStereo = 2;
                        atoms[k2].bHasStereoOrEquToStereo = 2;
                        continue;
                    }
                    let mut s2 = 0_usize;
                    m = 0;
                    while s2 < MAX_NUM_STEREO_BONDS as usize {
                        m = i32::from(atoms[k2].stereo_bond_neighbor[s2]);
                        if m == 0 || m - 1 == k1 as i32 {
                            break;
                        }
                        s2 += 1;
                    }
                    if m - 1 != k1 as i32 || s2 >= MAX_NUM_STEREO_BONDS as usize {
                        finish!(CT_STEREOCOUNT_ERR);
                    }
                    if atoms[k1].stereo_bond_parity[s1] != atoms[k2].stereo_bond_parity[s2] {
                        different_parities = 1;
                        continue;
                    }
                    let stereo_bond_parity2 =
                        i32::from(atoms[k1].stereo_bond_parity[s1]) & BITS_PARITY;
                    if stereo_bond_parity2 != stereo_bond_parity {
                        different_parities = 1;
                        continue;
                    }
                    if different_parities < 0 {
                        different_parities = 0;
                    }
                }
                j1 += 1;
            }

            if different_parities == 0 && (1..=4).contains(&stereo_bond_parity) {
                let mut j1 = 0_i32;
                while j1 <= i_max1 {
                    let k1 = usize::from(atom_numbers[(i_max1 - j1) as usize]);
                    if atom_rank1 != ranks[k1] {
                        break;
                    }
                    let mut s1 = 0_usize;
                    while s1 < MAX_NUM_STEREO_BONDS as usize {
                        let neighbor = atoms[k1].stereo_bond_neighbor[s1];
                        if neighbor == 0 {
                            break;
                        }
                        let k2 = usize::from(neighbor - 1);
                        if ranks[k2] == atom_rank2 {
                            let mut s2 = 0_usize;
                            m = 0;
                            while s2 < MAX_NUM_STEREO_BONDS as usize {
                                m = i32::from(atoms[k2].stereo_bond_neighbor[s2]);
                                if m == 0 || m - 1 == k1 as i32 {
                                    break;
                                }
                                s2 += 1;
                            }
                            if m - 1 != k1 as i32 || s2 >= MAX_NUM_STEREO_BONDS as usize {
                                finish!(CT_STEREOCOUNT_ERR);
                            }
                            let b1 =
                                i32::from(i32::from(atoms[k1].stereo_bond_parity[s1]) & 0x40 == 0);
                            let b2 =
                                i32::from(i32::from(atoms[k2].stereo_bond_parity[s2]) & 0x40 == 0);
                            if b1 + b2 == 2 {
                                atoms[k1].stereo_bond_parity[s1] |= 0x40;
                                atoms[k2].stereo_bond_parity[s2] |= 0x40;
                                num_set = num_set.wrapping_add(1);
                            } else if b1 != 0 || b2 != 0 {
                                finish!(CT_STEREOCOUNT_ERR);
                            }
                        }
                        s1 += 1;
                    }
                    j1 += 1;
                }
            }
            n1 += 1;
        }
    }

    finish!(num_set)
}

#[allow(non_snake_case)]
pub(crate) fn GetNextNeighborAndRank(
    heap: &SourceHeap,
    at: SourceConstPointer<sp_ATOM>,
    cur: AT_RANK,
    prev: AT_RANK,
    n: &mut AT_RANK,
    cr: &mut AT_RANK,
    nCanonRank: SourceConstPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1419 GetNextNeighborAndRank
    // INCHI✔️✔️: int GetNextNeighborAndRank( sp_ATOM       *at,
    // INCHI✔️✔️:                             AT_RANK       cur,
    // INCHI✔️✔️:                             AT_RANK       prev,
    // INCHI✔️✔️:                             AT_RANK       *n,
    // INCHI✔️✔️:                             AT_RANK       *cr,
    // INCHI✔️✔️:                             const AT_RANK *nCanonRank )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, val;
    // INCHI✔️✔️:     AT_RANK cr1 = MAX_ATOMS + 1, j, j1 = MAX_ATOMS + 1, crj;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0, val = at[(int) cur].valence; i < val; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (( j = at[cur].neighbor[i] ) != prev &&
    // INCHI✔️✔️:              cr1 > ( crj = nCanonRank[(int) j] ) && crj > *cr)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             cr1 = crj;
    // INCHI✔️✔️:             j1 = j;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (cr1 <= MAX_ATOMS)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *cr = cr1;
    // INCHI✔️✔️:         *n = (AT_RANK) j1;
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;  /*  program error */ /*   <BRKPT> */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetNextNeighborAndRank

    let atoms = heap.slice(at)?;
    let ranks = heap.slice(nCanonRank)?;
    let atom = atoms
        .get(usize::from(cur))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut cr1 = (MAX_ATOMS + 1) as AT_RANK;
    let mut j1 = (MAX_ATOMS + 1) as AT_RANK;
    for i in 0..i32::from(atom.valence) {
        let j = atom.neighbor[i as usize];
        let crj = *ranks
            .get(usize::from(j))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if j != prev && cr1 > crj && crj > *cr {
            cr1 = crj;
            j1 = j;
        }
    }
    if u32::from(cr1) <= MAX_ATOMS {
        *cr = cr1;
        *n = j1;
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetAndCheckNextNeighbors(
    heap: &SourceHeap,
    at: SourceConstPointer<sp_ATOM>,
    cur1: AT_RANK,
    prev1: AT_RANK,
    cur2: AT_RANK,
    prev2: AT_RANK,
    n1: &mut AT_RANK,
    n2: &mut AT_RANK,
    nVisited1: SourceConstPointer<AT_RANK>,
    nVisited2: SourceConstPointer<AT_RANK>,
    nRank: SourceConstPointer<AT_RANK>,
    nCanonRank: SourceConstPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1456 GetAndCheckNextNeighbors
    // INCHI✔️✔️: int GetAndCheckNextNeighbors( sp_ATOM       *at,
    // INCHI✔️✔️:                               AT_RANK       cur1,
    // INCHI✔️✔️:                               AT_RANK       prev1,
    // INCHI✔️✔️:                               AT_RANK       cur2,
    // INCHI✔️✔️:                               AT_RANK       prev2,
    // INCHI✔️✔️:                               AT_RANK       *n1,
    // INCHI✔️✔️:                               AT_RANK       *n2,
    // INCHI✔️✔️:                               AT_RANK       *nVisited1,
    // INCHI✔️✔️:                               AT_RANK       *nVisited2,
    // INCHI✔️✔️:                               const AT_RANK *nRank,
    // INCHI✔️✔️:                               const AT_RANK *nCanonRank )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_RANK cr1, cr2, s1, s2;
    // INCHI✔️✔️:     int     i1, i2, k1, k2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     cr1 = ( *n1 > MAX_ATOMS ) ? 0 : nCanonRank[(int) *n1];
    // INCHI✔️✔️:     cr2 = ( *n2 > MAX_ATOMS ) ? 0 : nCanonRank[(int) *n2];
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!GetNextNeighborAndRank( at, cur1, prev1, n1, &cr1, nCanonRank ) ||
    // INCHI✔️✔️:          !GetNextNeighborAndRank( at, cur2, prev2, n2, &cr2, nCanonRank ) ||
    // INCHI✔️✔️:          nRank[(int) *n1] != nRank[(int) *n2] || nVisited1[(int) *n1] != nVisited2[(int) *n2])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;  /*  program error; no breakpoint here */ /*   <BRKPT> */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  Even though the bond or cumulene might have already been checked, check it: this is */
    // INCHI✔️✔️:     /*  the only place we can check stereo bonds and cumulenes that are not edges of the DFS tree */
    // INCHI✔️✔️:     /*  The code works both for a stereo bond and a stereogenic cumulene. */
    // INCHI✔️✔️:     for (i1 = 0, k1 = 0;
    // INCHI✔️✔️:          i1 < MAX_NUM_STEREO_BONDS &&
    // INCHI✔️✔️:          ( s1 = at[cur1].stereo_bond_neighbor[i1] ) &&
    // INCHI✔️✔️:             !( k1 = ( at[cur1].neighbor[(int) at[cur1].stereo_bond_ord[i1]] == *n1 ) );
    // INCHI✔️✔️:          i1++) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i2 = 0, k2 = 0;
    // INCHI✔️✔️:          i2 < MAX_NUM_STEREO_BONDS &&
    // INCHI✔️✔️:          ( s2 = at[cur2].stereo_bond_neighbor[i2] ) &&
    // INCHI✔️✔️:             !( k2 = ( at[cur2].neighbor[(int) at[cur2].stereo_bond_ord[i2]] == *n2 ) );
    // INCHI✔️✔️:          i2++) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (k1 != k2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /*  possibly not an error: constit. equivalent atoms on a stereo bond and not on a stereo bond */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (k1 /* yes, it is a stereo bond */ &&
    // INCHI✔️✔️:         ( at[cur1].stereo_bond_parity[i1] != at[cur2].stereo_bond_parity[i2] ||
    // INCHI✔️✔️:           /* PARITY_KNOWN (at[cur1].stereo_bond_parity[i1] ) */  /*  replaced 08-13-2002 with the next: */
    // INCHI✔️✔️:             !PARITY_WELL_DEF( at[cur1].stereo_bond_parity[i1] ) /*  it suffices to check only one parity */
    // INCHI✔️✔️:             ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /*  different or (currently) unknown stereo bond parities */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 1; /*  stereo bonds have known parities */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetAndCheckNextNeighbors

    let canonical_ranks = heap.slice(nCanonRank)?;
    let mut cr1 = if u32::from(*n1) > MAX_ATOMS {
        0
    } else {
        canonical_ranks[usize::from(*n1)]
    };
    let mut cr2 = if u32::from(*n2) > MAX_ATOMS {
        0
    } else {
        canonical_ranks[usize::from(*n2)]
    };
    if GetNextNeighborAndRank(heap, at, cur1, prev1, n1, &mut cr1, nCanonRank)? == 0 {
        return Ok(0);
    }
    if GetNextNeighborAndRank(heap, at, cur2, prev2, n2, &mut cr2, nCanonRank)? == 0 {
        return Ok(0);
    }
    let ranks = heap.slice(nRank)?;
    let visited1 = heap.slice(nVisited1)?;
    let visited2 = heap.slice(nVisited2)?;
    if ranks[usize::from(*n1)] != ranks[usize::from(*n2)]
        || visited1[usize::from(*n1)] != visited2[usize::from(*n2)]
    {
        return Ok(0);
    }

    let atoms = heap.slice(at)?;
    let find_stereo = |cur: AT_RANK, next: AT_RANK| -> (bool, usize) {
        let atom = &atoms[usize::from(cur)];
        let mut index = 0_usize;
        while index < MAX_NUM_STEREO_BONDS as usize {
            if atom.stereo_bond_neighbor[index] == 0 {
                break;
            }
            let order = atom.stereo_bond_ord[index] as usize;
            if atom.neighbor[order] == next {
                return (true, index);
            }
            index += 1;
        }
        (false, index)
    };
    let (is_stereo1, stereo_index1) = find_stereo(cur1, *n1);
    let (is_stereo2, stereo_index2) = find_stereo(cur2, *n2);
    if is_stereo1 != is_stereo2 {
        return Ok(0);
    }
    if is_stereo1
        && (atoms[usize::from(cur1)].stereo_bond_parity[stereo_index1]
            != atoms[usize::from(cur2)].stereo_bond_parity[stereo_index2]
            || !parity_well_defined(i32::from(
                atoms[usize::from(cur1)].stereo_bond_parity[stereo_index1],
            )))
    {
        return Ok(0);
    }
    Ok(1)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn PathsHaveIdenticalKnownParities(
    heap: &mut SourceHeap,
    at: SourceConstPointer<sp_ATOM>,
    prev1: AT_RANK,
    cur1: AT_RANK,
    prev2: AT_RANK,
    cur2: AT_RANK,
    nVisited1: SourceMutPointer<AT_RANK>,
    nVisited2: SourceMutPointer<AT_RANK>,
    nRank: SourceConstPointer<AT_RANK>,
    nCanonRank: SourceConstPointer<AT_RANK>,
    mut nLength: AT_RANK,
) -> Result<AT_RANK, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1530 PathsHaveIdenticalKnownParities
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    AT_RANK PathsHaveIdenticalKnownParities( sp_ATOM       *at,
                                             AT_RANK       prev1,
                                             AT_RANK       cur1,
                                             AT_RANK       prev2,
                                             AT_RANK       cur2,
                                             AT_RANK       *nVisited1,
                                             AT_RANK       *nVisited2,
                                             const AT_RANK *nRank,
                                             const AT_RANK *nCanonRank,
                                             AT_RANK       nLength )
    {
        int k;
        AT_RANK n1, n2;

        nLength++;   /*  number of successfully traversed pairs of atoms */
        nVisited1[cur1] = nLength;
        nVisited2[cur2] = nLength;

        /*  the atoms must be either both stereogenic and have well-defined parities or non-stereogenic at all. */
        if (at[cur1].stereo_atom_parity != at[cur2].stereo_atom_parity ||
             (at[cur1].stereo_atom_parity && !PARITY_WELL_DEF( at[cur1].stereo_atom_parity )) ) /* djb-rwth: addressing LLVM warning */
        {
            return 0;  /*  Reject: Different or unknown in advance parities */
        }

        if (at[cur1].valence != at[cur2].valence)
        {
            return 0;  /*  program error */ /*   <BRKPT> */
        }

        if (at[cur1].valence == 1)
        {
            return nLength; /*  so far success */
        }


        for (k = 1, n1 = MAX_ATOMS + 1, n2 = MAX_ATOMS + 1; k < at[cur1].valence; k++)
        {
            /*  start from 1: since we do not go back, we have only (at[cur1].valence-1) bonds to try */

            if (!GetAndCheckNextNeighbors( at, cur1, prev1, cur2, prev2,
                &n1, &n2, nVisited1, nVisited2,
                nRank, nCanonRank ))
            {
                return 0; /*  different neighbors                       */
            }

            /*  In a DFS we do not traverse already visited atoms */
            if (!nVisited1[n1])
            {
                /*  recursion */
                if (!( nLength = PathsHaveIdenticalKnownParities( at, cur1, n1, cur2, n2, nVisited1, nVisited2, nRank, nCanonRank, nLength ) ))
                {
                    return 0;
                }
            }
        }

        /*  To be on a safe side, recheck after all nVisited[] have been set */
        for (k = 1, n1 = MAX_ATOMS + 1, n2 = MAX_ATOMS + 1; k < at[cur1].valence; k++)
        {
            /*  start from 1: since we do not go back, we have only (at[cur1].valence-1) bonds to try */
            if (!GetAndCheckNextNeighbors( at, cur1, prev1, cur2, prev2,
                &n1, &n2, nVisited1, nVisited2,
                nRank, nCanonRank ))
            {
                return 0; /*  different neighbors */
            }
        }

        return nLength;
    }
    */
    // END INCHI C FUNCTION: PathsHaveIdenticalKnownParities

    nLength = nLength.wrapping_add(1);
    *heap
        .slice_mut(nVisited1)?
        .get_mut(usize::from(cur1))
        .ok_or(SourceHeapError::PointerOutOfBounds)? = nLength;
    *heap
        .slice_mut(nVisited2)?
        .get_mut(usize::from(cur2))
        .ok_or(SourceHeapError::PointerOutOfBounds)? = nLength;

    let (parity1, parity2, valence1, valence2) = {
        let atoms = heap.slice(at)?;
        let atom1 = atoms
            .get(usize::from(cur1))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom2 = atoms
            .get(usize::from(cur2))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        (
            atom1.stereo_atom_parity,
            atom2.stereo_atom_parity,
            atom1.valence,
            atom2.valence,
        )
    };
    if parity1 != parity2 || (parity1 != 0 && !parity_well_defined(i32::from(parity1))) {
        return Ok(0);
    }
    if valence1 != valence2 {
        return Ok(0);
    }
    if valence1 == 1 {
        return Ok(nLength);
    }

    let mut n1 = (MAX_ATOMS + 1) as AT_RANK;
    let mut n2 = (MAX_ATOMS + 1) as AT_RANK;
    for _ in 1..i32::from(valence1) {
        if GetAndCheckNextNeighbors(
            heap,
            at,
            cur1,
            prev1,
            cur2,
            prev2,
            &mut n1,
            &mut n2,
            nVisited1.as_const(),
            nVisited2.as_const(),
            nRank,
            nCanonRank,
        )? == 0
        {
            return Ok(0);
        }
        let was_visited = *heap
            .slice(nVisited1.as_const())?
            .get(usize::from(n1))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if was_visited == 0 {
            nLength = PathsHaveIdenticalKnownParities(
                heap, at, cur1, n1, cur2, n2, nVisited1, nVisited2, nRank, nCanonRank, nLength,
            )?;
            if nLength == 0 {
                return Ok(0);
            }
        }
    }

    n1 = (MAX_ATOMS + 1) as AT_RANK;
    n2 = (MAX_ATOMS + 1) as AT_RANK;
    for _ in 1..i32::from(valence1) {
        if GetAndCheckNextNeighbors(
            heap,
            at,
            cur1,
            prev1,
            cur2,
            prev2,
            &mut n1,
            &mut n2,
            nVisited1.as_const(),
            nVisited2.as_const(),
            nRank,
            nCanonRank,
        )? == 0
        {
            return Ok(0);
        }
    }
    Ok(nLength)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveKnownNonStereoBondParities(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nCanonRank: SourceConstPointer<AT_RANK>,
    nRank: SourceConstPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1608 RemoveKnownNonStereoBondParities
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int RemoveKnownNonStereoBondParities( sp_ATOM       *at,
                                          int           num_atoms,
                                          const AT_RANK *nCanonRank,
                                          const AT_RANK *nRank,
                                          CANON_STAT    *pCS )
    {
        int j, n, m, ret;

        int i1, n1, s2;
        AT_RANK nAtomRank1, nAtomRank2, neigh[3], opposite_atom, *nVisited = NULL;
        ret = 0;
        for (i1 = 0; i1 < num_atoms; i1++)
        {
            if (at[i1].valence != 3 || !at[i1].stereo_bond_neighbor[0])
            {
                continue;
            }
            for (n1 = 0; n1 < MAX_NUM_STEREO_BONDS && ( s2 = at[i1].stereo_bond_neighbor[n1] ); n1++)
            {
                if (!PARITY_CALCULATE( at[i1].stereo_bond_parity[n1] ) && PARITY_WELL_DEF( at[i1].stereo_bond_parity[n1] ))
                {
                    continue;
                }
                opposite_atom = (AT_RANK) ( s2 - 1 );
                /* s2 = at[i1].neighbor[m=(int)at[i1].stereo_bond_ord[n1]]; */
                m = (int) at[i1].stereo_bond_ord[n1];
                for (j = 0, n = 0; j < at[i1].valence; j++)
                {
                    /* if ( at[i1].neighbor[j] != s2 ) */
                    if (j != m)
                    {
                        neigh[n++] = at[i1].neighbor[j];
                    }
                }
                if (n > 2)
                {
                    ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                    goto exit_function;
                }
                if (n != 2 || nRank[(int) neigh[0]] != nRank[(int) neigh[1]])
                {
                    continue; /*  may happen if another half-bond has not a defined parity */
                }
                if (at[i1].nRingSystem == at[(int) neigh[0]].nRingSystem)
                {
                    continue;  /*  no more ring system membership check is necessary because     */
                }              /*  the two neighbors are to be constitutionally equivalent atoms */
                if (!nVisited && !( nVisited = (AT_RANK*) inchi_malloc( sizeof( nVisited[0] )*num_atoms ) ))
                {
                    ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
                    goto exit_function;
                }
                memset( nVisited, 0, sizeof( nVisited[0] )*num_atoms ); /* djb-rwth: memset_s C11/Annex K variant? */
                nVisited[i1] = 1;
                if (PathsHaveIdenticalKnownParities( at, (AT_RANK) i1, neigh[0], (AT_RANK) i1, neigh[1], nVisited, nVisited, nRank, nCanonRank, 1 ))
                {
                    if (!RemoveOneStereoBond( at, i1, /* atom number*/ n1 /* stereo bond number*/ ))
                    {
                        ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                        goto exit_function;
                    }
                    n1--; /*  cycle counter may temporarily become negative */
                    /*  Remove from pCS */
                    nAtomRank1 = inchi_max( nCanonRank[i1], nCanonRank[opposite_atom] );
                    nAtomRank2 = inchi_min( nCanonRank[i1], nCanonRank[opposite_atom] );
                    for (n = 0, m = pCS->nLenLinearCTStereoDble - 1; n <= m; n++)
                    {
                        if (pCS->LinearCTStereoDble[n].at_num1 == nAtomRank1 &&
                             pCS->LinearCTStereoDble[n].at_num2 == nAtomRank2)
                        {
                            if (n < m)
                            { /*  remove pCS->LinearCTStereoDble[n] */
                                memmove( pCS->LinearCTStereoDble + n, pCS->LinearCTStereoDble + n + 1, ( (long long)m - (long long)n ) * sizeof( pCS->LinearCTStereoDble[0] ) ); /* djb-rwth: cast operators added */
                            }
                            pCS->nLenLinearCTStereoDble--;
    #if ( bRELEASE_VERSION == 0 )
                            pCS->bExtract |= EXTR_KNOWN_USED_TO_REMOVE_PARITY;
    #endif
                            m = -1;   /*  set flag "found" */
                            break;
                        }
                    }
                    if (m >= 0)
                    {
                        ret = CT_STEREOCOUNT_ERR;  /*  bond not found  <BRKPT> */
                        goto exit_function;
                    }
                    ret++;  /*  number of removed known in advance non-stereo bonds */
                }
            }
        }

    exit_function:
        if (nVisited)
        {
            inchi_free( nVisited );
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: RemoveKnownNonStereoBondParities

    let mut visited: Option<SourceMutPointer<AT_RANK>> = None;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut ret = 0_i32;
        for i1 in 0..num_atoms {
            let should_scan = {
                let atoms = heap.slice(at.as_const())?;
                let atom = atoms
                    .get(i1 as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                atom.valence == 3 && atom.stereo_bond_neighbor[0] != 0
            };
            if !should_scan {
                continue;
            }

            let mut n1 = 0_i32;
            while n1 < MAX_NUM_STEREO_BONDS as i32 {
                let (s2, bond_parity, stereo_order, valence) = {
                    let atoms = heap.slice(at.as_const())?;
                    let atom = &atoms[i1 as usize];
                    (
                        atom.stereo_bond_neighbor[n1 as usize],
                        atom.stereo_bond_parity[n1 as usize],
                        atom.stereo_bond_ord[n1 as usize],
                        atom.valence,
                    )
                };
                if s2 == 0 {
                    break;
                }
                if !parity_calculate(i32::from(bond_parity))
                    && parity_well_defined(i32::from(bond_parity))
                {
                    n1 += 1;
                    continue;
                }

                let opposite_atom = s2.wrapping_sub(1);
                let m = i32::from(stereo_order);
                let mut neigh = [0 as AT_RANK; 3];
                let mut n = 0_i32;
                for j in 0..i32::from(valence) {
                    if j != m {
                        neigh[n as usize] =
                            heap.slice(at.as_const())?[i1 as usize].neighbor[j as usize];
                        n += 1;
                    }
                }
                if n > 2 {
                    return Ok(CT_STEREOBOND_ERROR);
                }
                if n != 2 {
                    n1 += 1;
                    continue;
                }
                let ranks = heap.slice(nRank)?;
                if ranks[usize::from(neigh[0])] != ranks[usize::from(neigh[1])] {
                    n1 += 1;
                    continue;
                }
                let same_ring = {
                    let atoms = heap.slice(at.as_const())?;
                    atoms[i1 as usize].nRingSystem == atoms[usize::from(neigh[0])].nRingSystem
                };
                if same_ring {
                    n1 += 1;
                    continue;
                }

                let visited_pointer = if let Some(pointer) = visited {
                    pointer
                } else {
                    let atom_count = usize::try_from(num_atoms)
                        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
                    let pointer = match heap.allocate(vec![0 as AT_RANK; atom_count]) {
                        Ok(pointer) => pointer,
                        Err(SourceHeapError::AllocationFailed) => return Ok(CT_OUT_OF_RAM),
                        Err(error) => return Err(error),
                    };
                    visited = Some(pointer);
                    pointer
                };
                heap.slice_mut(visited_pointer)?.fill(0);
                heap.slice_mut(visited_pointer)?[i1 as usize] = 1;
                if PathsHaveIdenticalKnownParities(
                    heap,
                    at.as_const(),
                    i1 as AT_RANK,
                    neigh[0],
                    i1 as AT_RANK,
                    neigh[1],
                    visited_pointer,
                    visited_pointer,
                    nRank,
                    nCanonRank,
                    1,
                )? != 0
                {
                    if RemoveOneStereoBond(heap, at, i1, n1)? == 0 {
                        return Ok(CT_STEREOBOND_ERROR);
                    }
                    n1 -= 1;
                    let canonical_ranks = heap.slice(nCanonRank)?;
                    let first = canonical_ranks[i1 as usize];
                    let second = canonical_ranks[usize::from(opposite_atom)];
                    let atom_rank1 = first.max(second);
                    let atom_rank2 = first.min(second);
                    let mut m = pCS.nLenLinearCTStereoDble - 1;
                    let mut n = 0_i32;
                    while n <= m {
                        let record =
                            heap.slice(pCS.LinearCTStereoDble.as_const())?[n as usize].clone();
                        if record.at_num1 == atom_rank1 && record.at_num2 == atom_rank2 {
                            if n < m {
                                let records = heap.slice_mut(pCS.LinearCTStereoDble)?;
                                for index in n as usize..m as usize {
                                    records[index] = records[index + 1].clone();
                                }
                            }
                            pCS.nLenLinearCTStereoDble -= 1;
                            m = -1;
                            break;
                        }
                        n += 1;
                    }
                    if m >= 0 {
                        return Ok(CT_STEREOCOUNT_ERR);
                    }
                    ret = ret.wrapping_add(1);
                }
                n1 += 1;
            }
        }
        Ok(ret)
    })();

    let cleanup = if let Some(pointer) = visited {
        heap.free(pointer)
    } else {
        Ok(())
    };
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn SetKnownStereoCenterParities(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nCanonRank: SourceConstPointer<AT_RANK>,
    nRank: SourceConstPointer<AT_RANK>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1714 SetKnownStereoCenterParities
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int SetKnownStereoCenterParities( CANON_GLOBALS *pCG,
                                      sp_ATOM       *at,
                                      int           num_atoms,
                                      const AT_RANK *nCanonRank,
                                      const AT_RANK *nRank,
                                      const AT_RANK *nAtomNumber )
    {
        int i, j, n, m, j1, k, num_neigh, iMax, trans_i, trans_k, prev_trans, num_set;
        AT_RANK nAtomRank;
        AT_RANK nNeighRank[MAX_NUM_STEREO_ATOM_NEIGH];
        AT_RANK nNeighCanonRank[MAX_NUM_STEREO_ATOM_NEIGH];

        for (i = 0, num_set = 0; i < num_atoms; i++)
        {
            if (!at[i].parity || at[i].stereo_bond_neighbor[0])
            {
                continue;
            }
            if (at[i].stereo_atom_parity != AB_PARITY_CALC ||
                 !PARITY_WELL_DEF( at[i].parity ))
            {
                continue;
            }
            num_neigh = at[i].valence;
            for (j = 0; j < num_neigh; j++)
            {
                nNeighRank[j] = nRank[(int) at[i].neighbor[j]];
            }
            nAtomRank = nRank[i];
            if (num_neigh == 1)
            {
                /* other neighbors must be implicit H */
                at[i].stereo_atom_parity = at[i].parity;
                trans_i = 0;
            }
            else
            {
                /* sort constitutional equivalence ranks of the neighbors */
                trans_i = insertions_sort( pCG, nNeighRank, num_neigh, sizeof( nNeighRank[0] ), comp_AT_RANK );
                for (j = 1; j < num_neigh; j++)
                {
                    if (nNeighRank[j - 1] == nNeighRank[j])
                    {
                        break; /* at[i] has consitutionally identical neighbors */
                    }
                }
                if (j < num_neigh)
                {
                    /*  at least 2 neighbors are const. identical; parity cannot be calculated at this time */
                    continue; /*  try next stereo atom */
                }
            }
            prev_trans = -1;
            trans_k = 0;
            /*  find neighbors of constitutionally equivalent stereo centers */
            /*  and at[i] parities in case those centers are mapped on at[i] */
            for (iMax = (int) nAtomRank - 1, j1 = 0;
                    j1 <= iMax && nAtomRank == nRank[k = (int) nAtomNumber[iMax - j1]];
                        j1++)
            {
                /*  at[k] is constitutionally equivalent to at[i] */
                if ((int) at[k].valence != num_neigh)
                {
                    return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                }
                /* -- commented out to accept  non-stereogenic atoms since     --
                   -- they may participate in mapping stereocenters 12-16-2003 --
                if ( !PARITY_VAL(at[k].parity) ) {
                    continue; // not a stereogenic atom
                }
                */
                for (j = 0, m = 0; m < num_neigh; m++)
                {
                    for (n = 0; n < num_neigh; n++)
                    {
                        if (nRank[(int) at[k].neighbor[n]] == nNeighRank[m])
                        {
                            /* save canonical numbers (ranks) of the neighbors in
                             * order of increasing constit. equivalence ranks */
                            nNeighCanonRank[m] = nCanonRank[(int) at[k].neighbor[n]];
                            j++;
                            break;
                        }
                    }
                }
                if (j != num_neigh)
                {
                    return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                }
                trans_k = insertions_sort( pCG, nNeighCanonRank, num_neigh, sizeof( nNeighCanonRank[0] ), comp_AT_RANK );
                trans_k %= 2;
                if (prev_trans < 0)
                {
                    prev_trans = trans_k;
                }
                else if (trans_k != prev_trans)
                {
                    /*  different mappings may produce different parities. Cannot find the parity at this time */
                    /*  this may happen when a set of constit. equivalent atoms has non-contiguous canonical numbers */
                    break;
                }
            }
            if (trans_k == prev_trans)
            {
                at[i].stereo_atom_parity = 2 - ( at[i].parity + trans_i + prev_trans ) % 2;
                num_set++;
            }
        }

        return num_set;
    }
    */
    // END INCHI C FUNCTION: SetKnownStereoCenterParities

    fn sort_ranks(values: &mut [AT_RANK], count: usize) -> Result<i32, SourceHeapError> {
        let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(values);
        insertions_sort(
            bytes,
            count,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                Ok(i32::from(left) - i32::from(right))
            },
        )
    }

    let _ = pCG;
    let mut num_set = 0_i32;
    for i in 0..num_atoms {
        let (parity, stereo_atom_parity, has_stereo_bond, num_neigh, atom_neighbors) = {
            let atoms = heap.slice(at.as_const())?;
            let atom = atoms
                .get(i as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (
                atom.parity,
                atom.stereo_atom_parity,
                atom.stereo_bond_neighbor[0] != 0,
                i32::from(atom.valence),
                atom.neighbor,
            )
        };
        if parity == 0 || has_stereo_bond {
            continue;
        }
        if stereo_atom_parity != AB_PARITY_CALC as i8 || !parity_well_defined(i32::from(parity)) {
            continue;
        }
        let count = usize::try_from(num_neigh).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if count > MAX_NUM_STEREO_ATOM_NEIGH as usize {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let mut neighbor_ranks = [0 as AT_RANK; MAX_NUM_STEREO_ATOM_NEIGH as usize];
        for j in 0..count {
            neighbor_ranks[j] = heap.slice(nRank)?[usize::from(atom_neighbors[j])];
        }
        let atom_rank = heap.slice(nRank)?[i as usize];
        let trans_i = if num_neigh == 1 {
            heap.slice_mut(at)?[i as usize].stereo_atom_parity = parity;
            0
        } else {
            let transactions = sort_ranks(&mut neighbor_ranks, count)?;
            let mut j = 1_usize;
            while j < count && neighbor_ranks[j - 1] != neighbor_ranks[j] {
                j += 1;
            }
            if j < count {
                continue;
            }
            transactions
        };

        let mut previous_transactions = -1_i32;
        let mut transactions_for_equivalent = 0_i32;
        let maximum = i32::from(atom_rank) - 1;
        let mut j1 = 0_i32;
        while j1 <= maximum {
            let k = heap.slice(nAtomNumber)?[(maximum - j1) as usize];
            if heap.slice(nRank)?[usize::from(k)] != atom_rank {
                break;
            }
            let equivalent = heap
                .slice(at.as_const())?
                .get(usize::from(k))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(equivalent.valence) != num_neigh {
                return Ok(CT_STEREOCOUNT_ERR);
            }
            let equivalent_neighbors = equivalent.neighbor;
            let mut neighbor_canonical_ranks = [0 as AT_RANK; MAX_NUM_STEREO_ATOM_NEIGH as usize];
            let mut matched = 0_i32;
            for m in 0..count {
                for n in 0..count {
                    let neighbor = equivalent_neighbors[n];
                    if heap.slice(nRank)?[usize::from(neighbor)] == neighbor_ranks[m] {
                        neighbor_canonical_ranks[m] =
                            heap.slice(nCanonRank)?[usize::from(neighbor)];
                        matched += 1;
                        break;
                    }
                }
            }
            if matched != num_neigh {
                return Ok(CT_STEREOCOUNT_ERR);
            }
            transactions_for_equivalent = sort_ranks(&mut neighbor_canonical_ranks, count)? % 2;
            if previous_transactions < 0 {
                previous_transactions = transactions_for_equivalent;
            } else if transactions_for_equivalent != previous_transactions {
                break;
            }
            j1 += 1;
        }
        if transactions_for_equivalent == previous_transactions {
            heap.slice_mut(at)?[i as usize].stereo_atom_parity =
                (2 - (i32::from(parity) + trans_i + previous_transactions) % 2) as i8;
            num_set = num_set.wrapping_add(1);
        }
    }
    Ok(num_set)
}

#[allow(non_snake_case)]
pub(crate) fn RemoveKnownNonStereoCenterParities(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nCanonRank: SourceConstPointer<AT_RANK>,
    nRank: SourceConstPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1833 RemoveKnownNonStereoCenterParities
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int RemoveKnownNonStereoCenterParities( CANON_GLOBALS *pCG,
                                            sp_ATOM       *at,
                                            int           num_atoms,
                                            const AT_RANK *nCanonRank,
                                            const AT_RANK *nRank,
                                            CANON_STAT    *pCS )
    {
        int i, j, n, m, k, num_neigh, ret = 0;
        /*AT_RANK nAtomRank;*/
        AT_RANK nNeighRank[MAX_NUM_STEREO_ATOM_NEIGH], nNeighOrd[MAX_NUM_STEREO_ATOM_NEIGH];
        AT_RANK *nVisited = NULL;

        for (i = 0; i < num_atoms; i++)
        {
            if (!at[i].parity || at[i].stereo_bond_neighbor[0])
            {
                continue;
            }
            if (!PARITY_CALCULATE( at[i].stereo_atom_parity ) && PARITY_WELL_DEF( at[i].stereo_atom_parity ))
            {
                continue;
            }
            num_neigh = at[i].valence;
            for (j = 0; j < num_neigh; j++)
            {
                nNeighRank[j] = nRank[(int) at[i].neighbor[j]];
                nNeighOrd[j] = j;
            }
            /*nAtomRank = nRank[i];*/
            if (num_neigh == 1)
            {
                continue;
            }
            pCG->m_pn_RankForSort = nNeighRank;
            insertions_sort( pCG, nNeighOrd, num_neigh, sizeof( nNeighRank[0] ), CompRanksOrd );
            for (j = k = 1; k && j < num_neigh; j++)
            {
                if (at[i].nRingSystem != at[(int) at[i].neighbor[(int) nNeighOrd[j]]].nRingSystem &&
                     /*  no more ring system membership check is necessary because */
                     /*  the two neighbors are to be constitutionally equivalent atoms: */
                    nNeighRank[nNeighOrd[j - 1]] == nNeighRank[nNeighOrd[j]])
                {
                    k = j;
                    do
                    {
                        if (!nVisited && !( nVisited = (AT_RANK*) inchi_malloc( sizeof( nVisited[0] )*num_atoms ) ))
                        {
                            ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
                            goto exit_function;
                        }
                        memset( nVisited, 0, sizeof( nVisited[0] )*num_atoms ); /* djb-rwth: memset_s C11/Annex K variant? */
                        nVisited[i] = 1;
                        if (PathsHaveIdenticalKnownParities( at, (AT_RANK) i, at[i].neighbor[(int) nNeighOrd[j - 1]],
                            (AT_RANK) i, at[i].neighbor[(int) nNeighOrd[k]],
                            nVisited, nVisited, nRank, nCanonRank, 1 ))
                        {
                            at[i].parity = 0; /*  remove parity */
                            at[i].stereo_atom_parity = 0;
                            at[i].final_parity = 0;
                            /* at[i].bHasStereoOrEquToStereo = 0; */
                            for (n = 0, m = pCS->nLenLinearCTStereoCarb - 1; n <= m; n++)
                            {
                                if (pCS->LinearCTStereoCarb[n].at_num == nCanonRank[i])
                                {
                                    if (n < m)
                                    {    /*  remove pCS->LinearCTStereoCarb[n] */
                                        memmove( pCS->LinearCTStereoCarb + n, pCS->LinearCTStereoCarb + n + 1, ( (long long)m - (long long)n ) * sizeof( pCS->LinearCTStereoCarb[0] ) ); /* djb-rwth: cast operators added */
                                    }
                                    pCS->nLenLinearCTStereoCarb--;
                                    k = 0;

    #if ( bRELEASE_VERSION == 0 )
                                    pCS->bExtract |= EXTR_KNOWN_USED_TO_REMOVE_PARITY;
    #endif

                                    break;
                                }
                            }
                            if (k)
                            {
                                ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                                goto exit_function;
                            }
                            ret++;
                            break;
                        }
                    }
                    while (++k < num_neigh && nNeighRank[nNeighOrd[j - 1]] == nNeighRank[nNeighOrd[k]]);
                }
            }
        }

    exit_function:
        if (nVisited)
        {
            inchi_free( nVisited );
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: RemoveKnownNonStereoCenterParities

    fn sort_orders(
        orders: &mut [AT_RANK],
        ranks: &[AT_RANK],
        count: usize,
    ) -> Result<(), SourceHeapError> {
        let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(orders);
        insertions_sort(
            bytes,
            count,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                let difference =
                    i32::from(ranks[usize::from(left)]) - i32::from(ranks[usize::from(right)]);
                Ok(if difference == 0 {
                    i32::from(left) - i32::from(right)
                } else {
                    difference
                })
            },
        )?;
        Ok(())
    }

    let _ = pCG;
    let mut visited: Option<SourceMutPointer<AT_RANK>> = None;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut ret = 0_i32;
        for i in 0..num_atoms {
            let atom_snapshot = heap
                .slice(at.as_const())?
                .get(i as usize)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom_snapshot.parity == 0 || atom_snapshot.stereo_bond_neighbor[0] != 0 {
                continue;
            }
            let center_parity = i32::from(atom_snapshot.stereo_atom_parity);
            if !parity_calculate(center_parity) && parity_well_defined(center_parity) {
                continue;
            }
            let count = usize::try_from(atom_snapshot.valence)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if count > MAX_NUM_STEREO_ATOM_NEIGH as usize {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            let mut neighbor_ranks = [0 as AT_RANK; MAX_NUM_STEREO_ATOM_NEIGH as usize];
            let mut neighbor_orders = [0 as AT_RANK; MAX_NUM_STEREO_ATOM_NEIGH as usize];
            for j in 0..count {
                neighbor_ranks[j] = heap.slice(nRank)?[usize::from(atom_snapshot.neighbor[j])];
                neighbor_orders[j] = j as AT_RANK;
            }
            if count == 1 {
                continue;
            }
            sort_orders(&mut neighbor_orders, &neighbor_ranks, count)?;
            let mut j = 1_usize;
            let mut k = 1_usize;
            while k != 0 && j < count {
                let current_order = usize::from(neighbor_orders[j]);
                let previous_order = usize::from(neighbor_orders[j - 1]);
                let current_neighbor = atom_snapshot.neighbor[current_order];
                if atom_snapshot.nRingSystem
                    != heap.slice(at.as_const())?[usize::from(current_neighbor)].nRingSystem
                    && neighbor_ranks[previous_order] == neighbor_ranks[current_order]
                {
                    k = j;
                    loop {
                        let visited_pointer = if let Some(pointer) = visited {
                            pointer
                        } else {
                            let pointer = match heap
                                .allocate(vec![0 as AT_RANK; num_atoms as usize])
                            {
                                Ok(pointer) => pointer,
                                Err(SourceHeapError::AllocationFailed) => return Ok(CT_OUT_OF_RAM),
                                Err(error) => return Err(error),
                            };
                            visited = Some(pointer);
                            pointer
                        };
                        heap.slice_mut(visited_pointer)?.fill(0);
                        heap.slice_mut(visited_pointer)?[i as usize] = 1;
                        if PathsHaveIdenticalKnownParities(
                            heap,
                            at.as_const(),
                            i as AT_RANK,
                            atom_snapshot.neighbor[previous_order],
                            i as AT_RANK,
                            atom_snapshot.neighbor[usize::from(neighbor_orders[k])],
                            visited_pointer,
                            visited_pointer,
                            nRank,
                            nCanonRank,
                            1,
                        )? != 0
                        {
                            let atom = &mut heap.slice_mut(at)?[i as usize];
                            atom.parity = 0;
                            atom.stereo_atom_parity = 0;
                            atom.final_parity = 0;
                            let mut found = false;
                            let last = pCS.nLenLinearCTStereoCarb - 1;
                            let mut n = 0_i32;
                            while n <= last {
                                if heap.slice(pCS.LinearCTStereoCarb.as_const())?[n as usize].at_num
                                    == heap.slice(nCanonRank)?[i as usize]
                                {
                                    if n < last {
                                        let records = heap.slice_mut(pCS.LinearCTStereoCarb)?;
                                        for index in n as usize..last as usize {
                                            records[index] = records[index + 1].clone();
                                        }
                                    }
                                    pCS.nLenLinearCTStereoCarb -= 1;
                                    k = 0;
                                    found = true;
                                    break;
                                }
                                n += 1;
                            }
                            if !found {
                                return Ok(CT_STEREOCOUNT_ERR);
                            }
                            ret = ret.wrapping_add(1);
                            break;
                        }
                        k += 1;
                        if k >= count
                            || neighbor_ranks[previous_order]
                                != neighbor_ranks[usize::from(neighbor_orders[k])]
                        {
                            break;
                        }
                    }
                }
                j += 1;
            }
        }
        Ok(ret)
    })();
    let cleanup = if let Some(pointer) = visited {
        heap.free(pointer)
    } else {
        Ok(())
    };
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[allow(non_snake_case)]
pub(crate) fn MarkKnownEqualStereoCenterParities(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nRank: SourceConstPointer<AT_RANK>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1942 MarkKnownEqualStereoCenterParities
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int MarkKnownEqualStereoCenterParities( sp_ATOM       *at,
                                            int           num_atoms,
                                            const AT_RANK *nRank,
                                            const AT_RANK *nAtomNumber )
    {
        int i, j1, k, num_centers, iMax, bDifferentParities;
        AT_RANK nAtomRank;
        int  parity, parity_k;

        num_centers = 0;
        for (i = 0; i < num_atoms; i++)
        {
            if (!at[i].parity || at[i].stereo_bond_neighbor[0])
            {
                continue;
            }
            if (at[i].bHasStereoOrEquToStereo)
            {
                continue; /* already marked */
            }
            if ( /*!PARITY_KNOWN(at[i].stereo_atom_parity) ||*/ ( at[i].stereo_atom_parity & KNOWN_PARITIES_EQL ))
            {
                continue;
            }
            parity = PARITY_VAL( at[i].stereo_atom_parity );
            if (parity == AB_PARITY_NONE)
            {
                continue;
            }
            nAtomRank = nRank[i];
            bDifferentParities = -1;
            /*  find constitutionally equivalent stereo centers and compare their known at this time parities */
            for (iMax = (int) nAtomRank - 1, j1 = 0; j1 <= iMax && nAtomRank == nRank[k = (int) nAtomNumber[iMax - j1]]; j1++)
            {
                /*  at[k] is constitutionally equivalent to at[i] */
                parity_k = PARITY_VAL( at[k].stereo_atom_parity );
                if (parity_k != parity)
                {
                    bDifferentParities = 1;
                }
                else if (parity_k == parity && bDifferentParities < 0)
                {
                    bDifferentParities = 0;
                }
                if (!parity_k)
                {
                    at[k].bHasStereoOrEquToStereo = 2;
                }
                else if (!at[k].bHasStereoOrEquToStereo)
                {
                    at[k].bHasStereoOrEquToStereo = 1;
                }
            }
            if (0 == bDifferentParities && PARITY_KNOWN( parity ))
            {
                for (iMax = (int) nAtomRank - 1, j1 = 0; j1 <= iMax && nAtomRank == nRank[k = (int) nAtomNumber[iMax - j1]]; j1++)
                {
                    /*  at[k] is constitutionally equivalent to at[i] */
                    at[k].stereo_atom_parity |= KNOWN_PARITIES_EQL;
                    num_centers++;
                }
            }
        }

        return num_centers;
    }
    */
    // END INCHI C FUNCTION: MarkKnownEqualStereoCenterParities

    let mut num_centers = 0_i32;
    for i in 0..num_atoms {
        let atom = heap
            .slice(at.as_const())?
            .get(i as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.parity == 0
            || atom.stereo_bond_neighbor[0] != 0
            || atom.bHasStereoOrEquToStereo != 0
            || i32::from(atom.stereo_atom_parity) & 0x40 != 0
        {
            continue;
        }
        let parity = i32::from(atom.stereo_atom_parity) & BITS_PARITY;
        if parity == 0 {
            continue;
        }
        let atom_rank = heap.slice(nRank)?[i as usize];
        let maximum = i32::from(atom_rank) - 1;
        let mut different_parities = -1_i32;
        let mut j1 = 0_i32;
        while j1 <= maximum {
            let k = heap.slice(nAtomNumber)?[(maximum - j1) as usize];
            if heap.slice(nRank)?[usize::from(k)] != atom_rank {
                break;
            }
            let parity_k = i32::from(heap.slice(at.as_const())?[usize::from(k)].stereo_atom_parity)
                & BITS_PARITY;
            if parity_k != parity {
                different_parities = 1;
            } else if different_parities < 0 {
                different_parities = 0;
            }
            let equivalent = &mut heap.slice_mut(at)?[usize::from(k)];
            if parity_k == 0 {
                equivalent.bHasStereoOrEquToStereo = 2;
            } else if equivalent.bHasStereoOrEquToStereo == 0 {
                equivalent.bHasStereoOrEquToStereo = 1;
            }
            j1 += 1;
        }
        if different_parities == 0 && parity_known(parity) {
            let mut j1 = 0_i32;
            while j1 <= maximum {
                let k = heap.slice(nAtomNumber)?[(maximum - j1) as usize];
                if heap.slice(nRank)?[usize::from(k)] != atom_rank {
                    break;
                }
                heap.slice_mut(at)?[usize::from(k)].stereo_atom_parity |= 0x40;
                num_centers = num_centers.wrapping_add(1);
                j1 += 1;
            }
        }
    }
    Ok(num_centers)
}

#[allow(non_snake_case)]
pub(crate) fn UnmarkNonStereo(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nRank: SourceConstPointer<AT_RANK>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
    bIsotopic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:322 UnmarkNonStereo
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int UnmarkNonStereo( CANON_GLOBALS *pCG,
                         sp_ATOM       *at,
                         int           num_atoms,
                         const AT_RANK *nRank,
                         const AT_RANK *nAtomNumber,
                         int           bIsotopic )
    {
        int i, i1, i2, j, k, k1, k2, kn /* neigh*/, val, ic/* center*/, jc, num_implicit_H;
        int num_neighbors_with_parity, num_no_parity_atoms, num_removed_parities = -1, num_removed_parities0;
        AT_RANK nNeighborNumber[MAX_NUM_STEREO_ATOM_NEIGH];
        AT_RANK nPrevAtomRank, nPrevNeighRank;
    #ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
        int num_in_same_ring_system = 1, nRingSystem, num_with_eq_neigh_in_same_ring_system = 0; /* djb-rwth: although unlikely to ever occur, uninitialised num_in_same_ring_system variable can lead to garbage value, including 0 which leads to various errors and inconsistency with 1.06 outputs -- function rewriting and discussion required */
    #endif


        S_CHAR *visited = (S_CHAR *) inchi_malloc( num_atoms * sizeof( visited[0] ) );

        if (!visited)
        {
            goto exit_function;
        }

        num_removed_parities = 0;
        num_no_parity_atoms = 0;

        do
        {
            num_removed_parities0 = num_removed_parities;

            for (i = i1 = 0, nPrevAtomRank = 0; i <= num_atoms; i++)
            {
                /*  bounds violation check (i!=num_atoms) added 6-21-2002 */
                if (i == num_atoms || nPrevAtomRank != nRank[j = nAtomNumber[i]]
                     /* at[j].parity && 1 < at[j].valence && at[j].valence < MAX_NUM_STEREO_ATOM_NEIGH*/)
                {
                    /*  end of constitutionally equivalent atoms sequence */
                    /* nPrevRank = nRank[j]; */
                    i2 = i;
                    if (i2 - i1 > num_no_parity_atoms /*&& at[jc = nAtomNumber[i1]].parity*/)
                    {
                        /*  at[nAtomNumber[i1]]..at[nAtomNumber[i2-1]] are constitutionally equivalent and some of them have parity */
                        jc = nAtomNumber[i1];
                        num_no_parity_atoms = 0;
                        val = at[jc].valence; /*  all equivalent atoms have equal valences, etc. (except parities) */
                        num_implicit_H = at[jc].endpoint ? 0 : at[jc].num_H;
                        /*  Only atoms with valence <= MAX_NUM_STEREO_ATOM_NEIGH may have parity. However, check: */
                        if (val + num_implicit_H > MAX_NUM_STEREO_ATOM_NEIGH)
                        {
                            continue;  /*  program error ??? */ /*   <BRKPT> */
                        }
                        for (k = 0; k < val; k++)
                        {
                            nNeighborNumber[k] = k; /*  initialize an array of indexes for sorting */
                        }
                        /*  check parities */
                        for (ic = i1; ic < i2; ic++)
                        {
                            jc = nAtomNumber[ic];
                            /*  sort neighbors according to their canon. equivalence ranks */
                            pCG->m_pNeighborsForSort = at[jc].neighbor;
                            pCG->m_pn_RankForSort = nRank;
                            insertions_sort( pCG, nNeighborNumber, val, sizeof( nNeighborNumber[0] ), CompNeighborsAT_NUMBER );
                            num_neighbors_with_parity = -1; /*  non-zero */
                            for (k = k1 = 0, nPrevNeighRank = 0; k <= val; k++)
                            {
                                if (k == val || nPrevNeighRank != nRank[at[jc].neighbor[nNeighborNumber[k]]])
                                {
                                    k2 = k;
                                    if (k2 - k1 > 1)
                                    {
                                        /*  found 2 or more constitutionally equivalent neighbors */
                                        /*  Check if they have only non-stereogenic neighbors */
    #ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                                        num_in_same_ring_system = nRingSystem = 0;
                                        for (kn = k1; kn < k2; kn++) /* djb-rwth: removing redundant code */
                                        {
                                            int nCurNeighRingSystem = at[(int) at[jc].neighbor[nNeighborNumber[kn]]].nRingSystem;
                                            if (!nRingSystem)
                                            {
                                                nRingSystem = nCurNeighRingSystem;
                                            }
                                            else
                                            {
                                                num_in_same_ring_system += ( nRingSystem == nCurNeighRingSystem );
                                            }
                                        }
    #endif

                                        for (kn = k1, num_neighbors_with_parity = 0; kn < k2; kn++)
                                        {
                                            memset( visited, 0, num_atoms * sizeof( visited[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                                            visited[jc] = 1; /*  starting point; the only atom with parity */
                                            num_neighbors_with_parity +=
                                                find_atoms_with_parity( at, visited, jc, (int) at[jc].neighbor[nNeighborNumber[kn]] );
                                        }
                                    }
                                    /* if ( !num_neighbors_with_parity ) */
    #ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                                    if (!num_neighbors_with_parity && !num_in_same_ring_system)
    #else
                                    if (!num_neighbors_with_parity)
    #endif
                                    {
                                        break; /*  at[jc] cannot have defined parity */
                                    }
                                    if (k + 1 < val)
                                    {
                                        k1 = k; /*  at least 2 more neighbors left */
                                        nPrevNeighRank = nRank[at[jc].neighbor[nNeighborNumber[k]]];
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                            }
                            if (num_implicit_H > 1)
                            {
                                if ((bIsotopic && ( at[jc].num_iso_H[0] > 1 ||
                                    at[jc].num_iso_H[1] > 1 ||
                                    at[jc].num_iso_H[2] > 1 )) ||
                                      num_implicit_H > NUM_H_ISOTOPES ||
                                      !bIsotopic) /* djb-rwth: addressing LLVM warning */
                                {
                                    num_neighbors_with_parity = 0;
                                }
                            }
                            /*  increment if: */
                            /*  (a) constitutionally equivalent neighbors do exist, and */
                            /*  (b) all constitutionally equivalent neighbors do not have parity, and */
                            /*  (c) all constitutionally equivalent neighbors are not connected to atoms with parity */
                            num_no_parity_atoms += !num_neighbors_with_parity;
    #ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                            num_with_eq_neigh_in_same_ring_system += ( num_in_same_ring_system != 0 ); /* djb-rwth: initialisation of num_in_same_ring_system is required to avoid garbage value */
    #endif
                        }
    #ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                        if (num_no_parity_atoms == i2 - i1 && num_with_eq_neigh_in_same_ring_system != i2 - i1)
    #else
                        if (num_no_parity_atoms == i2 - i1)
    #endif

                        {
                            /*  all atoms at[nAtomNumber[i1]]..at[nAtomNumber[i2-1]] cannot be */
                            /*  stereo centers or belong to stereo bonds */
                            for (ic = i1; ic < i2; ic++)
                            {
                                int jn;
                                jc = nAtomNumber[ic];
                                at[jc].parity = 0; /*  remove parity */
                                at[jc].stereo_atom_parity = 0;
                                at[jc].final_parity = 0;
                                at[jc].bHasStereoOrEquToStereo = 0;
                                /*  remove stereo bonds */
                                for (k = 0; k < MAX_NUM_STEREO_BOND_NEIGH && ( jn = at[jc].stereo_bond_neighbor[k] ); k++)
                                {
                                    jn--; /*  stereo bond neighbor */
                                    /*  opposite end */
                                    for (k1 = 0; k1 < MAX_NUM_STEREO_BOND_NEIGH && ( kn = at[jn].stereo_bond_neighbor[k1] ); k1++)
                                    {
                                        if (kn - 1 == jc)
                                        {
                                            RemoveHalfStereoBond( at, jn, k1 );
                                            break;
                                        }
                                    }
                                    /*  at at[jc] stereo bond end; since references to all at[jc] */
                                    /*  stereo bond neighbors are to be removed, do not shift them */
                                    at[jc].stereo_bond_neighbor[k] = 0;
                                    at[jc].stereo_bond_ord[k] = 0;
                                    at[jc].stereo_bond_z_prod[k] = 0;
                                    at[jc].stereo_bond_parity[k] = 0;
                                }
                            }
                            num_removed_parities += num_no_parity_atoms;
                        }
                    }
                    if (i < num_atoms)
                    {
                        nPrevAtomRank = nRank[j];
                        i1 = i;
                    }
                    num_no_parity_atoms = 0;
                }
                num_no_parity_atoms += ( i < num_atoms && !at[j].parity );
            }
        }
        while (num_removed_parities != num_removed_parities0);

    exit_function:
        if (visited)
        {
            inchi_free( visited );
        }

        return num_removed_parities;
    }
    */
    // END INCHI C FUNCTION: UnmarkNonStereo

    fn sort_neighbor_numbers(
        values: &mut [AT_RANK],
        count: usize,
        neighbors: &[AT_RANK],
        ranks: &[AT_RANK],
    ) -> Result<(), SourceHeapError> {
        let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(values);
        insertions_sort(
            bytes,
            count,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                let left = usize::from(AT_RANK::from_ne_bytes([left[0], left[1]]));
                let right = usize::from(AT_RANK::from_ne_bytes([right[0], right[1]]));
                Ok(i32::from(ranks[usize::from(neighbors[left])])
                    - i32::from(ranks[usize::from(neighbors[right])]))
            },
        )?;
        Ok(())
    }

    let _ = pCG;
    let atom_count = usize::try_from(num_atoms)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let visited = match heap.allocate(vec![0 as S_CHAR; atom_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(-1),
        Err(error) => return Err(error),
    };
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut num_removed_parities = 0_i32;
        let mut num_no_parity_atoms = 0_i32;
        let mut num_in_same_ring_system = 1_i32;
        let mut num_with_eq_neigh_in_same_ring_system = 0_i32;
        loop {
            let removed_before = num_removed_parities;
            let mut i = 0_i32;
            let mut i1 = 0_i32;
            let mut previous_atom_rank: AT_RANK = 0;
            while i <= num_atoms {
                let current = if i < num_atoms {
                    Some(heap.slice(nAtomNumber)?[i as usize])
                } else {
                    None
                };
                let boundary = match current {
                    None => true,
                    Some(current) => previous_atom_rank != heap.slice(nRank)?[usize::from(current)],
                };
                if boundary {
                    let i2 = i;
                    if i2 - i1 > num_no_parity_atoms {
                        let first_center = heap.slice(nAtomNumber)?[i1 as usize];
                        num_no_parity_atoms = 0;
                        let first = heap.slice(at.as_const())?[usize::from(first_center)].clone();
                        let valence = i32::from(first.valence);
                        let implicit_h = if first.endpoint != 0 {
                            0
                        } else {
                            i32::from(first.num_H)
                        };
                        if valence + implicit_h > MAX_NUM_STEREO_ATOM_NEIGH as i32 {
                            i += 1;
                            continue;
                        }
                        {
                            let count = valence as usize;
                            let mut neighbor_numbers =
                                [0 as AT_RANK; MAX_NUM_STEREO_ATOM_NEIGH as usize];
                            for (index, value) in
                                neighbor_numbers.iter_mut().take(count).enumerate()
                            {
                                *value = index as AT_RANK;
                            }
                            for ic in i1..i2 {
                                let center = heap.slice(nAtomNumber)?[ic as usize];
                                let atom_snapshot =
                                    heap.slice(at.as_const())?[usize::from(center)].clone();
                                sort_neighbor_numbers(
                                    &mut neighbor_numbers,
                                    count,
                                    &atom_snapshot.neighbor,
                                    heap.slice(nRank)?,
                                )?;
                                let mut neighbors_with_parity = -1_i32;
                                let mut k = 0_usize;
                                let mut k1 = 0_usize;
                                let mut previous_neighbor_rank: AT_RANK = 0;
                                while k <= count {
                                    let rank_changed = k == count
                                        || previous_neighbor_rank
                                            != heap.slice(nRank)?[usize::from(
                                                atom_snapshot.neighbor
                                                    [usize::from(neighbor_numbers[k])],
                                            )];
                                    if rank_changed {
                                        let k2 = k;
                                        if k2 - k1 > 1 {
                                            num_in_same_ring_system = 0;
                                            let mut ring_system = 0_u16;
                                            for kn in k1..k2 {
                                                let neighbor = atom_snapshot.neighbor
                                                    [usize::from(neighbor_numbers[kn])];
                                                let current_ring = heap.slice(at.as_const())?
                                                    [usize::from(neighbor)]
                                                .nRingSystem;
                                                if ring_system == 0 {
                                                    ring_system = current_ring;
                                                } else {
                                                    num_in_same_ring_system +=
                                                        i32::from(ring_system == current_ring);
                                                }
                                            }
                                            neighbors_with_parity = 0;
                                            for kn in k1..k2 {
                                                heap.slice_mut(visited)?.fill(0);
                                                heap.slice_mut(visited)?[usize::from(center)] = 1;
                                                let neighbor = atom_snapshot.neighbor
                                                    [usize::from(neighbor_numbers[kn])];
                                                neighbors_with_parity += find_atoms_with_parity(
                                                    heap,
                                                    at.as_const(),
                                                    visited,
                                                    i32::from(center),
                                                    i32::from(neighbor),
                                                )?;
                                            }
                                        }
                                        if neighbors_with_parity == 0
                                            && num_in_same_ring_system == 0
                                        {
                                            break;
                                        }
                                        if k + 1 < count {
                                            k1 = k;
                                            previous_neighbor_rank = heap.slice(nRank)?
                                                [usize::from(
                                                    atom_snapshot.neighbor
                                                        [usize::from(neighbor_numbers[k])],
                                                )];
                                        } else {
                                            break;
                                        }
                                    }
                                    k += 1;
                                }
                                if implicit_h > 1
                                    && ((bIsotopic != 0
                                        && atom_snapshot.num_iso_H.iter().any(|&value| value > 1))
                                        || implicit_h > NUM_H_ISOTOPES as i32
                                        || bIsotopic == 0)
                                {
                                    neighbors_with_parity = 0;
                                }
                                num_no_parity_atoms += i32::from(neighbors_with_parity == 0);
                                num_with_eq_neigh_in_same_ring_system +=
                                    i32::from(num_in_same_ring_system != 0);
                            }
                            if num_no_parity_atoms == i2 - i1
                                && num_with_eq_neigh_in_same_ring_system != i2 - i1
                            {
                                for ic in i1..i2 {
                                    let center = heap.slice(nAtomNumber)?[ic as usize];
                                    {
                                        let atom = &mut heap.slice_mut(at)?[usize::from(center)];
                                        atom.parity = 0;
                                        atom.stereo_atom_parity = 0;
                                        atom.final_parity = 0;
                                        atom.bHasStereoOrEquToStereo = 0;
                                    }
                                    for slot in 0..MAX_NUM_STEREO_BOND_NEIGH as usize {
                                        let opposite = heap.slice(at.as_const())?
                                            [usize::from(center)]
                                        .stereo_bond_neighbor[slot];
                                        if opposite == 0 {
                                            break;
                                        }
                                        let opposite_index = opposite - 1;
                                        for reverse in 0..MAX_NUM_STEREO_BOND_NEIGH as usize {
                                            let reverse_neighbor = heap.slice(at.as_const())?
                                                [usize::from(opposite_index)]
                                            .stereo_bond_neighbor[reverse];
                                            if reverse_neighbor == 0 {
                                                break;
                                            }
                                            if reverse_neighbor - 1 == center {
                                                RemoveHalfStereoBond(
                                                    heap,
                                                    at,
                                                    i32::from(opposite_index),
                                                    reverse as i32,
                                                )?;
                                                break;
                                            }
                                        }
                                        let atom = &mut heap.slice_mut(at)?[usize::from(center)];
                                        atom.stereo_bond_neighbor[slot] = 0;
                                        atom.stereo_bond_ord[slot] = 0;
                                        atom.stereo_bond_z_prod[slot] = 0;
                                        atom.stereo_bond_parity[slot] = 0;
                                    }
                                }
                                num_removed_parities =
                                    num_removed_parities.wrapping_add(num_no_parity_atoms);
                            }
                        }
                    }
                    if let Some(current) = current {
                        previous_atom_rank = heap.slice(nRank)?[usize::from(current)];
                        i1 = i;
                    }
                    num_no_parity_atoms = 0;
                }
                if let Some(current) = current {
                    num_no_parity_atoms +=
                        i32::from(heap.slice(at.as_const())?[usize::from(current)].parity == 0);
                }
                i += 1;
            }
            if num_removed_parities == removed_before {
                break;
            }
        }
        Ok(num_removed_parities)
    })();
    let cleanup = heap.free(visited);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillOutStereoParities(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_RANK>,
    nRank: SourceConstPointer<AT_RANK>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
    pCG: &mut CANON_GLOBALS,
    bIsotopic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2127 FillOutStereoParities
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int FillOutStereoParities( sp_ATOM       *at,
                               int           num_atoms,
                               const AT_RANK *nCanonRank,
                               const AT_RANK *nAtomNumberCanon,
                               const AT_RANK *nRank,
                               const AT_RANK *nAtomNumber,
                               CANON_STAT    *pCS,
                               CANON_GLOBALS *pCG,
                               int           bIsotopic )
    {
        int ret;
        /*  unmark atoms with 2 or more constitutionally equivalent neighbors */
        /*  such that there is no path through them to an atom with parity */
        ret = UnmarkNonStereo( pCG, at, num_atoms, nRank, nAtomNumber, bIsotopic );
        if (ret < 0)
            return ret;  /*  program error? */ /*   <BRKPT> */
        ret = FillAllStereoDescriptors( pCG, at, num_atoms, nCanonRank, nAtomNumberCanon, pCS ); /*  ret<0: error */
        if (!ret)
        {
            ret = pCS->nLenLinearCTStereoCarb + pCS->nLenLinearCTStereoDble;
        }
        if (ret < 0)
        {
            return ret; /*  program error? */ /*   <BRKPT> */
        }

        if (ret >= 0)
        {
            int ret2;
            ret2 = SetKnownStereoCenterParities( pCG, at, num_atoms, nCanonRank, nRank, nAtomNumber );
            if (ret2 >= 0)
            {
                ret2 = MarkKnownEqualStereoCenterParities( at, num_atoms, nRank, nAtomNumber );
            }
            if (ret2 >= 0)
            {
                ret2 = SetKnownStereoBondParities( pCG, at, num_atoms, nCanonRank, nRank, nAtomNumber );
                if (ret2 >= 0)
                {
                    ret2 = MarkKnownEqualStereoBondParities( at, num_atoms, nRank, nAtomNumber );
                }
            }
    #if ( REMOVE_KNOWN_NONSTEREO == 1 ) /* { */
            if (ret2 >= 0)
            {
                int ret3;
                do
                {
                    ret2 = RemoveKnownNonStereoCenterParities( pCG, at, num_atoms, nCanonRank, nRank, pCS );
                    if (ret2 >= 0)
                    {
                        ret3 = RemoveKnownNonStereoBondParities( at, num_atoms, nCanonRank, nRank, pCS );
                        ret2 = ret3 >= 0 ? ret2 + ret3 : ret3;
                    }
                }
                while (ret2 > 0);
            }
            if (RETURNED_ERROR( ret2 ))
            {
                ret = ret2;
            }
    #endif /* } REMOVE_KNOWN_NONSTEREO */
        }

        return ret; /*  non-zero means error */
    }
    */
    // END INCHI C FUNCTION: FillOutStereoParities

    let mut ret = UnmarkNonStereo(heap, pCG, at, num_atoms, nRank, nAtomNumber, bIsotopic)?;
    if ret < 0 {
        return Ok(ret);
    }
    ret = FillAllStereoDescriptors(heap, pCG, at, num_atoms, nCanonRank, nAtomNumberCanon, pCS)?;
    if ret == 0 {
        ret = pCS
            .nLenLinearCTStereoCarb
            .wrapping_add(pCS.nLenLinearCTStereoDble);
    }
    if ret < 0 {
        return Ok(ret);
    }

    let mut ret2 = SetKnownStereoCenterParities(
        heap,
        pCG,
        at,
        num_atoms,
        nCanonRank.as_const(),
        nRank,
        nAtomNumber,
    )?;
    if ret2 >= 0 {
        ret2 = MarkKnownEqualStereoCenterParities(heap, at, num_atoms, nRank, nAtomNumber)?;
    }
    if ret2 >= 0 {
        ret2 = SetKnownStereoBondParities(
            heap,
            pCG,
            at,
            num_atoms,
            nCanonRank.as_const(),
            nRank,
            nAtomNumber,
        )?;
        if ret2 >= 0 {
            ret2 = MarkKnownEqualStereoBondParities(heap, at, num_atoms, nRank, nAtomNumber)?;
        }
    }
    if ret2 >= 0 {
        loop {
            ret2 = RemoveKnownNonStereoCenterParities(
                heap,
                pCG,
                at,
                num_atoms,
                nCanonRank.as_const(),
                nRank,
                pCS,
            )?;
            if ret2 >= 0 {
                let ret3 = RemoveKnownNonStereoBondParities(
                    heap,
                    at,
                    num_atoms,
                    nCanonRank.as_const(),
                    nRank,
                    pCS,
                )?;
                ret2 = if ret3 >= 0 {
                    ret2.wrapping_add(ret3)
                } else {
                    ret3
                };
            }
            if ret2 <= 0 {
                break;
            }
        }
    }
    if (CT_ERR_MIN..=CT_ERR_MAX).contains(&ret2) {
        ret = ret2;
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn GetStereoNeighborPos(
    heap: &SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    iAt1: i32,
    iAt2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2197 GetStereoNeighborPos
    // INCHI✔️✔️: int GetStereoNeighborPos( sp_ATOM *at,
    // INCHI✔️✔️:                           int     iAt1,
    // INCHI✔️✔️:                           int     iAt2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k1;
    // INCHI✔️✔️:     AT_RANK sNeigh = (AT_RANK) ( iAt2 + 1 );
    // INCHI✔️✔️:     AT_RANK s;
    // INCHI✔️✔️:     for (k1 = 0; k1 < MAX_NUM_STEREO_BONDS && ( s = at[iAt1].stereo_bond_neighbor[k1] ); k1++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (s == sNeigh)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return k1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return -1; /*  neighbor not found */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetStereoNeighborPos

    let atom = heap
        .slice(at.as_const())?
        .get(usize::try_from(iAt1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let sought = iAt2.wrapping_add(1) as AT_RANK;
    let mut k1 = 0_usize;
    while k1 < MAX_NUM_STEREO_BONDS as usize {
        let neighbor = atom.stereo_bond_neighbor[k1];
        if neighbor == 0 {
            break;
        }
        if neighbor == sought {
            return Ok(k1 as i32);
        }
        k1 += 1;
    }
    Ok(-1)
}

#[allow(non_snake_case)]
pub(crate) fn GetStereoBondParity(
    heap: &SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    i: i32,
    n: i32,
    nRank: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2219 GetStereoBondParity
    // INCHI✔️✔️: int GetStereoBondParity( sp_ATOM *at,
    // INCHI✔️✔️:                          int     i,
    // INCHI✔️✔️:                          int     n,
    // INCHI✔️✔️:                          AT_RANK *nRank )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k1, k2, s, parity;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at[i].stereo_bond_neighbor[0])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (k1 = 0; k1 < MAX_NUM_STEREO_BONDS && ( s = (int) at[i].stereo_bond_neighbor[k1] ); k1++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (--s == n)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 goto neigh1_found;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return -1; /*  error: not a stereo neighbor */
    // INCHI✔️✔️:     neigh1_found:
    // INCHI✔️✔️:         if (PARITY_KNOWN( at[i].stereo_bond_parity[k1] ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return PARITY_VAL( at[i].stereo_bond_parity[k1] );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         for (k2 = 0; k2 < MAX_NUM_STEREO_BONDS && ( s = (int) at[n].stereo_bond_neighbor[k2] ); k2++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (--s == i)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 goto neigh2_found;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return -1; /*  error: not a stereo neighbor */
    // INCHI✔️✔️:     neigh2_found:;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /*  error: not a stereo bond */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (ATOM_PARITY_WELL_DEF( at[i].parity ) &&
    // INCHI✔️✔️:          ATOM_PARITY_WELL_DEF( at[n].parity ) &&
    // INCHI✔️✔️:          MIN_DOT_PROD <= abs( at[i].stereo_bond_z_prod[k1] ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  bond parity can be calculated */
    // INCHI✔️✔️:         int half_parity1, half_parity2;
    // INCHI✔️✔️:         /*  check whether all neighbors are defined */
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:         half_parity1 = HalfStereoBondParity( at, i, k1, nRank );
    // INCHI✔️✔️:         half_parity2 = HalfStereoBondParity( at, n, k2, nRank );
    // INCHI✔️✔️:         if (!half_parity1 || !half_parity2)
    // INCHI✔️✔️:             return 0; /*  ranks undefined or not a stereo bond */
    // INCHI✔️✔️:         if (ATOM_PARITY_WELL_DEF( half_parity1 ) &&
    // INCHI✔️✔️:              ATOM_PARITY_WELL_DEF( half_parity2 ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             parity = 2 - ( half_parity1 + half_parity2
    // INCHI✔️✔️:                      + ( at[i].stereo_bond_z_prod[k1] < 0 ) ) % 2;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  parity cannot be calculated: not enough info or 'unknown' */
    // INCHI✔️✔️:         if (AB_PARITY_NONE != ( parity = inchi_max( at[i].parity, at[n].parity ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             parity = AB_PARITY_UNDF; /*  should not happen */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return parity;
    // INCHI✔️✔️: }
    // INCHI✔️✔️:
    // END INCHI C FUNCTION: GetStereoBondParity

    let atoms = heap.slice(at.as_const())?;
    let atom_i = atoms
        .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom_i.stereo_bond_neighbor[0] == 0 {
        return Ok(-1);
    }

    let mut k1 = 0_usize;
    let mut found1 = false;
    while k1 < MAX_NUM_STEREO_BONDS as usize {
        let s = i32::from(atom_i.stereo_bond_neighbor[k1]);
        if s == 0 {
            break;
        }
        if s.wrapping_sub(1) == n {
            found1 = true;
            break;
        }
        k1 += 1;
    }
    if !found1 {
        return Ok(-1);
    }

    let source_parity = i32::from(atom_i.stereo_bond_parity[k1]);
    if parity_known(source_parity) {
        return Ok(source_parity & BITS_PARITY);
    }

    let atom_n = atoms
        .get(usize::try_from(n).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut k2 = 0_usize;
    let mut found2 = false;
    while k2 < MAX_NUM_STEREO_BONDS as usize {
        let s = i32::from(atom_n.stereo_bond_neighbor[k2]);
        if s == 0 {
            break;
        }
        if s.wrapping_sub(1) == i {
            found2 = true;
            break;
        }
        k2 += 1;
    }
    if !found2 {
        return Ok(-1);
    }

    let parity = if atom_parity_well_defined(i32::from(atom_i.parity))
        && atom_parity_well_defined(i32::from(atom_n.parity))
        && MIN_DOT_PROD as i32 <= i32::from(atom_i.stereo_bond_z_prod[k1]).wrapping_abs()
    {
        let half_parity1 = HalfStereoBondParity(heap, at, i, k1 as i32, nRank)?;
        let half_parity2 = HalfStereoBondParity(heap, at, n, k2 as i32, nRank)?;
        if half_parity1 == 0 || half_parity2 == 0 {
            return Ok(0);
        }
        if atom_parity_well_defined(half_parity1) && atom_parity_well_defined(half_parity2) {
            2_i32.wrapping_sub(
                half_parity1
                    .wrapping_add(half_parity2)
                    .wrapping_add(i32::from(atom_i.stereo_bond_z_prod[k1] < 0))
                    .wrapping_rem(2),
            )
        } else {
            return Ok(CT_STEREOBOND_ERROR);
        }
    } else {
        let maximum = i32::from(atom_i.parity).max(i32::from(atom_n.parity));
        if maximum != AB_PARITY_NONE as i32 {
            AB_PARITY_UNDF as i32
        } else {
            maximum
        }
    };
    Ok(parity)
}

#[allow(non_snake_case)]
pub(crate) fn GetPermutationParity(
    heap: &SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    nAvoidNeighbor: AT_RANK,
    nCanonRank: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2296 GetPermutationParity
    // INCHI✔️✔️: int GetPermutationParity( CANON_GLOBALS *pCG,
    // INCHI✔️✔️:                           sp_ATOM       *at,
    // INCHI✔️✔️:                           AT_RANK       nAvoidNeighbor,
    // INCHI✔️✔️:                           AT_RANK       *nCanonRank )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_RANK nNeighRank[MAX_NUM_STEREO_ATOM_NEIGH];
    // INCHI✔️✔️:     int     j, k, parity;
    // INCHI✔️✔️:     if (at->valence > MAX_NUM_STEREO_ATOM_NEIGH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         parity = -1; /*  error */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (j = k = 0; j < at->valence; j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at->neighbor[j] != nAvoidNeighbor)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nNeighRank[k++] = nCanonRank[(int) at->neighbor[j]];
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (k)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             parity = insertions_sort( pCG, nNeighRank, k, sizeof( nNeighRank[0] ), comp_AT_RANK );
    // INCHI✔️✔️:             if (nNeighRank[0])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 parity = 2 - parity % 2;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 parity = 0; /*  not all ranks are known */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* special case: HX= with implicit H */
    // INCHI✔️✔️:             parity = 2;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return parity;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetPermutationParity

    let atom = heap
        .slice(at.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if i32::from(atom.valence) > MAX_NUM_STEREO_ATOM_NEIGH as i32 {
        return Ok(-1);
    }
    let ranks = heap.slice(nCanonRank.as_const())?;
    let mut neighbor_ranks = [0_u16; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut j = 0_i32;
    let mut k = 0_usize;
    while j < i32::from(atom.valence) {
        let neighbor = atom.neighbor[j as usize];
        if neighbor != nAvoidNeighbor {
            neighbor_ranks[k] = *ranks
                .get(usize::from(neighbor))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            k += 1;
        }
        j = j.wrapping_add(1);
    }
    let parity = if k != 0 {
        let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(&mut neighbor_ranks);
        let exchanges = insertions_sort(
            bytes,
            k,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                Ok(comp_AT_RANK(
                    AT_RANK::from_ne_bytes([left[0], left[1]]),
                    AT_RANK::from_ne_bytes([right[0], right[1]]),
                ))
            },
        )?;
        if neighbor_ranks[0] != 0 {
            2_i32.wrapping_sub(exchanges.wrapping_rem(2))
        } else {
            0
        }
    } else {
        2
    };
    let _ = pCG;
    Ok(parity)
}

#[allow(non_snake_case)]
pub(crate) fn GetStereoCenterParity(
    heap: &SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    i: i32,
    nRank: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2340 GetStereoCenterParity
    // INCHI✔️✔️: int GetStereoCenterParity( CANON_GLOBALS *pCG,
    // INCHI✔️✔️:                            sp_ATOM       *at,
    // INCHI✔️✔️:                            int           i,
    // INCHI✔️✔️:                            AT_RANK       *nRank )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_NUMB  nNeighborNumber2[MAXVAL];
    // INCHI✔️✔️:     int      parity;
    // INCHI✔️✔️:     int      k, num_trans;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!at[i].parity)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;   /*  not a stereo center                     */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (at[i].stereo_bond_neighbor[0])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;  /*  a stereo bond atom, not a stereo center */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (ATOM_PARITY_WELL_DEF( at[i].parity ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  number of neighbors transpositions to the sorted order is unknown. Find it. */
    // INCHI✔️✔️:         /*  If parity is not well-defined then doing this is a waste of time            */
    // INCHI✔️✔️:         int num_neigh = at[i].valence;
    // INCHI✔️✔️:         for (k = 0; k < num_neigh; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (!nRank[(int) at[i].neighbor[k]])
    // INCHI✔️✔️:                 return 0; /*  stereo at[i] does not belong to the traversed part of the structure */
    // INCHI✔️✔️:             nNeighborNumber2[k] = k;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         pCG->m_pNeighborsForSort = at[i].neighbor;
    // INCHI✔️✔️:         pCG->m_pn_RankForSort = nRank;
    // INCHI✔️✔️:         num_trans = insertions_sort( pCG, nNeighborNumber2, num_neigh, sizeof( nNeighborNumber2[0] ), CompNeighborsAT_NUMBER );
    // INCHI✔️✔️: #ifndef CT_NEIGH_INCREASE
    // INCHI✔️✔️:         num_trans += ( ( num_neigh*( num_neigh - 1 ) ) / 2 ) % 2;  /*  get correct parity for ascending order */
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         parity = 2 - ( at[i].parity + num_trans ) % 2;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         parity = at[i].parity;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return parity;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetStereoCenterParity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:850
    // INCHI✔️✔️: #define CT_NEIGH_INCREASE
    // END INCHI ACTIVE MACRO CONFIGURATION

    let atom = heap
        .slice(at.as_const())?
        .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.parity == 0 {
        return Ok(0);
    }
    if atom.stereo_bond_neighbor[0] != 0 {
        return Ok(-1);
    }
    if !atom_parity_well_defined(i32::from(atom.parity)) {
        return Ok(i32::from(atom.parity));
    }

    let count =
        usize::try_from(atom.valence).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if count > MAXVAL as usize {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let ranks = heap.slice(nRank.as_const())?;
    let mut neighbor_numbers = [0_u16; MAXVAL as usize];
    for (k, slot) in neighbor_numbers.iter_mut().take(count).enumerate() {
        let neighbor = atom.neighbor[k];
        if *ranks
            .get(usize::from(neighbor))
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            == 0
        {
            return Ok(0);
        }
        *slot = k as AT_NUMB;
    }
    pCG.m_pn_RankForSort = nRank.as_const();
    let bytes = bytemuck::cast_slice_mut::<AT_NUMB, u8>(&mut neighbor_numbers);
    let num_trans = insertions_sort(
        bytes,
        count,
        std::mem::size_of::<AT_NUMB>(),
        &mut |left, right| {
            CompNeighborsAT_NUMBER(
                AT_NUMB::from_ne_bytes([left[0], left[1]]),
                AT_NUMB::from_ne_bytes([right[0], right[1]]),
                CompNeighborsATNumberContext::Slices {
                    neighbors: &atom.neighbor,
                    ranks,
                },
            )
        },
    )?;
    Ok(2_i32.wrapping_sub(
        i32::from(atom.parity)
            .wrapping_add(num_trans)
            .wrapping_rem(2),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichicans__removeonestereocenter__line_300() {
        let mut heap = SourceHeap::default();
        let mut marked = sp_ATOM::default();
        marked.parity = 2;
        marked.stereo_atom_parity = 3;
        marked.final_parity = 4;
        let mut unmarked = sp_ATOM::default();
        unmarked.stereo_atom_parity = 5;
        unmarked.final_parity = 6;
        let atoms = heap.allocate_model_storage(vec![marked, unmarked]).unwrap();

        assert_eq!(RemoveOneStereoCenter(&mut heap, atoms, 0), Ok(1));
        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            (
                values[0].parity,
                values[0].stereo_atom_parity,
                values[0].final_parity
            ),
            (0, 0, 0)
        );
        assert_eq!(RemoveOneStereoCenter(&mut heap, atoms, 1), Ok(0));
        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            (
                values[1].parity,
                values[1].stereo_atom_parity,
                values[1].final_parity
            ),
            (0, 5, 6)
        );
        assert_eq!(
            RemoveOneStereoCenter(&mut heap, atoms, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }
    use crate::source_types::{AB_PARITY_ODD, AT_STEREO_CARB, AT_STEREO_DBLE};

    fn fixture(
        atoms: Vec<sp_ATOM>,
        carb: Vec<AT_STEREO_CARB>,
        double: Vec<AT_STEREO_DBLE>,
    ) -> (
        SourceHeap,
        SourceMutPointer<sp_ATOM>,
        SourceMutPointer<AT_RANK>,
        SourceMutPointer<AT_RANK>,
        CANON_STAT,
    ) {
        let mut heap = SourceHeap::default();
        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let ranks = heap.allocate_model_storage(vec![2, 1, 3, 4]).unwrap();
        let canonical = heap.allocate_model_storage(vec![99; 4]).unwrap();
        let carb_pointer = heap.allocate_model_storage(carb).unwrap();
        let double_pointer = heap.allocate_model_storage(double).unwrap();
        let stat = CANON_STAT {
            LinearCTStereoCarb: carb_pointer,
            nLenLinearCTStereoCarb: heap.slice(carb_pointer.as_const()).unwrap().len() as i32,
            LinearCTStereoDble: double_pointer,
            nLenLinearCTStereoDble: heap.slice(double_pointer.as_const()).unwrap().len() as i32,
            ..CANON_STAT::default()
        };
        (heap, atom_pointer, ranks, canonical, stat)
    }

    fn atom(parity: i8) -> sp_ATOM {
        sp_ATOM {
            parity,
            ..sp_ATOM::default()
        }
    }

    fn stereo_pair(parity: i8) -> Vec<sp_ATOM> {
        let mut atoms = vec![atom(1), atom(1)];
        atoms[0].valence = 1;
        atoms[0].neighbor[0] = 1;
        atoms[0].stereo_bond_neighbor[0] = 2;
        atoms[0].stereo_bond_parity[0] = parity;
        atoms[1].valence = 1;
        atoms[1].neighbor[0] = 0;
        atoms[1].stereo_bond_neighbor[0] = 1;
        atoms[1].stereo_bond_parity[0] = parity;
        atoms
    }

    fn run_mark_known(
        atoms: Vec<sp_ATOM>,
        ranks: Vec<AT_RANK>,
        atom_numbers: Vec<AT_RANK>,
    ) -> (i32, Vec<sp_ATOM>) {
        let mut heap = SourceHeap::default();
        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let rank_pointer = heap.allocate_model_storage(ranks).unwrap();
        let number_pointer = heap.allocate_model_storage(atom_numbers).unwrap();
        let count = heap.slice(atom_pointer.as_const()).unwrap().len() as i32;
        let result = MarkKnownEqualStereoBondParities(
            &mut heap,
            atom_pointer,
            count,
            rank_pointer.as_const(),
            number_pointer.as_const(),
        )
        .unwrap();
        (
            result,
            heap.slice(atom_pointer.as_const()).unwrap().to_vec(),
        )
    }

    #[test]
    fn source_port__ichicans__fillsinglestereodescriptors__line_525() {
        #[allow(clippy::too_many_arguments)]
        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            atom_index: i32,
            num_trans: i32,
            with_carb: bool,
            mut carb_len: i32,
            carb_max: i32,
            with_double: bool,
            mut double_len: i32,
            double_max: i32,
            allene: i32,
        ) -> (
            Result<i32, SourceHeapError>,
            Vec<AT_STEREO_CARB>,
            Vec<AT_STEREO_DBLE>,
            i32,
            i32,
            CANON_GLOBALS,
        ) {
            let mut heap = SourceHeap::default();
            let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
            let rank_pointer = heap.allocate_model_storage(ranks).unwrap();
            let carb_pointer = if with_carb {
                heap.allocate_model_storage(vec![
                    AT_STEREO_CARB {
                        at_num: 91,
                        parity: 92
                    };
                    3
                ])
                .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let double_pointer = if with_double {
                heap.allocate_model_storage(vec![
                    AT_STEREO_DBLE {
                        at_num1: 81,
                        at_num2: 82,
                        parity: 83,
                    };
                    3
                ])
                .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let mut globals = CANON_GLOBALS::default();
            let result = FillSingleStereoDescriptors(
                &mut heap,
                &mut globals,
                atom_pointer,
                atom_index,
                num_trans,
                rank_pointer,
                carb_pointer,
                &mut carb_len,
                carb_max,
                double_pointer,
                &mut double_len,
                double_max,
                allene,
            );
            let carb = if with_carb {
                heap.slice(carb_pointer.as_const()).unwrap().to_vec()
            } else {
                Vec::new()
            };
            let double = if with_double {
                heap.slice(double_pointer.as_const()).unwrap().to_vec()
            } else {
                Vec::new()
            };
            (result, carb, double, carb_len, double_len, globals)
        }

        fn center(parity: i8) -> Vec<sp_ATOM> {
            let mut atoms = vec![sp_ATOM::default(); 4];
            atoms[0].parity = parity;
            atoms[0].valence = 3;
            atoms[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
            atoms
        }

        fn stereo_bond(encoded_parity: i8, z_product: i8) -> Vec<sp_ATOM> {
            let mut atoms = center(1);
            atoms[0].stereo_bond_neighbor[0] = 2;
            atoms[0].stereo_bond_ord[0] = 0;
            atoms[0].stereo_bond_parity[0] = encoded_parity;
            atoms[0].stereo_bond_z_prod[0] = z_product;
            atoms[1].parity = 2;
            atoms[1].valence = 3;
            atoms[1].neighbor[..3].copy_from_slice(&[0, 2, 3]);
            atoms[1].stereo_bond_neighbor[0] = 1;
            atoms[1].stereo_bond_ord[0] = 0;
            atoms[1].stereo_bond_parity[0] = encoded_parity;
            atoms[1].stereo_bond_z_prod[0] = z_product;
            atoms
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty_globals = CANON_GLOBALS::default();
        let mut carb_len = 7;
        let mut double_len = 8;
        assert_eq!(
            FillSingleStereoDescriptors(
                &mut empty_heap,
                &mut empty_globals,
                SourceMutPointer::null(),
                i32::MIN,
                -1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut carb_len,
                0,
                SourceMutPointer::null(),
                &mut double_len,
                0,
                0,
            ),
            Ok(0)
        );
        assert_eq!((carb_len, double_len), (7, 8));

        let result = run(
            center(1),
            vec![4, 3, 1, 2],
            0,
            -1,
            true,
            0,
            3,
            false,
            0,
            0,
            0,
        );
        assert_eq!(result.0, Ok(0));
        assert_eq!(
            result.1[0],
            AT_STEREO_CARB {
                at_num: 4,
                parity: 1
            }
        );
        assert_eq!(result.3, 1);
        assert!(!result.5.m_pn_RankForSort.is_null());

        let manual = run(
            center(1),
            vec![4, 3, 1, 2],
            0,
            1,
            true,
            0,
            3,
            false,
            0,
            0,
            0,
        );
        assert_eq!(manual.1[0].parity, 2);
        assert!(manual.5.m_pn_RankForSort.is_null());

        let ill_defined = run(
            center(3),
            vec![4, 3, 1, 2],
            0,
            -99,
            true,
            0,
            3,
            false,
            0,
            0,
            0,
        );
        assert_eq!(ill_defined.1[0].parity, 3);

        let overflow = run(
            center(1),
            vec![4, 3, 1, 2],
            0,
            0,
            true,
            0,
            0,
            false,
            0,
            0,
            0,
        );
        assert_eq!(overflow.0, Ok(CT_OVERFLOW));
        assert_eq!(overflow.1[0].at_num, 91);
        assert_eq!(overflow.3, 0);

        let allene_center = run(
            center(1),
            vec![4, 3, 1, 2],
            0,
            0,
            true,
            0,
            3,
            false,
            0,
            0,
            1,
        );
        assert_eq!(allene_center.3, 0);
        assert_eq!(allene_center.1[0].at_num, 91);

        let known = run(
            stereo_bond(2, 0),
            vec![4, 3, 2, 1],
            0,
            -1,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(known.0, Ok(0));
        assert_eq!(
            known.2[0],
            AT_STEREO_DBLE {
                at_num1: 4,
                at_num2: 3,
                parity: 2,
            }
        );
        assert_eq!(known.4, 1);

        let higher_rank = run(
            stereo_bond(2, 0),
            vec![2, 3, 1, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(higher_rank.4, 0);
        assert_eq!(higher_rank.2[0].at_num1, 81);

        let no_parity = run(
            stereo_bond(0, 0),
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(no_parity.4, 0);

        let fallback = run(
            stereo_bond(AB_PARITY_CALC as i8, 49),
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(fallback.2[0].parity, AB_PARITY_UNDF as u8);

        let mut no_fallback_atoms = stereo_bond(AB_PARITY_CALC as i8, 49);
        no_fallback_atoms[0].parity = 0;
        no_fallback_atoms[1].parity = 0;
        let no_fallback = run(
            no_fallback_atoms,
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(no_fallback.4, 0);

        let calculated = run(
            stereo_bond(AB_PARITY_CALC as i8, 50),
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(calculated.0, Ok(0));
        assert_eq!(calculated.2[0].parity, 1);

        let calculated_negative = run(
            stereo_bond(AB_PARITY_CALC as i8, -50),
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(calculated_negative.2[0].parity, 2);

        let mut missing_reverse = stereo_bond(AB_PARITY_CALC as i8, 50);
        missing_reverse[1].stereo_bond_neighbor[0] = 0;
        let missing_reverse = run(
            missing_reverse,
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(missing_reverse.0, Ok(CT_STEREOBOND_ERROR));
        assert_eq!(missing_reverse.4, 0);

        let tied_half = run(
            stereo_bond(AB_PARITY_CALC as i8, 50),
            vec![4, 3, 2, 2],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(tied_half.0, Ok(CT_STEREOBOND_ERROR));

        let double_overflow = run(
            stereo_bond(2, 0),
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            0,
            0,
        );
        assert_eq!(double_overflow.0, Ok(CT_OVERFLOW));
        assert_eq!(double_overflow.2[0].at_num1, 81);
        assert_eq!(double_overflow.4, 0);

        let allene_atoms = stereo_bond((MULT_STEREOBOND + 1) as i8, 0);
        let skipped_allene = run(
            allene_atoms.clone(),
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            0,
        );
        assert_eq!(skipped_allene.4, 0);
        let selected_allene = run(
            allene_atoms,
            vec![4, 3, 2, 1],
            0,
            0,
            false,
            0,
            0,
            true,
            0,
            3,
            1,
        );
        assert_eq!(selected_allene.4, 1);
        assert_eq!(selected_allene.2[0].parity, 1);
    }

    #[test]
    fn source_port__ichicans__setcttononisotopicstereo__line_746() {
        let mut heap = SourceHeap::default();
        let double = heap
            .allocate_model_storage(vec![AT_STEREO_DBLE::default()])
            .unwrap();
        let carb = heap
            .allocate_model_storage(vec![AT_STEREO_CARB::default()])
            .unwrap();
        let double_inv = heap
            .allocate_model_storage(vec![AT_STEREO_DBLE::default(); 2])
            .unwrap();
        let carb_inv = heap
            .allocate_model_storage(vec![AT_STEREO_CARB::default(); 2])
            .unwrap();
        let source = CANON_STAT {
            LinearCTStereoDble: double,
            LinearCTStereoCarb: carb,
            LinearCTStereoDbleInv: double_inv,
            LinearCTStereoCarbInv: carb_inv,
            nMaxLenLinearCTStereoDble: 111,
            nMaxLenLinearCTStereoCarb: 112,
            nLenLinearCTStereoDble: 113,
            nLenLinearCTStereoCarb: 114,
            nLenLinearCTIsotopicStereoDble: 115,
            nLenLinearCTIsotopicStereoCarb: 116,
            ..CANON_STAT::default()
        };
        let mut target = CANON_STAT {
            nMaxLenLinearCTStereoDble: -1,
            nMaxLenLinearCTStereoCarb: -2,
            nLenLinearCTStereoDble: -3,
            nLenLinearCTStereoCarb: -4,
            nLenLinearCTIsotopicStereoDble: -5,
            nLenLinearCTIsotopicStereoCarb: -6,
            bFirstCT: 211,
            ..CANON_STAT::default()
        };

        SetCtToNonIsotopicStereo(&mut target, &source);

        assert_eq!(target.LinearCTStereoDble, double);
        assert_eq!(target.LinearCTStereoCarb, carb);
        assert_eq!(target.LinearCTStereoDbleInv, double_inv);
        assert_eq!(target.LinearCTStereoCarbInv, carb_inv);
        assert_eq!(target.nMaxLenLinearCTStereoDble, 111);
        assert_eq!(target.nMaxLenLinearCTStereoCarb, 112);
        assert_eq!(target.nLenLinearCTStereoDble, 113);
        assert_eq!(target.nLenLinearCTStereoCarb, 114);
        assert_eq!(target.nLenLinearCTIsotopicStereoDble, 115);
        assert_eq!(target.nLenLinearCTIsotopicStereoCarb, 116);
        assert_eq!(target.bFirstCT, 211);
    }

    #[test]
    fn source_port__ichicans__setcttoisotopicstereo__line_729() {
        let mut heap = SourceHeap::default();
        let isotopic_double = heap
            .allocate_model_storage(vec![AT_STEREO_DBLE::default()])
            .unwrap();
        let isotopic_carb = heap
            .allocate_model_storage(vec![AT_STEREO_CARB::default()])
            .unwrap();
        let isotopic_double_inv = heap
            .allocate_model_storage(vec![AT_STEREO_DBLE::default(); 2])
            .unwrap();
        let isotopic_carb_inv = heap
            .allocate_model_storage(vec![AT_STEREO_CARB::default(); 2])
            .unwrap();
        let source = CANON_STAT {
            LinearCTIsotopicStereoDble: isotopic_double,
            LinearCTIsotopicStereoCarb: isotopic_carb,
            LinearCTIsotopicStereoDbleInv: isotopic_double_inv,
            LinearCTIsotopicStereoCarbInv: isotopic_carb_inv,
            nMaxLenLinearCTIsotopicStereoDble: 101,
            nMaxLenLinearCTIsotopicStereoCarb: 102,
            nLenLinearCTIsotopicStereoDble: 103,
            nLenLinearCTIsotopicStereoCarb: 104,
            ..CANON_STAT::default()
        };
        let mut target = CANON_STAT {
            nLenLinearCTStereoDble: -1,
            nLenLinearCTStereoCarb: -2,
            nMaxLenLinearCTStereoDble: -3,
            nMaxLenLinearCTStereoCarb: -4,
            nLenLinearCTIsotopicStereoDble: 201,
            nLenLinearCTIsotopicStereoCarb: 202,
            bFirstCT: 203,
            ..CANON_STAT::default()
        };

        SetCtToIsotopicStereo(&mut target, &source);

        assert_eq!(target.LinearCTStereoDble, isotopic_double);
        assert_eq!(target.LinearCTStereoCarb, isotopic_carb);
        assert_eq!(target.LinearCTStereoDbleInv, isotopic_double_inv);
        assert_eq!(target.LinearCTStereoCarbInv, isotopic_carb_inv);
        assert_eq!(target.nMaxLenLinearCTStereoDble, 101);
        assert_eq!(target.nMaxLenLinearCTStereoCarb, 102);
        assert_eq!(target.nLenLinearCTStereoDble, 103);
        assert_eq!(target.nLenLinearCTStereoCarb, 104);
        assert_eq!(target.nLenLinearCTIsotopicStereoDble, 201);
        assert_eq!(target.nLenLinearCTIsotopicStereoCarb, 202);
        assert_eq!(target.bFirstCT, 203);
    }

    #[test]
    fn source_port__ichicans__switchatomstereoandisotopicstereo__line_705() {
        let mut first = sp_ATOM {
            parity: 1,
            parity2: 2,
            final_parity: 3,
            final_parity2: 4,
            stereo_atom_parity: 5,
            stereo_atom_parity2: 6,
            bHasStereoOrEquToStereo: 7,
            bHasStereoOrEquToStereo2: 8,
            ..sp_ATOM::default()
        };
        first.stereo_bond_neighbor = [11, 12, 13];
        first.stereo_bond_neighbor2 = [21, 22, 23];
        first.stereo_bond_ord = [31, 32, 33];
        first.stereo_bond_ord2 = [41, 42, 43];
        first.stereo_bond_z_prod = [51, 52, 53];
        first.stereo_bond_z_prod2 = [61, 62, 63];
        first.stereo_bond_parity = [71, 72, 73];
        first.stereo_bond_parity2 = [81, 82, 83];
        let second = sp_ATOM {
            parity: -1,
            parity2: -2,
            final_parity: -3,
            final_parity2: -4,
            ..sp_ATOM::default()
        };
        let untouched = sp_ATOM {
            parity: 91,
            parity2: 92,
            ..sp_ATOM::default()
        };
        let original = vec![first, second, untouched];
        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(original.clone()).unwrap();
        let mut switched = 0;
        assert_eq!(
            SwitchAtomStereoAndIsotopicStereo(&mut heap, atoms, 2, &mut switched),
            Ok(())
        );
        assert_eq!(switched, 1);
        let current = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            (
                current[0].parity,
                current[0].parity2,
                current[0].final_parity,
                current[0].final_parity2,
                current[0].stereo_atom_parity,
                current[0].stereo_atom_parity2,
                current[0].bHasStereoOrEquToStereo,
                current[0].bHasStereoOrEquToStereo2,
            ),
            (2, 1, 4, 3, 6, 5, 8, 7)
        );
        assert_eq!(current[0].stereo_bond_neighbor, [21, 22, 23]);
        assert_eq!(current[0].stereo_bond_neighbor2, [11, 12, 13]);
        assert_eq!(current[0].stereo_bond_ord, [41, 42, 43]);
        assert_eq!(current[0].stereo_bond_ord2, [31, 32, 33]);
        assert_eq!(current[0].stereo_bond_z_prod, [61, 62, 63]);
        assert_eq!(current[0].stereo_bond_z_prod2, [51, 52, 53]);
        assert_eq!(current[0].stereo_bond_parity, [81, 82, 83]);
        assert_eq!(current[0].stereo_bond_parity2, [71, 72, 73]);
        assert_eq!(current[1].parity, -2);
        assert_eq!(current[1].parity2, -1);
        assert_eq!(current[1].final_parity, -4);
        assert_eq!(current[1].final_parity2, -3);
        assert_eq!(current[2], original[2]);

        assert_eq!(
            SwitchAtomStereoAndIsotopicStereo(&mut heap, atoms, 2, &mut switched),
            Ok(())
        );
        assert_eq!(switched, 0);
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), original.as_slice());

        switched = -7;
        assert_eq!(
            SwitchAtomStereoAndIsotopicStereo(&mut heap, atoms, 0, &mut switched),
            Ok(())
        );
        assert_eq!(switched, 0);
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), original.as_slice());
    }

    #[test]
    fn source_port__ichicans__fillallstereodescriptors__line_766() {
        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            canonical_order: Vec<AT_RANK>,
            carb_max: i32,
            double_max: i32,
        ) -> (
            Result<i32, SourceHeapError>,
            CANON_STAT,
            Vec<AT_STEREO_CARB>,
            Vec<AT_STEREO_DBLE>,
        ) {
            let mut heap = SourceHeap::default();
            let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
            let rank_pointer = heap.allocate_model_storage(ranks).unwrap();
            let order_pointer = heap.allocate_model_storage(canonical_order).unwrap();
            let carb_pointer = heap
                .allocate_model_storage(vec![AT_STEREO_CARB::default(); 4])
                .unwrap();
            let double_pointer = heap
                .allocate_model_storage(vec![AT_STEREO_DBLE::default(); 4])
                .unwrap();
            let mut stat = CANON_STAT {
                LinearCTStereoCarb: carb_pointer,
                LinearCTStereoDble: double_pointer,
                nLenLinearCTStereoCarb: 91,
                nLenLinearCTStereoDble: 92,
                nMaxLenLinearCTStereoCarb: carb_max,
                nMaxLenLinearCTStereoDble: double_max,
                ..CANON_STAT::default()
            };
            let mut globals = CANON_GLOBALS::default();
            let atom_count = heap.slice(order_pointer.as_const()).unwrap().len() as i32;
            let result = FillAllStereoDescriptors(
                &mut heap,
                &mut globals,
                atom_pointer,
                atom_count,
                rank_pointer,
                order_pointer,
                &mut stat,
            );
            (
                result,
                stat,
                heap.slice(carb_pointer.as_const()).unwrap().to_vec(),
                heap.slice(double_pointer.as_const()).unwrap().to_vec(),
            )
        }

        let mut atoms = vec![sp_ATOM::default(); 5];
        atoms[0].parity = 1;
        atoms[0].valence = 3;
        atoms[0].neighbor[..3].copy_from_slice(&[2, 3, 4]);

        atoms[1].parity = 1;
        atoms[1].valence = 3;
        atoms[1].neighbor[..3].copy_from_slice(&[2, 3, 4]);
        atoms[1].stereo_bond_neighbor[0] = 3;
        atoms[1].stereo_bond_ord[0] = 0;
        atoms[1].stereo_bond_parity[0] = (MULT_STEREOBOND + 1) as i8;
        atoms[2].parity = 2;
        atoms[2].valence = 3;
        atoms[2].neighbor[..3].copy_from_slice(&[1, 3, 4]);
        atoms[2].stereo_bond_neighbor[0] = 2;
        atoms[2].stereo_bond_ord[0] = 0;
        atoms[2].stereo_bond_parity[0] = (MULT_STEREOBOND + 1) as i8;

        let completed = run(atoms.clone(), vec![2, 4, 3, 1, 1], vec![1, 0, 2], 4, 4);
        assert_eq!(completed.0, Ok(0));
        assert_eq!(completed.1.nLenLinearCTStereoCarb, 1);
        assert_eq!(completed.1.nLenLinearCTStereoDble, 1);
        assert_eq!(completed.2[0].at_num, 2);
        assert_eq!(completed.3[0].at_num1, 4);
        assert_eq!(completed.3[0].at_num2, 3);
        assert_eq!(completed.3[0].parity, 1);

        let reversed = run(atoms.clone(), vec![2, 4, 3, 1, 1], vec![0, 2, 1], 4, 4);
        assert_eq!(reversed.0, Ok(0));
        assert_eq!(reversed.1.nLenLinearCTStereoCarb, 1);
        assert_eq!(reversed.1.nLenLinearCTStereoDble, 1);

        let stopped = run(atoms, vec![2, 4, 3, 1, 1], vec![0, 1, 2], 0, 4);
        assert_eq!(stopped.0, Ok(CT_OVERFLOW));
        assert_eq!(stopped.1.nLenLinearCTStereoCarb, 0);
        assert_eq!(stopped.1.nLenLinearCTStereoDble, 0);

        let mut heap = SourceHeap::default();
        let mut globals = CANON_GLOBALS::default();
        let mut empty_stat = CANON_STAT {
            nLenLinearCTStereoCarb: 17,
            nLenLinearCTStereoDble: 18,
            ..CANON_STAT::default()
        };
        assert_eq!(
            FillAllStereoDescriptors(
                &mut heap,
                &mut globals,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty_stat,
            ),
            Ok(0)
        );
        assert_eq!(
            (
                empty_stat.nLenLinearCTStereoCarb,
                empty_stat.nLenLinearCTStereoDble
            ),
            (0, 0)
        );
    }

    #[test]
    fn source_port__ichicans__setknownstereobondparities__line_801() {
        fn direct(encoded_parity: i8, z_product: i8) -> Vec<sp_ATOM> {
            let mut atoms = vec![sp_ATOM::default(); 4];
            atoms[0].parity = 2;
            atoms[0].valence = 3;
            atoms[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
            atoms[0].stereo_bond_neighbor[0] = 2;
            atoms[0].stereo_bond_ord[0] = 0;
            atoms[0].stereo_bond_parity[0] = encoded_parity;
            atoms[0].stereo_bond_z_prod[0] = z_product;
            atoms[1].parity = 1;
            atoms[1].valence = 3;
            atoms[1].neighbor[..3].copy_from_slice(&[0, 2, 3]);
            atoms[1].stereo_bond_neighbor[0] = 1;
            atoms[1].stereo_bond_ord[0] = 0;
            atoms[1].stereo_bond_parity[0] = encoded_parity;
            atoms[1].stereo_bond_z_prod[0] = z_product;
            atoms
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            canonical_ranks: Vec<AT_RANK>,
            atom_numbers: Vec<AT_RANK>,
        ) -> (i32, Vec<sp_ATOM>) {
            let mut heap = SourceHeap::default();
            let atom_count = atoms.len() as i32;
            let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
            let rank_pointer = heap.allocate_model_storage(ranks).unwrap();
            let canonical_pointer = heap.allocate_model_storage(canonical_ranks).unwrap();
            let number_pointer = heap.allocate_model_storage(atom_numbers).unwrap();
            let mut globals = CANON_GLOBALS::default();
            let result = SetKnownStereoBondParities(
                &mut heap,
                &mut globals,
                atom_pointer,
                atom_count,
                canonical_pointer.as_const(),
                rank_pointer.as_const(),
                number_pointer.as_const(),
            )
            .unwrap();
            (
                result,
                heap.slice(atom_pointer.as_const()).unwrap().to_vec(),
            )
        }

        let ranks = vec![1, 2, 3, 4];
        let canonical = vec![1, 2, 3, 4];
        let numbers = vec![0, 1, 2, 3];

        let (count, positive) = run(
            direct(AB_PARITY_CALC as i8, 50),
            ranks.clone(),
            canonical.clone(),
            numbers.clone(),
        );
        assert_eq!(count, 1);
        assert_eq!(positive[0].stereo_bond_parity[0], 1);
        assert_eq!(positive[1].stereo_bond_parity[0], 1);

        let mut transposed = direct(AB_PARITY_CALC as i8, 50);
        transposed[1].neighbor[1..3].swap(0, 1);
        let (count, transposed) = run(
            transposed,
            ranks.clone(),
            canonical.clone(),
            numbers.clone(),
        );
        assert_eq!(count, 1);
        assert_eq!(transposed[0].stereo_bond_parity[0], 2);
        assert_eq!(transposed[1].stereo_bond_parity[0], 2);

        let (count, negative) = run(
            direct(AB_PARITY_CALC as i8, -50),
            ranks.clone(),
            canonical.clone(),
            numbers.clone(),
        );
        assert_eq!(count, 1);
        assert_eq!(negative[0].stereo_bond_parity[0], 2);
        assert_eq!(negative[1].stereo_bond_parity[0], 2);

        let (count, undefined) = run(
            direct(AB_PARITY_CALC as i8, 49),
            ranks.clone(),
            canonical.clone(),
            numbers.clone(),
        );
        assert_eq!(count, 1);
        assert_eq!(undefined[0].stereo_bond_parity[0], AB_PARITY_UNDF as i8);
        assert_eq!(undefined[1].stereo_bond_parity[0], AB_PARITY_UNDF as i8);

        let (count, cumulene_bits) = run(
            direct((MULT_STEREOBOND + AB_PARITY_CALC as i32) as i8, 50),
            ranks.clone(),
            canonical.clone(),
            numbers.clone(),
        );
        assert_eq!(count, 0);
        assert_eq!(
            cumulene_bits[0].stereo_bond_parity[0],
            (MULT_STEREOBOND + AB_PARITY_CALC as i32) as i8
        );

        for parity in [AB_PARITY_NONE as i8, 2_i8, 5_i8] {
            let (count, skipped) = run(
                direct(parity, 50),
                ranks.clone(),
                canonical.clone(),
                numbers.clone(),
            );
            assert_eq!(count, 0);
            assert_eq!(skipped[0].stereo_bond_parity[0], parity);
            assert_eq!(skipped[1].stereo_bond_parity[0], parity);
        }

        let (count, tied) = run(
            direct(AB_PARITY_CALC as i8, 50),
            vec![1, 2, 3, 3],
            canonical.clone(),
            numbers.clone(),
        );
        assert_eq!(count, 0);
        assert_eq!(tied[0].stereo_bond_parity[0], AB_PARITY_CALC as i8);

        let mut ill_endpoint = direct(AB_PARITY_CALC as i8, 50);
        ill_endpoint[0].parity = 3;
        let (count, ill_endpoint) = run(
            ill_endpoint,
            ranks.clone(),
            canonical.clone(),
            numbers.clone(),
        );
        assert_eq!(count, 0);
        assert_eq!(ill_endpoint[0].stereo_bond_parity[0], AB_PARITY_CALC as i8);

        let mut mismatch = direct(AB_PARITY_CALC as i8, 50);
        mismatch[0].stereo_bond_parity[0] = 14;
        let (error, mismatch) = run(mismatch, ranks.clone(), canonical.clone(), numbers.clone());
        assert_eq!(error, CT_STEREOCOUNT_ERR);
        assert_eq!(mismatch[1].stereo_bond_parity[0], AB_PARITY_CALC as i8);

        let mut cumulene = vec![sp_ATOM::default(); 5];
        cumulene[0].parity = 2;
        cumulene[0].valence = 3;
        cumulene[0].neighbor[..3].copy_from_slice(&[1, 3, 4]);
        cumulene[0].stereo_bond_neighbor[0] = 3;
        cumulene[0].stereo_bond_parity[0] = 14;
        cumulene[0].stereo_bond_z_prod[0] = 50;
        cumulene[1].valence = 2;
        cumulene[1].neighbor[..2].copy_from_slice(&[0, 2]);
        cumulene[2].parity = 1;
        cumulene[2].valence = 3;
        cumulene[2].neighbor[..3].copy_from_slice(&[1, 3, 4]);
        cumulene[2].stereo_bond_neighbor[0] = 1;
        cumulene[2].stereo_bond_parity[0] = 14;
        cumulene[2].stereo_bond_z_prod[0] = 50;
        let (count, cumulene) = run(
            cumulene,
            vec![1, 5, 2, 3, 4],
            vec![1, 5, 2, 3, 4],
            vec![0, 2, 3, 4, 1],
        );
        assert_eq!(count, 1);
        assert_eq!(cumulene[0].stereo_bond_parity[0], 9);
        assert_eq!(cumulene[2].stereo_bond_parity[0], 9);
    }

    #[test]
    fn source_port__ichicans__markknownequalstereobondparities__line_1093() {
        let (result, atoms) = run_mark_known(stereo_pair(1), vec![1, 2], vec![0, 1]);
        assert_eq!(result, 1);
        assert_eq!(atoms[0].stereo_bond_parity[0], 0x41);
        assert_eq!(atoms[1].stereo_bond_parity[0], 0x41);
        assert_eq!(
            (
                atoms[0].bHasStereoOrEquToStereo,
                atoms[1].bHasStereoOrEquToStereo
            ),
            (1, 1)
        );

        let (result, atoms) = run_mark_known(stereo_pair(0x41), vec![1, 2], vec![0, 1]);
        assert_eq!(result, 0);
        assert_eq!(atoms[0].stereo_bond_parity[0], 0x41);

        let mut missing_reverse = stereo_pair(1);
        missing_reverse[1].stereo_bond_neighbor[0] = 0;
        let (result, atoms) = run_mark_known(missing_reverse, vec![1, 2], vec![0, 1]);
        assert_eq!(result, CT_STEREOCOUNT_ERR);
        assert_eq!(atoms[0].stereo_bond_parity[0], 1);

        let mut different_reverse = stereo_pair(1);
        different_reverse[1].stereo_bond_parity[0] = 2;
        assert_eq!(
            run_mark_known(different_reverse, vec![1, 2], vec![0, 1]).0,
            CT_STEREOCOUNT_ERR
        );

        let mut equivalent = stereo_pair(1);
        equivalent.extend([atom(0), atom(0)]);
        equivalent[2].valence = 1;
        equivalent[2].neighbor[0] = 3;
        equivalent[3].valence = 1;
        equivalent[3].neighbor[0] = 2;
        let (result, atoms) =
            run_mark_known(equivalent.clone(), vec![2, 4, 2, 4], vec![0, 2, 1, 3]);
        assert_eq!(result, 0);
        assert_eq!(
            (
                atoms[2].bHasStereoOrEquToStereo,
                atoms[3].bHasStereoOrEquToStereo
            ),
            (2, 2)
        );
        assert_eq!(atoms[0].stereo_bond_parity[0], 1);

        equivalent[2].valence = 2;
        let (result, atoms) = run_mark_known(equivalent, vec![2, 4, 2, 4], vec![0, 2, 1, 3]);
        assert_eq!(result, CT_STEREOCOUNT_ERR);
        assert_eq!(atoms[2].bHasStereoOrEquToStereo, 0);

        let mut cumulene = vec![atom(1), atom(0), atom(1)];
        cumulene[0].valence = 1;
        cumulene[0].neighbor[0] = 1;
        cumulene[0].stereo_bond_neighbor[0] = 3;
        cumulene[0].stereo_bond_parity[0] = 0x08 | 4;
        cumulene[1].valence = 2;
        cumulene[1].neighbor[..2].copy_from_slice(&[0, 2]);
        cumulene[2].valence = 1;
        cumulene[2].neighbor[0] = 1;
        cumulene[2].stereo_bond_neighbor[0] = 1;
        cumulene[2].stereo_bond_parity[0] = 0x08 | 4;
        let (result, atoms) = run_mark_known(cumulene.clone(), vec![1, 2, 3], vec![0, 1, 2]);
        assert_eq!(result, 1);
        assert_eq!(atoms[0].stereo_bond_parity[0], 0x48 | 4);
        assert_eq!(atoms[2].stereo_bond_parity[0], 0x48 | 4);

        cumulene[1].num_H = 1;
        let (result, atoms) = run_mark_known(cumulene, vec![1, 2, 3], vec![0, 1, 2]);
        assert_eq!(result, 0);
        assert_eq!(atoms[0].stereo_bond_parity[0], 0x08 | 4);
        assert_eq!(atoms[0].bHasStereoOrEquToStereo, 1);

        let mut heap = SourceHeap::default();
        assert_eq!(
            MarkKnownEqualStereoBondParities(
                &mut heap,
                SourceMutPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null()
            ),
            Ok(0)
        );
        assert_eq!(
            MarkKnownEqualStereoBondParities(
                &mut heap,
                SourceMutPointer::null(),
                -1,
                SourceConstPointer::null(),
                SourceConstPointer::null()
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichicans__getnextneighborandrank__line_1419() {
        let mut heap = SourceHeap::default();
        let mut center = atom(0);
        center.valence = 5;
        center.neighbor[..5].copy_from_slice(&[1, 2, 3, 4, 5]);
        let atoms = heap
            .allocate_model_storage(vec![center, atom(0), atom(0), atom(0), atom(0), atom(0)])
            .unwrap();
        let ranks = heap.allocate_model_storage(vec![0, 8, 3, 5, 5, 9]).unwrap();
        let mut next = 77;
        let mut rank = 2;
        assert_eq!(
            GetNextNeighborAndRank(
                &heap,
                atoms.as_const(),
                0,
                2,
                &mut next,
                &mut rank,
                ranks.as_const()
            ),
            Ok(1)
        );
        assert_eq!((next, rank), (3, 5));
        assert_eq!(
            GetNextNeighborAndRank(
                &heap,
                atoms.as_const(),
                0,
                2,
                &mut next,
                &mut rank,
                ranks.as_const()
            ),
            Ok(1)
        );
        assert_eq!((next, rank), (1, 8));
        assert_eq!(
            GetNextNeighborAndRank(
                &heap,
                atoms.as_const(),
                0,
                2,
                &mut next,
                &mut rank,
                ranks.as_const()
            ),
            Ok(1)
        );
        assert_eq!((next, rank), (5, 9));
        assert_eq!(
            GetNextNeighborAndRank(
                &heap,
                atoms.as_const(),
                0,
                2,
                &mut next,
                &mut rank,
                ranks.as_const()
            ),
            Ok(0)
        );
        assert_eq!((next, rank), (5, 9));

        let boundary_ranks = heap
            .allocate_model_storage(vec![0, MAX_ATOMS as AT_RANK, (MAX_ATOMS + 1) as AT_RANK])
            .unwrap();
        let mut boundary_center = atom(0);
        boundary_center.valence = 2;
        boundary_center.neighbor[..2].copy_from_slice(&[1, 2]);
        let boundary_atoms = heap
            .allocate_model_storage(vec![boundary_center, atom(0), atom(0)])
            .unwrap();
        let mut next = 41;
        let mut rank = 0;
        assert_eq!(
            GetNextNeighborAndRank(
                &heap,
                boundary_atoms.as_const(),
                0,
                9,
                &mut next,
                &mut rank,
                boundary_ranks.as_const()
            ),
            Ok(1)
        );
        assert_eq!((next, rank), (1, MAX_ATOMS as AT_RANK));
        assert_eq!(
            GetNextNeighborAndRank(
                &heap,
                boundary_atoms.as_const(),
                0,
                9,
                &mut next,
                &mut rank,
                boundary_ranks.as_const()
            ),
            Ok(0)
        );
        assert_eq!((next, rank), (1, MAX_ATOMS as AT_RANK));

        let empty_atoms = heap.allocate_model_storage(vec![atom(0)]).unwrap();
        let empty_ranks = heap.allocate_model_storage(vec![0]).unwrap();
        let mut next = 12;
        let mut rank = 11;
        assert_eq!(
            GetNextNeighborAndRank(
                &heap,
                empty_atoms.as_const(),
                0,
                0,
                &mut next,
                &mut rank,
                empty_ranks.as_const()
            ),
            Ok(0)
        );
        assert_eq!((next, rank), (12, 11));
    }

    #[test]
    fn source_port__ichicans__getandchecknextneighbors__line_1456() {
        fn paired_atoms() -> Vec<sp_ATOM> {
            let mut atoms = vec![atom(0); 4];
            atoms[0].valence = 1;
            atoms[0].neighbor[0] = 1;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 3;
            atoms
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            visited1: Vec<AT_RANK>,
            visited2: Vec<AT_RANK>,
        ) -> (i32, AT_RANK, AT_RANK) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let canonical = heap.allocate_model_storage(vec![0, 5, 0, 5]).unwrap();
            let visited1 = heap.allocate_model_storage(visited1).unwrap();
            let visited2 = heap.allocate_model_storage(visited2).unwrap();
            let mut n1 = (MAX_ATOMS + 1) as AT_RANK;
            let mut n2 = (MAX_ATOMS + 1) as AT_RANK;
            let result = GetAndCheckNextNeighbors(
                &heap,
                atoms.as_const(),
                0,
                9,
                2,
                9,
                &mut n1,
                &mut n2,
                visited1.as_const(),
                visited2.as_const(),
                ranks.as_const(),
                canonical.as_const(),
            )
            .unwrap();
            (result, n1, n2)
        }

        let equal_visits = vec![0, 4, 0, 4];
        assert_eq!(
            run(
                paired_atoms(),
                vec![0, 7, 0, 7],
                equal_visits.clone(),
                equal_visits.clone()
            ),
            (1, 1, 3)
        );
        assert_eq!(
            run(
                paired_atoms(),
                vec![0, 7, 0, 8],
                equal_visits.clone(),
                equal_visits.clone()
            ),
            (0, 1, 3)
        );
        assert_eq!(
            run(
                paired_atoms(),
                vec![0, 7, 0, 7],
                equal_visits.clone(),
                vec![0, 4, 0, 5]
            ),
            (0, 1, 3)
        );

        let mut no_second = paired_atoms();
        no_second[2].valence = 0;
        assert_eq!(
            run(
                no_second,
                vec![0, 7, 0, 7],
                equal_visits.clone(),
                equal_visits.clone()
            ),
            (0, 1, (MAX_ATOMS + 1) as AT_RANK)
        );

        let mut stereo = paired_atoms();
        stereo[0].stereo_bond_neighbor[0] = 4;
        stereo[0].stereo_bond_parity[0] = 1;
        stereo[2].stereo_bond_neighbor[0] = 2;
        stereo[2].stereo_bond_parity[0] = 1;
        assert_eq!(
            run(
                stereo.clone(),
                vec![0, 7, 0, 7],
                equal_visits.clone(),
                equal_visits.clone()
            ),
            (1, 1, 3)
        );

        let mut one_stereo = stereo.clone();
        one_stereo[2].stereo_bond_neighbor[0] = 0;
        assert_eq!(
            run(
                one_stereo,
                vec![0, 7, 0, 7],
                equal_visits.clone(),
                equal_visits.clone()
            )
            .0,
            0
        );
        let mut different = stereo.clone();
        different[2].stereo_bond_parity[0] = 2;
        assert_eq!(
            run(
                different,
                vec![0, 7, 0, 7],
                equal_visits.clone(),
                equal_visits.clone()
            )
            .0,
            0
        );
        stereo[0].stereo_bond_parity[0] = 3;
        stereo[2].stereo_bond_parity[0] = 3;
        assert_eq!(
            run(stereo, vec![0, 7, 0, 7], equal_visits.clone(), equal_visits).0,
            0
        );
    }

    #[test]
    fn source_port__ichicans__pathshaveidenticalknownparities__line_1530() {
        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            canonical: Vec<AT_RANK>,
            visited: Option<Vec<AT_RANK>>,
            length: AT_RANK,
        ) -> (AT_RANK, Vec<AT_RANK>, Vec<AT_RANK>) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let canonical = heap.allocate_model_storage(canonical).unwrap();
            let visited1 = heap
                .allocate_model_storage(visited.clone().unwrap_or_else(|| vec![0; 6]))
                .unwrap();
            let visited2 = if let Some(values) = visited {
                heap.allocate_model_storage(values).unwrap()
            } else {
                visited1
            };
            let result = PathsHaveIdenticalKnownParities(
                &mut heap,
                atoms.as_const(),
                0,
                1,
                3,
                4,
                visited1,
                visited2,
                ranks.as_const(),
                canonical.as_const(),
                length,
            )
            .unwrap();
            let first = heap.slice(visited1.as_const()).unwrap().to_vec();
            let second = heap.slice(visited2.as_const()).unwrap().to_vec();
            (result, first, second)
        }

        let mut leaves = vec![atom(0); 6];
        leaves[1].valence = 1;
        leaves[4].valence = 1;
        assert_eq!(
            run(leaves.clone(), vec![0; 6], vec![0; 6], Some(vec![0; 6]), 7),
            (8, vec![0, 8, 0, 0, 0, 0], vec![0, 0, 0, 0, 8, 0])
        );

        let mut mismatch = leaves.clone();
        mismatch[1].stereo_atom_parity = 1;
        let (result, first, second) = run(mismatch, vec![0; 6], vec![0; 6], Some(vec![0; 6]), 2);
        assert_eq!(result, 0);
        assert_eq!((first[1], second[4]), (3, 3));

        let mut unknown = leaves.clone();
        unknown[1].stereo_atom_parity = 3;
        unknown[4].stereo_atom_parity = 3;
        assert_eq!(
            run(unknown, vec![0; 6], vec![0; 6], Some(vec![0; 6]), 0).0,
            0
        );

        let mut valence_mismatch = leaves.clone();
        valence_mismatch[1].valence = 2;
        assert_eq!(
            run(
                valence_mismatch,
                vec![0; 6],
                vec![0; 6],
                Some(vec![0; 6]),
                0
            )
            .0,
            0
        );

        let mut paths = vec![atom(0); 6];
        paths[1].valence = 2;
        paths[1].neighbor[..2].copy_from_slice(&[0, 2]);
        paths[2].valence = 1;
        paths[2].neighbor[0] = 1;
        paths[4].valence = 2;
        paths[4].neighbor[..2].copy_from_slice(&[3, 5]);
        paths[5].valence = 1;
        paths[5].neighbor[0] = 4;
        let (result, first, second) = run(
            paths.clone(),
            vec![0, 4, 8, 0, 4, 8],
            vec![0, 1, 2, 0, 1, 2],
            Some(vec![0; 6]),
            1,
        );
        assert_eq!(result, 3);
        assert_eq!((first[1], first[2], second[4], second[5]), (2, 3, 2, 3));

        let mut previsited = vec![0; 6];
        previsited[2] = 9;
        previsited[5] = 9;
        let (result, first, second) = run(
            paths.clone(),
            vec![0, 4, 8, 0, 4, 8],
            vec![0, 1, 2, 0, 1, 2],
            Some(previsited),
            1,
        );
        assert_eq!(result, 2);
        assert_eq!((first[2], second[5]), (9, 9));

        let (result, shared, second) = run(
            paths,
            vec![0, 4, 8, 0, 4, 8],
            vec![0, 1, 2, 0, 1, 2],
            None,
            1,
        );
        assert_eq!(result, 3);
        assert_eq!(shared, second);
        assert_eq!((shared[1], shared[2], shared[4], shared[5]), (2, 3, 2, 3));
    }

    #[test]
    fn source_port__ichicans__removehalfstereobond__line_208() {
        fn full_atom() -> sp_ATOM {
            sp_ATOM {
                stereo_bond_neighbor: [11, 22, 33],
                stereo_bond_ord: [1, 2, 3],
                stereo_bond_z_prod: [-4, 5, -6],
                stereo_bond_parity: [7, 8, 9],
                parity: 2,
                stereo_atom_parity: 1,
                final_parity: 2,
                bHasStereoOrEquToStereo: 7,
                ..sp_ATOM::default()
            }
        }

        fn run(atom: sp_ATOM, index: i32) -> (Result<i32, SourceHeapError>, sp_ATOM) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
            let result = RemoveHalfStereoBond(&mut heap, atoms, 0, index);
            (result, heap.slice(atoms.as_const()).unwrap()[0].clone())
        }

        let (result, atom) = run(full_atom(), 1);
        assert_eq!(result, Ok(1));
        assert_eq!(atom.stereo_bond_neighbor, [11, 33, 0]);
        assert_eq!(atom.stereo_bond_ord, [1, 3, 0]);
        assert_eq!(atom.stereo_bond_z_prod, [-4, -6, 0]);
        assert_eq!(atom.stereo_bond_parity, [7, 9, 0]);
        assert_eq!(
            (
                atom.parity,
                atom.stereo_atom_parity,
                atom.final_parity,
                atom.bHasStereoOrEquToStereo
            ),
            (2, 1, 2, 7)
        );

        let (result, atom) = run(full_atom(), 2);
        assert_eq!(result, Ok(1));
        assert_eq!(atom.stereo_bond_neighbor, [11, 22, 0]);
        assert_eq!(atom.stereo_bond_ord, [1, 2, 0]);
        assert_eq!(atom.stereo_bond_z_prod, [-4, 5, 0]);
        assert_eq!(atom.stereo_bond_parity, [7, 8, 0]);

        let mut only = full_atom();
        only.stereo_bond_neighbor = [11, 0, 0];
        let (result, atom) = run(only, 0);
        assert_eq!(result, Ok(1));
        assert_eq!(atom.stereo_bond_neighbor, [0, 0, 0]);
        assert_eq!(
            (
                atom.parity,
                atom.stereo_atom_parity,
                atom.final_parity,
                atom.bHasStereoOrEquToStereo
            ),
            (0, 0, 0, 7)
        );

        let original = full_atom();
        assert_eq!(run(original.clone(), 3), (Ok(0), original.clone()));
        assert_eq!(run(original.clone(), i32::MAX), (Ok(0), original.clone()));
        assert_eq!(
            run(original.clone(), -1).0,
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut empty = original.clone();
        empty.stereo_bond_neighbor[1] = 0;
        assert_eq!(run(empty.clone(), 1), (Ok(0), empty));

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![original]).unwrap();
        assert_eq!(
            RemoveHalfStereoBond(&mut heap, atoms, 1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichicans__removeonestereobond__line_271() {
        fn full_stereo_test_atom() -> sp_ATOM {
            sp_ATOM {
                stereo_bond_neighbor: [11, 22, 33],
                stereo_bond_ord: [1, 2, 3],
                stereo_bond_z_prod: [-4, 5, -6],
                stereo_bond_parity: [7, 8, 9],
                parity: 2,
                stereo_atom_parity: 1,
                final_parity: 2,
                ..sp_ATOM::default()
            }
        }

        fn run(atoms: Vec<sp_ATOM>, current: i32, slot: i32) -> (i32, Vec<sp_ATOM>) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let result = RemoveOneStereoBond(&mut heap, atoms, current, slot).unwrap();
            (result, heap.slice(atoms.as_const()).unwrap().to_vec())
        }

        let mut pair = vec![full_stereo_test_atom(), full_stereo_test_atom()];
        pair[0].stereo_bond_neighbor = [2, 0, 0];
        pair[1].stereo_bond_neighbor = [1, 0, 0];
        let (result, atoms) = run(pair, 0, 0);
        assert_eq!(result, 1);
        for atom in &atoms {
            assert_eq!(atom.stereo_bond_neighbor, [0, 0, 0]);
            assert_eq!(
                (atom.parity, atom.stereo_atom_parity, atom.final_parity),
                (0, 0, 0)
            );
        }

        let mut shifted = vec![full_stereo_test_atom(), full_stereo_test_atom()];
        shifted[0].stereo_bond_neighbor = [2, 0, 0];
        shifted[1].stereo_bond_neighbor = [9, 1, 7];
        shifted[1].stereo_bond_ord = [3, 4, 5];
        shifted[1].stereo_bond_z_prod = [6, 7, 8];
        shifted[1].stereo_bond_parity = [1, 2, 3];
        let (result, atoms) = run(shifted, 0, 0);
        assert_eq!(result, 1);
        assert_eq!(atoms[1].stereo_bond_neighbor, [9, 7, 0]);
        assert_eq!(atoms[1].stereo_bond_ord, [3, 5, 0]);
        assert_eq!(atoms[1].stereo_bond_z_prod, [6, 8, 0]);
        assert_eq!(atoms[1].stereo_bond_parity, [1, 3, 0]);

        let mut missing = vec![full_stereo_test_atom(), full_stereo_test_atom()];
        missing[0].stereo_bond_neighbor = [2, 0, 0];
        missing[1].stereo_bond_neighbor = [9, 0, 1];
        let original = missing.clone();
        let (result, atoms) = run(missing, 0, 0);
        assert_eq!(result, 0);
        assert_eq!(atoms, original);

        let mut self_bond = full_stereo_test_atom();
        self_bond.stereo_bond_neighbor = [1, 1, 0];
        self_bond.stereo_bond_ord = [4, 5, 0];
        let (result, atoms) = run(vec![self_bond], 0, 1);
        assert_eq!(result, 0);
        assert_eq!(atoms[0].stereo_bond_neighbor, [1, 0, 0]);
        assert_eq!(atoms[0].stereo_bond_ord, [5, 0, 0]);
    }

    #[test]
    fn source_port__ichicans__removeknownnonstereobondparities__line_1608() {
        fn candidate() -> (Vec<sp_ATOM>, Vec<AT_RANK>, Vec<AT_RANK>) {
            let mut atoms = vec![atom(0); 4];
            atoms[0].valence = 3;
            atoms[0].neighbor[..3].copy_from_slice(&[3, 1, 2]);
            atoms[0].nRingSystem = 1;
            atoms[0].stereo_bond_neighbor[0] = 4;
            atoms[0].stereo_bond_ord[0] = 0;
            atoms[0].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].nRingSystem = 2;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].nRingSystem = 2;
            atoms[3].valence = 1;
            atoms[3].neighbor[0] = 0;
            atoms[3].stereo_bond_neighbor[0] = 1;
            atoms[3].stereo_bond_ord[0] = 0;
            atoms[3].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
            (atoms, vec![0, 7, 7, 0], vec![2, 5, 6, 4])
        }

        fn descriptor(first: AT_RANK, second: AT_RANK, parity: u8) -> AT_STEREO_DBLE {
            AT_STEREO_DBLE {
                at_num1: first,
                at_num2: second,
                parity,
            }
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            canonical: Vec<AT_RANK>,
            records: Vec<AT_STEREO_DBLE>,
            fail_allocation: bool,
        ) -> (i32, Vec<sp_ATOM>, Vec<AT_STEREO_DBLE>, i32, usize) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let canonical = heap.allocate_model_storage(canonical).unwrap();
            let records = heap.allocate_model_storage(records).unwrap();
            let mut stat = CANON_STAT {
                LinearCTStereoDble: records,
                nLenLinearCTStereoDble: heap.slice(records.as_const()).unwrap().len() as i32,
                ..CANON_STAT::default()
            };
            if fail_allocation {
                heap.fail_after_allocations(0);
            }
            let result = RemoveKnownNonStereoBondParities(
                &mut heap,
                atoms,
                4,
                canonical.as_const(),
                ranks.as_const(),
                &mut stat,
            )
            .unwrap();
            (
                result,
                heap.slice(atoms.as_const()).unwrap().to_vec(),
                heap.slice(records.as_const()).unwrap().to_vec(),
                stat.nLenLinearCTStereoDble,
                heap.live_allocation_count(),
            )
        }

        let (atoms, ranks, canonical) = candidate();
        let records = vec![
            descriptor(9, 8, 1),
            descriptor(4, 2, 2),
            descriptor(7, 6, 3),
        ];
        let (result, atoms, records, length, allocations) =
            run(atoms, ranks, canonical, records, false);
        assert_eq!(result, 1);
        assert_eq!(length, 2);
        assert_eq!(records[0], descriptor(9, 8, 1));
        assert_eq!(records[1], descriptor(7, 6, 3));
        assert_eq!(records[2], descriptor(7, 6, 3));
        assert_eq!(atoms[0].stereo_bond_neighbor, [0, 0, 0]);
        assert_eq!(atoms[3].stereo_bond_neighbor, [0, 0, 0]);
        assert_eq!(allocations, 4);

        let (mut atoms, ranks, canonical) = candidate();
        atoms[0].stereo_bond_parity[0] = 1;
        let original = atoms.clone();
        let result = run(atoms, ranks, canonical, vec![descriptor(4, 2, 1)], false);
        assert_eq!(result.0, 0);
        assert_eq!(result.1, original);

        let (mut atoms, mut ranks, canonical) = candidate();
        ranks[2] = 8;
        assert_eq!(
            run(
                atoms.clone(),
                ranks,
                canonical.clone(),
                vec![descriptor(4, 2, 1)],
                false
            )
            .0,
            0
        );
        atoms[1].nRingSystem = 1;
        assert_eq!(
            run(
                atoms,
                vec![0, 7, 7, 0],
                canonical,
                vec![descriptor(4, 2, 1)],
                false
            )
            .0,
            0
        );

        let (mut atoms, ranks, canonical) = candidate();
        atoms[0].stereo_bond_ord[0] = -1;
        assert_eq!(
            run(atoms, ranks, canonical, vec![descriptor(4, 2, 1)], false).0,
            CT_STEREOBOND_ERROR
        );

        let (mut atoms, ranks, canonical) = candidate();
        atoms[1].stereo_atom_parity = 1;
        let result = run(atoms, ranks, canonical, vec![descriptor(4, 2, 1)], false);
        assert_eq!(result.0, 0);
        assert_eq!(result.4, 4);

        let (mut atoms, ranks, canonical) = candidate();
        atoms[3].stereo_bond_neighbor[0] = 0;
        let original = atoms.clone();
        let result = run(atoms, ranks, canonical, vec![descriptor(4, 2, 1)], false);
        assert_eq!(result.0, CT_STEREOBOND_ERROR);
        assert_eq!(result.1, original);
        assert_eq!(result.4, 4);

        let (atoms, ranks, canonical) = candidate();
        let result = run(atoms, ranks, canonical, vec![descriptor(9, 8, 1)], false);
        assert_eq!(result.0, CT_STEREOCOUNT_ERR);
        assert_eq!(result.3, 1);
        assert_eq!(result.1[0].stereo_bond_neighbor, [0, 0, 0]);
        assert_eq!(result.4, 4);

        let (atoms, ranks, canonical) = candidate();
        let result = run(atoms, ranks, canonical, vec![], false);
        assert_eq!(result.0, 1);
        assert_eq!(result.3, 0);
        assert_eq!(result.1[0].stereo_bond_neighbor, [0, 0, 0]);

        let (atoms, ranks, canonical) = candidate();
        let original = atoms.clone();
        let result = run(atoms, ranks, canonical, vec![descriptor(4, 2, 1)], true);
        assert_eq!(result.0, CT_OUT_OF_RAM);
        assert_eq!(result.1, original);
        assert_eq!(result.4, 4);

        let mut heap = SourceHeap::default();
        assert_eq!(
            RemoveKnownNonStereoBondParities(
                &mut heap,
                SourceMutPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &mut CANON_STAT::default()
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichicans__setknownstereocenterparities__line_1714() {
        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            canonical: Vec<AT_RANK>,
            atom_numbers: Vec<AT_RANK>,
        ) -> (i32, Vec<sp_ATOM>) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let canonical = heap.allocate_model_storage(canonical).unwrap();
            let atom_numbers = heap.allocate_model_storage(atom_numbers).unwrap();
            let count = heap.slice(atoms.as_const()).unwrap().len() as i32;
            let result = SetKnownStereoCenterParities(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                count,
                canonical.as_const(),
                ranks.as_const(),
                atom_numbers.as_const(),
            )
            .unwrap();
            (result, heap.slice(atoms.as_const()).unwrap().to_vec())
        }

        let mut leaf = vec![atom(0); 2];
        leaf[0].parity = 1;
        leaf[0].stereo_atom_parity = AB_PARITY_CALC as i8;
        leaf[0].valence = 1;
        leaf[0].neighbor[0] = 1;
        let (result, atoms) = run(leaf, vec![1, 2], vec![1, 2], vec![0, 1]);
        assert_eq!(result, 1);
        assert_eq!(atoms[0].stereo_atom_parity, 1);

        let mut transposed = vec![atom(0); 4];
        transposed[0].parity = 1;
        transposed[0].stereo_atom_parity = AB_PARITY_CALC as i8;
        transposed[0].valence = 3;
        transposed[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
        let (result, atoms) = run(
            transposed,
            vec![1, 4, 3, 5],
            vec![1, 20, 10, 30],
            vec![0, 1, 2, 3],
        );
        assert_eq!(result, 1);
        assert_eq!(atoms[0].stereo_atom_parity, 2);

        let mut duplicate = vec![atom(0); 3];
        duplicate[0].parity = 2;
        duplicate[0].stereo_atom_parity = AB_PARITY_CALC as i8;
        duplicate[0].valence = 2;
        duplicate[0].neighbor[..2].copy_from_slice(&[1, 2]);
        let (result, atoms) = run(duplicate, vec![1, 3, 3], vec![1, 2, 3], vec![0, 1, 2]);
        assert_eq!(result, 0);
        assert_eq!(atoms[0].stereo_atom_parity, AB_PARITY_CALC as i8);

        let mut equivalents = vec![atom(0); 6];
        for center in [0_usize, 3] {
            equivalents[center].parity = 1;
            equivalents[center].stereo_atom_parity = AB_PARITY_CALC as i8;
            equivalents[center].valence = 2;
        }
        equivalents[0].neighbor[..2].copy_from_slice(&[1, 2]);
        equivalents[3].neighbor[..2].copy_from_slice(&[4, 5]);
        let (result, atoms) = run(
            equivalents.clone(),
            vec![2, 3, 4, 2, 3, 4],
            vec![1, 10, 20, 2, 20, 10],
            vec![0, 3, 1, 4, 2, 5],
        );
        assert_eq!(result, 0);
        assert_eq!(atoms[0].stereo_atom_parity, AB_PARITY_CALC as i8);
        assert_eq!(atoms[3].stereo_atom_parity, AB_PARITY_CALC as i8);

        let mut wrong_valence = equivalents.clone();
        wrong_valence[3].valence = 1;
        let (result, atoms) = run(
            wrong_valence,
            vec![2, 3, 4, 2, 3, 4],
            vec![1, 10, 20, 2, 20, 10],
            vec![0, 3, 1, 4, 2, 5],
        );
        assert_eq!(result, CT_STEREOCOUNT_ERR);
        assert_eq!(atoms[0].stereo_atom_parity, AB_PARITY_CALC as i8);

        let mut partial = vec![atom(0); 4];
        partial[0].parity = 1;
        partial[0].stereo_atom_parity = AB_PARITY_CALC as i8;
        partial[0].valence = 1;
        partial[0].neighbor[0] = 1;
        partial[3].valence = 2;
        let (result, atoms) = run(
            partial,
            vec![2, 3, 4, 2],
            vec![1, 2, 3, 4],
            vec![0, 3, 1, 2],
        );
        assert_eq!(result, CT_STEREOCOUNT_ERR);
        assert_eq!(atoms[0].stereo_atom_parity, 1);

        let mut missing_rank = equivalents;
        missing_rank[3].neighbor[..2].copy_from_slice(&[4, 4]);
        assert_eq!(
            run(
                missing_rank,
                vec![2, 3, 4, 2, 3, 4],
                vec![1, 10, 20, 2, 20, 10],
                vec![0, 3, 1, 4, 2, 5]
            )
            .0,
            CT_STEREOCOUNT_ERR
        );

        for mode in 0..4 {
            let mut skipped = vec![atom(0); 2];
            skipped[0].parity = 1;
            skipped[0].stereo_atom_parity = AB_PARITY_CALC as i8;
            skipped[0].valence = 1;
            skipped[0].neighbor[0] = 1;
            match mode {
                0 => skipped[0].parity = 0,
                1 => skipped[0].stereo_bond_neighbor[0] = 2,
                2 => skipped[0].stereo_atom_parity = 1,
                3 => skipped[0].parity = 3,
                _ => unreachable!(),
            }
            let (result, atoms) = run(skipped, vec![1, 2], vec![1, 2], vec![0, 1]);
            assert_eq!(result, 0);
            assert_ne!(atoms[0].stereo_atom_parity, 2);
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            SetKnownStereoCenterParities(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null()
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichicans__removeknownnonstereocenterparities__line_1833() {
        fn candidate() -> (Vec<sp_ATOM>, Vec<AT_RANK>, Vec<AT_RANK>) {
            let mut atoms = vec![atom(0); 4];
            atoms[0].parity = 1;
            atoms[0].stereo_atom_parity = AB_PARITY_CALC as i8;
            atoms[0].final_parity = 2;
            atoms[0].valence = 3;
            atoms[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
            atoms[0].nRingSystem = 1;
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].nRingSystem = 2;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].nRingSystem = 2;
            atoms[3].valence = 1;
            atoms[3].neighbor[0] = 0;
            (atoms, vec![1, 5, 5, 9], vec![7, 2, 3, 4])
        }

        fn carb(atom: AT_RANK, parity: u8) -> AT_STEREO_CARB {
            AT_STEREO_CARB {
                at_num: atom,
                parity,
            }
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            canonical: Vec<AT_RANK>,
            records: Vec<AT_STEREO_CARB>,
            fail: bool,
        ) -> (i32, Vec<sp_ATOM>, Vec<AT_STEREO_CARB>, i32, usize) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let atom_count = heap.slice(atoms.as_const()).unwrap().len() as i32;
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let canonical = heap.allocate_model_storage(canonical).unwrap();
            let records = heap.allocate_model_storage(records).unwrap();
            let mut stat = CANON_STAT {
                LinearCTStereoCarb: records,
                nLenLinearCTStereoCarb: heap.slice(records.as_const()).unwrap().len() as i32,
                ..CANON_STAT::default()
            };
            if fail {
                heap.fail_after_allocations(0);
            }
            let result = RemoveKnownNonStereoCenterParities(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                atom_count,
                canonical.as_const(),
                ranks.as_const(),
                &mut stat,
            )
            .unwrap();
            (
                result,
                heap.slice(atoms.as_const()).unwrap().to_vec(),
                heap.slice(records.as_const()).unwrap().to_vec(),
                stat.nLenLinearCTStereoCarb,
                heap.live_allocation_count(),
            )
        }

        let (atoms, ranks, canonical) = candidate();
        let result = run(
            atoms,
            ranks,
            canonical,
            vec![carb(9, 1), carb(7, 2), carb(8, 3)],
            false,
        );
        assert_eq!(result.0, 1);
        assert_eq!(
            (
                result.1[0].parity,
                result.1[0].stereo_atom_parity,
                result.1[0].final_parity
            ),
            (0, 0, 0)
        );
        assert_eq!(result.3, 2);
        assert_eq!(result.2, vec![carb(9, 1), carb(8, 3), carb(8, 3)]);
        assert_eq!(result.4, 4);

        let (atoms, ranks, canonical) = candidate();
        let result = run(atoms, ranks, canonical, vec![carb(9, 1)], false);
        assert_eq!(result.0, CT_STEREOCOUNT_ERR);
        assert_eq!(result.1[0].parity, 0);
        assert_eq!(result.3, 1);
        assert_eq!(result.4, 4);

        let (mut atoms, ranks, canonical) = candidate();
        atoms[1].stereo_atom_parity = 1;
        let result = run(atoms, ranks, canonical, vec![carb(7, 1)], false);
        assert_eq!(result.0, 0);
        assert_eq!(result.1[0].parity, 1);
        assert_eq!(result.4, 4);

        let (atoms, ranks, canonical) = candidate();
        let original = atoms.clone();
        let result = run(atoms, ranks, canonical, vec![carb(7, 1)], true);
        assert_eq!(result.0, CT_OUT_OF_RAM);
        assert_eq!(result.1, original);
        assert_eq!(result.4, 4);

        let mut one = vec![atom(0); 2];
        one[0].parity = 1;
        one[0].stereo_atom_parity = AB_PARITY_CALC as i8;
        one[0].valence = 1;
        one[0].neighbor[0] = 1;
        assert_eq!(
            run(one, vec![1, 2], vec![7, 2], vec![carb(7, 1)], false).0,
            0
        );
    }

    #[test]
    fn source_port__ichicans__markknownequalstereocenterparities__line_1942() {
        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            atom_numbers: Vec<AT_RANK>,
        ) -> (i32, Vec<sp_ATOM>) {
            let mut heap = SourceHeap::default();
            let count = atoms.len() as i32;
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let atom_numbers = heap.allocate_model_storage(atom_numbers).unwrap();
            let result = MarkKnownEqualStereoCenterParities(
                &mut heap,
                atoms,
                count,
                ranks.as_const(),
                atom_numbers.as_const(),
            )
            .unwrap();
            (result, heap.slice(atoms.as_const()).unwrap().to_vec())
        }

        fn centers(parities: [i8; 3]) -> Vec<sp_ATOM> {
            parities
                .into_iter()
                .map(|stereo_atom_parity| sp_ATOM {
                    parity: 1,
                    stereo_atom_parity,
                    ..sp_ATOM::default()
                })
                .collect()
        }

        let (count, atoms) = run(centers([1, 1, 1]), vec![3, 3, 3], vec![0, 1, 2]);
        assert_eq!(count, 3);
        for atom in atoms {
            assert_eq!(atom.stereo_atom_parity, 0x41);
            assert_eq!(atom.bHasStereoOrEquToStereo, 1);
        }

        let (count, atoms) = run(centers([1, 2, 1]), vec![3, 3, 3], vec![0, 1, 2]);
        assert_eq!(count, 0);
        assert_eq!(
            atoms
                .iter()
                .map(|atom| atom.stereo_atom_parity)
                .collect::<Vec<_>>(),
            vec![1, 2, 1]
        );
        assert!(atoms.iter().all(|atom| atom.bHasStereoOrEquToStereo == 1));

        let (count, atoms) = run(centers([1, 0, 1]), vec![3, 3, 3], vec![0, 1, 2]);
        assert_eq!(count, 0);
        assert_eq!(atoms[1].bHasStereoOrEquToStereo, 2);
        assert_eq!(atoms[0].bHasStereoOrEquToStereo, 1);

        let (count, atoms) = run(centers([6, 6, 6]), vec![3, 3, 3], vec![0, 1, 2]);
        assert_eq!(count, 0);
        assert!(atoms.iter().all(|atom| atom.stereo_atom_parity == 6));
        assert!(atoms.iter().all(|atom| atom.bHasStereoOrEquToStereo == 1));

        let mut skipped = centers([1, 1, 1]);
        skipped.push(sp_ATOM {
            parity: 1,
            stereo_atom_parity: 0,
            ..sp_ATOM::default()
        });
        skipped[0].parity = 0;
        skipped[1].stereo_bond_neighbor[0] = 1;
        skipped[2].bHasStereoOrEquToStereo = 1;
        skipped[3].stereo_atom_parity = 0x41;
        let original = skipped.clone();
        let (count, atoms) = run(skipped, vec![1, 2, 3, 4], vec![0, 1, 2, 3]);
        assert_eq!(count, 0);
        assert_eq!(atoms, original);

        let lone = vec![sp_ATOM {
            parity: 1,
            stereo_atom_parity: 1,
            ..sp_ATOM::default()
        }];
        let (count, atoms) = run(lone, vec![0], vec![0]);
        assert_eq!(count, 0);
        assert_eq!(atoms[0].bHasStereoOrEquToStereo, 0);

        let mut heap = SourceHeap::default();
        assert_eq!(
            MarkKnownEqualStereoCenterParities(
                &mut heap,
                SourceMutPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null()
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichicans__find_atoms_with_parity__line_155() {
        fn run(
            atoms: Vec<sp_ATOM>,
            visited: Vec<S_CHAR>,
            from: i32,
            current: i32,
        ) -> (Result<i32, SourceHeapError>, Vec<S_CHAR>) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let visited = heap.allocate_model_storage(visited).unwrap();
            let result =
                find_atoms_with_parity(&mut heap, atoms.as_const(), visited, from, current);
            (result, heap.slice(visited.as_const()).unwrap().to_vec())
        }

        let mut parity_atom = atom(1);
        parity_atom.valence = 0;
        assert_eq!(
            run(vec![parity_atom.clone()], vec![0], -1, 0),
            (Ok(1), vec![0])
        );
        assert_eq!(run(vec![parity_atom], vec![1], -1, 0), (Ok(0), vec![1]));

        let mut path = vec![atom(0); 4];
        path[0].valence = 1;
        path[0].neighbor[0] = 1;
        path[1].valence = 2;
        path[1].neighbor[..2].copy_from_slice(&[0, 2]);
        path[2].valence = 2;
        path[2].neighbor[..2].copy_from_slice(&[1, 3]);
        path[2].parity = 1;
        path[3].parity = 1;
        assert_eq!(
            run(path.clone(), vec![0; 4], -1, 0),
            (Ok(1), vec![1, 1, 0, 0])
        );

        path[2].parity = 0;
        assert_eq!(
            run(path.clone(), vec![0; 4], 2, 1),
            (Ok(0), vec![1, 1, 0, 0])
        );

        path[2].valence = 2;
        path[2].neighbor[..2].copy_from_slice(&[1, 0]);
        path[3].parity = 0;
        assert_eq!(run(path, vec![0; 4], -1, 0), (Ok(0), vec![1, 1, 1, 0]));

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![atom(0)]).unwrap();
        let visited = heap.allocate_model_storage(vec![0 as S_CHAR]).unwrap();
        assert_eq!(
            find_atoms_with_parity(&mut heap, atoms.as_const(), visited, -1, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichicans__sethalfstereobondilldefpariy__line_189() {
        for (slot, encoded, new_parity, expected) in [
            (0_usize, 0x29_i8, AB_PARITY_UNDF as i32, 0x2c_i8),
            (1, 0x7e, AB_PARITY_ODD as i32, 0x79),
            (2, -1_i8, AB_PARITY_EVEN as i32, -6_i8),
        ] {
            let mut heap = SourceHeap::default();
            let mut atom = sp_ATOM::default();
            atom.parity = 7;
            atom.stereo_bond_neighbor[slot] = 1;
            atom.stereo_bond_parity[slot] = encoded;
            let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
            assert_eq!(
                SetHalfStereoBondIllDefPariy(&mut heap, atoms, 0, slot as i32, new_parity,),
                Ok(1)
            );
            let atom = &heap.slice(atoms.as_const()).unwrap()[0];
            assert_eq!(atom.stereo_bond_parity[slot], expected);
            assert_eq!(atom.parity, (new_parity & BITS_PARITY) as i8);
        }

        let mut heap = SourceHeap::default();
        let mut atom = sp_ATOM::default();
        atom.parity = 5;
        atom.stereo_bond_parity[1] = 0x31;
        let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
        assert_eq!(
            SetHalfStereoBondIllDefPariy(&mut heap, atoms, 0, 1, i32::MAX),
            Ok(0)
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap()[0].parity, 5);
        assert_eq!(
            heap.slice(atoms.as_const()).unwrap()[0].stereo_bond_parity[1],
            0x31
        );
        assert_eq!(
            SetHalfStereoBondIllDefPariy(
                &mut heap,
                SourceMutPointer::null(),
                i32::MIN,
                MAX_NUM_STEREO_BOND_NEIGH as i32,
                i32::MIN,
            ),
            Ok(0)
        );
        assert_eq!(
            SetHalfStereoBondIllDefPariy(&mut heap, atoms, 0, -1, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            SetHalfStereoBondIllDefPariy(&mut heap, atoms, 1, 0, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichicans__setonestereobondilldefparity__line_242() {
        let mut heap = SourceHeap::default();
        let mut atoms = vec![sp_ATOM::default(); 2];
        atoms[0].parity = 1;
        atoms[0].stereo_bond_neighbor[1] = 2;
        atoms[0].stereo_bond_parity[1] = 0x29;
        atoms[1].parity = 2;
        atoms[1].stereo_bond_neighbor[0] = 2;
        atoms[1].stereo_bond_neighbor[1] = 1;
        atoms[1].stereo_bond_parity[1] = 0x32;
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        assert_eq!(
            SetOneStereoBondIllDefParity(&mut heap, atoms, 0, 1, AB_PARITY_UNDF as i32,),
            Ok(1)
        );
        let changed = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(changed[0].parity, AB_PARITY_UNDF as i8);
        assert_eq!(changed[0].stereo_bond_parity[1], 0x2c);
        assert_eq!(changed[1].parity, AB_PARITY_UNDF as i8);
        assert_eq!(changed[1].stereo_bond_parity[1], 0x34);
        assert_eq!(changed[1].stereo_bond_parity[0], 0);

        let mut no_reverse_heap = SourceHeap::default();
        let mut no_reverse = vec![sp_ATOM::default(); 2];
        no_reverse[0].parity = 1;
        no_reverse[0].stereo_bond_neighbor[0] = 2;
        no_reverse[0].stereo_bond_parity[0] = 0x29;
        no_reverse[1].parity = 2;
        no_reverse[1].stereo_bond_neighbor[0] = 2;
        no_reverse[1].stereo_bond_parity[0] = 0x32;
        let original = no_reverse.clone();
        let no_reverse = no_reverse_heap.allocate_model_storage(no_reverse).unwrap();
        assert_eq!(
            SetOneStereoBondIllDefParity(&mut no_reverse_heap, no_reverse, 0, 0, 4),
            Ok(0)
        );
        assert_eq!(
            no_reverse_heap.slice(no_reverse.as_const()).unwrap(),
            original
        );

        let mut terminated_heap = SourceHeap::default();
        let mut terminated = vec![sp_ATOM::default(); 2];
        terminated[0].stereo_bond_neighbor[0] = 2;
        terminated[0].stereo_bond_parity[0] = 0x11;
        terminated[1].stereo_bond_neighbor[1] = 1;
        terminated[1].stereo_bond_parity[1] = 0x22;
        let original = terminated.clone();
        let terminated = terminated_heap.allocate_model_storage(terminated).unwrap();
        assert_eq!(
            SetOneStereoBondIllDefParity(&mut terminated_heap, terminated, 0, 0, 3),
            Ok(0)
        );
        assert_eq!(
            terminated_heap.slice(terminated.as_const()).unwrap(),
            original
        );

        let mut self_heap = SourceHeap::default();
        let mut self_atom = sp_ATOM::default();
        self_atom.stereo_bond_neighbor[0] = 1;
        self_atom.stereo_bond_parity[0] = 0x39;
        let self_pointer = self_heap.allocate_model_storage(vec![self_atom]).unwrap();
        assert_eq!(
            SetOneStereoBondIllDefParity(&mut self_heap, self_pointer, 0, 0, 2),
            Ok(1)
        );
        assert_eq!(
            self_heap.slice(self_pointer.as_const()).unwrap()[0].stereo_bond_parity[0],
            0x3a
        );

        assert_eq!(
            SetOneStereoBondIllDefParity(&mut self_heap, self_pointer, 0, -1, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            SetOneStereoBondIllDefParity(
                &mut self_heap,
                self_pointer,
                0,
                MAX_NUM_STEREO_BOND_NEIGH as i32,
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let zero_neighbor = self_heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        assert_eq!(
            SetOneStereoBondIllDefParity(&mut self_heap, zero_neighbor, 0, 0, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichicans__unmarknonstereo__line_322() {
        fn candidate() -> (Vec<sp_ATOM>, Vec<AT_RANK>, Vec<AT_RANK>) {
            let mut atoms = vec![atom(0); 4];
            atoms[0].parity = 1;
            atoms[0].stereo_atom_parity = 1;
            atoms[0].final_parity = 2;
            atoms[0].bHasStereoOrEquToStereo = 1;
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[0].stereo_bond_neighbor[0] = 4;
            atoms[0].stereo_bond_ord[0] = 1;
            atoms[0].stereo_bond_z_prod[0] = -7;
            atoms[0].stereo_bond_parity[0] = 2;
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].nRingSystem = 2;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].nRingSystem = 3;
            atoms[3].stereo_bond_neighbor[0] = 1;
            atoms[3].stereo_bond_ord[0] = 2;
            atoms[3].stereo_bond_z_prod[0] = 8;
            atoms[3].stereo_bond_parity[0] = 2;
            (atoms, vec![1, 2, 2, 3], vec![0, 1, 2, 3])
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            ranks: Vec<AT_RANK>,
            numbers: Vec<AT_RANK>,
            isotopic: i32,
            fail: bool,
        ) -> (i32, Vec<sp_ATOM>, usize) {
            let mut heap = SourceHeap::default();
            let count = atoms.len() as i32;
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let numbers = heap.allocate_model_storage(numbers).unwrap();
            if fail {
                heap.fail_after_allocations(0);
            }
            let result = UnmarkNonStereo(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                count,
                ranks.as_const(),
                numbers.as_const(),
                isotopic,
            )
            .unwrap();
            (
                result,
                heap.slice(atoms.as_const()).unwrap().to_vec(),
                heap.live_allocation_count(),
            )
        }

        let (atoms, ranks, numbers) = candidate();
        let result = run(atoms, ranks, numbers, 0, false);
        assert_eq!(result.0, 1);
        assert_eq!(
            (
                result.1[0].parity,
                result.1[0].stereo_atom_parity,
                result.1[0].final_parity,
                result.1[0].bHasStereoOrEquToStereo
            ),
            (0, 0, 0, 0)
        );
        assert_eq!(result.1[0].stereo_bond_neighbor, [0, 0, 0]);
        assert_eq!(result.1[0].stereo_bond_ord, [0, 0, 0]);
        assert_eq!(result.1[0].stereo_bond_z_prod, [0, 0, 0]);
        assert_eq!(result.1[0].stereo_bond_parity, [0, 0, 0]);
        assert_eq!(result.1[3].stereo_bond_neighbor, [0, 0, 0]);
        assert_eq!(result.2, 3);

        let (mut atoms, ranks, numbers) = candidate();
        atoms[2].nRingSystem = atoms[1].nRingSystem;
        let original = atoms.clone();
        let result = run(atoms, ranks, numbers, 0, false);
        assert_eq!(result.0, 0);
        assert_eq!(result.1, original);

        let (mut atoms, ranks, numbers) = candidate();
        atoms[1].parity = 1;
        let result = run(atoms, ranks, numbers, 0, false);
        assert_eq!(result.0, 0);
        assert_eq!(result.1[0].parity, 1);

        let (atoms, ranks, numbers) = candidate();
        let original = atoms.clone();
        let result = run(atoms, ranks, numbers, 0, true);
        assert_eq!(result.0, -1);
        assert_eq!(result.1, original);
        assert_eq!(result.2, 3);

        let (mut atoms, ranks, numbers) = candidate();
        atoms[0].num_H = 3;
        let original = atoms.clone();
        let result = run(atoms, ranks, numbers, 0, false);
        assert_eq!(result.0, 0);
        assert_eq!(result.1, original);

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(Vec::<sp_ATOM>::new()).unwrap();
        let ranks = heap.allocate_model_storage(Vec::<AT_RANK>::new()).unwrap();
        let numbers = heap.allocate_model_storage(Vec::<AT_RANK>::new()).unwrap();
        assert_eq!(
            UnmarkNonStereo(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                0,
                ranks.as_const(),
                numbers.as_const(),
                1
            ),
            Ok(0)
        );
        assert_eq!(heap.live_allocation_count(), 3);
    }

    #[test]
    fn source_port__ichicans__filloutstereoparities__line_2127() {
        #[allow(clippy::too_many_arguments)]
        fn run(
            atoms: Vec<sp_ATOM>,
            canonical_ranks: Vec<AT_RANK>,
            canonical_order: Vec<AT_RANK>,
            ranks: Vec<AT_RANK>,
            atom_numbers: Vec<AT_RANK>,
            carb_max: i32,
            double_max: i32,
            isotopic: i32,
            fail_first_work_allocation: bool,
        ) -> (i32, Vec<sp_ATOM>, CANON_STAT, Vec<AT_STEREO_CARB>) {
            let mut heap = SourceHeap::default();
            let atom_count = atoms.len() as i32;
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let canonical_ranks = heap.allocate_model_storage(canonical_ranks).unwrap();
            let canonical_order = heap.allocate_model_storage(canonical_order).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let atom_numbers = heap.allocate_model_storage(atom_numbers).unwrap();
            let carb = heap
                .allocate_model_storage(vec![AT_STEREO_CARB::default(); 8])
                .unwrap();
            let double = heap
                .allocate_model_storage(vec![AT_STEREO_DBLE::default(); 8])
                .unwrap();
            let mut stat = CANON_STAT {
                LinearCTStereoCarb: carb,
                LinearCTStereoDble: double,
                nLenLinearCTStereoCarb: 71,
                nLenLinearCTStereoDble: 72,
                nMaxLenLinearCTStereoCarb: carb_max,
                nMaxLenLinearCTStereoDble: double_max,
                ..CANON_STAT::default()
            };
            if fail_first_work_allocation {
                heap.fail_after_allocations(0);
            }
            let result = FillOutStereoParities(
                &mut heap,
                atoms,
                atom_count,
                canonical_ranks,
                canonical_order,
                ranks.as_const(),
                atom_numbers.as_const(),
                &mut stat,
                &mut CANON_GLOBALS::default(),
                isotopic,
            )
            .unwrap();
            (
                result,
                heap.slice(atoms.as_const()).unwrap().to_vec(),
                stat,
                heap.slice(carb.as_const()).unwrap().to_vec(),
            )
        }

        let empty = run(vec![], vec![], vec![], vec![], vec![], 0, 0, 0, false);
        assert_eq!(empty.0, 0);
        assert_eq!(empty.2.nLenLinearCTStereoCarb, 0);
        assert_eq!(empty.2.nLenLinearCTStereoDble, 0);

        let center = sp_ATOM {
            parity: AB_PARITY_ODD as i8,
            stereo_atom_parity: AB_PARITY_ODD as i8,
            ..sp_ATOM::default()
        };
        let descriptor = run(
            vec![center.clone()],
            vec![1],
            vec![0],
            vec![1],
            vec![0],
            8,
            8,
            0,
            false,
        );
        assert_eq!(descriptor.0, 1);
        assert_eq!(descriptor.2.nLenLinearCTStereoCarb, 1);
        assert_eq!(descriptor.2.nLenLinearCTStereoDble, 0);
        assert_eq!(descriptor.3[0].at_num, 1);
        assert_eq!(descriptor.3[0].parity, AB_PARITY_ODD as u8);

        let overflow = run(
            vec![center.clone()],
            vec![1],
            vec![0],
            vec![1],
            vec![0],
            0,
            8,
            1,
            false,
        );
        assert_eq!(overflow.0, CT_OVERFLOW);
        assert_eq!(overflow.2.nLenLinearCTStereoCarb, 0);

        let allocation_failure = run(
            vec![center],
            vec![1],
            vec![0],
            vec![1],
            vec![0],
            8,
            8,
            0,
            true,
        );
        assert_eq!(allocation_failure.0, -1);
        assert_eq!(allocation_failure.2.nLenLinearCTStereoCarb, 71);
        assert_eq!(allocation_failure.2.nLenLinearCTStereoDble, 72);

        let (mut removable_atoms, removable_ranks, removable_numbers) = {
            let mut atoms = vec![atom(0); 4];
            atoms[0].parity = 1;
            atoms[0].stereo_atom_parity = 1;
            atoms[0].final_parity = 2;
            atoms[0].bHasStereoOrEquToStereo = 1;
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].nRingSystem = 2;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].nRingSystem = 3;
            (atoms, vec![1, 2, 2, 3], vec![0, 1, 2, 3])
        };
        removable_atoms[0].num_iso_H[0] = 1;
        let removed = run(
            removable_atoms,
            vec![1, 2, 3, 4],
            vec![0, 1, 2, 3],
            removable_ranks,
            removable_numbers,
            8,
            8,
            1,
            false,
        );
        assert_eq!(removed.0, 0);
        assert_eq!(removed.1[0].parity, 0);
        assert_eq!(removed.1[0].stereo_atom_parity, 0);
        assert_eq!(removed.1[0].final_parity, 0);
        assert_eq!(removed.1[0].bHasStereoOrEquToStereo, 0);

        let mut inconsistent = vec![atom(0); 6];
        for center in [0_usize, 3] {
            inconsistent[center].parity = 1;
            inconsistent[center].stereo_atom_parity = AB_PARITY_CALC as i8;
            inconsistent[center].valence = 2;
        }
        inconsistent[0].neighbor[..2].copy_from_slice(&[1, 2]);
        inconsistent[3].neighbor[0] = 4;
        inconsistent[3].valence = 1;
        let error = run(
            inconsistent,
            vec![1, 10, 20, 2, 20, 10],
            vec![0, 3, 1, 5, 2, 4],
            vec![2, 3, 4, 2, 3, 4],
            vec![0, 3, 1, 4, 2, 5],
            8,
            8,
            0,
            false,
        );
        assert_eq!(error.0, CT_STEREOCOUNT_ERR);
    }

    #[test]
    fn source_port__ichicans__invertstereo__line_2016() {
        let mut atoms = vec![atom(AB_PARITY_ODD as i8), atom(AB_PARITY_EVEN as i8)];
        atoms[0].stereo_atom_parity = (0x40 | AB_PARITY_EVEN) as i8;
        atoms[0].final_parity = 4;
        let carb = vec![
            AT_STEREO_CARB {
                at_num: 2,
                parity: AB_PARITY_ODD as u8,
            },
            AT_STEREO_CARB {
                at_num: 1,
                parity: 3,
            },
        ];
        let (mut heap, atoms, ranks, canonical, mut stat) = fixture(atoms, carb, vec![]);
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                1
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(canonical.as_const()).unwrap()[..2], [1, 0]);
        assert_eq!(
            heap.slice(atoms.as_const()).unwrap()[0].parity,
            AB_PARITY_EVEN as i8
        );
        assert_eq!(
            heap.slice(atoms.as_const()).unwrap()[0].stereo_atom_parity,
            0x41
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap()[0].final_parity, 4);
        assert_eq!(
            heap.slice(stat.LinearCTStereoCarb.as_const()).unwrap()[0].parity,
            AB_PARITY_EVEN as u8
        );

        let mut atoms = vec![atom(1), atom(2)];
        atoms[0].stereo_bond_parity[0] = 0x08 | 1;
        atoms[1].stereo_bond_parity[0] = 0x08 | 2;
        atoms[0].stereo_bond_neighbor[0] = 2;
        atoms[1].stereo_bond_neighbor[0] = 1;
        let double = vec![
            AT_STEREO_DBLE {
                at_num1: 2,
                at_num2: 1,
                parity: 2,
            },
            AT_STEREO_DBLE {
                at_num1: 1,
                at_num2: 2,
                parity: 4,
            },
        ];
        let (mut heap, atoms, ranks, canonical, mut stat) = fixture(atoms, vec![], double);
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                1
            ),
            Ok(1)
        );
        let result = heap.slice(atoms.as_const()).unwrap();
        assert_eq!((result[0].parity, result[1].parity), (2, 2));
        assert_eq!(result[0].stereo_bond_parity[0], 0x08 | 2);
        assert_eq!(result[1].stereo_bond_parity[0], 0x08 | 1);
        assert_eq!(
            heap.slice(stat.LinearCTStereoDble.as_const()).unwrap()[0].parity,
            1
        );

        let mut atoms = vec![atom(1), atom(2)];
        atoms[0].stereo_bond_parity[0] = 0x10 | 1;
        atoms[1].stereo_bond_parity[0] = 0x10 | 2;
        let (mut heap, atoms, ranks, canonical, mut stat) = fixture(
            atoms.clone(),
            vec![],
            vec![AT_STEREO_DBLE {
                at_num1: 2,
                at_num2: 1,
                parity: 1,
            }],
        );
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                1
            ),
            Ok(0)
        );
        let result = heap.slice(atoms.as_const()).unwrap();
        assert_eq!((result[0].parity, result[1].parity), (1, 2));
        assert_eq!(result[0].stereo_bond_parity[0], 0x10 | 1);
        assert_eq!(result[1].stereo_bond_parity[0], 0x10 | 2);
        assert_eq!(
            heap.slice(stat.LinearCTStereoDble.as_const()).unwrap()[0].parity,
            1
        );

        let error_case = |mut atoms: Vec<sp_ATOM>| {
            atoms[0].stereo_bond_parity[0] = 0x08 | 1;
            atoms[1].stereo_bond_parity[0] = 0x08 | 2;
            atoms[0].stereo_bond_neighbor[0] = 2;
            atoms[1].stereo_bond_neighbor[0] = 1;
            fixture(
                atoms,
                vec![],
                vec![AT_STEREO_DBLE {
                    at_num1: 2,
                    at_num2: 1,
                    parity: 1,
                }],
            )
        };

        let (mut heap, atoms, ranks, canonical, mut stat) = error_case(vec![atom(1), atom(2)]);
        heap.slice_mut(atoms).unwrap()[0].stereo_bond_neighbor[1] = 3;
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                0
            ),
            Ok(CT_STEREOCOUNT_ERR)
        );

        let (mut heap, atoms, ranks, canonical, mut stat) = error_case(vec![atom(1), atom(2)]);
        heap.slice_mut(atoms).unwrap()[1].stereo_bond_parity[0] = 0x10 | 2;
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                0
            ),
            Ok(CT_STEREOCOUNT_ERR)
        );

        let (mut heap, atoms, ranks, canonical, mut stat) = error_case(vec![atom(1), atom(2)]);
        heap.slice_mut(atoms).unwrap()[1].stereo_bond_neighbor[0] = 2;
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                0
            ),
            Ok(CT_STEREOCOUNT_ERR)
        );

        let (mut heap, atoms, ranks, canonical, mut stat) = error_case(vec![atom(4), atom(2)]);
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                0
            ),
            Ok(CT_STEREOCOUNT_ERR)
        );

        let (mut heap, atoms, ranks, canonical, mut stat) = fixture(
            vec![atom(4), atom(2)],
            vec![AT_STEREO_CARB {
                at_num: 2,
                parity: 1,
            }],
            vec![],
        );
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                2,
                ranks.as_const(),
                canonical,
                &mut stat,
                1
            ),
            Ok(CT_STEREOCOUNT_ERR)
        );
        assert_eq!(
            heap.slice(stat.LinearCTStereoCarb.as_const()).unwrap()[0].parity,
            1
        );

        let (mut heap, atoms, ranks, canonical, mut stat) = fixture(vec![], vec![], vec![]);
        assert_eq!(
            InvertStereo(
                &mut heap,
                atoms,
                -1,
                ranks.as_const(),
                canonical,
                &mut stat,
                0
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(canonical.as_const()).unwrap(), &[99; 4]);
    }

    #[test]
    fn source_port__ichicans__getstereoneighborpos__line_2197() {
        let mut heap = SourceHeap::default();
        let mut atoms = vec![sp_ATOM::default()];
        atoms[0].stereo_bond_neighbor = [4, 8, 12];
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        assert_eq!(GetStereoNeighborPos(&heap, atoms, 0, 3), Ok(0));
        assert_eq!(GetStereoNeighborPos(&heap, atoms, 0, 7), Ok(1));
        assert_eq!(GetStereoNeighborPos(&heap, atoms, 0, 11), Ok(2));
        assert_eq!(GetStereoNeighborPos(&heap, atoms, 0, 9), Ok(-1));

        heap.slice_mut(atoms).unwrap()[0].stereo_bond_neighbor = [4, 0, 12];
        assert_eq!(GetStereoNeighborPos(&heap, atoms, 0, 11), Ok(-1));
        assert_eq!(GetStereoNeighborPos(&heap, atoms, 0, -1), Ok(-1));
        assert_eq!(
            GetStereoNeighborPos(&heap, atoms, -1, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            GetStereoNeighborPos(&heap, SourceMutPointer::null(), 0, 3),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichicans__getstereobondparity__line_2219() {
        fn fixture(
            z_product: i8,
            ranks: Vec<AT_RANK>,
        ) -> (
            SourceHeap,
            SourceMutPointer<sp_ATOM>,
            SourceMutPointer<AT_RANK>,
        ) {
            let mut atoms = vec![sp_ATOM::default(); 6];
            atoms[0].valence = 3;
            atoms[0].parity = 1;
            atoms[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
            atoms[0].stereo_bond_neighbor[0] = 2;
            atoms[0].stereo_bond_ord[0] = 0;
            atoms[0].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
            atoms[0].stereo_bond_z_prod[0] = z_product;
            atoms[1].valence = 3;
            atoms[1].parity = 1;
            atoms[1].neighbor[..3].copy_from_slice(&[0, 4, 5]);
            atoms[1].stereo_bond_neighbor[0] = 1;
            atoms[1].stereo_bond_ord[0] = 0;
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            (heap, atoms, ranks)
        }

        let (heap, atoms, ranks) = fixture(MIN_DOT_PROD as i8, vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(GetStereoBondParity(&heap, atoms, 0, 1, ranks), Ok(2));

        let (heap, atoms, ranks) = fixture(-(MIN_DOT_PROD as i8), vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(GetStereoBondParity(&heap, atoms, 0, 1, ranks), Ok(1));

        let (heap, atoms, ranks) = fixture(MIN_DOT_PROD as i8, vec![1, 2, 0, 4, 5, 6]);
        assert_eq!(GetStereoBondParity(&heap, atoms, 0, 1, ranks), Ok(0));

        let (heap, atoms, ranks) = fixture(MIN_DOT_PROD as i8, vec![1, 2, 3, 3, 5, 5]);
        assert_eq!(
            GetStereoBondParity(&heap, atoms, 0, 1, ranks),
            Ok(CT_STEREOBOND_ERROR)
        );

        let (mut heap, atoms, ranks) = fixture(0, vec![1, 2, 3, 4, 5, 6]);
        heap.slice_mut(atoms).unwrap()[0].parity = 0;
        heap.slice_mut(atoms).unwrap()[1].parity = 0;
        assert_eq!(GetStereoBondParity(&heap, atoms, 0, 1, ranks), Ok(0));
        heap.slice_mut(atoms).unwrap()[0].parity = 1;
        assert_eq!(
            GetStereoBondParity(&heap, atoms, 0, 1, ranks),
            Ok(AB_PARITY_UNDF as i32)
        );

        heap.slice_mut(atoms).unwrap()[0].stereo_bond_parity[0] = 0x09;
        heap.slice_mut(atoms).unwrap()[1].stereo_bond_neighbor[0] = 0;
        assert_eq!(GetStereoBondParity(&heap, atoms, 0, 1, ranks), Ok(1));
        heap.slice_mut(atoms).unwrap()[0].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
        assert_eq!(GetStereoBondParity(&heap, atoms, 0, 1, ranks), Ok(-1));
        heap.slice_mut(atoms).unwrap()[0].stereo_bond_neighbor[0] = 0;
        assert_eq!(GetStereoBondParity(&heap, atoms, 0, 1, ranks), Ok(-1));
        assert_eq!(
            GetStereoBondParity(&heap, SourceMutPointer::null(), 0, 1, ranks),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichicans__getpermutationparity__line_2296() {
        fn run(atom: sp_ATOM, ranks: Vec<AT_RANK>, avoid: AT_RANK) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atom = heap.allocate_model_storage(vec![atom]).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            GetPermutationParity(&heap, &mut CANON_GLOBALS::default(), atom, avoid, ranks)
        }

        let mut atom = sp_ATOM::default();
        atom.valence = 4;
        atom.neighbor[..4].copy_from_slice(&[0, 1, 2, 3]);
        assert_eq!(run(atom.clone(), vec![1, 2, 3, 4], MAX_ATOMS as u16), Ok(2));
        assert_eq!(run(atom.clone(), vec![2, 1, 3, 4], MAX_ATOMS as u16), Ok(1));
        assert_eq!(run(atom.clone(), vec![4, 3, 2, 1], MAX_ATOMS as u16), Ok(2));
        assert_eq!(run(atom.clone(), vec![0, 3, 2, 1], MAX_ATOMS as u16), Ok(0));
        assert_eq!(run(atom.clone(), vec![4, 3, 2, 1], 1), Ok(1));

        atom.valence = 0;
        assert_eq!(run(atom.clone(), vec![], 0), Ok(2));
        atom.valence = -1;
        assert_eq!(run(atom.clone(), vec![], 0), Ok(2));
        atom.valence = 5;
        assert_eq!(run(atom.clone(), vec![1; 5], 0), Ok(-1));

        atom.valence = 1;
        atom.neighbor[0] = 9;
        assert_eq!(
            run(atom, vec![1], MAX_ATOMS as u16),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichicans__getstereocenterparity__line_2340() {
        fn run(
            atom: sp_ATOM,
            ranks: Vec<AT_RANK>,
        ) -> (Result<i32, SourceHeapError>, CANON_GLOBALS) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
            let ranks = heap.allocate_model_storage(ranks).unwrap();
            let mut globals = CANON_GLOBALS::default();
            let result = GetStereoCenterParity(&heap, &mut globals, atoms, 0, ranks);
            (result, globals)
        }

        assert_eq!(run(sp_ATOM::default(), vec![]).0, Ok(0));
        let mut bond_atom = sp_ATOM {
            parity: 1,
            ..sp_ATOM::default()
        };
        bond_atom.stereo_bond_neighbor[0] = 2;
        assert_eq!(run(bond_atom, vec![]).0, Ok(-1));

        let undefined = sp_ATOM {
            parity: 3,
            ..sp_ATOM::default()
        };
        assert_eq!(run(undefined, vec![]).0, Ok(3));

        let mut center = sp_ATOM {
            valence: 4,
            parity: 1,
            ..sp_ATOM::default()
        };
        center.neighbor[..4].copy_from_slice(&[0, 1, 2, 3]);
        let (ascending, globals) = run(center.clone(), vec![1, 2, 3, 4]);
        assert_eq!(ascending, Ok(1));
        assert!(!globals.m_pn_RankForSort.is_null());
        assert_eq!(run(center.clone(), vec![2, 1, 3, 4]).0, Ok(2));
        assert_eq!(run(center.clone(), vec![4, 3, 2, 1]).0, Ok(1));
        assert_eq!(run(center.clone(), vec![1, 0, 3, 4]).0, Ok(0));

        center.parity = 2;
        assert_eq!(run(center.clone(), vec![1, 2, 3, 4]).0, Ok(2));
        center.valence = -1;
        assert_eq!(
            run(center, vec![1, 2, 3, 4]).0,
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }
}
