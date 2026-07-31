use crate::source::base::ichicans::{
    GetPermutationParity, GetStereoBondParity, GetStereoCenterParity, GetStereoNeighborPos,
    RemoveOneStereoBond, RemoveOneStereoCenter, SetOneStereoBondIllDefParity,
};
use crate::source::base::ichisort::{
    CompNeighListRanks, CompNeighListRanksOrd, CompNeighLists, CompNeighListsUpToMaxRank,
    CompNeighborsRanksCountEql, CompRank, CompRanksInvOrd, CompRanksOrd, CompareNeighListLex,
    CompareNeighListLexUpToMaxRank, CreateNeighList, FreeNeighList, inchi_qsort, insertions_sort,
    insertions_sort_AT_NUMBERS, insertions_sort_NeighList_AT_NUMBERS,
    insertions_sort_NeighList_AT_NUMBERS3, insertions_sort_NeighListBySymmAndCanonRank,
};
use crate::source::base::util::inchi_free;
use crate::source_types::local_ichimap2::{
    CHECKING_STEREOBOND, CHECKING_STEREOCENTER, COMP_STEREO_SUCCESS, MAP_MODE_C2, MAP_MODE_C2v,
    MAP_MODE_S4, MAP_MODE_STD, NEIGH_MODE_CHAIN, NEIGH_MODE_RING, NOT_WELL_DEF_UNDF,
    NOT_WELL_DEF_UNKN, PARITY_IMPOSSIBLE,
};
use crate::source_types::{
    AB_MAX_ILL_DEFINED_PARITY, AB_MAX_KNOWN_PARITY, AB_MAX_WELL_DEFINED_PARITY,
    AB_MIN_ILL_DEFINED_PARITY, AB_MIN_KNOWN_PARITY, AB_MIN_WELL_DEFINED_PARITY, AB_PARITY_CALC,
    AB_PARITY_UNDF, AT_NUMB, AT_RANK, BITS_PARITY, CANON_GLOBALS, CANON_STAT, CT_ERR_MAX,
    CT_ERR_MIN, CT_MAPCOUNT_ERR, CT_OUT_OF_RAM, CT_REMOVE_STEREO_ERR, CT_STEREOBOND_ERROR,
    CT_STEREOCOUNT_ERR, EQ_NEIGH, KNOWN_PARITIES_EQL, MASK_CUMULENE_LEN, MAX_ATOMS,
    MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BOND_NEIGH, MAX_NUM_STEREO_BONDS,
    MIN_NUM_STEREO_BOND_NEIGH, MULT_STEREOBOND, NEIGH_LIST, SourceHeap, SourceHeapError,
    SourceMutPointer, ppAT_RANK, sp_ATOM,
};

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn NumberOfTies(
    heap: &mut SourceHeap,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    length: i32,
    at_no1: i32,
    at_no2: i32,
    nNewRank: &mut AT_RANK,
    bAddStack: &mut i32,
    bMapped1: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:680 NumberOfTies
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int NumberOfTies( AT_RANK **pRankStack1,
                      AT_RANK **pRankStack2,
                      int length,
                      int at_no1,
                      int at_no2,
                      AT_RANK *nNewRank,
                      int *bAddStack,
                      int *bMapped1 )
    {

        AT_RANK *nRank1 = *pRankStack1++;
        AT_RANK *nAtomNumber1 = *pRankStack1++;  /*  ranks for mapping "1", "from" */

        AT_RANK *nRank2 = *pRankStack2++;
        AT_RANK *nAtomNumber2 = *pRankStack2++;  /*  ranks for mapping "2", "to" */

        AT_RANK r, *pTempArray;

        int iMax, i, i1, i2;

        *bAddStack = 0;
        *bMapped1 = 0;
        *nNewRank = 0;
        r = nRank1[at_no1];

        if (r != nRank2[at_no2])
        {
            return CT_MAPCOUNT_ERR; /*  atoms cannot be mapped onto each other: they have different ranks */ /*   <BRKPT> */
        }
        iMax = r - 1;
        /*  find i1 and i2 = numbers of ranks in nRank1[] and nRank2[] equal to r:  */
        for (i1 = 1; i1 <= iMax && r == nRank1[nAtomNumber1[iMax - i1]]; i1++)
        {
            ;
        }
        for (i2 = 1; i2 <= iMax && r == nRank2[nAtomNumber2[iMax - i2]]; i2++)
        {
            ;
        }
        if (i2 != i1)
            return CT_MAPCOUNT_ERR; /*  program error: must be identical number of equal ranks */ /*   <BRKPT> */
        /*  found i1 equal rank(s); preceding (smaller) non-equal rank is r-i1 */
        /*  To break the tie we have to reduce the rank r to r-i1+1 */

        /************ Note *******************************
         * IF ( i=r-1 && 0 <= i && i < num_atoms AND
         *      nRank[nAtomNumber1[i]] == r )
         * THEN:
         * nRank[nAtomNumber1[i+1]] >  r; (if i+1 < num_atoms)
         * nRank[nAtomNumber1[i-1]] <= r; (if i > 0)
         *
         * IF r = nRank[i] THEN
         * nRank[nAtomNumber1[r-1]] == r
         * nRank[nAtomNumber1[r-i-1]] <= nRank[nAtomNumber1[r-i]] (for 1 <= i < r )
         */
        if (i1 > 1)
        {
            /* int bAtFromHasAlreadyBeenMapped = 0; */
            *nNewRank = r - i1 + 1;
            /*  grab an existing or allocate a new array */
            /*  we need 4 arrays: 2 for ranks + 2 for numbers */
            for (i = 0; i < 4; i++)
            {
                if (i < 2)
                {
                    pTempArray = *pRankStack1;
                    *bMapped1 += ( pTempArray && pTempArray[0] );
                }
                else
                {
                    pTempArray = *pRankStack2;
                }
                if (!pTempArray && !( pTempArray = (AT_RANK *) inchi_malloc( length ) ))
                    return CT_OUT_OF_RAM;  /*  out of RAM */ /*   <BRKPT> */
                /*  copy "to" contents */
                switch (i)
                {
                    case 2:
                        memcpy(pTempArray, nRank2, length);
                        break;
                    case 3:
                        memcpy(pTempArray, nAtomNumber2, length);
                        break;
                }
                if (i < 2)
                    *pRankStack1++ = pTempArray;
                else
                {
                    *pRankStack2++ = pTempArray;
                }
            }
            *bAddStack = 2; /*  to break the tie we added 2 more arrays to pRankStack1 and pRankStack2 */
        }

        return i1;
    }
    */
    // END INCHI C FUNCTION: NumberOfTies

    let nRank1 = source_get(heap, pRankStack1, 0)?;
    let nAtomNumber1 = source_get(heap, pRankStack1, 1)?;
    let nRank2 = source_get(heap, pRankStack2, 0)?;
    let nAtomNumber2 = source_get(heap, pRankStack2, 1)?;

    *bAddStack = 0;
    *bMapped1 = 0;
    *nNewRank = 0;
    let r = source_get(heap, nRank1, at_no1)?;
    if r != source_get(heap, nRank2, at_no2)? {
        return Ok(CT_MAPCOUNT_ERR);
    }
    let i_max = i32::from(r).wrapping_sub(1);
    let mut i1 = 1_i32;
    while i1 <= i_max {
        let atom = source_get(heap, nAtomNumber1, i_max.wrapping_sub(i1))?;
        if r != source_get(heap, nRank1, i32::from(atom))? {
            break;
        }
        i1 = i1.wrapping_add(1);
    }
    let mut i2 = 1_i32;
    while i2 <= i_max {
        let atom = source_get(heap, nAtomNumber2, i_max.wrapping_sub(i2))?;
        if r != source_get(heap, nRank2, i32::from(atom))? {
            break;
        }
        i2 = i2.wrapping_add(1);
    }
    if i2 != i1 {
        return Ok(CT_MAPCOUNT_ERR);
    }
    if i1 > 1 {
        *nNewRank = i32::from(r).wrapping_sub(i1).wrapping_add(1) as AT_RANK;
        if length < 0 || length % 2 != 0 {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        let element_count = usize::try_from(length / 2)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        for i in 0..4_i32 {
            let (stack, slot) = if i < 2 {
                (pRankStack1, i.wrapping_add(2))
            } else {
                (pRankStack2, i)
            };
            let mut temporary = source_get(heap, stack, slot)?;
            if i < 2 && !temporary.is_null() {
                *bMapped1 = bMapped1.wrapping_add(i32::from(source_get(heap, temporary, 0)? != 0));
            }
            if temporary.is_null() {
                let mut values = Vec::new();
                if values.try_reserve_exact(element_count).is_err() {
                    return Ok(CT_OUT_OF_RAM);
                }
                values.resize(element_count, 0_u16);
                temporary = match heap.allocate(values) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed)
                    | Err(SourceHeapError::AllocationElementCountOutOfRange)
                    | Err(SourceHeapError::AllocationSizeOverflow) => return Ok(CT_OUT_OF_RAM),
                    Err(error) => return Err(error),
                };
            }
            if i == 2 || i == 3 {
                let source = if i == 2 { nRank2 } else { nAtomNumber2 };
                if temporary != source {
                    heap.with_slice_mut_and_heap(temporary, |target, heap| {
                        let source = heap.slice(source.as_const())?;
                        let source = source
                            .get(..element_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let target = target
                            .get_mut(..element_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        target.copy_from_slice(source);
                        Ok(())
                    })?;
                } else if heap.slice(source.as_const())?.len() < element_count {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
            }
            source_set(heap, stack, slot, temporary)?;
        }
        *bAddStack = 2;
    }
    Ok(i1)
}

#[allow(non_snake_case)]
pub(crate) fn ClearPreviousMappings(
    heap: &mut SourceHeap,
    pRankStack1: ppAT_RANK,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1334 ClearPreviousMappings
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int ClearPreviousMappings( AT_RANK **pRankStack1 )
    {
        int i;
        for (i = 0; pRankStack1[i]; i++)
        {
            pRankStack1[i][0] = 0;
        }

        return i;
    }
    */
    // END INCHI C FUNCTION: ClearPreviousMappings

    let mut i = 0_i32;
    loop {
        let row = source_get(heap, pRankStack1, i)?;
        if row.is_null() {
            return Ok(i);
        }
        source_set(heap, row, 0, 0)?;
        i = i.wrapping_add(1);
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn map_an_atom2(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    num_max: i32,
    at_no1: i32,
    at_no2: i32,
    nTempRank: SourceMutPointer<AT_RANK>,
    nNumMappedRanks: i32,
    pnNewNumMappedRanks: &mut i32,
    pCS: &mut CANON_STAT,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    bAddStack: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1350 map_an_atom2
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int map_an_atom2( CANON_GLOBALS *pCG,
                      int num_atoms,
                      int num_max,
                      int at_no1/*from*/,
                      int at_no2/*to*/,
                      AT_RANK *nTempRank,
                      int nNumMappedRanks,
                      int *pnNewNumMappedRanks,
                      CANON_STAT *pCS,
                      NEIGH_LIST    *NeighList,
                      AT_RANK  **pRankStack1,
                      AT_RANK  **pRankStack2,
                      int *bAddStack )
    {
        AT_RANK *nRank1, *nAtomNumber1;  /*  ranks for mapping "1", "from" */
        AT_RANK *nRank2, *nAtomNumber2;  /*  ranks for mapping "2", "to" */
        AT_RANK *nNewRank1 = NULL, *nNewAtomNumber1 = NULL;  /*  ranks for mapping "1", "from" */
        AT_RANK *nNewRank2 = NULL, *nNewAtomNumber2 = NULL;  /*  ranks for mapping "2", "to" */
        int     length = num_max * sizeof( AT_RANK );
        int     nNewNumRanks2, nNewNumRanks1;
        int     i, bAtFromHasAlreadyBeenMapped, nNumTies;
        AT_RANK nNewRank;

        nNumTies = NumberOfTies( pRankStack1, pRankStack2, length, at_no1, at_no2, &nNewRank, bAddStack, &bAtFromHasAlreadyBeenMapped );

        if (RETURNED_ERROR( nNumTies ))
        {
            return nNumTies;  /*  error */
        }

        nRank1 = *pRankStack1++;
        nAtomNumber1 = *pRankStack1++;  /*  ranks for mapping "1", "from" */

        nRank2 = *pRankStack2++;
        nAtomNumber2 = *pRankStack2++;  /*  ranks for mapping "2", "to" */

        if (nNumTies > 1)
        {

            nNewRank1 = *pRankStack1++;
            nNewAtomNumber1 = *pRankStack1++;  /*  ranks for mapping "1", "from" */

            nNewRank2 = *pRankStack2++;
            nNewAtomNumber2 = *pRankStack2++;  /*  ranks for mapping "2", "to" */
            /*  break a tie for "to" */
            memcpy(nNewRank2, nRank2, length);
            memcpy(nNewAtomNumber2, nAtomNumber2, length);
            nNewRank2[at_no2] = nNewRank;
            nNewNumRanks2 = DifferentiateRanks2( pCG, num_atoms, NeighList,
                                             nNumMappedRanks, nNewRank2, nTempRank,
                                             nNewAtomNumber2, &pCS->lNumNeighListIter, 1 );
            pCS->lNumBreakTies++;

            /*  Check whether the old mapping can be reused */
            if (2 == bAtFromHasAlreadyBeenMapped && nNewRank == nNewRank1[at_no1])
            {
                for (i = 0; i < num_atoms; i++)
                {
                    if (nNewRank1[nNewAtomNumber1[i]] != nNewRank2[nNewAtomNumber2[i]])
                    {
                        bAtFromHasAlreadyBeenMapped = 0; /*  It cannot. */
                        break;
                    }
                }
            }
            else
            {
                bAtFromHasAlreadyBeenMapped = 0;
            }
            if (2 != bAtFromHasAlreadyBeenMapped)
            {
                /*  break a tie for "from" */
                for (i = 0; pRankStack1[i]; i++)
                {
                    pRankStack1[i][0] = 0;
                }
                memcpy(nNewRank1, nRank1, length);
                memcpy(nNewAtomNumber1, nAtomNumber1, length);  /* GPF: bad nAtomNumber1 */
                nNewRank1[at_no1] = nNewRank;
                nNewNumRanks1 = DifferentiateRanks2( pCG, num_atoms, NeighList,
                                                 nNumMappedRanks, nNewRank1, nTempRank,
                                                 nNewAtomNumber1, &pCS->lNumNeighListIter, 1 );
                pCS->lNumBreakTies++;
            }
            else
            {
                nNewNumRanks1 = nNewNumRanks2;
            }

            if (nNewNumRanks1 != nNewNumRanks2)
                return CT_MAPCOUNT_ERR; /*  program error */ /*   <BRKPT> */
            *pnNewNumMappedRanks = nNewNumRanks2;
            /*  debug only */
            for (i = 0; i < num_atoms; i++)
            {
                if (nNewRank1[nNewAtomNumber1[i]] != nNewRank2[nNewAtomNumber2[i]])
                {
                    return CT_MAPCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                }
            }
        }
        else
        {
            *pnNewNumMappedRanks = nNumMappedRanks;
        }

        return ( nNewRank1 ) ? nNewRank1[at_no1] : nRank1[at_no1]; /*  mapping rank value */
    }
    */
    // END INCHI C FUNCTION: map_an_atom2

    let length = num_max.wrapping_mul(std::mem::size_of::<AT_RANK>() as i32);
    let mut nNewRank = 0_u16;
    let mut bAtFromHasAlreadyBeenMapped = 0_i32;
    let nNumTies = NumberOfTies(
        heap,
        pRankStack1,
        pRankStack2,
        length,
        at_no1,
        at_no2,
        &mut nNewRank,
        bAddStack,
        &mut bAtFromHasAlreadyBeenMapped,
    )?;
    if (CT_ERR_MIN..=CT_ERR_MAX).contains(&nNumTies) {
        return Ok(nNumTies);
    }

    let nRank1 = source_get(heap, pRankStack1, 0)?;
    let nAtomNumber1 = source_get(heap, pRankStack1, 1)?;
    let nRank2 = source_get(heap, pRankStack2, 0)?;
    let nAtomNumber2 = source_get(heap, pRankStack2, 1)?;

    let mut nNewRank1 = SourceMutPointer::null();
    if nNumTies > 1 {
        nNewRank1 = source_get(heap, pRankStack1, 2)?;
        let nNewAtomNumber1 = source_get(heap, pRankStack1, 3)?;
        let nNewRank2 = source_get(heap, pRankStack2, 2)?;
        let nNewAtomNumber2 = source_get(heap, pRankStack2, 3)?;
        let count = usize::try_from(num_max)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;

        copy_rank_prefix(heap, nNewRank2, nRank2, count)?;
        copy_rank_prefix(heap, nNewAtomNumber2, nAtomNumber2, count)?;
        source_set(heap, nNewRank2, at_no2, nNewRank)?;
        let nNewNumRanks2 = DifferentiateRanks2(
            heap,
            pCG,
            num_atoms,
            NeighList,
            nNumMappedRanks,
            nNewRank2,
            nTempRank,
            nNewAtomNumber2,
            &mut pCS.lNumNeighListIter,
            1,
        )?;
        pCS.lNumBreakTies = pCS.lNumBreakTies.wrapping_add(1);

        if bAtFromHasAlreadyBeenMapped == 2 && nNewRank == source_get(heap, nNewRank1, at_no1)? {
            let mut i = 0_i32;
            while i < num_atoms {
                let from_atom = source_get(heap, nNewAtomNumber1, i)?;
                let to_atom = source_get(heap, nNewAtomNumber2, i)?;
                if source_get(heap, nNewRank1, i32::from(from_atom))?
                    != source_get(heap, nNewRank2, i32::from(to_atom))?
                {
                    bAtFromHasAlreadyBeenMapped = 0;
                    break;
                }
                i = i.wrapping_add(1);
            }
        } else {
            bAtFromHasAlreadyBeenMapped = 0;
        }

        let nNewNumRanks1 = if bAtFromHasAlreadyBeenMapped != 2 {
            ClearPreviousMappings(heap, pRankStack1.offset(4)?)?;
            copy_rank_prefix(heap, nNewRank1, nRank1, count)?;
            copy_rank_prefix(heap, nNewAtomNumber1, nAtomNumber1, count)?;
            source_set(heap, nNewRank1, at_no1, nNewRank)?;
            let differentiated = DifferentiateRanks2(
                heap,
                pCG,
                num_atoms,
                NeighList,
                nNumMappedRanks,
                nNewRank1,
                nTempRank,
                nNewAtomNumber1,
                &mut pCS.lNumNeighListIter,
                1,
            )?;
            pCS.lNumBreakTies = pCS.lNumBreakTies.wrapping_add(1);
            differentiated
        } else {
            nNewNumRanks2
        };

        if nNewNumRanks1 != nNewNumRanks2 {
            return Ok(CT_MAPCOUNT_ERR);
        }
        *pnNewNumMappedRanks = nNewNumRanks2;
        let mut i = 0_i32;
        while i < num_atoms {
            let from_atom = source_get(heap, nNewAtomNumber1, i)?;
            let to_atom = source_get(heap, nNewAtomNumber2, i)?;
            if source_get(heap, nNewRank1, i32::from(from_atom))?
                != source_get(heap, nNewRank2, i32::from(to_atom))?
            {
                return Ok(CT_MAPCOUNT_ERR);
            }
            i = i.wrapping_add(1);
        }
    } else {
        *pnNewNumMappedRanks = nNumMappedRanks;
    }

    source_get(
        heap,
        if nNewRank1.is_null() {
            nRank1
        } else {
            nNewRank1
        },
        at_no1,
    )
    .map(i32::from)
}

fn copy_rank_prefix(
    heap: &mut SourceHeap,
    destination: SourceMutPointer<AT_RANK>,
    source: SourceMutPointer<AT_RANK>,
    count: usize,
) -> Result<(), SourceHeapError> {
    if destination == source {
        heap.slice(source.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        return Ok(());
    }
    heap.with_slice_mut_and_heap(destination, |target, heap| {
        let source = heap
            .slice(source.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let target = target
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        target.copy_from_slice(source);
        Ok(())
    })
}

#[allow(non_snake_case)]
pub(crate) fn SortedEquInfoToRanks(
    heap: &mut SourceHeap,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    num_atoms: i32,
    bChanged: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:148 SortedEquInfoToRanks
    // INCHI✔️✔️: int SortedEquInfoToRanks( const AT_RANK* nSymmRank, AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms, int *bChanged )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* v. 1.05 - changed declaration of nNumDiffRanks as suggested by Burt Leland
    // INCHI✔️✔️:                  to avoid the problem arising on compilation with VS2015
    // INCHI✔️✔️:
    // INCHI✔️✔️:     AT_RANK        rNew, rOld, nNumDiffRanks;*/
    // INCHI✔️✔️:     AT_RANK        rNew, rOld;
    // INCHI✔️✔️:     int            nNumDiffRanks = 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int            i, j, nNumChanges = 0, i_init;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     i_init = num_atoms - 1;
    // INCHI✔️✔️:     /* djb-rwth: fixing oss-fuzz issue #69965 */
    // INCHI✔️✔️:     if (i_init >= 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         j = (int)nAtomNumber[i_init];
    // INCHI✔️✔️:         rOld = nSymmRank[j];
    // INCHI✔️✔️:         rNew = (AT_RANK)num_atoms;
    // INCHI✔️✔️:         nRank[j] = (AT_RANK)num_atoms;
    // INCHI✔️✔️:         nNumDiffRanks = 1;
    // INCHI✔️✔️:         for (i = i_init; i > 0; i--)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             j = (int)nAtomNumber[i - 1];
    // INCHI✔️✔️:
    // INCHI✔️✔️:             if (nSymmRank[j] != rOld)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nNumDiffRanks++;
    // INCHI✔️✔️:                 rNew = (AT_RANK)i;
    // INCHI✔️✔️:                 nNumChanges += (rOld != rNew + 1);
    // INCHI✔️✔️:                 rOld = nSymmRank[j];
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:
    // INCHI✔️✔️:             nRank[j] = rNew;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (bChanged)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *bChanged = ( 0 != nNumChanges );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nNumDiffRanks;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SortedEquInfoToRanks

    let mut nNumDiffRanks = 1_i32;
    let mut nNumChanges = 0_i32;
    let i_init = num_atoms.wrapping_sub(1);
    if i_init >= 0 {
        let mut j = i32::from(source_get(heap, nAtomNumber, i_init)?);
        let mut rOld = source_get(heap, nSymmRank, j)?;
        let mut rNew = num_atoms as AT_RANK;
        source_set(heap, nRank, j, num_atoms as AT_RANK)?;
        nNumDiffRanks = 1;
        let mut i = i_init;
        while i > 0 {
            j = i32::from(source_get(heap, nAtomNumber, i.wrapping_sub(1))?);
            let current_symmetry_rank = source_get(heap, nSymmRank, j)?;
            if current_symmetry_rank != rOld {
                nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
                rNew = i as AT_RANK;
                nNumChanges = nNumChanges.wrapping_add(i32::from(
                    i32::from(rOld) != i32::from(rNew).wrapping_add(1),
                ));
                rOld = current_symmetry_rank;
            }
            source_set(heap, nRank, j, rNew)?;
            i = i.wrapping_sub(1);
        }
    }
    if let Some(changed) = bChanged {
        *changed = i32::from(nNumChanges != 0);
    }
    Ok(nNumDiffRanks)
}

#[allow(non_snake_case)]
pub(crate) fn SortedRanksToEquInfo(
    heap: &mut SourceHeap,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:199 SortedRanksToEquInfo
    // INCHI✔️✔️: int SortedRanksToEquInfo( AT_RANK* nSymmRank, const AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* v. 1.05 - changed declaration of nNumDiffRanks as suggested by Burt Leland
    // INCHI✔️✔️:                  to avoid the problem arising on compilation with VS2015
    // INCHI✔️✔️:
    // INCHI✔️✔️:     AT_RANK        rNew, rOld, nNumDiffRanks;*/
    // INCHI✔️✔️:     AT_RANK        rNew, rOld;
    // INCHI✔️✔️:     int nNumDiffRanks = 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int            i, j;
    // INCHI✔️✔️:     for (i = 1, j = (int) nAtomNumber[0],
    // INCHI✔️✔️:           rOld = nRank[j], rNew = nSymmRank[j] = 1,
    // INCHI✔️✔️:           nNumDiffRanks = 1;
    // INCHI✔️✔️:             i < num_atoms;
    // INCHI✔️✔️:                         i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         j = (int) nAtomNumber[i];
    // INCHI✔️✔️:         if (nRank[j] != rOld)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nNumDiffRanks++;
    // INCHI✔️✔️:             rNew = (AT_RANK) ( i + 1 );
    // INCHI✔️✔️:             rOld = nRank[j];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         nSymmRank[j] = rNew;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nNumDiffRanks;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SortedRanksToEquInfo

    let mut j = i32::from(source_get(heap, nAtomNumber, 0)?);
    let mut rOld = source_get(heap, nRank, j)?;
    let mut rNew = 1_u16;
    source_set(heap, nSymmRank, j, rNew)?;
    let mut nNumDiffRanks = 1_i32;
    let mut i = 1_i32;
    while i < num_atoms {
        j = i32::from(source_get(heap, nAtomNumber, i)?);
        let current_rank = source_get(heap, nRank, j)?;
        if current_rank != rOld {
            nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
            rNew = i.wrapping_add(1) as AT_RANK;
            rOld = current_rank;
        }
        source_set(heap, nSymmRank, j, rNew)?;
        i = i.wrapping_add(1);
    }
    Ok(nNumDiffRanks)
}

pub(crate) fn switch_ptrs(p1: &mut SourceMutPointer<AT_RANK>, p2: &mut SourceMutPointer<AT_RANK>) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:230 switch_ptrs
    // INCHI✔️✔️: void switch_ptrs( AT_RANK **p1, AT_RANK **p2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_RANK *tmp = *p1;
    // INCHI✔️✔️:     *p1 = *p2;
    // INCHI✔️✔️:     *p2 = tmp;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: switch_ptrs

    let tmp = *p1;
    *p1 = *p2;
    *p2 = tmp;
}

fn source_get<T: Copy + 'static>(
    heap: &SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
) -> Result<T, SourceHeapError> {
    heap.slice(pointer.as_const())?
        .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

fn source_set<T: Copy + 'static>(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
    value: T,
) -> Result<(), SourceHeapError> {
    heap.slice_mut(pointer)?
        .get_mut(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)
        .map(|slot| *slot = value)
}

#[allow(non_snake_case)]
pub(crate) fn SetNewRanksFromNeighLists3(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nRank: SourceMutPointer<AT_RANK>,
    nNewRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:241 SetNewRanksFromNeighLists3
    // INCHI✔️❌: int SetNewRanksFromNeighLists3( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                 int num_atoms,
    // INCHI✔️❌:                                 NEIGH_LIST *NeighList,
    // INCHI✔️❌:                                 AT_RANK *nRank,
    // INCHI✔️❌:                                 AT_RANK *nNewRank,
    // INCHI✔️❌:                                 AT_RANK *nAtomNumber )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int     i, j, k, nNumDiffRanks, nNumNewRanks;
    // INCHI✔️❌:     AT_RANK r1, r2;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  -- nAtomNumber[] is already properly set --
    // INCHI✔️❌:     for ( i = 0; i < num_atoms; i++ ) {
    // INCHI✔️❌:         nAtomNumber[i] = (AT_RANK)i;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:     /*  set globals for qsort */
    // INCHI✔️❌:
    // INCHI✔️❌:     pCG->m_pNeighList_RankForSort = NeighList;
    // INCHI✔️❌:     pCG->m_pn_RankForSort = nRank;
    // INCHI✔️❌:     nNumDiffRanks = 0;
    // INCHI✔️❌:     nNumNewRanks = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( nNewRank, 0, num_atoms * sizeof( nNewRank[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  sorting */
    // INCHI✔️❌:     for (i = 0, r1 = 1; i < num_atoms; r1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = (int)nAtomNumber[i];
    // INCHI✔️❌:         r2 = nRank[j];
    // INCHI✔️❌:         if (r1 == r2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nNewRank[j] = r2;
    // INCHI✔️❌:             nNumDiffRanks++;
    // INCHI✔️❌:             i++;
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         r1 = r2;
    // INCHI✔️❌:         insertions_sort_AT_NUMBERS( pCG, nAtomNumber + i, (int) r2 - i, CompNeighLists );
    // INCHI✔️❌:         /*insertions_sort( nAtomNumber+i, r2-i, sizeof( nAtomNumber[0] ), CompNeighLists );*/
    // INCHI✔️❌:         j = r2 - 1;
    // INCHI✔️❌:         k = (int)nAtomNumber[j];
    // INCHI✔️❌:         nNewRank[k] = r2;
    // INCHI✔️❌:         nNumDiffRanks++;
    // INCHI✔️❌:         while (j > i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (CompareNeighListLex( NeighList[(int) nAtomNumber[j - 1]],
    // INCHI✔️❌:                 NeighList[(int) nAtomNumber[j]], nRank ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 r2 = j;
    // INCHI✔️❌:                 nNumDiffRanks++;
    // INCHI✔️❌:                 nNumNewRanks++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             j--;
    // INCHI✔️❌:             nNewRank[(int) nAtomNumber[j]] = r2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         i = r1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumNewRanks ? -nNumDiffRanks : nNumDiffRanks;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SetNewRanksFromNeighLists3

    pCG.m_pNeighList_RankForSort = NeighList.as_const();
    pCG.m_pn_RankForSort = nRank.as_const();
    let mut nNumDiffRanks = 0_i32;
    let mut nNumNewRanks = 0_i32;

    let mut clear_index = 0_i32;
    while clear_index < num_atoms {
        source_set(heap, nNewRank, clear_index, 0)?;
        clear_index = clear_index.wrapping_add(1);
    }

    let mut i = 0_i32;
    let mut r1: AT_RANK = 1;
    while i < num_atoms {
        let mut j = i32::from(source_get(heap, nAtomNumber, i)?);
        let mut r2 = source_get(heap, nRank, j)?;
        if r1 == r2 {
            source_set(heap, nNewRank, j, r2)?;
            nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
            i = i.wrapping_add(1);
            r1 = r1.wrapping_add(1);
            continue;
        }

        r1 = r2;
        insertions_sort_AT_NUMBERS(
            heap,
            nAtomNumber.offset(i64::from(i))?,
            i32::from(r2).wrapping_sub(i),
            &mut |heap, left: AT_NUMB, right: AT_NUMB| CompNeighLists(heap, left, right, pCG),
        )?;
        j = i32::from(r2).wrapping_sub(1);
        let k = i32::from(source_get(heap, nAtomNumber, j)?);
        source_set(heap, nNewRank, k, r2)?;
        nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
        while j > i {
            let previous_atom = source_get(heap, nAtomNumber, j.wrapping_sub(1))?;
            let current_atom = source_get(heap, nAtomNumber, j)?;
            let previous_list = source_get(heap, NeighList, i32::from(previous_atom))?;
            let current_list = source_get(heap, NeighList, i32::from(current_atom))?;
            if CompareNeighListLex(heap, previous_list, current_list, nRank)? != 0 {
                r2 = j as AT_RANK;
                nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
                nNumNewRanks = nNumNewRanks.wrapping_add(1);
            }
            j = j.wrapping_sub(1);
            let atom = i32::from(source_get(heap, nAtomNumber, j)?);
            source_set(heap, nNewRank, atom, r2)?;
        }
        i = i32::from(r1);
        r1 = r1.wrapping_add(1);
    }

    Ok(if nNumNewRanks != 0 {
        nNumDiffRanks.wrapping_neg()
    } else {
        nNumDiffRanks
    })
}

#[allow(non_snake_case)]
pub(crate) fn SetNewRanksFromNeighLists4(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nRank: SourceMutPointer<AT_RANK>,
    nNewRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    nMaxAtRank: AT_RANK,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:308 SetNewRanksFromNeighLists4
    // INCHI✔️❌: int SetNewRanksFromNeighLists4( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                 int num_atoms,
    // INCHI✔️❌:                                 NEIGH_LIST *NeighList,
    // INCHI✔️❌:                                 AT_RANK *nRank,
    // INCHI✔️❌:                                 AT_RANK *nNewRank,
    // INCHI✔️❌:                                 AT_RANK *nAtomNumber,
    // INCHI✔️❌:                                 AT_RANK nMaxAtRank )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int     i, j, nNumDiffRanks, nNumNewRanks;
    // INCHI✔️❌:     AT_RANK r1, r2;
    // INCHI✔️❌:     /*  -- nAtomNumber[] is already properly set --
    // INCHI✔️❌:     for ( i = 0; i < num_atoms; i++ ) {
    // INCHI✔️❌:         nAtomNumber[i] = (AT_RANK)i;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  set globals for CompNeighListsUpToMaxRank */
    // INCHI✔️❌:     pCG->m_pNeighList_RankForSort = NeighList;
    // INCHI✔️❌:     pCG->m_pn_RankForSort = nRank;
    // INCHI✔️❌:     nNumDiffRanks = 0;
    // INCHI✔️❌:     nNumNewRanks = 0;
    // INCHI✔️❌:     pCG->m_nMaxAtNeighRankForSort = nMaxAtRank;
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( nNewRank, 0, num_atoms * sizeof( nNewRank[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  sorting */
    // INCHI✔️❌:     for (i = 0, r1 = 1; i < num_atoms; r1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (r1 == ( r2 = nRank[j = (int) nAtomNumber[i]] ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* non-tied rank: singleton */
    // INCHI✔️❌:             nNewRank[j] = r2;
    // INCHI✔️❌:             nNumDiffRanks++;
    // INCHI✔️❌:             i++;
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* tied rank r2
    // INCHI✔️❌:            r2-i atoms have rank r2
    // INCHI✔️❌:            next atom after them is in position r2
    // INCHI✔️❌:         */
    // INCHI✔️❌:         r1 = r2;
    // INCHI✔️❌:
    // INCHI✔️❌:         insertions_sort_AT_NUMBERS( pCG, nAtomNumber + i,
    // INCHI✔️❌:             (int) r2 - i, CompNeighListsUpToMaxRank );
    // INCHI✔️❌:         /*insertions_sort( nAtomNumber+i, r2-i, sizeof( nAtomNumber[0] ),  CompNeighListsUpToMaxRank );*/
    // INCHI✔️❌:
    // INCHI✔️❌:         j = r2 - 1; /* prepare cycle backward, from j to i step -1 */
    // INCHI✔️❌:         nNewRank[(int) nAtomNumber[j]] = r2;
    // INCHI✔️❌:         nNumDiffRanks++;
    // INCHI✔️❌:         while (j > i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (CompareNeighListLexUpToMaxRank( NeighList[nAtomNumber[j - 1]],
    // INCHI✔️❌:                 NeighList[nAtomNumber[j]], nRank, nMaxAtRank ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 r2 = j;
    // INCHI✔️❌:                 nNumDiffRanks++;
    // INCHI✔️❌:                 nNumNewRanks++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             j--;
    // INCHI✔️❌:             nNewRank[(int) nAtomNumber[j]] = r2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         i = r1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumNewRanks ? -nNumDiffRanks : nNumDiffRanks;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SetNewRanksFromNeighLists4

    pCG.m_pNeighList_RankForSort = NeighList.as_const();
    pCG.m_pn_RankForSort = nRank.as_const();
    pCG.m_nMaxAtNeighRankForSort = nMaxAtRank;
    let mut nNumDiffRanks = 0_i32;
    let mut nNumNewRanks = 0_i32;

    let mut clear_index = 0_i32;
    while clear_index < num_atoms {
        source_set(heap, nNewRank, clear_index, 0)?;
        clear_index = clear_index.wrapping_add(1);
    }

    let mut i = 0_i32;
    let mut r1: AT_RANK = 1;
    while i < num_atoms {
        let mut j = i32::from(source_get(heap, nAtomNumber, i)?);
        let mut r2 = source_get(heap, nRank, j)?;
        if r1 == r2 {
            source_set(heap, nNewRank, j, r2)?;
            nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
            i = i.wrapping_add(1);
            r1 = r1.wrapping_add(1);
            continue;
        }

        r1 = r2;
        insertions_sort_AT_NUMBERS(
            heap,
            nAtomNumber.offset(i64::from(i))?,
            i32::from(r2).wrapping_sub(i),
            &mut |heap, left: AT_NUMB, right: AT_NUMB| {
                CompNeighListsUpToMaxRank(heap, left, right, pCG)
            },
        )?;
        j = i32::from(r2).wrapping_sub(1);
        let atom = i32::from(source_get(heap, nAtomNumber, j)?);
        source_set(heap, nNewRank, atom, r2)?;
        nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
        while j > i {
            let previous_atom = source_get(heap, nAtomNumber, j.wrapping_sub(1))?;
            let current_atom = source_get(heap, nAtomNumber, j)?;
            let previous_list = source_get(heap, NeighList, i32::from(previous_atom))?;
            let current_list = source_get(heap, NeighList, i32::from(current_atom))?;
            if CompareNeighListLexUpToMaxRank(heap, previous_list, current_list, nRank, nMaxAtRank)?
                != 0
            {
                r2 = j as AT_RANK;
                nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
                nNumNewRanks = nNumNewRanks.wrapping_add(1);
            }
            j = j.wrapping_sub(1);
            let atom = i32::from(source_get(heap, nAtomNumber, j)?);
            source_set(heap, nNewRank, atom, r2)?;
        }
        i = i32::from(r1);
        r1 = r1.wrapping_add(1);
    }

    Ok(if nNumNewRanks != 0 {
        nNumDiffRanks.wrapping_neg()
    } else {
        nNumDiffRanks
    })
}

#[allow(non_snake_case)]
pub(crate) fn SetNewRanksFromNeighLists(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nRank: SourceMutPointer<AT_RANK>,
    nNewRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    bUseAltSort: i32,
    comp: &mut dyn FnMut(
        &SourceHeap,
        AT_RANK,
        AT_RANK,
        &CANON_GLOBALS,
    ) -> Result<i32, SourceHeapError>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:380 SetNewRanksFromNeighLists
    // INCHI✔️❌: int SetNewRanksFromNeighLists( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                int num_atoms,
    // INCHI✔️❌:                                NEIGH_LIST *NeighList,
    // INCHI✔️❌:                                AT_RANK *nRank,
    // INCHI✔️❌:                                AT_RANK *nNewRank,
    // INCHI✔️❌:                                AT_RANK *nAtomNumber,
    // INCHI✔️❌:                                int bUseAltSort,
    // INCHI✔️❌:                                int( *comp )( const void *, const void *, void * ) )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int     i, nNumDiffRanks, j;
    // INCHI✔️❌:     AT_RANK nCurrentRank;
    // INCHI✔️❌:     /*  -- nAtomNumber[] is already properly set --
    // INCHI✔️❌:     for ( i = 0; i < num_atoms; i++ ) {
    // INCHI✔️❌:         nAtomNumber[i] = (AT_RANK)i;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  set globals for qsort */
    // INCHI✔️❌:     pCG->m_pNeighList_RankForSort = NeighList;
    // INCHI✔️❌:     pCG->m_pn_RankForSort = nRank;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  sorting */
    // INCHI✔️❌:     if (bUseAltSort & 1)
    // INCHI✔️❌:         tsort( pCG, nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), comp /*CompNeighListRanksOrd*/ );
    // INCHI✔️❌:     else
    // INCHI✔️❌:         inchi_qsort( pCG, nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), comp /*CompNeighListRanksOrd*/ );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #69315 */
    // INCHI✔️❌:     nNumDiffRanks = 1;
    // INCHI✔️❌:     if (num_atoms > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nCurrentRank = (AT_RANK)num_atoms;
    // INCHI✔️❌:         j = num_atoms - 1;
    // INCHI✔️❌:         nNewRank[(int)nAtomNumber[j]] = nCurrentRank;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = num_atoms - 1; i > 0; i--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  Note: CompNeighListRanks() in following line implicitly reads nRank pointed by pn_RankForSort */
    // INCHI✔️❌:             if (CompNeighListRanks(&nAtomNumber[i - 1], &nAtomNumber[i], pCG))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nNumDiffRanks++;
    // INCHI✔️❌:                 nCurrentRank = (AT_RANK)i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nNewRank[(int)nAtomNumber[i - 1]] = nCurrentRank;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumDiffRanks;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SetNewRanksFromNeighLists
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicomn.h:53 tsort
    // INCHI✔️❌: #define tsort insertions_sort
    // END INCHI ACTIVE HEADER MACRO: tsort

    pCG.m_pNeighList_RankForSort = NeighList.as_const();
    pCG.m_pn_RankForSort = nRank.as_const();

    if num_atoms < 0 {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    if num_atoms > 1 && bUseAltSort & 1 != 0 {
        insertions_sort_AT_NUMBERS(heap, nAtomNumber, num_atoms, &mut |heap, left, right| {
            comp(heap, left, right, pCG)
        })?;
    } else if num_atoms > 1 {
        let count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let mut atoms = heap
            .slice(nAtomNumber.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        let sort_result = {
            let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(&mut atoms);
            inchi_qsort(
                bytes,
                count,
                std::mem::size_of::<AT_RANK>(),
                &mut |left, right| {
                    let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                    let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                    comp(heap, left, right, pCG)
                },
            )
        };
        for (index, atom) in atoms.into_iter().enumerate() {
            source_set(heap, nAtomNumber, index as i32, atom)?;
        }
        sort_result?;
    }

    let mut nNumDiffRanks = 1_i32;
    if num_atoms > 0 {
        let mut nCurrentRank = num_atoms as AT_RANK;
        let mut j = num_atoms.wrapping_sub(1);
        let atom = i32::from(source_get(heap, nAtomNumber, j)?);
        source_set(heap, nNewRank, atom, nCurrentRank)?;

        let mut i = num_atoms.wrapping_sub(1);
        while i > 0 {
            let previous = source_get(heap, nAtomNumber, i.wrapping_sub(1))?;
            let current = source_get(heap, nAtomNumber, i)?;
            if CompNeighListRanks(heap, previous, current, pCG)? != 0 {
                nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
                nCurrentRank = i as AT_RANK;
            }
            j = i.wrapping_sub(1);
            let atom = i32::from(source_get(heap, nAtomNumber, j)?);
            source_set(heap, nNewRank, atom, nCurrentRank)?;
            i = i.wrapping_sub(1);
        }
    }

    Ok(nNumDiffRanks)
}

#[allow(non_snake_case)]
pub(crate) fn SortNeighListsBySymmAndCanonRank(
    heap: &mut SourceHeap,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:434 SortNeighListsBySymmAndCanonRank
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    void SortNeighListsBySymmAndCanonRank( int num_atoms,
                                           NEIGH_LIST *NeighList,
                                           const AT_RANK *nSymmRank,
                                           const AT_RANK *nCanonRank )
    {
        int i;
        for (i = 0; i < num_atoms; i++)
        {
            insertions_sort_NeighListBySymmAndCanonRank( NeighList[i], nSymmRank, nCanonRank );
        }
    }
    */
    // END INCHI C FUNCTION: SortNeighListsBySymmAndCanonRank

    let mut i = 0_i32;
    while i < num_atoms {
        let list = source_get(heap, NeighList, i)?;
        insertions_sort_NeighListBySymmAndCanonRank(heap, list, nSymmRank, nCanonRank)?;
        i = i.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn SortNeighLists2(
    heap: &mut SourceHeap,
    num_atoms: i32,
    nRank: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:448 SortNeighLists2
    // INCHI✔️✔️: int SortNeighLists2( int num_atoms,
    // INCHI✔️✔️:                      AT_RANK *nRank,
    // INCHI✔️✔️:                      NEIGH_LIST *NeighList,
    // INCHI✔️✔️:                      AT_RANK *nAtomNumber )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k, i;
    // INCHI✔️✔️:     AT_RANK nPrevRank = 0; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:      * on entry nRank[nAtomNumber[k]] <= nRank[nAtomNumber[k+1]]  ( k < num_atoms-1 )
    // INCHI✔️✔️:      *          nRank[nAtomNumber[k]] >= k+1                      ( k < num_atoms )
    // INCHI✔️✔️:      *          nRank[nAtomNumber[k]] == k+1 if this nRank value is not tied OR if
    // INCHI✔️✔️:      *                nRank[nAtomNumber[k]] < nRank[nAtomNumber[k+1]] OR if k = num_atoms-1.
    // INCHI✔️✔️:      *
    // INCHI✔️✔️:      */
    // INCHI✔️✔️:     for (k = 0; k < num_atoms; k++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         i = nAtomNumber[k];
    // INCHI✔️✔️: #ifdef FIX_STEREOCOUNT_ERR
    // INCHI✔️✔️:         if (NeighList[i][0] > 1)
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:         if (( nRank[i] != k + 1 || nRank[i] == nPrevRank ) && NeighList[i][0] > 1)
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /*  nRank[i] is tied (duplicated) */
    // INCHI✔️✔️:             insertions_sort_NeighList_AT_NUMBERS( NeighList[i], nRank );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         nPrevRank = nRank[i]; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SortNeighLists2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:366 FIX_STEREOCOUNT_ERR
    // INCHI✔️✔️: #define FIX_STEREOCOUNT_ERR           1 /* (2018-01-09) Supplied by DT                              */
    // END INCHI ACTIVE MACRO CONFIGURATION: FIX_STEREOCOUNT_ERR

    let mut nPrevRank: AT_RANK = 0;
    let mut k = 0_i32;
    while k < num_atoms {
        let i = i32::from(source_get(heap, nAtomNumber, k)?);
        let list = source_get(heap, NeighList, i)?;
        if source_get(heap, list, 0)? > 1 {
            insertions_sort_NeighList_AT_NUMBERS(heap, list, nRank)?;
        }
        nPrevRank = source_get(heap, nRank, i)?;
        k = k.wrapping_add(1);
    }
    let _ = nPrevRank;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn SortNeighLists3(
    heap: &mut SourceHeap,
    num_atoms: i32,
    nRank: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:482 SortNeighLists3
    // INCHI✔️❌: int SortNeighLists3( int num_atoms,
    // INCHI✔️❌:                      AT_RANK *nRank,
    // INCHI✔️❌:                      NEIGH_LIST *NeighList,
    // INCHI✔️❌:                      AT_RANK *nAtomNumber )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int k, i;
    // INCHI✔️❌:     AT_RANK nPrevRank = 0;
    // INCHI✔️❌:     /*
    // INCHI✔️❌:      * on entry nRank[nAtomNumber[k]] <= nRank[nAtomNumber[k+1]]  ( k < num_atoms-1 )
    // INCHI✔️❌:      *          nRank[nAtomNumber[k]] >= k+1                      ( k < num_atoms )
    // INCHI✔️❌:      *          nRank[nAtomNumber[k]] == k+1 if this nRank value is not tied OR if
    // INCHI✔️❌:      *                nRank[nAtomNumber[k]] < nRank[nAtomNumber[k+1]] OR if k = num_atoms-1.
    // INCHI✔️❌:      *
    // INCHI✔️❌:      */
    // INCHI✔️❌:     for (k = 0; k < num_atoms; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         i = nAtomNumber[k];
    // INCHI✔️❌:         if (( nRank[i] != k + 1 || nRank[i] == nPrevRank ) && NeighList[i][0] > 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  nRank[i] is tied (duplicated) */
    // INCHI✔️❌:             insertions_sort_NeighList_AT_NUMBERS3( NeighList[i], nRank );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nPrevRank = nRank[i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SortNeighLists3

    let mut nPrevRank: AT_RANK = 0;
    let mut k = 0_i32;
    while k < num_atoms {
        let i = i32::from(source_get(heap, nAtomNumber, k)?);
        let rank = source_get(heap, nRank, i)?;
        if rank != k.wrapping_add(1) as AT_RANK || rank == nPrevRank {
            let list = source_get(heap, NeighList, i)?;
            if source_get(heap, list, 0)? > 1 {
                let _ = insertions_sort_NeighList_AT_NUMBERS3(heap, list, nRank)?;
            }
        }
        nPrevRank = rank;
        k = k.wrapping_add(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn DifferentiateRanks2(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    mut nNumCurrRanks: i32,
    mut pnCurrRank: SourceMutPointer<AT_RANK>,
    mut pnPrevRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    lNumIter: &mut i64,
    bUseAltSort: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:518 DifferentiateRanks2
    // INCHI✔️❌: int  DifferentiateRanks2( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                           int num_atoms,
    // INCHI✔️❌:                           NEIGH_LIST *NeighList,
    // INCHI✔️❌:                           int nNumCurrRanks,
    // INCHI✔️❌:                           AT_RANK *pnCurrRank,
    // INCHI✔️❌:                           AT_RANK *pnPrevRank,
    // INCHI✔️❌:                           AT_RANK *nAtomNumber,
    // INCHI✔️❌:                           long *lNumIter,
    // INCHI✔️❌:                           int bUseAltSort )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*int nNumPrevRanks;*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  SortNeighLists2 needs sorted ranks */
    // INCHI✔️❌:     pCG->m_pn_RankForSort = pnCurrRank;
    // INCHI✔️❌:     if (bUseAltSort & 1)
    // INCHI✔️❌:         tsort( pCG, nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), CompRank /* CompRanksOrd*/ );
    // INCHI✔️❌:     else
    // INCHI✔️❌:         inchi_qsort( pCG, nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), CompRanksOrd );
    // INCHI✔️❌:
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *lNumIter += 1;
    // INCHI✔️❌:         /*nNumPrevRanks = nNumCurrRanks;*/
    // INCHI✔️❌:         switch_ptrs( &pnCurrRank, &pnPrevRank );
    // INCHI✔️❌:         SortNeighLists2( num_atoms, pnPrevRank, NeighList, nAtomNumber );
    // INCHI✔️❌:         /*  the following call creates pnCurrRank out of pnPrevRank */
    // INCHI✔️❌:         nNumCurrRanks = SetNewRanksFromNeighLists( pCG, num_atoms, NeighList, pnPrevRank, pnCurrRank, nAtomNumber,
    // INCHI✔️❌:                                                  1, CompNeighListRanksOrd );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     while ( /*nNumPrevRanks != nNumCurrRanks ||*/ memcmp( pnPrevRank, pnCurrRank, num_atoms * sizeof( pnCurrRank[0] ) ));
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumCurrRanks;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DifferentiateRanks2
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicomn.h:53 tsort
    // INCHI✔️❌: #define tsort insertions_sort
    // END INCHI ACTIVE HEADER MACRO: tsort

    pCG.m_pn_RankForSort = pnCurrRank.as_const();
    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut atoms = if count == 0 {
        Vec::new()
    } else {
        heap.slice(nAtomNumber.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let sort_result = {
        let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(&mut atoms);
        if bUseAltSort & 1 != 0 {
            insertions_sort(
                bytes,
                count,
                std::mem::size_of::<AT_RANK>(),
                &mut |left, right| {
                    let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                    let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                    CompRank(heap, left, right, pCG)
                },
            )
            .map(|_| ())
        } else {
            inchi_qsort(
                bytes,
                count,
                std::mem::size_of::<AT_RANK>(),
                &mut |left, right| {
                    let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                    let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                    CompRanksOrd(heap, left, right, pCG)
                },
            )
        }
    };
    for (index, atom) in atoms.into_iter().enumerate() {
        source_set(heap, nAtomNumber, index as i32, atom)?;
    }
    sort_result?;

    loop {
        *lNumIter = lNumIter.wrapping_add(1);
        switch_ptrs(&mut pnCurrRank, &mut pnPrevRank);
        SortNeighLists2(heap, num_atoms, pnPrevRank, NeighList, nAtomNumber)?;
        nNumCurrRanks = SetNewRanksFromNeighLists(
            heap,
            pCG,
            num_atoms,
            NeighList,
            pnPrevRank,
            pnCurrRank,
            nAtomNumber,
            1,
            &mut CompNeighListRanksOrd,
        )?;
        let ranks_match = if count == 0 {
            true
        } else {
            let previous = heap
                .slice(pnPrevRank.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let current = heap
                .slice(pnCurrRank.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            previous == current
        };
        if ranks_match {
            break;
        }
    }

    Ok(nNumCurrRanks)
}

#[allow(non_snake_case)]
pub(crate) fn DifferentiateRanks3(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    mut nNumCurrRanks: i32,
    mut pnCurrRank: SourceMutPointer<AT_RANK>,
    mut pnPrevRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    lNumIter: &mut i64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:561 DifferentiateRanks3
    // INCHI✔️❌: int  DifferentiateRanks3( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                           int num_atoms,
    // INCHI✔️❌:                           NEIGH_LIST *NeighList,
    // INCHI✔️❌:                           int nNumCurrRanks,
    // INCHI✔️❌:                           AT_RANK *pnCurrRank,
    // INCHI✔️❌:                           AT_RANK *pnPrevRank,
    // INCHI✔️❌:                           AT_RANK *nAtomNumber,
    // INCHI✔️❌:                           long *lNumIter )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         static long count = 0;
    // INCHI✔️❌:         count ++;
    // INCHI✔️❌:         if ( count == 103 ) {
    // INCHI✔️❌:             int stop=1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  SortNeighLists3 needs sorted ranks: ranks/atnumbers must have been already sorted */
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *lNumIter += 1;
    // INCHI✔️❌:         switch_ptrs( &pnCurrRank, &pnPrevRank );
    // INCHI✔️❌:         SortNeighLists3( num_atoms, pnPrevRank, NeighList, nAtomNumber );
    // INCHI✔️❌:         /*  the following call creates pnCurrRank out of pnPrevRank */
    // INCHI✔️❌:         nNumCurrRanks = SetNewRanksFromNeighLists3( pCG, num_atoms, NeighList, pnPrevRank,
    // INCHI✔️❌:                                                     pnCurrRank, nAtomNumber );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     while (nNumCurrRanks < 0 /* memcmp( pnPrevRank, pnCurrRank, num_atoms*sizeof(pnCurrRank[0]) )*/);
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumCurrRanks;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DifferentiateRanks3

    loop {
        *lNumIter = lNumIter.wrapping_add(1);
        switch_ptrs(&mut pnCurrRank, &mut pnPrevRank);
        SortNeighLists3(heap, num_atoms, pnPrevRank, NeighList, nAtomNumber)?;
        nNumCurrRanks = SetNewRanksFromNeighLists3(
            heap,
            pCG,
            num_atoms,
            NeighList,
            pnPrevRank,
            pnCurrRank,
            nAtomNumber,
        )?;
        if nNumCurrRanks >= 0 {
            break;
        }
    }
    Ok(nNumCurrRanks)
}

#[allow(non_snake_case)]
pub(crate) fn DifferentiateRanks4(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    mut nNumCurrRanks: i32,
    mut pnCurrRank: SourceMutPointer<AT_RANK>,
    mut pnPrevRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    nMaxAtRank: AT_RANK,
    lNumIter: &mut i64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:602 DifferentiateRanks4
    // INCHI✔️❌: int  DifferentiateRanks4( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList,
    // INCHI✔️❌:                           int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
    // INCHI✔️❌:                           AT_RANK *nAtomNumber, AT_RANK nMaxAtRank, long *lNumIter )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         static long count = 0;
    // INCHI✔️❌:         count ++;
    // INCHI✔️❌:         if ( count == 103 ) {
    // INCHI✔️❌:             int stop=1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     */
    // INCHI✔️❌:     /*  SortNeighLists4 needs sorted ranks: ranks/atnumbers must have been already sorted */
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *lNumIter += 1;
    // INCHI✔️❌:         switch_ptrs( &pnCurrRank, &pnPrevRank );
    // INCHI✔️❌:         SortNeighLists3( num_atoms, pnPrevRank, NeighList, nAtomNumber );
    // INCHI✔️❌:         /*  the following call creates pnCurrRank out of pnPrevRank */
    // INCHI✔️❌:         nNumCurrRanks = SetNewRanksFromNeighLists4( pCG, num_atoms, NeighList, pnPrevRank,
    // INCHI✔️❌:                                                     pnCurrRank, nAtomNumber, nMaxAtRank );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     while (nNumCurrRanks < 0 /* memcmp( pnPrevRank, pnCurrRank, num_atoms*sizeof(pnCurrRank[0]) )*/);
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumCurrRanks;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DifferentiateRanks4

    loop {
        *lNumIter = lNumIter.wrapping_add(1);
        switch_ptrs(&mut pnCurrRank, &mut pnPrevRank);
        SortNeighLists3(heap, num_atoms, pnPrevRank, NeighList, nAtomNumber)?;
        nNumCurrRanks = SetNewRanksFromNeighLists4(
            heap,
            pCG,
            num_atoms,
            NeighList,
            pnPrevRank,
            pnCurrRank,
            nAtomNumber,
            nMaxAtRank,
        )?;
        if nNumCurrRanks >= 0 {
            break;
        }
    }
    Ok(nNumCurrRanks)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn DifferentiateRanksBasic(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    mut nNumCurrRanks: i32,
    mut pnCurrRank: SourceMutPointer<AT_RANK>,
    mut pnPrevRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    lNumIter: &mut i64,
    bUseAltSort: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:637 DifferentiateRanksBasic
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int  DifferentiateRanksBasic( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList,
                                     int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                     AT_RANK *nAtomNumber, long *lNumIter, int bUseAltSort )
    {
        int nNumPrevRanks;

        /*  SortNeighLists2 needs sorted ranks */
        pCG->m_pn_RankForSort = pnCurrRank;
        if (bUseAltSort & 1)
        {
            tsort( pCG, nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), CompRank );
        }
        else
        {
            inchi_qsort( pCG, nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), CompRank );
        }

        do
        {
            *lNumIter += 1;
            nNumPrevRanks = nNumCurrRanks;
            switch_ptrs( &pnCurrRank, &pnPrevRank );
            SortNeighLists2( num_atoms, pnPrevRank, NeighList, nAtomNumber );
            /*  the following call creates pnCurrRank out of pnPrevRank */
            nNumCurrRanks = SetNewRanksFromNeighLists( pCG, num_atoms, NeighList, pnPrevRank, pnCurrRank, nAtomNumber, bUseAltSort, CompNeighListRanks );
        }
        while (nNumPrevRanks != nNumCurrRanks || memcmp( pnPrevRank, pnCurrRank, num_atoms * sizeof( pnCurrRank[0] ) ));
        return nNumCurrRanks;
    }
    */
    // END INCHI C FUNCTION: DifferentiateRanksBasic
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicomn.h:53 tsort
    // INCHI✔️❌: #define tsort insertions_sort
    // END INCHI ACTIVE HEADER MACRO: tsort

    pCG.m_pn_RankForSort = pnCurrRank.as_const();
    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut atoms = if count == 0 {
        Vec::new()
    } else {
        heap.slice(nAtomNumber.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let sort_result = {
        let bytes = bytemuck::cast_slice_mut::<AT_RANK, u8>(&mut atoms);
        if bUseAltSort & 1 != 0 {
            insertions_sort(
                bytes,
                count,
                std::mem::size_of::<AT_RANK>(),
                &mut |left, right| {
                    let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                    let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                    CompRank(heap, left, right, pCG)
                },
            )
            .map(|_| ())
        } else {
            inchi_qsort(
                bytes,
                count,
                std::mem::size_of::<AT_RANK>(),
                &mut |left, right| {
                    let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                    let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                    CompRank(heap, left, right, pCG)
                },
            )
        }
    };
    for (index, atom) in atoms.into_iter().enumerate() {
        source_set(heap, nAtomNumber, index as i32, atom)?;
    }
    sort_result?;

    loop {
        *lNumIter = lNumIter.wrapping_add(1);
        let nNumPrevRanks = nNumCurrRanks;
        switch_ptrs(&mut pnCurrRank, &mut pnPrevRank);
        SortNeighLists2(heap, num_atoms, pnPrevRank, NeighList, nAtomNumber)?;
        nNumCurrRanks = SetNewRanksFromNeighLists(
            heap,
            pCG,
            num_atoms,
            NeighList,
            pnPrevRank,
            pnCurrRank,
            nAtomNumber,
            bUseAltSort,
            &mut CompNeighListRanks,
        )?;
        if nNumPrevRanks == nNumCurrRanks {
            let ranks_match = if count == 0 {
                true
            } else {
                let previous = heap
                    .slice(pnPrevRank.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let current = heap
                    .slice(pnCurrRank.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                previous == current
            };
            if ranks_match {
                break;
            }
        }
    }
    Ok(nNumCurrRanks)
}

#[allow(non_snake_case)]
pub(crate) fn HalfStereoBondParity(
    heap: &SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    at_no1: i32,
    i_sb_neigh: i32,
    nRank: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:802 HalfStereoBondParity
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int HalfStereoBondParity( sp_ATOM *at,
                              int at_no1,
                              int i_sb_neigh,
                              const AT_RANK *nRank )
    {
        /*
           Suppose neighbors #0,#1,#2 have ranks a, b, c. Remove rank of the neighbor connected
           by the stereogenic bond (NCSB) from the a, b, c list and denote the two left as r[0], r[1],
           in the same order. Let iNCSB be an ordering number (0,1,or 2) of the NCSB.
           Assume the neighbor connected by the stereogenic bond has infinite positive rank.
           Position the half-bond so that the stereogenic bond neighbor is to the right from the atom (see below)

           Definition.
           ===========
                           if rank(X) != rank(Y) then Half-bond parity = (rank(X) > rank(Y)), that is,
            Y
             \             if ( rank(X) < rank(Y) ) then Half-bond parity is Even
              C==NCSB      if ( rank(X) > rank(Y) ) then Half-bond parity is Odd
             /             if ( rank(X) = rank(Y) ) then Half-bond parity cannot be defined
            X

            1                          2             1
             \                          \             \
              C==NCSB       C==NCSB      C==NCSB       C==NCSB
             /             /            /
            2             1            1

            Parity = 1    Parity = 1   Parity = 2    Parity = 2
            (Odd)         (Odd)       (Even) or 0   (Even) or 0

           Half-bond parity =  (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

           Consider the following cases to prove the formula:

           Case 1: 3 explicit neighbors
           ============================
           If  (1) atom's geometric parity = even (which means neighbors #0, #1, #2 are located clockwise),
           and (2) neighbors other than NCSB have different ranks, then,
           assuming that NCSB always has the largest (infinite) rank (this is consistent with
           the assumption that implicit hydrogens have smallest ranks), we have 3 possibilities:

                                     c         a          b
                                      \         \          \
                                       C==a      C==b       C==c
                                      /         /          /
                                     b         c          a

                    iNCSB      =      0          1          2
               Half-bond parity =     b>c        a<c        a>b     (0=even, 1=odd)
                                   r[0]>r[1]  r[0]<r[1]  r[0]>r[1]
               Half-bond parity
               for all 3 cases      =    (iNCSB + (r[0] > r[1]))%2

               The following slight modification will work for both odd and even geometric parity:

               Half-bond parity     =    (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

               even parity (0) => atom above the bond has lower rank than the atom below the bond.


           Case 2: 2 explicit neighbors
           ============================
           One implicit hydrogen atom H or hydrogen isotope (implicit rank=0). Assume r[1]=0

                                     H         a            Note. The same method
                                      \         \                 works for
                                       C==a      C==b
                                      /         /             N==a   and   a
                                     b         H             /              \
                                                            b                N==b
                    iNCSB       =      0         1
               Half-bond parity =     b>0       a<0
               (r[1]=0, r[0]>0)    r[0]>r[1]  r[0]<r[1]

               Half-bond parity =  (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

           Case 3: 1 explicit neighbor (NCSB)
           ==================================
           Two implicit hydrogens, (number of neighbors on non-streogenic bonds)==0:

           Atom C geometric parity:  Even               Odd          Note. The same method
                                                                           works for
                                     D                  H
                                      \                  \           Even   and   Odd
                                       C==a               C==a
                                      /                  /           H               N==a
                                     H                  D             \             /
                                                                       N==a        H
                    iNCSB =           0                0
               Half-bond parity =    (0<0)=0         (0<0)+1 = 1
               (r[1]=0, r[0]=0)    r[1]<r[0]         (r[1]<r[0])+atom_parity

               Half-parity
               for this case  =    (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

        */
        int i, j, k, iNeigh, parity, at1_parity, at_no2;
        AT_RANK r[MAX_NUM_STEREO_BOND_NEIGH];

        if (at[at_no1].valence > MAX_NUM_STEREO_BOND_NEIGH || ( at1_parity = at[at_no1].parity ) <= 0)
        {
            return 0;
        }
        if (!PARITY_WELL_DEF( at1_parity ))
        {
            if (PARITY_KNOWN( at1_parity ))
            {
                return at1_parity;
            }
            return -at1_parity;
        }
        if (0 > i_sb_neigh || i_sb_neigh >= MAX_NUM_STEREO_BOND_NEIGH)
        {
            return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
        }
        for (i = 0; i <= i_sb_neigh; i++)
        {
            if (!at[at_no1].stereo_bond_neighbor[i])
            {
                return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
            }
        }
        at_no2 = at[at_no1].neighbor[(int) at[at_no1].stereo_bond_ord[i_sb_neigh]];
        memset( r, 0, sizeof( r ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        for (i = j = 0, iNeigh = -1; i < at[at_no1].valence; i++)
        {
            if (( k = (int) at[at_no1].neighbor[i] ) == at_no2)
            {
                iNeigh = i;
            }
            else
            {
                r[j++] = nRank[k];
            }
        }
        if (iNeigh < 0 || iNeigh != at[at_no1].stereo_bond_ord[i_sb_neigh])
        {
            return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
        }
        if ((j > 0 && !r[0]) || (j > 1 && !r[1])) /* djb-rwth: addressing LLVM warning */
            return 0; /*  undefined ranks */

        if ((j == 2 && r[0] == r[1]) || iNeigh < 0) /* djb-rwth: addressing LLVM warning */
        {
            parity = AB_PARITY_CALC;  /*  cannot calculate bond parity without additional breaking ties. */
        }
        else
        {
            parity = 2 - ( at[at_no1].parity + iNeigh + ( r[1] < r[0] ) ) % 2;
        }

        return parity;
    }
    */
    // END INCHI C FUNCTION: HalfStereoBondParity

    let atoms = heap.slice(at.as_const())?;
    let atom = atoms
        .get(usize::try_from(at_no1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let at1_parity = i32::from(atom.parity);
    if i32::from(atom.valence) > MAX_NUM_STEREO_BOND_NEIGH as i32 || at1_parity <= 0 {
        return Ok(0);
    }

    let parity_value = at1_parity & BITS_PARITY as i32;
    if !(AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
        .contains(&parity_value)
    {
        if (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&parity_value) {
            return Ok(at1_parity);
        }
        return Ok(at1_parity.wrapping_neg());
    }

    if i_sb_neigh < 0 || i_sb_neigh >= MAX_NUM_STEREO_BOND_NEIGH as i32 {
        return Ok(CT_STEREOBOND_ERROR);
    }
    for index in 0..=i_sb_neigh as usize {
        if atom.stereo_bond_neighbor[index] == 0 {
            return Ok(CT_STEREOBOND_ERROR);
        }
    }

    let stereo_bond_order = i32::from(atom.stereo_bond_ord[i_sb_neigh as usize]);
    let at_no2 = *atom
        .neighbor
        .get(usize::try_from(stereo_bond_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ranks = heap.slice(nRank.as_const())?;
    let mut r = [0_u16; MAX_NUM_STEREO_BOND_NEIGH as usize];
    let mut i = 0_i32;
    let mut j = 0_i32;
    let mut iNeigh = -1_i32;
    while i < i32::from(atom.valence) {
        let k = i32::from(
            *atom
                .neighbor
                .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if k == i32::from(at_no2) {
            iNeigh = i;
        } else {
            let rank = *ranks
                .get(usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *r.get_mut(usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = rank;
            j = j.wrapping_add(1);
        }
        i = i.wrapping_add(1);
    }

    if iNeigh < 0 || iNeigh != stereo_bond_order {
        return Ok(CT_STEREOBOND_ERROR);
    }
    if (j > 0 && r[0] == 0) || (j > 1 && r[1] == 0) {
        return Ok(0);
    }

    if j == 2 && r[0] == r[1] {
        Ok(AB_PARITY_CALC as i32)
    } else {
        Ok(2_i32.wrapping_sub(
            at1_parity
                .wrapping_add(iNeigh)
                .wrapping_add(i32::from(r[1] < r[0]))
                % 2,
        ))
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn parity_of_mapped_half_bond(
    heap: &SourceHeap,
    from_at: i32,
    to_at: i32,
    from_neigh: i32,
    to_neigh: i32,
    at: SourceMutPointer<sp_ATOM>,
    mut pEN: Option<&mut EQ_NEIGH>,
    nCanonRankFrom: SourceMutPointer<AT_RANK>,
    nRankFrom: SourceMutPointer<AT_RANK>,
    nRankTo: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:958 parity_of_mapped_half_bond
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int parity_of_mapped_half_bond( int from_at,
                                    int to_at,
                                    int from_neigh,
                                    int to_neigh,
                                    sp_ATOM *at,
                                    EQ_NEIGH *pEN,
                                    const AT_RANK *nCanonRankFrom,
                                    const AT_RANK *nRankFrom,
                                    const AT_RANK *nRankTo )
    {
        int     i, j, k, num_neigh;
        int     to_sb_neigh_ord, from_sb_neigh_ord, parity;
        AT_RANK r_to[MAX_NUM_STEREO_BOND_NEIGH], at_no_to[MAX_NUM_STEREO_BOND_NEIGH];
        AT_RANK r_canon_from[MAX_NUM_STEREO_BOND_NEIGH], at_no_from[MAX_NUM_STEREO_BOND_NEIGH];
        AT_RANK r, r_sb_neigh;

        for (i = 0; i < MAX_NUM_STEREO_BOND_NEIGH; i++)
        {
            r_to[i] = r_canon_from[i] = 0;
        }

        if (pEN)
        {
            memset( pEN, 0, sizeof( *pEN ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }

        /*  for debug only */
        if (nRankFrom[from_at] != nRankTo[to_at] ||
             nRankFrom[from_neigh] != nRankTo[to_neigh] ||
             at[to_at].valence != at[from_at].valence)
        {
            return 0;  /*  program error: both atoms must be mapped */ /*   <BRKPT> */
        }

        parity = PARITY_VAL( at[to_at].parity );
        num_neigh = at[to_at].valence;

        if (num_neigh > MAX_NUM_STEREO_BOND_NEIGH || num_neigh < MIN_NUM_STEREO_BOND_NEIGH)
        {
            /*  2 neighbors are possible in case of stereo bond with implicit H */
            /*  or a stereocenter -CHD- with an implicit H */
            if (num_neigh == 1 && at[to_at].stereo_bond_neighbor[0])
            {
                /*  1 neighbor can happen in case of a terminal =CHD */
                if (PARITY_WELL_DEF( parity ))
                {
                    return 2 - parity % 2;
                }
                else
                {
                    if (parity)
                    {
                        return parity;
                    }
                    else
                    {
                        return AB_PARITY_UNDF; /*  undefined parity */
                    }
                }
            }
            return 0;  /*  program error */ /*   <BRKPT> */
        }

        if (ATOM_PARITY_KNOWN( parity ))
        {
            if (!ATOM_PARITY_WELL_DEF( parity ))
            {
                return parity;
            }
        }
        else
        {
            if (parity)
            {
                return 0; /* parity; */
            }
            else
            {
                return 0; /* AB_PARITY_UNDF; */ /*  possibly program error: undefined parity */
            }
        }

        /*  locate at[to_at].stereo_bond_neighbor[] ordering numbers */
        for (i = 0, to_sb_neigh_ord = -1; i < MAX_NUM_STEREO_BONDS && ( k = (int) at[to_at].stereo_bond_neighbor[i] ); i++)
        {
            if (k == to_neigh + 1)
            {
                to_sb_neigh_ord = i;
                break;
            }
        }
        if (to_sb_neigh_ord < 0)
        {
            return 0;  /*  program error: not a stereo bond */ /*   <BRKPT> */
        }
        to_sb_neigh_ord = (int) at[to_at].stereo_bond_ord[to_sb_neigh_ord];
        r_sb_neigh = nRankTo[(int) at[to_at].neighbor[to_sb_neigh_ord]];
        for (i = j = 0; i < num_neigh; i++)
        {
            if (i != to_sb_neigh_ord)
            {
                r_to[j] = nRankTo[(int) ( at_no_to[j] = at[to_at].neighbor[i] )];
                if (r_sb_neigh == r_to[j])
                {
                    return 0; /*  stereo bond atoms are not fully mapped */
                }
                j++;
            }
        }
        if (j + 1 != num_neigh)
        {
            return 0; /*  program error */ /*   <BRKPT> */
        }
        if (j == 1)
        {
            /*  only one neighbor; no mapping needed */
            return 2 - ( parity + 1 + to_sb_neigh_ord ) % 2;
        }
        if (j != 2)
        {
            return 0; /*  program error: j can be only 0, 1, or 2 */ /*   <BRKPT> */ /* djb-rwth: addressing coverity ID #499526 -- refer to the first comment in this line */
        }

        if (r_to[0] == r_to[1])
        {
            /*  double bond neighbors need to be mapped */
            j = 0;
            from_sb_neigh_ord = -1;
            for (i = 0; i < num_neigh; i++)
            {
                k = at[from_at].neighbor[i];
                r = nRankFrom[k];
                if (r == r_sb_neigh)
                {
                    from_sb_neigh_ord = i;   /*  we need this value only for error-checking */
                }
                else
                {
                    if (r == r_to[0])
                    {
                        r_canon_from[j] = nCanonRankFrom[k];
                        at_no_from[j] = (AT_RANK) k;
                        j++;
                    }
                    else
                    {
                        return 0; /*  program error: unexpected rank, not fully mapped adjacent to the stereo bond atoms */ /*   <BRKPT> */
                    }
                }
            }
            if (from_sb_neigh_ord < 0 || j != 2)
            {
                return 0; /*  program error: rank of a neighbor not found */ /*   <BRKPT> */
            }
            if (pEN)
            {
                /*  j == 2 */
                pEN->to_at[0] = at_no_to[0];
                pEN->to_at[1] = at_no_to[1];
                pEN->num_to = 2;           /*  number of stored in pEN->to_at[] central atom neighbors */
                pEN->rank = r_to[0];     /*  mapping rank of the tied neighbors */
                 /*  i := index of the smaller out of r_canon_from[1] and r_canon_from[0] */
                i = ( r_canon_from[1] < r_canon_from[0] );
                pEN->from_at = at_no_from[i];
                pEN->canon_rank = r_canon_from[i];
            }
            return -( (int) r_to[0] );
        }
        /*  double bond neighbors a mapped: r_to[0] != r_to[1] */
        from_sb_neigh_ord = -1;
        for (i = 0; i < num_neigh; i++)
        {
            k = at[from_at].neighbor[i];
            r = nRankFrom[k];
            if (r == r_sb_neigh)
            {
                from_sb_neigh_ord = i;  /*  we need this value only for error-checking */
            }
            else
            {
                if (r == r_to[0])
                {
                    r_canon_from[0] = nCanonRankFrom[k];
                    /* at_no_from[0]   = (AT_RANK)k; */
                }
                else
                {
                    if (r == r_to[1])
                    {
                        r_canon_from[1] = nCanonRankFrom[k];
                        /* at_no_from[1]   = (AT_RANK)k; */
                    }
                    else
                    {
                        return 0; /*  program error: unexpected rank, not fully mapped adjacent to the stereo bond atoms */ /*   <BRKPT> */
                    }
                }
            }
        }

        if (!r_canon_from[0] || !r_canon_from[1] || from_sb_neigh_ord < 0)
        {
            return 0; /*  program error: neighbor rank not found */ /*   <BRKPT> */
        }

        return 2 - ( parity + to_sb_neigh_ord + ( r_canon_from[1] < r_canon_from[0] ) ) % 2;
    }
    */
    // END INCHI C FUNCTION: parity_of_mapped_half_bond

    if let Some(en) = pEN.as_deref_mut() {
        *en = EQ_NEIGH::default();
    }
    let atoms = heap.slice(at.as_const())?;
    let from_atom = atoms
        .get(usize::try_from(from_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let to_atom = atoms
        .get(usize::try_from(to_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let canon_from = heap.slice(nCanonRankFrom.as_const())?;
    let rank_from = heap.slice(nRankFrom.as_const())?;
    let rank_to = heap.slice(nRankTo.as_const())?;
    let rank_at = |ranks: &[AT_RANK], index: i32| -> Result<AT_RANK, SourceHeapError> {
        ranks
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };

    if rank_at(rank_from, from_at)? != rank_at(rank_to, to_at)?
        || rank_at(rank_from, from_neigh)? != rank_at(rank_to, to_neigh)?
        || to_atom.valence != from_atom.valence
    {
        return Ok(0);
    }

    let parity = i32::from(to_atom.parity) & BITS_PARITY as i32;
    let num_neigh = i32::from(to_atom.valence);
    if num_neigh > MAX_NUM_STEREO_BOND_NEIGH as i32 || num_neigh < MIN_NUM_STEREO_BOND_NEIGH as i32
    {
        if num_neigh == 1 && to_atom.stereo_bond_neighbor[0] != 0 {
            if (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
                .contains(&parity)
            {
                return Ok(2_i32.wrapping_sub(parity % 2));
            }
            return Ok(if parity != 0 {
                parity
            } else {
                AB_PARITY_UNDF as i32
            });
        }
        return Ok(0);
    }
    if (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&parity) {
        if !(AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
            .contains(&parity)
        {
            return Ok(parity);
        }
    } else {
        return Ok(0);
    }

    let mut to_sb_neigh_ord = -1_i32;
    let mut i = 0_i32;
    while i < MAX_NUM_STEREO_BONDS as i32 {
        let k = i32::from(to_atom.stereo_bond_neighbor[i as usize]);
        if k == 0 {
            break;
        }
        if k == to_neigh.wrapping_add(1) {
            to_sb_neigh_ord = i;
            break;
        }
        i = i.wrapping_add(1);
    }
    if to_sb_neigh_ord < 0 {
        return Ok(0);
    }
    to_sb_neigh_ord = i32::from(to_atom.stereo_bond_ord[to_sb_neigh_ord as usize]);
    let stereo_atom = *to_atom
        .neighbor
        .get(usize::try_from(to_sb_neigh_ord).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let r_sb_neigh = rank_at(rank_to, i32::from(stereo_atom))?;
    let mut r_to = [0_u16; MAX_NUM_STEREO_BOND_NEIGH as usize];
    let mut at_no_to = [0_u16; MAX_NUM_STEREO_BOND_NEIGH as usize];
    let mut j = 0_i32;
    i = 0;
    while i < num_neigh {
        if i != to_sb_neigh_ord {
            let neighbour = *to_atom
                .neighbor
                .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            at_no_to[j as usize] = neighbour;
            r_to[j as usize] = rank_at(rank_to, i32::from(neighbour))?;
            if r_sb_neigh == r_to[j as usize] {
                return Ok(0);
            }
            j = j.wrapping_add(1);
        }
        i = i.wrapping_add(1);
    }
    if j.wrapping_add(1) != num_neigh {
        return Ok(0);
    }
    if j == 1 {
        return Ok(2_i32.wrapping_sub(parity.wrapping_add(1).wrapping_add(to_sb_neigh_ord) % 2));
    }
    if j != 2 {
        return Ok(0);
    }

    let mut r_canon_from = [0_u16; MAX_NUM_STEREO_BOND_NEIGH as usize];
    if r_to[0] == r_to[1] {
        let mut at_no_from = [0_u16; MAX_NUM_STEREO_BOND_NEIGH as usize];
        j = 0;
        let mut from_sb_neigh_ord = -1_i32;
        i = 0;
        while i < num_neigh {
            let k = *from_atom
                .neighbor
                .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let r = rank_at(rank_from, i32::from(k))?;
            if r == r_sb_neigh {
                from_sb_neigh_ord = i;
            } else if r == r_to[0] {
                r_canon_from[j as usize] = rank_at(canon_from, i32::from(k))?;
                at_no_from[j as usize] = k;
                j = j.wrapping_add(1);
            } else {
                return Ok(0);
            }
            i = i.wrapping_add(1);
        }
        if from_sb_neigh_ord < 0 || j != 2 {
            return Ok(0);
        }
        if let Some(en) = pEN {
            en.to_at[0] = at_no_to[0];
            en.to_at[1] = at_no_to[1];
            en.num_to = 2;
            en.rank = r_to[0];
            let smaller = usize::from(r_canon_from[1] < r_canon_from[0]);
            en.from_at = at_no_from[smaller];
            en.canon_rank = r_canon_from[smaller];
        }
        return Ok(i32::from(r_to[0]).wrapping_neg());
    }

    let mut from_sb_neigh_ord = -1_i32;
    i = 0;
    while i < num_neigh {
        let k = *from_atom
            .neighbor
            .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let r = rank_at(rank_from, i32::from(k))?;
        if r == r_sb_neigh {
            from_sb_neigh_ord = i;
        } else if r == r_to[0] {
            r_canon_from[0] = rank_at(canon_from, i32::from(k))?;
        } else if r == r_to[1] {
            r_canon_from[1] = rank_at(canon_from, i32::from(k))?;
        } else {
            return Ok(0);
        }
        i = i.wrapping_add(1);
    }
    if r_canon_from[0] == 0 || r_canon_from[1] == 0 || from_sb_neigh_ord < 0 {
        return Ok(0);
    }
    Ok(2_i32.wrapping_sub(
        parity
            .wrapping_add(to_sb_neigh_ord)
            .wrapping_add(i32::from(r_canon_from[1] < r_canon_from[0]))
            % 2,
    ))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn parity_of_mapped_atom2(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    from_at: i32,
    to_at: i32,
    at: SourceMutPointer<sp_ATOM>,
    mut pEN: Option<&mut EQ_NEIGH>,
    nCanonRankFrom: SourceMutPointer<AT_RANK>,
    nRankFrom: SourceMutPointer<AT_RANK>,
    nRankTo: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1168 parity_of_mapped_atom2
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int parity_of_mapped_atom2( CANON_GLOBALS *pCG,
                                int from_at,
                                int to_at,
                                const sp_ATOM *at,
                                EQ_NEIGH *pEN,
                                const AT_RANK *nCanonRankFrom,
                                const AT_RANK *nRankFrom,
                                const AT_RANK *nRankTo )
    {
        AT_RANK nNeighRankFrom[4], nNeighNumberFrom[4], nNeighRankTo[4], nNeighNumberTo[4];
        AT_RANK nNeighRankFromCanon[4], nNeighRankToCanon[4];
        int     i, j, k, num_neigh;
        int     r1, r2, r, r_canon_from_min, neigh_canon_from_min = 0, r_canon_from; /* djb-rwth: proper initialisation required */
        int     num_trans_to, num_trans_from, neigh1, neigh2; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

        num_neigh = at[to_at].valence;

        if (pEN)
        {
            memset( pEN, 0, sizeof( *pEN ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }

        /*  for debug only */
        if (nRankFrom[from_at] != nRankTo[to_at])
        {
            return 0;  /*  program error */ /*   <BRKPT> */
        }
        if (num_neigh > MAX_NUM_STEREO_ATOM_NEIGH || num_neigh < 2)
        {
            /*  2 neighbors are possible in case of stereo bond with implicit H */
            /*  or a stereocenter >CHD with two implicit H */
            if (num_neigh == 1)
            {
                /*  1 neighbor can happen in case of a terminal -CHDT or =CHD */
                if (at[to_at].parity)
                {
                    return at[to_at].parity;
                }
                else
                {
                    return AB_PARITY_UNDF; /*  undefined parity */
                }
            }
            return 0;  /*  program error */ /*   <BRKPT> */
        }
        for (i = 0; i < num_neigh; i++)
        { /*  initialization of locals */
            nNeighNumberTo[i] =
                nNeighNumberFrom[i] = i;
            nNeighRankTo[i] = nRankTo[(int) at[to_at].neighbor[i]];       /* mapping rank */
            nNeighRankFrom[i] = nRankFrom[j = (int) at[from_at].neighbor[i]]; /* mapping rank */
            nNeighRankFromCanon[i] = nCanonRankFrom[j];                     /* canonical number */
        }

        pCG->m_pn_RankForSort = nNeighRankFrom;
        pCG->m_nNumCompNeighborsRanksCountEql = 0; /*  sort mapping ranks-from */
        num_trans_from = insertions_sort( pCG, nNeighNumberFrom, num_neigh, sizeof( nNeighNumberFrom[0] ), CompNeighborsRanksCountEql ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

        if (pCG->m_nNumCompNeighborsRanksCountEql)
        {
            /*  At least 2 neighbors have equal mapping ranks (are tied). */
            /*  Find tied from-neighbors with minimal canonical rank (nCanonRankFrom[]) */
            r_canon_from_min = MAX_ATOMS + 1; /*  max possible rank + 1 */
            for (i = 1, r = 0, r1 = nNeighRankFrom[neigh1 = nNeighNumberFrom[0]]; i < num_neigh; i++, r1 = r2, neigh1 = neigh2)
            {
                r2 = nNeighRankFrom[neigh2 = nNeighNumberFrom[i]];
                if (r2 == r1)
                {
                    /*  found neighbors with tied ranks */
                    if (r != r2)
                    {
                         /*  the 1st pair of neighbor with this rank */
                        r = r2;
                        if (( r_canon_from = nNeighRankFromCanon[neigh1] ) < r_canon_from_min)
                        {
                            r_canon_from_min = r_canon_from; /*  min canon rank */
                            neigh_canon_from_min = neigh1;       /*  neighbor number */
                        }
                    }
                    if (( r_canon_from = nNeighRankFromCanon[neigh2] ) < r_canon_from_min)
                    {
                        r_canon_from_min = r_canon_from;
                        neigh_canon_from_min = neigh2;
                    }
                }
            }
            if (r)
            {
                /*  neighbors with tied ranks have been found => parity cannot be determined without additional mapping */
                /*  find to-neighbors on which neigh_canon_from_min can be mapped */
                r1 = nNeighRankFrom[neigh_canon_from_min];
                if (pEN)
                {
                    for (i = j = 0; i < num_neigh; i++)
                    {
                        if (r1 == nNeighRankTo[i])
                        {
                            pEN->to_at[j++] = at[to_at].neighbor[i];
                        }
                    }
                    insertions_sort( pCG, pEN->to_at, j, sizeof( pEN->to_at[0] ), CompRanksInvOrd );
                    pEN->num_to = j;  /*  number of stored in pEN->to_at[] central atom neighbors */
                    pEN->from_at = at[from_at].neighbor[neigh_canon_from_min]; /*  neighbor with min. canon number */
                    pEN->rank = r1; /*  mapping rank of the tied neighbors */
                    pEN->canon_rank = r_canon_from_min;  /*  canon. rank of the pEN->from_at */
                }
                else
                {
                    /*  debug only */
                    for (i = j = 0; i < num_neigh; i++)
                    {
                        if (r1 == nNeighRankTo[i])
                        {
                            j++;
                        }
                    }
                }
                /*  debug only */
                if (j <= 1 || !r1 || r_canon_from_min > MAX_ATOMS)
                {
                    return 0; /*  program error */ /*   <BRKPT> */
                }
                return -r; /*  means parity cannot be determined */
            }
            return 0; /* program error */
        }

        /*  All neighbors have different mapping ranks; */
        /*  therefore no additional mapping of the neighbors is necessary */
        if (!ATOM_PARITY_WELL_DEF( at[to_at].parity ))
        {
            return at[to_at].parity; /*  unknown parity or cannot be determined */
        }

        pCG->m_pn_RankForSort = nNeighRankTo;
        num_trans_to = insertions_sort( pCG, nNeighNumberTo, num_neigh, sizeof( nNeighNumberTo[0] ), CompNeighborsRanksCountEql );

        /*  Map canonical ranks of neighbors. Mapped on each other "to" and "from" atoms have equal mapping ranks */
        for (i = 0; i < num_neigh; i++)
        {
            if (nNeighRankTo[j = nNeighNumberTo[i]] != nNeighRankFrom[k = nNeighNumberFrom[i]])
            {
                return 0; /*  program error: mapping ranks not equal, from_at neigborhood cannot be mapped on to_at neighbood. */ /*   <BRKPT> */
            }
            nNeighRankToCanon[j] = nNeighRankFromCanon[k]; /*  potential problem: other atom(s) may have same mapping rank and */
                                                           /*  different canon. rank(s). */
            /*  we may save some memory by eliminating nNeighRankFromCanon[]: */
            /*  nNeighRankToCanon[j] = nCanonRankFrom[at[from_at].neighbor[k]] */
        }

        pCG->m_pn_RankForSort = nNeighRankToCanon;
        num_trans_to += insertions_sort( pCG, nNeighNumberTo, num_neigh, sizeof( nNeighNumberTo[0] ), CompNeighborsRanksCountEql );
    #ifndef CT_NEIGH_INCREASE
        num_trans_to += ( ( num_neigh*( num_neigh - 1 ) ) / 2 ) % 2;  /*  get correct parity for ascending order of canon. numbers */
    #endif

        return 2 - ( num_trans_to + at[to_at].parity ) % 2;
    }
    */
    // END INCHI C FUNCTION: parity_of_mapped_atom2

    let to_atom = heap
        .slice(at.as_const())?
        .get(usize::try_from(to_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let num_neigh = i32::from(to_atom.valence);
    if let Some(en) = pEN.as_deref_mut() {
        *en = EQ_NEIGH::default();
    }
    let rank_from = heap.slice(nRankFrom.as_const())?.to_vec();
    let rank_to = heap.slice(nRankTo.as_const())?.to_vec();
    let rank_at = |ranks: &[AT_RANK], index: i32| -> Result<AT_RANK, SourceHeapError> {
        ranks
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };
    if rank_at(&rank_from, from_at)? != rank_at(&rank_to, to_at)? {
        return Ok(0);
    }
    if num_neigh > MAX_NUM_STEREO_ATOM_NEIGH as i32 || num_neigh < 2 {
        if num_neigh == 1 {
            return Ok(if to_atom.parity != 0 {
                i32::from(to_atom.parity)
            } else {
                AB_PARITY_UNDF as i32
            });
        }
        return Ok(0);
    }

    let from_atom = heap
        .slice(at.as_const())?
        .get(usize::try_from(from_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let canon_from = heap.slice(nCanonRankFrom.as_const())?.to_vec();
    let count = usize::try_from(num_neigh).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut neigh_rank_from = [0_u16; 4];
    let mut neigh_number_from = [0_u16; 4];
    let mut neigh_rank_to = [0_u16; 4];
    let mut neigh_number_to = [0_u16; 4];
    let mut neigh_rank_from_canon = [0_u16; 4];
    let mut neigh_rank_to_canon = [0_u16; 4];
    for i in 0..count {
        neigh_number_from[i] = i as AT_RANK;
        neigh_number_to[i] = i as AT_RANK;
        let to_neighbour = *to_atom
            .neighbor
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let from_neighbour = *from_atom
            .neighbor
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        neigh_rank_to[i] = rank_at(&rank_to, i32::from(to_neighbour))?;
        neigh_rank_from[i] = rank_at(&rank_from, i32::from(from_neighbour))?;
        neigh_rank_from_canon[i] = rank_at(&canon_from, i32::from(from_neighbour))?;
    }

    let from_stack = heap.allocate_model_storage(neigh_rank_from.to_vec())?;
    let to_stack = heap.allocate_model_storage(neigh_rank_to.to_vec())?;
    let canon_stack = heap.allocate_model_storage(neigh_rank_to_canon.to_vec())?;
    let result = (|| -> Result<i32, SourceHeapError> {
        pCG.m_pn_RankForSort = from_stack.as_const();
        pCG.m_nNumCompNeighborsRanksCountEql = 0;
        let _num_trans_from = insertions_sort(
            bytemuck::cast_slice_mut(&mut neigh_number_from),
            count,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                CompNeighborsRanksCountEql(left, right, &neigh_rank_from, pCG)
            },
        )?;

        if pCG.m_nNumCompNeighborsRanksCountEql != 0 {
            let mut r_canon_from_min = MAX_ATOMS as i32 + 1;
            let mut neigh_canon_from_min = 0_usize;
            let mut r = 0_i32;
            let mut neigh1 = usize::from(neigh_number_from[0]);
            let mut r1 = i32::from(neigh_rank_from[neigh1]);
            let mut i = 1_usize;
            while i < count {
                let neigh2 = usize::from(neigh_number_from[i]);
                let r2 = i32::from(neigh_rank_from[neigh2]);
                if r2 == r1 {
                    if r != r2 {
                        r = r2;
                        let r_canon_from = i32::from(neigh_rank_from_canon[neigh1]);
                        if r_canon_from < r_canon_from_min {
                            r_canon_from_min = r_canon_from;
                            neigh_canon_from_min = neigh1;
                        }
                    }
                    let r_canon_from = i32::from(neigh_rank_from_canon[neigh2]);
                    if r_canon_from < r_canon_from_min {
                        r_canon_from_min = r_canon_from;
                        neigh_canon_from_min = neigh2;
                    }
                }
                r1 = r2;
                neigh1 = neigh2;
                i += 1;
            }
            if r != 0 {
                r1 = i32::from(neigh_rank_from[neigh_canon_from_min]);
                let mut j = 0_usize;
                if let Some(en) = pEN.as_deref_mut() {
                    for i in 0..count {
                        if r1 == i32::from(neigh_rank_to[i]) {
                            en.to_at[j] = to_atom.neighbor[i];
                            j += 1;
                        }
                    }
                    let _ = insertions_sort(
                        bytemuck::cast_slice_mut(&mut en.to_at),
                        j,
                        std::mem::size_of::<AT_RANK>(),
                        &mut |left, right| {
                            let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                            let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                            Ok(CompRanksInvOrd(left, right))
                        },
                    )?;
                    en.num_to = j as i32;
                    en.from_at = from_atom.neighbor[neigh_canon_from_min];
                    en.rank = r1 as AT_RANK;
                    en.canon_rank = r_canon_from_min as AT_RANK;
                } else {
                    for &rank in neigh_rank_to.iter().take(count) {
                        if r1 == i32::from(rank) {
                            j += 1;
                        }
                    }
                }
                if j <= 1 || r1 == 0 || r_canon_from_min > MAX_ATOMS as i32 {
                    return Ok(0);
                }
                return Ok(r.wrapping_neg());
            }
            return Ok(0);
        }

        if !(AB_MIN_WELL_DEFINED_PARITY as i8..=AB_MAX_WELL_DEFINED_PARITY as i8)
            .contains(&to_atom.parity)
        {
            return Ok(i32::from(to_atom.parity));
        }

        pCG.m_pn_RankForSort = to_stack.as_const();
        let mut num_trans_to = insertions_sort(
            bytemuck::cast_slice_mut(&mut neigh_number_to),
            count,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                CompNeighborsRanksCountEql(left, right, &neigh_rank_to, pCG)
            },
        )?;
        for i in 0..count {
            let j = usize::from(neigh_number_to[i]);
            let k = usize::from(neigh_number_from[i]);
            if neigh_rank_to[j] != neigh_rank_from[k] {
                return Ok(0);
            }
            neigh_rank_to_canon[j] = neigh_rank_from_canon[k];
        }
        heap.slice_mut(canon_stack)?
            .get_mut(..4)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .copy_from_slice(&neigh_rank_to_canon);
        pCG.m_pn_RankForSort = canon_stack.as_const();
        num_trans_to = num_trans_to.wrapping_add(insertions_sort(
            bytemuck::cast_slice_mut(&mut neigh_number_to),
            count,
            std::mem::size_of::<AT_RANK>(),
            &mut |left, right| {
                let left = AT_RANK::from_ne_bytes([left[0], left[1]]);
                let right = AT_RANK::from_ne_bytes([right[0], right[1]]);
                CompNeighborsRanksCountEql(left, right, &neigh_rank_to_canon, pCG)
            },
        )?);
        Ok(2_i32.wrapping_sub(num_trans_to.wrapping_add(i32::from(to_atom.parity)) % 2))
    })();

    let free_from = heap.free(from_stack);
    let free_to = heap.free(to_stack);
    let free_canon = heap.free(canon_stack);
    let cleanup = free_from.and(free_to).and(free_canon);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[allow(non_snake_case)]
pub(crate) fn might_change_other_atom_parity(
    heap: &SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    at_no: i32,
    nRank2: SourceMutPointer<AT_RANK>,
    nRank1: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1461 might_change_other_atom_parity
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int might_change_other_atom_parity( sp_ATOM *at,
                                        int num_atoms,
                                        int at_no,
                                        AT_RANK *nRank2,
                                        AT_RANK *nRank1 )
    {
        int     i, j, neighbor_no;
        for (i = 0; i < num_atoms; i++)
        {
            if (nRank2[i] != nRank1[i])
            {
                if (i != at_no /*&& ATOM_PARITY_WELL_DEF(at[i].parity)*/
                    && at[i].bHasStereoOrEquToStereo
                    && !( at[i].stereo_atom_parity & KNOWN_PARITIES_EQL )
                    && !at[i].stereo_bond_neighbor[0]
                    )
                {

                    return 1; /*  may have changed stereo atoms order */
                }
                for (j = 0; j < at[i].valence; j++)
                {
                    neighbor_no = at[i].neighbor[j];
                    if (neighbor_no != at_no
                         /*&& ATOM_PARITY_WELL_DEF(at[neighbor_no].parity)*/
                         && at[neighbor_no].bHasStereoOrEquToStereo
                         && !( at[neighbor_no].stereo_atom_parity & KNOWN_PARITIES_EQL )
                         && !at[neighbor_no].stereo_bond_neighbor[0]
                       )
                    {
                        return 1; /*  may have changed stereo atom parity */
                    }
                }
            }
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: might_change_other_atom_parity

    if num_atoms <= 0 {
        return Ok(0);
    }
    let rank2 = heap.slice(nRank2.as_const())?;
    let rank1 = heap.slice(nRank1.as_const())?;
    let mut i = 0_i32;
    while i < num_atoms {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if rank2
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != rank1
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            let atoms = heap.slice(at.as_const())?;
            let atom = atoms
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i != at_no
                && atom.bHasStereoOrEquToStereo != 0
                && (i32::from(atom.stereo_atom_parity) & KNOWN_PARITIES_EQL as i32) == 0
                && atom.stereo_bond_neighbor[0] == 0
            {
                return Ok(1);
            }
            let mut j = 0_i32;
            while j < i32::from(atom.valence) {
                let neighbour = *atom
                    .neighbor
                    .get(usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let neighbour_atom = atoms
                    .get(usize::from(neighbour))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if i32::from(neighbour) != at_no
                    && neighbour_atom.bHasStereoOrEquToStereo != 0
                    && (i32::from(neighbour_atom.stereo_atom_parity) & KNOWN_PARITIES_EQL as i32)
                        == 0
                    && neighbour_atom.stereo_bond_neighbor[0] == 0
                {
                    return Ok(1);
                }
                j = j.wrapping_add(1);
            }
        }
        i = i.wrapping_add(1);
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn DeAllocateForNonStereoRemoval(
    heap: &mut SourceHeap,
    nAtomNumberCanon1: &mut SourceMutPointer<AT_RANK>,
    nAtomNumberCanon2: &mut SourceMutPointer<AT_RANK>,
    nl: &mut SourceMutPointer<NEIGH_LIST>,
    nl1: &mut SourceMutPointer<NEIGH_LIST>,
    nl2: &mut SourceMutPointer<NEIGH_LIST>,
    nVisited1: &mut SourceMutPointer<AT_RANK>,
    nVisited2: &mut SourceMutPointer<AT_RANK>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1504 DeAllocateForNonStereoRemoval
    // INCHI✔️❌: complete active REMOVE_CALC_NONSTEREO == 1 source frame follows verbatim.
    /*
    void DeAllocateForNonStereoRemoval( AT_RANK **nAtomNumberCanon1,
                                        AT_RANK **nAtomNumberCanon2,
                                        NEIGH_LIST **nl,
                                        NEIGH_LIST **nl1,
                                        NEIGH_LIST **nl2,
                                        AT_RANK **nVisited1,
                                        AT_RANK **nVisited2 )
    {
        if (*nAtomNumberCanon1)
        {
            inchi_free( *nAtomNumberCanon1 );
            *nAtomNumberCanon1 = NULL;
        }
        if (*nAtomNumberCanon2)
        {
            inchi_free( *nAtomNumberCanon2 );
            *nAtomNumberCanon2 = NULL;
        }
        if (*nl)
        {
            FreeNeighList( *nl );
            *nl = 0;
        }
        if (*nl1)
        {
            FreeNeighList( *nl1 );
            *nl1 = 0;
        }
        if (*nl2)
        {
            FreeNeighList( *nl2 );
            *nl2 = 0;
        }
        if (*nVisited1)
        {
            inchi_free( *nVisited1 );
            *nVisited1 = NULL;
        }
        if (*nVisited2)
        {
            inchi_free( *nVisited2 );
            *nVisited2 = NULL;
        }
    }
    */
    // END INCHI C FUNCTION: DeAllocateForNonStereoRemoval

    if !nAtomNumberCanon1.is_null() {
        inchi_free(heap, *nAtomNumberCanon1)?;
        *nAtomNumberCanon1 = SourceMutPointer::null();
    }
    if !nAtomNumberCanon2.is_null() {
        inchi_free(heap, *nAtomNumberCanon2)?;
        *nAtomNumberCanon2 = SourceMutPointer::null();
    }
    if !nl.is_null() {
        FreeNeighList(heap, *nl)?;
        *nl = SourceMutPointer::null();
    }
    if !nl1.is_null() {
        FreeNeighList(heap, *nl1)?;
        *nl1 = SourceMutPointer::null();
    }
    if !nl2.is_null() {
        FreeNeighList(heap, *nl2)?;
        *nl2 = SourceMutPointer::null();
    }
    if !nVisited1.is_null() {
        inchi_free(heap, *nVisited1)?;
        *nVisited1 = SourceMutPointer::null();
    }
    if !nVisited2.is_null() {
        inchi_free(heap, *nVisited2)?;
        *nVisited2 = SourceMutPointer::null();
    }
    Ok(())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AllocateForNonStereoRemoval(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon1: &mut SourceMutPointer<AT_RANK>,
    nAtomNumberCanon2: &mut SourceMutPointer<AT_RANK>,
    nl: &mut SourceMutPointer<NEIGH_LIST>,
    nl1: &mut SourceMutPointer<NEIGH_LIST>,
    nl2: &mut SourceMutPointer<NEIGH_LIST>,
    nVisited1: &mut SourceMutPointer<AT_RANK>,
    nVisited2: &mut SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1551 AllocateForNonStereoRemoval
    // INCHI✔️❌: complete active REMOVE_CALC_NONSTEREO == 1 source frame follows verbatim.
    /*
    int AllocateForNonStereoRemoval( sp_ATOM *at,
                                     int num_atoms,
                                     const AT_RANK *nSymmRank,
                                     AT_RANK *nCanonRank,
                                     AT_RANK **nAtomNumberCanon1,
                                     AT_RANK **nAtomNumberCanon2,
                                     NEIGH_LIST **nl,
                                     NEIGH_LIST **nl1,
                                     NEIGH_LIST **nl2,
                                     AT_RANK **nVisited1,
                                     AT_RANK **nVisited2 )
    {
        DeAllocateForNonStereoRemoval( nAtomNumberCanon1, nAtomNumberCanon2, nl, nl1, nl2, nVisited1, nVisited2 );
        *nAtomNumberCanon1 = (AT_RANK *) inchi_malloc( num_atoms * sizeof( **nAtomNumberCanon1 ) );
        *nAtomNumberCanon2 = (AT_RANK *) inchi_malloc( num_atoms * sizeof( **nAtomNumberCanon2 ) );
        *nl = CreateNeighList( num_atoms, num_atoms, at, 0, NULL );
        *nl1 = CreateNeighList( num_atoms, num_atoms, at, 0, NULL );
        *nl2 = CreateNeighList( num_atoms, num_atoms, at, 0, NULL );
        *nVisited1 = (AT_RANK *) inchi_malloc( num_atoms * sizeof( **nVisited1 ) );
        *nVisited2 = (AT_RANK *) inchi_malloc( num_atoms * sizeof( **nVisited2 ) );

        if (!*nl || !*nl1 || !*nl2 || !*nVisited1 || !*nVisited2 || !*nAtomNumberCanon1 || !*nAtomNumberCanon2)
        {
            DeAllocateForNonStereoRemoval( nAtomNumberCanon1, nAtomNumberCanon2, nl, nl1, nl2, nVisited1, nVisited2 );
            return 0;
        }
        /*  Sort neighbors according to symm. ranks (primary key) and canon. ranks (secondary key), in descending order */
        SortNeighListsBySymmAndCanonRank( num_atoms, *nl, nSymmRank, nCanonRank );
        SortNeighListsBySymmAndCanonRank( num_atoms, *nl1, nSymmRank, nCanonRank );
        SortNeighListsBySymmAndCanonRank( num_atoms, *nl2, nSymmRank, nCanonRank );

        return 1;
    }
    */
    // END INCHI C FUNCTION: AllocateForNonStereoRemoval

    DeAllocateForNonStereoRemoval(
        heap,
        nAtomNumberCanon1,
        nAtomNumberCanon2,
        nl,
        nl1,
        nl2,
        nVisited1,
        nVisited2,
    )?;
    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut allocate_rank = |heap: &mut SourceHeap| match heap.allocate(vec![0_u16; count]) {
        Ok(pointer) => Ok(pointer),
        Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
        Err(error) => Err(error),
    };
    *nAtomNumberCanon1 = allocate_rank(heap)?;
    *nAtomNumberCanon2 = allocate_rank(heap)?;
    *nl = CreateNeighList(heap, num_atoms, num_atoms, at.as_const(), 0, None)?;
    *nl1 = CreateNeighList(heap, num_atoms, num_atoms, at.as_const(), 0, None)?;
    *nl2 = CreateNeighList(heap, num_atoms, num_atoms, at.as_const(), 0, None)?;
    *nVisited1 = allocate_rank(heap)?;
    *nVisited2 = allocate_rank(heap)?;

    if nl.is_null()
        || nl1.is_null()
        || nl2.is_null()
        || nVisited1.is_null()
        || nVisited2.is_null()
        || nAtomNumberCanon1.is_null()
        || nAtomNumberCanon2.is_null()
    {
        DeAllocateForNonStereoRemoval(
            heap,
            nAtomNumberCanon1,
            nAtomNumberCanon2,
            nl,
            nl1,
            nl2,
            nVisited1,
            nVisited2,
        )?;
        return Ok(0);
    }
    SortNeighListsBySymmAndCanonRank(heap, num_atoms, *nl, nSymmRank, nCanonRank)?;
    SortNeighListsBySymmAndCanonRank(heap, num_atoms, *nl1, nSymmRank, nCanonRank)?;
    SortNeighListsBySymmAndCanonRank(heap, num_atoms, *nl2, nSymmRank, nCanonRank)?;
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn GetMinNewRank(
    heap: &SourceHeap,
    nAtomRank: SourceMutPointer<AT_RANK>,
    nAtomNumb: SourceMutPointer<AT_RANK>,
    nRank1: AT_RANK,
) -> Result<AT_RANK, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1587 GetMinNewRank
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    AT_RANK GetMinNewRank( AT_RANK *nAtomRank, AT_RANK *nAtomNumb, AT_RANK nRank1 )
    {
        int i;
        AT_RANK nRank2 = 0;
        for (i = (int) nRank1 - 1; 0 <= i && nRank1 == ( nRank2 = nAtomRank[(int) nAtomNumb[i]] ); i--)
        {
            ;
        }
        if (i >= 0)
        {
            nRank2++;
        }
        else
        {
            nRank2 = 1;
        }

        return nRank2;
    }
    */
    // END INCHI C FUNCTION: GetMinNewRank

    let mut i = i32::from(nRank1).wrapping_sub(1);
    let mut nRank2 = 0_u16;
    while i >= 0 {
        let atom = source_get(heap, nAtomNumb, i)?;
        nRank2 = source_get(heap, nAtomRank, i32::from(atom))?;
        if nRank1 != nRank2 {
            break;
        }
        i = i.wrapping_sub(1);
    }
    if i >= 0 {
        nRank2 = nRank2.wrapping_add(1);
    } else {
        nRank2 = 1;
    }
    Ok(nRank2)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn BreakNeighborsTie(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    num_at_tg: i32,
    ib: i32,
    ia: i32,
    neigh_num: SourceMutPointer<AT_RANK>,
    in1: i32,
    in2: i32,
    mode: i32,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nTempRank: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nl1: SourceMutPointer<NEIGH_LIST>,
    nl2: SourceMutPointer<NEIGH_LIST>,
    lNumIter: &mut i64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1609 BreakNeighborsTie
    // INCHI✔️❌: int BreakNeighborsTie( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                        sp_ATOM *at,
    // INCHI✔️❌:                        int num_atoms,
    // INCHI✔️❌:                        int num_at_tg,
    // INCHI✔️❌:                        int ib,
    // INCHI✔️❌:                        int ia,
    // INCHI✔️❌:                        AT_RANK *neigh_num,
    // INCHI✔️❌:                        int in1,
    // INCHI✔️❌:                        int in2,
    // INCHI✔️❌:                        int mode,
    // INCHI✔️❌:                        AT_RANK **pRankStack1,
    // INCHI✔️❌:                        AT_RANK **pRankStack2,
    // INCHI✔️❌:                        AT_RANK *nTempRank,
    // INCHI✔️❌:                        NEIGH_LIST *NeighList,
    // INCHI✔️❌:                        const AT_RANK *nSymmRank,
    // INCHI✔️❌:                        AT_RANK *nCanonRank,
    // INCHI✔️❌:                        NEIGH_LIST *nl1,
    // INCHI✔️❌:                        NEIGH_LIST *nl2,
    // INCHI✔️❌:                        long *lNumIter )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_RANK nRank1, nRank2;
    // INCHI✔️❌:     int     nNumDiffRanks, nNumDiffRanks1, nNumDiffRanks2, i;
    // INCHI✔️❌:     int n1 = (int) neigh_num[in1];
    // INCHI✔️❌:     int n2 = (int) neigh_num[in2];
    // INCHI✔️❌:     int other_neigh[2] = {0}, other_neig_ord[2] = {0}, num_other_neigh; /* djb-rwth: initialisations added */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  asymmetric calculation */
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((mode == MAP_MODE_S4  && in1) || /* for S4 we need only (in1,in2) = (0,1) (0,2) (0,3) pairs of neighbors */
    // INCHI✔️❌:          (mode != MAP_MODE_STD && at[ia].valence != MAX_NUM_STEREO_ATOM_NEIGH) ||
    // INCHI✔️❌:          (mode != MAP_MODE_STD && nSymmRank[n1] != nSymmRank[n2])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  1. Create initial ranks from equivalence information stored in nSymmRank */
    // INCHI✔️❌:     memcpy(pRankStack1[0], nSymmRank, num_at_tg * sizeof(pRankStack1[0][0]));
    // INCHI✔️❌:     pCG->m_pn_RankForSort = pRankStack1[0];
    // INCHI✔️❌:     tsort( pCG, pRankStack1[1], num_at_tg, sizeof( pRankStack1[1][0] ), CompRanksOrd );
    // INCHI✔️❌:     nNumDiffRanks = SortedEquInfoToRanks( pRankStack1[0]/*inp*/, pRankStack1[0]/*out*/, pRankStack1[1], num_at_tg, NULL );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* other neighbors */
    // INCHI✔️❌:     num_other_neigh = 0;
    // INCHI✔️❌:     if (at[ia].valence <= MAX_NUM_STEREO_ATOM_NEIGH && mode)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < at[ia].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i != in1 && i != in2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 other_neigh[num_other_neigh] = (int) neigh_num[i];
    // INCHI✔️❌:                 other_neig_ord[num_other_neigh] = i;
    // INCHI✔️❌:                 num_other_neigh++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((mode != MAP_MODE_STD && nSymmRank[other_neigh[0]] != nSymmRank[other_neigh[1]]) ||
    // INCHI✔️❌:          (mode == MAP_MODE_S4 && nSymmRank[n1] != nSymmRank[other_neigh[1]])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  2. Fix at[ia] */
    // INCHI✔️❌:     if (pRankStack1[0][ia] != nSymmRank[ia])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  at[ia] is constitutionally equivalent to some other atom. Fix at[ia]. */
    // INCHI✔️❌:         pRankStack1[0][ia] = nSymmRank[ia];
    // INCHI✔️❌:         nNumDiffRanks = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                      nNumDiffRanks, pRankStack1[0], nTempRank,
    // INCHI✔️❌:                                      pRankStack1[1], lNumIter, 1 );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*  3. In case of a double bond/cumulene only: */
    // INCHI✔️❌:     /*     fix at[ib] -- the opposite double bond/cumulene atom */
    // INCHI✔️❌:     if (ib < num_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  find the smallest possible rank */
    // INCHI✔️❌:         nRank1 = pRankStack1[0][ib];
    // INCHI✔️❌:         nRank2 = GetMinNewRank( pRankStack1[0], pRankStack1[1], nRank1 );
    // INCHI✔️❌:         /*  if the rank is smaller than pRankStack1[0][ib] then fix at[ib] */
    // INCHI✔️❌:         if (nRank2 != nRank1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pRankStack1[0][ib] = nRank2;
    // INCHI✔️❌:             nNumDiffRanks = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                          nNumDiffRanks, pRankStack1[0], nTempRank,
    // INCHI✔️❌:                                          pRankStack1[1], lNumIter, 1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**************************************************************************************
    // INCHI✔️❌:      * Note: It may (or may not?) make sense to fix "other neighbors":
    // INCHI✔️❌:      *       in case of a stereo center fix neighbors other than n1, n2
    // INCHI✔️❌:      *       in case of a double bond/cumulene fix the opposite atom neighbors
    // INCHI✔️❌:      *       The ranks assigned to the other neighbors in case of their equivalence
    // INCHI✔️❌:      *       should be in the ascending order of their canonical ranks ????
    // INCHI✔️❌:      *       *** For now we do not fix other neighbors ***
    // INCHI✔️❌:      **************************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  4. Check whether the neighbors still have equal ranks */
    // INCHI✔️❌:     if (pRankStack1[0][n1] != pRankStack1[0][n2])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /*  the two neighbors are not constitutionally equivalent */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*  5. Find new smallest possible rank for n1 and n2 */
    // INCHI✔️❌:     nRank1 = pRankStack1[0][n1];
    // INCHI✔️❌:     nRank2 = GetMinNewRank( pRankStack1[0], pRankStack1[1], nRank1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  6. Copy the results to the 2nd eq. rank arrays */
    // INCHI✔️❌:     memcpy(pRankStack2[0], pRankStack1[0], num_at_tg * sizeof(pRankStack2[0][0]));
    // INCHI✔️❌:     memcpy(pRankStack2[1], pRankStack1[1], num_at_tg * sizeof(pRankStack2[0][0]));
    // INCHI✔️❌:     /*  7. Break neighbor tie: map n1(1) <--> n2(2) */
    // INCHI✔️❌:     /*  7. Break neighbor tie: map n1(1) <--> n2(2) */
    // INCHI✔️❌:     pRankStack1[0][n1] = nRank2;
    // INCHI✔️❌:     nNumDiffRanks1 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                  nNumDiffRanks, pRankStack1[0], nTempRank,
    // INCHI✔️❌:                                  pRankStack1[1], lNumIter, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     pRankStack2[0][n2] = nRank2;
    // INCHI✔️❌:     nNumDiffRanks2 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                  nNumDiffRanks, pRankStack2[0], nTempRank,
    // INCHI✔️❌:                                  pRankStack2[1], lNumIter, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumDiffRanks1 != nNumDiffRanks2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  <BRKPT> */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (mode == MAP_MODE_C2v || mode == MAP_MODE_C2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Check for C2v reflection leading to parity inversion (mode=1) or C2 rotation (mode=2) */
    // INCHI✔️❌:         AT_RANK nRank10, nRank20;
    // INCHI✔️❌:         int     nn1, nn2;
    // INCHI✔️❌:         /*
    // INCHI✔️❌:          * C2v & C2: map
    // INCHI✔️❌:          * n1(1) <--> n2(2) -- at this point already done
    // INCHI✔️❌:          * n1(2) <--> n2(1) --> do at i = 0
    // INCHI✔️❌:          *
    // INCHI✔️❌:          * C2v: other neighbors must be unmoved: map
    // INCHI✔️❌:          * other_neigh[0](1) <--> other_neigh[0](2)
    // INCHI✔️❌:          * other_neigh[1](1) <--> other_neigh[1](2)
    // INCHI✔️❌:          *
    // INCHI✔️❌:          * C2:  other neighbors should be mapped on each other
    // INCHI✔️❌:          * other_neigh[0](1) <--> other_neigh[1](2)
    // INCHI✔️❌:          * other_neigh[1](1) <--> other_neigh[0](2)
    // INCHI✔️❌:          */
    // INCHI✔️❌:         for (i = 0; i <= 2; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* C2v & C2. Map n2(1) <--> n1(2) */
    // INCHI✔️❌:                 nn1 = n2;
    // INCHI✔️❌:                 nn2 = n1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (mode == MAP_MODE_C2v)
    // INCHI✔️❌:                 {   /* was '=', pointed by WDI */
    // INCHI✔️❌:                     /* i = 1 or 2
    // INCHI✔️❌:                      * C2v. Other neighbors must be unmoved: map
    // INCHI✔️❌:                      * i=1: other_neigh[0](1) <--> other_neigh[0](2)
    // INCHI✔️❌:                      * i=2: other_neigh[1](1) <--> other_neigh[1](2)
    // INCHI✔️❌:                      */
    // INCHI✔️❌:                     nn1 = other_neigh[i - 1]; /* 0 or 1 */
    // INCHI✔️❌:                     nn2 = other_neigh[i - 1]; /* 0 or 1 */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (mode == MAP_MODE_C2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* was '=', pointed by WDI */
    // INCHI✔️❌:                         /* i = 1 or 2
    // INCHI✔️❌:                          * C2.  Other neighbors should be mapped on each other
    // INCHI✔️❌:                          * i=1: other_neigh[0](1) <--> other_neigh[1](2)
    // INCHI✔️❌:                          * i=2: other_neigh[1](1) <--> other_neigh[0](2)
    // INCHI✔️❌:                          */
    // INCHI✔️❌:                         nn1 = other_neigh[i - 1]; /* 0 or 1 */
    // INCHI✔️❌:                         nn2 = other_neigh[2 - i]; /* 1 or 0 */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1; /* program error */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* map nn1(1) <--> nn2(2) */
    // INCHI✔️❌:             nRank10 = pRankStack1[0][nn1];
    // INCHI✔️❌:             nRank20 = pRankStack2[0][nn2];
    // INCHI✔️❌:             nRank1 = GetMinNewRank( pRankStack1[0], pRankStack1[1], nRank10 );
    // INCHI✔️❌:             nRank2 = GetMinNewRank( pRankStack2[0], pRankStack2[1], nRank20 );
    // INCHI✔️❌:             if (nRank10 == nRank20 && nRank1 == nRank2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nRank10 == nRank1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ; /* atoms are already mapped */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* need additional mapping: ranks are not fixed yet */
    // INCHI✔️❌:                     pRankStack1[0][nn1] = nRank1;
    // INCHI✔️❌:                     nNumDiffRanks1 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                                  nNumDiffRanks, pRankStack1[0], nTempRank,
    // INCHI✔️❌:                                                  pRankStack1[1], lNumIter, 1 );
    // INCHI✔️❌:                     pRankStack2[0][nn2] = nRank2;
    // INCHI✔️❌:                     nNumDiffRanks2 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                                  nNumDiffRanks, pRankStack2[0], nTempRank,
    // INCHI✔️❌:                                                  pRankStack2[1], lNumIter, 1 );
    // INCHI✔️❌:                     if (nNumDiffRanks1 != nNumDiffRanks2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1; /*  <BRKPT> */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;  /* mapping is not possible */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (mode == MAP_MODE_S4)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*
    // INCHI✔️❌:          *  Check for S4 reflection/rotation leading to parity inversion (mode=3)
    // INCHI✔️❌:          *
    // INCHI✔️❌:          * At this point n1(1) <--> n2(2) have been mapped and n1 has index in1 = 0
    // INCHI✔️❌:          * Below indexes in neigh_num[] are in brackets; [i] means neigh_num[i].
    // INCHI✔️❌:          * Numbers (#) in parentheses refer to pRankStack#
    // INCHI✔️❌:          *
    // INCHI✔️❌:          * in2=1: [0](1) <--> [1](2)  mapping has been done; add more mappings:
    // INCHI✔️❌:          *        [1](1) <--> [2](2)  [x]=[2]
    // INCHI✔️❌:          *        [2](1) <--> [3](2)  [y]=[3]
    // INCHI✔️❌:          *        [3](1) <--> [0](2)
    // INCHI✔️❌:          *        this will succeed if C2 axis crosses middle of [0]-[2] and [1]-[3] lines
    // INCHI✔️❌:          *
    // INCHI✔️❌:          * in2=2: [0](1) <--> [2](2) mapping has been done; add more mappings:
    // INCHI✔️❌:          *        [2](1) <--> [3](2)  [x]=[3]
    // INCHI✔️❌:          *        [3](1) <--> [1](2)  [y]=[1]
    // INCHI✔️❌:          *        [1](1) <--> [0](2)
    // INCHI✔️❌:          *        this will succeed if C2 axis crosses middle of [0]-[3] and [1]-[2] lines
    // INCHI✔️❌:          *
    // INCHI✔️❌:          * in2=3: [0](1) <--> [3](2) mapping has been done; add more mappings:
    // INCHI✔️❌:          *        [3](1) <--> [1](2)  [x]=[1]
    // INCHI✔️❌:          *        [1](1) <--> [2](2)  [y]=[2]
    // INCHI✔️❌:          *        [2](1) <--> [0](2)
    // INCHI✔️❌:          *        this will succeed if C2 axis crosses middle of [0]-[1] and [2]-[3] lines
    // INCHI✔️❌:          *
    // INCHI✔️❌:          * In general:
    // INCHI✔️❌:          *        [in1](1) <--> [in2](2)
    // INCHI✔️❌:          *        [in2](1) <--> [x]  (2)  i=0
    // INCHI✔️❌:          *        [x]  (1) <--> [y]  (2)  i=1
    // INCHI✔️❌:          *        [y]  (1) <--> [in1](2)  i=2
    // INCHI✔️❌:          *
    // INCHI✔️❌:          *    in1=0    always
    // INCHI✔️❌:          *    ===== how to find x, y from in2 ====
    // INCHI✔️❌:          *    in2=1 => x,y = 2, 3  or [x] = other_neigh[0], [y] = other_neigh[1]
    // INCHI✔️❌:          *    in2=2 => x,y = 3, 1  or [x] = other_neigh[1], [y] = other_neigh[0]
    // INCHI✔️❌:          *    in2=3 => x,y = 1, 2  or [x] = other_neigh[0], [y] = other_neigh[1]
    // INCHI✔️❌:          *    ====================================
    // INCHI✔️❌:          */
    // INCHI✔️❌:         AT_RANK nRank10, nRank20;
    // INCHI✔️❌:         int     nn1, nn2;
    // INCHI✔️❌:         for (i = 0; i <= 2; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             switch (i)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case 0:  /* [in2](1) <--> [x](2);  */
    // INCHI✔️❌:                     nn1 = n2;                    /* [in2] */
    // INCHI✔️❌:                     nn2 = other_neigh[1 - in2 % 2];  /* [x]   */
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case 1:  /* [x](1) <--> [y](2) */
    // INCHI✔️❌:                     nn1 = other_neigh[1 - in2 % 2];  /* [x]   */
    // INCHI✔️❌:                     nn2 = other_neigh[in2 % 2];  /* [y]   */
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case 2:
    // INCHI✔️❌:                     nn1 = other_neigh[in2 % 2];  /* [y]   */
    // INCHI✔️❌:                     nn2 = n1;                    /* [in1] */
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 default:
    // INCHI✔️❌:                     return -1; /* program error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* map nn1(1) <--> nn2(2) */
    // INCHI✔️❌:             nRank10 = pRankStack1[0][nn1];
    // INCHI✔️❌:             nRank20 = pRankStack2[0][nn2];
    // INCHI✔️❌:             nRank1 = GetMinNewRank( pRankStack1[0], pRankStack1[1], nRank10 );
    // INCHI✔️❌:             nRank2 = GetMinNewRank( pRankStack2[0], pRankStack2[1], nRank20 );
    // INCHI✔️❌:             if (nRank10 == nRank20 && nRank1 == nRank2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nRank10 == nRank1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ;/* atoms are already mapped */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                              /* need additional mapping: ranks are not fixed yet */
    // INCHI✔️❌:                     pRankStack1[0][nn1] = nRank1;
    // INCHI✔️❌:                     nNumDiffRanks1 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                                  nNumDiffRanks, pRankStack1[0], nTempRank,
    // INCHI✔️❌:                                                  pRankStack1[1], lNumIter, 1 );
    // INCHI✔️❌:                     pRankStack2[0][nn2] = nRank2;
    // INCHI✔️❌:                     nNumDiffRanks2 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                                  nNumDiffRanks, pRankStack2[0], nTempRank,
    // INCHI✔️❌:                                                  pRankStack2[1], lNumIter, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (nNumDiffRanks1 != nNumDiffRanks2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1; /*  <BRKPT> */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;  /* mapping is not possible */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( BREAK_ONE_MORE_SC_TIE == 1 ) /* { */
    // INCHI✔️❌:     /* Check for a very highly symmetrical stereo center 12-06-2002 */
    // INCHI✔️❌:     if (ib >= num_atoms && at[ia].valence == MAX_NUM_STEREO_ATOM_NEIGH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int num_eq;
    // INCHI✔️❌:         nRank1 = pRankStack1[0][n2];
    // INCHI✔️❌:         for (i = 0, num_eq = 0; i < at[ia].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_eq += ( nRank1 == pRankStack1[0][at[ia].neighbor[i]] );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (num_eq == MAX_NUM_STEREO_ATOM_NEIGH - 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = (int) nRank1 - 1; 0 <= i && nRank1 == ( nRank2 = pRankStack1[0][(int) pRankStack1[1][i]] ); i--)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nRank2++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nRank2 = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /*  7a. Break another neighbor tie */
    // INCHI✔️❌:             nNumDiffRanks = nNumDiffRanks1;
    // INCHI✔️❌:
    // INCHI✔️❌:             pRankStack1[0][n2] = nRank2;
    // INCHI✔️❌:             nNumDiffRanks1 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                          nNumDiffRanks, pRankStack1[0], nTempRank,
    // INCHI✔️❌:                                          pRankStack1[1], lNumIter, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:             pRankStack2[0][n1] = nRank2;
    // INCHI✔️❌:             nNumDiffRanks2 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                          nNumDiffRanks, pRankStack2[0], nTempRank,
    // INCHI✔️❌:                                          pRankStack2[1], lNumIter, 1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumDiffRanks1 != nNumDiffRanks2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  <BRKPT> */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /* } BREAK_ONE_MORE_SC_TIE */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( BREAK_ALSO_NEIGH_TIE == 1 )
    // INCHI✔️❌:     /* check whether neighbor's neighbors are tied and untie them */
    // INCHI✔️❌:     if (at[n1].nRingSystem == at[n2].nRingSystem &&  ib >= num_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AT_RANK NeighNeighList[MAX_NUM_STEREO_ATOM_NEIGH + 1];
    // INCHI✔️❌:         int m, neigh1 = -1, neigh2 = -1;
    // INCHI✔️❌:         nRank1 = nRank2 = 0;
    // INCHI✔️❌:         /* n1 */
    // INCHI✔️❌:         NeighNeighList[0] = at[n1].valence - 1; /* for insertions_sort_NeighListBySymmAndCanonRank() */
    // INCHI✔️❌:         for (i = 0, m = 1; i < at[n1].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int neigh = at[n1].neighbor[i];
    // INCHI✔️❌:             if (neigh != ia)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 NeighNeighList[m++] = neigh;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         insertions_sort_NeighListBySymmAndCanonRank( NeighNeighList, pRankStack1[0], nCanonRank );
    // INCHI✔️❌:         for (m = 2; m < at[n1].valence; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pRankStack1[0][NeighNeighList[m]] == pRankStack1[0][NeighNeighList[m - 1]])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 neigh1 = NeighNeighList[m - 1];
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* n2 */
    // INCHI✔️❌:         NeighNeighList[0] = at[n2].valence - 1; /* for insertions_sort_NeighListBySymmAndCanonRank() */
    // INCHI✔️❌:         for (i = 0, m = 1; i < at[n2].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int neigh = at[n2].neighbor[i];
    // INCHI✔️❌:             if (neigh != ia)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 NeighNeighList[m++] = neigh;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         insertions_sort_NeighListBySymmAndCanonRank( NeighNeighList, pRankStack2[0], nCanonRank );
    // INCHI✔️❌:         for (m = 2; m < at[n2].valence; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pRankStack2[0][NeighNeighList[m]] == pRankStack2[0][NeighNeighList[m - 1]])
    // INCHI✔️❌:             {
    // INCHI✔️❌: #if ( BREAK_ALSO_NEIGH_TIE_ROTATE == 1 )
    // INCHI✔️❌:                 neigh2 = NeighNeighList[m];    /* [m] to obtain same axis orientation  around ia<neigh */
    // INCHI✔️❌: #else
    // INCHI✔️❌:                 neigh2 = NeighNeighList[m - 1];  /* [m-1] to obtain reflection ??? */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (neigh1 >= 0 && neigh2 >= 0 && pRankStack1[0][neigh1] == pRankStack2[0][neigh2])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* neighbors' neighbors are tied */
    // INCHI✔️❌:             nRank1 = pRankStack1[0][neigh1];
    // INCHI✔️❌:             nRank2 = GetMinNewRank( pRankStack1[0], pRankStack1[1], nRank1 );
    // INCHI✔️❌:
    // INCHI✔️❌:             /*  Break neighbor's neighbor tie */
    // INCHI✔️❌:
    // INCHI✔️❌:             nNumDiffRanks = nNumDiffRanks1;
    // INCHI✔️❌:
    // INCHI✔️❌:             pRankStack1[0][neigh1] = nRank2;
    // INCHI✔️❌:             nNumDiffRanks1 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                          nNumDiffRanks, pRankStack1[0], nTempRank,
    // INCHI✔️❌:                                          pRankStack1[1], lNumIter, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:             pRankStack2[0][neigh2] = nRank2;
    // INCHI✔️❌:             nNumDiffRanks2 = DifferentiateRanksBasic( pCG, num_at_tg, NeighList,
    // INCHI✔️❌:                                          nNumDiffRanks, pRankStack2[0], nTempRank,
    // INCHI✔️❌:                                          pRankStack2[1], lNumIter, 1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  for debug only */
    // INCHI✔️❌:     for (i = 0; i < num_at_tg; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pRankStack1[0][(int) pRankStack1[1][i]] != pRankStack2[0][(int) pRankStack2[1][i]])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return -1;  /*  <BRKPT> */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Resort lists of  neighbors */
    // INCHI✔️❌:     SortNeighListsBySymmAndCanonRank( num_atoms, nl1, pRankStack1[0], nCanonRank );
    // INCHI✔️❌:     SortNeighListsBySymmAndCanonRank( num_atoms, nl2, pRankStack2[0], nCanonRank );
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumDiffRanks1 + 1;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: BreakNeighborsTie
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:638-640
    // INCHI✔️❌: #define BREAK_ONE_MORE_SC_TIE       1
    // INCHI✔️❌: #define BREAK_ALSO_NEIGH_TIE        0
    // END INCHI ACTIVE MACRO CONFIGURATION

    let n1 = i32::from(source_get(heap, neigh_num, in1)?);
    let n2 = i32::from(source_get(heap, neigh_num, in2)?);
    let (ia_valence, ia_neighbors) = {
        let atoms = heap.slice(at.as_const())?;
        let atom = atoms
            .get(usize::try_from(ia).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        (i32::from(atom.valence), atom.neighbor)
    };
    let mode_std = MAP_MODE_STD as i32;
    let mode_c2v = MAP_MODE_C2v as i32;
    let mode_c2 = MAP_MODE_C2 as i32;
    let mode_s4 = MAP_MODE_S4 as i32;

    if (mode == mode_s4 && in1 != 0)
        || (mode != mode_std && ia_valence != MAX_NUM_STEREO_ATOM_NEIGH as i32)
        || (mode != mode_std
            && source_get(heap, nSymmRank, n1)? != source_get(heap, nSymmRank, n2)?)
    {
        return Ok(0);
    }

    let count = usize::try_from(num_at_tg).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let rank10 = source_get(heap, pRankStack1, 0)?;
    let order10 = source_get(heap, pRankStack1, 1)?;
    copy_rank_prefix(heap, rank10, nSymmRank, count)?;
    pCG.m_pn_RankForSort = rank10.as_const();
    insertions_sort_AT_NUMBERS(heap, order10, num_at_tg, &mut |heap, left, right| {
        CompRanksOrd(heap, left, right, pCG)
    })?;
    let mut nNumDiffRanks = SortedEquInfoToRanks(heap, rank10, rank10, order10, num_at_tg, None)?;

    let mut other_neigh = [0_i32; 2];
    let mut other_neig_ord = [0_i32; 2];
    let mut num_other_neigh = 0_i32;
    if ia_valence <= MAX_NUM_STEREO_ATOM_NEIGH as i32 && mode != 0 {
        let mut i = 0_i32;
        while i < ia_valence {
            if i != in1 && i != in2 {
                let slot = usize::try_from(num_other_neigh)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                *other_neigh
                    .get_mut(slot)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? =
                    i32::from(source_get(heap, neigh_num, i)?);
                *other_neig_ord
                    .get_mut(slot)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = i;
                num_other_neigh = num_other_neigh.wrapping_add(1);
            }
            i = i.wrapping_add(1);
        }
    }
    let _ = (other_neig_ord, num_other_neigh);
    if (mode != mode_std
        && source_get(heap, nSymmRank, other_neigh[0])?
            != source_get(heap, nSymmRank, other_neigh[1])?)
        || (mode == mode_s4
            && source_get(heap, nSymmRank, n1)? != source_get(heap, nSymmRank, other_neigh[1])?)
    {
        return Ok(0);
    }

    if source_get(heap, rank10, ia)? != source_get(heap, nSymmRank, ia)? {
        let rank = source_get(heap, nSymmRank, ia)?;
        source_set(heap, rank10, ia, rank)?;
        nNumDiffRanks = DifferentiateRanksBasic(
            heap,
            pCG,
            num_at_tg,
            NeighList,
            nNumDiffRanks,
            rank10,
            nTempRank,
            order10,
            lNumIter,
            1,
        )?;
    }
    if ib < num_atoms {
        let old_rank = source_get(heap, rank10, ib)?;
        let new_rank = GetMinNewRank(heap, rank10, order10, old_rank)?;
        if new_rank != old_rank {
            source_set(heap, rank10, ib, new_rank)?;
            nNumDiffRanks = DifferentiateRanksBasic(
                heap,
                pCG,
                num_at_tg,
                NeighList,
                nNumDiffRanks,
                rank10,
                nTempRank,
                order10,
                lNumIter,
                1,
            )?;
        }
    }

    if source_get(heap, rank10, n1)? != source_get(heap, rank10, n2)? {
        return Ok(0);
    }
    let tied_rank = source_get(heap, rank10, n1)?;
    let new_rank = GetMinNewRank(heap, rank10, order10, tied_rank)?;

    let rank20 = source_get(heap, pRankStack2, 0)?;
    let order20 = source_get(heap, pRankStack2, 1)?;
    copy_rank_prefix(heap, rank20, rank10, count)?;
    copy_rank_prefix(heap, order20, order10, count)?;

    source_set(heap, rank10, n1, new_rank)?;
    let mut nNumDiffRanks1 = DifferentiateRanksBasic(
        heap,
        pCG,
        num_at_tg,
        NeighList,
        nNumDiffRanks,
        rank10,
        nTempRank,
        order10,
        lNumIter,
        1,
    )?;
    source_set(heap, rank20, n2, new_rank)?;
    let mut nNumDiffRanks2 = DifferentiateRanksBasic(
        heap,
        pCG,
        num_at_tg,
        NeighList,
        nNumDiffRanks,
        rank20,
        nTempRank,
        order20,
        lNumIter,
        1,
    )?;
    if nNumDiffRanks1 != nNumDiffRanks2 {
        return Ok(-1);
    }

    if mode == mode_c2v || mode == mode_c2 {
        let mut i = 0_i32;
        while i <= 2 {
            let (nn1, nn2) = if i == 0 {
                (n2, n1)
            } else if mode == mode_c2v {
                let index = usize::try_from(i.wrapping_sub(1))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                (other_neigh[index], other_neigh[index])
            } else if mode == mode_c2 {
                let first = usize::try_from(i.wrapping_sub(1))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let second = usize::try_from(2_i32.wrapping_sub(i))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                (other_neigh[first], other_neigh[second])
            } else {
                return Ok(-1);
            };
            let old1 = source_get(heap, rank10, nn1)?;
            let old2 = source_get(heap, rank20, nn2)?;
            let next1 = GetMinNewRank(heap, rank10, order10, old1)?;
            let next2 = GetMinNewRank(heap, rank20, order20, old2)?;
            if old1 == old2 && next1 == next2 {
                if old1 != next1 {
                    source_set(heap, rank10, nn1, next1)?;
                    nNumDiffRanks1 = DifferentiateRanksBasic(
                        heap,
                        pCG,
                        num_at_tg,
                        NeighList,
                        nNumDiffRanks,
                        rank10,
                        nTempRank,
                        order10,
                        lNumIter,
                        1,
                    )?;
                    source_set(heap, rank20, nn2, next2)?;
                    nNumDiffRanks2 = DifferentiateRanksBasic(
                        heap,
                        pCG,
                        num_at_tg,
                        NeighList,
                        nNumDiffRanks,
                        rank20,
                        nTempRank,
                        order20,
                        lNumIter,
                        1,
                    )?;
                    if nNumDiffRanks1 != nNumDiffRanks2 {
                        return Ok(-1);
                    }
                }
            } else {
                return Ok(0);
            }
            i = i.wrapping_add(1);
        }
    }

    if mode == mode_s4 {
        let mut i = 0_i32;
        while i <= 2 {
            let odd = usize::try_from(1_i32.wrapping_sub(in2 % 2))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let even = usize::try_from(in2 % 2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let (nn1, nn2) = match i {
                0 => (n2, other_neigh[odd]),
                1 => (other_neigh[odd], other_neigh[even]),
                2 => (other_neigh[even], n1),
                _ => return Ok(-1),
            };
            let old1 = source_get(heap, rank10, nn1)?;
            let old2 = source_get(heap, rank20, nn2)?;
            let next1 = GetMinNewRank(heap, rank10, order10, old1)?;
            let next2 = GetMinNewRank(heap, rank20, order20, old2)?;
            if old1 == old2 && next1 == next2 {
                if old1 != next1 {
                    source_set(heap, rank10, nn1, next1)?;
                    nNumDiffRanks1 = DifferentiateRanksBasic(
                        heap,
                        pCG,
                        num_at_tg,
                        NeighList,
                        nNumDiffRanks,
                        rank10,
                        nTempRank,
                        order10,
                        lNumIter,
                        1,
                    )?;
                    source_set(heap, rank20, nn2, next2)?;
                    nNumDiffRanks2 = DifferentiateRanksBasic(
                        heap,
                        pCG,
                        num_at_tg,
                        NeighList,
                        nNumDiffRanks,
                        rank20,
                        nTempRank,
                        order20,
                        lNumIter,
                        1,
                    )?;
                    if nNumDiffRanks1 != nNumDiffRanks2 {
                        return Ok(-1);
                    }
                }
            } else {
                return Ok(0);
            }
            i = i.wrapping_add(1);
        }
    }

    if ib >= num_atoms && ia_valence == MAX_NUM_STEREO_ATOM_NEIGH as i32 {
        let tied = source_get(heap, rank10, n2)?;
        let mut num_eq = 0_i32;
        let mut i = 0_i32;
        while i < ia_valence {
            let neighbor = *ia_neighbors
                .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            num_eq = num_eq.wrapping_add(i32::from(
                tied == source_get(heap, rank10, i32::from(neighbor))?,
            ));
            i = i.wrapping_add(1);
        }
        if num_eq == MAX_NUM_STEREO_ATOM_NEIGH as i32 - 1 {
            let next = GetMinNewRank(heap, rank10, order10, tied)?;
            nNumDiffRanks = nNumDiffRanks1;
            source_set(heap, rank10, n2, next)?;
            nNumDiffRanks1 = DifferentiateRanksBasic(
                heap,
                pCG,
                num_at_tg,
                NeighList,
                nNumDiffRanks,
                rank10,
                nTempRank,
                order10,
                lNumIter,
                1,
            )?;
            source_set(heap, rank20, n1, next)?;
            nNumDiffRanks2 = DifferentiateRanksBasic(
                heap,
                pCG,
                num_at_tg,
                NeighList,
                nNumDiffRanks,
                rank20,
                nTempRank,
                order20,
                lNumIter,
                1,
            )?;
        }
    }
    if nNumDiffRanks1 != nNumDiffRanks2 {
        return Ok(-1);
    }

    let mut i = 0_i32;
    while i < num_at_tg {
        let atom1 = source_get(heap, order10, i)?;
        let atom2 = source_get(heap, order20, i)?;
        if source_get(heap, rank10, i32::from(atom1))?
            != source_get(heap, rank20, i32::from(atom2))?
        {
            return Ok(-1);
        }
        i = i.wrapping_add(1);
    }

    SortNeighListsBySymmAndCanonRank(heap, num_atoms, nl1, rank10, nCanonRank)?;
    SortNeighListsBySymmAndCanonRank(heap, num_atoms, nl2, rank20, nCanonRank)?;
    Ok(nNumDiffRanks1.wrapping_add(1))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CheckNextSymmNeighborsAndBonds(
    heap: &SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    cur1: AT_RANK,
    cur2: AT_RANK,
    n1: AT_RANK,
    n2: AT_RANK,
    nAvoidCheckAtom: SourceMutPointer<AT_RANK>,
    nVisited1: SourceMutPointer<AT_RANK>,
    nVisited2: SourceMutPointer<AT_RANK>,
    nVisitOrd1: SourceMutPointer<AT_RANK>,
    nVisitOrd2: SourceMutPointer<AT_RANK>,
    nRank1: SourceMutPointer<AT_RANK>,
    nRank2: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2060 CheckNextSymmNeighborsAndBonds
    // INCHI✔️✔️: int CheckNextSymmNeighborsAndBonds( sp_ATOM *at,
    // INCHI✔️✔️:                                     AT_RANK cur1,
    // INCHI✔️✔️:                                     AT_RANK cur2,
    // INCHI✔️✔️:                                     AT_RANK n1,
    // INCHI✔️✔️:                                     AT_RANK n2,
    // INCHI✔️✔️:                                     AT_RANK *nAvoidCheckAtom,
    // INCHI✔️✔️:                                     AT_RANK *nVisited1,
    // INCHI✔️✔️:                                     AT_RANK *nVisited2,
    // INCHI✔️✔️:                                     AT_RANK *nVisitOrd1,
    // INCHI✔️✔️:                                     AT_RANK *nVisitOrd2,
    // INCHI✔️✔️:                                     const AT_RANK *nRank1,
    // INCHI✔️✔️:                                     const AT_RANK *nRank2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_RANK s1 = 0, s2 = 0;
    // INCHI✔️✔️:     int     i1, i2, k1, k2;
    // INCHI✔️✔️:     if (nRank1[n1] != nRank2[n2])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /*  parallel traversal in stereo removal failed */ /*   <BRKPT> */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     switch (!nVisited1[n1] + !nVisited2[n2])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         case 0:
    // INCHI✔️✔️:             if (nVisited1[n1] != n2 + 1 || nVisited2[n2] != n1 + 1)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return -1; /*  0; */ /*  possibly error???: we have come to an alreardy traversed pair and */
    // INCHI✔️✔️:                            /*  found that the pair previously has not been traversed synchroneously. */
    // INCHI✔️✔️:             }              /*  -- Happens in C60. */
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         case 1:
    // INCHI✔️✔️:             return -1; /*  0; */ /*  possibly error: one is zero, another is not a zero. Happens in C60 */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /*  case 2: */
    // INCHI✔️✔️:             /* both are zero, OK. */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (nVisitOrd1[n1] != nVisitOrd2[n2])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /*  0; */ /*  different DFS trees */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     /*  at[n1] and at[n2] are next to at[cur1] and at[cur2] respectively */
    // INCHI✔️✔️:     /*  Even though the bond might have already been checked, check whether */
    // INCHI✔️✔️:     /*  it is a stereo bond/cumulene. If it is, check the bond/cumulene parity. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  Even though the bond or cumulene might have already been checked, check it: this is */
    // INCHI✔️✔️:     /*  the only place we can check stereo bonds and cumulenes that are not edges of the DFS tree */
    // INCHI✔️✔️:     /*  The code works both for a stereo bond and a stereogenic cumulene. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i1 = 0, k1 = 0; i1 < MAX_NUM_STEREO_BONDS &&
    // INCHI✔️✔️:         ( s1 = at[cur1].stereo_bond_neighbor[i1] ) &&
    // INCHI✔️✔️:                          !( k1 = ( at[cur1].neighbor[(int) at[cur1].stereo_bond_ord[i1]] == n1 ) ); i1++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i2 = 0, k2 = 0; i2 < MAX_NUM_STEREO_BONDS &&
    // INCHI✔️✔️:         ( s2 = at[cur2].stereo_bond_neighbor[i2] ) &&
    // INCHI✔️✔️:                          !( k2 = ( at[cur2].neighbor[(int) at[cur2].stereo_bond_ord[i2]] == n2 ) ); i2++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* -- this does not work in case of cumulenes --
    // INCHI✔️✔️:     for ( i1 = 0, k1 = 0; i1 < MAX_NUM_STEREO_BONDS && (s1=at[cur1].stereo_bond_neighbor[i1]) && !(k1=(s1-1 == n1)); i1 ++ )
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     for ( i2 = 0, k2 = 0; i2 < MAX_NUM_STEREO_BONDS && (s2=at[cur2].stereo_bond_neighbor[i2]) && !(k2=(s2-1 == n2)); i2 ++ )
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (k1 != k2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /*  not an error: a stereo bond and not a stereo bond */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (k1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* here k1 == k2 */
    // INCHI✔️✔️:         int bCheckBond1, bCheckBond2;
    // INCHI✔️✔️:         s1--;
    // INCHI✔️✔️:         s2--;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         bCheckBond1 = ( cur1 != nAvoidCheckAtom[0] || s1 != nAvoidCheckAtom[1] ) &&
    // INCHI✔️✔️:             ( cur1 != nAvoidCheckAtom[1] || s1 != nAvoidCheckAtom[0] );
    // INCHI✔️✔️:         bCheckBond2 = ( cur2 != nAvoidCheckAtom[0] || s2 != nAvoidCheckAtom[1] ) &&
    // INCHI✔️✔️:             ( cur2 != nAvoidCheckAtom[1] || s2 != nAvoidCheckAtom[0] );
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (bCheckBond1 != bCheckBond2)
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!bCheckBond1 && !bCheckBond2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 1; /*  do not go any further in this direction */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (at[cur1].stereo_bond_parity[i1] != at[cur2].stereo_bond_parity[i2])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /*  different values of  at[].stereo_bond_parity: definitely different bonds */
    // INCHI✔️✔️:             /*  known parities */
    // INCHI✔️✔️:             if (PARITY_KNOWN( at[cur1].stereo_bond_parity[i1] ) &&
    // INCHI✔️✔️:                  PARITY_KNOWN( at[cur2].stereo_bond_parity[i2] ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return 0; /*  different currently known stereo bond parities */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️: #if ( PROPAGATE_ILL_DEF_STEREO != 1 )
    // INCHI✔️✔️:             /*  well defined and to be calculated from the ranks */
    // INCHI✔️✔️:             if (!( PARITY_CALCULATE( at[cur1].stereo_bond_parity[i1] ) && PARITY_WELL_DEF( at[cur2].stereo_bond_parity[i2] ) ||
    // INCHI✔️✔️:                 PARITY_WELL_DEF( at[cur1].stereo_bond_parity[i1] ) && PARITY_CALCULATE( at[cur2].stereo_bond_parity[i2] ) ||
    // INCHI✔️✔️:                 PARITY_CALCULATE( at[cur1].stereo_bond_parity[i1] ) && PARITY_CALCULATE( at[cur2].stereo_bond_parity[i2] ) ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️: /*  do not reject if: "well defined" and "calculate" or "calculate" and "calculate" */
    // INCHI✔️✔️:                 return 0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️: #if ( PROPAGATE_ILL_DEF_STEREO != 1 )
    // INCHI✔️✔️:         if (( cur1 != cur2 || s1 != s2 ) && ( cur1 != s2 || cur2 != s1 ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️: /*  two different stereo bonds */
    // INCHI✔️✔️:             if (PARITY_ILL_DEF( at[cur1].stereo_bond_parity[i1] ) ||
    // INCHI✔️✔️:                  PARITY_ILL_DEF( at[cur2].stereo_bond_parity[i2] ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return 0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 1; /*  stereo bonds to n1 and n2 have same known parities or are not stereo bonds */
    // INCHI✔️✔️: }
    // INCHI✔️✔️:
    // END INCHI C FUNCTION: CheckNextSymmNeighborsAndBonds
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:629
    // INCHI✔️✔️: #define PROPAGATE_ILL_DEF_STEREO    1
    // END INCHI ACTIVE MACRO CONFIGURATION

    if source_get(heap, nRank1, i32::from(n1))? != source_get(heap, nRank2, i32::from(n2))? {
        return Ok(-1);
    }
    let visited1 = source_get(heap, nVisited1, i32::from(n1))?;
    let visited2 = source_get(heap, nVisited2, i32::from(n2))?;
    match i32::from(visited1 == 0).wrapping_add(i32::from(visited2 == 0)) {
        0 => {
            if visited1 != n2.wrapping_add(1) || visited2 != n1.wrapping_add(1) {
                return Ok(-1);
            }
        }
        1 => return Ok(-1),
        _ => {}
    }
    if source_get(heap, nVisitOrd1, i32::from(n1))? != source_get(heap, nVisitOrd2, i32::from(n2))?
    {
        return Ok(-1);
    }

    let atoms = heap.slice(at.as_const())?;
    let atom1 = atoms
        .get(usize::from(cur1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atom2 = atoms
        .get(usize::from(cur2))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut i1 = 0_usize;
    let mut k1 = false;
    let mut s1 = 0_u16;
    while i1 < MAX_NUM_STEREO_BONDS as usize {
        s1 = atom1.stereo_bond_neighbor[i1];
        if s1 == 0 {
            break;
        }
        let order = usize::try_from(atom1.stereo_bond_ord[i1])
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        k1 = *atom1
            .neighbor
            .get(order)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            == n1;
        if k1 {
            break;
        }
        i1 += 1;
    }
    let mut i2 = 0_usize;
    let mut k2 = false;
    let mut s2 = 0_u16;
    while i2 < MAX_NUM_STEREO_BONDS as usize {
        s2 = atom2.stereo_bond_neighbor[i2];
        if s2 == 0 {
            break;
        }
        let order = usize::try_from(atom2.stereo_bond_ord[i2])
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        k2 = *atom2
            .neighbor
            .get(order)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            == n2;
        if k2 {
            break;
        }
        i2 += 1;
    }

    if k1 != k2 {
        return Ok(0);
    }
    if k1 {
        s1 = s1.wrapping_sub(1);
        s2 = s2.wrapping_sub(1);
        let avoid0 = source_get(heap, nAvoidCheckAtom, 0)?;
        let avoid1 = source_get(heap, nAvoidCheckAtom, 1)?;
        let check1 = (cur1 != avoid0 || s1 != avoid1) && (cur1 != avoid1 || s1 != avoid0);
        let check2 = (cur2 != avoid0 || s2 != avoid1) && (cur2 != avoid1 || s2 != avoid0);
        if check1 != check2 {
            return Ok(0);
        }
        if !check1 && !check2 {
            return Ok(1);
        }
        let parity1 = atom1.stereo_bond_parity[i1];
        let parity2 = atom2.stereo_bond_parity[i2];
        if parity1 != parity2 {
            let value1 = i32::from(parity1) & BITS_PARITY as i32;
            let value2 = i32::from(parity2) & BITS_PARITY as i32;
            let known1 =
                (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&value1);
            let known2 =
                (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&value2);
            if known1 && known2 {
                return Ok(0);
            }
        }
    }
    Ok(1)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CreateCheckSymmPaths(
    heap: &mut SourceHeap,
    at: SourceMutPointer<sp_ATOM>,
    prev1: AT_RANK,
    cur1: AT_RANK,
    prev2: AT_RANK,
    cur2: AT_RANK,
    nAvoidCheckAtom: SourceMutPointer<AT_RANK>,
    nVisited1: SourceMutPointer<AT_RANK>,
    nVisited2: SourceMutPointer<AT_RANK>,
    nVisitOrd1: SourceMutPointer<AT_RANK>,
    nVisitOrd2: SourceMutPointer<AT_RANK>,
    nl1: SourceMutPointer<NEIGH_LIST>,
    nl2: SourceMutPointer<NEIGH_LIST>,
    nRank1: SourceMutPointer<AT_RANK>,
    nRank2: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nLength: &mut AT_RANK,
    bParitiesInverted: &mut i32,
    mode: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2191 CreateCheckSymmPaths
    // INCHI✔️❌: int CreateCheckSymmPaths( sp_ATOM *at,
    // INCHI✔️❌:                           AT_RANK prev1,
    // INCHI✔️❌:                           AT_RANK cur1,
    // INCHI✔️❌:                           AT_RANK prev2,
    // INCHI✔️❌:                           AT_RANK cur2,
    // INCHI✔️❌:                           AT_RANK *nAvoidCheckAtom,
    // INCHI✔️❌:                           AT_RANK *nVisited1,
    // INCHI✔️❌:                           AT_RANK *nVisited2,
    // INCHI✔️❌:                           AT_RANK *nVisitOrd1,
    // INCHI✔️❌:                           AT_RANK *nVisitOrd2,
    // INCHI✔️❌:                           NEIGH_LIST *nl1,
    // INCHI✔️❌:                           NEIGH_LIST *nl2,
    // INCHI✔️❌:                           const AT_RANK *nRank1,
    // INCHI✔️❌:                           const AT_RANK *nRank2,
    // INCHI✔️❌:                           AT_RANK *nCanonRank,
    // INCHI✔️❌:                           AT_RANK *nLength,
    // INCHI✔️❌:                           int *bParitiesInverted,
    // INCHI✔️❌:                           int mode )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int k, k1, k2, ret = 0, bParitiesInvertedZero = 0, *pbParitiesInverted;
    // INCHI✔️❌:     AT_RANK n1, n2;
    // INCHI✔️❌:
    // INCHI✔️❌:     nVisited1[cur1] = cur2 + 1;  /*  symmetrically exchange atom numbers */
    // INCHI✔️❌:     nVisited2[cur2] = cur1 + 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     ( *nLength )++;
    // INCHI✔️❌:
    // INCHI✔️❌:     nVisitOrd1[cur1] = *nLength; /*  save DFS visit order */
    // INCHI✔️❌:     nVisitOrd2[cur2] = *nLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* new version allows all inverted parities */
    // INCHI✔️❌:     if (PARITY_WELL_DEF( at[cur1].stereo_atom_parity ) &&
    // INCHI✔️❌:          PARITY_WELL_DEF( at[cur2].stereo_atom_parity ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (*bParitiesInverted < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *bParitiesInverted = ( at[cur1].stereo_atom_parity + at[cur2].stereo_atom_parity ) % 2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (*bParitiesInverted != ( at[cur1].stereo_atom_parity + at[cur2].stereo_atom_parity ) % 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /*  Different known in advance parities have wrong "inverted" relation */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (PARITY_KNOWN( at[cur1].stereo_atom_parity ) &&
    // INCHI✔️❌:             PARITY_KNOWN( at[cur2].stereo_atom_parity ) &&
    // INCHI✔️❌:             at[cur1].stereo_atom_parity != at[cur2].stereo_atom_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;  /*  Different known in advance parities */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (cur1 != cur2 &&
    // INCHI✔️❌:          !at[cur1].stereo_bond_neighbor[0] && !at[cur2].stereo_bond_neighbor[0] &&
    // INCHI✔️❌:          PARITY_KNOWN( at[cur1].parity ) != PARITY_KNOWN( at[cur2].parity ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /*  one atom is stereogenic, another (presumably equivalent) is not. 9-11-2002 */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #if ( PROPAGATE_ILL_DEF_STEREO != 1 )
    // INCHI✔️❌:     if (cur1 != cur2 &&
    // INCHI✔️❌:         ( PARITY_ILL_DEF( at[cur1].stereo_atom_parity ) ||
    // INCHI✔️❌:             PARITY_ILL_DEF( at[cur2].stereo_atom_parity ) )
    // INCHI✔️❌:        )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;  /*  Cannot detect whether the paths are same or different */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[cur1].valence != at[cur2].valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return CT_REMOVE_STEREO_ERR; /*  program error */ /*   <BRKPT> */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (at[cur1].valence == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /*  so far success */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nl1[(int) cur1][0] != nl2[(int) cur2][0] || nl1[(int) cur1][0] != at[cur1].valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return CT_REMOVE_STEREO_ERR; /*  error: different valences */ /*   <BRKPT> */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k = 1, k1 = 1, k2 = 1; k < at[cur1].valence; k++, k1++, k2++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( n1 = nl1[(int) cur1][k1] ) == prev1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             n1 = nl1[(int) cur1][++k1]; /*  don't go back */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (( n2 = nl2[(int) cur2][k2] ) == prev2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             n2 = nl2[(int) cur2][++k2]; /*  don't go back */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (0 >= ( ret = CheckNextSymmNeighborsAndBonds( at, cur1, cur2, n1, n2, nAvoidCheckAtom,
    // INCHI✔️❌:             nVisited1, nVisited2, nVisitOrd1, nVisitOrd2, nRank1, nRank2 ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret; /*  different neighbors or bonds                       */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!nVisited1[n1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  recursion */
    // INCHI✔️❌:             /* allow all inverted parities only inside a single ring system containing the starting point */
    // INCHI✔️❌:             pbParitiesInverted = ( at[cur1].nRingSystem == at[n1].nRingSystem ) ? bParitiesInverted : &bParitiesInvertedZero;
    // INCHI✔️❌:             if (0 >= ( ret = CreateCheckSymmPaths( at, cur1, n1, cur2, n2, nAvoidCheckAtom,
    // INCHI✔️❌:                 nVisited1, nVisited2, nVisitOrd1, nVisitOrd2,
    // INCHI✔️❌:                 nl1, nl2, nRank1, nRank2, nCanonRank, nLength, pbParitiesInverted, mode ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 1; /*  Success */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateCheckSymmPaths
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:629
    // INCHI✔️❌: #define PROPAGATE_ILL_DEF_STEREO    1
    // END INCHI ACTIVE MACRO CONFIGURATION

    source_set(heap, nVisited1, i32::from(cur1), cur2.wrapping_add(1))?;
    source_set(heap, nVisited2, i32::from(cur2), cur1.wrapping_add(1))?;
    *nLength = nLength.wrapping_add(1);
    source_set(heap, nVisitOrd1, i32::from(cur1), *nLength)?;
    source_set(heap, nVisitOrd2, i32::from(cur2), *nLength)?;

    let (
        parity1,
        parity2,
        atom_parity1,
        atom_parity2,
        stereo_neighbor1,
        stereo_neighbor2,
        valence1,
        valence2,
    ) = {
        let atoms = heap.slice(at.as_const())?;
        let atom1 = atoms
            .get(usize::from(cur1))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom2 = atoms
            .get(usize::from(cur2))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        (
            atom1.stereo_atom_parity,
            atom2.stereo_atom_parity,
            atom1.parity,
            atom2.parity,
            atom1.stereo_bond_neighbor[0],
            atom2.stereo_bond_neighbor[0],
            i32::from(atom1.valence),
            i32::from(atom2.valence),
        )
    };
    let parity_value = |value: i8| i32::from(value) & BITS_PARITY as i32;
    let well_defined = |value: i8| {
        (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
            .contains(&parity_value(value))
    };
    let known = |value: i8| {
        (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&parity_value(value))
    };

    if well_defined(parity1) && well_defined(parity2) {
        let inverted = i32::from(parity1)
            .wrapping_add(i32::from(parity2))
            .wrapping_rem(2);
        if *bParitiesInverted < 0 {
            *bParitiesInverted = inverted;
        } else if *bParitiesInverted != inverted {
            return Ok(0);
        }
    } else if known(parity1) && known(parity2) && parity1 != parity2 {
        return Ok(0);
    }

    if cur1 != cur2
        && stereo_neighbor1 == 0
        && stereo_neighbor2 == 0
        && known(atom_parity1) != known(atom_parity2)
    {
        return Ok(0);
    }
    if valence1 != valence2 {
        return Ok(CT_REMOVE_STEREO_ERR);
    }
    if valence1 == 1 {
        return Ok(1);
    }

    let list1 = source_get(heap, nl1, i32::from(cur1))?;
    let list2 = source_get(heap, nl2, i32::from(cur2))?;
    let list_length1 = i32::from(source_get(heap, list1, 0)?);
    let list_length2 = i32::from(source_get(heap, list2, 0)?);
    if list_length1 != list_length2 || list_length1 != valence1 {
        return Ok(CT_REMOVE_STEREO_ERR);
    }

    let mut bParitiesInvertedZero = 0_i32;
    let mut k = 1_i32;
    let mut k1 = 1_i32;
    let mut k2 = 1_i32;
    while k < valence1 {
        let mut n1 = source_get(heap, list1, k1)?;
        if n1 == prev1 {
            k1 = k1.wrapping_add(1);
            n1 = source_get(heap, list1, k1)?;
        }
        let mut n2 = source_get(heap, list2, k2)?;
        if n2 == prev2 {
            k2 = k2.wrapping_add(1);
            n2 = source_get(heap, list2, k2)?;
        }
        let ret = CheckNextSymmNeighborsAndBonds(
            heap,
            at,
            cur1,
            cur2,
            n1,
            n2,
            nAvoidCheckAtom,
            nVisited1,
            nVisited2,
            nVisitOrd1,
            nVisitOrd2,
            nRank1,
            nRank2,
        )?;
        if ret <= 0 {
            return Ok(ret);
        }

        if source_get(heap, nVisited1, i32::from(n1))? == 0 {
            let same_ring = {
                let atoms = heap.slice(at.as_const())?;
                atoms
                    .get(usize::from(cur1))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nRingSystem
                    == atoms
                        .get(usize::from(n1))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nRingSystem
            };
            let parities = if same_ring {
                &mut *bParitiesInverted
            } else {
                &mut bParitiesInvertedZero
            };
            let ret = CreateCheckSymmPaths(
                heap,
                at,
                cur1,
                n1,
                cur2,
                n2,
                nAvoidCheckAtom,
                nVisited1,
                nVisited2,
                nVisitOrd1,
                nVisitOrd2,
                nl1,
                nl2,
                nRank1,
                nRank2,
                nCanonRank,
                nLength,
                parities,
                mode,
            )?;
            if ret <= 0 {
                return Ok(ret);
            }
        }
        k = k.wrapping_add(1);
        k1 = k1.wrapping_add(1);
        k2 = k2.wrapping_add(1);
    }
    let _ = (nCanonRank, mode);
    Ok(1)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CalculatedPathsParitiesAreIdentical(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon1: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon2: SourceMutPointer<AT_RANK>,
    nVisited1: SourceMutPointer<AT_RANK>,
    nVisited2: SourceMutPointer<AT_RANK>,
    prev_sb_neigh: AT_RANK,
    cur: AT_RANK,
    next1: AT_RANK,
    next2: AT_RANK,
    nNeighMode: i32,
    mut bParitiesInverted: i32,
    mode: i32,
    pCS: &mut CANON_STAT,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2399 CalculatedPathsParitiesAreIdentical
    // INCHI✔️❌: int CalculatedPathsParitiesAreIdentical( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                          sp_ATOM *at, int num_atoms,
    // INCHI✔️❌:                                          const AT_RANK *nSymmRank,
    // INCHI✔️❌:                                          AT_RANK *nCanonRank,
    // INCHI✔️❌:                                          AT_RANK *nAtomNumberCanon,
    // INCHI✔️❌:                                          AT_RANK *nAtomNumberCanon1,
    // INCHI✔️❌:                                          AT_RANK *nAtomNumberCanon2,
    // INCHI✔️❌:                                          AT_RANK *nVisited1,
    // INCHI✔️❌:                                          AT_RANK *nVisited2,
    // INCHI✔️❌:                                          AT_RANK prev_sb_neigh,
    // INCHI✔️❌:                                          AT_RANK cur,
    // INCHI✔️❌:                                          AT_RANK next1,
    // INCHI✔️❌:                                          AT_RANK next2,
    // INCHI✔️❌:                                          int nNeighMode,
    // INCHI✔️❌:                                          int bParitiesInverted,
    // INCHI✔️❌:                                          int mode,
    // INCHI✔️❌:                                          CANON_STAT *pCS,
    // INCHI✔️❌:                                          int vABParityUnknown )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, i01, i02, i11, i12, i21, i22, k, parity, parity1, parity2, parity12, num_other_neigh;
    // INCHI✔️❌:     int nNumEqStereogenic, nCheckingMode, not_well_def_parities;
    // INCHI✔️❌:     AT_RANK other_neigh[MAX_NUM_STEREO_ATOM_NEIGH], neigh, r1, r2;
    // INCHI✔️❌:     int  nNumComparedCenters = 0, nNumComparedBonds = 0, bCurParityInv1 = 0 /*, bCurParityInv2=0*/;
    // INCHI✔️❌:     int  bCurRotated = 0, nNumDiff = 0, nNumInv = 0;
    // INCHI✔️❌:     int  s1, s2;
    // INCHI✔️❌:
    // INCHI✔️❌:     nCheckingMode = ( prev_sb_neigh < num_atoms ) ? CHECKING_STEREOBOND : CHECKING_STEREOCENTER;
    // INCHI✔️❌:     not_well_def_parities = 0;
    // INCHI✔️❌:     nNumEqStereogenic = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((nNeighMode != NEIGH_MODE_RING &&
    // INCHI✔️❌:          bParitiesInverted != 0) || abs( bParitiesInverted ) != 1) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bParitiesInverted = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bParitiesInverted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0, i11 = i22 = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* count number of atoms that have not been visited */
    // INCHI✔️❌:             i11 += !nVisited1[i];
    // INCHI✔️❌:             i22 += !nVisited2[i];
    // INCHI✔️❌:             nAtomNumberCanon1[i] = MAX_ATOMS + 1;  /*  mark unchanged */
    // INCHI✔️❌:             nAtomNumberCanon2[i] = MAX_ATOMS + 1;  /*  mark unchanged */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i11 || i22)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (bParitiesInverted == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /* only a part of the structure has been inverted */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bParitiesInverted = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nAtomNumberCanon1[i] = MAX_ATOMS + 1;  /*  mark unchanged */
    // INCHI✔️❌:             nAtomNumberCanon2[i] = MAX_ATOMS + 1;  /*  mark unchanged */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((bParitiesInverted > 0 && !( mode == MAP_MODE_C2v || mode == MAP_MODE_S4 )) ||
    // INCHI✔️❌:          (bParitiesInverted == 0 && !( mode == MAP_MODE_C2 || mode == MAP_MODE_STD ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /**************************************************************************************
    // INCHI✔️❌:      *    The following discussion assumes that the canonical numbers are
    // INCHI✔️❌:      *    switched for some pairs of constitutionally identical atoms
    // INCHI✔️❌:      *    in such a way that the new numbering is an equivalent to the
    // INCHI✔️❌:      *    nCanonRank[] canonical numbering (the transposition belongs to the
    // INCHI✔️❌:      *    automorphism group of the chemical structure having no stereo).
    // INCHI✔️❌:      *    At this point non-zero elements of nVisited1[] and nVisited2[]
    // INCHI✔️❌:      *    together contain transposition P of the atom numbers.
    // INCHI✔️❌:      *    and P2 respectively of the ordering atom numbers: nVisitedi[k] = Pi(k)+1;
    // INCHI✔️❌:      *    In this implementation:
    // INCHI✔️❌:      *       P1(k)=k for all k
    // INCHI✔️❌:      *       P2(cur)=cur, P2(next1)=next2, P2(next2)=next1
    // INCHI✔️❌:      *
    // INCHI✔️❌:      *    Below we call one of the numberings "old", another "new".
    // INCHI✔️❌:      *
    // INCHI✔️❌:      *    *IF* the old and the new canonical numberings produce same parities for stereogenic
    // INCHI✔️❌:      *    elements for the same canonical number(s)
    // INCHI✔️❌:      *    (that is, old_parity(canon_number) == new_parity(canon_number)
    // INCHI✔️❌:      *    *except* the currently being tested stereocenter at[cur] or stereobond/cumulene
    // INCHI✔️❌:      *    at[cur]=at[prev_sb_neigh], whose parity MUST be inverted
    // INCHI✔️❌:      *
    // INCHI✔️❌:      *    *THEN* the stereocenter or stereobond/cumulene is not stereogenic with one
    // INCHI✔️❌:      *
    // INCHI✔️❌:      *    *EXCEPTION* If the currently tested stereogenic element is constitutionally
    // INCHI✔️❌:      *    equivalent to two or more other stereogenic elements that have been
    // INCHI✔️❌:      *    permuted then the currently tested one is still stereogenic.
    // INCHI✔️❌:      **************************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:      /*
    // INCHI✔️❌:      * 1. replace the assigned in each of the parallel traversals atom numbers
    // INCHI✔️❌:      *    with the canon. ranks corresponding to the atom numbers in the
    // INCHI✔️❌:      *    currently numbered atoms at[].
    // INCHI✔️❌:      *    One of obtained this way canonical numberings (probably nVisited1[])
    // INCHI✔️❌:      *    is same as the nCanonRank[] because usually nVisited1[i] = i+1 or 0
    // INCHI✔️❌:      */
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (nVisited1[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* canonical number of the atom mapped on atom #i in 'left' path */
    // INCHI✔️❌:             nVisited1[i] = nCanonRank[(int) nVisited1[i] - 1];
    // INCHI✔️❌:             /* reverse: atom # from the mapped canonical rank in 'left' path */
    // INCHI✔️❌:             nAtomNumberCanon1[nVisited1[i] - 1] = i;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nVisited2[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* canonical number of the atom mapped on atom #i in 'right' path */
    // INCHI✔️❌:             nVisited2[i] = nCanonRank[(int) nVisited2[i] - 1];
    // INCHI✔️❌:             /* reverse: atom # from the mapped canonical rank in 'right' path */
    // INCHI✔️❌:             nAtomNumberCanon2[nVisited2[i] - 1] = i;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* if 'left' and 'right' path do not have atoms in common except the
    // INCHI✔️❌:            starting atom (and in case of stereobond, the end atom) some of
    // INCHI✔️❌:            nVisitedi[i] elements may be zero.
    // INCHI✔️❌:         */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:      * if started with a stereobond then check whether its parity has changed.
    // INCHI✔️❌:      * If yes then continue, otherwise parities are different
    // INCHI✔️❌:      *
    // INCHI✔️❌:      * if started with a stereo center then prev_sb_neigh = MAX_ATOMS+1
    // INCHI✔️❌:      *
    // INCHI✔️❌:      * If the transposition of next1 and next2 changes only the parity of the starting stereo atom or stereo bond
    // INCHI✔️❌:      * then the stereo bond or stereo atom is not stereogenic
    // INCHI✔️❌:      *
    // INCHI✔️❌:      * The exception: the stereogenic elememt in question is equivalent
    // INCHI✔️❌:      *    to two or more traversed other stereogenic elememts
    // INCHI✔️❌:      *    (see nNumEqStereogenic below, case similar to trimethylcyclopropane:
    // INCHI✔️❌:      *     3 or more constitutionally equivalent stereogenic elements)
    // INCHI✔️❌:      */
    // INCHI✔️❌:     if (nCheckingMode == CHECKING_STEREOBOND)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /******************************************************************************
    // INCHI✔️❌:          *
    // INCHI✔️❌:          *  Possibly stereogenic starting bond or cumulene at[cur]-at[prev_sb_neigh]
    // INCHI✔️❌:          *
    // INCHI✔️❌:          *******************************************************************************/
    // INCHI✔️❌:         /*  checking the starting stereo bond */
    // INCHI✔️❌:         if (nVisited1[prev_sb_neigh] || nVisited2[prev_sb_neigh])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  the bond or cumulene is in the ring and the opposite atom has been visited */
    // INCHI✔️❌:             if (nVisited1[prev_sb_neigh] != nVisited2[prev_sb_neigh] ||
    // INCHI✔️❌:                  nCanonRank[prev_sb_neigh] != nVisited2[prev_sb_neigh])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /*  error: we came back to the same bond/cumulene and */ /*   <BRKPT> */
    // INCHI✔️❌:                           /*  assigned different canon. ranks to the opposite atom. */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (at[prev_sb_neigh].valence + at[prev_sb_neigh].num_H > 3)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;   /*  at[prev_sb_neigh] atom can not be adjacent to a stereo bond/cumulene */
    // INCHI✔️❌:                             /*  or does not have 3 attachments (hydrogens are not considered here) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (i = 0, k = 0; i < MAX_NUM_STEREO_BONDS &&
    // INCHI✔️❌:                 ( neigh = at[prev_sb_neigh].stereo_bond_neighbor[i] ) && !( k = ( neigh - 1 == cur ) );
    // INCHI✔️❌:                 i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!k)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  program error: could not locate stereogenic bond mark on the opposite atom */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             k = (int) at[prev_sb_neigh].stereo_bond_ord[i]; /*  seq. number of the double or cumulene bond on at[prev_sb_neigh] */
    // INCHI✔️❌:
    // INCHI✔️❌:             for (i = 0, num_other_neigh = 0; i < at[prev_sb_neigh].valence && num_other_neigh <= MAX_OTHER_NEIGH; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i != k)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  do not include the double or cumulene bond */
    // INCHI✔️❌:                     other_neigh[num_other_neigh++] = at[prev_sb_neigh].neighbor[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num_other_neigh + at[prev_sb_neigh].num_H > MAX_OTHER_NEIGH)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (i = 0; i < num_other_neigh; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k = (int) other_neigh[i];
    // INCHI✔️❌:                 if (nVisited1[k] && nVisited1[k] != nCanonRank[k])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0; /*  parity of the statring stereo bond/cumulene has not changed. */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (nVisited2[k] && nVisited2[k] != nCanonRank[k])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0; /*  parity of the statring stereo bond/cumulene has not changed. */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nCheckingMode == CHECKING_STEREOCENTER)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /**************************************************
    // INCHI✔️❌:          *
    // INCHI✔️❌:          *  Possibly stereogenic starting atom at[cur]
    // INCHI✔️❌:          *
    // INCHI✔️❌:          **************************************************/
    // INCHI✔️❌:         /*  checking the starting stereo center */
    // INCHI✔️❌:         for (i = 0, num_other_neigh = 0; i < at[cur].valence && num_other_neigh <= MAX_OTHER_NEIGH; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[cur].neighbor[i];
    // INCHI✔️❌:             if (neigh != next1 && neigh != next2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 other_neigh[num_other_neigh++] = neigh;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (num_other_neigh + at[cur].num_H > MAX_OTHER_NEIGH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         if ( bParitiesInverted && at[cur].valence == MAX_NUM_STEREO_ATOM_NEIGH ) {
    // INCHI✔️❌:             if ( nVisited1[other_neigh[0]] == nCanonRank[other_neigh[0]] ||
    // INCHI✔️❌:                  nVisited2[other_neigh[0]] == nCanonRank[other_neigh[0]] ||
    // INCHI✔️❌:                  nVisited1[other_neigh[1]] == nCanonRank[other_neigh[1]] ||
    // INCHI✔️❌:                  nVisited2[other_neigh[1]] == nCanonRank[other_neigh[1]] ) {
    // INCHI✔️❌:                 bParitiesInverted = 0;
    // INCHI✔️❌:                 bCurRotated = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         */
    // INCHI✔️❌:         /* bParitiesInverted = -1 means no predefined stereocenter has been checked */
    // INCHI✔️❌:         if (bParitiesInverted && at[cur].valence == MAX_NUM_STEREO_ATOM_NEIGH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* special case: 4 canonically eq. neighbors */
    // INCHI✔️❌:             int canon_parity, parity_vis_1, parity_vis_2;
    // INCHI✔️❌:             canon_parity = GetPermutationParity( pCG, at + cur, MAX_ATOMS + 1, nCanonRank );
    // INCHI✔️❌:             parity_vis_1 = GetPermutationParity( pCG, at + cur, MAX_ATOMS + 1, nVisited1 );
    // INCHI✔️❌:             parity_vis_2 = GetPermutationParity( pCG, at + cur, MAX_ATOMS + 1, nVisited2 );
    // INCHI✔️❌:             if (parity_vis_1 != parity_vis_2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (bParitiesInverted == 1 && parity_vis_1 == canon_parity)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /* not a typical case of inversion during the mapping of D4h stereocenter */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bParitiesInverted == -1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (parity_vis_1 == canon_parity)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bParitiesInverted = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bParitiesInverted = 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* at this point bParitiesInverted >= 0 */
    // INCHI✔️❌:         if (!bParitiesInverted && !bCurRotated)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_other_neigh; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k = (int) other_neigh[i];
    // INCHI✔️❌:                 if (nVisited1[k] && nVisited1[k] != nCanonRank[k])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0; /*  parity of the statring stereo center has not changed. */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (nVisited2[k] && nVisited2[k] != nCanonRank[k])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0; /*  parity of the statring stereo center has not changed. */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*****************************************************
    // INCHI✔️❌:      * Check other (non-starting) stereo centers
    // INCHI✔️❌:      ******************************************************/
    // INCHI✔️❌:     for (i = 0; i < pCS->nLenLinearCTStereoCarb; i++, nNumComparedCenters += ( k > 0 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         r1 = pCS->LinearCTStereoCarb[i].at_num;
    // INCHI✔️❌:         i01 = nAtomNumberCanon[r1 - 1]; /*  ord. number of the atom that has canon rank r1 */
    // INCHI✔️❌:
    // INCHI✔️❌:         i11 = nAtomNumberCanon1[r1 - 1]; /*  = (MAX_ATOMS+1) > num_atoms if the atom has not been traversed */
    // INCHI✔️❌:         i12 = nAtomNumberCanon2[r1 - 1]; /*  = otherwise < num_atoms */
    // INCHI✔️❌:
    // INCHI✔️❌:         s1 = ( i11 < num_atoms ); /*  1 => the center was traversed on path #1 */
    // INCHI✔️❌:         s2 = ( i12 < num_atoms ); /*  1 => the center was traversed on path #2 */
    // INCHI✔️❌:
    // INCHI✔️❌:         bCurParityInv1 = ( bParitiesInverted &&
    // INCHI✔️❌:                           at[cur].nRingSystem == at[i11].nRingSystem &&
    // INCHI✔️❌:                           at[cur].nRingSystem == at[i12].nRingSystem );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         k = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  check whether the two stereo centers (they can be one and the same atom) have been traversed */
    // INCHI✔️❌:         if (!s1 && !s2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;  /*  Both stereo centers have not been traversed; check the next pair. */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nCheckingMode == CHECKING_STEREOCENTER)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  check whether the stereocenters are the starting stereocenter */
    // INCHI✔️❌:             switch (( cur == i11 ) + ( cur == i12 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case 2:
    // INCHI✔️❌:                     continue; /*  do not recheck the starting atom */
    // INCHI✔️❌:                 case 1:
    // INCHI✔️❌:                     return -1; /*  possibly program error */ /*   <BRKPT> */
    // INCHI✔️❌:                 /* case 0: */
    // INCHI✔️❌:                 /*     break;  */  /*  the stereo centers are not the sarting stereo center */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (cur == i01)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;  /*  program error: in this case at least one of the i11, i12 must be == cur */ /*   <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNeighMode == NEIGH_MODE_RING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i11 != i12 && !bCurParityInv1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  failed: the two stereo atoms have not been traversed synchronously */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!at[i11].parity || !at[i12].parity)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /*  another atom does not have parity (it might have been removed) 9-11-2002 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nNeighMode == NEIGH_MODE_CHAIN)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (s1 + s2 != 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  program error: only one out of s1 and s2 must be 1, another must be 0. */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((s1 && !at[i11].parity) || (s2 && !at[i12].parity)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /*  another atom does not have parity (it might have been removed) 9-11-2002 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         parity = pCS->LinearCTStereoCarb[i].parity;
    // INCHI✔️❌:         if (nNeighMode == (NEIGH_MODE_RING && ( i11 != i01 ) && ( i12 != i01 )) ||
    // INCHI✔️❌:              /*  in NEIGH_MODE_RING case we know that i11 == i12 except bCurParityInv1 == 1 */
    // INCHI✔️❌:              nNeighMode == NEIGH_MODE_CHAIN
    // INCHI✔️❌:              /*  in NEIGH_MODE_CHAIN case here we always have 2 different atoms */
    // INCHI✔️❌:            ) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /****************************************************************
    // INCHI✔️❌:             * Case of two transposed atoms or a circular permutation in D4h
    // INCHI✔️❌:             */
    // INCHI✔️❌:             parity1 = s1 ? GetStereoCenterParity( pCG, at, i11, nVisited1 ) : PARITY_IMPOSSIBLE;
    // INCHI✔️❌:             parity2 = s2 ? GetStereoCenterParity( pCG, at, i12, nVisited2 ) : PARITY_IMPOSSIBLE;
    // INCHI✔️❌:             if (!ATOM_PARITY_KNOWN( parity1 ) && !ATOM_PARITY_KNOWN( parity2 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  should not happen: must have been detected at the time of the traversal */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (s1 && s2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bCurParityInv1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int parity1orig = GetStereoCenterParity( pCG, at, i11, nCanonRank );
    // INCHI✔️❌:                     int parity2orig = GetStereoCenterParity( pCG, at, i12, nCanonRank );
    // INCHI✔️❌:                     if (i11 == i12 ||
    // INCHI✔️❌:                         (( parity1 == parity1orig || parity2 == parity2orig || parity1 != parity2 ) &&
    // INCHI✔️❌:                          ATOM_PARITY_WELL_DEF( parity1 )) ||
    // INCHI✔️❌:                          ((parity1 != parity2 && ( !ATOM_PARITY_WELL_DEF( parity1 ))) ||
    // INCHI✔️❌:                              !ATOM_PARITY_WELL_DEF( parity2 ) )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         /*return -1; */ /* should be different atoms with inverted parities */
    // INCHI✔️❌:                         nNumDiff++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i11 != i12 || parity1 != parity2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1; /*  program error: must be the same atom */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             parity12 = s1 ? parity1 : parity2;
    // INCHI✔️❌:
    // INCHI✔️❌:             if (ATOM_PARITY_WELL_DEF( parity ) && parity == parity12)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*  symmetrical neighbors have well-defined equal parities */
    // INCHI✔️❌:                 k++;
    // INCHI✔️❌:                 if (nCheckingMode == CHECKING_STEREOCENTER && nNeighMode == NEIGH_MODE_RING)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  all 3: cur, i01, i11 are different atoms (here i11==i12) */
    // INCHI✔️❌:                     /*  here nSymmRank[i01]==nSymmRank[i11] due to the parallel traversal */
    // INCHI✔️❌:                     if (nSymmRank[cur] == nSymmRank[i01])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nNumEqStereogenic++;  /*  all 3 are equ */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (ATOM_PARITY_WELL_DEF( parity ) && ATOM_PARITY_WELL_DEF( parity12 ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  apparently different well-defined parities */
    // INCHI✔️❌:                     if (!bCurParityInv1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nNumInv++;
    // INCHI✔️❌:                         /* return 0; */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌: #if ( PROPAGATE_ILL_DEF_STEREO == 1 )
    // INCHI✔️❌:                     /*  at least one parity is ill-defined. Use parity1 and parity2 to temporarily save bitmaps */
    // INCHI✔️❌:                     parity1 = ( parity == vABParityUnknown /*AB_PARITY_UNKN*/ )
    // INCHI✔️❌:                         ? NOT_WELL_DEF_UNKN
    // INCHI✔️❌:                         : ( parity == AB_PARITY_UNDF ) ? NOT_WELL_DEF_UNDF : 0;
    // INCHI✔️❌:                     parity2 = ( parity12 == vABParityUnknown /*AB_PARITY_UNKN*/ )
    // INCHI✔️❌:                         ? NOT_WELL_DEF_UNKN : ( parity12 == AB_PARITY_UNDF ) ? NOT_WELL_DEF_UNDF : 0;
    // INCHI✔️❌:                     if (parity1 | parity2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         not_well_def_parities |= ( parity1 | parity2 );
    // INCHI✔️❌:                         k++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1;  /*  program error */ /*   <BRKPT> */
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #else
    // INCHI✔️❌:                     return 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i11 == i01 && i12 == i01)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /********************************************************************/
    // INCHI✔️❌:                 /*  i11 == i12 are same atom as i01, nNeighMode == NEIGH_MODE_RING */
    // INCHI✔️❌:                 if (!s1 || !s2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /*  the parity of the new neighbors permutation must be same as the old one */
    // INCHI✔️❌:                 /*  this must work for well-defined and ill-defined parities. */
    // INCHI✔️❌:                 /*  actual parity (that includes the geometry) is not important here. */
    // INCHI✔️❌:                 /*  old permutation */
    // INCHI✔️❌:                 parity = GetPermutationParity( pCG, at + i01, MAX_ATOMS + 1, nCanonRank );
    // INCHI✔️❌:                 /*  new parmutation */
    // INCHI✔️❌:                 parity1 = GetPermutationParity( pCG, at + i01, MAX_ATOMS + 1, nVisited1 );
    // INCHI✔️❌:                 parity2 = GetPermutationParity( pCG, at + i01, MAX_ATOMS + 1, nVisited2 );
    // INCHI✔️❌:                 if (parity != parity1 || parity != parity2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 k++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* nNeighMode == NEIGH_MODE_RING and only one out of the two (i11 == i01) (i12 == i01) is true */
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* nNumComparedCenters += (k > 0); */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (bCurRotated || nNumDiff || nNumInv)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:      /* !!!! Add here bParitiesInverted == 1 case !!!! */
    // INCHI✔️❌:     /******************************************************/
    // INCHI✔️❌:     /*  Check other (non-starting) stereo bonds/cumulenes */
    // INCHI✔️❌:     /******************************************************/
    // INCHI✔️❌:     for (i = 0; i < pCS->nLenLinearCTStereoDble; i++, nNumComparedBonds += ( k > 0 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         r1 = pCS->LinearCTStereoDble[i].at_num1;
    // INCHI✔️❌:         r2 = pCS->LinearCTStereoDble[i].at_num2;
    // INCHI✔️❌:         i01 = nAtomNumberCanon[r1 - 1];  /*  ord. number of the atom that originally has canon rank r1 */
    // INCHI✔️❌:         i02 = nAtomNumberCanon[r2 - 1];  /*  ord. number of the atom that originally has canon rank r2 */
    // INCHI✔️❌:
    // INCHI✔️❌:         i11 = nAtomNumberCanon1[r1 - 1]; /*  ord. number of the atom that got canon rank r1 during the parallel traversal */
    // INCHI✔️❌:         i12 = nAtomNumberCanon1[r2 - 1]; /*  ord. number of the atom that got canon rank r2 during the parallel traversal */
    // INCHI✔️❌:
    // INCHI✔️❌:         i21 = nAtomNumberCanon2[r1 - 1];
    // INCHI✔️❌:         i22 = nAtomNumberCanon2[r2 - 1];
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         s1 = ( i11 < num_atoms && i12 < num_atoms );
    // INCHI✔️❌:         s2 = ( i21 < num_atoms && i22 < num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:         k = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  check whether the two stereo bonds/allenes (they can be one and the same) have been traversed */
    // INCHI✔️❌:         if (!s1 && !s2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /*  Both stereo bonds/cumulenes have not been traversed; check the next pair. */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nCheckingMode == CHECKING_STEREOBOND)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             switch (( (i11 == cur && i12 == prev_sb_neigh) || (i12 == cur && i11 == prev_sb_neigh) ) +
    // INCHI✔️❌:                 ( (i21 == cur && i22 == prev_sb_neigh) || (i22 == cur && i21 == prev_sb_neigh) )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case 2:
    // INCHI✔️❌:                     continue; /*  do not recheck the starting bond/cumulene */
    // INCHI✔️❌:                 case 1:
    // INCHI✔️❌:                     return -1; /*  possibly program error  */ /*   <BRKPT> */
    // INCHI✔️❌:                 /* case 0: */
    // INCHI✔️❌:                 /*     break; */   /*  the stereo centers are not the sarting stereo center */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (( i01 == cur && i02 == prev_sb_neigh ) || ( i02 == cur && i01 == prev_sb_neigh ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;  /*  program error: in this case at least one of the i1x, i2x must be == cur */ /*   <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNeighMode == NEIGH_MODE_RING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (( i11 != i21 || i12 != i22 ) && ( i11 != i22 || i12 != i21 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  failed: the two bonds/cumulenes have not been traversed synchronously */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (0 > GetStereoNeighborPos( at, i11, i12 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /*  another bond is not stereo (the stereo might have been removed) 9-11-2002 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nNeighMode == NEIGH_MODE_CHAIN)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (s1 + s2 != 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  program error: only one out of s1 and s2 must be 1, another must be 0. */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((s1 && 0 > GetStereoNeighborPos( at, i11, i12 )) ||
    // INCHI✔️❌:                  (s2 && 0 > GetStereoNeighborPos( at, i21, i22 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /*  another bond is not stereo (the stereo might have been removed) 9-11-2002 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         parity = pCS->LinearCTStereoDble[i].parity;
    // INCHI✔️❌:         /* bMustBeIdentical  = ATOM_PARITY_ILL_DEF(parity); */
    // INCHI✔️❌:         /* nNumEqStereogenic = 0; */
    // INCHI✔️❌:
    // INCHI✔️❌:         if ((nNeighMode == NEIGH_MODE_RING && ( i11 != i01 || i12 != i02 ) && ( i11 != i02 || i12 != i01 )) ||
    // INCHI✔️❌:              nNeighMode == NEIGH_MODE_CHAIN                    /*  in NEIGH_MODE_CHAIN case here we always have 2 different atoms */
    // INCHI✔️❌:         ) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:          /*******************************************/
    // INCHI✔️❌:          /*  case of two transposed bonds/cumulenes */
    // INCHI✔️❌:             parity1 = s1 ? GetStereoBondParity( at, i11, i12, nVisited1 ) : PARITY_IMPOSSIBLE;
    // INCHI✔️❌:             parity2 = s2 ? GetStereoBondParity( at, i21, i22, nVisited2 ) : PARITY_IMPOSSIBLE;
    // INCHI✔️❌:             if (!ATOM_PARITY_KNOWN( parity1 ) && !ATOM_PARITY_KNOWN( parity2 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  should not happen: must have been detected at the time of traversal */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (s1 && s2 && ( (( i11 != i21 || i12 != i22 ) && ( i11 != i22 || i12 != i21 )) || parity1 != parity2 )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /*  program error: must be the same bond/cumulene */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             parity12 = s1 ? parity1 : parity2;
    // INCHI✔️❌:             if (ATOM_PARITY_WELL_DEF( parity ) && parity == parity12)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*  symmetrical neighbors have well-defined equal parities */
    // INCHI✔️❌:                 k++;
    // INCHI✔️❌:                 if (nCheckingMode == CHECKING_STEREOBOND && nNeighMode == NEIGH_MODE_RING)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  all 3 bonds: cur-prev_sb_neigh, i01-i02, i11-i12 are different */
    // INCHI✔️❌:                     /*  (here <i11,i12>==<i21,i22> compared as unordered pairs) */
    // INCHI✔️❌:                     if ((nSymmRank[cur] == nSymmRank[i01] && nSymmRank[prev_sb_neigh] == nSymmRank[i02]) ||
    // INCHI✔️❌:                          (nSymmRank[cur] == nSymmRank[i02] && nSymmRank[prev_sb_neigh] == nSymmRank[i01])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nNumEqStereogenic++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (ATOM_PARITY_WELL_DEF( parity ) && ATOM_PARITY_WELL_DEF( parity12 ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  apparently different well-defined parities */
    // INCHI✔️❌:                     return 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  at least one parity is ill-defined. Use parity1 and parity2 to temporarily save bitmaps */
    // INCHI✔️❌: #if ( PROPAGATE_ILL_DEF_STEREO == 1 )
    // INCHI✔️❌:                     parity1 = ( parity == vABParityUnknown /*AB_PARITY_UNKN*/ ) ? NOT_WELL_DEF_UNKN :
    // INCHI✔️❌:                         ( parity == AB_PARITY_UNDF ) ? NOT_WELL_DEF_UNDF : 0;
    // INCHI✔️❌:                     parity2 = ( parity12 == vABParityUnknown /*AB_PARITY_UNKN*/ ) ? NOT_WELL_DEF_UNKN :
    // INCHI✔️❌:                         ( parity12 == AB_PARITY_UNDF ) ? NOT_WELL_DEF_UNDF : 0;
    // INCHI✔️❌:                     if (parity1 | parity2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         not_well_def_parities |= ( parity1 | parity2 );
    // INCHI✔️❌:                         k++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1;  /*  program error */
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #else
    // INCHI✔️❌:                     return 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*****************************************************************************************/
    // INCHI✔️❌:             /*  i11-i12 and i21-i22 are same as i01-i02 bond/cumulene, nNeighMode == NEIGH_MODE_RING */
    // INCHI✔️❌:             AT_NUMB n1, n2;
    // INCHI✔️❌:             int       j;
    // INCHI✔️❌:             if (!s1 || !s2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /*  find neighbors along the stereo bond/cumulene */
    // INCHI✔️❌:             for (j = 0, n1 = MAX_ATOMS + 1; j < MAX_NUM_STEREO_BOND_NEIGH && at[i01].stereo_bond_neighbor[j]; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if ((int) at[i01].stereo_bond_neighbor[j] == i02 + 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     n1 = at[i01].neighbor[(int) at[i01].stereo_bond_ord[j]];
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (j = 0, n2 = MAX_ATOMS + 1; j < MAX_NUM_STEREO_BOND_NEIGH && at[i02].stereo_bond_neighbor[j]; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if ((int) at[i02].stereo_bond_neighbor[j] == i01 + 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     n2 = at[i02].neighbor[(int) at[i02].stereo_bond_ord[j]];
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (n1 > MAX_ATOMS || n2 > MAX_ATOMS)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return CT_REMOVE_STEREO_ERR;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /*  the parity of the new neighbors permutation must be same as the old one */
    // INCHI✔️❌:             /*  this must work for well-defined and ill-defined parities. */
    // INCHI✔️❌:             /*  actual parity (that includes the geometry) is not important here. */
    // INCHI✔️❌:             /*  old permutation */
    // INCHI✔️❌:             parity = GetPermutationParity( pCG, at + i01, n1, nCanonRank ) + GetPermutationParity( pCG, at + i02, n2, nCanonRank );
    // INCHI✔️❌:             /*  new parmutation */
    // INCHI✔️❌:             parity1 = GetPermutationParity( pCG, at + i01, n1, nVisited1 ) + GetPermutationParity( pCG, at + i02, n2, nVisited1 );
    // INCHI✔️❌:             parity2 = GetPermutationParity( pCG, at + i01, n1, nVisited2 ) + GetPermutationParity( pCG, at + i02, n2, nVisited2 );
    // INCHI✔️❌:             if (parity % 2 != parity1 % 2 || parity1 % 2 != parity2 % 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             k++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* nNumComparedBonds += ( k > 0 ); */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumEqStereogenic > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌: /*  case similar to trimethylcyclopropane: 3 constitutionally equivalent stereogenic elements */
    // INCHI✔️❌: /*  the transposition does not change the parities */
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:         pCS->bExtract |= EXTR_2EQL2CENTER_TO_REMOVE_PARITY;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /* =========================================================================================
    // INCHI✔️❌:     Note
    // INCHI✔️❌:     ====
    // INCHI✔️❌:     At this point the comparison is complete and no difference sufficient to establish
    // INCHI✔️❌:     absence of stereo parity has been found.
    // INCHI✔️❌:     However, non-zero not_well_def_parities means that an ill-defined parity was
    // INCHI✔️❌:     compared to an ill-defined or well-defined parity. This means that the parity
    // INCHI✔️❌:     of the atom or bond being checked cannot be well-defined anymore.
    // INCHI✔️❌: ========================================================================================*/
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     not_well_def_parities |= COMP_STEREO_SUCCESS;
    // INCHI✔️❌:
    // INCHI✔️❌:     return not_well_def_parities;
    // INCHI✔️❌:
    // INCHI✔️❌:    /*  Add 1 to indicate success. The stereogenic elements might have been */
    // INCHI✔️❌:    /*  removed while checking existence of the previous atom/bond stereo */
    // INCHI✔️❌:    /* return (nNumComparedCenters + nNumComparedBonds + 1);  */
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /********************************************************************************/
    // INCHI✔️❌: /*  Remove stereo marks from the bonds that are calculated to be non-stereo     */
    // END INCHI C FUNCTION: CalculatedPathsParitiesAreIdentical
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION
    // INCHI✔️❌: #define PROPAGATE_ILL_DEF_STEREO 1
    // INCHI✔️❌: #define bRELEASE_VERSION 1
    // END INCHI ACTIVE MACRO CONFIGURATION

    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let checking_mode = if i32::from(prev_sb_neigh) < num_atoms {
        CHECKING_STEREOBOND as i32
    } else {
        CHECKING_STEREOCENTER as i32
    };
    let mut not_well_def_parities = 0_i32;
    let mut n_num_eq_stereogenic = 0_i32;
    let mut n_num_compared_centers = 0_i32;
    let mut n_num_compared_bonds = 0_i32;
    let mut b_cur_rotated = 0_i32;
    let mut n_num_diff = 0_i32;
    let mut n_num_inv = 0_i32;

    if (nNeighMode != NEIGH_MODE_RING as i32 && bParitiesInverted != 0)
        || bParitiesInverted.wrapping_abs() != 1
    {
        bParitiesInverted = 0;
    }

    let sentinel = (MAX_ATOMS + 1) as AT_RANK;
    if bParitiesInverted != 0 {
        let mut missing1 = 0_i32;
        let mut missing2 = 0_i32;
        for index in 0..count {
            missing1 =
                missing1.wrapping_add(i32::from(source_get(heap, nVisited1, index as i32)? == 0));
            missing2 =
                missing2.wrapping_add(i32::from(source_get(heap, nVisited2, index as i32)? == 0));
            source_set(heap, nAtomNumberCanon1, index as i32, sentinel)?;
            source_set(heap, nAtomNumberCanon2, index as i32, sentinel)?;
        }
        if missing1 != 0 || missing2 != 0 {
            if bParitiesInverted == 1 {
                return Ok(0);
            }
            bParitiesInverted = 0;
        }
    } else {
        for index in 0..count {
            source_set(heap, nAtomNumberCanon1, index as i32, sentinel)?;
            source_set(heap, nAtomNumberCanon2, index as i32, sentinel)?;
        }
    }

    if (bParitiesInverted > 0 && !(mode == MAP_MODE_C2v as i32 || mode == MAP_MODE_S4 as i32))
        || (bParitiesInverted == 0 && !(mode == MAP_MODE_C2 as i32 || mode == MAP_MODE_STD as i32))
    {
        return Ok(0);
    }

    for index in 0..count {
        let visited1 = source_get(heap, nVisited1, index as i32)?;
        if visited1 != 0 {
            let canonical = source_get(heap, nCanonRank, i32::from(visited1) - 1)?;
            source_set(heap, nVisited1, index as i32, canonical)?;
            source_set(
                heap,
                nAtomNumberCanon1,
                i32::from(canonical) - 1,
                index as AT_RANK,
            )?;
        }
        let visited2 = source_get(heap, nVisited2, index as i32)?;
        if visited2 != 0 {
            let canonical = source_get(heap, nCanonRank, i32::from(visited2) - 1)?;
            source_set(heap, nVisited2, index as i32, canonical)?;
            source_set(
                heap,
                nAtomNumberCanon2,
                i32::from(canonical) - 1,
                index as AT_RANK,
            )?;
        }
    }

    let atoms = heap
        .slice(at.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let symm = heap
        .slice(nSymmRank.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let canonical = heap
        .slice(nCanonRank.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let atom_by_canon = heap
        .slice(nAtomNumberCanon.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let atom_by_canon1 = heap
        .slice(nAtomNumberCanon1.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let atom_by_canon2 = heap
        .slice(nAtomNumberCanon2.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let visited1 = heap
        .slice(nVisited1.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let visited2 = heap
        .slice(nVisited2.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let atom_index = |value: AT_RANK| -> Result<usize, SourceHeapError> {
        let index = usize::from(value);
        if index < count {
            Ok(index)
        } else {
            Err(SourceHeapError::PointerOutOfBounds)
        }
    };
    let rank_index = |rank: AT_RANK| -> Result<usize, SourceHeapError> {
        usize::from(rank)
            .checked_sub(1)
            .filter(|index| *index < count)
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };
    let atom_known =
        |parity: i32| (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&parity);
    let atom_well_defined = |parity: i32| {
        (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32).contains(&parity)
    };

    let cur_index = atom_index(cur)?;
    if checking_mode == CHECKING_STEREOBOND as i32 {
        let opposite = atom_index(prev_sb_neigh)?;
        if visited1[opposite] != 0 || visited2[opposite] != 0 {
            if visited1[opposite] != visited2[opposite] || canonical[opposite] != visited2[opposite]
            {
                return Ok(0);
            }
            if i32::from(atoms[opposite].valence).wrapping_add(i32::from(atoms[opposite].num_H)) > 3
            {
                return Ok(0);
            }
            let mut slot = 0_usize;
            let mut found = false;
            while slot < MAX_NUM_STEREO_BONDS as usize {
                let neighbor = atoms[opposite].stereo_bond_neighbor[slot];
                if neighbor == 0 {
                    break;
                }
                if neighbor.wrapping_sub(1) == cur {
                    found = true;
                    break;
                }
                slot += 1;
            }
            if !found {
                return Ok(-1);
            }
            let bond_order = usize::try_from(atoms[opposite].stereo_bond_ord[slot])
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut other = [0_u16; 2];
            let mut other_count = 0_usize;
            let mut index = 0_i32;
            while index < i32::from(atoms[opposite].valence) && other_count <= 2 {
                if index as usize != bond_order {
                    if other_count >= other.len() {
                        return Ok(CT_STEREOCOUNT_ERR);
                    }
                    other[other_count] = atoms[opposite].neighbor[index as usize];
                    other_count += 1;
                }
                index = index.wrapping_add(1);
            }
            if other_count as i32 + i32::from(atoms[opposite].num_H) > 2 {
                return Ok(CT_STEREOCOUNT_ERR);
            }
            for neighbor in other.into_iter().take(other_count) {
                let index = atom_index(neighbor)?;
                if (visited1[index] != 0 && visited1[index] != canonical[index])
                    || (visited2[index] != 0 && visited2[index] != canonical[index])
                {
                    return Ok(0);
                }
            }
        }
    }

    if checking_mode == CHECKING_STEREOCENTER as i32 {
        let mut other = [0_u16; 2];
        let mut other_count = 0_usize;
        let mut index = 0_i32;
        while index < i32::from(atoms[cur_index].valence) && other_count <= 2 {
            let neighbor = atoms[cur_index].neighbor[index as usize];
            if neighbor != next1 && neighbor != next2 {
                if other_count >= other.len() {
                    return Ok(CT_STEREOCOUNT_ERR);
                }
                other[other_count] = neighbor;
                other_count += 1;
            }
            index = index.wrapping_add(1);
        }
        if other_count as i32 + i32::from(atoms[cur_index].num_H) > 2 {
            return Ok(CT_STEREOCOUNT_ERR);
        }
        if bParitiesInverted != 0
            && i32::from(atoms[cur_index].valence) == MAX_NUM_STEREO_ATOM_NEIGH as i32
        {
            let atom_pointer = at.offset(i64::from(cur))?;
            let canon_parity = GetPermutationParity(heap, pCG, atom_pointer, sentinel, nCanonRank)?;
            let parity1 = GetPermutationParity(heap, pCG, atom_pointer, sentinel, nVisited1)?;
            let parity2 = GetPermutationParity(heap, pCG, atom_pointer, sentinel, nVisited2)?;
            if parity1 != parity2 {
                return Ok(0);
            }
            if bParitiesInverted == 1 && parity1 == canon_parity {
                return Ok(0);
            }
            if bParitiesInverted == -1 {
                bParitiesInverted = i32::from(parity1 != canon_parity);
            }
        }
        if bParitiesInverted == 0 && b_cur_rotated == 0 {
            for neighbor in other.into_iter().take(other_count) {
                let index = atom_index(neighbor)?;
                if (visited1[index] != 0 && visited1[index] != canonical[index])
                    || (visited2[index] != 0 && visited2[index] != canonical[index])
                {
                    return Ok(0);
                }
            }
        }
    }

    let center_count = usize::try_from(pCS.nLenLinearCTStereoCarb)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let centers = if center_count == 0 {
        Vec::new()
    } else {
        heap.slice(pCS.LinearCTStereoCarb.as_const())?
            .get(..center_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    for center in centers {
        let rank_slot = rank_index(center.at_num)?;
        let original = usize::from(atom_by_canon[rank_slot]);
        let mapped1 = usize::from(atom_by_canon1[rank_slot]);
        let mapped2 = usize::from(atom_by_canon2[rank_slot]);
        let traversed1 = mapped1 < count;
        let traversed2 = mapped2 < count;
        let parity_inverted_here = bParitiesInverted != 0
            && atoms[cur_index].nRingSystem
                == atoms
                    .get(mapped1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nRingSystem
            && atoms[cur_index].nRingSystem
                == atoms
                    .get(mapped2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nRingSystem;
        let mut compared = false;
        if !traversed1 && !traversed2 {
            continue;
        }
        if checking_mode == CHECKING_STEREOCENTER as i32 {
            match i32::from(cur_index == mapped1) + i32::from(cur_index == mapped2) {
                2 => continue,
                1 => return Ok(-1),
                _ => {}
            }
            if cur_index == original {
                return Ok(-1);
            }
        }
        if nNeighMode == NEIGH_MODE_RING as i32 {
            if mapped1 != mapped2 && !parity_inverted_here {
                return Ok(-1);
            }
            if atoms
                .get(mapped1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .parity
                == 0
                || atoms
                    .get(mapped2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .parity
                    == 0
            {
                return Ok(0);
            }
        }
        if nNeighMode == NEIGH_MODE_CHAIN as i32 {
            if i32::from(traversed1) + i32::from(traversed2) != 1 {
                return Ok(-1);
            }
            if (traversed1 && atoms[mapped1].parity == 0)
                || (traversed2 && atoms[mapped2].parity == 0)
            {
                return Ok(0);
            }
        }

        let stored_parity = i32::from(center.parity);
        let transposed_case = nNeighMode
            == i32::from(NEIGH_MODE_RING != 0 && mapped1 != original && mapped2 != original)
            || nNeighMode == NEIGH_MODE_CHAIN as i32;
        if transposed_case {
            let parity1 = if traversed1 {
                GetStereoCenterParity(heap, pCG, at, mapped1 as i32, nVisited1)?
            } else {
                PARITY_IMPOSSIBLE as i32
            };
            let parity2 = if traversed2 {
                GetStereoCenterParity(heap, pCG, at, mapped2 as i32, nVisited2)?
            } else {
                PARITY_IMPOSSIBLE as i32
            };
            if !atom_known(parity1) && !atom_known(parity2) {
                return Ok(-1);
            }
            if traversed1 && traversed2 {
                if parity_inverted_here {
                    let original1 =
                        GetStereoCenterParity(heap, pCG, at, mapped1 as i32, nCanonRank)?;
                    let original2 =
                        GetStereoCenterParity(heap, pCG, at, mapped2 as i32, nCanonRank)?;
                    if mapped1 == mapped2
                        || ((parity1 == original1 || parity2 == original2 || parity1 != parity2)
                            && atom_well_defined(parity1))
                        || (parity1 != parity2 && !atom_well_defined(parity1))
                        || !atom_well_defined(parity2)
                    {
                        n_num_diff = n_num_diff.wrapping_add(1);
                    }
                } else if mapped1 != mapped2 || parity1 != parity2 {
                    return Ok(-1);
                }
            }
            let mapped_parity = if traversed1 { parity1 } else { parity2 };
            if atom_well_defined(stored_parity) && stored_parity == mapped_parity {
                compared = true;
                if checking_mode == CHECKING_STEREOCENTER as i32
                    && nNeighMode == NEIGH_MODE_RING as i32
                    && symm[cur_index] == symm[original]
                {
                    n_num_eq_stereogenic = n_num_eq_stereogenic.wrapping_add(1);
                }
            } else if atom_well_defined(stored_parity) && atom_well_defined(mapped_parity) {
                if !parity_inverted_here {
                    n_num_inv = n_num_inv.wrapping_add(1);
                }
            } else {
                let bitmap1 = if stored_parity == vABParityUnknown {
                    NOT_WELL_DEF_UNKN as i32
                } else if stored_parity == AB_PARITY_UNDF as i32 {
                    NOT_WELL_DEF_UNDF as i32
                } else {
                    0
                };
                let bitmap2 = if mapped_parity == vABParityUnknown {
                    NOT_WELL_DEF_UNKN as i32
                } else if mapped_parity == AB_PARITY_UNDF as i32 {
                    NOT_WELL_DEF_UNDF as i32
                } else {
                    0
                };
                if bitmap1 | bitmap2 != 0 {
                    not_well_def_parities |= bitmap1 | bitmap2;
                    compared = true;
                } else {
                    return Ok(-1);
                }
            }
        } else if mapped1 == original && mapped2 == original {
            if !traversed1 || !traversed2 {
                return Ok(-1);
            }
            let atom_pointer = at.offset(original as i64)?;
            let old = GetPermutationParity(heap, pCG, atom_pointer, sentinel, nCanonRank)?;
            let new1 = GetPermutationParity(heap, pCG, atom_pointer, sentinel, nVisited1)?;
            let new2 = GetPermutationParity(heap, pCG, atom_pointer, sentinel, nVisited2)?;
            if old != new1 || old != new2 {
                return Ok(0);
            }
            compared = true;
        } else {
            return Ok(-1);
        }
        n_num_compared_centers = n_num_compared_centers.wrapping_add(i32::from(compared));
    }
    if b_cur_rotated != 0 || n_num_diff != 0 || n_num_inv != 0 {
        return Ok(0);
    }

    let bond_count = usize::try_from(pCS.nLenLinearCTStereoDble)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let bonds = if bond_count == 0 {
        Vec::new()
    } else {
        heap.slice(pCS.LinearCTStereoDble.as_const())?
            .get(..bond_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    for bond in bonds {
        let rank1_slot = rank_index(bond.at_num1)?;
        let rank2_slot = rank_index(bond.at_num2)?;
        let original1 = usize::from(atom_by_canon[rank1_slot]);
        let original2 = usize::from(atom_by_canon[rank2_slot]);
        let mapped11 = usize::from(atom_by_canon1[rank1_slot]);
        let mapped12 = usize::from(atom_by_canon1[rank2_slot]);
        let mapped21 = usize::from(atom_by_canon2[rank1_slot]);
        let mapped22 = usize::from(atom_by_canon2[rank2_slot]);
        let traversed1 = mapped11 < count && mapped12 < count;
        let traversed2 = mapped21 < count && mapped22 < count;
        let mut compared = false;
        if !traversed1 && !traversed2 {
            continue;
        }
        if checking_mode == CHECKING_STEREOBOND as i32 {
            let first_start = (mapped11 == cur_index && mapped12 == usize::from(prev_sb_neigh))
                || (mapped12 == cur_index && mapped11 == usize::from(prev_sb_neigh));
            let second_start = (mapped21 == cur_index && mapped22 == usize::from(prev_sb_neigh))
                || (mapped22 == cur_index && mapped21 == usize::from(prev_sb_neigh));
            match i32::from(first_start) + i32::from(second_start) {
                2 => continue,
                1 => return Ok(-1),
                _ => {}
            }
            if (original1 == cur_index && original2 == usize::from(prev_sb_neigh))
                || (original2 == cur_index && original1 == usize::from(prev_sb_neigh))
            {
                return Ok(-1);
            }
        }
        if nNeighMode == NEIGH_MODE_RING as i32 {
            if (mapped11 != mapped21 || mapped12 != mapped22)
                && (mapped11 != mapped22 || mapped12 != mapped21)
            {
                return Ok(-1);
            }
            if GetStereoNeighborPos(heap, at, mapped11 as i32, mapped12 as i32)? < 0 {
                return Ok(0);
            }
        }
        if nNeighMode == NEIGH_MODE_CHAIN as i32 {
            if i32::from(traversed1) + i32::from(traversed2) != 1 {
                return Ok(-1);
            }
            if (traversed1 && GetStereoNeighborPos(heap, at, mapped11 as i32, mapped12 as i32)? < 0)
                || (traversed2
                    && GetStereoNeighborPos(heap, at, mapped21 as i32, mapped22 as i32)? < 0)
            {
                return Ok(0);
            }
        }

        let stored_parity = i32::from(bond.parity);
        if (nNeighMode == NEIGH_MODE_RING as i32
            && (mapped11 != original1 || mapped12 != original2)
            && (mapped11 != original2 || mapped12 != original1))
            || nNeighMode == NEIGH_MODE_CHAIN as i32
        {
            let parity1 = if traversed1 {
                GetStereoBondParity(heap, at, mapped11 as i32, mapped12 as i32, nVisited1)?
            } else {
                PARITY_IMPOSSIBLE as i32
            };
            let parity2 = if traversed2 {
                GetStereoBondParity(heap, at, mapped21 as i32, mapped22 as i32, nVisited2)?
            } else {
                PARITY_IMPOSSIBLE as i32
            };
            if !atom_known(parity1) && !atom_known(parity2) {
                return Ok(-1);
            }
            if traversed1
                && traversed2
                && (((mapped11 != mapped21 || mapped12 != mapped22)
                    && (mapped11 != mapped22 || mapped12 != mapped21))
                    || parity1 != parity2)
            {
                return Ok(-1);
            }
            let mapped_parity = if traversed1 { parity1 } else { parity2 };
            if atom_well_defined(stored_parity) && stored_parity == mapped_parity {
                compared = true;
                if checking_mode == CHECKING_STEREOBOND as i32
                    && nNeighMode == NEIGH_MODE_RING as i32
                    && ((symm[cur_index] == symm[original1]
                        && symm[usize::from(prev_sb_neigh)] == symm[original2])
                        || (symm[cur_index] == symm[original2]
                            && symm[usize::from(prev_sb_neigh)] == symm[original1]))
                {
                    n_num_eq_stereogenic = n_num_eq_stereogenic.wrapping_add(1);
                }
            } else if atom_well_defined(stored_parity) && atom_well_defined(mapped_parity) {
                return Ok(0);
            } else {
                let bitmap1 = if stored_parity == vABParityUnknown {
                    NOT_WELL_DEF_UNKN as i32
                } else if stored_parity == AB_PARITY_UNDF as i32 {
                    NOT_WELL_DEF_UNDF as i32
                } else {
                    0
                };
                let bitmap2 = if mapped_parity == vABParityUnknown {
                    NOT_WELL_DEF_UNKN as i32
                } else if mapped_parity == AB_PARITY_UNDF as i32 {
                    NOT_WELL_DEF_UNDF as i32
                } else {
                    0
                };
                if bitmap1 | bitmap2 != 0 {
                    not_well_def_parities |= bitmap1 | bitmap2;
                    compared = true;
                } else {
                    return Ok(-1);
                }
            }
        } else {
            if !traversed1 || !traversed2 {
                return Ok(-1);
            }
            let mut avoid1 = sentinel;
            let mut avoid2 = sentinel;
            for slot in 0..MAX_NUM_STEREO_BOND_NEIGH as usize {
                let neighbor = atoms[original1].stereo_bond_neighbor[slot];
                if neighbor == 0 {
                    break;
                }
                if usize::from(neighbor.wrapping_sub(1)) == original2 {
                    let order = usize::try_from(atoms[original1].stereo_bond_ord[slot])
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    avoid1 = atoms[original1].neighbor[order];
                    break;
                }
            }
            for slot in 0..MAX_NUM_STEREO_BOND_NEIGH as usize {
                let neighbor = atoms[original2].stereo_bond_neighbor[slot];
                if neighbor == 0 {
                    break;
                }
                if usize::from(neighbor.wrapping_sub(1)) == original1 {
                    let order = usize::try_from(atoms[original2].stereo_bond_ord[slot])
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    avoid2 = atoms[original2].neighbor[order];
                    break;
                }
            }
            if avoid1 > MAX_ATOMS as AT_RANK || avoid2 > MAX_ATOMS as AT_RANK {
                return Ok(CT_REMOVE_STEREO_ERR);
            }
            let atom1 = at.offset(original1 as i64)?;
            let atom2 = at.offset(original2 as i64)?;
            let old = GetPermutationParity(heap, pCG, atom1, avoid1, nCanonRank)?
                .wrapping_add(GetPermutationParity(heap, pCG, atom2, avoid2, nCanonRank)?);
            let new1 = GetPermutationParity(heap, pCG, atom1, avoid1, nVisited1)?
                .wrapping_add(GetPermutationParity(heap, pCG, atom2, avoid2, nVisited1)?);
            let new2 = GetPermutationParity(heap, pCG, atom1, avoid1, nVisited2)?
                .wrapping_add(GetPermutationParity(heap, pCG, atom2, avoid2, nVisited2)?);
            if old.wrapping_rem(2) != new1.wrapping_rem(2)
                || new1.wrapping_rem(2) != new2.wrapping_rem(2)
            {
                return Ok(0);
            }
            compared = true;
        }
        n_num_compared_bonds = n_num_compared_bonds.wrapping_add(i32::from(compared));
    }

    if n_num_eq_stereogenic > 0 {
        return Ok(0);
    }
    not_well_def_parities |= COMP_STEREO_SUCCESS as i32;
    let _ = (n_num_compared_centers, n_num_compared_bonds);
    Ok(not_well_def_parities)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RemoveCalculatedNonStereoBondParities(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    num_at_tg: i32,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nTempRank: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon1: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon2: SourceMutPointer<AT_RANK>,
    nl: SourceMutPointer<NEIGH_LIST>,
    nl1: SourceMutPointer<NEIGH_LIST>,
    nl2: SourceMutPointer<NEIGH_LIST>,
    nVisited1: SourceMutPointer<AT_RANK>,
    nVisited2: SourceMutPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3099 RemoveCalculatedNonStereoBondParities
    // INCHI✔️❌: int RemoveCalculatedNonStereoBondParities( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                            sp_ATOM *at, int num_atoms,
    // INCHI✔️❌:                                            int num_at_tg,
    // INCHI✔️❌:                                            AT_RANK **pRankStack1,
    // INCHI✔️❌:                                            AT_RANK **pRankStack2,
    // INCHI✔️❌:                                            AT_RANK *nTempRank,
    // INCHI✔️❌:                                            NEIGH_LIST *NeighList,
    // INCHI✔️❌:                                            AT_RANK *nCanonRank,
    // INCHI✔️❌:                                            const AT_RANK *nSymmRank,
    // INCHI✔️❌:                                            AT_RANK *nAtomNumberCanon,
    // INCHI✔️❌:                                            AT_RANK *nAtomNumberCanon1,
    // INCHI✔️❌:                                            AT_RANK *nAtomNumberCanon2,
    // INCHI✔️❌:                                            NEIGH_LIST *nl,
    // INCHI✔️❌:                                            NEIGH_LIST *nl1,
    // INCHI✔️❌:                                            NEIGH_LIST *nl2,
    // INCHI✔️❌:                                            AT_RANK *nVisited1,
    // INCHI✔️❌:                                            AT_RANK *nVisited2,
    // INCHI✔️❌:                                            CANON_STAT *pCS,
    // INCHI✔️❌:                                            int vABParityUnknown )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, n, m, ret, ret1, ret2, ret_failed = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     int i1, n1, s2;  /*  n1 must be SIGNED integer */
    // INCHI✔️❌:     AT_RANK nAtomRank1, nAtomRank2, neigh[3] = { 0 }, nAvoidCheckAtom[2], opposite_atom, nLength; /* djb-rwth: initialisation of neigh required to avoid undefined array subscript */
    // INCHI✔️❌:     int         nNeighMode = NEIGH_MODE_CHAIN;
    // INCHI✔️❌:     int         nNumEqRingNeigh = 0, bRingNeigh, bSymmNeigh, bParitiesInverted;
    // INCHI✔️❌:     NEIGH_LIST *nl01, *nl02;
    // INCHI✔️❌:     const AT_RANK    *nSymmRank1, *nSymmRank2;
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: second_pass:
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i1 = 0; i1 < num_atoms && !RETURNED_ERROR(ret_failed); i1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[i1].valence != 3 || !at[i1].stereo_bond_neighbor[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (n1 = 0; n1 < MAX_NUM_STEREO_BONDS && !RETURNED_ERROR(ret_failed) && (s2 = at[i1].stereo_bond_neighbor[n1]); n1++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!PARITY_CALCULATE(at[i1].stereo_bond_parity[n1]) && PARITY_WELL_DEF(at[i1].stereo_bond_parity[n1]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             opposite_atom = (AT_RANK)(s2 - 1);
    // INCHI✔️❌:             s2 = at[i1].neighbor[(int)at[i1].stereo_bond_ord[n1]]; /*  different from opposite_atom in case of a cumulene */
    // INCHI✔️❌:             for (j = 1, n = 0; j <= (int)at[i1].valence; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nl[i1][j] != s2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     neigh[n++] = nl[i1][j]; /*  sorting guarantees that canon. rank of neigh[0] is greater or equal */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (n != 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nSymmRank[(int)neigh[0]] != nSymmRank[(int)neigh[1]])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue; /*  may happen if another half-bond has not a defined parity */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             bRingNeigh = (at[(int)neigh[0]].nRingSystem == at[(int)neigh[1]].nRingSystem);
    // INCHI✔️❌:             switch (nNeighMode)
    // INCHI✔️❌:             {
    // INCHI✔️❌:             case NEIGH_MODE_CHAIN:
    // INCHI✔️❌:                 if (bRingNeigh)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nNumEqRingNeigh++;
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 nl01 = nl;
    // INCHI✔️❌:                 nl02 = nl;
    // INCHI✔️❌:                 nSymmRank1 = nSymmRank;
    // INCHI✔️❌:                 nSymmRank2 = nSymmRank;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             case NEIGH_MODE_RING:
    // INCHI✔️❌:                 if (!bRingNeigh)
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 /*  break a tie between the two contitutionally equivalent neighbors, */
    // INCHI✔️❌:                 /*  refine the two partitions, sort neighbors lists nl1, nl2 */
    // INCHI✔️❌:
    // INCHI✔️❌:                 bSymmNeigh = BreakNeighborsTie(pCG,
    // INCHI✔️❌:                     at, num_atoms, num_at_tg,
    // INCHI✔️❌:                     opposite_atom, i1,
    // INCHI✔️❌:                     neigh, 0, 1, 0,
    // INCHI✔️❌:                     pRankStack1, pRankStack2, nTempRank,
    // INCHI✔️❌:                     NeighList, nSymmRank, nCanonRank,
    // INCHI✔️❌:                     nl1, nl2, &pCS->lNumNeighListIter);
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (bSymmNeigh <= 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ret_failed > bSymmNeigh)
    // INCHI✔️❌:                         ret_failed = bSymmNeigh;
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 nl01 = nl1;
    // INCHI✔️❌:                 nl02 = nl2;
    // INCHI✔️❌:                 nSymmRank1 = pRankStack1[0];
    // INCHI✔️❌:                 nSymmRank2 = pRankStack2[0];
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 return CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /*  initialize arrays */
    // INCHI✔️❌:             memset(nVisited1, 0, sizeof(nVisited1[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             memset(nVisited2, 0, sizeof(nVisited2[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             memset(nAtomNumberCanon1, 0, sizeof(nAtomNumberCanon1[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             memset(nAtomNumberCanon2, 0, sizeof(nAtomNumberCanon2[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             nLength = 1;
    // INCHI✔️❌:             nVisited1[i1] = i1 + 1;   /*  start atoms are the same */
    // INCHI✔️❌:             nVisited2[i1] = i1 + 1;
    // INCHI✔️❌:             nAtomNumberCanon1[i1] = nLength;
    // INCHI✔️❌:             nAtomNumberCanon2[i1] = nLength;
    // INCHI✔️❌:             nAvoidCheckAtom[0] = i1;
    // INCHI✔️❌:             nAvoidCheckAtom[1] = opposite_atom;
    // INCHI✔️❌:             bParitiesInverted = (nNeighMode == NEIGH_MODE_RING &&
    // INCHI✔️❌:                 IS_ALLENE_CHAIN(at[i1].stereo_bond_parity[n1]) &&
    // INCHI✔️❌:                 PARITY_CALCULATE(at[i1].stereo_bond_parity[n1]) &&
    // INCHI✔️❌:                 at[i1].nRingSystem == at[opposite_atom].nRingSystem &&
    // INCHI✔️❌:                 at[opposite_atom].valence == MAX_NUM_STEREO_BONDS) ? -1 : 0;
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:             /* djb-rwth: restructuring code to avoid garbage values -- discussion required */
    // INCHI✔️❌:
    // INCHI✔️❌:             ret1 = CreateCheckSymmPaths(at, (AT_RANK)i1, neigh[0], (AT_RANK)i1, neigh[1], nAvoidCheckAtom,
    // INCHI✔️❌:                 nVisited1, nVisited2, nAtomNumberCanon1, nAtomNumberCanon2,
    // INCHI✔️❌:                 nl01, nl02, nSymmRank1, nSymmRank2, nCanonRank, &nLength, &bParitiesInverted, 0);
    // INCHI✔️❌:             ret2 = CalculatedPathsParitiesAreIdentical(pCG, at, num_atoms, nSymmRank,
    // INCHI✔️❌:                 nCanonRank, nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
    // INCHI✔️❌:                 nVisited1, nVisited2, opposite_atom, (AT_RANK)i1,
    // INCHI✔️❌:                 neigh[0], neigh[1], nNeighMode, bParitiesInverted, 0,
    // INCHI✔️❌:                 pCS, vABParityUnknown);
    // INCHI✔️❌:
    // INCHI✔️❌:             if (0 < ret1) /* djb-rwth: or is this (0 < ret1) && (0 < ret2) ? */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (0 < ret2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ret2 & (NOT_WELL_DEF_UNKN | NOT_WELL_DEF_UNDF))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  possibly change the parity to unknown or undefined */
    // INCHI✔️❌:                             int new_parity = (ret2 & NOT_WELL_DEF_UNKN) ? vABParityUnknown /*AB_PARITY_UNKN*/ : AB_PARITY_UNDF;
    // INCHI✔️❌:                             if ((PARITY_ILL_DEF(at[i1].stereo_bond_parity[n1]) && PARITY_VAL(at[i1].stereo_bond_parity[n1]) > new_parity) ||
    // INCHI✔️❌:                                 PARITY_CALCULATE(at[i1].stereo_bond_parity[n1])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /*  set new unknown or undefined parity */
    // INCHI✔️❌:                                 SetOneStereoBondIllDefParity(at, i1, /* atom number*/ n1 /* stereo bond ord. number*/, new_parity);
    // INCHI✔️❌:                                 /*  change in pCS */
    // INCHI✔️❌:                                 nAtomRank1 = inchi_max(nCanonRank[i1], nCanonRank[opposite_atom]);
    // INCHI✔️❌:                                 nAtomRank2 = inchi_min(nCanonRank[i1], nCanonRank[opposite_atom]);
    // INCHI✔️❌:                                 for (n = 0, m = pCS->nLenLinearCTStereoDble - 1; n <= m; n++)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (pCS->LinearCTStereoDble[n].at_num1 == nAtomRank1 &&
    // INCHI✔️❌:                                         pCS->LinearCTStereoDble[n].at_num2 == nAtomRank2)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         pCS->LinearCTStereoDble[n].parity = new_parity;
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:                                         pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                         m = -1;
    // INCHI✔️❌:                                         break;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (m >= 0)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
    // INCHI✔️❌:                                     goto exit_function;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 ret++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  remove the parity */
    // INCHI✔️❌:                             if (!RemoveOneStereoBond(at, i1, /* atom number*/ n1 /* stereo bond ord. number*/))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
    // INCHI✔️❌:                                 goto exit_function;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             n1--;  /*  cycle counter may temporarily become negative */
    // INCHI✔️❌:                             /*  Remove from the pCS */
    // INCHI✔️❌:                             nAtomRank1 = inchi_max(nCanonRank[i1], nCanonRank[opposite_atom]);
    // INCHI✔️❌:                             nAtomRank2 = inchi_min(nCanonRank[i1], nCanonRank[opposite_atom]);
    // INCHI✔️❌:                             for (n = 0, m = pCS->nLenLinearCTStereoDble - 1; n <= m; n++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (pCS->LinearCTStereoDble[n].at_num1 == nAtomRank1 &&
    // INCHI✔️❌:                                     pCS->LinearCTStereoDble[n].at_num2 == nAtomRank2)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (n < m)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         /*  remove pCS->LinearCTStereoDble[n] */
    // INCHI✔️❌:                                         memmove(pCS->LinearCTStereoDble + n,
    // INCHI✔️❌:                                             pCS->LinearCTStereoDble + n + 1,
    // INCHI✔️❌:                                             ((long long)m - (long long)n) * sizeof(pCS->LinearCTStereoDble[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     pCS->nLenLinearCTStereoDble--;
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:                                     pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                     m = -1;
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (m >= 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
    // INCHI✔️❌:                                 goto exit_function;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             ret++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!ret_failed)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret_failed = (ret1 < 0) ? ret1 : (ret2 < 0) ? ret2 : 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (!RETURNED_ERROR(ret_failed))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (RETURNED_ERROR(ret1))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ret_failed = ret1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (RETURNED_ERROR(ret2))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ret_failed = ret2;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNeighMode == NEIGH_MODE_CHAIN && nNumEqRingNeigh && !RETURNED_ERROR( ret_failed ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nNeighMode = NEIGH_MODE_RING;
    // INCHI✔️❌:         goto second_pass;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return RETURNED_ERROR( ret_failed ) ? ret_failed : ret_failed ? -( ret_failed + 1 ) : ret;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /****************************************************************************/
    // END INCHI C FUNCTION: RemoveCalculatedNonStereoBondParities
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION
    // INCHI✔️❌: #define REMOVE_CALC_NONSTEREO 1
    // INCHI✔️❌: #define bRELEASE_VERSION 1
    // END INCHI ACTIVE MACRO CONFIGURATION

    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let returned_error = |value: i32| CT_ERR_MIN <= value && value <= CT_ERR_MAX;
    let parity_value = |value: i8| i32::from(value) & BITS_PARITY as i32;
    let parity_calculate = |value: i8| parity_value(value) == AB_PARITY_CALC as i32;
    let parity_well_defined = |value: i8| {
        (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
            .contains(&parity_value(value))
    };
    let parity_ill_defined = |value: i8| {
        (AB_MIN_ILL_DEFINED_PARITY as i32..=AB_MAX_ILL_DEFINED_PARITY as i32)
            .contains(&parity_value(value))
    };
    let is_allene_chain =
        |value: i8| ((u32::from(value as u8) & MASK_CUMULENE_LEN) / MULT_STEREOBOND) % 2 != 0;

    let mut ret = 0_i32;
    let mut ret_failed = 0_i32;
    let mut neigh_mode = NEIGH_MODE_CHAIN as i32;
    let mut num_equal_ring_neighbors = 0_i32;

    'passes: loop {
        let mut atom_number = 0_i32;
        while atom_number < num_atoms && !returned_error(ret_failed) {
            let atom_index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom_snapshot = heap
                .slice(at.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if atom_snapshot.valence != 3 || atom_snapshot.stereo_bond_neighbor[0] == 0 {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }

            let mut stereo_index = 0_i32;
            while stereo_index < MAX_NUM_STEREO_BONDS as i32 && !returned_error(ret_failed) {
                let atom_snapshot = heap
                    .slice(at.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let slot = usize::try_from(stereo_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let stereo_neighbor = atom_snapshot.stereo_bond_neighbor[slot];
                if stereo_neighbor == 0 {
                    break;
                }
                let bond_parity = atom_snapshot.stereo_bond_parity[slot];
                if !parity_calculate(bond_parity) && parity_well_defined(bond_parity) {
                    stereo_index = stereo_index.wrapping_add(1);
                    continue;
                }

                let opposite_atom = stereo_neighbor.wrapping_sub(1);
                let bond_order = usize::try_from(atom_snapshot.stereo_bond_ord[slot])
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let excluded_neighbor = *atom_snapshot
                    .neighbor
                    .get(bond_order)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let row = source_get(heap, nl, atom_number)?;
                let mut neighbors = [0_u16; 3];
                let mut neighbor_count = 0_usize;
                let mut position = 1_i32;
                while position <= i32::from(atom_snapshot.valence) {
                    let neighbor = source_get(heap, row, position)?;
                    if neighbor != excluded_neighbor {
                        if neighbor_count >= neighbors.len() {
                            ret = CT_STEREOBOND_ERROR;
                            break 'passes;
                        }
                        neighbors[neighbor_count] = neighbor;
                        neighbor_count += 1;
                    }
                    position = position.wrapping_add(1);
                }
                if neighbor_count != 2 {
                    ret = CT_STEREOBOND_ERROR;
                    break 'passes;
                }
                if source_get(heap, nSymmRank, i32::from(neighbors[0]))?
                    != source_get(heap, nSymmRank, i32::from(neighbors[1]))?
                {
                    stereo_index = stereo_index.wrapping_add(1);
                    continue;
                }

                let ring_neighbors = {
                    let atoms = heap.slice(at.as_const())?;
                    let first = atoms
                        .get(usize::from(neighbors[0]))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let second = atoms
                        .get(usize::from(neighbors[1]))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    first.nRingSystem == second.nRingSystem
                };

                let (path_list1, path_list2, symm_rank1, symm_rank2) =
                    if neigh_mode == NEIGH_MODE_CHAIN as i32 {
                        if ring_neighbors {
                            num_equal_ring_neighbors = num_equal_ring_neighbors.wrapping_add(1);
                            stereo_index = stereo_index.wrapping_add(1);
                            continue;
                        }
                        (nl, nl, nSymmRank, nSymmRank)
                    } else if neigh_mode == NEIGH_MODE_RING as i32 {
                        if !ring_neighbors {
                            stereo_index = stereo_index.wrapping_add(1);
                            continue;
                        }
                        let neighbor_storage = heap.allocate_model_storage(neighbors.to_vec())?;
                        let symmetric = BreakNeighborsTie(
                            heap,
                            pCG,
                            at,
                            num_atoms,
                            num_at_tg,
                            i32::from(opposite_atom),
                            atom_number,
                            neighbor_storage,
                            0,
                            1,
                            0,
                            pRankStack1,
                            pRankStack2,
                            nTempRank,
                            NeighList,
                            nSymmRank,
                            nCanonRank,
                            nl1,
                            nl2,
                            &mut pCS.lNumNeighListIter,
                        )?;
                        if symmetric <= 0 {
                            if ret_failed > symmetric {
                                ret_failed = symmetric;
                            }
                            stereo_index = stereo_index.wrapping_add(1);
                            continue;
                        }
                        (
                            nl1,
                            nl2,
                            source_get(heap, pRankStack1, 0)?,
                            source_get(heap, pRankStack2, 0)?,
                        )
                    } else {
                        return Ok(CT_STEREOCOUNT_ERR);
                    };

                heap.slice_mut(nVisited1)?
                    .get_mut(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .fill(0);
                heap.slice_mut(nVisited2)?
                    .get_mut(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .fill(0);
                heap.slice_mut(nAtomNumberCanon1)?
                    .get_mut(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .fill(0);
                heap.slice_mut(nAtomNumberCanon2)?
                    .get_mut(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .fill(0);

                let mut length = 1_u16;
                source_set(
                    heap,
                    nVisited1,
                    atom_number,
                    atom_number.wrapping_add(1) as AT_RANK,
                )?;
                source_set(
                    heap,
                    nVisited2,
                    atom_number,
                    atom_number.wrapping_add(1) as AT_RANK,
                )?;
                source_set(heap, nAtomNumberCanon1, atom_number, length)?;
                source_set(heap, nAtomNumberCanon2, atom_number, length)?;
                let avoid =
                    heap.allocate_model_storage(vec![atom_number as AT_RANK, opposite_atom])?;
                let opposite_index = usize::from(opposite_atom);
                let opposite_snapshot = heap
                    .slice(at.as_const())?
                    .get(opposite_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let mut parities_inverted = i32::from(
                    neigh_mode == NEIGH_MODE_RING as i32
                        && is_allene_chain(bond_parity)
                        && parity_calculate(bond_parity)
                        && atom_snapshot.nRingSystem == opposite_snapshot.nRingSystem
                        && i32::from(opposite_snapshot.valence) == MAX_NUM_STEREO_BONDS as i32,
                );
                if parities_inverted != 0 {
                    parities_inverted = -1;
                }

                let ret1 = CreateCheckSymmPaths(
                    heap,
                    at,
                    atom_number as AT_RANK,
                    neighbors[0],
                    atom_number as AT_RANK,
                    neighbors[1],
                    avoid,
                    nVisited1,
                    nVisited2,
                    nAtomNumberCanon1,
                    nAtomNumberCanon2,
                    path_list1,
                    path_list2,
                    symm_rank1,
                    symm_rank2,
                    nCanonRank,
                    &mut length,
                    &mut parities_inverted,
                    0,
                )?;
                let ret2 = CalculatedPathsParitiesAreIdentical(
                    heap,
                    pCG,
                    at,
                    num_atoms,
                    nSymmRank,
                    nCanonRank,
                    nAtomNumberCanon,
                    nAtomNumberCanon1,
                    nAtomNumberCanon2,
                    nVisited1,
                    nVisited2,
                    opposite_atom,
                    atom_number as AT_RANK,
                    neighbors[0],
                    neighbors[1],
                    neigh_mode,
                    parities_inverted,
                    MAP_MODE_STD as i32,
                    pCS,
                    vABParityUnknown,
                )?;

                if ret1 > 0 {
                    if ret2 > 0 {
                        if ret2 & (NOT_WELL_DEF_UNKN | NOT_WELL_DEF_UNDF) as i32 != 0 {
                            let new_parity = if ret2 & NOT_WELL_DEF_UNKN as i32 != 0 {
                                vABParityUnknown
                            } else {
                                AB_PARITY_UNDF as i32
                            };
                            let current_parity =
                                heap.slice(at.as_const())?[atom_index].stereo_bond_parity[slot];
                            if (parity_ill_defined(current_parity)
                                && parity_value(current_parity) > new_parity)
                                || parity_calculate(current_parity)
                            {
                                SetOneStereoBondIllDefParity(
                                    heap,
                                    at,
                                    atom_number,
                                    stereo_index,
                                    new_parity,
                                )?;
                                let rank_a = source_get(heap, nCanonRank, atom_number)?;
                                let rank_b =
                                    source_get(heap, nCanonRank, i32::from(opposite_atom))?;
                                let rank1 = rank_a.max(rank_b);
                                let rank2 = rank_a.min(rank_b);
                                let stereo_count = usize::try_from(pCS.nLenLinearCTStereoDble)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                let stereo_entries = heap.slice_mut(pCS.LinearCTStereoDble)?;
                                let mut found = false;
                                for entry in stereo_entries.iter_mut().take(stereo_count) {
                                    if entry.at_num1 == rank1 && entry.at_num2 == rank2 {
                                        entry.parity = new_parity as u8;
                                        found = true;
                                        break;
                                    }
                                }
                                if !found {
                                    ret = CT_STEREOCOUNT_ERR;
                                    break 'passes;
                                }
                                ret = ret.wrapping_add(1);
                            }
                        } else {
                            if RemoveOneStereoBond(heap, at, atom_number, stereo_index)? == 0 {
                                ret = CT_STEREOBOND_ERROR;
                                break 'passes;
                            }
                            stereo_index = stereo_index.wrapping_sub(1);
                            let rank_a = source_get(heap, nCanonRank, atom_number)?;
                            let rank_b = source_get(heap, nCanonRank, i32::from(opposite_atom))?;
                            let rank1 = rank_a.max(rank_b);
                            let rank2 = rank_a.min(rank_b);
                            let stereo_count = usize::try_from(pCS.nLenLinearCTStereoDble)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                            let stereo_entries = heap.slice_mut(pCS.LinearCTStereoDble)?;
                            let found = stereo_entries
                                .iter()
                                .take(stereo_count)
                                .position(|entry| entry.at_num1 == rank1 && entry.at_num2 == rank2);
                            if let Some(index) = found {
                                if index + 1 < stereo_count {
                                    for source in index + 1..stereo_count {
                                        stereo_entries[source - 1] = stereo_entries[source].clone();
                                    }
                                }
                                pCS.nLenLinearCTStereoDble =
                                    pCS.nLenLinearCTStereoDble.wrapping_sub(1);
                                ret = ret.wrapping_add(1);
                            } else {
                                ret = CT_STEREOCOUNT_ERR;
                                break 'passes;
                            }
                        }
                    } else {
                        if ret_failed == 0 {
                            ret_failed = if ret1 < 0 {
                                ret1
                            } else if ret2 < 0 {
                                ret2
                            } else {
                                0
                            };
                        }
                        if !returned_error(ret_failed) {
                            if returned_error(ret1) {
                                ret_failed = ret1;
                            } else if returned_error(ret2) {
                                ret_failed = ret2;
                            }
                        }
                    }
                }
                stereo_index = stereo_index.wrapping_add(1);
            }
            atom_number = atom_number.wrapping_add(1);
        }

        if neigh_mode == NEIGH_MODE_CHAIN as i32
            && num_equal_ring_neighbors != 0
            && !returned_error(ret_failed)
        {
            neigh_mode = NEIGH_MODE_RING as i32;
            continue;
        }
        break;
    }

    Ok(if returned_error(ret_failed) {
        ret_failed
    } else if ret_failed != 0 {
        ret_failed.wrapping_add(1).wrapping_neg()
    } else {
        ret
    })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RemoveCalculatedNonStereoCenterParities(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    num_at_tg: i32,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nTempRank: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon1: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon2: SourceMutPointer<AT_RANK>,
    nl: SourceMutPointer<NEIGH_LIST>,
    nl1: SourceMutPointer<NEIGH_LIST>,
    nl2: SourceMutPointer<NEIGH_LIST>,
    nVisited1: SourceMutPointer<AT_RANK>,
    nVisited2: SourceMutPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3356 RemoveCalculatedNonStereoCenterParities
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int RemoveCalculatedNonStereoCenterParities( CANON_GLOBALS *pCG,
                                                 sp_ATOM *at,
                                                 int num_atoms,
                                                 int num_at_tg,
                                                 AT_RANK **pRankStack1,
                                                 AT_RANK **pRankStack2,
                                                 AT_RANK *nTempRank,
                                                 NEIGH_LIST *NeighList,
                                                 AT_RANK *nCanonRank,
                                                 const AT_RANK *nSymmRank,
                                                 AT_RANK *nAtomNumberCanon,
                                                 AT_RANK *nAtomNumberCanon1,
                                                 AT_RANK *nAtomNumberCanon2,
                                                 NEIGH_LIST *nl,
                                                 NEIGH_LIST *nl1,
                                                 NEIGH_LIST *nl2,
                                                 AT_RANK *nVisited1,
                                                 AT_RANK *nVisited2,
                                                 CANON_STAT *pCS,
                                                 int vABParityUnknown )
    {
        int j, n, m, ret;

        int i, k, ret1, ret2, ret_failed = 0, mode, max_mode;
        AT_RANK nAtomRank1, neigh[MAX_NUM_STEREO_ATOM_NEIGH] = { 0 }, nAvoidCheckAtom[2], nLength; /* djb-rwth: initialisation of neigh required to avoid undefined array subscript */
        int         nNeighMode = NEIGH_MODE_CHAIN;
        int         nNumEqRingNeigh = 0, bRingNeigh, bSymmNeigh, bParitiesInverted;
        NEIGH_LIST *nl01, *nl02;
        const AT_RANK    *nSymmRank1, *nSymmRank2;

        ret = 0;

    second_pass:

        for (i = 0; i < num_atoms && !RETURNED_ERROR( ret_failed ); i++)
        {
            if (!at[i].parity || at[i].stereo_bond_neighbor[0])
            {
                continue;
            }

            if (at[i].valence > MAX_NUM_STEREO_ATOM_NEIGH)
            {
                continue; /*  error: stereo center cannot have more than 4 neighbors */ /*   <BRKPT> */
            }

            /*  at[i1] is a stereo center */
            if (!PARITY_CALCULATE( at[i].stereo_atom_parity ) && !PARITY_ILL_DEF( at[i].stereo_atom_parity ))
            {
                continue;
            }

            /* neighbors sorted according to symm. ranks (primary key) and canon. ranks (secondary key), in descending order */
            /* sorting guarantees that for two constit. equ. neighbors canon. ranks of the first is greater */
            /* !!! previously (but not anymore) the canon. rank of neigh[0] was greater than the others !!! */
            for (j = 0; j < at[i].valence; j++)
            {
                neigh[j] = nl[i][j + 1]; /*  sorting does NOT guarantee that canon. rank of neigh[0] is greater than others */
            }

            /*
             *  mode = 0 => Standard approach: switch 2 neighbors
             *         1 => Check for C2v reflection leading to parity inversion
             *         2 => Check for C2 rotation preserving parities
             *         3 => Check for S4 rotation/reflection leading to parity inversion
             */

    #if ( CHECK_C2v_S4_SYMM == 1 )
            if (nNeighMode = NEIGH_MODE_RING && at[i].valence == 4 &&
                 nSymmRank[(int) neigh[0]] == nSymmRank[(int) neigh[1]] &&
                 nSymmRank[(int) neigh[2]] == nSymmRank[(int) neigh[3]] &&
                 !at[i].bCutVertex
               )
            {
                if (nSymmRank[(int) neigh[1]] == nSymmRank[(int) neigh[2]])
                {
                    max_mode = MAP_MODE_S4;
                }
                else
                {
                    max_mode = inchi_max( MAP_MODE_C2v, MAP_MODE_C2 );
                }
            }
            else
            {
                max_mode = MAP_MODE_STD;
            }
    #else
            max_mode = MAP_MODE_STD;
    #endif

            for (j = 0; j < at[i].valence && at[i].parity && !RETURNED_ERROR( ret_failed ); j++)
            {
                for (k = j + 1; k < at[i].valence && at[i].parity && !RETURNED_ERROR( ret_failed ); k++)
                {
                    for (mode = 0; mode <= max_mode && at[i].parity && !RETURNED_ERROR(ret_failed); mode++)
                    {
                        if (nSymmRank[(int)neigh[j]] != nSymmRank[(int)neigh[k]])
                        {
                            continue; /*  the two neighbors are not constitutionally identical */
                        }
                        bRingNeigh = (at[(int)neigh[j]].nRingSystem == at[(int)neigh[k]].nRingSystem);
                        switch (nNeighMode)
                        {
                        case NEIGH_MODE_CHAIN:
                            if (bRingNeigh)
                            {
                                nNumEqRingNeigh++;
                                continue;
                            }
                            nl01 = nl;
                            nl02 = nl;
                            nSymmRank1 = nSymmRank;
                            nSymmRank2 = nSymmRank;
                            break;
                        case NEIGH_MODE_RING:
                            if (!bRingNeigh)
                                continue;
                            /*  break a tie between the two contitutionally equivalent neighbors, */
                            /*  refine the two partitions, sort neighbors lists nl1, nl2 */
                            bSymmNeigh = BreakNeighborsTie(pCG, at, num_atoms, num_at_tg, MAX_ATOMS + 1, i,
                                neigh, j, k, mode,
                                pRankStack1, pRankStack2, nTempRank, NeighList, nSymmRank, nCanonRank,
                                nl1, nl2, &pCS->lNumNeighListIter);
                            if (bSymmNeigh <= 0)
                            {
                                if (ret_failed > bSymmNeigh)
                                    ret_failed = bSymmNeigh;
                                continue;
                            }
                            nl01 = nl1;
                            nl02 = nl2;
                            nSymmRank1 = pRankStack1[0];
                            nSymmRank2 = pRankStack2[0];
                            break;
                        default:
                            return CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                        }

                        /*  initialize arrays */
                        memset(nVisited1, 0, sizeof(nVisited1[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
                        memset(nVisited2, 0, sizeof(nVisited2[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
                        memset(nAtomNumberCanon1, 0, sizeof(nAtomNumberCanon1[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
                        memset(nAtomNumberCanon2, 0, sizeof(nAtomNumberCanon2[0]) * num_atoms); /* djb-rwth: memset_s C11/Annex K variant? */
                        nLength = 1;
                        nVisited1[i] = i + 1;   /*  start atom is same */
                        nVisited2[i] = i + 1;
                        nAtomNumberCanon1[i] = nLength;
                        nAtomNumberCanon2[i] = nLength;
                        nAvoidCheckAtom[0] = i;
                        nAvoidCheckAtom[1] = MAX_ATOMS + 1;

                        bParitiesInverted = (mode == MAP_MODE_C2v || mode == MAP_MODE_S4) ? -1 : 0;

                        /*
                        if (nNeighMode==NEIGH_MODE_RING && at[i].valence==MAX_NUM_STEREO_ATOM_NEIGH) {
                            AT_RANK other_neigh[2];
                            int     n;
                            for ( m = n = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m ++ ) {
                                if ( at[i].neighbor[m] != neigh[j] && at[i].neighbor[m] != neigh[k] )
                                    other_neigh[n++] = at[i].neighbor[m];
                            }
                            if ( nSymmRank[(int)other_neigh[0]] == nSymmRank[(int)other_neigh[1]] )
                                bParitiesInverted = -1;
                        }
                        */

                        /* allow matching inverted centers only in case all equivalent neighbors in same ring system */

                        /* djb-rwth: restructuring code to avoid undefined array subscripts; removing unnecessary code; discussion required */

                        ret1 = CreateCheckSymmPaths(at, (AT_RANK)i, neigh[j], (AT_RANK)i, neigh[k], nAvoidCheckAtom,
                            nVisited1, nVisited2, nAtomNumberCanon1, nAtomNumberCanon2,
                            nl01, nl02, nSymmRank1, nSymmRank2, nCanonRank, &nLength,
                            &bParitiesInverted, mode);
                        ret2 = CalculatedPathsParitiesAreIdentical(pCG, at, num_atoms, nSymmRank,
                            nCanonRank, nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
                            nVisited1, nVisited2, (AT_RANK)MAX_ATOMS, (AT_RANK)i,
                            neigh[j], neigh[k], nNeighMode,
                            bParitiesInverted, mode, pCS, vABParityUnknown);

                        if (0 < ret1) /* djb-rwth: or is this (0 < ret1) && (0 < ret2) ? */
                        {
                            if (0 < ret2)
                            {
                                {
                                    if (ret2 & (NOT_WELL_DEF_UNKN | NOT_WELL_DEF_UNDF))
                                    {
                                        /*  possibly change the parity to unknown or undefined */
                                        int new_parity = (ret2 & NOT_WELL_DEF_UNKN) ? vABParityUnknown /*AB_PARITY_UNKN*/ : AB_PARITY_UNDF;
                                        if ((PARITY_ILL_DEF(at[i].stereo_atom_parity) &&
                                            PARITY_VAL(at[i].stereo_atom_parity) > new_parity) ||
                                            PARITY_CALCULATE(at[i].stereo_atom_parity)) /* djb-rwth: addressing LLVM warning */
                                        {
                                            /*  set new unknown or undefined parity */
                                            at[i].stereo_atom_parity = (at[i].stereo_atom_parity ^ PARITY_VAL(at[i].stereo_atom_parity)) | PARITY_VAL(new_parity);
                                            at[i].parity = PARITY_VAL(new_parity);
                                            /*  Remove from pCS */
                                            nAtomRank1 = nCanonRank[i];
                                            for (n = 0, m = pCS->nLenLinearCTStereoCarb - 1; n <= m; n++)
                                            {
                                                if (pCS->LinearCTStereoCarb[n].at_num == nAtomRank1)
                                                {
                                                    pCS->LinearCTStereoCarb[n].parity = PARITY_VAL(new_parity);
    #if ( bRELEASE_VERSION == 0 )
                                                    pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
    #endif
                                                    m = -1;
                                                    break;
                                                }
                                            }
                                            if (m >= 0)
                                            {
                                                ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                                                goto exit_function;
                                            }
                                            ret++; /*  number of removed or set unknown/undefined parities */
                                        }
                                    }
                                    else
                                    {
    #ifdef FIX_STEREOCOUNT_ERR
                                        if (at[i].stereo_atom_parity & KNOWN_PARITIES_EQL)
                                        {
                                            int jj;
                                            AT_RANK EqRank = pCS->nSymmRank[i];
                                            for (jj = 0; jj < num_atoms; jj++)
                                            {
                                                if (pCS->nSymmRank[jj] == EqRank)
                                                {
                                                    at[jj].stereo_atom_parity &= ~KNOWN_PARITIES_EQL;
                                                }
                                            }
                                        }
    #endif
                                        RemoveOneStereoCenter(at, i /* atom number*/);
                                        /*  Remove from pCS */
                                        nAtomRank1 = nCanonRank[i];
                                        for (n = 0, m = pCS->nLenLinearCTStereoCarb - 1; n <= m; n++)
                                        {
                                            if (pCS->LinearCTStereoCarb[n].at_num == nAtomRank1)
                                            {
                                                if (n < m)
                                                {
                                                    /*  remove pCS->LinearCTStereoDble[n] */
                                                    memmove(pCS->LinearCTStereoCarb + n,
                                                        pCS->LinearCTStereoCarb + n + 1,
                                                        ((long long)m - (long long)n) * sizeof(pCS->LinearCTStereoCarb[0])); /* djb-rwth: cast operators added */
                                                }
                                                pCS->nLenLinearCTStereoCarb--;
    #if ( bRELEASE_VERSION == 0 )
                                                pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
    #endif
                                                m = -1;
                                                break;
                                            }
                                        }
                                        if (m >= 0)
                                        {
                                            ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                                            goto exit_function;
                                        }
                                        ret++;  /*  number of removed or set unknown/undefined parities */
                                    }
                                }
                            }
                            else
                            {
                                if (!ret_failed)
                                {
                                    if (ret1 < 0)
                                    {
                                        ret_failed = ret1;
                                    }
                                    else
                                    {
                                        if (ret2 < 0)
                                        {
                                            ret_failed = ret2;
                                        }
                                    }
                                }
                                if (!RETURNED_ERROR(ret_failed))
                                {
                                    if (RETURNED_ERROR(ret1))
                                    {
                                        ret_failed = ret1;
                                    }
                                    else
                                    {
                                        if (RETURNED_ERROR(ret2))
                                        {
                                            ret_failed = ret2;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (nNeighMode == NEIGH_MODE_CHAIN && nNumEqRingNeigh && !RETURNED_ERROR( ret_failed ))
        {
            nNeighMode = NEIGH_MODE_RING;
            goto second_pass;
        }

    exit_function:

        return RETURNED_ERROR( ret_failed ) ? ret_failed : ret_failed ? -( ret + 1 ) : ret;
    }


    /****************************************************************************/
    */
    // END INCHI C FUNCTION: RemoveCalculatedNonStereoCenterParities
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION
    // INCHI✔️❌: #define REMOVE_CALC_NONSTEREO 1
    // INCHI✔️❌: #define CHECK_C2v_S4_SYMM 0
    // INCHI✔️❌: #define FIX_STEREOCOUNT_ERR 1
    // INCHI✔️❌: #define bRELEASE_VERSION 1
    // END INCHI ACTIVE MACRO CONFIGURATION

    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let returned_error = |value: i32| CT_ERR_MIN <= value && value <= CT_ERR_MAX;
    let parity_value = |value: i8| i32::from(value) & BITS_PARITY as i32;
    let parity_calculate = |value: i8| parity_value(value) == AB_PARITY_CALC as i32;
    let parity_ill_defined = |value: i8| {
        (AB_MIN_ILL_DEFINED_PARITY as i32..=AB_MAX_ILL_DEFINED_PARITY as i32)
            .contains(&parity_value(value))
    };

    let mut ret = 0_i32;
    let mut ret_failed = 0_i32;
    let mut neigh_mode = NEIGH_MODE_CHAIN as i32;
    let mut num_equal_ring_neighbors = 0_i32;

    'passes: loop {
        let mut atom_number = 0_i32;
        while atom_number < num_atoms && !returned_error(ret_failed) {
            let atom_index =
                usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom_snapshot = heap
                .slice(at.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if atom_snapshot.parity == 0 || atom_snapshot.stereo_bond_neighbor[0] != 0 {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }
            let valence = i32::from(atom_snapshot.valence);
            if valence > MAX_NUM_STEREO_ATOM_NEIGH as i32 {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }
            if !parity_calculate(atom_snapshot.stereo_atom_parity)
                && !parity_ill_defined(atom_snapshot.stereo_atom_parity)
            {
                atom_number = atom_number.wrapping_add(1);
                continue;
            }

            let mut neighbors = [0_u16; MAX_NUM_STEREO_ATOM_NEIGH as usize];
            let neighbor_count = if valence > 0 {
                usize::try_from(valence).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            } else {
                0
            };
            if neighbor_count != 0 {
                let row = source_get(heap, nl, atom_number)?;
                for position in 0..neighbor_count {
                    neighbors[position] = source_get(
                        heap,
                        row,
                        i32::try_from(position + 1)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                    )?;
                }
            }

            let max_mode = MAP_MODE_STD as i32;
            let mut j = 0_i32;
            while j < valence
                && heap.slice(at.as_const())?[atom_index].parity != 0
                && !returned_error(ret_failed)
            {
                let mut k = j.wrapping_add(1);
                while k < valence
                    && heap.slice(at.as_const())?[atom_index].parity != 0
                    && !returned_error(ret_failed)
                {
                    let mut mode = 0_i32;
                    while mode <= max_mode
                        && heap.slice(at.as_const())?[atom_index].parity != 0
                        && !returned_error(ret_failed)
                    {
                        let first = neighbors[usize::try_from(j)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                        let second = neighbors[usize::try_from(k)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                        if source_get(heap, nSymmRank, i32::from(first))?
                            != source_get(heap, nSymmRank, i32::from(second))?
                        {
                            mode = mode.wrapping_add(1);
                            continue;
                        }

                        let ring_neighbors = {
                            let atoms = heap.slice(at.as_const())?;
                            atoms
                                .get(usize::from(first))
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nRingSystem
                                == atoms
                                    .get(usize::from(second))
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .nRingSystem
                        };

                        let (path_list1, path_list2, symm_rank1, symm_rank2) = if neigh_mode
                            == NEIGH_MODE_CHAIN as i32
                        {
                            if ring_neighbors {
                                num_equal_ring_neighbors = num_equal_ring_neighbors.wrapping_add(1);
                                mode = mode.wrapping_add(1);
                                continue;
                            }
                            (nl, nl, nSymmRank, nSymmRank)
                        } else if neigh_mode == NEIGH_MODE_RING as i32 {
                            if !ring_neighbors {
                                mode = mode.wrapping_add(1);
                                continue;
                            }
                            let neighbor_storage =
                                heap.allocate_model_storage(neighbors.to_vec())?;
                            let symmetric = BreakNeighborsTie(
                                heap,
                                pCG,
                                at,
                                num_atoms,
                                num_at_tg,
                                i32::try_from(MAX_ATOMS.wrapping_add(1))
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                                atom_number,
                                neighbor_storage,
                                j,
                                k,
                                mode,
                                pRankStack1,
                                pRankStack2,
                                nTempRank,
                                NeighList,
                                nSymmRank,
                                nCanonRank,
                                nl1,
                                nl2,
                                &mut pCS.lNumNeighListIter,
                            );
                            let cleanup = heap.free(neighbor_storage);
                            let symmetric = match (symmetric, cleanup) {
                                (Err(error), _) | (Ok(_), Err(error)) => return Err(error),
                                (Ok(value), Ok(())) => value,
                            };
                            if symmetric <= 0 {
                                if ret_failed > symmetric {
                                    ret_failed = symmetric;
                                }
                                mode = mode.wrapping_add(1);
                                continue;
                            }
                            (
                                nl1,
                                nl2,
                                source_get(heap, pRankStack1, 0)?,
                                source_get(heap, pRankStack2, 0)?,
                            )
                        } else {
                            return Ok(CT_STEREOCOUNT_ERR);
                        };

                        heap.slice_mut(nVisited1)?
                            .get_mut(..count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .fill(0);
                        heap.slice_mut(nVisited2)?
                            .get_mut(..count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .fill(0);
                        heap.slice_mut(nAtomNumberCanon1)?
                            .get_mut(..count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .fill(0);
                        heap.slice_mut(nAtomNumberCanon2)?
                            .get_mut(..count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .fill(0);

                        let mut length = 1_u16;
                        let atom_rank = AT_RANK::try_from(atom_number)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                        source_set(heap, nVisited1, atom_number, atom_rank.wrapping_add(1))?;
                        source_set(heap, nVisited2, atom_number, atom_rank.wrapping_add(1))?;
                        source_set(heap, nAtomNumberCanon1, atom_number, length)?;
                        source_set(heap, nAtomNumberCanon2, atom_number, length)?;
                        let avoid = heap.allocate_model_storage(vec![
                            atom_rank,
                            AT_RANK::try_from(MAX_ATOMS.wrapping_add(1))
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                        ])?;
                        let mut parities_inverted =
                            if mode == MAP_MODE_C2v as i32 || mode == MAP_MODE_S4 as i32 {
                                -1
                            } else {
                                0
                            };

                        let checks = (|| -> Result<(i32, i32), SourceHeapError> {
                            let ret1 = CreateCheckSymmPaths(
                                heap,
                                at,
                                atom_rank,
                                first,
                                atom_rank,
                                second,
                                avoid,
                                nVisited1,
                                nVisited2,
                                nAtomNumberCanon1,
                                nAtomNumberCanon2,
                                path_list1,
                                path_list2,
                                symm_rank1,
                                symm_rank2,
                                nCanonRank,
                                &mut length,
                                &mut parities_inverted,
                                mode,
                            )?;
                            let ret2 = CalculatedPathsParitiesAreIdentical(
                                heap,
                                pCG,
                                at,
                                num_atoms,
                                nSymmRank,
                                nCanonRank,
                                nAtomNumberCanon,
                                nAtomNumberCanon1,
                                nAtomNumberCanon2,
                                nVisited1,
                                nVisited2,
                                MAX_ATOMS as AT_RANK,
                                atom_rank,
                                first,
                                second,
                                neigh_mode,
                                parities_inverted,
                                mode,
                                pCS,
                                vABParityUnknown,
                            )?;
                            Ok((ret1, ret2))
                        })();
                        let cleanup = heap.free(avoid);
                        let (ret1, ret2) = match (checks, cleanup) {
                            (Err(error), _) | (Ok(_), Err(error)) => return Err(error),
                            (Ok(values), Ok(())) => values,
                        };

                        if ret1 > 0 {
                            if ret2 > 0 {
                                if ret2 & (NOT_WELL_DEF_UNKN | NOT_WELL_DEF_UNDF) as i32 != 0 {
                                    let new_parity = if ret2 & NOT_WELL_DEF_UNKN as i32 != 0 {
                                        vABParityUnknown
                                    } else {
                                        AB_PARITY_UNDF as i32
                                    };
                                    let current =
                                        heap.slice(at.as_const())?[atom_index].stereo_atom_parity;
                                    if (parity_ill_defined(current)
                                        && parity_value(current) > new_parity)
                                        || parity_calculate(current)
                                    {
                                        let atom = &mut heap.slice_mut(at)?[atom_index];
                                        atom.stereo_atom_parity =
                                            (i32::from(atom.stereo_atom_parity)
                                                ^ parity_value(atom.stereo_atom_parity)
                                                | (new_parity & BITS_PARITY as i32))
                                                as i8;
                                        atom.parity = (new_parity & BITS_PARITY as i32) as i8;

                                        let canonical_rank =
                                            source_get(heap, nCanonRank, atom_number)?;
                                        let stereo_count = usize::try_from(
                                            pCS.nLenLinearCTStereoCarb,
                                        )
                                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                        let entries = heap.slice_mut(pCS.LinearCTStereoCarb)?;
                                        let mut found = false;
                                        for entry in entries.iter_mut().take(stereo_count) {
                                            if entry.at_num == canonical_rank {
                                                entry.parity =
                                                    (new_parity & BITS_PARITY as i32) as u8;
                                                found = true;
                                                break;
                                            }
                                        }
                                        if !found {
                                            ret = CT_STEREOCOUNT_ERR;
                                            break 'passes;
                                        }
                                        ret = ret.wrapping_add(1);
                                    }
                                } else {
                                    let current =
                                        heap.slice(at.as_const())?[atom_index].stereo_atom_parity;
                                    if current & KNOWN_PARITIES_EQL as i8 != 0 {
                                        let equivalent_rank =
                                            source_get(heap, pCS.nSymmRank, atom_number)?;
                                        let mut equivalent = 0_i32;
                                        while equivalent < num_atoms {
                                            if source_get(heap, pCS.nSymmRank, equivalent)?
                                                == equivalent_rank
                                            {
                                                let index =
                                                    usize::try_from(equivalent).map_err(|_| {
                                                        SourceHeapError::PointerOutOfBounds
                                                    })?;
                                                heap.slice_mut(at)?[index].stereo_atom_parity &=
                                                    !(KNOWN_PARITIES_EQL as i8);
                                            }
                                            equivalent = equivalent.wrapping_add(1);
                                        }
                                    }
                                    RemoveOneStereoCenter(heap, at, atom_number)?;

                                    let canonical_rank = source_get(heap, nCanonRank, atom_number)?;
                                    let stereo_count = usize::try_from(pCS.nLenLinearCTStereoCarb)
                                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                    let entries = heap.slice_mut(pCS.LinearCTStereoCarb)?;
                                    let found = entries
                                        .iter()
                                        .take(stereo_count)
                                        .position(|entry| entry.at_num == canonical_rank);
                                    if let Some(index) = found {
                                        for source in index + 1..stereo_count {
                                            entries[source - 1] = entries[source].clone();
                                        }
                                        pCS.nLenLinearCTStereoCarb =
                                            pCS.nLenLinearCTStereoCarb.wrapping_sub(1);
                                        ret = ret.wrapping_add(1);
                                    } else {
                                        ret = CT_STEREOCOUNT_ERR;
                                        break 'passes;
                                    }
                                }
                            } else {
                                if ret_failed == 0 {
                                    if ret1 < 0 {
                                        ret_failed = ret1;
                                    } else if ret2 < 0 {
                                        ret_failed = ret2;
                                    }
                                }
                                if !returned_error(ret_failed) {
                                    if returned_error(ret1) {
                                        ret_failed = ret1;
                                    } else if returned_error(ret2) {
                                        ret_failed = ret2;
                                    }
                                }
                            }
                        }

                        mode = mode.wrapping_add(1);
                    }
                    k = k.wrapping_add(1);
                }
                j = j.wrapping_add(1);
            }
            atom_number = atom_number.wrapping_add(1);
        }

        if neigh_mode == NEIGH_MODE_CHAIN as i32
            && num_equal_ring_neighbors != 0
            && !returned_error(ret_failed)
        {
            neigh_mode = NEIGH_MODE_RING as i32;
            continue;
        }
        break;
    }

    Ok(if returned_error(ret_failed) {
        ret_failed
    } else if ret_failed != 0 {
        ret.wrapping_add(1).wrapping_neg()
    } else {
        ret
    })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn RemoveCalculatedNonStereo(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    num_at_tg: i32,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nTempRank: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3672 RemoveCalculatedNonStereo
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int RemoveCalculatedNonStereo( CANON_GLOBALS *pCG,
                                   sp_ATOM *at,
                                   int num_atoms,
                                   int num_at_tg,
                                   AT_RANK **pRankStack1,
                                   AT_RANK **pRankStack2,
                                   AT_RANK *nTempRank,
                                   NEIGH_LIST *NeighList,
                                   const AT_RANK *nSymmRank,
                                   AT_RANK *nCanonRank,
                                   AT_RANK *nAtomNumberCanon,
                                   CANON_STAT *pCS,
                                   int vABParityUnknown )
    {
        NEIGH_LIST *nl = NULL, *nl1 = NULL, *nl2 = NULL;
        AT_RANK    *nVisited1 = NULL, *nVisited2 = NULL, *nAtomNumberCanon1 = NULL, *nAtomNumberCanon2 = NULL;
        int        nNumRemoved = 0, nTotRemoved = 0, ret = 0, ret1 = 0, ret2 = 0;

        if (!AllocateForNonStereoRemoval( at, num_atoms, nSymmRank, nCanonRank,
            &nAtomNumberCanon1, &nAtomNumberCanon2,
            &nl, &nl1, &nl2, &nVisited1, &nVisited2 ))
        {
            return CT_OUT_OF_RAM;  /*   <BRKPT> */
        }

        do
        {
            nNumRemoved = 0;
            /*  bonds */
            ret = RemoveCalculatedNonStereoBondParities( pCG, at, num_atoms, num_at_tg,
                                                  pRankStack1, pRankStack2, nTempRank, NeighList,
                                                  nCanonRank, nSymmRank,
                                                  nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
                                                  nl, nl1, nl2, nVisited1, nVisited2, pCS,
                                                  vABParityUnknown );
            if (RETURNED_ERROR( ret ))
            {
                goto exit_function;
            }
            if (ret < 0)
            {
                if (ret < ret1)
                {  /*   <BRKPT> */
                    ret1 = ret;
                }
                ret = -( ret + 1 ); /*  number of removed */
            }
            nNumRemoved += ret;

            /*  centers */
            ret = RemoveCalculatedNonStereoCenterParities( pCG, at, num_atoms, num_at_tg,
                                                  pRankStack1, pRankStack2, nTempRank, NeighList,
                                                  nCanonRank, nSymmRank,
                                                  nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
                                                  nl, nl1, nl2, nVisited1, nVisited2, pCS,
                                                  vABParityUnknown );
            if (RETURNED_ERROR( ret ))
            {
                goto exit_function;
            }
            if (ret < 0)
            {
                if (ret < ret2)
                {  /*   <BRKPT> */
                    ret2 = ret;
                }
                ret = -( ret + 1 ); /*  number of removed */
            }
            nNumRemoved += ret;

            nTotRemoved += nNumRemoved;
        }
        while (nNumRemoved);

        if (!RETURNED_ERROR( ret1 ) && !RETURNED_ERROR( ret2 ))
        {
            ret = inchi_min( ret1, ret2 );
            ret = ( ret >= 0 ) ? nTotRemoved : -( 1 + nTotRemoved );
        }

    exit_function:

        DeAllocateForNonStereoRemoval( &nAtomNumberCanon1, &nAtomNumberCanon2, &nl, &nl1, &nl2, &nVisited1, &nVisited2 );

        return ret;
    }

    #endif /* } REMOVE_CALC_NONSTEREO */
    */
    // END INCHI C FUNCTION: RemoveCalculatedNonStereo
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION
    // INCHI✔️❌: #define REMOVE_CALC_NONSTEREO 1
    // INCHI✔️❌: #define CHECK_C2v_S4_SYMM 0
    // INCHI✔️❌: #define FIX_STEREOCOUNT_ERR 1
    // INCHI✔️❌: #define bRELEASE_VERSION 1
    // END INCHI ACTIVE MACRO CONFIGURATION

    let mut atom_number_canon1 = SourceMutPointer::null();
    let mut atom_number_canon2 = SourceMutPointer::null();
    let mut nl = SourceMutPointer::null();
    let mut nl1 = SourceMutPointer::null();
    let mut nl2 = SourceMutPointer::null();
    let mut visited1 = SourceMutPointer::null();
    let mut visited2 = SourceMutPointer::null();

    if AllocateForNonStereoRemoval(
        heap,
        at,
        num_atoms,
        nSymmRank,
        nCanonRank,
        &mut atom_number_canon1,
        &mut atom_number_canon2,
        &mut nl,
        &mut nl1,
        &mut nl2,
        &mut visited1,
        &mut visited2,
    )? == 0
    {
        return Ok(CT_OUT_OF_RAM);
    }

    let returned_error = |value: i32| CT_ERR_MIN <= value && value <= CT_ERR_MAX;
    let operation = (|| -> Result<i32, SourceHeapError> {
        let mut total_removed = 0_i32;
        let mut bond_failure = 0_i32;
        let mut center_failure = 0_i32;

        let ret = loop {
            let mut removed = 0_i32;
            let mut ret = RemoveCalculatedNonStereoBondParities(
                heap,
                pCG,
                at,
                num_atoms,
                num_at_tg,
                pRankStack1,
                pRankStack2,
                nTempRank,
                NeighList,
                nCanonRank,
                nSymmRank,
                nAtomNumberCanon,
                atom_number_canon1,
                atom_number_canon2,
                nl,
                nl1,
                nl2,
                visited1,
                visited2,
                pCS,
                vABParityUnknown,
            )?;
            if returned_error(ret) {
                break ret;
            }
            if ret < 0 {
                if ret < bond_failure {
                    bond_failure = ret;
                }
                ret = ret.wrapping_add(1).wrapping_neg();
            }
            removed = removed.wrapping_add(ret);

            ret = RemoveCalculatedNonStereoCenterParities(
                heap,
                pCG,
                at,
                num_atoms,
                num_at_tg,
                pRankStack1,
                pRankStack2,
                nTempRank,
                NeighList,
                nCanonRank,
                nSymmRank,
                nAtomNumberCanon,
                atom_number_canon1,
                atom_number_canon2,
                nl,
                nl1,
                nl2,
                visited1,
                visited2,
                pCS,
                vABParityUnknown,
            )?;
            if returned_error(ret) {
                break ret;
            }
            if ret < 0 {
                if ret < center_failure {
                    center_failure = ret;
                }
                ret = ret.wrapping_add(1).wrapping_neg();
            }
            removed = removed.wrapping_add(ret);
            total_removed = total_removed.wrapping_add(removed);

            if removed == 0 {
                let failure = bond_failure.min(center_failure);
                break if failure >= 0 {
                    total_removed
                } else {
                    total_removed.wrapping_add(1).wrapping_neg()
                };
            }
        };
        Ok(ret)
    })();

    let cleanup = DeAllocateForNonStereoRemoval(
        heap,
        &mut atom_number_canon1,
        &mut atom_number_canon2,
        &mut nl,
        &mut nl1,
        &mut nl2,
        &mut visited1,
        &mut visited2,
    );
    match (operation, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(ret), Ok(())) => Ok(ret),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichimap2__removecalculatednonstereobondparities__line_3099() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        let mut globals = CANON_GLOBALS::default();
        let mut stats = CANON_STAT::default();
        assert_eq!(
            RemoveCalculatedNonStereoBondParities(
                &mut heap,
                &mut globals,
                atoms,
                1,
                1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(0)
        );

        let mut skip_heap = SourceHeap::default();
        let mut skip_atom = sp_ATOM::default();
        skip_atom.valence = 3;
        skip_atom.stereo_bond_neighbor[0] = 1;
        skip_atom.stereo_bond_parity[0] = crate::source_types::AB_PARITY_ODD as i8;
        let skip_atoms = skip_heap.allocate_model_storage(vec![skip_atom]).unwrap();
        let mut skip_globals = CANON_GLOBALS::default();
        let mut skip_stats = CANON_STAT::default();
        assert_eq!(
            RemoveCalculatedNonStereoBondParities(
                &mut skip_heap,
                &mut skip_globals,
                skip_atoms,
                1,
                1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut skip_stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(0)
        );

        let mut error_heap = SourceHeap::default();
        let mut error_atoms = vec![sp_ATOM::default(); 2];
        error_atoms[0].valence = 3;
        error_atoms[0].neighbor[0] = 1;
        error_atoms[0].stereo_bond_neighbor[0] = 2;
        error_atoms[0].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
        let error_atoms = error_heap.allocate_model_storage(error_atoms).unwrap();
        let error_row = error_heap
            .allocate_model_storage(vec![3_u16, 1, 1, 1])
            .unwrap();
        let error_nl = error_heap
            .allocate_model_storage(vec![error_row, error_row])
            .unwrap();
        let mut error_globals = CANON_GLOBALS::default();
        let mut error_stats = CANON_STAT::default();
        assert_eq!(
            RemoveCalculatedNonStereoBondParities(
                &mut error_heap,
                &mut error_globals,
                error_atoms,
                2,
                2,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                error_nl,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut error_stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(CT_STEREOBOND_ERROR)
        );

        let mut remove_heap = SourceHeap::default();
        let mut atoms = vec![sp_ATOM::default(); 4];
        atoms[0].valence = 3;
        atoms[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
        atoms[0].stereo_bond_neighbor[0] = 2;
        atoms[0].stereo_bond_ord[0] = 0;
        atoms[0].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
        atoms[1].valence = 1;
        atoms[1].neighbor[0] = 0;
        atoms[1].stereo_bond_neighbor[0] = 1;
        atoms[1].stereo_bond_ord[0] = 0;
        atoms[1].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
        atoms[2].valence = 1;
        atoms[2].neighbor[0] = 0;
        atoms[2].nRingSystem = 1;
        atoms[2].stereo_atom_parity = crate::source_types::AB_PARITY_ODD as i8;
        atoms[3].valence = 1;
        atoms[3].neighbor[0] = 0;
        atoms[3].nRingSystem = 2;
        atoms[3].stereo_atom_parity = crate::source_types::AB_PARITY_ODD as i8;
        let atoms = remove_heap.allocate_model_storage(atoms).unwrap();
        let row0 = remove_heap
            .allocate_model_storage(vec![3_u16, 3, 2, 1])
            .unwrap();
        let row1 = remove_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let row2 = remove_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let row3 = remove_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let nl = remove_heap
            .allocate_model_storage(vec![row0, row1, row2, row3])
            .unwrap();
        let canon = remove_heap
            .allocate_model_storage(vec![1_u16, 2, 3, 4])
            .unwrap();
        let symm = remove_heap
            .allocate_model_storage(vec![1_u16, 2, 3, 3])
            .unwrap();
        let atom_by_canon = remove_heap
            .allocate_model_storage(vec![0_u16, 1, 2, 3])
            .unwrap();
        let atom_by_canon1 = remove_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let atom_by_canon2 = remove_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let visited1 = remove_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let visited2 = remove_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let stereo_entries = remove_heap
            .allocate_model_storage(vec![crate::source_types::AT_STEREO_DBLE {
                at_num1: 2,
                at_num2: 1,
                parity: AB_PARITY_CALC as u8,
            }])
            .unwrap();
        let mut stats = CANON_STAT::default();
        stats.LinearCTStereoDble = stereo_entries;
        stats.nLenLinearCTStereoDble = 1;
        let mut globals = CANON_GLOBALS::default();
        assert_eq!(
            RemoveCalculatedNonStereoBondParities(
                &mut remove_heap,
                &mut globals,
                atoms,
                4,
                4,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                canon,
                symm,
                atom_by_canon,
                atom_by_canon1,
                atom_by_canon2,
                nl,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                visited1,
                visited2,
                &mut stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(1)
        );
        assert_eq!(stats.nLenLinearCTStereoDble, 0);
        assert_eq!(
            remove_heap.slice(atoms.as_const()).unwrap()[0].stereo_bond_neighbor[0],
            0
        );
        assert_eq!(
            remove_heap.slice(atoms.as_const()).unwrap()[1].stereo_bond_neighbor[0],
            0
        );
    }

    #[test]
    fn source_port__ichimap2__calculatedpathsparitiesareidentical__line_2399() {
        let sentinel = (MAX_ATOMS + 1) as AT_RANK;

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default(); 3])
            .unwrap();
        let symm = heap.allocate_model_storage(vec![1_u16, 1, 1]).unwrap();
        let canon = heap.allocate_model_storage(vec![2_u16, 3, 1]).unwrap();
        let atom_by_canon = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let atom_by_canon1 = heap.allocate_model_storage(vec![99_u16; 3]).unwrap();
        let atom_by_canon2 = heap.allocate_model_storage(vec![88_u16; 3]).unwrap();
        let visited1 = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let visited2 = heap.allocate_model_storage(vec![3_u16, 1, 2]).unwrap();
        let mut globals = CANON_GLOBALS::default();
        let mut stats = CANON_STAT::default();
        assert_eq!(
            CalculatedPathsParitiesAreIdentical(
                &mut heap,
                &mut globals,
                atoms,
                3,
                symm,
                canon,
                atom_by_canon,
                atom_by_canon1,
                atom_by_canon2,
                visited1,
                visited2,
                sentinel,
                0,
                1,
                2,
                NEIGH_MODE_RING as i32,
                0,
                MAP_MODE_STD as i32,
                &mut stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(COMP_STEREO_SUCCESS as i32)
        );
        assert_eq!(heap.slice(visited1.as_const()).unwrap(), &[2, 3, 1]);
        assert_eq!(heap.slice(visited2.as_const()).unwrap(), &[1, 2, 3]);
        assert_eq!(heap.slice(atom_by_canon1.as_const()).unwrap(), &[2, 0, 1]);
        assert_eq!(heap.slice(atom_by_canon2.as_const()).unwrap(), &[0, 1, 2]);

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default(); 3])
            .unwrap();
        let symm = heap.allocate_model_storage(vec![1_u16; 3]).unwrap();
        let canon = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let atom_by_canon = heap.allocate_model_storage(vec![0_u16, 1, 2]).unwrap();
        let atom_by_canon1 = heap.allocate_model_storage(vec![77_u16; 3]).unwrap();
        let atom_by_canon2 = heap.allocate_model_storage(vec![66_u16; 3]).unwrap();
        let visited1 = heap.allocate_model_storage(vec![1_u16, 0, 3]).unwrap();
        let visited2 = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let mut globals = CANON_GLOBALS::default();
        let mut stats = CANON_STAT::default();
        assert_eq!(
            CalculatedPathsParitiesAreIdentical(
                &mut heap,
                &mut globals,
                atoms,
                3,
                symm,
                canon,
                atom_by_canon,
                atom_by_canon1,
                atom_by_canon2,
                visited1,
                visited2,
                sentinel,
                0,
                1,
                2,
                NEIGH_MODE_RING as i32,
                1,
                MAP_MODE_C2v as i32,
                &mut stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(visited1.as_const()).unwrap(), &[1, 0, 3]);
        assert_eq!(
            heap.slice(atom_by_canon1.as_const()).unwrap(),
            &[sentinel; 3]
        );
        assert_eq!(
            heap.slice(atom_by_canon2.as_const()).unwrap(),
            &[sentinel; 3]
        );

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default(); 3])
            .unwrap();
        let symm = heap.allocate_model_storage(vec![1_u16; 3]).unwrap();
        let canon = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let atom_by_canon = heap.allocate_model_storage(vec![0_u16, 1, 2]).unwrap();
        let atom_by_canon1 = heap.allocate_model_storage(vec![77_u16; 3]).unwrap();
        let atom_by_canon2 = heap.allocate_model_storage(vec![66_u16; 3]).unwrap();
        let visited1 = heap.allocate_model_storage(vec![1_u16, 0, 3]).unwrap();
        let visited2 = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let mut globals = CANON_GLOBALS::default();
        let mut stats = CANON_STAT::default();
        assert_eq!(
            CalculatedPathsParitiesAreIdentical(
                &mut heap,
                &mut globals,
                atoms,
                3,
                symm,
                canon,
                atom_by_canon,
                atom_by_canon1,
                atom_by_canon2,
                visited1,
                visited2,
                sentinel,
                0,
                1,
                2,
                NEIGH_MODE_RING as i32,
                -1,
                MAP_MODE_STD as i32,
                &mut stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(COMP_STEREO_SUCCESS as i32)
        );
        assert_eq!(heap.slice(visited1.as_const()).unwrap(), &[1, 0, 3]);
        assert_eq!(
            heap.slice(atom_by_canon1.as_const()).unwrap(),
            &[0, sentinel, 2]
        );

        let mut heap = SourceHeap::default();
        let mut atom_values = vec![sp_ATOM::default(); 4];
        atom_values[0].valence = 3;
        atom_values[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
        let atoms = heap.allocate_model_storage(atom_values).unwrap();
        let symm = heap.allocate_model_storage(vec![1_u16; 4]).unwrap();
        let canon = heap.allocate_model_storage(vec![1_u16, 2, 3, 4]).unwrap();
        let atom_by_canon = heap.allocate_model_storage(vec![0_u16, 1, 2, 3]).unwrap();
        let atom_by_canon1 = heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let atom_by_canon2 = heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let visited1 = heap.allocate_model_storage(vec![1_u16, 2, 3, 4]).unwrap();
        let visited2 = heap.allocate_model_storage(vec![1_u16, 2, 3, 4]).unwrap();
        let mut globals = CANON_GLOBALS::default();
        let mut stats = CANON_STAT::default();
        assert_eq!(
            CalculatedPathsParitiesAreIdentical(
                &mut heap,
                &mut globals,
                atoms,
                4,
                symm,
                canon,
                atom_by_canon,
                atom_by_canon1,
                atom_by_canon2,
                visited1,
                visited2,
                sentinel,
                0,
                sentinel,
                sentinel,
                NEIGH_MODE_RING as i32,
                0,
                MAP_MODE_STD as i32,
                &mut stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(CT_STEREOCOUNT_ERR)
        );

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default(); 2])
            .unwrap();
        let symm = heap.allocate_model_storage(vec![1_u16; 2]).unwrap();
        let canon = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let atom_by_canon = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let atom_by_canon1 = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let atom_by_canon2 = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let visited1 = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let visited2 = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let mut globals = CANON_GLOBALS::default();
        let mut stats = CANON_STAT::default();
        assert_eq!(
            CalculatedPathsParitiesAreIdentical(
                &mut heap,
                &mut globals,
                atoms,
                2,
                symm,
                canon,
                atom_by_canon,
                atom_by_canon1,
                atom_by_canon2,
                visited1,
                visited2,
                1,
                0,
                0,
                0,
                NEIGH_MODE_RING as i32,
                0,
                MAP_MODE_STD as i32,
                &mut stats,
                crate::source_types::AB_PARITY_UNKN as i32,
            ),
            Ok(-1)
        );
    }

    #[allow(clippy::type_complexity)]
    fn tie_fixture(
        heap: &mut SourceHeap,
        ranks1: Vec<AT_RANK>,
        order1: Vec<AT_RANK>,
        ranks2: Vec<AT_RANK>,
        order2: Vec<AT_RANK>,
        extra1: [SourceMutPointer<AT_RANK>; 2],
        extra2: [SourceMutPointer<AT_RANK>; 2],
    ) -> (ppAT_RANK, ppAT_RANK) {
        let ranks1 = heap.allocate_model_storage(ranks1).unwrap();
        let order1 = heap.allocate_model_storage(order1).unwrap();
        let ranks2 = heap.allocate_model_storage(ranks2).unwrap();
        let order2 = heap.allocate_model_storage(order2).unwrap();
        let stack1 = heap
            .allocate_model_storage(vec![ranks1, order1, extra1[0], extra1[1]])
            .unwrap();
        let stack2 = heap
            .allocate_model_storage(vec![ranks2, order2, extra2[0], extra2[1]])
            .unwrap();
        (stack1, stack2)
    }

    #[test]
    fn source_port__ichimap2__numberofties__line_680() {
        let mut heap = SourceHeap::default();
        let mapped = heap.allocate_model_storage(vec![5_u16, 9, 9]).unwrap();
        let unmapped = heap.allocate_model_storage(vec![0_u16, 9, 9]).unwrap();
        let copy_rank = heap.allocate_model_storage(vec![99_u16; 3]).unwrap();
        let copy_order = heap.allocate_model_storage(vec![88_u16; 3]).unwrap();
        let (stack1, stack2) = tie_fixture(
            &mut heap,
            vec![2, 2, 3],
            vec![0, 1, 2],
            vec![2, 2, 3],
            vec![1, 0, 2],
            [mapped, unmapped],
            [copy_rank, copy_order],
        );
        let mut new_rank = 77;
        let mut add_stack = 77;
        let mut mapped_count = 77;
        assert_eq!(
            NumberOfTies(
                &mut heap,
                stack1,
                stack2,
                6,
                0,
                1,
                &mut new_rank,
                &mut add_stack,
                &mut mapped_count
            ),
            Ok(2)
        );
        assert_eq!((new_rank, add_stack, mapped_count), (1, 2, 1));
        assert_eq!(heap.slice(mapped.as_const()).unwrap(), &[5, 9, 9]);
        assert_eq!(heap.slice(unmapped.as_const()).unwrap(), &[0, 9, 9]);
        assert_eq!(heap.slice(copy_rank.as_const()).unwrap(), &[2, 2, 3]);
        assert_eq!(heap.slice(copy_order.as_const()).unwrap(), &[1, 0, 2]);

        let mut no_tie_heap = SourceHeap::default();
        let (no_tie1, no_tie2) = tie_fixture(
            &mut no_tie_heap,
            vec![1, 2],
            vec![0, 1],
            vec![1, 2],
            vec![0, 1],
            [SourceMutPointer::null(); 2],
            [SourceMutPointer::null(); 2],
        );
        new_rank = 9;
        add_stack = 9;
        mapped_count = 9;
        assert_eq!(
            NumberOfTies(
                &mut no_tie_heap,
                no_tie1,
                no_tie2,
                4,
                0,
                0,
                &mut new_rank,
                &mut add_stack,
                &mut mapped_count
            ),
            Ok(1)
        );
        assert_eq!((new_rank, add_stack, mapped_count), (0, 0, 0));

        let mut mismatch_heap = SourceHeap::default();
        let (mismatch1, mismatch2) = tie_fixture(
            &mut mismatch_heap,
            vec![2, 2],
            vec![0, 1],
            vec![1, 2],
            vec![0, 1],
            [SourceMutPointer::null(); 2],
            [SourceMutPointer::null(); 2],
        );
        assert_eq!(
            NumberOfTies(
                &mut mismatch_heap,
                mismatch1,
                mismatch2,
                4,
                0,
                0,
                &mut new_rank,
                &mut add_stack,
                &mut mapped_count
            ),
            Ok(CT_MAPCOUNT_ERR)
        );
        assert_eq!((new_rank, add_stack, mapped_count), (0, 0, 0));

        let mut count_mismatch_heap = SourceHeap::default();
        let (count1, count2) = tie_fixture(
            &mut count_mismatch_heap,
            vec![2, 2],
            vec![0, 1],
            vec![2, 1],
            vec![1, 0],
            [SourceMutPointer::null(); 2],
            [SourceMutPointer::null(); 2],
        );
        assert_eq!(
            NumberOfTies(
                &mut count_mismatch_heap,
                count1,
                count2,
                4,
                0,
                0,
                &mut new_rank,
                &mut add_stack,
                &mut mapped_count
            ),
            Ok(CT_MAPCOUNT_ERR)
        );

        let mut allocation_heap = SourceHeap::default();
        let (allocation1, allocation2) = tie_fixture(
            &mut allocation_heap,
            vec![2, 2],
            vec![0, 1],
            vec![2, 2],
            vec![0, 1],
            [SourceMutPointer::null(); 2],
            [SourceMutPointer::null(); 2],
        );
        allocation_heap.fail_after_allocations(1);
        assert_eq!(
            NumberOfTies(
                &mut allocation_heap,
                allocation1,
                allocation2,
                4,
                0,
                0,
                &mut new_rank,
                &mut add_stack,
                &mut mapped_count
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert_eq!((new_rank, add_stack, mapped_count), (1, 0, 0));
        let allocation_slots = allocation_heap.slice(allocation1.as_const()).unwrap();
        assert!(!allocation_slots[2].is_null());
        assert!(allocation_slots[3].is_null());

        assert_eq!(
            NumberOfTies(
                &mut allocation_heap,
                allocation1,
                allocation2,
                3,
                0,
                0,
                &mut new_rank,
                &mut add_stack,
                &mut mapped_count
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
    }

    #[test]
    fn source_port__ichimap2__map_an_atom2__line_1350() {
        fn empty_neighbours(heap: &mut SourceHeap) -> SourceMutPointer<NEIGH_LIST> {
            let first = heap.allocate_model_storage(vec![0_u16]).unwrap();
            let second = heap.allocate_model_storage(vec![0_u16]).unwrap();
            heap.allocate_model_storage(vec![first, second]).unwrap()
        }

        let mut no_tie_heap = SourceHeap::default();
        let (no_tie1, no_tie2) = tie_fixture(
            &mut no_tie_heap,
            vec![1, 2],
            vec![0, 1],
            vec![1, 2],
            vec![0, 1],
            [SourceMutPointer::null(); 2],
            [SourceMutPointer::null(); 2],
        );
        let no_tie_temp = no_tie_heap
            .allocate_model_storage(vec![91_u16, 92])
            .unwrap();
        let no_tie_neighbours = empty_neighbours(&mut no_tie_heap);
        let mut no_tie_new_count = 77;
        let mut no_tie_add_stack = 77;
        let mut no_tie_stat = CANON_STAT {
            lNumBreakTies: 5,
            lNumNeighListIter: 6,
            ..CANON_STAT::default()
        };
        assert_eq!(
            map_an_atom2(
                &mut no_tie_heap,
                &mut CANON_GLOBALS::default(),
                2,
                2,
                0,
                0,
                no_tie_temp,
                2,
                &mut no_tie_new_count,
                &mut no_tie_stat,
                no_tie_neighbours,
                no_tie1,
                no_tie2,
                &mut no_tie_add_stack,
            ),
            Ok(1)
        );
        assert_eq!((no_tie_new_count, no_tie_add_stack), (2, 0));
        assert_eq!(
            (no_tie_stat.lNumBreakTies, no_tie_stat.lNumNeighListIter),
            (5, 6)
        );
        assert_eq!(
            no_tie_heap.slice(no_tie_temp.as_const()).unwrap(),
            &[91, 92]
        );

        let mut mismatch_heap = SourceHeap::default();
        let (mismatch1, mismatch2) = tie_fixture(
            &mut mismatch_heap,
            vec![2, 2],
            vec![0, 1],
            vec![1, 2],
            vec![0, 1],
            [SourceMutPointer::null(); 2],
            [SourceMutPointer::null(); 2],
        );
        let mismatch_temp = mismatch_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        let mismatch_neighbours = empty_neighbours(&mut mismatch_heap);
        let mut mismatch_new_count = 83;
        let mut mismatch_add_stack = 84;
        let mut mismatch_stat = CANON_STAT::default();
        assert_eq!(
            map_an_atom2(
                &mut mismatch_heap,
                &mut CANON_GLOBALS::default(),
                2,
                2,
                0,
                0,
                mismatch_temp,
                1,
                &mut mismatch_new_count,
                &mut mismatch_stat,
                mismatch_neighbours,
                mismatch1,
                mismatch2,
                &mut mismatch_add_stack,
            ),
            Ok(CT_MAPCOUNT_ERR)
        );
        assert_eq!((mismatch_new_count, mismatch_add_stack), (83, 0));
        assert_eq!(mismatch_stat, CANON_STAT::default());

        let mut rebuild_heap = SourceHeap::default();
        let rebuild_rank1 = rebuild_heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let rebuild_order1 = rebuild_heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let rebuild_rank2 = rebuild_heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let rebuild_order2 = rebuild_heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let previous_mapping = rebuild_heap
            .allocate_model_storage(vec![19_u16, 20])
            .unwrap();
        let rebuild_stack1 = rebuild_heap
            .allocate_model_storage(vec![
                rebuild_rank1,
                rebuild_order1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                previous_mapping,
                SourceMutPointer::null(),
            ])
            .unwrap();
        let rebuild_stack2 = rebuild_heap
            .allocate_model_storage(vec![
                rebuild_rank2,
                rebuild_order2,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ])
            .unwrap();
        let rebuild_temp = rebuild_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let rebuild_neighbours = empty_neighbours(&mut rebuild_heap);
        let mut rebuild_new_count = 71;
        let mut rebuild_add_stack = 72;
        let mut rebuild_stat = CANON_STAT {
            lNumBreakTies: i64::MAX,
            lNumNeighListIter: i64::MAX - 1,
            ..CANON_STAT::default()
        };
        assert_eq!(
            map_an_atom2(
                &mut rebuild_heap,
                &mut CANON_GLOBALS::default(),
                2,
                2,
                0,
                0,
                rebuild_temp,
                1,
                &mut rebuild_new_count,
                &mut rebuild_stat,
                rebuild_neighbours,
                rebuild_stack1,
                rebuild_stack2,
                &mut rebuild_add_stack,
            ),
            Ok(1)
        );
        assert_eq!((rebuild_new_count, rebuild_add_stack), (2, 2));
        assert_eq!(
            (rebuild_stat.lNumBreakTies, rebuild_stat.lNumNeighListIter),
            (i64::MIN + 1, i64::MIN)
        );
        assert_eq!(
            rebuild_heap.slice(previous_mapping.as_const()).unwrap(),
            &[0, 20]
        );
        let rebuild_slots = rebuild_heap.slice(rebuild_stack1.as_const()).unwrap();
        assert_eq!(
            rebuild_heap.slice(rebuild_slots[2].as_const()).unwrap(),
            &[1, 2]
        );
        assert_eq!(
            rebuild_heap.slice(rebuild_slots[3].as_const()).unwrap(),
            &[0, 1]
        );

        let mut reuse_heap = SourceHeap::default();
        let reuse_rank1 = reuse_heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let reuse_order1 = reuse_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let reuse_rank2 = reuse_heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let reuse_order2 = reuse_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let mapped_rank1 = reuse_heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let mapped_order1 = reuse_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let reuse_stack1 = reuse_heap
            .allocate_model_storage(vec![
                reuse_rank1,
                reuse_order1,
                mapped_rank1,
                mapped_order1,
                SourceMutPointer::null(),
            ])
            .unwrap();
        let reuse_stack2 = reuse_heap
            .allocate_model_storage(vec![
                reuse_rank2,
                reuse_order2,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ])
            .unwrap();
        let reuse_temp = reuse_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let reuse_neighbours = empty_neighbours(&mut reuse_heap);
        let mut reuse_new_count = 0;
        let mut reuse_add_stack = 0;
        let mut reuse_stat = CANON_STAT::default();
        assert_eq!(
            map_an_atom2(
                &mut reuse_heap,
                &mut CANON_GLOBALS::default(),
                2,
                2,
                1,
                1,
                reuse_temp,
                1,
                &mut reuse_new_count,
                &mut reuse_stat,
                reuse_neighbours,
                reuse_stack1,
                reuse_stack2,
                &mut reuse_add_stack,
            ),
            Ok(1)
        );
        assert_eq!((reuse_new_count, reuse_add_stack), (2, 2));
        assert_eq!(
            (reuse_stat.lNumBreakTies, reuse_stat.lNumNeighListIter),
            (1, 1)
        );
        assert_eq!(reuse_heap.slice(mapped_rank1.as_const()).unwrap(), &[2, 1]);
        assert_eq!(reuse_heap.slice(mapped_order1.as_const()).unwrap(), &[1, 0]);

        let mut allocation_heap = SourceHeap::default();
        let (allocation1, allocation2) = tie_fixture(
            &mut allocation_heap,
            vec![2, 2],
            vec![0, 1],
            vec![2, 2],
            vec![0, 1],
            [SourceMutPointer::null(); 2],
            [SourceMutPointer::null(); 2],
        );
        let allocation_temp = allocation_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        let allocation_neighbours = empty_neighbours(&mut allocation_heap);
        allocation_heap.fail_after_allocations(0);
        let mut allocation_new_count = 31;
        let mut allocation_add_stack = 32;
        let mut allocation_stat = CANON_STAT::default();
        assert_eq!(
            map_an_atom2(
                &mut allocation_heap,
                &mut CANON_GLOBALS::default(),
                2,
                2,
                0,
                0,
                allocation_temp,
                1,
                &mut allocation_new_count,
                &mut allocation_stat,
                allocation_neighbours,
                allocation1,
                allocation2,
                &mut allocation_add_stack,
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert_eq!((allocation_new_count, allocation_add_stack), (31, 0));
        assert_eq!(allocation_stat, CANON_STAT::default());
    }

    #[test]
    fn source_port__ichimap2__clearpreviousmappings__line_1334() {
        let mut heap = SourceHeap::default();
        let first = heap.allocate_model_storage(vec![7_u16, 8]).unwrap();
        let second = heap.allocate_model_storage(vec![9_u16]).unwrap();
        let stack = heap
            .allocate_model_storage(vec![first, second, SourceMutPointer::null()])
            .unwrap();
        assert_eq!(ClearPreviousMappings(&mut heap, stack), Ok(2));
        assert_eq!(heap.slice(first.as_const()).unwrap(), &[0, 8]);
        assert_eq!(heap.slice(second.as_const()).unwrap(), &[0]);

        let empty = heap
            .allocate_model_storage(vec![SourceMutPointer::<AT_RANK>::null()])
            .unwrap();
        assert_eq!(ClearPreviousMappings(&mut heap, empty), Ok(0));
        assert_eq!(
            ClearPreviousMappings(&mut heap, SourceMutPointer::null()),
            Err(SourceHeapError::NullPointer)
        );

        let prefix = heap.allocate_model_storage(vec![5_u16]).unwrap();
        let missing_terminator = heap.allocate_model_storage(vec![prefix]).unwrap();
        assert_eq!(
            ClearPreviousMappings(&mut heap, missing_terminator),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(prefix.as_const()).unwrap(), &[0]);

        let short = heap.allocate_model_storage(Vec::<AT_RANK>::new()).unwrap();
        let invalid_row = heap
            .allocate_model_storage(vec![short, SourceMutPointer::null()])
            .unwrap();
        assert_eq!(
            ClearPreviousMappings(&mut heap, invalid_row),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap2__halfstereobondparity__line_802() {
        fn run(center: sp_ATOM, ranks: Vec<AT_RANK>, slot: i32) -> i32 {
            let mut atoms = vec![sp_ATOM::default(); ranks.len()];
            atoms[0] = center.clone();
            let mut heap = SourceHeap::default();
            let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
            let rank_pointer = heap.allocate_model_storage(ranks).unwrap();
            HalfStereoBondParity(&heap, atom_pointer, 0, slot, rank_pointer).unwrap()
        }

        fn center(valence: i8, parity: i8, order: i8) -> sp_ATOM {
            let mut atom = sp_ATOM {
                valence,
                parity,
                neighbor: [1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                ..sp_ATOM::default()
            };
            atom.stereo_bond_neighbor[0] = 2;
            atom.stereo_bond_ord[0] = order;
            atom
        }

        assert_eq!(run(center(4, 1, 0), vec![0, 1, 2, 3], 0), 0);
        assert_eq!(run(center(3, 0, 0), vec![0, 1, 2, 3], 0), 0);
        assert_eq!(run(center(3, 3, 0), vec![0, 1, 2, 3], 0), 3);
        assert_eq!(run(center(3, 0x43, 0), vec![0, 1, 2, 3], 0), 0x43);
        assert_eq!(
            run(center(3, AB_PARITY_CALC as i8, 0), vec![0, 1, 2, 3], 0),
            -6
        );

        assert_eq!(
            run(center(3, 1, 0), vec![0, 1, 2, 3], -1),
            CT_STEREOBOND_ERROR
        );
        assert_eq!(
            run(center(3, 1, 0), vec![0, 1, 2, 3], 3),
            CT_STEREOBOND_ERROR
        );

        let mut missing = center(3, 1, 0);
        missing.stereo_bond_neighbor[1] = 0;
        assert_eq!(run(missing, vec![0, 1, 2, 3], 1), CT_STEREOBOND_ERROR);

        let mut second_slot = center(3, 1, 1);
        second_slot.stereo_bond_neighbor[1] = 3;
        second_slot.stereo_bond_ord[1] = 1;
        assert_eq!(run(second_slot, vec![0, 1, 2, 3], 1), 2);

        assert_eq!(run(center(3, 2, 1), vec![0, 1, 2, 3], 0), 1);
        assert_eq!(run(center(3, 2, 1), vec![0, 3, 2, 1], 0), 2);
        assert_eq!(
            run(center(3, 2, 1), vec![0, 4, 2, 4], 0),
            AB_PARITY_CALC as i32
        );
        assert_eq!(run(center(3, 2, 1), vec![0, 0, 2, 3], 0), 0);
        assert_eq!(run(center(3, 2, 1), vec![0, 1, 2, 0], 0), 0);

        assert_eq!(run(center(2, 2, 1), vec![0, 1, 2, 3], 0), 2);
        assert_eq!(run(center(1, 1, 0), vec![0, 1, 2, 3], 0), 1);

        let mut absent = center(3, 2, 1);
        absent.neighbor = [1; 20];
        assert_eq!(run(absent, vec![0, 1, 2, 3], 0), CT_STEREOBOND_ERROR);

        let mut duplicate = center(3, 2, 1);
        duplicate.neighbor[2] = 2;
        assert_eq!(run(duplicate, vec![0, 1, 2, 3], 0), CT_STEREOBOND_ERROR);
    }

    #[test]
    fn source_port__ichimap2__parity_of_mapped_half_bond__line_958() {
        fn fixture() -> (Vec<sp_ATOM>, Vec<AT_RANK>, Vec<AT_RANK>, Vec<AT_RANK>) {
            let mut atoms = vec![sp_ATOM::default(); 8];
            atoms[0].valence = 3;
            atoms[0].parity = 1;
            atoms[0].neighbor[..3].copy_from_slice(&[2, 4, 5]);
            atoms[0].stereo_bond_neighbor[0] = 3;
            atoms[0].stereo_bond_ord[0] = 0;
            atoms[1].valence = 3;
            atoms[1].neighbor[..3].copy_from_slice(&[3, 6, 7]);
            let mut canon_from = vec![0_u16; 8];
            canon_from[6] = 100;
            canon_from[7] = 200;
            let mut rank_from = vec![0_u16; 8];
            rank_from[1] = 10;
            rank_from[3] = 20;
            rank_from[6] = 30;
            rank_from[7] = 40;
            let mut rank_to = vec![0_u16; 8];
            rank_to[0] = 10;
            rank_to[2] = 20;
            rank_to[4] = 30;
            rank_to[5] = 40;
            (atoms, canon_from, rank_from, rank_to)
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            canon_from: Vec<AT_RANK>,
            rank_from: Vec<AT_RANK>,
            rank_to: Vec<AT_RANK>,
        ) -> (Result<i32, SourceHeapError>, EQ_NEIGH) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let canon_from = heap.allocate_model_storage(canon_from).unwrap();
            let rank_from = heap.allocate_model_storage(rank_from).unwrap();
            let rank_to = heap.allocate_model_storage(rank_to).unwrap();
            let mut en = EQ_NEIGH {
                num_to: 99,
                to_at: [99; 4],
                from_at: 99,
                rank: 99,
                canon_rank: 99,
            };
            let result = parity_of_mapped_half_bond(
                &heap,
                1,
                0,
                3,
                2,
                atoms,
                Some(&mut en),
                canon_from,
                rank_from,
                rank_to,
            );
            (result, en)
        }

        let (atoms, canon_from, rank_from, rank_to) = fixture();
        assert_eq!(
            run(
                atoms.clone(),
                canon_from.clone(),
                rank_from.clone(),
                rank_to.clone()
            ),
            (Ok(1), EQ_NEIGH::default())
        );

        let mut reversed_canon = canon_from.clone();
        reversed_canon[6] = 200;
        reversed_canon[7] = 100;
        assert_eq!(
            run(
                atoms.clone(),
                reversed_canon,
                rank_from.clone(),
                rank_to.clone()
            )
            .0,
            Ok(2)
        );

        let mut tied_rank_from = rank_from.clone();
        tied_rank_from[7] = 30;
        let mut tied_rank_to = rank_to.clone();
        tied_rank_to[5] = 30;
        let (tied_result, tied_en) = run(
            atoms.clone(),
            canon_from.clone(),
            tied_rank_from,
            tied_rank_to,
        );
        assert_eq!(tied_result, Ok(-30));
        assert_eq!(
            tied_en,
            EQ_NEIGH {
                num_to: 2,
                to_at: [4, 5, 0, 0],
                from_at: 6,
                rank: 30,
                canon_rank: 100,
            }
        );

        let mut two_neighbour_atoms = atoms.clone();
        two_neighbour_atoms[0].valence = 2;
        two_neighbour_atoms[1].valence = 2;
        assert_eq!(
            run(
                two_neighbour_atoms,
                canon_from.clone(),
                rank_from.clone(),
                rank_to.clone()
            )
            .0,
            Ok(2)
        );

        for (parity, expected) in [(1_i8, 1), (2, 2), (3, 3), (0, AB_PARITY_UNDF as i32)] {
            let mut terminal_atoms = atoms.clone();
            terminal_atoms[0].valence = 1;
            terminal_atoms[1].valence = 1;
            terminal_atoms[0].parity = parity;
            assert_eq!(
                run(
                    terminal_atoms,
                    canon_from.clone(),
                    rank_from.clone(),
                    rank_to.clone()
                )
                .0,
                Ok(expected)
            );
        }

        let mut mismatch_rank_to = rank_to.clone();
        mismatch_rank_to[0] = 11;
        let (mismatch_result, mismatch_en) = run(
            atoms.clone(),
            canon_from.clone(),
            rank_from.clone(),
            mismatch_rank_to,
        );
        assert_eq!((mismatch_result, mismatch_en), (Ok(0), EQ_NEIGH::default()));

        let mut invalid_valence = atoms.clone();
        invalid_valence[0].valence = 4;
        invalid_valence[1].valence = 4;
        assert_eq!(
            run(
                invalid_valence,
                canon_from.clone(),
                rank_from.clone(),
                rank_to.clone()
            )
            .0,
            Ok(0)
        );
        let mut terminal_without_bond = atoms.clone();
        terminal_without_bond[0].valence = 1;
        terminal_without_bond[1].valence = 1;
        terminal_without_bond[0].stereo_bond_neighbor[0] = 0;
        assert_eq!(
            run(
                terminal_without_bond,
                canon_from.clone(),
                rank_from.clone(),
                rank_to.clone()
            )
            .0,
            Ok(0)
        );

        for parity in [3_i8, 5, 0] {
            let mut parity_atoms = atoms.clone();
            parity_atoms[0].parity = parity;
            assert_eq!(
                run(
                    parity_atoms,
                    canon_from.clone(),
                    rank_from.clone(),
                    rank_to.clone()
                )
                .0,
                Ok(if parity == 3 { 3 } else { 0 })
            );
        }

        let mut missing_stereo_bond = atoms.clone();
        missing_stereo_bond[0].stereo_bond_neighbor[0] = 8;
        assert_eq!(
            run(
                missing_stereo_bond,
                canon_from.clone(),
                rank_from.clone(),
                rank_to.clone()
            )
            .0,
            Ok(0)
        );
        let mut tied_to_stereo_rank = rank_to.clone();
        tied_to_stereo_rank[4] = 20;
        assert_eq!(
            run(
                atoms.clone(),
                canon_from.clone(),
                rank_from.clone(),
                tied_to_stereo_rank
            )
            .0,
            Ok(0)
        );

        let mut unexpected_from_rank = rank_from.clone();
        unexpected_from_rank[7] = 41;
        assert_eq!(
            run(
                atoms.clone(),
                canon_from.clone(),
                unexpected_from_rank,
                rank_to.clone()
            )
            .0,
            Ok(0)
        );
        let mut missing_canon = canon_from.clone();
        missing_canon[7] = 0;
        assert_eq!(run(atoms, missing_canon, rank_from, rank_to).0, Ok(0));
    }

    #[test]
    fn source_port__ichimap2__parity_of_mapped_atom2__line_1168() {
        fn fixture() -> (Vec<sp_ATOM>, Vec<AT_RANK>, Vec<AT_RANK>, Vec<AT_RANK>) {
            let mut atoms = vec![sp_ATOM::default(); 8];
            atoms[0].valence = 3;
            atoms[0].neighbor[..3].copy_from_slice(&[5, 6, 7]);
            atoms[1].valence = 3;
            atoms[1].parity = 1;
            atoms[1].neighbor[..3].copy_from_slice(&[2, 3, 4]);
            let mut canon = vec![0_u16; 8];
            canon[5] = 300;
            canon[6] = 100;
            canon[7] = 200;
            let mut from = vec![0_u16; 8];
            from[0] = 10;
            from[5] = 40;
            from[6] = 20;
            from[7] = 30;
            let mut to = vec![0_u16; 8];
            to[1] = 10;
            to[2] = 20;
            to[3] = 30;
            to[4] = 40;
            (atoms, canon, from, to)
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            canon: Vec<AT_RANK>,
            from: Vec<AT_RANK>,
            to: Vec<AT_RANK>,
            with_en: bool,
        ) -> (Result<i32, SourceHeapError>, EQ_NEIGH, i32, bool, usize) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let canon = heap.allocate_model_storage(canon).unwrap();
            let from = heap.allocate_model_storage(from).unwrap();
            let to = heap.allocate_model_storage(to).unwrap();
            let baseline = heap.live_allocation_count();
            let mut en = EQ_NEIGH {
                num_to: 99,
                to_at: [99; 4],
                from_at: 99,
                rank: 99,
                canon_rank: 99,
            };
            let mut globals = CANON_GLOBALS::default();
            let result = parity_of_mapped_atom2(
                &mut heap,
                &mut globals,
                0,
                1,
                atoms,
                with_en.then_some(&mut en),
                canon,
                from,
                to,
            );
            let dangling = matches!(
                heap.slice(globals.m_pn_RankForSort),
                Err(SourceHeapError::MissingAllocation)
            );
            (
                result,
                en,
                globals.m_nNumCompNeighborsRanksCountEql,
                dangling,
                heap.live_allocation_count() - baseline,
            )
        }

        let (atoms, canon, from, to) = fixture();
        let distinct = run(atoms.clone(), canon.clone(), from.clone(), to.clone(), true);
        assert_eq!(distinct.0, Ok(1));
        assert_eq!(distinct.1, EQ_NEIGH::default());
        assert!(distinct.2 >= 0);
        assert!(distinct.3);
        assert_eq!(distinct.4, 0);

        let mut tied_from = from.clone();
        tied_from[5] = 20;
        let mut tied_to = to.clone();
        tied_to[3] = 20;
        let tied = run(
            atoms.clone(),
            canon.clone(),
            tied_from.clone(),
            tied_to.clone(),
            true,
        );
        assert_eq!(tied.0, Ok(-20));
        assert_eq!(
            tied.1,
            EQ_NEIGH {
                num_to: 2,
                to_at: [3, 2, 0, 0],
                from_at: 6,
                rank: 20,
                canon_rank: 100,
            }
        );
        assert!(tied.2 > 0);
        assert!(tied.3);
        assert_eq!(tied.4, 0);
        assert_eq!(
            run(atoms.clone(), canon.clone(), tied_from, tied_to, false).0,
            Ok(-20)
        );

        let mut unknown_atoms = atoms.clone();
        unknown_atoms[1].parity = 3;
        assert_eq!(
            run(unknown_atoms, canon.clone(), from.clone(), to.clone(), true).0,
            Ok(3)
        );
        for (parity, expected) in [(0_i8, AB_PARITY_UNDF as i32), (-1, -1)] {
            let mut terminal_atoms = atoms.clone();
            terminal_atoms[1].valence = 1;
            terminal_atoms[1].parity = parity;
            assert_eq!(
                run(
                    terminal_atoms,
                    canon.clone(),
                    from.clone(),
                    to.clone(),
                    true
                )
                .0,
                Ok(expected)
            );
        }

        let mut invalid_valence = atoms.clone();
        invalid_valence[1].valence = 5;
        assert_eq!(
            run(
                invalid_valence,
                canon.clone(),
                from.clone(),
                to.clone(),
                true
            )
            .0,
            Ok(0)
        );
        let mut central_mismatch = to.clone();
        central_mismatch[1] = 11;
        let mismatch = run(
            atoms.clone(),
            canon.clone(),
            from.clone(),
            central_mismatch,
            true,
        );
        assert_eq!(mismatch.0, Ok(0));
        assert_eq!(mismatch.1, EQ_NEIGH::default());
        assert!(!mismatch.3);
        assert_eq!(mismatch.4, 0);

        let mut neighbourhood_mismatch = to.clone();
        neighbourhood_mismatch[4] = 41;
        assert_eq!(
            run(
                atoms.clone(),
                canon.clone(),
                from.clone(),
                neighbourhood_mismatch,
                true
            )
            .0,
            Ok(0)
        );
        let mut zero_tied_from = from.clone();
        zero_tied_from[5] = 0;
        zero_tied_from[6] = 0;
        let mut zero_tied_to = to;
        zero_tied_to[2] = 0;
        zero_tied_to[3] = 0;
        assert_eq!(
            run(atoms, canon, zero_tied_from, zero_tied_to, true).0,
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimap2__might_change_other_atom_parity__line_1461() {
        fn run(
            atoms: Vec<sp_ATOM>,
            rank2: Vec<AT_RANK>,
            rank1: Vec<AT_RANK>,
            at_no: i32,
        ) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let rank2 = heap.allocate_model_storage(rank2).unwrap();
            let rank1 = heap.allocate_model_storage(rank1).unwrap();
            might_change_other_atom_parity(&heap, atoms, 3, at_no, rank2, rank1)
        }

        let atoms = vec![sp_ATOM::default(); 3];
        assert_eq!(run(atoms.clone(), vec![1, 2, 3], vec![1, 2, 3], 0), Ok(0));

        let mut direct = atoms.clone();
        direct[1].bHasStereoOrEquToStereo = 1;
        assert_eq!(run(direct.clone(), vec![1, 9, 3], vec![1, 2, 3], 0), Ok(1));
        direct[1].stereo_atom_parity = KNOWN_PARITIES_EQL as i8;
        assert_eq!(run(direct.clone(), vec![1, 9, 3], vec![1, 2, 3], 0), Ok(0));
        direct[1].stereo_atom_parity = 0;
        direct[1].stereo_bond_neighbor[0] = 1;
        assert_eq!(run(direct, vec![1, 9, 3], vec![1, 2, 3], 0), Ok(0));

        let mut neighbour = atoms.clone();
        neighbour[0].valence = 2;
        neighbour[0].neighbor[..2].copy_from_slice(&[1, 2]);
        neighbour[1].bHasStereoOrEquToStereo = 1;
        assert_eq!(
            run(neighbour.clone(), vec![9, 2, 3], vec![1, 2, 3], 0),
            Ok(1)
        );
        neighbour[1].stereo_atom_parity = KNOWN_PARITIES_EQL as i8;
        assert_eq!(
            run(neighbour.clone(), vec![9, 2, 3], vec![1, 2, 3], 0),
            Ok(0)
        );
        neighbour[1].stereo_atom_parity = 0;
        neighbour[1].stereo_bond_neighbor[0] = 2;
        assert_eq!(run(neighbour, vec![9, 2, 3], vec![1, 2, 3], 0), Ok(0));

        let mut negative_valence = atoms.clone();
        negative_valence[0].valence = -1;
        assert_eq!(
            run(negative_valence, vec![9, 2, 3], vec![1, 2, 3], 0),
            Ok(0)
        );

        let mut no_access_heap = SourceHeap::default();
        let equal2 = no_access_heap
            .allocate_model_storage(vec![1_u16, 2])
            .unwrap();
        let equal1 = no_access_heap
            .allocate_model_storage(vec![1_u16, 2])
            .unwrap();
        assert_eq!(
            might_change_other_atom_parity(
                &no_access_heap,
                SourceMutPointer::null(),
                2,
                0,
                equal2,
                equal1,
            ),
            Ok(0)
        );
        assert_eq!(
            might_change_other_atom_parity(
                &no_access_heap,
                SourceMutPointer::null(),
                -1,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        let short = no_access_heap.allocate_model_storage(vec![1_u16]).unwrap();
        assert_eq!(
            might_change_other_atom_parity(
                &no_access_heap,
                SourceMutPointer::null(),
                2,
                0,
                short,
                equal1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap2__deallocatefornonstereoremoval__line_1504() {
        fn neighbour_list(heap: &mut SourceHeap) -> SourceMutPointer<NEIGH_LIST> {
            let atoms = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
            heap.allocate_model_storage(vec![atoms]).unwrap()
        }

        let mut heap = SourceHeap::default();
        let mut canon1 = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let mut canon2 = heap.allocate_model_storage(vec![2_u16]).unwrap();
        let mut nl = neighbour_list(&mut heap);
        let mut nl1 = neighbour_list(&mut heap);
        let mut nl2 = neighbour_list(&mut heap);
        let mut visited1 = heap.allocate_model_storage(vec![3_u16]).unwrap();
        let mut visited2 = heap.allocate_model_storage(vec![4_u16]).unwrap();
        assert_eq!(heap.live_allocation_count(), 10);
        assert_eq!(
            DeAllocateForNonStereoRemoval(
                &mut heap,
                &mut canon1,
                &mut canon2,
                &mut nl,
                &mut nl1,
                &mut nl2,
                &mut visited1,
                &mut visited2,
            ),
            Ok(())
        );
        assert!(canon1.is_null());
        assert!(canon2.is_null());
        assert!(nl.is_null());
        assert!(nl1.is_null());
        assert!(nl2.is_null());
        assert!(visited1.is_null());
        assert!(visited2.is_null());
        assert_eq!(heap.live_allocation_count(), 0);
        assert_eq!(
            DeAllocateForNonStereoRemoval(
                &mut heap,
                &mut canon1,
                &mut canon2,
                &mut nl,
                &mut nl1,
                &mut nl2,
                &mut visited1,
                &mut visited2,
            ),
            Ok(())
        );

        let mut failed_heap = SourceHeap::default();
        let mut stale = failed_heap.allocate_model_storage(vec![1_u16]).unwrap();
        failed_heap.free(stale).unwrap();
        let mut later_rank = failed_heap.allocate_model_storage(vec![2_u16]).unwrap();
        let mut later_nl = neighbour_list(&mut failed_heap);
        let mut null_nl1 = SourceMutPointer::null();
        let mut null_nl2 = SourceMutPointer::null();
        let mut null_rank1 = SourceMutPointer::null();
        let mut null_rank2 = SourceMutPointer::null();
        assert_eq!(
            DeAllocateForNonStereoRemoval(
                &mut failed_heap,
                &mut stale,
                &mut later_rank,
                &mut later_nl,
                &mut null_nl1,
                &mut null_nl2,
                &mut null_rank1,
                &mut null_rank2,
            ),
            Err(SourceHeapError::MissingAllocation)
        );
        assert!(!stale.is_null());
        assert!(!later_rank.is_null());
        assert!(!later_nl.is_null());

        let mut partial_heap = SourceHeap::default();
        let row = partial_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut broken_nl = partial_heap.allocate_model_storage(vec![row]).unwrap();
        partial_heap.free(row).unwrap();
        let mut later_visited = partial_heap.allocate_model_storage(vec![7_u16]).unwrap();
        let mut null_rank1 = SourceMutPointer::null();
        let mut null_rank2 = SourceMutPointer::null();
        let mut null_rank3 = SourceMutPointer::null();
        let mut null_list1 = SourceMutPointer::null();
        let mut null_list2 = SourceMutPointer::null();
        assert_eq!(
            DeAllocateForNonStereoRemoval(
                &mut partial_heap,
                &mut null_rank1,
                &mut null_rank2,
                &mut broken_nl,
                &mut null_list1,
                &mut null_list2,
                &mut later_visited,
                &mut null_rank3,
            ),
            Err(SourceHeapError::MissingAllocation)
        );
        assert!(!broken_nl.is_null());
        assert!(!later_visited.is_null());
    }

    #[test]
    fn source_port__ichimap2__allocatefornonstereoremoval__line_1551() {
        fn inputs(
            heap: &mut SourceHeap,
        ) -> (
            SourceMutPointer<sp_ATOM>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
        ) {
            let mut atoms = vec![sp_ATOM::default(); 3];
            atoms[0].valence = 2;
            atoms[0].neighbor[0] = 1;
            atoms[0].neighbor[1] = 2;
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let symmetry = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
            let canonical = heap.allocate_model_storage(vec![3_u16, 2, 1]).unwrap();
            (atoms, symmetry, canonical)
        }

        let mut heap = SourceHeap::default();
        let (atoms, symmetry, canonical) = inputs(&mut heap);
        let mut canon1 = heap.allocate_model_storage(vec![91_u16]).unwrap();
        let mut canon2 = heap.allocate_model_storage(vec![92_u16]).unwrap();
        let old_row = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut lists = heap.allocate_model_storage(vec![old_row]).unwrap();
        let old_row1 = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut lists1 = heap.allocate_model_storage(vec![old_row1]).unwrap();
        let old_row2 = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut lists2 = heap.allocate_model_storage(vec![old_row2]).unwrap();
        let mut visited1 = heap.allocate_model_storage(vec![93_u16]).unwrap();
        let mut visited2 = heap.allocate_model_storage(vec![94_u16]).unwrap();
        assert_eq!(heap.live_allocation_count(), 13);
        assert_eq!(
            AllocateForNonStereoRemoval(
                &mut heap,
                atoms,
                3,
                symmetry,
                canonical,
                &mut canon1,
                &mut canon2,
                &mut lists,
                &mut lists1,
                &mut lists2,
                &mut visited1,
                &mut visited2,
            ),
            Ok(1)
        );
        assert_eq!(heap.live_allocation_count(), 13);
        assert_eq!(heap.slice(canon1.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(canon2.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(visited1.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(visited2.as_const()).unwrap(), &[0, 0, 0]);
        for pointer_list in [lists, lists1, lists2] {
            let rows = heap.slice(pointer_list.as_const()).unwrap();
            assert_eq!(rows.len(), 4);
            assert!(rows[3].is_null());
            assert_eq!(
                heap.slice(rows[0].as_const()).unwrap(),
                &[2, 2, 1, 1, 0, 1, 0, 0]
            );
            assert_eq!(&heap.slice(rows[1].as_const()).unwrap()[..2], &[1, 0]);
            assert_eq!(&heap.slice(rows[2].as_const()).unwrap()[..2], &[1, 0]);
        }

        for successful_allocations in 0..10 {
            let mut failed_heap = SourceHeap::default();
            let (atoms, symmetry, canonical) = inputs(&mut failed_heap);
            let baseline = failed_heap.live_allocation_count();
            failed_heap.fail_after_allocations(successful_allocations);
            let mut canon1 = SourceMutPointer::null();
            let mut canon2 = SourceMutPointer::null();
            let mut lists = SourceMutPointer::null();
            let mut lists1 = SourceMutPointer::null();
            let mut lists2 = SourceMutPointer::null();
            let mut visited1 = SourceMutPointer::null();
            let mut visited2 = SourceMutPointer::null();
            assert_eq!(
                AllocateForNonStereoRemoval(
                    &mut failed_heap,
                    atoms,
                    3,
                    symmetry,
                    canonical,
                    &mut canon1,
                    &mut canon2,
                    &mut lists,
                    &mut lists1,
                    &mut lists2,
                    &mut visited1,
                    &mut visited2,
                ),
                Ok(0),
                "failure after {successful_allocations} successful allocations"
            );
            assert!(canon1.is_null());
            assert!(canon2.is_null());
            assert!(lists.is_null());
            assert!(lists1.is_null());
            assert!(lists2.is_null());
            assert!(visited1.is_null());
            assert!(visited2.is_null());
            assert_eq!(failed_heap.live_allocation_count(), baseline);
        }
    }

    #[test]
    fn source_port__ichimap2__getminnewrank__line_1587() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![2_u16, 4, 4, 7, 7])
            .unwrap();
        let atoms = heap
            .allocate_model_storage(vec![0_u16, 2, 1, 4, 3])
            .unwrap();
        assert_eq!(GetMinNewRank(&heap, ranks, atoms, 0), Ok(1));
        assert_eq!(GetMinNewRank(&heap, ranks, atoms, 1), Ok(3));
        assert_eq!(GetMinNewRank(&heap, ranks, atoms, 2), Ok(5));
        assert_eq!(GetMinNewRank(&heap, ranks, atoms, 3), Ok(5));
        assert_eq!(GetMinNewRank(&heap, ranks, atoms, 4), Ok(8));
        assert_eq!(GetMinNewRank(&heap, ranks, atoms, 5), Ok(8));

        let suffix_ranks = heap.allocate_model_storage(vec![1_u16, 2, 4, 4]).unwrap();
        let suffix_atoms = heap.allocate_model_storage(vec![0_u16, 1, 2, 3]).unwrap();
        assert_eq!(GetMinNewRank(&heap, suffix_ranks, suffix_atoms, 4), Ok(3));

        let all_tied_ranks = heap.allocate_model_storage(vec![3_u16; 3]).unwrap();
        let all_tied_atoms = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        assert_eq!(
            GetMinNewRank(&heap, all_tied_ranks, all_tied_atoms, 3),
            Ok(1)
        );

        let wrapping_ranks = heap.allocate_model_storage(vec![u16::MAX, 9_u16]).unwrap();
        let wrapping_atoms = heap.allocate_model_storage(vec![0_u16]).unwrap();
        assert_eq!(
            GetMinNewRank(&heap, wrapping_ranks, wrapping_atoms, 1),
            Ok(0)
        );

        let bad_atoms = heap.allocate_model_storage(vec![9_u16]).unwrap();
        assert_eq!(
            GetMinNewRank(&heap, ranks, bad_atoms, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            GetMinNewRank(&heap, SourceMutPointer::null(), SourceMutPointer::null(), 0,),
            Ok(1)
        );
    }

    #[test]
    fn source_port__ichimap2__sortedequinfotoranks__line_148() {
        let mut empty_heap = SourceHeap::default();
        let mut empty_changed = 7;
        assert_eq!(
            SortedEquInfoToRanks(
                &mut empty_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                Some(&mut empty_changed),
            ),
            Ok(1)
        );
        assert_eq!(empty_changed, 0);

        let mut heap = SourceHeap::default();
        let symmetry = heap
            .allocate_model_storage(vec![10_u16, 30, 10, 20])
            .unwrap();
        let ranks = heap.allocate_model_storage(vec![99_u16; 4]).unwrap();
        let atom_numbers = heap.allocate_model_storage(vec![2_u16, 0, 3, 1]).unwrap();
        let mut changed = 0;
        assert_eq!(
            SortedEquInfoToRanks(
                &mut heap,
                symmetry,
                ranks,
                atom_numbers,
                4,
                Some(&mut changed),
            ),
            Ok(3)
        );
        assert_eq!(heap.slice(ranks.as_const()).unwrap(), &[2, 4, 2, 3]);
        assert_eq!(changed, 1);

        let alias = heap.allocate_model_storage(vec![2_u16, 4, 2, 3]).unwrap();
        let mut alias_changed = 9;
        assert_eq!(
            SortedEquInfoToRanks(
                &mut heap,
                alias,
                alias,
                atom_numbers,
                4,
                Some(&mut alias_changed),
            ),
            Ok(3)
        );
        assert_eq!(heap.slice(alias.as_const()).unwrap(), &[2, 4, 2, 3]);
        assert_eq!(alias_changed, 0);

        let bad_order = heap.allocate_model_storage(vec![2_u16, 0, 9, 1]).unwrap();
        let partial = heap.allocate_model_storage(vec![88_u16; 4]).unwrap();
        let mut unchanged_on_error = 5;
        assert_eq!(
            SortedEquInfoToRanks(
                &mut heap,
                symmetry,
                partial,
                bad_order,
                4,
                Some(&mut unchanged_on_error),
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), &[88, 4, 88, 88]);
        assert_eq!(unchanged_on_error, 5);

        assert_eq!(
            SortedEquInfoToRanks(&mut heap, symmetry, ranks, atom_numbers, -1, None,),
            Ok(1)
        );
    }

    #[test]
    fn source_port__ichimap2__sortedrankstoequinfo__line_199() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate_model_storage(vec![2_u16, 4, 2, 3]).unwrap();
        let symmetry = heap.allocate_model_storage(vec![99_u16; 4]).unwrap();
        let atom_numbers = heap.allocate_model_storage(vec![2_u16, 0, 3, 1]).unwrap();
        assert_eq!(
            SortedRanksToEquInfo(&mut heap, symmetry, ranks, atom_numbers, 4),
            Ok(3)
        );
        assert_eq!(heap.slice(symmetry.as_const()).unwrap(), &[1, 4, 1, 3]);

        let alias = heap.allocate_model_storage(vec![2_u16, 4, 2, 3]).unwrap();
        assert_eq!(
            SortedRanksToEquInfo(&mut heap, alias, alias, atom_numbers, 4),
            Ok(3)
        );
        assert_eq!(heap.slice(alias.as_const()).unwrap(), &[1, 4, 1, 3]);

        let one_rank = heap.allocate_model_storage(vec![77_u16]).unwrap();
        let one_output = heap.allocate_model_storage(vec![88_u16]).unwrap();
        let one_number = heap.allocate_model_storage(vec![0_u16]).unwrap();
        assert_eq!(
            SortedRanksToEquInfo(&mut heap, one_output, one_rank, one_number, 1),
            Ok(1)
        );
        assert_eq!(heap.slice(one_output.as_const()).unwrap(), &[1]);

        let bad_order = heap.allocate_model_storage(vec![2_u16, 0, 9, 1]).unwrap();
        let partial = heap.allocate_model_storage(vec![88_u16; 4]).unwrap();
        assert_eq!(
            SortedRanksToEquInfo(&mut heap, partial, ranks, bad_order, 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), &[1, 88, 1, 88]);

        let mut empty_heap = SourceHeap::default();
        assert_eq!(
            SortedRanksToEquInfo(
                &mut empty_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimap2__switch_ptrs__line_230() {
        let mut heap = SourceHeap::default();
        let first_base = heap.allocate(vec![10_u16, 11, 12]).unwrap();
        let second_base = heap.allocate(vec![20_u16, 21, 22, 23]).unwrap();
        let mut first = first_base.offset(1).unwrap();
        let mut second = second_base.offset(2).unwrap();

        switch_ptrs(&mut first, &mut second);
        assert_eq!(first, second_base.offset(2).unwrap());
        assert_eq!(second, first_base.offset(1).unwrap());
        assert_eq!(heap.slice(first.as_const()).unwrap(), &[22, 23]);
        assert_eq!(heap.slice(second.as_const()).unwrap(), &[11, 12]);
        assert_eq!(heap.slice(first_base.as_const()).unwrap(), &[10, 11, 12]);
        assert_eq!(
            heap.slice(second_base.as_const()).unwrap(),
            &[20, 21, 22, 23]
        );

        let mut null = SourceMutPointer::<AT_RANK>::null();
        switch_ptrs(&mut first, &mut null);
        assert!(first.is_null());
        assert_eq!(null, second_base.offset(2).unwrap());
    }

    #[test]
    fn source_port__ichimap2__setnewranksfromneighlists3__line_241() {
        let mut empty_heap = SourceHeap::default();
        let mut empty_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetNewRanksFromNeighLists3(
                &mut empty_heap,
                &mut empty_globals,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert!(empty_globals.m_pNeighList_RankForSort.is_null());
        assert!(empty_globals.m_pn_RankForSort.is_null());

        let mut singleton_heap = SourceHeap::default();
        let singleton_ranks = singleton_heap.allocate(vec![2_u16, 3, 1]).unwrap();
        let singleton_new = singleton_heap.allocate(vec![9_u16; 3]).unwrap();
        let singleton_atoms = singleton_heap.allocate(vec![2_u16, 0, 1]).unwrap();
        let mut singleton_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetNewRanksFromNeighLists3(
                &mut singleton_heap,
                &mut singleton_globals,
                3,
                SourceMutPointer::null(),
                singleton_ranks,
                singleton_new,
                singleton_atoms,
            ),
            Ok(3)
        );
        assert_eq!(
            singleton_heap.slice(singleton_new.as_const()).unwrap(),
            &[2, 3, 1]
        );
        assert_eq!(
            singleton_heap.slice(singleton_atoms.as_const()).unwrap(),
            &[2, 0, 1]
        );

        let mut tied_heap = SourceHeap::default();
        let list0 = tied_heap.allocate(vec![0_u16]).unwrap();
        let list1 = tied_heap.allocate(vec![1_u16, 0]).unwrap();
        let list2 = tied_heap.allocate(vec![2_u16, 0, 0]).unwrap();
        let list3 = tied_heap.allocate(vec![1_u16, 0]).unwrap();
        let lists = tied_heap
            .allocate(vec![list0, list1, list2, list3])
            .unwrap();
        let tied_ranks = tied_heap.allocate(vec![4_u16; 4]).unwrap();
        let tied_new = tied_heap.allocate(vec![9_u16; 4]).unwrap();
        let tied_atoms = tied_heap.allocate(vec![2_u16, 0, 3, 1]).unwrap();
        let mut tied_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetNewRanksFromNeighLists3(
                &mut tied_heap,
                &mut tied_globals,
                4,
                lists,
                tied_ranks,
                tied_new,
                tied_atoms,
            ),
            Ok(-3)
        );
        assert_eq!(
            tied_heap.slice(tied_atoms.as_const()).unwrap(),
            &[0, 3, 1, 2]
        );
        assert_eq!(tied_heap.slice(tied_new.as_const()).unwrap(), &[1, 3, 4, 3]);
        assert_eq!(tied_globals.m_pNeighList_RankForSort, lists.as_const());
        assert_eq!(tied_globals.m_pn_RankForSort, tied_ranks.as_const());

        let equal_lists = tied_heap
            .allocate(vec![list1, list1, list1, list1])
            .unwrap();
        let equal_new = tied_heap.allocate(vec![8_u16; 4]).unwrap();
        let equal_atoms = tied_heap.allocate(vec![3_u16, 2, 1, 0]).unwrap();
        assert_eq!(
            SetNewRanksFromNeighLists3(
                &mut tied_heap,
                &mut tied_globals,
                4,
                equal_lists,
                tied_ranks,
                equal_new,
                equal_atoms,
            ),
            Ok(1)
        );
        assert_eq!(
            tied_heap.slice(equal_atoms.as_const()).unwrap(),
            &[3, 2, 1, 0]
        );
        assert_eq!(
            tied_heap.slice(equal_new.as_const()).unwrap(),
            &[4, 4, 4, 4]
        );

        let mut partial_heap = SourceHeap::default();
        let partial_new = partial_heap.allocate(vec![7_u16, 8]).unwrap();
        let partial_lists = partial_heap.allocate(Vec::<NEIGH_LIST>::new()).unwrap();
        let partial_ranks = partial_heap.allocate(Vec::<AT_RANK>::new()).unwrap();
        let partial_atoms = partial_heap.allocate(Vec::<AT_RANK>::new()).unwrap();
        let mut partial_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetNewRanksFromNeighLists3(
                &mut partial_heap,
                &mut partial_globals,
                3,
                partial_lists,
                partial_ranks,
                partial_new,
                partial_atoms,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(partial_heap.slice(partial_new.as_const()).unwrap(), &[0, 0]);
        assert_eq!(
            partial_globals.m_pNeighList_RankForSort,
            partial_lists.as_const()
        );
        assert_eq!(partial_globals.m_pn_RankForSort, partial_ranks.as_const());
    }

    #[test]
    fn source_port__ichimap2__setnewranksfromneighlists4__line_308() {
        let mut heap = SourceHeap::default();
        let list0 = heap.allocate(vec![0_u16]).unwrap();
        let list1 = heap.allocate(vec![1_u16, 0]).unwrap();
        let list2 = heap.allocate(vec![2_u16, 0, 0]).unwrap();
        let list3 = heap.allocate(vec![1_u16, 0]).unwrap();
        let lists = heap.allocate(vec![list0, list1, list2, list3]).unwrap();
        let ranks = heap.allocate(vec![4_u16; 4]).unwrap();
        let new_ranks = heap.allocate(vec![9_u16; 4]).unwrap();
        let atoms = heap.allocate(vec![2_u16, 0, 3, 1]).unwrap();
        let mut globals = CANON_GLOBALS::default();

        assert_eq!(
            SetNewRanksFromNeighLists4(
                &mut heap,
                &mut globals,
                4,
                lists,
                ranks,
                new_ranks,
                atoms,
                2,
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), &[2, 0, 3, 1]);
        assert_eq!(heap.slice(new_ranks.as_const()).unwrap(), &[4, 4, 4, 4]);
        assert_eq!(globals.m_nMaxAtNeighRankForSort, 2);
        assert_eq!(globals.m_pNeighList_RankForSort, lists.as_const());
        assert_eq!(globals.m_pn_RankForSort, ranks.as_const());

        heap.slice_mut(new_ranks).unwrap().fill(9);
        assert_eq!(
            SetNewRanksFromNeighLists4(
                &mut heap,
                &mut globals,
                4,
                lists,
                ranks,
                new_ranks,
                atoms,
                4,
            ),
            Ok(-3)
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), &[0, 3, 1, 2]);
        assert_eq!(heap.slice(new_ranks.as_const()).unwrap(), &[1, 3, 4, 3]);
        assert_eq!(globals.m_nMaxAtNeighRankForSort, 4);

        let singleton_ranks = heap.allocate(vec![1_u16, 2]).unwrap();
        let singleton_new = heap.allocate(vec![8_u16, 8]).unwrap();
        let singleton_atoms = heap.allocate(vec![0_u16, 1]).unwrap();
        assert_eq!(
            SetNewRanksFromNeighLists4(
                &mut heap,
                &mut globals,
                2,
                SourceMutPointer::null(),
                singleton_ranks,
                singleton_new,
                singleton_atoms,
                u16::MAX,
            ),
            Ok(2)
        );
        assert_eq!(heap.slice(singleton_new.as_const()).unwrap(), &[1, 2]);

        let mut partial_heap = SourceHeap::default();
        let partial_new = partial_heap.allocate(vec![7_u16, 8]).unwrap();
        let partial_lists = partial_heap.allocate(Vec::<NEIGH_LIST>::new()).unwrap();
        let partial_ranks = partial_heap.allocate(Vec::<AT_RANK>::new()).unwrap();
        let partial_atoms = partial_heap.allocate(Vec::<AT_RANK>::new()).unwrap();
        let mut partial_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetNewRanksFromNeighLists4(
                &mut partial_heap,
                &mut partial_globals,
                3,
                partial_lists,
                partial_ranks,
                partial_new,
                partial_atoms,
                17,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(partial_heap.slice(partial_new.as_const()).unwrap(), &[0, 0]);
        assert_eq!(partial_globals.m_nMaxAtNeighRankForSort, 17);
    }

    #[test]
    fn source_port__ichimap2__setnewranksfromneighlists__line_380() {
        let mut empty_heap = SourceHeap::default();
        let mut empty_globals = CANON_GLOBALS::default();
        let mut empty_calls = 0;
        assert_eq!(
            SetNewRanksFromNeighLists(
                &mut empty_heap,
                &mut empty_globals,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                &mut |_, _, _, _| {
                    empty_calls += 1;
                    Ok(0)
                },
            ),
            Ok(1)
        );
        assert_eq!(empty_calls, 0);
        assert!(empty_globals.m_pNeighList_RankForSort.is_null());
        assert!(empty_globals.m_pn_RankForSort.is_null());

        for alternate_sort in [0, 1, 3] {
            let mut heap = SourceHeap::default();
            let list0 = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
            let list1 = heap.allocate_model_storage(vec![0_u16]).unwrap();
            let list2 = heap.allocate_model_storage(vec![1_u16, 3]).unwrap();
            let list3 = heap.allocate_model_storage(vec![0_u16]).unwrap();
            let lists = heap
                .allocate_model_storage(vec![list0, list1, list2, list3])
                .unwrap();
            let ranks = heap.allocate_model_storage(vec![3_u16, 1, 3, 4]).unwrap();
            let new_ranks = heap.allocate_model_storage(vec![99_u16; 4]).unwrap();
            let atoms = heap.allocate_model_storage(vec![3_u16, 2, 1, 0]).unwrap();
            let mut globals = CANON_GLOBALS::default();
            let mut calls = 0_i32;
            assert_eq!(
                SetNewRanksFromNeighLists(
                    &mut heap,
                    &mut globals,
                    4,
                    lists,
                    ranks,
                    new_ranks,
                    atoms,
                    alternate_sort,
                    &mut |heap, left, right, globals| {
                        calls = calls.wrapping_add(1);
                        CompNeighListRanks(heap, left, right, globals)
                    },
                ),
                Ok(4)
            );
            assert!(calls > 0);
            assert_eq!(heap.slice(atoms.as_const()).unwrap(), &[1, 0, 2, 3]);
            assert_eq!(heap.slice(new_ranks.as_const()).unwrap(), &[2, 1, 3, 4]);
            assert_eq!(globals.m_pNeighList_RankForSort, lists.as_const());
            assert_eq!(globals.m_pn_RankForSort, ranks.as_const());
        }

        let mut tied_heap = SourceHeap::default();
        let tied_list = tied_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let tied_lists = tied_heap
            .allocate_model_storage(vec![tied_list, tied_list])
            .unwrap();
        let tied_ranks = tied_heap.allocate_model_storage(vec![2_u16; 2]).unwrap();
        let tied_new = tied_heap.allocate_model_storage(vec![9_u16; 2]).unwrap();
        let tied_atoms = tied_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let mut tied_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetNewRanksFromNeighLists(
                &mut tied_heap,
                &mut tied_globals,
                2,
                tied_lists,
                tied_ranks,
                tied_new,
                tied_atoms,
                1,
                &mut CompNeighListRanks,
            ),
            Ok(1)
        );
        assert_eq!(tied_heap.slice(tied_new.as_const()).unwrap(), &[2, 2]);

        assert_eq!(
            SetNewRanksFromNeighLists(
                &mut tied_heap,
                &mut tied_globals,
                -1,
                tied_lists,
                tied_ranks,
                tied_new,
                tied_atoms,
                0,
                &mut CompNeighListRanks,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__ichimap2__sortneighlistsbysymmandcanonrank__line_434() {
        let mut heap = SourceHeap::default();
        let symmetry = heap.allocate_model_storage(vec![1_u16, 3, 3]).unwrap();
        let canonical = heap.allocate_model_storage(vec![9_u16, 1, 5]).unwrap();
        let first = heap
            .allocate_model_storage(vec![3_u16, 0, 1, 2, 77])
            .unwrap();
        let second = heap.allocate_model_storage(vec![2_u16, 0, 2, 88]).unwrap();
        let empty = heap.allocate_model_storage(vec![0_u16, 99]).unwrap();
        let lists = heap
            .allocate_model_storage(vec![first, second, empty])
            .unwrap();

        assert_eq!(
            SortNeighListsBySymmAndCanonRank(&mut heap, 3, lists, symmetry, canonical,),
            Ok(())
        );
        assert_eq!(heap.slice(first.as_const()).unwrap(), &[3, 2, 1, 0, 77]);
        assert_eq!(heap.slice(second.as_const()).unwrap(), &[2, 2, 0, 88]);
        assert_eq!(heap.slice(empty.as_const()).unwrap(), &[0, 99]);

        assert_eq!(
            SortNeighListsBySymmAndCanonRank(
                &mut heap,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(())
        );
        assert_eq!(
            SortNeighListsBySymmAndCanonRank(
                &mut heap,
                -1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(())
        );

        let partial_symmetry = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let partial_canonical = heap.allocate_model_storage(vec![0_u16, 0]).unwrap();
        let partial = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let partial_lists = heap
            .allocate_model_storage(vec![partial, SourceMutPointer::null()])
            .unwrap();
        assert_eq!(
            SortNeighListsBySymmAndCanonRank(
                &mut heap,
                2,
                partial_lists,
                partial_symmetry,
                partial_canonical,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), &[2, 1, 0]);
    }

    #[test]
    fn source_port__ichimap2__sortneighlists2__line_448() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate(vec![1_u16, 2, 4, 4, 5]).unwrap();
        let untied = heap.allocate(vec![3_u16, 2, 1, 4]).unwrap();
        let singleton = heap.allocate(vec![1_u16, 4]).unwrap();
        let equal_ranks = heap.allocate(vec![2_u16, 3, 2]).unwrap();
        let tied = heap.allocate(vec![3_u16, 4, 2, 1]).unwrap();
        let empty = heap.allocate(vec![0_u16]).unwrap();
        let lists = heap
            .allocate(vec![untied, singleton, equal_ranks, tied, empty])
            .unwrap();
        let atoms = heap.allocate(vec![0_u16, 1, 3, 2, 4]).unwrap();

        assert_eq!(SortNeighLists2(&mut heap, 5, ranks, lists, atoms), Ok(0));
        assert_eq!(heap.slice(untied.as_const()).unwrap(), &[3, 1, 2, 4]);
        assert_eq!(heap.slice(singleton.as_const()).unwrap(), &[1, 4]);
        assert_eq!(heap.slice(equal_ranks.as_const()).unwrap(), &[2, 3, 2]);
        assert_eq!(heap.slice(tied.as_const()).unwrap(), &[3, 1, 2, 4]);
        assert_eq!(heap.slice(empty.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), &[0, 1, 3, 2, 4]);

        assert_eq!(
            SortNeighLists2(
                &mut heap,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            SortNeighLists2(
                &mut heap,
                -1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let partial_ranks = heap.allocate(vec![3_u16, 1, 2]).unwrap();
        let partial_list = heap.allocate(vec![2_u16, 0, 1]).unwrap();
        let partial_lists = heap.allocate(vec![partial_list]).unwrap();
        let partial_atoms = heap.allocate(vec![0_u16, 1]).unwrap();
        assert_eq!(
            SortNeighLists2(&mut heap, 2, partial_ranks, partial_lists, partial_atoms,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial_list.as_const()).unwrap(), &[2, 1, 0]);
    }

    #[test]
    fn source_port__ichimap2__sortneighlists3__line_482() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate(vec![2_u16, 2, 3, 5, 5]).unwrap();
        let list0 = heap.allocate(vec![3_u16, 2, 0, 1]).unwrap();
        let list1 = heap.allocate(vec![1_u16, 2]).unwrap();
        let list3 = heap.allocate(vec![2_u16, 2, 0]).unwrap();
        let list4 = heap.allocate(vec![0_u16]).unwrap();
        let lists = heap
            .allocate(vec![list0, list1, SourceMutPointer::null(), list3, list4])
            .unwrap();
        let atoms = heap.allocate(vec![0_u16, 1, 2, 3, 4]).unwrap();

        assert_eq!(SortNeighLists3(&mut heap, 5, ranks, lists, atoms), Ok(0));
        assert_eq!(heap.slice(list0.as_const()).unwrap(), &[3, 0, 1, 2]);
        assert_eq!(heap.slice(list1.as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(list3.as_const()).unwrap(), &[2, 0, 2]);
        assert_eq!(heap.slice(list4.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), &[0, 1, 2, 3, 4]);

        assert_eq!(
            SortNeighLists3(
                &mut heap,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            SortNeighLists3(
                &mut heap,
                -1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let partial_ranks = heap.allocate(vec![2_u16, 2, 1]).unwrap();
        let partial_list = heap.allocate(vec![2_u16, 0, 2]).unwrap();
        let partial_lists = heap.allocate(vec![partial_list]).unwrap();
        let partial_atoms = heap.allocate(vec![0_u16, 1]).unwrap();
        assert_eq!(
            SortNeighLists3(&mut heap, 2, partial_ranks, partial_lists, partial_atoms,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial_list.as_const()).unwrap(), &[2, 2, 0]);
    }

    #[test]
    fn source_port__ichimap2__differentiateranks2__line_518() {
        fn fixture(
            heap: &mut SourceHeap,
        ) -> (
            SourceMutPointer<NEIGH_LIST>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
        ) {
            let list0 = heap.allocate(vec![0_u16]).unwrap();
            let list1 = heap.allocate(vec![1_u16, 0]).unwrap();
            let list2 = heap.allocate(vec![2_u16, 0, 0]).unwrap();
            let list3 = heap.allocate(vec![1_u16, 0]).unwrap();
            (
                heap.allocate(vec![list0, list1, list2, list3]).unwrap(),
                heap.allocate(vec![4_u16; 4]).unwrap(),
                heap.allocate(vec![0_u16; 4]).unwrap(),
                heap.allocate(vec![2_u16, 0, 3, 1]).unwrap(),
            )
        }

        for alternate_sort in [0, 1, 3] {
            let mut heap = SourceHeap::default();
            let (lists, current, previous, atoms) = fixture(&mut heap);
            let mut globals = CANON_GLOBALS::default();
            let mut iterations = -2_i64;
            assert_eq!(
                DifferentiateRanks2(
                    &mut heap,
                    &mut globals,
                    4,
                    lists,
                    99,
                    current,
                    previous,
                    atoms,
                    &mut iterations,
                    alternate_sort,
                ),
                Ok(3)
            );
            assert_eq!(iterations, 0);
            assert_eq!(heap.slice(current.as_const()).unwrap(), &[1, 3, 4, 3]);
            assert_eq!(heap.slice(previous.as_const()).unwrap(), &[1, 3, 4, 3]);
            assert_eq!(heap.slice(atoms.as_const()).unwrap(), &[0, 1, 3, 2]);
            assert_eq!(globals.m_pNeighList_RankForSort, lists.as_const());
            assert_eq!(globals.m_pn_RankForSort, previous.as_const());
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty_globals = CANON_GLOBALS::default();
        let mut empty_iterations = 10_i64;
        assert_eq!(
            DifferentiateRanks2(
                &mut empty_heap,
                &mut empty_globals,
                0,
                SourceMutPointer::null(),
                -7,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty_iterations,
                0,
            ),
            Ok(1)
        );
        assert_eq!(empty_iterations, 11);

        let mut invalid_heap = SourceHeap::default();
        let invalid_atoms = invalid_heap.allocate(vec![0_u16]).unwrap();
        let mut invalid_globals = CANON_GLOBALS::default();
        let mut invalid_iterations = 7_i64;
        assert_eq!(
            DifferentiateRanks2(
                &mut invalid_heap,
                &mut invalid_globals,
                1,
                SourceMutPointer::null(),
                1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                invalid_atoms,
                &mut invalid_iterations,
                1,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(invalid_iterations, 8);
        assert!(invalid_globals.m_pn_RankForSort.is_null());

        assert_eq!(
            DifferentiateRanks2(
                &mut invalid_heap,
                &mut invalid_globals,
                -1,
                SourceMutPointer::null(),
                1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut invalid_iterations,
                0,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(invalid_iterations, 8);
    }

    #[test]
    fn source_port__ichimap2__differentiateranks3__line_561() {
        let mut singleton_heap = SourceHeap::default();
        let singleton_current = singleton_heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let singleton_previous = singleton_heap.allocate(vec![9_u16; 3]).unwrap();
        let singleton_atoms = singleton_heap.allocate(vec![0_u16, 1, 2]).unwrap();
        let mut singleton_globals = CANON_GLOBALS::default();
        let mut singleton_iterations = 7_i64;
        assert_eq!(
            DifferentiateRanks3(
                &mut singleton_heap,
                &mut singleton_globals,
                3,
                SourceMutPointer::null(),
                99,
                singleton_current,
                singleton_previous,
                singleton_atoms,
                &mut singleton_iterations,
            ),
            Ok(3)
        );
        assert_eq!(singleton_iterations, 8);
        assert_eq!(
            singleton_heap.slice(singleton_current.as_const()).unwrap(),
            &[1, 2, 3]
        );
        assert_eq!(
            singleton_heap.slice(singleton_previous.as_const()).unwrap(),
            &[1, 2, 3]
        );

        let mut tied_heap = SourceHeap::default();
        let list0 = tied_heap.allocate(vec![0_u16]).unwrap();
        let list1 = tied_heap.allocate(vec![1_u16, 0]).unwrap();
        let list2 = tied_heap.allocate(vec![2_u16, 0, 0]).unwrap();
        let list3 = tied_heap.allocate(vec![1_u16, 0]).unwrap();
        let lists = tied_heap
            .allocate(vec![list0, list1, list2, list3])
            .unwrap();
        let tied_current = tied_heap.allocate(vec![4_u16; 4]).unwrap();
        let tied_previous = tied_heap.allocate(vec![0_u16; 4]).unwrap();
        let tied_atoms = tied_heap.allocate(vec![2_u16, 0, 3, 1]).unwrap();
        let mut tied_globals = CANON_GLOBALS::default();
        let mut tied_iterations = -2_i64;
        assert_eq!(
            DifferentiateRanks3(
                &mut tied_heap,
                &mut tied_globals,
                4,
                lists,
                1,
                tied_current,
                tied_previous,
                tied_atoms,
                &mut tied_iterations,
            ),
            Ok(3)
        );
        assert_eq!(tied_iterations, 0);
        assert_eq!(
            tied_heap.slice(tied_current.as_const()).unwrap(),
            &[1, 3, 4, 3]
        );
        assert_eq!(
            tied_heap.slice(tied_previous.as_const()).unwrap(),
            &[1, 3, 4, 3]
        );
        assert_eq!(
            tied_heap.slice(tied_atoms.as_const()).unwrap(),
            &[0, 3, 1, 2]
        );

        let mut empty_iterations = 10_i64;
        assert_eq!(
            DifferentiateRanks3(
                &mut tied_heap,
                &mut tied_globals,
                0,
                SourceMutPointer::null(),
                -7,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty_iterations,
            ),
            Ok(0)
        );
        assert_eq!(empty_iterations, 11);

        let mut error_heap = SourceHeap::default();
        let error_atoms = error_heap.allocate(vec![0_u16]).unwrap();
        let error_scratch = error_heap.allocate(vec![8_u16]).unwrap();
        let mut error_globals = CANON_GLOBALS::default();
        let mut error_iterations = 4_i64;
        assert_eq!(
            DifferentiateRanks3(
                &mut error_heap,
                &mut error_globals,
                1,
                SourceMutPointer::null(),
                1,
                SourceMutPointer::null(),
                error_scratch,
                error_atoms,
                &mut error_iterations,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(error_iterations, 5);
        assert_eq!(error_heap.slice(error_scratch.as_const()).unwrap(), &[8]);
    }

    #[test]
    fn source_port__ichimap2__differentiateranks4__line_602() {
        fn fixture(
            heap: &mut SourceHeap,
        ) -> (
            SourceMutPointer<NEIGH_LIST>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
        ) {
            let list0 = heap.allocate(vec![0_u16]).unwrap();
            let list1 = heap.allocate(vec![1_u16, 0]).unwrap();
            let list2 = heap.allocate(vec![2_u16, 0, 0]).unwrap();
            let list3 = heap.allocate(vec![1_u16, 0]).unwrap();
            (
                heap.allocate(vec![list0, list1, list2, list3]).unwrap(),
                heap.allocate(vec![4_u16; 4]).unwrap(),
                heap.allocate(vec![0_u16; 4]).unwrap(),
                heap.allocate(vec![2_u16, 0, 3, 1]).unwrap(),
            )
        }

        let mut low_heap = SourceHeap::default();
        let (low_lists, low_current, low_previous, low_atoms) = fixture(&mut low_heap);
        let mut low_globals = CANON_GLOBALS::default();
        let mut low_iterations = 0_i64;
        assert_eq!(
            DifferentiateRanks4(
                &mut low_heap,
                &mut low_globals,
                4,
                low_lists,
                4,
                low_current,
                low_previous,
                low_atoms,
                2,
                &mut low_iterations,
            ),
            Ok(1)
        );
        assert_eq!(low_iterations, 1);
        assert_eq!(low_heap.slice(low_current.as_const()).unwrap(), &[4; 4]);
        assert_eq!(low_heap.slice(low_previous.as_const()).unwrap(), &[4; 4]);
        assert_eq!(low_heap.slice(low_atoms.as_const()).unwrap(), &[2, 0, 3, 1]);
        assert_eq!(low_globals.m_nMaxAtNeighRankForSort, 2);

        let mut full_heap = SourceHeap::default();
        let (full_lists, full_current, full_previous, full_atoms) = fixture(&mut full_heap);
        let mut full_globals = CANON_GLOBALS::default();
        let mut full_iterations = 11_i64;
        assert_eq!(
            DifferentiateRanks4(
                &mut full_heap,
                &mut full_globals,
                4,
                full_lists,
                -3,
                full_current,
                full_previous,
                full_atoms,
                4,
                &mut full_iterations,
            ),
            Ok(3)
        );
        assert_eq!(full_iterations, 13);
        assert_eq!(
            full_heap.slice(full_current.as_const()).unwrap(),
            &[1, 3, 4, 3]
        );
        assert_eq!(
            full_heap.slice(full_previous.as_const()).unwrap(),
            &[1, 3, 4, 3]
        );
        assert_eq!(
            full_heap.slice(full_atoms.as_const()).unwrap(),
            &[0, 3, 1, 2]
        );
        assert_eq!(full_globals.m_nMaxAtNeighRankForSort, 4);

        let mut error_heap = SourceHeap::default();
        let error_atoms = error_heap.allocate(vec![0_u16]).unwrap();
        let error_scratch = error_heap.allocate(vec![8_u16]).unwrap();
        let mut error_globals = CANON_GLOBALS::default();
        let mut error_iterations = -1_i64;
        assert_eq!(
            DifferentiateRanks4(
                &mut error_heap,
                &mut error_globals,
                1,
                SourceMutPointer::null(),
                1,
                SourceMutPointer::null(),
                error_scratch,
                error_atoms,
                9,
                &mut error_iterations,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(error_iterations, 0);
        assert_eq!(error_heap.slice(error_scratch.as_const()).unwrap(), &[8]);
    }

    #[test]
    fn source_port__ichimap2__differentiateranksbasic__line_637() {
        fn fixture(
            heap: &mut SourceHeap,
        ) -> (
            SourceMutPointer<NEIGH_LIST>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
        ) {
            let list0 = heap.allocate_model_storage(vec![0_u16]).unwrap();
            let list1 = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
            let list2 = heap.allocate_model_storage(vec![2_u16, 0, 0]).unwrap();
            let list3 = heap.allocate_model_storage(vec![3_u16, 0, 0, 0]).unwrap();
            (
                heap.allocate_model_storage(vec![list0, list1, list2, list3])
                    .unwrap(),
                heap.allocate_model_storage(vec![4_u16; 4]).unwrap(),
                heap.allocate_model_storage(vec![0_u16; 4]).unwrap(),
                heap.allocate_model_storage(vec![2_u16, 0, 3, 1]).unwrap(),
            )
        }

        for alternate_sort in [0, 1, 3] {
            let mut heap = SourceHeap::default();
            let (lists, current, previous, atoms) = fixture(&mut heap);
            let mut globals = CANON_GLOBALS::default();
            let mut iterations = -2_i64;
            assert_eq!(
                DifferentiateRanksBasic(
                    &mut heap,
                    &mut globals,
                    4,
                    lists,
                    1,
                    current,
                    previous,
                    atoms,
                    &mut iterations,
                    alternate_sort,
                ),
                Ok(4)
            );
            assert_eq!(iterations, 0);
            assert_eq!(heap.slice(current.as_const()).unwrap(), &[1, 2, 3, 4]);
            assert_eq!(heap.slice(previous.as_const()).unwrap(), &[1, 2, 3, 4]);
            assert_eq!(heap.slice(atoms.as_const()).unwrap(), &[0, 1, 2, 3]);
            assert_eq!(globals.m_pNeighList_RankForSort, lists.as_const());
            assert_eq!(globals.m_pn_RankForSort, previous.as_const());
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty_globals = CANON_GLOBALS::default();
        let mut empty_iterations = 10_i64;
        assert_eq!(
            DifferentiateRanksBasic(
                &mut empty_heap,
                &mut empty_globals,
                0,
                SourceMutPointer::null(),
                -7,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty_iterations,
                0,
            ),
            Ok(1)
        );
        assert_eq!(empty_iterations, 12);

        let mut error_heap = SourceHeap::default();
        let error_atoms = error_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut error_globals = CANON_GLOBALS::default();
        let mut error_iterations = 7_i64;
        assert_eq!(
            DifferentiateRanksBasic(
                &mut error_heap,
                &mut error_globals,
                1,
                SourceMutPointer::null(),
                1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                error_atoms,
                &mut error_iterations,
                1,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(error_iterations, 8);
        assert!(error_globals.m_pn_RankForSort.is_null());

        assert_eq!(
            DifferentiateRanksBasic(
                &mut error_heap,
                &mut error_globals,
                -1,
                SourceMutPointer::null(),
                1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut error_iterations,
                0,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(error_iterations, 8);
    }

    #[test]
    fn source_port__ichimap2__breakneighborstie__line_1609() {
        fn empty_lists(heap: &mut SourceHeap, count: usize) -> SourceMutPointer<NEIGH_LIST> {
            let mut lists = Vec::with_capacity(count);
            for _ in 0..count {
                lists.push(heap.allocate_model_storage(vec![0_u16]).unwrap());
            }
            heap.allocate_model_storage(lists).unwrap()
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default(), sp_ATOM::default()])
            .unwrap();
        let neighbors = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let symmetry = heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let canonical = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let rank1 = heap.allocate_model_storage(vec![99_u16; 2]).unwrap();
        let order1 = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let rank2 = heap.allocate_model_storage(vec![88_u16; 2]).unwrap();
        let order2 = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let stack1 = heap.allocate_model_storage(vec![rank1, order1]).unwrap();
        let stack2 = heap.allocate_model_storage(vec![rank2, order2]).unwrap();
        let scratch = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let graph = empty_lists(&mut heap, 2);
        let lists1 = empty_lists(&mut heap, 2);
        let lists2 = empty_lists(&mut heap, 2);
        let mut globals = CANON_GLOBALS::default();
        let mut iterations = 0_i64;
        assert_eq!(
            BreakNeighborsTie(
                &mut heap,
                &mut globals,
                atoms,
                2,
                2,
                2,
                0,
                neighbors,
                0,
                1,
                MAP_MODE_STD as i32,
                stack1,
                stack2,
                scratch,
                graph,
                symmetry,
                canonical,
                lists1,
                lists2,
                &mut iterations,
            ),
            Ok(3)
        );
        assert_eq!(iterations, 4);
        assert_eq!(heap.slice(rank1.as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(rank2.as_const()).unwrap(), &[2, 1]);
        assert_eq!(heap.slice(order1.as_const()).unwrap(), &[0, 1]);
        assert_eq!(heap.slice(order2.as_const()).unwrap(), &[1, 0]);
        assert_eq!(globals.m_pn_RankForSort, scratch.as_const());

        let mut early_heap = SourceHeap::default();
        let early_atoms = early_heap
            .allocate_model_storage(vec![sp_ATOM::default(), sp_ATOM::default()])
            .unwrap();
        let early_neighbors = early_heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let early_symmetry = early_heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let mut early_globals = CANON_GLOBALS::default();
        let mut early_iterations = 17_i64;
        assert_eq!(
            BreakNeighborsTie(
                &mut early_heap,
                &mut early_globals,
                early_atoms,
                2,
                2,
                2,
                0,
                early_neighbors,
                1,
                0,
                MAP_MODE_S4 as i32,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                early_symmetry,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut early_iterations,
            ),
            Ok(0)
        );
        assert_eq!(early_iterations, 17);
        assert!(early_globals.m_pn_RankForSort.is_null());

        assert_eq!(
            BreakNeighborsTie(
                &mut early_heap,
                &mut early_globals,
                early_atoms,
                2,
                2,
                2,
                0,
                early_neighbors,
                0,
                1,
                MAP_MODE_C2 as i32,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                early_symmetry,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut early_iterations,
            ),
            Ok(0)
        );

        let unequal_symmetry = early_heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let unequal_rank1 = early_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let unequal_order1 = early_heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let unequal_stack1 = early_heap
            .allocate_model_storage(vec![unequal_rank1, unequal_order1])
            .unwrap();
        let unequal_graph = empty_lists(&mut early_heap, 2);
        assert_eq!(
            BreakNeighborsTie(
                &mut early_heap,
                &mut early_globals,
                early_atoms,
                2,
                2,
                2,
                0,
                early_neighbors,
                0,
                1,
                MAP_MODE_STD as i32,
                unequal_stack1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                unequal_graph,
                unequal_symmetry,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut early_iterations,
            ),
            Ok(0)
        );
        assert_eq!(early_heap.slice(unequal_rank1.as_const()).unwrap(), &[1, 2]);
    }

    #[test]
    fn source_port__ichimap2__checknextsymmneighborsandbonds__line_2060() {
        let mut heap = SourceHeap::default();
        let mut atom_values = vec![sp_ATOM::default(); 5];
        atom_values[0].neighbor[0] = 2;
        atom_values[0].stereo_bond_neighbor[0] = 4;
        atom_values[0].stereo_bond_ord[0] = 0;
        atom_values[0].stereo_bond_parity[0] = 1;
        atom_values[1].neighbor[0] = 3;
        atom_values[1].stereo_bond_neighbor[0] = 5;
        atom_values[1].stereo_bond_ord[0] = 0;
        atom_values[1].stereo_bond_parity[0] = 1;
        let atoms = heap.allocate_model_storage(atom_values).unwrap();
        let avoid = heap.allocate_model_storage(vec![9_u16, 9]).unwrap();
        let visited1 = heap.allocate_model_storage(vec![0_u16; 5]).unwrap();
        let visited2 = heap.allocate_model_storage(vec![0_u16; 5]).unwrap();
        let order1 = heap
            .allocate_model_storage(vec![0_u16, 0, 7, 0, 0])
            .unwrap();
        let order2 = heap
            .allocate_model_storage(vec![0_u16, 0, 0, 7, 0])
            .unwrap();
        let rank1 = heap
            .allocate_model_storage(vec![0_u16, 0, 6, 0, 0])
            .unwrap();
        let rank2 = heap
            .allocate_model_storage(vec![0_u16, 0, 0, 6, 0])
            .unwrap();

        macro_rules! check {
            () => {
                CheckNextSymmNeighborsAndBonds(
                    &heap, atoms, 0, 1, 2, 3, avoid, visited1, visited2, order1, order2, rank1,
                    rank2,
                )
            };
        }
        assert_eq!(check!(), Ok(1));

        source_set(&mut heap, rank2, 3, 5).unwrap();
        assert_eq!(check!(), Ok(-1));
        source_set(&mut heap, rank2, 3, 6).unwrap();

        source_set(&mut heap, visited1, 2, 4).unwrap();
        assert_eq!(check!(), Ok(-1));
        source_set(&mut heap, visited2, 3, 3).unwrap();
        assert_eq!(check!(), Ok(1));
        source_set(&mut heap, visited2, 3, 9).unwrap();
        assert_eq!(check!(), Ok(-1));
        source_set(&mut heap, visited1, 2, 0).unwrap();
        source_set(&mut heap, visited2, 3, 0).unwrap();

        source_set(&mut heap, order2, 3, 8).unwrap();
        assert_eq!(check!(), Ok(-1));
        source_set(&mut heap, order2, 3, 7).unwrap();

        heap.slice_mut(atoms).unwrap()[1].stereo_bond_neighbor[0] = 0;
        assert_eq!(check!(), Ok(0));
        heap.slice_mut(atoms).unwrap()[1].stereo_bond_neighbor[0] = 5;

        source_set(&mut heap, avoid, 0, 0).unwrap();
        source_set(&mut heap, avoid, 1, 3).unwrap();
        assert_eq!(check!(), Ok(0));
        source_set(&mut heap, avoid, 0, 9).unwrap();
        source_set(&mut heap, avoid, 1, 9).unwrap();

        heap.slice_mut(atoms).unwrap()[1].stereo_bond_parity[0] = 2;
        assert_eq!(check!(), Ok(0));
        heap.slice_mut(atoms).unwrap()[0].stereo_bond_parity[0] = 3;
        heap.slice_mut(atoms).unwrap()[1].stereo_bond_parity[0] = 5;
        assert_eq!(check!(), Ok(1));

        assert_eq!(
            CheckNextSymmNeighborsAndBonds(
                &heap,
                atoms,
                0,
                1,
                2,
                3,
                SourceMutPointer::null(),
                visited1,
                visited2,
                order1,
                order2,
                rank1,
                rank2,
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimap2__createchecksymmpaths__line_2191() {
        let mut terminal_heap = SourceHeap::default();
        let mut terminal_atoms = vec![sp_ATOM::default(); 2];
        terminal_atoms[0].valence = 1;
        terminal_atoms[1].valence = 1;
        terminal_atoms[0].stereo_atom_parity = 1;
        terminal_atoms[1].stereo_atom_parity = 2;
        let terminal_atoms = terminal_heap
            .allocate_model_storage(terminal_atoms)
            .unwrap();
        let visited1 = terminal_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        let visited2 = terminal_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        let order1 = terminal_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        let order2 = terminal_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        let mut length = 0_u16;
        let mut inverted = -1_i32;
        assert_eq!(
            CreateCheckSymmPaths(
                &mut terminal_heap,
                terminal_atoms,
                0,
                0,
                0,
                1,
                SourceMutPointer::null(),
                visited1,
                visited2,
                order1,
                order2,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut length,
                &mut inverted,
                0,
            ),
            Ok(1)
        );
        assert_eq!((length, inverted), (1, 1));
        assert_eq!(terminal_heap.slice(visited1.as_const()).unwrap(), &[2, 0]);
        assert_eq!(terminal_heap.slice(visited2.as_const()).unwrap(), &[0, 1]);
        assert_eq!(terminal_heap.slice(order1.as_const()).unwrap(), &[1, 0]);
        assert_eq!(terminal_heap.slice(order2.as_const()).unwrap(), &[0, 1]);

        source_set(&mut terminal_heap, visited1, 0, 0).unwrap();
        source_set(&mut terminal_heap, visited2, 1, 0).unwrap();
        length = 7;
        inverted = 0;
        assert_eq!(
            CreateCheckSymmPaths(
                &mut terminal_heap,
                terminal_atoms,
                0,
                0,
                0,
                1,
                SourceMutPointer::null(),
                visited1,
                visited2,
                order1,
                order2,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut length,
                &mut inverted,
                0,
            ),
            Ok(0)
        );
        assert_eq!(length, 8);
        assert_eq!(terminal_heap.slice(visited1.as_const()).unwrap()[0], 2);

        terminal_heap.slice_mut(terminal_atoms).unwrap()[1].valence = 2;
        inverted = -1;
        assert_eq!(
            CreateCheckSymmPaths(
                &mut terminal_heap,
                terminal_atoms,
                0,
                0,
                0,
                1,
                SourceMutPointer::null(),
                visited1,
                visited2,
                order1,
                order2,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut length,
                &mut inverted,
                0,
            ),
            Ok(CT_REMOVE_STEREO_ERR)
        );

        let mut recursive_heap = SourceHeap::default();
        let mut atoms = vec![sp_ATOM::default(); 5];
        atoms[0].valence = 2;
        atoms[1].valence = 2;
        atoms[2].valence = 1;
        atoms[3].valence = 1;
        atoms[0].nRingSystem = 7;
        atoms[2].nRingSystem = 7;
        atoms[1].nRingSystem = 8;
        atoms[3].nRingSystem = 8;
        let atoms = recursive_heap.allocate_model_storage(atoms).unwrap();
        let list0 = recursive_heap
            .allocate_model_storage(vec![2_u16, 4, 2])
            .unwrap();
        let list1 = recursive_heap
            .allocate_model_storage(vec![2_u16, 4, 3])
            .unwrap();
        let dummy = recursive_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let lists1 = recursive_heap
            .allocate_model_storage(vec![list0, dummy, dummy, dummy, dummy])
            .unwrap();
        let lists2 = recursive_heap
            .allocate_model_storage(vec![dummy, list1, dummy, dummy, dummy])
            .unwrap();
        let avoid = recursive_heap
            .allocate_model_storage(vec![9_u16, 9])
            .unwrap();
        let visited1 = recursive_heap
            .allocate_model_storage(vec![0_u16; 5])
            .unwrap();
        let visited2 = recursive_heap
            .allocate_model_storage(vec![0_u16; 5])
            .unwrap();
        let order1 = recursive_heap
            .allocate_model_storage(vec![0_u16; 5])
            .unwrap();
        let order2 = recursive_heap
            .allocate_model_storage(vec![0_u16; 5])
            .unwrap();
        let ranks1 = recursive_heap
            .allocate_model_storage(vec![1_u16, 1, 5, 0, 0])
            .unwrap();
        let ranks2 = recursive_heap
            .allocate_model_storage(vec![1_u16, 1, 0, 5, 0])
            .unwrap();
        let mut recursive_length = 0_u16;
        let mut recursive_inverted = -1_i32;
        assert_eq!(
            CreateCheckSymmPaths(
                &mut recursive_heap,
                atoms,
                4,
                0,
                4,
                1,
                avoid,
                visited1,
                visited2,
                order1,
                order2,
                lists1,
                lists2,
                ranks1,
                ranks2,
                SourceMutPointer::null(),
                &mut recursive_length,
                &mut recursive_inverted,
                0,
            ),
            Ok(1)
        );
        assert_eq!(recursive_length, 2);
        assert_eq!(
            recursive_heap.slice(visited1.as_const()).unwrap(),
            &[2, 0, 4, 0, 0]
        );
        assert_eq!(
            recursive_heap.slice(visited2.as_const()).unwrap(),
            &[0, 1, 0, 3, 0]
        );
        assert_eq!(
            recursive_heap.slice(order1.as_const()).unwrap(),
            &[1, 0, 2, 0, 0]
        );
        assert_eq!(
            recursive_heap.slice(order2.as_const()).unwrap(),
            &[0, 1, 0, 2, 0]
        );
    }

    #[test]
    fn source_port__ichimap2__removecalculatednonstereocenterparities__line_3356() {
        use crate::source_types::{AB_PARITY_ODD, AB_PARITY_UNKN, AT_STEREO_CARB};

        fn candidate(extra_atom: bool) -> Vec<sp_ATOM> {
            let mut atoms = vec![sp_ATOM::default(); if extra_atom { 4 } else { 3 }];
            atoms[0].parity = AB_PARITY_CALC as i8;
            atoms[0].stereo_atom_parity = AB_PARITY_CALC as i8;
            atoms[0].final_parity = AB_PARITY_CALC as i8;
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].nRingSystem = 1;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].nRingSystem = 2;
            atoms
        }

        fn carb(at_num: AT_RANK, parity: i32) -> AT_STEREO_CARB {
            AT_STEREO_CARB {
                at_num,
                parity: parity as u8,
            }
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            symmetry: Vec<AT_RANK>,
            canonical: Vec<AT_RANK>,
            atom_by_canonical: Vec<AT_RANK>,
            carb_entries: Vec<AT_STEREO_CARB>,
        ) -> (i32, Vec<sp_ATOM>, Vec<AT_STEREO_CARB>, i32) {
            let mut heap = SourceHeap::default();
            let count = atoms.len();
            let atoms_pointer = heap.allocate_model_storage(atoms).unwrap();
            let symmetry_pointer = heap.allocate_model_storage(symmetry).unwrap();
            let canonical_pointer = heap.allocate_model_storage(canonical).unwrap();
            let atom_by_canonical_pointer = heap.allocate_model_storage(atom_by_canonical).unwrap();
            let rank1 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let order1 = heap
                .allocate_model_storage((0..count as AT_RANK).collect())
                .unwrap();
            let rank2 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let order2 = heap
                .allocate_model_storage((0..count as AT_RANK).collect())
                .unwrap();
            let rank_stack1 = heap.allocate_model_storage(vec![rank1, order1]).unwrap();
            let rank_stack2 = heap.allocate_model_storage(vec![rank2, order2]).unwrap();
            let temporary_rank = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let mut atom_by_canonical1 = SourceMutPointer::null();
            let mut atom_by_canonical2 = SourceMutPointer::null();
            let mut lists = SourceMutPointer::null();
            let mut lists1 = SourceMutPointer::null();
            let mut lists2 = SourceMutPointer::null();
            let mut visited1 = SourceMutPointer::null();
            let mut visited2 = SourceMutPointer::null();
            assert_eq!(
                AllocateForNonStereoRemoval(
                    &mut heap,
                    atoms_pointer,
                    count as i32,
                    symmetry_pointer,
                    canonical_pointer,
                    &mut atom_by_canonical1,
                    &mut atom_by_canonical2,
                    &mut lists,
                    &mut lists1,
                    &mut lists2,
                    &mut visited1,
                    &mut visited2,
                ),
                Ok(1)
            );
            let carb_pointer = heap.allocate_model_storage(carb_entries).unwrap();
            let mut stats = CANON_STAT::default();
            stats.LinearCTStereoCarb = carb_pointer;
            stats.nLenLinearCTStereoCarb =
                heap.slice(carb_pointer.as_const()).unwrap().len() as i32;
            stats.nSymmRank = symmetry_pointer;
            let mut globals = CANON_GLOBALS::default();
            let result = RemoveCalculatedNonStereoCenterParities(
                &mut heap,
                &mut globals,
                atoms_pointer,
                count as i32,
                count as i32,
                rank_stack1,
                rank_stack2,
                temporary_rank,
                lists,
                canonical_pointer,
                symmetry_pointer,
                atom_by_canonical_pointer,
                atom_by_canonical1,
                atom_by_canonical2,
                lists,
                lists1,
                lists2,
                visited1,
                visited2,
                &mut stats,
                AB_PARITY_UNKN as i32,
            )
            .unwrap();
            (
                result,
                heap.slice(atoms_pointer.as_const()).unwrap().to_vec(),
                heap.slice(carb_pointer.as_const()).unwrap().to_vec(),
                stats.nLenLinearCTStereoCarb,
            )
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        let mut globals = CANON_GLOBALS::default();
        let mut stats = CANON_STAT::default();
        let call_without_arrays =
            |heap: &mut SourceHeap, globals: &mut CANON_GLOBALS, stats: &mut CANON_STAT| {
                RemoveCalculatedNonStereoCenterParities(
                    heap,
                    globals,
                    atoms,
                    1,
                    1,
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    stats,
                    AB_PARITY_UNKN as i32,
                )
            };
        assert_eq!(
            call_without_arrays(&mut heap, &mut globals, &mut stats),
            Ok(0)
        );

        {
            let atom = &mut heap.slice_mut(atoms).unwrap()[0];
            atom.parity = 1;
            atom.stereo_atom_parity = AB_PARITY_CALC as i8;
            atom.stereo_bond_neighbor[0] = 1;
        }
        assert_eq!(
            call_without_arrays(&mut heap, &mut globals, &mut stats),
            Ok(0)
        );

        {
            let atom = &mut heap.slice_mut(atoms).unwrap()[0];
            atom.stereo_bond_neighbor[0] = 0;
            atom.valence = 5;
        }
        assert_eq!(
            call_without_arrays(&mut heap, &mut globals, &mut stats),
            Ok(0)
        );

        {
            let atom = &mut heap.slice_mut(atoms).unwrap()[0];
            atom.valence = -1;
        }
        assert_eq!(
            call_without_arrays(&mut heap, &mut globals, &mut stats),
            Ok(0)
        );

        {
            let atom = &mut heap.slice_mut(atoms).unwrap()[0];
            atom.valence = 2;
            atom.stereo_atom_parity = AB_PARITY_ODD as i8;
        }
        assert_eq!(
            call_without_arrays(&mut heap, &mut globals, &mut stats),
            Ok(0)
        );

        let removed = run(
            candidate(true),
            vec![1, 2, 2, 4],
            vec![1, 2, 3, 4],
            vec![0, 1, 2, 3],
            vec![
                carb(1, AB_PARITY_CALC as i32),
                carb(4, AB_PARITY_ODD as i32),
            ],
        );
        assert_eq!(removed.0, 1);
        assert_eq!(removed.1[0].parity, 0);
        assert_eq!(removed.1[0].stereo_atom_parity, 0);
        assert_eq!(removed.1[0].final_parity, 0);
        assert_eq!(removed.3, 1);
        assert_eq!(removed.2[0], carb(4, AB_PARITY_ODD as i32));

        let mut ring_atoms = candidate(false);
        ring_atoms[1].nRingSystem = 5;
        ring_atoms[2].nRingSystem = 5;
        let ring = run(
            ring_atoms,
            vec![1, 2, 2],
            vec![1, 2, 3],
            vec![0, 1, 2],
            vec![carb(1, AB_PARITY_CALC as i32)],
        );
        assert_eq!(ring.0, 1);
        assert_eq!(ring.1[0].parity, 0);
        assert_eq!(ring.3, 0);

        let missing = run(
            candidate(false),
            vec![1, 2, 2],
            vec![1, 2, 3],
            vec![0, 1, 2],
            vec![],
        );
        assert_eq!(missing.0, CT_STEREOCOUNT_ERR);
        assert_eq!(missing.1[0].parity, 0);
        assert_eq!(missing.3, 0);

        let mut known = candidate(true);
        known[0].stereo_atom_parity |= KNOWN_PARITIES_EQL as i8;
        known[3].stereo_atom_parity = KNOWN_PARITIES_EQL as i8;
        let known = run(
            known,
            vec![7, 2, 2, 7],
            vec![1, 2, 3, 4],
            vec![0, 1, 2, 3],
            vec![carb(1, AB_PARITY_CALC as i32)],
        );
        assert_eq!(known.0, 1);
        assert_eq!(known.1[0].stereo_atom_parity, 0);
        assert_eq!(known.1[3].stereo_atom_parity & KNOWN_PARITIES_EQL as i8, 0);

        let mut propagated = candidate(false);
        propagated[1].parity = AB_PARITY_UNKN as i8;
        propagated[1].stereo_atom_parity = AB_PARITY_UNKN as i8;
        propagated[2].parity = AB_PARITY_UNDF as i8;
        propagated[2].stereo_atom_parity = AB_PARITY_UNDF as i8;
        let propagated = run(
            propagated,
            vec![1, 2, 2],
            vec![1, 2, 3],
            vec![0, 1, 2],
            vec![
                carb(1, AB_PARITY_CALC as i32),
                carb(2, AB_PARITY_UNKN as i32),
            ],
        );
        assert_eq!(propagated.0, 0);
        assert_eq!(propagated.1[0].parity, AB_PARITY_CALC as i8);
        assert_eq!(propagated.1[0].stereo_atom_parity, AB_PARITY_CALC as i8);
        assert_eq!(propagated.2[0].parity, AB_PARITY_CALC as u8);
        assert_eq!(propagated.3, 2);
    }

    #[test]
    fn source_port__ichimap2__removecalculatednonstereo__line_3672() {
        use crate::source_types::{AB_PARITY_ODD, AB_PARITY_UNKN, AT_STEREO_CARB};

        fn candidate() -> Vec<sp_ATOM> {
            let mut atoms = vec![sp_ATOM::default(); 4];
            atoms[0].parity = AB_PARITY_CALC as i8;
            atoms[0].stereo_atom_parity = AB_PARITY_CALC as i8;
            atoms[0].final_parity = AB_PARITY_CALC as i8;
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].nRingSystem = 1;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].nRingSystem = 2;
            atoms
        }

        fn run(
            atoms: Vec<sp_ATOM>,
            symmetry: Vec<AT_RANK>,
            canonical: Vec<AT_RANK>,
            atom_by_canonical: Vec<AT_RANK>,
            carb_entries: Vec<AT_STEREO_CARB>,
        ) -> (i32, Vec<sp_ATOM>, i32, usize, usize) {
            let mut heap = SourceHeap::default();
            let count = atoms.len();
            let atoms_pointer = heap.allocate_model_storage(atoms).unwrap();
            let symmetry_pointer = heap.allocate_model_storage(symmetry).unwrap();
            let canonical_pointer = heap.allocate_model_storage(canonical).unwrap();
            let atom_by_canonical_pointer = heap.allocate_model_storage(atom_by_canonical).unwrap();
            let rank1 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let order1 = heap
                .allocate_model_storage((0..count as AT_RANK).collect())
                .unwrap();
            let rank2 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let order2 = heap
                .allocate_model_storage((0..count as AT_RANK).collect())
                .unwrap();
            let rank_stack1 = heap.allocate_model_storage(vec![rank1, order1]).unwrap();
            let rank_stack2 = heap.allocate_model_storage(vec![rank2, order2]).unwrap();
            let temporary_rank = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let neighbor_list = CreateNeighList(
                &mut heap,
                count as i32,
                count as i32,
                atoms_pointer.as_const(),
                0,
                None,
            )
            .unwrap();
            let carb_pointer = heap.allocate_model_storage(carb_entries).unwrap();
            let mut stats = CANON_STAT::default();
            stats.LinearCTStereoCarb = carb_pointer;
            stats.nLenLinearCTStereoCarb =
                heap.slice(carb_pointer.as_const()).unwrap().len() as i32;
            stats.nSymmRank = symmetry_pointer;
            let baseline = heap.live_allocation_count();
            let result = RemoveCalculatedNonStereo(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                atoms_pointer,
                count as i32,
                count as i32,
                rank_stack1,
                rank_stack2,
                temporary_rank,
                neighbor_list,
                symmetry_pointer,
                canonical_pointer,
                atom_by_canonical_pointer,
                &mut stats,
                AB_PARITY_UNKN as i32,
            )
            .unwrap();
            let after = heap.live_allocation_count();
            (
                result,
                heap.slice(atoms_pointer.as_const()).unwrap().to_vec(),
                stats.nLenLinearCTStereoCarb,
                baseline,
                after,
            )
        }

        let unchanged = run(vec![sp_ATOM::default()], vec![1], vec![1], vec![0], vec![]);
        assert_eq!(unchanged.0, 0);
        assert_eq!(unchanged.3, unchanged.4);

        let removed = run(
            candidate(),
            vec![1, 2, 2, 4],
            vec![1, 2, 3, 4],
            vec![0, 1, 2, 3],
            vec![AT_STEREO_CARB {
                at_num: 1,
                parity: AB_PARITY_CALC as u8,
            }],
        );
        assert_eq!(removed.0, 1);
        assert_eq!(removed.1[0].parity, 0);
        assert_eq!(removed.1[0].stereo_atom_parity, 0);
        assert_eq!(removed.1[0].final_parity, 0);
        assert_eq!(removed.2, 0);
        assert_eq!(removed.3, removed.4);

        let missing = run(
            candidate(),
            vec![1, 2, 2, 4],
            vec![1, 2, 3, 4],
            vec![0, 1, 2, 3],
            vec![AT_STEREO_CARB {
                at_num: 4,
                parity: AB_PARITY_ODD as u8,
            }],
        );
        assert_eq!(missing.0, CT_STEREOCOUNT_ERR);
        assert_eq!(missing.1[0].parity, 0);
        assert_eq!(missing.2, 1);
        assert_eq!(missing.3, missing.4);

        let mut failure_heap = SourceHeap::default();
        let atoms = failure_heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        let symmetry = failure_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let canonical = failure_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let atom_by_canonical = failure_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let baseline = failure_heap.live_allocation_count();
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            RemoveCalculatedNonStereo(
                &mut failure_heap,
                &mut CANON_GLOBALS::default(),
                atoms,
                1,
                1,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                symmetry,
                canonical,
                atom_by_canonical,
                &mut CANON_STAT::default(),
                AB_PARITY_UNKN as i32,
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert_eq!(failure_heap.live_allocation_count(), baseline);
    }
}
