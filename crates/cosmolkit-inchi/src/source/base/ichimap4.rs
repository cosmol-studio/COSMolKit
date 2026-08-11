use crate::source::base::ichicano::bInchiTimeIsOver;
use crate::source::base::ichimap1::{
    All_SB_Same, All_SC_Same, CompareLinCtStereoAtomToValues, CompareLinCtStereoDoubleToValues,
    CurTreeAddAtom, CurTreeAddRank, CurTreeGetPos, CurTreeIsLastAtomEqu, CurTreeIsLastRank,
    CurTreeKeepLastAtomsOnly, CurTreeRemoveIfLastAtom, CurTreeRemoveLastRank,
    CurTreeRemoveLastRankIfNoAtoms, CurTreeSetPos, Next_SB_At_CanonRanks2, Next_SC_At_CanonRank2,
    NextStereoParity2Test, SetUseAtomForStereo, bUniqueAtNbrFromMappingRank, nGetMcr, nJoin2Mcrs,
};
use crate::source::base::ichimap2::{
    ClearPreviousMappings, RemoveCalculatedNonStereo, map_an_atom2, might_change_other_atom_parity,
    parity_of_mapped_atom2, parity_of_mapped_half_bond,
};
use crate::source::base::ichisort::BreakAllTies;
use crate::source_types::{
    AB_MAX_KNOWN_PARITY, AB_MIN_KNOWN_PARITY, AB_PARITY_CALC, AB_PARITY_NONE, AB_PARITY_UNDF,
    AB_PARITY_UNKN, AT_RANK, AT_STEREO_CARB, AT_STEREO_DBLE, BEST_PARITY, BITS_PARITY,
    CANON_GLOBALS, CANON_STAT, CMODE_REDNDNT_STEREO, CT_ERR_MAX, CT_ERR_MIN, CT_MAPCOUNT_ERR,
    CT_OVERFLOW, CT_STEREO_CANON_ERR, CT_STEREOBOND_ERROR, CT_STEREOCOUNT_ERR, CT_TIMEOUT_ERR,
    CT_USER_QUIT_ERR, CUR_TREE, EQ_NEIGH, INCHI_CLOCK, KNOWN_PARITIES_EQL, MAX_NUM_STEREO_BONDS,
    NEIGH_LIST, STEREO_AT_MARK, SourceHeap, SourceHeapError, SourceMutPointer, WORSE_PARITY,
    inchiTime, ppAT_RANK, sp_ATOM,
};

#[inline(always)]
fn source_get<T: Clone + 'static>(
    heap: &SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
) -> Result<T, SourceHeapError> {
    heap.slice(pointer.as_const())?
        .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

#[inline(always)]
fn source_set<T: Clone + 'static>(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
    value: T,
) -> Result<(), SourceHeapError> {
    *heap
        .slice_mut(pointer)?
        .get_mut(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
    Ok(())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn map_stereo_bonds4(
    heap: &mut SourceHeap,
    ic: &mut INCHI_CLOCK,
    clock_result: i64,
    user_action: Option<fn() -> i32>,
    console_quit: Option<fn() -> i32>,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    num_at_tg: i32,
    num_max: i32,
    mut bAllene: i32,
    nCanonRankFrom: SourceMutPointer<AT_RANK>,
    nAtomNumberCanonFrom: SourceMutPointer<AT_RANK>,
    nCanonRankTo: SourceMutPointer<AT_RANK>,
    nSymmRank: SourceMutPointer<AT_RANK>,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nTempRank: SourceMutPointer<AT_RANK>,
    nNumMappedRanksInput: i32,
    nSymmStereo: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    pCS: &mut CANON_STAT,
    mut cur_tree: Option<&mut CUR_TREE>,
    nNumMappedBonds: i32,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c:83 map_stereo_bonds4
    // INCHI✔️❌: complete active source frame follows verbatim.
    /*
    int map_stereo_bonds4( struct tagINCHI_CLOCK *ic,
                            CANON_GLOBALS *pCG,
                            sp_ATOM *at,
                            int num_atoms,
                            int num_at_tg,
                            int num_max,
                            int bAllene,
                            const AT_RANK *nCanonRankFrom,
                            const AT_RANK *nAtomNumberCanonFrom, /*  non-stereo canon ranking */
                            AT_RANK *nCanonRankTo,  /* output canonical stereo numbering*/
                            const AT_RANK *nSymmRank,
                            AT_RANK   **pRankStack1 /* from */,
                            AT_RANK **pRankStack2   /* to */,
                            AT_RANK *nTempRank,
                            int nNumMappedRanksInput,
                            AT_RANK *nSymmStereo,
                            NEIGH_LIST *NeighList,
                            CANON_STAT *pCS,
                            CUR_TREE *cur_tree,
                            int nNumMappedBonds,
                            int vABParityUnknown )
    {
        int nTotSuccess = 0; /* 1=>full mapping has been completed;
                              * 2=>obtained a better stereo;
                              * 4=>restart (stereo bond or atom removed from the stereo CT)
                              */
        int tpos1 = 0;
        AT_STEREO_DBLE prevBond;
        memset( &prevBond, 0, sizeof( prevBond ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        tpos1 = CurTreeGetPos( cur_tree );

    total_restart:

        if (!nNumMappedBonds)
        {

            memset( pCS->bRankUsedForStereo, 0, sizeof( pCS->bRankUsedForStereo[0] )*num_atoms ); /* djb-rwth: memset_s C11/Annex K variant? */
            SetUseAtomForStereo( pCS->bAtomUsedForStereo, at, num_atoms );

            if (pCS->bFirstCT && nSymmStereo && !pCS->bKeepSymmRank)
            {
                int i;
                for (i = 0; i < num_at_tg; i++)
                {
                    /*  nSymmStereo[i] = min. {k | at[k] stereo eq. to at[i]} */
                    nSymmStereo[i] = i; /*  for union-join to keep track of stereo-equivalent atoms */
                }
            }
    #ifdef FIX_STEREOCOUNT_ERR
            if (pCS->bFirstCT)
            {
                CurTreeSetPos( cur_tree, tpos1 = 0 );
            }
    #endif
        }

        if (nNumMappedBonds < pCS->nLenLinearCTStereoDble)
        {

            int at_rank1, at_rank2, bStereoIsBetterWasSetHere;
            /* AT_RANK *nRankFrom=*pRankStack1++,  AT_RANK *nAtomNumberFrom=pRankStack1++; */
            /* AT_RANK *nRankTo  =*pRankStack2++,  AT_RANK *nAtomNumberTo  =pRankStack2++; */
            AT_RANK canon_min1, canon_min2;
            int bFirstCanonRank;
            int i, j, j1, j2, at_from1, at_from2, at_to1, at_to2, iMax, c;
            int nStackPtr[SB_DEPTH], nNumMappedRanks[SB_DEPTH], LastMappedTo1;
            int istk, istk2, istk3, bAddStack, nNumAtTo1Success;
            int ret1, ret2, parity1, parity2;

            AT_RANK at_rank_canon1; /*  = pCS->LinearCTStereoDble[nNumMappedBonds].at_num1; */ /*  canonical numbers of atoms */
            AT_RANK at_rank_canon2; /*  = pCS->LinearCTStereoDble[nNumMappedBonds].at_num2; */ /*  adjacent to the stereogenic bond */
            int nNumChoices, nNumUnkn, nNumUndf, nNumBest, nNumWorse, nNumCalc, sb_parity_calc;
            int stereo_bond_parity = 0, prev_stereo_bond_parity, pass, bAllParitiesIdentical, bAllParitiesIdentical2;
            AT_STEREO_DBLE prevBond2;

            prevBond = pCS->LinearCTStereoDble[nNumMappedBonds];
            bFirstCanonRank = 1;
            canon_min1 = canon_min2 = 0;
            /*
                    // find candidates for atom_from1, atom_to1; they must have identical mapping ranks
                    at_rank1=pRankStack1[0][at_from1=nAtomNumberCanonFrom[(int)at_rank_canon1 - 1]]; // rank "from" for mapping
                    at_rank2=pRankStack1[0][at_from2=nAtomNumberCanonFrom[(int)at_rank_canon2 - 1]]; // rank "from" for mapping
            */
            if (nNumMappedBonds)
            {
                at_rank_canon1 = pCS->LinearCTStereoDble[nNumMappedBonds - 1].at_num1;
                at_rank_canon2 = pCS->LinearCTStereoDble[nNumMappedBonds - 1].at_num2;
            }
            else
            {
                at_rank_canon1 = 0;
                at_rank_canon2 = 0;
            }
            goto bypass_next_canon_ranks_check;

        next_canon_ranks:

                /*  Save time: avoid calling Next_SB_At_CanonRanks2() */
            if ((!pCS->bStereoIsBetter /* ??? && !pCS->bFirstCT ???*/ &&
                  at_rank_canon1 > pCS->LinearCTStereoDble[nNumMappedBonds].at_num1) ||
                  (at_rank_canon1 == pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 &&
                  at_rank_canon2 >= pCS->LinearCTStereoDble[nNumMappedBonds].at_num2)) /* djb-rwth: addressing LLVM warning */
            {

                if (!nTotSuccess)
                {
                    pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
                }
                CurTreeSetPos( cur_tree, tpos1 );
                return nTotSuccess;
            }

        bypass_next_canon_ranks_check:

            CurTreeSetPos( cur_tree, tpos1 );

            /*  find next available canon. numbers for a stereogenic bond pair of atoms */
            /*  process allenes AFTER all double bonds and odd-number-of-double-bonds cumulenes */
            if (!( ret1 = Next_SB_At_CanonRanks2( &at_rank_canon1, &at_rank_canon2, /*  canonical numbers */
                &canon_min1, &canon_min2,
                &bFirstCanonRank, pCS->bAtomUsedForStereo,
                pRankStack1, pRankStack2,
                nCanonRankFrom, nAtomNumberCanonFrom,
                at, num_atoms, bAllene ) ))
            {
                /* failed to find next stereo bond to assign parity */
                if (!bAllene && bFirstCanonRank)
                {
                    /* all stereobond have been processed; try to find allene to continue */
                    AT_RANK at_rank_canon1_Allene = 0, canon_min1_Allene = 0;
                    AT_RANK at_rank_canon2_Allene = 0, canon_min2_Allene = 0;
                    if ((ret1 = Next_SB_At_CanonRanks2( &at_rank_canon1_Allene, &at_rank_canon2_Allene,
                        &canon_min1_Allene, &canon_min2_Allene,
                        &bFirstCanonRank, pCS->bAtomUsedForStereo,
                        pRankStack1, pRankStack2,
                        nCanonRankFrom, nAtomNumberCanonFrom,
                        at, num_atoms, 1 ))) /* djb-rwth: addressing LLVM warning */
                    {
                        at_rank_canon1 = at_rank_canon1_Allene;
                        at_rank_canon2 = at_rank_canon2_Allene;
                        canon_min1 = canon_min1_Allene;
                        canon_min2 = canon_min2_Allene;
                        bAllene = 1; /* switch to allenes */
                    }
                }
            }

            if (!ret1 || (!pCS->bStereoIsBetter &&
                ( at_rank_canon1 > pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 ||
                    (at_rank_canon1 == pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 &&
                    at_rank_canon2 > pCS->LinearCTStereoDble[nNumMappedBonds].at_num2) ))) /* djb-rwth: addressing LLVM warnings */
            {
                /* new ranks provide greater pCS->LinearCTStereoDble[nNumMappedBonds] and therefore rejected */
                if (!nTotSuccess)
                {
                    pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond; /* restore stereo bond CT for the current bond */
                }
                return nTotSuccess;
            }

            /* current stereo bond initialization */
            nNumChoices = 0;
            nNumUnkn = 0;
            nNumUndf = 0;
            nNumBest = 0;
            nNumWorse = 0;
            nNumCalc = 0;
            pass = 0;
            prev_stereo_bond_parity = 0;

            at_rank1 = pRankStack1[0][at_from1 = nAtomNumberCanonFrom[(int) at_rank_canon1 - 1]]; /* atom 1 rank "from" for mapping */
            at_rank2 = pRankStack1[0][at_from2 = nAtomNumberCanonFrom[(int) at_rank_canon2 - 1]]; /* atom 2 rank "from" for mapping */
            /* we are going to map bond (at[at_from1], at[at_from2]) and
               canonical ranks of its atoms (at_rank_canon1, at_rank_canon2)
               onto a stereogenic bond (at[at_to1], at[at_to2])
             */
            iMax = at_rank1 - 1;
            /*  test correctness: sorted pRankStack2[0][] and pRankStack1[0][] should have same ranks for both atoms */
            if (at_rank1 != pRankStack2[0][pRankStack2[1][at_rank1 - 1]] ||
                 at_rank2 != pRankStack2[0][pRankStack2[1][at_rank2 - 1]])
            {
                /* program error: "from" and "to" mapping ranks are not equal */
                return CT_STEREOCOUNT_ERR; /*   <BRKPT> */
            }
            /* -- do not check stereo features of "from" atoms:
               -- in case of "bond/charge isomerism" they may be missing.
            if ( !at[at_from1].stereo_bond_neighbor[0] ||
                 !at[at_from2].stereo_bond_neighbor[0] )
                return CT_STEREOCOUNT_ERR; // program error
            */

            /*  find out if we have a choice in mapping: check all possible pairs (at_to1, at_to2)
                such that at_from1 is possibly constitutionally equivalent to at_to1, at_from2 to at_to2 */
            for (j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1 = pRankStack2[1][iMax - j1]]; j1++)
            {
                if (!at[at_to1].stereo_bond_neighbor[0])
                    continue; /*  at_to1 does not belong to a stereo bond */
                for (j2 = 0; j2 < MAX_NUM_STEREO_BONDS &&
                    ( at_to2 = at[at_to1].stereo_bond_neighbor[j2] ); j2++)
                {
                    at_to2--;
                    if (pRankStack1[0][at_from2] != pRankStack2[0][at_to2])
                        continue; /*  at_from2 cannot be mapped on at_to2 */
                    stereo_bond_parity = PARITY_VAL( at[at_to1].stereo_bond_parity[j2] );
                    i = 0;
                    switch (stereo_bond_parity)
                    {

                        case AB_PARITY_UNDF: nNumUndf++;
                            break; /*  4 */
                        case AB_PARITY_UNKN: nNumUnkn++;
                            break; /*  3 (occurs if forced different to UNDF)*/

                        case BEST_PARITY:    nNumBest++; break; /*  1 */
                        case WORSE_PARITY:   nNumWorse++; break; /*  2 */
                        case AB_PARITY_CALC: nNumCalc++; break; /*  6 */
                        case AB_PARITY_NONE: i++;         break; /*  0 */
                    }
                    nNumChoices += !i;
                }
            }
            if (nNumChoices != nNumCalc + nNumUndf + nNumUnkn + nNumBest + nNumWorse)
            {
                return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
            }
            if (!nNumChoices)
            {
                goto next_canon_ranks;
            }
            /*  Determine the first parity to search */
            sb_parity_calc = ( nNumCalc > 0 ) ? BEST_PARITY : 0;

            /*  ==============================================================
                Search sequence:           sb_parity_calc    stereo_bond_parity
                ==============================================================
                BEST_PARITY   (calc)       BEST_PARITY     BEST_PARITY
                BEST_PARITY   (known)      BEST_PARITY     WORSE_PARITY  or 0
                WORSE_PARITY  (calc)       WORSE_PARITY    WORSE_PARITY
                WORSE_PARITY  (known)      WORSE_PARITY    0
                AB_PARITY_UNKN(known)      AB_PARITY_UNKN  0
                AB_PARITY_UNDF(known)      AB_PARITY_UNDF  0

                if (sb_parity_calc==stereo_bond_parity) then "calc" else "known"
             */

        repeat_all:

            if (!nNumMappedBonds)
                pCS->bStereoIsBetter = 0;  /*  the first stereo feature in the canonical CT; moved here 7-13-2002 */

            if (!pass++)
            {
                /*  select the smallest (best) parity to search */
                if (sb_parity_calc)
                {
                    stereo_bond_parity = BEST_PARITY;
                }
                else
                {
                    stereo_bond_parity = nNumBest ? BEST_PARITY :
                        nNumWorse ? WORSE_PARITY :
                        nNumUnkn ? AB_PARITY_UNKN :
                        nNumUndf ? AB_PARITY_UNDF : AB_PARITY_NONE;
                }
            }
            else
            {
                 /* second pass: since the first pass failed, search for a worse result */
                prev_stereo_bond_parity = stereo_bond_parity;
                i = NextStereoParity2Test( &stereo_bond_parity, &sb_parity_calc,
                                         nNumBest, nNumWorse, nNumUnkn, nNumUndf, nNumCalc, vABParityUnknown );
                switch (i)
                {
                    case 0:
                        break; /* obtained next parity to test */
                    case 1:
                        goto next_canon_ranks;
                    default:
                        return i; /* program error */
                }
            }
            if (stereo_bond_parity == AB_PARITY_NONE)
            {
                /*  error? */
                return CT_STEREOCOUNT_ERR;                   /*   <BRKPT> */
            }
            /*  check if the new requested parity is good (small) enough */
            if (!pCS->bStereoIsBetter)
            {
                c = CompareLinCtStereoDoubleToValues( nTotSuccess ? pCS->LinearCTStereoDble + nNumMappedBonds : &prevBond,
                                  at_rank_canon1, at_rank_canon2, (U_CHAR) stereo_bond_parity );
                if (c < 0)
                {
                    if (!nTotSuccess)
                    {
                        pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
                    }
                    CurTreeSetPos( cur_tree, tpos1 );
                    return nTotSuccess;
                }
            }

            bAllParitiesIdentical = 0;
            /* djb-rwth: removing redundant code */
            LastMappedTo1 = -1;
            bStereoIsBetterWasSetHere = 0;
            /* djb-rwth: removing redundant code */

            if (!nNumMappedBonds && prev_stereo_bond_parity != stereo_bond_parity)
                pCS->bStereoIsBetter = 0;  /*  the first stereo feature in the canonical CT; moved here 5-24-2002 */

            if (prev_stereo_bond_parity != stereo_bond_parity)
            {
                CurTreeSetPos( cur_tree, tpos1 );  /*  start over */
            }

            /* Mapping: here at_rank1 = nRankTo, at_to1 = nAtomNumberTo */
            for (j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1 = pRankStack2[1][iMax - j1]]; j1++)
            {
                nNumAtTo1Success = 0;
                if (!at[at_to1].stereo_bond_neighbor[0])
                    continue; /*  at_to1 does not belong to a stereo bond */
                if (tpos1 < CurTreeGetPos( cur_tree ) &&
                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) &&
                     1 == CurTreeIsLastAtomEqu( cur_tree, at_to1, nSymmStereo ))
                {
                    /* at_to1 is known to be stereogenically equivalent to another atom tried with at_rank_canon1 */
                    continue;
                }
                bAllParitiesIdentical2 = 0;
                for (j2 = 0; j2 < MAX_NUM_STEREO_BONDS && ( at_to2 = at[at_to1].stereo_bond_neighbor[j2] ); j2++)
                {
                    EQ_NEIGH  EN1[2], EN2[2];
                    int bond_parity, num1, num2;
                    AT_RANK at_rank_canon_n1, at_rank_canon_n2;

                    at_to2--;
                    if (pRankStack1[0][at_from2] != pRankStack2[0][at_to2])
                        continue; /*  at_from2 cannot be mapped on at_to2 even without mapping at_from1 to at_to1 */

                    /*  check whether the bond parity corresponds to the requested bond parity */
                    if (PARITY_KNOWN( at[at_to1].stereo_bond_parity[j2] ))
                    {
                        if (stereo_bond_parity == sb_parity_calc)
                        {
                            continue;  /*  requested parity to be calculated, found known parity */
                        }
                        if (stereo_bond_parity != PARITY_VAL( at[at_to1].stereo_bond_parity[j2] ))
                        {
                            continue;  /*  parity differs from the requested parity */
                        }
                    }
                    else
                        if (PARITY_CALCULATE( at[at_to1].stereo_bond_parity[j2] ))
                        {
                            if (stereo_bond_parity != sb_parity_calc)
                            {
                                continue;  /*  requested known parity, found parity to be calculated */
                            }
                        }
                        else
                        {
                            return CT_STEREOCOUNT_ERR;  /*  unknown parity type */ /*   <BRKPT> */
                        }
                        /*  initialize stack pointer nStackPtr[istk] for "hand-made" recursion */
                        /*  stacks are pRankStack1[], pRankStack2[], nNumMappedRanks[] */
                    istk = 0;
                    nStackPtr[0] = 0;
                    nNumMappedRanks[0] = nNumMappedRanksInput;
                    bAddStack = 0;
                    bAllParitiesIdentical = ( ( at[at_to1].stereo_bond_parity[j2] & KNOWN_PARITIES_EQL ) ) &&
                        PARITY_KNOWN( at[at_to1].stereo_bond_parity[j2] );

                    if (!bAllParitiesIdentical && !nNumCalc &&
                        ( !nNumUndf + !nNumUnkn + !nNumBest + !nNumWorse ) == 3)
                    {
                        /* only one kind of bond parity is present; check whether all parities are really same */
                        bAllParitiesIdentical = All_SB_Same( at_rank_canon1, at_rank_canon2, /*  canonical numbers */
                                                              pRankStack1, pRankStack2,
                                                              nAtomNumberCanonFrom, at );
                        if (bAllParitiesIdentical < 0)
                        {
                            return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                        }
                    }

                    /*****************************************************************
                     * do the mapping only if parities are not same
                     */
                    if (!bAllParitiesIdentical)
                    {
                        /*  map atom 1 or reuse previous mapping */
                        if (LastMappedTo1 != at_to1)
                        {
                            /*  avoid repetitve mapping to the same first at_to1 using LastMappedTo1 variable */
                            /*  map atom 1 */
                            ret1 = map_an_atom2( pCG, num_at_tg, num_max, at_from1, at_to1,
                                                nTempRank, nNumMappedRanks[istk], &nNumMappedRanks[istk + 1], pCS,
                                                NeighList, pRankStack1 + nStackPtr[istk], pRankStack2 + nStackPtr[istk],
                                                &bAddStack );
                            if (RETURNED_ERROR( ret1 ))
                            {
                                return ret1; /*  error */
                            }
                            nStackPtr[istk + 1] = nStackPtr[istk] + bAddStack;
                            LastMappedTo1 = at_to1;
                            if (bAddStack)
                            {
                                if (tpos1 == CurTreeGetPos( cur_tree ) ||
                                     0 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
                                {
                                    CurTreeAddRank( cur_tree, at_rank_canon1 );
                                }
                                CurTreeAddAtom( cur_tree, at_to1 );
                            }
                        }
                        istk++; /*  = 1 */
                        /*  check if we can map atom 2 */
                        if (pRankStack1[nStackPtr[istk]][at_from2] != pRankStack2[nStackPtr[istk]][at_to2])
                        {
                            /*
                             * This may happen when:
                             * A) Charge/bond isomerism, for example cyclopentadiene(-), or
                             * B) possibly stereogenic bond in an alternating ring has heighbors
                             * in 2 symmetrically attached rings.
                             * Such an alternating bond cannot be mapped on possibly stereogenic bond
                             * that has neighbors belonging to 1 of the symmetrically attached rings only.
                             * For example:
                             *   A---B---C---D  If all atoms are Carbons then B, C, F, G are constitutionally
                             *  ||  ||  ||  ||  equivalent. However, bonds B-C, F-G are not equivalent to
                             *  ||  ||  ||  ||  B-F and C-G and cannot be mapped on them.
                             *   E---F---G---H  If at_from1=B, at_from2=F, at_to1=B, then at_from2 cannot be mapped on at_to2=C
                             *                  If at_from1=B, at_from2=F, at_to1=C, then at_from2 cannot be mapped on at_to2=B
                             *                  etc.
                             */
                            if (sb_parity_calc != stereo_bond_parity)
                            {
                                /* can be passed only once for each bond */
                                nNumChoices--;
                                nNumUndf -= ( stereo_bond_parity == AB_PARITY_UNDF );
                                nNumUnkn -= ( stereo_bond_parity == AB_PARITY_UNKN );
                                nNumBest -= ( stereo_bond_parity == BEST_PARITY );
                                nNumWorse -= ( stereo_bond_parity == WORSE_PARITY );
                                /* nNumCalc  = nNumChoices - (nNumUndf + nNumUnkn + nNumBest + nNumWorse); */
                            }
                            else
                                if (sb_parity_calc == BEST_PARITY)
                                {
                                    /* can be passed 2 times: for BEST_PARITY and WORSE_PARITY in this order */
                                    nNumChoices--; /*  do not repeate for WORSE_PARITY */
                                    nNumCalc--;
                                }
                            continue;  /*  Happens for ID=80036,80253,91354,95532,101532,103788 */
                        }
                        if (nStackPtr[istk] > nStackPtr[istk - 1])
                        {
                            bAllParitiesIdentical2 = All_SB_Same( at_rank_canon1, at_rank_canon2,
                                                      pRankStack1 + nStackPtr[istk], pRankStack2 + nStackPtr[istk],
                                                      nAtomNumberCanonFrom, at );
                            if (bAllParitiesIdentical2 < 0)
                            {
                                return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                            }
                        }
                        else
                        {
                            bAllParitiesIdentical2 = 0;
                        }
                        if (bAllParitiesIdentical2)
                        {
                            /*  do no mapping when all equivalent bonds have same parity */
                            /*  stereo_bond_parity = PARITY_VAL(at[at_to1].stereo_bond_parity[j2]); */
                            ClearPreviousMappings( pRankStack1 + nStackPtr[istk] + 2 );
                        }
                        else
                        {
                            if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ) &&
                                 1 == CurTreeIsLastAtomEqu( cur_tree, at_to2, nSymmStereo ))
                            {
                                continue;
                            }
                            /*  map atom 2 */
                            ret2 = map_an_atom2( pCG, num_at_tg, num_max, at_from2, at_to2,
                                                nTempRank, nNumMappedRanks[istk], &nNumMappedRanks[istk + 1], pCS,
                                                NeighList, pRankStack1 + nStackPtr[istk], pRankStack2 + nStackPtr[istk],
                                                &bAddStack );
                            if (RETURNED_ERROR( ret2 ))
                            {
                                return ret2; /*  program error */
                            }
                            nStackPtr[istk + 1] = nStackPtr[istk] + bAddStack;
                            istk++; /*  = 2 */
                            if (bAddStack)
                            {
                                if (tpos1 == CurTreeGetPos( cur_tree ) ||
                                     0 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ))
                                {
                                    CurTreeAddRank( cur_tree, at_rank_canon2 );
                                }
                                CurTreeAddAtom( cur_tree, at_to2 );
                            }
                        }
                    }
                    else
                    {
                        /*  do no mapping when all equivalent bonds have same parity */
                        /*  stereo_bond_parity = PARITY_VAL(at[at_to1].stereo_bond_parity[j2]); */
                        ClearPreviousMappings( pRankStack1 + 2 );
                    }

                    /*  we have a precalculated (known) bond parity */

                    /************************************************************
                     *
                     *   Known Bond Parity case: do not map stereo bond neighbors
                     */
                    if (stereo_bond_parity != sb_parity_calc) /*  parity is known */
                    {
                        /*  accept bond parity and do not map the neighbors */
                        bond_parity = stereo_bond_parity;
                        /*  same code as under " make a decision to accept current mapping" comment below */
                        /*  with one exception: istk instead of istk3 */
                        c = CompareLinCtStereoDoubleToValues( pCS->LinearCTStereoDble + nNumMappedBonds,
                                                  at_rank_canon1, at_rank_canon2, (U_CHAR) bond_parity );
                        if (c < 0 && !pCS->bStereoIsBetter)
                        {
                            /*  reject */

                            pCS->lNumRejectedCT++;
                            /*  remove failed atom2 from the tree */
                            if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ))
                            {
                                CurTreeRemoveIfLastAtom( cur_tree, at_to2 );
                                CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                            }
                            continue;  /*  to next at_to2; Reject this at_to2: not a minimal CT. */
                        }
                        else
                        {
                            /*  accept */

                            if (c > 0 && !pCS->bStereoIsBetter)
                            {
                                /*  bond entry is less than the previusly found */
                                pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                                prevBond2 = pCS->LinearCTStereoDble[nNumMappedBonds];
                            }
                            pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 = at_rank_canon1;
                            pCS->LinearCTStereoDble[nNumMappedBonds].at_num2 = at_rank_canon2;
                            pCS->LinearCTStereoDble[nNumMappedBonds].parity = bond_parity;
                            /*  recursive call */
                            pCS->bRankUsedForStereo[at_from1] ++;
                            pCS->bRankUsedForStereo[at_from2] ++;
                            if (!bAllParitiesIdentical)
                            {
                                pCS->bAtomUsedForStereo[at_to1] --;
                                pCS->bAtomUsedForStereo[at_to2] --;
                            }
                            ret2 = map_stereo_bonds4( ic, pCG, at, num_atoms, num_at_tg, num_max, bAllene, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                                                      nSymmRank, pRankStack1 + nStackPtr[istk], pRankStack2 + nStackPtr[istk],
                                                      nTempRank, nNumMappedRanks[istk], nSymmStereo, NeighList,
                                                      pCS, cur_tree, nNumMappedBonds + 1,
                                                      vABParityUnknown );
                            if (!bAllParitiesIdentical)
                            {
                                pCS->bAtomUsedForStereo[at_to1] ++;
                                pCS->bAtomUsedForStereo[at_to2] ++;
                            }
                            pCS->bRankUsedForStereo[at_from1] --;
                            pCS->bRankUsedForStereo[at_from2] --;
                            if (ret2 == 4)
                            {
                                if (nNumMappedBonds)
                                {
                                    return ret2;
                                }
                                else
                                {
                                    pCS->bFirstCT = 1;
                                    goto total_restart;
                                }
                            }

                            if (RETURNED_ERROR( ret2 ))
                            {
                                if (ret2 == CT_TIMEOUT_ERR)
                                    return ret2;
                                else
                                    return ret2; /*  program error */
                            }
                            if (ret2 > 0)
                            {
                                nTotSuccess |= 1;
                                nNumAtTo1Success++;
                                if (bStereoIsBetterWasSetHere || ( ret2 & 2 ))
                                {
                                    CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                                    nTotSuccess |= 2; /*  Obtained a smaller CT */
                                }
                            }
                            else
                            {
                                if (bStereoIsBetterWasSetHere)
                                { /*  rollback */
                                    pCS->bStereoIsBetter = 0;
                                    pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond2;
                                }
                                /*  remove failed atom2 from the tree */
                                if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ))
                                {
                                    CurTreeRemoveIfLastAtom( cur_tree, at_to2 );
                                    CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                                }
                                /*
                                if ( 1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                                    CurTreeRemoveLastAtom( cur_tree, at_to1 );
                                    CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                                }
                                */
                            }
                            bStereoIsBetterWasSetHere = 0;
                        }
                        if (bAllParitiesIdentical || bAllParitiesIdentical2)
                        {
                            break; /* j2 cycle, at_to2 (no need to repeat) */
                        }
                        continue; /*  to next at_to2 */
                    }
                    /***************************************************************************
                     *
                     *   Unknown Bond Parity case: may need to map stereo bond neighbors
                     *
                     ****************************************************************************
                     * Ranks are not known in advance
                     * check if at_from1/at_to1 half-bond has neighbors with equal mapping ranks
                     */

                    parity1 = parity_of_mapped_half_bond( at_from1, at_to1, at_from2, at_to2, at, &EN1[0],
                               nCanonRankFrom, pRankStack1[nStackPtr[istk]], pRankStack2[nStackPtr[istk]] );
                    /* old approach -- before E/Z parities
                    parity1 = parity_of_mapped_atom2( pCG, at_from1, at_to1, at, &EN1[0],
                                   nCanonRankFrom, pRankStack1[nStackPtr[istk]], pRankStack2[nStackPtr[istk]] );
                     */
                    /*  the following commented out statement is not needed here. */
                    /*  parity2 = parity_of_mapped_atom2( pCG, at_from2, at_to2, at, &EN2[0],
                                                         nCanonRankFrom, pRankStack1[nStackPtr[istk]],
                                                         pRankStack2[nStackPtr[istk]] );
                     */
                    if (!parity1)
                    {
                        return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                    }
                    num1 = parity1 > 0 ? 1 : 2; /*  parity < 0 means additional mapping is needed to set parity */

                    /*  --- try all possible mappings of the stereo bond ending atoms' neighbors --- */
                    at_rank_canon_n1 = 0;
                    /* djb-rwth: removing redundant code */
                    for (i = 0; i < num1; i++)
                    {
                        int at_from_n1, at_to_n1, at_no_n1_num_success = 0;
                        istk2 = istk;
                        if (num1 == 2)
                        {
                            at_rank_canon_n1 = nCanonRankFrom[EN1[0].from_at];
                            /*  an additional neighbor mapping is necessary; */
                            /*  we need to map only one at_from1 neighbor to make all neighbors have different ranks */

                            at_from_n1 = EN1[0].from_at;
                            at_to_n1 = EN1[0].to_at[i];

                            if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ) &&
                                 1 == CurTreeIsLastAtomEqu( cur_tree, at_to_n1, nSymmStereo ))
                                continue;
                            /*
                            if ( nSymmStereo && !pCS->bFirstCT ) {
                                if ( i && nSymmStereo[at_to_n1] == nSymmStereo[(int)EN1[0].to_at[0]] ) {
                                    continue; // do not test stereo equivalent atoms except the first one
                                }
                            }
                            */
                            /*  neighbors are tied. Untie them by breaking a tie on ONE of them. */
                            ret1 = map_an_atom2( pCG, num_at_tg, num_max, at_from_n1, at_to_n1,
                                                nTempRank, nNumMappedRanks[istk2], &nNumMappedRanks[istk2 + 1], pCS,
                                                NeighList, pRankStack1 + nStackPtr[istk2], pRankStack2 + nStackPtr[istk2],
                                                &bAddStack );
                            if (RETURNED_ERROR( ret1 ))
                            {
                                return ret1; /*  program error */ /*   <BRKPT> */
                            }
                            nStackPtr[istk2 + 1] = nStackPtr[istk2] + bAddStack;
                            istk2++;  /*  <= 3 */
                            /*  debug */
                            if (istk2 >= SB_DEPTH)
                            {
                                return CT_OVERFLOW; /*  program error */ /*   <BRKPT> */
                            }
                            if (bAddStack)
                            {
                                if (tpos1 == CurTreeGetPos( cur_tree ) ||
                                     0 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ))
                                {
                                    CurTreeAddRank( cur_tree, at_rank_canon_n1 );
                                }
                                CurTreeAddAtom( cur_tree, at_to_n1 );
                            }


                            /*  now that all at_from1 neighbors have been mapped the parity must be defined */
                            parity1 = parity_of_mapped_half_bond( at_from1, at_to1, at_from2, at_to2, at, &EN1[1],
                               nCanonRankFrom, pRankStack1[nStackPtr[istk2]], pRankStack2[nStackPtr[istk2]] );
                            if (parity1 <= 0)
                                return CT_STEREOCOUNT_ERR;  /*  program error */ /*   <BRKPT> */
                        }
                        else
                        {
                            nNumMappedRanks[istk2 + 1] = nNumMappedRanks[istk2];
                            nStackPtr[istk2 + 1] = nStackPtr[istk2];
                            istk2++;  /*  <= 3 */
                        }

                        /*  check if at_from2/at_to2 half-bond has neighbors with equal mapping ranks */
                        parity2 = parity_of_mapped_half_bond( at_from2, at_to2, at_from1, at_to1, at, &EN2[0],
                               nCanonRankFrom, pRankStack1[nStackPtr[istk2]], pRankStack2[nStackPtr[istk2]] );
                        if (!parity2)
                        {
                            return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                        }
                        num2 = parity2 > 0 ? 1 : 2;
                        at_rank_canon_n2 = 0;
                        for (j = 0; j < num2; j++)
                        {
                            int at_from_n2, at_to_n2 = 0;
                            istk3 = istk2;
                            if (num2 == 2)
                            {
                                at_rank_canon_n2 = nCanonRankFrom[EN2[0].from_at];
                                /*  we need to map only one at_from2 neighbor to make its neighbors have different ranks */
                                at_from_n2 = EN2[0].from_at;
                                at_to_n2 = EN2[0].to_at[j];

                                if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ) &&
                                     1 == CurTreeIsLastAtomEqu( cur_tree, at_to_n2, nSymmStereo ))
                                    continue;

                                /*
                                if ( nSymmStereo && !pCS->bFirstCT ) {
                                    if ( j && nSymmStereo[at_to_n2] == nSymmStereo[(int)EN2[0].to_at[0]] ) {
                                        continue; // do not test stereo equivalent atoms except the first one
                                    }
                                }
                                */
                                /*  neighbors are tied. Untie them by breaking a tie on ONE of them. */
                                ret1 = map_an_atom2( pCG, num_at_tg, num_max, at_from_n2, at_to_n2,
                                                    nTempRank, nNumMappedRanks[istk3], &nNumMappedRanks[istk3 + 1], pCS,
                                                    NeighList, pRankStack1 + nStackPtr[istk3],
                                                    pRankStack2 + nStackPtr[istk3],
                                                    &bAddStack );
                                if (RETURNED_ERROR( ret1 ))
                                {
                                    return ret1; /*  program error */
                                }
                                nStackPtr[istk3 + 1] = nStackPtr[istk3] + bAddStack;
                                istk3++;  /*  <= 4 */

                                if (bAddStack)
                                {
                                    if (tpos1 == CurTreeGetPos( cur_tree ) ||
                                         0 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ))
                                    {
                                        CurTreeAddRank( cur_tree, at_rank_canon_n2 );
                                    }
                                    CurTreeAddAtom( cur_tree, at_to_n2 );
                                }

                                parity2 = parity_of_mapped_half_bond( at_from2, at_to2, at_from1, at_to1, at, &EN2[1],
                                         nCanonRankFrom, pRankStack1[nStackPtr[istk3]], pRankStack2[nStackPtr[istk3]] );
                                if (parity2 <= 0)
                                {
                                    return CT_STEREOCOUNT_ERR;  /*  program error */ /*   <BRKPT> */
                                }
                            }
                            else
                            {
                                                 /*  no additional mapping is needed to set atom's parity */
                                nNumMappedRanks[istk3 + 1] = nNumMappedRanks[istk3];
                                nStackPtr[istk3 + 1] = nStackPtr[istk3];
                                istk3++;  /*  <= 4 */
                            }

                            /*******************************************************************
                             * at this point the stereo bond is fully mapped to find its parity
                             *******************************************************************/

                            if (parity1 <= 0 || parity2 <= 0)
                            {
                                return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                            }

                            /*  find current bond parity  AB_PARITY_ODD */
                            if (ATOM_PARITY_WELL_DEF( parity1 ) && ATOM_PARITY_WELL_DEF( parity2 ))
                            {
                                bond_parity = 2 - ( parity1 + parity2 ) % 2;
                            }
                            else
                            {
                                bond_parity = inchi_max( parity1, parity2 );
                            }
                            if (ATOM_PARITY_WELL_DEF( bond_parity ) && at[at_to1].stereo_bond_z_prod[j2] < 0)
                                bond_parity = 2 - ( bond_parity + 1 ) % 2; /*  invert the bond parity */


                            /********************************************************
                             * make a decision whether to accept the current mapping
                             */
                            c = CompareLinCtStereoDoubleToValues( pCS->LinearCTStereoDble + nNumMappedBonds,
                                                      at_rank_canon1, at_rank_canon2, (U_CHAR) bond_parity );
                            if (sb_parity_calc != bond_parity ||
                                 (c < 0 && !pCS->bStereoIsBetter)) /* djb-rwth: addressing LLVM warning */
                            {
                                /*  reject */
                                pCS->lNumRejectedCT++;
                                /*  remove failed atom2 from the tree */
                                if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ))
                                {
                                    CurTreeRemoveIfLastAtom( cur_tree, at_to_n2 );
                                    CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                                }
                                continue;  /*  Reject: not a minimal CT. */
                            }
                            else
                            {

                                                /*  try to accept */

                                if (c > 0 && !pCS->bStereoIsBetter)
                                {
                                    /*  bond_parity is less than the previusly found */
                                    pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                                    prevBond2 = pCS->LinearCTStereoDble[nNumMappedBonds];
                                }
                                /*  accept */
                                pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 = at_rank_canon1;
                                pCS->LinearCTStereoDble[nNumMappedBonds].at_num2 = at_rank_canon2;
                                pCS->LinearCTStereoDble[nNumMappedBonds].parity = bond_parity;
                                /*  recursive call */
                                pCS->bRankUsedForStereo[at_from1] ++;
                                pCS->bRankUsedForStereo[at_from2] ++;
                                pCS->bAtomUsedForStereo[at_to1] --;
                                pCS->bAtomUsedForStereo[at_to2] --;
                                ret2 = map_stereo_bonds4( ic, pCG, at, num_atoms, num_at_tg, num_max, bAllene, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                                                          nSymmRank, pRankStack1 + nStackPtr[istk3], pRankStack2 + nStackPtr[istk3],
                                                          nTempRank, nNumMappedRanks[istk3], nSymmStereo, NeighList,
                                                          pCS, cur_tree, nNumMappedBonds + 1,
                                                          vABParityUnknown );
                                pCS->bRankUsedForStereo[at_from1] --;
                                pCS->bRankUsedForStereo[at_from2] --;
                                pCS->bAtomUsedForStereo[at_to1] ++;
                                pCS->bAtomUsedForStereo[at_to2] ++;
                                if (ret2 == 4)
                                {
                                    if (nNumMappedBonds)
                                    {
                                        return ret2;
                                    }
                                    else
                                    {
                                        pCS->bFirstCT = 1;
                                        goto total_restart;
                                    }
                                }
                                if (RETURNED_ERROR( ret2 ))
                                {
                                    return ret2; /*  program error */ /* djb-rwth: fixing coverity ID #499569 */
                                }
                                if (ret2 > 0)
                                {
                                    nTotSuccess |= 1;
                                    nNumAtTo1Success++;
                                    if (bStereoIsBetterWasSetHere || ( ret2 & 2 ))
                                    {
                                        CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                                        nTotSuccess |= 2; /*  Obtained a smaller CT */
                                    }
                                    at_no_n1_num_success++;
                                }
                                else
                                {
                                    if (bStereoIsBetterWasSetHere)
                                    {  /*  rollback */
                                        pCS->bStereoIsBetter = 0;
                                        pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond2;
                                    }
                                    if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                         1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ))
                                    {
                                        CurTreeRemoveIfLastAtom( cur_tree, at_to_n2 );
                                        CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                                    }
                                }
                                bStereoIsBetterWasSetHere = 0;
                            }
                        } /*  end choices in mapping neighbors of the 2nd half-bond */
                        if (tpos1 < CurTreeGetPos( cur_tree ) &&
                             1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ))
                        {
                            CurTreeRemoveLastRank( cur_tree );
                        }
                        /* added 2006-07-20 */
                        if (!at_no_n1_num_success && tpos1 < CurTreeGetPos( cur_tree ) &&
                            1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ))
                        {
                            CurTreeRemoveIfLastAtom( cur_tree, at_to_n1 );
                        }
                    } /*  end choices in mapping neighbors of the 1st half-bond */
                    if (tpos1 < CurTreeGetPos( cur_tree ) &&
                          1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ))
                    {
                        CurTreeRemoveLastRank( cur_tree );
                    }
                } /*  end of choices in mapping at_from2 */
                if (tpos1 < CurTreeGetPos( cur_tree ) &&
                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ))
                {
                    CurTreeRemoveLastRank( cur_tree );
                }
                if (!nNumAtTo1Success)
                {
                    if (tpos1 < CurTreeGetPos( cur_tree ) &&
                         1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
                    {
                        CurTreeRemoveIfLastAtom( cur_tree, at_to1 );
                        CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                    }
                }
                if (bAllParitiesIdentical /*&& !nSymmStereo*/)
                {
                    break;
                }
            } /*  end of choices in mapping at_from1 */

            if (tpos1 < CurTreeGetPos( cur_tree ) &&
                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
            {
                CurTreeRemoveLastRank( cur_tree );
            }
            else
            {
                /*  CurTree consistecy check (debug only) */
                if (tpos1 != CurTreeGetPos( cur_tree ))
                {
                    return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                }
            }

            if (!nTotSuccess || stereo_bond_parity == sb_parity_calc)
            {
                goto repeat_all; /*  repeat with next parity if no success or with the same parity, now known */
            }

            /*  Previously the control flow never came here... */
            if (!nTotSuccess)
            {
                pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
                CurTreeSetPos( cur_tree, tpos1 );
                /*  Occurs when atoms are not really equvalent ( -O= without positive charge in "aromatic" ring) */
                return 0; /* Happens for ID=92439,100318,100319 when EXCL_ALL_AROM_BOND_PARITY=0 and
                           * nNumChoices=0.
                           * Results from impossible previous mapping of symmetric relatively
                           * to a central ring aromatic circles while central ring is not symmetrical due to
                           * alternate bonds (in the central ring number of pi-electrons, atoms and bonds
                           * are symmetrical).
                           * Does not happen when alternate bonds of the central ring
                           * are treated as aromatic by attaching a (+) charge to the oxygen.
                           */
            }
        }
        else
        {
            int ret;

            if (!nNumMappedBonds)
            {
                pCS->bStereoIsBetter = 0;  /*  the first stereo feature in the canonical CT has not been processed yet */
            }

            if (nNumMappedBonds < pCS->nLenLinearCTStereoDble)
            {
                prevBond = pCS->LinearCTStereoDble[nNumMappedBonds];
            }

            /*deep_map_stereo_atoms4=0;*/

            /*  all stereo bonds have been mapped; now start processing stereo atoms... */
            ret = map_stereo_atoms4( ic, pCG, at, num_atoms, num_at_tg, num_max, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                            nSymmRank, pRankStack1, pRankStack2, nTempRank, nNumMappedRanksInput,
                            nSymmStereo, NeighList, pCS, cur_tree, 0, vABParityUnknown );

            if (ret == 4)
            {
                if (nNumMappedBonds)
                {
                    return ret;
                }
                else
                {
                    pCS->bFirstCT = 1;
                    goto total_restart;
                }
            }
            if (RETURNED_ERROR( ret ))
            {
                if (ret == CT_TIMEOUT_ERR)
                    return ret;
                else
                    return ret; /*  program error */
            }
            if (ret > 0)
            {
                nTotSuccess |= 1;
                if (ret & 2)
                {
                    CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                    nTotSuccess |= 2; /*  Obtained a smaller CT */
                }
            }
        }
        if (!nTotSuccess && pCS->nLenLinearCTStereoDble &&
             nNumMappedBonds < pCS->nLenLinearCTStereoDble)
        {
            pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
        }

        return nTotSuccess;  /*  ok */
    }
    */
    // END INCHI C FUNCTION: map_stereo_bonds4
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION
    // INCHI✔️❌: FIX_STEREOCOUNT_ERR=1 for the selected production target.
    // INCHI✔️❌: bRELEASE_VERSION=1; the inactive debug-only branches are excluded.
    // END INCHI ACTIVE MACRO CONFIGURATION

    let returned_error = |value: i32| (CT_ERR_MIN..=CT_ERR_MAX).contains(&value);
    let parity_value = |value: i8| i32::from(value) & BITS_PARITY as i32;
    let parity_known = |value: i8| {
        (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&parity_value(value))
    };
    let parity_calculate = |value: i8| parity_value(value) == AB_PARITY_CALC as i32;
    let atom_parity_well_defined =
        |value: i32| (AB_MIN_KNOWN_PARITY as i32..=WORSE_PARITY as i32).contains(&value);
    let mut tpos1 = CurTreeGetPos(cur_tree.as_deref());

    'total_restart: loop {
        let mut total_success = 0_i32;
        let mut previous_bond = AT_STEREO_DBLE::default();

        if nNumMappedBonds == 0 {
            let count =
                usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(pCS.bRankUsedForStereo)?
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .fill(0);
            SetUseAtomForStereo(heap, pCS.bAtomUsedForStereo, at.as_const(), num_atoms)?;
            if pCS.bFirstCT != 0 && !nSymmStereo.is_null() && pCS.bKeepSymmRank == 0 {
                let mut index = 0_i32;
                while index < num_at_tg {
                    source_set(heap, nSymmStereo, index, index as AT_RANK)?;
                    index = index.wrapping_add(1);
                }
            }
            if pCS.bFirstCT != 0 {
                tpos1 = 0;
                CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
            }
        }

        if nNumMappedBonds >= pCS.nLenLinearCTStereoDble {
            if nNumMappedBonds == 0 {
                pCS.bStereoIsBetter = 0;
            }
            let result = map_stereo_atoms4(
                heap,
                ic,
                clock_result,
                user_action,
                console_quit,
                pCG,
                at,
                num_atoms,
                num_at_tg,
                num_max,
                nCanonRankFrom,
                nAtomNumberCanonFrom,
                nCanonRankTo,
                nSymmRank,
                pRankStack1,
                pRankStack2,
                nTempRank,
                nNumMappedRanksInput,
                nSymmStereo,
                NeighList,
                pCS,
                cur_tree.as_deref_mut(),
                0,
                vABParityUnknown,
            )?;
            if result == 4 {
                if nNumMappedBonds != 0 {
                    return Ok(result);
                }
                pCS.bFirstCT = 1;
                continue 'total_restart;
            }
            if returned_error(result) {
                return Ok(result);
            }
            if result > 0 {
                total_success |= 1;
                if result & 2 != 0 {
                    CurTreeKeepLastAtomsOnly(heap, cur_tree.as_deref_mut(), tpos1, 1)?;
                    total_success |= 2;
                }
            }
            return Ok(total_success);
        }

        previous_bond = source_get(heap, pCS.LinearCTStereoDble, nNumMappedBonds)?;
        let mut canon_rank1 = if nNumMappedBonds != 0 {
            source_get(
                heap,
                pCS.LinearCTStereoDble,
                nNumMappedBonds.wrapping_sub(1),
            )?
            .at_num1
        } else {
            0
        };
        let mut canon_rank2 = if nNumMappedBonds != 0 {
            source_get(
                heap,
                pCS.LinearCTStereoDble,
                nNumMappedBonds.wrapping_sub(1),
            )?
            .at_num2
        } else {
            0
        };
        let mut canon_min1 = 0_u16;
        let mut canon_min2 = 0_u16;
        let mut first_canon_rank = 1_i32;
        let mut bypass_limit = true;

        'canon_ranks: loop {
            if !bypass_limit {
                let current = source_get(heap, pCS.LinearCTStereoDble, nNumMappedBonds)?;
                if (pCS.bStereoIsBetter == 0 && canon_rank1 > current.at_num1)
                    || (canon_rank1 == current.at_num1 && canon_rank2 >= current.at_num2)
                {
                    if total_success == 0 {
                        source_set(heap, pCS.LinearCTStereoDble, nNumMappedBonds, previous_bond)?;
                    }
                    CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
                    return Ok(total_success);
                }
            }
            bypass_limit = false;
            CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);

            let mut found = Next_SB_At_CanonRanks2(
                heap,
                &mut canon_rank1,
                &mut canon_rank2,
                &mut canon_min1,
                &mut canon_min2,
                &mut first_canon_rank,
                pCS.bAtomUsedForStereo.as_const(),
                pRankStack1,
                pRankStack2,
                nCanonRankFrom.as_const(),
                nAtomNumberCanonFrom.as_const(),
                at.as_const(),
                num_atoms,
                bAllene,
            )?;
            if found == 0 && bAllene == 0 && first_canon_rank != 0 {
                let mut allene_rank1 = 0_u16;
                let mut allene_rank2 = 0_u16;
                let mut allene_min1 = 0_u16;
                let mut allene_min2 = 0_u16;
                found = Next_SB_At_CanonRanks2(
                    heap,
                    &mut allene_rank1,
                    &mut allene_rank2,
                    &mut allene_min1,
                    &mut allene_min2,
                    &mut first_canon_rank,
                    pCS.bAtomUsedForStereo.as_const(),
                    pRankStack1,
                    pRankStack2,
                    nCanonRankFrom.as_const(),
                    nAtomNumberCanonFrom.as_const(),
                    at.as_const(),
                    num_atoms,
                    1,
                )?;
                if found != 0 {
                    canon_rank1 = allene_rank1;
                    canon_rank2 = allene_rank2;
                    canon_min1 = allene_min1;
                    canon_min2 = allene_min2;
                    bAllene = 1;
                }
            }
            let current = source_get(heap, pCS.LinearCTStereoDble, nNumMappedBonds)?;
            if found == 0
                || (pCS.bStereoIsBetter == 0
                    && (canon_rank1 > current.at_num1
                        || (canon_rank1 == current.at_num1 && canon_rank2 > current.at_num2)))
            {
                if total_success == 0 {
                    source_set(heap, pCS.LinearCTStereoDble, nNumMappedBonds, previous_bond)?;
                }
                return Ok(total_success);
            }

            let rank1_pointer = source_get(heap, pRankStack1, 0)?;
            let order2_pointer = source_get(heap, pRankStack2, 1)?;
            let rank2_pointer = source_get(heap, pRankStack2, 0)?;
            let from1 = source_get(
                heap,
                nAtomNumberCanonFrom,
                i32::from(canon_rank1).wrapping_sub(1),
            )?;
            let from2 = source_get(
                heap,
                nAtomNumberCanonFrom,
                i32::from(canon_rank2).wrapping_sub(1),
            )?;
            let mapping_rank1 = source_get(heap, rank1_pointer, i32::from(from1))?;
            let mapping_rank2 = source_get(heap, rank1_pointer, i32::from(from2))?;
            let max_index = i32::from(mapping_rank1).wrapping_sub(1);
            let rank1_check_atom = source_get(
                heap,
                order2_pointer,
                i32::from(mapping_rank1).wrapping_sub(1),
            )?;
            let rank2_check_atom = source_get(
                heap,
                order2_pointer,
                i32::from(mapping_rank2).wrapping_sub(1),
            )?;
            if mapping_rank1 != source_get(heap, rank2_pointer, i32::from(rank1_check_atom))?
                || mapping_rank2 != source_get(heap, rank2_pointer, i32::from(rank2_check_atom))?
            {
                return Ok(CT_STEREOCOUNT_ERR);
            }

            let mut number_choices = 0_i32;
            let mut number_unknown = 0_i32;
            let mut number_undefined = 0_i32;
            let mut number_best = 0_i32;
            let mut number_worse = 0_i32;
            let mut number_calculate = 0_i32;
            let mut candidate1 = 0_i32;
            while candidate1 <= max_index {
                let to1 = source_get(heap, order2_pointer, max_index.wrapping_sub(candidate1))?;
                if mapping_rank1 != source_get(heap, rank2_pointer, i32::from(to1))? {
                    break;
                }
                let atom1 = source_get(heap, at, i32::from(to1))?;
                if atom1.stereo_bond_neighbor[0] != 0 {
                    for slot in 0..MAX_NUM_STEREO_BONDS as usize {
                        let neighbor = atom1.stereo_bond_neighbor[slot];
                        if neighbor == 0 {
                            break;
                        }
                        let to2 = neighbor.wrapping_sub(1);
                        if mapping_rank2 != source_get(heap, rank2_pointer, i32::from(to2))? {
                            continue;
                        }
                        match parity_value(atom1.stereo_bond_parity[slot]) {
                            value if value == AB_PARITY_UNDF as i32 => number_undefined += 1,
                            value if value == AB_PARITY_UNKN as i32 => number_unknown += 1,
                            value if value == BEST_PARITY as i32 => number_best += 1,
                            value if value == WORSE_PARITY as i32 => number_worse += 1,
                            value if value == AB_PARITY_CALC as i32 => number_calculate += 1,
                            value if value == AB_PARITY_NONE as i32 => {
                                continue;
                            }
                            _ => {}
                        }
                        number_choices += 1;
                    }
                }
                candidate1 = candidate1.wrapping_add(1);
            }
            if number_choices
                != number_calculate + number_undefined + number_unknown + number_best + number_worse
            {
                return Ok(CT_STEREOCOUNT_ERR);
            }
            if number_choices == 0 {
                continue 'canon_ranks;
            }

            let mut calculated_parity = if number_calculate > 0 {
                BEST_PARITY as i32
            } else {
                0
            };
            let mut requested_parity = 0_i32;
            let mut previous_requested_parity = 0_i32;
            let mut pass = 0_i32;

            'parity_pass: loop {
                if nNumMappedBonds == 0 {
                    pCS.bStereoIsBetter = 0;
                }
                if pass == 0 {
                    requested_parity = if calculated_parity != 0 {
                        BEST_PARITY as i32
                    } else if number_best != 0 {
                        BEST_PARITY as i32
                    } else if number_worse != 0 {
                        WORSE_PARITY as i32
                    } else if number_unknown != 0 {
                        AB_PARITY_UNKN as i32
                    } else if number_undefined != 0 {
                        AB_PARITY_UNDF as i32
                    } else {
                        AB_PARITY_NONE as i32
                    };
                } else {
                    previous_requested_parity = requested_parity;
                    let next = NextStereoParity2Test(
                        &mut requested_parity,
                        &mut calculated_parity,
                        number_best,
                        number_worse,
                        number_unknown,
                        number_undefined,
                        number_calculate,
                        vABParityUnknown,
                    );
                    if next == 1 {
                        continue 'canon_ranks;
                    }
                    if next != 0 {
                        return Ok(next);
                    }
                }
                pass = pass.wrapping_add(1);
                if requested_parity == AB_PARITY_NONE as i32 {
                    return Ok(CT_STEREOCOUNT_ERR);
                }
                if pCS.bStereoIsBetter == 0 {
                    let compare = if total_success != 0 {
                        CompareLinCtStereoDoubleToValues(
                            heap,
                            pCS.LinearCTStereoDble
                                .offset(i64::from(nNumMappedBonds))?
                                .as_const(),
                            canon_rank1,
                            canon_rank2,
                            requested_parity as u8,
                        )?
                    } else if previous_bond.at_num1 > canon_rank1 {
                        1
                    } else if previous_bond.at_num1 != canon_rank1 {
                        -1
                    } else if previous_bond.at_num2 > canon_rank2 {
                        1
                    } else if previous_bond.at_num2 != canon_rank2 {
                        -1
                    } else if previous_bond.parity > requested_parity as u8 {
                        1
                    } else if previous_bond.parity != requested_parity as u8 {
                        -1
                    } else {
                        0
                    };
                    if compare < 0 {
                        if total_success == 0 {
                            source_set(
                                heap,
                                pCS.LinearCTStereoDble,
                                nNumMappedBonds,
                                previous_bond,
                            )?;
                        }
                        CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
                        return Ok(total_success);
                    }
                }

                let mut all_identical = false;
                let mut last_mapped_to1 = None;
                let mut better_set_here = false;
                if nNumMappedBonds == 0 && previous_requested_parity != requested_parity {
                    pCS.bStereoIsBetter = 0;
                }
                if previous_requested_parity != requested_parity {
                    CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
                }

                let mut to1_scan = 0_i32;
                while to1_scan <= max_index {
                    let to1 = source_get(heap, order2_pointer, max_index.wrapping_sub(to1_scan))?;
                    if mapping_rank1 != source_get(heap, rank2_pointer, i32::from(to1))? {
                        break;
                    }
                    to1_scan = to1_scan.wrapping_add(1);
                    let atom1 = source_get(heap, at, i32::from(to1))?;
                    if atom1.stereo_bond_neighbor[0] == 0 {
                        continue;
                    }
                    if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                        && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_rank1)? == 1
                        && CurTreeIsLastAtomEqu(
                            heap,
                            cur_tree.as_deref(),
                            i32::from(to1),
                            nSymmStereo.as_const(),
                        )? == 1
                    {
                        continue;
                    }
                    let mut candidate_success = 0_i32;
                    let mut all_identical_after_first = false;
                    let mut stack_ptr = [0_i32; 6];
                    let mut mapped_ranks = [0_i32; 6];

                    for slot in 0..MAX_NUM_STEREO_BONDS as usize {
                        let neighbor = atom1.stereo_bond_neighbor[slot];
                        if neighbor == 0 {
                            break;
                        }
                        let to2 = neighbor.wrapping_sub(1);
                        if mapping_rank2 != source_get(heap, rank2_pointer, i32::from(to2))? {
                            continue;
                        }
                        let encoded_parity = atom1.stereo_bond_parity[slot];
                        if parity_known(encoded_parity) {
                            if requested_parity == calculated_parity
                                || requested_parity != parity_value(encoded_parity)
                            {
                                continue;
                            }
                        } else if parity_calculate(encoded_parity) {
                            if requested_parity != calculated_parity {
                                continue;
                            }
                        } else {
                            return Ok(CT_STEREOCOUNT_ERR);
                        }

                        mapped_ranks[0] = nNumMappedRanksInput;
                        let mut stack_index = 0_usize;
                        let mut add_stack = 0_i32;
                        all_identical = encoded_parity & KNOWN_PARITIES_EQL as i8 != 0
                            && parity_known(encoded_parity);
                        if !all_identical
                            && number_calculate == 0
                            && i32::from(number_undefined == 0)
                                + i32::from(number_unknown == 0)
                                + i32::from(number_best == 0)
                                + i32::from(number_worse == 0)
                                == 3
                        {
                            let same = All_SB_Same(
                                heap,
                                canon_rank1,
                                canon_rank2,
                                pRankStack1,
                                pRankStack2,
                                nAtomNumberCanonFrom.as_const(),
                                at.as_const(),
                            )?;
                            if same < 0 {
                                return Ok(CT_STEREOCOUNT_ERR);
                            }
                            all_identical = same != 0;
                        }

                        if !all_identical {
                            if last_mapped_to1 != Some(to1) {
                                let mut next_mapped = 0_i32;
                                let result = map_an_atom2(
                                    heap,
                                    pCG,
                                    num_at_tg,
                                    num_max,
                                    i32::from(from1),
                                    i32::from(to1),
                                    nTempRank,
                                    mapped_ranks[0],
                                    &mut next_mapped,
                                    pCS,
                                    NeighList,
                                    pRankStack1,
                                    pRankStack2,
                                    &mut add_stack,
                                )?;
                                if returned_error(result) {
                                    return Ok(result);
                                }
                                mapped_ranks[1] = next_mapped;
                                stack_ptr[1] = add_stack;
                                last_mapped_to1 = Some(to1);
                                if add_stack != 0 {
                                    if tpos1 == CurTreeGetPos(cur_tree.as_deref())
                                        || CurTreeIsLastRank(
                                            heap,
                                            cur_tree.as_deref(),
                                            canon_rank1,
                                        )? == 0
                                    {
                                        let _ = CurTreeAddRank(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            canon_rank1,
                                        )?;
                                    }
                                    let _ = CurTreeAddAtom(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                        i32::from(to1),
                                    )?;
                                }
                            }
                            stack_index = 1;
                            let mapped_rank1_pointer = source_get(heap, pRankStack1, stack_ptr[1])?;
                            let mapped_rank2_pointer = source_get(heap, pRankStack2, stack_ptr[1])?;
                            if source_get(heap, mapped_rank1_pointer, i32::from(from2))?
                                != source_get(heap, mapped_rank2_pointer, i32::from(to2))?
                            {
                                if calculated_parity != requested_parity {
                                    number_choices = number_choices.wrapping_sub(1);
                                    number_undefined = number_undefined.wrapping_sub(i32::from(
                                        requested_parity == AB_PARITY_UNDF as i32,
                                    ));
                                    number_unknown = number_unknown.wrapping_sub(i32::from(
                                        requested_parity == AB_PARITY_UNKN as i32,
                                    ));
                                    number_best = number_best.wrapping_sub(i32::from(
                                        requested_parity == BEST_PARITY as i32,
                                    ));
                                    number_worse = number_worse.wrapping_sub(i32::from(
                                        requested_parity == WORSE_PARITY as i32,
                                    ));
                                } else if calculated_parity == BEST_PARITY as i32 {
                                    number_choices = number_choices.wrapping_sub(1);
                                    number_calculate = number_calculate.wrapping_sub(1);
                                }
                                continue;
                            }
                            if stack_ptr[1] > stack_ptr[0] {
                                let same = All_SB_Same(
                                    heap,
                                    canon_rank1,
                                    canon_rank2,
                                    pRankStack1.offset(i64::from(stack_ptr[1]))?,
                                    pRankStack2.offset(i64::from(stack_ptr[1]))?,
                                    nAtomNumberCanonFrom.as_const(),
                                    at.as_const(),
                                )?;
                                if same < 0 {
                                    return Ok(CT_STEREOBOND_ERROR);
                                }
                                all_identical_after_first = same != 0;
                            } else {
                                all_identical_after_first = false;
                            }
                            if all_identical_after_first {
                                ClearPreviousMappings(
                                    heap,
                                    pRankStack1.offset(i64::from(stack_ptr[1].wrapping_add(2)))?,
                                )?;
                            } else {
                                if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                    && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_rank2)?
                                        == 1
                                    && CurTreeIsLastAtomEqu(
                                        heap,
                                        cur_tree.as_deref(),
                                        i32::from(to2),
                                        nSymmStereo.as_const(),
                                    )? == 1
                                {
                                    continue;
                                }
                                let mut next_mapped = 0_i32;
                                let result = map_an_atom2(
                                    heap,
                                    pCG,
                                    num_at_tg,
                                    num_max,
                                    i32::from(from2),
                                    i32::from(to2),
                                    nTempRank,
                                    mapped_ranks[1],
                                    &mut next_mapped,
                                    pCS,
                                    NeighList,
                                    pRankStack1.offset(i64::from(stack_ptr[1]))?,
                                    pRankStack2.offset(i64::from(stack_ptr[1]))?,
                                    &mut add_stack,
                                )?;
                                if returned_error(result) {
                                    return Ok(result);
                                }
                                mapped_ranks[2] = next_mapped;
                                stack_ptr[2] = stack_ptr[1].wrapping_add(add_stack);
                                stack_index = 2;
                                if add_stack != 0 {
                                    if tpos1 == CurTreeGetPos(cur_tree.as_deref())
                                        || CurTreeIsLastRank(
                                            heap,
                                            cur_tree.as_deref(),
                                            canon_rank2,
                                        )? == 0
                                    {
                                        let _ = CurTreeAddRank(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            canon_rank2,
                                        )?;
                                    }
                                    let _ = CurTreeAddAtom(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                        i32::from(to2),
                                    )?;
                                }
                            }
                        } else {
                            ClearPreviousMappings(heap, pRankStack1.offset(2)?)?;
                        }

                        if requested_parity != calculated_parity {
                            let bond_parity = requested_parity;
                            let compare = CompareLinCtStereoDoubleToValues(
                                heap,
                                pCS.LinearCTStereoDble
                                    .offset(i64::from(nNumMappedBonds))?
                                    .as_const(),
                                canon_rank1,
                                canon_rank2,
                                bond_parity as u8,
                            )?;
                            if compare < 0 && pCS.bStereoIsBetter == 0 {
                                pCS.lNumRejectedCT = pCS.lNumRejectedCT.wrapping_add(1);
                                if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                    && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_rank2)?
                                        == 1
                                {
                                    let _ = CurTreeRemoveIfLastAtom(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                        i32::from(to2),
                                    )?;
                                    let _ = CurTreeRemoveLastRankIfNoAtoms(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                    )?;
                                }
                                continue;
                            }

                            let saved = if compare > 0 && pCS.bStereoIsBetter == 0 {
                                pCS.bStereoIsBetter = 1;
                                better_set_here = true;
                                source_get(heap, pCS.LinearCTStereoDble, nNumMappedBonds)?
                            } else {
                                AT_STEREO_DBLE::default()
                            };
                            source_set(
                                heap,
                                pCS.LinearCTStereoDble,
                                nNumMappedBonds,
                                AT_STEREO_DBLE {
                                    at_num1: canon_rank1,
                                    at_num2: canon_rank2,
                                    parity: bond_parity as u8,
                                },
                            )?;
                            let used_rank1 =
                                source_get(heap, pCS.bRankUsedForStereo, i32::from(from1))?;
                            let used_rank2 =
                                source_get(heap, pCS.bRankUsedForStereo, i32::from(from2))?;
                            source_set(
                                heap,
                                pCS.bRankUsedForStereo,
                                i32::from(from1),
                                used_rank1.wrapping_add(1),
                            )?;
                            source_set(
                                heap,
                                pCS.bRankUsedForStereo,
                                i32::from(from2),
                                used_rank2.wrapping_add(1),
                            )?;
                            let used_atom1 =
                                source_get(heap, pCS.bAtomUsedForStereo, i32::from(to1))?;
                            let used_atom2 =
                                source_get(heap, pCS.bAtomUsedForStereo, i32::from(to2))?;
                            if !all_identical {
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to1),
                                    used_atom1.wrapping_sub(1),
                                )?;
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to2),
                                    used_atom2.wrapping_sub(1),
                                )?;
                            }
                            let recursive = map_stereo_bonds4(
                                heap,
                                ic,
                                clock_result,
                                user_action,
                                console_quit,
                                pCG,
                                at,
                                num_atoms,
                                num_at_tg,
                                num_max,
                                bAllene,
                                nCanonRankFrom,
                                nAtomNumberCanonFrom,
                                nCanonRankTo,
                                nSymmRank,
                                pRankStack1.offset(i64::from(stack_ptr[stack_index]))?,
                                pRankStack2.offset(i64::from(stack_ptr[stack_index]))?,
                                nTempRank,
                                mapped_ranks[stack_index],
                                nSymmStereo,
                                NeighList,
                                pCS,
                                cur_tree.as_deref_mut(),
                                nNumMappedBonds.wrapping_add(1),
                                vABParityUnknown,
                            )?;
                            if !all_identical {
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to1),
                                    used_atom1,
                                )?;
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to2),
                                    used_atom2,
                                )?;
                            }
                            source_set(heap, pCS.bRankUsedForStereo, i32::from(from1), used_rank1)?;
                            source_set(heap, pCS.bRankUsedForStereo, i32::from(from2), used_rank2)?;
                            if recursive == 4 {
                                if nNumMappedBonds != 0 {
                                    return Ok(recursive);
                                }
                                pCS.bFirstCT = 1;
                                continue 'total_restart;
                            }
                            if returned_error(recursive) {
                                return Ok(recursive);
                            }
                            if recursive > 0 {
                                total_success |= 1;
                                candidate_success = candidate_success.wrapping_add(1);
                                if better_set_here || recursive & 2 != 0 {
                                    CurTreeKeepLastAtomsOnly(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                        tpos1,
                                        1,
                                    )?;
                                    total_success |= 2;
                                }
                            } else {
                                if better_set_here {
                                    pCS.bStereoIsBetter = 0;
                                    source_set(
                                        heap,
                                        pCS.LinearCTStereoDble,
                                        nNumMappedBonds,
                                        saved,
                                    )?;
                                }
                                if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                    && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_rank2)?
                                        == 1
                                {
                                    let _ = CurTreeRemoveIfLastAtom(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                        i32::from(to2),
                                    )?;
                                    let _ = CurTreeRemoveLastRankIfNoAtoms(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                    )?;
                                }
                            }
                            better_set_here = false;
                            if all_identical || all_identical_after_first {
                                break;
                            }
                            continue;
                        }

                        // Calculated parity maps tied half-bond neighbors below.
                        let rank_from = source_get(heap, pRankStack1, stack_ptr[stack_index])?;
                        let rank_to = source_get(heap, pRankStack2, stack_ptr[stack_index])?;
                        let mut equivalent1 = [EQ_NEIGH::default(), EQ_NEIGH::default()];
                        let mut equivalent2 = [EQ_NEIGH::default(), EQ_NEIGH::default()];
                        let mut parity1 = parity_of_mapped_half_bond(
                            heap,
                            i32::from(from1),
                            i32::from(to1),
                            i32::from(from2),
                            i32::from(to2),
                            at,
                            Some(&mut equivalent1[0]),
                            nCanonRankFrom,
                            rank_from,
                            rank_to,
                        )?;
                        if parity1 == 0 {
                            return Ok(CT_STEREOCOUNT_ERR);
                        }
                        let number1 = if parity1 > 0 { 1 } else { 2 };
                        let mut canon_neighbor1 = 0_u16;

                        for choice1 in 0..number1 {
                            let mut stack_index2 = stack_index;
                            let mut to_neighbor1 = 0_u16;
                            let mut neighbor1_success = 0_i32;
                            if number1 == 2 {
                                canon_neighbor1 = source_get(
                                    heap,
                                    nCanonRankFrom,
                                    i32::from(equivalent1[0].from_at),
                                )?;
                                let from_neighbor1 = equivalent1[0].from_at;
                                to_neighbor1 = equivalent1[0].to_at[choice1 as usize];
                                if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                    && CurTreeIsLastRank(
                                        heap,
                                        cur_tree.as_deref(),
                                        canon_neighbor1,
                                    )? == 1
                                    && CurTreeIsLastAtomEqu(
                                        heap,
                                        cur_tree.as_deref(),
                                        i32::from(to_neighbor1),
                                        nSymmStereo.as_const(),
                                    )? == 1
                                {
                                    continue;
                                }
                                let mut next_mapped = 0_i32;
                                let result = map_an_atom2(
                                    heap,
                                    pCG,
                                    num_at_tg,
                                    num_max,
                                    i32::from(from_neighbor1),
                                    i32::from(to_neighbor1),
                                    nTempRank,
                                    mapped_ranks[stack_index2],
                                    &mut next_mapped,
                                    pCS,
                                    NeighList,
                                    pRankStack1.offset(i64::from(stack_ptr[stack_index2]))?,
                                    pRankStack2.offset(i64::from(stack_ptr[stack_index2]))?,
                                    &mut add_stack,
                                )?;
                                if returned_error(result) {
                                    return Ok(result);
                                }
                                mapped_ranks[stack_index2 + 1] = next_mapped;
                                stack_ptr[stack_index2 + 1] =
                                    stack_ptr[stack_index2].wrapping_add(add_stack);
                                stack_index2 += 1;
                                if stack_index2 >= 6 {
                                    return Ok(CT_OVERFLOW);
                                }
                                if add_stack != 0 {
                                    if tpos1 == CurTreeGetPos(cur_tree.as_deref())
                                        || CurTreeIsLastRank(
                                            heap,
                                            cur_tree.as_deref(),
                                            canon_neighbor1,
                                        )? == 0
                                    {
                                        let _ = CurTreeAddRank(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            canon_neighbor1,
                                        )?;
                                    }
                                    let _ = CurTreeAddAtom(
                                        heap,
                                        cur_tree.as_deref_mut(),
                                        i32::from(to_neighbor1),
                                    )?;
                                }
                                parity1 = parity_of_mapped_half_bond(
                                    heap,
                                    i32::from(from1),
                                    i32::from(to1),
                                    i32::from(from2),
                                    i32::from(to2),
                                    at,
                                    Some(&mut equivalent1[1]),
                                    nCanonRankFrom,
                                    source_get(heap, pRankStack1, stack_ptr[stack_index2])?,
                                    source_get(heap, pRankStack2, stack_ptr[stack_index2])?,
                                )?;
                                if parity1 <= 0 {
                                    return Ok(CT_STEREOCOUNT_ERR);
                                }
                            } else {
                                mapped_ranks[stack_index2 + 1] = mapped_ranks[stack_index2];
                                stack_ptr[stack_index2 + 1] = stack_ptr[stack_index2];
                                stack_index2 += 1;
                            }

                            let mut parity2 = parity_of_mapped_half_bond(
                                heap,
                                i32::from(from2),
                                i32::from(to2),
                                i32::from(from1),
                                i32::from(to1),
                                at,
                                Some(&mut equivalent2[0]),
                                nCanonRankFrom,
                                source_get(heap, pRankStack1, stack_ptr[stack_index2])?,
                                source_get(heap, pRankStack2, stack_ptr[stack_index2])?,
                            )?;
                            if parity2 == 0 {
                                return Ok(CT_STEREOCOUNT_ERR);
                            }
                            let number2 = if parity2 > 0 { 1 } else { 2 };
                            let mut canon_neighbor2 = 0_u16;

                            for choice2 in 0..number2 {
                                let mut stack_index3 = stack_index2;
                                let mut to_neighbor2 = 0_u16;
                                if number2 == 2 {
                                    canon_neighbor2 = source_get(
                                        heap,
                                        nCanonRankFrom,
                                        i32::from(equivalent2[0].from_at),
                                    )?;
                                    let from_neighbor2 = equivalent2[0].from_at;
                                    to_neighbor2 = equivalent2[0].to_at[choice2 as usize];
                                    if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                        && CurTreeIsLastRank(
                                            heap,
                                            cur_tree.as_deref(),
                                            canon_neighbor2,
                                        )? == 1
                                        && CurTreeIsLastAtomEqu(
                                            heap,
                                            cur_tree.as_deref(),
                                            i32::from(to_neighbor2),
                                            nSymmStereo.as_const(),
                                        )? == 1
                                    {
                                        continue;
                                    }
                                    let mut next_mapped = 0_i32;
                                    let result = map_an_atom2(
                                        heap,
                                        pCG,
                                        num_at_tg,
                                        num_max,
                                        i32::from(from_neighbor2),
                                        i32::from(to_neighbor2),
                                        nTempRank,
                                        mapped_ranks[stack_index3],
                                        &mut next_mapped,
                                        pCS,
                                        NeighList,
                                        pRankStack1.offset(i64::from(stack_ptr[stack_index3]))?,
                                        pRankStack2.offset(i64::from(stack_ptr[stack_index3]))?,
                                        &mut add_stack,
                                    )?;
                                    if returned_error(result) {
                                        return Ok(result);
                                    }
                                    mapped_ranks[stack_index3 + 1] = next_mapped;
                                    stack_ptr[stack_index3 + 1] =
                                        stack_ptr[stack_index3].wrapping_add(add_stack);
                                    stack_index3 += 1;
                                    if add_stack != 0 {
                                        if tpos1 == CurTreeGetPos(cur_tree.as_deref())
                                            || CurTreeIsLastRank(
                                                heap,
                                                cur_tree.as_deref(),
                                                canon_neighbor2,
                                            )? == 0
                                        {
                                            let _ = CurTreeAddRank(
                                                heap,
                                                cur_tree.as_deref_mut(),
                                                canon_neighbor2,
                                            )?;
                                        }
                                        let _ = CurTreeAddAtom(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            i32::from(to_neighbor2),
                                        )?;
                                    }
                                    parity2 = parity_of_mapped_half_bond(
                                        heap,
                                        i32::from(from2),
                                        i32::from(to2),
                                        i32::from(from1),
                                        i32::from(to1),
                                        at,
                                        Some(&mut equivalent2[1]),
                                        nCanonRankFrom,
                                        source_get(heap, pRankStack1, stack_ptr[stack_index3])?,
                                        source_get(heap, pRankStack2, stack_ptr[stack_index3])?,
                                    )?;
                                    if parity2 <= 0 {
                                        return Ok(CT_STEREOCOUNT_ERR);
                                    }
                                } else {
                                    mapped_ranks[stack_index3 + 1] = mapped_ranks[stack_index3];
                                    stack_ptr[stack_index3 + 1] = stack_ptr[stack_index3];
                                    stack_index3 += 1;
                                }

                                if parity1 <= 0 || parity2 <= 0 {
                                    return Ok(CT_STEREOCOUNT_ERR);
                                }
                                let mut bond_parity = if atom_parity_well_defined(parity1)
                                    && atom_parity_well_defined(parity2)
                                {
                                    2 - (parity1 + parity2) % 2
                                } else {
                                    parity1.max(parity2)
                                };
                                if atom_parity_well_defined(bond_parity)
                                    && atom1.stereo_bond_z_prod[slot] < 0
                                {
                                    bond_parity = 2 - (bond_parity + 1) % 2;
                                }
                                let compare = CompareLinCtStereoDoubleToValues(
                                    heap,
                                    pCS.LinearCTStereoDble
                                        .offset(i64::from(nNumMappedBonds))?
                                        .as_const(),
                                    canon_rank1,
                                    canon_rank2,
                                    bond_parity as u8,
                                )?;
                                if calculated_parity != bond_parity
                                    || (compare < 0 && pCS.bStereoIsBetter == 0)
                                {
                                    pCS.lNumRejectedCT = pCS.lNumRejectedCT.wrapping_add(1);
                                    if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                        && CurTreeIsLastRank(
                                            heap,
                                            cur_tree.as_deref(),
                                            canon_neighbor2,
                                        )? == 1
                                    {
                                        let _ = CurTreeRemoveIfLastAtom(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            i32::from(to_neighbor2),
                                        )?;
                                        let _ = CurTreeRemoveLastRankIfNoAtoms(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                        )?;
                                    }
                                    continue;
                                }

                                let saved = if compare > 0 && pCS.bStereoIsBetter == 0 {
                                    pCS.bStereoIsBetter = 1;
                                    better_set_here = true;
                                    source_get(heap, pCS.LinearCTStereoDble, nNumMappedBonds)?
                                } else {
                                    AT_STEREO_DBLE::default()
                                };
                                source_set(
                                    heap,
                                    pCS.LinearCTStereoDble,
                                    nNumMappedBonds,
                                    AT_STEREO_DBLE {
                                        at_num1: canon_rank1,
                                        at_num2: canon_rank2,
                                        parity: bond_parity as u8,
                                    },
                                )?;
                                let used_rank1 =
                                    source_get(heap, pCS.bRankUsedForStereo, i32::from(from1))?;
                                let used_rank2 =
                                    source_get(heap, pCS.bRankUsedForStereo, i32::from(from2))?;
                                let used_atom1 =
                                    source_get(heap, pCS.bAtomUsedForStereo, i32::from(to1))?;
                                let used_atom2 =
                                    source_get(heap, pCS.bAtomUsedForStereo, i32::from(to2))?;
                                source_set(
                                    heap,
                                    pCS.bRankUsedForStereo,
                                    i32::from(from1),
                                    used_rank1.wrapping_add(1),
                                )?;
                                source_set(
                                    heap,
                                    pCS.bRankUsedForStereo,
                                    i32::from(from2),
                                    used_rank2.wrapping_add(1),
                                )?;
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to1),
                                    used_atom1.wrapping_sub(1),
                                )?;
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to2),
                                    used_atom2.wrapping_sub(1),
                                )?;
                                let recursive = map_stereo_bonds4(
                                    heap,
                                    ic,
                                    clock_result,
                                    user_action,
                                    console_quit,
                                    pCG,
                                    at,
                                    num_atoms,
                                    num_at_tg,
                                    num_max,
                                    bAllene,
                                    nCanonRankFrom,
                                    nAtomNumberCanonFrom,
                                    nCanonRankTo,
                                    nSymmRank,
                                    pRankStack1.offset(i64::from(stack_ptr[stack_index3]))?,
                                    pRankStack2.offset(i64::from(stack_ptr[stack_index3]))?,
                                    nTempRank,
                                    mapped_ranks[stack_index3],
                                    nSymmStereo,
                                    NeighList,
                                    pCS,
                                    cur_tree.as_deref_mut(),
                                    nNumMappedBonds.wrapping_add(1),
                                    vABParityUnknown,
                                )?;
                                source_set(
                                    heap,
                                    pCS.bRankUsedForStereo,
                                    i32::from(from1),
                                    used_rank1,
                                )?;
                                source_set(
                                    heap,
                                    pCS.bRankUsedForStereo,
                                    i32::from(from2),
                                    used_rank2,
                                )?;
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to1),
                                    used_atom1,
                                )?;
                                source_set(
                                    heap,
                                    pCS.bAtomUsedForStereo,
                                    i32::from(to2),
                                    used_atom2,
                                )?;
                                if recursive == 4 {
                                    if nNumMappedBonds != 0 {
                                        return Ok(recursive);
                                    }
                                    pCS.bFirstCT = 1;
                                    continue 'total_restart;
                                }
                                if returned_error(recursive) {
                                    return Ok(recursive);
                                }
                                if recursive > 0 {
                                    total_success |= 1;
                                    candidate_success = candidate_success.wrapping_add(1);
                                    neighbor1_success = neighbor1_success.wrapping_add(1);
                                    if better_set_here || recursive & 2 != 0 {
                                        CurTreeKeepLastAtomsOnly(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            tpos1,
                                            1,
                                        )?;
                                        total_success |= 2;
                                    }
                                } else {
                                    if better_set_here {
                                        pCS.bStereoIsBetter = 0;
                                        source_set(
                                            heap,
                                            pCS.LinearCTStereoDble,
                                            nNumMappedBonds,
                                            saved,
                                        )?;
                                    }
                                    if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                        && CurTreeIsLastRank(
                                            heap,
                                            cur_tree.as_deref(),
                                            canon_neighbor2,
                                        )? == 1
                                    {
                                        let _ = CurTreeRemoveIfLastAtom(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            i32::from(to_neighbor2),
                                        )?;
                                        let _ = CurTreeRemoveLastRankIfNoAtoms(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                        )?;
                                    }
                                }
                                better_set_here = false;
                            }
                            if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_neighbor2)?
                                    == 1
                            {
                                let _ = CurTreeRemoveLastRank(heap, cur_tree.as_deref_mut())?;
                            }
                            if neighbor1_success == 0
                                && tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_neighbor1)?
                                    == 1
                            {
                                let _ = CurTreeRemoveIfLastAtom(
                                    heap,
                                    cur_tree.as_deref_mut(),
                                    i32::from(to_neighbor1),
                                )?;
                            }
                        }
                        if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                            && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_neighbor1)? == 1
                        {
                            let _ = CurTreeRemoveLastRank(heap, cur_tree.as_deref_mut())?;
                        }
                    }
                    if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                        && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_rank2)? == 1
                    {
                        let _ = CurTreeRemoveLastRank(heap, cur_tree.as_deref_mut())?;
                    }
                    if candidate_success == 0
                        && tpos1 < CurTreeGetPos(cur_tree.as_deref())
                        && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_rank1)? == 1
                    {
                        let _ =
                            CurTreeRemoveIfLastAtom(heap, cur_tree.as_deref_mut(), i32::from(to1))?;
                        let _ = CurTreeRemoveLastRankIfNoAtoms(heap, cur_tree.as_deref_mut())?;
                    }
                    if all_identical {
                        break;
                    }
                }

                if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                    && CurTreeIsLastRank(heap, cur_tree.as_deref(), canon_rank1)? == 1
                {
                    let _ = CurTreeRemoveLastRank(heap, cur_tree.as_deref_mut())?;
                } else if tpos1 != CurTreeGetPos(cur_tree.as_deref()) {
                    return Ok(CT_STEREOCOUNT_ERR);
                }
                if total_success == 0 || requested_parity == calculated_parity {
                    continue 'parity_pass;
                }
                return Ok(total_success);
            }
        }
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn map_stereo_atoms4(
    heap: &mut SourceHeap,
    ic: &mut INCHI_CLOCK,
    clock_result: i64,
    user_action: Option<fn() -> i32>,
    console_quit: Option<fn() -> i32>,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<sp_ATOM>,
    num_atoms: i32,
    num_at_tg: i32,
    num_max: i32,
    nCanonRankFrom: SourceMutPointer<AT_RANK>,
    nAtomNumberCanonFrom: SourceMutPointer<AT_RANK>,
    nCanonRankTo: SourceMutPointer<AT_RANK>,
    nSymmRank: SourceMutPointer<AT_RANK>,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nTempRank: SourceMutPointer<AT_RANK>,
    nNumMappedRanksInput: i32,
    nSymmStereo: SourceMutPointer<AT_RANK>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    pCS: &mut CANON_STAT,
    mut cur_tree: Option<&mut CUR_TREE>,
    nNumMappedAtoms: i32,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c:1126 map_stereo_atoms4
    // INCHI✔️❌: complete active source frame follows verbatim.
    /*
    int map_stereo_atoms4( struct tagINCHI_CLOCK *ic,
                           CANON_GLOBALS *pCG,
                           sp_ATOM *at,
                           int num_atoms,
                           int num_at_tg,
                           int num_max,
                           const AT_RANK *nCanonRankFrom,
                           const AT_RANK *nAtomNumberCanonFrom,
                           AT_RANK *nCanonRankTo, /*  canonical numbering to be mapped */
                           const AT_RANK *nSymmRank,
                           AT_RANK **pRankStack1/*from*/,
                           AT_RANK **pRankStack2/*to*/,
                           AT_RANK *nTempRank,
                           int nNumMappedRanksInput,
                           AT_RANK *nSymmStereo,
                           NEIGH_LIST *NeighList,
                           CANON_STAT *pCS,
                           CUR_TREE *cur_tree,
                           int nNumMappedAtoms,
                           int vABParityUnknown )
    {
        /*
         *   Do not check whether "from" atoms have any stereo features.
         */
        int nTotSuccess = 0;
        AT_STEREO_CARB prevAtom;
        int tpos1;

        memset( &prevAtom, 0, sizeof( prevAtom ) );

        tpos1 = CurTreeGetPos( cur_tree );

        if (nNumMappedAtoms < pCS->nLenLinearCTStereoCarb)
        {
            /* AT_RANK *nRankFrom=*pRankStack1++,  AT_RANK *nAtomNumberFrom=pRankStack1++; */
            /* AT_RANK *nRankTo  =*pRankStack2++,  AT_RANK *nAtomNumberTo  =pRankStack2++; */
            int j1, at_from1, at_to1, /*at_from2, at_to2,*/ iMax, lvl, bStereoIsBetterWasSetHere;
            int istk, bAddStack, nNumAtTo1Success, c, bFirstTime = 1, bAllParitiesIdentical; /* djb-rwth: removing redundant variables */
            EQ_NEIGH EN[5], *pEN;
            int nStackPtr[5], nMappedRanks[5], j[5], *nSP, *nMR, bLastLvlFailed;

            AT_RANK at_rank_canon1, cr[5], at_to[5];
            AT_RANK canon_rank1_min = 0;
            int at_rank1; /*  rank for mapping */
            int nNumChoices, nNumUnkn, nNumUndf, nNumWorse, nNumBest, nNumCalc;
            int stereo_center_parity = 0, sb_parity_calc, pass; /* djb-rwth: removing redundant variable */
            AT_STEREO_CARB prevAtom2;

            prevAtom = pCS->LinearCTStereoCarb[nNumMappedAtoms]; /*  save to restore in case of failure */
            at_rank_canon1 = nNumMappedAtoms ? pCS->LinearCTStereoCarb[nNumMappedAtoms - 1].at_num : 0;

            goto bypass_next_canon_rank_check;

        next_canon_rank:

            if (!pCS->bStereoIsBetter /*??? && !pCS->bFirstCT ???*/ &&
                  at_rank_canon1 >= pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num)
            {
                /*  cannot find next available canonical number */
                if (!nTotSuccess)
                {
                    pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom; /*  restore because of failure */
                }
                CurTreeSetPos( cur_tree, tpos1 );
                return nTotSuccess;
            }

    bypass_next_canon_rank_check:

            CurTreeSetPos( cur_tree, tpos1 );

            /*  find next available canon. number for a stereogenic atom */
            if (!Next_SC_At_CanonRank2( &at_rank_canon1, &canon_rank1_min, &bFirstTime,
                pCS->bAtomUsedForStereo, pRankStack1, pRankStack2,
                nAtomNumberCanonFrom, num_atoms ) ||
                  (!pCS->bStereoIsBetter &&
                  at_rank_canon1 > pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num)) /* djb-rwth: addressing LLVM warning */
            {
                /*  cannot find next available canonical number */
                if (!nTotSuccess)
                {
                    pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom; /*  restore because of failure */
                }
                return nTotSuccess;
            }

            nNumChoices = 0;
            nNumUnkn = 0;
            nNumUndf = 0;
            nNumBest = 0;
            nNumWorse = 0;
            nNumCalc = 0;
            pass = 0;
            /* djb-rwth: removing redundant code */

            /*  get mapping rank for the canon. number */
            at_rank1 = pRankStack1[0][at_from1 = (int) nAtomNumberCanonFrom[at_rank_canon1 - 1]];
            iMax = at_rank1 - 1;
            /*  for debug only */
            if (at_rank1 != pRankStack2[0][pRankStack2[1][at_rank1 - 1]])
            {
                return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
            }

            /*  count special parities of the not mapped yet "to" atoms */
            for (j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1 = pRankStack2[1][iMax - j1]]; j1++)
            {
                if (!at[at_to1].stereo_bond_neighbor[0] && pCS->bAtomUsedForStereo[at_to1] == STEREO_AT_MARK)
                {
                    int no_choice = 0;
                    stereo_center_parity = PARITY_VAL( at[at_to1].stereo_atom_parity );
                    switch (stereo_center_parity)
                    {

                        case AB_PARITY_UNDF: nNumUndf++; break; /*  4 */

                        case AB_PARITY_UNKN: nNumUnkn++;
                            break; /*  3 */

                        case BEST_PARITY:    nNumBest++; break; /*  1 */
                        case WORSE_PARITY:   nNumWorse++; break; /*  2 */
                        case AB_PARITY_CALC: nNumCalc++; break;
                        case AB_PARITY_NONE: no_choice++; break; /*  0 */
                    }
                    nNumChoices += !no_choice;
                }
            }

            if (nNumChoices != nNumCalc + nNumUndf + nNumUnkn + nNumBest + nNumWorse)
            {
                return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
            }
            if (!nNumChoices)
            {
                goto next_canon_rank;
            }
            /*  Determine the first parity to search */
            sb_parity_calc = ( nNumCalc > 0 ) ? BEST_PARITY : 0;

            /*  ==============================================================
                Search sequence:           sb_parity_calc    stereo_center_parity
                ==============================================================
                BEST_PARITY   (calc)       BEST_PARITY     BEST_PARITY
                BEST_PARITY   (known)      BEST_PARITY     WORSE_PARITY  or 0
                WORSE_PARITY  (calc)       WORSE_PARITY    WORSE_PARITY
                WORSE_PARITY  (known)      WORSE_PARITY    0
                AB_PARITY_UNKN(known)      AB_PARITY_UNKN  0
                AB_PARITY_UNDF(known)      AB_PARITY_UNDF  0

                if (sb_parity_calc==stereo_center_parity) then "calc" else "known"
             */

        repeat_all:

            if (!pass++)
            {
                /*  select the smallest parity to search */
                if (sb_parity_calc)
                {
                    stereo_center_parity = BEST_PARITY;
                }
                else
                {
                    stereo_center_parity = nNumBest ? BEST_PARITY :
                        nNumWorse ? WORSE_PARITY :
                        nNumUnkn ? AB_PARITY_UNKN :
                        nNumUndf ? AB_PARITY_UNDF : AB_PARITY_NONE;
                }
            }
            else
            {
                /* djb-rwth: removing redundant code */
                j1 = NextStereoParity2Test( &stereo_center_parity, &sb_parity_calc,
                                         nNumBest, nNumWorse, nNumUnkn, nNumUndf, nNumCalc,
                                         vABParityUnknown );
                switch (j1)
                {
                    case 0:
                        break; /* obtained next parity to test */
                    case 1:
                        goto next_canon_rank;
                    default:
                        return j1; /* program error */
                }
            }

            if (stereo_center_parity == AB_PARITY_NONE)
            {
                /*  error? */
                return CT_STEREOCOUNT_ERR;                  /*   <BRKPT> */
            }

            /*  check if the new requested parity is small enough */
            if (!pCS->bStereoIsBetter)
            {
                c = CompareLinCtStereoAtomToValues( nTotSuccess ? pCS->LinearCTStereoCarb + nNumMappedAtoms : &prevAtom,
                                  at_rank_canon1, (U_CHAR) stereo_center_parity );
                if (c < 0)
                {
                    if (!nTotSuccess)
                    {
                        pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom;
                    }
                    CurTreeSetPos( cur_tree, tpos1 );
                    return nTotSuccess;
                }
            }

            /* djb-rwth: removing redundant code */
            bStereoIsBetterWasSetHere = 0;
            /* djb-rwth: removing redundant code */
            CurTreeSetPos( cur_tree, tpos1 );  /*  start over */

            /*
            if ( prev_stereo_center_parity != stereo_center_parity ) {
                CurTreeSetPos( cur_tree, tpos1 );
            }
            */
            /*  nRankTo                 nAtomNumberTo */
            for (j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1 = pRankStack2[1][iMax - j1]]; j1++)
            {
                int ret, ret1, ret2, parity1;
                nNumAtTo1Success = 0;
                /*
                if ( !(at[at_to1].stereo_atom_parity && !at[at_to1].stereo_bond_neighbor[0] &&
                       pCS->bAtomUsedForStereo[at_to1] == STEREO_AT_MARK ) )
                */
                if (!at[at_to1].stereo_atom_parity || at[at_to1].stereo_bond_neighbor[0] ||
                      pCS->bAtomUsedForStereo[at_to1] != STEREO_AT_MARK) /* simplify 12-17-2003 */
                {
                    continue;
                }
                /* Do not map on non-stereogenic atom constitutionally
                * equivalent to a steregenic atom. Here
                * at[at_to1] is not a sterereo center;  |       |
                * bonds tautomerism is a usual cause.  -P(+)-CH=P-
                * For example, consider a fragment:     |       |
                * The two atoms P may be constitutionally
                * equivalent, P(+) may be seen as a stereocenter
                * while another P has a double bond (Now such a P(V) IS a stereocenter).
                */

                /*  check whether the stereocenter parity corresponds to the requested stereocenter parity */
                if (PARITY_KNOWN( at[at_to1].stereo_atom_parity ))
                {
                    if (stereo_center_parity == sb_parity_calc)
                    {
                        continue;  /*  requested parity to be calculated, found known parity */
                    }
                    if (stereo_center_parity != PARITY_VAL( at[at_to1].stereo_atom_parity ))
                    {
                        continue;  /*  parity differs from the requested parity */
                    }
                }
                else
                {
                    if (PARITY_CALCULATE( at[at_to1].stereo_atom_parity ))
                    {
                        if (stereo_center_parity != sb_parity_calc)
                        {
                            continue;  /*  requested known parity, found patity to be calculated */
                        }
                    }
                    else
                    {
                        return CT_STEREOCOUNT_ERR;  /*  unknown parity type */
                    }
                }

                bAllParitiesIdentical = ( ( at[at_to1].stereo_atom_parity & KNOWN_PARITIES_EQL ) &&
                                         PARITY_KNOWN( at[at_to1].stereo_atom_parity ) );

                if (!bAllParitiesIdentical && !nNumCalc &&
                    ( !nNumUndf + !nNumUnkn + !nNumBest + !nNumWorse ) == 3)
                {
                    /* only one kind of stereocenter parity is present; check whether all parities are really same */
                    bAllParitiesIdentical = All_SC_Same( at_rank_canon1, /*  canonical number */
                                                          pRankStack1, pRankStack2,
                                                          nAtomNumberCanonFrom, at );
                    if (bAllParitiesIdentical < 0)
                    {
                        return CT_STEREOCOUNT_ERR;
                    }
                }
                if (tpos1 < CurTreeGetPos( cur_tree ) &&
                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) &&
                     1 == CurTreeIsLastAtomEqu( cur_tree, at_to1, nSymmStereo ))
                {
                    continue;
                }


                /*  initialize stack pointer nStackPtr[istk] for "hand-made" recursion */
                /*  stacks are pRankStack1[], pRankStack2[], nNumMappedRanks[] */
                istk = 0;
                nStackPtr[istk] = 0;
                nMappedRanks[istk] = nNumMappedRanksInput;
                bAddStack = 0;
                /*  if all equivalent atoms have same known parity, do not map any of them here */
                if (!bAllParitiesIdentical)
                {
                    /*  map the central atom */
                    /*  this mapping is always possible */
                    ret1 = map_an_atom2( pCG, num_at_tg, num_max, at_from1, at_to1,
                                        nTempRank, nMappedRanks[istk], &nMappedRanks[istk + 1], pCS,
                                        NeighList, pRankStack1 + nStackPtr[istk], pRankStack2 + nStackPtr[istk],
                                        &bAddStack );
                    if (RETURNED_ERROR( ret1 ))
                    {
                        return ret1; /*  error */
                    }
                    nStackPtr[istk + 1] = nStackPtr[istk] + bAddStack;
                    istk++;
                }
                else
                {
                    ClearPreviousMappings( pRankStack1 + 2 ); /*  precaution */
                }

                /*********************************************************************************
                 *
                 *   Unknown Stereocenter Parity case: possibly need to map stereo center neighbors
                 */
                if (stereo_center_parity == sb_parity_calc)
                {
                    /*  find out the parity */
                    parity1 = parity_of_mapped_atom2( pCG, at_from1, at_to1, at, &EN[istk],
                                                     nCanonRankFrom, pRankStack1[nStackPtr[istk]],
                                                     pRankStack2[nStackPtr[istk]] );
                    /*  if parity is well-defined then returned EN[istk].num_to=0 */
                    if (!parity1)
                    {
                        return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                    }
                    if (!EN[istk].num_to && parity1 != sb_parity_calc)
                    {
                        continue; /*  looking for the parity value = sb_parity_calc */
                    }
                }
                else
                {
                    /*  Known parity */
                    parity1 = stereo_center_parity;
                    EN[istk].num_to = 0;
                }

                /***********************************************************************
                 * no need to map the neighbors: parity is known or has been calculated
                 */
                if ((stereo_center_parity == sb_parity_calc && !EN[istk].num_to) ||
                     /*  now well-defined, but unknown in advance atom parity  OR   */
                     stereo_center_parity != sb_parity_calc) /* djb-rwth: addressing LLVM warning */
                     /*  known in advance parity = stereo_center_parity */
                {
                    /*  do not need to map the neighbors */
                    c = CompareLinCtStereoAtomToValues( pCS->LinearCTStereoCarb + nNumMappedAtoms,
                                                at_rank_canon1, (U_CHAR) parity1 );
                    if (c < 0 && !pCS->bStereoIsBetter)
                    {
                        /*  reject */
                        pCS->lNumRejectedCT++;
                        continue;  /*  Reject: not a minimal CT. Should not happen */
                    }
                    else
                    {
                        /*  accept */
                        if (bAddStack)
                        {
                            if (tpos1 == CurTreeGetPos( cur_tree ) ||
                                 0 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
                            {
                                CurTreeAddRank( cur_tree, at_rank_canon1 );
                            }
                            CurTreeAddAtom( cur_tree, at_to1 );
                        }

                        if (c > 0 && !pCS->bStereoIsBetter)
                        {
                            /*  stereo center entry is less than the previusly found */
                            pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                            prevAtom2 = pCS->LinearCTStereoCarb[nNumMappedAtoms];
                        }
                        pCS->LinearCTStereoCarb[nNumMappedAtoms].parity = parity1;
                        pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num = at_rank_canon1;
                        pCS->bRankUsedForStereo[at_from1] = 3;
    #if ( FIX_ChCh_STEREO_CANON_BUG == 1 )
                        if (!bAllParitiesIdentical)
    #endif
                            pCS->bAtomUsedForStereo[at_to1] -= STEREO_AT_MARK;


                        ret = map_stereo_atoms4( ic, pCG,
                                                  at, num_atoms, num_at_tg, num_max,
                                                  nCanonRankFrom,
                                                  nAtomNumberCanonFrom,
                                                  nCanonRankTo,
                                                  nSymmRank,
                                                  pRankStack1 + nStackPtr[istk],
                                                  pRankStack2 + nStackPtr[istk],
                                                  nTempRank,
                                                  nMappedRanks[istk],
                                                  nSymmStereo,
                                                  NeighList,
                                                  pCS,
                                                  cur_tree,
                                                  nNumMappedAtoms + 1,
                                                  vABParityUnknown );

                        pCS->bRankUsedForStereo[at_from1] = 0;
    #if ( FIX_ChCh_STEREO_CANON_BUG == 1 )
                        if (!bAllParitiesIdentical)
    #endif
                            pCS->bAtomUsedForStereo[at_to1] += STEREO_AT_MARK;

                        if (ret == 4)
                        {
                            return ret;
                        }

                        if (RETURNED_ERROR( ret ))
                        {
                            return ret; /*  program error */ /* djb-rwth: fixing coverity ID #499567 */
                        }

                        if (ret > 0)
                        {
                            nTotSuccess |= 1;
                            nNumAtTo1Success++;
                            if (bStereoIsBetterWasSetHere || ( ret & 2 ))
                            {
                                CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                                nTotSuccess |= 2; /*  Obtained a smaller CT */
                            }
                        }
                        else
                        {
                            if (bStereoIsBetterWasSetHere)
                            {
                                pCS->bStereoIsBetter = 0;
                                pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom2;
                            }
                            /*  remove failed atom1 from the tree */

                            if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
                            {
                                CurTreeRemoveIfLastAtom( cur_tree, at_to1 );
                                CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                            }
                        }
                        bStereoIsBetterWasSetHere = 0;
                    }

                    /*
                    if ( (at[at_to1].stereo_atom_parity & KNOWN_PARITIES_EQL ) &&
                         ATOM_PARITY_KNOWN(stereo_center_parity) && !nSymmStereo ) { // ??? add && !nSymmStereo ???
                        break; // do not repeat for the same kind of stereo atom with the parity known in advance
                    }
                    */
                    if (bAllParitiesIdentical)
                    {
                        break;  /*  do not repeat for the same kind of stereo atom with the parity known in advance */
                    }
                    continue;
                }

                /***************************************************
                 *
                 * Need to map the neighbors
                 */
                if (stereo_center_parity != sb_parity_calc)
                {
                    return CT_STEREOCOUNT_ERR;  /*  program error */ /*   <BRKPT> */
                }
                /* -- has already been calculated --
                parity1 = parity_of_mapped_atom2( pCG, at_from1, at_to1, at, &EN[istk],
                                                 nCanonRankFrom, pRankStack1[nStackPtr[istk]], pRankStack2[nStackPtr[istk]] );
                */
                if (!parity1)
                {
                    return CT_STEREOCOUNT_ERR; /*  1/25/2002 */ /*   <BRKPT> */
                }

                if (bAddStack)
                {
                    if (tpos1 == CurTreeGetPos( cur_tree ) ||
                         0 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
                    {
                        CurTreeAddRank( cur_tree, at_rank_canon1 );
                    }
                    CurTreeAddAtom( cur_tree, at_to1 );
                }
                /******************************************************
                 * Need to fix the neighbors to define the atom parity
                 ******************************************************/
                /*  a recursion replaced with the hand-made stack */
                lvl = 0;              /*  the "recursion" depth level */
                nSP = &nStackPtr[istk];
                nMR = &nMappedRanks[istk];
                pEN = &EN[istk];
                bLastLvlFailed = 0;

                /*  entering "recursion" depth level lvl */
            next_lvl:
                if (pEN[lvl].num_to)
                {
                    /* Found tied neighbors. Try all transpositions of the tied neighbors.
                     * j is a number of the "to" tied neighbor in the pEN[lvl].to_at[*] to
                     * which the pEN[lvl].from_at "from" neighbor's canonical number is mapped
                     */
                    j[lvl] = 0;
                next_j:
                    cr[lvl] = nCanonRankFrom[pEN[lvl].from_at];
                    at_to[lvl] = pEN[lvl].to_at[j[lvl]];
                    if (j[lvl] && tpos1 < CurTreeGetPos( cur_tree ) &&
                         1 == CurTreeIsLastRank( cur_tree, cr[lvl] ) &&
                         1 == CurTreeIsLastAtomEqu( cur_tree, at_to[lvl], nSymmStereo ))
                    {
                        lvl++;
                        bLastLvlFailed = 0;
                        goto backup; /*  do not test stereo equivalent atoms except the first one */
                    }

                    ret2 = map_an_atom2( pCG, num_at_tg, num_max,
                                        pEN[lvl].from_at,        /* from */
                                        pEN[lvl].to_at[j[lvl]],  /* to */
                                        nTempRank, nMR[lvl], &nMR[lvl + 1], pCS,
                                        NeighList, pRankStack1 + nSP[lvl], pRankStack2 + nSP[lvl],
                                        &bAddStack );

                    if (RETURNED_ERROR( ret2 ))
                    {
                        return ret2; /*  program error */
                    }

                    /*  next recursion depth level */
                    if (bAddStack)
                    {
                        if (tpos1 == CurTreeGetPos( cur_tree ) ||
                             0 == CurTreeIsLastRank( cur_tree, cr[lvl] ))
                        {
                            CurTreeAddRank( cur_tree, cr[lvl] );
                        }
                        CurTreeAddAtom( cur_tree, at_to[lvl] );
                    }
                    nSP[lvl + 1] = nSP[lvl] + bAddStack;
                    lvl++; /*  upon increment lvl = number of additionally mapped neighbors
                             *  (entering next recursion level) */

                    /*  check if the mapping has defined the parity */
                    parity1 = parity_of_mapped_atom2( pCG, at_from1, at_to1, at, &pEN[lvl],
                                                 nCanonRankFrom, pRankStack1[nSP[lvl]], pRankStack2[nSP[lvl]] );

                    if (!parity1)
                    {
                        return CT_STEREOCOUNT_ERR; /*  1/25/2002 */ /*   <BRKPT> */
                    }
                    if (parity1 < 0)
                    {
                        goto next_lvl; /*  we need at least one more mapping to define the parity */
                    }

                    /**********************************************************
                     *
                     *  Check the parity
                     *
                     **********************************************************
                     *  make a decision whether to accept the current mapping */

                    c = CompareLinCtStereoAtomToValues( pCS->LinearCTStereoCarb + nNumMappedAtoms,
                                                at_rank_canon1, (U_CHAR) parity1 );
                    if (sb_parity_calc != parity1 ||
                         (c < 0 && !pCS->bStereoIsBetter)) /* djb-rwth: addressing LLVM warning */
                    {
                        pCS->lNumRejectedCT++;
                        bLastLvlFailed = 1;
                    }
                    else
                       /*  the parity has been defined (all neighbors have untied ranks) */
                       /*  if ( bAcceptAllParities || parity1 == BEST_PARITY ) */
                    {
                        /*********************************************************************
                         *
                         * Process the parity here. We are at the top of the recursion stack.
                         *
                         *********************************************************************/
                        /*  try to accept current neighbors mapping */
                        if (c > 0 && !pCS->bStereoIsBetter)
                        {
                            pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                            prevAtom2 = pCS->LinearCTStereoCarb[nNumMappedAtoms];
                        }
                        pCS->LinearCTStereoCarb[nNumMappedAtoms].parity = parity1;
                        pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num = at_rank_canon1;
                        pCS->bRankUsedForStereo[at_from1] = 3;
                        pCS->bAtomUsedForStereo[at_to1] -= STEREO_AT_MARK;

                        ret = map_stereo_atoms4( ic, pCG, at, num_atoms, num_at_tg, num_max, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                                           nSymmRank, pRankStack1 + nSP[lvl], pRankStack2 + nSP[lvl],
                                           nTempRank, nMR[lvl], nSymmStereo, NeighList,
                                           pCS, cur_tree, nNumMappedAtoms + 1,
                                           vABParityUnknown );

                        pCS->bRankUsedForStereo[at_from1] = 0;
                        pCS->bAtomUsedForStereo[at_to1] += STEREO_AT_MARK;
                        if (ret == 4)
                        {
                            return ret;
                        }
                        if (RETURNED_ERROR( ret ))
                        {
                            if (ret == CT_TIMEOUT_ERR)
                                return ret;
                            else
                                return ret; /*  program error */
                        }
                        if (ret > 0)
                        {
                            nTotSuccess |= 1;
                            nNumAtTo1Success++;
                            if (bStereoIsBetterWasSetHere || ( ret & 2 ))
                            {
                                CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                                nTotSuccess |= 2; /*  Obtained a smaller CT */
                            }
                        }
                        else
                        {
                            if (bStereoIsBetterWasSetHere)
                            {
                                pCS->bStereoIsBetter = 0;
                                pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom2;
                            }
                            bLastLvlFailed = 1;
                        }
                        bStereoIsBetterWasSetHere = 0;

                        /*  avoid redundant repetitions: */
                        /*  check if neighbors mappings have altered another stereo center parity */
                        if (!nSymmStereo && !might_change_other_atom_parity( at, num_atoms, at_to1,
                            pRankStack2[nSP[lvl]] /* ranks after neigbors mapping */,
                            pRankStack2[nStackPtr[istk]] /* ranks before the mapping neighbors */ ))
                        {
                            goto done;
                        }
                    }
                    /*  Continue the cycle. Go to the previous "recursion" level */
                backup:
                    while (lvl-- > 0)
                    {

                        j[lvl] ++; /*  next neighbor at this level */
                        if (j[lvl] < pEN[lvl].num_to)
                        {
                            if (bLastLvlFailed)
                            {
                                if (tpos1 < CurTreeGetPos( cur_tree ) &&
                                     1 == CurTreeIsLastRank( cur_tree, cr[lvl] ))
                                {
                                    CurTreeRemoveIfLastAtom( cur_tree, at_to[lvl] );
                                    CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                                }
                                bLastLvlFailed = 0;
                            }
                            /*  Done with this level. Go back one level */
                            goto next_j;
                        }
                        /*  remove failed atom from the tree */
                        if (tpos1 < CurTreeGetPos( cur_tree ) &&
                             1 == CurTreeIsLastRank( cur_tree, cr[lvl] ))
                        {
                            CurTreeRemoveLastRank( cur_tree );
                        }
                    }
                    goto done;
                }
                else
                {
                    cr[lvl] = 0;
                }

    done:;      /*  at this point lvl=0. */
                if (!nNumAtTo1Success)
                {
                    if (tpos1 < CurTreeGetPos( cur_tree ) &&
                         1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
                    {
                        CurTreeRemoveIfLastAtom( cur_tree, at_to1 );
                        CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                    }
                }
            } /*  end of stereo atom mapping cycle */

            if (tpos1 < CurTreeGetPos( cur_tree ) &&
                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ))
            {
                CurTreeRemoveLastRank( cur_tree );
            }
            else
            {
                /*  CurTree consistency check (debug only) */
                if (tpos1 != CurTreeGetPos( cur_tree ))
                {
                    return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                }
            }

            if (!nTotSuccess || stereo_center_parity == sb_parity_calc)
            {
                goto repeat_all; /*  repeat with next parity if no success or with the same parity, now known */
            }
        }

        else
        {

            /****************************************************
             *
             *  All stereogenic atoms and bonds have been mapped
             *
             ****************************************************/

            if ((UserAction && USER_ACTION_QUIT == ( *UserAction )( )) ||
                 (ConsoleQuit && ( *ConsoleQuit )( ))) /* djb-rwth: addressing LLVM warning */
            {
                return CT_USER_QUIT_ERR;
            }

            if (pCS->bStereoIsBetter || pCS->bFirstCT)
            {
                /* All stereo atoms have been mapped. Current stereo name is better than all previous.
                 * Create new numbering for the new CT
                 * break all remaining "from" ties
                 */
                int i1, ret;
                AT_RANK rc, n1, n2;

                ret = BreakAllTies( pCG, num_at_tg, num_max, pRankStack1,
                                  NeighList, nTempRank, pCS );

                if (RETURNED_ERROR( ret ))
                {
                    return ret;
                }

                /*  break all remaining "from" ties */
                ret = BreakAllTies( pCG, num_at_tg, num_max, pRankStack2,
                                  NeighList, nTempRank, pCS );

                if (RETURNED_ERROR( ret ))
                {
                    return ret;
                }
                /*  move stack pointers to the "nAtomNumber[*]" after all ties are broken */
                pRankStack1 += 2;
                pRankStack2 += 2;
                /* Now final mapping ranks of "to" atom (*pRankStack2)[i] and "from" atom (*pRankStack1)[i]
                 * are equal and all ranks are different, that is, we have a full mapping
                 * Copy so far best canonical numbering from "from" to "to".
                 */
                memset( pCS->nPrevAtomNumber, 0, num_at_tg * sizeof( pCS->nPrevAtomNumber[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                for (i1 = 0; i1 < num_at_tg; i1++)
                {
                    n1 = pRankStack1[1][i1];
                    rc = nCanonRankFrom[n1]; /*  new canon. rank */
                    n2 = pRankStack2[1][i1];                  /*  orig. atom number */
                    nCanonRankTo[n2] = rc;                    /*  assign new canon. number to the atom */
                    /*  use this array to find stereo-equivalent atoms */
                    pCS->nPrevAtomNumber[rc - 1] = n2; /*  ord. number of the atom having canon. rank = rc */
                    nSymmStereo[i1] = i1;            /*  restart search for stereo equivalent atoms */
                    /* check mapping correctness */
                    if (pRankStack1[0][n1] != pRankStack2[0][n2] ||
                         nSymmRank[n1] != nSymmRank[n2])
                    {
                        return CT_STEREO_CANON_ERR; /* stereo mapping error */
                    }
                }
                /*  statistics */
                pCS->lNumTotCT++;
                pCS->lNumEqualCT = 1;
                pCS->lNumDecreasedCT++;
                pCS->bStereoIsBetter = 0; /*  prepare to start over */
                nTotSuccess = 1;
                pCS->bFirstCT = 0;

    #if ( REMOVE_CALC_NONSTEREO == 1 ) /* { */
                if (!( pCS->nMode & CMODE_REDNDNT_STEREO ))
                {
                    i1 = RemoveCalculatedNonStereo( pCG, at, num_atoms, num_at_tg,
                                      pRankStack1, pRankStack2, nTempRank, NeighList,
                                      nSymmRank, nCanonRankTo, pCS->nPrevAtomNumber, pCS,
                                      vABParityUnknown );
                    if (RETURNED_ERROR( i1 ))
                    {
                        return i1;
                    }
                    if (i1 < 0)
                    {
    #if ( bRELEASE_VERSION == 0 )
                        pCS->bExtract |= EXTR_REMOVE_PARITY_WARNING;
    #endif
                        i1 = -( 1 + i1 );
                    }
                    if (i1 > 0)
                    {
                        return 4; /*  total restart: due to newly found stereo equivalence */
                                  /*  the length of the stereo CT has changed */
                    }
                }
    #endif /* } REMOVE_CALC_NONSTEREO */

                /* djb-rwth: removing redundant code */
            }
            else
            {
                /*  current stereo name is same as previous. We do not need a full mapping. */
                if (nSymmStereo)
                {
                    int num_changes = 0;
                    AT_RANK r, n1, n2, r_max, cr;
                    r_max = (AT_RANK) num_at_tg;
                    for (r = 1; r <= r_max; r++)
                    {
                        if (bUniqueAtNbrFromMappingRank( pRankStack1, r, &n1 ))
                        {
                            if (bUniqueAtNbrFromMappingRank( pRankStack2, r, &n2 ))
                            {
                                /*  atoms at[n1], at[n2] have identical untied mapping rank r */
                                cr = nCanonRankFrom[(int) n1] - 1; /*  (new at[n2] canonical rank)-1 */
                                /*  pCS->nPrevAtomNumber[(int)cr] = */
                                /*    previous ordering number of an atom with the canon. rank = cr+1; */
                                /*    make this atom equivalent to atom at[n2]: */
                                num_changes += nJoin2Mcrs( nSymmStereo, pCS->nPrevAtomNumber[(int) cr], n2 );
                            }
                            else
                            {
                                return CT_MAPCOUNT_ERR; /*  mapping ranks must be either both tied or untied. */ /*   <BRKPT> */
                            }
                        }
                    }
                    if (num_changes)
                    { /*  compress trees to stars */
                        for (r = r_max - 1; r; r--)
                        {
                            nGetMcr( nSymmStereo, r );
                        }
                    }
                }
                /*  statistics */
                pCS->lNumEqualCT++;
                pCS->lNumTotCT++;
                nTotSuccess = 1;
            }
            if (bInchiTimeIsOver( ic, pCS->ulTimeOutTime ))
            {
                return CT_TIMEOUT_ERR;
            }
        }
        if (!nTotSuccess && nNumMappedAtoms < pCS->nLenLinearCTStereoCarb)
        {
            pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom;
            CurTreeSetPos( cur_tree, tpos1 );
        }

        return nTotSuccess;  /*  return to the previous level of the recursion. */
    }
        */
    // END INCHI C FUNCTION: map_stereo_atoms4
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION
    // INCHI✔️❌: #define FIX_ChCh_STEREO_CANON_BUG 1
    // INCHI✔️❌: #define REMOVE_CALC_NONSTEREO 1
    // INCHI✔️❌: #define bRELEASE_VERSION 1
    // END INCHI ACTIVE MACRO CONFIGURATION

    let returned_error = |value: i32| (CT_ERR_MIN..=CT_ERR_MAX).contains(&value);
    let parity_value = |value: i8| i32::from(value) & BITS_PARITY as i32;
    let parity_known = |value: i8| {
        (AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&parity_value(value))
    };
    let parity_calculate = |value: i8| parity_value(value) == AB_PARITY_CALC as i32;
    let tpos1 = CurTreeGetPos(cur_tree.as_deref());
    let mut total_success = 0_i32;
    let mut previous_atom = AT_STEREO_CARB::default();

    if nNumMappedAtoms < pCS.nLenLinearCTStereoCarb {
        previous_atom = source_get(heap, pCS.LinearCTStereoCarb, nNumMappedAtoms)?;
        let mut atom_rank_canon = if nNumMappedAtoms != 0 {
            source_get(
                heap,
                pCS.LinearCTStereoCarb,
                nNumMappedAtoms.wrapping_sub(1),
            )?
            .at_num
        } else {
            0
        };
        let mut canon_rank_min = 0_u16;
        let mut first_time = 1_i32;
        let mut bypass_limit_check = true;

        'canon_rank: loop {
            if !bypass_limit_check
                && pCS.bStereoIsBetter == 0
                && atom_rank_canon
                    >= source_get(heap, pCS.LinearCTStereoCarb, nNumMappedAtoms)?.at_num
            {
                if total_success == 0 {
                    source_set(
                        heap,
                        pCS.LinearCTStereoCarb,
                        nNumMappedAtoms,
                        previous_atom.clone(),
                    )?;
                }
                CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
                return Ok(total_success);
            }
            bypass_limit_check = false;
            CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);

            if Next_SC_At_CanonRank2(
                heap,
                &mut atom_rank_canon,
                &mut canon_rank_min,
                &mut first_time,
                pCS.bAtomUsedForStereo.as_const(),
                pRankStack1,
                pRankStack2,
                nAtomNumberCanonFrom.as_const(),
                num_atoms,
            )? == 0
                || (pCS.bStereoIsBetter == 0
                    && atom_rank_canon
                        > source_get(heap, pCS.LinearCTStereoCarb, nNumMappedAtoms)?.at_num)
            {
                if total_success == 0 {
                    source_set(
                        heap,
                        pCS.LinearCTStereoCarb,
                        nNumMappedAtoms,
                        previous_atom.clone(),
                    )?;
                }
                return Ok(total_success);
            }

            // These C arrays are fixed allocations for the lifetime of the
            // canonicalization workspace. Helpers below may update elements in
            // other rank rows and the stereo marker row, but they neither free
            // nor resize these buffers. Stable views therefore retain C pointer
            // identity without retaining element references across helper calls.
            let canon_from = unsafe { heap.stable_slice(nAtomNumberCanonFrom.as_const())? };
            let canon_index = usize::try_from(i32::from(atom_rank_canon).wrapping_sub(1))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:         at_rank1 = pRankStack1[0][at_from1 = (int) nAtomNumberCanonFrom[at_rank_canon1 - 1]];
            let from_atom = *canon_from.get(canon_index)?;
            let rank_stack1 = unsafe { heap.stable_slice(pRankStack1.as_const())? };
            let rank1_pointer = *rank_stack1.get(0)?;
            let rank_stack2 = unsafe { heap.stable_slice(pRankStack2.as_const())? };
            let rank2_pointer = *rank_stack2.get(0)?;
            let order2_pointer = *rank_stack2.get(1)?;
            let rank1 = unsafe { heap.stable_slice(rank1_pointer.as_const())? };
            let rank2 = unsafe { heap.stable_slice(rank2_pointer.as_const())? };
            let order2 = unsafe { heap.stable_slice(order2_pointer.as_const())? };
            let mapping_rank = *rank1.get(usize::from(from_atom))?;
            // INCHI✔️✔️:         iMax = at_rank1 - 1;
            let max_index = i32::from(mapping_rank).wrapping_sub(1);
            let max_index_usize =
                usize::try_from(max_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let check_to = *order2.get(max_index_usize)?;
            // INCHI✔️✔️:         if (at_rank1 != pRankStack2[0][pRankStack2[1][at_rank1 - 1]])
            if mapping_rank != *rank2.get(usize::from(check_to))? {
                // INCHI✔️✔️:             return CT_STEREOCOUNT_ERR;
                return Ok(CT_STEREOCOUNT_ERR);
            }

            let mut number_choices = 0_i32;
            let mut number_unknown = 0_i32;
            let mut number_undefined = 0_i32;
            let mut number_worse = 0_i32;
            let mut number_best = 0_i32;
            let mut number_calculate = 0_i32;
            let atoms = unsafe { heap.stable_slice(at.as_const())? };
            let used_atoms = unsafe { heap.stable_slice(pCS.bAtomUsedForStereo.as_const())? };
            // INCHI✔️✔️:         for (j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1 = pRankStack2[1][iMax - j1]]; j1++)
            let mut scan = 0_i32;
            while scan <= max_index {
                let order_index = usize::try_from(max_index.wrapping_sub(scan))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let to_atom = *order2.get(order_index)?;
                if mapping_rank != *rank2.get(usize::from(to_atom))? {
                    break;
                }
                let atom_value = atoms.get(usize::from(to_atom))?;
                // INCHI✔️✔️:             if (!at[at_to1].stereo_bond_neighbor[0] && pCS->bAtomUsedForStereo[at_to1] == STEREO_AT_MARK)
                if atom_value.stereo_bond_neighbor[0] == 0
                    && *used_atoms.get(usize::from(to_atom))? == STEREO_AT_MARK as i8
                {
                    // INCHI✔️✔️:                 stereo_center_parity = PARITY_VAL( at[at_to1].stereo_atom_parity );
                    match parity_value(atom_value.stereo_atom_parity) {
                        // INCHI✔️✔️:                     case AB_PARITY_UNDF: nNumUndf++; break;
                        value if value == AB_PARITY_UNDF as i32 => number_undefined += 1,
                        // INCHI✔️✔️:                     case AB_PARITY_UNKN: nNumUnkn++; break;
                        value if value == AB_PARITY_UNKN as i32 => number_unknown += 1,
                        // INCHI✔️✔️:                     case BEST_PARITY: nNumBest++; break;
                        value if value == BEST_PARITY as i32 => number_best += 1,
                        // INCHI✔️✔️:                     case WORSE_PARITY: nNumWorse++; break;
                        value if value == WORSE_PARITY as i32 => number_worse += 1,
                        // INCHI✔️✔️:                     case AB_PARITY_CALC: nNumCalc++; break;
                        value if value == AB_PARITY_CALC as i32 => number_calculate += 1,
                        // INCHI✔️✔️:                     case AB_PARITY_NONE: no_choice++; break;
                        value if value == AB_PARITY_NONE as i32 => {
                            scan = scan.wrapping_add(1);
                            continue;
                        }
                        _ => {}
                    }
                    // INCHI✔️✔️:                 nNumChoices += !no_choice;
                    number_choices = number_choices.wrapping_add(1);
                }
                scan = scan.wrapping_add(1);
            }
            // INCHI✔️✔️:         if (nNumChoices != nNumCalc + nNumUndf + nNumUnkn + nNumBest + nNumWorse)
            if number_choices
                != number_calculate + number_undefined + number_unknown + number_best + number_worse
            {
                // INCHI✔️✔️:             return CT_STEREOCOUNT_ERR;
                return Ok(CT_STEREOCOUNT_ERR);
            }
            // INCHI✔️✔️:         if (!nNumChoices)
            if number_choices == 0 {
                // INCHI✔️✔️:             goto next_canon_rank;
                continue 'canon_rank;
            }

            let mut calculated_parity = if number_calculate > 0 {
                BEST_PARITY as i32
            } else {
                0
            };
            let mut stereo_parity = 0_i32;
            let mut pass = 0_i32;

            'parity: loop {
                if pass == 0 {
                    stereo_parity = if calculated_parity != 0 {
                        BEST_PARITY as i32
                    } else if number_best != 0 {
                        BEST_PARITY as i32
                    } else if number_worse != 0 {
                        WORSE_PARITY as i32
                    } else if number_unknown != 0 {
                        AB_PARITY_UNKN as i32
                    } else if number_undefined != 0 {
                        AB_PARITY_UNDF as i32
                    } else {
                        AB_PARITY_NONE as i32
                    };
                } else {
                    let next = NextStereoParity2Test(
                        &mut stereo_parity,
                        &mut calculated_parity,
                        number_best,
                        number_worse,
                        number_unknown,
                        number_undefined,
                        number_calculate,
                        vABParityUnknown,
                    );
                    if next == 1 {
                        continue 'canon_rank;
                    }
                    if next != 0 {
                        return Ok(next);
                    }
                }
                pass = pass.wrapping_add(1);
                if stereo_parity == AB_PARITY_NONE as i32 {
                    return Ok(CT_STEREOCOUNT_ERR);
                }

                if pCS.bStereoIsBetter == 0 {
                    let compare_pointer = if total_success != 0 {
                        pCS.LinearCTStereoCarb.offset(i64::from(nNumMappedAtoms))?
                    } else {
                        let value = if previous_atom.at_num > atom_rank_canon {
                            1
                        } else if previous_atom.at_num != atom_rank_canon {
                            -1
                        } else if previous_atom.parity > stereo_parity as u8 {
                            1
                        } else if previous_atom.parity != stereo_parity as u8 {
                            -1
                        } else {
                            0
                        };
                        if value < 0 {
                            if total_success == 0 {
                                source_set(
                                    heap,
                                    pCS.LinearCTStereoCarb,
                                    nNumMappedAtoms,
                                    previous_atom.clone(),
                                )?;
                            }
                            CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
                            return Ok(total_success);
                        }
                        SourceMutPointer::null()
                    };
                    if !compare_pointer.is_null()
                        && CompareLinCtStereoAtomToValues(
                            heap,
                            compare_pointer.as_const(),
                            atom_rank_canon,
                            stereo_parity as u8,
                        )? < 0
                    {
                        if total_success == 0 {
                            source_set(
                                heap,
                                pCS.LinearCTStereoCarb,
                                nNumMappedAtoms,
                                previous_atom.clone(),
                            )?;
                        }
                        CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
                        return Ok(total_success);
                    }
                }
                CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);

                // INCHI✔️✔️:         for (j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1 = pRankStack2[1][iMax - j1]]; j1++)
                let mut candidate_index = 0_i32;
                while candidate_index <= max_index {
                    let order_index = usize::try_from(max_index.wrapping_sub(candidate_index))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let to_atom = *order2.get(order_index)?;
                    if mapping_rank != *rank2.get(usize::from(to_atom))? {
                        break;
                    }
                    candidate_index = candidate_index.wrapping_add(1);
                    let atom_value = atoms.get(usize::from(to_atom))?;
                    let atom_parity = atom_value.stereo_atom_parity;
                    let has_stereo_bond = atom_value.stereo_bond_neighbor[0] != 0;
                    // End the element borrow before helpers can update other
                    // source allocations; the stable handle itself is a raw
                    // fixed-buffer view and carries no Rust reference alias.
                    let _ = atom_value;
                    // INCHI✔️✔️:             if (!at[at_to1].stereo_atom_parity || at[at_to1].stereo_bond_neighbor[0] || pCS->bAtomUsedForStereo[at_to1] != STEREO_AT_MARK)
                    if atom_parity == 0
                        || has_stereo_bond
                        || *used_atoms.get(usize::from(to_atom))? != STEREO_AT_MARK as i8
                    {
                        continue;
                    }
                    // INCHI✔️✔️:             if (PARITY_KNOWN( at[at_to1].stereo_atom_parity ))
                    if parity_known(atom_parity) {
                        if stereo_parity == calculated_parity
                            || stereo_parity != parity_value(atom_parity)
                        {
                            continue;
                        }
                    } else if parity_calculate(atom_parity) {
                        if stereo_parity != calculated_parity {
                            continue;
                        }
                    } else {
                        return Ok(CT_STEREOCOUNT_ERR);
                    }

                    // INCHI✔️✔️:             bAllParitiesIdentical = ( ( at[at_to1].stereo_atom_parity & KNOWN_PARITIES_EQL ) && PARITY_KNOWN( at[at_to1].stereo_atom_parity ) );
                    let mut all_identical =
                        atom_parity & KNOWN_PARITIES_EQL as i8 != 0 && parity_known(atom_parity);
                    if !all_identical
                        && number_calculate == 0
                        && i32::from(number_undefined == 0)
                            + i32::from(number_unknown == 0)
                            + i32::from(number_best == 0)
                            + i32::from(number_worse == 0)
                            == 3
                    {
                        let same = All_SC_Same(
                            heap,
                            atom_rank_canon,
                            pRankStack1,
                            pRankStack2,
                            nAtomNumberCanonFrom.as_const(),
                            at.as_const(),
                        )?;
                        if same < 0 {
                            return Ok(CT_STEREOCOUNT_ERR);
                        }
                        all_identical = same != 0;
                    }
                    if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                        && CurTreeIsLastRank(heap, cur_tree.as_deref(), atom_rank_canon)? == 1
                        && CurTreeIsLastAtomEqu(
                            heap,
                            cur_tree.as_deref(),
                            i32::from(to_atom),
                            nSymmStereo.as_const(),
                        )? == 1
                    {
                        continue;
                    }

                    let mut stack_ptr = [0_i32; 5];
                    let mut mapped_ranks = [0_i32; 5];
                    mapped_ranks[0] = nNumMappedRanksInput;
                    let mut stack_index = 0_usize;
                    let mut add_stack = 0_i32;
                    if !all_identical {
                        let mut next_mapped = 0_i32;
                        let map_result = map_an_atom2(
                            heap,
                            pCG,
                            num_at_tg,
                            num_max,
                            i32::from(from_atom),
                            i32::from(to_atom),
                            nTempRank,
                            mapped_ranks[0],
                            &mut next_mapped,
                            pCS,
                            NeighList,
                            pRankStack1,
                            pRankStack2,
                            &mut add_stack,
                        )?;
                        if returned_error(map_result) {
                            return Ok(map_result);
                        }
                        mapped_ranks[1] = next_mapped;
                        stack_ptr[1] = add_stack;
                        stack_index = 1;
                    } else {
                        ClearPreviousMappings(heap, pRankStack1.offset(2)?)?;
                    }

                    let rank_from = source_get(heap, pRankStack1, stack_ptr[stack_index])?;
                    let rank_to = source_get(heap, pRankStack2, stack_ptr[stack_index])?;
                    let mut equivalent = std::array::from_fn::<_, 5, _>(|_| EQ_NEIGH::default());
                    let mut parity = if stereo_parity == calculated_parity {
                        parity_of_mapped_atom2(
                            heap,
                            pCG,
                            i32::from(from_atom),
                            i32::from(to_atom),
                            at,
                            Some(&mut equivalent[stack_index]),
                            nCanonRankFrom,
                            rank_from,
                            rank_to,
                        )?
                    } else {
                        equivalent[stack_index].num_to = 0;
                        stereo_parity
                    };
                    if parity == 0 {
                        return Ok(CT_STEREOCOUNT_ERR);
                    }
                    if stereo_parity == calculated_parity
                        && equivalent[stack_index].num_to == 0
                        && parity != calculated_parity
                    {
                        continue;
                    }

                    let direct = (stereo_parity == calculated_parity
                        && equivalent[stack_index].num_to == 0)
                        || stereo_parity != calculated_parity;
                    let mut candidate_success = 0_i32;
                    if direct {
                        let compare = CompareLinCtStereoAtomToValues(
                            heap,
                            pCS.LinearCTStereoCarb
                                .offset(i64::from(nNumMappedAtoms))?
                                .as_const(),
                            atom_rank_canon,
                            parity as u8,
                        )?;
                        if compare < 0 && pCS.bStereoIsBetter == 0 {
                            pCS.lNumRejectedCT = pCS.lNumRejectedCT.wrapping_add(1);
                            continue;
                        }
                        if add_stack != 0 {
                            if tpos1 == CurTreeGetPos(cur_tree.as_deref())
                                || CurTreeIsLastRank(heap, cur_tree.as_deref(), atom_rank_canon)?
                                    == 0
                            {
                                let _ =
                                    CurTreeAddRank(heap, cur_tree.as_deref_mut(), atom_rank_canon)?;
                            }
                            let _ =
                                CurTreeAddAtom(heap, cur_tree.as_deref_mut(), i32::from(to_atom))?;
                        }
                        let mut better_here = false;
                        let saved = if compare > 0 && pCS.bStereoIsBetter == 0 {
                            pCS.bStereoIsBetter = 1;
                            better_here = true;
                            source_get(heap, pCS.LinearCTStereoCarb, nNumMappedAtoms)?
                        } else {
                            AT_STEREO_CARB::default()
                        };
                        source_set(
                            heap,
                            pCS.LinearCTStereoCarb,
                            nNumMappedAtoms,
                            AT_STEREO_CARB {
                                at_num: atom_rank_canon,
                                parity: parity as u8,
                            },
                        )?;
                        source_set(heap, pCS.bRankUsedForStereo, i32::from(from_atom), 3)?;
                        if !all_identical {
                            let value =
                                source_get(heap, pCS.bAtomUsedForStereo, i32::from(to_atom))?;
                            source_set(
                                heap,
                                pCS.bAtomUsedForStereo,
                                i32::from(to_atom),
                                value.wrapping_sub(STEREO_AT_MARK as i8),
                            )?;
                        }
                        let recursive = map_stereo_atoms4(
                            heap,
                            ic,
                            clock_result,
                            user_action,
                            console_quit,
                            pCG,
                            at,
                            num_atoms,
                            num_at_tg,
                            num_max,
                            nCanonRankFrom,
                            nAtomNumberCanonFrom,
                            nCanonRankTo,
                            nSymmRank,
                            pRankStack1.offset(i64::from(stack_ptr[stack_index]))?,
                            pRankStack2.offset(i64::from(stack_ptr[stack_index]))?,
                            nTempRank,
                            mapped_ranks[stack_index],
                            nSymmStereo,
                            NeighList,
                            pCS,
                            cur_tree.as_deref_mut(),
                            nNumMappedAtoms.wrapping_add(1),
                            vABParityUnknown,
                        )?;
                        source_set(heap, pCS.bRankUsedForStereo, i32::from(from_atom), 0)?;
                        if !all_identical {
                            let value =
                                source_get(heap, pCS.bAtomUsedForStereo, i32::from(to_atom))?;
                            source_set(
                                heap,
                                pCS.bAtomUsedForStereo,
                                i32::from(to_atom),
                                value.wrapping_add(STEREO_AT_MARK as i8),
                            )?;
                        }
                        if recursive == 4 {
                            return Ok(4);
                        }
                        if returned_error(recursive) {
                            return Ok(recursive);
                        }
                        if recursive > 0 {
                            candidate_success += 1;
                            total_success |= 1;
                            if better_here || recursive & 2 != 0 {
                                CurTreeKeepLastAtomsOnly(heap, cur_tree.as_deref_mut(), tpos1, 1)?;
                                total_success |= 2;
                            }
                        } else if better_here {
                            pCS.bStereoIsBetter = 0;
                            source_set(heap, pCS.LinearCTStereoCarb, nNumMappedAtoms, saved)?;
                        } else if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                            && CurTreeIsLastRank(heap, cur_tree.as_deref(), atom_rank_canon)? == 1
                        {
                            let _ = CurTreeRemoveIfLastAtom(
                                heap,
                                cur_tree.as_deref_mut(),
                                i32::from(to_atom),
                            )?;
                            let _ = CurTreeRemoveLastRankIfNoAtoms(heap, cur_tree.as_deref_mut())?;
                        }
                    } else {
                        if stereo_parity != calculated_parity || parity == 0 {
                            return Ok(CT_STEREOCOUNT_ERR);
                        }
                        if add_stack != 0 {
                            if tpos1 == CurTreeGetPos(cur_tree.as_deref())
                                || CurTreeIsLastRank(heap, cur_tree.as_deref(), atom_rank_canon)?
                                    == 0
                            {
                                let _ =
                                    CurTreeAddRank(heap, cur_tree.as_deref_mut(), atom_rank_canon)?;
                            }
                            let _ =
                                CurTreeAddAtom(heap, cur_tree.as_deref_mut(), i32::from(to_atom))?;
                        }

                        let mut choices = [0_i32; 5];
                        let mut canonical = [0_u16; 5];
                        let mut chosen_atoms = [0_u16; 5];
                        let base_stack_index = stack_index;
                        let mut level = 0_usize;
                        let mut last_failed = false;
                        'neighbor_search: loop {
                            if equivalent[base_stack_index + level].num_to == 0 {
                                break;
                            }
                            choices[level] = 0;
                            'choice: loop {
                                let en_index = base_stack_index + level;
                                canonical[level] = source_get(
                                    heap,
                                    nCanonRankFrom,
                                    i32::from(equivalent[en_index].from_at),
                                )?;
                                let choice = usize::try_from(choices[level])
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                chosen_atoms[level] = *equivalent[en_index]
                                    .to_at
                                    .get(choice)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if choices[level] != 0
                                    && tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                    && CurTreeIsLastRank(
                                        heap,
                                        cur_tree.as_deref(),
                                        canonical[level],
                                    )? == 1
                                    && CurTreeIsLastAtomEqu(
                                        heap,
                                        cur_tree.as_deref(),
                                        i32::from(chosen_atoms[level]),
                                        nSymmStereo.as_const(),
                                    )? == 1
                                {
                                    last_failed = false;
                                    level += 1;
                                } else {
                                    let current_stack = stack_ptr[base_stack_index + level];
                                    let mut next_mapped = 0_i32;
                                    let map_result = map_an_atom2(
                                        heap,
                                        pCG,
                                        num_at_tg,
                                        num_max,
                                        i32::from(equivalent[en_index].from_at),
                                        i32::from(chosen_atoms[level]),
                                        nTempRank,
                                        mapped_ranks[base_stack_index + level],
                                        &mut next_mapped,
                                        pCS,
                                        NeighList,
                                        pRankStack1.offset(i64::from(current_stack))?,
                                        pRankStack2.offset(i64::from(current_stack))?,
                                        &mut add_stack,
                                    )?;
                                    if returned_error(map_result) {
                                        return Ok(map_result);
                                    }
                                    if add_stack != 0 {
                                        if tpos1 == CurTreeGetPos(cur_tree.as_deref())
                                            || CurTreeIsLastRank(
                                                heap,
                                                cur_tree.as_deref(),
                                                canonical[level],
                                            )? == 0
                                        {
                                            let _ = CurTreeAddRank(
                                                heap,
                                                cur_tree.as_deref_mut(),
                                                canonical[level],
                                            )?;
                                        }
                                        let _ = CurTreeAddAtom(
                                            heap,
                                            cur_tree.as_deref_mut(),
                                            i32::from(chosen_atoms[level]),
                                        )?;
                                    }
                                    stack_ptr[base_stack_index + level + 1] =
                                        current_stack.wrapping_add(add_stack);
                                    mapped_ranks[base_stack_index + level + 1] = next_mapped;
                                    level += 1;
                                    let new_index = base_stack_index + level;
                                    parity = parity_of_mapped_atom2(
                                        heap,
                                        pCG,
                                        i32::from(from_atom),
                                        i32::from(to_atom),
                                        at,
                                        Some(&mut equivalent[new_index]),
                                        nCanonRankFrom,
                                        source_get(heap, pRankStack1, stack_ptr[new_index])?,
                                        source_get(heap, pRankStack2, stack_ptr[new_index])?,
                                    )?;
                                    if parity == 0 {
                                        return Ok(CT_STEREOCOUNT_ERR);
                                    }
                                    if parity < 0 {
                                        continue 'neighbor_search;
                                    }

                                    let compare = CompareLinCtStereoAtomToValues(
                                        heap,
                                        pCS.LinearCTStereoCarb
                                            .offset(i64::from(nNumMappedAtoms))?
                                            .as_const(),
                                        atom_rank_canon,
                                        parity as u8,
                                    )?;
                                    if calculated_parity != parity
                                        || (compare < 0 && pCS.bStereoIsBetter == 0)
                                    {
                                        pCS.lNumRejectedCT = pCS.lNumRejectedCT.wrapping_add(1);
                                        last_failed = true;
                                    } else {
                                        let mut better_here = false;
                                        let saved = if compare > 0 && pCS.bStereoIsBetter == 0 {
                                            pCS.bStereoIsBetter = 1;
                                            better_here = true;
                                            source_get(
                                                heap,
                                                pCS.LinearCTStereoCarb,
                                                nNumMappedAtoms,
                                            )?
                                        } else {
                                            AT_STEREO_CARB::default()
                                        };
                                        source_set(
                                            heap,
                                            pCS.LinearCTStereoCarb,
                                            nNumMappedAtoms,
                                            AT_STEREO_CARB {
                                                at_num: atom_rank_canon,
                                                parity: parity as u8,
                                            },
                                        )?;
                                        source_set(
                                            heap,
                                            pCS.bRankUsedForStereo,
                                            i32::from(from_atom),
                                            3,
                                        )?;
                                        let used = source_get(
                                            heap,
                                            pCS.bAtomUsedForStereo,
                                            i32::from(to_atom),
                                        )?;
                                        source_set(
                                            heap,
                                            pCS.bAtomUsedForStereo,
                                            i32::from(to_atom),
                                            used.wrapping_sub(STEREO_AT_MARK as i8),
                                        )?;
                                        let recursive = map_stereo_atoms4(
                                            heap,
                                            ic,
                                            clock_result,
                                            user_action,
                                            console_quit,
                                            pCG,
                                            at,
                                            num_atoms,
                                            num_at_tg,
                                            num_max,
                                            nCanonRankFrom,
                                            nAtomNumberCanonFrom,
                                            nCanonRankTo,
                                            nSymmRank,
                                            pRankStack1.offset(i64::from(stack_ptr[new_index]))?,
                                            pRankStack2.offset(i64::from(stack_ptr[new_index]))?,
                                            nTempRank,
                                            mapped_ranks[new_index],
                                            nSymmStereo,
                                            NeighList,
                                            pCS,
                                            cur_tree.as_deref_mut(),
                                            nNumMappedAtoms.wrapping_add(1),
                                            vABParityUnknown,
                                        )?;
                                        source_set(
                                            heap,
                                            pCS.bRankUsedForStereo,
                                            i32::from(from_atom),
                                            0,
                                        )?;
                                        source_set(
                                            heap,
                                            pCS.bAtomUsedForStereo,
                                            i32::from(to_atom),
                                            used,
                                        )?;
                                        if recursive == 4 {
                                            return Ok(4);
                                        }
                                        if returned_error(recursive) {
                                            return Ok(recursive);
                                        }
                                        if recursive > 0 {
                                            total_success |= 1;
                                            candidate_success += 1;
                                            if better_here || recursive & 2 != 0 {
                                                CurTreeKeepLastAtomsOnly(
                                                    heap,
                                                    cur_tree.as_deref_mut(),
                                                    tpos1,
                                                    1,
                                                )?;
                                                total_success |= 2;
                                            }
                                        } else {
                                            if better_here {
                                                pCS.bStereoIsBetter = 0;
                                                source_set(
                                                    heap,
                                                    pCS.LinearCTStereoCarb,
                                                    nNumMappedAtoms,
                                                    saved,
                                                )?;
                                            }
                                            last_failed = true;
                                        }
                                        if nSymmStereo.is_null()
                                            && might_change_other_atom_parity(
                                                heap,
                                                at,
                                                num_atoms,
                                                i32::from(to_atom),
                                                source_get(
                                                    heap,
                                                    pRankStack2,
                                                    stack_ptr[new_index],
                                                )?,
                                                source_get(
                                                    heap,
                                                    pRankStack2,
                                                    stack_ptr[base_stack_index],
                                                )?,
                                            )? == 0
                                        {
                                            break 'neighbor_search;
                                        }
                                    }
                                }

                                loop {
                                    if level == 0 {
                                        break 'neighbor_search;
                                    }
                                    level -= 1;
                                    choices[level] = choices[level].wrapping_add(1);
                                    let en_index = base_stack_index + level;
                                    if choices[level] < equivalent[en_index].num_to {
                                        if last_failed {
                                            if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                                && CurTreeIsLastRank(
                                                    heap,
                                                    cur_tree.as_deref(),
                                                    canonical[level],
                                                )? == 1
                                            {
                                                let _ = CurTreeRemoveIfLastAtom(
                                                    heap,
                                                    cur_tree.as_deref_mut(),
                                                    i32::from(chosen_atoms[level]),
                                                )?;
                                                let _ = CurTreeRemoveLastRankIfNoAtoms(
                                                    heap,
                                                    cur_tree.as_deref_mut(),
                                                )?;
                                            }
                                            last_failed = false;
                                        }
                                        continue 'choice;
                                    }
                                    if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                                        && CurTreeIsLastRank(
                                            heap,
                                            cur_tree.as_deref(),
                                            canonical[level],
                                        )? == 1
                                    {
                                        let _ =
                                            CurTreeRemoveLastRank(heap, cur_tree.as_deref_mut())?;
                                    }
                                }
                            }
                        }
                    }

                    if candidate_success == 0
                        && tpos1 < CurTreeGetPos(cur_tree.as_deref())
                        && CurTreeIsLastRank(heap, cur_tree.as_deref(), atom_rank_canon)? == 1
                    {
                        let _ = CurTreeRemoveIfLastAtom(
                            heap,
                            cur_tree.as_deref_mut(),
                            i32::from(to_atom),
                        )?;
                        let _ = CurTreeRemoveLastRankIfNoAtoms(heap, cur_tree.as_deref_mut())?;
                    }
                    if all_identical {
                        break;
                    }
                }

                if tpos1 < CurTreeGetPos(cur_tree.as_deref())
                    && CurTreeIsLastRank(heap, cur_tree.as_deref(), atom_rank_canon)? == 1
                {
                    let _ = CurTreeRemoveLastRank(heap, cur_tree.as_deref_mut())?;
                } else if tpos1 != CurTreeGetPos(cur_tree.as_deref()) {
                    return Ok(CT_STEREOCOUNT_ERR);
                }

                if total_success == 0 || stereo_parity == calculated_parity {
                    continue 'parity;
                }
                return Ok(total_success);
            }
        }
    } else {
        if user_action.is_some_and(|callback| callback() == 1)
            || console_quit.is_some_and(|callback| callback() != 0)
        {
            return Ok(CT_USER_QUIT_ERR);
        }

        if pCS.bStereoIsBetter != 0 || pCS.bFirstCT != 0 {
            let first_break = BreakAllTies(
                heap,
                pCG,
                num_at_tg,
                num_max,
                pRankStack1,
                NeighList,
                nTempRank,
                pCS,
            )?;
            if returned_error(first_break) {
                return Ok(first_break);
            }
            let second_break = BreakAllTies(
                heap,
                pCG,
                num_at_tg,
                num_max,
                pRankStack2,
                NeighList,
                nTempRank,
                pCS,
            )?;
            if returned_error(second_break) {
                return Ok(second_break);
            }
            let final_stack1 = pRankStack1.offset(2)?;
            let final_stack2 = pRankStack2.offset(2)?;
            let rank1 = source_get(heap, final_stack1, 0)?;
            let order1 = source_get(heap, final_stack1, 1)?;
            let rank2 = source_get(heap, final_stack2, 0)?;
            let order2 = source_get(heap, final_stack2, 1)?;
            heap.slice_mut(pCS.nPrevAtomNumber)?
                .get_mut(
                    ..usize::try_from(num_at_tg)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .fill(0);
            let mut index = 0_i32;
            while index < num_at_tg {
                let from = source_get(heap, order1, index)?;
                let canonical = source_get(heap, nCanonRankFrom, i32::from(from))?;
                let to = source_get(heap, order2, index)?;
                source_set(heap, nCanonRankTo, i32::from(to), canonical)?;
                source_set(
                    heap,
                    pCS.nPrevAtomNumber,
                    i32::from(canonical).wrapping_sub(1),
                    to,
                )?;
                source_set(heap, nSymmStereo, index, index as AT_RANK)?;
                if source_get(heap, rank1, i32::from(from))?
                    != source_get(heap, rank2, i32::from(to))?
                    || source_get(heap, nSymmRank, i32::from(from))?
                        != source_get(heap, nSymmRank, i32::from(to))?
                {
                    return Ok(CT_STEREO_CANON_ERR);
                }
                index = index.wrapping_add(1);
            }
            pCS.lNumTotCT = pCS.lNumTotCT.wrapping_add(1);
            pCS.lNumEqualCT = 1;
            pCS.lNumDecreasedCT = pCS.lNumDecreasedCT.wrapping_add(1);
            pCS.bStereoIsBetter = 0;
            total_success = 1;
            pCS.bFirstCT = 0;

            if pCS.nMode & u64::from(CMODE_REDNDNT_STEREO) == 0 {
                let removed = RemoveCalculatedNonStereo(
                    heap,
                    pCG,
                    at,
                    num_atoms,
                    num_at_tg,
                    final_stack1,
                    final_stack2,
                    nTempRank,
                    NeighList,
                    nSymmRank,
                    nCanonRankTo,
                    pCS.nPrevAtomNumber,
                    pCS,
                    vABParityUnknown,
                )?;
                if returned_error(removed) {
                    return Ok(removed);
                }
                let removed = if removed < 0 {
                    removed.wrapping_add(1).wrapping_neg()
                } else {
                    removed
                };
                if removed > 0 {
                    return Ok(4);
                }
            }
        } else {
            if !nSymmStereo.is_null() {
                let mut changes = 0_i32;
                let mut rank = 1_u16;
                let maximum = num_at_tg as AT_RANK;
                while rank <= maximum {
                    let mut from = 0_u16;
                    if bUniqueAtNbrFromMappingRank(heap, pRankStack1, rank, &mut from)? != 0 {
                        let mut to = 0_u16;
                        if bUniqueAtNbrFromMappingRank(heap, pRankStack2, rank, &mut to)? == 0 {
                            return Ok(CT_MAPCOUNT_ERR);
                        }
                        let canonical =
                            source_get(heap, nCanonRankFrom, i32::from(from))?.wrapping_sub(1);
                        changes = changes.wrapping_add(nJoin2Mcrs(
                            heap,
                            nSymmStereo,
                            source_get(heap, pCS.nPrevAtomNumber, i32::from(canonical))?,
                            to,
                        )?);
                    }
                    rank = rank.wrapping_add(1);
                }
                if changes != 0 {
                    let mut rank = maximum.wrapping_sub(1);
                    while rank != 0 {
                        let _ = nGetMcr(heap, nSymmStereo, rank)?;
                        rank = rank.wrapping_sub(1);
                    }
                }
            }
            pCS.lNumEqualCT = pCS.lNumEqualCT.wrapping_add(1);
            pCS.lNumTotCT = pCS.lNumTotCT.wrapping_add(1);
            total_success = 1;
        }

        let timeout: Option<inchiTime> = if pCS.ulTimeOutTime.is_null() {
            None
        } else {
            Some(source_get(heap, pCS.ulTimeOutTime, 0)?)
        };
        if bInchiTimeIsOver(ic, timeout.as_ref(), clock_result) != 0 {
            return Ok(CT_TIMEOUT_ERR);
        }
    }

    if total_success == 0 && nNumMappedAtoms < pCS.nLenLinearCTStereoCarb {
        source_set(heap, pCS.LinearCTStereoCarb, nNumMappedAtoms, previous_atom)?;
        CurTreeSetPos(cur_tree.as_deref_mut(), tpos1);
    }
    Ok(total_success)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::ichisort::CreateNeighList;
    use crate::source_types::{CT_OUT_OF_RAM, tagInchiTime};

    struct Fixture {
        heap: SourceHeap,
        atoms: SourceMutPointer<sp_ATOM>,
        canon_from: SourceMutPointer<AT_RANK>,
        atom_by_canon: SourceMutPointer<AT_RANK>,
        canon_to: SourceMutPointer<AT_RANK>,
        symm: SourceMutPointer<AT_RANK>,
        stack1: ppAT_RANK,
        stack2: ppAT_RANK,
        temporary: SourceMutPointer<AT_RANK>,
        symm_stereo: SourceMutPointer<AT_RANK>,
        neighbors: SourceMutPointer<NEIGH_LIST>,
        stats: CANON_STAT,
    }

    impl Fixture {
        fn new(
            atoms: Vec<sp_ATOM>,
            rank1: Vec<AT_RANK>,
            order1: Vec<AT_RANK>,
            rank2: Vec<AT_RANK>,
            order2: Vec<AT_RANK>,
        ) -> Self {
            let mut heap = SourceHeap::default();
            let count = atoms.len();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let canon_from = heap
                .allocate_model_storage((1..=count as AT_RANK).collect())
                .unwrap();
            let atom_by_canon = heap
                .allocate_model_storage((0..count as AT_RANK).collect())
                .unwrap();
            let canon_to = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let symm = heap
                .allocate_model_storage((1..=count as AT_RANK).collect())
                .unwrap();
            let rank1 = heap.allocate_model_storage(rank1).unwrap();
            let order1 = heap.allocate_model_storage(order1).unwrap();
            let output_rank1 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let output_order1 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let stack1 = heap
                .allocate_model_storage(vec![
                    rank1,
                    order1,
                    output_rank1,
                    output_order1,
                    SourceMutPointer::null(),
                ])
                .unwrap();
            let rank2 = heap.allocate_model_storage(rank2).unwrap();
            let order2 = heap.allocate_model_storage(order2).unwrap();
            let output_rank2 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let output_order2 = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let stack2 = heap
                .allocate_model_storage(vec![
                    rank2,
                    order2,
                    output_rank2,
                    output_order2,
                    SourceMutPointer::null(),
                ])
                .unwrap();
            let temporary = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let symm_stereo = heap
                .allocate_model_storage((0..count as AT_RANK).collect())
                .unwrap();
            let neighbors = CreateNeighList(
                &mut heap,
                count as i32,
                count as i32,
                atoms.as_const(),
                0,
                None,
            )
            .unwrap();
            let b_rank = heap.allocate_model_storage(vec![0_i8; count]).unwrap();
            let b_atom = heap.allocate_model_storage(vec![0_i8; count]).unwrap();
            let previous = heap.allocate_model_storage(vec![0_u16; count]).unwrap();
            let stereo = heap
                .allocate_model_storage(vec![AT_STEREO_CARB::default(); count.max(1)])
                .unwrap();
            let stereo_bonds = heap
                .allocate_model_storage(vec![AT_STEREO_DBLE::default(); count.max(1)])
                .unwrap();
            let stats = CANON_STAT {
                bRankUsedForStereo: b_rank,
                bAtomUsedForStereo: b_atom,
                nPrevAtomNumber: previous,
                LinearCTStereoDble: stereo_bonds,
                LinearCTStereoCarb: stereo,
                nMode: u64::from(CMODE_REDNDNT_STEREO),
                ..CANON_STAT::default()
            };
            Self {
                heap,
                atoms,
                canon_from,
                atom_by_canon,
                canon_to,
                symm,
                stack1,
                stack2,
                temporary,
                symm_stereo,
                neighbors,
                stats,
            }
        }

        fn call(
            &mut self,
            mapped_atoms: i32,
            user_action: Option<fn() -> i32>,
            clock_result: i64,
        ) -> Result<i32, SourceHeapError> {
            let count = self.heap.slice(self.atoms.as_const())?.len() as i32;
            map_stereo_atoms4(
                &mut self.heap,
                &mut INCHI_CLOCK::default(),
                clock_result,
                user_action,
                None,
                &mut CANON_GLOBALS::default(),
                self.atoms,
                count,
                count,
                count,
                self.canon_from,
                self.atom_by_canon,
                self.canon_to,
                self.symm,
                self.stack1,
                self.stack2,
                self.temporary,
                count,
                self.symm_stereo,
                self.neighbors,
                &mut self.stats,
                None,
                mapped_atoms,
                AB_PARITY_UNKN as i32,
            )
        }

        fn call_bonds(&mut self, mapped_bonds: i32) -> Result<i32, SourceHeapError> {
            let count = self.heap.slice(self.atoms.as_const())?.len() as i32;
            let tree_capacity = usize::try_from(count.max(1))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                .saturating_mul(8);
            let tree_pointer = self
                .heap
                .allocate_model_storage(vec![0_u16; tree_capacity])?;
            let mut tree = CUR_TREE {
                tree: tree_pointer,
                max_len: tree_capacity as i32,
                incr_len: count.max(1),
                ..CUR_TREE::default()
            };
            map_stereo_bonds4(
                &mut self.heap,
                &mut INCHI_CLOCK::default(),
                0,
                None,
                None,
                &mut CANON_GLOBALS::default(),
                self.atoms,
                count,
                count,
                count,
                0,
                self.canon_from,
                self.atom_by_canon,
                self.canon_to,
                self.symm,
                self.stack1,
                self.stack2,
                self.temporary,
                count,
                self.symm_stereo,
                self.neighbors,
                &mut self.stats,
                Some(&mut tree),
                mapped_bonds,
                AB_PARITY_UNKN as i32,
            )
        }
    }

    fn stereo_bond_atom(neighbor: u16, parity: i8) -> sp_ATOM {
        let mut atom = sp_ATOM {
            valence: 1,
            parity,
            ..sp_ATOM::default()
        };
        atom.neighbor[0] = neighbor;
        atom.stereo_bond_neighbor[0] = neighbor.wrapping_add(1);
        atom.stereo_bond_parity[0] = parity;
        atom
    }

    #[test]
    fn source_port__ichimap4__map_stereo_bonds4__line_83() {
        let mut handoff =
            Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        handoff.stats.bFirstCT = 1;
        assert_eq!(handoff.call_bonds(0), Ok(1));
        assert_eq!(handoff.stats.bFirstCT, 0);

        let mut no_candidate =
            Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        no_candidate.stats.nLenLinearCTStereoDble = 1;
        no_candidate
            .heap
            .slice_mut(no_candidate.stats.LinearCTStereoDble)
            .unwrap()[0] = AT_STEREO_DBLE {
            at_num1: 7,
            at_num2: 6,
            parity: 2,
        };
        assert_eq!(no_candidate.call_bonds(0), Ok(0));
        assert_eq!(
            no_candidate
                .heap
                .slice(no_candidate.stats.LinearCTStereoDble.as_const())
                .unwrap()[0],
            AT_STEREO_DBLE {
                at_num1: 7,
                at_num2: 6,
                parity: 2,
            }
        );

        let invalid_atoms = vec![stereo_bond_atom(1, 7), stereo_bond_atom(0, 7)];
        let mut invalid = Fixture::new(
            invalid_atoms,
            vec![1, 2],
            vec![0, 1],
            vec![1, 2],
            vec![0, 1],
        );
        invalid.stats.nLenLinearCTStereoDble = 1;
        invalid
            .heap
            .slice_mut(invalid.stats.LinearCTStereoDble)
            .unwrap()[0] = AT_STEREO_DBLE {
            at_num1: 3,
            at_num2: 2,
            parity: 1,
        };
        assert_eq!(invalid.call_bonds(0), Ok(CT_STEREOCOUNT_ERR));

        let known_atoms = vec![
            stereo_bond_atom(1, BEST_PARITY as i8),
            stereo_bond_atom(0, BEST_PARITY as i8),
        ];
        let mut known = Fixture::new(known_atoms, vec![1, 2], vec![0, 1], vec![1, 2], vec![0, 1]);
        known.stats.nLenLinearCTStereoDble = 1;
        known.stats.bFirstCT = 1;
        known
            .heap
            .slice_mut(known.stats.LinearCTStereoDble)
            .unwrap()[0] = AT_STEREO_DBLE {
            at_num1: 3,
            at_num2: 2,
            parity: 2,
        };
        assert_eq!(known.call_bonds(0), Ok(3));
        assert_eq!(
            known
                .heap
                .slice(known.stats.LinearCTStereoDble.as_const())
                .unwrap()[0],
            AT_STEREO_DBLE {
                at_num1: 2,
                at_num2: 1,
                parity: BEST_PARITY as u8,
            }
        );
        assert_eq!(
            known
                .heap
                .slice(known.stats.bRankUsedForStereo.as_const())
                .unwrap(),
            &[0, 0]
        );
        assert_eq!(
            known
                .heap
                .slice(known.stats.bAtomUsedForStereo.as_const())
                .unwrap(),
            &[1, 1]
        );

        let mut calculated_atoms = vec![sp_ATOM::default(); 6];
        calculated_atoms[0].valence = 1;
        calculated_atoms[0].neighbor[0] = 4;
        calculated_atoms[1].valence = 1;
        calculated_atoms[1].neighbor[0] = 4;
        calculated_atoms[2].valence = 1;
        calculated_atoms[2].neighbor[0] = 5;
        calculated_atoms[3].valence = 1;
        calculated_atoms[3].neighbor[0] = 5;
        calculated_atoms[4].valence = 3;
        calculated_atoms[4].neighbor[..3].copy_from_slice(&[5, 0, 1]);
        calculated_atoms[4].parity = BEST_PARITY as i8;
        calculated_atoms[4].stereo_bond_neighbor[0] = 6;
        calculated_atoms[4].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
        calculated_atoms[5].valence = 3;
        calculated_atoms[5].neighbor[..3].copy_from_slice(&[4, 2, 3]);
        calculated_atoms[5].parity = BEST_PARITY as i8;
        calculated_atoms[5].stereo_bond_neighbor[0] = 5;
        calculated_atoms[5].stereo_bond_parity[0] = AB_PARITY_CALC as i8;
        let mut calculated = Fixture::new(
            calculated_atoms,
            vec![1, 2, 3, 4, 5, 6],
            vec![0, 1, 2, 3, 4, 5],
            vec![1, 2, 3, 4, 5, 6],
            vec![0, 1, 2, 3, 4, 5],
        );
        calculated.stats.nLenLinearCTStereoDble = 1;
        calculated.stats.bFirstCT = 1;
        calculated
            .heap
            .slice_mut(calculated.stats.LinearCTStereoDble)
            .unwrap()[0] = AT_STEREO_DBLE {
            at_num1: 7,
            at_num2: 6,
            parity: 2,
        };
        let calculated_result = calculated.call_bonds(0).unwrap();
        assert!(
            calculated_result > 0,
            "calculated result: {calculated_result}"
        );
        let calculated_ct = &calculated
            .heap
            .slice(calculated.stats.LinearCTStereoDble.as_const())
            .unwrap()[0];
        assert_eq!((calculated_ct.at_num1, calculated_ct.at_num2), (6, 5));
        assert!((BEST_PARITY as u8..=WORSE_PARITY as u8).contains(&calculated_ct.parity));
    }

    #[test]
    fn source_port__ichimap4__map_stereo_atoms4__line_1126() {
        let mut no_candidate =
            Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        no_candidate.stats.nLenLinearCTStereoCarb = 1;
        no_candidate
            .heap
            .slice_mut(no_candidate.stats.LinearCTStereoCarb)
            .unwrap()[0] = AT_STEREO_CARB {
            at_num: 7,
            parity: 2,
        };
        assert_eq!(no_candidate.call(0, None, 0), Ok(0));
        assert_eq!(
            no_candidate
                .heap
                .slice(no_candidate.stats.LinearCTStereoCarb.as_const())
                .unwrap()[0],
            AT_STEREO_CARB {
                at_num: 7,
                parity: 2
            }
        );

        let mut invalid = Fixture::new(
            vec![sp_ATOM {
                stereo_atom_parity: 7,
                parity: 1,
                ..sp_ATOM::default()
            }],
            vec![1],
            vec![0],
            vec![1],
            vec![0],
        );
        invalid.stats.nLenLinearCTStereoCarb = 1;
        invalid
            .heap
            .slice_mut(invalid.stats.bAtomUsedForStereo)
            .unwrap()[0] = STEREO_AT_MARK as i8;
        invalid
            .heap
            .slice_mut(invalid.stats.LinearCTStereoCarb)
            .unwrap()[0] = AT_STEREO_CARB {
            at_num: 2,
            parity: 1,
        };
        assert_eq!(invalid.call(0, None, 0), Ok(CT_STEREOCOUNT_ERR));

        let mut recursive = Fixture::new(
            vec![sp_ATOM {
                stereo_atom_parity: BEST_PARITY as i8,
                parity: BEST_PARITY as i8,
                ..sp_ATOM::default()
            }],
            vec![1],
            vec![0],
            vec![1],
            vec![0],
        );
        recursive.stats.nLenLinearCTStereoCarb = 1;
        recursive.stats.bFirstCT = 1;
        recursive
            .heap
            .slice_mut(recursive.stats.bAtomUsedForStereo)
            .unwrap()[0] = STEREO_AT_MARK as i8;
        recursive
            .heap
            .slice_mut(recursive.stats.LinearCTStereoCarb)
            .unwrap()[0] = AT_STEREO_CARB {
            at_num: 2,
            parity: 2,
        };
        assert_eq!(recursive.call(0, None, 0), Ok(3));
        assert_eq!(
            recursive
                .heap
                .slice(recursive.stats.LinearCTStereoCarb.as_const())
                .unwrap()[0],
            AT_STEREO_CARB {
                at_num: 1,
                parity: BEST_PARITY as u8
            }
        );
        assert_eq!(
            recursive
                .heap
                .slice(recursive.stats.bAtomUsedForStereo.as_const())
                .unwrap()[0],
            STEREO_AT_MARK as i8
        );
        assert_eq!(
            recursive
                .heap
                .slice(recursive.stats.bRankUsedForStereo.as_const())
                .unwrap()[0],
            0
        );
        assert_eq!(
            (recursive.stats.lNumTotCT, recursive.stats.lNumDecreasedCT),
            (1, 1)
        );

        let mut first_ct =
            Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        first_ct.stats.bFirstCT = 1;
        first_ct.stats.lNumTotCT = i64::MAX;
        assert_eq!(first_ct.call(0, None, 0), Ok(1));
        assert_eq!(
            first_ct.heap.slice(first_ct.canon_to.as_const()).unwrap(),
            &[1]
        );
        assert_eq!(
            first_ct
                .heap
                .slice(first_ct.stats.nPrevAtomNumber.as_const())
                .unwrap(),
            &[0]
        );
        assert_eq!(first_ct.stats.lNumTotCT, i64::MIN);
        assert_eq!(first_ct.stats.bFirstCT, 0);

        let mut equal = Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        equal.stats.lNumEqualCT = i64::MAX;
        equal.stats.lNumTotCT = 4;
        assert_eq!(equal.call(0, None, 0), Ok(1));
        assert_eq!(
            (equal.stats.lNumEqualCT, equal.stats.lNumTotCT),
            (i64::MIN, 5)
        );

        let mut mismatch = Fixture::new(
            vec![sp_ATOM::default(); 2],
            vec![1, 2],
            vec![0, 1],
            vec![2, 1],
            vec![1, 0],
        );
        mismatch.stats.bFirstCT = 1;
        assert_eq!(mismatch.call(0, None, 0), Ok(CT_STEREO_CANON_ERR));

        fn quit() -> i32 {
            1
        }
        let mut quit_fixture =
            Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        assert_eq!(quit_fixture.call(0, Some(quit), 0), Ok(CT_USER_QUIT_ERR));

        let mut timeout =
            Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        timeout.stats.ulTimeOutTime = timeout
            .heap
            .allocate_model_storage(vec![tagInchiTime { clockTime: 0 }])
            .unwrap();
        assert_eq!(timeout.call(0, None, 1), Ok(CT_TIMEOUT_ERR));

        let mut oom = Fixture::new(vec![sp_ATOM::default()], vec![1], vec![0], vec![1], vec![0]);
        oom.stats.bFirstCT = 1;
        oom.stats.nMode = 0;
        let baseline = oom.heap.live_allocation_count();
        oom.heap.fail_after_allocations(0);
        assert_eq!(oom.call(0, None, 0), Ok(CT_OUT_OF_RAM));
        assert_eq!(oom.heap.live_allocation_count(), baseline);
    }
}
