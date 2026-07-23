use crate::source::base::util::{inchi_calloc, inchi_free, is_el_a_metal};
use crate::source_types::{
    AB_MAX_WELL_DEFINED_PARITY, AB_MIN_WELL_DEFINED_PARITY, AB_PARITY_UNKN, BOND_TYPE_DOUBLE,
    MAX_NUM_STEREO_BONDS, SB_PARITY_FLAG, SB_PARITY_MASK, SB_PARITY_SHFT, STEREO_DBLE_EITHER,
    SourceHeap, SourceHeapError, SourceMutPointer, inp_ATOM,
};

#[allow(non_snake_case)]
pub(crate) fn FixUnkn0DStereoBonds(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_at: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2006 FixUnkn0DStereoBonds
    // INCHI✔️❌: int FixUnkn0DStereoBonds( inp_ATOM *at, int num_at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, m, num = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Add usual Unknown stereobond descriptors to each Unknown bond */
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (AB_PARITY_UNKN == at[i].sb_parity[m])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[i].bond_stereo[(int) at[i].sb_ord[m]] = STEREO_DBLE_EITHER;
    // INCHI✔️❌:                 num++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FixUnkn0DStereoBonds

    let mut num = 0_i32;
    for i in 0..num_at {
        let i = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        for m in 0..MAX_NUM_STEREO_BONDS as usize {
            let (parity, bond_ordinal) = {
                let atom = heap
                    .slice(atoms.as_const())?
                    .get(i)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                (atom.sb_parity[m], atom.sb_ord[m])
            };
            if parity == 0 {
                break;
            }
            if parity == AB_PARITY_UNKN as i8 {
                let bond_ordinal = usize::try_from(bond_ordinal)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice_mut(atoms)?
                    .get_mut(i)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let stereo = atom
                    .bond_stereo
                    .get_mut(bond_ordinal)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                *stereo = STEREO_DBLE_EITHER as i8;
                num = num
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }
    }
    Ok(num)
}

#[allow(non_snake_case)]
pub(crate) fn ReconcileAllCmlBondParities(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    b_disconnected: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4663 ReconcileAllCmlBondParities
    // INCHI✔️❌: int ReconcileAllCmlBondParities( inp_ATOM *at,
    // INCHI✔️❌:                                  int num_atoms,
    // INCHI✔️❌:                                  int bDisconnected )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, ret = 0;
    // INCHI✔️❌:     S_CHAR *visited = (S_CHAR*) inchi_calloc( num_atoms, sizeof( *visited ) );
    // INCHI✔️❌:     if (!visited)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /* out of RAM */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[i].sb_parity[0] && !visited[i] && !( bDisconnected && is_el_a_metal( at[i].el_number ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = ReconcileCmlIncidentBondParities( at, i, -1, visited, bDisconnected ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break; /* error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_free( visited );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ReconcileAllCmlBondParities

    let visited = match inchi_calloc::<i8>(heap, num_atoms as u64, 1) {
        Ok(visited) => visited,
        Err(_) => return Ok(-1),
    };
    let result = (|| {
        let mut ret = 0_i32;
        for i in 0..num_atoms {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = heap
                .slice(atoms.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let was_visited = *heap
                .slice(visited.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom.sb_parity[0] != 0
                && was_visited == 0
                && !(b_disconnected != 0 && is_el_a_metal(i32::from(atom.el_number))? != 0)
            {
                ret =
                    ReconcileCmlIncidentBondParities(heap, atoms, i, -1, visited, b_disconnected)?;
                if ret != 0 {
                    break;
                }
            }
        }
        Ok(ret)
    })();
    inchi_free(heap, visited)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn ReconcileCmlIncidentBondParities(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    cur_atom: i32,
    prev_atom: i32,
    visited: SourceMutPointer<i8>,
    b_disconnected: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4690 ReconcileCmlIncidentBondParities
    // INCHI✔️❌: int ReconcileCmlIncidentBondParities( inp_ATOM *at,
    // INCHI✔️❌:                                       int cur_atom,
    // INCHI✔️❌:                                       int prev_atom,
    // INCHI✔️❌:                                       S_CHAR *visited,
    // INCHI✔️❌:                                       int bDisconnected )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* visited = 0 or parity => atom has not been visited
    // INCHI✔️❌:                  10 + parity => currently is on the stack + its final parity
    // INCHI✔️❌:                  20 + parity => has been visited; is not on the stack anymore + its final parity */
    // INCHI✔️❌:     int i, j, nxt_atom, ret = 0, len;
    // INCHI✔️❌:     int icur2nxt, icur2neigh;   /* cur atom neighbors */
    // INCHI✔️❌:     int inxt2cur, inxt2neigh;   /* next atom neighbors */
    // INCHI✔️❌:     int cur_parity, nxt_parity;
    // INCHI✔️❌:     int cur_order_parity, nxt_order_parity, cur_sb_parity, nxt_sb_parity, bCurMask, bNxtMask;
    // INCHI✔️❌:     /* !(bDisconnected && is_el_a_metal(at[i].el_number) */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[cur_atom].valence > MAX_NUM_STEREO_BONDS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* ignore */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!at[cur_atom].sb_parity[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (visited[cur_atom] >= 10)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 2; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     cur_parity = visited[cur_atom] % 10;
    // INCHI✔️❌:
    // INCHI✔️❌:     visited[cur_atom] += 10;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < MAX_NUM_STEREO_BONDS && at[cur_atom].sb_parity[i]; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         icur2nxt = (int) at[cur_atom].sb_ord[i];
    // INCHI✔️❌:         len = get_opposite_sb_atom( at, cur_atom, icur2nxt, &nxt_atom, &inxt2cur, &j );
    // INCHI✔️❌:         if (!len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 4; /* could not find the opposite atom: bond parity data error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nxt_atom == prev_atom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (visited[nxt_atom] >= 20)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* back edge, second visit: ignore */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[nxt_atom].valence > MAX_NUM_STEREO_BONDS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* may be treated only after metal disconnection */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (bDisconnected && ( at[cur_atom].sb_parity[i] & SB_PARITY_FLAG ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cur_sb_parity = ( at[cur_atom].sb_parity[i] >> SB_PARITY_SHFT );
    // INCHI✔️❌:             bCurMask = 3 << SB_PARITY_SHFT;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cur_sb_parity = ( at[cur_atom].sb_parity[i] & SB_PARITY_MASK );
    // INCHI✔️❌:             bCurMask = 3;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (bDisconnected && ( at[nxt_atom].sb_parity[j] & SB_PARITY_FLAG ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nxt_sb_parity = ( at[nxt_atom].sb_parity[j] >> SB_PARITY_SHFT );
    // INCHI✔️❌:             bNxtMask = 3 << SB_PARITY_SHFT;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nxt_sb_parity = ( at[nxt_atom].sb_parity[j] & SB_PARITY_MASK );
    // INCHI✔️❌:             bNxtMask = 3;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!ATOM_PARITY_WELL_DEF( cur_sb_parity ) ||
    // INCHI✔️❌:              !ATOM_PARITY_WELL_DEF( nxt_sb_parity ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_sb_parity == nxt_sb_parity)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:                 /* goto move_forward;       */
    // INCHI✔️❌:                 /* bypass unknown/undefined */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return 3; /* sb parities do not match: bond parity data error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         icur2neigh = (int) at[cur_atom].sn_ord[i];
    // INCHI✔️❌:         inxt2neigh = (int) at[nxt_atom].sn_ord[j];
    // INCHI✔️❌:         /* Parity of at[cur_atom].neighbor[] premutation to reach this order: { next_atom, neigh_atom, ...} */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* 1. move next_atom  from position=icur2nxt to position=0 =>
    // INCHI✔️❌:          *         icur2nxt permutations
    // INCHI✔️❌:          * 2. move neigh_atom from position=inxt2neigh+(inxt2cur > inxt2neigh) to position=1 =>
    // INCHI✔️❌:          *         inxt2neigh+(inxt2cur > inxt2neigh)-1 permutations.
    // INCHI✔️❌:          * Note if (inxt2cur > inxt2neigh) then move #1 increments neigh_atom position
    // INCHI✔️❌:          * Note add 4 because icur2neigh may be negative due to isotopic H removal
    // INCHI✔️❌:          */
    // INCHI✔️❌:         cur_order_parity = ( 4 + icur2nxt + icur2neigh + ( icur2neigh > icur2nxt ) ) % 2;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Same for next atom: */
    // INCHI✔️❌:         /* parity of at[nxt_atom].neighbor[] premutation to reach this order: { cur_atom, neigh_atom, ...} */
    // INCHI✔️❌:         nxt_order_parity = ( 4 + inxt2cur + inxt2neigh + ( inxt2neigh > inxt2cur ) ) % 2;
    // INCHI✔️❌:         nxt_parity = visited[nxt_atom] % 10;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!cur_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cur_parity = 2 - ( cur_order_parity + cur_sb_parity ) % 2;
    // INCHI✔️❌:             visited[cur_atom] += cur_parity;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_parity != 2 - ( cur_order_parity + cur_sb_parity ) % 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /***** Reconcile bond parities *****/
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Each bond parity is split into two values located at the end atoms.
    // INCHI✔️❌:                    For T (trans) the values are (1,1) or (2,2)
    // INCHI✔️❌:                    For C (cis)   the values are (1,2) or (2,1)
    // INCHI✔️❌:                    The fact that one pair = another with inverted parities, namely
    // INCHI✔️❌:                    Inv(1,1) = (2,2) and Inv(1,2) = (2,1), allows to
    // INCHI✔️❌:                    simultaneouly invert parities of the current bond end atoms
    // INCHI✔️❌:                    (at[cur_atom].sb_parity[i], at[nxt_atom].sb_parity[j])
    // INCHI✔️❌:                    so that the final current atom parity cur_parity
    // INCHI✔️❌:                    calculated later in stereochemical canonicalization for
    // INCHI✔️❌:                    each stereobond incident with the current atomis same.
    // INCHI✔️❌:                    Achieving this is called here RECONCILIATION.
    // INCHI✔️❌:                    If at the closure of an aromatic circuit the parities of
    // INCHI✔️❌:                    next atom cannot be reconciled with already calculated then
    // INCHI✔️❌:                    this function returns 5 (error).
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 at[cur_atom].sb_parity[i] ^= bCurMask;
    // INCHI✔️❌:                 at[nxt_atom].sb_parity[j] ^= bNxtMask;
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                 nxt_sb_parity ^= 3;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!nxt_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nxt_parity = 2 - ( nxt_order_parity + nxt_sb_parity ) % 2;
    // INCHI✔️❌:             visited[nxt_atom] += nxt_parity;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nxt_parity != 2 - ( nxt_order_parity + nxt_sb_parity ) % 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 5; /* algorithm does not work for Mebius-like structures */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Move_forward: */
    // INCHI✔️❌:         if (visited[nxt_atom] < 10)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = ReconcileCmlIncidentBondParities( at, nxt_atom, cur_atom, visited, bDisconnected );
    // INCHI✔️❌:             if (ret)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     visited[cur_atom] += 10; /* all bonds incident to the current atom have
    // INCHI✔️❌:                                 been processed or an error occurred. */
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ReconcileCmlIncidentBondParities

    let cur_index = usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(atoms.as_const())?
        .get(cur_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if i32::from(current.valence) > MAX_NUM_STEREO_BONDS as i32 {
        return Ok(0);
    }
    if current.sb_parity[0] == 0 {
        return Ok(1);
    }
    let visited_value = *heap
        .slice(visited.as_const())?
        .get(cur_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if visited_value >= 10 {
        return Ok(2);
    }
    let mut cur_parity = i32::from(visited_value % 10);
    heap.slice_mut(visited)?[cur_index] = visited_value.wrapping_add(10);

    let mut ret = 0_i32;
    for i in 0..MAX_NUM_STEREO_BONDS as usize {
        let current_sb_parity = heap.slice(atoms.as_const())?[cur_index].sb_parity[i];
        if current_sb_parity == 0 {
            break;
        }
        let current_to_next = i32::from(heap.slice(atoms.as_const())?[cur_index].sb_ord[i]);
        let (mut next_atom, mut next_to_current, mut j) = (0_i32, 0_i32, 0_i32);
        let len = get_opposite_sb_atom(
            heap,
            atoms,
            cur_atom,
            current_to_next,
            Some(&mut next_atom),
            Some(&mut next_to_current),
            Some(&mut j),
        )?;
        if len == 0 {
            return Ok(4);
        }
        if next_atom == prev_atom {
            continue;
        }
        let next_index =
            usize::try_from(next_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if heap.slice(visited.as_const())?[next_index] >= 20 {
            continue;
        }
        if i32::from(heap.slice(atoms.as_const())?[next_index].valence)
            > MAX_NUM_STEREO_BONDS as i32
        {
            continue;
        }
        let j = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let current_parity_bits = i32::from(heap.slice(atoms.as_const())?[cur_index].sb_parity[i]);
        let next_parity_bits = i32::from(heap.slice(atoms.as_const())?[next_index].sb_parity[j]);
        let (cur_sb_parity, current_mask) =
            if b_disconnected != 0 && current_parity_bits & SB_PARITY_FLAG as i32 != 0 {
                (
                    current_parity_bits >> SB_PARITY_SHFT,
                    3_i32 << SB_PARITY_SHFT,
                )
            } else {
                (current_parity_bits & SB_PARITY_MASK as i32, 3)
            };
        let (mut next_sb_parity, next_mask) =
            if b_disconnected != 0 && next_parity_bits & SB_PARITY_FLAG as i32 != 0 {
                (next_parity_bits >> SB_PARITY_SHFT, 3_i32 << SB_PARITY_SHFT)
            } else {
                (next_parity_bits & SB_PARITY_MASK as i32, 3)
            };
        let well_defined = |parity: i32| {
            AB_MIN_WELL_DEFINED_PARITY as i32 <= parity
                && parity <= AB_MAX_WELL_DEFINED_PARITY as i32
        };
        if !well_defined(cur_sb_parity) || !well_defined(next_sb_parity) {
            if cur_sb_parity == next_sb_parity {
                continue;
            }
            return Ok(3);
        }

        let current_to_neighbor = i32::from(heap.slice(atoms.as_const())?[cur_index].sn_ord[i]);
        let next_to_neighbor = i32::from(heap.slice(atoms.as_const())?[next_index].sn_ord[j]);
        let current_order_parity = (4
            + current_to_next
            + current_to_neighbor
            + i32::from(current_to_neighbor > current_to_next))
            % 2;
        let next_order_parity = (4
            + next_to_current
            + next_to_neighbor
            + i32::from(next_to_neighbor > next_to_current))
            % 2;
        let next_visited = heap.slice(visited.as_const())?[next_index];
        let next_parity = i32::from(next_visited % 10);

        if cur_parity == 0 {
            cur_parity = 2 - (current_order_parity + cur_sb_parity) % 2;
            let current_visited = heap.slice(visited.as_const())?[cur_index];
            heap.slice_mut(visited)?[cur_index] = current_visited.wrapping_add(cur_parity as i8);
        } else if cur_parity != 2 - (current_order_parity + cur_sb_parity) % 2 {
            heap.slice_mut(atoms)?[cur_index].sb_parity[i] ^= current_mask as i8;
            heap.slice_mut(atoms)?[next_index].sb_parity[j] ^= next_mask as i8;
            next_sb_parity ^= 3;
        }

        if next_parity == 0 {
            let calculated = 2 - (next_order_parity + next_sb_parity) % 2;
            heap.slice_mut(visited)?[next_index] = next_visited.wrapping_add(calculated as i8);
        } else if next_parity != 2 - (next_order_parity + next_sb_parity) % 2 {
            return Ok(5);
        }

        if heap.slice(visited.as_const())?[next_index] < 10 {
            ret = ReconcileCmlIncidentBondParities(
                heap,
                atoms,
                next_atom,
                cur_atom,
                visited,
                b_disconnected,
            )?;
            if ret != 0 {
                break;
            }
        }
    }
    let value = heap.slice(visited.as_const())?[cur_index];
    heap.slice_mut(visited)?[cur_index] = value.wrapping_add(10);
    Ok(ret)
}

pub(crate) fn get_opposite_sb_atom(
    heap: &SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    cur_atom: i32,
    icur2nxt: i32,
    next_atom_out: Option<&mut i32>,
    next_to_current_out: Option<&mut i32>,
    next_parity_ordinal_out: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4861 get_opposite_sb_atom
    // INCHI✔️❌: int get_opposite_sb_atom( inp_ATOM *at,
    // INCHI✔️❌:                           int cur_atom,
    // INCHI✔️❌:                           int icur2nxt,
    // INCHI✔️❌:                           int *pnxt_atom,
    // INCHI✔️❌:                           int *pinxt2cur,
    // INCHI✔️❌:                           int *pinxt_sb_parity_ord )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB nxt_atom;
    // INCHI✔️❌:     int     j, len;
    // INCHI✔️❌:
    // INCHI✔️❌:     len = 0;
    // INCHI✔️❌:     while (len++ < 20)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Arbitrarily set cumulene length limit to avoid infinite loop */
    // INCHI✔️❌:         nxt_atom = at[cur_atom].neighbor[icur2nxt];
    // INCHI✔️❌:         for (j = 0; j < MAX_NUM_STEREO_BONDS && at[nxt_atom].sb_parity[j]; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_atom == at[nxt_atom].neighbor[(int) at[nxt_atom].sb_ord[j]])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Found the opposite atom */
    // INCHI✔️❌:                 *pnxt_atom = nxt_atom;
    // INCHI✔️❌:                 *pinxt2cur = at[nxt_atom].sb_ord[j];
    // INCHI✔️❌:                 *pinxt_sb_parity_ord = j;
    // INCHI✔️❌:                 return len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (j)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /* reached atom(s) with stereobond (sb) parity, the opposite atom has not been found */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[nxt_atom].valence == 2 && 2 * BOND_TYPE_DOUBLE == at[nxt_atom].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Follow cumulene =X= path */
    // INCHI✔️❌:             icur2nxt = ( at[nxt_atom].neighbor[0] == cur_atom );
    // INCHI✔️❌:             cur_atom = nxt_atom;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /* neither atom with a sb parity not middle cumulene could be reached */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0; /* too long chain of cumulene was found */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: get_opposite_sb_atom

    get_opposite_sb_atom_slice(
        heap.slice(atoms.as_const())?,
        cur_atom,
        icur2nxt,
        next_atom_out,
        next_to_current_out,
        next_parity_ordinal_out,
    )
}

pub(crate) fn get_opposite_sb_atom_slice(
    atoms: &[inp_ATOM],
    mut cur_atom: i32,
    mut icur2nxt: i32,
    mut next_atom_out: Option<&mut i32>,
    mut next_to_current_out: Option<&mut i32>,
    mut next_parity_ordinal_out: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    for len in 1..=20_i32 {
        let current_index =
            usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next_ordinal =
            usize::try_from(icur2nxt).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next_atom = i32::from(
            *atoms
                .get(current_index)
                .and_then(|atom| atom.neighbor.get(next_ordinal))
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let next_index =
            usize::try_from(next_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next = atoms
            .get(next_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut j = 0_usize;
        while j < MAX_NUM_STEREO_BONDS as usize && next.sb_parity[j] != 0 {
            let reverse_ordinal =
                usize::try_from(next.sb_ord[j]).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let reverse_atom = *next
                .neighbor
                .get(reverse_ordinal)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if cur_atom == i32::from(reverse_atom) {
                *next_atom_out
                    .as_deref_mut()
                    .ok_or(SourceHeapError::NullPointer)? = next_atom;
                *next_to_current_out
                    .as_deref_mut()
                    .ok_or(SourceHeapError::NullPointer)? = i32::from(next.sb_ord[j]);
                *next_parity_ordinal_out
                    .as_deref_mut()
                    .ok_or(SourceHeapError::NullPointer)? = j as i32;
                return Ok(len);
            }
            j += 1;
        }
        if j != 0 {
            return Ok(0);
        }
        if next.valence == 2 && (2 * BOND_TYPE_DOUBLE as i32) == i32::from(next.chem_bonds_valence)
        {
            icur2nxt = i32::from(next.neighbor[0] == cur_atom as u16);
            cur_atom = next_atom;
        } else {
            return Ok(0);
        }
    }
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichister__fixunkn0dstereobonds__line_2006() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            FixUnkn0DStereoBonds(&mut heap, SourceMutPointer::null(), 0).unwrap(),
            0
        );
        assert_eq!(
            FixUnkn0DStereoBonds(&mut heap, SourceMutPointer::null(), -1).unwrap(),
            0
        );

        let mut first = inp_ATOM::default();
        first.sb_parity = [AB_PARITY_UNKN as i8, 0, AB_PARITY_UNKN as i8];
        first.sb_ord = [2, 3, 4];
        let mut second = inp_ATOM::default();
        second.sb_parity = [1, AB_PARITY_UNKN as i8, AB_PARITY_UNKN as i8];
        second.sb_ord = [0, 19, 7];
        let atoms = heap.allocate(vec![first, second]).unwrap();
        assert_eq!(FixUnkn0DStereoBonds(&mut heap, atoms, 2).unwrap(), 3);

        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(values[0].bond_stereo[2], STEREO_DBLE_EITHER as i8);
        assert_eq!(values[0].bond_stereo[3], 0);
        assert_eq!(values[0].bond_stereo[4], 0);
        assert_eq!(values[1].bond_stereo[0], 0);
        assert_eq!(values[1].bond_stereo[19], STEREO_DBLE_EITHER as i8);
        assert_eq!(values[1].bond_stereo[7], STEREO_DBLE_EITHER as i8);
    }

    fn stereo_pair(current_parity: i8, next_parity: i8) -> Vec<inp_ATOM> {
        let mut current = inp_ATOM::default();
        current.valence = 2;
        current.neighbor[0] = 1;
        current.sb_ord[0] = 0;
        current.sn_ord[0] = 0;
        current.sb_parity[0] = current_parity;
        let mut next = inp_ATOM::default();
        next.valence = 2;
        next.neighbor[0] = 0;
        next.sb_ord[0] = 0;
        next.sn_ord[0] = 0;
        next.sb_parity[0] = next_parity;
        vec![current, next]
    }

    #[test]
    fn source_port__ichister__reconcilecmlincidentbondparities__line_4690() {
        let mut heap = SourceHeap::default();

        let mut ignored = inp_ATOM::default();
        ignored.valence = (MAX_NUM_STEREO_BONDS + 1) as i8;
        ignored.sb_parity[0] = 1;
        let ignored = heap.allocate(vec![ignored]).unwrap();
        let ignored_visited = heap.allocate(vec![0_i8]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, ignored, 0, -1, ignored_visited, 0,)
                .unwrap(),
            0
        );

        let no_parity = heap.allocate(vec![inp_ATOM::default()]).unwrap();
        let no_parity_visited = heap.allocate(vec![0_i8]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, no_parity, 0, -1, no_parity_visited, 0,)
                .unwrap(),
            1
        );

        let already = heap.allocate(stereo_pair(1, 1)).unwrap();
        let already_visited = heap.allocate(vec![10_i8, 0]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, already, 0, -1, already_visited, 0,)
                .unwrap(),
            2
        );

        let mut broken = inp_ATOM::default();
        broken.valence = 1;
        broken.neighbor[0] = 1;
        broken.sb_parity[0] = 1;
        let cannot_find = heap.allocate(vec![broken, inp_ATOM::default()]).unwrap();
        let cannot_find_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(
                &mut heap,
                cannot_find,
                0,
                -1,
                cannot_find_visited,
                0,
            )
            .unwrap(),
            4
        );
        assert_eq!(heap.slice(cannot_find_visited.as_const()).unwrap()[0], 10);

        let mismatch = heap.allocate(stereo_pair(3, 4)).unwrap();
        let mismatch_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, mismatch, 0, -1, mismatch_visited, 0,)
                .unwrap(),
            3
        );

        let bypass = heap.allocate(stereo_pair(3, 3)).unwrap();
        let bypass_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, bypass, 0, -1, bypass_visited, 0,).unwrap(),
            0
        );
        assert_eq!(heap.slice(bypass_visited.as_const()).unwrap(), &[20, 0]);

        let normal = heap.allocate(stereo_pair(1, 1)).unwrap();
        let normal_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, normal, 0, -1, normal_visited, 0,).unwrap(),
            0
        );
        assert_eq!(heap.slice(normal_visited.as_const()).unwrap(), &[21, 21]);

        let reconcile = heap.allocate(stereo_pair(1, 1)).unwrap();
        let reconcile_visited = heap.allocate(vec![2_i8, 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, reconcile, 0, -1, reconcile_visited, 0,)
                .unwrap(),
            0
        );
        assert_eq!(heap.slice(reconcile.as_const()).unwrap()[0].sb_parity[0], 2);
        assert_eq!(heap.slice(reconcile.as_const()).unwrap()[1].sb_parity[0], 2);
        assert_eq!(heap.slice(reconcile_visited.as_const()).unwrap(), &[22, 22]);

        let mobius = heap.allocate(stereo_pair(1, 1)).unwrap();
        let mobius_visited = heap.allocate(vec![0_i8, 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, mobius, 0, -1, mobius_visited, 0,).unwrap(),
            5
        );

        let disconnected = heap
            .allocate(stereo_pair(1 << SB_PARITY_SHFT, 1 << SB_PARITY_SHFT))
            .unwrap();
        let disconnected_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(
                &mut heap,
                disconnected,
                0,
                -1,
                disconnected_visited,
                1,
            )
            .unwrap(),
            0
        );
        assert_eq!(
            heap.slice(disconnected_visited.as_const()).unwrap(),
            &[21, 21]
        );

        let completed = heap.allocate(stereo_pair(1, 2)).unwrap();
        let completed_visited = heap.allocate(vec![0_i8, 20]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, completed, 0, -1, completed_visited, 0,)
                .unwrap(),
            0
        );
        assert_eq!(heap.slice(completed_visited.as_const()).unwrap(), &[20, 20]);

        let mut high_valence_pair = stereo_pair(1, 2);
        high_valence_pair[1].valence = (MAX_NUM_STEREO_BONDS + 1) as i8;
        let high_valence = heap.allocate(high_valence_pair).unwrap();
        let high_valence_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(
                &mut heap,
                high_valence,
                0,
                -1,
                high_valence_visited,
                0,
            )
            .unwrap(),
            0
        );
        assert_eq!(
            heap.slice(high_valence_visited.as_const()).unwrap(),
            &[20, 0]
        );
    }

    #[test]
    fn source_port__ichister__reconcileallcmlbondparities__line_4663() {
        let mut heap = SourceHeap::default();
        let empty = heap.allocate(Vec::<inp_ATOM>::new()).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, empty, 0, 0).unwrap(),
            0
        );

        let normal = heap.allocate(stereo_pair(1, 1)).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, normal, 2, 0).unwrap(),
            0
        );

        let mismatch = heap.allocate(stereo_pair(3, 4)).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, mismatch, 2, 0).unwrap(),
            3
        );

        let mut metal_pair = stereo_pair(3, 4);
        metal_pair[0].el_number = 25;
        metal_pair[1].el_number = 25;
        let metal = heap.allocate(metal_pair).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, metal, 2, 1).unwrap(),
            0
        );

        let malformed = {
            let mut atom = inp_ATOM::default();
            atom.valence = 1;
            atom.neighbor[0] = 1;
            atom.sb_parity[0] = 1;
            heap.allocate(vec![atom, inp_ATOM::default()]).unwrap()
        };
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, malformed, 2, 0).unwrap(),
            4
        );

        heap.fail_after_allocations(0);
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, normal, 2, 0).unwrap(),
            -1
        );
    }

    #[test]
    fn source_port__ichister__get_opposite_sb_atom__line_4861() {
        let mut heap = SourceHeap::default();
        let mut start = inp_ATOM::default();
        start.neighbor[0] = 1;
        let mut opposite = inp_ATOM::default();
        opposite.sb_parity[0] = 1;
        opposite.sb_ord[0] = 2;
        opposite.neighbor[2] = 0;
        let direct = heap.allocate(vec![start.clone(), opposite]).unwrap();
        let (mut next, mut reverse, mut parity) = (-1, -1, -1);
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                direct,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            1
        );
        assert_eq!((next, reverse, parity), (1, 2, 0));

        let mut wrong = inp_ATOM::default();
        wrong.sb_parity[0] = 1;
        wrong.sb_ord[0] = 0;
        wrong.neighbor[0] = 7;
        let no_match = heap.allocate(vec![start.clone(), wrong]).unwrap();
        let (mut next, mut reverse, mut parity) = (-2, -3, -4);
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                no_match,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            0
        );
        assert_eq!((next, reverse, parity), (-2, -3, -4));

        let mut middle = inp_ATOM::default();
        middle.neighbor[0] = 0;
        middle.neighbor[1] = 2;
        middle.valence = 2;
        middle.chem_bonds_valence = (2 * BOND_TYPE_DOUBLE) as i8;
        let mut far = inp_ATOM::default();
        far.sb_parity[0] = 2;
        far.sb_ord[0] = 1;
        far.neighbor[1] = 1;
        let chain = heap.allocate(vec![start, middle, far]).unwrap();
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                chain,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            2
        );
        assert_eq!((next, reverse, parity), (2, 1, 0));

        let mut loop_atom = inp_ATOM::default();
        loop_atom.neighbor[0] = 1;
        loop_atom.neighbor[1] = 0;
        loop_atom.valence = 2;
        loop_atom.chem_bonds_valence = (2 * BOND_TYPE_DOUBLE) as i8;
        let cycle = heap.allocate(vec![loop_atom.clone(), loop_atom]).unwrap();
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                cycle,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            0
        );
    }
}
