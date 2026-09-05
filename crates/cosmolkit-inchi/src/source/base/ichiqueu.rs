use crate::source::base::{
    ichi_bns::bExistsAnyAltPath,
    ichitaut::{AddAtom2DA, AddAtom2num, is_centerpoint_elem_strict, nGetEndpointInfo},
    util::get_endpoint_valence,
};
use crate::source_types::{
    ALT_PATH_MODE_TAUTOM, AT_RANK, BOND_ALT_13, BOND_ALT12NS, BOND_ALTERN, BOND_DOUBLE, BOND_MARK_ALL, BOND_SINGLE,
    BOND_TAUTOM, BOND_TRIPLE, BalancedNetworkData, BalancedNetworkStructure, CANON_GLOBALS, DFS_PATH, ENDPOINT_INFO,
    SourceHeap, SourceHeapError, SourceMutPointer, T_BONDPOS, T_ENDPOINT, clock_t, inp_ATOM,
};

pub(crate) fn are_alt_bonds(bonds: &[u8], len: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:778 are_alt_bonds
    // INCHI✔️✔️: int are_alt_bonds( U_CHAR *bonds, int len )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     U_CHAR next_bond;
    // INCHI✔️✔️:     int           i, bAnyBond, bTautBondPresent = BOND_ALTERN;
    // INCHI✔️✔️:     if (len < 2 || bonds[0] == BOND_TRIPLE || bonds[0] == BOND_ALT_13)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     next_bond = bonds[0] == BOND_SINGLE ? BOND_DOUBLE : bonds[0] == BOND_DOUBLE ? BOND_SINGLE : 0; /* djb-rwth: removing redundant code; ignoring LLVM warning: possible presence of global variables */
    // INCHI✔️✔️:     if (bonds[0] == BOND_TAUTOM)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         bTautBondPresent = BOND_TAUTOM;
    // INCHI✔️✔️:         next_bond = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         next_bond = bonds[0] == BOND_SINGLE ? BOND_DOUBLE : bonds[0] == BOND_DOUBLE ? BOND_SINGLE : 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 1; i < len; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (bonds[i] == BOND_TAUTOM)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bTautBondPresent = BOND_TAUTOM;
    // INCHI✔️✔️:             bAnyBond = 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bAnyBond = ( bonds[i] == BOND_ALTERN || bonds[i] == BOND_ALT12NS );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (next_bond)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (bonds[i] == next_bond || bAnyBond)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 next_bond = ( next_bond == BOND_SINGLE ) ? BOND_DOUBLE : BOND_SINGLE;
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (bonds[i] == BOND_SINGLE)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 next_bond = BOND_DOUBLE;
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (bonds[i] == BOND_DOUBLE)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     next_bond = BOND_SINGLE;
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (!bAnyBond)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         return 0;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return !next_bond ? bTautBondPresent : ( next_bond == BOND_SINGLE )
    // INCHI✔️✔️:                             ? BOND_DOUBLE : BOND_SINGLE; /* bond to the end atom */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: are_alt_bonds

    if len < 2 {
        return Ok(0);
    }
    let len = usize::try_from(len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let bonds = bonds.get(..len).ok_or(SourceHeapError::PointerOutOfBounds)?;
    if bonds[0] == BOND_TRIPLE as u8 || bonds[0] == BOND_ALT_13 as u8 {
        return Ok(0);
    }

    let mut bTautBondPresent = BOND_ALTERN as i32;
    let mut next_bond = if bonds[0] == BOND_SINGLE as u8 {
        BOND_DOUBLE as u8
    } else if bonds[0] == BOND_DOUBLE as u8 {
        BOND_SINGLE as u8
    } else {
        0
    };
    if bonds[0] == BOND_TAUTOM as u8 {
        bTautBondPresent = BOND_TAUTOM as i32;
        next_bond = 0;
    } else {
        next_bond = if bonds[0] == BOND_SINGLE as u8 {
            BOND_DOUBLE as u8
        } else if bonds[0] == BOND_DOUBLE as u8 {
            BOND_SINGLE as u8
        } else {
            0
        };
    }

    for &bond in &bonds[1..] {
        let bAnyBond = if bond == BOND_TAUTOM as u8 {
            bTautBondPresent = BOND_TAUTOM as i32;
            true
        } else {
            bond == BOND_ALTERN as u8 || bond == BOND_ALT12NS as u8
        };
        if next_bond != 0 {
            if bond == next_bond || bAnyBond {
                next_bond = if next_bond == BOND_SINGLE as u8 {
                    BOND_DOUBLE as u8
                } else {
                    BOND_SINGLE as u8
                };
                continue;
            }
            return Ok(0);
        }
        if bond == BOND_SINGLE as u8 {
            next_bond = BOND_DOUBLE as u8;
            continue;
        }
        if bond == BOND_DOUBLE as u8 {
            next_bond = BOND_SINGLE as u8;
            continue;
        }
        if !bAnyBond {
            return Ok(0);
        }
    }

    Ok(if next_bond == 0 {
        bTautBondPresent
    } else if next_bond == BOND_SINGLE as u8 {
        BOND_DOUBLE as i32
    } else {
        BOND_SINGLE as i32
    })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AddBondsPos(
    atom: &[inp_ATOM],
    BondPosTmp: &mut [T_BONDPOS],
    nNumBondPosTmp: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    mut nNumBondPos: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:846 AddBondsPos
    // INCHI✔️✔️: int AddBondsPos( inp_ATOM *atom,
    // INCHI✔️✔️:                  T_BONDPOS *BondPosTmp,
    // INCHI✔️✔️:                  int nNumBondPosTmp,
    // INCHI✔️✔️:                  T_BONDPOS *BondPos,
    // INCHI✔️✔️:                  int nMaxNumBondPos,
    // INCHI✔️✔️:                  int nNumBondPos )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j, k, cur_at, nxt_at;
    // INCHI✔️✔️:     /*  add opposite direction bonds to BondPosTmp */
    // INCHI✔️✔️:     for (j = 0; j < nNumBondPosTmp; j += 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         cur_at = BondPosTmp[j].nAtomNumber;
    // INCHI✔️✔️:         nxt_at = atom[cur_at].neighbor[(int) BondPosTmp[j].neighbor_index];
    // INCHI✔️✔️:         for (k = 0; k < atom[nxt_at].valence; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (cur_at == atom[nxt_at].neighbor[k])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 BondPosTmp[j + 1].nAtomNumber = nxt_at;
    // INCHI✔️✔️:                 BondPosTmp[j + 1].neighbor_index = k;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     /*  add new tautomeric bonds */
    // INCHI✔️✔️:     for (j = 0; j < nNumBondPosTmp; j += 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0; i < nNumBondPos; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if ((BondPos[i].nAtomNumber == BondPosTmp[j].nAtomNumber &&
    // INCHI✔️✔️:                 BondPos[i].neighbor_index == BondPosTmp[j].neighbor_index) ||
    // INCHI✔️✔️:                 (BondPos[i].nAtomNumber == BondPosTmp[j + 1].nAtomNumber &&
    // INCHI✔️✔️:                 BondPos[i].neighbor_index == BondPosTmp[j + 1].neighbor_index))  /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 break; /*  bond has already been added */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (i == nNumBondPos)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (i > nMaxNumBondPos)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return -1; /*  overflow */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             BondPos[nNumBondPos++] = BondPosTmp[j];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nNumBondPos;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddBondsPos

    let temporary_count = usize::try_from(nNumBondPosTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if temporary_count > BondPosTmp.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut j = 0_usize;
    while j < temporary_count {
        let first = BondPosTmp.get(j).cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let cur_at = usize::from(first.nAtomNumber);
        let neighbor_index = usize::from(first.neighbor_index);
        let current = atom.get(cur_at).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let nxt_at = usize::from(
            *current
                .neighbor
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let next = atom.get(nxt_at).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = usize::try_from(next.valence).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        for k in 0..valence {
            if cur_at == usize::from(*next.neighbor.get(k).ok_or(SourceHeapError::PointerOutOfBounds)?) {
                let opposite = BondPosTmp.get_mut(j + 1).ok_or(SourceHeapError::PointerOutOfBounds)?;
                opposite.nAtomNumber = nxt_at as AT_RANK;
                opposite.neighbor_index = k as AT_RANK;
                break;
            }
        }
        j = j.wrapping_add(2);
    }

    let mut j = 0_usize;
    while j < temporary_count {
        let existing_count = usize::try_from(nNumBondPos).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if existing_count > BondPos.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let first = BondPosTmp.get(j).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let opposite = BondPosTmp.get(j + 1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut i = 0_usize;
        while i < existing_count {
            let existing = BondPos.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if (existing.nAtomNumber == first.nAtomNumber && existing.neighbor_index == first.neighbor_index)
                || (existing.nAtomNumber == opposite.nAtomNumber && existing.neighbor_index == opposite.neighbor_index)
            {
                break;
            }
            i += 1;
        }
        if i == existing_count {
            if i32::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)? > nMaxNumBondPos {
                return Ok(-1);
            }
            *BondPos
                .get_mut(existing_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = first.clone();
            nNumBondPos = nNumBondPos.wrapping_add(1);
        }
        j = j.wrapping_add(2);
    }

    Ok(nNumBondPos)
}

#[allow(non_snake_case)]
pub(crate) fn AddEndPoints(
    EndPointTmp: &[T_ENDPOINT],
    nNumNewEndPoint: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    mut nNumEndPoint: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:897 AddEndPoints
    // INCHI✔️✔️: int AddEndPoints( T_ENDPOINT *EndPointTmp,
    // INCHI✔️✔️:                   int nNumNewEndPoint,
    // INCHI✔️✔️:                   T_ENDPOINT *EndPoint,
    // INCHI✔️✔️:                   int nMaxNumEndPoint,
    // INCHI✔️✔️:                   int nNumEndPoint )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j;
    // INCHI✔️✔️:     /*  add new endpoints */
    // INCHI✔️✔️:     for (j = 0; j < nNumNewEndPoint; j++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0; i < nNumEndPoint; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (EndPoint[i].nAtomNumber == EndPointTmp[j].nAtomNumber)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (i == nNumEndPoint)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (i > nMaxNumEndPoint)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return -1; /*  overflow */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             EndPoint[nNumEndPoint++] = EndPointTmp[j];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nNumEndPoint;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddEndPoints

    let new_count = usize::try_from(nNumNewEndPoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let new_endpoints = EndPointTmp
        .get(..new_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for candidate in new_endpoints {
        let existing_count = usize::try_from(nNumEndPoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if existing_count > EndPoint.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let mut i = 0_usize;
        while i < existing_count {
            if EndPoint.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?.nAtomNumber == candidate.nAtomNumber {
                break;
            }
            i += 1;
        }
        if i == existing_count {
            if i32::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)? > nMaxNumEndPoint {
                return Ok(-1);
            }
            *EndPoint
                .get_mut(existing_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = candidate.clone();
            nNumEndPoint = nNumEndPoint.wrapping_add(1);
        }
    }
    Ok(nNumEndPoint)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Check7MembTautRing(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &[DFS_PATH],
    nLenDfsPath: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:951 Check7MembTautRing
    // INCHI✔️❌: int Check7MembTautRing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                         inp_ATOM *atom,
    // INCHI✔️❌:                         DFS_PATH *DfsPath,
    // INCHI✔️❌:                         int nLenDfsPath,
    // INCHI✔️❌:                         int nStartAtomNeighbor,
    // INCHI✔️❌:                         int nStartAtomNeighbor2,
    // INCHI✔️❌:                         int nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                         T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                         int nMaxNumEndPoint,
    // INCHI✔️❌:                         T_BONDPOS  *BondPos,
    // INCHI✔️❌:                         int nMaxNumBondPos,
    // INCHI✔️❌:                         int *pnNumEndPoint,
    // INCHI✔️❌:                         int *pnNumBondPos,
    // INCHI✔️❌:                         struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                         struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                         int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌: #define PATH_LEN 8
    // INCHI✔️❌:
    // INCHI✔️❌:     int i, j, k, /*m,*/ nNumEndPoint, nNumEndPointTmp, nNumBondPos, nNumBondPosTmp;
    // INCHI✔️❌:     int endpoint, /*nMobile, nMobile1, nMobile2,*/ o1_at, o2_at;
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     U_CHAR path_bonds[PATH_LEN + 1], bond_type;
    // INCHI✔️❌:     T_ENDPOINT EndPointTmp[2];
    // INCHI✔️❌:     T_BONDPOS  BondPosTmp[2 * PATH_LEN];
    // INCHI✔️❌:     ENDPOINT_INFO eif1, eif2;
    // INCHI✔️❌:     int nErr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nLenDfsPath + 2 > PATH_LEN)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  too long path */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nLenDfsPath != 6 && nLenDfsPath != 4)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumBondPos = *pnNumBondPos;
    // INCHI✔️❌:     nNumEndPoint = *pnNumEndPoint;
    // INCHI✔️❌:     nNumBondPosTmp = 0;
    // INCHI✔️❌:     nNumEndPointTmp = 0;
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     o1_at = atom[(int) DfsPath[1].at_no].neighbor[nStartAtomNeighborNeighbor];
    // INCHI✔️❌:     o2_at = atom[(int) DfsPath[0].at_no].neighbor[nStartAtomNeighbor2];
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     nMobile1 = (atom[o1_at].charge == -1) + atom[o1_at].num_H;
    // INCHI✔️❌:     nMobile2 = (atom[o2_at].charge == -1) + atom[o2_at].num_H;
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (!nGetEndpointInfo( atom, o1_at, &eif1 ) ||
    // INCHI✔️❌:         !nGetEndpointInfo( atom, o2_at, &eif2 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  save endpoints */
    // INCHI✔️❌:     for (j = 0; j < 2; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         endpoint = j ? o2_at : o1_at;
    // INCHI✔️❌:         if (!atom[endpoint].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
    // INCHI✔️❌:             AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
    // INCHI✔️❌:             /*
    // INCHI✔️❌:                    nMobile  = j? nMobile2 : nMobile1;
    // INCHI✔️❌:                } else {
    // INCHI✔️❌:                    nMobile  = 0;
    // INCHI✔️❌:                }
    // INCHI✔️❌:                if ( nMobile ) {
    // INCHI✔️❌:                    EndPointTmp[nNumEndPointTmp].num[1] = (atom[endpoint].charge == -1);
    // INCHI✔️❌:                    EndPointTmp[nNumEndPointTmp].num[0] = nMobile;
    // INCHI✔️❌:                    for ( m = 0; m < T_NUM_ISOTOPIC; m ++ ) {
    // INCHI✔️❌:                        EndPointTmp[nNumEndPointTmp].num[T_NUM_NO_ISOTOPIC+m] = atom[endpoint].num_iso_H[NUM_H_ISOTOPES-m-1];
    // INCHI✔️❌:                    }
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nEquNumber = 0;
    // INCHI✔️❌:         nNumEndPointTmp++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  extract bonds */
    // INCHI✔️❌:     k = (int) DfsPath[1].at_no;
    // INCHI✔️❌:     bond_type = ( atom[k].bond_type[nStartAtomNeighborNeighbor] & ~BOND_MARK_ALL );
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:     bond_type = ACTUAL_ORDER( pBNS, k, nStartAtomNeighborNeighbor, bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:     path_bonds[0] = bond_type;
    // INCHI✔️❌:     if (REPLACE_THE_BOND( bond_type ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         BondPosTmp[nNumBondPosTmp].nAtomNumber = k;
    // INCHI✔️❌:         BondPosTmp[nNumBondPosTmp].neighbor_index = nStartAtomNeighborNeighbor;
    // INCHI✔️❌:         nNumBondPosTmp += 2;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 1; i <= nLenDfsPath; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bond_type = DfsPath[i].bond_type;
    // INCHI✔️❌:         path_bonds[i] = bond_type;
    // INCHI✔️❌:         if (REPLACE_THE_BOND( bond_type ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[i].at_no;
    // INCHI✔️❌:             BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[i].bond_pos;
    // INCHI✔️❌:             nNumBondPosTmp += 2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     bond_type = ( atom[(int) DfsPath[0].at_no].bond_type[nStartAtomNeighbor2] & ~BOND_MARK_ALL );
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:     bond_type = ACTUAL_ORDER( pBNS, (int) DfsPath[0].at_no, nStartAtomNeighbor2, bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:     path_bonds[i++] = bond_type;
    // INCHI✔️❌:     if (REPLACE_THE_BOND( bond_type ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[0].at_no;
    // INCHI✔️❌:         BondPosTmp[nNumBondPosTmp].neighbor_index = nStartAtomNeighbor2;
    // INCHI✔️❌:         nNumBondPosTmp += 2;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!are_alt_bonds( path_bonds, i ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* path_bonds is from at_n1 to at_n2 */
    // INCHI✔️❌:     if (!( j = are_alt_bonds( path_bonds, i ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* j is a bond type of the last bond to o2_at, the first bond from o1_at is 2-j if j=1 or 2 */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* single bond at o2_at: it should have a mobile atom, o1_at should not */
    // INCHI✔️❌:     if ((j == BOND_SINGLE && ( (!atom[o2_at].endpoint && !eif2.cDonor) || (!atom[o1_at].endpoint && !eif1.cAcceptor) )) ||
    // INCHI✔️❌:         /* double bond at o2_at: it should not have a mobile atom, o1_at should */
    // INCHI✔️❌:         (j == BOND_DOUBLE && ( (!atom[o2_at].endpoint && !eif2.cAcceptor) || (!atom[o1_at].endpoint && !eif1.cDonor) ))) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* bond pattern does not fit */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    // INCHI✔️❌:     nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *pnNumBondPos = nNumBondPos;
    // INCHI✔️❌:             *pnNumEndPoint = nNumEndPoint;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* finally check whether the bonds allow moving the hydrogens */
    // INCHI✔️❌:         if (( atom[o1_at].endpoint != atom[o2_at].endpoint || !atom[o1_at].endpoint ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌:             nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, o1_at, o2_at, ALT_PATH_MODE_TAUTOM );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (nErr <= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return nErr;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: #undef PATH_LEN
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Check7MembTautRing
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: Check7MembTautRing
    // INCHI✔️❌: #define PATH_LEN 8
    // INCHI✔️❌: #define FIX_BOND23_IN_TAUT 0
    // INCHI✔️❌: #define REPLACE_ALT_WITH_TAUT 1
    // INCHI✔️❌: #define REPLACE_THE_BOND(X) ( (X) == BOND_SINGLE || (X) == BOND_DOUBLE || (X) == BOND_ALTERN || (X) == BOND_ALT12NS )
    // INCHI✔️❌: #define BOND_MARK_ALL 0xf0
    // END INCHI ACTIVE MACRO CONFIGURATION: Check7MembTautRing

    const PATH_LEN: i32 = 8;
    if nLenDfsPath.wrapping_add(2) > PATH_LEN {
        return Ok(-1);
    }
    if nLenDfsPath != 6 && nLenDfsPath != 4 {
        return Ok(-1);
    }

    let mut nNumBondPos = *pnNumBondPos;
    let mut nNumEndPoint = *pnNumEndPoint;
    let mut nNumBondPosTmp = 0_i32;
    let mut nNumEndPointTmp = 0_i32;
    let mut ret = 0_i32;
    let mut o1_at = 0_i32;
    let mut o2_at = 0_i32;

    {
        let atoms = heap.slice(atom.as_const())?;
        let path0 = DfsPath.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let path1 = DfsPath.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let path1_atom = atoms
            .get(usize::from(path1.at_no))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        o1_at = i32::from(
            *path1_atom
                .neighbor
                .get(usize::try_from(nStartAtomNeighborNeighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let path0_atom = atoms
            .get(usize::from(path0.at_no))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        o2_at = i32::from(
            *path0_atom
                .neighbor
                .get(usize::try_from(nStartAtomNeighbor2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );

        let mut eif1 = ENDPOINT_INFO::default();
        let mut eif2 = ENDPOINT_INFO::default();
        if nGetEndpointInfo(atoms, o1_at, &mut eif1) == 0 || nGetEndpointInfo(atoms, o2_at, &mut eif2) == 0 {
            return Ok(0);
        }

        let mut EndPointTmp: [T_ENDPOINT; 2] = std::array::from_fn(|_| T_ENDPOINT::default());
        for endpoint in [o1_at, o2_at] {
            let temporary_index = usize::try_from(nNumEndPointTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let endpoint_atom = atoms
                .get(usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if endpoint_atom.endpoint == 0 {
                AddAtom2num(&mut EndPointTmp[temporary_index].num, atoms, endpoint, 2)?;
                AddAtom2DA(&mut EndPointTmp[temporary_index].num_DA, atoms, endpoint, 2)?;
            } else {
                EndPointTmp[temporary_index] = T_ENDPOINT::default();
            }
            EndPointTmp[temporary_index].nAtomNumber = endpoint as AT_RANK;
            EndPointTmp[temporary_index].nGroupNumber = endpoint_atom.endpoint;
            EndPointTmp[temporary_index].nEquNumber = 0;
            nNumEndPointTmp = nNumEndPointTmp.wrapping_add(1);
        }

        let mut path_bonds = [0_u8; PATH_LEN as usize + 1];
        let mut BondPosTmp: [T_BONDPOS; 2 * PATH_LEN as usize] = std::array::from_fn(|_| T_BONDPOS::default());
        let k = i32::from(path1.at_no);
        let first_neighbor_index =
            usize::try_from(nStartAtomNeighborNeighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut bond_type = *path1_atom
            .bond_type
            .get(first_neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            & !(BOND_MARK_ALL as u8);
        path_bonds[0] = bond_type;
        if matches!(
            bond_type,
            x if x == BOND_SINGLE as u8
                || x == BOND_DOUBLE as u8
                || x == BOND_ALTERN as u8
                || x == BOND_ALT12NS as u8
        ) {
            let index = usize::try_from(nNumBondPosTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            BondPosTmp[index].nAtomNumber = k as AT_RANK;
            BondPosTmp[index].neighbor_index = nStartAtomNeighborNeighbor as AT_RANK;
            nNumBondPosTmp = nNumBondPosTmp.wrapping_add(2);
        }

        let mut i = 1_i32;
        while i <= nLenDfsPath {
            let path = DfsPath
                .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            bond_type = path.bond_type;
            path_bonds[usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = bond_type;
            if matches!(
                bond_type,
                x if x == BOND_SINGLE as u8
                    || x == BOND_DOUBLE as u8
                    || x == BOND_ALTERN as u8
                    || x == BOND_ALT12NS as u8
            ) {
                let index = usize::try_from(nNumBondPosTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                BondPosTmp[index].nAtomNumber = path.at_no;
                BondPosTmp[index].neighbor_index = path.bond_pos as AT_RANK;
                nNumBondPosTmp = nNumBondPosTmp.wrapping_add(2);
            }
            i = i.wrapping_add(1);
        }

        bond_type = *path0_atom
            .bond_type
            .get(usize::try_from(nStartAtomNeighbor2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            & !(BOND_MARK_ALL as u8);
        path_bonds[usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = bond_type;
        i = i.wrapping_add(1);
        if matches!(
            bond_type,
            x if x == BOND_SINGLE as u8
                || x == BOND_DOUBLE as u8
                || x == BOND_ALTERN as u8
                || x == BOND_ALT12NS as u8
        ) {
            let index = usize::try_from(nNumBondPosTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            BondPosTmp[index].nAtomNumber = path0.at_no;
            BondPosTmp[index].neighbor_index = nStartAtomNeighbor2 as AT_RANK;
            nNumBondPosTmp = nNumBondPosTmp.wrapping_add(2);
        }

        if are_alt_bonds(&path_bonds, i)? == 0 {
            return Ok(0);
        }
        let j = are_alt_bonds(&path_bonds, i)?;
        if j == 0 {
            return Ok(0);
        }

        let o1 = atoms
            .get(usize::try_from(o1_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let o2 = atoms
            .get(usize::try_from(o2_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if (j == BOND_SINGLE as i32
            && ((o2.endpoint == 0 && eif2.cDonor == 0) || (o1.endpoint == 0 && eif1.cAcceptor == 0)))
            || (j == BOND_DOUBLE as i32
                && ((o2.endpoint == 0 && eif2.cAcceptor == 0) || (o1.endpoint == 0 && eif1.cDonor == 0)))
        {
            return Ok(0);
        }

        nNumBondPos = AddBondsPos(
            atoms,
            &mut BondPosTmp,
            nNumBondPosTmp,
            BondPos,
            nMaxNumBondPos,
            nNumBondPos,
        )?;
        nNumEndPoint = AddEndPoints(&EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint)?;

        if nNumBondPos >= 0 && nNumEndPoint >= 0 && (nNumBondPos > *pnNumBondPos || nNumEndPoint > *pnNumEndPoint) {
            ret = 1;
            *pnNumBondPos = nNumBondPos;
            *pnNumEndPoint = nNumEndPoint;
        }
    }

    if ret != 0 {
        let atoms = heap.slice(atom.as_const())?;
        let o1_endpoint = atoms
            .get(usize::try_from(o1_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint;
        let o2_endpoint = atoms
            .get(usize::try_from(o2_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint;
        if o1_endpoint != o2_endpoint || o1_endpoint == 0 {
            let nErr = bExistsAnyAltPath(
                heap,
                pCG,
                pBNS,
                pBD,
                atom,
                num_atoms,
                o1_at,
                o2_at,
                ALT_PATH_MODE_TAUTOM as i32,
                clock_result,
            )?;
            if nErr <= 0 {
                return Ok(nErr);
            }
        }
    }

    Ok(ret)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Check6MembTautRing(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &[DFS_PATH],
    nLenDfsPath: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1142 Check6MembTautRing
    // INCHI✔️❌: int Check6MembTautRing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                         inp_ATOM *atom,
    // INCHI✔️❌:                         DFS_PATH *DfsPath,
    // INCHI✔️❌:                         int nLenDfsPath,
    // INCHI✔️❌:                         int nStartAtomNeighbor,
    // INCHI✔️❌:                         int nStartAtomNeighbor2,
    // INCHI✔️❌:                         int nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                         T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                         int nMaxNumEndPoint,
    // INCHI✔️❌:                         T_BONDPOS  *BondPos,
    // INCHI✔️❌:                         int nMaxNumBondPos,
    // INCHI✔️❌:                         int *pnNumEndPoint,
    // INCHI✔️❌:                         int *pnNumBondPos,
    // INCHI✔️❌:                         struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                         struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                         int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌: #define PATH_LEN 4
    // INCHI✔️❌:     int i, j, k, /*m,*/ nNumBondPos, nNumEndPoint, ept, eptn;
    // INCHI✔️❌:     int nNumEndPointTmp, nNumBondPosTmp, o_at = 0, ret;
    // INCHI✔️❌:     /* int num_taut_endpoints, num_H; */
    // INCHI✔️❌:     int middle_pos;
    // INCHI✔️❌:     int nMobile, endpoint, endpoint_valence, chem_bonds_valence;
    // INCHI✔️❌:     int nMobile1, endpoint_valence1;  /*  o_at */
    // INCHI✔️❌:     int nMobile2, endpoint_valence2;  /*  n_at */
    // INCHI✔️❌:     int nxt_at;
    // INCHI✔️❌:     int n_at;
    // INCHI✔️❌:     U_CHAR path_bonds[2][PATH_LEN + 1], bond_type;
    // INCHI✔️❌:     T_ENDPOINT EndPointTmp[2];
    // INCHI✔️❌:     T_BONDPOS  BondPosTmp[4 * PATH_LEN];
    // INCHI✔️❌:     ENDPOINT_INFO eif1, eif2;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nStartAtomNeighbor >= 0 || nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nLenDfsPath != 5)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumBondPos = *pnNumBondPos;
    // INCHI✔️❌:     nNumEndPoint = *pnNumEndPoint;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     nNumEndPointTmp = 0;
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:     for (ept = 0; ept < 2; ept++) /* djb-rwth: initialisation needed for num array */
    // INCHI✔️❌:         for (eptn = 0; eptn < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; eptn++)
    // INCHI✔️❌:             EndPointTmp[ept].num[eptn] = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     n_at = (int) DfsPath[0].at_no;   /*  -N= or -NH- atom */
    // INCHI✔️❌:     nxt_at = DfsPath[middle_pos = ( nLenDfsPath + 1 ) / 2].at_no;  /*  must have tautomeric neighbor -OH or =O or -NH2 or =NH */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (atom[nxt_at].valence != 3
    // INCHI✔️❌: #if ( TAUT_RINGS_ATTACH_CHAIN == 1 )
    // INCHI✔️❌:         || !atom[nxt_at].bCutVertex
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < atom[nxt_at].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         o_at = atom[nxt_at].neighbor[i];
    // INCHI✔️❌:         if (o_at != DfsPath[middle_pos - 1].at_no && o_at != DfsPath[middle_pos + 1].at_no)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break; /*  >=O or />-OH has been found */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i == atom[nxt_at].valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /*  no neighboring atom >=O or />-OH */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     bond_type = ( atom[nxt_at].bond_type[i] & ~BOND_MARK_ALL );
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:     bond_type = ACTUAL_ORDER( pBNS, nxt_at, i, bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:     if (bond_type != BOND_SINGLE &&
    // INCHI✔️❌:         bond_type != BOND_DOUBLE &&
    // INCHI✔️❌:         bond_type != BOND_TAUTOM &&
    // INCHI✔️❌:         bond_type != BOND_ALT12NS &&
    // INCHI✔️❌:         bond_type != BOND_ALTERN)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  check whether the two atoms already belong to one tautomeric group */
    // INCHI❌❌: #if ( TAUT_IGNORE_EQL_ENDPOINTS == 1 )
    // INCHI❌❌:     if (atom[n_at].endpoint && atom[n_at].endpoint == atom[o_at].endpoint)
    // INCHI❌❌:     {
    // INCHI❌❌:         return 0;
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:     /*  check =O valence; must be 2 for O, S, Se or 3 for N */
    // INCHI✔️❌:     if (!( endpoint_valence1 = nGetEndpointInfo( atom, o_at, &eif1 ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /*  n_at has been checked in MarkTautomerGroups(...) */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         if ( 2 != endpoint_valence1 )
    // INCHI✔️❌:             return 0; // accept only O, S, Se
    // INCHI✔️❌:     */
    // INCHI✔️❌:     /*  check hydrogens/endpoints */
    // INCHI✔️❌:     nMobile1 = atom[o_at].num_H + ( atom[o_at].charge == -1 );
    // INCHI✔️❌:     if (bond_type == BOND_SINGLE && !eif1.cDonor && !atom[o_at].endpoint)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* not needed since nGetEndpointInfo returned non-zero
    // INCHI✔️❌:     if ( nMobile1 + atom[o_at].chem_bonds_valence != endpoint_valence1 )
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!( endpoint_valence2 = nGetEndpointInfo( atom, n_at, &eif2 ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* should not happen here */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     nMobile2 = atom[n_at].num_H + ( atom[n_at].charge == -1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     nMobile = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  can mobile group move from o_at to n_at? */
    // INCHI✔️❌:     nMobile += ( atom[o_at].endpoint || eif1.cDonor ) &&  /*  from o_at */
    // INCHI✔️❌:         bond_type != BOND_DOUBLE &&
    // INCHI✔️❌:         ( atom[n_at].endpoint ||                   /*  to n_at */
    // INCHI✔️❌:             eif2.cNeutralBondsValence > atom[n_at].valence );
    // INCHI✔️❌:     /*  can mobile group move from n_at to o_at? */
    // INCHI✔️❌:     nMobile += ( atom[n_at].endpoint || eif2.cDonor ) && /*  from n_at */
    // INCHI✔️❌:         ( atom[o_at].endpoint ||          /*  to o_at */
    // INCHI✔️❌:             eif1.cNeutralBondsValence > atom[o_at].valence ) &&
    // INCHI✔️❌:         bond_type != BOND_SINGLE;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!nMobile)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     num_H = atom[n_at].num_H + atom[o_at].num_H;
    // INCHI✔️❌:     num_taut_endpoints = (0!=atom[n_at].endpoint) + (0!=atom[o_at].endpoint); // if O, N already are endpoints
    // INCHI✔️❌:     if ( num_H != 1 && num_taut_endpoints != 2 && !(num_H==2 && num_taut_endpoints >= 1) ) {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:     /*  extract -OH bond */
    // INCHI✔️❌:     nNumBondPosTmp = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     path_bonds[0][0] = path_bonds[1][0] = bond_type;
    // INCHI✔️❌:     if (REPLACE_THE_BOND( bond_type ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         BondPosTmp[nNumBondPosTmp].nAtomNumber = nxt_at; /*  accumulate bonds to be */
    // INCHI✔️❌:         BondPosTmp[nNumBondPosTmp].neighbor_index = i;   /*  marked as tautomeric */
    // INCHI✔️❌:         nNumBondPosTmp += 2; /*  leave room for the same bond in the opposite direction */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  extract other bonds */
    // INCHI✔️❌:     /* path_bonds[] contents:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                     O              OH            OH
    // INCHI✔️❌:                     ||             |             |
    // INCHI✔️❌:                    /  \          //  \          /  \\
    // INCHI✔️❌:                   ||   ||  <-->  |   ||  <-->  ||   |
    // INCHI✔️❌:                    \  /          \\  /          \  //
    // INCHI✔️❌:                     NH             N             N
    // INCHI✔️❌:
    // INCHI✔️❌:         path[0]:  O=NH-=-      OH-N...         OH.N...
    // INCHI✔️❌:         path[1]   O=NH-=-      OH-N...         OH.N...
    // INCHI✔️❌:                  bonds are    all bonds       all bonds
    // INCHI✔️❌:                  single and   are either      are either
    // INCHI✔️❌:                  double       alt or taut     alt or taut
    // INCHI✔️❌:     */
    // INCHI✔️❌:     for (j = 0; j < middle_pos; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < 2; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* k = i? j : middle_pos-1-j; */
    // INCHI✔️❌:             k = i ? middle_pos + j : middle_pos - 1 - j;
    // INCHI✔️❌:             /*  i=0: from O neighbor i=0: down to N, i=1: up to N */
    // INCHI✔️❌:             bond_type = DfsPath[k].bond_type;
    // INCHI✔️❌:
    // INCHI✔️❌:             path_bonds[i][j + 1] = bond_type;
    // INCHI✔️❌:             if (REPLACE_THE_BOND( bond_type ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[k].at_no;       /*  accumulate bonds to be */
    // INCHI✔️❌:                 BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[k].bond_pos; /*  marked as tautomeric */
    // INCHI✔️❌:                 nNumBondPosTmp += 2;   /*  leave room for the same bond in the opposite direction */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!are_alt_bonds( path_bonds[0], middle_pos + 1 ) || !are_alt_bonds( path_bonds[1], middle_pos + 1 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* finally check whether the bonds allow moving the hydrogens */
    // INCHI✔️❌:     if (( atom[o_at].endpoint != atom[n_at].endpoint || !atom[o_at].endpoint ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int nErr;
    // INCHI✔️❌:         nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, n_at, o_at, ALT_PATH_MODE_TAUTOM );
    // INCHI✔️❌:         if (nErr <= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return nErr;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  save endpoints */
    // INCHI✔️❌:     for (j = 0; j < 2; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         endpoint = j ? n_at :      /*  =N-  2 */
    // INCHI✔️❌:             o_at;       /*  -OH  1 */
    // INCHI✔️❌:         if (!atom[endpoint].endpoint)
    // INCHI✔️❌:         { /* not a known endpoint */
    // INCHI✔️❌:             endpoint_valence = j ? endpoint_valence2 : endpoint_valence1;
    // INCHI✔️❌:             chem_bonds_valence = j ? eif2.cNeutralBondsValence : eif1.cNeutralBondsValence;
    // INCHI✔️❌:             /* endpoint_valence = get_endpoint_valence( atom[endpoint].el_number ); */
    // INCHI✔️❌:             nMobile = j ? nMobile2 : nMobile1;
    // INCHI✔️❌:             /* nMobile  = (atom[endpoint].charge == -1) + atom[endpoint].num_H; */
    // INCHI✔️❌:             /* if ( nMobile + atom[endpoint].chem_bonds_valence != endpoint_valence ) -- fixed 02-06-2003*/
    // INCHI✔️❌:             if (nMobile + chem_bonds_valence != endpoint_valence)
    // INCHI✔️❌:                 return 0; /*  abnormal endpoint valence; ignore. */
    // INCHI✔️❌:             AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
    // INCHI✔️❌:             AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
    // INCHI✔️❌:             /*
    // INCHI✔️❌:                         EndPointTmp[nNumEndPointTmp].num[1] = (atom[endpoint].charge == -1);
    // INCHI✔️❌:                         EndPointTmp[nNumEndPointTmp].num[0] = nMobile;
    // INCHI✔️❌:                         for ( m = 0; m < T_NUM_ISOTOPIC; m ++ ) {
    // INCHI✔️❌:                             EndPointTmp[nNumEndPointTmp].num[T_NUM_NO_ISOTOPIC+m] = atom[endpoint].num_iso_H[NUM_H_ISOTOPES-m-1];
    // INCHI✔️❌:                         }
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         { /* already an endpoint */ /* **now it is wrong:** no mobile atom/charge at this endpoint */
    // INCHI✔️❌:             memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nEquNumber = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         nNumEndPointTmp++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*  add collected tautomeric bonds and endpoints to the input/output data */
    // INCHI✔️❌:     nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    // INCHI✔️❌:     nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *pnNumBondPos = nNumBondPos;
    // INCHI✔️❌:             *pnNumEndPoint = nNumEndPoint;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: #undef PATH_LEN
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Check6MembTautRing
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: Check6MembTautRing
    // INCHI✔️❌: #define PATH_LEN 4
    // INCHI✔️❌: #define TAUT_RINGS_ATTACH_CHAIN 1
    // INCHI✔️❌: #define FIX_BOND23_IN_TAUT 0
    // INCHI✔️❌: #define TAUT_IGNORE_EQL_ENDPOINTS 0
    // INCHI✔️❌: #define REPLACE_ALT_WITH_TAUT 1
    // INCHI✔️❌: #define REPLACE_THE_BOND(X) ( (X) == BOND_SINGLE || (X) == BOND_DOUBLE || (X) == BOND_ALTERN || (X) == BOND_ALT12NS )
    // END INCHI ACTIVE MACRO CONFIGURATION: Check6MembTautRing
    if nStartAtomNeighbor >= 0 || nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0 {
        return Ok(-1);
    }
    if nLenDfsPath != 5 {
        return Ok(-1);
    }

    let middle_pos = nLenDfsPath.wrapping_add(1) / 2;
    let n_at = i32::from(DfsPath.first().ok_or(SourceHeapError::PointerOutOfBounds)?.at_no);
    let nxt_at = i32::from(
        DfsPath
            .get(usize::try_from(middle_pos).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at_no,
    );

    let (o_at, external_bond_pos, external_bond, eif1, eif2, endpoint_valence1, endpoint_valence2, nMobile1, nMobile2) = {
        let atoms = heap.slice(atom.as_const())?;
        let center = atoms
            .get(usize::try_from(nxt_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if center.valence != 3 || center.bCutVertex == 0 {
            return Ok(0);
        }
        let previous = DfsPath
            .get(usize::try_from(middle_pos.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at_no;
        let next = DfsPath
            .get(usize::try_from(middle_pos.wrapping_add(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at_no;
        let valence = usize::try_from(center.valence).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut external = None;
        for i in 0..valence {
            let neighbor = *center.neighbor.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor != previous && neighbor != next {
                external = Some((i, i32::from(neighbor)));
                break;
            }
        }
        let Some((external_bond_pos, o_at)) = external else {
            return Ok(0);
        };
        let external_bond = *center
            .bond_type
            .get(external_bond_pos)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            & !(BOND_MARK_ALL as u8);
        if !matches!(
            external_bond,
            x if x == BOND_SINGLE as u8
                || x == BOND_DOUBLE as u8
                || x == BOND_TAUTOM as u8
                || x == BOND_ALT12NS as u8
                || x == BOND_ALTERN as u8
        ) {
            return Ok(0);
        }
        let mut eif1 = ENDPOINT_INFO::default();
        let endpoint_valence1 = nGetEndpointInfo(atoms, o_at, &mut eif1);
        if endpoint_valence1 == 0 {
            return Ok(0);
        }
        let oxygen = atoms
            .get(usize::try_from(o_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let nMobile1 = i32::from(oxygen.num_H) + i32::from(oxygen.charge == -1);
        if external_bond == BOND_SINGLE as u8 && eif1.cDonor == 0 && oxygen.endpoint == 0 {
            return Ok(0);
        }
        let mut eif2 = ENDPOINT_INFO::default();
        let endpoint_valence2 = nGetEndpointInfo(atoms, n_at, &mut eif2);
        if endpoint_valence2 == 0 {
            return Ok(0);
        }
        let nitrogen = atoms
            .get(usize::try_from(n_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let nMobile2 = i32::from(nitrogen.num_H) + i32::from(nitrogen.charge == -1);
        let mobile_from_o = (oxygen.endpoint != 0 || eif1.cDonor != 0)
            && external_bond != BOND_DOUBLE as u8
            && (nitrogen.endpoint != 0 || eif2.cNeutralBondsValence > nitrogen.valence);
        let mobile_from_n = (nitrogen.endpoint != 0 || eif2.cDonor != 0)
            && (oxygen.endpoint != 0 || eif1.cNeutralBondsValence > oxygen.valence)
            && external_bond != BOND_SINGLE as u8;
        if !mobile_from_o && !mobile_from_n {
            return Ok(0);
        }
        (
            o_at,
            external_bond_pos,
            external_bond,
            eif1,
            eif2,
            endpoint_valence1,
            endpoint_valence2,
            nMobile1,
            nMobile2,
        )
    };

    const PATH_LEN: usize = 4;
    let mut path_bonds = [[0_u8; PATH_LEN + 1]; 2];
    path_bonds[0][0] = external_bond;
    path_bonds[1][0] = external_bond;
    let mut BondPosTmp: [T_BONDPOS; 4 * PATH_LEN] = std::array::from_fn(|_| T_BONDPOS::default());
    let mut nNumBondPosTmp = 0_i32;
    if matches!(
        external_bond,
        x if x == BOND_SINGLE as u8
            || x == BOND_DOUBLE as u8
            || x == BOND_ALTERN as u8
            || x == BOND_ALT12NS as u8
    ) {
        BondPosTmp[0].nAtomNumber = nxt_at as AT_RANK;
        BondPosTmp[0].neighbor_index = external_bond_pos as AT_RANK;
        nNumBondPosTmp = 2;
    }
    for j in 0..middle_pos {
        for side in 0..2_i32 {
            let k = if side != 0 {
                middle_pos.wrapping_add(j)
            } else {
                middle_pos.wrapping_sub(1).wrapping_sub(j)
            };
            let path = DfsPath
                .get(usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let bond_type = path.bond_type;
            path_bonds[usize::try_from(side).unwrap()][usize::try_from(j.wrapping_add(1)).unwrap()] = bond_type;
            if matches!(
                bond_type,
                x if x == BOND_SINGLE as u8
                    || x == BOND_DOUBLE as u8
                    || x == BOND_ALTERN as u8
                    || x == BOND_ALT12NS as u8
            ) {
                let index = usize::try_from(nNumBondPosTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                BondPosTmp[index].nAtomNumber = path.at_no;
                BondPosTmp[index].neighbor_index = path.bond_pos as AT_RANK;
                nNumBondPosTmp = nNumBondPosTmp.wrapping_add(2);
            }
        }
    }
    if are_alt_bonds(&path_bonds[0], middle_pos.wrapping_add(1))? == 0
        || are_alt_bonds(&path_bonds[1], middle_pos.wrapping_add(1))? == 0
    {
        return Ok(0);
    }

    let (endpoint_o, endpoint_n) = {
        let atoms = heap.slice(atom.as_const())?;
        let endpoint_o = atoms[usize::try_from(o_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?].endpoint;
        let endpoint_n = atoms[usize::try_from(n_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?].endpoint;
        (endpoint_o, endpoint_n)
    };
    if endpoint_o != endpoint_n || endpoint_o == 0 {
        let nErr = bExistsAnyAltPath(
            heap,
            pCG,
            pBNS,
            pBD,
            atom,
            num_atoms,
            n_at,
            o_at,
            ALT_PATH_MODE_TAUTOM as i32,
            clock_result,
        )?;
        if nErr <= 0 {
            return Ok(nErr);
        }
    }

    let atoms = heap.slice(atom.as_const())?;
    let mut EndPointTmp: [T_ENDPOINT; 2] = std::array::from_fn(|_| T_ENDPOINT::default());
    let mut nNumEndPointTmp = 0_i32;
    for (j, endpoint) in [o_at, n_at].into_iter().enumerate() {
        let endpoint_atom = atoms
            .get(usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if endpoint_atom.endpoint == 0 {
            let endpoint_valence = if j != 0 { endpoint_valence2 } else { endpoint_valence1 };
            let chem_bonds_valence = if j != 0 {
                i32::from(eif2.cNeutralBondsValence)
            } else {
                i32::from(eif1.cNeutralBondsValence)
            };
            let mobile = if j != 0 { nMobile2 } else { nMobile1 };
            if mobile.wrapping_add(chem_bonds_valence) != endpoint_valence {
                return Ok(0);
            }
            AddAtom2num(&mut EndPointTmp[j].num, atoms, endpoint, 2)?;
            AddAtom2DA(&mut EndPointTmp[j].num_DA, atoms, endpoint, 2)?;
        }
        EndPointTmp[j].nAtomNumber = endpoint as AT_RANK;
        EndPointTmp[j].nGroupNumber = endpoint_atom.endpoint;
        EndPointTmp[j].nEquNumber = 0;
        nNumEndPointTmp = nNumEndPointTmp.wrapping_add(1);
    }

    let nNumBondPos = AddBondsPos(
        atoms,
        &mut BondPosTmp,
        nNumBondPosTmp,
        BondPos,
        nMaxNumBondPos,
        *pnNumBondPos,
    )?;
    let nNumEndPoint = AddEndPoints(&EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, *pnNumEndPoint)?;
    let mut ret = 0_i32;
    if nNumBondPos >= 0 && nNumEndPoint >= 0 && (nNumBondPos > *pnNumBondPos || nNumEndPoint > *pnNumEndPoint) {
        ret = 1;
        *pnNumBondPos = nNumBondPos;
        *pnNumEndPoint = nNumEndPoint;
    }
    Ok(ret)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Check5MembTautRing(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &[DFS_PATH],
    nLenDfsPath: i32,
    _nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1666 Check5MembTautRing
    // INCHI✔️❌: int Check5MembTautRing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                         inp_ATOM *atom,
    // INCHI✔️❌:                         DFS_PATH *DfsPath,
    // INCHI✔️❌:                         int nLenDfsPath,
    // INCHI✔️❌:                         int nStartAtomNeighbor,
    // INCHI✔️❌:                         int nStartAtomNeighbor2,
    // INCHI✔️❌:                         int nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                         T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                         int nMaxNumEndPoint,
    // INCHI✔️❌:                         T_BONDPOS  *BondPos,
    // INCHI✔️❌:                         int nMaxNumBondPos,
    // INCHI✔️❌:                         int *pnNumEndPoint,
    // INCHI✔️❌:                         int *pnNumBondPos,
    // INCHI✔️❌:                         struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                         struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                         int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌: #define PATH_LEN 4
    // INCHI✔️❌:     int i, j, /*m,*/ nMobile, nMobile1, nMobile2, ept, eptn;
    // INCHI✔️❌:     int num_taut_endpoints, nNumBondPos, nNumBondPosTmp, nNumEndPoint, nNumEndPointTmp, ret;
    // INCHI✔️❌:     int endpoint;
    // INCHI✔️❌:     int n1_at = (int) DfsPath[0].at_no;
    // INCHI✔️❌:     int n2_at = (int) DfsPath[1].at_no;
    // INCHI✔️❌:     U_CHAR path_bonds[PATH_LEN + 1], bond_type;
    // INCHI✔️❌:     T_ENDPOINT EndPointTmp[2];
    // INCHI✔️❌:     T_BONDPOS  BondPosTmp[2 * PATH_LEN];
    // INCHI✔️❌:     ENDPOINT_INFO eif1, eif2;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  the two root atoms (atom[n1_at] and atom[n2_at]) cannot belong */
    // INCHI✔️❌:     /*  to one and the same tautomeric group: it has been verified in MarkTautomerGroups() */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  check hydrogens/endpoints */
    // INCHI✔️❌:     if (nLenDfsPath != 4)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /*  program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /*  program error: wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumBondPos = *pnNumBondPos;
    // INCHI✔️❌:     nNumEndPoint = *pnNumEndPoint;
    // INCHI✔️❌:     nNumEndPointTmp = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:     for (ept = 0; ept < 2; ept++) /* djb-rwth: initialisation needed for num array */
    // INCHI✔️❌:         for (eptn = 0; eptn < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; eptn++)
    // INCHI✔️❌:             EndPointTmp[ept].num[eptn] = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!nGetEndpointInfo( atom, n1_at, &eif1 ) ||
    // INCHI✔️❌:         !nGetEndpointInfo( atom, n2_at, &eif2 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nMobile1 = atom[n1_at].num_H + ( atom[n1_at].charge == -1 );
    // INCHI✔️❌:     nMobile2 = atom[n2_at].num_H + ( atom[n2_at].charge == -1 );
    // INCHI✔️❌:     nMobile = nMobile1 + nMobile2;
    // INCHI✔️❌:     num_taut_endpoints = ( 0 != atom[n1_at].endpoint ) + ( 0 != atom[n2_at].endpoint ); /*  if both N atoms already are endpoints */
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( !(nMobile == 1 || num_taut_endpoints == 2) && !(nMobile>1 && num_taut_endpoints >= 1) ) {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (num_taut_endpoints == 0 && nMobile != 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* finally check whether the bonds allow moving the hydrogens */
    // INCHI✔️❌:     if (( atom[n1_at].endpoint != atom[n2_at].endpoint || !atom[n1_at].endpoint ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int nErr;
    // INCHI✔️❌:         nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, n1_at, n2_at, ALT_PATH_MODE_TAUTOM );
    // INCHI✔️❌:         if (nErr <= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return nErr;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  save endpoints */
    // INCHI✔️❌:     for (j = 0; j < 2; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         endpoint = j ? n1_at : n2_at;
    // INCHI✔️❌:         if (!atom[endpoint].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* not a known endpoint */
    // INCHI✔️❌:             /*
    // INCHI✔️❌:                         nMobile  = (atom[endpoint].charge == -1) + atom[endpoint].num_H;
    // INCHI✔️❌:                     } else {
    // INCHI✔️❌:                         nMobile  = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ( nMobile ) {
    // INCHI✔️❌:             */
    // INCHI✔️❌:             AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
    // INCHI✔️❌:             AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             EndPointTmp[nNumEndPointTmp].num[1] = (atom[endpoint].charge == -1);
    // INCHI✔️❌:             EndPointTmp[nNumEndPointTmp].num[0] = nMobile;
    // INCHI✔️❌:             for ( m = 0; m < T_NUM_ISOTOPIC; m ++ ) {
    // INCHI✔️❌:                 EndPointTmp[nNumEndPointTmp].num[T_NUM_NO_ISOTOPIC+m] = atom[endpoint].num_iso_H[NUM_H_ISOTOPES-m-1];
    // INCHI✔️❌:             }
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nEquNumber = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         nNumEndPointTmp++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  extract bonds */
    // INCHI✔️❌:     nNumBondPosTmp = 0;
    // INCHI✔️❌:     for (i = 1; i <= nLenDfsPath; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bond_type = DfsPath[i].bond_type;
    // INCHI✔️❌:         path_bonds[i - 1] = bond_type;
    // INCHI✔️❌:         if (REPLACE_THE_BOND( bond_type ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[i].at_no;
    // INCHI✔️❌:             BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[i].bond_pos;
    // INCHI✔️❌:             nNumBondPosTmp += 2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* path_bonds is from at_n2 to at_n1 */
    // INCHI✔️❌:     if (!( i = are_alt_bonds( path_bonds, nLenDfsPath ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* i is a bond type of the last bond to at_n1, the first bond from at_n2 is 2-i if i=1 or 2 */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* single bond at n1_at: it should have a mobile atom, n2_at should not */
    // INCHI✔️❌:     if ((i == BOND_SINGLE && ( (!atom[n1_at].endpoint && !eif1.cDonor) || (!atom[n2_at].endpoint && !eif2.cAcceptor) )) ||
    // INCHI✔️❌:         /* double bond at n1_at: it should not have a mobile atom, n2_at should */
    // INCHI✔️❌:         (i == BOND_DOUBLE && ( (!atom[n1_at].endpoint && !eif1.cAcceptor) || (!atom[n2_at].endpoint && !eif2.cDonor) ))) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* bond pattern does not fit */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    // INCHI✔️❌:     nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *pnNumBondPos = nNumBondPos;
    // INCHI✔️❌:             *pnNumEndPoint = nNumEndPoint;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌:
    // INCHI✔️❌: #undef PATH_LEN
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Check5MembTautRing
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: Check5MembTautRing
    // INCHI✔️❌: #define PATH_LEN 4
    // INCHI✔️❌: #define REPLACE_ALT_WITH_TAUT 1
    // INCHI✔️❌: #define REPLACE_THE_BOND(X) ( (X) == BOND_SINGLE || (X) == BOND_DOUBLE || (X) == BOND_ALTERN || (X) == BOND_ALT12NS )
    // END INCHI ACTIVE MACRO CONFIGURATION: Check5MembTautRing

    const PATH_LEN: usize = 4;
    if nLenDfsPath != PATH_LEN as i32 {
        return Ok(0);
    }
    if nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0 {
        return Ok(0);
    }

    let path0 = DfsPath.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
    let path1 = DfsPath.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let n1_at = i32::from(path0.at_no);
    let n2_at = i32::from(path1.at_no);
    let (eif1, eif2, endpoint1, endpoint2) = {
        let atoms = heap.slice(atom.as_const())?;
        let mut eif1 = ENDPOINT_INFO::default();
        let mut eif2 = ENDPOINT_INFO::default();
        if nGetEndpointInfo(atoms, n1_at, &mut eif1) == 0 || nGetEndpointInfo(atoms, n2_at, &mut eif2) == 0 {
            return Ok(0);
        }
        let atom1 = atoms
            .get(usize::try_from(n1_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom2 = atoms
            .get(usize::try_from(n2_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let nMobile1 = i32::from(atom1.num_H) + i32::from(atom1.charge == -1);
        let nMobile2 = i32::from(atom2.num_H) + i32::from(atom2.charge == -1);
        let nMobile = nMobile1.wrapping_add(nMobile2);
        let num_taut_endpoints = i32::from(atom1.endpoint != 0) + i32::from(atom2.endpoint != 0);
        if num_taut_endpoints == 0 && nMobile != 1 {
            return Ok(0);
        }
        (eif1, eif2, atom1.endpoint, atom2.endpoint)
    };

    if endpoint1 != endpoint2 || endpoint1 == 0 {
        let nErr = bExistsAnyAltPath(
            heap,
            pCG,
            pBNS,
            pBD,
            atom,
            num_atoms,
            n1_at,
            n2_at,
            ALT_PATH_MODE_TAUTOM as i32,
            clock_result,
        )?;
        if nErr <= 0 {
            return Ok(nErr);
        }
    }

    let atoms = heap.slice(atom.as_const())?;
    let mut nNumBondPos = *pnNumBondPos;
    let mut nNumEndPoint = *pnNumEndPoint;
    let mut nNumEndPointTmp = 0_i32;
    let mut EndPointTmp: [T_ENDPOINT; 2] = std::array::from_fn(|_| T_ENDPOINT::default());
    for endpoint in [n2_at, n1_at] {
        let temporary_index = usize::try_from(nNumEndPointTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let endpoint_atom = atoms
            .get(usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if endpoint_atom.endpoint == 0 {
            AddAtom2num(&mut EndPointTmp[temporary_index].num, atoms, endpoint, 2)?;
            AddAtom2DA(&mut EndPointTmp[temporary_index].num_DA, atoms, endpoint, 2)?;
        } else {
            EndPointTmp[temporary_index] = T_ENDPOINT::default();
        }
        EndPointTmp[temporary_index].nAtomNumber = endpoint as AT_RANK;
        EndPointTmp[temporary_index].nGroupNumber = endpoint_atom.endpoint;
        EndPointTmp[temporary_index].nEquNumber = 0;
        nNumEndPointTmp = nNumEndPointTmp.wrapping_add(1);
    }

    let mut path_bonds = [0_u8; PATH_LEN + 1];
    let mut BondPosTmp: [T_BONDPOS; 2 * PATH_LEN] = std::array::from_fn(|_| T_BONDPOS::default());
    let mut nNumBondPosTmp = 0_i32;
    let mut i = 1_i32;
    while i <= nLenDfsPath {
        let path = DfsPath
            .get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bond_type = path.bond_type;
        path_bonds[usize::try_from(i.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = bond_type;
        if matches!(
            bond_type,
            x if x == BOND_SINGLE as u8
                || x == BOND_DOUBLE as u8
                || x == BOND_ALTERN as u8
                || x == BOND_ALT12NS as u8
        ) {
            let temporary_index = usize::try_from(nNumBondPosTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            BondPosTmp[temporary_index].nAtomNumber = path.at_no;
            BondPosTmp[temporary_index].neighbor_index = path.bond_pos as AT_RANK;
            nNumBondPosTmp = nNumBondPosTmp.wrapping_add(2);
        }
        i = i.wrapping_add(1);
    }
    let terminal_bond = are_alt_bonds(&path_bonds, nLenDfsPath)?;
    if terminal_bond == 0 {
        return Ok(0);
    }

    let atom1 = atoms
        .get(usize::try_from(n1_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atom2 = atoms
        .get(usize::try_from(n2_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if (terminal_bond == BOND_SINGLE as i32
        && ((atom1.endpoint == 0 && eif1.cDonor == 0) || (atom2.endpoint == 0 && eif2.cAcceptor == 0)))
        || (terminal_bond == BOND_DOUBLE as i32
            && ((atom1.endpoint == 0 && eif1.cAcceptor == 0) || (atom2.endpoint == 0 && eif2.cDonor == 0)))
    {
        return Ok(0);
    }

    nNumBondPos = AddBondsPos(
        atoms,
        &mut BondPosTmp,
        nNumBondPosTmp,
        BondPos,
        nMaxNumBondPos,
        nNumBondPos,
    )?;
    nNumEndPoint = AddEndPoints(&EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint)?;

    let mut ret = 0_i32;
    if nNumBondPos >= 0 && nNumEndPoint >= 0 && (nNumBondPos > *pnNumBondPos || nNumEndPoint > *pnNumEndPoint) {
        ret = 1;
        *pnNumBondPos = nNumBondPos;
        *pnNumEndPoint = nNumEndPoint;
    }
    Ok(ret)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn nGet14TautIn7MembAltRing(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    nStartAtom: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighborEndpoint: i32,
    nStartAtomNeighborNeighborEndpoint: i32,
    nDfsPathPos: SourceMutPointer<AT_RANK>,
    DfsPath: &mut [DFS_PATH],
    nMaxLenDfsPath: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:240 nGet14TautIn7MembAltRing
    // INCHI✔️❌: int nGet14TautIn7MembAltRing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                               inp_ATOM *atom,
    // INCHI✔️❌:                               int nStartAtom,
    // INCHI✔️❌:                               int nStartAtomNeighbor,
    // INCHI✔️❌:                               int nStartAtomNeighborEndpoint,
    // INCHI✔️❌:                               int nStartAtomNeighborNeighborEndpoint,
    // INCHI✔️❌:                               AT_RANK  *nDfsPathPos,
    // INCHI✔️❌:                               DFS_PATH *DfsPath,
    // INCHI✔️❌:                               int nMaxLenDfsPath,
    // INCHI✔️❌:                               T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                               int nMaxNumEndPoint,
    // INCHI✔️❌:                               T_BONDPOS  *BondPos,
    // INCHI✔️❌:                               int nMaxNumBondPos,
    // INCHI✔️❌:                               int *pnNumEndPoint,
    // INCHI✔️❌:                               int *pnNumBondPos,
    // INCHI✔️❌:                               struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                               struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                               int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet;
    // INCHI✔️❌:
    // INCHI✔️❌:     *pnNumEndPoint = 0;
    // INCHI✔️❌:     *pnNumBondPos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nMaxLenDfsPath <= 7)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  path is too short */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRet = DFS_FindTautInARing( pCG, atom, nStartAtom, nStartAtomNeighbor,
    // INCHI✔️❌:                                 nStartAtomNeighborEndpoint,
    // INCHI✔️❌:                                 nStartAtomNeighborNeighborEndpoint, 7,
    // INCHI✔️❌:                                 nDfsPathPos, DfsPath, Check7MembTautRing,
    // INCHI✔️❌:                                 bIsCenterPointStrict, EndPoint, nMaxNumEndPoint,
    // INCHI✔️❌:                                 BondPos, nMaxNumBondPos, pnNumEndPoint,
    // INCHI✔️❌:                                 pnNumBondPos, pBNS, pBD, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: nGet14TautIn7MembAltRing
    // SourceHeap preserves allocation identity but retains the DFS/BNS lookup performance gap.

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;
    if nMaxLenDfsPath <= 7 {
        return Ok(-1);
    }

    let mut check_ring = |args: DfsRingCallbackArgs<'_>| {
        let pBNS = args.pBNS.ok_or(SourceHeapError::NullPointer)?;
        let pBD = args.pBD.ok_or(SourceHeapError::NullPointer)?;
        Check7MembTautRing(
            args.heap,
            args.pCG,
            args.atom,
            args.DfsPath,
            args.nLenDfsPath,
            args.nStartAtomNeighbor,
            args.nStartAtomNeighbor2,
            args.nStartAtomNeighborNeighbor,
            args.EndPoint,
            args.nMaxNumEndPoint,
            args.BondPos,
            args.nMaxNumBondPos,
            args.pnNumEndPoint,
            args.pnNumBondPos,
            pBNS,
            pBD,
            args.num_atoms,
            clock_result,
        )
    };
    let mut check_center =
        |heap: &SourceHeap, atom: SourceMutPointer<inp_ATOM>, iat: i32| bIsCenterPointStrict(heap, atom, iat);

    DFS_FindTautInARing(
        heap,
        pCG,
        atom,
        nStartAtom,
        nStartAtomNeighbor,
        nStartAtomNeighborEndpoint,
        nStartAtomNeighborNeighborEndpoint,
        7,
        nDfsPathPos,
        DfsPath,
        &mut check_ring,
        &mut check_center,
        EndPoint,
        nMaxNumEndPoint,
        BondPos,
        nMaxNumBondPos,
        pnNumEndPoint,
        pnNumBondPos,
        Some(pBNS),
        Some(pBD),
        num_atoms,
    )
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn nGet14TautIn5MembAltRing(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    nStartAtom: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighborEndpoint: i32,
    nStartAtomNeighborNeighborEndpoint: i32,
    nDfsPathPos: SourceMutPointer<AT_RANK>,
    DfsPath: &mut [DFS_PATH],
    nMaxLenDfsPath: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:283 nGet14TautIn5MembAltRing
    // INCHI✔️❌: int nGet14TautIn5MembAltRing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                               inp_ATOM *atom,
    // INCHI✔️❌:                               int nStartAtom,
    // INCHI✔️❌:                               int nStartAtomNeighbor,
    // INCHI✔️❌:                               int nStartAtomNeighborEndpoint,
    // INCHI✔️❌:                               int nStartAtomNeighborNeighborEndpoint,
    // INCHI✔️❌:                               AT_RANK  *nDfsPathPos,
    // INCHI✔️❌:                               DFS_PATH *DfsPath,
    // INCHI✔️❌:                               int nMaxLenDfsPath,
    // INCHI✔️❌:                               T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                               int nMaxNumEndPoint,
    // INCHI✔️❌:                               T_BONDPOS  *BondPos,
    // INCHI✔️❌:                               int nMaxNumBondPos,
    // INCHI✔️❌:                               int *pnNumEndPoint,
    // INCHI✔️❌:                               int *pnNumBondPos,
    // INCHI✔️❌:                               struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                               struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                               int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet;
    // INCHI✔️❌:
    // INCHI✔️❌:     *pnNumEndPoint = 0;
    // INCHI✔️❌:     *pnNumBondPos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nMaxLenDfsPath <= 5)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  path is too short */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRet = DFS_FindTautInARing( pCG, atom, nStartAtom, nStartAtomNeighbor,
    // INCHI✔️❌:         nStartAtomNeighborEndpoint, nStartAtomNeighborNeighborEndpoint, 5,
    // INCHI✔️❌:         nDfsPathPos, DfsPath, Check7MembTautRing, bIsCenterPointStrict,
    // INCHI✔️❌:         EndPoint, nMaxNumEndPoint, BondPos, nMaxNumBondPos,
    // INCHI✔️❌:         pnNumEndPoint, pnNumBondPos, pBNS, pBD, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: nGet14TautIn5MembAltRing
    // SourceHeap preserves allocation identity but retains the DFS/BNS lookup performance gap.

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;
    if nMaxLenDfsPath <= 5 {
        return Ok(-1);
    }

    let mut check_ring = |args: DfsRingCallbackArgs<'_>| {
        let pBNS = args.pBNS.ok_or(SourceHeapError::NullPointer)?;
        let pBD = args.pBD.ok_or(SourceHeapError::NullPointer)?;
        Check7MembTautRing(
            args.heap,
            args.pCG,
            args.atom,
            args.DfsPath,
            args.nLenDfsPath,
            args.nStartAtomNeighbor,
            args.nStartAtomNeighbor2,
            args.nStartAtomNeighborNeighbor,
            args.EndPoint,
            args.nMaxNumEndPoint,
            args.BondPos,
            args.nMaxNumBondPos,
            args.pnNumEndPoint,
            args.pnNumBondPos,
            pBNS,
            pBD,
            args.num_atoms,
            clock_result,
        )
    };
    let mut check_center =
        |heap: &SourceHeap, atom: SourceMutPointer<inp_ATOM>, iat: i32| bIsCenterPointStrict(heap, atom, iat);

    DFS_FindTautInARing(
        heap,
        pCG,
        atom,
        nStartAtom,
        nStartAtomNeighbor,
        nStartAtomNeighborEndpoint,
        nStartAtomNeighborNeighborEndpoint,
        5,
        nDfsPathPos,
        DfsPath,
        &mut check_ring,
        &mut check_center,
        EndPoint,
        nMaxNumEndPoint,
        BondPos,
        nMaxNumBondPos,
        pnNumEndPoint,
        pnNumBondPos,
        Some(pBNS),
        Some(pBD),
        num_atoms,
    )
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn nGet12TautIn5MembAltRing(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    nStartAtom: i32,
    nStartAtomNeighbor: i32,
    nDfsPathPos: SourceMutPointer<AT_RANK>,
    DfsPath: &mut [DFS_PATH],
    nMaxLenDfsPath: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:322 nGet12TautIn5MembAltRing
    // INCHI✔️❌: int nGet12TautIn5MembAltRing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                               inp_ATOM *atom,
    // INCHI✔️❌:                               int nStartAtom,
    // INCHI✔️❌:                               int nStartAtomNeighbor,
    // INCHI✔️❌:                               AT_RANK  *nDfsPathPos,
    // INCHI✔️❌:                               DFS_PATH *DfsPath,
    // INCHI✔️❌:                               int nMaxLenDfsPath,
    // INCHI✔️❌:                               T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                               int nMaxNumEndPoint,
    // INCHI✔️❌:                               T_BONDPOS  *BondPos,
    // INCHI✔️❌:                               int nMaxNumBondPos,
    // INCHI✔️❌:                               int *pnNumEndPoint,
    // INCHI✔️❌:                               int *pnNumBondPos,
    // INCHI✔️❌:                               struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                               struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                               int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet;
    // INCHI✔️❌:
    // INCHI✔️❌:     *pnNumEndPoint = 0;
    // INCHI✔️❌:     *pnNumBondPos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nMaxLenDfsPath <= 5)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  path is too short */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRet = DFS_FindTautInARing( pCG, atom, nStartAtom, nStartAtomNeighbor,
    // INCHI✔️❌:                                 -1, -1, 5,
    // INCHI✔️❌:                                 nDfsPathPos, DfsPath, Check5MembTautRing,
    // INCHI✔️❌:                                 bIsCenterPointStrict, EndPoint, nMaxNumEndPoint,
    // INCHI✔️❌:                                 BondPos, nMaxNumBondPos, pnNumEndPoint,
    // INCHI✔️❌:                                 pnNumBondPos, pBNS, pBD, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: nGet12TautIn5MembAltRing
    // SourceHeap preserves allocation identity but retains the DFS/BNS lookup performance gap.

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;
    if nMaxLenDfsPath <= 5 {
        return Ok(-1);
    }

    let mut check_ring = |args: DfsRingCallbackArgs<'_>| {
        let pBNS = args.pBNS.ok_or(SourceHeapError::NullPointer)?;
        let pBD = args.pBD.ok_or(SourceHeapError::NullPointer)?;
        Check5MembTautRing(
            args.heap,
            args.pCG,
            args.atom,
            args.DfsPath,
            args.nLenDfsPath,
            args.nStartAtomNeighbor,
            args.nStartAtomNeighbor2,
            args.nStartAtomNeighborNeighbor,
            args.EndPoint,
            args.nMaxNumEndPoint,
            args.BondPos,
            args.nMaxNumBondPos,
            args.pnNumEndPoint,
            args.pnNumBondPos,
            pBNS,
            pBD,
            args.num_atoms,
            clock_result,
        )
    };
    let mut check_center =
        |heap: &SourceHeap, atom: SourceMutPointer<inp_ATOM>, iat: i32| bIsCenterPointStrict(heap, atom, iat);

    DFS_FindTautInARing(
        heap,
        pCG,
        atom,
        nStartAtom,
        nStartAtomNeighbor,
        -1,
        -1,
        5,
        nDfsPathPos,
        DfsPath,
        &mut check_ring,
        &mut check_center,
        EndPoint,
        nMaxNumEndPoint,
        BondPos,
        nMaxNumBondPos,
        pnNumEndPoint,
        pnNumBondPos,
        Some(pBNS),
        Some(pBD),
        num_atoms,
    )
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn nGet15TautIn6MembAltRing(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    nStartAtom: i32,
    nDfsPathPos: SourceMutPointer<AT_RANK>,
    DfsPath: &mut [DFS_PATH],
    nMaxLenDfsPath: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:361 nGet15TautIn6MembAltRing
    // INCHI✔️❌: int nGet15TautIn6MembAltRing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                               inp_ATOM *atom,
    // INCHI✔️❌:                               int nStartAtom,
    // INCHI✔️❌:                               AT_RANK  *nDfsPathPos,
    // INCHI✔️❌:                               DFS_PATH *DfsPath,
    // INCHI✔️❌:                               int nMaxLenDfsPath,
    // INCHI✔️❌:                               T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                               int nMaxNumEndPoint,
    // INCHI✔️❌:                               T_BONDPOS  *BondPos,
    // INCHI✔️❌:                               int nMaxNumBondPos,
    // INCHI✔️❌:                               int *pnNumEndPoint,
    // INCHI✔️❌:                               int *pnNumBondPos,
    // INCHI✔️❌:                               struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                               struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                               int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet;
    // INCHI✔️❌:
    // INCHI✔️❌:     *pnNumEndPoint = 0;
    // INCHI✔️❌:     *pnNumBondPos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nMaxLenDfsPath <= 7)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  path is too short */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRet = DFS_FindTautInARing( pCG, atom, nStartAtom,
    // INCHI✔️❌:                                 -1/*nStartAtomNeighbor*/,
    // INCHI✔️❌:                                 -1/*nStartAtomNeighbor2*/,
    // INCHI✔️❌:                                 -1/*nStartAtomNeighborNeighbor*/,
    // INCHI✔️❌:                                 6 /* nCycleLen*/,
    // INCHI✔️❌:                                 nDfsPathPos, DfsPath,
    // INCHI✔️❌:                                 Check6MembTautRing, bIsCenterPointStrict,
    // INCHI✔️❌:                                 EndPoint, nMaxNumEndPoint,
    // INCHI✔️❌:                                 BondPos, nMaxNumBondPos,
    // INCHI✔️❌:                                 pnNumEndPoint, pnNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: nGet15TautIn6MembAltRing
    // SourceHeap preserves allocation identity but retains the DFS/BNS lookup performance gap.

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;
    if nMaxLenDfsPath <= 7 {
        return Ok(-1);
    }
    let mut check_ring = |args: DfsRingCallbackArgs<'_>| {
        let pBNS = args.pBNS.ok_or(SourceHeapError::NullPointer)?;
        let pBD = args.pBD.ok_or(SourceHeapError::NullPointer)?;
        Check6MembTautRing(
            args.heap,
            args.pCG,
            args.atom,
            args.DfsPath,
            args.nLenDfsPath,
            args.nStartAtomNeighbor,
            args.nStartAtomNeighbor2,
            args.nStartAtomNeighborNeighbor,
            args.EndPoint,
            args.nMaxNumEndPoint,
            args.BondPos,
            args.nMaxNumBondPos,
            args.pnNumEndPoint,
            args.pnNumBondPos,
            pBNS,
            pBD,
            args.num_atoms,
            clock_result,
        )
    };
    let mut check_center =
        |heap: &SourceHeap, atom: SourceMutPointer<inp_ATOM>, iat: i32| bIsCenterPointStrict(heap, atom, iat);
    DFS_FindTautInARing(
        heap,
        pCG,
        atom,
        nStartAtom,
        -1,
        -1,
        -1,
        6,
        nDfsPathPos,
        DfsPath,
        &mut check_ring,
        &mut check_center,
        EndPoint,
        nMaxNumEndPoint,
        BondPos,
        nMaxNumBondPos,
        pnNumEndPoint,
        pnNumBondPos,
        Some(pBNS),
        Some(pBD),
        num_atoms,
    )
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn nGet15TautInAltPath(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    nStartAtom: i32,
    nDfsPathPos: SourceMutPointer<AT_RANK>,
    DfsPath: &mut [DFS_PATH],
    nMaxLenDfsPath: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:407 nGet15TautInAltPath
    // INCHI✔️❌: int nGet15TautInAltPath( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                          inp_ATOM *atom,
    // INCHI✔️❌:                          int nStartAtom,
    // INCHI✔️❌:                          AT_RANK  *nDfsPathPos,
    // INCHI✔️❌:                          DFS_PATH *DfsPath,
    // INCHI✔️❌:                          int nMaxLenDfsPath,
    // INCHI✔️❌:                          T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                          int nMaxNumEndPoint,
    // INCHI✔️❌:                          T_BONDPOS  *BondPos,
    // INCHI✔️❌:                          int nMaxNumBondPos,
    // INCHI✔️❌:                          int *pnNumEndPoint,
    // INCHI✔️❌:                          int *pnNumBondPos,
    // INCHI✔️❌:                          struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                          struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                          int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet;
    // INCHI✔️❌:
    // INCHI✔️❌:     *pnNumEndPoint = 0;
    // INCHI✔️❌:     *pnNumBondPos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nMaxLenDfsPath <= 7)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  path is too short */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRet = DFS_FindTautAltPath( pCG, atom, nStartAtom,
    // INCHI✔️❌:                                 -1/*nStartAtomNeighbor*/,
    // INCHI✔️❌:                                 -1/*nStartAtomNeighbor2*/,
    // INCHI✔️❌:                                 -1/*nStartAtomNeighborNeighbor*/,
    // INCHI✔️❌:                                 4 /* nCycleLen*/,
    // INCHI✔️❌:                                 nDfsPathPos, DfsPath,
    // INCHI✔️❌:                                 Check15TautPath, Check15TautPathCenterpoint,
    // INCHI✔️❌:                                 EndPoint, nMaxNumEndPoint,
    // INCHI✔️❌:                                 BondPos, nMaxNumBondPos,
    // INCHI✔️❌:                                 pnNumEndPoint, pnNumBondPos,
    // INCHI✔️❌:                                 pBNS, pBD, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: nGet15TautInAltPath
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nGet15TautInAltPath
    // INCHI✔️❌: #define TAUT_15_NON_RING 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nGet15TautInAltPath
    // SourceHeap preserves allocation identity but retains the DFS/BNS lookup performance gap.

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;
    if nMaxLenDfsPath <= 7 {
        return Ok(-1);
    }

    let mut check_path = |args: DfsPathCallbackArgs<'_>| {
        let pBNS = args.pBNS.ok_or(SourceHeapError::NullPointer)?;
        let pBD = args.pBD.ok_or(SourceHeapError::NullPointer)?;
        Check15TautPath(
            args.heap,
            args.pCG,
            args.atom,
            args.DfsPath,
            args.nLenDfsPath,
            args.jNxtNeigh,
            args.nStartAtomNeighbor,
            args.nStartAtomNeighbor2,
            args.nStartAtomNeighborNeighbor,
            args.EndPoint,
            args.nMaxNumEndPoint,
            args.BondPos,
            args.nMaxNumBondPos,
            args.pnNumEndPoint,
            args.pnNumBondPos,
            pBNS,
            pBD,
            args.num_atoms,
            clock_result,
        )
    };
    let mut check_center = |args: DfsPathCenterCallbackArgs<'_>| {
        Check15TautPathCenterpoint(
            args.heap,
            args.atom,
            args.DfsPath,
            args.nLenDfsPath,
            args.jNxtNeigh,
            args.pBNS,
            args.pBD,
            args.num_atoms,
        )
    };

    DFS_FindTautAltPath(
        heap,
        pCG,
        atom,
        nStartAtom,
        -1,
        -1,
        -1,
        4,
        nDfsPathPos,
        DfsPath,
        &mut check_path,
        &mut check_center,
        EndPoint,
        nMaxNumEndPoint,
        BondPos,
        nMaxNumBondPos,
        pnNumEndPoint,
        pnNumBondPos,
        Some(pBNS),
        Some(pBD),
        num_atoms,
    )
}

#[allow(non_snake_case)]
pub(crate) fn bIsCenterPointStrict(
    heap: &SourceHeap,
    atom: SourceMutPointer<inp_ATOM>,
    iat: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:215 bIsCenterPointStrict
    // INCHI✔️❌: int bIsCenterPointStrict( inp_ATOM *atom, int iat )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (atom[iat].valence == atom[iat].chem_bonds_valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int endpoint_valence = get_endpoint_valence( atom[iat].el_number );
    // INCHI✔️❌:         if (endpoint_valence && ( (endpoint_valence > atom[iat].valence && /* added a check for negative charge or H 3-31-03 */
    // INCHI✔️❌:             ( atom[iat].num_H || atom[iat].charge == -1 )) ||
    // INCHI✔️❌:             (!atom[iat].charge && atom[iat].c_point) )) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1; /*  may appear to be tautomeric or chargable
    // INCHI✔️❌:                           (this increases chem_bonds_valence), should be explored */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (atom[iat].valence + 1 == atom[iat].chem_bonds_valence &&
    // INCHI✔️❌:         is_centerpoint_elem_strict( atom[iat].el_number ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: bIsCenterPointStrict
    // BEGIN INCHI ACTIVE HEADER CONFIGURATION: bIsCenterPointStrict
    // INCHI✔️❌: typedef signed char S_CHAR;
    // END INCHI ACTIVE HEADER CONFIGURATION: bIsCenterPointStrict
    // SourceHeap preserves C pointer identity but adds a BTreeMap allocation lookup.

    let atom = heap
        .slice(atom.as_const())?
        .get(usize::try_from(iat).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.valence == atom.chem_bonds_valence {
        let endpoint_valence = get_endpoint_valence(atom.el_number);
        if endpoint_valence != 0
            && ((endpoint_valence > i32::from(atom.valence) && (atom.num_H != 0 || atom.charge == -1))
                || (atom.charge == 0 && atom.c_point != 0))
        {
            return Ok(1);
        }
        return Ok(0);
    }
    if i32::from(atom.valence) + 1 == i32::from(atom.chem_bonds_valence)
        && is_centerpoint_elem_strict(atom.el_number) != 0
    {
        return Ok(1);
    }

    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Check15TautPathCenterpoint(
    heap: &SourceHeap,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &[DFS_PATH],
    nLenDfsPath: i32,
    jNxtNeigh: i32,
    _pBNS: Option<&mut BalancedNetworkStructure>,
    _pBD: Option<&mut BalancedNetworkData>,
    _num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1409 Check15TautPathCenterpoint
    // INCHI✔️❌: int Check15TautPathCenterpoint( inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int jNxtNeigh,
    // INCHI✔️❌:     struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:     struct BalancedNetworkData *pBD, int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nxt_at = atom[DfsPath[nLenDfsPath].at_no].neighbor[jNxtNeigh];
    // INCHI✔️❌:     /* atom[nxt_at].endpoint below allows for keto-enol -CH< or -CH2- endpoints  */
    // INCHI✔️❌:     return atom[nxt_at].endpoint || bIsCenterPointStrict( atom, nxt_at );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Check15TautPathCenterpoint
    // BEGIN INCHI ACTIVE TYPE CONFIGURATION: Check15TautPathCenterpoint
    // INCHI✔️❌: typedef unsigned short AT_RANK;
    // INCHI✔️❌: typedef signed char S_CHAR;
    // END INCHI ACTIVE TYPE CONFIGURATION: Check15TautPathCenterpoint
    // SourceHeap preserves C pointer identity but adds BTreeMap lookups for atom access.

    let path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current_atom_index = usize::from(
        DfsPath
            .get(path_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at_no,
    );
    let neighbor_index = usize::try_from(jNxtNeigh).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let nxt_at = i32::from(
        *heap
            .slice(atom.as_const())?
            .get(current_atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .neighbor
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let next_atom_index = usize::try_from(nxt_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let endpoint = heap
        .slice(atom.as_const())?
        .get(next_atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .endpoint;
    if endpoint != 0 {
        return Ok(1);
    }
    bIsCenterPointStrict(heap, atom, nxt_at)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Check15TautPath(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &mut [DFS_PATH],
    mut nLenDfsPath: i32,
    jNxtNeigh: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    pBNS: &mut BalancedNetworkStructure,
    pBD: &mut BalancedNetworkData,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1434 Check15TautPath
    // INCHI✔️❌: int Check15TautPath( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                      inp_ATOM *atom, DFS_PATH *DfsPath,
    // INCHI✔️❌:                      int nLenDfsPath,
    // INCHI✔️❌:                      int jNxtNeigh,
    // INCHI✔️❌:                      int nStartAtomNeighbor,
    // INCHI✔️❌:                      int nStartAtomNeighbor2,
    // INCHI✔️❌:                      int nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                      T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                      int nMaxNumEndPoint,
    // INCHI✔️❌:                      T_BONDPOS  *BondPos,
    // INCHI✔️❌:                      int nMaxNumBondPos,
    // INCHI✔️❌:                      int *pnNumEndPoint,
    // INCHI✔️❌:                      int *pnNumBondPos,
    // INCHI✔️❌:                      struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                      struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                      int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌: #define PATH_LEN 4
    // INCHI✔️❌:     int i, j, k, /*m,*/ nNumBondPos, nNumEndPoint, cur_at, prv_at, at1, at2 /*, at3, step_at*/;
    // INCHI✔️❌:     int nNumEndPointTmp, nNumBondPosTmp, ret;
    // INCHI✔️❌:     /* int num_taut_endpoints, num_H; */
    // INCHI✔️❌:     int nMobile, endpoint, endpoint_valence, chem_bonds_valence;
    // INCHI✔️❌:     int nMobile1, endpoint_valence1;  /* start atom, at1 */
    // INCHI✔️❌:     int nMobile2, endpoint_valence2;  /* end atom,   at2 */
    // INCHI✔️❌:     /*int nMobile3, endpoint_valence3=-1;*/  /* middle atom, at3 */
    // INCHI✔️❌:     /*int nxt_at;*/
    // INCHI✔️❌:     int alt_bonds[2];
    // INCHI✔️❌:     U_CHAR /*path_bonds[2][PATH_LEN+1],*/ bond_type;
    // INCHI✔️❌:     T_ENDPOINT EndPointTmp[2];
    // INCHI✔️❌:     T_BONDPOS  BondPosTmp[4 * PATH_LEN];
    // INCHI✔️❌:     ENDPOINT_INFO eif1, eif2/*, eif3*/;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nStartAtomNeighbor >= 0 || nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nLenDfsPath != 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /*  wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumBondPos = *pnNumBondPos;
    // INCHI✔️❌:     nNumEndPoint = *pnNumEndPoint;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     nNumEndPointTmp = 0;
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*-------add the last atom, nLenDfsPath=4 --*/
    // INCHI✔️❌:     j = jNxtNeigh;
    // INCHI✔️❌:     prv_at = DfsPath[nLenDfsPath].at_no;
    // INCHI✔️❌:     cur_at = atom[prv_at].neighbor[j];
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_type = ( atom[prv_at].bond_type[j] & ~BOND_MARK_ALL );
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:     DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, prv_at, j, DfsPath[nLenDfsPath].bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_pos = j; /* fix index of the bond to the next atom */
    // INCHI✔️❌:
    // INCHI✔️❌:     nLenDfsPath++;
    // INCHI✔️❌:
    // INCHI✔️❌:     DfsPath[nLenDfsPath].at_no = cur_at;
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_type = 0;
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_pos = -1;
    // INCHI✔️❌:     /*nDfsPathPos[cur_at]                = nLenDfsPath+1;*/ /* mark with distance + 1 */
    // INCHI✔️❌: /*------------------------------------------*/
    // INCHI✔️❌:     at1 = (int) DfsPath[0].at_no;
    // INCHI✔️❌:     at2 = (int) DfsPath[nLenDfsPath].at_no;
    // INCHI✔️❌:     /*at3 = (int)DfsPath[2].at_no;*/
    // INCHI✔️❌:     if (atom[at1].endpoint && atom[at1].endpoint == atom[at2].endpoint)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* start & end already belong to the same taut group */
    // INCHI✔️❌:         goto exit_function;  /* nothing to do */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check bond types along alt path */
    // INCHI✔️❌:     alt_bonds[0] = alt_bonds[1] = 0;
    // INCHI✔️❌:     for (i = 0; i < nLenDfsPath; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         alt_bonds[i % 2] |= IS_ALT_OR_DBLBOND( DfsPath[i].bond_type );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (( alt_bonds[0] & alt_bonds[1] & ( BOND_SINGLE | BOND_DOUBLE ) ) ||
    // INCHI✔️❌:         ( alt_bonds[0] & BOND_WRONG ) || ( alt_bonds[1] & BOND_WRONG ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function; /* incompatible with alt path or wrong bonds */\
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* check possibly tautomeric endpoints at the ends */
    // INCHI✔️❌:     endpoint_valence1 = nGetEndpointInfo( atom, at1, &eif1 );
    // INCHI✔️❌:     endpoint_valence2 = nGetEndpointInfo( atom, at2, &eif2 );
    // INCHI❌❌: #ifdef NEVER   /* do not use C-endpoint of keto-enol tautomer to find 1,5 the taut path */
    // INCHI❌❌:     if (!endpoint_valence1 && !atom[at1].endpoint ||
    // INCHI❌❌:         !endpoint_valence2 && !atom[at2].endpoint)
    // INCHI❌❌:         goto exit_function; /* at least one of the end atoms cannot be an endpoint */
    // INCHI❌❌: #endif
    // INCHI✔️❌:     if (!endpoint_valence1 || !endpoint_valence2)
    // INCHI✔️❌:         goto exit_function;  /* require both endpoints be heteroatoms */
    // INCHI✔️❌:     /*  check hydrogens/endpoints */
    // INCHI✔️❌:     nMobile1 = atom[at1].num_H + ( atom[at1].charge == -1 );
    // INCHI✔️❌:     if (!atom[at1].endpoint)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( alt_bonds[0] & BOND_SINGLE ) && !eif1.cDonor)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (( alt_bonds[0] & BOND_DOUBLE ) && !eif1.cAcceptor)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     nMobile2 = atom[at2].num_H + ( atom[at2].charge == -1 );
    // INCHI✔️❌:     if (!atom[at2].endpoint)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( alt_bonds[1] & BOND_SINGLE ) && !eif2.cDonor)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (( alt_bonds[1] & BOND_DOUBLE ) && !eif2.cAcceptor)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nMobile = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  can mobile group move from at1=o_at to at2=n_at? */
    // INCHI✔️❌:     nMobile += ( atom[at1].endpoint || eif1.cDonor ) &&  /*  from o_at */
    // INCHI✔️❌:         !( alt_bonds[0] & BOND_DOUBLE ) &&
    // INCHI✔️❌:         ( atom[at2].endpoint ||                   /*  to n_at */
    // INCHI✔️❌:             eif2.cNeutralBondsValence > atom[at2].valence );
    // INCHI✔️❌:     /*  can mobile group move from at2=n_at to at1=o_at? */
    // INCHI✔️❌:     nMobile += ( atom[at2].endpoint || eif2.cDonor ) && /*  from n_at */
    // INCHI✔️❌:         !( alt_bonds[1] & BOND_DOUBLE ) &&
    // INCHI✔️❌:         ( atom[at1].endpoint ||          /*  to o_at */
    // INCHI✔️❌:             eif1.cNeutralBondsValence > atom[at1].valence );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!nMobile)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check whether the bonds allow moving the hydrogens between at1 and at2 */
    // INCHI✔️❌:     if (( atom[at1].endpoint != atom[at2].endpoint || !atom[at1].endpoint ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int nErr;
    // INCHI✔️❌:         nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, at1, at2, ALT_PATH_MODE_TAUTOM );
    // INCHI✔️❌:         if (nErr <= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = nErr;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* save tautomeric bonds */
    // INCHI✔️❌:     nNumBondPosTmp = 0;
    // INCHI✔️❌:     for (k = 0; k < nLenDfsPath; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bond_type = DfsPath[k].bond_type;
    // INCHI✔️❌:         if (REPLACE_THE_BOND( bond_type ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[k].at_no;     /*  accumulate bonds to be */
    // INCHI✔️❌:             BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[k].bond_pos;  /*  marked as tautomeric */
    // INCHI✔️❌:             nNumBondPosTmp += 2; /*  leave room for the same bond in opposite direction */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*  save endpoints */
    // INCHI✔️❌:     for (j = 0; j < 2; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         endpoint = j ? at2 : at1;
    // INCHI✔️❌:         if (!atom[endpoint].endpoint)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* not a known endpoint */
    // INCHI✔️❌:             endpoint_valence = j ? endpoint_valence2 : endpoint_valence1;
    // INCHI✔️❌:             chem_bonds_valence = j ? eif2.cNeutralBondsValence : eif1.cNeutralBondsValence;
    // INCHI✔️❌:             /* endpoint_valence = get_endpoint_valence( atom[endpoint].el_number ); */
    // INCHI✔️❌:             nMobile = j ? nMobile2 : nMobile1;
    // INCHI✔️❌:             /* nMobile  = (atom[endpoint].charge == -1) + atom[endpoint].num_H; */
    // INCHI✔️❌:             /* if ( nMobile + atom[endpoint].chem_bonds_valence != endpoint_valence ) -- fixed 02-06-2003*/
    // INCHI✔️❌:             if (nMobile + chem_bonds_valence != endpoint_valence)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_function; /*  abnormal endpoint valence; ignore. */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
    // INCHI✔️❌:             AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* already an endpoint */ /* **now it is wrong:** no mobile atom/charge at this endpoint */
    // INCHI✔️❌:             memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
    // INCHI✔️❌:         EndPointTmp[nNumEndPointTmp].nEquNumber = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         nNumEndPointTmp++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*  add collected tautomeric bonds and endpoints to the input/output data */
    // INCHI✔️❌:     nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    // INCHI✔️❌:     nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *pnNumBondPos = nNumBondPos;
    // INCHI✔️❌:             *pnNumEndPoint = nNumEndPoint;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     /*nDfsPathPos[DfsPath[nLenDfsPath].at_no] = 0;*/
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌:
    // INCHI✔️❌: #undef PATH_LEN
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Check15TautPath
    // BEGIN INCHI ACTIVE MACRO AND TYPE CONFIGURATION: Check15TautPath
    // INCHI✔️❌: #define PATH_LEN 4
    // INCHI✔️❌: #define TAUT_15_NON_RING 1
    // INCHI✔️❌: #define FIX_BOND23_IN_TAUT 0
    // INCHI✔️❌: #define REPLACE_ALT_WITH_TAUT 1
    // INCHI✔️❌: #define BOND_WRONG 64
    // INCHI✔️❌: #define IS_ALT_OR_DBLBOND(X) (((X) == BOND_SINGLE || (X) == BOND_DOUBLE)? (X) : ((X) == BOND_ALTERN || (X) == BOND_TAUTOM || (X) == BOND_ALT12NS)? BOND_ALTERN : BOND_WRONG)
    // INCHI✔️❌: #define REPLACE_THE_BOND(X) ( (X) == BOND_SINGLE || (X) == BOND_DOUBLE || (X) == BOND_ALTERN || (X) == BOND_ALT12NS )
    // INCHI✔️❌: #define BOND_MARK_ALL 0xf0
    // INCHI✔️❌: typedef unsigned short AT_RANK;
    // INCHI✔️❌: typedef signed char S_CHAR;
    // END INCHI ACTIVE MACRO AND TYPE CONFIGURATION: Check15TautPath
    // SourceHeap preserves C pointer identity but adds BTreeMap lookups for atom and BNS access.

    const PATH_LEN: usize = 4;
    const BOND_WRONG: i32 = 64;
    if nStartAtomNeighbor >= 0 || nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0 {
        return Ok(-1);
    }
    if nLenDfsPath != 3 {
        return Ok(-1);
    }

    let mut nNumBondPos = *pnNumBondPos;
    let mut nNumEndPoint = *pnNumEndPoint;
    let mut ret = 0_i32;

    let path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let prv_at = i32::from(
        DfsPath
            .get(path_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at_no,
    );
    let neighbor_index = usize::try_from(jNxtNeigh).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let (cur_at, bond_type) = {
        let previous = heap
            .slice(atom.as_const())?
            .get(usize::try_from(prv_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        (
            i32::from(
                *previous
                    .neighbor
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            ),
            *previous
                .bond_type
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                & !(BOND_MARK_ALL as u8),
        )
    };
    let path = DfsPath.get_mut(path_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    path.bond_type = bond_type;
    path.bond_pos = jNxtNeigh as i8;
    nLenDfsPath = nLenDfsPath.wrapping_add(1);
    let end_path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let end_path = DfsPath
        .get_mut(end_path_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    end_path.at_no = cur_at as AT_RANK;
    end_path.bond_type = 0;
    end_path.bond_pos = -1;

    let at1 = i32::from(DfsPath.first().ok_or(SourceHeapError::PointerOutOfBounds)?.at_no);
    let at2 = i32::from(
        DfsPath
            .get(end_path_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .at_no,
    );
    let at1_index = usize::try_from(at1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let at2_index = usize::try_from(at2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    {
        let atoms = heap.slice(atom.as_const())?;
        let endpoint1 = atoms
            .get(at1_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint;
        let endpoint2 = atoms
            .get(at2_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint;
        if endpoint1 != 0 && endpoint1 == endpoint2 {
            return Ok(ret);
        }
    }

    let classify_bond = |bond: u8| -> i32 {
        if bond == BOND_SINGLE as u8 || bond == BOND_DOUBLE as u8 {
            i32::from(bond)
        } else if bond == BOND_ALTERN as u8 || bond == BOND_TAUTOM as u8 || bond == BOND_ALT12NS as u8 {
            BOND_ALTERN as i32
        } else {
            BOND_WRONG
        }
    };
    let mut alt_bonds = [0_i32; 2];
    for i in 0..end_path_index {
        let bond = DfsPath.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?.bond_type;
        alt_bonds[i % 2] |= classify_bond(bond);
    }
    if (alt_bonds[0] & alt_bonds[1] & (BOND_SINGLE as i32 | BOND_DOUBLE as i32)) != 0
        || (alt_bonds[0] & BOND_WRONG) != 0
        || (alt_bonds[1] & BOND_WRONG) != 0
    {
        return Ok(ret);
    }

    let mut eif1 = ENDPOINT_INFO::default();
    let mut eif2 = ENDPOINT_INFO::default();
    let (endpoint_valence1, endpoint_valence2) = {
        let atoms = heap.slice(atom.as_const())?;
        (
            nGetEndpointInfo(atoms, at1, &mut eif1),
            nGetEndpointInfo(atoms, at2, &mut eif2),
        )
    };
    if endpoint_valence1 == 0 || endpoint_valence2 == 0 {
        return Ok(ret);
    }

    let (nMobile1, nMobile2, endpoint1, endpoint2, valence1, valence2) = {
        let atoms = heap.slice(atom.as_const())?;
        let atom1 = atoms.get(at1_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom2 = atoms.get(at2_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        (
            i32::from(atom1.num_H) + i32::from(atom1.charge == -1),
            i32::from(atom2.num_H) + i32::from(atom2.charge == -1),
            atom1.endpoint,
            atom2.endpoint,
            i32::from(atom1.valence),
            i32::from(atom2.valence),
        )
    };
    if endpoint1 == 0 {
        if (alt_bonds[0] & BOND_SINGLE as i32) != 0 && eif1.cDonor == 0 {
            return Ok(ret);
        }
        if (alt_bonds[0] & BOND_DOUBLE as i32) != 0 && eif1.cAcceptor == 0 {
            return Ok(ret);
        }
    }
    if endpoint2 == 0 {
        if (alt_bonds[1] & BOND_SINGLE as i32) != 0 && eif2.cDonor == 0 {
            return Ok(ret);
        }
        if (alt_bonds[1] & BOND_DOUBLE as i32) != 0 && eif2.cAcceptor == 0 {
            return Ok(ret);
        }
    }

    let mut nMobile = 0_i32;
    nMobile += i32::from(
        (endpoint1 != 0 || eif1.cDonor != 0)
            && (alt_bonds[0] & BOND_DOUBLE as i32) == 0
            && (endpoint2 != 0 || i32::from(eif2.cNeutralBondsValence) > valence2),
    );
    nMobile += i32::from(
        (endpoint2 != 0 || eif2.cDonor != 0)
            && (alt_bonds[1] & BOND_DOUBLE as i32) == 0
            && (endpoint1 != 0 || i32::from(eif1.cNeutralBondsValence) > valence1),
    );
    if nMobile == 0 {
        return Ok(ret);
    }

    if endpoint1 != endpoint2 || endpoint1 == 0 {
        let nErr = bExistsAnyAltPath(
            heap,
            pCG,
            pBNS,
            pBD,
            atom,
            num_atoms,
            at1,
            at2,
            ALT_PATH_MODE_TAUTOM as i32,
            clock_result,
        )?;
        if nErr <= 0 {
            return Ok(nErr);
        }
    }

    let mut BondPosTmp: [T_BONDPOS; 4 * PATH_LEN] = std::array::from_fn(|_| T_BONDPOS::default());
    let mut nNumBondPosTmp = 0_i32;
    for k in 0..end_path_index {
        let path = DfsPath.get(k).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bond_type = path.bond_type;
        if bond_type == BOND_SINGLE as u8
            || bond_type == BOND_DOUBLE as u8
            || bond_type == BOND_ALTERN as u8
            || bond_type == BOND_ALT12NS as u8
        {
            let index = usize::try_from(nNumBondPosTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            BondPosTmp[index].nAtomNumber = path.at_no;
            BondPosTmp[index].neighbor_index = path.bond_pos as AT_RANK;
            nNumBondPosTmp = nNumBondPosTmp.wrapping_add(2);
        }
    }

    let mut EndPointTmp: [T_ENDPOINT; 2] = std::array::from_fn(|_| T_ENDPOINT::default());
    let mut nNumEndPointTmp = 0_i32;
    for (j, endpoint) in [at1, at2].into_iter().enumerate() {
        let temporary_index = usize::try_from(nNumEndPointTmp).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let endpoint_index = usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let known_endpoint = heap
            .slice(atom.as_const())?
            .get(endpoint_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint;
        if known_endpoint == 0 {
            let endpoint_valence = if j != 0 { endpoint_valence2 } else { endpoint_valence1 };
            let chem_bonds_valence = if j != 0 {
                i32::from(eif2.cNeutralBondsValence)
            } else {
                i32::from(eif1.cNeutralBondsValence)
            };
            let mobile = if j != 0 { nMobile2 } else { nMobile1 };
            if mobile + chem_bonds_valence != endpoint_valence {
                return Ok(ret);
            }
            let atoms = heap.slice(atom.as_const())?;
            AddAtom2num(&mut EndPointTmp[temporary_index].num, atoms, endpoint, 2)?;
            AddAtom2DA(&mut EndPointTmp[temporary_index].num_DA, atoms, endpoint, 2)?;
        } else {
            EndPointTmp[temporary_index] = T_ENDPOINT::default();
        }
        EndPointTmp[temporary_index].nAtomNumber = endpoint as AT_RANK;
        EndPointTmp[temporary_index].nGroupNumber = known_endpoint;
        EndPointTmp[temporary_index].nEquNumber = 0;
        nNumEndPointTmp = nNumEndPointTmp.wrapping_add(1);
    }

    {
        let atoms = heap.slice(atom.as_const())?;
        nNumBondPos = AddBondsPos(
            atoms,
            &mut BondPosTmp,
            nNumBondPosTmp,
            BondPos,
            nMaxNumBondPos,
            nNumBondPos,
        )?;
        nNumEndPoint = AddEndPoints(&EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint)?;
    }
    if nNumBondPos >= 0 && nNumEndPoint >= 0 {
        ret = i32::from(nNumBondPos > *pnNumBondPos || nNumEndPoint > *pnNumEndPoint);
        if ret != 0 {
            *pnNumBondPos = nNumBondPos;
            *pnNumEndPoint = nNumEndPoint;
        }
    }

    Ok(ret)
}

pub(crate) struct DfsRingCallbackArgs<'a> {
    heap: &'a mut SourceHeap,
    pCG: &'a mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &'a mut [DFS_PATH],
    nLenDfsPath: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    EndPoint: &'a mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &'a mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &'a mut i32,
    pnNumBondPos: &'a mut i32,
    pBNS: Option<&'a mut BalancedNetworkStructure>,
    pBD: Option<&'a mut BalancedNetworkData>,
    num_atoms: i32,
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn DFS_FindTautInARing<CheckRing, CheckCenter>(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    nStartAtom: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    mut nCycleLen: i32,
    nDfsPathPos: SourceMutPointer<AT_RANK>,
    DfsPath: &mut [DFS_PATH],
    CheckDfsRing: &mut CheckRing,
    CheckCenterPoint: &mut CheckCenter,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    mut pBNS: Option<&mut BalancedNetworkStructure>,
    mut pBD: Option<&mut BalancedNetworkData>,
    num_atoms: i32,
) -> Result<i32, SourceHeapError>
where
    CheckRing: for<'a> FnMut(DfsRingCallbackArgs<'a>) -> Result<i32, SourceHeapError>,
    CheckCenter: FnMut(&SourceHeap, SourceMutPointer<inp_ATOM>, i32) -> Result<i32, SourceHeapError>,
{
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:458 DFS_FindTautInARing
    // INCHI✔️❌: int DFS_FindTautInARing( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                          inp_ATOM *atom,
    // INCHI✔️❌:                          int nStartAtom,
    // INCHI✔️❌:                          int nStartAtomNeighbor,
    // INCHI✔️❌:                          int nStartAtomNeighbor2,
    // INCHI✔️❌:                          int nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                          int nCycleLen,
    // INCHI✔️❌:                          AT_RANK  *nDfsPathPos,
    // INCHI✔️❌:                          DFS_PATH *DfsPath,
    // INCHI✔️❌:                          CHECK_DFS_RING *CheckDfsRing,
    // INCHI✔️❌:                          CHECK_CENTERPOINT *CheckCenterPoint,
    // INCHI✔️❌:                          T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                          int nMaxNumEndPoint,
    // INCHI✔️❌:                          T_BONDPOS  *BondPos,
    // INCHI✔️❌:                          int nMaxNumBondPos,
    // INCHI✔️❌:                          int *pnNumEndPoint,
    // INCHI✔️❌:                          int *pnNumBondPos,
    // INCHI✔️❌:                          struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                          struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                          int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  Depth First Search */
    // INCHI✔️❌:     /*  Ignore all atoms not belonging to the current ring system (=biconnected component) */
    // INCHI✔️❌:     AT_RANK      nMinLenDfsPath;
    // INCHI✔️❌:     int          j, cur_at, nxt_at, prv_at;
    // INCHI✔️❌:     int          nLenDfsPath, nNumFound, ret;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int          nDoNotTouchAtom1 = -1, nDoNotTouchAtom2 = -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLenDfsPath = 0;
    // INCHI✔️❌:     nNumFound = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     nCycleLen--;
    // INCHI✔️❌:
    // INCHI✔️❌:     DfsPath[nLenDfsPath].at_no = cur_at = nStartAtom;
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_type = 0;
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_pos = -1;
    // INCHI✔️❌:     nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark */
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     nMinLenDfsPath = 0;
    // INCHI✔️❌:     if (nStartAtomNeighbor2 >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nDoNotTouchAtom1 = (int) atom[cur_at].neighbor[nStartAtomNeighbor2];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add the first neighbor to the 2nd tree position if required */
    // INCHI✔️❌:     if (nStartAtomNeighbor >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = nStartAtomNeighbor;
    // INCHI✔️❌:         prv_at = cur_at;
    // INCHI✔️❌:         cur_at = atom[prv_at].neighbor[j];
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_type = ( atom[prv_at].bond_type[j] & ~BOND_MARK_ALL );
    // INCHI✔️❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, prv_at, j, DfsPath[nLenDfsPath].bond_type );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_pos = j;
    // INCHI✔️❌:
    // INCHI✔️❌:         nLenDfsPath++;
    // INCHI✔️❌:
    // INCHI✔️❌:         DfsPath[nLenDfsPath].at_no = cur_at;
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_type = 0;
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_pos = -1;
    // INCHI✔️❌:         nDfsPathPos[cur_at] = nLenDfsPath + 1;
    // INCHI✔️❌:         nMinLenDfsPath++;
    // INCHI✔️❌:         if (nStartAtomNeighborNeighbor >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nDoNotTouchAtom2 = (int) atom[cur_at].neighbor[nStartAtomNeighborNeighbor];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  MAIN DFS CYCLE: may find one and the same t-group 2 times; saves only one instance */
    // INCHI✔️❌:     /*  traverse *all* paths starting at atom[nStartAtom]; max. path length = (nCycleLen+1)  */
    // INCHI✔️❌:     while (nLenDfsPath >= nMinLenDfsPath)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = ++DfsPath[nLenDfsPath].bond_pos;
    // INCHI✔️❌:         if (j < atom[cur_at = (int) DfsPath[nLenDfsPath].at_no].valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             DfsPath[nLenDfsPath].bond_type = ( atom[cur_at].bond_type[j] & ~BOND_MARK_ALL );
    // INCHI✔️❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI✔️❌:             DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, cur_at, j, DfsPath[nLenDfsPath].bond_type );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             nxt_at = (int) atom[cur_at].neighbor[j];
    // INCHI✔️❌:             if (nxt_at == nDoNotTouchAtom1 ||
    // INCHI✔️❌:                 nxt_at == nDoNotTouchAtom2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ; /*  ignore */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nDfsPathPos[nxt_at])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  found a ring closure or a step backwards */
    // INCHI✔️❌:                     if (1 == nDfsPathPos[nxt_at] && nLenDfsPath == nCycleLen)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*  we have found the cycle; check it */
    // INCHI✔️❌:                         ret = ( *CheckDfsRing )( pCG,
    // INCHI✔️❌:                             atom,
    // INCHI✔️❌:                             DfsPath, nLenDfsPath,
    // INCHI✔️❌:                             nStartAtomNeighbor,
    // INCHI✔️❌:                             nStartAtomNeighbor2,
    // INCHI✔️❌:                             nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                             EndPoint, nMaxNumEndPoint,
    // INCHI✔️❌:                             BondPos, nMaxNumBondPos,
    // INCHI✔️❌:                             pnNumEndPoint, pnNumBondPos,
    // INCHI✔️❌:                             pBNS, pBD, num_atoms );
    // INCHI✔️❌:                         if (ret < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nNumFound = ret;
    // INCHI✔️❌:                             goto clear_path;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         nNumFound += ret;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!( *CheckCenterPoint )( atom, nxt_at ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ; /*  cannot advance to a non-centerpoint; ignore */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nLenDfsPath < nCycleLen)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  advance */
    // INCHI✔️❌:                             nLenDfsPath++;
    // INCHI✔️❌:                             cur_at = nxt_at;
    // INCHI✔️❌:                             DfsPath[nLenDfsPath].at_no = cur_at;
    // INCHI✔️❌:                             DfsPath[nLenDfsPath].bond_type = 0;
    // INCHI✔️❌:                             DfsPath[nLenDfsPath].bond_pos = -1;
    // INCHI✔️❌:                             nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  retract */
    // INCHI✔️❌:             nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
    // INCHI✔️❌:             nLenDfsPath--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: clear_path:
    // INCHI✔️❌:     while (0 <= nLenDfsPath)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
    // INCHI✔️❌:         nLenDfsPath--;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumFound;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DFS_FindTautInARing
    // BEGIN INCHI ACTIVE MACRO AND TYPE CONFIGURATION: DFS_FindTautInARing
    // INCHI✔️❌: #define FIND_RING_SYSTEMS 1
    // INCHI✔️❌: #define FIX_BOND23_IN_TAUT 0
    // INCHI✔️❌: #define BOND_MARK_ALL 0xf0
    // INCHI✔️❌: typedef unsigned short AT_RANK;
    // INCHI✔️❌: typedef signed char S_CHAR;
    // INCHI✔️❌: The two ACTUAL_ORDER blocks are inactive after preprocessing.
    // END INCHI ACTIVE MACRO AND TYPE CONFIGURATION: DFS_FindTautInARing
    // SourceHeap preserves C pointer identity but adds BTreeMap lookups for atom and path marks.

    let mut nLenDfsPath = 0_i32;
    let mut nNumFound = 0_i32;
    let mut nDoNotTouchAtom1 = -1_i32;
    let mut nDoNotTouchAtom2 = -1_i32;

    nCycleLen = nCycleLen.wrapping_sub(1);

    let mut cur_at = nStartAtom;
    let path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let path = DfsPath.get_mut(path_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    path.at_no = nStartAtom as AT_RANK;
    path.bond_type = 0;
    path.bond_pos = -1;
    let atom_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(nDfsPathPos)?
        .get_mut(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = nLenDfsPath.wrapping_add(1) as AT_RANK;

    let mut nMinLenDfsPath = 0_u16;
    if nStartAtomNeighbor2 >= 0 {
        let neighbor_index = usize::try_from(nStartAtomNeighbor2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        nDoNotTouchAtom1 = i32::from(
            *heap
                .slice(atom.as_const())?
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .neighbor
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
    }

    if nStartAtomNeighbor >= 0 {
        let j = usize::try_from(nStartAtomNeighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let prv_at = cur_at;
        let previous_index = usize::try_from(prv_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let (next_atom, bond_type) = {
            let previous = heap
                .slice(atom.as_const())?
                .get(previous_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (
                *previous.neighbor.get(j).ok_or(SourceHeapError::PointerOutOfBounds)?,
                *previous.bond_type.get(j).ok_or(SourceHeapError::PointerOutOfBounds)? & !0xf0_u8,
            )
        };
        cur_at = i32::from(next_atom);
        let path = DfsPath.get_mut(path_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        path.bond_type = bond_type;
        path.bond_pos = nStartAtomNeighbor as i8;

        nLenDfsPath = nLenDfsPath.wrapping_add(1);
        let path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let path = DfsPath.get_mut(path_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        path.at_no = cur_at as AT_RANK;
        path.bond_type = 0;
        path.bond_pos = -1;
        let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(nDfsPathPos)?
            .get_mut(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = nLenDfsPath.wrapping_add(1) as AT_RANK;
        nMinLenDfsPath = nMinLenDfsPath.wrapping_add(1);
        if nStartAtomNeighborNeighbor >= 0 {
            let neighbor_index =
                usize::try_from(nStartAtomNeighborNeighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            nDoNotTouchAtom2 = i32::from(
                *heap
                    .slice(atom.as_const())?
                    .get(current_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .neighbor
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
        }
    }

    while nLenDfsPath >= i32::from(nMinLenDfsPath) {
        let path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let path = DfsPath.get_mut(path_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        path.bond_pos = path.bond_pos.wrapping_add(1);
        let j = i32::from(path.bond_pos);
        cur_at = i32::from(path.at_no);
        let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let valence = i32::from(
            heap.slice(atom.as_const())?
                .get(current_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .valence,
        );
        if j < valence {
            let neighbor_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let (bond_type, next_atom) = {
                let current = heap
                    .slice(atom.as_const())?
                    .get(current_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                (
                    *current
                        .bond_type
                        .get(neighbor_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        & !0xf0_u8,
                    *current
                        .neighbor
                        .get(neighbor_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let nxt_at = i32::from(next_atom);
            DfsPath
                .get_mut(path_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bond_type = bond_type;

            if nxt_at != nDoNotTouchAtom1 && nxt_at != nDoNotTouchAtom2 {
                let next_index = usize::try_from(nxt_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let next_path_position = *heap
                    .slice(nDfsPathPos.as_const())?
                    .get(next_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if next_path_position != 0 {
                    if next_path_position == 1 && nLenDfsPath == nCycleLen {
                        let ret = CheckDfsRing(DfsRingCallbackArgs {
                            heap,
                            pCG,
                            atom,
                            DfsPath,
                            nLenDfsPath,
                            nStartAtomNeighbor,
                            nStartAtomNeighbor2,
                            nStartAtomNeighborNeighbor,
                            EndPoint,
                            nMaxNumEndPoint,
                            BondPos,
                            nMaxNumBondPos,
                            pnNumEndPoint,
                            pnNumBondPos,
                            pBNS: pBNS.as_deref_mut(),
                            pBD: pBD.as_deref_mut(),
                            num_atoms,
                        })?;
                        if ret < 0 {
                            nNumFound = ret;
                            while nLenDfsPath >= 0 {
                                let clear_index =
                                    usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let atom_to_clear = i32::from(
                                    DfsPath
                                        .get(clear_index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                                        .at_no,
                                );
                                let atom_to_clear =
                                    usize::try_from(atom_to_clear).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                *heap
                                    .slice_mut(nDfsPathPos)?
                                    .get_mut(atom_to_clear)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                                nLenDfsPath = nLenDfsPath.wrapping_sub(1);
                            }
                            return Ok(nNumFound);
                        }
                        nNumFound = nNumFound.wrapping_add(ret);
                    }
                } else if CheckCenterPoint(heap, atom, nxt_at)? != 0 && nLenDfsPath < nCycleLen {
                    nLenDfsPath = nLenDfsPath.wrapping_add(1);
                    cur_at = nxt_at;
                    let next_path_index =
                        usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let next_path = DfsPath
                        .get_mut(next_path_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    next_path.at_no = cur_at as AT_RANK;
                    next_path.bond_type = 0;
                    next_path.bond_pos = -1;
                    let next_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    *heap
                        .slice_mut(nDfsPathPos)?
                        .get_mut(next_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = nLenDfsPath.wrapping_add(1) as AT_RANK;
                }
            }
        } else {
            let atom_to_clear = i32::from(
                DfsPath
                    .get(path_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at_no,
            );
            let atom_to_clear = usize::try_from(atom_to_clear).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(nDfsPathPos)?
                .get_mut(atom_to_clear)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            nLenDfsPath = nLenDfsPath.wrapping_sub(1);
        }
    }

    while nLenDfsPath >= 0 {
        let clear_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom_to_clear = i32::from(
            DfsPath
                .get(clear_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .at_no,
        );
        let atom_to_clear = usize::try_from(atom_to_clear).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(nDfsPathPos)?
            .get_mut(atom_to_clear)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        nLenDfsPath = nLenDfsPath.wrapping_sub(1);
    }

    Ok(nNumFound)
}

pub(crate) struct DfsPathCallbackArgs<'a> {
    heap: &'a mut SourceHeap,
    pCG: &'a mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &'a mut [DFS_PATH],
    nLenDfsPath: i32,
    jNxtNeigh: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    EndPoint: &'a mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &'a mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &'a mut i32,
    pnNumBondPos: &'a mut i32,
    pBNS: Option<&'a mut BalancedNetworkStructure>,
    pBD: Option<&'a mut BalancedNetworkData>,
    num_atoms: i32,
}

pub(crate) struct DfsPathCenterCallbackArgs<'a> {
    heap: &'a mut SourceHeap,
    atom: SourceMutPointer<inp_ATOM>,
    DfsPath: &'a mut [DFS_PATH],
    nLenDfsPath: i32,
    jNxtNeigh: i32,
    pBNS: Option<&'a mut BalancedNetworkStructure>,
    pBD: Option<&'a mut BalancedNetworkData>,
    num_atoms: i32,
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn DFS_FindTautAltPath<CheckPath, CheckCenter>(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    atom: SourceMutPointer<inp_ATOM>,
    nStartAtom: i32,
    nStartAtomNeighbor: i32,
    nStartAtomNeighbor2: i32,
    nStartAtomNeighborNeighbor: i32,
    mut nCycleLen: i32,
    nDfsPathPos: SourceMutPointer<AT_RANK>,
    DfsPath: &mut [DFS_PATH],
    CheckDfsPath: &mut CheckPath,
    CheckCenterPointPath: &mut CheckCenter,
    EndPoint: &mut [T_ENDPOINT],
    nMaxNumEndPoint: i32,
    BondPos: &mut [T_BONDPOS],
    nMaxNumBondPos: i32,
    pnNumEndPoint: &mut i32,
    pnNumBondPos: &mut i32,
    mut pBNS: Option<&mut BalancedNetworkStructure>,
    mut pBD: Option<&mut BalancedNetworkData>,
    num_atoms: i32,
) -> Result<i32, SourceHeapError>
where
    CheckPath: for<'a> FnMut(DfsPathCallbackArgs<'a>) -> Result<i32, SourceHeapError>,
    CheckCenter: for<'a> FnMut(DfsPathCenterCallbackArgs<'a>) -> Result<i32, SourceHeapError>,
{
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:616 DFS_FindTautAltPath
    // INCHI✔️❌: int DFS_FindTautAltPath( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                          inp_ATOM *atom,
    // INCHI✔️❌:                          int nStartAtom,
    // INCHI✔️❌:                          int nStartAtomNeighbor,
    // INCHI✔️❌:                          int nStartAtomNeighbor2,
    // INCHI✔️❌:                          int nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                          int nCycleLen,
    // INCHI✔️❌:                          AT_RANK  *nDfsPathPos,
    // INCHI✔️❌:                          DFS_PATH *DfsPath,
    // INCHI✔️❌:                          CHECK_DFS_PATH *CheckDfsPath,
    // INCHI✔️❌:                          CHECK_DFS_CENTERPOINT *CheckCenterPointPath,
    // INCHI✔️❌:                          T_ENDPOINT *EndPoint,
    // INCHI✔️❌:                          int nMaxNumEndPoint,
    // INCHI✔️❌:                          T_BONDPOS  *BondPos,
    // INCHI✔️❌:                          int nMaxNumBondPos,
    // INCHI✔️❌:                          int *pnNumEndPoint,
    // INCHI✔️❌:                          int *pnNumBondPos,
    // INCHI✔️❌:                          struct BalancedNetworkStructure *pBNS,
    // INCHI✔️❌:                          struct BalancedNetworkData *pBD,
    // INCHI✔️❌:                          int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  Naive Depth First Search: same atom may be approached along different alt paths */
    // INCHI✔️❌:     /*  Ignore all atoms not belonging to the current ring system (=biconnected component) */
    // INCHI✔️❌:     AT_RANK      nMinLenDfsPath;
    // INCHI✔️❌:     int          j, cur_at, nxt_at, prv_at;
    // INCHI✔️❌:     int          nLenDfsPath, nNumFound, ret;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int          nDoNotTouchAtom1 = -1, nDoNotTouchAtom2 = -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     nLenDfsPath = 0;
    // INCHI✔️❌:     nNumFound = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     nCycleLen--; /* indef of the last atom in the alt path, statring from 0 */
    // INCHI✔️❌:
    // INCHI✔️❌:     DfsPath[nLenDfsPath].at_no = cur_at = nStartAtom;
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_type = 0;
    // INCHI✔️❌:     DfsPath[nLenDfsPath].bond_pos = -1;  /* initialize index of the bond to the next atom */
    // INCHI✔️❌:     nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark with distance + 1 */
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables/code */
    // INCHI✔️❌:     nMinLenDfsPath = 0;  /* allow to restart from nStartAtom */
    // INCHI✔️❌:     if (nStartAtomNeighbor2 >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nDoNotTouchAtom1 = (int) atom[cur_at].neighbor[nStartAtomNeighbor2];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  add the first neighbor to the 2nd tree position if required */
    // INCHI✔️❌:     if (nStartAtomNeighbor >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = nStartAtomNeighbor;
    // INCHI✔️❌:         prv_at = cur_at;
    // INCHI✔️❌:         cur_at = atom[prv_at].neighbor[j];
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_type = ( atom[prv_at].bond_type[j] & ~BOND_MARK_ALL );
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:         DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, prv_at, j, DfsPath[nLenDfsPath].bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_pos = j; /* fix index of the bond to the next atom */
    // INCHI✔️❌:
    // INCHI✔️❌:         nLenDfsPath++;
    // INCHI✔️❌:
    // INCHI✔️❌:         DfsPath[nLenDfsPath].at_no = cur_at;
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_type = 0;
    // INCHI✔️❌:         DfsPath[nLenDfsPath].bond_pos = -1;
    // INCHI✔️❌:         nDfsPathPos[cur_at] = nLenDfsPath + 1; /* mark with distance + 1 */
    // INCHI✔️❌:         nMinLenDfsPath++;                 /* allow to restart from nStartAtom's neighbor */
    // INCHI✔️❌:         if (nStartAtomNeighborNeighbor >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nDoNotTouchAtom2 = (int) atom[cur_at].neighbor[nStartAtomNeighborNeighbor];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  MAIN DFS CYCLE: may find one and the same t-group 2 times; saves only one instance */
    // INCHI✔️❌:     /*  traverse *all* paths starting at atom[nStartAtom]; max. path length = (nCycleLen+1)  */
    // INCHI✔️❌:     while (nLenDfsPath >= nMinLenDfsPath)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = ++DfsPath[nLenDfsPath].bond_pos;
    // INCHI✔️❌:         if (j < atom[cur_at = (int) DfsPath[nLenDfsPath].at_no].valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             DfsPath[nLenDfsPath].bond_type = ( atom[cur_at].bond_type[j] & ~BOND_MARK_ALL );
    // INCHI❌❌: #if ( FIX_BOND23_IN_TAUT == 1 )
    // INCHI❌❌:             DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, cur_at, j, DfsPath[nLenDfsPath].bond_type );
    // INCHI❌❌: #endif
    // INCHI✔️❌:             nxt_at = (int) atom[cur_at].neighbor[j];
    // INCHI✔️❌:             if (nxt_at == nDoNotTouchAtom1 || /* forbidden */
    // INCHI✔️❌:                 nxt_at == nDoNotTouchAtom2 || /* forbidden */
    // INCHI✔️❌:                 nDfsPathPos[nxt_at] || /* ring closure */
    // INCHI✔️❌:                 (nLenDfsPath && nxt_at == (int) DfsPath[nLenDfsPath - 1].at_no) /* step backwards */
    // INCHI✔️❌:                 ) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ; /* ignore nxt_at */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nLenDfsPath == nCycleLen &&
    // INCHI✔️❌:                        /* 1,5 and at least one of the endpoints is not in a ring */
    // INCHI✔️❌:                     ( atom[nxt_at].nNumAtInRingSystem == 1 || atom[nStartAtom].nNumAtInRingSystem == 1 ) &&
    // INCHI✔️❌:                        /*  we have found the alt path of the requested length; check it */
    // INCHI✔️❌:                        /* calling Check15TautPath() */
    // INCHI✔️❌:                        ( ret = ( *CheckDfsPath )( pCG,
    // INCHI✔️❌:                                                   atom,
    // INCHI✔️❌:                                                   DfsPath, nLenDfsPath,
    // INCHI✔️❌:                                                   j, nStartAtomNeighbor,
    // INCHI✔️❌:                                                   nStartAtomNeighbor2,
    // INCHI✔️❌:                                                   nStartAtomNeighborNeighbor,
    // INCHI✔️❌:                                                   EndPoint, nMaxNumEndPoint, BondPos, nMaxNumBondPos,
    // INCHI✔️❌:                                                   pnNumEndPoint, pnNumBondPos,
    // INCHI✔️❌:                                                   pBNS, pBD, num_atoms ) ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ret < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nNumFound = ret;
    // INCHI✔️❌:                         goto clear_path; /* program error */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     nNumFound += ret; /* success */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* calling Check15TautPathCenterpoint() */
    // INCHI✔️❌:                     if (!( *CheckCenterPointPath )( atom, DfsPath, nLenDfsPath, j,
    // INCHI✔️❌:                                                     pBNS, pBD, num_atoms ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ; /*  cannot advance to a non-centerpoint; ignore */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nLenDfsPath < nCycleLen)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  advance */
    // INCHI✔️❌:                             nLenDfsPath++;
    // INCHI✔️❌:                             cur_at = nxt_at;
    // INCHI✔️❌:                             DfsPath[nLenDfsPath].at_no = cur_at;
    // INCHI✔️❌:                             DfsPath[nLenDfsPath].bond_type = 0;
    // INCHI✔️❌:                             DfsPath[nLenDfsPath].bond_pos = -1;
    // INCHI✔️❌:                             nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  retract */
    // INCHI✔️❌:             nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
    // INCHI✔️❌:             nLenDfsPath--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: clear_path:
    // INCHI✔️❌:     while (0 <= nLenDfsPath)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
    // INCHI✔️❌:         nLenDfsPath--;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumFound;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DFS_FindTautAltPath
    // BEGIN INCHI ACTIVE MACRO AND TYPE CONFIGURATION: DFS_FindTautAltPath
    // INCHI✔️❌: #define TAUT_15_NON_RING 1
    // INCHI✔️❌: #define FIX_BOND23_IN_TAUT 0
    // INCHI✔️❌: #define BOND_MARK_ALL 0xf0
    // INCHI✔️❌: typedef unsigned short AT_RANK;
    // INCHI✔️❌: typedef signed char S_CHAR;
    // INCHI✔️❌: The two ACTUAL_ORDER blocks are inactive after preprocessing.
    // END INCHI ACTIVE MACRO AND TYPE CONFIGURATION: DFS_FindTautAltPath
    // SourceHeap preserves C pointer identity but adds BTreeMap lookups for atoms and path marks.

    let mut nLenDfsPath = 0_i32;
    let mut nNumFound = 0_i32;
    let mut nDoNotTouchAtom1 = -1_i32;
    let mut nDoNotTouchAtom2 = -1_i32;

    nCycleLen = nCycleLen.wrapping_sub(1);

    let mut cur_at = nStartAtom;
    let path = DfsPath.get_mut(0).ok_or(SourceHeapError::PointerOutOfBounds)?;
    path.at_no = nStartAtom as AT_RANK;
    path.bond_type = 0;
    path.bond_pos = -1;
    let start_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(nDfsPathPos)?
        .get_mut(start_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 1;

    let mut nMinLenDfsPath = 0_u16;
    if nStartAtomNeighbor2 >= 0 {
        let neighbor_index = usize::try_from(nStartAtomNeighbor2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        nDoNotTouchAtom1 = i32::from(
            *heap
                .slice(atom.as_const())?
                .get(start_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .neighbor
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
    }

    if nStartAtomNeighbor >= 0 {
        let j = usize::try_from(nStartAtomNeighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let prv_at = cur_at;
        let previous_index = usize::try_from(prv_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let (next_atom, bond_type) = {
            let previous = heap
                .slice(atom.as_const())?
                .get(previous_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (
                *previous.neighbor.get(j).ok_or(SourceHeapError::PointerOutOfBounds)?,
                *previous.bond_type.get(j).ok_or(SourceHeapError::PointerOutOfBounds)? & !(BOND_MARK_ALL as u8),
            )
        };
        cur_at = i32::from(next_atom);
        let path = DfsPath.get_mut(0).ok_or(SourceHeapError::PointerOutOfBounds)?;
        path.bond_type = bond_type;
        path.bond_pos = nStartAtomNeighbor as i8;

        nLenDfsPath = nLenDfsPath.wrapping_add(1);
        let path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let path = DfsPath.get_mut(path_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        path.at_no = cur_at as AT_RANK;
        path.bond_type = 0;
        path.bond_pos = -1;
        let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(nDfsPathPos)?
            .get_mut(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = nLenDfsPath.wrapping_add(1) as AT_RANK;
        nMinLenDfsPath = nMinLenDfsPath.wrapping_add(1);
        if nStartAtomNeighborNeighbor >= 0 {
            let neighbor_index =
                usize::try_from(nStartAtomNeighborNeighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            nDoNotTouchAtom2 = i32::from(
                *heap
                    .slice(atom.as_const())?
                    .get(current_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .neighbor
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
        }
    }

    while nLenDfsPath >= i32::from(nMinLenDfsPath) {
        let path_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let path = DfsPath.get_mut(path_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        path.bond_pos = path.bond_pos.wrapping_add(1);
        let j = i32::from(path.bond_pos);
        cur_at = i32::from(path.at_no);
        let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let valence = i32::from(
            heap.slice(atom.as_const())?
                .get(current_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .valence,
        );
        if j < valence {
            let neighbor_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let (bond_type, next_atom) = {
                let current = heap
                    .slice(atom.as_const())?
                    .get(current_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                (
                    *current
                        .bond_type
                        .get(neighbor_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        & !(BOND_MARK_ALL as u8),
                    *current
                        .neighbor
                        .get(neighbor_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let nxt_at = i32::from(next_atom);
            DfsPath
                .get_mut(path_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bond_type = bond_type;

            let next_index = usize::try_from(nxt_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let next_path_position = *heap
                .slice(nDfsPathPos.as_const())?
                .get(next_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let is_step_backwards = if nLenDfsPath != 0 {
                let previous_path_index =
                    usize::try_from(nLenDfsPath.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                nxt_at
                    == i32::from(
                        DfsPath
                            .get(previous_path_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .at_no,
                    )
            } else {
                false
            };
            if nxt_at == nDoNotTouchAtom1 || nxt_at == nDoNotTouchAtom2 || next_path_position != 0 || is_step_backwards
            {
                continue;
            }

            let endpoint_is_outside_ring = {
                let atoms = heap.slice(atom.as_const())?;
                atoms
                    .get(next_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNumAtInRingSystem
                    == 1
                    || atoms
                        .get(start_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNumAtInRingSystem
                        == 1
            };
            let ret = if nLenDfsPath == nCycleLen && endpoint_is_outside_ring {
                CheckDfsPath(DfsPathCallbackArgs {
                    heap,
                    pCG,
                    atom,
                    DfsPath,
                    nLenDfsPath,
                    jNxtNeigh: j,
                    nStartAtomNeighbor,
                    nStartAtomNeighbor2,
                    nStartAtomNeighborNeighbor,
                    EndPoint,
                    nMaxNumEndPoint,
                    BondPos,
                    nMaxNumBondPos,
                    pnNumEndPoint,
                    pnNumBondPos,
                    pBNS: pBNS.as_deref_mut(),
                    pBD: pBD.as_deref_mut(),
                    num_atoms,
                })?
            } else {
                0
            };
            if ret != 0 {
                if ret < 0 {
                    nNumFound = ret;
                    while nLenDfsPath >= 0 {
                        let clear_index =
                            usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let atom_to_clear = usize::from(
                            DfsPath
                                .get(clear_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .at_no,
                        );
                        *heap
                            .slice_mut(nDfsPathPos)?
                            .get_mut(atom_to_clear)
                            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                        nLenDfsPath = nLenDfsPath.wrapping_sub(1);
                    }
                    return Ok(nNumFound);
                }
                nNumFound = nNumFound.wrapping_add(ret);
            } else {
                let center = CheckCenterPointPath(DfsPathCenterCallbackArgs {
                    heap,
                    atom,
                    DfsPath,
                    nLenDfsPath,
                    jNxtNeigh: j,
                    pBNS: pBNS.as_deref_mut(),
                    pBD: pBD.as_deref_mut(),
                    num_atoms,
                })?;
                if center != 0 && nLenDfsPath < nCycleLen {
                    nLenDfsPath = nLenDfsPath.wrapping_add(1);
                    cur_at = nxt_at;
                    let next_path_index =
                        usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let next_path = DfsPath
                        .get_mut(next_path_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    next_path.at_no = cur_at as AT_RANK;
                    next_path.bond_type = 0;
                    next_path.bond_pos = -1;
                    *heap
                        .slice_mut(nDfsPathPos)?
                        .get_mut(next_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = nLenDfsPath.wrapping_add(1) as AT_RANK;
                }
            }
        } else {
            let atom_to_clear = usize::from(
                DfsPath
                    .get(path_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at_no,
            );
            *heap
                .slice_mut(nDfsPathPos)?
                .get_mut(atom_to_clear)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            nLenDfsPath = nLenDfsPath.wrapping_sub(1);
        }
    }

    while nLenDfsPath >= 0 {
        let clear_index = usize::try_from(nLenDfsPath).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom_to_clear = usize::from(
            DfsPath
                .get(clear_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .at_no,
        );
        *heap
            .slice_mut(nDfsPathPos)?
            .get_mut(atom_to_clear)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        nLenDfsPath = nLenDfsPath.wrapping_sub(1);
    }

    Ok(nNumFound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::ichi_bns::{AllocateAndInitBnData, AllocateAndInitBnStruct};
    use crate::source_types::{BNS_ALT_PATH, BNS_EDGE, BNS_VERTEX, INCHI_CLOCK, RADICAL_DOUBLET};

    fn atom(element: u8, valence: i8, chem_bonds_valence: i8, num_h: i8, charge: i8, c_point: u16) -> inp_ATOM {
        inp_ATOM {
            el_number: element,
            valence,
            chem_bonds_valence,
            num_H: num_h,
            charge,
            c_point,
            ..inp_ATOM::default()
        }
    }

    #[test]
    fn source_port__ichiqueu__biscenterpointstrict__line_215() {
        fn check(atoms: Vec<inp_ATOM>, iat: i32) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            bIsCenterPointStrict(&heap, atoms, iat)
        }

        assert_eq!(check(vec![atom(7, 2, 2, 1, 0, 0)], 0), Ok(1));
        assert_eq!(check(vec![atom(7, 2, 2, 0, -1, 0)], 0), Ok(1));
        assert_eq!(check(vec![atom(7, 3, 3, 0, 0, 9)], 0), Ok(1));
        assert_eq!(check(vec![atom(7, 3, 3, 0, 0, 0)], 0), Ok(0));
        assert_eq!(check(vec![atom(0, 0, 0, 1, 0, 1)], 0), Ok(0));

        for element in [6, 7, 15, 33, 51] {
            assert_eq!(
                check(vec![atom(element, 2, 3, 0, 0, 0)], 0),
                Ok(1),
                "strict element {element}",
            );
        }
        for element in [0, 8, 16, 17, 35, 53, u8::MAX] {
            assert_eq!(
                check(vec![atom(element, 2, 3, 0, 0, 0)], 0),
                Ok(0),
                "non-strict element {element}",
            );
        }
        assert_eq!(check(vec![atom(6, 2, 4, 0, 0, 0)], 0), Ok(0));
        assert_eq!(check(vec![atom(6, i8::MAX, i8::MIN, 0, 0, 0)], 0), Ok(0));

        assert_eq!(
            check(vec![atom(6, 3, 3, 1, 0, 0)], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            check(vec![atom(6, 3, 3, 1, 0, 0)], 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiqueu__check15tautpathcenterpoint__line_1409() {
        fn check(next: inp_ATOM) -> Result<i32, SourceHeapError> {
            let mut atoms = vec![inp_ATOM::default(), next];
            atoms[0].neighbor[0] = 1;
            atoms[0].valence = 1;
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let path = vec![DFS_PATH {
                at_no: 0,
                bond_type: 0,
                bond_pos: -1,
            }];
            let mut bns = BalancedNetworkStructure::default();
            let mut bd = BalancedNetworkData::default();
            Check15TautPathCenterpoint(&heap, atoms, &path, 0, 0, Some(&mut bns), Some(&mut bd), 2)
        }

        let mut known_endpoint = atom(0, 0, 0, 0, 0, 0);
        known_endpoint.endpoint = 7;
        assert_eq!(check(known_endpoint), Ok(1));
        assert_eq!(check(atom(6, 2, 3, 0, 0, 0)), Ok(1));
        assert_eq!(check(atom(8, 2, 3, 0, 0, 0)), Ok(0));
        assert_eq!(check(atom(7, 2, 2, 1, 0, 0)), Ok(1));
        assert_eq!(check(atom(7, 3, 3, 0, 0, 0)), Ok(0));

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![inp_ATOM::default(); 2]).unwrap();
        let path = vec![DFS_PATH {
            at_no: 0,
            bond_type: 0,
            bond_pos: -1,
        }];
        assert_eq!(
            Check15TautPathCenterpoint(&heap, atoms, &path, -1, 0, None, None, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            Check15TautPathCenterpoint(&heap, atoms, &path, 1, 0, None, None, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            Check15TautPathCenterpoint(&heap, atoms, &path, 0, -1, None, None, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            Check15TautPathCenterpoint(&heap, atoms, &path, 0, 99, None, None, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let invalid_path = vec![DFS_PATH {
            at_no: 9,
            bond_type: 0,
            bond_pos: -1,
        }];
        assert_eq!(
            Check15TautPathCenterpoint(&heap, atoms, &invalid_path, 0, 0, None, None, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    fn add_bond(atoms: &mut [inp_ATOM], first: usize, second: usize, bond_type: u8) {
        let first_pos = usize::try_from(atoms[first].valence).unwrap();
        atoms[first].neighbor[first_pos] = second as u16;
        atoms[first].bond_type[first_pos] = bond_type;
        atoms[first].valence += 1;

        let second_pos = usize::try_from(atoms[second].valence).unwrap();
        atoms[second].neighbor[second_pos] = first as u16;
        atoms[second].bond_type[second_pos] = bond_type;
        atoms[second].valence += 1;
    }

    fn add_bond_with_positions(atoms: &mut [inp_ATOM], first: usize, second: usize, bond_type: u8) -> (i32, i32) {
        let first_pos = i32::from(atoms[first].valence);
        let second_pos = i32::from(atoms[second].valence);
        add_bond(atoms, first, second, bond_type);
        (first_pos, second_pos)
    }

    fn check7_bns(heap: &mut SourceHeap, atom_count: usize, left: usize, right: usize) -> BalancedNetworkStructure {
        let edge = heap
            .allocate_model_storage(vec![BNS_EDGE {
                neighbor1: left.min(right) as AT_RANK,
                neighbor12: (left as AT_RANK) ^ (right as AT_RANK),
                neigh_ord: [0, 0],
                ..BNS_EDGE::default()
            }])
            .unwrap();
        let mut vertices = Vec::with_capacity(atom_count);
        for index in 0..atom_count {
            let adjacency = if index == left || index == right {
                vec![0_i32]
            } else {
                Vec::new()
            };
            let num_adj_edges = adjacency.len() as AT_RANK;
            let iedge = heap.allocate_model_storage(adjacency).unwrap();
            vertices.push(BNS_VERTEX {
                num_adj_edges,
                iedge,
                ..BNS_VERTEX::default()
            });
        }
        let vert = heap.allocate_model_storage(vertices).unwrap();
        let mut bns = BalancedNetworkStructure {
            num_atoms: atom_count as i32,
            num_vertices: atom_count as i32,
            num_bonds: 1,
            num_edges: 1,
            max_vertices: atom_count as i32,
            max_edges: 1,
            max_iedges: 2,
            max_altp: 1,
            vert,
            edge,
            pbTautFlags: heap.allocate_model_storage(vec![0_u64]).unwrap(),
            ..BalancedNetworkStructure::default()
        };
        bns.altp[0] = heap.allocate_model_storage(vec![BNS_ALT_PATH::default(); 8]).unwrap();
        bns
    }

    struct Check15Fixture {
        heap: SourceHeap,
        atoms: SourceMutPointer<inp_ATOM>,
        positions: SourceMutPointer<AT_RANK>,
        path: Vec<DFS_PATH>,
        bns: BalancedNetworkStructure,
        bd: BalancedNetworkData,
    }

    fn check15_fixture(endpoint1: AT_RANK, endpoint2: AT_RANK) -> Check15Fixture {
        let mut atoms = vec![inp_ATOM::default(); 5];
        let bond_types = [
            BOND_SINGLE as u8,
            BOND_DOUBLE as u8,
            BOND_SINGLE as u8,
            BOND_DOUBLE as u8,
        ];
        let mut path = vec![DFS_PATH::default(); 5];
        for index in 0..4 {
            let (bond_pos, _) = add_bond_with_positions(&mut atoms, index, index + 1, bond_types[index]);
            if index < 3 {
                path[index] = DFS_PATH {
                    at_no: index as AT_RANK,
                    bond_type: bond_types[index],
                    bond_pos: bond_pos as i8,
                };
            }
        }
        path[3].at_no = 3;
        path[3].bond_pos = -1;

        atoms[0].el_number = 8;
        atoms[0].chem_bonds_valence = 1;
        atoms[0].num_H = 1;
        atoms[0].endpoint = endpoint1;
        atoms[0].nNumAtInRingSystem = 1;
        atoms[4].el_number = 8;
        atoms[4].chem_bonds_valence = 2;
        atoms[4].endpoint = endpoint2;
        atoms[4].nNumAtInRingSystem = 1;
        for atom in &mut atoms[1..4] {
            atom.el_number = 6;
            atom.chem_bonds_valence = 3;
        }

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let mut changed_bonds = 0;
        let bns_pointer = AllocateAndInitBnStruct(&mut heap, atoms, 5, 0, 0, 1, &mut changed_bonds).unwrap();
        assert_eq!(changed_bonds, 0);
        let mut bns = heap.slice(bns_pointer.as_const()).unwrap()[0].clone();
        bns.pbTautFlags = heap.allocate_model_storage(vec![0_u64]).unwrap();
        bns.ic = heap.allocate_model_storage(vec![INCHI_CLOCK::default()]).unwrap();
        let bd_pointer = AllocateAndInitBnData(&mut heap, bns.max_vertices).unwrap();
        let bd = heap.slice(bd_pointer.as_const()).unwrap()[0].clone();
        Check15Fixture {
            heap,
            atoms,
            positions,
            path,
            bns,
            bd,
        }
    }

    #[test]
    fn source_port__ichiqueu__check15tautpath__line_1434() {
        fn run(
            fixture: &mut Check15Fixture,
            path_len: i32,
            start_neighbor: i32,
            start_neighbor2: i32,
            start_neighbor_neighbor: i32,
            endpoints: &mut [T_ENDPOINT],
            max_endpoints: i32,
            bonds: &mut [T_BONDPOS],
            max_bonds: i32,
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            Check15TautPath(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                &mut fixture.path,
                path_len,
                1,
                start_neighbor,
                start_neighbor2,
                start_neighbor_neighbor,
                endpoints,
                max_endpoints,
                bonds,
                max_bonds,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                5,
                0,
            )
        }

        let mut accepted = check15_fixture(1, 2);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 4];
        let mut endpoint_count = 0;
        let mut bond_count = 0;
        assert_eq!(
            run(
                &mut accepted,
                3,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 4));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (
                    endpoint.nAtomNumber,
                    endpoint.nGroupNumber,
                    endpoint.nEquNumber,
                    endpoint.num,
                    endpoint.num_DA,
                ))
                .collect::<Vec<_>>(),
            vec![(0, 1, 0, [0; 5], [0; 6]), (4, 2, 0, [0; 5], [0; 6])]
        );
        assert_eq!(
            bonds
                .iter()
                .map(|bond| (bond.nAtomNumber, bond.neighbor_index))
                .collect::<Vec<_>>(),
            vec![(0, 0), (1, 1), (2, 1), (3, 1)]
        );
        assert_eq!(
            accepted.path[3],
            DFS_PATH {
                at_no: 3,
                bond_type: BOND_DOUBLE as u8,
                bond_pos: 1,
            }
        );
        assert_eq!(accepted.path[4].at_no, 4);

        let mut unknown = check15_fixture(0, 0);
        let mut unknown_endpoints = vec![T_ENDPOINT::default(); 2];
        let mut unknown_bonds = vec![T_BONDPOS::default(); 4];
        let mut unknown_endpoint_count = 0;
        let mut unknown_bond_count = 0;
        assert_eq!(
            run(
                &mut unknown,
                3,
                -1,
                -1,
                -1,
                &mut unknown_endpoints,
                2,
                &mut unknown_bonds,
                4,
                &mut unknown_endpoint_count,
                &mut unknown_bond_count,
            ),
            Ok(1)
        );
        assert_eq!((unknown_endpoint_count, unknown_bond_count), (2, 4));
        assert_eq!(unknown_endpoints[0].nAtomNumber, 0);
        assert_eq!(unknown_endpoints[0].nGroupNumber, 0);
        assert_eq!(unknown_endpoints[0].num[0], 1);
        assert_eq!(unknown_endpoints[1].nAtomNumber, 4);
        assert_eq!(unknown_endpoints[1].nGroupNumber, 0);

        for (path_len, n1, n2, n3) in [(2, -1, -1, -1), (3, 0, -1, -1), (3, -1, 0, -1), (3, -1, -1, 0)] {
            let mut invalid = check15_fixture(1, 2);
            let original_path = invalid.path.clone();
            let mut count1 = 7;
            let mut count2 = 9;
            assert_eq!(
                run(
                    &mut invalid,
                    path_len,
                    n1,
                    n2,
                    n3,
                    &mut endpoints,
                    2,
                    &mut bonds,
                    4,
                    &mut count1,
                    &mut count2,
                ),
                Ok(-1)
            );
            assert_eq!(invalid.path, original_path);
            assert_eq!((count1, count2), (7, 9));
        }

        let mut same_group = check15_fixture(3, 3);
        let mut same_endpoint_count = 0;
        let mut same_bond_count = 0;
        assert_eq!(
            run(
                &mut same_group,
                3,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut same_endpoint_count,
                &mut same_bond_count,
            ),
            Ok(0)
        );
        assert_eq!(same_group.path[3].bond_type, BOND_DOUBLE as u8);
        assert_eq!(same_group.path[4].at_no, 4);
        assert_eq!((same_endpoint_count, same_bond_count), (0, 0));

        for mutation in 0..4 {
            let mut rejected = check15_fixture(1, 2);
            match mutation {
                0 => rejected.path[0].bond_type = BOND_TRIPLE as u8,
                1 => rejected.path[1].bond_type = BOND_SINGLE as u8,
                2 => {
                    rejected.heap.slice_mut(rejected.atoms).unwrap()[0].el_number = 6;
                }
                3 => {
                    let first = &mut rejected.heap.slice_mut(rejected.atoms).unwrap()[0];
                    first.num_H = 0;
                    first.chem_bonds_valence = 2;
                    first.endpoint = 0;
                }
                _ => unreachable!(),
            }
            let mut count1 = 0;
            let mut count2 = 0;
            assert_eq!(
                run(
                    &mut rejected,
                    3,
                    -1,
                    -1,
                    -1,
                    &mut endpoints,
                    2,
                    &mut bonds,
                    4,
                    &mut count1,
                    &mut count2,
                ),
                Ok(0),
                "rejection branch {mutation}",
            );
            assert_eq!((count1, count2), (0, 0));
        }

        let mut no_alt_path = check15_fixture(0, 0);
        {
            let vertices = no_alt_path.heap.slice_mut(no_alt_path.bns.vert).unwrap();
            vertices[0].st_edge.flow = 0;
            vertices[4].st_edge.flow = 0;
        }
        let mut no_alt_endpoint_count = 0;
        let mut no_alt_bond_count = 0;
        assert_eq!(
            run(
                &mut no_alt_path,
                3,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut no_alt_endpoint_count,
                &mut no_alt_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((no_alt_endpoint_count, no_alt_bond_count), (0, 0));

        let mut overflow = check15_fixture(1, 2);
        let mut overflow_endpoints = vec![T_ENDPOINT::default(); 2];
        let mut overflow_bonds = vec![T_BONDPOS::default(); 4];
        let mut overflow_endpoint_count = 0;
        let mut overflow_bond_count = 0;
        assert_eq!(
            run(
                &mut overflow,
                3,
                -1,
                -1,
                -1,
                &mut overflow_endpoints,
                0,
                &mut overflow_bonds,
                0,
                &mut overflow_endpoint_count,
                &mut overflow_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((overflow_endpoint_count, overflow_bond_count), (0, 0));
        assert_eq!(overflow_endpoints[0].nAtomNumber, 0);
        assert_eq!(overflow_bonds[0].nAtomNumber, 0);

        let mut missing_path = check15_fixture(1, 2);
        missing_path.path.truncate(4);
        assert_eq!(
            run(
                &mut missing_path,
                3,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiqueu__nget15tautinaltpath__line_407() {
        fn run(
            fixture: &mut Check15Fixture,
            start_atom: i32,
            max_path: i32,
            endpoints: &mut [T_ENDPOINT],
            bonds: &mut [T_BONDPOS],
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            nGet15TautInAltPath(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                start_atom,
                fixture.positions,
                &mut fixture.path,
                max_path,
                endpoints,
                2,
                bonds,
                4,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                5,
                0,
            )
        }

        let mut accepted = check15_fixture(1, 2);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 4];
        let mut endpoint_count = 77;
        let mut bond_count = 88;
        assert_eq!(
            run(
                &mut accepted,
                0,
                8,
                &mut endpoints,
                &mut bonds,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 4));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (endpoint.nAtomNumber, endpoint.nGroupNumber))
                .collect::<Vec<_>>(),
            vec![(0, 1), (4, 2)]
        );
        assert_eq!(
            bonds
                .iter()
                .map(|bond| (bond.nAtomNumber, bond.neighbor_index))
                .collect::<Vec<_>>(),
            vec![(0, 0), (1, 1), (2, 1), (3, 1)]
        );
        assert_eq!(accepted.heap.slice(accepted.positions.as_const()).unwrap(), &[0; 5]);

        let mut too_short = check15_fixture(1, 2);
        too_short.heap.slice_mut(too_short.positions).unwrap().fill(13);
        let mut short_endpoint_count = 77;
        let mut short_bond_count = 88;
        assert_eq!(
            run(
                &mut too_short,
                0,
                7,
                &mut endpoints,
                &mut bonds,
                &mut short_endpoint_count,
                &mut short_bond_count,
            ),
            Ok(-1)
        );
        assert_eq!((short_endpoint_count, short_bond_count), (0, 0));
        assert_eq!(too_short.heap.slice(too_short.positions.as_const()).unwrap(), &[13; 5]);

        let mut no_path = check15_fixture(1, 2);
        no_path.heap.slice_mut(no_path.atoms).unwrap()[2].chem_bonds_valence = 2;
        let mut no_path_endpoint_count = 77;
        let mut no_path_bond_count = 88;
        assert_eq!(
            run(
                &mut no_path,
                0,
                8,
                &mut endpoints,
                &mut bonds,
                &mut no_path_endpoint_count,
                &mut no_path_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((no_path_endpoint_count, no_path_bond_count), (0, 0));
        assert_eq!(no_path.heap.slice(no_path.positions.as_const()).unwrap(), &[0; 5]);

        let mut invalid_start = check15_fixture(1, 2);
        let mut invalid_endpoint_count = 77;
        let mut invalid_bond_count = 88;
        assert_eq!(
            run(
                &mut invalid_start,
                5,
                8,
                &mut endpoints,
                &mut bonds,
                &mut invalid_endpoint_count,
                &mut invalid_bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((invalid_endpoint_count, invalid_bond_count), (0, 0));
    }

    struct Check5Fixture {
        heap: SourceHeap,
        atoms: SourceMutPointer<inp_ATOM>,
        positions: SourceMutPointer<AT_RANK>,
        path: Vec<DFS_PATH>,
        start_neighbor: i32,
        bns: BalancedNetworkStructure,
        bd: BalancedNetworkData,
    }

    fn check5_fixture(endpoint_group: AT_RANK) -> Check5Fixture {
        let mut atoms = vec![inp_ATOM::default(); 5];
        let (start_neighbor, _) = add_bond_with_positions(&mut atoms, 0, 1, BOND_SINGLE as u8);
        let mut path = vec![DFS_PATH::default(); 5];
        path[0].at_no = 0;
        for index in 1..=4 {
            let next = if index == 4 { 0 } else { index + 1 };
            let bond_type = if index % 2 == 1 {
                BOND_SINGLE as u8
            } else {
                BOND_DOUBLE as u8
            };
            let (bond_pos, _) = add_bond_with_positions(&mut atoms, index, next, bond_type);
            path[index] = DFS_PATH {
                at_no: index as AT_RANK,
                bond_type,
                bond_pos: bond_pos as i8,
            };
        }

        atoms[0].el_number = 7;
        atoms[0].chem_bonds_valence = 3;
        atoms[0].endpoint = endpoint_group;
        atoms[1].el_number = 7;
        atoms[1].chem_bonds_valence = 2;
        atoms[1].num_H = 1;
        atoms[1].endpoint = endpoint_group;
        for atom in &mut atoms[2..] {
            atom.el_number = 6;
            atom.chem_bonds_valence = 3;
        }

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let bns = check7_bns(&mut heap, 5, 0, 1);
        Check5Fixture {
            heap,
            atoms,
            positions,
            path,
            start_neighbor,
            bns,
            bd: BalancedNetworkData::default(),
        }
    }

    #[test]
    fn source_port__ichiqueu__check5membtautring__line_1666() {
        fn run(
            fixture: &mut Check5Fixture,
            path_len: i32,
            start_neighbor2: i32,
            start_neighbor_neighbor: i32,
            endpoints: &mut [T_ENDPOINT],
            max_endpoints: i32,
            bonds: &mut [T_BONDPOS],
            max_bonds: i32,
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            Check5MembTautRing(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                &fixture.path,
                path_len,
                fixture.start_neighbor,
                start_neighbor2,
                start_neighbor_neighbor,
                endpoints,
                max_endpoints,
                bonds,
                max_bonds,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                5,
                0,
            )
        }

        let mut accepted = check5_fixture(9);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 4];
        let mut endpoint_count = 0;
        let mut bond_count = 0;
        assert_eq!(
            run(
                &mut accepted,
                4,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 4));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (
                    endpoint.nAtomNumber,
                    endpoint.nGroupNumber,
                    endpoint.nEquNumber,
                    endpoint.num,
                    endpoint.num_DA,
                ))
                .collect::<Vec<_>>(),
            vec![(1, 9, 0, [0; 5], [0; 6]), (0, 9, 0, [0; 5], [0; 6])]
        );
        assert_eq!(
            bonds
                .iter()
                .map(|bond| (bond.nAtomNumber, bond.neighbor_index))
                .collect::<Vec<_>>(),
            vec![(1, 1), (2, 1), (3, 1), (4, 1)]
        );
        assert_eq!(
            run(
                &mut accepted,
                4,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(0)
        );
        assert_eq!((endpoint_count, bond_count), (2, 4));

        for (path_len, start_neighbor2, start_neighbor_neighbor) in [(3, -1, -1), (4, 0, -1), (4, -1, 0)] {
            let mut invalid = check5_fixture(9);
            let mut invalid_endpoint_count = 7;
            let mut invalid_bond_count = 8;
            assert_eq!(
                run(
                    &mut invalid,
                    path_len,
                    start_neighbor2,
                    start_neighbor_neighbor,
                    &mut endpoints,
                    2,
                    &mut bonds,
                    4,
                    &mut invalid_endpoint_count,
                    &mut invalid_bond_count,
                ),
                Ok(0)
            );
            assert_eq!((invalid_endpoint_count, invalid_bond_count), (7, 8));
        }

        let mut invalid_endpoint = check5_fixture(9);
        invalid_endpoint.heap.slice_mut(invalid_endpoint.atoms).unwrap()[0].radical = RADICAL_DOUBLET as i8;
        let mut rejected_endpoint_count = 0;
        let mut rejected_bond_count = 0;
        assert_eq!(
            run(
                &mut invalid_endpoint,
                4,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut invalid_mobile = check5_fixture(0);
        {
            let atoms = invalid_mobile.heap.slice_mut(invalid_mobile.atoms).unwrap();
            atoms[0].chem_bonds_valence = 2;
            atoms[0].num_H = 1;
        }
        assert_eq!(
            run(
                &mut invalid_mobile,
                4,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut rejected_alt_path = check5_fixture(0);
        assert_eq!(
            run(
                &mut rejected_alt_path,
                4,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((rejected_endpoint_count, rejected_bond_count), (0, 0));

        let mut invalid_bond = check5_fixture(9);
        invalid_bond.path[2].bond_type = BOND_TRIPLE as u8;
        assert_eq!(
            run(
                &mut invalid_bond,
                4,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut overflow = check5_fixture(9);
        let mut overflow_endpoints = vec![T_ENDPOINT::default(); 2];
        let mut overflow_bonds = vec![T_BONDPOS::default(); 4];
        let mut overflow_endpoint_count = 0;
        let mut overflow_bond_count = 0;
        assert_eq!(
            run(
                &mut overflow,
                4,
                -1,
                -1,
                &mut overflow_endpoints,
                0,
                &mut overflow_bonds,
                0,
                &mut overflow_endpoint_count,
                &mut overflow_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((overflow_endpoint_count, overflow_bond_count), (0, 0));
        assert_eq!(overflow_endpoints[0].nAtomNumber, 1);
        assert_eq!(overflow_bonds[0].nAtomNumber, 1);

        let mut missing_path = check5_fixture(9);
        missing_path.path.clear();
        assert_eq!(
            run(
                &mut missing_path,
                4,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                4,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiqueu__nget12tautin5membaltring__line_322() {
        fn run(
            fixture: &mut Check5Fixture,
            max_path: i32,
            endpoints: &mut [T_ENDPOINT],
            bonds: &mut [T_BONDPOS],
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            nGet12TautIn5MembAltRing(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                0,
                fixture.start_neighbor,
                fixture.positions,
                &mut fixture.path,
                max_path,
                endpoints,
                2,
                bonds,
                4,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                5,
                0,
            )
        }

        let mut accepted = check5_fixture(9);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 4];
        let mut endpoint_count = 77;
        let mut bond_count = 88;
        assert_eq!(
            run(
                &mut accepted,
                6,
                &mut endpoints,
                &mut bonds,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 4));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (endpoint.nAtomNumber, endpoint.nGroupNumber))
                .collect::<Vec<_>>(),
            vec![(1, 9), (0, 9)]
        );
        assert_eq!(
            bonds
                .iter()
                .map(|bond| (bond.nAtomNumber, bond.neighbor_index))
                .collect::<Vec<_>>(),
            vec![(1, 1), (2, 1), (3, 1), (4, 1)]
        );
        assert_eq!(accepted.heap.slice(accepted.positions.as_const()).unwrap(), &[0; 5]);
        assert_eq!(accepted.path[0].at_no, 0);
        assert_eq!(accepted.path[1].at_no, 1);

        let mut too_short = check5_fixture(9);
        too_short.heap.slice_mut(too_short.positions).unwrap().fill(13);
        let mut too_short_endpoint_count = 77;
        let mut too_short_bond_count = 88;
        assert_eq!(
            run(
                &mut too_short,
                5,
                &mut endpoints,
                &mut bonds,
                &mut too_short_endpoint_count,
                &mut too_short_bond_count,
            ),
            Ok(-1)
        );
        assert_eq!((too_short_endpoint_count, too_short_bond_count), (0, 0));
        assert_eq!(too_short.heap.slice(too_short.positions.as_const()).unwrap(), &[13; 5]);

        let mut no_ring = check5_fixture(9);
        {
            let atoms = no_ring.heap.slice_mut(no_ring.atoms).unwrap();
            atoms[2].el_number = 0;
            atoms[2].chem_bonds_valence = 0;
        }
        let mut no_ring_endpoint_count = 77;
        let mut no_ring_bond_count = 88;
        assert_eq!(
            run(
                &mut no_ring,
                6,
                &mut endpoints,
                &mut bonds,
                &mut no_ring_endpoint_count,
                &mut no_ring_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((no_ring_endpoint_count, no_ring_bond_count), (0, 0));
        assert_eq!(no_ring.heap.slice(no_ring.positions.as_const()).unwrap(), &[0; 5]);

        let mut missing_path = check5_fixture(9);
        missing_path.path.clear();
        let mut missing_endpoint_count = 77;
        let mut missing_bond_count = 88;
        assert_eq!(
            run(
                &mut missing_path,
                6,
                &mut endpoints,
                &mut bonds,
                &mut missing_endpoint_count,
                &mut missing_bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((missing_endpoint_count, missing_bond_count), (0, 0));
    }

    struct Check6Fixture {
        heap: SourceHeap,
        atoms: SourceMutPointer<inp_ATOM>,
        positions: SourceMutPointer<AT_RANK>,
        path: Vec<DFS_PATH>,
        bns: BalancedNetworkStructure,
        bd: BalancedNetworkData,
    }

    fn check6_fixture(endpoint_group: AT_RANK) -> Check6Fixture {
        let mut atoms = vec![inp_ATOM::default(); 7];
        let bond_types = [
            BOND_SINGLE as u8,
            BOND_DOUBLE as u8,
            BOND_SINGLE as u8,
            BOND_SINGLE as u8,
            BOND_DOUBLE as u8,
            BOND_SINGLE as u8,
        ];
        let mut path = vec![DFS_PATH::default(); 6];
        for index in 0..6 {
            let next = (index + 1) % 6;
            let (bond_pos, _) = add_bond_with_positions(&mut atoms, index, next, bond_types[index]);
            path[index] = DFS_PATH {
                at_no: index as AT_RANK,
                bond_type: bond_types[index],
                bond_pos: bond_pos as i8,
            };
        }
        add_bond_with_positions(&mut atoms, 3, 6, BOND_DOUBLE as u8);
        for atom in &mut atoms[..6] {
            atom.el_number = 6;
            atom.chem_bonds_valence = 3;
        }
        atoms[0].el_number = 7;
        atoms[0].chem_bonds_valence = 2;
        atoms[0].num_H = 1;
        atoms[0].endpoint = endpoint_group;
        atoms[3].chem_bonds_valence = 4;
        atoms[3].bCutVertex = 1;
        atoms[6].el_number = 8;
        atoms[6].chem_bonds_valence = 2;
        atoms[6].endpoint = endpoint_group;

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let positions = heap.allocate_model_storage(vec![0; 7]).unwrap();
        let bns = check7_bns(&mut heap, 7, 0, 6);
        Check6Fixture {
            heap,
            atoms,
            positions,
            path,
            bns,
            bd: BalancedNetworkData::default(),
        }
    }

    #[test]
    fn source_port__ichiqueu__check6membtautring__line_1142() {
        fn run(
            fixture: &mut Check6Fixture,
            path_len: i32,
            start_neighbor: i32,
            start_neighbor2: i32,
            start_neighbor_neighbor: i32,
            endpoints: &mut [T_ENDPOINT],
            max_endpoints: i32,
            bonds: &mut [T_BONDPOS],
            max_bonds: i32,
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            Check6MembTautRing(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                &fixture.path,
                path_len,
                start_neighbor,
                start_neighbor2,
                start_neighbor_neighbor,
                endpoints,
                max_endpoints,
                bonds,
                max_bonds,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                7,
                0,
            )
        }

        let mut accepted = check6_fixture(9);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 7];
        let mut endpoint_count = 0;
        let mut bond_count = 0;
        assert_eq!(
            run(
                &mut accepted,
                5,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                7,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 7));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (endpoint.nAtomNumber, endpoint.nGroupNumber))
                .collect::<Vec<_>>(),
            vec![(6, 9), (0, 9)]
        );
        assert_eq!(
            bonds.iter().map(|bond| bond.nAtomNumber).collect::<Vec<_>>(),
            vec![3, 2, 3, 1, 4, 0, 5]
        );

        for (path_len, start_neighbor, start_neighbor2, start_neighbor_neighbor) in
            [(4, -1, -1, -1), (5, 0, -1, -1), (5, -1, 0, -1), (5, -1, -1, 0)]
        {
            let mut invalid = check6_fixture(9);
            let mut invalid_endpoint_count = 7;
            let mut invalid_bond_count = 8;
            assert_eq!(
                run(
                    &mut invalid,
                    path_len,
                    start_neighbor,
                    start_neighbor2,
                    start_neighbor_neighbor,
                    &mut endpoints,
                    2,
                    &mut bonds,
                    7,
                    &mut invalid_endpoint_count,
                    &mut invalid_bond_count,
                ),
                Ok(-1)
            );
            assert_eq!((invalid_endpoint_count, invalid_bond_count), (7, 8));
        }

        let mut not_cut = check6_fixture(9);
        not_cut.heap.slice_mut(not_cut.atoms).unwrap()[3].bCutVertex = 0;
        let mut rejected_endpoint_count = 0;
        let mut rejected_bond_count = 0;
        assert_eq!(
            run(
                &mut not_cut,
                5,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                7,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut invalid_external_bond = check6_fixture(9);
        invalid_external_bond
            .heap
            .slice_mut(invalid_external_bond.atoms)
            .unwrap()[3]
            .bond_type[2] = BOND_TRIPLE as u8;
        assert_eq!(
            run(
                &mut invalid_external_bond,
                5,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                7,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut invalid_endpoint = check6_fixture(9);
        invalid_endpoint.heap.slice_mut(invalid_endpoint.atoms).unwrap()[6].radical = RADICAL_DOUBLET as i8;
        assert_eq!(
            run(
                &mut invalid_endpoint,
                5,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                7,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut rejected_alt_path = check6_fixture(0);
        assert_eq!(
            run(
                &mut rejected_alt_path,
                5,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                7,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut invalid_path_bond = check6_fixture(9);
        invalid_path_bond.path[1].bond_type = BOND_TRIPLE as u8;
        assert_eq!(
            run(
                &mut invalid_path_bond,
                5,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                7,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut overflow = check6_fixture(9);
        let mut overflow_endpoints = vec![T_ENDPOINT::default(); 2];
        let mut overflow_bonds = vec![T_BONDPOS::default(); 7];
        let mut overflow_endpoint_count = 0;
        let mut overflow_bond_count = 0;
        assert_eq!(
            run(
                &mut overflow,
                5,
                -1,
                -1,
                -1,
                &mut overflow_endpoints,
                0,
                &mut overflow_bonds,
                0,
                &mut overflow_endpoint_count,
                &mut overflow_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((overflow_endpoint_count, overflow_bond_count), (0, 0));
        assert_eq!(overflow_endpoints[0].nAtomNumber, 6);
        assert_eq!(overflow_bonds[0].nAtomNumber, 3);

        let mut missing_path = check6_fixture(9);
        missing_path.path.clear();
        assert_eq!(
            run(
                &mut missing_path,
                5,
                -1,
                -1,
                -1,
                &mut endpoints,
                2,
                &mut bonds,
                7,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiqueu__nget15tautin6membaltring__line_361() {
        fn run(
            fixture: &mut Check6Fixture,
            max_path: i32,
            endpoints: &mut [T_ENDPOINT],
            bonds: &mut [T_BONDPOS],
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            nGet15TautIn6MembAltRing(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                0,
                fixture.positions,
                &mut fixture.path,
                max_path,
                endpoints,
                2,
                bonds,
                7,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                7,
                0,
            )
        }

        let mut accepted = check6_fixture(9);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 7];
        let mut endpoint_count = 77;
        let mut bond_count = 88;
        assert_eq!(
            run(
                &mut accepted,
                8,
                &mut endpoints,
                &mut bonds,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 7));
        assert_eq!(accepted.heap.slice(accepted.positions.as_const()).unwrap(), &[0; 7]);

        let mut too_short = check6_fixture(9);
        too_short.heap.slice_mut(too_short.positions).unwrap().fill(13);
        let mut too_short_endpoint_count = 77;
        let mut too_short_bond_count = 88;
        assert_eq!(
            run(
                &mut too_short,
                7,
                &mut endpoints,
                &mut bonds,
                &mut too_short_endpoint_count,
                &mut too_short_bond_count,
            ),
            Ok(-1)
        );
        assert_eq!((too_short_endpoint_count, too_short_bond_count), (0, 0));
        assert_eq!(too_short.heap.slice(too_short.positions.as_const()).unwrap(), &[13; 7]);

        let mut no_ring = check6_fixture(9);
        no_ring.heap.slice_mut(no_ring.atoms).unwrap()[3].bCutVertex = 0;
        let mut no_ring_endpoint_count = 77;
        let mut no_ring_bond_count = 88;
        assert_eq!(
            run(
                &mut no_ring,
                8,
                &mut endpoints,
                &mut bonds,
                &mut no_ring_endpoint_count,
                &mut no_ring_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((no_ring_endpoint_count, no_ring_bond_count), (0, 0));
        assert_eq!(no_ring.heap.slice(no_ring.positions.as_const()).unwrap(), &[0; 7]);

        let mut missing_path = check6_fixture(9);
        missing_path.path.clear();
        let mut missing_endpoint_count = 77;
        let mut missing_bond_count = 88;
        assert_eq!(
            run(
                &mut missing_path,
                8,
                &mut endpoints,
                &mut bonds,
                &mut missing_endpoint_count,
                &mut missing_bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((missing_endpoint_count, missing_bond_count), (0, 0));
    }

    struct Check7Fixture {
        heap: SourceHeap,
        atoms: SourceMutPointer<inp_ATOM>,
        path: Vec<DFS_PATH>,
        start_neighbor2: i32,
        start_neighbor_neighbor: i32,
        first_endpoint: usize,
        second_endpoint: usize,
        bns: BalancedNetworkStructure,
        bd: BalancedNetworkData,
    }

    fn check7_fixture(path_len: usize, endpoint_group: AT_RANK) -> Check7Fixture {
        let first_endpoint = path_len + 1;
        let second_endpoint = path_len + 2;
        let atom_count = path_len + 3;
        let mut atoms = vec![inp_ATOM::default(); atom_count];

        let (start_neighbor_neighbor, _) = add_bond_with_positions(&mut atoms, 1, first_endpoint, BOND_DOUBLE as u8);
        let (start_neighbor2, _) = add_bond_with_positions(&mut atoms, 0, second_endpoint, BOND_SINGLE as u8);

        let mut path = vec![DFS_PATH::default(); path_len + 1];
        path[0].at_no = 0;
        for index in 1..=path_len {
            let next = if index == path_len { 0 } else { index + 1 };
            let bond_type = if index % 2 == 1 {
                BOND_SINGLE as u8
            } else {
                BOND_DOUBLE as u8
            };
            let (bond_pos, _) = add_bond_with_positions(&mut atoms, index, next, bond_type);
            path[index] = DFS_PATH {
                at_no: index as AT_RANK,
                bond_type,
                bond_pos: bond_pos as i8,
            };
        }

        atoms[first_endpoint].el_number = 8;
        atoms[first_endpoint].chem_bonds_valence = 2;
        atoms[first_endpoint].endpoint = endpoint_group;
        atoms[second_endpoint].el_number = 8;
        atoms[second_endpoint].chem_bonds_valence = 1;
        atoms[second_endpoint].num_H = 1;
        atoms[second_endpoint].endpoint = endpoint_group;

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let bns = check7_bns(&mut heap, atom_count, first_endpoint, second_endpoint);
        Check7Fixture {
            heap,
            atoms,
            path,
            start_neighbor2,
            start_neighbor_neighbor,
            first_endpoint,
            second_endpoint,
            bns,
            bd: BalancedNetworkData::default(),
        }
    }

    #[test]
    fn source_port__ichiqueu__check7membtautring__line_951() {
        fn run(
            fixture: &mut Check7Fixture,
            path_len: i32,
            endpoints: &mut [T_ENDPOINT],
            max_endpoints: i32,
            bonds: &mut [T_BONDPOS],
            max_bonds: i32,
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            Check7MembTautRing(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                &fixture.path,
                path_len,
                -1,
                fixture.start_neighbor2,
                fixture.start_neighbor_neighbor,
                endpoints,
                max_endpoints,
                bonds,
                max_bonds,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                (path_len + 3).max(0),
                0,
            )
        }

        let mut accepted = check7_fixture(4, 9);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 6];
        let mut endpoint_count = 0;
        let mut bond_count = 0;
        assert_eq!(
            run(
                &mut accepted,
                4,
                &mut endpoints,
                2,
                &mut bonds,
                6,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 6));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (endpoint.nAtomNumber, endpoint.nGroupNumber, endpoint.nEquNumber))
                .collect::<Vec<_>>(),
            vec![(5, 9, 0), (6, 9, 0)]
        );
        assert!(endpoints.iter().all(|endpoint| endpoint.num == [0; 5]));
        assert!(endpoints.iter().all(|endpoint| endpoint.num_DA == [0; 6]));
        assert_eq!(
            bonds
                .iter()
                .map(|bond| (bond.nAtomNumber, bond.neighbor_index))
                .collect::<Vec<_>>(),
            vec![(1, 0), (1, 1), (2, 1), (3, 1), (4, 1), (0, 0)]
        );
        assert_eq!(
            run(
                &mut accepted,
                4,
                &mut endpoints,
                2,
                &mut bonds,
                6,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(0)
        );
        assert_eq!((endpoint_count, bond_count), (2, 6));

        let mut accepted_six = check7_fixture(6, 11);
        let mut endpoints_six = vec![T_ENDPOINT::default(); 2];
        let mut bonds_six = vec![T_BONDPOS::default(); 8];
        let mut endpoint_count_six = 0;
        let mut bond_count_six = 0;
        assert_eq!(
            run(
                &mut accepted_six,
                6,
                &mut endpoints_six,
                2,
                &mut bonds_six,
                8,
                &mut endpoint_count_six,
                &mut bond_count_six,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count_six, bond_count_six), (2, 8));

        let mut invalid_call = check7_fixture(4, 9);
        for invalid_length in [7, 5] {
            let mut invalid_endpoints = vec![T_ENDPOINT::default(); 2];
            let mut invalid_bonds = vec![T_BONDPOS::default(); 6];
            let mut invalid_endpoint_count = 0;
            let mut invalid_bond_count = 0;
            assert_eq!(
                run(
                    &mut invalid_call,
                    invalid_length,
                    &mut invalid_endpoints,
                    2,
                    &mut invalid_bonds,
                    6,
                    &mut invalid_endpoint_count,
                    &mut invalid_bond_count,
                ),
                Ok(-1)
            );
            assert_eq!((invalid_endpoint_count, invalid_bond_count), (0, 0));
        }

        let mut invalid_endpoint = check7_fixture(4, 9);
        invalid_endpoint.heap.slice_mut(invalid_endpoint.atoms).unwrap()[invalid_endpoint.first_endpoint].radical =
            RADICAL_DOUBLET as i8;
        let mut rejected_endpoints = vec![T_ENDPOINT::default(); 2];
        let mut rejected_bonds = vec![T_BONDPOS::default(); 6];
        let mut rejected_endpoint_count = 0;
        let mut rejected_bond_count = 0;
        assert_eq!(
            run(
                &mut invalid_endpoint,
                4,
                &mut rejected_endpoints,
                2,
                &mut rejected_bonds,
                6,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((rejected_endpoint_count, rejected_bond_count), (0, 0));

        let mut invalid_bond = check7_fixture(4, 9);
        invalid_bond.path[2].bond_type = BOND_TRIPLE as u8;
        assert_eq!(
            run(
                &mut invalid_bond,
                4,
                &mut rejected_endpoints,
                2,
                &mut rejected_bonds,
                6,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut wrong_roles = check7_fixture(4, 0);
        {
            let atoms = wrong_roles.heap.slice_mut(wrong_roles.atoms).unwrap();
            atoms[wrong_roles.second_endpoint].num_H = 0;
            atoms[wrong_roles.second_endpoint].chem_bonds_valence = 2;
        }
        assert_eq!(
            run(
                &mut wrong_roles,
                4,
                &mut rejected_endpoints,
                2,
                &mut rejected_bonds,
                6,
                &mut rejected_endpoint_count,
                &mut rejected_bond_count,
            ),
            Ok(0)
        );

        let mut overflow = check7_fixture(4, 9);
        let mut overflow_endpoints = vec![T_ENDPOINT::default(); 2];
        let mut overflow_bonds = vec![T_BONDPOS::default(); 6];
        let mut overflow_endpoint_count = 0;
        let mut overflow_bond_count = 0;
        assert_eq!(
            run(
                &mut overflow,
                4,
                &mut overflow_endpoints,
                0,
                &mut overflow_bonds,
                0,
                &mut overflow_endpoint_count,
                &mut overflow_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((overflow_endpoint_count, overflow_bond_count), (0, 0));
        assert_eq!(overflow_endpoints[0].nAtomNumber, 5);
        assert_eq!(overflow_bonds[0].nAtomNumber, 1);

        let mut rejected_path = check7_fixture(4, 0);
        let mut path_endpoints = vec![T_ENDPOINT::default(); 2];
        let mut path_bonds = vec![T_BONDPOS::default(); 6];
        let mut path_endpoint_count = 0;
        let mut path_bond_count = 0;
        assert_eq!(
            run(
                &mut rejected_path,
                4,
                &mut path_endpoints,
                2,
                &mut path_bonds,
                6,
                &mut path_endpoint_count,
                &mut path_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((path_endpoint_count, path_bond_count), (2, 6));
        assert_eq!(path_endpoints[0].nAtomNumber, 5);
        assert_eq!(path_bonds[0].nAtomNumber, 1);
    }

    struct NGet14Fixture {
        heap: SourceHeap,
        atoms: SourceMutPointer<inp_ATOM>,
        positions: SourceMutPointer<AT_RANK>,
        path: Vec<DFS_PATH>,
        start_neighbor: i32,
        start_endpoint: i32,
        neighbor_endpoint: i32,
        num_atoms: i32,
        bns: BalancedNetworkStructure,
        bd: BalancedNetworkData,
    }

    fn nget14_fixture(path_len: usize) -> NGet14Fixture {
        let first_endpoint = path_len + 1;
        let second_endpoint = path_len + 2;
        let atom_count = path_len + 3;
        let mut atoms = vec![inp_ATOM::default(); atom_count];
        let (start_neighbor, _) = add_bond_with_positions(&mut atoms, 0, 1, BOND_SINGLE as u8);
        let (neighbor_endpoint, _) = add_bond_with_positions(&mut atoms, 1, first_endpoint, BOND_DOUBLE as u8);
        let (start_endpoint, _) = add_bond_with_positions(&mut atoms, 0, second_endpoint, BOND_SINGLE as u8);
        for index in 1..=path_len {
            let next = if index == path_len { 0 } else { index + 1 };
            let bond_type = if index % 2 == 1 {
                BOND_SINGLE as u8
            } else {
                BOND_DOUBLE as u8
            };
            add_bond_with_positions(&mut atoms, index, next, bond_type);
        }
        for atom in &mut atoms[..=path_len] {
            atom.el_number = 6;
            atom.chem_bonds_valence = atom.valence.wrapping_add(1);
        }
        atoms[first_endpoint].el_number = 8;
        atoms[first_endpoint].chem_bonds_valence = 2;
        atoms[first_endpoint].endpoint = 9;
        atoms[second_endpoint].el_number = 8;
        atoms[second_endpoint].chem_bonds_valence = 1;
        atoms[second_endpoint].num_H = 1;
        atoms[second_endpoint].endpoint = 9;

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let positions = heap.allocate_model_storage(vec![0; atom_count]).unwrap();
        let bns = check7_bns(&mut heap, atom_count, first_endpoint, second_endpoint);
        NGet14Fixture {
            heap,
            atoms,
            positions,
            path: vec![DFS_PATH::default(); 8],
            start_neighbor,
            start_endpoint,
            neighbor_endpoint,
            num_atoms: atom_count as i32,
            bns,
            bd: BalancedNetworkData::default(),
        }
    }

    #[test]
    fn source_port__ichiqueu__nget14tautin7membaltring__line_240() {
        fn run(
            fixture: &mut NGet14Fixture,
            max_path: i32,
            endpoints: &mut [T_ENDPOINT],
            bonds: &mut [T_BONDPOS],
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            nGet14TautIn7MembAltRing(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                0,
                fixture.start_neighbor,
                fixture.start_endpoint,
                fixture.neighbor_endpoint,
                fixture.positions,
                &mut fixture.path,
                max_path,
                endpoints,
                2,
                bonds,
                8,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                fixture.num_atoms,
                0,
            )
        }

        let mut fixture = nget14_fixture(6);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 8];
        let mut endpoint_count = 77;
        let mut bond_count = 88;
        assert_eq!(
            run(
                &mut fixture,
                8,
                &mut endpoints,
                &mut bonds,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 8));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (endpoint.nAtomNumber, endpoint.nGroupNumber))
                .collect::<Vec<_>>(),
            vec![(7, 9), (8, 9)]
        );
        assert_eq!(fixture.heap.slice(fixture.positions.as_const()).unwrap(), &[0; 9]);
        assert_eq!(fixture.path[0].at_no, 0);
        assert_eq!(fixture.path[1].at_no, 1);

        let mut too_short = nget14_fixture(6);
        too_short.heap.slice_mut(too_short.positions).unwrap().fill(13);
        let mut too_short_endpoint_count = 77;
        let mut too_short_bond_count = 88;
        assert_eq!(
            run(
                &mut too_short,
                7,
                &mut endpoints,
                &mut bonds,
                &mut too_short_endpoint_count,
                &mut too_short_bond_count,
            ),
            Ok(-1)
        );
        assert_eq!((too_short_endpoint_count, too_short_bond_count), (0, 0));
        assert_eq!(too_short.heap.slice(too_short.positions.as_const()).unwrap(), &[13; 9]);

        let mut no_ring = nget14_fixture(6);
        {
            let atoms = no_ring.heap.slice_mut(no_ring.atoms).unwrap();
            atoms[2].el_number = 0;
            atoms[2].chem_bonds_valence = 0;
        }
        let mut no_ring_endpoint_count = 77;
        let mut no_ring_bond_count = 88;
        assert_eq!(
            run(
                &mut no_ring,
                8,
                &mut endpoints,
                &mut bonds,
                &mut no_ring_endpoint_count,
                &mut no_ring_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((no_ring_endpoint_count, no_ring_bond_count), (0, 0));
        assert_eq!(no_ring.heap.slice(no_ring.positions.as_const()).unwrap(), &[0; 9]);

        let mut missing_path = nget14_fixture(6);
        missing_path.path.clear();
        let mut missing_endpoint_count = 77;
        let mut missing_bond_count = 88;
        assert_eq!(
            run(
                &mut missing_path,
                8,
                &mut endpoints,
                &mut bonds,
                &mut missing_endpoint_count,
                &mut missing_bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((missing_endpoint_count, missing_bond_count), (0, 0));
    }

    #[test]
    fn source_port__ichiqueu__nget14tautin5membaltring__line_283() {
        fn run(
            fixture: &mut NGet14Fixture,
            max_path: i32,
            endpoints: &mut [T_ENDPOINT],
            bonds: &mut [T_BONDPOS],
            endpoint_count: &mut i32,
            bond_count: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            nGet14TautIn5MembAltRing(
                &mut fixture.heap,
                &mut CANON_GLOBALS::default(),
                fixture.atoms,
                0,
                fixture.start_neighbor,
                fixture.start_endpoint,
                fixture.neighbor_endpoint,
                fixture.positions,
                &mut fixture.path,
                max_path,
                endpoints,
                2,
                bonds,
                6,
                endpoint_count,
                bond_count,
                &mut fixture.bns,
                &mut fixture.bd,
                fixture.num_atoms,
                0,
            )
        }

        let mut fixture = nget14_fixture(4);
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 6];
        let mut endpoint_count = 77;
        let mut bond_count = 88;
        assert_eq!(
            run(
                &mut fixture,
                6,
                &mut endpoints,
                &mut bonds,
                &mut endpoint_count,
                &mut bond_count,
            ),
            Ok(1)
        );
        assert_eq!((endpoint_count, bond_count), (2, 6));
        assert_eq!(
            endpoints
                .iter()
                .map(|endpoint| (endpoint.nAtomNumber, endpoint.nGroupNumber))
                .collect::<Vec<_>>(),
            vec![(5, 9), (6, 9)]
        );
        assert_eq!(fixture.heap.slice(fixture.positions.as_const()).unwrap(), &[0; 7]);
        assert_eq!(fixture.path[0].at_no, 0);
        assert_eq!(fixture.path[1].at_no, 1);

        let mut too_short = nget14_fixture(4);
        too_short.heap.slice_mut(too_short.positions).unwrap().fill(13);
        let mut too_short_endpoint_count = 77;
        let mut too_short_bond_count = 88;
        assert_eq!(
            run(
                &mut too_short,
                5,
                &mut endpoints,
                &mut bonds,
                &mut too_short_endpoint_count,
                &mut too_short_bond_count,
            ),
            Ok(-1)
        );
        assert_eq!((too_short_endpoint_count, too_short_bond_count), (0, 0));
        assert_eq!(too_short.heap.slice(too_short.positions.as_const()).unwrap(), &[13; 7]);

        let mut no_ring = nget14_fixture(4);
        {
            let atoms = no_ring.heap.slice_mut(no_ring.atoms).unwrap();
            atoms[2].el_number = 0;
            atoms[2].chem_bonds_valence = 0;
        }
        let mut no_ring_endpoint_count = 77;
        let mut no_ring_bond_count = 88;
        assert_eq!(
            run(
                &mut no_ring,
                6,
                &mut endpoints,
                &mut bonds,
                &mut no_ring_endpoint_count,
                &mut no_ring_bond_count,
            ),
            Ok(0)
        );
        assert_eq!((no_ring_endpoint_count, no_ring_bond_count), (0, 0));
        assert_eq!(no_ring.heap.slice(no_ring.positions.as_const()).unwrap(), &[0; 7]);

        let mut missing_path = nget14_fixture(4);
        missing_path.path.clear();
        let mut missing_endpoint_count = 77;
        let mut missing_bond_count = 88;
        assert_eq!(
            run(
                &mut missing_path,
                6,
                &mut endpoints,
                &mut bonds,
                &mut missing_endpoint_count,
                &mut missing_bond_count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((missing_endpoint_count, missing_bond_count), (0, 0));
    }

    #[test]
    fn source_port__ichiqueu__dfs_findtautinaring__line_458() {
        let mut atom_values = vec![inp_ATOM::default(); 4];
        add_bond(&mut atom_values, 0, 1, 0xf1);
        add_bond(&mut atom_values, 1, 2, 0xf2);
        add_bond(&mut atom_values, 2, 3, 0xf3);
        add_bond(&mut atom_values, 3, 0, 0xf4);

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(atom_values).unwrap();
        let mut canon = CANON_GLOBALS::default();
        let mut bns = BalancedNetworkStructure::default();
        let mut bd = BalancedNetworkData::default();
        let positions = heap.allocate_model_storage(vec![0; 4]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 3];
        let mut endpoint_count = 0;
        let mut bond_count = 0;
        let mut calls = Vec::new();
        let mut ring = |args: DfsRingCallbackArgs<'_>| {
            assert_eq!(args.nLenDfsPath, 3);
            assert_eq!(args.nStartAtomNeighbor, -1);
            assert_eq!(args.nStartAtomNeighbor2, -1);
            assert_eq!(args.nStartAtomNeighborNeighbor, -1);
            assert_eq!(args.nMaxNumEndPoint, 2);
            assert_eq!(args.nMaxNumBondPos, 3);
            assert_eq!(args.num_atoms, 4);
            assert_eq!(args.heap.slice(args.atom.as_const()).unwrap().len(), 4);
            assert_eq!(args.EndPoint.len(), 2);
            assert_eq!(args.BondPos.len(), 3);
            assert_eq!(args.pBNS.as_ref().unwrap().num_atoms, 0);
            assert_eq!(args.pBD.as_ref().unwrap().QSize, 0);
            args.pCG.m_num_bit += 1;
            *args.pnNumEndPoint += 1;
            *args.pnNumBondPos += 2;
            args.EndPoint[0].nAtomNumber = 17;
            args.BondPos[0].nAtomNumber = 19;
            calls.push(
                args.DfsPath[..=usize::try_from(args.nLenDfsPath).unwrap()]
                    .iter()
                    .map(|entry| (entry.at_no, entry.bond_type))
                    .collect::<Vec<_>>(),
            );
            Ok(2)
        };
        let mut center = |_heap: &SourceHeap, _atoms: SourceMutPointer<inp_ATOM>, _iat: i32| Ok(1);
        assert_eq!(
            DFS_FindTautInARing(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut ring,
                &mut center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                Some(&mut bns),
                Some(&mut bd),
                4,
            ),
            Ok(4)
        );
        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0], vec![(0, 1), (1, 2), (2, 3), (3, 4)]);
        assert_eq!(calls[1], vec![(0, 4), (3, 3), (2, 2), (1, 1)]);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 4]);
        assert_eq!(canon.m_num_bit, 2);
        assert_eq!(endpoint_count, 2);
        assert_eq!(bond_count, 4);
        assert_eq!(endpoints[0].nAtomNumber, 17);
        assert_eq!(bonds[0].nAtomNumber, 19);

        let positions = heap.allocate_model_storage(vec![0; 4]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut calls = 0;
        let mut negative_ring = |_args: DfsRingCallbackArgs<'_>| {
            calls += 1;
            Ok(-71)
        };
        assert_eq!(
            DFS_FindTautInARing(
                &mut heap,
                &mut canon,
                atoms,
                0,
                0,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut negative_ring,
                &mut center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                4,
            ),
            Ok(-71)
        );
        assert_eq!(calls, 1);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 4]);

        for (neighbor2, neighbor_neighbor, center_result) in [(1, -1, 1), (-1, 1, 1), (-1, -1, 0)] {
            let positions = heap.allocate_model_storage(vec![0; 4]).unwrap();
            let mut paths = vec![DFS_PATH::default(); 8];
            let mut ring_calls = 0;
            let mut no_ring = |_args: DfsRingCallbackArgs<'_>| {
                ring_calls += 1;
                Ok(1)
            };
            let mut selective_center =
                move |_heap: &SourceHeap, _atoms: SourceMutPointer<inp_ATOM>, _iat: i32| Ok(center_result);
            assert_eq!(
                DFS_FindTautInARing(
                    &mut heap,
                    &mut canon,
                    atoms,
                    0,
                    if neighbor_neighbor >= 0 { 0 } else { -1 },
                    neighbor2,
                    neighbor_neighbor,
                    4,
                    positions,
                    &mut paths,
                    &mut no_ring,
                    &mut selective_center,
                    &mut endpoints,
                    2,
                    &mut bonds,
                    3,
                    &mut endpoint_count,
                    &mut bond_count,
                    None,
                    None,
                    4,
                ),
                Ok(0)
            );
            assert_eq!(ring_calls, 0);
            assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 4]);
        }

        let positions = heap.allocate_model_storage(vec![0; 4]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut no_ring = |_args: DfsRingCallbackArgs<'_>| Ok(1);
        assert_eq!(
            DFS_FindTautInARing(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                3,
                positions,
                &mut paths,
                &mut no_ring,
                &mut center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                4,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 4]);

        let mut too_short_paths = Vec::new();
        let too_short_positions = heap.allocate_model_storage(vec![0; 4]).unwrap();
        assert_eq!(
            DFS_FindTautInARing(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                4,
                too_short_positions,
                &mut too_short_paths,
                &mut no_ring,
                &mut center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                4,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiqueu__dfs_findtautaltpath__line_616() {
        fn chain_atoms(outside_ring: bool) -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 5];
            add_bond(&mut atoms, 0, 1, 0xf1);
            add_bond(&mut atoms, 1, 2, 0xf2);
            add_bond(&mut atoms, 2, 3, 0xf3);
            add_bond(&mut atoms, 3, 4, 0xf4);
            atoms[0].nNumAtInRingSystem = if outside_ring { 1 } else { 2 };
            atoms[4].nNumAtInRingSystem = if outside_ring { 1 } else { 2 };
            atoms
        }

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(chain_atoms(true)).unwrap();
        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut canon = CANON_GLOBALS::default();
        let mut bns = BalancedNetworkStructure::default();
        let mut bd = BalancedNetworkData::default();
        let mut endpoints = vec![T_ENDPOINT::default(); 2];
        let mut bonds = vec![T_BONDPOS::default(); 3];
        let mut endpoint_count = 0;
        let mut bond_count = 0;
        let mut path_calls = Vec::new();
        let mut center_calls = Vec::new();
        let mut check_path = |args: DfsPathCallbackArgs<'_>| {
            assert_eq!(args.nLenDfsPath, 3);
            assert_eq!(args.jNxtNeigh, 1);
            assert_eq!(args.nStartAtomNeighbor, -1);
            assert_eq!(args.nStartAtomNeighbor2, -1);
            assert_eq!(args.nStartAtomNeighborNeighbor, -1);
            assert_eq!(args.nMaxNumEndPoint, 2);
            assert_eq!(args.nMaxNumBondPos, 3);
            assert_eq!(args.num_atoms, 5);
            assert_eq!(args.heap.slice(args.atom.as_const()).unwrap().len(), 5);
            assert_eq!(args.pBNS.as_ref().unwrap().num_atoms, 0);
            assert_eq!(args.pBD.as_ref().unwrap().QSize, 0);
            args.pCG.m_num_bit = args.pCG.m_num_bit.wrapping_add(1);
            *args.pnNumEndPoint = args.pnNumEndPoint.wrapping_add(1);
            *args.pnNumBondPos = args.pnNumBondPos.wrapping_add(2);
            args.EndPoint[0].nAtomNumber = 17;
            args.BondPos[0].nAtomNumber = 19;
            path_calls.push(
                args.DfsPath[..=usize::try_from(args.nLenDfsPath).unwrap()]
                    .iter()
                    .map(|entry| (entry.at_no, entry.bond_type, entry.bond_pos))
                    .collect::<Vec<_>>(),
            );
            Ok(2)
        };
        let mut check_center = |args: DfsPathCenterCallbackArgs<'_>| {
            let current = usize::from(args.DfsPath[usize::try_from(args.nLenDfsPath).unwrap()].at_no);
            let next = usize::from(
                args.heap.slice(args.atom.as_const()).unwrap()[current].neighbor
                    [usize::try_from(args.jNxtNeigh).unwrap()],
            );
            assert_eq!(args.num_atoms, 5);
            assert!(args.pBNS.is_some());
            assert!(args.pBD.is_some());
            center_calls.push((args.nLenDfsPath, current, next));
            Ok(1)
        };
        assert_eq!(
            DFS_FindTautAltPath(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut check_path,
                &mut check_center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                Some(&mut bns),
                Some(&mut bd),
                5,
            ),
            Ok(2)
        );
        assert_eq!(path_calls, vec![vec![(0, 1, 0), (1, 2, 1), (2, 3, 1), (3, 4, 1)]]);
        assert_eq!(center_calls, vec![(0, 0, 1), (1, 1, 2), (2, 2, 3)]);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 5]);
        assert_eq!(canon.m_num_bit, 1);
        assert_eq!((endpoint_count, bond_count), (1, 2));
        assert_eq!(endpoints[0].nAtomNumber, 17);
        assert_eq!(bonds[0].nAtomNumber, 19);

        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut negative_calls = 0;
        let mut negative_path = |_args: DfsPathCallbackArgs<'_>| {
            negative_calls += 1;
            Ok(-71)
        };
        let mut all_centers = |_args: DfsPathCenterCallbackArgs<'_>| Ok(1);
        assert_eq!(
            DFS_FindTautAltPath(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut negative_path,
                &mut all_centers,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                5,
            ),
            Ok(-71)
        );
        assert_eq!(negative_calls, 1);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 5]);

        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut zero_path_calls = 0;
        let mut terminal_center_calls = Vec::new();
        let mut zero_path = |_args: DfsPathCallbackArgs<'_>| {
            zero_path_calls += 1;
            Ok(0)
        };
        let mut terminal_center = |args: DfsPathCenterCallbackArgs<'_>| {
            terminal_center_calls.push(args.nLenDfsPath);
            Ok(1)
        };
        assert_eq!(
            DFS_FindTautAltPath(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut zero_path,
                &mut terminal_center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                5,
            ),
            Ok(0)
        );
        assert_eq!(zero_path_calls, 1);
        assert_eq!(terminal_center_calls, vec![0, 1, 2, 3]);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 5]);

        let ring_atoms = heap.allocate_model_storage(chain_atoms(false)).unwrap();
        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut gated_path_calls = 0;
        let mut gated_center_calls = 0;
        let mut gated_path = |_args: DfsPathCallbackArgs<'_>| {
            gated_path_calls += 1;
            Ok(1)
        };
        let mut gated_center = |_args: DfsPathCenterCallbackArgs<'_>| {
            gated_center_calls += 1;
            Ok(1)
        };
        assert_eq!(
            DFS_FindTautAltPath(
                &mut heap,
                &mut canon,
                ring_atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut gated_path,
                &mut gated_center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                5,
            ),
            Ok(0)
        );
        assert_eq!(gated_path_calls, 0);
        assert_eq!(gated_center_calls, 4);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 5]);

        for (first_neighbor, forbidden_from_start, forbidden_from_neighbor, expected_calls) in
            [(-1, 0, -1, 0), (0, -1, 1, 0), (0, -1, -1, 1)]
        {
            let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
            let mut paths = vec![DFS_PATH::default(); 8];
            let mut calls = 0;
            let mut path = |_args: DfsPathCallbackArgs<'_>| {
                calls += 1;
                Ok(1)
            };
            let mut center = |_args: DfsPathCenterCallbackArgs<'_>| Ok(1);
            assert_eq!(
                DFS_FindTautAltPath(
                    &mut heap,
                    &mut canon,
                    atoms,
                    0,
                    first_neighbor,
                    forbidden_from_start,
                    forbidden_from_neighbor,
                    4,
                    positions,
                    &mut paths,
                    &mut path,
                    &mut center,
                    &mut endpoints,
                    2,
                    &mut bonds,
                    3,
                    &mut endpoint_count,
                    &mut bond_count,
                    None,
                    None,
                    5,
                ),
                Ok(expected_calls)
            );
            assert_eq!(calls, expected_calls);
            assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 5]);
        }

        let mut branched = vec![inp_ATOM::default(); 6];
        for (first, second, bond_type) in [(0, 1, 1), (1, 2, 2), (0, 4, 1), (4, 2, 2), (2, 3, 3), (3, 5, 4)] {
            add_bond(&mut branched, first, second, bond_type);
        }
        branched[0].nNumAtInRingSystem = 1;
        let branched_atoms = heap.allocate_model_storage(branched).unwrap();
        let positions = heap.allocate_model_storage(vec![0; 6]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut repeated_paths = Vec::new();
        let mut path = |args: DfsPathCallbackArgs<'_>| {
            repeated_paths.push(
                args.DfsPath[..=usize::try_from(args.nLenDfsPath).unwrap()]
                    .iter()
                    .map(|entry| entry.at_no)
                    .collect::<Vec<_>>(),
            );
            Ok(1)
        };
        let mut center = |_args: DfsPathCenterCallbackArgs<'_>| Ok(1);
        assert_eq!(
            DFS_FindTautAltPath(
                &mut heap,
                &mut canon,
                branched_atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut path,
                &mut center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                6,
            ),
            Ok(2)
        );
        assert_eq!(repeated_paths, vec![vec![0, 1, 2, 3], vec![0, 4, 2, 3]]);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 6]);

        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let mut no_paths = Vec::new();
        let mut path = |_args: DfsPathCallbackArgs<'_>| Ok(1);
        let mut center = |_args: DfsPathCenterCallbackArgs<'_>| Ok(1);
        assert_eq!(
            DFS_FindTautAltPath(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut no_paths,
                &mut path,
                &mut center,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                5,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let positions = heap.allocate_model_storage(vec![0; 5]).unwrap();
        let mut paths = vec![DFS_PATH::default(); 8];
        let mut no_centers = |_args: DfsPathCenterCallbackArgs<'_>| Ok(0);
        assert_eq!(
            DFS_FindTautAltPath(
                &mut heap,
                &mut canon,
                atoms,
                0,
                -1,
                -1,
                -1,
                4,
                positions,
                &mut paths,
                &mut path,
                &mut no_centers,
                &mut endpoints,
                2,
                &mut bonds,
                3,
                &mut endpoint_count,
                &mut bond_count,
                None,
                None,
                5,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[0; 5]);
    }

    #[test]
    fn source_port__ichiqueu__are_alt_bonds__line_778() {
        assert_eq!(are_alt_bonds(&[], -1), Ok(0));
        assert_eq!(are_alt_bonds(&[], 0), Ok(0));
        assert_eq!(are_alt_bonds(&[BOND_SINGLE as u8], 1), Ok(0));
        assert_eq!(are_alt_bonds(&[BOND_TRIPLE as u8, 1], 2), Ok(0));
        assert_eq!(are_alt_bonds(&[BOND_ALT_13 as u8, 1], 2), Ok(0));

        assert_eq!(are_alt_bonds(&[1, 2], 2), Ok(BOND_DOUBLE as i32));
        assert_eq!(are_alt_bonds(&[2, 1], 2), Ok(BOND_SINGLE as i32));
        assert_eq!(are_alt_bonds(&[1, 2, 1], 3), Ok(BOND_SINGLE as i32));
        assert_eq!(are_alt_bonds(&[2, 1, 2], 3), Ok(BOND_DOUBLE as i32));
        assert_eq!(are_alt_bonds(&[1, 1], 2), Ok(0));
        assert_eq!(are_alt_bonds(&[2, 2], 2), Ok(0));

        assert_eq!(are_alt_bonds(&[4, 4], 2), Ok(BOND_ALTERN as i32));
        assert_eq!(are_alt_bonds(&[9, 8], 2), Ok(BOND_TAUTOM as i32));
        assert_eq!(are_alt_bonds(&[8, 4], 2), Ok(BOND_TAUTOM as i32));
        assert_eq!(are_alt_bonds(&[4, 1], 2), Ok(BOND_SINGLE as i32));
        assert_eq!(are_alt_bonds(&[4, 2], 2), Ok(BOND_DOUBLE as i32));
        assert_eq!(are_alt_bonds(&[0, 4], 2), Ok(BOND_ALTERN as i32));
        assert_eq!(are_alt_bonds(&[0, 8], 2), Ok(BOND_TAUTOM as i32));
        assert_eq!(are_alt_bonds(&[0, 0], 2), Ok(0));
        assert_eq!(are_alt_bonds(&[1, 4, 1], 3), Ok(BOND_SINGLE as i32));
        assert_eq!(are_alt_bonds(&[2, 9, 2], 3), Ok(BOND_DOUBLE as i32));

        assert_eq!(
            are_alt_bonds(&[BOND_SINGLE as u8], 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiqueu__addbondspos__line_846() {
        let mut atoms = vec![inp_ATOM::default(); 3];
        add_bond(&mut atoms, 0, 1, 1);
        add_bond(&mut atoms, 1, 2, 2);

        let mut temporary = vec![T_BONDPOS::default(); 4];
        temporary[0].nAtomNumber = 0;
        temporary[0].neighbor_index = 0;
        temporary[1].nAtomNumber = 77;
        temporary[1].neighbor_index = 78;
        temporary[2].nAtomNumber = 1;
        temporary[2].neighbor_index = 1;
        temporary[3].nAtomNumber = 79;
        temporary[3].neighbor_index = 80;
        let mut output = vec![T_BONDPOS::default(); 4];
        assert_eq!(AddBondsPos(&atoms, &mut temporary, 4, &mut output, 3, 0), Ok(2));
        assert_eq!((temporary[1].nAtomNumber, temporary[1].neighbor_index), (1, 0));
        assert_eq!((temporary[3].nAtomNumber, temporary[3].neighbor_index), (2, 0));
        assert_eq!((output[0].nAtomNumber, output[0].neighbor_index), (0, 0));
        assert_eq!((output[1].nAtomNumber, output[1].neighbor_index), (1, 1));

        let mut reverse_duplicate = vec![T_BONDPOS::default(); 2];
        reverse_duplicate[0].nAtomNumber = 1;
        reverse_duplicate[0].neighbor_index = 0;
        let preserved = output.clone();
        assert_eq!(
            AddBondsPos(&atoms, &mut reverse_duplicate, 2, &mut output, 3, 2,),
            Ok(2)
        );
        assert_eq!(output, preserved);

        let mut equality = vec![T_BONDPOS::default(); 2];
        equality[0].nAtomNumber = 1;
        equality[0].neighbor_index = 1;
        let mut equality_output = vec![T_BONDPOS::default(); 2];
        equality_output[0].nAtomNumber = 9;
        assert_eq!(AddBondsPos(&atoms, &mut equality, 2, &mut equality_output, 1, 1), Ok(2));
        assert_eq!(equality_output[1], equality[0]);

        let mut overflow = vec![T_BONDPOS::default(); 2];
        overflow[0].nAtomNumber = 1;
        overflow[0].neighbor_index = 1;
        let mut full = vec![T_BONDPOS::default(); 3];
        full[0].nAtomNumber = 7;
        full[1].nAtomNumber = 8;
        let before = full.clone();
        assert_eq!(AddBondsPos(&atoms, &mut overflow, 2, &mut full, 1, 2), Ok(-1));
        assert_eq!(full, before);
        assert_eq!((overflow[1].nAtomNumber, overflow[1].neighbor_index), (2, 0));

        let mut odd = vec![T_BONDPOS::default(); 1];
        assert_eq!(
            AddBondsPos(&atoms, &mut odd, 1, &mut output, 3, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            AddBondsPos(&atoms, &mut temporary, -1, &mut output, 3, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    fn endpoint(atom_number: u16, group: u16, marker: u16) -> T_ENDPOINT {
        T_ENDPOINT {
            nAtomNumber: atom_number,
            nGroupNumber: group,
            num: [marker; 5],
            ..T_ENDPOINT::default()
        }
    }

    #[test]
    fn source_port__ichiqueu__addendpoints__line_897() {
        let temporary = vec![endpoint(2, 20, 200), endpoint(3, 30, 300)];
        let mut output = vec![T_ENDPOINT::default(); 4];
        output[0] = endpoint(1, 10, 100);
        assert_eq!(AddEndPoints(&temporary, 2, &mut output, 3, 1), Ok(3));
        assert_eq!(output[1], temporary[0]);
        assert_eq!(output[2], temporary[1]);

        let duplicates = vec![endpoint(3, 99, 999), endpoint(1, 88, 888)];
        let preserved = output.clone();
        assert_eq!(AddEndPoints(&duplicates, 2, &mut output, 3, 3), Ok(3));
        assert_eq!(output, preserved);

        let equality = vec![endpoint(5, 50, 500)];
        let mut equality_output = vec![T_ENDPOINT::default(); 2];
        equality_output[0] = endpoint(4, 40, 400);
        assert_eq!(AddEndPoints(&equality, 1, &mut equality_output, 1, 1), Ok(2));
        assert_eq!(equality_output[1], equality[0]);

        let partial = vec![endpoint(6, 60, 600), endpoint(7, 70, 700)];
        let mut partial_output = vec![T_ENDPOINT::default(); 3];
        partial_output[0] = endpoint(4, 40, 400);
        assert_eq!(AddEndPoints(&partial, 2, &mut partial_output, 1, 1), Ok(-1));
        assert_eq!(partial_output[1], partial[0]);
        assert_eq!(partial_output[2], T_ENDPOINT::default());

        assert_eq!(
            AddEndPoints(&temporary, -1, &mut output, 3, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            AddEndPoints(&temporary[..1], 2, &mut output, 3, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }
}
