use crate::source::base::ichi_bns::mark_alt_bonds_and_taut_groups;
use crate::source::base::mol2atom::CreateInpAtomData;
use crate::source::base::runichi3::UnMarkRingSystemsInp;
use crate::source::base::strutil::{
    DisconnectInpAtBond, ExtractConnectedComponent, MarkDisconnectedComponents, add_DT_to_num_H,
    remove_terminal_HDT,
};
use crate::source::base::util::{inchi_calloc, inchi_free, is_el_a_metal, is_in_the_list};
use crate::source_types::local_ichinorm::{
    BIT_UNDERIV, CFLAG_MARK_BLOCK, CFLAG_MARK_BRANCH, DERIV_AMINE_tN, DERIV_AT, DERIV_BRIDGE_NH,
    DERIV_BRIDGE_O, DERIV_DANSYL, DERIV_DUPLIC, DERIV_NOT, DERIV_REPL_N_WITH_O,
    DERIV_REPL_N_WITH_OH, DERIV_RING_DMOX_DEOX, DERIV_RING_DMOX_DEOX_N, DERIV_RING_DMOX_DEOX_O,
    DERIV_RING_NH_OUTSIDE_PRECURSOR, DERIV_RING_O_OUTSIDE_PRECURSOR, DERIV_RING_OUTSIDE_PRECURSOR,
    DERIV_RING2_OUTSIDE_PRECUR, DERIV_RING2_PPRDN_OUTSIDE_PRECUR,
    DERIV_RING2_PRRLDD_OUTSIDE_PRECUR, DERIV_RO_COX, DERIV_UNEXPADABLE, DERIV_UNMARK,
    DERIV_X_OXIME, MAX_AT_DERIV, MIN_AT_LEFT_DERIV, NOT_AT_DERIV, R2C_AT, R2C_ATPAIR, R2C_EMPTY,
    tagDerivBit_DERIV_BIT_Acetate as DERIV_BIT_ACETATE,
    tagDerivBit_DERIV_BIT_BenzOX as DERIV_BIT_BENZOX,
    tagDerivBit_DERIV_BIT_Benzoate as DERIV_BIT_BENZOATE,
    tagDerivBit_DERIV_BIT_DEOX as DERIV_BIT_DEOX, tagDerivBit_DERIV_BIT_DMOX as DERIV_BIT_DMOX,
    tagDerivBit_DERIV_BIT_Dansyl as DERIV_BIT_DANSYL, tagDerivBit_DERIV_BIT_EtOX as DERIV_BIT_ETOX,
    tagDerivBit_DERIV_BIT_HFB as DERIV_BIT_HFB, tagDerivBit_DERIV_BIT_MOX as DERIV_BIT_MOX,
    tagDerivBit_DERIV_BIT_PFP as DERIV_BIT_PFP,
    tagDerivBit_DERIV_BIT_Piperidine as DERIV_BIT_PIPERIDINE,
    tagDerivBit_DERIV_BIT_Pyrrolidide as DERIV_BIT_PYRROLIDIDE,
    tagDerivBit_DERIV_BIT_TBDMS as DERIV_BIT_TBDMS, tagDerivBit_DERIV_BIT_TFA as DERIV_BIT_TFA,
    tagDerivBit_DERIV_BIT_TMS as DERIV_BIT_TMS, tagDerivBit_DERIV_BIT_Unknown as DERIV_BIT_UNKNOWN,
};
use crate::source_types::{
    BOND_DOUBLE, BOND_SINGLE, BOND_TRIPLE, BOND_TYPE_ALTERN, CANON_GLOBALS, CT_OUT_OF_RAM,
    CT_OVERFLOW, INCHI_CLOCK, INP_ATOM_DATA, MAX_SDF_VALUE, ORIG_ATOM_DATA, SourceConstPointer,
    SourceHeap, SourceHeapError, SourceMutPointer, clock_t, copy_inp_atom_gcc_lp64_byte_prefix,
    inp_ATOM,
};

const EL_NUMBER_C: u8 = 6;
const EL_NUMBER_B: u8 = 5;
const EL_NUMBER_F: u8 = 9;
const EL_NUMBER_N: u8 = 7;
const EL_NUMBER_O: u8 = 8;
const EL_NUMBER_P: u8 = 15;
const EL_NUMBER_SI: u8 = 14;
const EL_NUMBER_S: u8 = 16;
const EL_NUMBER_H: u8 = 1;

#[allow(non_snake_case)]
pub(crate) fn mark_arom_bonds(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: &mut CANON_GLOBALS,
    at: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:784 mark_arom_bonds
    // INCHI✔️❌: int mark_arom_bonds( struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     INCHI_MODE bTautFlags = 0, bTautFlagsDone = 0;
    // INCHI✔️❌:     inp_ATOM *at_fixed_bonds_out = NULL;
    // INCHI✔️❌:     T_GROUP_INFO *t_group_info = NULL;
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = mark_alt_bonds_and_taut_groups( ic, pCG, at, at_fixed_bonds_out, num_atoms,
    // INCHI✔️❌:                                           NULL,
    // INCHI✔️❌:                                           t_group_info, &bTautFlags, &bTautFlagsDone, 0, NULL );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: mark_arom_bonds

    let taut_flags = heap.allocate_model_storage(vec![0_u64])?;
    let taut_flags_done = match heap.allocate_model_storage(vec![0_u64]) {
        Ok(pointer) => pointer,
        Err(error) => {
            heap.free(taut_flags)?;
            return Err(error);
        }
    };
    let result = mark_alt_bonds_and_taut_groups(
        heap,
        ic,
        pCG,
        at,
        SourceMutPointer::null(),
        num_atoms,
        SourceMutPointer::null(),
        SourceMutPointer::null(),
        taut_flags,
        taut_flags_done,
        0,
        SourceMutPointer::null(),
        clock_result,
    );
    let free_flags = heap.free(taut_flags);
    let free_flags_done = heap.free(taut_flags_done);
    free_flags?;
    free_flags_done?;
    result
}

pub(crate) fn cmp_r2c_atpair(pair1: &R2C_ATPAIR, pair2: &R2C_ATPAIR) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:800 cmp_r2c_atpair
    // INCHI✔️❌: int cmp_r2c_atpair( const void *p1, const void *p2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     const R2C_ATPAIR *ap1 = (const R2C_ATPAIR *) p1;
    // INCHI✔️❌:     const R2C_ATPAIR *ap2 = (const R2C_ATPAIR *) p2;
    // INCHI✔️❌:     int diff = (int) ap1->at[0] - (int) ap2->at[0];
    // INCHI✔️❌:     if (!diff)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         diff = (int) ap1->at[1] - (int) ap2->at[1];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return diff;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: cmp_r2c_atpair

    let difference = i32::from(pair1.at[0]) - i32::from(pair2.at[0]);
    if difference != 0 {
        difference
    } else {
        i32::from(pair1.at[1]) - i32::from(pair2.at[1])
    }
}

pub(crate) fn has_atom_pair_seq(
    atom_pairs: &[R2C_ATPAIR],
    num_atom_pairs: i32,
    atom1: u16,
    atom2: u16,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:815 has_atom_pair_seq
    // INCHI✔️❌: int has_atom_pair_seq( R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     R2C_ATPAIR ap1;
    // INCHI✔️❌:     int i1;
    // INCHI✔️❌:     int n = at1 > at2;
    // INCHI✔️❌:
    // INCHI✔️❌:     ap1.at[n] = at1;
    // INCHI✔️❌:     ap1.at[1 - n] = at2;
    // INCHI✔️❌:     for (i1 = 0; i1 < num_ap; i1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ap[i1].at[0] == ap1.at[0] &&
    // INCHI✔️❌:              ap[i1].at[1] == ap1.at[1])
    // INCHI✔️❌:             return i1 + 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0; /* not found */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: has_atom_pair_seq

    let position = usize::from(atom1 > atom2);
    let mut sought = [0_u16; 2];
    sought[position] = atom1;
    sought[1 - position] = atom2;
    let mut index = 0_i32;
    while index < num_atom_pairs {
        let pair = atom_pairs
            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if pair.at[0] == sought[0] && pair.at[1] == sought[1] {
            return Ok(index + 1);
        }
        index += 1;
    }
    Ok(0)
}

pub(crate) fn mark_atoms_ap(
    atoms: &mut [inp_ATOM],
    start: u16,
    atom_pairs: &[R2C_ATPAIR],
    num_atom_pairs: i32,
    mut num: i32,
    flags: u16,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:871 mark_atoms_ap
    // INCHI✔️❌: int mark_atoms_ap( inp_ATOM *at,
    // INCHI✔️❌:                    AT_NUMB start,
    // INCHI✔️❌:                    R2C_ATPAIR *ap,
    // INCHI✔️❌:                    int num_ap, int num,
    // INCHI✔️❌:                    AT_NUMB cFlags )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (!at[start].at_type)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i;
    // INCHI✔️❌:         AT_NUMB neigh;
    // INCHI✔️❌:         at[start].at_type = cFlags;
    // INCHI✔️❌:         num++;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[start].neighbor[i];
    // INCHI✔️❌:             if (has_atom_pair_seq( ap, num_ap, start, neigh ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             num = mark_atoms_ap( at, neigh, ap, num_ap, num, cFlags );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num; /* number of atoms traversed forward from at[start] */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: mark_atoms_ap

    let start_index = usize::from(start);
    let start_atom = atoms
        .get(start_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if start_atom.at_type != 0 {
        return Ok(num);
    }
    atoms[start_index].at_type = flags;
    num = num.wrapping_add(1);
    let valence = i32::from(atoms[start_index].valence);
    let mut neighbor_index = 0_i32;
    while neighbor_index < valence {
        let neighbor = *atoms[start_index]
            .neighbor
            .get(usize::try_from(neighbor_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if has_atom_pair_seq(atom_pairs, num_atom_pairs, start, neighbor)? == 0 {
            num = mark_atoms_ap(atoms, neighbor, atom_pairs, num_atom_pairs, num, flags)?;
        }
        neighbor_index += 1;
    }
    Ok(num)
}

pub(crate) fn mark_deriv_agents(
    atoms: &mut [inp_ATOM],
    derivatives: &[DERIV_AT],
    num_atoms: i32,
    atom_pairs: &[R2C_ATPAIR],
    num_cuts_to_check: i32,
    num_components: &mut u16,
    current_num_atoms: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5743 mark_deriv_agents
    // INCHI✔️❌: int mark_deriv_agents( inp_ATOM *at,
    // INCHI✔️❌:                        DERIV_AT *da,
    // INCHI✔️❌:                        int num_atoms,
    // INCHI✔️❌:                        R2C_ATPAIR *ap,
    // INCHI✔️❌:                        int num_cuts_to_check,
    // INCHI✔️❌:                        AT_NUMB *pnum_comp,
    // INCHI✔️❌:                        int *pcur_num_at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* mark components to be disconnected */
    // INCHI✔️❌:     int comp_num = 0;   /* number of components */
    // INCHI✔️❌:     int cur_num_at = 0; /* number of atoms left after disconnecting the derivatizing agent */
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     int i, j, k = -1, n;
    // INCHI✔️❌:     *pnum_comp = 0;
    // INCHI✔️❌:     *pcur_num_at = 0;
    // INCHI✔️❌:     UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:     for (i = 0; i < num_cuts_to_check; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         n = 0;
    // INCHI✔️❌:         for (j = 0; j < 2; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da[(int) ap[i].at[j]].typ[0])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k = j;
    // INCHI✔️❌:                 n++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (n != 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = -3;
    // INCHI✔️❌:             goto exit_r2c_num; /* wrong atom pair */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         n = ap[i].at[k]; /* marked atom */
    // INCHI✔️❌:         if (( da[n].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             n = ap[i].at[1 - k];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* at[n] belongs to the precursor */
    // INCHI✔️❌:         if (!at[n].at_type)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             comp_num++;
    // INCHI✔️❌:             cur_num_at = mark_atoms_ap( at, n, ap, num_cuts_to_check, cur_num_at, comp_num );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *pnum_comp = comp_num;
    // INCHI✔️❌:     *pcur_num_at = cur_num_at;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_r2c_num:
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: mark_deriv_agents

    *num_components = 0;
    *current_num_atoms = 0;
    UnMarkOtherIndicators(atoms, num_atoms)?;
    let mut component_count = 0_i32;
    let mut atom_count = 0_i32;
    let mut pair_index = 0_i32;
    while pair_index < num_cuts_to_check {
        let pair = atom_pairs
            .get(usize::try_from(pair_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut marked_position = 0_usize;
        let mut marked_count = 0_i32;
        for position in 0..2_usize {
            let derivative = derivatives
                .get(usize::from(pair.at[position]))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if derivative.typ[0] != 0 {
                marked_position = position;
                marked_count += 1;
            }
        }
        if marked_count != 1 {
            return Ok(-3);
        }
        let marked_atom = pair.at[marked_position];
        let marked_derivative = derivatives
            .get(usize::from(marked_atom))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let precursor =
            if i32::from(marked_derivative.typ[0]) & DERIV_RING_OUTSIDE_PRECURSOR as i32 != 0 {
                pair.at[1 - marked_position]
            } else {
                marked_atom
            };
        let precursor_atom = atoms
            .get(usize::from(precursor))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if precursor_atom.at_type == 0 {
            component_count = component_count.wrapping_add(1);
            atom_count = mark_atoms_ap(
                atoms,
                precursor,
                atom_pairs,
                num_cuts_to_check,
                atom_count,
                component_count as u16,
            )?;
        }
        pair_index += 1;
    }
    *num_components = component_count as u16;
    *current_num_atoms = atom_count;
    Ok(0)
}

pub(crate) fn replace_arom_bonds(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
    original_atoms: &[inp_ATOM],
    num_original_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5809 replace_arom_bonds
    // INCHI✔️❌: int replace_arom_bonds( inp_ATOM *at,
    // INCHI✔️❌:                         int num_atoms,
    // INCHI✔️❌:                         inp_ATOM *at2,
    // INCHI✔️❌:                         int num_atoms2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, num_err = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[i].bond_type[j] > BOND_TRIPLE)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* find pairs of atoms using orig. atom numbers */
    // INCHI✔️❌:                 int i1, i2;
    // INCHI✔️❌:                 char bSuccess = 0;
    // INCHI✔️❌:                 int neigh = at[i].neighbor[j];
    // INCHI✔️❌:                 AT_NUMB orig_no1 = at[i].orig_at_number;
    // INCHI✔️❌:                 AT_NUMB orig_no2 = at[neigh].orig_at_number;
    // INCHI✔️❌:                 for (i1 = 0; i1 < num_atoms2 && at2[i1].orig_at_number != orig_no1; i1++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 for (i2 = 0; i2 < num_atoms2 && at2[i2].orig_at_number != orig_no2; i2++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (i1 < num_atoms2 && i2 < num_atoms2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     AT_NUMB *p1 = is_in_the_list( at2[i1].neighbor, (AT_NUMB) i2, at[i1].valence );
    // INCHI✔️❌:                     AT_NUMB *pneigh = is_in_the_list( at[neigh].neighbor, (AT_NUMB) i, at[neigh].valence );
    // INCHI✔️❌:                     if (p1 && pneigh)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int n1 = p1 - at2[i1].neighbor;
    // INCHI✔️❌:                         int nneigh = pneigh - at[neigh].neighbor;
    // INCHI✔️❌:                         at[i].bond_type[j] = at[neigh].bond_type[nneigh] = at2[i1].bond_type[n1];
    // INCHI✔️❌:                         bSuccess = 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (!bSuccess)
    // INCHI✔️❌:                 {
    // INCHI❌❌: #ifdef _DEBUG
    // INCHI❌❌:                     int stop_here = 1;
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     num_err++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_err;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: replace_arom_bonds

    let mut errors = 0_i32;
    let mut atom_index = 0_i32;
    while atom_index < num_atoms {
        let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let valence = i32::from(
            atoms
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .valence,
        );
        let mut bond_index = 0_i32;
        while bond_index < valence {
            let bond =
                usize::try_from(bond_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let bond_type = *atoms
                .get(index)
                .and_then(|atom| atom.bond_type.get(bond))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(bond_type) > BOND_TRIPLE as i32 {
                let neighbor = usize::from(
                    *atoms
                        .get(index)
                        .and_then(|atom| atom.neighbor.get(bond))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let first_number = atoms
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .orig_at_number;
                let second_number = atoms
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .orig_at_number;
                let mut first_original = 0_i32;
                while first_original < num_original_atoms {
                    let original = original_atoms
                        .get(
                            usize::try_from(first_original)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if original.orig_at_number == first_number {
                        break;
                    }
                    first_original += 1;
                }
                let mut second_original = 0_i32;
                while second_original < num_original_atoms {
                    let original = original_atoms
                        .get(
                            usize::try_from(second_original)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if original.orig_at_number == second_number {
                        break;
                    }
                    second_original += 1;
                }
                let mut success = false;
                if first_original < num_original_atoms && second_original < num_original_atoms {
                    let first_original_index = usize::try_from(first_original)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let second_original_atom = second_original as u16;
                    let source_search_len = i32::from(
                        atoms
                            .get(first_original_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .valence,
                    );
                    let source_bond = is_in_the_list(
                        Some(
                            &original_atoms
                                .get(first_original_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .neighbor,
                        ),
                        second_original_atom,
                        source_search_len,
                    )?;
                    let reverse_bond = is_in_the_list(
                        Some(
                            &atoms
                                .get(neighbor)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .neighbor,
                        ),
                        atom_index as u16,
                        i32::from(
                            atoms
                                .get(neighbor)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .valence,
                        ),
                    )?;
                    if let (Some(source_bond), Some(reverse_bond)) = (source_bond, reverse_bond) {
                        let replacement = *original_atoms[first_original_index]
                            .bond_type
                            .get(source_bond)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        atoms[neighbor].bond_type[reverse_bond] = replacement;
                        atoms[index].bond_type[bond] = replacement;
                        success = true;
                    }
                }
                if !success {
                    errors = errors.wrapping_add(1);
                }
            }
            bond_index += 1;
        }
        atom_index += 1;
    }
    Ok(errors)
}

#[allow(non_snake_case)]
pub(crate) fn add_explicit_H(
    heap: &mut SourceHeap,
    input: &mut INP_ATOM_DATA,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5868 add_explicit_H
    // INCHI✔️❌: int add_explicit_H( INP_ATOM_DATA *inp_cur_data )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* do not care about stereo parities for now */
    // INCHI✔️❌:     int curRemovedH, num_added_explicit_H, iat, num_removed_H, m, num_H, num_atoms = inp_cur_data->num_at;
    // INCHI✔️❌:     if (( num_removed_H = inp_cur_data->num_removed_H ) > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inp_ATOM *at_H = (inp_ATOM *) inchi_calloc( num_removed_H, sizeof( inp_ATOM ) );
    // INCHI✔️❌:         inp_ATOM *at = inp_cur_data->at;
    // INCHI✔️❌:         for (curRemovedH = num_atoms, num_added_explicit_H = 0; curRemovedH < num_atoms + num_removed_H; curRemovedH++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[curRemovedH].el_number == EL_NUMBER_H &&
    // INCHI✔️❌:                  1 == at[curRemovedH].valence &&
    // INCHI✔️❌:                  ( iat = (int) at[curRemovedH].neighbor[0] ) < num_atoms &&
    // INCHI✔️❌:                  0 <= ( m = at[curRemovedH].iso_atw_diff ) &&
    // INCHI✔️❌:                  m <= NUM_H_ISOTOPES)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* num_H is the total number of all implicit H including isotopic H */
    // INCHI✔️❌:                 if (at_H && 0 < at[iat].num_H && 0 <= ( num_H = at[iat].num_H - NUM_ISO_H( at, iat ) ) && ( (!m && num_H) || (m && at[iat].num_iso_H[m - 1]) )) /* djb-rwth: addressing LLVM warning; fixing a NULL pointer dereference */
    // INCHI✔️❌:                 { /* number of implicit H > 0 */
    // INCHI✔️❌:                     int val = at[iat].valence;
    // INCHI✔️❌:                     /* set hydrogen atom */
    // INCHI✔️❌:                     at_H[num_added_explicit_H] = at[curRemovedH];
    // INCHI✔️❌:                     at_H[num_added_explicit_H].neighbor[0] = iat;
    // INCHI✔️❌:                     at_H[num_added_explicit_H].valence = 1;
    // INCHI✔️❌:                     at_H[num_added_explicit_H].chem_bonds_valence = at_H[num_added_explicit_H].bond_type[0];
    // INCHI✔️❌:                     /* set heavy atom */
    // INCHI✔️❌:                     at[iat].neighbor[val] = num_atoms + num_added_explicit_H;
    // INCHI✔️❌:                     at[iat].bond_type[val] = at_H[num_added_explicit_H].bond_type[0];
    // INCHI✔️❌:                     at[iat].bond_stereo[val] = -at_H[num_added_explicit_H].bond_stereo[0];
    // INCHI✔️❌:                     at[iat].valence++;
    // INCHI✔️❌:                     if (BOND_SINGLE <= at[iat].bond_type[val] && at[iat].bond_type[val] <= BOND_TRIPLE)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at[iat].chem_bonds_valence += at[iat].bond_type[val];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* should not happen */
    // INCHI✔️❌:                         at_H[num_added_explicit_H].bond_type[0] = at[iat].bond_type[val] = BOND_SINGLE;
    // INCHI✔️❌:                         at[iat].chem_bonds_valence += BOND_SINGLE;
    // INCHI✔️❌:                         at_H[num_added_explicit_H].chem_bonds_valence = BOND_SINGLE;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     at[iat].num_H--;
    // INCHI✔️❌:                     if (m)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at[iat].num_iso_H[m - 1] --;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     num_added_explicit_H++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (0 < num_added_explicit_H && num_added_explicit_H <= num_removed_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memcpy(at + num_atoms, at_H, num_added_explicit_H * sizeof(at));
    // INCHI✔️❌:             inp_cur_data->num_removed_H = 0;
    // INCHI✔️❌:             inp_cur_data->num_at = ( num_atoms += num_added_explicit_H );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_free( at_H );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return num_atoms;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: add_explicit_H

    let mut num_atoms = input.num_at;
    let num_removed_h = input.num_removed_H;
    if num_removed_h <= 0 {
        return Ok(num_atoms);
    }
    let temporary = inchi_calloc::<inp_ATOM>(heap, num_removed_h as u64, 176)
        .unwrap_or_else(|_| SourceMutPointer::null());
    let processing = (|| -> Result<i32, SourceHeapError> {
        let end = num_atoms
            .checked_add(num_removed_h)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let mut removed_index = num_atoms;
        let mut added = 0_i32;
        while removed_index < end {
            let removed = heap
                .slice(input.at.as_const().offset(i64::from(removed_index))?)?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let heavy_index = i32::from(removed.neighbor[0]);
            let isotope = i32::from(removed.iso_atw_diff);
            if removed.el_number == EL_NUMBER_H
                && removed.valence == 1
                && heavy_index < num_atoms
                && isotope >= 0
                && isotope <= 3
                && !temporary.is_null()
            {
                let heavy = heap
                    .slice(input.at.as_const().offset(i64::from(heavy_index))?)?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let ordinary_h = i32::from(heavy.num_H)
                    - heavy
                        .num_iso_H
                        .iter()
                        .map(|&value| i32::from(value))
                        .sum::<i32>();
                let isotope_available = isotope != 0
                    && *heavy
                        .num_iso_H
                        .get(
                            usize::try_from(isotope - 1)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0;
                if heavy.num_H > 0
                    && ordinary_h >= 0
                    && ((isotope == 0 && ordinary_h != 0) || isotope_available)
                {
                    let mut hydrogen = removed;
                    hydrogen.neighbor[0] = heavy_index as u16;
                    hydrogen.valence = 1;
                    hydrogen.chem_bonds_valence = hydrogen.bond_type[0] as i8;
                    *heap
                        .slice_mut(temporary.offset(i64::from(added))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = hydrogen.clone();

                    let invalid_bond = {
                        let heavy_atom = heap
                            .slice_mut(input.at.offset(i64::from(heavy_index))?)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let valence = usize::try_from(i32::from(heavy_atom.valence))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        *heavy_atom
                            .neighbor
                            .get_mut(valence)
                            .ok_or(SourceHeapError::PointerOutOfBounds)? =
                            num_atoms.wrapping_add(added) as u16;
                        *heavy_atom
                            .bond_type
                            .get_mut(valence)
                            .ok_or(SourceHeapError::PointerOutOfBounds)? = hydrogen.bond_type[0];
                        *heavy_atom
                            .bond_stereo
                            .get_mut(valence)
                            .ok_or(SourceHeapError::PointerOutOfBounds)? =
                            hydrogen.bond_stereo[0].wrapping_neg();
                        heavy_atom.valence = heavy_atom.valence.wrapping_add(1);
                        let bond_type = heavy_atom.bond_type[valence];
                        let invalid =
                            bond_type < BOND_SINGLE as u8 || bond_type > BOND_TRIPLE as u8;
                        if invalid {
                            heavy_atom.bond_type[valence] = BOND_SINGLE as u8;
                            heavy_atom.chem_bonds_valence = heavy_atom
                                .chem_bonds_valence
                                .wrapping_add(BOND_SINGLE as i8);
                        } else {
                            heavy_atom.chem_bonds_valence =
                                heavy_atom.chem_bonds_valence.wrapping_add(bond_type as i8);
                        }
                        invalid
                    };
                    if invalid_bond {
                        let temporary_h = heap
                            .slice_mut(temporary.offset(i64::from(added))?)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        temporary_h.bond_type[0] = BOND_SINGLE as u8;
                        temporary_h.chem_bonds_valence = BOND_SINGLE as i8;
                    }
                    let heavy_atom = heap
                        .slice_mut(input.at.offset(i64::from(heavy_index))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    heavy_atom.num_H = heavy_atom.num_H.wrapping_sub(1);
                    if isotope != 0 {
                        let isotope_count = heavy_atom
                            .num_iso_H
                            .get_mut(
                                usize::try_from(isotope - 1)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        *isotope_count = isotope_count.wrapping_sub(1);
                    }
                    added = added.wrapping_add(1);
                }
            }
            removed_index = removed_index.wrapping_add(1);
        }
        if added > 0 && added <= num_removed_h {
            let temporary_atoms = heap
                .slice(temporary.as_const())?
                .get(..usize::try_from(added).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let destination = heap.slice_mut(input.at.offset(i64::from(num_atoms))?)?;
            copy_inp_atom_gcc_lp64_byte_prefix(
                destination,
                &temporary_atoms,
                usize::try_from(added)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    .checked_mul(8)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?,
            )?;
            input.num_removed_H = 0;
            num_atoms = num_atoms.wrapping_add(added);
            input.num_at = num_atoms;
        }
        Ok(num_atoms)
    })();
    inchi_free(heap, temporary)?;
    processing
}

pub(crate) fn fill_out_bond_cuts(
    atoms: &[inp_ATOM],
    derivatives: &[DERIV_AT],
    num_atoms: i32,
    atom_pairs: &mut [R2C_ATPAIR],
    num_cuts_to_check: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5628 fill_out_bond_cuts
    // INCHI✔️❌: int fill_out_bond_cuts( inp_ATOM *at,
    // INCHI✔️❌:                         DERIV_AT *da,
    // INCHI✔️❌:                         int num_atoms,
    // INCHI✔️❌:                         R2C_ATPAIR *ap,
    // INCHI✔️❌:                         int num_cuts_to_check )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, n;
    // INCHI✔️❌:     AT_NUMB at1, at2;
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     /* fill out the array of bonds to be cut */
    // INCHI✔️❌:     for (i = j = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR ) &&
    // INCHI✔️❌:              da[i].num[0] <= MAX_AT_DERIV && da[i].num[1] <= MAX_AT_DERIV)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (j + 1 >= num_cuts_to_check)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = -2;
    // INCHI✔️❌:                 goto exit_r2c_num; /* wrong number of cuts = num */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (k = 0; k < 2; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at1 = i;
    // INCHI✔️❌:                 at2 = at[at1].neighbor[(int) da[at1].ord[k]];
    // INCHI✔️❌:                 n = ( at1 > at2 );
    // INCHI✔️❌:                 ap[j].at[n] = at1;
    // INCHI✔️❌:                 ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
    // INCHI✔️❌:                 ap[j].atno = i;
    // INCHI✔️❌:                 j++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (0 < cmp_r2c_atpair( ap + j - 2, ap + j - 1 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 R2C_ATPAIR ap1 = ap[j - 2];
    // INCHI✔️❌:                 ap[j - 2] = ap[j - 1];
    // INCHI✔️❌:                 ap[j - 1] = ap1; /* sort each pair */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* 2013-12-04 DT */
    // INCHI✔️❌:             if (da[i].typ[0] & DERIV_RING_DMOX_DEOX)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int other_atom = (int) da[i].other_atom - 1;
    // INCHI✔️❌:                 if (da[i].typ[1] || other_atom < 0 || i == other_atom || da[other_atom].other_atom != i + 1 ||
    // INCHI✔️❌:                      !( da[other_atom].typ[0] & DERIV_RING_DMOX_DEOX ) || da[other_atom].typ[1])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = -3;
    // INCHI✔️❌:                     goto exit_r2c_num;
    // INCHI✔️❌:                     /* no other cut may be at the atom in addition to DERIV_RING_DMOX_DEOX
    // INCHI✔️❌:                     or no other_atom or other_atom has wrong deriv. type */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* make sure the ap[] for the two cuts in the same ring are adjacent */
    // INCHI✔️❌:                 if (other_atom > i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (j + 1 >= num_cuts_to_check)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = -2;
    // INCHI✔️❌:                         goto exit_r2c_num; /* wrong number of cuts = num */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* cut #1 */
    // INCHI✔️❌:                     at1 = i;
    // INCHI✔️❌:                     at2 = at[at1].neighbor[(int) da[at1].ord[0]];
    // INCHI✔️❌:                     n = ( at1 > at2 );
    // INCHI✔️❌:                     ap[j].at[n] = at1;
    // INCHI✔️❌:                     ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
    // INCHI✔️❌:                     ap[j].atno = i;
    // INCHI✔️❌:                     j++;
    // INCHI✔️❌:                     /* cut #2 */
    // INCHI✔️❌:                     at1 = other_atom;
    // INCHI✔️❌:                     at2 = at[at1].neighbor[(int) da[at1].ord[0]];
    // INCHI✔️❌:                     n = ( at1 > at2 );
    // INCHI✔️❌:                     ap[j].at[n] = at1;
    // INCHI✔️❌:                     ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
    // INCHI✔️❌:                     ap[j].atno = i;
    // INCHI✔️❌:                     j++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* add each pair of cuts only once */
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #if ( COUNT_ALL_NOT_DERIV == 1 )
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (k = 0; k < DERIV_AT_LEN && da[i].typ[k]; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (j >= num_cuts_to_check || ( da[i].typ[k] & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = -2;
    // INCHI✔️❌:                         goto exit_r2c_num; /* wrong number of cuts = num or wrong type */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     at1 = i;
    // INCHI✔️❌:                     at2 = at[i].neighbor[(int) da[i].ord[k]];
    // INCHI✔️❌:                     n = ( at1 > at2 );
    // INCHI✔️❌:                     /* pair of atoms possibly to be disconnected */
    // INCHI✔️❌:                     ap[j].at[n] = at1;
    // INCHI✔️❌:                     ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
    // INCHI✔️❌:                                                /* precursor's atom */
    // INCHI✔️❌:                     ap[j].atno = i;
    // INCHI✔️❌:                     j++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return j;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_r2c_num:
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: fill_out_bond_cuts

    let mut atom_index = 0_i32;
    let mut pair_count = 0_i32;
    while atom_index < num_atoms {
        let derivative_index =
            usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let derivative = derivatives
            .get(derivative_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let first_type = i32::from(derivative.typ[0]);
        let second_type = i32::from(derivative.typ[1]);

        if first_type & DERIV_RING_OUTSIDE_PRECURSOR as i32 != 0
            && second_type & DERIV_RING_OUTSIDE_PRECURSOR as i32 != 0
            && i32::from(derivative.num[0]) <= MAX_AT_DERIV as i32
            && i32::from(derivative.num[1]) <= MAX_AT_DERIV as i32
        {
            if pair_count + 1 >= num_cuts_to_check {
                return Ok(-2);
            }
            for cut_index in 0..2_usize {
                let first_atom = atom_index as u16;
                let order = usize::try_from(i32::from(derivative.ord[cut_index]))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let second_atom = *atoms
                    .get(usize::from(first_atom))
                    .and_then(|atom| atom.neighbor.get(order))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let position = usize::from(first_atom > second_atom);
                let pair = atom_pairs
                    .get_mut(
                        usize::try_from(pair_count)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                pair.at[position] = first_atom;
                pair.at[1 - position] = second_atom;
                pair.atno = atom_index as u16;
                pair_count += 1;
            }
            let first_pair =
                usize::try_from(pair_count - 2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let second_pair =
                usize::try_from(pair_count - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if cmp_r2c_atpair(&atom_pairs[first_pair], &atom_pairs[second_pair]) > 0 {
                atom_pairs.swap(first_pair, second_pair);
            }
        } else if first_type & DERIV_RING_DMOX_DEOX as i32 != 0 {
            let other_atom = i32::from(derivative.other_atom) - 1;
            if second_type != 0 || other_atom < 0 || atom_index == other_atom {
                return Ok(-3);
            }
            let other_index =
                usize::try_from(other_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let other_derivative = derivatives
                .get(other_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(other_derivative.other_atom) != atom_index + 1
                || i32::from(other_derivative.typ[0]) & DERIV_RING_DMOX_DEOX as i32 == 0
                || other_derivative.typ[1] != 0
            {
                return Ok(-3);
            }
            if other_atom > atom_index {
                if pair_count + 1 >= num_cuts_to_check {
                    return Ok(-2);
                }
                for source_atom in [atom_index, other_atom] {
                    let source_index = usize::try_from(source_atom)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let source_derivative = derivatives
                        .get(source_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let order = usize::try_from(i32::from(source_derivative.ord[0]))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let first_atom = source_atom as u16;
                    let second_atom = *atoms
                        .get(usize::from(first_atom))
                        .and_then(|atom| atom.neighbor.get(order))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let position = usize::from(first_atom > second_atom);
                    let pair = atom_pairs
                        .get_mut(
                            usize::try_from(pair_count)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    pair.at[position] = first_atom;
                    pair.at[1 - position] = second_atom;
                    pair.atno = atom_index as u16;
                    pair_count += 1;
                }
            }
        } else {
            for cut_index in 0..derivative.typ.len() {
                let type_ = i32::from(derivative.typ[cut_index]);
                if type_ == 0 {
                    break;
                }
                if pair_count >= num_cuts_to_check
                    || type_ & DERIV_RING_OUTSIDE_PRECURSOR as i32 != 0
                {
                    return Ok(-2);
                }
                let order = usize::try_from(i32::from(derivative.ord[cut_index]))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let first_atom = atom_index as u16;
                let second_atom = *atoms
                    .get(derivative_index)
                    .and_then(|atom| atom.neighbor.get(order))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let position = usize::from(first_atom > second_atom);
                let pair = atom_pairs
                    .get_mut(
                        usize::try_from(pair_count)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                pair.at[position] = first_atom;
                pair.at[1 - position] = second_atom;
                pair.atno = atom_index as u16;
                pair_count += 1;
            }
        }
        atom_index += 1;
    }
    Ok(pair_count)
}

#[allow(non_snake_case)]
pub(crate) fn UnMarkOtherIndicators(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:681 UnMarkOtherIndicators
    // INCHI✔️❌: int UnMarkOtherIndicators( inp_ATOM *at, int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[i].at_type = 0;
    // INCHI✔️❌:         at[i].cFlags = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: UnMarkOtherIndicators

    let mut index = 0_i32;
    while index < num_atoms {
        let atom = atoms
            .get_mut(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.at_type = 0;
        atom.cFlags = 0;
        index += 1;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn UnMarkDisconnectedComponents(
    heap: &mut SourceHeap,
    original: &mut ORIG_ATOM_DATA,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:652 UnMarkDisconnectedComponents
    // INCHI✔️❌: int UnMarkDisconnectedComponents( ORIG_ATOM_DATA *orig_inp_data )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < orig_inp_data->num_inp_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         orig_inp_data->at[i].orig_compt_at_numb = 0;
    // INCHI✔️❌:         orig_inp_data->at[i].component = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_inp_data->nCurAtLen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( orig_inp_data->nCurAtLen );
    // INCHI✔️❌:         orig_inp_data->nCurAtLen = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_inp_data->nOldCompNumber)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( orig_inp_data->nOldCompNumber );
    // INCHI✔️❌:         orig_inp_data->nOldCompNumber = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_inp_data->num_components = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: UnMarkDisconnectedComponents

    let mut index = 0_i32;
    while index < original.num_inp_atoms {
        let atom = heap
            .slice_mut(original.at.offset(i64::from(index))?)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.orig_compt_at_numb = 0;
        atom.component = 0;
        index += 1;
    }
    if !original.nCurAtLen.is_null() {
        inchi_free(heap, original.nCurAtLen)?;
        original.nCurAtLen = SourceMutPointer::null();
    }
    if !original.nOldCompNumber.is_null() {
        inchi_free(heap, original.nOldCompNumber)?;
        original.nOldCompNumber = SourceMutPointer::null();
    }
    original.num_components = 0;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn UnMarkOneComponent(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:698 UnMarkOneComponent
    // INCHI✔️❌: int UnMarkOneComponent( inp_ATOM *at, int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[i].orig_compt_at_numb = 0;
    // INCHI✔️❌:         at[i].component = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: UnMarkOneComponent

    let mut index = 0_i32;
    while index < num_atoms {
        let atom = atoms
            .get_mut(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.orig_compt_at_numb = 0;
        atom.component = 0;
        index += 1;
    }
    Ok(0)
}

pub(crate) fn subtract_DT_from_num_H(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:734 subtract_DT_from_num_H
    // INCHI✔️❌: int subtract_DT_from_num_H( int num_atoms, inp_ATOM *at )
    // INCHI✔️❌: /*  assume num_1H, num_D and num_T are included in num_H */
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j;
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = 0; j < NUM_H_ISOTOPES; j++)
    // INCHI✔️❌:             at[i].num_H -= at[i].num_iso_H[j];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: subtract_DT_from_num_H

    let mut atom_index = 0_i32;
    while atom_index < num_atoms {
        let atom = atoms
            .get_mut(usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for isotope_count in atom.num_iso_H {
            atom.num_H = atom.num_H.wrapping_sub(isotope_count);
        }
        atom_index += 1;
    }
    Ok(0)
}

pub(crate) fn add_inp_ATOM(
    atoms: &mut [inp_ATOM],
    len_at: i32,
    len_cur: i32,
    add: &[inp_ATOM],
    len_add: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:749 add_inp_ATOM
    // INCHI✔️❌: int add_inp_ATOM( inp_ATOM *at,
    // INCHI✔️❌:                   int len_at,
    // INCHI✔️❌:                   int len_cur,
    // INCHI✔️❌:                   inp_ATOM *add,
    // INCHI✔️❌:                   int len_add )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j;
    // INCHI✔️❌:     inp_ATOM *a;
    // INCHI✔️❌:     /* chack correctness */
    // INCHI✔️❌:     if (len_cur < 0)
    // INCHI✔️❌:         return len_cur;
    // INCHI✔️❌:     if (len_add < 0)
    // INCHI✔️❌:         return len_add;
    // INCHI✔️❌:     if (len_cur + len_add > len_at)
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     /* copy */
    // INCHI✔️❌:     memcpy(at + len_cur, add, len_add * sizeof(at[0]));
    // INCHI✔️❌:     /* modify */
    // INCHI✔️❌:     if (len_cur)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         a = at + len_cur;
    // INCHI✔️❌:         for (i = 0; i < len_add; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = 0; j < a[i].valence; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 a[i].neighbor[j] += len_cur;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return len_cur + len_add;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: add_inp_ATOM

    if len_cur < 0 {
        return Ok(len_cur);
    }
    if len_add < 0 {
        return Ok(len_add);
    }
    let result_len = len_cur.wrapping_add(len_add);
    if result_len > len_at {
        return Ok(-1);
    }

    let destination_start =
        usize::try_from(len_cur).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let added_len = usize::try_from(len_add).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let destination_end = destination_start
        .checked_add(added_len)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let source = add
        .get(..added_len)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    atoms
        .get_mut(destination_start..destination_end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone_from_slice(source);

    if len_cur != 0 {
        let offset = len_cur as u16;
        for atom in atoms
            .get_mut(destination_start..destination_end)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            let valence = usize::try_from(i32::from(atom.valence).max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            for neighbor in atom
                .neighbor
                .get_mut(..valence)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                *neighbor = neighbor.wrapping_add(offset);
            }
        }
    }
    Ok(result_len)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Edit_MergeComponentsAndRecreateOAD(
    heap: &mut SourceHeap,
    original: &mut ORIG_ATOM_DATA,
    current: &mut [INP_ATOM_DATA],
    num_components: i32,
    error_code: &mut i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7563 OAD_Edit_MergeComponentsAndRecreateOAD
    // INCHI✔️❌: void OAD_Edit_MergeComponentsAndRecreateOAD( ORIG_ATOM_DATA *orig_OrigAtomData,
    // INCHI✔️❌:                                              INP_ATOM_DATA *curr_InpAtomData,
    // INCHI✔️❌:                                              int num_components,
    // INCHI✔️❌:                                              int *errcode )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, num_atoms = 0, cur_num_at = 0;
    // INCHI✔️❌:     inp_ATOM *at;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_components <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *errcode = -999; /* num atoms mismatch */
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Merge kept components into 'at' */
    // INCHI✔️❌:     for (i = 0; i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_atoms += curr_InpAtomData[i].num_at;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     at = (inp_ATOM *) inchi_calloc( num_atoms, sizeof( at[0] ) );
    // INCHI✔️❌:     cur_num_at = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Clean and prepare */
    // INCHI✔️❌:         UnMarkRingSystemsInp( curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );
    // INCHI✔️❌:         UnMarkOtherIndicators( curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );
    // INCHI✔️❌:         UnMarkOneComponent( curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );
    // INCHI✔️❌:
    // INCHI✔️❌:         subtract_DT_from_num_H( curr_InpAtomData[i].num_at, curr_InpAtomData[i].at );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Merge one by one */
    // INCHI✔️❌:         cur_num_at = add_inp_ATOM( at, num_atoms, cur_num_at, curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Replace original OrigAtomData structure */
    // INCHI✔️❌:     if (cur_num_at == num_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( orig_OrigAtomData->at );
    // INCHI✔️❌:         orig_OrigAtomData->at = at;
    // INCHI✔️❌:
    // INCHI✔️❌:         orig_OrigAtomData->num_inp_atoms = cur_num_at;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Destroy original coordinates as we destroyed a part of original input structure */
    // INCHI✔️❌:         if (orig_OrigAtomData->szCoord)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( orig_OrigAtomData->szCoord );
    // INCHI✔️❌:             orig_OrigAtomData->szCoord = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         UnMarkDisconnectedComponents( orig_OrigAtomData );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Create copy error! */
    // INCHI✔️❌:         if (at)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( at );
    // INCHI✔️❌:             at = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         *errcode = -999; /* num atoms mismatch */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_Edit_MergeComponentsAndRecreateOAD

    if num_components <= 0 {
        *error_code = -999;
        return Ok(());
    }
    let component_count =
        usize::try_from(num_components).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let components = current
        .get_mut(..component_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut num_atoms = 0_i32;
    for component in components.iter() {
        num_atoms = num_atoms.wrapping_add(component.num_at);
    }
    let atom_count =
        u64::try_from(num_atoms).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let merged = inchi_calloc::<inp_ATOM>(heap, atom_count, 176)?;
    let merge_result = (|| -> Result<i32, SourceHeapError> {
        let mut current_count = 0_i32;
        for component in components.iter_mut() {
            UnMarkRingSystemsInp(heap, component.at, component.num_at)?;
            {
                let atoms = heap.slice_mut(component.at)?;
                UnMarkOtherIndicators(atoms, component.num_at)?;
                UnMarkOneComponent(atoms, component.num_at)?;
                subtract_DT_from_num_H(atoms, component.num_at)?;
            }
            let component_len = usize::try_from(component.num_at)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let component_atoms = heap
                .slice(component.at.as_const())?
                .get(..component_len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            current_count = add_inp_ATOM(
                heap.slice_mut(merged)?,
                num_atoms,
                current_count,
                &component_atoms,
                component.num_at,
            )?;
        }
        Ok(current_count)
    })();
    let current_count = match merge_result {
        Ok(value) => value,
        Err(error) => {
            inchi_free(heap, merged)?;
            return Err(error);
        }
    };

    if current_count == num_atoms {
        inchi_free(heap, original.at)?;
        original.at = merged;
        original.num_inp_atoms = current_count;
        if !original.szCoord.is_null() {
            inchi_free(heap, original.szCoord)?;
            original.szCoord = SourceMutPointer::null();
        }
        UnMarkDisconnectedComponents(heap, original)?;
    } else {
        inchi_free(heap, merged)?;
        *error_code = -999;
    }
    Ok(())
}

fn oad_allocate_atom_pairs(
    heap: &mut SourceHeap,
    pointer: &mut SourceMutPointer<R2C_ATPAIR>,
    allocated: &mut i32,
    required: i32,
) -> Result<bool, SourceHeapError> {
    if required > 0 && (*allocated < required || pointer.is_null()) {
        if !pointer.is_null() {
            inchi_free(heap, *pointer)?;
            *pointer = SourceMutPointer::null();
        }
        let count = usize::try_from(required)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        match heap.allocate(vec![R2C_ATPAIR::default(); count]) {
            Ok(allocation) => {
                *pointer = allocation;
                *allocated = required;
            }
            Err(SourceHeapError::AllocationFailed) => return Ok(false),
            Err(error) => return Err(error),
        }
    }
    Ok(true)
}

fn oad_prepare_derivative_cuts(
    atoms: &[inp_ATOM],
    derivatives: &mut [DERIV_AT],
    num_atoms: i32,
) -> Result<Result<(i32, i32, i32), i32>, SourceHeapError> {
    let atom_count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if atoms.len() < atom_count || derivatives.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut num_cuts = 0_i32;
    let mut num_ring_cuts = 0_i32;
    let mut num_cut_pieces = 0_i32;

    for atom_index in 0..atom_count {
        let length = derivatives[atom_index]
            .typ
            .iter()
            .take_while(|type_| **type_ != 0)
            .count();
        match length {
            0 => continue,
            1 => {
                if derivatives[atom_index].typ[0] & DERIV_RING_DMOX_DEOX as i16 != 0 {
                    let other = usize::from(derivatives[atom_index].other_atom.wrapping_sub(1));
                    if derivatives[atom_index].other_atom == 0
                        || other >= atom_count
                        || (derivatives[other].typ[0] ^ DERIV_RING_DMOX_DEOX as i16)
                            != derivatives[atom_index].typ[0]
                        || usize::from(derivatives[other].other_atom.wrapping_sub(1)) != atom_index
                    {
                        derivatives[atom_index].num[0] = NOT_AT_DERIV as i8;
                    } else {
                        num_cuts += 1;
                        num_ring_cuts += 1;
                        num_cut_pieces += i32::from(atom_index > other);
                    }
                } else {
                    num_cuts += 1;
                    num_cut_pieces += 1;
                }
            }
            2 => {
                let first_type = derivatives[atom_index].typ[0];
                let second_type = derivatives[atom_index].typ[1];
                if first_type & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0
                    && second_type & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0
                {
                    num_ring_cuts += 2;
                    num_cuts += 2;
                    num_cut_pieces += 1;
                } else if first_type != 0
                    && first_type & DERIV_RING2_OUTSIDE_PRECUR as i16 == first_type
                    && second_type == first_type
                {
                    num_ring_cuts += 2;
                    num_cuts += 2;
                    num_cut_pieces += 1;
                } else if first_type == DERIV_AMINE_tN as i16
                    && second_type == DERIV_AMINE_tN as i16
                {
                    num_cuts += 2;
                    num_cut_pieces += 2;
                } else if first_type == DERIV_RO_COX as i16 && second_type == DERIV_RO_COX as i16 {
                    if derivatives[atom_index].num[0] == derivatives[atom_index].num[1] {
                        derivatives[atom_index] = DERIV_AT::default();
                    } else {
                        let first_num = derivatives[atom_index].num[0];
                        let second_num = derivatives[atom_index].num[1];
                        let selected = if first_num != 0 && second_num != 0 {
                            const PRIORITY: [i8; 6] = [12, 9, 6, 13, 3, 8];
                            let first = PRIORITY.iter().position(|value| *value == first_num);
                            let second = PRIORITY.iter().position(|value| *value == second_num);
                            match (first, second) {
                                (Some(first), Some(second)) => i32::from(second < first),
                                (Some(_), None) => 0,
                                (None, Some(_)) => 1,
                                (None, None) => -1,
                            }
                        } else if first_num != 0 {
                            0
                        } else if second_num != 0 {
                            1
                        } else {
                            -1
                        };
                        match selected {
                            1 => {
                                derivatives[atom_index].num[0] = second_num;
                                derivatives[atom_index].ord[0] = derivatives[atom_index].ord[1];
                                derivatives[atom_index].typ[0] = second_type;
                                derivatives[atom_index].typ[1] = 0;
                                derivatives[atom_index].num[1] = 0;
                                derivatives[atom_index].ord[1] = 0;
                                num_cuts += 1;
                                num_cut_pieces += 1;
                            }
                            0 => {
                                derivatives[atom_index].typ[1] = 0;
                                derivatives[atom_index].num[1] = 0;
                                derivatives[atom_index].ord[1] = 0;
                                num_cuts += 1;
                                num_cut_pieces += 1;
                            }
                            _ => derivatives[atom_index] = DERIV_AT::default(),
                        }
                    }
                } else if first_type == second_type {
                    let mut silyl = [0_i32; 2];
                    for derivative_index in 0..2 {
                        if is_deriv_chain(
                            atoms,
                            atom_index as i32,
                            num_atoms,
                            &derivatives[atom_index],
                            derivative_index as i32,
                            None,
                            0,
                            None,
                            0,
                            None,
                        )? == 0
                        {
                            derivatives[atom_index].num[derivative_index] = NOT_AT_DERIV as i8;
                        } else {
                            let ordinal = i32::from(derivatives[atom_index].ord[derivative_index])
                                - i32::from(b'0');
                            let neighbor = usize::from(
                                atoms[atom_index].neighbor[usize::try_from(ordinal)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                            );
                            silyl[derivative_index] =
                                is_silyl2(atoms, neighbor as i32, atom_index as i32)?;
                        }
                    }
                    let first_num = derivatives[atom_index].num[0];
                    let second_num = derivatives[atom_index].num[1];
                    if (silyl[1] != 0 && (silyl[0] == 0 || silyl[1] < silyl[0]))
                        || (silyl == [0, 0] && first_num > second_num)
                    {
                        derivatives[atom_index].num[0] = second_num;
                        derivatives[atom_index].ord[0] = derivatives[atom_index].ord[1];
                        derivatives[atom_index].typ[0] = second_type;
                        derivatives[atom_index].typ[1] = 0;
                        num_cuts += 1;
                        num_cut_pieces += 1;
                    } else if (silyl[0] != 0 && (silyl[1] == 0 || silyl[0] < silyl[1]))
                        || (silyl == [0, 0] && first_num < second_num)
                    {
                        derivatives[atom_index].typ[1] = 0;
                        num_cuts += 1;
                        num_cut_pieces += 1;
                    } else {
                        derivatives[atom_index].typ[0] = 0;
                        derivatives[atom_index].typ[1] = 0;
                    }
                } else if first_type & (DERIV_RO_COX | DERIV_BRIDGE_O) as i16 != 0
                    && second_type & (DERIV_RO_COX | DERIV_BRIDGE_O) as i16 != 0
                {
                    if first_type == DERIV_BRIDGE_O as i16 {
                        derivatives[atom_index].typ[0] = second_type;
                        derivatives[atom_index].num[0] = derivatives[atom_index].num[1];
                        derivatives[atom_index].ord[0] = derivatives[atom_index].ord[1];
                    }
                    derivatives[atom_index].typ[1] = 0;
                    derivatives[atom_index].num[1] = 0;
                    derivatives[atom_index].ord[1] = 0;
                    num_cuts += 1;
                    num_cut_pieces += 1;
                } else {
                    return Ok(Err(-88));
                }
            }
            3 => {
                if derivatives[atom_index].typ[0] != DERIV_AMINE_tN as i16
                    || derivatives[atom_index].typ[1] != derivatives[atom_index].typ[0]
                    || derivatives[atom_index].typ[2] != derivatives[atom_index].typ[0]
                {
                    return Ok(Err(-88));
                }
                let mut silyl = [0_i32; 3];
                for derivative_index in 0..3 {
                    if is_deriv_chain(
                        atoms,
                        atom_index as i32,
                        num_atoms,
                        &derivatives[atom_index],
                        derivative_index as i32,
                        None,
                        0,
                        None,
                        0,
                        None,
                    )? == 0
                    {
                        derivatives[atom_index].num[derivative_index] = NOT_AT_DERIV as i8;
                    } else {
                        let ordinal = i32::from(derivatives[atom_index].ord[derivative_index])
                            - i32::from(b'0');
                        let neighbor = usize::from(
                            atoms[atom_index].neighbor[usize::try_from(ordinal)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                        );
                        silyl[derivative_index] =
                            is_silyl2(atoms, neighbor as i32, atom_index as i32)?;
                    }
                }
                let less = |left: usize, right: usize, derivative: &DERIV_AT| {
                    (silyl[left] != 0 && (silyl[right] == 0 || silyl[left] < silyl[right]))
                        || (silyl[left] == 0
                            && silyl[right] == 0
                            && derivative.num[left] < derivative.num[right])
                };
                let mut smallest = if less(0, 1, &derivatives[atom_index]) {
                    0
                } else {
                    1
                };
                let mut largest = usize::from(smallest == 0);
                smallest = if less(smallest, 2, &derivatives[atom_index]) {
                    smallest
                } else {
                    2
                };
                largest = if less(smallest, 2, &derivatives[atom_index]) {
                    2
                } else {
                    largest
                };
                let middle = ((smallest + 1) ^ (largest + 1)) - 1;
                if derivatives[atom_index].num[smallest] == derivatives[atom_index].num[largest]
                    && silyl[smallest] == silyl[largest]
                {
                    derivatives[atom_index].typ[..3].fill(0);
                } else if (derivatives[atom_index].num[smallest]
                    == derivatives[atom_index].num[middle]
                    && silyl[smallest] == silyl[middle])
                    || (silyl[smallest] != 0 && silyl[middle] != 0 && silyl[largest] == 0)
                {
                    match largest {
                        0 => {
                            derivatives[atom_index].num[0] = derivatives[atom_index].num[1];
                            derivatives[atom_index].ord[0] = derivatives[atom_index].ord[1];
                            derivatives[atom_index].typ[0] = derivatives[atom_index].typ[1];
                            derivatives[atom_index].num[1] = derivatives[atom_index].num[2];
                            derivatives[atom_index].ord[1] = derivatives[atom_index].ord[2];
                            derivatives[atom_index].typ[1] = derivatives[atom_index].typ[2];
                        }
                        1 => {
                            derivatives[atom_index].num[1] = derivatives[atom_index].num[2];
                            derivatives[atom_index].ord[1] = derivatives[atom_index].ord[2];
                            derivatives[atom_index].typ[1] = derivatives[atom_index].typ[2];
                        }
                        _ => {}
                    }
                    derivatives[atom_index].typ[2] = 0;
                    num_cuts += 2;
                    num_cut_pieces += 2;
                } else {
                    if smallest != 0 {
                        derivatives[atom_index].num[0] = derivatives[atom_index].num[smallest];
                        derivatives[atom_index].ord[0] = derivatives[atom_index].ord[smallest];
                        derivatives[atom_index].typ[0] = derivatives[atom_index].typ[smallest];
                    }
                    derivatives[atom_index].typ[1] = 0;
                    derivatives[atom_index].typ[2] = 0;
                    num_cuts += 1;
                    num_cut_pieces += 1;
                }
            }
            4 => {
                let all_ring = derivatives[atom_index]
                    .typ
                    .iter()
                    .all(|type_| *type_ & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0);
                if !all_ring {
                    return Ok(Err(-88));
                }
                let first_size = derivatives[atom_index].num[0].max(derivatives[atom_index].num[1]);
                let second_size =
                    derivatives[atom_index].num[2].max(derivatives[atom_index].num[3]);
                if first_size < second_size
                    && is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                        atoms,
                        atom_index as i32,
                        num_atoms,
                        &derivatives[atom_index],
                        0,
                        None,
                        0,
                        None,
                        0,
                        None,
                    )? > 0
                {
                    derivatives[atom_index].typ[2] = 0;
                    derivatives[atom_index].typ[3] = 0;
                    num_cuts += 2;
                    num_ring_cuts += 2;
                    num_cut_pieces += 1;
                } else if first_size > second_size
                    && is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                        atoms,
                        atom_index as i32,
                        num_atoms,
                        &derivatives[atom_index],
                        2,
                        None,
                        0,
                        None,
                        0,
                        None,
                    )? > 0
                {
                    derivatives[atom_index].num[0] = derivatives[atom_index].num[2];
                    derivatives[atom_index].ord[0] = derivatives[atom_index].ord[2];
                    derivatives[atom_index].typ[0] = derivatives[atom_index].typ[2];
                    derivatives[atom_index].num[1] = derivatives[atom_index].num[3];
                    derivatives[atom_index].ord[1] = derivatives[atom_index].ord[3];
                    derivatives[atom_index].typ[1] = derivatives[atom_index].typ[3];
                    derivatives[atom_index].typ[2] = 0;
                    derivatives[atom_index].typ[3] = 0;
                    num_cuts += 2;
                    num_ring_cuts += 2;
                    num_cut_pieces += 1;
                } else {
                    derivatives[atom_index].typ.fill(0);
                }
            }
            _ => return Ok(Err(-88)),
        }
    }

    let mut found = 0_i32;
    for atom_index in 0..atom_count {
        let mut derivative_index = 0_usize;
        while derivative_index < derivatives[atom_index].typ.len()
            && derivatives[atom_index].typ[derivative_index] != 0
        {
            if derivatives[atom_index].typ[derivative_index] & DERIV_DUPLIC as i16 != 0 {
                derivative_index += 1;
                continue;
            }
            let ordinal = usize::try_from(derivatives[atom_index].ord[derivative_index])
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let neighbor = usize::from(atoms[atom_index].neighbor[ordinal]);
            if neighbor < atom_index {
                derivative_index += 1;
                continue;
            }
            if neighbor < atom_count {
                let mut neighbor_derivative = 0_usize;
                while neighbor_derivative < derivatives[neighbor].typ.len()
                    && derivatives[neighbor].typ[neighbor_derivative] != 0
                {
                    if derivatives[neighbor].typ[neighbor_derivative] & DERIV_DUPLIC as i16 == 0 {
                        let neighbor_ordinal =
                            usize::try_from(derivatives[neighbor].ord[neighbor_derivative])
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if usize::from(atoms[neighbor].neighbor[neighbor_ordinal]) == atom_index {
                            let mut first_index = derivative_index as i32;
                            let mut second_index = neighbor_derivative as i32;
                            let first_result = is_deriv_chain_or_ring(
                                atoms,
                                atom_index as i32,
                                num_atoms,
                                &derivatives[atom_index],
                                &mut first_index,
                            )?;
                            let second_result = is_deriv_chain_or_ring(
                                atoms,
                                neighbor as i32,
                                num_atoms,
                                &derivatives[neighbor],
                                &mut second_index,
                            )?;
                            if first_result < 0 {
                                return Ok(Err(first_result));
                            }
                            if second_result < 0 {
                                return Ok(Err(second_result));
                            }
                            if first_result == 0 || (first_result != 0 && second_result != 0) {
                                let selected = usize::try_from(first_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                if derivatives[atom_index].typ[selected]
                                    & DERIV_RING_OUTSIDE_PRECURSOR as i16
                                    != 0
                                {
                                    num_cuts -= 2;
                                    num_ring_cuts -= 2;
                                } else {
                                    num_cuts -= 1;
                                }
                                num_cut_pieces -= 1;
                                let result =
                                    remove_deriv_mark(&mut derivatives[atom_index], first_index)?;
                                if result != 0 {
                                    return Ok(Err(result));
                                }
                                found += 1;
                            }
                            if second_result == 0 || (first_result != 0 && second_result != 0) {
                                let selected = usize::try_from(second_index)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                if derivatives[neighbor].typ[selected]
                                    & DERIV_RING_OUTSIDE_PRECURSOR as i16
                                    != 0
                                {
                                    num_cuts -= 2;
                                    num_ring_cuts -= 2;
                                } else {
                                    num_cuts -= 1;
                                }
                                num_cut_pieces -= 1;
                                let result =
                                    remove_deriv_mark(&mut derivatives[neighbor], second_index)?;
                                if result != 0 {
                                    return Ok(Err(result));
                                }
                                found += 1;
                            }
                        }
                    }
                    neighbor_derivative += 1;
                }
            }
            derivative_index += 1;
        }
    }
    if found != 0 {
        for derivative in derivatives.iter_mut().take(atom_count) {
            let mut index = 0_usize;
            while index < derivative.typ.len() && derivative.typ[index] != 0 {
                if derivative.typ[index] & DERIV_DUPLIC as i16 != 0 {
                    let result = remove_deriv(derivative, index as i32)?;
                    if result != 0 {
                        return Ok(Err(result));
                    }
                } else {
                    index += 1;
                }
            }
        }
    }
    Ok(Ok((num_cuts, num_ring_cuts, num_cut_pieces)))
}

fn oad_locate_cut_owner(
    atoms: &[inp_ATOM],
    derivatives: &[DERIV_AT],
    pair: &R2C_ATPAIR,
    require_atom_bound: Option<usize>,
) -> Result<Result<(usize, usize, usize), i32>, SourceHeapError> {
    let mut found = 0_i32;
    let mut pair_side = 0_usize;
    let mut owner = 0_usize;
    let mut derivative_index = 0_usize;
    for side in 0..2 {
        let candidate = usize::from(pair.at[side]);
        if require_atom_bound.is_some_and(|bound| candidate >= bound) {
            continue;
        }
        let candidate_derivative = derivatives
            .get(candidate)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for selected in 0..2 {
            if candidate_derivative.typ[selected] == 0 {
                continue;
            }
            let ordinal = usize::try_from(candidate_derivative.ord[selected])
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if atoms
                .get(candidate)
                .and_then(|atom| atom.neighbor.get(ordinal))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                == pair.at[1 - side]
            {
                pair_side = side;
                owner = candidate;
                derivative_index = selected;
                found += 1;
            }
        }
    }
    if found == 1 {
        Ok(Ok((pair_side, owner, derivative_index)))
    } else {
        Ok(Err(-3))
    }
}

fn oad_cut_span(type_: i16) -> usize {
    let mut span = 1_usize;
    if type_ & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
        return 2;
    }
    span += usize::from(type_ != 0 && type_ & DERIV_RING_DMOX_DEOX as i16 == type_);
    span += usize::from(type_ != 0 && type_ & DERIV_RING2_OUTSIDE_PRECUR as i16 == type_);
    span
}

fn oad_validate_derivative_agents(
    atoms: &mut [inp_ATOM],
    derivatives: &mut [DERIV_AT],
    num_atoms: i32,
    atom_pairs: &mut [R2C_ATPAIR],
    num_cuts: &mut i32,
    num_ring_cuts: &mut i32,
    num_cut_pieces: &mut i32,
    num_cuts_to_check: &mut i32,
    current_num_atoms: &mut i32,
) -> Result<i32, SourceHeapError> {
    if *num_cuts_to_check < 2 {
        return Ok(0);
    }
    loop {
        let expected =
            usize::try_from(*num_cuts_to_check).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if atom_pairs.len() < expected {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let filled = fill_out_bond_cuts(
            atoms,
            derivatives,
            num_atoms,
            atom_pairs,
            *num_cuts_to_check,
        )?;
        if filled < 0 {
            UnMarkOtherIndicators(atoms, num_atoms)?;
            return Ok(filled);
        }
        if filled != *num_cuts_to_check {
            UnMarkOtherIndicators(atoms, num_atoms)?;
            return Ok(-3);
        }

        let mut components = 0_u16;
        let mut index = 0_usize;
        while index < expected {
            let (precursor_side, owner, derivative_index) =
                match oad_locate_cut_owner(atoms, derivatives, &atom_pairs[index], None)? {
                    Ok(value) => value,
                    Err(value) => {
                        UnMarkOtherIndicators(atoms, num_atoms)?;
                        return Ok(value);
                    }
                };
            let type_ = derivatives[owner].typ[derivative_index];
            let span = oad_cut_span(type_);
            let start = if type_ & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
                usize::from(atom_pairs[index].at[precursor_side])
            } else {
                usize::from(atom_pairs[index].at[1 - precursor_side])
            };
            if index + span > expected
                || (span == 2
                    && usize::from(atom_pairs[index + 1].at[0]) != start
                    && usize::from(atom_pairs[index + 1].at[1]) != start)
            {
                UnMarkOtherIndicators(atoms, num_atoms)?;
                return Ok(-3);
            }
            *current_num_atoms = mark_atoms_ap(
                atoms,
                start as u16,
                &atom_pairs[index..],
                span as i32,
                0,
                CFLAG_MARK_BRANCH as u16,
            )?;
            for candidate in 0..expected {
                if candidate == index || candidate == index + span - 1 {
                    continue;
                }
                let pair = &atom_pairs[candidate];
                if atoms[usize::from(pair.at[0])].at_type != 0
                    || atoms[usize::from(pair.at[1])].at_type != 0
                {
                    if type_ & DERIV_UNEXPADABLE as i16 == type_ {
                        let candidate_owner = usize::from(pair.atno);
                        if candidate_owner != owner {
                            let mut marked = 0_usize;
                            while marked < 2 && derivatives[candidate_owner].typ[marked] != 0 {
                                derivatives[candidate_owner].typ[marked] |= DERIV_UNMARK as i16;
                                marked += 1;
                            }
                            *num_cuts -= 1;
                            *num_cut_pieces -= 1;
                            if marked == 2 {
                                *num_cuts -= 1;
                                *num_ring_cuts -= 2;
                            }
                            components = components.wrapping_add(1);
                        }
                    } else {
                        derivatives[owner].typ[derivative_index] |= DERIV_UNMARK as i16;
                        *num_cuts -= 1;
                        *num_cut_pieces -= 1;
                        if span == 2 {
                            derivatives[owner].typ[1 - derivative_index] |= DERIV_UNMARK as i16;
                            *num_cuts -= 1;
                            *num_ring_cuts -= 2;
                        }
                        components = components.wrapping_add(1);
                    }
                    break;
                }
            }
            UnMarkOtherIndicators(atoms, num_atoms)?;
            index += span;
        }

        if components != 0 {
            for derivative in derivatives
                .iter_mut()
                .take(usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            {
                let second = if derivative.typ[0] & DERIV_UNMARK as i16 != 0 {
                    derivative.num[0] = derivative.num[1];
                    derivative.ord[0] = derivative.ord[1];
                    derivative.typ[0] = derivative.typ[1];
                    derivative.typ[1] = 0;
                    0
                } else {
                    1
                };
                if derivative.typ[second] & DERIV_UNMARK as i16 != 0 {
                    derivative.typ[second] = 0;
                }
            }
            *num_cuts_to_check = *num_cuts;
            if *num_cuts < 0 || *num_ring_cuts < 0 || *num_cut_pieces < 0 {
                UnMarkOtherIndicators(atoms, num_atoms)?;
                return Ok(-3);
            }
            if *num_cuts_to_check > 0 {
                continue;
            }
        }

        components = 0;
        let result = mark_deriv_agents(
            atoms,
            derivatives,
            num_atoms,
            atom_pairs,
            *num_cuts_to_check,
            &mut components,
            current_num_atoms,
        )?;
        if result != 0 {
            UnMarkOtherIndicators(atoms, num_atoms)?;
            return Ok(result);
        }
        if components > 1 {
            if *num_ring_cuts <= 2 {
                UnMarkOtherIndicators(atoms, num_atoms)?;
                return Ok(-99);
            }
            let mut eliminated = 0_i32;
            for atom_index in
                0..usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?
            {
                if derivatives[atom_index].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0
                    && derivatives[atom_index].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0
                {
                    let first = usize::from(
                        atoms[atom_index].neighbor[usize::try_from(derivatives[atom_index].ord[0])
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                    );
                    let second = usize::from(
                        atoms[atom_index].neighbor[usize::try_from(derivatives[atom_index].ord[1])
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                    );
                    if atoms[first].at_type != atoms[second].at_type {
                        derivatives[atom_index].typ[0] = 0;
                        derivatives[atom_index].typ[1] = 0;
                        eliminated += 1;
                        *num_cuts_to_check -= 2;
                        *num_cuts -= 2;
                        *num_ring_cuts -= 2;
                        *num_cut_pieces -= 1;
                    }
                }
            }
            if eliminated > 0 && *num_cuts_to_check > 2 {
                UnMarkOtherIndicators(atoms, num_atoms)?;
                continue;
            }
        }
        UnMarkOtherIndicators(atoms, num_atoms)?;
        return Ok(0);
    }
}

fn oad_check_final_precursor(
    atoms: &mut [inp_ATOM],
    derivatives: &[DERIV_AT],
    num_atoms: i32,
    atom_pairs: &mut [R2C_ATPAIR],
    num_cuts_to_check: i32,
    current_num_atoms: &mut i32,
) -> Result<i32, SourceHeapError> {
    let expected =
        usize::try_from(num_cuts_to_check).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let filled = fill_out_bond_cuts(atoms, derivatives, num_atoms, atom_pairs, num_cuts_to_check)?;
    if filled < 0 {
        return Ok(filled);
    }
    if filled != num_cuts_to_check {
        return Ok(-3);
    }
    let atom_count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut index = 0_usize;
    while index < expected {
        let (precursor_side, owner, derivative_index) =
            match oad_locate_cut_owner(atoms, derivatives, &atom_pairs[index], Some(atom_count))? {
                Ok(value) => value,
                Err(value) => return Ok(value),
            };
        let type_ = derivatives[owner].typ[derivative_index];
        let span = oad_cut_span(type_);
        let start = if type_ & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
            usize::from(atom_pairs[index].at[precursor_side])
        } else {
            usize::from(atom_pairs[index].at[1 - precursor_side])
        };
        if index + span > expected
            || (span == 2
                && usize::from(atom_pairs[index + 1].at[0]) != start
                && usize::from(atom_pairs[index + 1].at[1]) != start)
        {
            return Ok(-3);
        }
        *current_num_atoms = mark_atoms_ap(
            atoms,
            start as u16,
            &atom_pairs[index..],
            span as i32,
            0,
            CFLAG_MARK_BRANCH as u16,
        )?;
        UnMarkOtherIndicators(atoms, num_atoms)?;
        index += span;
    }
    let mut components = 0_u16;
    let result = mark_deriv_agents(
        atoms,
        derivatives,
        num_atoms,
        atom_pairs,
        num_cuts_to_check,
        &mut components,
        current_num_atoms,
    )?;
    UnMarkOtherIndicators(atoms, num_atoms)?;
    Ok(result)
}

fn oad_make_selected_cuts(
    atoms: &mut [inp_ATOM],
    derivatives: &mut [DERIV_AT],
    num_atoms: i32,
) -> Result<Result<i32, i32>, SourceHeapError> {
    let mut num_cuts = 0_i32;
    for atom_index in
        0..usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?
    {
        let length = derivatives[atom_index]
            .typ
            .iter()
            .take_while(|type_| **type_ != 0)
            .count();
        match length {
            0 => continue,
            1 => {
                make_single_cut(atoms, derivatives, atom_index as i32, 0)?;
                num_cuts += 1;
            }
            2 => {
                let first = derivatives[atom_index].typ[0];
                let second = derivatives[atom_index].typ[1];
                if (first & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0
                    && second & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0)
                    || (first == DERIV_AMINE_tN as i16 && second == DERIV_AMINE_tN as i16)
                    || (first != 0
                        && first & DERIV_RING2_OUTSIDE_PRECUR as i16 == first
                        && second == first)
                {
                    make_single_cut(atoms, derivatives, atom_index as i32, 1)?;
                    make_single_cut(atoms, derivatives, atom_index as i32, 0)?;
                    num_cuts += 1;
                } else if first == second {
                    if derivatives[atom_index].num[0] > derivatives[atom_index].num[1] {
                        make_single_cut(atoms, derivatives, atom_index as i32, 1)?;
                        num_cuts += 1;
                    } else if derivatives[atom_index].num[0] < derivatives[atom_index].num[1] {
                        make_single_cut(atoms, derivatives, atom_index as i32, 0)?;
                        num_cuts += 1;
                    }
                } else {
                    return Ok(Err(-88));
                }
            }
            3 => {
                if derivatives[atom_index].typ[0] != DERIV_AMINE_tN as i16
                    || derivatives[atom_index].typ[1] != derivatives[atom_index].typ[0]
                    || derivatives[atom_index].typ[2] != derivatives[atom_index].typ[0]
                {
                    return Ok(Err(-88));
                }
                let mut smallest =
                    if derivatives[atom_index].num[0] < derivatives[atom_index].num[1] {
                        0
                    } else {
                        1
                    };
                smallest = if derivatives[atom_index].num[smallest] < derivatives[atom_index].num[2]
                {
                    smallest
                } else {
                    2
                };
                let mut largest = if derivatives[atom_index].num[0] < derivatives[atom_index].num[1]
                {
                    1
                } else {
                    0
                };
                largest = if derivatives[atom_index].num[smallest] < derivatives[atom_index].num[2]
                {
                    2
                } else {
                    largest
                };
                let mut middle = ((smallest + 1) ^ (largest + 1)) - 1;
                if derivatives[atom_index].num[smallest] == derivatives[atom_index].num[largest] {
                    continue;
                }
                if derivatives[atom_index].num[smallest] == derivatives[atom_index].num[middle]
                    && smallest < middle
                {
                    std::mem::swap(&mut smallest, &mut middle);
                }
                make_single_cut(atoms, derivatives, atom_index as i32, smallest as i32)?;
                num_cuts += 1;
                if derivatives[atom_index].num[smallest] == derivatives[atom_index].num[middle] {
                    make_single_cut(atoms, derivatives, atom_index as i32, middle as i32)?;
                    num_cuts += 1;
                }
            }
            4 => {
                if !derivatives[atom_index]
                    .typ
                    .iter()
                    .all(|type_| *type_ & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0)
                {
                    return Ok(Err(-88));
                }
                let first = derivatives[atom_index].num[0].max(derivatives[atom_index].num[1]);
                let second = derivatives[atom_index].num[2].max(derivatives[atom_index].num[3]);
                if first < second {
                    make_single_cut(atoms, derivatives, atom_index as i32, 1)?;
                    make_single_cut(atoms, derivatives, atom_index as i32, 0)?;
                    num_cuts += 1;
                } else if first > second {
                    make_single_cut(atoms, derivatives, atom_index as i32, 3)?;
                    make_single_cut(atoms, derivatives, atom_index as i32, 2)?;
                    num_cuts += 1;
                }
            }
            _ => return Ok(Err(-88)),
        }
    }
    Ok(Ok(num_cuts))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_Edit_Underivatize(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canon_globals: &mut CANON_GLOBALS,
    original: &mut ORIG_ATOM_DATA,
    output_sdf: i32,
    output_report: i32,
    sdf_value: SourceMutPointer<i8>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5934 OAD_Edit_Underivatize
    // INCHI✔️❌: int OAD_Edit_Underivatize( struct tagINCHI_CLOCK *ic,
    // INCHI✔️❌:                            struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                            ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                            int bOutputSdf,
    // INCHI✔️❌:                            int bOutputReport,
    // INCHI✔️❌:                            char *pSdfValue )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌: #define ALLOC_AP \
    // INCHI✔️❌:         if ( 0 < num_cuts_to_check && (lenAllocated_ap < num_cuts_to_check || !ap) ) {\
    // INCHI✔️❌:             if ( ap )\
    // INCHI✔️❌:                 inchi_free( ap );\
    // INCHI✔️❌:             ap = (R2C_ATPAIR *) inchi_malloc( num_cuts_to_check * sizeof(ap[0]) );\
    // INCHI✔️❌:             if ( !ap ) {\
    // INCHI✔️❌:                 ret = -1;\
    // INCHI✔️❌:                 goto exit_function; /* malloc failure */\
    // INCHI✔️❌:             }\
    // INCHI✔️❌:             lenAllocated_ap = num_cuts_to_check;\
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     int ret = 0, i, j, k, m, n, num_atoms, num_components, i_component, nFound, num, cur_num_at = 0, len, ind1, ind2, ind3; /* djb-rwth: adding variables for char -> int conversion of subscripts; initialisation added */
    // INCHI✔️❌:     int num_cuts, num_ring_cuts, num_cut_pieces, num_cuts_to_check;
    // INCHI✔️❌:     inp_ATOM *at = orig_inp_data->at; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔️❌:     INP_ATOM_DATA *inp_cur_data = NULL;
    // INCHI✔️❌:     DERIV_AT      *da = NULL;
    // INCHI✔️❌:     R2C_ATPAIR    *ap = NULL;
    // INCHI✔️❌:     int            lenAllocated_ap = 0;
    // INCHI✔️❌:     int  nTotNumCuts = 0;
    // INCHI✔️❌:     int  num_removed_H = 0; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌: #ifdef FIX_UNDERIV_TO_SDF
    // INCHI✔️❌:     inp_ATOM *at2 = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if ( UNDERIVATIZE_REPORT == 1 )
    // INCHI✔️❌: #define UNDERIV_LIST_LEN 2048
    // INCHI✔️❌: #define UNDERIV_LIST_LEN2 2048
    // INCHI✔️❌:     char szUnderivList[UNDERIV_LIST_LEN] = "";
    // INCHI✔️❌:     char szUnderivList2[UNDERIV_LIST_LEN2] = "";
    // INCHI✔️❌:     char underivPrefix[] = "\tDeriv=";
    // INCHI✔️❌:     char underivPostfix[] = "";
    // INCHI✔️❌:     char underivPrefix2[] = "\tDeriv2=";
    // INCHI✔️❌:     char underivPostfix2[] = "";
    // INCHI✔️❌:     char underivPrefix3[] = "\tDerivBits=";
    // INCHI✔️❌:     char underivPostfix3[] = "";
    // INCHI✔️❌:     char cDerivSeparator = ',';
    // INCHI✔️❌:     int numUnderiv = 0, numUnderiv2 = 0, numUnderiv3 = 0; /* djb-rwth: ignoring LLVM warning: variables used to store function return value */
    // INCHI✔️❌:     BIT_UNDERIV bitUnderivList = 0;
    // INCHI✔️❌:     char szbitUnderivList[16]; /* int32 has at most 32/4=8 hexadecimal digits + 0x prefix + zero termination = 8+2+1=11 */
    // INCHI✔️❌: #else
    // INCHI✔️❌: #define UNDERIV_LIST_LEN 0
    // INCHI✔️❌: #define UNDERIV_LIST_LEN2 0
    // INCHI✔️❌: #define UNDERIV_MAX_NUM  0   /*max. number of records in szUnderivList */
    // INCHI✔️❌:     char *szUnderivList = NULL;
    // INCHI✔️❌:     char *szUnderivList2 = NULL;
    // INCHI✔️❌:     char **pszUnderiv = NULL;
    // INCHI✔️❌:     BIT_UNDERIV bitUnderivList = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Prepare */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*set_R2C_el_numbers( );*/
    // INCHI✔️❌:
    // INCHI✔️❌: #ifndef UNDERIV_ADD_EXPLICIT_H
    // INCHI✔️❌:     num_atoms = remove_terminal_HDT( orig_inp_data->num_inp_atoms, at, 1 );
    // INCHI✔️❌:     /*^^^^^ always accomodate accomodate FIX_TERM_H_CHRG_BUG - IPl, July 2008*/
    // INCHI✔️❌:     num_removed_H = orig_inp_data->num_inp_atoms - num_atoms;
    // INCHI✔️❌:     orig_inp_data->num_inp_atoms = num_atoms;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Initialize */
    // INCHI✔️❌:     UnMarkDisconnectedComponents( orig_inp_data );
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Mark */
    // INCHI✔️❌:     num_components = MarkDisconnectedComponents( orig_inp_data, 0 );
    // INCHI✔️❌:     inp_cur_data = (INP_ATOM_DATA *) inchi_calloc( num_components, sizeof( inp_cur_data[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i_component = 0; i_component < num_components; i_component++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         CreateInpAtomData(inp_cur_data + i_component, orig_inp_data->nCurAtLen[i_component], 0);
    // INCHI✔️❌:
    // INCHI✔️❌:         inp_cur_data[i_component].num_at = ExtractConnectedComponent(orig_inp_data->at, orig_inp_data->num_inp_atoms, i_component + 1, inp_cur_data[i_component].at);
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  error processing */
    // INCHI✔️❌:         if (inp_cur_data[i_component].num_at <= 0 || orig_inp_data->nCurAtLen[i_component] != inp_cur_data[i_component].num_at)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = -(i_component + 1); /* severe error */
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef UNDERIV_ADD_EXPLICIT_H
    // INCHI✔️❌:         num_atoms = remove_terminal_HDT(inp_cur_data[i_component].num_at, inp_cur_data[i_component].at, 1);
    // INCHI✔️❌:         inp_cur_data[i_component].num_removed_H = inp_cur_data[i_component].num_at - num_atoms;
    // INCHI✔️❌:         inp_cur_data[i_component].num_at = num_atoms;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Initialize */
    // INCHI✔️❌:         num_atoms = inp_cur_data[i_component].num_at;
    // INCHI✔️❌:         at = inp_cur_data[i_component].at;
    // INCHI✔️❌:         add_DT_to_num_H(num_atoms, at);
    // INCHI✔️❌:
    // INCHI✔️❌:         UnMarkRingSystemsInp(at, num_atoms);
    // INCHI✔️❌:         UnMarkOtherIndicators(at, num_atoms);
    // INCHI✔️❌:         UnMarkOneComponent(at, num_atoms);
    // INCHI✔️❌:         MarkRingSystemsInp(at, num_atoms, 0);
    // INCHI✔️❌: #ifdef FIX_UNDERIV_TO_SDF
    // INCHI✔️❌:         if (bOutputSdf)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* save orig. bond types to restore them after replacing them with aromatic */
    // INCHI✔️❌:             if (at2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free(at2);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((at2 = (inp_ATOM*)inchi_malloc(num_atoms * sizeof(inp_ATOM)))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memcpy(at2, at, num_atoms * sizeof(inp_ATOM));
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Mark aromatic bonds */
    // INCHI✔️❌:         ret = mark_arom_bonds(ic, pCG, at, num_atoms);
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (da)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free(da);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         da = (DERIV_AT*)inchi_calloc(num_atoms, sizeof(da[0]));
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Detect derivatives */
    // INCHI✔️❌:         nFound = 0;
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[i].bCutVertex && !da[i].typ[inchi_min(at[i].valence, DERIV_AT_LEN) - 1])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (k = 0; k < at[i].valence; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num = count_one_bond_atoms(at, da, i, k, CFLAG_MARK_BRANCH, &nFound);
    // INCHI✔️❌:                     UnMarkOtherIndicators(at, num_atoms);
    // INCHI✔️❌:                     if (num < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = num; /* severe error */
    // INCHI✔️❌:                         goto exit_function;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Prepare cuts: remove cuts that are not to be done */
    // INCHI✔️❌:         /* in addition, count ring cuts DERIV_RING_OUTSIDE_PRECURSOR */
    // INCHI✔️❌:         num_ring_cuts = 0;
    // INCHI✔️❌:         num_cuts = 0;
    // INCHI✔️❌:         num_cut_pieces = 0;
    // INCHI✔️❌:         if (da) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = num = 0; i < num_atoms; i++) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*for ( len = 0; len < MAX_AT_DERIV && da[i].typ[len]; len ++ ) -- bug fixed 2013-11-07 DCh */
    // INCHI✔️❌:                 for (len = 0; len < DERIV_AT_LEN && da[i].typ[len]; len++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 switch (len)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:
    // INCHI✔️❌:                 case 0:
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                 case 1:
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:                     if (da[i].typ[0] & DERIV_RING_DMOX_DEOX)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (!da[i].other_atom || (j = da[i].other_atom - 1) >= num_atoms ||
    // INCHI✔️❌:                             (da[j].typ[0] ^ DERIV_RING_DMOX_DEOX) != da[i].typ[0] ||
    // INCHI✔️❌:                             da[j].other_atom - 1 != i)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* program error */
    // INCHI✔️❌:                             da[i].num[0] = NOT_AT_DERIV;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* one of 2 ring cuts */
    // INCHI✔️❌:                             num_cuts += 1;
    // INCHI✔️❌:                             num_ring_cuts += 1;
    // INCHI✔️❌:                             num_cut_pieces += (i > j); /* count pieces only once */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* single cut: unconditional */
    // INCHI✔️❌:                         num_cuts += 1;
    // INCHI✔️❌:                         num_cut_pieces += 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:
    // INCHI✔️❌:                 case 2:
    // INCHI✔️❌:                     if ((da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR) && (da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* double cut, unconditional */
    // INCHI✔️❌:                         num_ring_cuts += 2;
    // INCHI✔️❌:                         num_cuts += 2;
    // INCHI✔️❌:                         num_cut_pieces += 1;
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌: #ifdef DERIV_RING2_OUTSIDE_PRECUR
    // INCHI✔️❌:                         if (da[i].typ[0] && (da[i].typ[0] & DERIV_RING2_OUTSIDE_PRECUR) == da[i].typ[0] && da[i].typ[1] == da[i].typ[0])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* double cut, unconditional */
    // INCHI✔️❌:                             num_ring_cuts += 2;
    // INCHI✔️❌:                             num_cuts += 2;
    // INCHI✔️❌:                             num_cut_pieces += 1;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                             if (da[i].typ[0] == DERIV_AMINE_tN && da[i].typ[1] == DERIV_AMINE_tN)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* double cut, unconditional */
    // INCHI✔️❌:                                 num_cuts += 2;
    // INCHI✔️❌:                                 num_cut_pieces += 2;
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌: #ifdef DERIV_RO_COX
    // INCHI✔️❌:                                 if (da[i].typ[0] == DERIV_RO_COX && da[i].typ[1] == DERIV_RO_COX)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (da[i].num[0] == da[i].num[1])
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         memset(da + i, 0, sizeof(da[0])); /* don't remove if the two agents are identical */ /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         if (da[i].num[0] && da[i].num[1])
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             static char pref_RO_COX[] = {/*likely deriv.agent*/12, 9, 6, 13, 3, 8/*likely precursor*/, 0 };
    // INCHI✔️❌:                                             char* p0 = strchr(pref_RO_COX, da[i].num[0]);
    // INCHI✔️❌:                                             char* p1 = strchr(pref_RO_COX, da[i].num[1]);
    // INCHI✔️❌:                                             if (p0 && p1)
    // INCHI✔️❌:                                             {
    // INCHI✔️❌:                                                 j = p1 < p0; /* j=1 => deriv. agent num[1] has higher priority */
    // INCHI✔️❌:                                             }
    // INCHI✔️❌:                                             else
    // INCHI✔️❌:                                             {
    // INCHI✔️❌:                                                 j = p0 ? 0 : p1 ? 1 : -1; /* we are here if there is a program error */
    // INCHI✔️❌:                                             }
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         else
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             j = da[i].num[0] ? 0 : da[i].num[1] ? 1 : -1; /* we are here if there is a program error */
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                                     switch (j)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                     case 1:
    // INCHI✔️❌:                                         da[i].num[0] = da[i].num[1];
    // INCHI✔️❌:                                         da[i].ord[0] = da[i].ord[1];
    // INCHI✔️❌:                                         /* fall through */
    // INCHI✔️❌:                                     case 0:
    // INCHI✔️❌:                                         da[i].typ[1] = 0;
    // INCHI✔️❌:                                         da[i].num[1] = 0;
    // INCHI✔️❌:                                         da[i].ord[1] = 0;
    // INCHI✔️❌:                                         num_cuts += 1;
    // INCHI✔️❌:                                         num_cut_pieces += 1;
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:                                     case -1:
    // INCHI✔️❌:                                         memset(da + i, 0, sizeof(da[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                                         break; /* will produce error */
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌: #endif  /* DERIV_RO_COX */
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (da[i].typ[0] == da[i].typ[1])
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         int sy0 = 0, sy1 = 0;
    // INCHI✔️❌:                                         /* DERIV_BRIDGE_O or DERIV_BRIDGE_NH ; cut off the smallest */
    // INCHI✔️❌:                                         if (0 == is_deriv_chain(at, i, num_atoms, da + i, 0, NULL, 0, NULL, 0, NULL))
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             da[i].num[0] = NOT_AT_DERIV;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         else
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             ind1 = da[i].ord[0] - '0'; /* djb-rwth: converting char to int for subscript use */
    // INCHI✔️❌:                                             sy0 = is_silyl2(at, at[i].neighbor[ind1], i);
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         if (0 == is_deriv_chain(at, i, num_atoms, da + i, 1, NULL, 0, NULL, 0, NULL))
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             da[i].num[1] = NOT_AT_DERIV;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         else
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             ind2 = da[i].ord[1] - '0'; /* djb-rwth: converting char to int for subscript use */
    // INCHI✔️❌:                                             sy1 = is_silyl2(at, at[i].neighbor[ind2], i);
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         if ((sy1 && (!sy0 || sy1 < sy0)) || (!(sy0 || sy1) && da[i].num[0] > da[i].num[1])) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             da[i].num[0] = da[i].num[1];
    // INCHI✔️❌:                                             da[i].ord[0] = da[i].ord[1];
    // INCHI✔️❌:                                             da[i].typ[0] = da[i].typ[1];
    // INCHI✔️❌:                                             da[i].typ[1] = 0;
    // INCHI✔️❌:                                             num_cuts += 1;
    // INCHI✔️❌:                                             num_cut_pieces += 1;
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         else
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             if ((sy0 && (!sy1 || sy0 < sy1)) || (!(sy0 || sy1) && da[i].num[0] < da[i].num[1])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                                             {
    // INCHI✔️❌:                                                 da[i].typ[1] = 0;
    // INCHI✔️❌:                                                 num_cuts += 1;
    // INCHI✔️❌:                                                 num_cut_pieces += 1;
    // INCHI✔️❌:                                             }
    // INCHI✔️❌:                                             else
    // INCHI✔️❌:                                             {
    // INCHI✔️❌:                                                 /* attachments have same size: ignore both */
    // INCHI✔️❌:                                                 /* ??? check for standard derivatizations ??? */
    // INCHI✔️❌:                                                 da[i].typ[0] = 0;
    // INCHI✔️❌:                                                 da[i].typ[1] = 0;
    // INCHI✔️❌:                                             }
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                         continue;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌: #ifdef DERIV_RO_COX
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         if ((da[i].typ[0] & (DERIV_RO_COX | DERIV_BRIDGE_O)) && (da[i].typ[1] & (DERIV_RO_COX | DERIV_BRIDGE_O)))
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             /*static char pref_RO_COX2[] = {13, 3, 8*/
    // INCHI✔️❌:                                             /*likely precursor*/
    // INCHI✔️❌:                                             /*, 0};*/
    // INCHI✔️❌:                                             j = (da[i].typ[0] == DERIV_BRIDGE_O);  /* da[i].typ[j] == DERIV_RO_COX */
    // INCHI✔️❌:                                             /*------------------------------------------------------------------------*
    // INCHI✔️❌:                                             * discard DERIV_RO_COX only in case [CH3-C(=O)]-O-[CH3]                  *
    // INCHI✔️❌:                                             *                      precursor DERIV_BRIDGE_O   DERIV_RO_COX precursor *
    // INCHI✔️❌:                                             * has already been done in finding DERIV_RO_COX                          *
    // INCHI✔️❌:                                             *------------------------------------------------------------------------*/
    // INCHI✔️❌: #ifdef NEVER
    // INCHI✔️❌:                                             /* methyl/ethyl alcoholes have already been excluded in get_traversed_deriv_type(...)
    // INCHI✔️❌:                                             for etyl/methyl acetate/benzoate. Only acetate/benzoic acids may their precursors */
    // INCHI✔️❌:                                             if (da[i].num[j] == 3 &&  /* -O-[C(=O)-CH3]  : DERIV_RO_COX*/
    // INCHI✔️❌:                                                 da[i].num[1 - j] == 1)
    // INCHI✔️❌:                                             { /* [CH3]-O-C(=O)-R : DERIV_BRIDGE_O; derivatizin agent in [] */
    // INCHI✔️❌:                                                 j = 1 - j; /* remove da[i].xxx[j] */
    // INCHI✔️❌:                                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                             /*
    // INCHI✔️❌:                                             O
    // INCHI✔️❌:                                             ||
    // INCHI✔️❌:                                             R----O----C---X
    // INCHI✔️❌:
    // INCHI✔️❌:                                             R----O---OC---X
    // INCHI✔️❌:                                             DERIV_BRIDGE_O \___/\____precursor__/
    // INCHI✔️❌:
    // INCHI✔️❌:                                             R----O----COX
    // INCHI✔️❌:                                             \___precursor__/ \___/ DERIV_RO_COX
    // INCHI✔️❌:
    // INCHI✔️❌:                                             rule: If R is >Si< then select DERIV_BRIDGE_O
    // INCHI✔️❌:
    // INCHI✔️❌:                                             */
    // INCHI✔️❌:
    // INCHI✔️❌:                                             if (j /* choose DERIV_RO_COX */)
    // INCHI✔️❌:                                             {
    // INCHI✔️❌:                                                 /* da[i].typ[1]=DERIV_RO_COX is not likely a precursor */
    // INCHI✔️❌:                                                 da[i].typ[0] = da[i].typ[1];
    // INCHI✔️❌:                                                 da[i].num[0] = da[i].num[1];
    // INCHI✔️❌:                                                 da[i].ord[0] = da[i].ord[1];
    // INCHI✔️❌:                                             }
    // INCHI✔️❌:                                             da[i].typ[1] = 0;
    // INCHI✔️❌:                                             da[i].num[1] = 0;
    // INCHI✔️❌:                                             da[i].ord[1] = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                                             num_cuts += 1;
    // INCHI✔️❌:                                             num_cut_pieces += 1;
    // INCHI✔️❌:                                             continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     ret = -88;
    // INCHI✔️❌:                     goto exit_function; /* unexpected */
    // INCHI✔️❌:
    // INCHI✔️❌:                 case 3:
    // INCHI✔️❌:                     if (da[i].typ[0] == da[i].typ[1] &&
    // INCHI✔️❌:                         da[i].typ[0] == da[i].typ[2] &&
    // INCHI✔️❌:                         da[i].typ[0] == DERIV_AMINE_tN)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int x, y, z;
    // INCHI✔️❌:                         int sy[3] = { 0, 0, 0 }; /* silyl */
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (0 == is_deriv_chain(at, i, num_atoms, da + i, 0, NULL, 0, NULL, 0, NULL))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             da[i].num[0] = NOT_AT_DERIV;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ind1 = da[i].ord[0] - '0'; /* djb-rwth: converting char to int for subscript use */
    // INCHI✔️❌:                             sy[0] = is_silyl2(at, at[i].neighbor[ind1], i);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (0 == is_deriv_chain(at, i, num_atoms, da + i, 1, NULL, 0, NULL, 0, NULL))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             da[i].num[1] = NOT_AT_DERIV;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ind2 = da[i].ord[1] - '0'; /* djb-rwth: converting char to int for subscript use */
    // INCHI✔️❌:                             sy[1] = is_silyl2(at, at[i].neighbor[ind2], i);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (0 == is_deriv_chain(at, i, num_atoms, da + i, 2, NULL, 0, NULL, 0, NULL))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             da[i].num[2] = NOT_AT_DERIV;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ind3 = da[i].ord[2] - '0'; /* djb-rwth: converting char to int for subscript use */
    // INCHI✔️❌:                             sy[2] = is_silyl2(at, at[i].neighbor[ind3], i);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:                         x = ((sy[0] && (!sy[1] || sy[0] < sy[1])) || (!(sy[0] || sy[1]) && da[i].num[0] < da[i].num[1])) ? 0 : 1; /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         z = !x;
    // INCHI✔️❌:                         x = ((sy[x] && (!sy[2] || sy[x] < sy[2])) || (!(sy[x] || sy[2]) && da[i].num[x] < da[i].num[2])) ? x : 2; /* min */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         /*z = (da[i].num[0] < da[i].num[1])? 1 : 0;*/
    // INCHI✔️❌:                         /*z = (da[i].num[x] < da[i].num[2])? 2 : z;*/ /* max */
    // INCHI✔️❌:                         z = ((sy[x] && (!sy[2] || sy[x] < sy[2])) || (!(sy[x] || sy[2]) && da[i].num[x] < da[i].num[2])) ? 2 : z; /* max */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         y = ((x + 1) ^ (z + 1)) - 1;                      /* median */
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (da[i].num[x] == da[i].num[z] && sy[x] == sy[z])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* no cuts */
    // INCHI✔️❌:                             da[i].typ[0] = 0;
    // INCHI✔️❌:                             da[i].typ[1] = 0;
    // INCHI✔️❌:                             da[i].typ[2] = 0;
    // INCHI✔️❌:                             continue; /* all deriv. agents have same size, ignore */
    // INCHI✔️❌:                             /* ??? check for standard derivatizations ??? */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if ((da[i].num[x] == da[i].num[y] && sy[x] == sy[y]) || (sy[x] && sy[y] && !sy[z])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* two smallest */
    // INCHI✔️❌:                                 switch (z)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                 case 0:
    // INCHI✔️❌:                                     da[i].num[0] = da[i].num[1];
    // INCHI✔️❌:                                     da[i].ord[0] = da[i].ord[1];
    // INCHI✔️❌:                                     da[i].typ[0] = da[i].typ[1];
    // INCHI✔️❌:
    // INCHI✔️❌:                                     da[i].num[1] = da[i].num[2];
    // INCHI✔️❌:                                     da[i].ord[1] = da[i].ord[2];
    // INCHI✔️❌:                                     da[i].typ[1] = da[i].typ[2];
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 case 1:
    // INCHI✔️❌:                                     da[i].num[1] = da[i].num[2];
    // INCHI✔️❌:                                     da[i].ord[1] = da[i].ord[2];
    // INCHI✔️❌:                                     da[i].typ[1] = da[i].typ[2];
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 case 2:
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 da[i].typ[2] = 0;
    // INCHI✔️❌:                                 num_cuts += 2;
    // INCHI✔️❌:                                 num_cut_pieces += 2;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* one smallest */
    // INCHI✔️❌:                                 if (x)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     da[i].num[0] = da[i].num[x];
    // INCHI✔️❌:                                     da[i].ord[0] = da[i].ord[x];
    // INCHI✔️❌:                                     da[i].typ[0] = da[i].typ[x];
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 da[i].typ[1] = 0;
    // INCHI✔️❌:                                 da[i].typ[2] = 0;
    // INCHI✔️❌:                                 num_cuts += 1;
    // INCHI✔️❌:                                 num_cut_pieces += 1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ret = -88;
    // INCHI✔️❌:                     goto exit_function; /* unexpected */
    // INCHI✔️❌:
    // INCHI✔️❌:                 case 4:
    // INCHI✔️❌:                     if ((da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR) && (da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR) &&
    // INCHI✔️❌:                         (da[i].typ[2] & DERIV_RING_OUTSIDE_PRECURSOR) && (da[i].typ[3] & DERIV_RING_OUTSIDE_PRECURSOR))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int n01 = inchi_max(da[i].num[0], da[i].num[1]);
    // INCHI✔️❌:                         int n23 = inchi_max(da[i].num[2], da[i].num[3]);
    // INCHI✔️❌:                         if (n01 < n23 && 0 < is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(at, i, num_atoms, da + i, 0, NULL, 0, NULL, 0, NULL))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             da[i].typ[2] = 0;
    // INCHI✔️❌:                             da[i].typ[3] = 0;
    // INCHI✔️❌:                             num_cuts += 2;
    // INCHI✔️❌:                             num_ring_cuts += 2;
    // INCHI✔️❌:                             num_cut_pieces += 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (n01 > n23 && 0 < is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(at, i, num_atoms, da + i, 2, NULL, 0, NULL, 0, NULL))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 da[i].num[0] = da[i].num[2];
    // INCHI✔️❌:                                 da[i].ord[0] = da[i].ord[2];
    // INCHI✔️❌:                                 da[i].typ[0] = da[i].typ[2];
    // INCHI✔️❌:
    // INCHI✔️❌:                                 da[i].num[1] = da[i].num[3];
    // INCHI✔️❌:                                 da[i].ord[1] = da[i].ord[3];
    // INCHI✔️❌:                                 da[i].typ[1] = da[i].typ[3];
    // INCHI✔️❌:
    // INCHI✔️❌:                                 da[i].typ[2] = 0;
    // INCHI✔️❌:                                 da[i].typ[3] = 0;
    // INCHI✔️❌:                                 num_cuts += 2;
    // INCHI✔️❌:                                 num_ring_cuts += 2;
    // INCHI✔️❌:                                 num_cut_pieces += 1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 da[i].typ[0] = 0;
    // INCHI✔️❌:                                 da[i].typ[1] = 0;
    // INCHI✔️❌:                                 da[i].typ[2] = 0;
    // INCHI✔️❌:                                 da[i].typ[3] = 0;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ret = -88;
    // INCHI✔️❌:                     goto exit_function; /* unexpected */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             Eliminate cases when
    // INCHI✔️❌:             da[i1].typ[j1] && da[i2].typ[j2] &&
    // INCHI✔️❌:             at[i1].neighbor[da[i1].ord[j1]] == i2 && at[i2].neighbor[da[i2].ord[j2]] == i1
    // INCHI✔️❌:             that is, same bond is in the da[] elements of the adjacent atoms
    // INCHI✔️❌:             */
    // INCHI✔️❌:
    // INCHI✔️❌:             nFound = 0; /* number of cuts to remove */
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*for ( j = 0; j < MAX_AT_DERIV && da[i].typ[j]; j ++ ) -- bug fixed 2013-11-07 DCh */
    // INCHI✔️❌:                 for (j = 0; j < DERIV_AT_LEN && da[i].typ[j]; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (da[i].typ[j] & DERIV_DUPLIC)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     n = at[i].neighbor[(int)da[i].ord[j]];
    // INCHI✔️❌:                     if (n < i)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /*for ( k = 0; k < MAX_AT_DERIV && da[n].typ[k]; k ++ ) -- bug fixed 2013-11-07 DCh */
    // INCHI✔️❌:                     for (k = 0; k < DERIV_AT_LEN && n< num_atoms && da[n].typ[k]; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (da[n].typ[k] & DERIV_DUPLIC)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         m = at[n].neighbor[(int)da[n].ord[k]];
    // INCHI✔️❌:                         if (m == i)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* same bond in da[i].typ[j] and da[n].typ[k] */
    // INCHI✔️❌:                             /* check whether both derivatives are acceptable */
    // INCHI✔️❌:                             int k1 = k, j1 = j;
    // INCHI✔️❌:                             int ret_i = is_deriv_chain_or_ring(at, i, num_atoms, da + i, &j1);
    // INCHI✔️❌:                             int ret_n = is_deriv_chain_or_ring(at, n, num_atoms, da + n, &k1);
    // INCHI✔️❌:                             if (ret_i < 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ret = ret_i;
    // INCHI✔️❌:                                 goto exit_function;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (ret_n < 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ret = ret_n;
    // INCHI✔️❌:                                 goto exit_function;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!ret_i || (ret_i && ret_n)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (da[i].typ[j1] & DERIV_RING_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     num_cuts -= 2;
    // INCHI✔️❌:                                     num_ring_cuts -= 2;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     num_cuts -= 1;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 num_cut_pieces -= 1;
    // INCHI✔️❌:                                 if ((ret = remove_deriv_mark(da + i, j1))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     goto exit_function;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 nFound++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (!ret_n || (ret_i && ret_n)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (da[n].typ[k1] & DERIV_RING_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     num_cuts -= 2;
    // INCHI✔️❌:                                     num_ring_cuts -= 2;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     num_cuts -= 1;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 num_cut_pieces -= 1;
    // INCHI✔️❌:                                 if ((ret = remove_deriv_mark(da + n, k1))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     goto exit_function;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 nFound++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (nFound)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*for ( j = 0; j < MAX_AT_DERIV && da[i].typ[j]; j ++ ) -- bug fixed 2013-11-07 DCh */
    // INCHI✔️❌:                     for (j = 0; j < DERIV_AT_LEN && da[i].typ[j]; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* attn: j is changed inside the cycle body */
    // INCHI✔️❌:                         if (da[i].typ[j] & DERIV_DUPLIC)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if ((ret = remove_deriv(da + i, j))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 goto exit_function;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             j--;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* make sure DERIV_RING_OUTSIDE_PRECURSOR type disconnections increase */
    // INCHI✔️❌:         /* number of components by the number of disconnected derivateves */
    // INCHI✔️❌:         /* Avoid cases like these:
    // INCHI✔️❌:
    // INCHI✔️❌:         O--R--O             DO--R--OD
    // INCHI✔️❌:         /       \
    // INCHI✔️❌:         R--X         Y--R => R--XT2     T2Y--R
    // INCHI✔️❌:         \       /
    // INCHI✔️❌:         O--R--O             DO--R--OD
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         O--O                 DO--OD
    // INCHI✔️❌:         /    \
    // INCHI✔️❌:         R--X--O---Y--R    =>  R--X  OD2 Y--R
    // INCHI✔️❌:         T2     T2
    // INCHI✔️❌:
    // INCHI✔️❌:         */
    // INCHI✔️❌:         /* count DERIV_RING_OUTSIDE_PRECURSOR-type attachments */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( COUNT_ALL_NOT_DERIV == 1 )
    // INCHI✔️❌:         num_cuts_to_check = num_cuts; /* STOP HERE */
    // INCHI✔️❌: #else
    // INCHI✔️❌:         num_cuts_to_check = num_ring_cuts;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         ALLOC_AP
    // INCHI✔️❌:
    // INCHI✔️❌:             if (num_cuts_to_check >= 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* check */
    // INCHI✔️❌:                 AT_NUMB    comp_num;
    // INCHI✔️❌:                 int        /*n,*/ m_at, m_ord;
    // INCHI✔️❌:                 /*AT_NUMB at1, at2;*/
    // INCHI✔️❌:
    // INCHI✔️❌:             repeat_without_deriv_ring:
    // INCHI✔️❌:
    // INCHI✔️❌:                 ALLOC_AP
    // INCHI✔️❌:
    // INCHI✔️❌:                     comp_num = 0;
    // INCHI✔️❌:                 /* fill out the array of bonds to be cut */
    // INCHI✔️❌:                 j = fill_out_bond_cuts( at, da, num_atoms, ap, num_cuts_to_check );
    // INCHI✔️❌:                 if (j < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = j;
    // INCHI✔️❌:                     goto exit_r2c_num; /* wrong number of cuts = num */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (j != num_cuts_to_check)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = -3;
    // INCHI✔️❌:                     goto exit_r2c_num; /* wrong number of cuts = num */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* can't sort the bonds for subsequent searching by bisections in
    // INCHI✔️❌:                 mark_atoms_ap() -> has_atom_pair() 2013-08-48 DCh */
    // INCHI✔️❌:                 /* !!!!!!!! check that there are no derivatives inside a derivative */
    // INCHI✔️❌:                 comp_num = 0; /* here it is the number of removed cuts */
    // INCHI✔️❌:                 for (i = 0; i < num_cuts_to_check; i += j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (j = n = 0; j < 2; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int atj = (int) ap[i].at[j];
    // INCHI✔️❌:                         if (da[atj].typ[0] && at[atj].neighbor[(int) da[atj].ord[0]] == ap[i].at[1 - j])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = j;      /* ap[i].at[k] is precursor atom */
    // INCHI✔️❌:                             n++;
    // INCHI✔️❌:                             m_at = atj; /* precursor atom at[m_at], da[m_at] */
    // INCHI✔️❌:                             m_ord = 0;  /* da[m_at].typ[m_ord] - type of the deriv.bond to break  */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (da[atj].typ[1] && at[atj].neighbor[(int) da[atj].ord[1]] == ap[i].at[1 - j])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 k = j;
    // INCHI✔️❌:                                 n++;
    // INCHI✔️❌:                                 m_at = atj;
    // INCHI✔️❌:                                 m_ord = 1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (n != 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = -3;
    // INCHI✔️❌:                         goto exit_r2c_num; /* wrong atom pair */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (( da[m_at].typ[m_ord] & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = (int) ap[i].at[k];   /* atom inside the derivation attachment */
    // INCHI✔️❌:                         j = 2;             /* number of bonds to cut */
    // INCHI✔️❌:                         if (i + j > num_cuts_to_check || ((int) ap[i + 1].at[0] != n && (int) ap[i + 1].at[1] != n)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ret = -3;
    // INCHI✔️❌:                             goto exit_r2c_num; /* wrong atom pair */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = ap[i].at[1 - k]; /* atom inside the tentative derivation attachment */
    // INCHI✔️❌:                         j = 1;             /* number of bonds to cut */
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:                                            /*j += (0 != (da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX)); */
    // INCHI✔️❌:                                            /* these 2 cuts are always adjacent */
    // INCHI✔️❌:                         j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX ) ); /* these 2 cuts are always adjacent */
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_RING2_OUTSIDE_PRECUR
    // INCHI✔️❌:                         j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING2_OUTSIDE_PRECUR ) ); /* these 2 cuts are always adjacent */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* at[n] belongs to the derivation agent  */
    // INCHI✔️❌:                     cur_num_at = mark_atoms_ap( at, n, ap + i, j, 0, 1 );
    // INCHI✔️❌:                     for (k = 0; k < num_cuts_to_check; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (k == i || k == i + j - 1) /* not more than 2 cuts per derivatization agent */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue; /* skip current 1 or 2 cuts */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (at[(int) ap[k].at[0]].at_type || at[(int) ap[k].at[1]].at_type)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*#ifdef DERIV_X_OXIME*/
    // INCHI✔️❌:                             if (( da[m_at].typ[m_ord] & DERIV_UNEXPADABLE ) == da[m_at].typ[m_ord])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* this derivatization agent cannot be inside another, larger deriv. agent */
    // INCHI✔️❌:                                 /* disable cuts of the larger derivatization agent */
    // INCHI✔️❌:                                 if (ap[k].atno != m_at)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     for (j = 0; j < 2 && da[ap[k].atno].typ[j]; j++)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         da[ap[k].atno].typ[j] |= DERIV_UNMARK;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     num_cuts -= 1;
    // INCHI✔️❌:                                     num_cut_pieces -= 1;
    // INCHI✔️❌:                                     if (j == 2)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         num_cuts -= 1;
    // INCHI✔️❌:                                         num_ring_cuts -= 2;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                     comp_num++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 /* #endif */
    // INCHI✔️❌:                                 /* unmark the cut: found a cut inside the derivatizing agent */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 da[m_at].typ[m_ord] |= DERIV_UNMARK;
    // INCHI✔️❌:                                 num_cuts -= 1;
    // INCHI✔️❌:                                 num_cut_pieces -= 1;
    // INCHI✔️❌:                                 if (j == 2)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     da[m_at].typ[1 - m_ord] |= DERIV_UNMARK;
    // INCHI✔️❌:                                     num_cuts -= 1;
    // INCHI✔️❌:                                     num_ring_cuts -= 2;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 comp_num++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (comp_num)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (da[i].typ[0] & DERIV_UNMARK)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             da[i].num[0] = da[i].num[1];
    // INCHI✔️❌:                             da[i].ord[0] = da[i].ord[1];
    // INCHI✔️❌:                             da[i].typ[0] = da[i].typ[1];
    // INCHI✔️❌:                             da[i].typ[1] = 0;
    // INCHI✔️❌:                             j = 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             j = 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (da[i].typ[j] & DERIV_UNMARK)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             da[i].typ[j] = 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( COUNT_ALL_NOT_DERIV == 1 )
    // INCHI✔️❌:                     num_cuts_to_check = num_cuts;
    // INCHI✔️❌: #else
    // INCHI✔️❌:                     num_cuts_to_check = num_ring_cuts;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     if (num_cuts < 0 || num_ring_cuts < 0 || num_cut_pieces < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = -3;
    // INCHI✔️❌:                         goto exit_r2c_num; /* wrong number of cuts = num */
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #if( defined(DERIV_X_OXIME) || defined(DERIV_RO_COX) )
    // INCHI✔️❌:                     if (num_cuts_to_check > 0)
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto repeat_without_deriv_ring;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Sort the bonds for subsequent searching by bisections
    // INCHI✔️❌:                 -- disabled because DERIV_RING_DMOX_DEOX have to be adjacent in ap
    // INCHI✔️❌:                 if ( num_cuts_to_check > 1 ) {
    // INCHI✔️❌:                 qsort( ap, num_cuts_to_check, sizeof(ap[0]), cmp_r2c_atpair);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 if ((ret = mark_deriv_agents( at, da, num_atoms, ap, num_cuts_to_check, &comp_num, &cur_num_at ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_r2c_num; /* wrong atom pair */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (comp_num > 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* eliminate offending DERIV_RING_OUTSIDE_PRECURSOR type derivatives */
    // INCHI✔️❌:                     if (num_ring_cuts <= 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = -99;
    // INCHI✔️❌:                         goto exit_r2c_num;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     n = 0;
    // INCHI✔️❌:                     for (i = j = 0; i < num_atoms; i++) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (( da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int at1a = at[i].neighbor[(int) da[i].ord[0]];
    // INCHI✔️❌:                             int at2a = at[i].neighbor[(int) da[i].ord[1]];
    // INCHI✔️❌:                             if (at[at1a].at_type != at[at2a].at_type)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 da[i].typ[0] = 0; /* eliminate this cut */
    // INCHI✔️❌:                                 da[i].typ[1] = 0;
    // INCHI✔️❌:                                 n++;
    // INCHI✔️❌:                                 num_cuts_to_check -= 2;
    // INCHI✔️❌:                                 num_cuts -= 2;
    // INCHI✔️❌:                                 num_ring_cuts -= 2;
    // INCHI✔️❌:                                 num_cut_pieces -= 1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (n > 0 && num_cuts_to_check > 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:                         goto repeat_without_deriv_ring;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:             exit_r2c_num:
    // INCHI✔️❌:                 /*inchi_free( ap );*/
    // INCHI✔️❌:                 UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:                 /*if ( ret < 0 || num_cuts_to_check >= 2 && cur_num_at < MIN_AT_LEFT_DERIV ) */
    // INCHI✔️❌:                 /* -- bug: cur_num_at may include later rejected deriv. agents 2013-11-08 DCh */
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_function; /* unexpected  error or nothing left */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!num_cuts)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /*goto exit_function;*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Eliminate derivatives that are not in the list */
    // INCHI✔️❌:         num_cuts = eliminate_deriv_not_in_list( at, da, num_atoms, szUnderivList, UNDERIV_LIST_LEN, szUnderivList2, UNDERIV_LIST_LEN2, &bitUnderivList );
    // INCHI✔️❌:         if (num_cuts < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = num_cuts;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Check how many atoms was left in the precursor - begin - 2013-11-12 DT */
    // INCHI✔️❌:         if (( num_cuts_to_check = num_cuts ) >= 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             AT_NUMB comp_num = 0; /* here it is the number of removed cuts */
    // INCHI✔️❌:             int        /*n,*/ m_at, m_ord;
    // INCHI✔️❌:
    // INCHI✔️❌:             ALLOC_AP
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Fill out the array of bonds to be cut */
    // INCHI✔️❌:                 j = fill_out_bond_cuts( at, da, num_atoms, ap, num_cuts_to_check );
    // INCHI✔️❌:             if (j < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = j;
    // INCHI✔️❌:                 goto exit_r2c_num2; /* wrong number of cuts = num */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (j != num_cuts_to_check)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret = -3;
    // INCHI✔️❌:                 goto exit_r2c_num2; /* wrong number of cuts = num */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* The following code was copied from above to mark removed */
    // INCHI✔️❌:             for (i = 0; i < num_cuts_to_check; i += j)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (j = n = 0; j < 2; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int atj = (int) ap[i].at[j];
    // INCHI✔️❌:                     if (atj < num_atoms && da[atj].typ[0] && at[atj].neighbor[(int) da[atj].ord[0]] == ap[i].at[1 - j])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;      /* ap[i].at[k] is precursor atom */
    // INCHI✔️❌:                         n++;
    // INCHI✔️❌:                         m_at = atj; /* precursor atom at[m_at], da[m_at] */
    // INCHI✔️❌:                         m_ord = 0;  /* da[m_at].typ[m_ord] - type of the deriv.bond to break  */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (atj < num_atoms && da[atj].typ[1] && at[atj].neighbor[(int) da[atj].ord[1]] == ap[i].at[1 - j])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = j;
    // INCHI✔️❌:                             n++;
    // INCHI✔️❌:                             m_at = atj;
    // INCHI✔️❌:                             m_ord = 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (n != 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = -3;
    // INCHI✔️❌:                     goto exit_r2c_num2; /* wrong atom pair */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (( da[m_at].typ[m_ord] & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     n = (int) ap[i].at[k];   /* atom inside the derivation attachment */
    // INCHI✔️❌:                     j = 2;             /* number of bonds to cut */
    // INCHI✔️❌:                     if (i + j > num_cuts_to_check || ((int) ap[i + 1].at[0] != n && (int) ap[i + 1].at[1] != n)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = -3;
    // INCHI✔️❌:                         goto exit_r2c_num2; /* wrong atom pair */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     n = ap[i].at[1 - k];    /* atom inside the tentative derivation attachment */
    // INCHI✔️❌:                     j = 1;                  /* number of bonds to cut */
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:                                             /*j += (0 != (da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX));*/
    // INCHI✔️❌:                                             /* these 2 cuts are always adjacent */
    // INCHI✔️❌:                     j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX ) ); /* these 2 cuts are always adjacent */
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_RING2_OUTSIDE_PRECUR
    // INCHI✔️❌:                     j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING2_OUTSIDE_PRECUR ) ); /* these 2 cuts are always adjacent */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* at[n] belongs to the derivation agent  */
    // INCHI✔️❌:                 cur_num_at = mark_atoms_ap( at, n, ap + i, j, 0, 1 );
    // INCHI✔️❌:                 UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             ret = mark_deriv_agents( at, da, num_atoms, ap, num_cuts_to_check, &comp_num, &cur_num_at );
    // INCHI✔️❌:             UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:             if (ap)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( ap );
    // INCHI✔️❌:                 ap = NULL;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (ret)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_r2c_num2; /* wrong atom pair */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     exit_r2c_num2:
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ap)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( ap );
    // INCHI✔️❌:             ap = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (ret < 0 || (num_cuts_to_check >= 2 && cur_num_at < MIN_AT_LEFT_DERIV)) /* -- bug: cur_num_at may include later rejected deriv. agents 2013-11-08 DCh */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function; /* unexpected  error or nothing left */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Check how many atoms was left in the precursor - end - 2013-11-12 DT */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* make cuts */
    // INCHI✔️❌:         num_cuts = 0;
    // INCHI✔️❌:         for (i = num = 0; i < num_atoms; i++) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*for ( len = 0; len < MAX_AT_DERIV && da[i].typ[len]; len ++ ) -- bug fixed 2013-11-07 DCh */
    // INCHI✔️❌:             for (len = 0; len < DERIV_AT_LEN && da[i].typ[len]; len++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             switch (len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case 0:
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 case 1:
    // INCHI✔️❌:                     /* single cut: unconditional */
    // INCHI✔️❌:                     make_single_cut( at, da, i, 0 );
    // INCHI✔️❌:                     num_cuts += 1;
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 case 2:
    // INCHI✔️❌:                     if (( (da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR )) ||
    // INCHI✔️❌:                          (da[i].typ[0] == DERIV_AMINE_tN && da[i].typ[1] == DERIV_AMINE_tN)
    // INCHI✔️❌: #ifdef DERIV_RING2_OUTSIDE_PRECUR
    // INCHI✔️❌:                          || (da[i].typ[0] && da[i].typ[0] == ( da[i].typ[0] & DERIV_RING2_OUTSIDE_PRECUR ) &&
    // INCHI✔️❌:                          da[i].typ[1] == da[i].typ[0])
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                          ) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* double cut, unconditional */
    // INCHI✔️❌:                         make_single_cut( at, da, i, 1 );
    // INCHI✔️❌:                         make_single_cut( at, da, i, 0 );
    // INCHI✔️❌:                         num_cuts += 1;
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (da[i].typ[0] == da[i].typ[1])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* DERIV_BRIDGE_O or DERIV_BRIDGE_NH; cut off the smallest */
    // INCHI✔️❌:                         if (da[i].num[0] > da[i].num[1])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             make_single_cut( at, da, i, 1 );
    // INCHI✔️❌:                             num_cuts += 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (da[i].num[0] < da[i].num[1])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 make_single_cut( at, da, i, 0 );
    // INCHI✔️❌:                                 num_cuts += 1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ret = -88;
    // INCHI✔️❌:                     goto exit_function; /* unexpected */
    // INCHI✔️❌:                 case 3:
    // INCHI✔️❌:                     if (da[i].typ[0] == da[i].typ[1] &&
    // INCHI✔️❌:                          da[i].typ[0] == da[i].typ[2] &&
    // INCHI✔️❌:                          da[i].typ[0] == DERIV_AMINE_tN)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int x, y, z;
    // INCHI✔️❌:                         x = ( da[i].num[0] < da[i].num[1] ) ? 0 : 1;
    // INCHI✔️❌:                         x = ( da[i].num[x] < da[i].num[2] ) ? x : 2; /* min */
    // INCHI✔️❌:                         z = ( da[i].num[0] < da[i].num[1] ) ? 1 : 0;
    // INCHI✔️❌:                         z = ( da[i].num[x] < da[i].num[2] ) ? 2 : z; /* max */
    // INCHI✔️❌:                         y = ( ( x + 1 ) ^ ( z + 1 ) ) - 1;                      /* median */
    // INCHI✔️❌:                         if (da[i].num[x] == da[i].num[z])
    // INCHI✔️❌:                             continue; /* all deriv. agents have same size */
    // INCHI✔️❌:                                       /* two smallest */
    // INCHI✔️❌:                         if (da[i].num[x] == da[i].num[y] && x < y)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int t = x; /* first cut x > y */
    // INCHI✔️❌:                             x = y;
    // INCHI✔️❌:                             y = t;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         make_single_cut( at, da, i, x );
    // INCHI✔️❌:                         num_cuts += 1;
    // INCHI✔️❌:                         if (da[i].num[x] == da[i].num[y])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* equally small */
    // INCHI✔️❌:                             make_single_cut( at, da, i, y );
    // INCHI✔️❌:                             num_cuts += 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ret = -88;
    // INCHI✔️❌:                     goto exit_function; /* unexpected */
    // INCHI✔️❌:                 case 4:
    // INCHI✔️❌:                     if (( da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR ) &&
    // INCHI✔️❌:                         ( da[i].typ[2] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[3] & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int n01 = inchi_max( da[i].num[0], da[i].num[1] );
    // INCHI✔️❌:                         int n23 = inchi_max( da[i].num[2], da[i].num[3] );
    // INCHI✔️❌:                         if (n01 < n23)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             make_single_cut( at, da, i, 1 );
    // INCHI✔️❌:                             make_single_cut( at, da, i, 0 );
    // INCHI✔️❌:                             num_cuts += 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (n01 > n23)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 make_single_cut( at, da, i, 3 );
    // INCHI✔️❌:                                 make_single_cut( at, da, i, 2 );
    // INCHI✔️❌:                                 num_cuts += 1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         nTotNumCuts += num_cuts;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef FIX_UNDERIV_TO_SDF
    // INCHI✔️❌:         if (bOutputSdf && at2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* replace arom bonds with original */
    // INCHI✔️❌:             replace_arom_bonds( at, num_atoms, at2, num_atoms );
    // INCHI✔️❌:             if (at2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( at2 );
    // INCHI✔️❌:                 at2 = NULL;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef UNDERIV_ADD_EXPLICIT_H
    // INCHI✔️❌:         /**********  Add explicit hydrogens ************************/
    // INCHI✔️❌:         num_atoms = add_explicit_H( inp_cur_data + i_component );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         if (num_cuts)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             remove_cut_derivs( num_atoms, at, inp_cur_data, i_component, &ret );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     } /* for (i_component = 0; i_component < num_components; i_component++) */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nTotNumCuts)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         OAD_Edit_MergeComponentsAndRecreateOAD( orig_inp_data, inp_cur_data, num_components, &ret );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     free_underiv_temp_data( ap, da, at2, inp_cur_data, num_components );
    // INCHI✔️❌:
    // INCHI✔️❌: #if( UNDERIVATIZE_REPORT == 1 )
    // INCHI✔️❌:     if (!ret && nTotNumCuts && pSdfValue && bOutputReport)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         numUnderiv = sort_merge_underiv( pSdfValue, bOutputSdf, szUnderivList, cDerivSeparator, underivPrefix, underivPostfix ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:         numUnderiv2 = sort_merge_underiv( pSdfValue, bOutputSdf, szUnderivList2, cDerivSeparator, underivPrefix2, underivPostfix2 ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:         sprintf(szbitUnderivList, "0x%.8X", bitUnderivList); /* djb-rwth: addressing GCC warning about 0 flag being ignored */
    // INCHI✔️❌:         numUnderiv3 = sort_merge_underiv( pSdfValue, bOutputSdf, szbitUnderivList, cDerivSeparator, underivPrefix3, underivPostfix3 ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret ? ret : nTotNumCuts;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: #endif /* UNDERIVATIZE */
    // INCHI✔️❌:
    // END INCHI C FUNCTION: OAD_Edit_Underivatize
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Edit_Underivatize
    // INCHI❗❌: #define UNDERIVATIZE 1
    // INCHI❗❌: #define UNDERIVATIZE_REPORT 1
    // INCHI❗❌: #define FIX_UNDERIV_TO_SDF
    // INCHI❗❌: #define UNDERIV_ADD_EXPLICIT_H
    // INCHI❗❌: #define COUNT_ALL_NOT_DERIV 1
    // INCHI❗❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Edit_Underivatize

    let mut ret = 0_i32;
    let mut num_components = 0_i32;
    let mut total_num_cuts = 0_i32;
    let mut current_num_atoms = 0_i32;
    let mut current = SourceMutPointer::<INP_ATOM_DATA>::null();
    let mut derivatives = SourceMutPointer::<DERIV_AT>::null();
    let mut atom_pairs = SourceMutPointer::<R2C_ATPAIR>::null();
    let mut original_bonds = SourceMutPointer::<inp_ATOM>::null();
    let mut allocated_atom_pairs = 0_i32;
    let mut underiv_list = [0_i8; 2048];
    let mut underiv_list2 = [0_i8; 2048];
    let mut bit_underiv_list = 0_i32;

    let processing = (|| -> Result<(), SourceHeapError> {
        UnMarkDisconnectedComponents(heap, original)?;
        num_components = MarkDisconnectedComponents(heap, original, 0)?;
        current = inchi_calloc(
            heap,
            u64::try_from(num_components)
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
            std::mem::size_of::<INP_ATOM_DATA>() as u64,
        )?;

        'components: for component_index in 0..num_components {
            let component_offset = current.offset(i64::from(component_index))?;
            let mut component = INP_ATOM_DATA::default();
            let component_length = i32::from(
                *heap
                    .slice(original.nCurAtLen.as_const())?
                    .get(component_index as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            CreateInpAtomData(heap, &mut component, component_length, 0)?;
            *heap
                .slice_mut(component_offset)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component.clone();

            let source_atoms = heap
                .slice(original.at.as_const())?
                .get(
                    ..usize::try_from(original.num_inp_atoms)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let mut component_atoms = heap.slice(component.at.as_const())?.to_vec();
            component.num_at = ExtractConnectedComponent(
                heap,
                &source_atoms,
                original.num_inp_atoms,
                component_index + 1,
                &mut component_atoms,
            )?;
            heap.slice_mut(component.at)?
                .clone_from_slice(&component_atoms);
            if component.num_at <= 0 || component_length != component.num_at {
                ret = -(component_index + 1);
                *heap
                    .slice_mut(component_offset)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                break 'components;
            }

            let num_atoms = remove_terminal_HDT(heap, component.num_at, component.at, 1)?;
            component.num_removed_H = component.num_at.wrapping_sub(num_atoms);
            component.num_at = num_atoms;
            let mut num_atoms = component.num_at;
            {
                let atoms = heap.slice_mut(component.at)?;
                add_DT_to_num_H(num_atoms, atoms)?;
                UnMarkOtherIndicators(atoms, num_atoms)?;
                UnMarkOneComponent(atoms, num_atoms)?;
            }
            UnMarkRingSystemsInp(heap, component.at, num_atoms)?;
            MarkRingSystemsInp(heap, component.at, num_atoms, 0)?;

            if output_sdf != 0 {
                if !original_bonds.is_null() {
                    inchi_free(heap, original_bonds)?;
                    original_bonds = SourceMutPointer::null();
                }
                let atoms = heap
                    .slice(component.at.as_const())?
                    .get(
                        ..usize::try_from(num_atoms)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                original_bonds = match heap.allocate(atoms) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
            }

            ret = mark_arom_bonds(
                heap,
                clock,
                canon_globals,
                component.at,
                num_atoms,
                clock_result,
            )?;
            if ret < 0 {
                *heap
                    .slice_mut(component_offset)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                break 'components;
            }
            ret = 0;

            if !derivatives.is_null() {
                inchi_free(heap, derivatives)?;
            }
            derivatives = inchi_calloc(
                heap,
                u64::try_from(num_atoms)
                    .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
                std::mem::size_of::<DERIV_AT>() as u64,
            )?;
            let atom_count =
                usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut atoms = heap
                .slice(component.at.as_const())?
                .get(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let mut derivative_values = heap
                .slice(derivatives.as_const())?
                .get(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();

            let mut found = 0_i32;
            for atom_index in 0..atom_count {
                let valence = usize::try_from(atoms[atom_index].valence)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if atoms[atom_index].bCutVertex != 0
                    && derivative_values[atom_index].typ[valence.min(4) - 1] == 0
                {
                    for ordinal in 0..valence {
                        let count = count_one_bond_atoms(
                            &mut atoms,
                            &mut derivative_values,
                            atom_index as i32,
                            ordinal as i32,
                            CFLAG_MARK_BRANCH as i8,
                            &mut found,
                        )?;
                        UnMarkOtherIndicators(&mut atoms, num_atoms)?;
                        if count < 0 {
                            ret = count;
                            break;
                        }
                    }
                    if ret < 0 {
                        break;
                    }
                }
            }
            if ret < 0 {
                heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
                *heap
                    .slice_mut(component_offset)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                break 'components;
            }

            let (mut num_cuts, mut num_ring_cuts, mut num_cut_pieces) =
                match oad_prepare_derivative_cuts(&atoms, &mut derivative_values, num_atoms)? {
                    Ok(values) => values,
                    Err(value) => {
                        ret = value;
                        heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                        heap.slice_mut(derivatives)?[..atom_count]
                            .clone_from_slice(&derivative_values);
                        *heap
                            .slice_mut(component_offset)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                        break 'components;
                    }
                };
            let mut num_cuts_to_check = num_cuts;
            if !oad_allocate_atom_pairs(
                heap,
                &mut atom_pairs,
                &mut allocated_atom_pairs,
                num_cuts_to_check,
            )? {
                ret = -1;
                heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
                *heap
                    .slice_mut(component_offset)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                break 'components;
            }
            if num_cuts_to_check >= 2 {
                let mut pair_values = heap.slice(atom_pairs.as_const())?.to_vec();
                ret = oad_validate_derivative_agents(
                    &mut atoms,
                    &mut derivative_values,
                    num_atoms,
                    &mut pair_values,
                    &mut num_cuts,
                    &mut num_ring_cuts,
                    &mut num_cut_pieces,
                    &mut num_cuts_to_check,
                    &mut current_num_atoms,
                )?;
                heap.slice_mut(atom_pairs)?.clone_from_slice(&pair_values);
                if ret < 0 {
                    heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                    heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
                    *heap
                        .slice_mut(component_offset)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                    break 'components;
                }
            }
            if num_cuts == 0 {
                heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
                *heap
                    .slice_mut(component_offset)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                continue;
            }

            num_cuts = eliminate_deriv_not_in_list(
                &atoms,
                &mut derivative_values,
                num_atoms,
                Some(&mut underiv_list),
                2048,
                Some(&mut underiv_list2),
                2048,
                Some(&mut bit_underiv_list),
            )?;
            if num_cuts < 0 {
                ret = num_cuts;
                heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
                *heap
                    .slice_mut(component_offset)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                break 'components;
            }

            num_cuts_to_check = num_cuts;
            if num_cuts_to_check >= 1 {
                if !oad_allocate_atom_pairs(
                    heap,
                    &mut atom_pairs,
                    &mut allocated_atom_pairs,
                    num_cuts_to_check,
                )? {
                    ret = -1;
                } else {
                    let mut pair_values = heap.slice(atom_pairs.as_const())?.to_vec();
                    ret = oad_check_final_precursor(
                        &mut atoms,
                        &derivative_values,
                        num_atoms,
                        &mut pair_values,
                        num_cuts_to_check,
                        &mut current_num_atoms,
                    )?;
                    heap.slice_mut(atom_pairs)?.clone_from_slice(&pair_values);
                }
            }
            if !atom_pairs.is_null() {
                inchi_free(heap, atom_pairs)?;
                atom_pairs = SourceMutPointer::null();
            }
            if ret < 0 || (num_cuts_to_check >= 2 && current_num_atoms < MIN_AT_LEFT_DERIV as i32) {
                heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
                *heap
                    .slice_mut(component_offset)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                break 'components;
            }

            num_cuts = match oad_make_selected_cuts(&mut atoms, &mut derivative_values, num_atoms)?
            {
                Ok(value) => value,
                Err(value) => {
                    ret = value;
                    heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
                    heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
                    *heap
                        .slice_mut(component_offset)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
                    break 'components;
                }
            };
            total_num_cuts = total_num_cuts.wrapping_add(num_cuts);

            if output_sdf != 0 && !original_bonds.is_null() {
                let original_atoms = heap.slice(original_bonds.as_const())?.to_vec();
                replace_arom_bonds(&mut atoms, num_atoms, &original_atoms, num_atoms)?;
                inchi_free(heap, original_bonds)?;
                original_bonds = SourceMutPointer::null();
            }
            heap.slice_mut(component.at)?[..atom_count].clone_from_slice(&atoms);
            heap.slice_mut(derivatives)?[..atom_count].clone_from_slice(&derivative_values);
            *heap
                .slice_mut(component_offset)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component.clone();

            num_atoms = add_explicit_H(heap, &mut component)?;
            *heap
                .slice_mut(component_offset)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component.clone();
            if num_cuts != 0 {
                let mut components = heap.slice(current.as_const())?.to_vec();
                remove_cut_derivs(
                    heap,
                    num_atoms,
                    component.at,
                    &mut components,
                    component_index,
                    &mut ret,
                )?;
                heap.slice_mut(current)?.clone_from_slice(&components);
                component = heap
                    .slice(component_offset.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
            }
            *heap
                .slice_mut(component_offset)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
        }

        if total_num_cuts != 0 && ret == 0 {
            let mut components = heap.slice(current.as_const())?.to_vec();
            OAD_Edit_MergeComponentsAndRecreateOAD(
                heap,
                original,
                &mut components,
                num_components,
                &mut ret,
            )?;
            heap.slice_mut(current)?.clone_from_slice(&components);
        }
        Ok(())
    })();

    let cleanup = if current.is_null() {
        if !atom_pairs.is_null() {
            inchi_free(heap, atom_pairs)?;
        }
        if !derivatives.is_null() {
            inchi_free(heap, derivatives)?;
        }
        if !original_bonds.is_null() {
            inchi_free(heap, original_bonds)?;
        }
        Ok(())
    } else {
        free_underiv_temp_data(
            heap,
            atom_pairs,
            derivatives,
            original_bonds,
            current,
            num_components,
        )
    };
    match (processing, cleanup) {
        (Err(error), _) => return Err(error),
        (Ok(()), Err(error)) => return Err(error),
        (Ok(()), Ok(())) => {}
    }

    if ret == 0 && total_num_cuts != 0 && !sdf_value.is_null() && output_report != 0 {
        let target = heap.slice_mut(sdf_value)?;
        let _ = sort_merge_underiv(
            target,
            output_sdf,
            &mut underiv_list,
            b',' as i8,
            &b"\tDeriv=\0".map(|byte| byte as i8),
            &[0],
        )?;
        let _ = sort_merge_underiv(
            target,
            output_sdf,
            &mut underiv_list2,
            b',' as i8,
            &b"\tDeriv2=\0".map(|byte| byte as i8),
            &[0],
        )?;
        let text = format!("0x{:08X}", bit_underiv_list as u32);
        let mut bit_list = [0_i8; 16];
        for (target, source) in bit_list.iter_mut().zip(text.bytes()) {
            *target = source as i8;
        }
        let _ = sort_merge_underiv(
            target,
            output_sdf,
            &mut bit_list,
            b',' as i8,
            &b"\tDerivBits=\0".map(|byte| byte as i8),
            &[0],
        )?;
    }
    Ok(if ret != 0 { ret } else { total_num_cuts })
}

#[allow(non_snake_case)]
pub(crate) fn detect_r2c_Zatom(
    atoms: &[inp_ATOM],
    derivatives: &mut [R2C_AT],
    iZ: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7159 detect_r2c_Zatom
    // INCHI✔️❌: int detect_r2c_Zatom( inp_ATOM *at, R2C_AT *da, int iZ )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, neigh, neighneigh, nRingSystem, num_found;
    // INCHI✔️❌:     R2C_AT da1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[iZ].valence > 4)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (at[iZ].valence != at[iZ].chem_bonds_valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* approach limitation: no double bonds */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (at[iZ].el_number != EL_NUMBER_C)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* sugar-specific */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (at[iZ].nNumAtInRingSystem < 5)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* not in a suitable ring */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!at[iZ].bCutVertex)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;  /* recognize only type 1 for now */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRingSystem = at[iZ].nRingSystem;
    // INCHI✔️❌:     memset( &da1, R2C_EMPTY, sizeof( da1 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0, num_found = 0; i < at[iZ].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[iZ].neighbor[i];
    // INCHI✔️❌:         if (at[neigh].charge || at[neigh].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].el_number == EL_NUMBER_O &&
    // INCHI✔️❌:              at[neigh].valence == 1 &&
    // INCHI✔️❌:              at[neigh].chem_bonds_valence == 1 &&
    // INCHI✔️❌:              at[neigh].num_H == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* found Z-OH, i.e. Z-YH */
    // INCHI✔️❌:             if (da1.ordY == R2C_EMPTY)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1.ordY = i;
    // INCHI✔️❌:                 num_found++;
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].el_number == EL_NUMBER_O &&
    // INCHI✔️❌:              at[neigh].valence == 2 &&
    // INCHI✔️❌:              at[neigh].chem_bonds_valence == 2 &&
    // INCHI✔️❌:              at[neigh].num_H == 0 &&
    // INCHI✔️❌:              at[neigh].nRingSystem == nRingSystem)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* found Z-O-, i.e. Z-W- */
    // INCHI✔️❌:             if (da1.ordW == R2C_EMPTY)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* j = index of the oppozite to at[iZ] neighbor of at[neigh] */
    // INCHI✔️❌:                 j = ( at[neigh].neighbor[0] == iZ );
    // INCHI✔️❌:                 neighneigh = at[neigh].neighbor[j];
    // INCHI✔️❌:                 if (at[neighneigh].valence != at[neighneigh].chem_bonds_valence ||
    // INCHI✔️❌:                      at[neighneigh].el_number != EL_NUMBER_C)
    // INCHI✔️❌:                     return 0; /* sugar-specific */
    // INCHI✔️❌:                 da1.ordW = i;
    // INCHI✔️❌:                 num_found++;
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:              at[neigh].valence > 2 &&
    // INCHI✔️❌:              at[neigh].chem_bonds_valence == at[neigh].valence &&
    // INCHI✔️❌:              at[neigh].num_H <= 1 &&
    // INCHI✔️❌:              at[neigh].nRingSystem == nRingSystem)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* sugar-specfic: carbon in the ring should have -OH neighbor */
    // INCHI✔️❌:             int iOH;
    // INCHI✔️❌:             for (j = 0; j < at[neigh].valence; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 iOH = at[neigh].neighbor[j];
    // INCHI✔️❌:                 if (at[iOH].el_number == EL_NUMBER_O &&
    // INCHI✔️❌:                      at[iOH].valence == 1 &&
    // INCHI✔️❌:                      at[iOH].chem_bonds_valence == 1 &&
    // INCHI✔️❌:                      at[iOH].num_H == 1 &&
    // INCHI✔️❌:                      !at[iOH].charge && !at[iOH].radical)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (da1.ordC == R2C_EMPTY)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.ordC = i;
    // INCHI✔️❌:                         num_found++;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (j < at[neigh].valence)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:              at[neigh].chem_bonds_valence == at[neigh].valence &&
    // INCHI✔️❌:              at[neigh].nRingSystem != nRingSystem)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* extra carbon neighbor of Z */
    // INCHI✔️❌:             if (da1.ordCopt == R2C_EMPTY)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1.ordCopt = i;
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return 0; /* unexpectd neighbor */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_found == 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         da1.type = 1;
    // INCHI✔️❌:         da[iZ] = da1;
    // INCHI✔️❌:         return 1; /* disconnection found */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: detect_r2c_Zatom
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: detect_r2c_Zatom
    // INCHI✔️❌: #define R2C_EMPTY  127
    // INCHI✔️❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️❌: #define EL_NUMBER_O  ((U_CHAR)8)
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: detect_r2c_Zatom

    let center_index = usize::try_from(iZ).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let center = atoms
        .get(center_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if center.valence > 4
        || center.valence != center.chem_bonds_valence
        || center.el_number != EL_NUMBER_C
        || center.nNumAtInRingSystem < 5
        || center.bCutVertex == 0
    {
        return Ok(0);
    }

    let ring_system = center.nRingSystem;
    let empty = R2C_EMPTY as i8;
    let mut candidate = R2C_AT {
        type_: empty,
        ordW: empty,
        ordY: empty,
        ordC: empty,
        ordCopt: empty,
    };
    let mut index = 0_i32;
    let mut num_found = 0_i32;
    while index < i32::from(center.valence) {
        let ordinal = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let neighbor_index = usize::from(
            *center
                .neighbor
                .get(ordinal)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.charge != 0 || neighbor.radical != 0 {
            return Ok(0);
        }
        if neighbor.el_number == EL_NUMBER_O
            && neighbor.valence == 1
            && neighbor.chem_bonds_valence == 1
            && neighbor.num_H == 1
        {
            if candidate.ordY == empty {
                candidate.ordY = index as i8;
                num_found += 1;
                index += 1;
                continue;
            }
            return Ok(0);
        }
        if neighbor.el_number == EL_NUMBER_O
            && neighbor.valence == 2
            && neighbor.chem_bonds_valence == 2
            && neighbor.num_H == 0
            && neighbor.nRingSystem == ring_system
        {
            if candidate.ordW != empty {
                return Ok(0);
            }
            let opposite_ordinal = usize::from(i32::from(neighbor.neighbor[0]) == iZ);
            let opposite_index = usize::from(neighbor.neighbor[opposite_ordinal]);
            let opposite = atoms
                .get(opposite_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if opposite.valence != opposite.chem_bonds_valence || opposite.el_number != EL_NUMBER_C
            {
                return Ok(0);
            }
            candidate.ordW = index as i8;
            num_found += 1;
            index += 1;
            continue;
        }
        if neighbor.el_number == EL_NUMBER_C
            && neighbor.valence > 2
            && neighbor.chem_bonds_valence == neighbor.valence
            && neighbor.num_H <= 1
            && neighbor.nRingSystem == ring_system
        {
            let mut neighbor_ordinal = 0_i32;
            while neighbor_ordinal < i32::from(neighbor.valence) {
                let hydroxyl_index = usize::from(
                    neighbor.neighbor[usize::try_from(neighbor_ordinal)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?],
                );
                let hydroxyl = atoms
                    .get(hydroxyl_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if hydroxyl.el_number == EL_NUMBER_O
                    && hydroxyl.valence == 1
                    && hydroxyl.chem_bonds_valence == 1
                    && hydroxyl.num_H == 1
                    && hydroxyl.charge == 0
                    && hydroxyl.radical == 0
                {
                    if candidate.ordC != empty {
                        return Ok(0);
                    }
                    candidate.ordC = index as i8;
                    num_found += 1;
                    break;
                }
                neighbor_ordinal += 1;
            }
            if neighbor_ordinal < i32::from(neighbor.valence) {
                index += 1;
                continue;
            }
        }
        if neighbor.el_number == EL_NUMBER_C
            && neighbor.chem_bonds_valence == neighbor.valence
            && neighbor.nRingSystem != ring_system
        {
            if candidate.ordCopt == empty {
                candidate.ordCopt = index as i8;
                index += 1;
                continue;
            }
        }
        return Ok(0);
    }

    if num_found == 3 {
        candidate.type_ = 1;
        *derivatives
            .get_mut(center_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = candidate;
        Ok(1)
    } else {
        Ok(0)
    }
}

pub(crate) fn cut_ring_to_chain(
    atoms: &mut [inp_ATOM],
    derivatives: &[R2C_AT],
    iZ: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7296 cut_ring_to_chain
    // INCHI✔️❌: int cut_ring_to_chain( inp_ATOM *at, R2C_AT *da, int iZ )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = -1; /* error flag */
    // INCHI✔️❌:     int iordW = (int) da[iZ].ordW; /* ord of the bond in iZ */
    // INCHI✔️❌:     int iordY = (int) da[iZ].ordY; /* ord of the bond in iZ */
    // INCHI✔️❌:     int iordC = (int) da[iZ].ordC;
    // INCHI✔️❌:     int iW, iY, num_iso_H, i, jordZ;
    // INCHI✔️❌:     AT_NUMB *p;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (da[iZ].type != 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (0 > iordW || iordW >= at[iZ].valence ||
    // INCHI✔️❌:          0 > iordY || iordY >= at[iZ].valence ||
    // INCHI✔️❌:          0 > iordC || iordC >= at[iZ].valence /* suger-specific*/)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* find other da[] that affect at[iZ] */
    // INCHI✔️❌:     iW = at[iZ].neighbor[iordW];  /* opposite atom to disconnect and add H */
    // INCHI✔️❌:     iY = at[iZ].neighbor[iordY];  /* opposite atom to increment the bond and remove H*/
    // INCHI✔️❌:     if (!at[iY].num_H || at[iZ].bond_type[iordY] != BOND_TYPE_SINGLE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -2; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* increment at[iZ]--at[iY] bond */
    // INCHI✔️❌:     p = is_in_the_list( at[iY].neighbor, (AT_NUMB) iZ, at[iY].valence );
    // INCHI✔️❌:     if (!p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -3; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     jordZ = p - at[iY].neighbor;
    // INCHI✔️❌:     at[iZ].bond_type[iordY] ++;
    // INCHI✔️❌:     at[iZ].chem_bonds_valence++;
    // INCHI✔️❌:     at[iY].bond_type[jordZ] ++;
    // INCHI✔️❌:     at[iY].chem_bonds_valence++;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* disconnect at[iZ]--at[iW] bond */
    // INCHI✔️❌:     ret = DisconnectInpAtBond( at, NULL, iZ, iordW );
    // INCHI✔️❌:     if (ret != 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -4; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* disconnection succeeded */
    // INCHI✔️❌:     /* transfer H from at[iY] to at[iW] */
    // INCHI✔️❌:     num_iso_H = NUM_ISO_H( at, iY );
    // INCHI✔️❌:     if (at[iY].num_H == num_iso_H)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < NUM_H_ISOTOPES; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[iY].num_iso_H[i])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[iY].num_iso_H[i] --;
    // INCHI✔️❌:                 at[iW].num_iso_H[i] ++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     at[iY].num_H--;
    // INCHI✔️❌:     at[iW].num_H++;
    // INCHI✔️❌:
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: cut_ring_to_chain
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: cut_ring_to_chain
    // INCHI✔️❌: #define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
    // INCHI✔️❌: #define NUM_H_ISOTOPES        3  /* number of hydrogen isotopes: protium, deuterium, tritium */
    // INCHI✔️❌: #define BOND_TYPE_SINGLE    1
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: cut_ring_to_chain

    let center_index = usize::try_from(iZ).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let derivative = derivatives
        .get(center_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if derivative.type_ != 1 {
        return Ok(0);
    }
    let order_w = i32::from(derivative.ordW);
    let order_y = i32::from(derivative.ordY);
    let order_c = i32::from(derivative.ordC);
    let center_valence = i32::from(
        atoms
            .get(center_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .valence,
    );
    if order_w < 0
        || order_w >= center_valence
        || order_y < 0
        || order_y >= center_valence
        || order_c < 0
        || order_c >= center_valence
    {
        return Ok(-1);
    }
    let order_w_index =
        usize::try_from(order_w).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let order_y_index =
        usize::try_from(order_y).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let neighbor_w = usize::from(atoms[center_index].neighbor[order_w_index]);
    let neighbor_y = usize::from(atoms[center_index].neighbor[order_y_index]);
    if atoms
        .get(neighbor_y)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .num_H
        == 0
        || atoms[center_index].bond_type[order_y_index] != BOND_SINGLE as u8
    {
        return Ok(-2);
    }
    let reverse_y = is_in_the_list(
        Some(&atoms[neighbor_y].neighbor),
        iZ as u16,
        i32::from(atoms[neighbor_y].valence),
    )?;
    let Some(reverse_y) = reverse_y else {
        return Ok(-3);
    };

    atoms[center_index].bond_type[order_y_index] =
        atoms[center_index].bond_type[order_y_index].wrapping_add(1);
    atoms[center_index].chem_bonds_valence = atoms[center_index].chem_bonds_valence.wrapping_add(1);
    atoms[neighbor_y].bond_type[reverse_y] = atoms[neighbor_y].bond_type[reverse_y].wrapping_add(1);
    atoms[neighbor_y].chem_bonds_valence = atoms[neighbor_y].chem_bonds_valence.wrapping_add(1);

    if DisconnectInpAtBond(atoms, None, iZ, order_w)? != 1 {
        return Ok(-4);
    }

    let num_iso_h = i32::from(atoms[neighbor_y].num_iso_H[0])
        + i32::from(atoms[neighbor_y].num_iso_H[1])
        + i32::from(atoms[neighbor_y].num_iso_H[2]);
    if i32::from(atoms[neighbor_y].num_H) == num_iso_h {
        for isotope in 0..3 {
            if atoms[neighbor_y].num_iso_H[isotope] != 0 {
                atoms[neighbor_y].num_iso_H[isotope] =
                    atoms[neighbor_y].num_iso_H[isotope].wrapping_sub(1);
                atoms[neighbor_w].num_iso_H[isotope] =
                    atoms[neighbor_w].num_iso_H[isotope].wrapping_add(1);
            }
        }
    }
    atoms[neighbor_y].num_H = atoms[neighbor_y].num_H.wrapping_sub(1);
    atoms[neighbor_w].num_H = atoms[neighbor_w].num_H.wrapping_add(1);
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn Ring2Chain(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canon_globals: &mut CANON_GLOBALS,
    original: &mut ORIG_ATOM_DATA,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7362 Ring2Chain
    // INCHI✔️❌: int Ring2Chain( struct tagINCHI_CLOCK *ic,
    // INCHI✔️❌:                 struct tagCANON_GLOBALS *pCG,
    // INCHI✔️❌:                 ORIG_ATOM_DATA *orig_inp_data )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, i, j, n, num_atoms, num_components, num, num_cuts, iZ; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     inp_ATOM *at = orig_inp_data->at;
    // INCHI✔️❌:     INP_ATOM_DATA *inp_cur_data = NULL;
    // INCHI✔️❌:     R2C_AT        *da = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* prepare */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*set_R2C_el_numbers( );*/
    // INCHI✔️❌:
    // INCHI✔️❌:     num_atoms = remove_terminal_HDT( orig_inp_data->num_inp_atoms, at, 1 );
    // INCHI✔️❌:     /*^^^^^ always accomodate accomodate FIX_TERM_H_CHRG_BUG - IPl, July 2008*/
    // INCHI✔️❌:     orig_inp_data->num_inp_atoms = num_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* initialize */
    // INCHI✔️❌:     UnMarkDisconnectedComponents( orig_inp_data );
    // INCHI✔️❌:     num_cuts = 0;
    // INCHI✔️❌:     /* mark */
    // INCHI✔️❌:     num_components = MarkDisconnectedComponents( orig_inp_data, 0 );
    // INCHI✔️❌:     inp_cur_data = (INP_ATOM_DATA *) inchi_calloc( num_components, sizeof( inp_cur_data[0] ) );
    // INCHI✔️❌:     iZ = -1;
    // INCHI✔️❌:     for (j = 0; j < num_components; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         CreateInpAtomData( inp_cur_data + j, orig_inp_data->nCurAtLen[j], 0 );
    // INCHI✔️❌:         inp_cur_data[j].num_at = ExtractConnectedComponent( orig_inp_data->at, orig_inp_data->num_inp_atoms, j + 1, inp_cur_data[j].at );
    // INCHI✔️❌:         /*  error processing */
    // INCHI✔️❌:         if (inp_cur_data[j].num_at <= 0 || orig_inp_data->nCurAtLen[j] != inp_cur_data[j].num_at)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = -( j + 1 ); /* severe error */
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* initialize */
    // INCHI✔️❌:         num_atoms = inp_cur_data[j].num_at;
    // INCHI✔️❌:         at = inp_cur_data[j].at;
    // INCHI✔️❌:         add_DT_to_num_H( num_atoms, at );
    // INCHI✔️❌:
    // INCHI✔️❌:         UnMarkRingSystemsInp( at, num_atoms );
    // INCHI✔️❌:         UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:         UnMarkOneComponent( at, num_atoms );
    // INCHI✔️❌:         MarkRingSystemsInp( at, num_atoms, 0 );
    // INCHI✔️❌:         ret = mark_arom_bonds( ic, pCG, at, num_atoms );
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ret = 0;
    // INCHI✔️❌:         if (da)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( da );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         da = (R2C_AT *) inchi_calloc( num_atoms, sizeof( da[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* detect ring-to-chain possibilities */
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         for (i = 0, num = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[i].bCutVertex /* type 1 specific*/ && !da[i].type)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num += ( n = detect_r2c_Zatom( at, da, i ) );
    // INCHI✔️❌:                 if (n == 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iZ = i;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (num == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* convert ring to chain: make single cut */
    // INCHI✔️❌:             ret = cut_ring_to_chain( at, da, iZ );
    // INCHI✔️❌:             if (ret < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             num_cuts += ( ret == 1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (num)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* allocate an array of bonds to be cut */
    // INCHI✔️❌:                 R2C_ATPAIR *ap = (R2C_ATPAIR *) inchi_malloc( sizeof( ap[0] ) * num );
    // INCHI✔️❌:                 AT_NUMB    comp_num = 0;
    // INCHI✔️❌:                 if (!ap)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = -1; /* malloc failure */
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* fill out the array of bonds to be cut */
    // INCHI✔️❌:                 for (i = j = 0; i < num_atoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (da[i].type == 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         AT_NUMB at1 = i;
    // INCHI✔️❌:                         AT_NUMB at2 = at[i].neighbor[(int) da[i].ordW];
    // INCHI✔️❌:                         if (j >= num)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ret = -2;
    // INCHI✔️❌:                             goto exit_r2c_num; /* wrong number of cuts = num */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         n = ( at1 > at2 );
    // INCHI✔️❌:                         ap[j].at[n] = at1;
    // INCHI✔️❌:                         ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
    // INCHI✔️❌:                         j++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (j != num)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = -3;
    // INCHI✔️❌:                     goto exit_r2c_num; /* wrong number of cuts = num */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* sort the bonds for subsequent searching by bisections */
    // INCHI✔️❌:                 qsort( ap, num, sizeof( ap[0] ), cmp_r2c_atpair );
    // INCHI✔️❌:                 /* mark components to be disconnected */
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (j = 0; j < 2; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (!at[ap[i].at[j]].at_type)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             comp_num++;
    // INCHI✔️❌:                             mark_atoms_ap( at, (int) ap[i].at[j], ap, num, 0, comp_num );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* convert ring to chain */
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int i1 = ap[i].at[0];
    // INCHI✔️❌:                     int i2 = ap[i].at[1];
    // INCHI✔️❌:                     iZ = -1;
    // INCHI✔️❌:                     if (at[i1].at_type == at[i2].at_type)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* by definition, there are no adjacent iZ atoms; one iZ atom per bond */
    // INCHI✔️❌:                         if (da[i1].type == 1 && at[i1].neighbor[(int) da[i1].ordW] == i2)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             iZ = i1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                             if (da[i2].type == 1 && at[i2].neighbor[(int) da[i2].ordW] == i1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 iZ = i2;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ret = -3;
    // INCHI✔️❌:                                 goto exit_r2c_num;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         ret = cut_ring_to_chain( at, da, iZ );
    // INCHI✔️❌:                         if (ret < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto exit_r2c_num;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         num_cuts += ( ret == 1 );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 ret = 0;
    // INCHI✔️❌:             exit_r2c_num:
    // INCHI✔️❌:                 inchi_free( ap );
    // INCHI✔️❌:                 UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_cuts)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         OAD_Edit_MergeComponentsAndRecreateOAD( orig_inp_data, inp_cur_data, num_components, &ret );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (da)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( da );
    // INCHI✔️❌:         da = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (j = 0; j < num_components; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         FreeInpAtomData( inp_cur_data + j );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_free( inp_cur_data );
    // INCHI✔️❌:     inp_cur_data = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret ? ret : num_cuts;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Ring2Chain
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: Ring2Chain
    // INCHI✔️❌: #define RING2CHAIN                 1 /* open rings R-C(-OH)-O-R => R-C(=O) OH-R   */
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: Ring2Chain

    let mut ret = 0_i32;
    let mut num_components = 0_i32;
    let mut num_cuts = 0_i32;
    let mut current = SourceMutPointer::<INP_ATOM_DATA>::null();
    let mut derivatives = SourceMutPointer::<R2C_AT>::null();

    let processing = (|| -> Result<(), SourceHeapError> {
        let num_atoms = remove_terminal_HDT(heap, original.num_inp_atoms, original.at, 1)?;
        original.num_inp_atoms = num_atoms;
        UnMarkDisconnectedComponents(heap, original)?;
        num_components = MarkDisconnectedComponents(heap, original, 0)?;
        current = inchi_calloc(
            heap,
            u64::try_from(num_components)
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
            std::mem::size_of::<INP_ATOM_DATA>() as u64,
        )?;

        'components: for component_index in 0..num_components {
            let component_pointer = current.offset(i64::from(component_index))?;
            let component_length = i32::from(
                *heap
                    .slice(original.nCurAtLen.as_const())?
                    .get(component_index as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let mut component = INP_ATOM_DATA::default();
            if CreateInpAtomData(heap, &mut component, component_length, 0)? == 0
                || component.at.is_null()
            {
                return Err(SourceHeapError::AllocationFailed);
            }
            *heap
                .slice_mut(component_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component.clone();

            let source_atoms = heap
                .slice(original.at.as_const())?
                .get(
                    ..usize::try_from(original.num_inp_atoms)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let mut component_atoms = heap.slice(component.at.as_const())?.to_vec();
            component.num_at = ExtractConnectedComponent(
                heap,
                &source_atoms,
                original.num_inp_atoms,
                component_index + 1,
                &mut component_atoms,
            )?;
            heap.slice_mut(component.at)?
                .clone_from_slice(&component_atoms);
            *heap
                .slice_mut(component_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component.clone();
            if component.num_at <= 0 || component_length != component.num_at {
                ret = -(component_index + 1);
                break 'components;
            }

            let component_atom_count = usize::try_from(component.num_at)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            {
                let atoms = heap.slice_mut(component.at)?;
                add_DT_to_num_H(component.num_at, atoms)?;
                UnMarkOtherIndicators(atoms, component.num_at)?;
                UnMarkOneComponent(atoms, component.num_at)?;
            }
            UnMarkRingSystemsInp(heap, component.at, component.num_at)?;
            MarkRingSystemsInp(heap, component.at, component.num_at, 0)?;
            ret = mark_arom_bonds(
                heap,
                clock,
                canon_globals,
                component.at,
                component.num_at,
                clock_result,
            )?;
            if ret < 0 {
                break 'components;
            }
            ret = 0;

            if !derivatives.is_null() {
                inchi_free(heap, derivatives)?;
                derivatives = SourceMutPointer::null();
            }
            derivatives = inchi_calloc(
                heap,
                u64::try_from(component.num_at)
                    .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
                std::mem::size_of::<R2C_AT>() as u64,
            )?;
            let mut atoms = heap
                .slice(component.at.as_const())?
                .get(..component_atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let mut derivative_values = heap
                .slice(derivatives.as_const())?
                .get(..component_atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let mut num = 0_i32;
            let mut center = -1_i32;
            for atom_index in 0..component_atom_count {
                if atoms[atom_index].bCutVertex != 0 && derivative_values[atom_index].type_ == 0 {
                    let found =
                        detect_r2c_Zatom(&atoms, &mut derivative_values, atom_index as i32)?;
                    num = num.wrapping_add(found);
                    if found == 1 {
                        center = atom_index as i32;
                    }
                    UnMarkOtherIndicators(&mut atoms, component.num_at)?;
                }
            }

            if num == 1 {
                ret = cut_ring_to_chain(&mut atoms, &derivative_values, center)?;
                if ret < 0 {
                    heap.slice_mut(component.at)?[..component_atom_count].clone_from_slice(&atoms);
                    heap.slice_mut(derivatives)?[..component_atom_count]
                        .clone_from_slice(&derivative_values);
                    break 'components;
                }
                num_cuts = num_cuts.wrapping_add(i32::from(ret == 1));
            } else if num != 0 {
                let pair_count = usize::try_from(num)
                    .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
                let atom_pairs = match heap.allocate(vec![R2C_ATPAIR::default(); pair_count]) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        ret = -1;
                        heap.slice_mut(component.at)?[..component_atom_count]
                            .clone_from_slice(&atoms);
                        heap.slice_mut(derivatives)?[..component_atom_count]
                            .clone_from_slice(&derivative_values);
                        break 'components;
                    }
                    Err(error) => return Err(error),
                };
                let pair_processing = (|| -> Result<(), SourceHeapError> {
                    let pairs = heap.slice_mut(atom_pairs)?;
                    let mut pair_index = 0_usize;
                    for atom_index in 0..component_atom_count {
                        if derivative_values[atom_index].type_ == 1 {
                            if pair_index >= pair_count {
                                ret = -2;
                                return Ok(());
                            }
                            let ordinal =
                                usize::try_from(i32::from(derivative_values[atom_index].ordW))
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let first = atom_index as u16;
                            let second = atoms[atom_index].neighbor[ordinal];
                            let position = usize::from(first > second);
                            pairs[pair_index].at[position] = first;
                            pairs[pair_index].at[1 - position] = second;
                            pair_index += 1;
                        }
                    }
                    if pair_index != pair_count {
                        ret = -3;
                        return Ok(());
                    }
                    pairs.sort_unstable_by(|left, right| cmp_r2c_atpair(left, right).cmp(&0));
                    let pair_values = pairs.to_vec();
                    let mut component_number = 0_u16;
                    for pair in &pair_values {
                        for side in 0..2 {
                            let endpoint = usize::from(pair.at[side]);
                            if atoms
                                .get(endpoint)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .at_type
                                == 0
                            {
                                component_number = component_number.wrapping_add(1);
                                mark_atoms_ap(
                                    &mut atoms,
                                    pair.at[side],
                                    &pair_values,
                                    num,
                                    0,
                                    component_number,
                                )?;
                            }
                        }
                    }
                    for pair in &pair_values {
                        let first = usize::from(pair.at[0]);
                        let second = usize::from(pair.at[1]);
                        if atoms[first].at_type != atoms[second].at_type {
                            continue;
                        }
                        let selected = if derivative_values[first].type_ == 1 {
                            let ordinal = usize::try_from(i32::from(derivative_values[first].ordW))
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            (atoms[first].neighbor[ordinal] == pair.at[1]).then_some(first)
                        } else {
                            None
                        }
                        .or_else(|| {
                            if derivative_values[second].type_ != 1 {
                                return None;
                            }
                            let ordinal =
                                usize::try_from(i32::from(derivative_values[second].ordW)).ok()?;
                            (atoms[second].neighbor[ordinal] == pair.at[0]).then_some(second)
                        });
                        let Some(selected) = selected else {
                            ret = -3;
                            return Ok(());
                        };
                        ret = cut_ring_to_chain(&mut atoms, &derivative_values, selected as i32)?;
                        if ret < 0 {
                            return Ok(());
                        }
                        num_cuts = num_cuts.wrapping_add(i32::from(ret == 1));
                    }
                    ret = 0;
                    Ok(())
                })();
                inchi_free(heap, atom_pairs)?;
                pair_processing?;
                UnMarkOtherIndicators(&mut atoms, component.num_at)?;
                if ret < 0 {
                    heap.slice_mut(component.at)?[..component_atom_count].clone_from_slice(&atoms);
                    heap.slice_mut(derivatives)?[..component_atom_count]
                        .clone_from_slice(&derivative_values);
                    break 'components;
                }
            }
            heap.slice_mut(component.at)?[..component_atom_count].clone_from_slice(&atoms);
            heap.slice_mut(derivatives)?[..component_atom_count]
                .clone_from_slice(&derivative_values);
        }

        if num_cuts != 0 {
            let mut current_values = heap
                .slice(current.as_const())?
                .get(
                    ..usize::try_from(num_components)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            OAD_Edit_MergeComponentsAndRecreateOAD(
                heap,
                original,
                &mut current_values,
                num_components,
                &mut ret,
            )?;
            heap.slice_mut(current)?[..current_values.len()].clone_from_slice(&current_values);
        }
        Ok(())
    })();

    let cleanup = (|| -> Result<(), SourceHeapError> {
        if !derivatives.is_null() {
            inchi_free(heap, derivatives)?;
            derivatives = SourceMutPointer::null();
        }
        if !current.is_null() {
            for component_index in 0..num_components {
                let pointer = current.offset(i64::from(component_index))?;
                let mut component = heap
                    .slice(pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                FreeInpAtomData(heap, &mut component)?;
                *heap
                    .slice_mut(pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
            }
            inchi_free(heap, current)?;
            current = SourceMutPointer::null();
        }
        Ok(())
    })();
    processing?;
    cleanup?;
    Ok(if ret != 0 { ret } else { num_cuts })
}

pub(crate) fn FreeInpAtomData(
    heap: &mut SourceHeap,
    input: &mut INP_ATOM_DATA,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1058 FreeInpAtomData
    // INCHI✔️❌: void FreeInpAtomData(INP_ATOM_DATA *inp_at_data)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (inp_at_data)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         FreeInpAtom(&inp_at_data->at);
    // INCHI✔️❌:         FreeInpAtom(&inp_at_data->at_fixed_bonds);
    // INCHI✔️❌:         memset(inp_at_data, 0, sizeof(*inp_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeInpAtomData

    // BEGIN INCHI C HELPER FRAME: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1043 FreeInpAtom
    // INCHI✔️❌: void FreeInpAtom(inp_ATOM **at)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (at && *at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(*at);
    // INCHI✔️❌:         *at = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C HELPER FRAME: FreeInpAtom

    if !input.at.is_null() {
        inchi_free(heap, input.at)?;
        input.at = SourceMutPointer::null();
    }
    if !input.at_fixed_bonds.is_null() {
        inchi_free(heap, input.at_fixed_bonds)?;
        input.at_fixed_bonds = SourceMutPointer::null();
    }
    *input = INP_ATOM_DATA::default();
    Ok(())
}

pub(crate) fn free_underiv_temp_data(
    heap: &mut SourceHeap,
    atom_pairs: SourceMutPointer<R2C_ATPAIR>,
    derivatives: SourceMutPointer<DERIV_AT>,
    atoms2: SourceMutPointer<inp_ATOM>,
    current: SourceMutPointer<INP_ATOM_DATA>,
    num_components: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7632 free_underiv_temp_data
    // INCHI✔️❌: void free_underiv_temp_data( R2C_ATPAIR *ap,
    // INCHI✔️❌:                              DERIV_AT *da,
    // INCHI✔️❌:                              inp_ATOM *at2,
    // INCHI✔️❌:                              INP_ATOM_DATA *inp_cur_data,
    // INCHI✔️❌:                              int num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i_component;
    // INCHI✔️❌:     if (ap)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( ap );
    // INCHI✔️❌:         ap = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (da)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( da );
    // INCHI✔️❌:         da = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #ifdef FIX_UNDERIV_TO_SDF
    // INCHI✔️❌:     if (at2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( at2 );
    // INCHI✔️❌:         at2 = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     for (i_component = 0; i_component < num_components; i_component++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         FreeInpAtomData( inp_cur_data + i_component );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_free( inp_cur_data );
    // INCHI✔️❌:     inp_cur_data = NULL;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: free_underiv_temp_data

    if !atom_pairs.is_null() {
        inchi_free(heap, atom_pairs)?;
    }
    if !derivatives.is_null() {
        inchi_free(heap, derivatives)?;
    }
    if !atoms2.is_null() {
        inchi_free(heap, atoms2)?;
    }
    let mut component_index = 0_i32;
    while component_index < num_components {
        let component_pointer = current.offset(i64::from(component_index))?;
        let mut component = heap
            .slice(component_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        FreeInpAtomData(heap, &mut component)?;
        *heap
            .slice_mut(component_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
        component_index += 1;
    }
    if !current.is_null() {
        inchi_free(heap, current)?;
    }
    Ok(())
}

pub(crate) fn remove_cut_derivs(
    heap: &mut SourceHeap,
    mut num_atoms: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    current: &mut [INP_ATOM_DATA],
    component_index: i32,
    error_code: &mut i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7666 remove_cut_derivs
    // INCHI✔️❌: void remove_cut_derivs( int num_atoms,
    // INCHI✔️❌:                         inp_ATOM *at,
    // INCHI✔️❌:                         INP_ATOM_DATA *inp_cur_data,
    // INCHI✔️❌:                         int i_component,
    // INCHI✔️❌:                         int *errcode )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* Remove marked with Tritium disconnected derivative attachments */
    // INCHI✔️❌:     ORIG_ATOM_DATA Orig_inp_data1, *orig_inp_data1;
    // INCHI✔️❌:     INP_ATOM_DATA *inp_cur_data1 = NULL;
    // INCHI✔️❌:     int i, num_components1, i_component1, cur_num_at = 0; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_inp_data1 = &Orig_inp_data1;
    // INCHI✔️❌:     memset( orig_inp_data1, 0, sizeof( orig_inp_data1[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     UnMarkRingSystemsInp( at, num_atoms );
    // INCHI✔️❌:     UnMarkOtherIndicators( at, num_atoms );
    // INCHI✔️❌:     UnMarkOneComponent( at, num_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         orig_inp_data1->num_inp_bonds += at[i].valence;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     orig_inp_data1->num_inp_bonds /= 2;
    // INCHI✔️❌:     orig_inp_data1->num_inp_atoms = num_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_inp_data1->at = at; /* = from inp_cur_data[i_component].at */
    // INCHI✔️❌:
    // INCHI✔️❌:     num_components1 = MarkDisconnectedComponents( orig_inp_data1, 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:     inp_cur_data1 = (INP_ATOM_DATA *) inchi_calloc( num_components1, sizeof( inp_cur_data1[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (inp_cur_data1 && (num_components1 > 0)) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:     {
    // INCHI✔️❌:     /* Extract components and discard disconnected derivatizing agents */
    // INCHI✔️❌:     for (i_component1 = 0; i_component1 < num_components1; i_component1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         CreateInpAtomData( inp_cur_data1 + i_component1, orig_inp_data1->nCurAtLen[i_component1], 0 );
    // INCHI✔️❌:         inp_cur_data1[i_component1].num_at = ExtractConnectedComponent( orig_inp_data1->at, orig_inp_data1->num_inp_atoms,
    // INCHI✔️❌:                                                                         i_component1 + 1, inp_cur_data1[i_component1].at );
    // INCHI✔️❌:         /*  error processing */
    // INCHI✔️❌:         if (inp_cur_data1[i_component1].num_at <= 0 || orig_inp_data1->nCurAtLen[i_component1] != inp_cur_data1[i_component1].num_at)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *errcode = -( i_component1 + 1 ); /* severe error */
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* if the component has tritium then discard: it is a derivatizing agent */
    // INCHI✔️❌:         for (i = 0; i < inp_cur_data1[i_component1].num_at; i++)
    // INCHI✔️❌:         {
    // INCHI❌❌: #ifdef UNDERIV_ADD_D_TO_PRECURSOR
    // INCHI❌❌:             if (inp_cur_data1[i_component1].at[i].num_iso_H[1])
    // INCHI❌❌:             {
    // INCHI❌❌:                 inp_cur_data1[i_component1].at[i].num_iso_H[1] = 0; /* remove deuterium */
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:             if (inp_cur_data1[i_component1].at[i].num_iso_H[2])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 FreeInpAtomData( inp_cur_data1 + i_component1 );
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* Merge components into one -- must be only one */
    // INCHI✔️❌:     for (i_component1 = 0, num_atoms = 0; i_component1 < num_components1; i_component1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_atoms += inp_cur_data1[i_component1].num_at;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     at = (inp_ATOM *) inchi_calloc( num_atoms, sizeof( at[0] ) );
    // INCHI✔️❌:     cur_num_at = 0;
    // INCHI✔️❌:     for (i_component1 = 0; i_component1 < num_components1; i_component1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* clean and prepare */
    // INCHI✔️❌:         if (!inp_cur_data1[i_component1].num_at)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* removed derivatizing object */
    // INCHI✔️❌:                       /*UnMarkOneComponent( inp_cur_data1[i_component1].at, inp_cur_data1[i_component1].num_at );*/
    // INCHI✔️❌:                       /* merge one by one */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         cur_num_at = add_inp_ATOM( at, num_atoms, cur_num_at, inp_cur_data1[i_component1].at, inp_cur_data1[i_component1].num_at );
    // INCHI✔️❌:         FreeInpAtomData( inp_cur_data1 + i_component1 ); /* cleanup */
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Replace the component */
    // INCHI✔️❌:     /* Order of the following two statements is critically important */
    // INCHI✔️❌:
    // INCHI✔️❌:     UnMarkDisconnectedComponents( orig_inp_data1 ); /* orig_inp_data1->at is same as inp_cur_data[i_component].at */
    // INCHI✔️❌:     FreeInpAtomData( inp_cur_data + i_component ); /* cleanup the original component */
    // INCHI✔️❌:
    // INCHI✔️❌:     inp_cur_data[i_component].at = at;
    // INCHI✔️❌:     inp_cur_data[i_component].num_at = cur_num_at;
    // INCHI✔️❌:     inchi_free( inp_cur_data1 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: remove_cut_derivs

    let component_index =
        usize::try_from(component_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if component_index >= current.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let mut original = ORIG_ATOM_DATA {
        at: atoms,
        num_inp_atoms: num_atoms,
        ..ORIG_ATOM_DATA::default()
    };
    UnMarkRingSystemsInp(heap, atoms, num_atoms)?;
    {
        let source = heap.slice_mut(atoms)?;
        UnMarkOtherIndicators(source, num_atoms)?;
        UnMarkOneComponent(source, num_atoms)?;
        for atom in source
            .get(..usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            original.num_inp_bonds = original.num_inp_bonds.wrapping_add(i32::from(atom.valence));
        }
    }
    original.num_inp_bonds /= 2;
    let num_components = MarkDisconnectedComponents(heap, &mut original, 0)?;

    let component_data = match u64::try_from(num_components) {
        Ok(count) => match inchi_calloc::<INP_ATOM_DATA>(heap, count, 96) {
            Ok(pointer) => pointer,
            Err(
                SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange,
            ) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        },
        Err(_) => SourceMutPointer::null(),
    };

    let mut merged_atoms = atoms;
    let mut merged_count = 0_i32;
    if !component_data.is_null() && num_components > 0 {
        let mut component_index1 = 0_i32;
        while component_index1 < num_components {
            let component_data_pointer = component_data.offset(i64::from(component_index1))?;
            let component_len = heap
                .slice(original.nCurAtLen.as_const())?
                .get(
                    usize::try_from(component_index1)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut component = heap
                .slice(component_data_pointer.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let _ = CreateInpAtomData(heap, &mut component, i32::from(component_len), 0)?;
            *heap
                .slice_mut(component_data_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component;
            let output_pointer = heap
                .slice(component_data_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .at;
            if output_pointer.is_null() {
                return Err(SourceHeapError::AllocationFailed);
            }
            let source_atoms = heap
                .slice(original.at.as_const())?
                .get(
                    ..usize::try_from(num_atoms)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            let output_len = heap.slice(output_pointer.as_const())?.len();
            let mut output_atoms = vec![inp_ATOM::default(); output_len];
            let extracted = ExtractConnectedComponent(
                heap,
                &source_atoms,
                num_atoms,
                component_index1 + 1,
                &mut output_atoms,
            )?;
            heap.slice_mut(output_pointer)?[..output_len].clone_from_slice(&output_atoms);
            heap.slice_mut(component_data_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .num_at = extracted;
            if extracted <= 0 || i32::from(component_len) != extracted {
                *error_code = -(component_index1 + 1);
                break;
            }
            let mut atom_index = 0_i32;
            while atom_index < extracted {
                let atom = heap
                    .slice(output_pointer.as_const())?
                    .get(
                        usize::try_from(atom_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if atom.num_iso_H[2] != 0 {
                    let mut discarded = heap
                        .slice_mut(component_data_pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    FreeInpAtomData(heap, &mut discarded)?;
                    *heap
                        .slice_mut(component_data_pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = discarded;
                    break;
                }
                atom_index += 1;
            }
            component_index1 += 1;
        }

        num_atoms = 0;
        for component_index1 in 0..num_components {
            let component = heap
                .slice(
                    component_data
                        .offset(i64::from(component_index1))?
                        .as_const(),
                )?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            num_atoms = num_atoms.wrapping_add(component.num_at);
        }
        merged_atoms = match u64::try_from(num_atoms) {
            Ok(count) => match inchi_calloc::<inp_ATOM>(heap, count, 176) {
                Ok(pointer) => pointer,
                Err(
                    SourceHeapError::AllocationFailed
                    | SourceHeapError::AllocationSizeOverflow
                    | SourceHeapError::AllocationElementCountOutOfRange,
                ) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            },
            Err(_) => SourceMutPointer::null(),
        };
        let mut component_index1 = 0_i32;
        while component_index1 < num_components {
            let component_data_pointer = component_data.offset(i64::from(component_index1))?;
            let component = heap
                .slice(component_data_pointer.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if component.num_at != 0 {
                if merged_atoms.is_null() {
                    return Err(SourceHeapError::AllocationFailed);
                }
                let component_atoms = heap
                    .slice(component.at.as_const())?
                    .get(
                        ..usize::try_from(component.num_at)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                merged_count = add_inp_ATOM(
                    heap.slice_mut(merged_atoms)?,
                    num_atoms,
                    merged_count,
                    &component_atoms,
                    component.num_at,
                )?;
                let mut component_copy = component.clone();
                FreeInpAtomData(heap, &mut component_copy)?;
                *heap
                    .slice_mut(component_data_pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = component_copy;
            }
            component_index1 += 1;
        }
    }

    UnMarkDisconnectedComponents(heap, &mut original)?;
    let input = current
        .get_mut(component_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    FreeInpAtomData(heap, input)?;
    input.at = merged_atoms;
    input.num_at = merged_count;
    if !component_data.is_null() {
        inchi_free(heap, component_data)?;
    }
    Ok(())
}

fn work_get<T: Copy + 'static>(
    heap: &SourceHeap,
    pointer: SourceConstPointer<T>,
    index: i32,
) -> Result<T, SourceHeapError> {
    heap.slice(pointer.offset(i64::from(index))?)?
        .first()
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

fn work_set<T: 'static>(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
    value: T,
) -> Result<(), SourceHeapError> {
    *heap
        .slice_mut(pointer.offset(i64::from(index))?)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn MarkRingSystemsInp(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    start: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:59 MarkRingSystemsInp
    // INCHI✔️❌: int MarkRingSystemsInp( inp_ATOM *at, int num_atoms, int start )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB   *nStackAtom = NULL;
    // INCHI✔️❌:     int        nTopStackAtom = -1;
    // INCHI✔️❌:     AT_NUMB   *nRingStack = NULL;
    // INCHI✔️❌:     int        nTopRingStack = -1; /* was AT_NUMB */
    // INCHI✔️❌:     AT_NUMB   *nDfsNumber = NULL;
    // INCHI✔️❌:     AT_NUMB   *nLowNumber = NULL;
    // INCHI✔️❌:     S_CHAR    *cNeighNumb = NULL;
    // INCHI✔️❌:     AT_NUMB    nDfs;
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:     AT_NUMB    nRs, *nRsConnect = NULL;
    // INCHI❌❌:     int        k;
    // INCHI❌❌:     AT_NUMB   *tree = NULL;
    // INCHI❌❌:     int        nNumConnect, nMaxNumConnect, nLenConnect;
    // INCHI❌❌: #endif
    // INCHI✔️❌:     AT_NUMB    nNumAtInRingSystem;
    // INCHI✔️❌:     int        i, j, u, /*start,*/ nNumRingSystems, nNumStartChildren;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  allocate arrays */
    // INCHI✔️❌:     nStackAtom = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nStackAtom[0] ) );
    // INCHI✔️❌:     nRingStack = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nRingStack[0] ) );
    // INCHI✔️❌:     nDfsNumber = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nDfsNumber[0] ) );
    // INCHI✔️❌:     nLowNumber = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nLowNumber[0] ) );
    // INCHI✔️❌:     cNeighNumb = (S_CHAR  *) inchi_malloc( num_atoms * sizeof( cNeighNumb[0] ) );
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:     nRsConnect = (AT_NUMB *) inchi_calloc( 3 * num_atoms + 3, sizeof( nRsConnect[0] ) );
    // INCHI❌❌: #endif
    // INCHI✔️❌:     /*  check allocation */
    // INCHI✔️❌:     if (!nStackAtom || !nRingStack || !nDfsNumber || !nLowNumber || !cNeighNumb
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:          || !nRsConnect
    // INCHI❌❌: #endif
    // INCHI✔️❌:          )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nNumRingSystems = CT_OUT_OF_RAM;  /*  program error */ /*   <BRKPT> */
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /********************************************
    // INCHI✔️❌:     *
    // INCHI✔️❌:     * Find Cut-vertices & Blocks
    // INCHI✔️❌:     *
    // INCHI✔️❌:     ********************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  initiation */
    // INCHI✔️❌:     /*start           = 0;*/
    // INCHI✔️❌:     nNumRingSystems = 0;
    // INCHI✔️❌:     u = start; /*  start atom */
    // INCHI✔️❌:     nDfs = 0;
    // INCHI✔️❌:     nTopStackAtom = -1;
    // INCHI✔️❌:     nTopRingStack = -1;
    // INCHI✔️❌:     memset( nDfsNumber, 0, num_atoms * sizeof( nDfsNumber[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( cNeighNumb, 0, num_atoms * sizeof( cNeighNumb[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     /*  push the start atom on the stack */
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #66720 */
    // INCHI✔️❌:     if (u <= num_atoms - 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nLowNumber[u] = nDfsNumber[u] = ++nDfs;
    // INCHI✔️❌:         nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️❌:         nRingStack[++nTopRingStack] = (AT_NUMB)u;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nNumRingSystems = CT_OVERFLOW;  /*  program error */ /*   <BRKPT> */
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumStartChildren = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /* advance */
    // INCHI✔️❌:     advance_block:
    // INCHI✔️❌:
    // INCHI✔️❌:         /*if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) )*/
    // INCHI✔️❌:         /* replaced due to missing sequence point */
    // INCHI✔️❌:         if (i = (int) nStackAtom[nTopStackAtom], j = (int) cNeighNumb[i], (int) at[i].valence > j)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cNeighNumb[i] ++;
    // INCHI✔️❌:             u = (int) at[i].neighbor[j];
    // INCHI✔️❌:             if (!nDfsNumber[u])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* tree edge, 1st visit -- advance */
    // INCHI✔️❌:                 nStackAtom[++nTopStackAtom] = (AT_NUMB) u;
    // INCHI✔️❌:                 nRingStack[++nTopRingStack] = (AT_NUMB) u;
    // INCHI✔️❌:                 nLowNumber[u] = nDfsNumber[u] = ++nDfs;
    // INCHI✔️❌:                 nNumStartChildren += ( i == start );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!nTopStackAtom || u != (int) nStackAtom[nTopStackAtom - 1])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  may comment out ? */
    // INCHI✔️❌:                     /* back edge: u is not a predecessor of i */
    // INCHI✔️❌:                     if (nDfsNumber[u] < nDfsNumber[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* Back edge, 1st visit: u is an ancestor of i. Compare */
    // INCHI✔️❌:                         if (nLowNumber[i] > nDfsNumber[u])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nLowNumber[i] = nDfsNumber[u];
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 } /*  may comment out ? */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             goto advance_block;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cNeighNumb[i] = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* back up */
    // INCHI✔️❌:         if (i != start)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             u = (int) nStackAtom[nTopStackAtom - 1]; /* predecessor of i */
    // INCHI✔️❌:             if (nLowNumber[i] >= nDfsNumber[u])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* output the block */
    // INCHI✔️❌:                 nNumRingSystems++;
    // INCHI✔️❌:                 at[u].nBlockSystem = nNumRingSystems;
    // INCHI✔️❌:                 if (u != start || nNumStartChildren > 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at[u].bCutVertex += 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 while (nTopRingStack >= 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j = nRingStack[nTopRingStack--];
    // INCHI✔️❌:                     at[j].nBlockSystem = nNumRingSystems; /*  mark the atom */
    // INCHI✔️❌:                     if (i == j)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nLowNumber[u] > nLowNumber[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* inherit */
    // INCHI✔️❌:                     nLowNumber[u] = nLowNumber[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } while (--nTopStackAtom >= 0);
    // INCHI✔️❌:
    // INCHI✔️❌:     /****************************************************************************
    // INCHI✔️❌:     *
    // INCHI✔️❌:     * Find Ring Systems
    // INCHI✔️❌:     * Including chain atoms X: A-X-B, where the bonds (of any kind) are bridges.
    // INCHI✔️❌:     *
    // INCHI✔️❌:     ****************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  initiation */
    // INCHI✔️❌:     /* start           = 0;*/
    // INCHI✔️❌:     nNumRingSystems = 0;
    // INCHI✔️❌:     u = start; /*  start atom */
    // INCHI✔️❌:     nDfs = 0;
    // INCHI✔️❌:     nTopStackAtom = -1;
    // INCHI✔️❌:     nTopRingStack = -1;
    // INCHI✔️❌:     memset( nDfsNumber, 0, num_atoms * sizeof( nDfsNumber[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( cNeighNumb, 0, num_atoms * sizeof( cNeighNumb[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     /*  push the start atom on the stack */
    // INCHI✔️❌:     nLowNumber[u] = nDfsNumber[u] = ++nDfs;
    // INCHI✔️❌:     nStackAtom[++nTopStackAtom] = (AT_NUMB) u;
    // INCHI✔️❌:     nRingStack[++nTopRingStack] = (AT_NUMB) u;
    // INCHI✔️❌:
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:     nNumConnect = nLenConnect = nMaxNumConnect = 0;
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* advance */
    // INCHI✔️❌:     advance_ring:
    // INCHI✔️❌:         /*if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) )*/
    // INCHI✔️❌:         /* replaced due to missing sequence point */
    // INCHI✔️❌:         if (i = (int) nStackAtom[nTopStackAtom], j = (int) cNeighNumb[i], (int) at[i].valence > j)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cNeighNumb[i] ++;
    // INCHI✔️❌:             u = (int) at[i].neighbor[j];
    // INCHI✔️❌:             if (!nDfsNumber[u])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* tree edge, 1st visit -- advance */
    // INCHI✔️❌:                 nStackAtom[++nTopStackAtom] = (AT_NUMB) u;
    // INCHI✔️❌:                 nRingStack[++nTopRingStack] = (AT_NUMB) u;
    // INCHI✔️❌:                 nLowNumber[u] = nDfsNumber[u] = ++nDfs;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!nTopStackAtom || u != (int) nStackAtom[nTopStackAtom - 1])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* back edge: u is not a predecessor of i */
    // INCHI✔️❌:                     if (nDfsNumber[u] < nDfsNumber[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* Back edge, 1st visit: u is ancestor of i. Compare */
    // INCHI✔️❌:                         if (nLowNumber[i] > nDfsNumber[u])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nLowNumber[i] = nDfsNumber[u];
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             goto advance_ring;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cNeighNumb[i] = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* back up */
    // INCHI✔️❌:         if (nDfsNumber[i] == nLowNumber[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  found a ring system */
    // INCHI✔️❌:             nNumRingSystems++;
    // INCHI✔️❌:             /*  unwind nRingStack[] down to i */
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:             nNumConnect = 2;
    // INCHI❌❌:             /* data structure: for each ring system nRsConnect[] contains:
    // INCHI❌❌:             * 1) nNumConnect+1 = (number of already discovered neighboring "ring systems" + 1)+1
    // INCHI❌❌:             * 2) nNumAtInRingSystem
    // INCHI❌❌:             * 3) (nNumConnect-1) numbers (IDs) of neighboring ring systems.
    // INCHI❌❌:             * BFS guarantees that each neighboring ring system is encountered only one time
    // INCHI❌❌:             * Number of all neighboring ring systems = (nNumConnect-1)+1 = nNumConnect
    // INCHI❌❌:             * (One additional ring system is where the BFS retracts from the vertex #i,
    // INCHI❌❌:             * except when i=DFS root node. In the latter case there is/are only (nNumConnect-1)
    // INCHI❌❌:             * neighboring ring system(s).
    // INCHI❌❌:             */
    // INCHI❌❌: #endif
    // INCHI✔️❌:             /*  count atoms in a ring system */
    // INCHI✔️❌:             for (nNumAtInRingSystem = 0, j = nTopRingStack; 0 <= j; j--)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nNumAtInRingSystem++;
    // INCHI✔️❌:                 if (i == (int) nRingStack[j])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             while (nTopRingStack >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 j = (int) nRingStack[nTopRingStack--];
    // INCHI✔️❌:                 at[j].nRingSystem = (AT_NUMB) nNumRingSystems; /*  ring system id */
    // INCHI✔️❌:                 at[j].nNumAtInRingSystem = nNumAtInRingSystem;
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:                 for (k = 0; k < at[j].valence; k++)
    // INCHI❌❌:                 {
    // INCHI❌❌:                     if (( nRs = at[at[j].neighbor[k]].nRingSystem ) && (int) nRs != nNumRingSystems)
    // INCHI❌❌:                     {
    // INCHI❌❌:                         nRsConnect[nLenConnect + ( nNumConnect++ )] = nRs; /*  adjacent ring system id */
    // INCHI❌❌:                     }
    // INCHI❌❌:                 }
    // INCHI❌❌: #endif
    // INCHI✔️❌:                 if (i == j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  reached atom on the top of nStackAtom[] stack  */
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:             nRsConnect[nLenConnect] = nNumConnect;
    // INCHI❌❌:             nRsConnect[nLenConnect + 1] = nNumAtInRingSystem;
    // INCHI❌❌:             nLenConnect += nNumConnect;
    // INCHI❌❌:             if (nMaxNumConnect < nNumConnect)
    // INCHI❌❌:             {
    // INCHI❌❌:                 /*  max number of neighboring ring systems */
    // INCHI❌❌:                 nMaxNumConnect = nNumConnect;
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nTopStackAtom > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 j = (int) nStackAtom[nTopStackAtom - 1];
    // INCHI✔️❌:                 /* inherit nLowNumber */
    // INCHI✔️❌:                 if (nLowNumber[j] > nLowNumber[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nLowNumber[j] = nLowNumber[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } while (--nTopStackAtom >= 0);
    // INCHI✔️❌:
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 ) /*  normally disabled */
    // INCHI❌❌:     nMaxNumConnect++;
    // INCHI❌❌:     if (nNumRingSystems > 1)
    // INCHI❌❌:     {
    // INCHI❌❌:         int nCol = nMaxNumConnect + 1;
    // INCHI❌❌:         int nNumInSyst = nMaxNumConnect;
    // INCHI❌❌:         int nMaxNeigh = nMaxNumConnect - 1;
    // INCHI❌❌: #define T(a,b) tree[(a)*nCol+b]
    // INCHI❌❌:         if (tree = (AT_NUMB *) inchi_calloc( nCol * ( nNumRingSystems + 1 ), sizeof( tree[0] ) ))
    // INCHI❌❌:         {
    // INCHI❌❌:             int len, neigh;
    // INCHI❌❌:             /*  reuse previous allocations */
    // INCHI❌❌:             AT_NUMB *nNumVisitedNeighbors = nStackAtom;
    // INCHI❌❌:             AT_NUMB *nDistanceFromTerminal = nRingStack;
    // INCHI❌❌:             AT_NUMB *nCurrActiveRingSystem = nDfsNumber;
    // INCHI❌❌:             AT_NUMB *nNextActiveRingSystem = nLowNumber;
    // INCHI❌❌:             int        nNumCurrActiveRingSystems, nNumNextActiveRingSystems, pass;
    // INCHI❌❌:             /* build a "condensation graph (actually, a tree)" in which
    // INCHI❌❌:             * each vertex corresponds to a ring system T(row, col) = T(ring syst, neighbors)
    // INCHI❌❌:             * Number of rows = column length = max. number of ring system neighbors + 2
    // INCHI❌❌:             * Number of cols = row length    = number of ring systems + 1
    // INCHI❌❌:             * Neighboring ring systems are contiguously stored in a row
    // INCHI❌❌:             * T(i,0) = number of neighbors,  1 <= i <= nNumRingSystems;
    // INCHI❌❌:             * T(i,k) = number of a neighboring ring system, 1 <= k <= T(i,0)
    // INCHI❌❌:             * T(i,nCol-1) = number of atoms in the system #i
    // INCHI❌❌:             */
    // INCHI❌❌:             for (i = 1, j = 0; len = nRsConnect[j]; i++)
    // INCHI❌❌:             {
    // INCHI❌❌:                 T( i, nNumInSyst ) = nRsConnect[j + 1];
    // INCHI❌❌:                 for (k = 2; k < len; k++)
    // INCHI❌❌:                 {
    // INCHI❌❌:                     neigh = nRsConnect[j + k];
    // INCHI❌❌:                     if (T( i, 0 ) < nMaxNeigh && T( neigh, 0 ) < nMaxNeigh)
    // INCHI❌❌:                     {
    // INCHI❌❌:                         T( i, 0 )++;
    // INCHI❌❌:                         T( neigh, 0 )++;
    // INCHI❌❌:                         T( i, T( i, 0 ) ) = neigh;
    // INCHI❌❌:                         T( neigh, T( neigh, 0 ) ) = i;
    // INCHI❌❌:                     }
    // INCHI❌❌:                     else
    // INCHI❌❌:                     {
    // INCHI❌❌:                         nNumRingSystems = CT_OVERFLOW;  /*  program error */ /*   <BRKPT> */
    // INCHI❌❌:                         goto exit_function;
    // INCHI❌❌:                     }
    // INCHI❌❌:                 }
    // INCHI❌❌:                 j += len;
    // INCHI❌❌:             }
    // INCHI❌❌:             /*  clear memory */
    // INCHI❌❌:             memset( nNumVisitedNeighbors, 0, nNumRingSystems * sizeof( nNumVisitedNeighbors[0] ) );
    // INCHI❌❌:             memset( nDistanceFromTerminal, 0, nNumRingSystems * sizeof( nDistanceFromTerminal[0] ) );
    // INCHI❌❌:             memset( nCurrActiveRingSystem, 0, nNumRingSystems * sizeof( nCurrActiveRingSystem[0] ) );
    // INCHI❌❌:             memset( nNextActiveRingSystem, 0, nNumRingSystems * sizeof( nNextActiveRingSystem[0] ) );
    // INCHI❌❌:             nNumNextActiveRingSystems = 0;
    // INCHI❌❌:             for (i = 0; i < nNumRingSystems; i++)
    // INCHI❌❌:             {
    // INCHI❌❌:                 if (1 == T( i + 1, 0 ))
    // INCHI❌❌:                 {
    // INCHI❌❌:                     nNextActiveRingSystem[i] = 1; /*  number of traversed neighbors + 1 */
    // INCHI❌❌:                     nDistanceFromTerminal[i] = 1;
    // INCHI❌❌:                     nNumNextActiveRingSystems++;
    // INCHI❌❌:                 }
    // INCHI❌❌:                 else
    // INCHI❌❌:                 {
    // INCHI❌❌:                     nNextActiveRingSystem[i] = 0;
    // INCHI❌❌:                     nDistanceFromTerminal[i] = 0;
    // INCHI❌❌:                 }
    // INCHI❌❌:                 nNumVisitedNeighbors[i] = 0;
    // INCHI❌❌:             }
    // INCHI❌❌:
    // INCHI❌❌:             /* nCurrActiveRingSystem[i] = a sum of:
    // INCHI❌❌:             * 1) +1 if it is or was active
    // INCHI❌❌:             * 2) +(number of neighbors from which it was reached)
    // INCHI❌❌:             * 3) +1 if it was left and not active anymore
    // INCHI❌❌:             */
    // INCHI❌❌:             pass = 0;
    // INCHI❌❌:             do
    // INCHI❌❌:             {
    // INCHI❌❌:                 nNumCurrActiveRingSystems = nNumNextActiveRingSystems;
    // INCHI❌❌:                 nNumNextActiveRingSystems = 0;
    // INCHI❌❌:                 memcpy( nCurrActiveRingSystem, nNextActiveRingSystem,
    // INCHI❌❌:                         nNumRingSystems * sizeof( nNextActiveRingSystem[0] ) );
    // INCHI❌❌:                 for (i = 0; i < nNumRingSystems; i++)
    // INCHI❌❌:                 {
    // INCHI❌❌:                     if (T( i + 1, 0 ) == nCurrActiveRingSystem[i])
    // INCHI❌❌:                     {
    // INCHI❌❌:                         /* on the previous pass currently active ring system i+1 bas been reached
    // INCHI❌❌:                         * from all neighbors except one;
    // INCHI❌❌:                         * the neighbors from which it was reached have
    // INCHI❌❌:                         * T(neigh,0)+1 == nCurrActiveRingSystem[i]
    // INCHI❌❌:                         * this ring system has not been left yet
    // INCHI❌❌:                         */
    // INCHI❌❌:                         for (k = 1, len = T( i + 1, 0 ); k <= len; k++)
    // INCHI❌❌:                         {
    // INCHI❌❌:                             neigh = (int) T( i + 1, k );
    // INCHI❌❌:                             if (T( neigh, 0 ) >= nCurrActiveRingSystem[neigh - 1])
    // INCHI❌❌:                             {
    // INCHI❌❌:                                 if (0 == pass)
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     nDistanceFromTerminal[i] = 1;
    // INCHI❌❌:                                 }
    // INCHI❌❌:                                 break;
    // INCHI❌❌:                             }
    // INCHI❌❌:                         }
    // INCHI❌❌:                         if (k <= len)
    // INCHI❌❌:                         {
    // INCHI❌❌:                             /* neigh was not reached from at least 2 neighbors
    // INCHI❌❌:                             * walk along -R- chain (T(neigh,0)==2) up to
    // INCHI❌❌:                             * 1)  a terminal system, not including it or
    // INCHI❌❌:                             * 2)  a branching point.
    // INCHI❌❌:                             *
    // INCHI❌❌:                             * pass = 0: started from terminal systems:
    // INCHI❌❌:                             *     reach the branching point.
    // INCHI❌❌:                             * If chain system next to a terminal system has already been reached
    // INCHI❌❌:                             * then walk along it according to Note below
    // INCHI❌❌:                             *
    // INCHI❌❌:                             * pass > 0: started from branching points
    // INCHI❌❌:                             * 2a) If the branching point has not been reached from 2 or more neighbors,
    // INCHI❌❌:                             *     then include it
    // INCHI❌❌:                             * 2b) If the branching point has not been reached from 1 neighbor only,
    // INCHI❌❌:                             *     then do not include it: it will be a starting point later
    // INCHI❌❌:                             * Note: if a chain atom already has nDistanceFromTerminal[i] > 0, then
    // INCHI❌❌:                             *     the last atom should be the one such that
    // INCHI❌❌:                             *     its nDistanceFromTerminal[]+1>= nDistanceFromTerminal[] of the
    // INCHI❌❌:                             *     next in the chain
    // INCHI❌❌:                             */
    // INCHI❌❌:                             int bOk = 0;
    // INCHI❌❌:                             k = i + 1; /*  starting point */
    // INCHI❌❌:                             if (0 == pass && T( k, nNumInSyst ) > 1)
    // INCHI❌❌:                             {
    // INCHI❌❌:                                 nNumNextActiveRingSystems++; /*  request next pass */
    // INCHI❌❌:                                 continue; /*  stop a the terminal ring system */
    // INCHI❌❌:                             }
    // INCHI❌❌:                             while (2 == T( neigh, 0 ))
    // INCHI❌❌:                             {
    // INCHI❌❌:                                 /*  walk along a chain */
    // INCHI❌❌:                                 if (!nNextActiveRingSystem[neigh - 1])
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     nNextActiveRingSystem[neigh - 1] = 1; /*  make neighbor active */
    // INCHI❌❌:                                 }
    // INCHI❌❌:                                 else
    // INCHI❌❌:                                     if (nDistanceFromTerminal[k - 1] + 1 <= nDistanceFromTerminal[neigh - 1])
    // INCHI❌❌:                                     {
    // INCHI❌❌:                                         /*  walking along the chain; already have had a walk */
    // INCHI❌❌:                                         /*  in the opposite direction at this pass */
    // INCHI❌❌:                                     }
    // INCHI❌❌:                                     else
    // INCHI❌❌:                                     {
    // INCHI❌❌:                                         /*  k is the last; neigh (it is a bridge -X-) has not been reached */
    // INCHI❌❌:                                         bOk = 1;
    // INCHI❌❌:                                         break;
    // INCHI❌❌:                                     }
    // INCHI❌❌:                                 nNextActiveRingSystem[k - 1] ++; /*  leave system k */
    // INCHI❌❌:                                 if (nNextActiveRingSystem[neigh - 1] < T( neigh, 0 ))
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     nNextActiveRingSystem[neigh - 1] ++; /*  add one connection to neigh */
    // INCHI❌❌:                                 }
    // INCHI❌❌:                                 nDistanceFromTerminal[neigh - 1] = nDistanceFromTerminal[k - 1] + 1;
    // INCHI❌❌:                                 j = ( T( neigh, 1 ) == k ) ? 2 : 1;
    // INCHI❌❌:                                 k = neigh;
    // INCHI❌❌:                                 neigh = T( k, j ); /*  next in the chain */
    // INCHI❌❌:                                 nNumNextActiveRingSystems++;
    // INCHI❌❌:                                 if (T( k, nNumInSyst ) > 1)
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     bOk = 1;
    // INCHI❌❌:                                     break; /*  stop on a ring system */
    // INCHI❌❌:                                 }
    // INCHI❌❌:                             }
    // INCHI❌❌:                             /*  neigh is a terminal or a bridge or a branching point */
    // INCHI❌❌:                             if (2 > T( neigh, 0 ))
    // INCHI❌❌:                             {
    // INCHI❌❌:                                 /*  neighbor is a terminal atom */
    // INCHI❌❌:                                 if (1 < pass)
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     nNumRingSystems = CT_UNKNOWN_ERR; /*  error (debug only) */ /*   <BRKPT> */
    // INCHI❌❌:                                     goto exit_function;
    // INCHI❌❌:                                 }
    // INCHI❌❌:                                 continue;
    // INCHI❌❌:                             }
    // INCHI❌❌:                             if (2 == T( neigh, 0 ))
    // INCHI❌❌:                             {
    // INCHI❌❌:                                 /*  neighbor is a bridge */
    // INCHI❌❌:                                 continue;
    // INCHI❌❌:                             }
    // INCHI❌❌:                             /*  neighbor is a branching point */
    // INCHI❌❌:                             if (T( neigh, 0 ) > nCurrActiveRingSystem[neigh - 1])
    // INCHI❌❌:                             {
    // INCHI❌❌:                                 /*  move to the neigh (make neigh active): on previous pass it */
    // INCHI❌❌:                                 /*  has not been reached from 2 or more neighbors */
    // INCHI❌❌:                                 if (!nNextActiveRingSystem[neigh - 1])
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     nNextActiveRingSystem[neigh - 1] = 1;
    // INCHI❌❌:                                 }
    // INCHI❌❌:                                 if (nDistanceFromTerminal[neigh - 1] < nDistanceFromTerminal[k - 1] + 1)
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     nDistanceFromTerminal[neigh - 1] = nDistanceFromTerminal[k - 1] + 1;
    // INCHI❌❌:                                 }
    // INCHI❌❌:                                 nNextActiveRingSystem[k - 1] ++; /*  leave system k */
    // INCHI❌❌:                                 if (nNextActiveRingSystem[neigh - 1] < T( neigh, 0 ))
    // INCHI❌❌:                                 {
    // INCHI❌❌:                                     nNextActiveRingSystem[neigh - 1] ++; /*  add one connection to neigh */
    // INCHI❌❌:                                 }
    // INCHI❌❌:                                 nNumNextActiveRingSystems++;
    // INCHI❌❌:                             }
    // INCHI❌❌:                         }
    // INCHI❌❌:                     }
    // INCHI❌❌:                 }
    // INCHI❌❌:                 pass++;
    // INCHI❌❌:             } while (nNumNextActiveRingSystems);
    // INCHI❌❌:
    // INCHI❌❌:             for (i = 0; i < num_atoms; i++)
    // INCHI❌❌:             {
    // INCHI❌❌:                 at[i].nDistanceFromTerminal = nDistanceFromTerminal[(int) at[i].nRingSystem - 1];
    // INCHI❌❌:             }
    // INCHI❌❌:
    // INCHI❌❌:             inchi_free( tree );
    // INCHI❌❌:             tree = NULL;
    // INCHI❌❌: #undef T
    // INCHI❌❌:         }
    // INCHI❌❌:         else
    // INCHI❌❌:         {
    // INCHI❌❌:             nNumRingSystems = CT_OUT_OF_RAM; /*  error */ /*   <BRKPT> */
    // INCHI❌❌:             goto exit_function;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nStackAtom)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nStackAtom );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nRingStack)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nRingStack );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nDfsNumber)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nDfsNumber );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nLowNumber)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nLowNumber );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (cNeighNumb)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( cNeighNumb );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    // INCHI❌❌:     if (nRsConnect)
    // INCHI❌❌:         inchi_free( nRsConnect );
    // INCHI❌❌:     if (tree)
    // INCHI❌❌:         inchi_free( tree );
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumRingSystems;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkRingSystemsInp
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:706
    // INCHI❌❌: #define FIND_RING_SYSTEMS           1  /* 1 => find and mark ring systems, blocks, cut-vertices */
    // INCHI❌❌: #define FIND_RINS_SYSTEMS_DISTANCES 0  /* 1 => find ring system and atom distance from terminal */
    // END INCHI ACTIVE MACRO CONFIGURATION: MarkRingSystemsInp

    let count = match usize::try_from(num_atoms) {
        Ok(count) => count,
        Err(_) => return Ok(CT_OUT_OF_RAM),
    };
    let stack = match heap.allocate(vec![0_u16; count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(CT_OUT_OF_RAM),
        Err(error) => return Err(error),
    };
    let ring_stack = match heap.allocate(vec![0_u16; count]) {
        Ok(pointer) => pointer,
        Err(error) => {
            inchi_free(heap, stack)?;
            return match error {
                SourceHeapError::AllocationFailed => Ok(CT_OUT_OF_RAM),
                _ => Err(error),
            };
        }
    };
    let dfs_number = match heap.allocate(vec![0_u16; count]) {
        Ok(pointer) => pointer,
        Err(error) => {
            inchi_free(heap, stack)?;
            inchi_free(heap, ring_stack)?;
            return match error {
                SourceHeapError::AllocationFailed => Ok(CT_OUT_OF_RAM),
                _ => Err(error),
            };
        }
    };
    let low_number = match heap.allocate(vec![0_u16; count]) {
        Ok(pointer) => pointer,
        Err(error) => {
            inchi_free(heap, stack)?;
            inchi_free(heap, ring_stack)?;
            inchi_free(heap, dfs_number)?;
            return match error {
                SourceHeapError::AllocationFailed => Ok(CT_OUT_OF_RAM),
                _ => Err(error),
            };
        }
    };
    let neighbor_number = match heap.allocate(vec![0_i8; count]) {
        Ok(pointer) => pointer,
        Err(error) => {
            inchi_free(heap, stack)?;
            inchi_free(heap, ring_stack)?;
            inchi_free(heap, dfs_number)?;
            inchi_free(heap, low_number)?;
            return match error {
                SourceHeapError::AllocationFailed => Ok(CT_OUT_OF_RAM),
                _ => Err(error),
            };
        }
    };

    let result = (|| -> Result<i32, SourceHeapError> {
        let mut ring_systems = 0_i32;
        let mut u = start;
        let mut dfs = 0_u16;
        let mut stack_top = -1_i32;
        let mut ring_top = -1_i32;
        if u <= num_atoms.wrapping_sub(1) {
            dfs = dfs.wrapping_add(1);
            work_set(heap, low_number, u, dfs)?;
            work_set(heap, dfs_number, u, dfs)?;
            stack_top += 1;
            work_set(heap, stack, stack_top, u as u16)?;
            ring_top += 1;
            work_set(heap, ring_stack, ring_top, u as u16)?;
        } else {
            return Ok(CT_OVERFLOW);
        }
        let mut start_children = 0_i32;

        loop {
            let i = i32::from(work_get(heap, stack.as_const(), stack_top)?);
            let j = i32::from(work_get(heap, neighbor_number.as_const(), i)?);
            let atom = heap
                .slice(atoms.as_const().offset(i64::from(i))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if i32::from(atom.valence) > j {
                work_set(heap, neighbor_number, i, (j as i8).wrapping_add(1))?;
                let neighbor_index =
                    usize::try_from(j).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                u = i32::from(
                    *atom
                        .neighbor
                        .get(neighbor_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                if work_get(heap, dfs_number.as_const(), u)? == 0 {
                    stack_top += 1;
                    work_set(heap, stack, stack_top, u as u16)?;
                    ring_top += 1;
                    work_set(heap, ring_stack, ring_top, u as u16)?;
                    dfs = dfs.wrapping_add(1);
                    work_set(heap, low_number, u, dfs)?;
                    work_set(heap, dfs_number, u, dfs)?;
                    start_children = start_children.wrapping_add(i32::from(i == start));
                } else if stack_top == 0
                    || u != i32::from(work_get(heap, stack.as_const(), stack_top - 1)?)
                {
                    let neighbor_dfs = work_get(heap, dfs_number.as_const(), u)?;
                    let current_dfs = work_get(heap, dfs_number.as_const(), i)?;
                    if neighbor_dfs < current_dfs
                        && work_get(heap, low_number.as_const(), i)? > neighbor_dfs
                    {
                        work_set(heap, low_number, i, neighbor_dfs)?;
                    }
                }
                continue;
            }
            work_set(heap, neighbor_number, i, 0)?;
            if i != start {
                u = i32::from(work_get(heap, stack.as_const(), stack_top - 1)?);
                if work_get(heap, low_number.as_const(), i)?
                    >= work_get(heap, dfs_number.as_const(), u)?
                {
                    ring_systems = ring_systems.wrapping_add(1);
                    {
                        let atom = heap
                            .slice_mut(atoms.offset(i64::from(u))?)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        atom.nBlockSystem = ring_systems as u16;
                        if u != start || start_children > 1 {
                            atom.bCutVertex = atom.bCutVertex.wrapping_add(1);
                        }
                    }
                    while ring_top >= 0 {
                        let j = i32::from(work_get(heap, ring_stack.as_const(), ring_top)?);
                        ring_top -= 1;
                        heap.slice_mut(atoms.offset(i64::from(j))?)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nBlockSystem = ring_systems as u16;
                        if i == j {
                            break;
                        }
                    }
                } else {
                    let child_low = work_get(heap, low_number.as_const(), i)?;
                    if work_get(heap, low_number.as_const(), u)? > child_low {
                        work_set(heap, low_number, u, child_low)?;
                    }
                }
            }
            stack_top -= 1;
            if stack_top < 0 {
                break;
            }
        }

        ring_systems = 0;
        u = start;
        dfs = 1;
        stack_top = 0;
        ring_top = 0;
        heap.slice_mut(dfs_number)?.fill(0);
        heap.slice_mut(neighbor_number)?.fill(0);
        work_set(heap, low_number, u, dfs)?;
        work_set(heap, dfs_number, u, dfs)?;
        work_set(heap, stack, stack_top, u as u16)?;
        work_set(heap, ring_stack, ring_top, u as u16)?;

        loop {
            let i = i32::from(work_get(heap, stack.as_const(), stack_top)?);
            let j = i32::from(work_get(heap, neighbor_number.as_const(), i)?);
            let atom = heap
                .slice(atoms.as_const().offset(i64::from(i))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if i32::from(atom.valence) > j {
                work_set(heap, neighbor_number, i, (j as i8).wrapping_add(1))?;
                let neighbor_index =
                    usize::try_from(j).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                u = i32::from(
                    *atom
                        .neighbor
                        .get(neighbor_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                if work_get(heap, dfs_number.as_const(), u)? == 0 {
                    stack_top += 1;
                    work_set(heap, stack, stack_top, u as u16)?;
                    ring_top += 1;
                    work_set(heap, ring_stack, ring_top, u as u16)?;
                    dfs = dfs.wrapping_add(1);
                    work_set(heap, low_number, u, dfs)?;
                    work_set(heap, dfs_number, u, dfs)?;
                } else if stack_top == 0
                    || u != i32::from(work_get(heap, stack.as_const(), stack_top - 1)?)
                {
                    let neighbor_dfs = work_get(heap, dfs_number.as_const(), u)?;
                    let current_dfs = work_get(heap, dfs_number.as_const(), i)?;
                    if neighbor_dfs < current_dfs
                        && work_get(heap, low_number.as_const(), i)? > neighbor_dfs
                    {
                        work_set(heap, low_number, i, neighbor_dfs)?;
                    }
                }
                continue;
            }
            work_set(heap, neighbor_number, i, 0)?;
            if work_get(heap, dfs_number.as_const(), i)?
                == work_get(heap, low_number.as_const(), i)?
            {
                ring_systems = ring_systems.wrapping_add(1);
                let mut ring_size = 0_u16;
                let mut j = ring_top;
                while j >= 0 {
                    ring_size = ring_size.wrapping_add(1);
                    if i == i32::from(work_get(heap, ring_stack.as_const(), j)?) {
                        break;
                    }
                    j -= 1;
                }
                while ring_top >= 0 {
                    let j = i32::from(work_get(heap, ring_stack.as_const(), ring_top)?);
                    ring_top -= 1;
                    let atom = heap
                        .slice_mut(atoms.offset(i64::from(j))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    atom.nRingSystem = ring_systems as u16;
                    atom.nNumAtInRingSystem = ring_size;
                    if i == j {
                        break;
                    }
                }
            } else if stack_top > 0 {
                let predecessor = i32::from(work_get(heap, stack.as_const(), stack_top - 1)?);
                let child_low = work_get(heap, low_number.as_const(), i)?;
                if work_get(heap, low_number.as_const(), predecessor)? > child_low {
                    work_set(heap, low_number, predecessor, child_low)?;
                }
            }
            stack_top -= 1;
            if stack_top < 0 {
                break;
            }
        }
        Ok(ring_systems)
    })();

    inchi_free(heap, stack)?;
    inchi_free(heap, ring_stack)?;
    inchi_free(heap, dfs_number)?;
    inchi_free(heap, low_number)?;
    inchi_free(heap, neighbor_number)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn mark_atoms_cFlags(
    atoms: &mut [inp_ATOM],
    start: i32,
    mut num: i32,
    flags: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1368 mark_atoms_cFlags
    // INCHI✔️❌: int mark_atoms_cFlags( inp_ATOM *at, int start, int num, char cFlags )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (!( at[start].cFlags & cFlags ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i;
    // INCHI✔️❌:         at[start].cFlags |= cFlags;
    // INCHI✔️❌:         num++;
    // INCHI✔️❌:         for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = mark_atoms_cFlags( at, at[start].neighbor[i], num, cFlags );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num; /* number of atoms traversed forward from at[start] */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: mark_atoms_cFlags

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(start)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.cFlags & flags == 0 {
        let valence = usize::try_from(i32::from(atom.valence))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        atoms[start].cFlags |= flags;
        num = num.wrapping_add(1);
        for neighbor_order in 0..valence {
            let neighbor = i32::from(atoms[start].neighbor[neighbor_order]);
            num = mark_atoms_cFlags(atoms, neighbor, num, flags)?;
        }
    }
    Ok(num)
}

#[allow(non_snake_case)]
pub(crate) fn unmark_atoms_cFlags(
    atoms: &mut [inp_ATOM],
    start: i32,
    mut num: i32,
    flags: i8,
    inverse_flags: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1388 unmark_atoms_cFlags
    // INCHI✔️❌: int unmark_atoms_cFlags( inp_ATOM *at,
    // INCHI✔️❌:                          int start,
    // INCHI✔️❌:                          int num,
    // INCHI✔️❌:                          char cFlags,
    // INCHI✔️❌:                          char cInvFlags )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (at[start].cFlags & cFlags)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i;
    // INCHI✔️❌:         at[start].cFlags &= cInvFlags;
    // INCHI✔️❌:         num++;
    // INCHI✔️❌:         for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = unmark_atoms_cFlags( at, at[start].neighbor[i], num, cFlags, cInvFlags );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num; /* number of atoms traversed forward from at[start] */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: unmark_atoms_cFlags

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(start)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.cFlags & flags != 0 {
        let valence = usize::try_from(i32::from(atom.valence))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        atoms[start].cFlags &= inverse_flags;
        num = num.wrapping_add(1);
        for neighbor_order in 0..valence {
            let neighbor = i32::from(atoms[start].neighbor[neighbor_order]);
            num = unmark_atoms_cFlags(atoms, neighbor, num, flags, inverse_flags)?;
        }
    }
    Ok(num)
}

#[allow(non_snake_case)]
pub(crate) fn is_C_or_S_DB_O(atoms: &[inp_ATOM], atom_index: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1410 is_C_or_S_DB_O
    // INCHI✔️❌: int is_C_or_S_DB_O( inp_ATOM *at, int i )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, neigh;
    // INCHI✔️❌:     if ((at[i].el_number != EL_NUMBER_C &&
    // INCHI✔️❌:          at[i].el_number != EL_NUMBER_S) ||
    // INCHI✔️❌:          at[i].charge || at[i].radical) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[i].neighbor[j];
    // INCHI✔️❌:         if (( at[neigh].el_number == EL_NUMBER_O ||
    // INCHI✔️❌:               at[neigh].el_number == EL_NUMBER_S ) &&
    // INCHI✔️❌:              !at[neigh].num_H && 1 == at[neigh].valence &&
    // INCHI✔️❌:              2 == at[neigh].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_C_or_S_DB_O
    // BEGIN ACTIVE INCHI HEADER MACROS: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h
    // INCHI✔️❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️❌: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️❌: #define EL_NUMBER_S  ((U_CHAR)16)
    // END ACTIVE INCHI HEADER MACROS: EL_NUMBER_*

    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if (atom.el_number != EL_NUMBER_C as u8 && atom.el_number != EL_NUMBER_S as u8)
        || atom.charge != 0
        || atom.radical != 0
    {
        return Ok(0);
    }
    let valence = usize::try_from(i32::from(atom.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for neighbor_order in 0..valence {
        let neighbor = atoms
            .get(atom.neighbor[neighbor_order] as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if (neighbor.el_number == EL_NUMBER_O as u8 || neighbor.el_number == EL_NUMBER_S as u8)
            && neighbor.num_H == 0
            && neighbor.valence == 1
            && neighbor.chem_bonds_valence == 2
        {
            return Ok(1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn is_C_DB_O(atoms: &[inp_ATOM], atom_index: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1434 is_C_DB_O
    // INCHI✔️❌: int is_C_DB_O( inp_ATOM *at, int i )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, neigh;
    // INCHI✔️❌:     if (at[i].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:          at[i].charge || at[i].radical ||
    // INCHI✔️❌:          at[i].valence != 3 || at[i].chem_bonds_valence != 4)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[i].neighbor[j];
    // INCHI✔️❌:         if (( at[neigh].el_number == EL_NUMBER_O ) &&
    // INCHI✔️❌:              !at[neigh].num_H && 1 == at[neigh].valence &&
    // INCHI✔️❌:              2 == at[neigh].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return j + 1; /* =O ord */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_C_DB_O
    // BEGIN ACTIVE INCHI HEADER MACROS: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h
    // INCHI✔️❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️❌: #define EL_NUMBER_O  ((U_CHAR)8)
    // END ACTIVE INCHI HEADER MACROS: EL_NUMBER_C, EL_NUMBER_O

    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.el_number != EL_NUMBER_C
        || atom.charge != 0
        || atom.radical != 0
        || atom.valence != 3
        || atom.chem_bonds_valence != 4
    {
        return Ok(0);
    }
    let valence = usize::try_from(i32::from(atom.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for neighbor_order in 0..valence {
        let neighbor = atoms
            .get(atom.neighbor[neighbor_order] as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.el_number == EL_NUMBER_O
            && neighbor.num_H == 0
            && neighbor.valence == 1
            && neighbor.chem_bonds_valence == 2
        {
            return Ok(i32::try_from(neighbor_order).unwrap() + 1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn is_C_unsat_not_arom(
    atoms: &[inp_ATOM],
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1457 is_C_unsat_not_arom
    // INCHI✔️❌: int is_C_unsat_not_arom( inp_ATOM *at, int i )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, neigh, num_arom, num_DB;
    // INCHI✔️❌:     if (at[i].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:          at[i].valence == at[i].chem_bonds_valence || /* no double/triple bonds */
    // INCHI✔️❌:          at[i].valence + 1 < at[i].chem_bonds_valence || /* >1 double bond or >=1 triple bond */
    // INCHI✔️❌:          at[i].chem_bonds_valence + at[i].num_H != 4 || /* C has wrong valence */
    // INCHI✔️❌:          at[i].charge || at[i].radical)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     num_arom = num_DB = 0;
    // INCHI✔️❌:     for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[i].neighbor[j];
    // INCHI✔️❌:         num_arom += at[i].bond_type[j] == BOND_TYPE_ALTERN;
    // INCHI✔️❌:         if (( at[neigh].el_number == EL_NUMBER_O ||
    // INCHI✔️❌:               at[neigh].el_number == EL_NUMBER_S ) &&
    // INCHI✔️❌:              !at[neigh].num_H && 1 == at[neigh].valence &&
    // INCHI✔️❌:              2 == at[neigh].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* do not count double bonds to terminal =O or =S */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         num_DB += at[i].bond_type[j] == BOND_TYPE_DOUBLE;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_DB && !num_arom;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_C_unsat_not_arom

    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = i32::from(atom.valence);
    let bond_valence = i32::from(atom.chem_bonds_valence);
    if atom.el_number != EL_NUMBER_C
        || valence == bond_valence
        || valence + 1 < bond_valence
        || bond_valence + i32::from(atom.num_H) != 4
        || atom.charge != 0
        || atom.radical != 0
    {
        return Ok(0);
    }
    let mut aromatic = 0_i32;
    let mut double_bonds = 0_i32;
    for order in 0..usize::try_from(valence).map_err(|_| SourceHeapError::PointerOutOfBounds)? {
        let neighbor = atoms
            .get(usize::from(atom.neighbor[order]))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        aromatic += i32::from(atom.bond_type[order] == BOND_TYPE_ALTERN as u8);
        if matches!(neighbor.el_number, EL_NUMBER_O | EL_NUMBER_S)
            && neighbor.num_H == 0
            && neighbor.valence == 1
            && neighbor.chem_bonds_valence == 2
        {
            continue;
        }
        double_bonds += i32::from(atom.bond_type[order] == BOND_DOUBLE as u8);
    }
    Ok(i32::from(double_bonds != 0 && aromatic == 0))
}

#[allow(non_snake_case)]
pub(crate) fn is_Phenyl(
    atoms: &[inp_ATOM],
    outside_point: i32,
    attachment_point: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1559 is_Phenyl
    // INCHI✔️❌: int is_Phenyl( inp_ATOM *at, int outside_point, int attachment_point )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int iNext, iCur, iNewNext, k;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[attachment_point].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:          at[attachment_point].valence == 3 && /*at[attachment_point].chem_bonds_valence == 4 &&*/
    // INCHI✔️❌:          !at[attachment_point].num_H && !at[attachment_point].charge && !at[attachment_point].radical &&
    // INCHI✔️❌:          at[attachment_point].nRingSystem != at[outside_point].nRingSystem &&
    // INCHI✔️❌:          at[attachment_point].bCutVertex && at[attachment_point].nNumAtInRingSystem == 6)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         for (iNext = 0; iNext < at[attachment_point].valence; iNext++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[attachment_point].neighbor[iNext] != outside_point)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (iNext == at[attachment_point].valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /* program error*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:         iCur = attachment_point;
    // INCHI✔️❌:         iNext = at[attachment_point].neighbor[iNext];
    // INCHI✔️❌:         for (k = 0; k < 5; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* here we do not check bond type in the aromatic ring */
    // INCHI✔️❌:             if (at[iNext].el_number != EL_NUMBER_C || at[iNext].valence != 2 || at[iNext].num_H != 1 || at[iNext].charge || at[iNext].radical)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             iNewNext = at[iNext].neighbor[at[iNext].neighbor[0] == iCur];
    // INCHI✔️❌:             iCur = iNext;
    // INCHI✔️❌:             iNext = iNewNext;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return ( iNext == attachment_point );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_Phenyl
    // BEGIN ACTIVE INCHI HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h:48
    // INCHI✔️❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // END ACTIVE INCHI HEADER MACRO: EL_NUMBER_C

    let outside_point =
        usize::try_from(outside_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let attachment_point =
        usize::try_from(attachment_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let outside = atoms
        .get(outside_point)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let attachment = atoms
        .get(attachment_point)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    if attachment.el_number != EL_NUMBER_C
        || attachment.valence != 3
        || attachment.num_H != 0
        || attachment.charge != 0
        || attachment.radical != 0
        || attachment.nRingSystem == outside.nRingSystem
        || attachment.bCutVertex == 0
        || attachment.nNumAtInRingSystem != 6
    {
        return Ok(0);
    }

    let valence = usize::try_from(i32::from(attachment.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let Some(first_order) =
        (0..valence).find(|&order| usize::from(attachment.neighbor[order]) != outside_point)
    else {
        return Ok(0);
    };
    let mut current = attachment_point;
    let mut next = usize::from(attachment.neighbor[first_order]);
    for _ in 0..5 {
        let atom = atoms.get(next).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.el_number != EL_NUMBER_C
            || atom.valence != 2
            || atom.num_H != 1
            || atom.charge != 0
            || atom.radical != 0
        {
            return Ok(0);
        }
        let new_next =
            usize::from(atom.neighbor[usize::from(usize::from(atom.neighbor[0]) == current)]);
        current = next;
        next = new_next;
    }
    Ok(i32::from(next == attachment_point))
}

#[allow(non_snake_case)]
pub(crate) fn is_PentaFluoroPhenyl(
    atoms: &[inp_ATOM],
    outside_point: i32,
    attachment_point: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1615 is_PentaFluoroPhenyl
    // INCHI✔️❌: int is_PentaFluoroPhenyl( inp_ATOM *at,
    // INCHI✔️❌:                           int outside_point,
    // INCHI✔️❌:                           int attachment_point )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int iNext, iCur, iNewNext, nF, k, i, neigh;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[attachment_point].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:          at[attachment_point].valence == 3 && /*at[attachment_point].chem_bonds_valence == 4 &&*/
    // INCHI✔️❌:          !at[attachment_point].num_H && !at[attachment_point].charge && !at[attachment_point].radical &&
    // INCHI✔️❌:          at[attachment_point].nRingSystem != at[outside_point].nRingSystem &&
    // INCHI✔️❌:          at[attachment_point].bCutVertex && at[attachment_point].nNumAtInRingSystem == 6)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         for (iNext = 0; iNext < at[attachment_point].valence; iNext++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[attachment_point].neighbor[iNext] != outside_point)
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (iNext == at[attachment_point].valence)
    // INCHI✔️❌:             return 0; /* program error*/
    // INCHI✔️❌:         iCur = attachment_point;
    // INCHI✔️❌:         iNext = at[attachment_point].neighbor[iNext];
    // INCHI✔️❌:         for (k = 0; k < 5; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* here we do not check bond type in the aromatic ring */
    // INCHI✔️❌:             if (at[iNext].el_number != EL_NUMBER_C || at[iNext].valence != 3 || at[iNext].num_H != 0 || at[iNext].charge || at[iNext].radical)
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             for (i = 0, nF = 0, iNewNext = -1; i < at[iNext].valence; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 neigh = at[iNext].neighbor[i];
    // INCHI✔️❌:                 if (neigh == iCur)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                     if (at[neigh].el_number == EL_NUMBER_F && at[neigh].chem_bonds_valence == 1 && !at[neigh].charge && !at[neigh].radical && !at[neigh].num_H)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nF++; /* terminal flourine */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                         if (iNewNext == -1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             iNewNext = neigh; /* Carbon will be checked on the next pass */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (iNewNext == -1 || nF != 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             iCur = iNext;
    // INCHI✔️❌:             iNext = iNewNext;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return ( iNext == attachment_point );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_PentaFluoroPhenyl
    // BEGIN ACTIVE INCHI HEADER MACROS: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h
    // INCHI✔️❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️❌: #define EL_NUMBER_F  ((U_CHAR)9)
    // END ACTIVE INCHI HEADER MACROS: EL_NUMBER_C, EL_NUMBER_F

    let outside_point =
        usize::try_from(outside_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let attachment_point =
        usize::try_from(attachment_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let outside = atoms
        .get(outside_point)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let attachment = atoms
        .get(attachment_point)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if attachment.el_number != EL_NUMBER_C
        || attachment.valence != 3
        || attachment.num_H != 0
        || attachment.charge != 0
        || attachment.radical != 0
        || attachment.nRingSystem == outside.nRingSystem
        || attachment.bCutVertex == 0
        || attachment.nNumAtInRingSystem != 6
    {
        return Ok(0);
    }

    let valence = usize::try_from(i32::from(attachment.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let Some(first_order) =
        (0..valence).find(|&order| usize::from(attachment.neighbor[order]) != outside_point)
    else {
        return Ok(0);
    };
    let mut current = attachment_point;
    let mut next = usize::from(attachment.neighbor[first_order]);
    for _ in 0..5 {
        let atom = atoms.get(next).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.el_number != EL_NUMBER_C
            || atom.valence != 3
            || atom.num_H != 0
            || atom.charge != 0
            || atom.radical != 0
        {
            return Ok(0);
        }
        let valence = usize::try_from(i32::from(atom.valence))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut fluorines = 0_i32;
        let mut new_next = None;
        for order in 0..valence {
            let neighbor_index = usize::from(atom.neighbor[order]);
            if neighbor_index == current {
                continue;
            }
            let neighbor = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.el_number == EL_NUMBER_F
                && neighbor.chem_bonds_valence == 1
                && neighbor.charge == 0
                && neighbor.radical == 0
                && neighbor.num_H == 0
            {
                fluorines += 1;
            } else if new_next.is_none() {
                new_next = Some(neighbor_index);
            } else {
                return Ok(0);
            }
        }
        let Some(candidate) = new_next else {
            return Ok(0);
        };
        if fluorines != 1 {
            return Ok(0);
        }
        current = next;
        next = candidate;
    }
    Ok(i32::from(next == attachment_point))
}

#[allow(non_snake_case)]
pub(crate) fn is_Methyl(atoms: &[inp_ATOM], attachment_point: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1679 is_Methyl
    // INCHI✔️❌: int is_Methyl( inp_ATOM *at, int attachment_point )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (at[attachment_point].valence == 1 && at[attachment_point].chem_bonds_valence == 1 &&
    // INCHI✔️❌:          at[attachment_point].el_number == EL_NUMBER_C && at[attachment_point].num_H == 3 &&
    // INCHI✔️❌:          !at[attachment_point].charge && !at[attachment_point].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* methyl */
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_Methyl
    // BEGIN ACTIVE INCHI HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h:48
    // INCHI✔️❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // END ACTIVE INCHI HEADER MACRO: EL_NUMBER_C

    let attachment_point =
        usize::try_from(attachment_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(attachment_point)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(
        atom.valence == 1
            && atom.chem_bonds_valence == 1
            && atom.el_number == EL_NUMBER_C
            && atom.num_H == 3
            && atom.charge == 0
            && atom.radical == 0,
    ))
}

#[allow(non_snake_case)]
pub(crate) fn is_Ethyl(
    atoms: &[inp_ATOM],
    outside_point: i32,
    attachment_point: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1694 is_Ethyl
    // INCHI✔️❌: int is_Ethyl( inp_ATOM *at, int outside_point, int attachment_point )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (at[attachment_point].valence == 2 && at[attachment_point].chem_bonds_valence == 2 &&
    // INCHI✔️❌:          at[attachment_point].el_number == EL_NUMBER_C && at[attachment_point].num_H == 2 &&
    // INCHI✔️❌:          !at[attachment_point].charge && !at[attachment_point].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* methanediyl */
    // INCHI✔️❌:         int iat_methyl = at[attachment_point].neighbor[( at[attachment_point].neighbor[0] == outside_point )];
    // INCHI✔️❌:         return is_Methyl( at, iat_methyl );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_Ethyl
    // BEGIN ACTIVE INCHI HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h:48
    // INCHI✔️❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // END ACTIVE INCHI HEADER MACRO: EL_NUMBER_C

    let outside_point =
        usize::try_from(outside_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let attachment_point =
        usize::try_from(attachment_point).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let attachment = atoms
        .get(attachment_point)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if attachment.valence != 2
        || attachment.chem_bonds_valence != 2
        || attachment.el_number != EL_NUMBER_C
        || attachment.num_H != 2
        || attachment.charge != 0
        || attachment.radical != 0
    {
        return Ok(0);
    }
    let methyl_order = usize::from(attachment.neighbor[0] as usize == outside_point);
    let methyl = i32::from(attachment.neighbor[methyl_order]);
    is_Methyl(atoms, methyl)
}

#[allow(non_snake_case)]
pub(crate) fn is_Methyl_or_Etyl(
    atoms: &[inp_ATOM],
    outside_point: i32,
    attachment_point: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1709 is_Methyl_or_Etyl
    // INCHI✔️❌: int is_Methyl_or_Etyl( inp_ATOM *at, int outside_point, int attachment_point )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (is_Methyl( at, attachment_point ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (is_Ethyl( at, outside_point, attachment_point ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 2;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_Methyl_or_Etyl

    if is_Methyl(atoms, attachment_point)? != 0 {
        return Ok(1);
    }
    if is_Ethyl(atoms, outside_point, attachment_point)? != 0 {
        return Ok(2);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn is_Si_IV(atoms: &[inp_ATOM], atom_index: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1727 is_Si_IV
    // INCHI✔️❌: int is_Si_IV( inp_ATOM *at, int i )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (at[i].el_number != EL_NUMBER_SI ||
    // INCHI✔️❌:          at[i].charge || at[i].radical || at[i].valence != 4 || at[i].chem_bonds_valence != 4)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_Si_IV
    // BEGIN ACTIVE INCHI HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h:52
    // INCHI✔️❌: #define EL_NUMBER_SI ((U_CHAR)14)
    // END ACTIVE INCHI HEADER MACRO: EL_NUMBER_SI

    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(
        atom.el_number == EL_NUMBER_SI
            && atom.charge == 0
            && atom.radical == 0
            && atom.valence == 4
            && atom.chem_bonds_valence == 4,
    ))
}

#[allow(non_snake_case)]
pub(crate) fn is_DERIV_RING_DMOX_DEOX_O(
    atoms: &[inp_ATOM],
    cur_atom: i32,
    from_ord: i32,
    _da: Option<&mut DERIV_AT>,
    da1: Option<&mut DERIV_AT>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1802 is_DERIV_RING_DMOX_DEOX_O
    // INCHI✔️❌: int is_DERIV_RING_DMOX_DEOX_O( inp_ATOM *at,
    // INCHI✔️❌:                                int cur_atom,
    // INCHI✔️❌:                                int from_ord,
    // INCHI✔️❌:                                DERIV_AT *da,
    // INCHI✔️❌:                                DERIV_AT *da1 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*-----------------------
    // INCHI✔️❌:     <-
    // INCHI✔️❌:     #4 #3
    // INCHI✔️❌:     R--C==N   Me or Et
    // INCHI✔️❌:     |   \ /
    // INCHI✔️❌:     |    C #2
    // INCHI✔️❌:     |   / \
    // INCHI✔️❌:     at[k]:O--CH2 Me or ET
    // INCHI✔️❌:     #0  #1
    // INCHI✔️❌:     ->
    // INCHI✔️❌:     --------------------------*/
    // INCHI✔️❌:     /*            #0           #1           #2           #3           #4 */
    // INCHI✔️❌:     static const U_CHAR bond_type[OX_RING_SIZE] = { BOND_SINGLE, BOND_SINGLE, BOND_SINGLE, BOND_DOUBLE, BOND_SINGLE };
    // INCHI✔️❌:     static const S_CHAR valence[OX_RING_SIZE] = { 2,           2,           4,           2,           3 };
    // INCHI✔️❌:     static const S_CHAR bonds_valence[OX_RING_SIZE] = { 2,           2,           4,           3,           4 };
    // INCHI✔️❌:     static const S_CHAR num_H[OX_RING_SIZE] = { 0,           2,           0,           0,           0 };
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_NUMB from, curr, next, nRingSystem, at_no[OX_RING_SIZE];
    // INCHI✔️❌:     S_CHAR  bond_no[OX_RING_SIZE];
    // INCHI✔️❌:     int     i, n0, n1, attach1, attach2, neigh;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[cur_atom].el_number == EL_NUMBER_O && at[cur_atom].nNumAtInRingSystem == OX_RING_SIZE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AT_NUMB attype[OX_RING_SIZE] = { (AT_NUMB) EL_NUMBER_O, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_N, (AT_NUMB) EL_NUMBER_C };
    // INCHI✔️❌:
    // INCHI✔️❌:         curr = cur_atom;
    // INCHI✔️❌:         from = at[curr].neighbor[from_ord];
    // INCHI✔️❌:         nRingSystem = at[curr].nRingSystem;
    // INCHI✔️❌:         n0 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         do
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* find next atom in a simple ring */
    // INCHI✔️❌:             for (i = 0; i < at[curr].valence &&
    // INCHI✔️❌:                 ( from == ( next = at[curr].neighbor[i] ) ||
    // INCHI✔️❌:                   nRingSystem != at[next].nRingSystem ); i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i == at[curr].valence)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* check curr atom */
    // INCHI✔️❌:             if (at[curr].charge || at[curr].radical)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (at[curr].bond_type[i] != bond_type[n0] ||
    // INCHI✔️❌:                  at[curr].valence != valence[n0] ||
    // INCHI✔️❌:                  at[curr].chem_bonds_valence != bonds_valence[n0] ||
    // INCHI✔️❌:                  at[curr].num_H != num_H[n0] ||
    // INCHI✔️❌:                  at[curr].el_number != (U_CHAR) attype[n0])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* save current atom */
    // INCHI✔️❌:             at_no[n0] = curr;
    // INCHI✔️❌:             bond_no[n0] = i;
    // INCHI✔️❌:             /* prepare for the next */
    // INCHI✔️❌:             from = curr;
    // INCHI✔️❌:             curr = next;
    // INCHI✔️❌:             n0++;
    // INCHI✔️❌:         } while (n0 < OX_RING_SIZE && curr != cur_atom);
    // INCHI✔️❌:         /* check completion */
    // INCHI✔️❌:         if (OX_RING_SIZE != n0 || curr != cur_atom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* check if R is C */
    // INCHI✔️❌:         n1 = at_no[4];
    // INCHI✔️❌:         for (i = 0; i < at[n1].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[n1].neighbor[i];
    // INCHI✔️❌:             if (neigh != at_no[0] && neigh != at_no[3])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].el_number != EL_NUMBER_C)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto check_next_derivative2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     break; /* checked */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* check >C< attachments */
    // INCHI✔️❌:         n1 = at_no[2];
    // INCHI✔️❌:         attach1 = attach2 = 0;
    // INCHI✔️❌:         for (i = 0; i < at[n1].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[n1].neighbor[i] != at_no[1] &&
    // INCHI✔️❌:                  at[n1].neighbor[i] != at_no[3])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!attach1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     attach1 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!attach2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         attach2 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative2;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!attach2 || attach2 != attach1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* all checks are done */
    // INCHI✔️❌:         if ( /*da &&*/ da1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             short ord_O = bond_no[0];
    // INCHI✔️❌:             /*short ord_N = !bond_no[3];*/
    // INCHI✔️❌:             AT_NUMB iN = at_no[3];
    // INCHI✔️❌:             /*AT_NUMB iO  = at_no[0];*/
    // INCHI✔️❌:             char num_2remove = 2 + attach1 + attach2;
    // INCHI✔️❌:             da1->typ[0]    /* = da1->typ[1] */ = DERIV_RING_DMOX_DEOX_O;
    // INCHI✔️❌:             da1->ord[0]    /* = da1->ord[1] */ = ord_O;
    // INCHI✔️❌:             da1->num[0]    /* = da1->num[1] */ = num_2remove;
    // INCHI✔️❌:             da1->other_atom = iN + 1;
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             if ( da1->typ[0] ) {
    // INCHI✔️❌:             if ( da1->typ[0]     != DERIV_RING_DMOX_DEOX_O ||
    // INCHI✔️❌:             da1->ord[0]     != ord_O ||
    // INCHI✔️❌:             da1->num[0]     != num_2remove ||
    // INCHI✔️❌:             da1->other_atom != iN ) {
    // INCHI✔️❌:             goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             } else {
    // INCHI✔️❌:             da1->typ[0]     = DERIV_RING_DMOX_DEOX_O;
    // INCHI✔️❌:             da1->ord[0]     = ord_O;
    // INCHI✔️❌:             da1->num[0]     = num_2remove;
    // INCHI✔️❌:             da1->other_atom = iN;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return DERIV_RING_DMOX_DEOX_O;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: check_next_derivative2:;
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_DERIV_RING_DMOX_DEOX_O

    const OX_RING_SIZE: usize = 5;
    const BOND_TYPES: [u8; OX_RING_SIZE] = [
        BOND_SINGLE as u8,
        BOND_SINGLE as u8,
        BOND_SINGLE as u8,
        BOND_DOUBLE as u8,
        BOND_SINGLE as u8,
    ];
    const VALENCES: [i8; OX_RING_SIZE] = [2, 2, 4, 2, 3];
    const BOND_VALENCES: [i8; OX_RING_SIZE] = [2, 2, 4, 3, 4];
    const HYDROGENS: [i8; OX_RING_SIZE] = [0, 2, 0, 0, 0];
    const ATOM_TYPES: [u8; OX_RING_SIZE] = [
        EL_NUMBER_O,
        EL_NUMBER_C,
        EL_NUMBER_C,
        EL_NUMBER_N,
        EL_NUMBER_C,
    ];

    let cur_atom = usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let from_ord = usize::try_from(from_ord).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let start = atoms
        .get(cur_atom)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if start.el_number != EL_NUMBER_O || usize::from(start.nNumAtInRingSystem) != OX_RING_SIZE {
        return Ok(0);
    }

    let mut current = cur_atom;
    let mut from = usize::from(
        *start
            .neighbor
            .get(from_ord)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let ring_system = start.nRingSystem;
    let mut atom_numbers = [0_usize; OX_RING_SIZE];
    let mut bond_numbers = [0_i8; OX_RING_SIZE];
    let mut count = 0_usize;
    loop {
        let atom = atoms
            .get(current)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = usize::try_from(i32::from(atom.valence))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut selected = None;
        for order in 0..valence {
            let next = usize::from(atom.neighbor[order]);
            let next_atom = atoms.get(next).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if from != next && ring_system == next_atom.nRingSystem {
                selected = Some((order, next));
                break;
            }
        }
        let Some((order, next)) = selected else {
            return Ok(0);
        };
        if atom.charge != 0 || atom.radical != 0 {
            return Ok(0);
        }
        if atom.bond_type[order] != BOND_TYPES[count]
            || atom.valence != VALENCES[count]
            || atom.chem_bonds_valence != BOND_VALENCES[count]
            || atom.num_H != HYDROGENS[count]
            || atom.el_number != ATOM_TYPES[count]
        {
            return Ok(0);
        }
        atom_numbers[count] = current;
        bond_numbers[count] = i8::try_from(order).unwrap();
        from = current;
        current = next;
        count += 1;
        if count >= OX_RING_SIZE || current == cur_atom {
            break;
        }
    }
    if count != OX_RING_SIZE || current != cur_atom {
        return Ok(0);
    }

    let carbon_r = atom_numbers[4];
    let atom = atoms
        .get(carbon_r)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = usize::try_from(i32::from(atom.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for order in 0..valence {
        let neighbor = usize::from(atom.neighbor[order]);
        if neighbor != atom_numbers[0] && neighbor != atom_numbers[3] {
            if atoms
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .el_number
                != EL_NUMBER_C
            {
                return Ok(0);
            }
            break;
        }
    }

    let central = atom_numbers[2];
    let atom = atoms
        .get(central)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = usize::try_from(i32::from(atom.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut attach1 = 0_i32;
    let mut attach2 = 0_i32;
    for order in 0..valence {
        let neighbor = usize::from(atom.neighbor[order]);
        if neighbor != atom_numbers[1] && neighbor != atom_numbers[3] {
            if attach1 == 0 {
                attach1 = is_Methyl_or_Etyl(atoms, central as i32, neighbor as i32)?;
            } else if attach2 == 0 {
                attach2 = is_Methyl_or_Etyl(atoms, central as i32, neighbor as i32)?;
            } else {
                return Ok(0);
            }
        }
    }
    if attach2 == 0 || attach2 != attach1 {
        return Ok(0);
    }

    if let Some(output) = da1 {
        output.typ[0] = DERIV_RING_DMOX_DEOX_O as i16;
        output.ord[0] = bond_numbers[0];
        output.num[0] = i8::try_from(2 + attach1 + attach2).unwrap();
        output.other_atom = (atom_numbers[3] as u16).wrapping_add(1);
    }
    Ok(DERIV_RING_DMOX_DEOX_O as i32)
}

#[allow(non_snake_case)]
pub(crate) fn is_DERIV_RING_DMOX_DEOX_N(
    atoms: &[inp_ATOM],
    cur_atom: i32,
    from_ord: i32,
    _da: Option<&mut DERIV_AT>,
    da1: Option<&mut DERIV_AT>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1962 is_DERIV_RING_DMOX_DEOX_N
    // INCHI✔️❌: int is_DERIV_RING_DMOX_DEOX_N( inp_ATOM *at,
    // INCHI✔️❌:                                int cur_atom,
    // INCHI✔️❌:                                int from_ord,
    // INCHI✔️❌:                                DERIV_AT *da,
    // INCHI✔️❌:                                DERIV_AT *da1 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*-----------------------
    // INCHI✔️❌:     ->
    // INCHI✔️❌:     #4 #0
    // INCHI✔️❌:     R--C==N   Me or Et
    // INCHI✔️❌:     |   \ /
    // INCHI✔️❌:     |    C #1
    // INCHI✔️❌:     |   / \
    // INCHI✔️❌:     at[k]:O--CH2 Me or ET
    // INCHI✔️❌:     #3  #2
    // INCHI✔️❌:     <-
    // INCHI✔️❌:     --------------------------*/
    // INCHI✔️❌:     /*            #0           #1           #2           #3           #4 */
    // INCHI✔️❌:     static const U_CHAR bond_type[OX_RING_SIZE] = { BOND_SINGLE, BOND_SINGLE, BOND_SINGLE, BOND_SINGLE, BOND_DOUBLE };
    // INCHI✔️❌:     static const S_CHAR valence[OX_RING_SIZE] = { 2,           4,           2,           2,           3 };
    // INCHI✔️❌:     static const S_CHAR bonds_valence[OX_RING_SIZE] = { 3,           4,           2,           2,           4 };
    // INCHI✔️❌:     static const S_CHAR num_H[OX_RING_SIZE] = { 0,           0,           2,           0,           0 };
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_NUMB from, curr, next, nRingSystem, at_no[OX_RING_SIZE];
    // INCHI✔️❌:     S_CHAR  bond_no[OX_RING_SIZE];
    // INCHI✔️❌:     int     i, n0, n1, attach1, attach2, neigh;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[cur_atom].el_number == EL_NUMBER_N && at[cur_atom].nNumAtInRingSystem == OX_RING_SIZE &&
    // INCHI✔️❌:          at[cur_atom].valence == 2 && at[cur_atom].chem_bonds_valence == 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AT_NUMB attype[OX_RING_SIZE] = { (AT_NUMB) EL_NUMBER_N, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_O, (AT_NUMB) EL_NUMBER_C };
    // INCHI✔️❌:
    // INCHI✔️❌:         curr = cur_atom;
    // INCHI✔️❌:         from = at[curr].neighbor[from_ord];
    // INCHI✔️❌:         nRingSystem = at[curr].nRingSystem;
    // INCHI✔️❌:         n0 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         do
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* find next atom in a simple ring */
    // INCHI✔️❌:             for (i = 0; i < at[curr].valence &&
    // INCHI✔️❌:                 ( from == ( next = at[curr].neighbor[i] ) ||
    // INCHI✔️❌:                   nRingSystem != at[next].nRingSystem ); i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i == at[curr].valence)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* check curr atom */
    // INCHI✔️❌:             if (at[curr].charge || at[curr].radical)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (at[curr].bond_type[i] != bond_type[n0] ||
    // INCHI✔️❌:                  at[curr].valence != valence[n0] ||
    // INCHI✔️❌:                  at[curr].chem_bonds_valence != bonds_valence[n0] ||
    // INCHI✔️❌:                  at[curr].num_H != num_H[n0] ||
    // INCHI✔️❌:                  at[curr].el_number != (U_CHAR) attype[n0])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* save current atom */
    // INCHI✔️❌:             at_no[n0] = curr;
    // INCHI✔️❌:             bond_no[n0] = i;
    // INCHI✔️❌:             /* prepare for the next */
    // INCHI✔️❌:             from = curr;
    // INCHI✔️❌:             curr = next;
    // INCHI✔️❌:             n0++;
    // INCHI✔️❌:         } while (n0 < OX_RING_SIZE && curr != cur_atom);
    // INCHI✔️❌:
    // INCHI✔️❌:         /* check completion */
    // INCHI✔️❌:         if (OX_RING_SIZE != n0 || curr != cur_atom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* check if R is C */
    // INCHI✔️❌:         n1 = at_no[4];
    // INCHI✔️❌:         for (i = 0; i < at[n1].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[n1].neighbor[i];
    // INCHI✔️❌:             if (neigh != at_no[0] && neigh != at_no[3])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].el_number != EL_NUMBER_C)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto check_next_derivative2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     break; /* checked */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* check >C< attachments */
    // INCHI✔️❌:         n1 = at_no[1];
    // INCHI✔️❌:         attach1 = attach2 = 0;
    // INCHI✔️❌:         for (i = 0; i < at[n1].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[n1].neighbor[i] != at_no[0] &&
    // INCHI✔️❌:                  at[n1].neighbor[i] != at_no[2])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!attach1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     attach1 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!attach2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         attach2 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative2;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!attach2 || attach2 != attach1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* all checks are done */
    // INCHI✔️❌:         if ( /*da &&*/ da1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*short ord_O = !bond_no[3];*/
    // INCHI✔️❌:             short ord_N = bond_no[0];
    // INCHI✔️❌:             /*AT_NUMB iN  = at_no[0];*/
    // INCHI✔️❌:             AT_NUMB iO = at_no[3];
    // INCHI✔️❌:             char num_2remove = 2 + attach1 + attach2;
    // INCHI✔️❌:             da1->typ[0]    /* = da1->typ[1] */ = DERIV_RING_DMOX_DEOX_N;
    // INCHI✔️❌:             da1->ord[0]    /* = da1->ord[1] */ = ord_N;
    // INCHI✔️❌:             da1->num[0]    /* = da1->num[1] */ = num_2remove;
    // INCHI✔️❌:             da1->other_atom = iO + 1;
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             if ( da1->typ[0] ) {
    // INCHI✔️❌:             if ( da1->typ[0]     != DERIV_RING_DMOX_DEOX_O ||
    // INCHI✔️❌:             da1->ord[0]     != ord_O ||
    // INCHI✔️❌:             da1->num[0]     != num_2remove ||
    // INCHI✔️❌:             da1->other_atom != iN ) {
    // INCHI✔️❌:             goto check_next_derivative2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             } else {
    // INCHI✔️❌:             da1->typ[0]     = DERIV_RING_DMOX_DEOX_O;
    // INCHI✔️❌:             da1->ord[0]     = ord_O;
    // INCHI✔️❌:             da1->num[0]     = num_2remove;
    // INCHI✔️❌:             da1->other_atom = iN;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return DERIV_RING_DMOX_DEOX_N;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: check_next_derivative2:;
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_DERIV_RING_DMOX_DEOX_N

    const OX_RING_SIZE: usize = 5;
    const BOND_TYPES: [u8; OX_RING_SIZE] = [
        BOND_SINGLE as u8,
        BOND_SINGLE as u8,
        BOND_SINGLE as u8,
        BOND_SINGLE as u8,
        BOND_DOUBLE as u8,
    ];
    const VALENCES: [i8; OX_RING_SIZE] = [2, 4, 2, 2, 3];
    const BOND_VALENCES: [i8; OX_RING_SIZE] = [3, 4, 2, 2, 4];
    const HYDROGENS: [i8; OX_RING_SIZE] = [0, 0, 2, 0, 0];
    const ATOM_TYPES: [u8; OX_RING_SIZE] = [
        EL_NUMBER_N,
        EL_NUMBER_C,
        EL_NUMBER_C,
        EL_NUMBER_O,
        EL_NUMBER_C,
    ];

    let cur_atom = usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let from_ord = usize::try_from(from_ord).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let start = atoms
        .get(cur_atom)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if start.el_number != EL_NUMBER_N
        || usize::from(start.nNumAtInRingSystem) != OX_RING_SIZE
        || start.valence != 2
        || start.chem_bonds_valence != 3
    {
        return Ok(0);
    }

    let mut current = cur_atom;
    let mut from = usize::from(
        *start
            .neighbor
            .get(from_ord)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let ring_system = start.nRingSystem;
    let mut atom_numbers = [0_usize; OX_RING_SIZE];
    let mut bond_numbers = [0_i8; OX_RING_SIZE];
    let mut count = 0_usize;
    loop {
        let atom = atoms
            .get(current)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = usize::try_from(i32::from(atom.valence))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut selected = None;
        for order in 0..valence {
            let next = usize::from(atom.neighbor[order]);
            let next_atom = atoms.get(next).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if from != next && ring_system == next_atom.nRingSystem {
                selected = Some((order, next));
                break;
            }
        }
        let Some((order, next)) = selected else {
            return Ok(0);
        };
        if atom.charge != 0 || atom.radical != 0 {
            return Ok(0);
        }
        if atom.bond_type[order] != BOND_TYPES[count]
            || atom.valence != VALENCES[count]
            || atom.chem_bonds_valence != BOND_VALENCES[count]
            || atom.num_H != HYDROGENS[count]
            || atom.el_number != ATOM_TYPES[count]
        {
            return Ok(0);
        }
        atom_numbers[count] = current;
        bond_numbers[count] = i8::try_from(order).unwrap();
        from = current;
        current = next;
        count += 1;
        if count >= OX_RING_SIZE || current == cur_atom {
            break;
        }
    }
    if count != OX_RING_SIZE || current != cur_atom {
        return Ok(0);
    }

    let carbon_r = atom_numbers[4];
    let atom = atoms
        .get(carbon_r)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = usize::try_from(i32::from(atom.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for order in 0..valence {
        let neighbor = usize::from(atom.neighbor[order]);
        if neighbor != atom_numbers[0] && neighbor != atom_numbers[3] {
            if atoms
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .el_number
                != EL_NUMBER_C
            {
                return Ok(0);
            }
            break;
        }
    }

    let central = atom_numbers[1];
    let atom = atoms
        .get(central)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = usize::try_from(i32::from(atom.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut attach1 = 0_i32;
    let mut attach2 = 0_i32;
    for order in 0..valence {
        let neighbor = usize::from(atom.neighbor[order]);
        if neighbor != atom_numbers[0] && neighbor != atom_numbers[2] {
            if attach1 == 0 {
                attach1 = is_Methyl_or_Etyl(atoms, central as i32, neighbor as i32)?;
            } else if attach2 == 0 {
                attach2 = is_Methyl_or_Etyl(atoms, central as i32, neighbor as i32)?;
            } else {
                return Ok(0);
            }
        }
    }
    if attach2 == 0 || attach2 != attach1 {
        return Ok(0);
    }

    if let Some(output) = da1 {
        output.typ[0] = DERIV_RING_DMOX_DEOX_N as i16;
        output.ord[0] = bond_numbers[0];
        output.num[0] = i8::try_from(2 + attach1 + attach2).unwrap();
        output.other_atom = (atom_numbers[3] as u16).wrapping_add(1);
    }
    Ok(DERIV_RING_DMOX_DEOX_N as i32)
}

#[allow(non_snake_case)]
pub(crate) fn is_DERIV_RING2_PRRLDD_PPRDN(
    atoms: &[inp_ATOM],
    cur_atom: i32,
    from_ord: i32,
    _da: Option<&mut DERIV_AT>,
    da1: Option<&mut DERIV_AT>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2147 is_DERIV_RING2_PRRLDD_PPRDN
    // INCHI✔️❌: int is_DERIV_RING2_PRRLDD_PPRDN( inp_ATOM *at,
    // INCHI✔️❌:                                  int cur_atom,
    // INCHI✔️❌:                                  int from_ord,
    // INCHI✔️❌:                                  DERIV_AT *da,
    // INCHI✔️❌:                                  DERIV_AT *da1 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     #1   #2                                      #1   #2
    // INCHI✔️❌:     O   CH2--CH2                                 O   CH2--CH2            O
    // INCHI✔️❌:     ||  /     |                                  ||  /     |             ||
    // INCHI✔️❌:     R--C--N #0   |     N is cur_atom             R--C--N #0  CH2 #3  =>  R--C--OH
    // INCHI✔️❌:     from  \     |     C(IV) is from atom         from  \     |
    // INCHI✔️❌:     CH2--CH2                                     CH2--CH2
    // INCHI✔️❌:     #4   #3                                      #5   #4
    // INCHI✔️❌:
    // INCHI✔️❌:     Pyrrolidides: Replace -N< with -OH           Piperidines: Replace -N< with -OH
    // INCHI✔️❌:     DERIV_RING2_PRRLDD_OUTSIDE_PRECUR            DERIV_RING2_PPRDN_OUTSIDE_PRECUR
    // INCHI✔️❌:     */
    // INCHI✔️❌:     int iat_from, i, neigh, k;
    // INCHI✔️❌:     char ord[2];
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check cur atom */
    // INCHI✔️❌:     if (
    // INCHI✔️❌:         at[cur_atom].el_number == EL_NUMBER_N &&
    // INCHI✔️❌:         MIN_PRRLDD_PPRDN_RING_SIZE <= at[cur_atom].nNumAtInRingSystem &&
    // INCHI✔️❌:         at[cur_atom].nNumAtInRingSystem <= MAX_PRRLDD_PPRDN_RING_SIZE &&
    // INCHI✔️❌:         at[cur_atom].valence == 3 && at[cur_atom].chem_bonds_valence == 3 &&
    // INCHI✔️❌:         !at[cur_atom].charge && !at[cur_atom].radical && !at[cur_atom].num_H &&
    // INCHI✔️❌:         /* check the "from" atom (on the left from cur atom) */
    // INCHI✔️❌:         at[iat_from = at[cur_atom].neighbor[from_ord]].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:         at[iat_from].nNumAtInRingSystem == 1 &&
    // INCHI✔️❌:         at[iat_from].valence == 3 &&
    // INCHI✔️❌:         at[iat_from].chem_bonds_valence == 4 &&
    // INCHI✔️❌:         !at[iat_from].charge && !at[iat_from].radical && !at[iat_from].num_H)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* check neighbors of the "from" atom (on the left from cur atom) */
    // INCHI✔️❌:         for (i = 0; i < at[iat_from].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[iat_from].neighbor[i];
    // INCHI✔️❌:             if (neigh == cur_atom)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[iat_from].bond_type[i] == BOND_SINGLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[neigh].el_number != EL_NUMBER_C)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[iat_from].bond_type[i] == BOND_DOUBLE)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (at[neigh].el_number != EL_NUMBER_O ||
    // INCHI✔️❌:                              at[neigh].valence != 1 || at[neigh].chem_bonds_valence != 2 ||
    // INCHI✔️❌:                              at[neigh].charge || at[neigh].radical || at[neigh].num_H)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto check_next_derivative;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* check the ring */
    // INCHI✔️❌:         iat_from = cur_atom;
    // INCHI✔️❌:         neigh = at[cur_atom].neighbor[!from_ord]; /* any except "from" atom */
    // INCHI✔️❌:         i = 1;
    // INCHI✔️❌:         if (at[cur_atom].nRingSystem != at[neigh].nRingSystem)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         do
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[neigh].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:                  at[neigh].valence != 2 ||
    // INCHI✔️❌:                  at[neigh].chem_bonds_valence != 2 ||
    // INCHI✔️❌:                  at[neigh].num_H != 2 ||
    // INCHI✔️❌:                  at[neigh].charge || at[neigh].radical)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             i++; /* at[neigh] satisfied the conditions */
    // INCHI✔️❌:             k = at[neigh].neighbor[at[neigh].neighbor[0] == iat_from];
    // INCHI✔️❌:             iat_from = neigh;
    // INCHI✔️❌:             neigh = k;
    // INCHI✔️❌:         } while (neigh != cur_atom && i <= MAX_PRRLDD_PPRDN_RING_SIZE);
    // INCHI✔️❌:
    // INCHI✔️❌:         if (neigh == cur_atom &&
    // INCHI✔️❌:              MIN_PRRLDD_PPRDN_RING_SIZE <= i && i <= MAX_PRRLDD_PPRDN_RING_SIZE)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da1)
    // INCHI✔️❌:             {
    // INCHI✔️❌: #ifdef DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
    // INCHI✔️❌:                 if (i == PRRLDD_RING_SIZE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da1->typ[0] = da1->typ[1] = DERIV_RING2_PRRLDD_OUTSIDE_PRECUR;
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_RING2_PPRDN_OUTSIDE_PRECUR
    // INCHI✔️❌:                 if (i == PPRDN_RING_SIZE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da1->typ[0] = da1->typ[1] = DERIV_RING2_PPRDN_OUTSIDE_PRECUR;
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 ord[0] = !from_ord;
    // INCHI✔️❌:                 ord[1] = ( ( !from_ord + 1 ) ^ ( from_ord + 1 ) ) - 1; /* the 3rd possible index out of (0,1,2) */
    // INCHI✔️❌:                 k = ( ord[1] < ord[0] );
    // INCHI✔️❌:                 da1->ord[0] = ord[k];  /* smaller */
    // INCHI✔️❌:                 da1->ord[1] = ord[!k]; /* greater */
    // INCHI✔️❌:                 /*da1->num[0] = */da1->num[0] = i - 1; /* djb-rwth: unresolved issue -- revision required? / da1->num[1] = i-1? */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return i;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: check_next_derivative:
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_DERIV_RING2_PRRLDD_PPRDN

    const MIN_RING_SIZE: i32 = 5;
    const MAX_RING_SIZE: i32 = 6;
    let cur_atom = usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let from_ord_usize =
        usize::try_from(from_ord).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = atoms
        .get(cur_atom)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let from_atom = usize::from(
        *current
            .neighbor
            .get(from_ord_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let from = atoms
        .get(from_atom)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if current.el_number != EL_NUMBER_N
        || i32::from(current.nNumAtInRingSystem) < MIN_RING_SIZE
        || i32::from(current.nNumAtInRingSystem) > MAX_RING_SIZE
        || current.valence != 3
        || current.chem_bonds_valence != 3
        || current.charge != 0
        || current.radical != 0
        || current.num_H != 0
        || from.el_number != EL_NUMBER_C
        || from.nNumAtInRingSystem != 1
        || from.valence != 3
        || from.chem_bonds_valence != 4
        || from.charge != 0
        || from.radical != 0
        || from.num_H != 0
    {
        return Ok(0);
    }

    let from_valence = usize::try_from(i32::from(from.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for order in 0..from_valence {
        let neighbor_index = usize::from(from.neighbor[order]);
        if neighbor_index == cur_atom {
            continue;
        }
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if from.bond_type[order] == BOND_SINGLE as u8 {
            if neighbor.el_number != EL_NUMBER_C {
                return Ok(0);
            }
        } else if from.bond_type[order] == BOND_DOUBLE as u8 {
            if neighbor.el_number != EL_NUMBER_O
                || neighbor.valence != 1
                || neighbor.chem_bonds_valence != 2
                || neighbor.charge != 0
                || neighbor.radical != 0
                || neighbor.num_H != 0
            {
                return Ok(0);
            }
        } else {
            return Ok(0);
        }
    }

    let not_from_order = usize::from(from_ord == 0);
    let mut previous = cur_atom;
    let mut neighbor_index = usize::from(current.neighbor[not_from_order]);
    let first_ring_atom = atoms
        .get(neighbor_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut ring_size = 1_i32;
    if current.nRingSystem != first_ring_atom.nRingSystem {
        return Ok(0);
    }
    loop {
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.el_number != EL_NUMBER_C
            || neighbor.valence != 2
            || neighbor.chem_bonds_valence != 2
            || neighbor.num_H != 2
            || neighbor.charge != 0
            || neighbor.radical != 0
        {
            return Ok(0);
        }
        ring_size += 1;
        let next_order = usize::from(usize::from(neighbor.neighbor[0]) == previous);
        let next = usize::from(neighbor.neighbor[next_order]);
        previous = neighbor_index;
        neighbor_index = next;
        if neighbor_index == cur_atom || ring_size > MAX_RING_SIZE {
            break;
        }
    }

    if neighbor_index != cur_atom || !(MIN_RING_SIZE..=MAX_RING_SIZE).contains(&ring_size) {
        return Ok(0);
    }
    if let Some(output) = da1 {
        let derivative_type = if ring_size == 5 {
            DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
        } else {
            DERIV_RING2_PPRDN_OUTSIDE_PRECUR
        };
        output.typ[1] = derivative_type as i16;
        output.typ[0] = output.typ[1];
        let ord0 = i8::from(from_ord == 0);
        let ord1 = (((i32::from(ord0) + 1) ^ (from_ord + 1)) - 1) as i8;
        let smaller_is_second = i32::from(ord1 < ord0);
        let orders = [ord0, ord1];
        output.ord[0] = orders[smaller_is_second as usize];
        output.ord[1] = orders[usize::from(smaller_is_second == 0)];
        output.num[0] = i8::try_from(ring_size - 1).unwrap();
    }
    Ok(ring_size)
}

pub(crate) fn check_arom_chain(
    atoms: &[inp_ATOM],
    cur: i32,
    from: i32,
    last: i32,
    len: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2281 check_arom_chain
    // INCHI✔️❌: int check_arom_chain( inp_ATOM *at,
    // INCHI✔️❌:                       int cur /* first*/,
    // INCHI✔️❌:                       int from,
    // INCHI✔️❌:                       int last,
    // INCHI✔️❌:                       int len )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, num;
    // INCHI✔️❌:     num = 0;
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* check this on all except at[last], which is typically different */
    // INCHI✔️❌:         if (at[cur].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:              at[cur].valence != 2 ||
    // INCHI✔️❌:              at[cur].chem_bonds_valence != 3 ||
    // INCHI✔️❌:              at[cur].num_H != 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* bond to the next atom - check on all, cur..last, atoms */
    // INCHI✔️❌:         i = ( at[cur].neighbor[0] == from ); /* index of a bond to the next atom */
    // INCHI✔️❌:         if (at[cur].bond_type[i] != BOND_ALTERN)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         num++; /* checks are complete */
    // INCHI✔️❌:                /* prepare for the next atom */
    // INCHI✔️❌:         from = cur;
    // INCHI✔️❌:         cur = at[cur].neighbor[i];
    // INCHI✔️❌:     } while (cur != last && num < len);
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( cur == last && ++num == len );
    // INCHI✔️❌:
    // INCHI✔️❌: check_next_derivative:
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: check_arom_chain
    // BEGIN ACTIVE INCHI HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/extr_ct.h:215
    // INCHI✔️❌: #define BOND_ALTERN BOND_TYPE_ALTERN  /* 4 single/double                    */
    // END ACTIVE INCHI HEADER MACRO: BOND_ALTERN

    let mut cur = usize::try_from(cur).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut from = usize::try_from(from).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let last = usize::try_from(last).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut num = 0_i32;
    loop {
        let atom = atoms.get(cur).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.el_number != EL_NUMBER_C
            || atom.valence != 2
            || atom.chem_bonds_valence != 3
            || atom.num_H != 1
        {
            return Ok(0);
        }
        let next_order = usize::from(usize::from(atom.neighbor[0]) == from);
        if atom.bond_type[next_order] != BOND_TYPE_ALTERN as u8 {
            return Ok(0);
        }
        num = num.wrapping_add(1);
        from = cur;
        cur = usize::from(atom.neighbor[next_order]);
        if cur == last || num >= len {
            break;
        }
    }
    Ok(i32::from(cur == last && num.wrapping_add(1) == len))
}

#[allow(non_snake_case)]
pub(crate) fn is_Dansyl(
    atoms: &[inp_ATOM],
    cur_atom: i32,
    to_ord: i32,
    _da: Option<&mut DERIV_AT>,
    da1: Option<&mut DERIV_AT>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2334 is_Dansyl
    // INCHI✔️❌: int is_Dansyl( inp_ATOM *at,
    // INCHI✔️❌:                int cur_atom,
    // INCHI✔️❌:                int to_ord, DERIV_AT *da,
    // INCHI✔️❌:                DERIV_AT *da1 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, a, b, c, d, cj = -1, dj = -1, neigh, k, iS /* S */, iN /* N */;
    // INCHI✔️❌:     if (( (( at[cur_atom].el_number == EL_NUMBER_O || at[cur_atom].el_number == EL_NUMBER_S || at[cur_atom].el_number == EL_NUMBER_N ) &&
    // INCHI✔️❌:           at[cur_atom].valence == 2 && at[cur_atom].num_H == ( at[cur_atom].el_number == EL_NUMBER_N ) &&
    // INCHI✔️❌:           at[cur_atom].nNumAtInRingSystem == 1) ||
    // INCHI✔️❌:           (at[cur_atom].el_number == EL_NUMBER_N && at[cur_atom].valence == 3 && at[cur_atom].num_H == 0) ) &&
    // INCHI✔️❌:          at[cur_atom].valence == at[cur_atom].chem_bonds_valence &&
    // INCHI✔️❌:          at[iS /* S */ = at[cur_atom].neighbor[to_ord]].el_number == EL_NUMBER_S &&
    // INCHI✔️❌:          at[iS].valence == 4 && at[iS].chem_bonds_valence == 6) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* neighbors of S; 6=1+1+2+2, 1+1+1+3 only. Therefore, we do not need to count (=O) neighbors */
    // INCHI✔️❌:         for (i = 0, a = -1; i < at[iS].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_atom == ( neigh = at[iS].neighbor[i] ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[iS].bond_type[i] == BOND_DOUBLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[neigh].el_number != EL_NUMBER_O ||
    // INCHI✔️❌:                          at[neigh].valence != 1 ||
    // INCHI✔️❌:                          at[neigh].chem_bonds_valence != 2 ||
    // INCHI✔️❌:                          at[neigh].num_H || at[neigh].charge || at[neigh].radical)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[iS].bond_type[i] == BOND_SINGLE && a == -1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (at[neigh].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:                              at[neigh].nNumAtInRingSystem != 10 ||
    // INCHI✔️❌:                              at[neigh].valence != 3 ||
    // INCHI✔️❌:                              at[neigh].chem_bonds_valence != 4 ||
    // INCHI✔️❌:                              at[neigh].num_H || at[neigh].charge || at[neigh].radical ||
    // INCHI✔️❌:                              at[neigh].bond_type[at[neigh].neighbor[0] == iS] != BOND_ALTERN)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             goto check_next_derivative;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         a = neigh;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (a < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* at[a] - aromatic, j1 - index of its bond to its aromatic neighbor */
    // INCHI✔️❌:         /* find at[b] */
    // INCHI✔️❌:         for (i = k = 0, b = -1; i < at[a].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[a].neighbor[i];
    // INCHI✔️❌:             if (neigh == iS)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k += 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[a].bond_type[i] != BOND_ALTERN)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[neigh].valence == 3 && b == -1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     b = neigh;
    // INCHI✔️❌:                     k += 10;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[neigh].valence == 2 && at[neigh].chem_bonds_valence == 3)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k += 100;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* structure check: C[b] */
    // INCHI✔️❌:         if (k != 111)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[b].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:              at[b].valence != 3 ||
    // INCHI✔️❌:              at[b].chem_bonds_valence != 4 ||
    // INCHI✔️❌:              at[b].charge || at[b].radical || at[b].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* find at[c] */
    // INCHI✔️❌:         for (i = k = 0, c = -1; i < at[b].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[b].bond_type[i] != BOND_ALTERN)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             neigh = at[b].neighbor[i];
    // INCHI✔️❌:             if (neigh == a)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k += 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].valence == 3 && c == -1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     c = neigh;
    // INCHI✔️❌:                     k += 10;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[neigh].valence == 2 && at[neigh].chem_bonds_valence == 3)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k += 100;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (k != 111)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* structure check: C[c] */
    // INCHI✔️❌:         if (at[c].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:              at[c].valence != 3 ||
    // INCHI✔️❌:              at[c].chem_bonds_valence != 4 ||
    // INCHI✔️❌:              at[c].charge || at[c].radical || at[c].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* find at[d] */
    // INCHI✔️❌:         for (i = k = 0, d = -1; i < at[c].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[c].bond_type[i] != BOND_ALTERN)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             neigh = at[c].neighbor[i];
    // INCHI✔️❌:             if (neigh == b)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k += 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].valence == 3 && d == -1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     d = neigh;
    // INCHI✔️❌:                     k += 10;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[neigh].valence == 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         cj = i;
    // INCHI✔️❌:                         k += 100;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (k != 111)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* structure check: C[d] */
    // INCHI✔️❌:         if (at[d].el_number != EL_NUMBER_C ||
    // INCHI✔️❌:              at[d].valence != 3 ||
    // INCHI✔️❌:              at[d].chem_bonds_valence != 4 ||
    // INCHI✔️❌:              at[d].charge || at[d].radical || at[d].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* find at[iN]  */
    // INCHI✔️❌:         for (i = k = 0, iN = -1; i < at[d].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[d].neighbor[i];
    // INCHI✔️❌:             if (neigh == c)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k += 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].valence == 3 && iN == -1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iN = neigh;
    // INCHI✔️❌:                     k += 10;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[neigh].valence == 2 && at[neigh].chem_bonds_valence == 3 && at[d].bond_type[i] == BOND_ALTERN)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         dj = i;
    // INCHI✔️❌:                         k += 100;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* structure check: at[iN] */
    // INCHI✔️❌:         if (k != 111)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[iN].el_number != EL_NUMBER_N ||
    // INCHI✔️❌:              at[iN].valence != 3 ||
    // INCHI✔️❌:              at[iN].chem_bonds_valence != 3 ||
    // INCHI✔️❌:              at[iN].charge || at[d].radical || at[d].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* find attached to N 2 methyls */
    // INCHI✔️❌:         for (i = 0; i < at[iN].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (( neigh = at[iN].neighbor[i] ) != d &&
    // INCHI✔️❌:                  !is_Methyl( at, neigh ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative; /* not Dansyl */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* check aromatic chain d-dj-...b and c-cj-...a */
    // INCHI✔️❌:         if (check_arom_chain( at, at[d].neighbor[dj] /* first*/, d /*from*/, b /*to*/, 4 ) &&
    // INCHI✔️❌:              check_arom_chain( at, at[c].neighbor[cj] /* first*/, c /*from*/, a /*to*/, 4 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1->typ[0] = DERIV_DANSYL;
    // INCHI✔️❌:                 da1->ord[0] = to_ord;
    // INCHI✔️❌:                 da1->num[0] = 16;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return DERIV_DANSYL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: check_next_derivative:
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_Dansyl

    let cur_atom = usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let to_ord_index = usize::try_from(to_ord).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = atoms
        .get(cur_atom)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let first_shape = matches!(current.el_number, EL_NUMBER_O | EL_NUMBER_S | EL_NUMBER_N)
        && current.valence == 2
        && current.num_H == i8::from(current.el_number == EL_NUMBER_N)
        && current.nNumAtInRingSystem == 1;
    let second_shape =
        current.el_number == EL_NUMBER_N && current.valence == 3 && current.num_H == 0;
    let sulfur_index = usize::from(
        *current
            .neighbor
            .get(to_ord_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let sulfur = atoms
        .get(sulfur_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if !(first_shape || second_shape)
        || current.valence != current.chem_bonds_valence
        || sulfur.el_number != EL_NUMBER_S
        || sulfur.valence != 4
        || sulfur.chem_bonds_valence != 6
    {
        return Ok(0);
    }

    let mut a = None;
    for order in 0..usize::try_from(i32::from(sulfur.valence)).unwrap() {
        let neighbor_index = usize::from(sulfur.neighbor[order]);
        if neighbor_index == cur_atom {
            continue;
        }
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if sulfur.bond_type[order] == BOND_DOUBLE as u8 {
            if neighbor.el_number != EL_NUMBER_O
                || neighbor.valence != 1
                || neighbor.chem_bonds_valence != 2
                || neighbor.num_H != 0
                || neighbor.charge != 0
                || neighbor.radical != 0
            {
                return Ok(0);
            }
        } else if sulfur.bond_type[order] == BOND_SINGLE as u8 && a.is_none() {
            let aromatic_order = usize::from(usize::from(neighbor.neighbor[0]) == sulfur_index);
            if neighbor.el_number != EL_NUMBER_C
                || neighbor.nNumAtInRingSystem != 10
                || neighbor.valence != 3
                || neighbor.chem_bonds_valence != 4
                || neighbor.num_H != 0
                || neighbor.charge != 0
                || neighbor.radical != 0
                || neighbor.bond_type[aromatic_order] != BOND_TYPE_ALTERN as u8
            {
                return Ok(0);
            }
            a = Some(neighbor_index);
        } else {
            return Ok(0);
        }
    }
    let Some(a) = a else { return Ok(0) };

    let atom_a = atoms.get(a).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut k = 0_i32;
    let mut b = None;
    for order in 0..usize::try_from(i32::from(atom_a.valence)).unwrap() {
        let neighbor_index = usize::from(atom_a.neighbor[order]);
        if neighbor_index == sulfur_index {
            k += 1;
        } else {
            if atom_a.bond_type[order] != BOND_TYPE_ALTERN as u8 {
                return Ok(0);
            }
            let neighbor = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.valence == 3 && b.is_none() {
                b = Some(neighbor_index);
                k += 10;
            } else if neighbor.valence == 2 && neighbor.chem_bonds_valence == 3 {
                k += 100;
            } else {
                return Ok(0);
            }
        }
    }
    if k != 111 {
        return Ok(0);
    }
    let b = b.unwrap();
    let atom_b = atoms.get(b).ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom_b.el_number != EL_NUMBER_C
        || atom_b.valence != 3
        || atom_b.chem_bonds_valence != 4
        || atom_b.charge != 0
        || atom_b.radical != 0
        || atom_b.num_H != 0
    {
        return Ok(0);
    }

    k = 0;
    let mut c = None;
    for order in 0..usize::try_from(i32::from(atom_b.valence)).unwrap() {
        if atom_b.bond_type[order] != BOND_TYPE_ALTERN as u8 {
            return Ok(0);
        }
        let neighbor_index = usize::from(atom_b.neighbor[order]);
        if neighbor_index == a {
            k += 1;
        } else {
            let neighbor = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.valence == 3 && c.is_none() {
                c = Some(neighbor_index);
                k += 10;
            } else if neighbor.valence == 2 && neighbor.chem_bonds_valence == 3 {
                k += 100;
            } else {
                return Ok(0);
            }
        }
    }
    if k != 111 {
        return Ok(0);
    }
    let c = c.unwrap();
    let atom_c = atoms.get(c).ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom_c.el_number != EL_NUMBER_C
        || atom_c.valence != 3
        || atom_c.chem_bonds_valence != 4
        || atom_c.charge != 0
        || atom_c.radical != 0
        || atom_c.num_H != 0
    {
        return Ok(0);
    }

    k = 0;
    let mut d = None;
    let mut cj = None;
    for order in 0..usize::try_from(i32::from(atom_c.valence)).unwrap() {
        if atom_c.bond_type[order] != BOND_TYPE_ALTERN as u8 {
            return Ok(0);
        }
        let neighbor_index = usize::from(atom_c.neighbor[order]);
        if neighbor_index == b {
            k += 1;
        } else {
            let neighbor = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.valence == 3 && d.is_none() {
                d = Some(neighbor_index);
                k += 10;
            } else if neighbor.valence == 2 {
                cj = Some(order);
                k += 100;
            } else {
                return Ok(0);
            }
        }
    }
    if k != 111 {
        return Ok(0);
    }
    let d = d.unwrap();
    let cj = cj.unwrap();
    let atom_d = atoms.get(d).ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom_d.el_number != EL_NUMBER_C
        || atom_d.valence != 3
        || atom_d.chem_bonds_valence != 4
        || atom_d.charge != 0
        || atom_d.radical != 0
        || atom_d.num_H != 0
    {
        return Ok(0);
    }

    k = 0;
    let mut nitrogen = None;
    let mut dj = None;
    for order in 0..usize::try_from(i32::from(atom_d.valence)).unwrap() {
        let neighbor_index = usize::from(atom_d.neighbor[order]);
        if neighbor_index == c {
            k += 1;
        } else {
            let neighbor = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.valence == 3 && nitrogen.is_none() {
                nitrogen = Some(neighbor_index);
                k += 10;
            } else if neighbor.valence == 2
                && neighbor.chem_bonds_valence == 3
                && atom_d.bond_type[order] == BOND_TYPE_ALTERN as u8
            {
                dj = Some(order);
                k += 100;
            } else {
                return Ok(0);
            }
        }
    }
    if k != 111 {
        return Ok(0);
    }
    let nitrogen = nitrogen.unwrap();
    let dj = dj.unwrap();
    let atom_n = atoms
        .get(nitrogen)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom_n.el_number != EL_NUMBER_N
        || atom_n.valence != 3
        || atom_n.chem_bonds_valence != 3
        || atom_n.charge != 0
        || atom_d.radical != 0
        || atom_d.num_H != 0
    {
        return Ok(0);
    }
    for order in 0..usize::try_from(i32::from(atom_n.valence)).unwrap() {
        let neighbor = usize::from(atom_n.neighbor[order]);
        if neighbor != d && is_Methyl(atoms, neighbor as i32)? == 0 {
            return Ok(0);
        }
    }

    if check_arom_chain(atoms, i32::from(atom_d.neighbor[dj]), d as i32, b as i32, 4)? != 0
        && check_arom_chain(atoms, i32::from(atom_c.neighbor[cj]), c as i32, a as i32, 4)? != 0
    {
        if let Some(output) = da1 {
            output.typ[0] = DERIV_DANSYL as i16;
            output.ord[0] = to_ord as i8;
            output.num[0] = 16;
        }
        return Ok(DERIV_DANSYL as i32);
    }
    Ok(0)
}

pub(crate) fn is_silyl2(
    atoms: &[inp_ATOM],
    start: i32,
    from_atom: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3738 is_silyl2
    // INCHI✔️❌: int is_silyl2( inp_ATOM *at, int start, int from_at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, neigh, num_Me = 0, iC_IV = -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[start].el_number != EL_NUMBER_SI || at[start].valence != 4 ||
    // INCHI✔️❌:          at[start].valence != at[start].chem_bonds_valence ||
    // INCHI✔️❌:          at[start].charge || at[start].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[start].neighbor[i];
    // INCHI✔️❌:         if (neigh == from_at)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].charge || at[neigh].radical ||
    // INCHI✔️❌:              at[neigh].valence != at[neigh].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].valence == 4)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[neigh].el_number == EL_NUMBER_C && iC_IV < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 iC_IV = neigh;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[neigh].valence == 1 &&
    // INCHI✔️❌:                  at[neigh].valence == at[neigh].chem_bonds_valence &&
    // INCHI✔️❌:                  at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_Me++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (num_Me == 3 && iC_IV < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* Si(CH3)3 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* next C(IV) */
    // INCHI✔️❌:     if (num_Me != 2 || iC_IV < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     num_Me = 0;
    // INCHI✔️❌:     for (i = 0; i < at[iC_IV].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[iC_IV].neighbor[i];
    // INCHI✔️❌:         if (neigh == start)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].charge || at[neigh].radical ||
    // INCHI✔️❌:              at[neigh].valence != at[neigh].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].valence == 1 &&
    // INCHI✔️❌:              at[neigh].valence == at[neigh].chem_bonds_valence &&
    // INCHI✔️❌:              at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_Me++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (num_Me == 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 2; /* Si(CH3)2-Si/C(CH3)3, not Si(CH3)2-Si(CH3)3 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: is_silyl2

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(start)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.el_number != EL_NUMBER_SI
        || atom.valence != 4
        || atom.valence != atom.chem_bonds_valence
        || atom.charge != 0
        || atom.radical != 0
    {
        return Ok(0);
    }

    let mut methyls = 0_i32;
    let mut carbon_iv = None;
    for order in 0..usize::try_from(i32::from(atom.valence)).unwrap_or(0) {
        let neighbor_index = usize::from(atom.neighbor[order]);
        if neighbor_index as i32 == from_atom {
            continue;
        }
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.charge != 0
            || neighbor.radical != 0
            || neighbor.valence != neighbor.chem_bonds_valence
        {
            return Ok(0);
        }
        if neighbor.valence == 4 {
            if neighbor.el_number == EL_NUMBER_C && carbon_iv.is_none() {
                carbon_iv = Some(neighbor_index);
            } else {
                return Ok(0);
            }
        } else if neighbor.valence == 1
            && neighbor.valence == neighbor.chem_bonds_valence
            && neighbor.el_number == EL_NUMBER_C
            && neighbor.num_H == 3
        {
            methyls += 1;
        } else {
            return Ok(0);
        }
    }
    if methyls == 3 && carbon_iv.is_none() {
        return Ok(1);
    }
    if methyls != 2 {
        return Ok(0);
    }
    let Some(carbon_iv) = carbon_iv else {
        return Ok(0);
    };

    methyls = 0;
    let carbon = atoms
        .get(carbon_iv)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for order in 0..usize::try_from(i32::from(carbon.valence)).unwrap_or(0) {
        let neighbor_index = usize::from(carbon.neighbor[order]);
        if neighbor_index == start {
            continue;
        }
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.charge != 0
            || neighbor.radical != 0
            || neighbor.valence != neighbor.chem_bonds_valence
        {
            return Ok(0);
        }
        if neighbor.valence == 1
            && neighbor.valence == neighbor.chem_bonds_valence
            && neighbor.el_number == EL_NUMBER_C
            && neighbor.num_H == 3
        {
            methyls += 1;
        } else {
            return Ok(0);
        }
    }
    Ok(i32::from(methyls == 3) * 2)
}

pub(crate) fn is_silyl(
    atoms: &[inp_ATOM],
    start: i32,
    previous_order: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3633 is_silyl
    // INCHI✔️❌: int is_silyl( inp_ATOM *at, int start, int ord_prev )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, neigh, num_Me = 0, iC_IV = -1, iSi_IV = -1, i_C_or_Si_IV;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[start].el_number != EL_NUMBER_SI || at[start].valence != 4 ||
    // INCHI✔️❌:          at[start].valence != at[start].chem_bonds_valence ||
    // INCHI✔️❌:          at[start].charge || at[start].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i == ord_prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         neigh = at[start].neighbor[i];
    // INCHI✔️❌:         if (at[neigh].charge || at[neigh].radical ||
    // INCHI✔️❌:              at[neigh].valence != at[neigh].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].valence == 4)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[neigh].el_number == EL_NUMBER_C && iC_IV < 0 && iSi_IV < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 iC_IV = neigh;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].el_number == EL_NUMBER_SI && iC_IV < 0 && iSi_IV < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iSi_IV = neigh;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[neigh].valence == 1 &&
    // INCHI✔️❌:                  at[neigh].valence == at[neigh].chem_bonds_valence &&
    // INCHI✔️❌:                  at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_Me++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_Me == 3 && iC_IV < 0 && iSi_IV < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* Si(CH3)3 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* next C(IV) or Si(IV) */
    // INCHI✔️❌:     /* this is a fix requested by Anz. and suggested by Gary 09/21/2011
    // INCHI✔️❌:     it rejects -Si(CH3)2-Si(CH3)3 and allows only -Si(CH3)2-C(CH3)3
    // INCHI✔️❌:     */
    // INCHI✔️❌:     i_C_or_Si_IV = iC_IV >= 0 ? iC_IV : -1;
    // INCHI✔️❌:     if (num_Me != 2 || i_C_or_Si_IV < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     num_Me = 0;
    // INCHI✔️❌:     for (i = 0; i < at[i_C_or_Si_IV].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[i_C_or_Si_IV].neighbor[i];
    // INCHI✔️❌:         if (neigh == start)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].charge || at[neigh].radical ||
    // INCHI✔️❌:              at[neigh].valence != at[neigh].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[neigh].valence == 1 &&
    // INCHI✔️❌:              at[neigh].valence == at[neigh].chem_bonds_valence &&
    // INCHI✔️❌:              at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_Me++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (num_Me == 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 2; /* Si(CH3)2Si/C(CH3)3 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_silyl

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(start)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.el_number != EL_NUMBER_SI
        || atom.valence != 4
        || atom.valence != atom.chem_bonds_valence
        || atom.charge != 0
        || atom.radical != 0
    {
        return Ok(0);
    }
    let mut methyls = 0_i32;
    let mut carbon_iv = None;
    let mut silicon_iv = None;
    for order in 0..usize::try_from(i32::from(atom.valence)).unwrap() {
        if order as i32 == previous_order {
            continue;
        }
        let neighbor_index = usize::from(atom.neighbor[order]);
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.charge != 0
            || neighbor.radical != 0
            || neighbor.valence != neighbor.chem_bonds_valence
        {
            return Ok(0);
        }
        if neighbor.valence == 4 {
            if neighbor.el_number == EL_NUMBER_C && carbon_iv.is_none() && silicon_iv.is_none() {
                carbon_iv = Some(neighbor_index);
            } else if neighbor.el_number == EL_NUMBER_SI
                && carbon_iv.is_none()
                && silicon_iv.is_none()
            {
                silicon_iv = Some(neighbor_index);
            } else {
                return Ok(0);
            }
        } else if neighbor.valence == 1
            && neighbor.valence == neighbor.chem_bonds_valence
            && neighbor.el_number == EL_NUMBER_C
            && neighbor.num_H == 3
        {
            methyls += 1;
        } else {
            return Ok(0);
        }
    }
    if methyls == 3 && carbon_iv.is_none() && silicon_iv.is_none() {
        return Ok(1);
    }
    let Some(carbon_iv) = carbon_iv else {
        return Ok(0);
    };
    if methyls != 2 {
        return Ok(0);
    }
    methyls = 0;
    let carbon = atoms
        .get(carbon_iv)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for order in 0..usize::try_from(i32::from(carbon.valence)).unwrap() {
        let neighbor_index = usize::from(carbon.neighbor[order]);
        if neighbor_index == start {
            continue;
        }
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.charge != 0
            || neighbor.radical != 0
            || neighbor.valence != neighbor.chem_bonds_valence
        {
            return Ok(0);
        }
        if neighbor.valence == 1
            && neighbor.valence == neighbor.chem_bonds_valence
            && neighbor.el_number == EL_NUMBER_C
            && neighbor.num_H == 3
        {
            methyls += 1;
        } else {
            return Ok(0);
        }
    }
    Ok(if methyls == 3 { 2 } else { 0 })
}

#[allow(non_snake_case)]
pub(crate) fn is_CF3_or_linC3F7a(
    atoms: &[inp_ATOM],
    start: i32,
    previous_atom: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4000 is_CF3_or_linC3F7a
    // INCHI✔️❌: int is_CF3_or_linC3F7a( inp_ATOM *at, int start, int iat_prev )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (iat_prev == at[start].neighbor[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return is_CF3_or_linC3F7( at, start, i );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_CF3_or_linC3F7a

    let start_index = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(start_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for order in 0..usize::try_from(i32::from(atom.valence)).unwrap_or(0) {
        if i32::from(atom.neighbor[order]) == previous_atom {
            return is_CF3_or_linC3F7(atoms, start, order as i32);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn is_CF3_or_linC3F7(
    atoms: &[inp_ATOM],
    start: i32,
    previous_order: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4017 is_CF3_or_linC3F7
    // INCHI✔️❌: int is_CF3_or_linC3F7( inp_ATOM *at, int start, int ord_prev )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, num_C_IV, iC_IV, neigh, num_F, num_C = 0;
    // INCHI✔️❌:     AT_NUMB *p;
    // INCHI✔️❌:
    // INCHI✔️❌:     while (num_C < 4)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         if (at[start].el_number != EL_NUMBER_C || at[start].valence != 4 ||
    // INCHI✔️❌:              at[start].valence != at[start].chem_bonds_valence ||
    // INCHI✔️❌:              at[start].chem_bonds_valence + at[start].num_H != 4 ||
    // INCHI✔️❌:              at[start].charge || at[start].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         iC_IV = num_F = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = num_C_IV = 0; i < at[start].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i == ord_prev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             neigh = at[start].neighbor[i];
    // INCHI✔️❌:             if (at[neigh].valence != at[neigh].chem_bonds_valence ||
    // INCHI✔️❌:                  at[neigh].charge || at[neigh].radical)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (at[neigh].el_number == EL_NUMBER_F)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].chem_bonds_valence + at[neigh].num_H != 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 num_F++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[neigh].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:                      at[neigh].valence == 4 &&
    // INCHI✔️❌:                      !at[neigh].num_H && !at[neigh].charge && !at[neigh].radical && num_C_IV < 2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (num_C_IV)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     iC_IV = neigh;
    // INCHI✔️❌:                     num_C_IV++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (num_C_IV + num_F != 3)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         num_C++; /* found -CF2-C or -CF3 */
    // INCHI✔️❌:         if (!num_C_IV)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break; /* -CF3 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* treat next C(IV) as a new start atom */
    // INCHI✔️❌:         if ((p = is_in_the_list( at[iC_IV].neighbor, (AT_NUMB) start, at[iC_IV].valence ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             start = iC_IV;
    // INCHI✔️❌:             ord_prev = p - at[iC_IV].neighbor;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return -1; /* program error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Corrected by DT below - ? was - return num_C == 1 ? 1 : num_C == 3 ? 2 : 0;*/
    // INCHI✔️❌:     return num_C == 1 ? 1 : num_C == 2 ? 2 : num_C == 3 ? 3 : 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_CF3_or_linC3F7

    let mut start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut previous_order = previous_order;
    let mut carbon_count = 0_i32;
    while carbon_count < 4 {
        let atom = atoms
            .get(start)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.el_number != EL_NUMBER_C
            || atom.valence != 4
            || atom.valence != atom.chem_bonds_valence
            || i32::from(atom.chem_bonds_valence) + i32::from(atom.num_H) != 4
            || atom.charge != 0
            || atom.radical != 0
        {
            return Ok(0);
        }
        let mut next_carbon = 0_usize;
        let mut fluorines = 0_i32;
        let mut carbon_neighbors = 0_i32;
        for order in 0..usize::try_from(i32::from(atom.valence)).unwrap() {
            if order as i32 == previous_order {
                continue;
            }
            let neighbor_index = usize::from(atom.neighbor[order]);
            let neighbor = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.valence != neighbor.chem_bonds_valence
                || neighbor.charge != 0
                || neighbor.radical != 0
            {
                return Ok(0);
            }
            if neighbor.el_number == EL_NUMBER_F {
                if i32::from(neighbor.chem_bonds_valence) + i32::from(neighbor.num_H) != 1 {
                    return Ok(0);
                }
                fluorines += 1;
            } else if neighbor.el_number == EL_NUMBER_C
                && neighbor.valence == 4
                && neighbor.num_H == 0
                && neighbor.charge == 0
                && neighbor.radical == 0
                && carbon_neighbors < 2
            {
                if carbon_neighbors != 0 {
                    return Ok(0);
                }
                next_carbon = neighbor_index;
                carbon_neighbors += 1;
            }
        }
        if carbon_neighbors + fluorines != 3 {
            return Ok(0);
        }
        carbon_count += 1;
        if carbon_neighbors == 0 {
            break;
        }
        let next = atoms
            .get(next_carbon)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let Some(order) =
            is_in_the_list(Some(&next.neighbor), start as u16, i32::from(next.valence))?
        else {
            return Ok(-1);
        };
        start = next_carbon;
        previous_order = order as i32;
    }
    Ok(match carbon_count {
        1 => 1,
        2 => 2,
        3 => 3,
        _ => 0,
    })
}

pub(crate) fn underiv_buf_clear(buffer: Option<&mut [i8]>) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5147 underiv_buf_clear
    // INCHI✔️❌: void underiv_buf_clear( char *szUnderiv )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (NULL != szUnderiv)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szUnderiv[0] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: underiv_buf_clear
    if let Some(buffer) = buffer {
        *buffer
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    }
    Ok(())
}

pub(crate) fn underiv_list_add(
    list: Option<&mut [i8]>,
    list_capacity: i32,
    added: Option<&[i8]>,
    delimiter: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5156 underiv_list_add
    // INCHI✔️❌: int underiv_list_add( char *szUnderivList, int lenUnderivList, const char *szUnderiv, char cDelimiter )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (NULL != szUnderivList && lenUnderivList > 0 && NULL != szUnderiv && *szUnderiv)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int lenList = strlen( szUnderivList );
    // INCHI✔️❌:         int lenAdd = strlen( szUnderiv );
    // INCHI✔️❌:         int bDelimiter = cDelimiter ? 1 : 0;
    // INCHI✔️❌:         if (lenList + lenAdd + bDelimiter < lenUnderivList)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (lenList && bDelimiter)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szUnderivList[lenList++] = cDelimiter;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             memcpy(szUnderivList + lenList, szUnderiv, (long long)lenAdd + 1); /* +1 adds zero termination */ /* djb-rwth: cast operator added */
    // INCHI✔️❌:             return lenList + lenAdd;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: underiv_list_add

    let Some(list) = list else { return Ok(0) };
    if list_capacity <= 0 {
        return Ok(0);
    }
    let Some(added) = added else { return Ok(0) };
    if added.first().copied().unwrap_or(0) == 0 {
        return Ok(0);
    }
    let list_length = list
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let added_length = added
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let delimiter_length = usize::from(delimiter != 0);
    let required = list_length
        .checked_add(added_length)
        .and_then(|value| value.checked_add(delimiter_length))
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if i32::try_from(required).map_err(|_| SourceHeapError::SourceIntegerOverflow)? >= list_capacity
    {
        return Ok(0);
    }
    let mut write_at = list_length;
    if list_length != 0 && delimiter != 0 {
        *list
            .get_mut(write_at)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = delimiter;
        write_at += 1;
    }
    let write_end = write_at
        .checked_add(added_length + 1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let destination = list
        .get_mut(write_at..write_end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    destination.copy_from_slice(&added[..=added_length]);
    i32::try_from(write_at + added_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

pub(crate) fn underiv_list_get_last(
    list: Option<&[i8]>,
    delimiter: i8,
) -> Result<Option<usize>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5179 underiv_list_get_last
    // INCHI✔️❌: const char* underiv_list_get_last( const char *szUnderivList,
    // INCHI✔️❌:                                    char cDelimiter )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (szUnderivList && cDelimiter)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         const char *p = strrchr( szUnderivList, cDelimiter );
    // INCHI✔️❌:         if (NULL == p)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             p = szUnderivList;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return p;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return NULL;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: underiv_list_get_last

    let Some(list) = list else { return Ok(None) };
    if delimiter == 0 {
        return Ok(None);
    }
    let length = list
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(Some(
        list[..length]
            .iter()
            .rposition(|&byte| byte == delimiter)
            .unwrap_or(0),
    ))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn mark_atoms_deriv(
    atoms: &mut [inp_ATOM],
    derivatives: &mut [DERIV_AT],
    start: i32,
    mut num: i32,
    flags: i8,
    found: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3332 mark_atoms_deriv
    // INCHI✔️❌: int mark_atoms_deriv( inp_ATOM *at,
    // INCHI✔️❌:                       DERIV_AT *da,
    // INCHI✔️❌:                       int start,
    // INCHI✔️❌:                       int num,
    // INCHI✔️❌:                       char cFlags,
    // INCHI✔️❌:                       int *pbFound )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, nFound = 0, ret; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     DERIV_AT da1;
    // INCHI✔️❌:     int      ret2;   /* moved from below 2024-09-01 DT */
    // INCHI✔️❌:     DERIV_AT da2;    /* moved from below 2024-09-01 DT */
    // INCHI✔️❌:     da1.other_atom = 0; /* djb-rwth: initialisation needed for if conditons */
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:     /* djb-rwth: initialisation needed to avoid garbage values in add_to_da function call; fixing coverity ID #499492 */
    // INCHI✔️❌:     memset(da2.typ, 0, DERIV_AT_LEN * sizeof(da2.typ[0]));
    // INCHI✔️❌:     memset(da2.ord, '\0', DERIV_AT_LEN * sizeof(da2.ord[0]));
    // INCHI✔️❌:     memset(da2.num, '\0', DERIV_AT_LEN * sizeof(da2.num[0]));
    // INCHI✔️❌:     da2.other_atom = 0; /* djb-rwth: initialisation needed for if conditons */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     if (!( at[start].cFlags & cFlags ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (DERIV_NOT == ( ret = get_traversed_deriv_type( at, da, start, &da1, cFlags ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nFound++; /* at[start] cannot belong to a derivatizing agent */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         num++;
    // INCHI✔️❌:         at[start].cFlags |= cFlags;
    // INCHI✔️❌:         if (da1.typ[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* possibly a derivatization agent attachment point. */
    // INCHI✔️❌:             /* check neighbors that have not been traversed yet */
    // INCHI✔️❌:             int n1 = 0, n2 = 0, i1 = -1, i2 = -1, nFound1 = 0, nFound2 = 0;
    // INCHI✔️❌:             switch (da1.typ[0])
    // INCHI✔️❌:             {
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:                 case DERIV_RING_DMOX_DEOX_N:
    // INCHI✔️❌:                 case DERIV_RING_DMOX_DEOX_O:
    // INCHI✔️❌:                     ret2 = get_traversed_deriv_type( at, da, da1.other_atom - 1, &da2, cFlags );
    // INCHI✔️❌:                     if (ret != ( ret2 ^ DERIV_RING_DMOX_DEOX ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* bug */
    // INCHI✔️❌:                         nFound++; /* at[start] cannot belong to a derivatizing agent */
    // INCHI✔️❌:                         goto check_neighbors; /* bypass add_to_da( da+start, &da1 ) */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /*at[da1.other_atom-1].cFlags |= cFlags;*/
    // INCHI✔️❌:                     n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
    // INCHI✔️❌:                     nFound++;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_X_OXIME
    // INCHI✔️❌:                 case DERIV_X_OXIME:
    // INCHI✔️❌:                     n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
    // INCHI✔️❌:                     nFound++;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_RO_COX
    // INCHI✔️❌:                 case DERIV_RO_COX:
    // INCHI✔️❌:                     n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
    // INCHI✔️❌:                     nFound++;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) || defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) || defined(DERIV_DANSYL) )
    // INCHI✔️❌: #ifdef DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
    // INCHI✔️❌:                 case DERIV_RING2_PRRLDD_OUTSIDE_PRECUR:
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_RING2_PPRDN_OUTSIDE_PRECUR
    // INCHI✔️❌:                 case DERIV_RING2_PPRDN_OUTSIDE_PRECUR:
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_DANSYL
    // INCHI✔️❌:                 case DERIV_DANSYL:
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
    // INCHI✔️❌:                     nFound++;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                 case DERIV_BRIDGE_O:
    // INCHI✔️❌:                 case DERIV_BRIDGE_NH:
    // INCHI✔️❌:                     n1 = mark_atoms_deriv( at, da, at[start].neighbor[(int) da1.ord[0]], 0, cFlags, &nFound1 );
    // INCHI✔️❌:                     if (n1 > MAX_AT_DERIV || nFound1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.typ[0] = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.num[0] = n1;
    // INCHI✔️❌:                         nFound++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case DERIV_AMINE_tN:
    // INCHI✔️❌:                     n1 = mark_atoms_deriv( at, da, at[start].neighbor[i1 = da1.ord[0]], 0, cFlags, &nFound1 ); /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                     if (da1.typ[1])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n2 = mark_atoms_deriv( at, da, at[start].neighbor[i2 = da1.ord[1]], 0, cFlags, &nFound2 ); /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (0 < n1 && n1 <= MAX_AT_DERIV && !nFound1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.num[0] = n1;
    // INCHI✔️❌:                         i = 1;
    // INCHI✔️❌:                         nFound++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.ord[0] = da1.ord[1];
    // INCHI✔️❌:                         da1.num[0] = da1.num[1];
    // INCHI✔️❌:                         da1.typ[0] = da1.typ[1];
    // INCHI✔️❌:                         da1.typ[1] = 0;
    // INCHI✔️❌:                         i = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (0 < n2 && n2 <= MAX_AT_DERIV && !nFound2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.num[i] = n2;
    // INCHI✔️❌:                         nFound++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.typ[i] = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case DERIV_RING_O_OUTSIDE_PRECURSOR:
    // INCHI✔️❌:                 case DERIV_RING_NH_OUTSIDE_PRECURSOR:
    // INCHI✔️❌:                     for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (i != da1.ord[0] && i != da1.ord[1] && !( at[at[start].neighbor[i]].cFlags & cFlags ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (!n1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 n1 = mark_atoms_deriv( at, da, at[start].neighbor[i1 = i], 0, cFlags, &nFound1 ); /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 n2 = mark_atoms_deriv( at, da, at[start].neighbor[i2 = i], 0, cFlags, &nFound2 ); /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (n1 + n2 + 1 > MAX_AT_DERIV || nFound1 || nFound2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.typ[1] = da1.typ[0] = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1.num[0] = n1;
    // INCHI✔️❌:                         da1.num[1] = n2;
    // INCHI✔️❌:                         nFound++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* --------- end of switch( da1.typ[0] ) ------- */
    // INCHI✔️❌:             if (n1 < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return n1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (n2 < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return n2; /* errors */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:             if (da1.other_atom)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (da2.other_atom == start + 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if ((i = add_to_da( da + da1.other_atom - 1, &da2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return i;  /* error */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -4; /* error */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             if ((i = add_to_da( da + start, &da1 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return i;  /* error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             *pbFound += nFound1 + nFound2 + nFound;
    // INCHI✔️❌:             num += n1 + n2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *pbFound += nFound;
    // INCHI✔️❌:         }
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:         check_neighbors:
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                        for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:                        {
    // INCHI✔️❌:                            num = mark_atoms_deriv( at, da, at[start].neighbor[i], num, cFlags, pbFound );
    // INCHI✔️❌:                            if (num < 0)
    // INCHI✔️❌:                            {
    // INCHI✔️❌:                                return num;
    // INCHI✔️❌:                            }
    // INCHI✔️❌:                        }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* *pbFound =  number of derivatizer attachment points traversed forward from at[start] */
    // INCHI✔️❌:
    // INCHI✔️❌:     return num; /* number of atoms traversed forward from at[start] */
    // END INCHI C FUNCTION: mark_atoms_deriv

    let start_index = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = atoms
        .get(start_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if current.cFlags & flags != 0 {
        return Ok(num);
    }

    let mut da1 = DERIV_AT::default();
    let derivative_type = get_traversed_deriv_type(atoms, derivatives, start, &mut da1, flags)?;
    let mut found_here = i32::from(derivative_type == DERIV_NOT as i32);
    num = num.wrapping_add(1);
    atoms[start_index].cFlags |= flags;
    let mut da2 = DERIV_AT::default();
    if da1.typ[0] != 0 {
        let mut first_count = 0_i32;
        let mut second_count = 0_i32;
        let mut found_first = 0_i32;
        let mut found_second = 0_i32;
        let mut bypass_add = false;
        match i32::from(da1.typ[0]) {
            value
                if value == DERIV_RING_DMOX_DEOX_N as i32
                    || value == DERIV_RING_DMOX_DEOX_O as i32 =>
            {
                let other = usize::from(da1.other_atom)
                    .checked_sub(1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let ret2 =
                    get_traversed_deriv_type(atoms, derivatives, other as i32, &mut da2, flags)?;
                if derivative_type != (ret2 ^ DERIV_RING_DMOX_DEOX as i32) {
                    found_here = found_here.wrapping_add(1);
                    bypass_add = true;
                } else {
                    first_count = i32::from(da1.num[0]);
                    found_here = found_here.wrapping_add(1);
                }
            }
            value
                if value == DERIV_X_OXIME as i32
                    || value == DERIV_RO_COX as i32
                    || value == DERIV_RING2_PRRLDD_OUTSIDE_PRECUR as i32
                    || value == DERIV_RING2_PPRDN_OUTSIDE_PRECUR as i32
                    || value == DERIV_DANSYL as i32 =>
            {
                first_count = i32::from(da1.num[0]);
                found_here = found_here.wrapping_add(1);
            }
            value if value == DERIV_BRIDGE_O as i32 || value == DERIV_BRIDGE_NH as i32 => {
                let order = usize::try_from(i32::from(da1.ord[0]))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let neighbor = current
                    .neighbor
                    .get(order)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                first_count = mark_atoms_deriv(
                    atoms,
                    derivatives,
                    i32::from(neighbor),
                    0,
                    flags,
                    &mut found_first,
                )?;
                if first_count > MAX_AT_DERIV as i32 || found_first != 0 {
                    da1.typ[0] = 0;
                } else {
                    da1.num[0] = first_count as i8;
                    found_here = found_here.wrapping_add(1);
                }
            }
            value if value == DERIV_AMINE_tN as i32 => {
                let first_order = usize::try_from(i32::from(da1.ord[0]))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let first_neighbor = current
                    .neighbor
                    .get(first_order)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                first_count = mark_atoms_deriv(
                    atoms,
                    derivatives,
                    i32::from(first_neighbor),
                    0,
                    flags,
                    &mut found_first,
                )?;
                if da1.typ[1] != 0 {
                    let second_order = usize::try_from(i32::from(da1.ord[1]))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let second_neighbor = current
                        .neighbor
                        .get(second_order)
                        .copied()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    second_count = mark_atoms_deriv(
                        atoms,
                        derivatives,
                        i32::from(second_neighbor),
                        0,
                        flags,
                        &mut found_second,
                    )?;
                }
                let mut slot = 0_usize;
                if first_count > 0 && first_count <= MAX_AT_DERIV as i32 && found_first == 0 {
                    da1.num[0] = first_count as i8;
                    slot = 1;
                    found_here = found_here.wrapping_add(1);
                } else {
                    da1.ord[0] = da1.ord[1];
                    da1.num[0] = da1.num[1];
                    da1.typ[0] = da1.typ[1];
                    da1.typ[1] = 0;
                }
                if second_count > 0 && second_count <= MAX_AT_DERIV as i32 && found_second == 0 {
                    da1.num[slot] = second_count as i8;
                    found_here = found_here.wrapping_add(1);
                } else {
                    da1.typ[slot] = 0;
                }
            }
            value
                if value == DERIV_RING_O_OUTSIDE_PRECURSOR as i32
                    || value == DERIV_RING_NH_OUTSIDE_PRECURSOR as i32 =>
            {
                for order in 0..usize::try_from(i32::from(current.valence)).unwrap_or(0) {
                    if order != usize::try_from(i32::from(da1.ord[0])).unwrap_or(usize::MAX)
                        && order != usize::try_from(i32::from(da1.ord[1])).unwrap_or(usize::MAX)
                        && atoms
                            .get(usize::from(current.neighbor[order]))
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .cFlags
                            & flags
                            == 0
                    {
                        if first_count == 0 {
                            first_count = mark_atoms_deriv(
                                atoms,
                                derivatives,
                                i32::from(current.neighbor[order]),
                                0,
                                flags,
                                &mut found_first,
                            )?;
                        } else {
                            second_count = mark_atoms_deriv(
                                atoms,
                                derivatives,
                                i32::from(current.neighbor[order]),
                                0,
                                flags,
                                &mut found_second,
                            )?;
                        }
                    }
                }
                if first_count.wrapping_add(second_count).wrapping_add(1) > MAX_AT_DERIV as i32
                    || found_first != 0
                    || found_second != 0
                {
                    da1.typ[0] = 0;
                    da1.typ[1] = 0;
                } else {
                    da1.num[0] = first_count as i8;
                    da1.num[1] = second_count as i8;
                    found_here = found_here.wrapping_add(1);
                }
            }
            _ => {}
        }

        if !bypass_add {
            if first_count < 0 {
                return Ok(first_count);
            }
            if second_count < 0 {
                return Ok(second_count);
            }
            if da1.other_atom != 0 {
                let other = usize::from(da1.other_atom)
                    .checked_sub(1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if da2.other_atom == (start_index as u16).wrapping_add(1) {
                    let error = add_to_da(
                        derivatives
                            .get_mut(other)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        &da2,
                    );
                    if error != 0 {
                        return Ok(error);
                    }
                } else {
                    return Ok(-4);
                }
            }
            let error = add_to_da(
                derivatives
                    .get_mut(start_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                &da1,
            );
            if error != 0 {
                return Ok(error);
            }
            *found = found
                .wrapping_add(found_first)
                .wrapping_add(found_second)
                .wrapping_add(found_here);
            num = num.wrapping_add(first_count).wrapping_add(second_count);
        }
    } else {
        *found = found.wrapping_add(found_here);
    }

    let valence = usize::try_from(i32::from(current.valence)).unwrap_or(0);
    for order in 0..valence {
        num = mark_atoms_deriv(
            atoms,
            derivatives,
            i32::from(current.neighbor[order]),
            num,
            flags,
            found,
        )?;
        if num < 0 {
            return Ok(num);
        }
    }
    Ok(num)
}

pub(crate) fn count_one_bond_atoms(
    atoms: &mut [inp_ATOM],
    derivatives: &mut [DERIV_AT],
    start: i32,
    order: i32,
    flags: i8,
    found: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3536 count_one_bond_atoms
    // INCHI✔️❌: int count_one_bond_atoms( inp_ATOM *at,
    // INCHI✔️❌:                           DERIV_AT *da,
    // INCHI✔️❌:                           int start,
    // INCHI✔️❌:                           int ord,
    // INCHI✔️❌:                           char cFlags,
    // INCHI✔️❌:                           int *bFound )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int num = 0;
    // INCHI✔️❌:     if (!( at[at[start].neighbor[ord]].cFlags & cFlags ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[at[start].neighbor[ord]].cFlags |= cFlags;
    // INCHI✔️❌:         num++;
    // INCHI✔️❌:         num = mark_atoms_deriv( at, da, start, num, cFlags, bFound );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: count_one_bond_atoms

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let order = usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let neighbor = usize::from(
        *atoms
            .get(start)
            .and_then(|atom| atom.neighbor.get(order))
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let neighbor_atom = atoms
        .get_mut(neighbor)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if neighbor_atom.cFlags & flags != 0 {
        return Ok(0);
    }
    neighbor_atom.cFlags |= flags;
    mark_atoms_deriv(atoms, derivatives, start as i32, 1, flags, found)
}

#[allow(non_snake_case)]
pub(crate) fn is_nButyl(
    atoms: &[inp_ATOM],
    start: i32,
    previous_order: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3831 is_nButyl
    // INCHI✔️❌: int is_nButyl( inp_ATOM *at, int start, int ord_prev )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, next, curr = start;
    // INCHI✔️❌:     int prev = at[curr].neighbor[ord_prev];
    // INCHI✔️❌:     const int len = 4;
    // INCHI✔️❌:     for (i = 0; i < len; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[curr].el_number != EL_NUMBER_C || at[curr].valence > 2 ||
    // INCHI✔️❌:              at[curr].valence != at[curr].chem_bonds_valence ||
    // INCHI✔️❌:              at[curr].chem_bonds_valence + at[curr].num_H != 4 ||
    // INCHI✔️❌:              at[curr].charge || at[curr].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[curr].valence == 2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             next = at[curr].neighbor[prev == at[curr].neighbor[0]];
    // INCHI✔️❌:             prev = curr;
    // INCHI✔️❌:             curr = next;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return i + 1 == len;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_nButyl

    let mut current = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let previous_order =
        usize::try_from(previous_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut previous = usize::from(
        *atoms
            .get(current)
            .and_then(|atom| atom.neighbor.get(previous_order))
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    for index in 0..4 {
        let atom = atoms
            .get(current)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.el_number != EL_NUMBER_C
            || atom.valence > 2
            || atom.valence != atom.chem_bonds_valence
            || atom.chem_bonds_valence.wrapping_add(atom.num_H) != 4
            || atom.charge != 0
            || atom.radical != 0
        {
            return Ok(0);
        }
        if atom.valence == 2 {
            let next_order = usize::from(previous == usize::from(atom.neighbor[0]));
            let next = usize::from(atom.neighbor[next_order]);
            previous = current;
            current = next;
        } else {
            return Ok(i32::from(index + 1 == 4));
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn is_Me_or_Et(
    atoms: &[inp_ATOM],
    start: i32,
    previous_order: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3862 is_Me_or_Et
    // INCHI✔️❌: int is_Me_or_Et( inp_ATOM *at, int start, int ord_prev )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, neigh = -1;
    // INCHI✔️❌:     if (at[start].el_number != EL_NUMBER_C || at[start].valence > 2 ||
    // INCHI✔️❌:          at[start].valence != at[start].chem_bonds_valence ||
    // INCHI✔️❌:          at[start].chem_bonds_valence + at[start].num_H != 4 ||
    // INCHI✔️❌:          at[start].charge || at[start].radical)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     for (i = 0; i < at[start].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i == ord_prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (neigh >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         neigh = at[start].neighbor[i];
    // INCHI✔️❌:         if (at[neigh].el_number != EL_NUMBER_C || at[neigh].valence > 1 ||
    // INCHI✔️❌:              at[neigh].valence != at[neigh].chem_bonds_valence ||
    // INCHI✔️❌:              at[neigh].chem_bonds_valence + at[neigh].num_H != 4 ||
    // INCHI✔️❌:              at[neigh].charge || at[neigh].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 1 + ( neigh >= 0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_Me_or_Et

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(start)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.el_number != EL_NUMBER_C
        || atom.valence > 2
        || atom.valence != atom.chem_bonds_valence
        || atom.chem_bonds_valence.wrapping_add(atom.num_H) != 4
        || atom.charge != 0
        || atom.radical != 0
    {
        return Ok(0);
    }
    let mut neighbor = None;
    for order in 0..usize::try_from(i32::from(atom.valence)).unwrap_or(0) {
        if order as i32 == previous_order {
            continue;
        }
        if neighbor.is_some() {
            return Ok(0);
        }
        let candidate = atoms
            .get(usize::from(atom.neighbor[order]))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        neighbor = Some(candidate);
        if candidate.el_number != EL_NUMBER_C
            || candidate.valence > 1
            || candidate.valence != candidate.chem_bonds_valence
            || candidate.chem_bonds_valence.wrapping_add(candidate.num_H) != 4
            || candidate.charge != 0
            || candidate.radical != 0
        {
            return Ok(0);
        }
    }
    Ok(1 + i32::from(neighbor.is_some()))
}

#[allow(non_snake_case)]
pub(crate) fn is_phenyl(
    atoms: &[inp_ATOM],
    start: i32,
    previous_order: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4103 is_phenyl
    // INCHI✔️❌: int is_phenyl( inp_ATOM *at, int start, int ord_prev )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int k, neigh, cur_at, ord;
    // INCHI✔️❌:     if (at[start].el_number != EL_NUMBER_C || at[start].valence != 3 ||
    // INCHI✔️❌:          at[start].valence + 1 != at[start].chem_bonds_valence ||
    // INCHI✔️❌:          at[start].chem_bonds_valence + at[start].num_H != 4 ||
    // INCHI✔️❌:          at[start].charge || at[start].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ord = ( ord_prev + 1 ) % 3;
    // INCHI✔️❌:     cur_at = start;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k = 0; k < 5; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         neigh = at[cur_at].neighbor[ord];
    // INCHI✔️❌:         if (at[neigh].el_number != EL_NUMBER_C || at[neigh].valence != 2 ||
    // INCHI✔️❌:              at[neigh].valence + 1 != at[neigh].chem_bonds_valence ||
    // INCHI✔️❌:              at[neigh].chem_bonds_valence + at[neigh].num_H != 4 ||
    // INCHI✔️❌:              at[neigh].charge || at[neigh].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ord = ( at[neigh].neighbor[0] == cur_at );
    // INCHI✔️❌:         cur_at = neigh;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( at[cur_at].neighbor[ord] == start );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_phenyl

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let start_atom = atoms
        .get(start)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if start_atom.el_number != EL_NUMBER_C
        || start_atom.valence != 3
        || i16::from(start_atom.valence) + 1 != i16::from(start_atom.chem_bonds_valence)
        || i16::from(start_atom.chem_bonds_valence) + i16::from(start_atom.num_H) != 4
        || start_atom.charge != 0
        || start_atom.radical != 0
    {
        return Ok(0);
    }
    let direction = (previous_order + 1) % 3;
    let mut order = usize::try_from(direction).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut current = start;
    for _ in 0..5 {
        let current_atom = atoms
            .get(current)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let neighbor = usize::from(
            *current_atom
                .neighbor
                .get(order)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let neighbor_atom = atoms
            .get(neighbor)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor_atom.el_number != EL_NUMBER_C
            || neighbor_atom.valence != 2
            || i16::from(neighbor_atom.valence) + 1 != i16::from(neighbor_atom.chem_bonds_valence)
            || i16::from(neighbor_atom.chem_bonds_valence) + i16::from(neighbor_atom.num_H) != 4
            || neighbor_atom.charge != 0
            || neighbor_atom.radical != 0
        {
            return Ok(0);
        }
        order = usize::from(neighbor_atom.neighbor[0] == current as u16);
        current = neighbor;
    }
    Ok(i32::from(
        atoms
            .get(current)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .neighbor[order]
            == start as u16,
    ))
}

#[allow(clippy::too_many_arguments, non_snake_case)]
pub(crate) fn is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
    atoms: &[inp_ATOM],
    start: i32,
    _num_atoms: i32,
    derivative: &DERIV_AT,
    derivative_index: i32,
    mut underiv: Option<&mut [i8]>,
    len_underiv: i32,
    _underiv2: Option<&mut [i8]>,
    _len_underiv2: i32,
    _bit_underiv: Option<&mut BIT_UNDERIV>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4136 is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR
    // INCHI✔️❌: int is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( inp_ATOM *at,
    // INCHI✔️❌:                                              int start,
    // INCHI✔️❌:                                              int num_atoms,
    // INCHI✔️❌:                                              DERIV_AT *da1,
    // INCHI✔️❌:                                              int idrv,
    // INCHI✔️❌:                                              char *szUnderiv,
    // INCHI✔️❌:                                              int lenUnderiv,
    // INCHI✔️❌:                                              char *szUnderiv2,
    // INCHI✔️❌:                                              int lenUnderiv2,
    // INCHI✔️❌:                                              BIT_UNDERIV *bitUnderiv )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, neigh_at[2], prev_ord[2], neigh, is_B = 0, is_C = 0, n0, n1, n2, n3, n[4] = {0}, nFound, ind1, ind2; /* djb-rwth: adding variables for char -> int conversion of subscripts */
    // INCHI✔️❌:     AT_NUMB *p;
    // INCHI✔️❌:     const char *pUnk;
    // INCHI✔️❌:     char str[16] = { '\0' };
    // INCHI✔️❌: #if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:     char strO[8] = { '\0' };
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if( defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:     char strN[8] = { '\0' };
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (da1->typ[idrv] && ( da1->typ[idrv] & DERIV_RING_OUTSIDE_PRECURSOR ) == da1->typ[idrv] &&
    // INCHI✔️❌:          da1->typ[idrv + 1] && ( da1->typ[idrv + 1] & DERIV_RING_OUTSIDE_PRECURSOR ) == da1->typ[idrv + 1])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( (da1->typ[idrv] & DERIV_RING_O_OUTSIDE_PRECURSOR  || da1->typ[idrv+1] != DERIV_RING_O_OUTSIDE_PRECURSOR) &&
    // INCHI✔️❌:     (da1->typ[idrv] != DERIV_RING_NH_OUTSIDE_PRECURSOR || da1->typ[idrv+1] != DERIV_RING_NH_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (at[start].charge || at[start].radical || at[start].valence != at[start].chem_bonds_valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* check bond types start-n0 and start-n3 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check
    // INCHI✔️❌:     n1    n0
    // INCHI✔️❌:     R1---O
    // INCHI✔️❌:     |     \
    // INCHI✔️❌:     |      B [start]
    // INCHI✔️❌:     |     /
    // INCHI✔️❌:     R2---O
    // INCHI✔️❌:     n2    n3
    // INCHI✔️❌:
    // INCHI✔️❌:     All bond are single except n1-n2 (R1-R2), which may be either single or aromatic
    // INCHI✔️❌:
    // INCHI✔️❌:     */
    // INCHI✔️❌:     nFound = 0;
    // INCHI✔️❌:     ind1 = da1->ord[0] - '0'; /* djb-rwth: converting char to int for subscript use */
    // INCHI✔️❌:     ind2 = da1->ord[1] - '0'; /* djb-rwth: converting char to int for subscript use */
    // INCHI✔️❌:     n0 = at[start].neighbor[ind1];
    // INCHI✔️❌:     n3 = at[start].neighbor[ind2];
    // INCHI✔️❌:     /* search for i, j, k such that at[at[n1]neighbor[i]].neighbor[k]= at[n2]neighbor[j] */
    // INCHI✔️❌:     for (i = 0; i < at[n0].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( n1 = at[n0].neighbor[i] ) == start)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* don't go back */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (BOND_SINGLE != at[n0].bond_type[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* check bond type n0-n1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (j = 0; j < at[n1].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (( n2 = at[n1].neighbor[j] ) == n0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue; /* don't go back */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((p = is_in_the_list( at[n3].neighbor, (AT_NUMB) n2, at[n3].valence ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (( BOND_SINGLE == at[n1].bond_type[j] || BOND_ALTERN == at[n1].bond_type[j] ) && /* check bond type n1-n2 */
    // INCHI✔️❌:                      BOND_SINGLE == at[n3].bond_type[p - at[n3].neighbor] && /* check bond type n3-n2 */
    // INCHI✔️❌:                      !nFound++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     n[0] = n0;
    // INCHI✔️❌:                     n[1] = n1;
    // INCHI✔️❌:                     n[2] = n2;
    // INCHI✔️❌:                     n[3] = n3;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nFound != 1 || at[n[1]].el_number != EL_NUMBER_C || at[n[2]].el_number != EL_NUMBER_C)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* n[1] and n[2] cannot have 3 neighboring heteroatoms */
    // INCHI✔️❌:     for (i = 1; i <= 2; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[n1 = n[i]].valence > 3)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (k = 0, j = 0; j < at[n1].valence; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k += ( at[at[n1].neighbor[j]].el_number != EL_NUMBER_C );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (k >= 3)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     n0 = n[0];
    // INCHI✔️❌:     n3 = n[3];
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL != szUnderiv)
    // INCHI✔️❌:     {
    // INCHI✔️❌: #if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:         if (da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR && da1->typ[idrv + 1] == DERIV_RING_O_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[n0].el_number <= at[n3].el_number)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(strO, at[n0].elname);
    // INCHI✔️❌:                 strcat(strO, at[n3].elname);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(strO, at[n3].elname);
    // INCHI✔️❌:                 strcat(strO, at[n0].elname);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(strO, at[n0].elname);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (da1->typ[idrv + 1] == DERIV_RING_O_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     strcat(strO, at[n3].elname);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if( defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:         if (da1->typ[idrv] == DERIV_RING_NH_OUTSIDE_PRECURSOR && da1->typ[idrv + 1] == DERIV_RING_NH_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (1 == at[n0].num_H && 1 == at[n3].num_H)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(strN, "(NH)2");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (1 == at[n0].num_H || 1 == at[n3].num_H)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     strcat(strN, "(NH)N");
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     strcat(strN, "NN");
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da1->typ[idrv] == DERIV_RING_NH_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(strN, 1 == at[n0].num_H ? "(NH)" : "N");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (da1->typ[idrv + 1] == DERIV_RING_NH_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     strcat(strN, 1 == at[n3].num_H ? "(NH)" : "N");
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         switch (da1->typ[idrv] | da1->typ[idrv + 1])
    // INCHI✔️❌:         {
    // INCHI✔️❌: #if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) && defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:             case ( DERIV_RING_O_OUTSIDE_PRECURSOR | DERIV_RING_NH_OUTSIDE_PRECURSOR ):
    // INCHI✔️❌:                 strcat(str, strN);
    // INCHI✔️❌:                 strcat(str, strO); /* "(NH)O" or "(NH)S" */
    // INCHI✔️❌:                 break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:             case ( DERIV_RING_O_OUTSIDE_PRECURSOR ):
    // INCHI✔️❌:                 strcat(str, strO); /* "OO" or "OS" or "SS" */
    // INCHI✔️❌:                 break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if( defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
    // INCHI✔️❌:             case ( DERIV_RING_NH_OUTSIDE_PRECURSOR ):
    // INCHI✔️❌:                 strcat(str, strN);
    // INCHI✔️❌:                 break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 strcat(str, "???");
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         strcat(str, "-");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     underiv_list_add( szUnderiv, lenUnderiv, str, 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:     /*underiv_list_add( szUnderiv, lenUnderiv, da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR? "OO-" : "(NH)2-", 0 );*/
    // INCHI✔️❌:     if (at[start].el_number == EL_NUMBER_B && at[start].valence == 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         is_B = 1;
    // INCHI✔️❌:         underiv_list_add( szUnderiv, lenUnderiv, "B", 0 );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[start].el_number == EL_NUMBER_C && ( at[start].valence == 3 || at[start].valence == 4 ) &&
    // INCHI✔️❌:              at[start].chem_bonds_valence == at[start].valence &&
    // INCHI✔️❌:              at[start].num_H + at[start].chem_bonds_valence == 4)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             is_C = 1;
    // INCHI✔️❌:             underiv_list_add( szUnderiv, lenUnderiv, at[start].valence == 3 ? "CH" : "C", 0 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* locate bonds connecting >B- or >C< or >C- to the rest of the derivative */
    // INCHI✔️❌:     for (i = k = 0; i < at[start].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i != da1->ord[idrv] && i != da1->ord[idrv + 1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             neigh = at[start].neighbor[i];
    // INCHI✔️❌:             if (k < 2 && ( p = is_in_the_list( at[neigh].neighbor, (AT_NUMB) start, at[neigh].valence ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 neigh_at[k] = neigh;
    // INCHI✔️❌:                 prev_ord[k] = p - at[neigh].neighbor;
    // INCHI✔️❌:                 k++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /* program error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = n2 = 0; /* djb-rwth: ignoring LLVM warning: variables used */
    // INCHI✔️❌:     if (is_B && k == 1 && da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (0 < ( n1 = is_Me_or_Et( at, neigh_at[0], prev_ord[0] ) )
    // INCHI✔️❌: #ifdef UNDERIV_OOB_nButyl
    // INCHI✔️❌:              || 0 < ( n2 = is_nButyl( at, neigh_at[0], prev_ord[0] ) )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:              )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             underiv_list_add( szUnderiv, lenUnderiv, n1 == 1 ? "Me" : n1 == 2 ? "Et" : n2 ? "nButyl" : "???", 0 );
    // INCHI✔️❌:             /* is_B */
    // INCHI✔️❌:             underiv_list_add( szUnderiv2, lenUnderiv2, n1 == 1 ? pszDerivName[DERIV_ID_MeBorate] :
    // INCHI✔️❌:                               n1 == 2 ? pszDerivName[DERIV_ID_EtBorate] :
    // INCHI✔️❌:                               n2 ? pszDerivName[DERIV_ID_BuBorate] : "???", ' ' );
    // INCHI✔️❌:             *bitUnderiv |= n1 == 1 ? DERIV_BIT_MeBorate :
    // INCHI✔️❌:                 n1 == 2 ? DERIV_BIT_EtBorate :
    // INCHI✔️❌:                 n2 ? DERIV_BIT_BuBorate : DERIV_BIT_Unknown;
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (is_C)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (k == 1 && at[start].num_H == 1 && is_phenyl( at, neigh_at[0], prev_ord[0] ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv, "Phe", 0 );
    // INCHI✔️❌:                 if (strN[0])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if ((pUnk = underiv_list_get_last( szUnderiv, ' ' ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         underiv_list_add( szUnderiv2, lenUnderiv2, pUnk, ' ' );
    // INCHI✔️❌:                         *bitUnderiv |= DERIV_BIT_Unknown;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_Benzlidene], ' ' );
    // INCHI✔️❌:                     *bitUnderiv |= DERIV_BIT_Benzlidene;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (k == 2 && 0 < ( n1 = is_Me_or_Et( at, neigh_at[0], prev_ord[0] ) ) &&
    // INCHI✔️❌:                  0 < ( n2 = is_Me_or_Et( at, neigh_at[1], prev_ord[1] ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (n1 != n2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     underiv_list_add( szUnderiv, lenUnderiv, "MeEt", 0 );
    // INCHI✔️❌:                     if ((pUnk = underiv_list_get_last( szUnderiv, ' ' ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         underiv_list_add( szUnderiv2, lenUnderiv2, pUnk, ' ' );
    // INCHI✔️❌:                         *bitUnderiv |= DERIV_BIT_Unknown;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (n1 == 1 || n1 == 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         underiv_list_add( szUnderiv, lenUnderiv, n1 == 1 ? "Me2" : "Et2", 0 );
    // INCHI✔️❌:                         if (strN[0] || n1 != 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if ((pUnk = underiv_list_get_last( szUnderiv, ' ' ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 underiv_list_add( szUnderiv2, lenUnderiv2, pUnk, ' ' );
    // INCHI✔️❌:                                 *bitUnderiv |= DERIV_BIT_Unknown;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[n1 == 1 ? DERIV_ID_Acentonate : DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                             *bitUnderiv |= ( n1 == 1 ? DERIV_BIT_Acentonate : DERIV_BIT_Unknown );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( is_B && k == 1 && is_Me_or_Et( at, neigh_at[0], prev_ord[0]) )
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌:     if ( is_B && k == 1 && is_nButyl( at, neigh_at[0], prev_ord[0]) )
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌:     if ( is_C && k == 1 && at[start].num_H == 1 && is_phenyl( at, neigh_at[0], prev_ord[0]) )
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌:     if ( is_C && k == 2 && is_Me_or_Et( at, neigh_at[0], prev_ord[0]) &&
    // INCHI✔️❌:     is_Me_or_Et( at, neigh_at[1], prev_ord[1]) )
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // END INCHI C FUNCTION: is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR

    let start = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let derivative_index =
        usize::try_from(derivative_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let first_type = i32::from(
        *derivative
            .typ
            .get(derivative_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let second_type = i32::from(
        *derivative
            .typ
            .get(derivative_index + 1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    if first_type == 0
        || first_type & DERIV_RING_OUTSIDE_PRECURSOR as i32 != first_type
        || second_type == 0
        || second_type & DERIV_RING_OUTSIDE_PRECURSOR as i32 != second_type
    {
        return Ok(0);
    }
    let start_atom = atoms
        .get(start)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if start_atom.charge != 0
        || start_atom.radical != 0
        || start_atom.valence != start_atom.chem_bonds_valence
    {
        return Ok(0);
    }

    let first_order = usize::try_from(i32::from(derivative.ord[0]) - i32::from(b'0'))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let second_order = usize::try_from(i32::from(derivative.ord[1]) - i32::from(b'0'))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let n0 = usize::from(
        *start_atom
            .neighbor
            .get(first_order)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let n3 = usize::from(
        *start_atom
            .neighbor
            .get(second_order)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let mut found = 0_i32;
    let mut ring = [0_usize; 4];
    let n0_atom = atoms.get(n0).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let n3_atom = atoms.get(n3).ok_or(SourceHeapError::PointerOutOfBounds)?;
    for first_neighbor_order in 0..usize::try_from(i32::from(n0_atom.valence)).unwrap_or(0) {
        let n1 = usize::from(n0_atom.neighbor[first_neighbor_order]);
        if n1 == start || n0_atom.bond_type[first_neighbor_order] != BOND_SINGLE as u8 {
            continue;
        }
        let n1_atom = atoms.get(n1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        for second_neighbor_order in 0..usize::try_from(i32::from(n1_atom.valence)).unwrap_or(0) {
            let n2 = usize::from(n1_atom.neighbor[second_neighbor_order]);
            if n2 == n0 {
                continue;
            }
            let Some(n3_order) = is_in_the_list(
                Some(&n3_atom.neighbor),
                n2 as u16,
                i32::from(n3_atom.valence),
            )?
            else {
                continue;
            };
            if (n1_atom.bond_type[second_neighbor_order] == BOND_SINGLE as u8
                || n1_atom.bond_type[second_neighbor_order] == BOND_TYPE_ALTERN as u8)
                && n3_atom.bond_type[n3_order] == BOND_SINGLE as u8
            {
                if found == 0 {
                    ring = [n0, n1, n2, n3];
                }
                found = found.wrapping_add(1);
            }
        }
    }
    if found != 1
        || atoms
            .get(ring[1])
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .el_number
            != EL_NUMBER_C
        || atoms
            .get(ring[2])
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .el_number
            != EL_NUMBER_C
    {
        return Ok(0);
    }
    for &ring_index in &ring[1..=2] {
        let ring_atom = atoms
            .get(ring_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if ring_atom.valence > 3 {
            let mut heteroatoms = 0_i32;
            for order in 0..usize::try_from(i32::from(ring_atom.valence)).unwrap_or(0) {
                heteroatoms += i32::from(
                    atoms
                        .get(usize::from(ring_atom.neighbor[order]))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .el_number
                        != EL_NUMBER_C,
                );
            }
            if heteroatoms >= 3 {
                return Ok(0);
            }
        }
    }

    let mut oxygen = String::new();
    let mut nitrogen = String::new();
    if underiv.is_some() {
        let first_atom = atoms.get(n0).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let second_atom = atoms.get(n3).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if first_type == DERIV_RING_O_OUTSIDE_PRECURSOR as i32
            && second_type == DERIV_RING_O_OUTSIDE_PRECURSOR as i32
        {
            let first_name = atom_element_name(first_atom)?;
            let second_name = atom_element_name(second_atom)?;
            if first_atom.el_number <= second_atom.el_number {
                oxygen.push_str(&first_name);
                oxygen.push_str(&second_name);
            } else {
                oxygen.push_str(&second_name);
                oxygen.push_str(&first_name);
            }
        } else if first_type == DERIV_RING_O_OUTSIDE_PRECURSOR as i32 {
            oxygen.push_str(&atom_element_name(first_atom)?);
        } else if second_type == DERIV_RING_O_OUTSIDE_PRECURSOR as i32 {
            oxygen.push_str(&atom_element_name(second_atom)?);
        }
        if first_type == DERIV_RING_NH_OUTSIDE_PRECURSOR as i32
            && second_type == DERIV_RING_NH_OUTSIDE_PRECURSOR as i32
        {
            nitrogen.push_str(match (first_atom.num_H == 1, second_atom.num_H == 1) {
                (true, true) => "(NH)2",
                (true, false) | (false, true) => "(NH)N",
                (false, false) => "NN",
            });
        } else if first_type == DERIV_RING_NH_OUTSIDE_PRECURSOR as i32 {
            nitrogen.push_str(if first_atom.num_H == 1 { "(NH)" } else { "N" });
        } else if second_type == DERIV_RING_NH_OUTSIDE_PRECURSOR as i32 {
            nitrogen.push_str(if second_atom.num_H == 1 { "(NH)" } else { "N" });
        }
    }
    let mut prefix = match first_type | second_type {
        value
            if value
                == (DERIV_RING_O_OUTSIDE_PRECURSOR | DERIV_RING_NH_OUTSIDE_PRECURSOR) as i32 =>
        {
            format!("{nitrogen}{oxygen}")
        }
        value if value == DERIV_RING_O_OUTSIDE_PRECURSOR as i32 => oxygen,
        value if value == DERIV_RING_NH_OUTSIDE_PRECURSOR as i32 => nitrogen,
        _ => "???".to_owned(),
    };
    prefix.push('-');
    underiv_add_text(&mut underiv, len_underiv, &prefix, 0)?;

    if start_atom.el_number == EL_NUMBER_B && start_atom.valence == 3 {
        underiv_add_text(&mut underiv, len_underiv, "B", 0)?;
    } else if start_atom.el_number == EL_NUMBER_C
        && (start_atom.valence == 3 || start_atom.valence == 4)
        && start_atom.chem_bonds_valence == start_atom.valence
        && start_atom.num_H.wrapping_add(start_atom.chem_bonds_valence) == 4
    {
        underiv_add_text(
            &mut underiv,
            len_underiv,
            if start_atom.valence == 3 { "CH" } else { "C" },
            0,
        )?;
    } else {
        return Ok(0);
    }

    let mut connected = 0_usize;
    for order in 0..usize::try_from(i32::from(start_atom.valence)).unwrap_or(0) {
        if order as i8 != derivative.ord[derivative_index]
            && order as i8 != derivative.ord[derivative_index + 1]
        {
            let neighbor = atoms
                .get(usize::from(start_atom.neighbor[order]))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let back = is_in_the_list(
                Some(&neighbor.neighbor),
                start as u16,
                i32::from(neighbor.valence),
            )?;
            if connected < 2 && back.is_some() {
                connected += 1;
            } else {
                return Ok(-1);
            }
        }
    }

    // Numeric ord values enter C undefined behavior at the subtraction above; ASCII ord values
    // make every start bond enter this loop and therefore cannot reach the remaining success cases.
    Ok(0)
}

pub(crate) fn add_to_da(derivative: &mut DERIV_AT, added: &DERIV_AT) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3245 add_to_da
    // INCHI✔️❌: int add_to_da( DERIV_AT *da, DERIV_AT *add )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* if add has more than 1 element the elements are in ascending add.ord[] order */
    // INCHI✔️❌:     int i, j, len_da, len_add, numAddHiPri, numDaHiPri;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (len_da = 0, numDaHiPri = 0; len_da < DERIV_AT_LEN && da->typ[len_da]; len_da++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         numDaHiPri += ( 0 != ( da->typ[len_da] & DERIV_UNEXPADABLE ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (len_add = 0, numAddHiPri = 0; len_add < DERIV_AT_LEN && da->typ[len_add]; len_add++) /* djb-rwth: addressing coverity ID #499516 -- definitely not a copy-paste error */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         numAddHiPri += ( 0 != ( add->typ[len_add] & DERIV_UNEXPADABLE ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* HiPri replaces non-HiPri derivatives */
    // INCHI✔️❌:     if (numAddHiPri && !numDaHiPri)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* no harm if already len_da=0 */
    // INCHI✔️❌:         memset( da, 0, sizeof( *da ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         len_da = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* non-HiPri derivatives cannot be added to HiPri */
    // INCHI✔️❌:         if (!numAddHiPri && numDaHiPri)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0; j < DERIV_AT_LEN && add->typ[j]; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < len_da; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (add->ord[j] == da->ord[i])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (add->typ[j] != da->typ[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* error, should not happen */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (add->num[j] != da->num[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -2; /* error, should not happen */
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:                 if ((( len_da>1 || j ) && ( add->other_atom || da->other_atom )) || (1 == len_da && add->other_atom != da->other_atom)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -3; /* other_atom implies single bond to cut */
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (i == len_da)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* add->ord[j] has different ord values from all da->ord[]; add or replace */
    // INCHI✔️❌:             if (len_da < DERIV_AT_LEN)
    // INCHI✔️❌:             {
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:                 if (( i || j ) && add->other_atom)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -3; /* other_atom implies single bond to cut */
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 da->ord[i] = add->ord[j];
    // INCHI✔️❌:                 da->typ[i] = add->typ[j];
    // INCHI✔️❌:                 da->num[i] = add->num[j];
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:                 da->other_atom = add->other_atom;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 len_da++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -4; /* overflow, should not happen */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /****************************************************************************
    // END INCHI C FUNCTION: add_to_da

    let mut derivative_length = 0_usize;
    let mut derivative_high_priority = 0_i32;
    while derivative_length < derivative.typ.len() && derivative.typ[derivative_length] != 0 {
        derivative_high_priority +=
            i32::from(i32::from(derivative.typ[derivative_length]) & DERIV_UNEXPADABLE != 0);
        derivative_length += 1;
    }

    let mut added_length = 0_usize;
    let mut added_high_priority = 0_i32;
    while added_length < derivative.typ.len() && derivative.typ[added_length] != 0 {
        added_high_priority +=
            i32::from(i32::from(added.typ[added_length]) & DERIV_UNEXPADABLE != 0);
        added_length += 1;
    }

    if added_high_priority != 0 && derivative_high_priority == 0 {
        *derivative = DERIV_AT::default();
        derivative_length = 0;
    } else if added_high_priority == 0 && derivative_high_priority != 0 {
        return 0;
    }

    for added_index in 0..added.typ.len() {
        if added.typ[added_index] == 0 {
            break;
        }
        let mut index = 0_usize;
        while index < derivative_length {
            if added.ord[added_index] == derivative.ord[index] {
                if added.typ[added_index] != derivative.typ[index] {
                    return -1;
                }
                if added.num[added_index] != derivative.num[index] {
                    return -2;
                }
                if ((derivative_length > 1 || added_index != 0)
                    && (added.other_atom != 0 || derivative.other_atom != 0))
                    || (derivative_length == 1 && added.other_atom != derivative.other_atom)
                {
                    return -3;
                }
                break;
            }
            index += 1;
        }
        if index == derivative_length {
            if derivative_length >= derivative.typ.len() {
                return -4;
            }
            if (index != 0 || added_index != 0) && added.other_atom != 0 {
                return -3;
            }
            derivative.ord[index] = added.ord[added_index];
            derivative.typ[index] = added.typ[added_index];
            derivative.num[index] = added.num[added_index];
            derivative.other_atom = added.other_atom;
            derivative_length += 1;
        }
    }
    0
}

fn classify_ring_back_connector(
    atoms: &[inp_ATOM],
    center_index: usize,
    connector_index: usize,
) -> Result<Option<(i16, u16)>, SourceHeapError> {
    let connector = atoms
        .get(connector_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let outside_order = usize::from(connector.neighbor[0] as usize == center_index);
    let outside = atoms
        .get(usize::from(connector.neighbor[outside_order]))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let common = connector.valence == 2
        && connector.chem_bonds_valence == connector.valence
        && connector.nNumAtInRingSystem >= 5
        && outside.el_number == EL_NUMBER_C
        && connector.charge == 0
        && connector.radical == 0;
    if common
        && (connector.el_number == EL_NUMBER_O || connector.el_number == EL_NUMBER_S)
        && connector.num_H == 0
    {
        Ok(Some((
            crate::source_types::local_ichinorm::DERIV_RING_O_OUTSIDE_PRECURSOR as i16,
            connector.nBlockSystem,
        )))
    } else if common && connector.el_number == EL_NUMBER_N && connector.num_H == 1 {
        Ok(Some((
            crate::source_types::local_ichinorm::DERIV_RING_NH_OUTSIDE_PRECURSOR as i16,
            connector.nBlockSystem,
        )))
    } else {
        Ok(None)
    }
}

#[allow(clippy::too_many_lines)]
pub(crate) fn get_traversed_deriv_type(
    atoms: &mut [inp_ATOM],
    derivatives: &mut [DERIV_AT],
    atom_index: i32,
    output: &mut DERIV_AT,
    flags: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2677 get_traversed_deriv_type
    // INCHI✔️❌: int get_traversed_deriv_type( inp_ATOM *at,
    // INCHI✔️❌:                               DERIV_AT *da,
    // INCHI✔️❌:                               int k,
    // INCHI✔️❌:                               DERIV_AT *da1,
    // INCHI✔️❌:                               char cFlags )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* at[k] is attachment point of the precursor */
    // INCHI✔️❌:     /* at[(int)at[k].neighbor[m]] is inside precursor */
    // INCHI✔️❌:     /* at[(int)at[k].neighbor[!m]] is inside derivatizing agent */
    // INCHI✔️❌:     /* !!! Except DERIV_RING_O_OUTSIDE_PRECURSOR, DERIV_RING_NH_OUTSIDE_PRECURSOR !!! */
    // INCHI✔️❌:     /* when at[k] is B or C attached to two atoms of the precursor */
    // INCHI✔️❌:     int i, j, m, n1, nBlockSystemFrom, nOrdBack1, nOrdBack2, nOrdBack3, nBackType1, nBackType2;
    // INCHI✔️❌:
    // INCHI✔️❌: #if( defined(DERIV_X_OXIME) || defined(DERIV_RO_COX) || defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:     int n0, n2, n3;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( da1, 0, sizeof( *da1 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     if (at[k].cFlags & cFlags)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (m = 0; m < at[k].valence && !( at[(int) at[k].neighbor[m]].cFlags & cFlags ); m++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (m == at[k].valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;  /* error: at least one neighbor must have cFlags */
    // INCHI✔️❌:                     /* traversing at[k] from at[(int)at[k].neighbor[m]] */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (at[k].valence == 1 && at[k].num_H && (
    // INCHI✔️❌:         at[k].el_number == EL_NUMBER_O ||
    // INCHI✔️❌:         at[k].el_number == EL_NUMBER_N ||
    // INCHI✔️❌:         at[k].el_number == EL_NUMBER_S ||
    // INCHI✔️❌:         at[k].el_number == EL_NUMBER_P ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return DERIV_NOT;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (is_el_a_metal( at[k].el_number ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return DERIV_NOT;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #ifdef NEVER
    // INCHI✔️❌:     if (at[k].el_number == EL_NUMBER_N && at[k].valence >= 2 && at[k].chem_bonds_valence <= 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return DERIV_NOT; /* prohibit -N-, -N=, allow -N# as in isocyano -N#C or NO2 */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     /* m is the ord of the bond from which at[k] was reached first time */
    // INCHI✔️❌:     if (da[k].typ[0] && ( da[k].typ[0] & DERIV_UNEXPADABLE ) == da[k].typ[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #ifdef DERIV_X_OXIME
    // INCHI✔️❌:     if (at[k].nNumAtInRingSystem == 1 && at[k].el_number == EL_NUMBER_N &&
    // INCHI✔️❌:          at[k].valence == 2 && at[k].chem_bonds_valence == 3 &&
    // INCHI✔️❌:          !at[k].num_H && !at[k].charge && !at[k].radical &&
    // INCHI✔️❌:          at[n0 = at[k].neighbor[m]].el_number == EL_NUMBER_C &&  /* inside precursor */
    // INCHI✔️❌:          at[n1 = at[k].neighbor[!m]].el_number == EL_NUMBER_O && /* inside derivatizing agent */
    // INCHI✔️❌:          at[k].bond_type[m] == BOND_DOUBLE && at[k].bond_type[!m] == BOND_SINGLE &&
    // INCHI✔️❌:          at[n0].valence == 3 - at[n0].num_H && at[n0].chem_bonds_valence == 4 - at[n0].num_H && !at[n0].charge && !at[n0].radical && /* C */
    // INCHI✔️❌:          at[n1].valence == 2 && at[n1].chem_bonds_valence == 2 && !at[n1].charge && !at[n1].radical                   /* O */
    // INCHI✔️❌:          )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* found C==N--O */
    // INCHI✔️❌:         /* traversing from C to O; C(at[neighbor[m]]) has cFlag; N is at[k]; N-O bond is to be broken */
    // INCHI✔️❌:         /* m     !m
    // INCHI✔️❌:            n0 k  n1 n2       n2    n2  n3   n2  n3      n2  n3...
    // INCHI✔️❌:            C==N--O--R; -R: -CH3, -CH2-CH3, -Si(CH3)3, -CH2-C6H5; cut N-O and replace =N- with =O 2013-08-22 DT */
    // INCHI✔️❌:         /* check other neighbors of C: they should be C,H or C,C */
    // INCHI✔️❌:         for (i = 0; i < at[n0].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[n0].neighbor[i] != k && at[at[n0].neighbor[i]].el_number != EL_NUMBER_C)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* found    C
    // INCHI✔️❌:         |
    // INCHI✔️❌:         C==N==O
    // INCHI✔️❌:         |
    // INCHI✔️❌:         C
    // INCHI✔️❌:         */
    // INCHI✔️❌:         /* find other neighbor of O */
    // INCHI✔️❌:         n2 = at[n1].neighbor[at[n1].neighbor[0] == k];
    // INCHI✔️❌:         if (is_Si_IV( at, n2 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int n4 = -1;
    // INCHI✔️❌:             for (i = 0; i < at[n2].valence; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int n3 = at[n2].neighbor[i];
    // INCHI✔️❌:                 if (n3 == n1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue; /* atom O */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[n3].el_number != EL_NUMBER_C || at[n3].charge ||
    // INCHI✔️❌:                      at[n3].radical || at[n3].chem_bonds_valence != at[n3].valence)
    // INCHI✔️❌:                     goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌:                 if (n4 == -1 && at[n3].valence == 4 && !at[n3].num_H)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     n4 = n3; /* possibly tret-butyl */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[n3].chem_bonds_valence != 1 || at[n3].num_H != 3)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌:                                                     /* methyl identified */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (n4 == -1)
    // INCHI✔️❌:             {
    // INCHI✔️❌: #ifdef UNDERIV_X_OXIME_TMS
    // INCHI✔️❌:                 /* found    C
    // INCHI✔️❌:                 |        n2
    // INCHI✔️❌:                 C==N==O--Si(CH3)3
    // INCHI✔️❌:                 |
    // INCHI✔️❌:                 C
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 da1->ord[0] = !m;         /* ord of neighbor O, already checked */
    // INCHI✔️❌:                 da1->typ[0] = DERIV_X_OXIME;   /* type */
    // INCHI✔️❌:                 da1->num[0] = 5; /* 5 atoms: -O-Si(CH3)3 -O-TMS*/
    // INCHI✔️❌:                 return DERIV_X_OXIME;   /* >C=N-O-Si(CH3)3 */
    // INCHI✔️❌: #else
    // INCHI✔️❌:                 goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             }
    // INCHI✔️❌: #ifndef UNDERIV_X_OXIME_TBDMS
    // INCHI✔️❌:             goto check_next_derivative; /* do not include TBDMS */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             for (i = 0; i < at[n4].valence; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int n3 = at[n4].neighbor[i];
    // INCHI✔️❌:                 if (n3 == n2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue; /* atom Si */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[n3].el_number != EL_NUMBER_C || at[n3].charge ||
    // INCHI✔️❌:                      at[n3].radical || at[n3].chem_bonds_valence != at[n3].valence)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[n3].chem_bonds_valence != 1 || at[n3].num_H != 3)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌:                                                 /* methyl identified */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             da1->ord[0] = !m;         /* ord of neighbor O, already checked */
    // INCHI✔️❌:             da1->typ[0] = DERIV_X_OXIME;   /* type */
    // INCHI✔️❌:             da1->num[0] = 8; /* 8 atoms: -O-Si(CH3)2-C(CH3)3 -O-TBDMS  */
    // INCHI✔️❌:             return DERIV_X_OXIME;   /* >C=N-O-Si(CH3)2-C(CH3)3 */
    // INCHI✔️❌:
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[n2].el_number != EL_NUMBER_C)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[n2].chem_bonds_valence == 1 && at[n2].num_H == 3 && !at[n2].charge && !at[n2].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             da1->ord[0] = !m;         /* ord of neighbor O, already checked */
    // INCHI✔️❌:             da1->typ[0] = DERIV_X_OXIME;   /* type */
    // INCHI✔️❌:             da1->num[0] = 2; /* 2 atoms: -O-CH3 */
    // INCHI✔️❌:             return DERIV_X_OXIME;   /* >C=N-O-CH3 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[n2].valence != 2 || at[n2].chem_bonds_valence != 2 || at[n2].num_H != 2 || at[n2].charge || at[n2].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto check_next_derivative; /* wrong neighbor */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         n3 = at[n2].neighbor[at[n2].neighbor[0] == n1];
    // INCHI✔️❌:         if (at[n3].chem_bonds_valence == 1 && at[n3].num_H == 3 && !at[n3].charge && !at[n3].radical)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             da1->ord[0] = !m;         /* ord of neighbor O, already checked */
    // INCHI✔️❌:             da1->typ[0] = DERIV_X_OXIME;   /* type */
    // INCHI✔️❌:             da1->num[0] = 3; /* 3 atoms: -O-CH2-CH3 */
    // INCHI✔️❌:             return DERIV_X_OXIME;   /* >C=N-O-CH2-CH3 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[n3].valence == 3 && at[n3].chem_bonds_valence == 4 && !at[n3].num_H && !at[n3].charge && !at[n3].radical &&
    // INCHI✔️❌:              at[n3].nRingSystem != at[n2].nRingSystem && at[n3].bCutVertex && at[n3].nNumAtInRingSystem == 6)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (is_Phenyl( at, n2, n3 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1->ord[0] = !m;         /* ord of neighbor O, already checked */
    // INCHI✔️❌:                 da1->typ[0] = DERIV_X_OXIME;   /* type */
    // INCHI✔️❌:                 da1->num[0] = 8; /* 8 atoms: O--CH2-C6H5 */
    // INCHI✔️❌:                 return DERIV_X_OXIME;   /* >C=N-O-CH2-C6H5 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: check_next_derivative:
    // INCHI✔️❌: #endif  /* DERIV_X_OXIME */
    // INCHI✔️❌: #ifdef DERIV_DANSYL
    // INCHI✔️❌:     if (at[k].nNumAtInRingSystem == 1 &&
    // INCHI✔️❌:         ( (( at[k].el_number == EL_NUMBER_O || at[k].el_number == EL_NUMBER_S ) && at[k].valence == 2) ||
    // INCHI✔️❌:           (at[k].el_number == EL_NUMBER_N && at[k].valence == 2 && at[k].num_H == 1) || (at[k].valence == 3 && at[k].num_H == 2) )) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         DERIV_AT da2;
    // INCHI✔️❌:         memset( &da2, 0, sizeof( da2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         for (j = 0, n1 = 0; j < at[k].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (j == m)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (at[i = at[k].neighbor[j]].el_number == EL_NUMBER_S &&
    // INCHI✔️❌:                  at[i].valence == 4 && at[i].chem_bonds_valence == 6 &&
    // INCHI✔️❌:                  is_Dansyl( at, k, j, da, &da2 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 n1++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (n1 == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *da1 = da2;
    // INCHI✔️❌:             return DERIV_DANSYL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[k].nNumAtInRingSystem == 1 && ( at[k].el_number == EL_NUMBER_O || at[k].el_number == EL_NUMBER_S ) &&
    // INCHI✔️❌:          at[k].valence == 2 && at[k].chem_bonds_valence == 2 &&
    // INCHI✔️❌:          !at[k].num_H && !at[k].charge && !at[k].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*   at[k].neighbor[m] k n1==at[k].neighbor[!m]  */
    // INCHI✔️❌:         /*               -> A--O--B -> traversing from A to B; cut O--B */
    // INCHI✔️❌:         /* check for carboxy A(=O)-O-B and A--O--B(=O) */
    // INCHI✔️❌:         /*int has_A_CO   = is_C_or_S_DB_O( at, at[k].neighbor[m] );*/
    // INCHI✔️❌:         int has_B_CO = is_C_or_S_DB_O( at, n1 = at[k].neighbor[!m] );/* B is C(=o) or S(=O) */
    // INCHI✔️❌:         int is_A_Si_IV = is_Si_IV( at, at[k].neighbor[m] ); /* A is >Si< */
    // INCHI✔️❌:                                                             /* int is_B_Si_IV = is_Si_IV( at, at[k].neighbor[!m] );*/
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef DERIV_RO_COX
    // INCHI✔️❌:                                                             /*               n3
    // INCHI✔️❌:                                                             at[k].neighbor[m]  k n1 n2
    // INCHI✔️❌:                                                             R--O--C--X; -X = -CH3,  -Phenyl,  -C[n]F[2n+1] 0 < n < 4
    // INCHI✔️❌:                                                             ||        (acetate)(benzoate)
    // INCHI✔️❌:                                                             O
    // INCHI✔️❌:                                                             */
    // INCHI✔️❌:         n3 = at[k].neighbor[m];  /* R */
    // INCHI✔️❌:         if (has_B_CO && is_C_DB_O( at, n1 ) &&  /* B:n1 is >C=O  */
    // INCHI✔️❌:              !is_silyl2( at, n3, k ) &&
    // INCHI✔️❌:              !is_el_a_metal( at[n3].el_number ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = 0; j < at[n1].valence; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[n1].neighbor[j] != k && at[n1].bond_type[j] == BOND_SINGLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* the only suspected neighbor */
    // INCHI✔️❌:                     n2 = at[n1].neighbor[j]; /* X */
    // INCHI✔️❌:                     n0 = is_CF3_or_linC3F7a(at, n2, n1); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     if (n0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n0 = 2 + 3 * n0 + 1; /* (6,9,12 atoms) -C(=O)C[n]F[2n+1]; is_CF3_or_linC3F7a returns n */
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #ifdef UNDERIV_RO_COX_Me
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                         if (is_Methyl( at, n2 ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* methyl */
    // INCHI✔️❌:                             if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 n0 = 3; /* 3 atoms: -C(=O)-CH3 */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif /* UNDERIV_RO_COX_Me */
    // INCHI✔️❌: #ifdef UNDERIV_RO_COX_Et
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                             if (is_Ethyl( at, n1, n2 ))
    // INCHI✔️❌:                             { /* ethyl */
    // INCHI✔️❌:                                 if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     n0 = 4; /* 4 atoms: -C(=O)-CH2-CH3 */
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌: #endif /* UNDERIV_RO_COX_Et */
    // INCHI✔️❌: #ifdef UNDERIV_RO_COX_BENZOATES
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 if (is_Phenyl( at, n1, n2 ))
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         n0 = 2 + 6; /* 8  atoms -C(=O)-C6H5 */
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌: #endif /* UNDERIV_RO_COX_BENZOATES */
    // INCHI✔️❌: #ifdef UNDERIV_RO_COX_PENTAFLOUROBENZOATES
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (is_PentaFluoroPhenyl( at, n1, n2 ))
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
    // INCHI✔️❌:                                         {
    // INCHI✔️❌:                                             n0 = 13; /* 13  atoms -C(=O)-C6F5 */
    // INCHI✔️❌:                                         }
    // INCHI✔️❌:                                     }
    // INCHI✔️❌: #endif  /* UNDERIV_RO_COX_PENTAFLOUROBENZOATES */
    // INCHI✔️❌:                     if (n0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         da1->ord[0] = !m;         /* ord of at[k]'s, that is, -O-'s neighbor X in C(=O)-X */
    // INCHI✔️❌:                         da1->typ[0] = DERIV_RO_COX;   /* type */
    // INCHI✔️❌:                         da1->num[0] = n0; /* num atoms in derivatizing attachement */
    // INCHI✔️❌:                         return DERIV_RO_COX;   /* R--O--C(=O)--X; -X = -CH3, -C[n]F[2n+1] 0 < n < 4, -Phenyl, -C6F5 */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif   /* DERIV_RO_COX */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (is_A_Si_IV && has_B_CO)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*                                             precursor | deriv.agent */
    // INCHI✔️❌:             ; /* do not cut bond --- in A=>Si<, B(=O), B=C,S: Si(IV)-O---B(=O) */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* at[k] is precursor's attachment point;
    // INCHI✔️❌:             at[at[k].neighbor[!m]] belongs to drivatizing agent,
    // INCHI✔️❌:             at[at[k].neighbor[m]]  was marked (from_atom) */
    // INCHI✔️❌:             if (is_possibly_deriv_neigh( at, k, !m, DERIV_BRIDGE_O, cFlags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1->ord[0] = !m;         /* ord of neighbor B, not traversed yet */
    // INCHI✔️❌:                 da1->typ[0] = DERIV_BRIDGE_O;   /* type */
    // INCHI✔️❌:                 return DERIV_BRIDGE_O;   /* Representative: R-C(=O)-O(!m)--[D]  */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     R--C==N   Me or Et
    // INCHI✔️❌:     |   \ /
    // INCHI✔️❌:     |    C
    // INCHI✔️❌:     |   / \
    // INCHI✔️❌:     O--CH2 Me or ET
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (at[k].el_number == EL_NUMBER_O && at[k].nNumAtInRingSystem == 5 &&
    // INCHI✔️❌:          is_DERIV_RING_DMOX_DEOX_O( at, k, m, da, da1 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return DERIV_RING_DMOX_DEOX_O;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (at[k].el_number == EL_NUMBER_N && at[k].nNumAtInRingSystem == 5 &&
    // INCHI✔️❌:          at[k].valence == 2 && at[k].chem_bonds_valence == 3 &&
    // INCHI✔️❌:          is_DERIV_RING_DMOX_DEOX_N( at, k, m, da, da1 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return DERIV_RING_DMOX_DEOX_N;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* R1--NH--R2 */
    // INCHI✔️❌:     if (at[k].bCutVertex && at[k].el_number == EL_NUMBER_N &&
    // INCHI✔️❌:          at[k].valence == 2 && at[k].chem_bonds_valence == at[k].valence &&
    // INCHI✔️❌:          at[k].valence + at[k].num_H == 3 && !at[k].charge && !at[k].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         da1->ord[0] = !m;         /* ord of neighbor B, not traversed yet */
    // INCHI✔️❌:         da1->typ[0] = DERIV_BRIDGE_NH;   /* type */
    // INCHI✔️❌:         return DERIV_BRIDGE_NH;   /* R1-NH-R2  amine */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     R2
    // INCHI✔️❌:     /
    // INCHI✔️❌:     R1----N
    // INCHI✔️❌:     \
    // INCHI✔️❌:     R3
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[k].bCutVertex && at[k].el_number == EL_NUMBER_N &&
    // INCHI✔️❌:          at[k].valence == 3 && at[k].chem_bonds_valence == at[k].valence &&
    // INCHI✔️❌:          at[k].valence + at[k].num_H == 3 && !at[k].charge && !at[k].radical)
    // INCHI✔️❌:     {
    // INCHI✔️❌: #if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) || defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )
    // INCHI✔️❌:         if (at[j = at[k].neighbor[m]].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:              at[j].valence == 3 && at[j].chem_bonds_valence == 4 &&
    // INCHI✔️❌:              ( i = is_DERIV_RING2_PRRLDD_PPRDN( at, k, m, da, da1 ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return i;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int rm1 = ( at[at[k].neighbor[m]].nRingSystem == at[at[k].neighbor[( m + 1 ) % 3]].nRingSystem );
    // INCHI✔️❌:             int rm2 = ( at[at[k].neighbor[m]].nRingSystem == at[at[k].neighbor[( m + 2 ) % 3]].nRingSystem );
    // INCHI✔️❌:             int r12 = ( at[at[k].neighbor[( m + 1 ) % 3]].nRingSystem == at[at[k].neighbor[( m + 2 ) % 3]].nRingSystem );
    // INCHI✔️❌:             int ord[2] = { -1, -1 };
    // INCHI✔️❌:             i = 0; /* type: tertriary amine: DERIV_AMINE_tN */
    // INCHI✔️❌:             switch (rm1 + rm2 + r12)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case 0:
    // INCHI✔️❌:                     /* -N< has no ring bonds */
    // INCHI✔️❌:                     if (is_possibly_deriv_neigh( at, k, ( m + 1 ) % 3, DERIV_AMINE_tN, cFlags ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ord[i++] = ( m + 1 ) % 3; /* ord of a non-ring neighbor, not traversed yet */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (is_possibly_deriv_neigh( at, k, ( m + 2 ) % 3, DERIV_AMINE_tN, cFlags ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ord[i++] = ( m + 2 ) % 3; /* ord of another non-ring neighbor, not traversed yet */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (i == 2 && ord[0] > ord[1])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int tmp = ord[0];
    // INCHI✔️❌:                         ord[0] = ord[1];
    // INCHI✔️❌:                         ord[1] = tmp;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:
    // INCHI✔️❌:                 case 1:
    // INCHI✔️❌:                     /* -N< has one non-ring bond; do not consider [m] because it is "from" bond */
    // INCHI✔️❌:                     if (rm1 && is_possibly_deriv_neigh( at, k, ( m + 2 ) % 3, DERIV_AMINE_tN, cFlags ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ord[i++] = ( m + 2 ) % 3;   /* ord of a single non-ring neighbor, not traversed yet */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                         if (rm2 && is_possibly_deriv_neigh( at, k, ( m + 1 ) % 3, DERIV_AMINE_tN, cFlags ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ord[i++] = ( m + 1 ) % 3; /* ord of a single non-ring neighbor, not traversed yet */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     /* r12 != 0 <=> at[k]neighbor[m] is the only non-ring bond; ignore it because it is "from" bond */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (j = 0; j < i; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1->ord[j] = ord[j];
    // INCHI✔️❌:                 da1->typ[j] = DERIV_AMINE_tN;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return DERIV_AMINE_tN;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return 0; /* all neighbors or two untraversed edges are in one ring system */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*-----------------------------------------------------------------*/
    // INCHI✔️❌:     /* DERIV_RING_O_OUTSIDE_PRECURSOR, DERIV_RING_NH_OUTSIDE_PRECURSOR */
    // INCHI✔️❌:     if (at[k].bCutVertex && /* DD */
    // INCHI✔️❌:          at[k].valence == at[k].chem_bonds_valence &&
    // INCHI✔️❌:          ( !at[k].num_H || (at[k].el_number == EL_NUMBER_C && 1 == at[k].num_H) ) &&
    // INCHI✔️❌:          !at[k].charge && !at[k].radical &&
    // INCHI✔️❌:          ( (at[k].el_number == EL_NUMBER_C  && at[k].valence + at[k].num_H == 4) ||
    // INCHI✔️❌:            (at[k].el_number == EL_NUMBER_SI && at[k].valence == 4) ||
    // INCHI✔️❌:            (at[k].el_number == EL_NUMBER_B  && at[k].valence == 3) )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /*-->    j \        entering path: ->X--O--DD
    // INCHI✔️❌:         --X--O  v
    // INCHI✔️❌:         |    \ k /    DD = C, CH, Si, B
    // INCHI✔️❌:         |     DD
    // INCHI✔️❌:         |    /   \     O = O, S, NH  = at[j], going from DD
    // INCHI✔️❌:         --Y--O
    // INCHI✔️❌:         X, Y -- must be C
    // INCHI✔️❌:         */
    // INCHI✔️❌:         nBlockSystemFrom = 0;
    // INCHI✔️❌:         nBackType1 = nBackType2 = 0;
    // INCHI✔️❌:         nOrdBack1 = nOrdBack2 = nOrdBack3 = -1;
    // INCHI✔️❌:         j = (int) at[k].neighbor[m];
    // INCHI✔️❌:         ; /* X */
    // INCHI✔️❌:         if (( at[j].el_number == EL_NUMBER_O || at[j].el_number == EL_NUMBER_S ) && at[j].valence == 2 &&
    // INCHI✔️❌:              at[j].chem_bonds_valence == at[j].valence &&
    // INCHI✔️❌:              at[j].nNumAtInRingSystem >= 5 &&
    // INCHI✔️❌:              at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* X is C */
    // INCHI✔️❌:              !at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nBackType1 = DERIV_RING_O_OUTSIDE_PRECURSOR;
    // INCHI✔️❌:             nBlockSystemFrom = at[j].nBlockSystem;
    // INCHI✔️❌:             nOrdBack1 = m; /* predecessor */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:             if (at[j].el_number == EL_NUMBER_N && at[j].valence == 2 &&
    // INCHI✔️❌:                  at[j].chem_bonds_valence == at[j].valence &&
    // INCHI✔️❌:                  at[j].nNumAtInRingSystem >= 5 &&
    // INCHI✔️❌:                  at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* X is C */
    // INCHI✔️❌:                  1 == at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nBackType1 = DERIV_RING_NH_OUTSIDE_PRECURSOR;
    // INCHI✔️❌:                 nBlockSystemFrom = at[j].nBlockSystem;
    // INCHI✔️❌:                 nOrdBack1 = m; /* predecessor */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         if (nBlockSystemFrom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int num1, num2, bFound;
    // INCHI✔️❌:             at[k].cFlags |= CFLAG_MARK_BLOCK;
    // INCHI✔️❌:             /* mark precursor atoms + at[k] */
    // INCHI✔️❌:             num1 = mark_atoms_cFlags( at, at[k].neighbor[nOrdBack1], 1, CFLAG_MARK_BLOCK );
    // INCHI✔️❌:             for (i = 0; i < at[k].valence; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i == nOrdBack1)
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 j = (int) at[k].neighbor[i];
    // INCHI✔️❌:                 bFound = 0;
    // INCHI✔️❌:                 if (at[j].cFlags & CFLAG_MARK_BLOCK)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (( at[j].el_number == EL_NUMBER_O || at[j].el_number == EL_NUMBER_S ) && at[j].valence == 2 &&
    // INCHI✔️❌:                          at[j].chem_bonds_valence == at[j].valence &&
    // INCHI✔️❌:                          at[j].nNumAtInRingSystem >= 5 &&
    // INCHI✔️❌:                          at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* Y is C */
    // INCHI✔️❌:                          !at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bFound = 1;
    // INCHI✔️❌:                         if (nOrdBack2 < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nOrdBack2 = i; /* predecessor #2 */
    // INCHI✔️❌:                             nBackType2 = DERIV_RING_O_OUTSIDE_PRECURSOR;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nOrdBack3 = i; /* predecessor #3 -- should not happen */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (at[j].el_number == EL_NUMBER_N && at[j].valence == 2 &&
    // INCHI✔️❌:                          at[j].chem_bonds_valence == at[j].valence &&
    // INCHI✔️❌:                          at[j].nNumAtInRingSystem >= 5 &&
    // INCHI✔️❌:                          at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* Y is C */
    // INCHI✔️❌:                          1 == at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bFound = 1;
    // INCHI✔️❌:                         if (nOrdBack2 < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nOrdBack2 = i; /* predecessor #2 */
    // INCHI✔️❌:                             nBackType2 = DERIV_RING_NH_OUTSIDE_PRECURSOR;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nOrdBack3 = i; /* predecessor #3 -- should not happen */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (!bFound)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nOrdBack3 = 99; /* reject: wrong neighboring atom in the same block */
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             num2 = unmark_atoms_cFlags( at, k, 0, CFLAG_MARK_BLOCK, CFLAG_MARK_BLOCK_INV );
    // INCHI✔️❌:             if (num1 != num2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1; /* mark_atoms_cFlags() program error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nOrdBack2 >= 0 && nOrdBack3 < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* note: da1 refers to DD, which is a neighbor of 2 precursor atoms; ord point to precursor attachment points, O */
    // INCHI✔️❌:                 if (nOrdBack1 < nOrdBack2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da1->ord[0] = nOrdBack1;  /* ord of a ring neighbor, traversed */
    // INCHI✔️❌:                     da1->typ[0] = nBackType1;
    // INCHI✔️❌:                     da1->ord[1] = nOrdBack2;  /* ord of another ring neighbor, not traversed yet */
    // INCHI✔️❌:                     da1->typ[1] = nBackType2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da1->ord[0] = nOrdBack2;  /* ord of a ring neighbor, traversed */
    // INCHI✔️❌:                     da1->typ[0] = nBackType2;
    // INCHI✔️❌:                     da1->ord[1] = nOrdBack1;  /* ord of another ring neighbor, not traversed yet */
    // INCHI✔️❌:                     da1->typ[1] = nBackType1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 return nBackType1 | nBackType2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: get_traversed_deriv_type

    *output = DERIV_AT::default();
    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if atom.cFlags & flags != 0 {
        return Ok(0);
    }
    let valence = usize::try_from(i32::from(atom.valence)).unwrap_or(0);
    let mut from_order = None;
    for order in 0..valence {
        let neighbor = atoms
            .get(usize::from(atom.neighbor[order]))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.cFlags & flags != 0 {
            from_order = Some(order);
            break;
        }
    }
    let Some(from_order) = from_order else {
        return Ok(-1);
    };
    if atom.valence == 1
        && atom.num_H != 0
        && matches!(
            atom.el_number,
            EL_NUMBER_O | EL_NUMBER_N | EL_NUMBER_S | EL_NUMBER_P
        )
    {
        return Ok(DERIV_NOT as i32);
    }
    if is_el_a_metal(i32::from(atom.el_number))? != 0 {
        return Ok(DERIV_NOT as i32);
    }
    let derivative = derivatives
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if derivative.typ[0] != 0
        && i32::from(derivative.typ[0]) & DERIV_UNEXPADABLE == i32::from(derivative.typ[0])
    {
        return Ok(0);
    }

    let other_order = usize::from(from_order == 0);
    if atom.nNumAtInRingSystem == 1
        && atom.el_number == EL_NUMBER_N
        && atom.valence == 2
        && atom.chem_bonds_valence == 3
        && atom.num_H == 0
        && atom.charge == 0
        && atom.radical == 0
    {
        let precursor_index = usize::from(atom.neighbor[from_order]);
        let oxygen_index = usize::from(atom.neighbor[other_order]);
        let precursor = atoms
            .get(precursor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let oxygen = atoms
            .get(oxygen_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if precursor.el_number == EL_NUMBER_C
            && oxygen.el_number == EL_NUMBER_O
            && atom.bond_type[from_order] == BOND_DOUBLE as u8
            && atom.bond_type[other_order] == BOND_SINGLE as u8
            && precursor.valence == 3_i8.wrapping_sub(precursor.num_H)
            && precursor.chem_bonds_valence == 4_i8.wrapping_sub(precursor.num_H)
            && precursor.charge == 0
            && precursor.radical == 0
            && oxygen.valence == 2
            && oxygen.chem_bonds_valence == 2
            && oxygen.charge == 0
            && oxygen.radical == 0
        {
            let precursor_neighbors_are_carbon =
                (0..usize::try_from(i32::from(precursor.valence)).unwrap_or(0)).all(|order| {
                    let candidate = usize::from(precursor.neighbor[order]);
                    candidate == atom_index
                        || atoms
                            .get(candidate)
                            .is_some_and(|value| value.el_number == EL_NUMBER_C)
                });
            if precursor_neighbors_are_carbon {
                let substituent_order = usize::from(oxygen.neighbor[0] as usize == atom_index);
                let substituent_index = usize::from(oxygen.neighbor[substituent_order]);
                let substituent = atoms
                    .get(substituent_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if is_Si_IV(atoms, substituent_index as i32)? != 0 {
                    let mut carbon_iv = None;
                    let mut valid = true;
                    for order in 0..usize::try_from(i32::from(substituent.valence)).unwrap_or(0) {
                        let candidate_index = usize::from(substituent.neighbor[order]);
                        if candidate_index == oxygen_index {
                            continue;
                        }
                        let candidate = atoms
                            .get(candidate_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if candidate.el_number != EL_NUMBER_C
                            || candidate.charge != 0
                            || candidate.radical != 0
                            || candidate.chem_bonds_valence != candidate.valence
                        {
                            valid = false;
                            break;
                        }
                        if carbon_iv.is_none() && candidate.valence == 4 && candidate.num_H == 0 {
                            carbon_iv = Some(candidate_index);
                        } else if candidate.chem_bonds_valence != 1 || candidate.num_H != 3 {
                            valid = false;
                            break;
                        }
                    }
                    if valid {
                        if let Some(carbon_iv) = carbon_iv {
                            let carbon = atoms
                                .get(carbon_iv)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            for order in 0..usize::try_from(i32::from(carbon.valence)).unwrap_or(0)
                            {
                                let candidate_index = usize::from(carbon.neighbor[order]);
                                if candidate_index == substituent_index {
                                    continue;
                                }
                                let candidate = atoms
                                    .get(candidate_index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if candidate.el_number != EL_NUMBER_C
                                    || candidate.charge != 0
                                    || candidate.radical != 0
                                    || candidate.chem_bonds_valence != candidate.valence
                                    || candidate.chem_bonds_valence != 1
                                    || candidate.num_H != 3
                                {
                                    valid = false;
                                    break;
                                }
                            }
                            if valid {
                                output.ord[0] = other_order as i8;
                                output.typ[0] = DERIV_X_OXIME as i16;
                                output.num[0] = 8;
                                return Ok(DERIV_X_OXIME as i32);
                            }
                        } else {
                            output.ord[0] = other_order as i8;
                            output.typ[0] = DERIV_X_OXIME as i16;
                            output.num[0] = 5;
                            return Ok(DERIV_X_OXIME as i32);
                        }
                    }
                } else if substituent.el_number == EL_NUMBER_C {
                    if substituent.chem_bonds_valence == 1
                        && substituent.num_H == 3
                        && substituent.charge == 0
                        && substituent.radical == 0
                    {
                        output.ord[0] = other_order as i8;
                        output.typ[0] = DERIV_X_OXIME as i16;
                        output.num[0] = 2;
                        return Ok(DERIV_X_OXIME as i32);
                    }
                    if substituent.valence == 2
                        && substituent.chem_bonds_valence == 2
                        && substituent.num_H == 2
                        && substituent.charge == 0
                        && substituent.radical == 0
                    {
                        let next_order =
                            usize::from(substituent.neighbor[0] as usize == oxygen_index);
                        let next_index = usize::from(substituent.neighbor[next_order]);
                        let next = atoms
                            .get(next_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if next.chem_bonds_valence == 1
                            && next.num_H == 3
                            && next.charge == 0
                            && next.radical == 0
                        {
                            output.ord[0] = other_order as i8;
                            output.typ[0] = DERIV_X_OXIME as i16;
                            output.num[0] = 3;
                            return Ok(DERIV_X_OXIME as i32);
                        }
                        if next.valence == 3
                            && next.chem_bonds_valence == 4
                            && next.num_H == 0
                            && next.charge == 0
                            && next.radical == 0
                            && next.nRingSystem != substituent.nRingSystem
                            && next.bCutVertex != 0
                            && next.nNumAtInRingSystem == 6
                            && is_Phenyl(atoms, substituent_index as i32, next_index as i32)? != 0
                        {
                            output.ord[0] = other_order as i8;
                            output.typ[0] = DERIV_X_OXIME as i16;
                            output.num[0] = 8;
                            return Ok(DERIV_X_OXIME as i32);
                        }
                    }
                }
            }
        }
    }

    if atom.nNumAtInRingSystem == 1
        && (((atom.el_number == EL_NUMBER_O || atom.el_number == EL_NUMBER_S) && atom.valence == 2)
            || (atom.el_number == EL_NUMBER_N && atom.valence == 2 && atom.num_H == 1)
            || (atom.valence == 3 && atom.num_H == 2))
    {
        let mut found = 0_i32;
        let mut dansyl = DERIV_AT::default();
        for order in 0..valence {
            if order == from_order {
                continue;
            }
            let sulfur_index = usize::from(atom.neighbor[order]);
            let sulfur = atoms
                .get(sulfur_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if sulfur.el_number == EL_NUMBER_S
                && sulfur.valence == 4
                && sulfur.chem_bonds_valence == 6
            {
                let current_derivative = derivatives
                    .get_mut(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let mut candidate = DERIV_AT::default();
                if is_Dansyl(
                    atoms,
                    atom_index as i32,
                    order as i32,
                    Some(current_derivative),
                    Some(&mut candidate),
                )? != 0
                {
                    found += 1;
                    dansyl = candidate;
                }
            }
        }
        if found == 1 {
            *output = dansyl;
            return Ok(DERIV_DANSYL as i32);
        }
    }

    if atom.nNumAtInRingSystem == 1
        && (atom.el_number == EL_NUMBER_O || atom.el_number == EL_NUMBER_S)
        && atom.valence == 2
        && atom.chem_bonds_valence == 2
        && atom.num_H == 0
        && atom.charge == 0
        && atom.radical == 0
    {
        let agent_index = usize::from(atom.neighbor[other_order]);
        let precursor_index = usize::from(atom.neighbor[from_order]);
        let has_agent_carbonyl = is_C_or_S_DB_O(atoms, agent_index as i32)? != 0;
        let precursor_is_silicon = is_Si_IV(atoms, precursor_index as i32)? != 0;
        if has_agent_carbonyl
            && is_C_DB_O(atoms, agent_index as i32)? != 0
            && is_silyl2(atoms, precursor_index as i32, atom_index as i32)? == 0
            && is_el_a_metal(i32::from(
                atoms
                    .get(precursor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .el_number,
            ))? == 0
        {
            let agent = atoms
                .get(agent_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for order in 0..usize::try_from(i32::from(agent.valence)).unwrap_or(0) {
                let substituent_index = usize::from(agent.neighbor[order]);
                if substituent_index == atom_index || agent.bond_type[order] != BOND_SINGLE as u8 {
                    continue;
                }
                let mut count =
                    is_CF3_or_linC3F7a(atoms, substituent_index as i32, agent_index as i32)?;
                if count != 0 {
                    count = 2_i32
                        .wrapping_add(3_i32.wrapping_mul(count))
                        .wrapping_add(1);
                } else if is_Methyl(atoms, substituent_index as i32)? != 0 {
                    if is_Methyl(atoms, precursor_index as i32)? == 0
                        && is_Ethyl(atoms, atom_index as i32, precursor_index as i32)? == 0
                    {
                        count = 3;
                    }
                } else if is_Phenyl(atoms, agent_index as i32, substituent_index as i32)? != 0 {
                    if is_Methyl(atoms, precursor_index as i32)? == 0
                        && is_Ethyl(atoms, atom_index as i32, precursor_index as i32)? == 0
                    {
                        count = 8;
                    }
                } else if is_PentaFluoroPhenyl(atoms, agent_index as i32, substituent_index as i32)?
                    != 0
                    && is_Methyl(atoms, precursor_index as i32)? == 0
                    && is_Ethyl(atoms, atom_index as i32, precursor_index as i32)? == 0
                {
                    count = 13;
                }
                if count != 0 {
                    output.ord[0] = other_order as i8;
                    output.typ[0] = DERIV_RO_COX as i16;
                    output.num[0] = count as i8;
                    return Ok(DERIV_RO_COX as i32);
                }
                break;
            }
        }
        if !(precursor_is_silicon && has_agent_carbonyl)
            && is_possibly_deriv_neigh(
                atoms,
                atom_index as i32,
                other_order as i32,
                DERIV_BRIDGE_O as i32,
                flags,
            )? != 0
        {
            output.ord[0] = other_order as i8;
            output.typ[0] = DERIV_BRIDGE_O as i16;
            return Ok(DERIV_BRIDGE_O as i32);
        }
    }

    if atom.el_number == EL_NUMBER_O
        && atom.nNumAtInRingSystem == 5
        && is_DERIV_RING_DMOX_DEOX_O(
            atoms,
            atom_index as i32,
            from_order as i32,
            derivatives.get_mut(atom_index),
            Some(output),
        )? != 0
    {
        return Ok(DERIV_RING_DMOX_DEOX_O as i32);
    }
    if atom.el_number == EL_NUMBER_N
        && atom.nNumAtInRingSystem == 5
        && atom.valence == 2
        && atom.chem_bonds_valence == 3
        && is_DERIV_RING_DMOX_DEOX_N(
            atoms,
            atom_index as i32,
            from_order as i32,
            derivatives.get_mut(atom_index),
            Some(output),
        )? != 0
    {
        return Ok(DERIV_RING_DMOX_DEOX_N as i32);
    }

    if atom.bCutVertex != 0
        && atom.el_number == EL_NUMBER_N
        && atom.valence == 2
        && atom.chem_bonds_valence == atom.valence
        && i32::from(atom.valence) + i32::from(atom.num_H) == 3
        && atom.charge == 0
        && atom.radical == 0
    {
        output.ord[0] = other_order as i8;
        output.typ[0] = DERIV_BRIDGE_NH as i16;
        return Ok(DERIV_BRIDGE_NH as i32);
    }

    if atom.bCutVertex != 0
        && atom.el_number == EL_NUMBER_N
        && atom.valence == 3
        && atom.chem_bonds_valence == atom.valence
        && i32::from(atom.valence) + i32::from(atom.num_H) == 3
        && atom.charge == 0
        && atom.radical == 0
    {
        let from_index = usize::from(atom.neighbor[from_order]);
        let from = atoms
            .get(from_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if from.el_number == EL_NUMBER_C && from.valence == 3 && from.chem_bonds_valence == 4 {
            let ring_type = is_DERIV_RING2_PRRLDD_PPRDN(
                atoms,
                atom_index as i32,
                from_order as i32,
                derivatives.get_mut(atom_index),
                Some(output),
            )?;
            if ring_type != 0 {
                return Ok(ring_type);
            }
        }

        let neighbor_ring = |order: usize| -> Result<u16, SourceHeapError> {
            Ok(atoms
                .get(usize::from(atom.neighbor[order]))
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .nRingSystem)
        };
        let order1 = (from_order + 1) % 3;
        let order2 = (from_order + 2) % 3;
        let rm1 = neighbor_ring(from_order)? == neighbor_ring(order1)?;
        let rm2 = neighbor_ring(from_order)? == neighbor_ring(order2)?;
        let r12 = neighbor_ring(order1)? == neighbor_ring(order2)?;
        let mut orders = [-1_i32; 2];
        let mut count = 0_usize;
        match i32::from(rm1) + i32::from(rm2) + i32::from(r12) {
            0 => {
                for order in [order1, order2] {
                    if is_possibly_deriv_neigh(
                        atoms,
                        atom_index as i32,
                        order as i32,
                        DERIV_AMINE_tN as i32,
                        flags,
                    )? != 0
                    {
                        orders[count] = order as i32;
                        count += 1;
                    }
                }
                if count == 2 && orders[0] > orders[1] {
                    orders.swap(0, 1);
                }
            }
            1 => {
                let candidate = if rm1 {
                    Some(order2)
                } else if rm2 {
                    Some(order1)
                } else {
                    None
                };
                if let Some(order) = candidate
                    && is_possibly_deriv_neigh(
                        atoms,
                        atom_index as i32,
                        order as i32,
                        DERIV_AMINE_tN as i32,
                        flags,
                    )? != 0
                {
                    orders[0] = order as i32;
                    count = 1;
                }
            }
            _ => {}
        }
        for index in 0..count {
            output.ord[index] = orders[index] as i8;
            output.typ[index] = DERIV_AMINE_tN as i16;
        }
        if count != 0 {
            return Ok(DERIV_AMINE_tN as i32);
        }
        return Ok(0);
    }

    if atom.bCutVertex != 0
        && atom.valence == atom.chem_bonds_valence
        && (atom.num_H == 0 || (atom.el_number == EL_NUMBER_C && atom.num_H == 1))
        && atom.charge == 0
        && atom.radical == 0
        && ((atom.el_number == EL_NUMBER_C && i32::from(atom.valence) + i32::from(atom.num_H) == 4)
            || (atom.el_number == EL_NUMBER_SI && atom.valence == 4)
            || (atom.el_number == EL_NUMBER_B && atom.valence == 3))
    {
        let first_connector = usize::from(atom.neighbor[from_order]);
        if let Some((first_type, block_system)) =
            classify_ring_back_connector(atoms, atom_index, first_connector)?
            && block_system != 0
        {
            atoms[atom_index].cFlags |= CFLAG_MARK_BLOCK as i8;
            let first_neighbor = i32::from(atom.neighbor[from_order]);
            let marked = mark_atoms_cFlags(atoms, first_neighbor, 1, CFLAG_MARK_BLOCK as i8)?;
            let current = atoms
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let mut second = None;
            let mut third = false;
            for order in 0..valence {
                if order == from_order {
                    continue;
                }
                let connector = atoms
                    .get(usize::from(current.neighbor[order]))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if connector.cFlags & CFLAG_MARK_BLOCK as i8 == 0 {
                    continue;
                }
                if let Some((type_, _)) = classify_ring_back_connector(
                    atoms,
                    atom_index,
                    usize::from(current.neighbor[order]),
                )? {
                    if second.is_none() {
                        second = Some((order, type_));
                    } else {
                        third = true;
                    }
                } else {
                    third = true;
                    break;
                }
            }
            let unmarked = unmark_atoms_cFlags(
                atoms,
                atom_index as i32,
                0,
                CFLAG_MARK_BLOCK as i8,
                !(CFLAG_MARK_BLOCK as i8),
            )?;
            if marked != unmarked {
                return Ok(-1);
            }
            if let Some((second_order, second_type)) = second
                && !third
            {
                let mut entries = [(from_order, first_type), (second_order, second_type)];
                entries.sort_by_key(|entry| entry.0);
                output.ord[0] = entries[0].0 as i8;
                output.typ[0] = entries[0].1;
                output.ord[1] = entries[1].0 as i8;
                output.typ[1] = entries[1].1;
                return Ok(i32::from(first_type | second_type));
            }
        }
    }

    Ok(0)
}

pub(crate) fn is_possibly_deriv_neigh(
    atoms: &[inp_ATOM],
    atom_index: i32,
    order: i32,
    type_: i32,
    _flags: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2610 is_possibly_deriv_neigh
    // INCHI✔️❌: int is_possibly_deriv_neigh( inp_ATOM *at,
    // INCHI✔️❌:                              int iat,
    // INCHI✔️❌:                              int iord,
    // INCHI✔️❌:                              int type,
    // INCHI✔️❌:                              char cFlags )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int neigh = at[iat].neighbor[iord]; /* inside derivatizing agent */
    // INCHI✔️❌:     int neigh_from = -1;
    // INCHI✔️❌:     U_CHAR el = at[neigh].el_number;
    // INCHI✔️❌:     int    bOk = 0;
    // INCHI✔️❌:     switch (type)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         case DERIV_BRIDGE_O:
    // INCHI✔️❌:             neigh_from = at[iat].neighbor[!iord]; /* inside precursor
    // INCHI✔️❌:                                                   neigh_from  iat
    // INCHI✔️❌:                                                          -> A--O--B -> traversing from A(neigh_from) to B(neigh); may we cut O--B bond? */
    // INCHI✔️❌:                                                   /* do not cut bond "---" in A=Si(IV), B(=O), B=C: Si(IV)-O---B(=O) */
    // INCHI✔️❌:             if (!( is_C_or_S_DB_O( at, neigh ) && is_Si_IV( at, neigh_from ) ) &&
    // INCHI✔️❌:                  !is_C_unsat_not_arom( at, neigh ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bOk = ( el == EL_NUMBER_C ||
    // INCHI✔️❌:                         el == EL_NUMBER_SI ||
    // INCHI✔️❌:                         el == EL_NUMBER_S ||
    // INCHI✔️❌:                         el == EL_NUMBER_P ) &&
    // INCHI✔️❌:                     is_deriv_chain2( at, iat, DERIV_BRIDGE_O, -1, iord, 0, NULL, 0, NULL, 0, NULL );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef DERIV_RO_COX
    // INCHI✔️❌:         case DERIV_RO_COX:
    // INCHI✔️❌:             /*           iord
    // INCHI✔️❌:             -> R-O--[C(=O)-B]; -B: -CH3, C[n]F[2n+1] 0 < n < 4; may we cut O--C bond? */
    // INCHI✔️❌:             neigh_from = at[iat].neighbor[!iord];
    // INCHI✔️❌:             if (at[neigh_from].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:                  at[iat].el_number == EL_NUMBER_O &&
    // INCHI✔️❌:                  at[neigh].el_number == EL_NUMBER_C &&
    // INCHI✔️❌:                  is_C_or_S_DB_O( at, neigh ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bOk = 1; /*is_deriv_chain2( at, iat, DERIV_RO_COX, iord, 0 ); does nothing */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         case DERIV_BRIDGE_NH:
    // INCHI✔️❌:             /* -> A--NH--B -> traversing from A(neigh_from) to B(neigh); may we cut NH--B bond? */
    // INCHI✔️❌:             bOk = ( is_C_or_S_DB_O( at, neigh ) ||
    // INCHI✔️❌:                     /*is_C_Alk( at, neigh, cFlags ) ||*/
    // INCHI✔️❌:                     is_Si_IV( at, neigh ) /*||
    // INCHI✔️❌:                                           is_P_TB_N( at, neigh )*/ ) && !( is_C_unsat_not_arom( at, neigh ) ) &&
    // INCHI✔️❌:                 is_deriv_chain2( at, iat, DERIV_BRIDGE_NH, -1, iord, 0, NULL, 0, NULL, 0, NULL );
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case DERIV_AMINE_tN:
    // INCHI✔️❌:             bOk = ( is_C_or_S_DB_O( at, neigh ) ||
    // INCHI✔️❌:                     /*is_C_Alk( at, neigh, cFlags ) ||*/
    // INCHI✔️❌:                     is_Si_IV( at, neigh ) /*||
    // INCHI✔️❌:                                           is_P_TB_N( at, neigh )*/ ) && !( is_C_unsat_not_arom( at, neigh ) ) &&
    // INCHI✔️❌:                 is_deriv_chain2( at, iat, DERIV_AMINE_tN, -1, iord, 0, NULL, 0, NULL, 0, NULL );
    // INCHI✔️❌:             break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return bOk;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /****************************************************************************
    // END INCHI C FUNCTION: is_possibly_deriv_neigh

    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let order_index = usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let neighbor_index = usize::from(
        *atom
            .neighbor
            .get(order_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let neighbor = atoms
        .get(neighbor_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let result = match type_ {
        value if value == DERIV_BRIDGE_O as i32 => {
            let previous_index = usize::from(atom.neighbor[usize::from(order == 0)]);
            let blocked_carbonyl = is_C_or_S_DB_O(atoms, neighbor_index as i32)? != 0
                && is_Si_IV(atoms, previous_index as i32)? != 0;
            if !blocked_carbonyl
                && is_C_unsat_not_arom(atoms, neighbor_index as i32)? == 0
                && matches!(
                    neighbor.el_number,
                    EL_NUMBER_C | EL_NUMBER_SI | EL_NUMBER_S | 15
                )
            {
                i32::from(
                    is_deriv_chain2(
                        atoms,
                        atom_index as i32,
                        DERIV_BRIDGE_O as i32,
                        -1,
                        order,
                        0,
                        None,
                        0,
                        None,
                        0,
                        None,
                    )? != 0,
                )
            } else {
                0
            }
        }
        value if value == DERIV_RO_COX as i32 => {
            let previous_index = usize::from(atom.neighbor[usize::from(order == 0)]);
            let previous = atoms
                .get(previous_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            i32::from(
                previous.el_number == EL_NUMBER_C
                    && atom.el_number == EL_NUMBER_O
                    && neighbor.el_number == EL_NUMBER_C
                    && is_C_or_S_DB_O(atoms, neighbor_index as i32)? != 0,
            )
        }
        value if value == DERIV_BRIDGE_NH as i32 || value == DERIV_AMINE_tN as i32 => {
            if (is_C_or_S_DB_O(atoms, neighbor_index as i32)? != 0
                || is_Si_IV(atoms, neighbor_index as i32)? != 0)
                && is_C_unsat_not_arom(atoms, neighbor_index as i32)? == 0
            {
                i32::from(
                    is_deriv_chain2(
                        atoms,
                        atom_index as i32,
                        type_,
                        -1,
                        order,
                        0,
                        None,
                        0,
                        None,
                        0,
                        None,
                    )? != 0,
                )
            } else {
                0
            }
        }
        _ => 0,
    };
    Ok(result)
}

fn underiv_add_text(
    target: &mut Option<&mut [i8]>,
    capacity: i32,
    text: &str,
    delimiter: i8,
) -> Result<(), SourceHeapError> {
    let mut source: Vec<i8> = text.bytes().map(|byte| byte as i8).collect();
    source.push(0);
    underiv_list_add(target.as_deref_mut(), capacity, Some(&source), delimiter)?;
    Ok(())
}

fn underiv_add_c_string(
    target: &mut Option<&mut [i8]>,
    capacity: i32,
    source: &[i8],
    delimiter: i8,
) -> Result<(), SourceHeapError> {
    underiv_list_add(target.as_deref_mut(), capacity, Some(source), delimiter)?;
    Ok(())
}

fn atom_element_name(atom: &inp_ATOM) -> Result<String, SourceHeapError> {
    let length = atom
        .elname
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(atom.elname[..length]
        .iter()
        .map(|&byte| byte as u8 as char)
        .collect())
}

fn underiv_or_bit(
    bit_underiv: &mut Option<&mut BIT_UNDERIV>,
    value: u32,
) -> Result<(), SourceHeapError> {
    let output = bit_underiv
        .as_deref_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *output |= value as BIT_UNDERIV;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn is_deriv_chain2(
    atoms: &[inp_ATOM],
    start: i32,
    type_: i32,
    num: i32,
    ord: i32,
    idrv: i32,
    mut underiv: Option<&mut [i8]>,
    len_underiv: i32,
    mut underiv2: Option<&mut [i8]>,
    len_underiv2: i32,
    mut bit_underiv: Option<&mut BIT_UNDERIV>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4483 is_deriv_chain2
    // INCHI✔️❌: int is_deriv_chain2( inp_ATOM *at,
    // INCHI✔️❌:                      int start,
    // INCHI✔️❌:                      int type,
    // INCHI✔️❌:                      int num,
    // INCHI✔️❌:                      int ord,
    // INCHI✔️❌:                      int idrv,
    // INCHI✔️❌:                      char *szUnderiv,
    // INCHI✔️❌:                      int lenUnderiv,
    // INCHI✔️❌:                      char *szUnderiv2,
    // INCHI✔️❌:                      int lenUnderiv2,
    // INCHI✔️❌:                      BIT_UNDERIV *bitUnderiv )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, k, prev_ord, neigh, iC, iO /* O or N */, iNxt, n1 = 0, n2 = 0;
    // INCHI✔️❌:     AT_NUMB *p;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!type || ( type & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     #ifdef DERIV_RING2_OUTSIDE_PRECUR
    // INCHI✔️❌:     if (type & DERIV_RING2_OUTSIDE_PRECUR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     #endif
    // INCHI✔️❌:     */
    // INCHI✔️❌:     /* reject unexpected unsaturated */
    // INCHI✔️❌:     if (at[start].charge || at[start].radical || (at[start].valence != at[start].chem_bonds_valence
    // INCHI✔️❌: #ifdef DERIV_X_OXIME
    // INCHI✔️❌:          && type != DERIV_X_OXIME
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_RO_COX
    // INCHI✔️❌:          && type != DERIV_RO_COX
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef DERIV_RING_DMOX_DEOX_N
    // INCHI✔️❌:          && type != DERIV_RING_DMOX_DEOX_N
    // INCHI✔️❌: #endif
    // INCHI✔️❌:          )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     neigh = at[start].neighbor[(int) ord];
    // INCHI✔️❌:     p = is_in_the_list( at[neigh].neighbor, (AT_NUMB) start, at[neigh].valence );
    // INCHI✔️❌:     if (!p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     prev_ord = p - at[neigh].neighbor;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* eliminate silyl possibility */
    // INCHI✔️❌:     if (type == DERIV_BRIDGE_O || type == DERIV_BRIDGE_NH || type == DERIV_AMINE_tN)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((n1 = is_silyl( at, neigh, prev_ord ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[start].valence != 2 || /* amine ? */
    // INCHI✔️❌:                                           /*at[start].el_number == EL_NUMBER_O && */
    // INCHI✔️❌:                  at[at[start].neighbor[!ord]].el_number != EL_NUMBER_SI)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Gari's request: disconnect only from C-O->-Si..., not Si-O->-Si... (e.g.  CASr.n.= 141-63-9 ). */
    // INCHI✔️❌:                 /* ???? in case of type == DERIV_AMINE_tN why don't we check more neighbors ???? */
    // INCHI✔️❌:                 if (NULL != szUnderiv && 0 < lenUnderiv)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     char szPrecur[16] = "R";
    // INCHI✔️❌:                     switch (type)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         case DERIV_BRIDGE_O:
    // INCHI✔️❌:                             strcat(szPrecur, at[start].elname);
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case DERIV_BRIDGE_NH:
    // INCHI✔️❌:                             strcat(szPrecur, "NH");
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case DERIV_AMINE_tN:
    // INCHI✔️❌:                             strcat(szPrecur, "N");
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         default:
    // INCHI✔️❌:                             strcat(szPrecur, "??");
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     strcat(szPrecur, "-");
    // INCHI✔️❌:                     switch (n1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         case 1:
    // INCHI✔️❌:                             strcat(szPrecur, "TMS");
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case 2:
    // INCHI✔️❌:                             strcat(szPrecur, "TBDMS");
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         default:
    // INCHI✔️❌:                             strcat(szPrecur, "???");
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     underiv_list_add( szUnderiv, lenUnderiv, szPrecur, ' ' );
    // INCHI✔️❌:                     underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[n1 == 1 ? DERIV_ID_TMS : n1 == 2 ? DERIV_ID_TBDMS : DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                     *bitUnderiv |= n1 == 1 ? DERIV_BIT_TMS : n1 == 2 ? DERIV_BIT_TBDMS : DERIV_BIT_Unknown;
    // INCHI✔️❌:
    // INCHI✔️❌:                     /*
    // INCHI✔️❌:                     underiv_list_add( szUnderiv, lenUnderiv, type == DERIV_BRIDGE_O? "RO-" : type == DERIV_BRIDGE_NH? "RNH-" : "RN-", ' ' );
    // INCHI✔️❌:                     underiv_list_add( szUnderiv, lenUnderiv, n1==1?"TMS" : n1==2? "TBDMS" : "???", 0 );
    // INCHI✔️❌:                     */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #ifdef UNDERIV_SYLYL_ONLY
    // INCHI✔️❌:     return 0; /* if it is not Sylyl then it is not a derivative */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = n2 = 0; /* djb-rwth: ignoring LLVM warning: variables used */
    // INCHI✔️❌:     if (type == DERIV_BRIDGE_O)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* check acetyl */
    // INCHI✔️❌: #if ( !defined(UNDERIV_ACETATE_Me) && !defined(UNDERIV_ACETATE_Et) && !defined(UNDERIV_ACETATE_CnF2np1) )
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌: #else
    // INCHI✔️❌:         iC = at[start].neighbor[!ord];
    // INCHI✔️❌:         if (at[iC].charge || at[iC].radical || at[iC].num_H ||
    // INCHI✔️❌:              at[iC].el_number != EL_NUMBER_C || at[iC].valence != 3 ||
    // INCHI✔️❌:              at[iC].valence + 1 != at[iC].chem_bonds_valence)
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         for (i = k = 0; i < at[iC].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iO = at[iC].neighbor[i];
    // INCHI✔️❌:             if (at[iO].charge || at[iO].radical || at[iO].num_H ||
    // INCHI✔️❌:                  at[iO].el_number != EL_NUMBER_O || at[iO].valence != 1 ||
    // INCHI✔️❌:                  at[iO].valence + 1 != at[iO].chem_bonds_valence)
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             k++; /* number of =O */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (k != 1)
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         /* check derivative */
    // INCHI✔️❌: #if defined(UNDERIV_ACETATE_Et) || defined(UNDERIV_ACETATE_Me)
    // INCHI✔️❌:         n1 = is_Me_or_Et( at, neigh, prev_ord );
    // INCHI✔️❌:         if (
    // INCHI✔️❌: #if defined(UNDERIV_ACETATE_Et) && defined(UNDERIV_ACETATE_Me)
    // INCHI✔️❌:             0 !=
    // INCHI✔️❌: #elif defined(UNDERIV_ACETATE_Et)
    // INCHI✔️❌:             2 ==
    // INCHI✔️❌: #elif defined(UNDERIV_ACETATE_Me)
    // INCHI✔️❌:             1 ==
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             n1)
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         else
    // INCHI✔️❌:             n1 = 0;
    // INCHI✔️❌: #endif /* defined(UNDERIV_ACETATE_Et) || defined(UNDERIV_ACETATE_Me) */
    // INCHI✔️❌: #ifdef UNDERIV_ACETATE_CnF2np1
    // INCHI✔️❌:         if (n1 <= 0)
    // INCHI✔️❌:             n2 = is_CF3_or_linC3F7( at, neigh, prev_ord );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         if (n1 || n2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (szUnderiv)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv, "RCOO-", ' ' );
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv,
    // INCHI✔️❌:                                   n1 == 1 ? "Me" :
    // INCHI✔️❌:                                   n1 == 2 ? "Et" :
    // INCHI✔️❌:                                   n2 == 1 ? "CF3" :
    // INCHI✔️❌:                                   n2 == 2 ? "C2F5" :
    // INCHI✔️❌:                                   n2 == 3 ? "C3F7" :
    // INCHI✔️❌:                                   "C?F?", 0 );
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[
    // INCHI✔️❌: #if  defined(UNDERIV_ACETATE_Me)
    // INCHI✔️❌:                     n1 == 1 ? DERIV_ID_Methylation :
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if  defined(UNDERIV_ACETATE_Et)
    // INCHI✔️❌:                         n1 == 2 ? DERIV_ID_Ethylation :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         n2 == 1 ? DERIV_ID_TFA :
    // INCHI✔️❌:                         n2 == 2 ? DERIV_ID_PFP :
    // INCHI✔️❌:                         n2 == 3 ? DERIV_ID_HFB :
    // INCHI✔️❌:                         DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                 *bitUnderiv |=
    // INCHI✔️❌: #if  defined(UNDERIV_ACETATE_Me)
    // INCHI✔️❌:                     n1 == 1 ? DERIV_BIT_Methylation :
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if  defined(UNDERIV_ACETATE_Et)
    // INCHI✔️❌:                     n1 == 2 ? DERIV_BIT_Ethylation :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     n2 == 1 ? DERIV_BIT_TFA :
    // INCHI✔️❌:                     n2 == 2 ? DERIV_BIT_PFP :
    // INCHI✔️❌:                     n2 == 3 ? DERIV_BIT_HFB :
    // INCHI✔️❌:                     DERIV_BIT_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌: #endif /* !defined(UNDERIV_ACETATE_Me) && !defined(UNDERIV_ACETATE_Et) && !defined(UNDERIV_ACETATE_CnF2np1) */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = n2 = 0;
    // INCHI✔️❌:     if (type == DERIV_BRIDGE_NH || type == DERIV_AMINE_tN)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* check acetyl */
    // INCHI✔️❌:         iNxt = -1;
    // INCHI✔️❌:         iC = at[start].neighbor[(int) ord];
    // INCHI✔️❌:         if (at[iC].charge || at[iC].radical || at[iC].num_H ||
    // INCHI✔️❌:              at[iC].el_number != EL_NUMBER_C || at[iC].valence != 3 ||
    // INCHI✔️❌:              at[iC].valence + 1 != at[iC].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (i = k = 0; i < at[iC].valence; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iO = at[iC].neighbor[i];
    // INCHI✔️❌:             if (at[iO].charge || at[iO].radical || at[iO].num_H ||
    // INCHI✔️❌:                  at[iO].el_number != EL_NUMBER_O || at[iO].valence != 1 ||
    // INCHI✔️❌:                  at[iO].valence + 1 != at[iO].chem_bonds_valence)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (iO != start)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (iNxt < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iNxt = iO;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             k++; /* number of =O */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (k != 1 || iNxt < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* find bond from iNxt to iC */
    // INCHI✔️❌:         p = is_in_the_list( at[iNxt].neighbor, (AT_NUMB) iC, at[iNxt].valence );
    // INCHI✔️❌:         if (!p)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return -1; /* program error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         prev_ord = p - at[iNxt].neighbor;
    // INCHI✔️❌:         /* check derivative */
    // INCHI✔️❌: #if ( defined( UNDERIV_RN_AcEt ) || defined( UNDERIV_RNH_AcEt )  )
    // INCHI✔️❌:         n1 = is_Me_or_Et( at, iNxt, prev_ord );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         if (!n1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             n2 = is_CF3_or_linC3F7( at, iNxt, prev_ord );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (n1 == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌: #if ( !defined(UNDERIV_RN_AcMe) && defined(DERIV_BRIDGE_tN) )
    // INCHI✔️❌:                 if (type == DERIV_BRIDGE_tN) n1 = 0;
    // INCHI✔️❌: #elif ( !defined(UNDERIV_RNH_AcMe) && defined(DERIV_BRIDGE_NH) )
    // INCHI✔️❌:                 if (type == DERIV_BRIDGE_NH) n1 = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 ; /* keep C-compiler happy */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (n1 == 2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌: #if ( !defined(UNDERIV_RN_AcEt) && defined(DERIV_BRIDGE_tN) )
    // INCHI✔️❌:                     if (type == DERIV_BRIDGE_tN) n1 = 0;
    // INCHI✔️❌: #elif ( !defined(UNDERIV_RNH_AcEt) && defined(DERIV_BRIDGE_NH) )
    // INCHI✔️❌:                     if (type == DERIV_BRIDGE_NH) n1 = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     ; /* keep C-compiler happy */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     n1 = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (n1 || n2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (szUnderiv)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv, type == DERIV_AMINE_tN ? "RN-C(O)" : "RNH-C(O)", ' ' );
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv,
    // INCHI✔️❌:                                   n1 == 1 ? "Me" :
    // INCHI✔️❌:                                   n1 == 2 ? "Et" :
    // INCHI✔️❌:                                   n2 == 1 ? "CF3" :
    // INCHI✔️❌:                                   n2 == 2 ? "C2F5" :
    // INCHI✔️❌:                                   n2 == 3 ? "C3F7" :
    // INCHI✔️❌:                                   "C?F?", 0 );
    // INCHI✔️❌:                 /* djb-rwth: addressing coverity ID #499506 -- condition is correct for n1 != 1 */
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[
    // INCHI✔️❌: #if defined(UNDERIV_RN_AcMe) || defined(UNDERIV_RNH_AcMe)
    // INCHI✔️❌:                     n1 == 1 ? DERIV_ID_Acetate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if defined(UNDERIV_RN_AcEt) || defined(UNDERIV_RNH_AcEt)
    // INCHI✔️❌:                         n1 == 2 ? DERIV_ID_Propanoate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         n2 == 1 ? DERIV_ID_TFA :
    // INCHI✔️❌:                         n2 == 2 ? DERIV_ID_PFP :
    // INCHI✔️❌:                         n2 == 3 ? DERIV_ID_HFB :
    // INCHI✔️❌:                         DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                 *bitUnderiv |=
    // INCHI✔️❌: #if defined(UNDERIV_RN_AcMe) || defined(UNDERIV_RNH_AcMe)
    // INCHI✔️❌:                     n1 == 1 ? DERIV_BIT_Acetate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if defined(UNDERIV_RN_AcEt) || defined(UNDERIV_RNH_AcEt)
    // INCHI✔️❌:                     n1 == 2 ? DERIV_BIT_Propanoate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     n2 == 1 ? DERIV_BIT_TFA :
    // INCHI✔️❌:                     n2 == 2 ? DERIV_BIT_PFP :
    // INCHI✔️❌:                     n2 == 3 ? DERIV_BIT_HFB :
    // INCHI✔️❌:                     DERIV_BIT_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*underiv_buf_clear( szUnderiv );*/
    // INCHI✔️❌: #ifdef DERIV_X_OXIME
    // INCHI✔️❌:     if (type == DERIV_X_OXIME)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (szUnderiv)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iO = 0;
    // INCHI✔️❌:             if (num == 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv, "MOX", ' ' );
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_MOX], ' ' );
    // INCHI✔️❌:                 *bitUnderiv |= DERIV_BIT_MOX;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (num == 3)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     underiv_list_add( szUnderiv, lenUnderiv, "EtOX", ' ' );
    // INCHI✔️❌:                     underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_EtOX], ' ' );
    // INCHI✔️❌:                     *bitUnderiv |= DERIV_BIT_EtOX;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (num == 8)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         neigh = at[start].neighbor[ord];
    // INCHI✔️❌:                         iC = at[neigh].neighbor[start == at[neigh].neighbor[0]];
    // INCHI✔️❌:                         iO = at[iC].el_number == EL_NUMBER_SI;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     underiv_list_add( szUnderiv, lenUnderiv, "OX-", ' ' );
    // INCHI✔️❌:                     underiv_list_add( szUnderiv, lenUnderiv, num == 5 ? "TMS" : num == 8 && iO ? "TBDMS" : num == 8 && !iO ? "CH2Phe" : "???", 0 );
    // INCHI✔️❌:                     underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[num == 5 ? DERIV_ID_TMS : num == 8 && iO ? DERIV_ID_TBDMS : num == 8 && !iO ? DERIV_ID_BenzOX : DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                     *bitUnderiv |= num == 5 ? DERIV_BIT_TMS : num == 8 && iO ? DERIV_BIT_TBDMS : num == 8 && !iO ? DERIV_BIT_BenzOX : DERIV_BIT_Unknown;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef DERIV_RO_COX
    // INCHI✔️❌:     if (type == DERIV_RO_COX)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (szUnderiv)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             underiv_list_add( szUnderiv, lenUnderiv, at[start].el_number == EL_NUMBER_O ? "RO-C(O)" :
    // INCHI✔️❌:                               at[start].el_number == EL_NUMBER_S ? "RS-C(O)" : "R?-C(O)", ' ' );
    // INCHI✔️❌:             underiv_list_add( szUnderiv, lenUnderiv, num == 3 ? "Me" :
    // INCHI✔️❌:                               num == 4 ? "Et" :
    // INCHI✔️❌:                               num == 6 ? "CF3" :
    // INCHI✔️❌:                               num == 8 ? "Phe" :
    // INCHI✔️❌:                               num == 9 ? "C2F5" :
    // INCHI✔️❌:                               num == 12 ? "C3F7" :
    // INCHI✔️❌:                               num == 13 ? "PheF5" :
    // INCHI✔️❌:                               "???", 0 );
    // INCHI✔️❌:             if (num == 13)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, underiv_list_get_last( szUnderiv, ' ' ), ' ' );
    // INCHI✔️❌:                 *bitUnderiv |= DERIV_BIT_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[
    // INCHI✔️❌: #if  defined(UNDERIV_RO_COX_Me)
    // INCHI✔️❌:                     num == 3 ? DERIV_ID_Acetate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if  defined(UNDERIV_RO_COX_Et)
    // INCHI✔️❌:                         num == 4 ? DERIV_ID_Propanoate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         num == 6 ? DERIV_ID_TFA :
    // INCHI✔️❌: #if  defined(UNDERIV_RO_COX_BENZOATES)
    // INCHI✔️❌:                         num == 8 ? DERIV_ID_Benzoate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         num == 9 ? DERIV_ID_PFP :
    // INCHI✔️❌:                         num == 12 ? DERIV_ID_HFB :
    // INCHI✔️❌:                         num == 13 ? DERIV_ID_PFB :
    // INCHI✔️❌:                         DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                 *bitUnderiv |=
    // INCHI✔️❌: #if  defined(UNDERIV_RO_COX_Me)
    // INCHI✔️❌:                     num == 3 ? DERIV_BIT_Acetate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if  defined(UNDERIV_RO_COX_Et)
    // INCHI✔️❌:                     num == 4 ? DERIV_BIT_Propanoate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     num == 6 ? DERIV_BIT_TFA :
    // INCHI✔️❌: #if  defined(UNDERIV_RO_COX_BENZOATES)
    // INCHI✔️❌:                     num == 8 ? DERIV_BIT_Benzoate :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     num == 9 ? DERIV_BIT_PFP :
    // INCHI✔️❌:                     num == 12 ? DERIV_BIT_HFB :
    // INCHI✔️❌:                     num == 13 ? DERIV_BIT_PFB :
    // INCHI✔️❌:                     DERIV_ID_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    // INCHI✔️❌:     if (!idrv && ( type == DERIV_RING_DMOX_DEOX_N ||
    // INCHI✔️❌:                    type == DERIV_RING_DMOX_DEOX_O ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (szUnderiv && type == DERIV_RING_DMOX_DEOX_N)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* add only once; do not add upon DERIV_RING_DMOX_DEOX_O */
    // INCHI✔️❌:             underiv_list_add( szUnderiv, lenUnderiv, num == 4 ? "DMOX" : num == 6 ? "DEOX" : num ? "D?OX" : "???", ' ' );
    // INCHI✔️❌:             if (num == 4 || num == 6)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[num == 4 ? DERIV_ID_DMOX : num == 6 ? DERIV_ID_DEOX : DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                 *bitUnderiv |= num == 4 ? DERIV_BIT_DMOX : num == 6 ? DERIV_BIT_DEOX : DERIV_BIT_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, underiv_list_get_last( szUnderiv, ' ' ), ' ' );
    // INCHI✔️❌:                 *bitUnderiv |= DERIV_BIT_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef DERIV_RING2_OUTSIDE_PRECUR
    // INCHI✔️❌:     if (type && ( type & DERIV_RING2_OUTSIDE_PRECUR ) == type)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (szUnderiv && !idrv)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (num == 4 || num == 5)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv, num == 4 ? "Pyrrolidide" : num == 5 ? "Piperidine" : "???", ' ' );
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[num == 4 ? DERIV_ID_Pyrrolidide : num == 5 ? DERIV_ID_Piperidine : DERIV_ID_Unknown], ' ' ); /* djb-rwth: addressing coverity ID #499491 -- working correctly for num == 5 */
    // INCHI✔️❌:                 *bitUnderiv |= num == 4 ? DERIV_BIT_Pyrrolidide : num == 5 ? DERIV_BIT_Piperidine : DERIV_BIT_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌: #ifdef _DEBUG
    // INCHI✔️❌:                 int stop_DERIV_RING2_OUTSIDE_PRECUR = 1; /* debug only */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 underiv_list_add( szUnderiv, lenUnderiv, num ? "Pyrrol?Piper?" : "???", ' ' );
    // INCHI✔️❌:                 underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_Unknown], ' ' );
    // INCHI✔️❌:                 *bitUnderiv |= DERIV_BIT_Unknown;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef DERIV_DANSYL
    // INCHI✔️❌:     if (!idrv && type && ( type & DERIV_DANSYL ) == type)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (szUnderiv)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char szRO[16] = "R";
    // INCHI✔️❌:             strcat(szRO, at[start].elname);
    // INCHI✔️❌:             if (at[start].num_H == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(szRO, "H");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             strcat(szRO, "-Dansyl");
    // INCHI✔️❌:             underiv_list_add( szUnderiv, lenUnderiv, szRO, ' ' );
    // INCHI✔️❌:             underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_Dansyl], ' ' );
    // INCHI✔️❌:             *bitUnderiv |= DERIV_BIT_Dansyl;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: is_deriv_chain2

    let start_index = usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let start_atom = atoms
        .get(start_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if type_ == 0 || type_ & DERIV_RING_OUTSIDE_PRECURSOR as i32 != 0 {
        return Ok(0);
    }
    if start_atom.charge != 0
        || start_atom.radical != 0
        || (start_atom.valence != start_atom.chem_bonds_valence
            && type_ != DERIV_X_OXIME as i32
            && type_ != DERIV_RO_COX as i32
            && type_ != DERIV_RING_DMOX_DEOX_N as i32)
    {
        return Ok(0);
    }

    let order = usize::try_from(ord).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let neighbor_index = usize::from(
        *start_atom
            .neighbor
            .get(order)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let neighbor = atoms
        .get(neighbor_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let Some(mut previous_order) = is_in_the_list(
        Some(&neighbor.neighbor),
        start as u16,
        i32::from(neighbor.valence),
    )?
    else {
        return Ok(-1);
    };

    if type_ == DERIV_BRIDGE_O as i32
        || type_ == DERIV_BRIDGE_NH as i32
        || type_ == DERIV_AMINE_tN as i32
    {
        let silyl = is_silyl(atoms, neighbor_index as i32, previous_order as i32)?;
        if silyl != 0 {
            let opposite_order = usize::from(ord == 0);
            let opposite_is_silicon = start_atom
                .neighbor
                .get(opposite_order)
                .and_then(|&index| atoms.get(usize::from(index)))
                .is_some_and(|atom| atom.el_number == EL_NUMBER_SI);
            if start_atom.valence != 2 || !opposite_is_silicon {
                if underiv.is_some() && len_underiv > 0 {
                    let precursor = format!(
                        "R{}-{}",
                        match type_ {
                            value if value == DERIV_BRIDGE_O as i32 => {
                                atom_element_name(start_atom)?
                            }
                            value if value == DERIV_BRIDGE_NH as i32 => "NH".to_owned(),
                            value if value == DERIV_AMINE_tN as i32 => "N".to_owned(),
                            _ => "??".to_owned(),
                        },
                        match silyl {
                            1 => "TMS",
                            2 => "TBDMS",
                            _ => "???",
                        }
                    );
                    underiv_add_text(&mut underiv, len_underiv, &precursor, b' ' as i8)?;
                    let (name, bit) = match silyl {
                        1 => ("TMS", DERIV_BIT_TMS),
                        2 => ("TBDMS", DERIV_BIT_TBDMS),
                        _ => ("???", DERIV_BIT_UNKNOWN),
                    };
                    underiv_add_text(&mut underiv2, len_underiv2, name, b' ' as i8)?;
                    underiv_or_bit(&mut bit_underiv, bit)?;
                }
                return Ok(1);
            }
            return Ok(0);
        }
    }

    if type_ == DERIV_BRIDGE_O as i32 {
        return Ok(0);
    }

    if type_ == DERIV_BRIDGE_NH as i32 || type_ == DERIV_AMINE_tN as i32 {
        let carbon_index = neighbor_index;
        let carbon = atoms
            .get(carbon_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if carbon.charge != 0
            || carbon.radical != 0
            || carbon.num_H != 0
            || carbon.el_number != EL_NUMBER_C
            || carbon.valence != 3
            || i32::from(carbon.valence) + 1 != i32::from(carbon.chem_bonds_valence)
        {
            return Ok(0);
        }
        let mut next = None;
        let mut oxygens = 0_i32;
        for carbon_order in 0..usize::try_from(i32::from(carbon.valence)).unwrap_or(0) {
            let candidate_index = usize::from(carbon.neighbor[carbon_order]);
            let candidate = atoms
                .get(candidate_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let is_carbonyl_oxygen = candidate.charge == 0
                && candidate.radical == 0
                && candidate.num_H == 0
                && candidate.el_number == EL_NUMBER_O
                && candidate.valence == 1
                && i32::from(candidate.valence) + 1 == i32::from(candidate.chem_bonds_valence);
            if is_carbonyl_oxygen {
                oxygens += 1;
            } else if candidate_index != start_index {
                if next.replace(candidate_index).is_some() {
                    return Ok(0);
                }
            }
        }
        let Some(next_index) = next else {
            return Ok(0);
        };
        if oxygens != 1 {
            return Ok(0);
        }
        let next_atom = atoms
            .get(next_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let Some(found_order) = is_in_the_list(
            Some(&next_atom.neighbor),
            carbon_index as u16,
            i32::from(next_atom.valence),
        )?
        else {
            return Ok(-1);
        };
        previous_order = found_order;
        let perfluoro = is_CF3_or_linC3F7(atoms, next_index as i32, previous_order as i32)?;
        if perfluoro != 0 {
            if underiv.is_some() {
                underiv_add_text(
                    &mut underiv,
                    len_underiv,
                    if type_ == DERIV_AMINE_tN as i32 {
                        "RN-C(O)"
                    } else {
                        "RNH-C(O)"
                    },
                    b' ' as i8,
                )?;
                let suffix = match perfluoro {
                    1 => "CF3",
                    2 => "C2F5",
                    3 => "C3F7",
                    _ => "C?F?",
                };
                underiv_add_text(&mut underiv, len_underiv, suffix, 0)?;
                let (name, bit) = match perfluoro {
                    1 => ("TFA", DERIV_BIT_TFA),
                    2 => ("PFP", DERIV_BIT_PFP),
                    3 => ("HFB", DERIV_BIT_HFB),
                    _ => ("???", DERIV_BIT_UNKNOWN),
                };
                underiv_add_text(&mut underiv2, len_underiv2, name, b' ' as i8)?;
                underiv_or_bit(&mut bit_underiv, bit)?;
            }
            return Ok(1);
        }
        return Ok(0);
    }

    if type_ == DERIV_X_OXIME as i32 {
        if underiv.is_some() {
            let (primary, name, bit) = match num {
                2 => ("MOX", "MOX", DERIV_BIT_MOX),
                3 => ("EtOX", "EtOX", DERIV_BIT_ETOX),
                _ => {
                    let silicon = if num == 8 {
                        let oxime_neighbor = atoms
                            .get(neighbor_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let next_order = usize::from(start as u16 == oxime_neighbor.neighbor[0]);
                        let carbon = atoms
                            .get(usize::from(oxime_neighbor.neighbor[next_order]))
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        carbon.el_number == EL_NUMBER_SI
                    } else {
                        false
                    };
                    underiv_add_text(&mut underiv, len_underiv, "OX-", b' ' as i8)?;
                    match (num, silicon) {
                        (5, _) => ("TMS", "TMS", DERIV_BIT_TMS),
                        (8, true) => ("TBDMS", "TBDMS", DERIV_BIT_TBDMS),
                        (8, false) => ("CH2Phe", "BenzOX", DERIV_BIT_BENZOX),
                        _ => ("???", "???", DERIV_BIT_UNKNOWN),
                    }
                }
            };
            underiv_add_text(
                &mut underiv,
                len_underiv,
                primary,
                if num == 2 || num == 3 { b' ' as i8 } else { 0 },
            )?;
            underiv_add_text(&mut underiv2, len_underiv2, name, b' ' as i8)?;
            underiv_or_bit(&mut bit_underiv, bit)?;
        }
        return Ok(1);
    }

    if type_ == DERIV_RO_COX as i32 {
        if underiv.is_some() {
            let prefix = match start_atom.el_number {
                EL_NUMBER_O => "RO-C(O)",
                EL_NUMBER_S => "RS-C(O)",
                _ => "R?-C(O)",
            };
            underiv_add_text(&mut underiv, len_underiv, prefix, b' ' as i8)?;
            underiv_add_text(
                &mut underiv,
                len_underiv,
                match num {
                    3 => "Me",
                    4 => "Et",
                    6 => "CF3",
                    8 => "Phe",
                    9 => "C2F5",
                    12 => "C3F7",
                    13 => "PheF5",
                    _ => "???",
                },
                0,
            )?;
            if num == 13 {
                let primary = underiv
                    .as_deref()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let offset = underiv_list_get_last(Some(primary), b' ' as i8)?
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                underiv_add_c_string(
                    &mut underiv2,
                    len_underiv2,
                    primary
                        .get(offset..)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    b' ' as i8,
                )?;
                underiv_or_bit(&mut bit_underiv, DERIV_BIT_UNKNOWN)?;
            } else {
                let (name, bit) = match num {
                    3 => ("Acetate", DERIV_BIT_ACETATE),
                    6 => ("TFA", DERIV_BIT_TFA),
                    8 => ("Benzoate", DERIV_BIT_BENZOATE),
                    9 => ("PFP", DERIV_BIT_PFP),
                    12 => ("HFB", DERIV_BIT_HFB),
                    _ => ("???", 19),
                };
                underiv_add_text(&mut underiv2, len_underiv2, name, b' ' as i8)?;
                underiv_or_bit(&mut bit_underiv, bit)?;
            }
        }
        return Ok(1);
    }

    if idrv == 0
        && (type_ == DERIV_RING_DMOX_DEOX_N as i32 || type_ == DERIV_RING_DMOX_DEOX_O as i32)
    {
        if underiv.is_some() && type_ == DERIV_RING_DMOX_DEOX_N as i32 {
            let primary = match num {
                4 => "DMOX",
                6 => "DEOX",
                0 => "???",
                _ => "D?OX",
            };
            underiv_add_text(&mut underiv, len_underiv, primary, b' ' as i8)?;
            if num == 4 || num == 6 {
                let (name, bit) = if num == 4 {
                    ("DMOX", DERIV_BIT_DMOX)
                } else {
                    ("DEOX", DERIV_BIT_DEOX)
                };
                underiv_add_text(&mut underiv2, len_underiv2, name, b' ' as i8)?;
                underiv_or_bit(&mut bit_underiv, bit)?;
            } else {
                let primary = underiv
                    .as_deref()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let offset = underiv_list_get_last(Some(primary), b' ' as i8)?
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                underiv_add_c_string(
                    &mut underiv2,
                    len_underiv2,
                    primary
                        .get(offset..)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    b' ' as i8,
                )?;
                underiv_or_bit(&mut bit_underiv, DERIV_BIT_UNKNOWN)?;
            }
        }
        return Ok(1);
    }

    if type_ != 0 && type_ & DERIV_RING2_OUTSIDE_PRECUR as i32 == type_ {
        if underiv.is_some() && idrv == 0 {
            let (primary, name, bit) = match num {
                4 => ("Pyrrolidide", "Pyrrolidide", DERIV_BIT_PYRROLIDIDE),
                5 => ("Piperidine", "Piperidine", DERIV_BIT_PIPERIDINE),
                0 => ("???", "???", DERIV_BIT_UNKNOWN),
                _ => ("Pyrrol?Piper?", "???", DERIV_BIT_UNKNOWN),
            };
            underiv_add_text(&mut underiv, len_underiv, primary, b' ' as i8)?;
            underiv_add_text(&mut underiv2, len_underiv2, name, b' ' as i8)?;
            underiv_or_bit(&mut bit_underiv, bit)?;
        }
        return Ok(1);
    }

    if idrv == 0 && type_ != 0 && type_ & DERIV_DANSYL as i32 == type_ {
        if underiv.is_some() {
            let mut primary = format!("R{}", atom_element_name(start_atom)?);
            if start_atom.num_H == 1 {
                primary.push('H');
            }
            primary.push_str("-Dansyl");
            underiv_add_text(&mut underiv, len_underiv, &primary, b' ' as i8)?;
            underiv_add_text(&mut underiv2, len_underiv2, "Dansyl", b' ' as i8)?;
            underiv_or_bit(&mut bit_underiv, DERIV_BIT_DANSYL)?;
        }
        return Ok(1);
    }

    Ok(0)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn is_deriv_chain(
    atoms: &[inp_ATOM],
    start: i32,
    _num_atoms: i32,
    derivative: &DERIV_AT,
    derivative_index: i32,
    underiv: Option<&mut [i8]>,
    len_underiv: i32,
    underiv2: Option<&mut [i8]>,
    len_underiv2: i32,
    bit_underiv: Option<&mut BIT_UNDERIV>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4980 is_deriv_chain
    // INCHI✔️❌: int is_deriv_chain( inp_ATOM *at,
    // INCHI✔️❌:                     int start,
    // INCHI✔️❌:                     int num_atoms,
    // INCHI✔️❌:                     DERIV_AT *da1,
    // INCHI✔️❌:                     int idrv,
    // INCHI✔️❌:                     char *szUnderiv,
    // INCHI✔️❌:                     int lenUnderiv,
    // INCHI✔️❌:                     char *szUnderiv2,
    // INCHI✔️❌:                     int lenUnderiv2,
    // INCHI✔️❌:                     BIT_UNDERIV *bitUnderiv )
    // INCHI✔️❌: {
    // INCHI✔️❌:     return is_deriv_chain2( at, start, da1->typ[idrv], da1->num[idrv], da1->ord[idrv], idrv, szUnderiv, lenUnderiv, szUnderiv2, lenUnderiv2, bitUnderiv );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_deriv_chain
    let index =
        usize::try_from(derivative_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    is_deriv_chain2(
        atoms,
        start,
        i32::from(
            *derivative
                .typ
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ),
        i32::from(
            *derivative
                .num
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ),
        i32::from(
            *derivative
                .ord
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ),
        derivative_index,
        underiv,
        len_underiv,
        underiv2,
        len_underiv2,
        bit_underiv,
    )
}

pub(crate) fn is_deriv_chain_or_ring(
    atoms: &[inp_ATOM],
    start: i32,
    num_atoms: i32,
    derivative: &DERIV_AT,
    derivative_index: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4996 is_deriv_chain_or_ring
    // INCHI✔️❌: int is_deriv_chain_or_ring( inp_ATOM *at,
    // INCHI✔️❌:                             int start,
    // INCHI✔️❌:                             int num_atoms,
    // INCHI✔️❌:                             DERIV_AT *da1,
    // INCHI✔️❌:                             int *idrv )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, ret = -1;
    // INCHI✔️❌:     if (da1->typ[*idrv] & DERIV_RING_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* find the first ord of this derivative; ord of ring derivatives are in pairs */
    // INCHI✔️❌:         int j = -1;
    // INCHI✔️❌:         for (i = 0; i < DERIV_AT_LEN && da1->typ[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da1->typ[i] & DERIV_RING_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i == *idrv || i + 1 == *idrv)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *idrv = j = i;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 i++; /* bypass the second bond to the same derivatization agent */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* check consistency */
    // INCHI✔️❌:         if (j == -1 || j + 1 >= DERIV_AT_LEN ||
    // INCHI✔️❌:              !( da1->typ[j] & DERIV_RING_OUTSIDE_PRECURSOR ) || !( da1->typ[j + 1] & DERIV_RING_OUTSIDE_PRECURSOR ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = -1; /* program error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( at, start, num_atoms, da1, j, NULL, 0, NULL, 0, NULL );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (da1->typ[*idrv])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( at, start, num_atoms, da1, *idrv, NULL, 0, NULL, 0, NULL );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_deriv_chain_or_ring
    let requested =
        usize::try_from(*derivative_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let requested_type = *derivative
        .typ
        .get(requested)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if requested_type & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
        let mut found = None;
        let mut index = 0_usize;
        while index < derivative.typ.len() && derivative.typ[index] != 0 {
            if derivative.typ[index] & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
                if index == requested || index + 1 == requested {
                    *derivative_index = index as i32;
                    found = Some(index);
                    break;
                }
                index += 1;
            }
            index += 1;
        }
        let Some(found) = found else {
            return Ok(-1);
        };
        if found + 1 >= derivative.typ.len()
            || derivative.typ[found] & DERIV_RING_OUTSIDE_PRECURSOR as i16 == 0
            || derivative.typ[found + 1] & DERIV_RING_OUTSIDE_PRECURSOR as i16 == 0
        {
            return Ok(-1);
        }
        is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
            atoms,
            start,
            num_atoms,
            derivative,
            found as i32,
            None,
            0,
            None,
            0,
            None,
        )
    } else if requested_type != 0 {
        is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
            atoms,
            start,
            num_atoms,
            derivative,
            *derivative_index,
            None,
            0,
            None,
            0,
            None,
        )
    } else {
        Ok(-1)
    }
}

pub(crate) fn remove_deriv(
    derivative: &mut DERIV_AT,
    derivative_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5046 remove_deriv
    // INCHI✔️❌: int remove_deriv( DERIV_AT *da1, int idrv )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, ret = -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (da1->typ[idrv] & DERIV_RING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* find the first ord of this derivative; ord of ring derivatives are in pairs */
    // INCHI✔️❌:         j = -1;
    // INCHI✔️❌:         for (i = 0; i < DERIV_AT_LEN && da1->typ[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da1->typ[i] & DERIV_RING)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i == idrv || i + 1 == idrv)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j = i;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 i++; /* bypass the second bond to the same derivatization agent */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* delete if data are consistent */
    // INCHI✔️❌:         if (j >= 0 && j + 1 < DERIV_AT_LEN && ( da1->typ[j] & DERIV_RING ) && ( da1->typ[j + 1] & DERIV_RING ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (; j < DERIV_AT_LEN - 2 && da1->typ[j + 2]; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1->typ[j] = da1->typ[j + 2];
    // INCHI✔️❌:                 da1->num[j] = da1->num[j + 2];
    // INCHI✔️❌:                 da1->ord[j] = da1->ord[j + 2];
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (; j < DERIV_AT_LEN; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 da1->typ[j] = 0;
    // INCHI✔️❌:                 da1->num[j] = 0;
    // INCHI✔️❌:                 da1->ord[j] = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             ret = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = idrv;
    // INCHI✔️❌:         for (; j < DERIV_AT_LEN - 1 && da1->typ[j + 1]; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             da1->typ[j] = da1->typ[j + 1];
    // INCHI✔️❌:             da1->num[j] = da1->num[j + 1];
    // INCHI✔️❌:             da1->ord[j] = da1->ord[j + 1];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (; j < DERIV_AT_LEN; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             da1->typ[j] = 0;
    // INCHI✔️❌:             da1->num[j] = 0;
    // INCHI✔️❌:             da1->ord[j] = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ret = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: remove_deriv
    let requested =
        usize::try_from(derivative_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let requested_type = *derivative
        .typ
        .get(requested)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if requested_type & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
        let mut found = None;
        let mut index = 0_usize;
        while index < derivative.typ.len() && derivative.typ[index] != 0 {
            if derivative.typ[index] & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
                if index == requested || index + 1 == requested {
                    found = Some(index);
                    break;
                }
                index += 1;
            }
            index += 1;
        }
        let Some(mut index) = found else {
            return Ok(-1);
        };
        if index + 1 >= derivative.typ.len()
            || derivative.typ[index] & DERIV_RING_OUTSIDE_PRECURSOR as i16 == 0
            || derivative.typ[index + 1] & DERIV_RING_OUTSIDE_PRECURSOR as i16 == 0
        {
            return Ok(-1);
        }
        while index + 2 < derivative.typ.len() && derivative.typ[index + 2] != 0 {
            derivative.typ[index] = derivative.typ[index + 2];
            derivative.num[index] = derivative.num[index + 2];
            derivative.ord[index] = derivative.ord[index + 2];
            index += 1;
        }
        while index < derivative.typ.len() {
            derivative.typ[index] = 0;
            derivative.num[index] = 0;
            derivative.ord[index] = 0;
            index += 1;
        }
    } else {
        let mut index = requested;
        while index + 1 < derivative.typ.len() && derivative.typ[index + 1] != 0 {
            derivative.typ[index] = derivative.typ[index + 1];
            derivative.num[index] = derivative.num[index + 1];
            derivative.ord[index] = derivative.ord[index + 1];
            index += 1;
        }
        while index < derivative.typ.len() {
            derivative.typ[index] = 0;
            derivative.num[index] = 0;
            derivative.ord[index] = 0;
            index += 1;
        }
    }
    Ok(0)
}

pub(crate) fn remove_deriv_mark(
    derivative: &mut DERIV_AT,
    derivative_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5108 remove_deriv_mark
    // INCHI✔️❌: int remove_deriv_mark( DERIV_AT *da1, int idrv )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, ret = -1;
    // INCHI✔️❌:     if (da1->typ[idrv] & DERIV_RING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* find the first ord of this derivative; ord of ring derivatives are in pairs */
    // INCHI✔️❌:         j = -1;
    // INCHI✔️❌:         for (i = 0; i < DERIV_AT_LEN && da1->typ[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da1->typ[i] & DERIV_RING)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i == idrv || i + 1 == idrv)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j = i;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 i++; /* bypass the second bond to the same derivatization agent */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* delete if data are consistent */
    // INCHI✔️❌:         if (j >= 0 && j + 1 < DERIV_AT_LEN && ( da1->typ[j] & DERIV_RING ) && ( da1->typ[j + 1] & DERIV_RING ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             da1->typ[j] |= DERIV_DUPLIC;
    // INCHI✔️❌:             da1->typ[j + 1] |= DERIV_DUPLIC;
    // INCHI✔️❌:             ret = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = idrv;
    // INCHI✔️❌:         da1->typ[j] |= DERIV_DUPLIC;
    // INCHI✔️❌:         ret = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: remove_deriv_mark
    let requested =
        usize::try_from(derivative_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let requested_type = *derivative
        .typ
        .get(requested)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if requested_type & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
        let mut found = None;
        let mut index = 0_usize;
        while index < derivative.typ.len() && derivative.typ[index] != 0 {
            if derivative.typ[index] & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
                if index == requested || index + 1 == requested {
                    found = Some(index);
                    break;
                }
                index += 1;
            }
            index += 1;
        }
        let Some(index) = found else {
            return Ok(-1);
        };
        if index + 1 >= derivative.typ.len()
            || derivative.typ[index] & DERIV_RING_OUTSIDE_PRECURSOR as i16 == 0
            || derivative.typ[index + 1] & DERIV_RING_OUTSIDE_PRECURSOR as i16 == 0
        {
            return Ok(-1);
        }
        derivative.typ[index] |= DERIV_DUPLIC as i16;
        derivative.typ[index + 1] |= DERIV_DUPLIC as i16;
    } else {
        derivative.typ[requested] |= DERIV_DUPLIC as i16;
    }
    Ok(0)
}

pub(crate) fn underiv_compare(first: &[i8], second: &[i8]) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5197 underiv_compare
    // INCHI✔️❌: int underiv_compare( const void *p1, const void *p2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     return strcmp( *(const char **) p1, *(const char **) p2 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: underiv_compare
    let mut index = 0_usize;
    loop {
        let left = *first
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
        let right = *second
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
        if left != right || left == 0 {
            return Ok(i32::from(left) - i32::from(right));
        }
        index += 1;
    }
}

pub(crate) fn underiv_list_add_two_cuts(
    mut list: Option<&mut [i8]>,
    list_capacity: i32,
    underiv: &mut [i8],
    delimiter: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5204 underiv_list_add_two_cuts
    // INCHI✔️❌: int underiv_list_add_two_cuts( char *szUnderivList,
    // INCHI✔️❌:                                int lenUnderivList,
    // INCHI✔️❌:                                char *szUnderiv,
    // INCHI✔️❌:                                const char cDelim )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* may happen only in RN- case, DERIV_AMINE_tN, but not always  */
    // INCHI✔️❌:     const char szDelim[] = { cDelim, 0 };
    // INCHI✔️❌:     char *p1 = strtok( szUnderiv, szDelim );
    // INCHI✔️❌:     char *p2 = p1 ? strtok( NULL, szDelim ) : NULL;
    // INCHI✔️❌:     char *p1m = p1 ? strchr( p1, '-' ) : NULL;
    // INCHI✔️❌:     char *p2m = p2 ? strchr( p2, '-' ) : NULL;
    // INCHI✔️❌:     char *pm;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (p1m && p2m)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (p1m - p1 == p2m - p2 && !memcmp( p1, p2, p1m - p1 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* found common prefix */
    // INCHI✔️❌:             int diff;
    // INCHI✔️❌:             *p1m++ = '\0';
    // INCHI✔️❌:             *p2m++ = '\0';
    // INCHI✔️❌:             /* output
    // INCHI✔️❌:             [common prefix]-(suffix1)+(suffix2): (suffix1)<(suffix2)  or
    // INCHI✔️❌:             [common prefix]-2(suffix1):          (suffix1)==(suffix2)
    // INCHI✔️❌:             */
    // INCHI✔️❌:             underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim ); /* [common prefix] */
    // INCHI✔️❌:             underiv_list_add( szUnderivList, lenUnderivList, "-", 0 );  /* - */
    // INCHI✔️❌:             diff = strcmp( p1m, p2m );
    // INCHI✔️❌:             if (diff > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pm = p1m;
    // INCHI✔️❌:                 p1m = p2m; /* (suffix1) - smaller */
    // INCHI✔️❌:                 p2m = pm;  /* (suffix2) - greater */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (diff)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderivList, lenUnderivList, p1m, 0 );  /* (suffix1) */
    // INCHI✔️❌:                 underiv_list_add( szUnderivList, lenUnderivList, "+", 0 );  /* - */
    // INCHI✔️❌:                 underiv_list_add( szUnderivList, lenUnderivList, p2m, 0 );  /* (suffix2) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( szUnderivList, lenUnderivList, "2", 0 );  /* 2 */
    // INCHI✔️❌:                 underiv_list_add( szUnderivList, lenUnderivList, p1m, 0 );  /* (suffix1) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* should not happen */
    // INCHI✔️❌:             underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim );
    // INCHI✔️❌:             underiv_list_add( szUnderivList, lenUnderivList, p2, cDelim );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (p1 && p2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* should not happen */
    // INCHI✔️❌:             underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim );
    // INCHI✔️❌:             underiv_list_add( szUnderivList, lenUnderivList, p2, cDelim );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (p1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* happens in case only one num[] is not zero, namely, DERIV_RING2_OUTSIDE_PRECUR */
    // INCHI✔️❌:                 underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: underiv_list_add_two_cuts

    fn next_token(
        input: &mut [i8],
        cursor: &mut usize,
        delimiter: i8,
    ) -> Result<Option<(usize, usize)>, SourceHeapError> {
        if delimiter != 0 {
            loop {
                let byte = *input
                    .get(*cursor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if byte == 0 {
                    return Ok(None);
                }
                if byte != delimiter {
                    break;
                }
                *cursor += 1;
            }
        }
        let start = *cursor;
        if *input
            .get(start)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            == 0
        {
            return Ok(None);
        }
        loop {
            let byte = *input
                .get(*cursor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if byte == 0 {
                return Ok(Some((start, *cursor)));
            }
            if delimiter != 0 && byte == delimiter {
                let end = *cursor;
                input[end] = 0;
                *cursor += 1;
                return Ok(Some((start, end)));
            }
            *cursor += 1;
        }
    }

    fn add(
        list: &mut Option<&mut [i8]>,
        list_capacity: i32,
        value: &[i8],
        delimiter: i8,
    ) -> Result<(), SourceHeapError> {
        let _ = underiv_list_add(list.as_deref_mut(), list_capacity, Some(value), delimiter)?;
        Ok(())
    }

    let mut cursor = 0_usize;
    let Some((first_start, first_end)) = next_token(underiv, &mut cursor, delimiter)? else {
        return Ok(0);
    };
    let second = next_token(underiv, &mut cursor, delimiter)?;
    let first_hyphen = underiv[first_start..first_end]
        .iter()
        .position(|&byte| byte == b'-' as i8)
        .map(|offset| first_start + offset);
    let second_hyphen = second.and_then(|(start, end)| {
        underiv[start..end]
            .iter()
            .position(|&byte| byte == b'-' as i8)
            .map(|offset| start + offset)
    });

    if let (Some((second_start, second_end)), Some(first_minus), Some(second_minus)) =
        (second, first_hyphen, second_hyphen)
    {
        let first_prefix_length = first_minus - first_start;
        let second_prefix_length = second_minus - second_start;
        if first_prefix_length == second_prefix_length
            && underiv[first_start..first_minus] == underiv[second_start..second_minus]
        {
            underiv[first_minus] = 0;
            underiv[second_minus] = 0;
            add(
                &mut list,
                list_capacity,
                &underiv[first_start..=first_minus],
                delimiter,
            )?;
            const HYPHEN: [i8; 2] = [b'-' as i8, 0];
            add(&mut list, list_capacity, &HYPHEN, 0)?;
            let first_suffix = &underiv[first_minus + 1..=first_end];
            let second_suffix = &underiv[second_minus + 1..=second_end];
            let diff = underiv_compare(first_suffix, second_suffix)?;
            let (smaller, greater) = if diff > 0 {
                (second_suffix, first_suffix)
            } else {
                (first_suffix, second_suffix)
            };
            if diff != 0 {
                add(&mut list, list_capacity, smaller, 0)?;
                const PLUS: [i8; 2] = [b'+' as i8, 0];
                add(&mut list, list_capacity, &PLUS, 0)?;
                add(&mut list, list_capacity, greater, 0)?;
            } else {
                const TWO: [i8; 2] = [b'2' as i8, 0];
                add(&mut list, list_capacity, &TWO, 0)?;
                add(&mut list, list_capacity, first_suffix, 0)?;
            }
        } else {
            add(
                &mut list,
                list_capacity,
                &underiv[first_start..=first_end],
                delimiter,
            )?;
            add(
                &mut list,
                list_capacity,
                &underiv[second_start..=second_end],
                delimiter,
            )?;
        }
    } else if let Some((second_start, second_end)) = second {
        add(
            &mut list,
            list_capacity,
            &underiv[first_start..=first_end],
            delimiter,
        )?;
        add(
            &mut list,
            list_capacity,
            &underiv[second_start..=second_end],
            delimiter,
        )?;
    } else {
        add(
            &mut list,
            list_capacity,
            &underiv[first_start..=first_end],
            delimiter,
        )?;
    }

    Ok(0)
}

pub(crate) fn sort_merge_underiv(
    sdf_value: &mut [i8],
    output_sdf: i32,
    underiv_list: &mut [i8],
    derivative_separator: i8,
    prefix: &[i8],
    postfix: &[i8],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5281 sort_merge_underiv
    // INCHI✔️❌: int sort_merge_underiv( char *pSdfValue,
    // INCHI✔️❌:                         int bOutputSdf,
    // INCHI✔️❌:                         char *szUnderivList,
    // INCHI✔️❌:                         char cDerivSeparator,
    // INCHI✔️❌:                         const char *pszUnderivPrefix,
    // INCHI✔️❌:                         const char *pszUnderivPostfix )
    // INCHI✔️❌: {
    // INCHI✔️❌: #define UNDERIV_MAX_NUM  512   /*max. number of records in szUnderivList */
    // INCHI✔️❌:     int num, numUnderiv = 0, i, j, k, m, n;
    // INCHI✔️❌:     char *q;
    // INCHI✔️❌:     char coeff[16];
    // INCHI✔️❌:     char *pszUnderiv[UNDERIV_MAX_NUM];
    // INCHI✔️❌:
    // INCHI✔️❌:     num = strlen( pSdfValue );
    // INCHI✔️❌:     n = num + strlen( pszUnderivPrefix ) + 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (n + strlen( pszUnderivPostfix ) + 1 + strlen( szUnderivList ) < MAX_SDF_VALUE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (numUnderiv = 0, q = strtok( szUnderivList, " " ); numUnderiv < UNDERIV_MAX_NUM && q; q = strtok( NULL, " " ), numUnderiv++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pszUnderiv[numUnderiv] = q;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*if ( !bOutputSdf || num ) {*/
    // INCHI✔️❌:         n = underiv_list_add( pSdfValue, MAX_SDF_VALUE, pszUnderivPrefix, 0 );
    // INCHI✔️❌:         /*}*/
    // INCHI✔️❌:         if (numUnderiv > 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             qsort( pszUnderiv, numUnderiv, sizeof( pszUnderiv[0] ), underiv_compare );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (i = 0; i < numUnderiv; i = j)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = i + 1; j < numUnderiv && !underiv_compare( pszUnderiv + i, pszUnderiv + j ); j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ; /* find identical derivatives */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             k = strlen( pszUnderiv[i] );
    // INCHI✔️❌:             if (1 < ( num = j - i ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 k = sprintf(coeff, "%d", num);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((long long)n + (long long)k + sizeof( pszUnderivPostfix ) < MAX_SDF_VALUE) /* djb-rwth: cast operators added */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 m = i ? cDerivSeparator : 0;
    // INCHI✔️❌:                 if (num > 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     underiv_list_add( pSdfValue, MAX_SDF_VALUE, coeff, m );
    // INCHI✔️❌:                     m = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 n = underiv_list_add( pSdfValue, MAX_SDF_VALUE, pszUnderiv[i], m );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 underiv_list_add( pSdfValue, MAX_SDF_VALUE, "!", 0 ); /* overflow indicator */
    // INCHI✔️❌:                 numUnderiv = -( 1 + numUnderiv );
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!bOutputSdf || num)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             underiv_list_add( pSdfValue, MAX_SDF_VALUE, pszUnderivPostfix, 0 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return numUnderiv;
    // INCHI✔️❌: #undef UNDERIV_MAX_NUM
    // INCHI✔️❌: }
    // INCHI✔️❌: #endif  /* UNDERIVATIZE_REPORT == 1 */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // END INCHI C FUNCTION: sort_merge_underiv
    fn c_len(value: &[i8]) -> Result<usize, SourceHeapError> {
        value
            .iter()
            .position(|&byte| byte == 0)
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    let initial_length = c_len(sdf_value)?;
    let prefix_length = c_len(prefix)?;
    let postfix_length = c_len(postfix)?;
    let list_length = c_len(underiv_list)?;
    let mut reused_num = initial_length as i32;
    let initial_n = initial_length
        .checked_add(prefix_length)
        .and_then(|value| value.checked_add(1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if initial_n + postfix_length + 1 + list_length >= MAX_SDF_VALUE as usize {
        return Ok(0);
    }

    let mut tokens = Vec::<Vec<i8>>::new();
    let mut cursor = 0_usize;
    while tokens.len() < 512 {
        while cursor < underiv_list.len() && underiv_list[cursor] == b' ' as i8 {
            cursor += 1;
        }
        let byte = *underiv_list
            .get(cursor)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if byte == 0 {
            break;
        }
        let start = cursor;
        while cursor < underiv_list.len()
            && underiv_list[cursor] != 0
            && underiv_list[cursor] != b' ' as i8
        {
            cursor += 1;
        }
        let end_byte = *underiv_list
            .get(cursor)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut token = underiv_list[start..cursor].to_vec();
        token.push(0);
        tokens.push(token);
        if end_byte == 0 {
            break;
        }
        underiv_list[cursor] = 0;
        cursor += 1;
    }

    let mut n = underiv_list_add(Some(sdf_value), MAX_SDF_VALUE as i32, Some(prefix), 0)?;
    tokens.sort_by(|left, right| {
        let result = underiv_compare(left, right).expect("owned tokens are NUL-terminated");
        result.cmp(&0)
    });
    let token_count = tokens.len() as i32;
    let mut index = 0_usize;
    while index < tokens.len() {
        let mut next = index + 1;
        while next < tokens.len() && underiv_compare(&tokens[index], &tokens[next])? == 0 {
            next += 1;
        }
        reused_num = (next - index) as i32;
        let coefficient = reused_num.to_string();
        let item_length = if reused_num > 1 {
            coefficient.len()
        } else {
            c_len(&tokens[index])?
        };
        if i64::from(n) + item_length as i64 + 8 < i64::from(MAX_SDF_VALUE) {
            let mut delimiter = if index != 0 { derivative_separator } else { 0 };
            if reused_num > 1 {
                let mut coefficient_c: Vec<i8> =
                    coefficient.bytes().map(|byte| byte as i8).collect();
                coefficient_c.push(0);
                underiv_list_add(
                    Some(sdf_value),
                    MAX_SDF_VALUE as i32,
                    Some(&coefficient_c),
                    delimiter,
                )?;
                delimiter = 0;
            }
            n = underiv_list_add(
                Some(sdf_value),
                MAX_SDF_VALUE as i32,
                Some(&tokens[index]),
                delimiter,
            )?;
        } else {
            underiv_list_add(
                Some(sdf_value),
                MAX_SDF_VALUE as i32,
                Some(&[b'!' as i8, 0]),
                0,
            )?;
            return Ok(-(1 + token_count));
        }
        index = next;
    }
    if output_sdf == 0 || reused_num != 0 {
        underiv_list_add(Some(sdf_value), MAX_SDF_VALUE as i32, Some(postfix), 0)?;
    }
    Ok(token_count)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn eliminate_deriv_not_in_list(
    atoms: &[inp_ATOM],
    derivatives: &mut [DERIV_AT],
    num_atoms: i32,
    mut underiv_list: Option<&mut [i8]>,
    len_underiv_list: i32,
    mut underiv_list2: Option<&mut [i8]>,
    len_underiv_list2: i32,
    mut bit_underiv_list: Option<&mut BIT_UNDERIV>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5351 eliminate_deriv_not_in_list
    // INCHI✔️❌: int eliminate_deriv_not_in_list( inp_ATOM *at,
    // INCHI✔️❌:                                  DERIV_AT *da,
    // INCHI✔️❌:                                  int num_atoms,
    // INCHI✔️❌:                                  char *szUnderivList,
    // INCHI✔️❌:                                  int lenUnderivList,
    // INCHI✔️❌:                                  char *szUnderivList2,
    // INCHI✔️❌:                                  int lenUnderivList2,
    // INCHI✔️❌:                                  BIT_UNDERIV *bitUnderivList )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, num_da, num_cuts = 0, num_cuts_this_atom, ret = 0;
    // INCHI✔️❌: #if( UNDERIVATIZE_REPORT == 1 )
    // INCHI✔️❌: #define UNDERIV_LEN 128
    // INCHI✔️❌: #define UNDERIV_LEN2 128
    // INCHI✔️❌:     char szUnderiv[UNDERIV_LEN];
    // INCHI✔️❌:     char szUnderiv2[UNDERIV_LEN2];
    // INCHI✔️❌:     BIT_UNDERIV  bitUnderiv;
    // INCHI✔️❌: #else
    // INCHI✔️❌: #define UNDERIV_LEN 0
    // INCHI✔️❌: #define UNDERIV_LEN2 0
    // INCHI✔️❌:     char *szUnderiv = NULL;
    // INCHI✔️❌:     char *szUnderiv2 = NULL;
    // INCHI✔️❌:     BIT_UNDERIV  bitUnderiv;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!da[i].typ[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* count deriative attachments */
    // INCHI✔️❌:         for (num_da = 0; num_da < DERIV_AT_LEN && da[i].typ[num_da]; num_da++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (num_da > 2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return -1; /* should not happen */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (num_da == 2 && da[i].typ[0] != da[i].typ[1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             da[i].typ[0] = da[i].typ[1] = 0; /* do not allow */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         underiv_buf_clear( szUnderiv );
    // INCHI✔️❌:         underiv_buf_clear( szUnderiv2 );
    // INCHI✔️❌:         bitUnderiv = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         num_cuts_this_atom = 0;
    // INCHI✔️❌:         if (da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = 0;
    // INCHI✔️❌: #ifndef UNDERIV_SYLYL_ONLY
    // INCHI✔️❌:             if (num_da == 2 && 1 + da[i].num[0] + da[i].num[1] <= MAX_AT_DERIV &&
    // INCHI✔️❌:                  0 < ( ret = is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( at, i, num_atoms, da + i, 0, szUnderiv, UNDERIV_LEN, szUnderiv2, UNDERIV_LEN2, &bitUnderiv ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_cuts += 2;
    // INCHI✔️❌:                 underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );
    // INCHI✔️❌:                 underiv_list_add( szUnderivList2, lenUnderivList2, szUnderiv2, ' ' );
    // INCHI✔️❌:                 *bitUnderivList |= bitUnderiv;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da[i].typ[0] = da[i].typ[1] = 0; /* not a derivative */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = 0;
    // INCHI✔️❌:             /*if ( da[i].num[0] <= MAX_AT_DERIV && 0 < (ret = is_deriv_chain( at, i, num_atoms, da+i, 0 )) )*/
    // INCHI✔️❌:             if (IS_DA_NUM_LE( da + i, 0, MAX_AT_DERIV ) && 0 < ( ret = is_deriv_chain( at, i, num_atoms, da + i, 0, szUnderiv, UNDERIV_LEN, szUnderiv2, UNDERIV_LEN2, &bitUnderiv ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_cuts++;
    // INCHI✔️❌:                 num_cuts_this_atom++;
    // INCHI✔️❌:                 j = 1;
    // INCHI✔️❌:                 /*underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );*/
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da[i].ord[0] = da[i].ord[1];
    // INCHI✔️❌:                     da[i].num[0] = da[i].num[1];
    // INCHI✔️❌:                     da[i].typ[0] = da[i].typ[1];
    // INCHI✔️❌:                     da[i].typ[1] = 0;
    // INCHI✔️❌:                     j = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /*underiv_buf_clear( szUnderiv );*/
    // INCHI✔️❌:             /*if ( da[i].typ[j] && da[i].num[j] <= MAX_AT_DERIV &&*/
    // INCHI✔️❌:             if (IS_DA_NUM_LE( da + i, j, MAX_AT_DERIV ) &&
    // INCHI✔️❌:                  0 < ( ret = is_deriv_chain( at, i, num_atoms, da + i, j, szUnderiv, UNDERIV_LEN, szUnderiv2, UNDERIV_LEN2, &bitUnderiv ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_cuts++;
    // INCHI✔️❌:                 num_cuts_this_atom++;
    // INCHI✔️❌:                 /*underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );*/
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da[i].typ[j] = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num_cuts_this_atom == 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* may happen only in RN- case, DERIV_AMINE_tN, but not always  */
    // INCHI✔️❌:                 underiv_list_add_two_cuts( szUnderivList, lenUnderivList, szUnderiv, ' ' );
    // INCHI✔️❌:                 underiv_list_add( szUnderivList2, lenUnderivList2, szUnderiv2, ' ' );
    // INCHI✔️❌:                 *bitUnderivList |= bitUnderiv;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (num_cuts_this_atom == 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );
    // INCHI✔️❌:                     underiv_list_add( szUnderivList2, lenUnderivList2, szUnderiv2, ' ' );
    // INCHI✔️❌:                     *bitUnderivList |= bitUnderiv;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_cuts;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: eliminate_deriv_not_in_list

    fn is_da_num_le(derivative: &DERIV_AT, index: usize) -> Result<bool, SourceHeapError> {
        let type_ = i32::from(
            *derivative
                .typ
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let number = i32::from(
            *derivative
                .num
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        Ok((type_ != 0 && (type_ & DERIV_UNEXPADABLE) == type_) || number <= MAX_AT_DERIV as i32)
    }

    fn add(
        list: &mut Option<&mut [i8]>,
        capacity: i32,
        value: &[i8],
    ) -> Result<(), SourceHeapError> {
        let _ = underiv_list_add(list.as_deref_mut(), capacity, Some(value), b' ' as i8)?;
        Ok(())
    }

    fn merge_bits(
        target: &mut Option<&mut BIT_UNDERIV>,
        value: BIT_UNDERIV,
    ) -> Result<(), SourceHeapError> {
        let target = target
            .as_deref_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *target |= value;
        Ok(())
    }

    let mut num_cuts = 0_i32;
    let mut atom_index = 0_i32;
    while atom_index < num_atoms {
        let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let derivative = derivatives
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if derivative.typ[0] == 0 {
            atom_index += 1;
            continue;
        }
        let num_da = derivative
            .typ
            .iter()
            .position(|&type_| type_ == 0)
            .unwrap_or(derivative.typ.len());
        if num_da > 2 {
            return Ok(-1);
        }
        if num_da == 2 && derivative.typ[0] != derivative.typ[1] {
            let derivative = derivatives
                .get_mut(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            derivative.typ[0] = 0;
            derivative.typ[1] = 0;
            atom_index += 1;
            continue;
        }

        let mut underiv = [0_i8; 128];
        let mut underiv2 = [0_i8; 128];
        underiv_buf_clear(Some(&mut underiv))?;
        underiv_buf_clear(Some(&mut underiv2))?;
        let mut bit_underiv = 0_i32;
        let mut num_cuts_this_atom = 0_i32;

        if derivatives[index].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR as i16 != 0 {
            let mut ret = 0_i32;
            if num_da == 2
                && 1 + i32::from(derivatives[index].num[0]) + i32::from(derivatives[index].num[1])
                    <= MAX_AT_DERIV as i32
            {
                ret = is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                    atoms,
                    atom_index,
                    num_atoms,
                    &derivatives[index],
                    0,
                    Some(&mut underiv),
                    128,
                    Some(&mut underiv2),
                    128,
                    Some(&mut bit_underiv),
                )?;
            }
            if ret > 0 {
                num_cuts += 2;
                add(&mut underiv_list, len_underiv_list, &underiv)?;
                add(&mut underiv_list2, len_underiv_list2, &underiv2)?;
                merge_bits(&mut bit_underiv_list, bit_underiv)?;
            } else if ret < 0 {
                return Ok(ret);
            } else {
                derivatives[index].typ[0] = 0;
                derivatives[index].typ[1] = 0;
            }
        } else {
            let mut ret = 0_i32;
            let derivative0_eligible = is_da_num_le(&derivatives[index], 0)?;
            if derivative0_eligible {
                ret = is_deriv_chain(
                    atoms,
                    atom_index,
                    num_atoms,
                    &derivatives[index],
                    0,
                    Some(&mut underiv),
                    128,
                    Some(&mut underiv2),
                    128,
                    Some(&mut bit_underiv),
                )?;
            }
            let derivative_index = if derivative0_eligible && ret > 0 {
                num_cuts += 1;
                num_cuts_this_atom += 1;
                1_usize
            } else if ret < 0 {
                return Ok(ret);
            } else {
                derivatives[index].ord[0] = derivatives[index].ord[1];
                derivatives[index].num[0] = derivatives[index].num[1];
                derivatives[index].typ[0] = derivatives[index].typ[1];
                derivatives[index].typ[1] = 0;
                0_usize
            };

            ret = 0;
            let derivative1_eligible = is_da_num_le(&derivatives[index], derivative_index)?;
            if derivative1_eligible {
                ret = is_deriv_chain(
                    atoms,
                    atom_index,
                    num_atoms,
                    &derivatives[index],
                    derivative_index as i32,
                    Some(&mut underiv),
                    128,
                    Some(&mut underiv2),
                    128,
                    Some(&mut bit_underiv),
                )?;
            }
            if derivative1_eligible && ret > 0 {
                num_cuts += 1;
                num_cuts_this_atom += 1;
            } else if ret < 0 {
                return Ok(ret);
            } else {
                derivatives[index].typ[derivative_index] = 0;
            }

            if num_cuts_this_atom == 2 {
                underiv_list_add_two_cuts(
                    underiv_list.as_deref_mut(),
                    len_underiv_list,
                    &mut underiv,
                    b' ' as i8,
                )?;
                add(&mut underiv_list2, len_underiv_list2, &underiv2)?;
                merge_bits(&mut bit_underiv_list, bit_underiv)?;
            } else if num_cuts_this_atom == 1 {
                add(&mut underiv_list, len_underiv_list, &underiv)?;
                add(&mut underiv_list2, len_underiv_list2, &underiv2)?;
                merge_bits(&mut bit_underiv_list, bit_underiv)?;
            }
        }
        atom_index += 1;
    }

    Ok(num_cuts)
}

pub(crate) fn make_single_cut(
    atoms: &mut [inp_ATOM],
    derivatives: &mut [DERIV_AT],
    atom_index: i32,
    cut_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5495 make_single_cut
    // INCHI✔️❌: int make_single_cut( inp_ATOM *at, DERIV_AT *da, int iat, int icut )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = -1; /* error flag */
    // INCHI✔️❌:     int iord = (int) da[iat].ord[icut]; /* ord of the bond in iat */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (da[iat].typ[icut] & DERIV_DUPLIC)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (iord < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* find other da[] that affect at[iat] */
    // INCHI✔️❌:         int jat = at[iat].neighbor[iord];  /* opposite atom */
    // INCHI✔️❌:         AT_NUMB *p = is_in_the_list( at[jat].neighbor, (AT_NUMB) iat, at[jat].valence );
    // INCHI✔️❌:         int jord = p ? ( p - at[jat].neighbor ) : -1;
    // INCHI✔️❌:         int i;
    // INCHI✔️❌:         const int iT = 2;
    // INCHI✔️❌: #ifdef UNDERIV_ADD_D_TO_PRECURSOR
    // INCHI✔️❌:         const int iD = 1;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         if (jord < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return -1;  /* program error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = DisconnectInpAtBond( at, NULL, iat, iord );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ret == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (da[iat].typ[icut] & DERIV_RING_OUTSIDE_PRECURSOR)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* at[jat] belongs to the main structure */
    // INCHI✔️❌:                 at[jat].num_H++;
    // INCHI✔️❌: #ifdef UNDERIV_ADD_D_TO_PRECURSOR
    // INCHI✔️❌:                 at[jat].num_iso_H[iD] ++; /* add D to the main structure */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 at[iat].num_H++;        /* add T to the derivatizing fragment */
    // INCHI✔️❌:                 at[iat].num_iso_H[iT] ++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (da[iat].typ[icut] && ( da[iat].typ[icut] & DERIV_REPL_N_WITH_O ) == da[iat].typ[icut])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at[jat].num_H++;
    // INCHI✔️❌:                     at[jat].num_iso_H[iT] ++; /* add T to the derivatizing fragment ??? */
    // INCHI✔️❌:                                               /* replace R=N-DerivAgent with R=O H-DerivAgent */
    // INCHI✔️❌:                     at[iat].elname[0] = 'O';
    // INCHI✔️❌:                     at[iat].el_number = EL_NUMBER_O;   /* since N replaced with O, do not add H */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (!icut && da[iat].typ[icut] && ( da[iat].typ[icut] & DERIV_REPL_N_WITH_OH ) == da[iat].typ[icut])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* cut #0 is the second; in case of DERIV_RING2_OUTSIDE_PRECUR on the 1st cut H has already been added */
    // INCHI✔️❌:                         at[jat].num_H++;
    // INCHI✔️❌:                         at[jat].num_iso_H[iT] ++; /* add T to the derivatizing fragment ??? */
    // INCHI✔️❌:                                                   /* replace R=N-DerivAgent with R=O H-DerivAgent */
    // INCHI✔️❌:                         at[iat].elname[0] = 'O';
    // INCHI✔️❌:                         at[iat].el_number = EL_NUMBER_O;   /* since N replaced with O, do not add H */
    // INCHI✔️❌: #ifdef DERIV_RING2_OUTSIDE_PRECUR
    // INCHI✔️❌:                         if (!( da[iat].typ[icut] & DERIV_RING2_OUTSIDE_PRECUR ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             at[iat].num_H++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at[iat].num_H++;
    // INCHI✔️❌: #ifdef UNDERIV_ADD_D_TO_PRECURSOR
    // INCHI✔️❌:                         at[iat].num_iso_H[iD] ++; /* add D to the main structure */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         at[jat].num_H++;        /* add T to the derivatizing fragment */
    // INCHI✔️❌:                         at[jat].num_iso_H[iT] ++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* adjust ord for other bonds */
    // INCHI✔️❌:             for (i = 0; i < DERIV_AT_LEN && da[iat].typ[i]; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (da[iat].ord[i] == iord)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da[iat].ord[i] = -( 1 + da[iat].ord[i] ); /* mark the bond being disconnected */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (da[iat].ord[i] > iord)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da[iat].ord[i] --;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (i = 0; i < DERIV_AT_LEN && da[jat].typ[i]; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (da[jat].ord[i] == jord)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* opposite atom needs the same bond to be disconnected */
    // INCHI✔️❌: #ifdef NEVER        /* not needed */
    // INCHI✔️❌:                     if (da[iat].num[icut] == da[jat].num[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iD = 2; /* mark both as fragments */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                         if (da[iat].num[icut] > da[jat].num[i])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             iD = 2; /* greater as a main structure */
    // INCHI✔️❌:                             iT = 1; /* mark smaller as a derivatizing fragment */
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     da[jat].ord[i] = -( 1 + da[jat].ord[i] );
    // INCHI✔️❌:                     da[jat].typ[i] |= DERIV_DUPLIC; /* do not cut again */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (da[jat].ord[i] > jord)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     da[jat].ord[i] --;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: make_single_cut

    let iat = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let cut = usize::try_from(cut_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let derivative = derivatives
        .get(iat)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let order = i32::from(
        *derivative
            .ord
            .get(cut)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let type_ = i32::from(
        *derivative
            .typ
            .get(cut)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    if type_ & DERIV_DUPLIC as i32 != 0 {
        return Ok(0);
    }
    if order < 0 {
        return Ok(-1);
    }
    let order_index = usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let opposite = usize::from(
        *atoms
            .get(iat)
            .and_then(|atom| atom.neighbor.get(order_index))
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let opposite_atom = atoms
        .get(opposite)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let reverse = is_in_the_list(
        Some(&opposite_atom.neighbor),
        atom_index as u16,
        i32::from(opposite_atom.valence),
    )?;
    let Some(reverse) = reverse else {
        return Ok(-1);
    };
    let result = DisconnectInpAtBond(atoms, None, atom_index, order)?;
    if result != 1 {
        return Ok(result);
    }

    if type_ & DERIV_RING_OUTSIDE_PRECURSOR as i32 != 0 {
        atoms[opposite].num_H = atoms[opposite].num_H.wrapping_add(1);
        atoms[iat].num_H = atoms[iat].num_H.wrapping_add(1);
        atoms[iat].num_iso_H[2] = atoms[iat].num_iso_H[2].wrapping_add(1);
    } else if type_ != 0 && (type_ & DERIV_REPL_N_WITH_O) == type_ {
        atoms[opposite].num_H = atoms[opposite].num_H.wrapping_add(1);
        atoms[opposite].num_iso_H[2] = atoms[opposite].num_iso_H[2].wrapping_add(1);
        atoms[iat].elname[0] = b'O' as i8;
        atoms[iat].el_number = EL_NUMBER_O;
    } else if cut == 0 && type_ != 0 && (type_ & DERIV_REPL_N_WITH_OH) == type_ {
        atoms[opposite].num_H = atoms[opposite].num_H.wrapping_add(1);
        atoms[opposite].num_iso_H[2] = atoms[opposite].num_iso_H[2].wrapping_add(1);
        atoms[iat].elname[0] = b'O' as i8;
        atoms[iat].el_number = EL_NUMBER_O;
        if type_ & DERIV_RING2_OUTSIDE_PRECUR as i32 == 0 {
            atoms[iat].num_H = atoms[iat].num_H.wrapping_add(1);
        }
    } else {
        atoms[iat].num_H = atoms[iat].num_H.wrapping_add(1);
        atoms[opposite].num_H = atoms[opposite].num_H.wrapping_add(1);
        atoms[opposite].num_iso_H[2] = atoms[opposite].num_iso_H[2].wrapping_add(1);
    }

    for index in 0..derivatives[iat].typ.len() {
        if derivatives[iat].typ[index] == 0 {
            break;
        }
        let value = i32::from(derivatives[iat].ord[index]);
        if value == order {
            derivatives[iat].ord[index] = (-(1 + value)) as i8;
        } else if value > order {
            derivatives[iat].ord[index] = derivatives[iat].ord[index].wrapping_sub(1);
        }
    }
    for index in 0..derivatives[opposite].typ.len() {
        if derivatives[opposite].typ[index] == 0 {
            break;
        }
        let value = i32::from(derivatives[opposite].ord[index]);
        if value == reverse as i32 {
            derivatives[opposite].ord[index] = (-(1 + value)) as i8;
            derivatives[opposite].typ[index] |= DERIV_DUPLIC as i16;
        } else if value > reverse as i32 {
            derivatives[opposite].ord[index] = derivatives[opposite].ord[index].wrapping_sub(1);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn oad_atom(element: u8, original_number: u16) -> inp_ATOM {
        let mut atom = inp_ATOM {
            el_number: element,
            orig_at_number: original_number,
            ..inp_ATOM::default()
        };
        atom.elname[0] = match element {
            EL_NUMBER_C => b'C',
            EL_NUMBER_H => b'H',
            EL_NUMBER_O => b'O',
            _ => b'N',
        } as i8;
        atom
    }

    fn oad_original(heap: &mut SourceHeap, disconnected: bool) -> ORIG_ATOM_DATA {
        let mut atoms = vec![oad_atom(EL_NUMBER_C, 1), oad_atom(EL_NUMBER_C, 2)];
        if !disconnected {
            atoms[0].valence = 1;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].neighbor[0] = 1;
            atoms[0].bond_type[0] = BOND_SINGLE as u8;
            atoms[1].valence = 1;
            atoms[1].chem_bonds_valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].bond_type[0] = BOND_SINGLE as u8;
        }
        ORIG_ATOM_DATA {
            at: heap.allocate_model_storage(atoms).unwrap(),
            num_inp_atoms: 2,
            num_inp_bonds: i32::from(!disconnected),
            ..ORIG_ATOM_DATA::default()
        }
    }

    #[test]
    fn source_port__ichinorm__oad_edit_underivatize__line_5934() {
        for (disconnected, output_sdf) in [(false, 0), (false, 1), (true, 0)] {
            let mut heap = SourceHeap::default();
            let mut original = oad_original(&mut heap, disconnected);
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let report = heap.allocate_model_storage(vec![0_i8; 128]).unwrap();
            assert_eq!(
                OAD_Edit_Underivatize(
                    &mut heap,
                    clock,
                    &mut CANON_GLOBALS::default(),
                    &mut original,
                    output_sdf,
                    1,
                    report,
                    0,
                ),
                Ok(0)
            );
            assert_eq!(original.num_inp_atoms, 2);
            assert_eq!(
                heap.slice(original.at.as_const()).unwrap()[..2]
                    .iter()
                    .map(|atom| atom.orig_at_number)
                    .collect::<Vec<_>>(),
                vec![1, 2]
            );
            assert_eq!(heap.slice(report.as_const()).unwrap()[0], 0);
            assert_eq!(
                heap.live_allocations_of::<INP_ATOM_DATA>(),
                0,
                "component descriptor allocation must be cleaned"
            );
            assert_eq!(heap.live_allocations_of::<DERIV_AT>(), 0);
            assert_eq!(heap.live_allocations_of::<R2C_ATPAIR>(), 0);
        }

        let normal = DERIV_AT {
            typ: [DERIV_BRIDGE_O as i16, 0, 0, 0],
            ..DERIV_AT::default()
        };
        let mut normal_atoms = [oad_atom(EL_NUMBER_C, 1), oad_atom(EL_NUMBER_C, 2)];
        normal_atoms[0].neighbor[0] = 1;
        assert_eq!(
            oad_prepare_derivative_cuts(&normal_atoms, &mut [normal, DERIV_AT::default()], 2,),
            Ok(Ok((1, 0, 1)))
        );

        let mut ring_atoms = [oad_atom(EL_NUMBER_C, 1), oad_atom(EL_NUMBER_C, 2)];
        ring_atoms[0].neighbor[0] = 1;
        ring_atoms[0].neighbor[4] = 1;
        ring_atoms[0].neighbor[7] = 1;
        ring_atoms[1].neighbor[0] = 0;
        let ring_type = DERIV_RING_OUTSIDE_PRECURSOR as i16;
        let ring = DERIV_AT {
            typ: [ring_type, ring_type, 0, 0],
            ..DERIV_AT::default()
        };
        assert_eq!(
            oad_prepare_derivative_cuts(&ring_atoms, &mut [ring, DERIV_AT::default()], 2),
            Ok(Ok((2, 2, 1)))
        );

        let ring2_type = DERIV_RING2_OUTSIDE_PRECUR as i16;
        let ring2 = DERIV_AT {
            typ: [ring2_type, ring2_type, 0, 0],
            ..DERIV_AT::default()
        };
        assert_eq!(
            oad_prepare_derivative_cuts(&ring_atoms, &mut [ring2, DERIV_AT::default()], 2),
            Ok(Ok((2, 2, 1)))
        );

        let amine = DERIV_AT {
            typ: [DERIV_AMINE_tN as i16, DERIV_AMINE_tN as i16, 0, 0],
            ..DERIV_AT::default()
        };
        assert_eq!(
            oad_prepare_derivative_cuts(&ring_atoms, &mut [amine, DERIV_AT::default()], 2),
            Ok(Ok((2, 0, 2)))
        );

        let equal_ro_cox = DERIV_AT {
            typ: [DERIV_RO_COX as i16, DERIV_RO_COX as i16, 0, 0],
            num: [3, 3, 0, 0],
            ..DERIV_AT::default()
        };
        let mut equal_case = [equal_ro_cox, DERIV_AT::default()];
        assert_eq!(
            oad_prepare_derivative_cuts(&ring_atoms, &mut equal_case, 2),
            Ok(Ok((0, 0, 0)))
        );
        assert_eq!(equal_case[0], DERIV_AT::default());

        let mixed = DERIV_AT {
            typ: [DERIV_BRIDGE_O as i16, DERIV_RO_COX as i16, 0, 0],
            ord: [4, 7, 0, 0],
            num: [1, 8, 0, 0],
            ..DERIV_AT::default()
        };
        let mut mixed_case = [mixed, DERIV_AT::default()];
        assert_eq!(
            oad_prepare_derivative_cuts(&ring_atoms, &mut mixed_case, 2),
            Ok(Ok((1, 0, 1)))
        );
        assert_eq!(
            (mixed_case[0].typ, mixed_case[0].ord, mixed_case[0].num),
            ([DERIV_RO_COX as i16, 0, 0, 0], [7, 0, 0, 0], [8, 0, 0, 0])
        );

        let invalid = DERIV_AT {
            typ: [1, 2, 3, 0],
            ..DERIV_AT::default()
        };
        assert_eq!(
            oad_prepare_derivative_cuts(&ring_atoms, &mut [invalid, DERIV_AT::default()], 2),
            Ok(Err(-88))
        );

        let mut allocation_heap = SourceHeap::default();
        let mut pairs = SourceMutPointer::null();
        let mut allocated = 0;
        allocation_heap.fail_after_allocations(0);
        assert_eq!(
            oad_allocate_atom_pairs(&mut allocation_heap, &mut pairs, &mut allocated, 2),
            Ok(false)
        );
        assert!(pairs.is_null());
        assert_eq!(allocation_heap.source_allocation_calls(), 1);
        assert_eq!(allocation_heap.live_allocation_count(), 0);
    }

    fn r2c_fixture() -> Vec<inp_ATOM> {
        let mut atoms = vec![inp_ATOM::default(); 6];
        atoms[0] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 4,
            chem_bonds_valence: 4,
            bCutVertex: 1,
            nRingSystem: 1,
            nNumAtInRingSystem: 5,
            ..inp_ATOM::default()
        };
        atoms[0].neighbor[..4].copy_from_slice(&[1, 2, 3, 4]);
        atoms[1] = inp_ATOM {
            el_number: EL_NUMBER_O,
            valence: 1,
            chem_bonds_valence: 1,
            num_H: 1,
            ..inp_ATOM::default()
        };
        atoms[1].neighbor[0] = 0;
        atoms[2] = inp_ATOM {
            el_number: EL_NUMBER_O,
            valence: 2,
            chem_bonds_valence: 2,
            nRingSystem: 1,
            ..inp_ATOM::default()
        };
        atoms[2].neighbor[..2].copy_from_slice(&[0, 3]);
        atoms[3] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 3,
            chem_bonds_valence: 3,
            num_H: 1,
            nRingSystem: 1,
            ..inp_ATOM::default()
        };
        atoms[3].neighbor[..3].copy_from_slice(&[0, 2, 5]);
        atoms[4] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 1,
            chem_bonds_valence: 1,
            nRingSystem: 2,
            ..inp_ATOM::default()
        };
        atoms[4].neighbor[0] = 0;
        atoms[5] = inp_ATOM {
            el_number: EL_NUMBER_O,
            valence: 1,
            chem_bonds_valence: 1,
            num_H: 1,
            ..inp_ATOM::default()
        };
        atoms[5].neighbor[0] = 3;
        atoms
    }

    #[test]
    fn source_port__ichinorm__detect_r2c_zatom__line_7159() {
        let untouched = R2C_AT {
            type_: -11,
            ordW: -12,
            ordY: -13,
            ordC: -14,
            ordCopt: -15,
        };
        let run = |atoms: &[inp_ATOM]| {
            let mut derivatives = vec![untouched.clone(); atoms.len()];
            let result = detect_r2c_Zatom(atoms, &mut derivatives, 0).unwrap();
            (result, derivatives[0].clone())
        };

        let atoms = r2c_fixture();
        assert_eq!(
            run(&atoms),
            (
                1,
                R2C_AT {
                    type_: 1,
                    ordW: 1,
                    ordY: 0,
                    ordC: 2,
                    ordCopt: 3,
                }
            )
        );

        let mut no_optional = atoms.clone();
        no_optional[0].valence = 3;
        no_optional[0].chem_bonds_valence = 3;
        assert_eq!(
            run(&no_optional),
            (
                1,
                R2C_AT {
                    type_: 1,
                    ordW: 1,
                    ordY: 0,
                    ordC: 2,
                    ordCopt: R2C_EMPTY as i8,
                }
            )
        );

        let mut opposite_first = atoms.clone();
        opposite_first[2].neighbor[..2].copy_from_slice(&[3, 0]);
        assert_eq!(run(&opposite_first).0, 1);

        let assert_rejected = |candidate: Vec<inp_ATOM>| {
            assert_eq!(run(&candidate), (0, untouched.clone()));
        };
        for candidate in [
            {
                let mut value = atoms.clone();
                value[0].valence = 5;
                value
            },
            {
                let mut value = atoms.clone();
                value[0].chem_bonds_valence = 3;
                value
            },
            {
                let mut value = atoms.clone();
                value[0].el_number = EL_NUMBER_N;
                value
            },
            {
                let mut value = atoms.clone();
                value[0].nNumAtInRingSystem = 4;
                value
            },
            {
                let mut value = atoms.clone();
                value[0].bCutVertex = 0;
                value
            },
            {
                let mut value = atoms.clone();
                value[0].valence = -1;
                value[0].chem_bonds_valence = -1;
                value
            },
            {
                let mut value = atoms.clone();
                value[1].charge = 1;
                value
            },
            {
                let mut value = atoms.clone();
                value[1].radical = 1;
                value
            },
            {
                let mut value = atoms.clone();
                value[4] = value[1].clone();
                value[4].neighbor[0] = 0;
                value
            },
            {
                let mut value = atoms.clone();
                value[3].chem_bonds_valence = 2;
                value
            },
            {
                let mut value = atoms.clone();
                value[3].el_number = EL_NUMBER_N;
                value
            },
            {
                let mut value = atoms.clone();
                value[4] = value[2].clone();
                value[4].neighbor[..2].copy_from_slice(&[0, 3]);
                value
            },
            {
                let mut value = atoms.clone();
                value[5].charge = 1;
                value
            },
            {
                let mut value = atoms.clone();
                value[5].radical = 1;
                value
            },
            {
                let mut value = atoms.clone();
                value[5].num_H = 0;
                value
            },
            {
                let mut value = atoms.clone();
                value[4] = value[3].clone();
                value[4].neighbor[..3].copy_from_slice(&[0, 2, 5]);
                value
            },
            {
                let mut value = atoms.clone();
                value.push(value[4].clone());
                value[0].neighbor[..4].copy_from_slice(&[1, 2, 4, 6]);
                value
            },
            {
                let mut value = atoms.clone();
                value[4].el_number = EL_NUMBER_N;
                value
            },
            {
                let mut value = atoms.clone();
                value[0].valence = 2;
                value[0].chem_bonds_valence = 2;
                value
            },
        ] {
            assert_rejected(candidate);
        }

        let mut missing_neighbor = atoms.clone();
        missing_neighbor[0].neighbor[0] = u16::MAX;
        assert_eq!(
            detect_r2c_Zatom(
                &missing_neighbor,
                &mut vec![R2C_AT::default(); missing_neighbor.len()],
                0,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            detect_r2c_Zatom(&atoms, &mut [], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            detect_r2c_Zatom(&atoms, &mut [], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    fn r2c_cut_fixture() -> (Vec<inp_ATOM>, Vec<R2C_AT>) {
        let mut atoms = r2c_fixture();
        atoms[0].bond_type[..4].fill(BOND_SINGLE as u8);
        atoms[1].bond_type[0] = BOND_SINGLE as u8;
        atoms[2].bond_type[..2].fill(BOND_SINGLE as u8);
        atoms[3].bond_type[..3].fill(BOND_SINGLE as u8);
        atoms[4].bond_type[0] = BOND_SINGLE as u8;
        atoms[5].bond_type[0] = BOND_SINGLE as u8;
        let mut derivatives = vec![R2C_AT::default(); atoms.len()];
        assert_eq!(detect_r2c_Zatom(&atoms, &mut derivatives, 0), Ok(1));
        (atoms, derivatives)
    }

    #[test]
    fn source_port__ichinorm__cut_ring_to_chain__line_7296() {
        let (baseline, derivatives) = r2c_cut_fixture();

        let mut atoms = baseline.clone();
        assert_eq!(cut_ring_to_chain(&mut atoms, &derivatives, 0), Ok(1));
        let mut expected = baseline.clone();
        expected[0].neighbor[..4].copy_from_slice(&[1, 3, 4, 0]);
        expected[0].bond_type[..4].copy_from_slice(&[2, 1, 1, 0]);
        expected[0].valence = 3;
        expected[0].chem_bonds_valence = 4;
        expected[1].bond_type[0] = 2;
        expected[1].chem_bonds_valence = 2;
        expected[1].num_H = 0;
        expected[2].neighbor[..2].copy_from_slice(&[3, 0]);
        expected[2].bond_type[..2].copy_from_slice(&[1, 0]);
        expected[2].valence = 1;
        expected[2].chem_bonds_valence = 1;
        expected[2].num_H = 1;
        assert_eq!(atoms, expected);

        let mut inactive = baseline.clone();
        let mut inactive_derivatives = derivatives.clone();
        inactive_derivatives[0].type_ = 0;
        assert_eq!(
            cut_ring_to_chain(&mut inactive, &inactive_derivatives, 0),
            Ok(0)
        );
        assert_eq!(inactive, baseline);

        for (field, value) in [(0, -1), (0, 4), (1, -1), (1, 4), (2, -1), (2, 4)] {
            let mut atoms = baseline.clone();
            let mut invalid = derivatives.clone();
            match field {
                0 => invalid[0].ordW = value,
                1 => invalid[0].ordY = value,
                _ => invalid[0].ordC = value,
            }
            assert_eq!(cut_ring_to_chain(&mut atoms, &invalid, 0), Ok(-1));
            assert_eq!(atoms, baseline);
        }

        let mut no_hydrogen = baseline.clone();
        no_hydrogen[1].num_H = 0;
        let before = no_hydrogen.clone();
        assert_eq!(cut_ring_to_chain(&mut no_hydrogen, &derivatives, 0), Ok(-2));
        assert_eq!(no_hydrogen, before);

        let mut wrong_bond = baseline.clone();
        wrong_bond[0].bond_type[0] = BOND_DOUBLE as u8;
        let before = wrong_bond.clone();
        assert_eq!(cut_ring_to_chain(&mut wrong_bond, &derivatives, 0), Ok(-2));
        assert_eq!(wrong_bond, before);

        let mut no_reverse_y = baseline.clone();
        no_reverse_y[1].neighbor[0] = 5;
        let before = no_reverse_y.clone();
        assert_eq!(
            cut_ring_to_chain(&mut no_reverse_y, &derivatives, 0),
            Ok(-3)
        );
        assert_eq!(no_reverse_y, before);

        let mut no_reverse_w = baseline.clone();
        no_reverse_w[2].neighbor[..2].copy_from_slice(&[3, 3]);
        let before = no_reverse_w.clone();
        assert_eq!(
            cut_ring_to_chain(&mut no_reverse_w, &derivatives, 0),
            Ok(-4)
        );
        let mut partial = before.clone();
        partial[0].bond_type[0] = 2;
        partial[0].chem_bonds_valence = 5;
        partial[1].bond_type[0] = 2;
        partial[1].chem_bonds_valence = 2;
        assert_eq!(no_reverse_w, partial);

        let mut isotopic = baseline.clone();
        isotopic[1].num_H = 3;
        isotopic[1].num_iso_H = [1, 1, 1];
        isotopic[2].num_iso_H = [2, 3, 4];
        assert_eq!(cut_ring_to_chain(&mut isotopic, &derivatives, 0), Ok(1));
        assert_eq!(isotopic[1].num_H, 2);
        assert_eq!(isotopic[1].num_iso_H, [0, 0, 0]);
        assert_eq!(isotopic[2].num_H, 1);
        assert_eq!(isotopic[2].num_iso_H, [3, 4, 5]);

        let mut nonmatching_isotopes = baseline.clone();
        nonmatching_isotopes[1].num_H = 2;
        nonmatching_isotopes[1].num_iso_H = [1, 0, 0];
        assert_eq!(
            cut_ring_to_chain(&mut nonmatching_isotopes, &derivatives, 0),
            Ok(1)
        );
        assert_eq!(nonmatching_isotopes[1].num_H, 1);
        assert_eq!(nonmatching_isotopes[1].num_iso_H, [1, 0, 0]);
        assert_eq!(nonmatching_isotopes[2].num_iso_H, [0, 0, 0]);

        let mut signed_boundary = baseline.clone();
        signed_boundary[1].num_H = i8::MIN;
        signed_boundary[2].num_H = i8::MAX;
        assert_eq!(
            cut_ring_to_chain(&mut signed_boundary, &derivatives, 0),
            Ok(1)
        );
        assert_eq!(signed_boundary[1].num_H, i8::MAX);
        assert_eq!(signed_boundary[2].num_H, i8::MIN);

        let mut atoms = baseline.clone();
        assert_eq!(
            cut_ring_to_chain(&mut atoms, &[], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            cut_ring_to_chain(&mut atoms, &derivatives, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    fn connect_r2c_atoms(atoms: &mut [inp_ATOM], first: usize, second: usize) {
        let first_ordinal = usize::try_from(atoms[first].valence).unwrap();
        let second_ordinal = usize::try_from(atoms[second].valence).unwrap();
        atoms[first].neighbor[first_ordinal] = second as u16;
        atoms[first].bond_type[first_ordinal] = BOND_SINGLE as u8;
        atoms[first].valence += 1;
        atoms[first].chem_bonds_valence += 1;
        atoms[second].neighbor[second_ordinal] = first as u16;
        atoms[second].bond_type[second_ordinal] = BOND_SINGLE as u8;
        atoms[second].valence += 1;
        atoms[second].chem_bonds_valence += 1;
    }

    fn r2c_ring_atoms(original_number_offset: u16) -> Vec<inp_ATOM> {
        let mut atoms = vec![inp_ATOM::default(); 8];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.orig_at_number = original_number_offset + index as u16 + 1;
        }
        for index in [0, 3, 4, 6, 7] {
            atoms[index].el_number = EL_NUMBER_C;
            atoms[index].elname[0] = b'C' as i8;
        }
        for index in [1, 2, 5] {
            atoms[index].el_number = EL_NUMBER_O;
            atoms[index].elname[0] = b'O' as i8;
        }
        atoms[1].num_H = 1;
        atoms[3].num_H = 1;
        atoms[4].num_H = 3;
        atoms[5].num_H = 1;
        atoms[6].num_H = 2;
        atoms[7].num_H = 2;
        connect_r2c_atoms(&mut atoms, 0, 1);
        connect_r2c_atoms(&mut atoms, 0, 2);
        connect_r2c_atoms(&mut atoms, 2, 6);
        connect_r2c_atoms(&mut atoms, 6, 7);
        connect_r2c_atoms(&mut atoms, 7, 3);
        connect_r2c_atoms(&mut atoms, 3, 0);
        connect_r2c_atoms(&mut atoms, 0, 4);
        connect_r2c_atoms(&mut atoms, 3, 5);
        atoms
    }

    fn r2c_original(heap: &mut SourceHeap, atoms: Vec<inp_ATOM>) -> ORIG_ATOM_DATA {
        let atom_count = atoms.len() as i32;
        let bond_count = atoms
            .iter()
            .map(|atom| i32::from(atom.valence))
            .sum::<i32>()
            / 2;
        ORIG_ATOM_DATA {
            at: heap.allocate_model_storage(atoms).unwrap(),
            num_inp_atoms: atom_count,
            num_inp_bonds: bond_count,
            ..ORIG_ATOM_DATA::default()
        }
    }

    fn r2c_bond(atoms: &[inp_ATOM], first_number: u16, second_number: u16) -> Option<u8> {
        let first = atoms
            .iter()
            .position(|atom| atom.orig_at_number == first_number)?;
        let second = atoms
            .iter()
            .position(|atom| atom.orig_at_number == second_number)?;
        let atom = &atoms[first];
        (0..usize::try_from(atom.valence).ok()?)
            .find(|ordinal| usize::from(atom.neighbor[*ordinal]) == second)
            .map(|ordinal| atom.bond_type[ordinal])
    }

    #[test]
    fn source_port__ichinorm__ring2chain__line_7362() {
        let mut heap = SourceHeap::default();
        let mut original = r2c_original(&mut heap, r2c_ring_atoms(0));
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        assert_eq!(
            Ring2Chain(
                &mut heap,
                clock,
                &mut CANON_GLOBALS::default(),
                &mut original,
                0,
            ),
            Ok(1)
        );
        assert_eq!(original.num_inp_atoms, 8);
        assert!(original.nCurAtLen.is_null());
        let atoms = &heap.slice(original.at.as_const()).unwrap()[..8];
        assert_eq!(r2c_bond(atoms, 1, 2), Some(BOND_DOUBLE as u8));
        assert_eq!(r2c_bond(atoms, 1, 3), None);
        let hydroxyl = atoms.iter().find(|atom| atom.orig_at_number == 2).unwrap();
        let opened_oxygen = atoms.iter().find(|atom| atom.orig_at_number == 3).unwrap();
        assert_eq!(hydroxyl.num_H, 0);
        assert_eq!(opened_oxygen.num_H, 1);
        assert_eq!(heap.live_allocations_of::<R2C_AT>(), 0);
        assert_eq!(heap.live_allocations_of::<R2C_ATPAIR>(), 0);
        assert_eq!(heap.live_allocations_of::<INP_ATOM_DATA>(), 0);

        let mut multi_heap = SourceHeap::default();
        let mut multi_atoms = r2c_ring_atoms(0);
        let mut second = r2c_ring_atoms(8);
        for atom in &mut second {
            let valence = usize::try_from(atom.valence).unwrap();
            for neighbor in &mut atom.neighbor[..valence] {
                *neighbor += 8;
            }
        }
        multi_atoms.extend(second);
        connect_r2c_atoms(&mut multi_atoms, 4, 12);
        let mut multi_original = r2c_original(&mut multi_heap, multi_atoms);
        let multi_clock = multi_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        assert_eq!(
            Ring2Chain(
                &mut multi_heap,
                multi_clock,
                &mut CANON_GLOBALS::default(),
                &mut multi_original,
                0,
            ),
            Ok(2)
        );
        let multi_result = &multi_heap.slice(multi_original.at.as_const()).unwrap()[..16];
        assert_eq!(r2c_bond(multi_result, 1, 3), None);
        assert_eq!(r2c_bond(multi_result, 9, 11), None);
        assert_eq!(r2c_bond(multi_result, 1, 2), Some(BOND_DOUBLE as u8));
        assert_eq!(r2c_bond(multi_result, 9, 10), Some(BOND_DOUBLE as u8));
        assert_eq!(multi_heap.live_allocations_of::<R2C_AT>(), 0);
        assert_eq!(multi_heap.live_allocations_of::<R2C_ATPAIR>(), 0);
        assert_eq!(multi_heap.live_allocations_of::<INP_ATOM_DATA>(), 0);

        let mut empty_heap = SourceHeap::default();
        let mut carbon = inp_ATOM {
            el_number: EL_NUMBER_C,
            num_H: 4,
            orig_at_number: 1,
            ..inp_ATOM::default()
        };
        carbon.elname[0] = b'C' as i8;
        let mut no_candidate = r2c_original(&mut empty_heap, vec![carbon]);
        let empty_clock = empty_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        assert_eq!(
            Ring2Chain(
                &mut empty_heap,
                empty_clock,
                &mut CANON_GLOBALS::default(),
                &mut no_candidate,
                0,
            ),
            Ok(0)
        );
        assert_eq!(no_candidate.num_inp_atoms, 1);
        assert!(!no_candidate.nCurAtLen.is_null());
        assert_eq!(empty_heap.live_allocations_of::<R2C_AT>(), 0);
        assert_eq!(empty_heap.live_allocations_of::<R2C_ATPAIR>(), 0);
        assert_eq!(empty_heap.live_allocations_of::<INP_ATOM_DATA>(), 0);
    }

    #[test]
    fn source_port__ichinorm__mark_arom_bonds__line_784() {
        fn two_carbons(heap: &mut SourceHeap) -> SourceMutPointer<inp_ATOM> {
            let mut atoms = vec![
                inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 1,
                    chem_bonds_valence: 1,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 1,
                    chem_bonds_valence: 1,
                    ..inp_ATOM::default()
                },
            ];
            atoms[0].neighbor[0] = 1;
            atoms[0].bond_type[0] = BOND_SINGLE as u8;
            atoms[1].neighbor[0] = 0;
            atoms[1].bond_type[0] = BOND_SINGLE as u8;
            heap.allocate_model_storage(atoms).unwrap()
        }

        let mut heap = SourceHeap::default();
        let atoms = two_carbons(&mut heap);
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let baseline = heap.live_allocation_count();
        assert_eq!(
            mark_arom_bonds(&mut heap, clock, &mut CANON_GLOBALS::default(), atoms, 2, 0,),
            Ok(2)
        );
        assert_eq!(heap.live_allocation_count(), baseline);

        let mut failure_heap = SourceHeap::default();
        let failure_atoms = two_carbons(&mut failure_heap);
        let failure_clock = failure_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let failure_baseline = failure_heap.live_allocation_count();
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            mark_arom_bonds(
                &mut failure_heap,
                failure_clock,
                &mut CANON_GLOBALS::default(),
                failure_atoms,
                2,
                0,
            ),
            Ok(crate::source_types::BNS_OUT_OF_RAM)
        );
        assert_eq!(failure_heap.live_allocation_count(), failure_baseline);
    }

    fn atom(neighbors: &[u16]) -> inp_ATOM {
        let mut atom = inp_ATOM {
            valence: i8::try_from(neighbors.len()).unwrap(),
            ..inp_ATOM::default()
        };
        atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
        atom
    }

    #[test]
    fn source_port__ichinorm__cmp_r2c_atpair__line_800() {
        let pair = |first, second, atno| R2C_ATPAIR {
            at: [first, second],
            atno,
        };

        assert_eq!(cmp_r2c_atpair(&pair(1, 50, 3), &pair(2, 0, 3)), -1);
        assert_eq!(cmp_r2c_atpair(&pair(2, 0, 3), &pair(1, 50, 3)), 1);
        assert_eq!(cmp_r2c_atpair(&pair(7, 8, 3), &pair(7, 10, 3)), -2);
        assert_eq!(cmp_r2c_atpair(&pair(7, 10, 3), &pair(7, 8, 3)), 2);
        assert_eq!(cmp_r2c_atpair(&pair(7, 8, 0), &pair(7, 8, u16::MAX)), 0);
        assert_eq!(
            cmp_r2c_atpair(&pair(0, u16::MAX, 0), &pair(u16::MAX, 0, 0)),
            -65_535
        );
        assert_eq!(
            cmp_r2c_atpair(&pair(u16::MAX, 0, 0), &pair(0, u16::MAX, 0)),
            65_535
        );
        assert_eq!(
            cmp_r2c_atpair(&pair(0, 0, 0), &pair(0, u16::MAX, 0)),
            -65_535
        );
    }

    #[test]
    fn source_port__ichinorm__has_atom_pair_seq__line_815() {
        let pairs = [
            R2C_ATPAIR {
                at: [1, 4],
                atno: 99,
            },
            R2C_ATPAIR {
                at: [2, 7],
                atno: 5,
            },
            R2C_ATPAIR {
                at: [2, 7],
                atno: u16::MAX,
            },
            R2C_ATPAIR {
                at: [0, u16::MAX],
                atno: 3,
            },
        ];

        assert_eq!(has_atom_pair_seq(&pairs, -1, 1, 4), Ok(0));
        assert_eq!(has_atom_pair_seq(&pairs, 0, 1, 4), Ok(0));
        assert_eq!(has_atom_pair_seq(&pairs, 4, 1, 4), Ok(1));
        assert_eq!(has_atom_pair_seq(&pairs, 4, 4, 1), Ok(1));
        assert_eq!(has_atom_pair_seq(&pairs, 4, 7, 2), Ok(2));
        assert_eq!(has_atom_pair_seq(&pairs, 4, 0, u16::MAX), Ok(4));
        assert_eq!(has_atom_pair_seq(&pairs, 4, 9, 9), Ok(0));
        assert_eq!(has_atom_pair_seq(&pairs, 1, 2, 7), Ok(0));
        assert_eq!(
            has_atom_pair_seq(&pairs, 5, 9, 9),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__mark_atoms_ap__line_871() {
        let blocked = [R2C_ATPAIR {
            at: [1, 2],
            atno: 77,
        }];
        let mut chain = [atom(&[1]), atom(&[0, 2]), atom(&[1])];
        assert_eq!(
            mark_atoms_ap(&mut chain, 0, &blocked, 1, 10, 0x1234),
            Ok(12)
        );
        assert_eq!(chain[0].at_type, 0x1234);
        assert_eq!(chain[1].at_type, 0x1234);
        assert_eq!(chain[2].at_type, 0);

        let mut premarked = [atom(&[1]), atom(&[0])];
        premarked[0].at_type = 9;
        assert_eq!(
            mark_atoms_ap(&mut premarked, 0, &[], i32::MAX, -7, 42),
            Ok(-7)
        );
        assert_eq!(premarked[1].at_type, 0);

        let mut cycle = [atom(&[1, 2]), atom(&[0, 2]), atom(&[0, 1])];
        assert_eq!(mark_atoms_ap(&mut cycle, 0, &[], 0, 0, u16::MAX), Ok(3));
        assert!(cycle.iter().all(|atom| atom.at_type == u16::MAX));

        let mut stopped = [atom(&[1]), atom(&[0, 2]), atom(&[1])];
        stopped[1].at_type = 88;
        assert_eq!(mark_atoms_ap(&mut stopped, 0, &[], 0, 5, 7), Ok(6));
        assert_eq!(stopped[0].at_type, 7);
        assert_eq!(stopped[1].at_type, 88);
        assert_eq!(stopped[2].at_type, 0);

        let mut negative_valence = [inp_ATOM {
            valence: -1,
            ..inp_ATOM::default()
        }];
        assert_eq!(
            mark_atoms_ap(&mut negative_valence, 0, &[], 0, i32::MAX, 1),
            Ok(i32::MIN)
        );
        assert_eq!(negative_valence[0].at_type, 1);

        assert_eq!(
            mark_atoms_ap(&mut [], 0, &[], 0, 0, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut malformed = [atom(&[1])];
        assert_eq!(
            mark_atoms_ap(&mut malformed, 0, &[], 0, 0, 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(malformed[0].at_type, 5);

        let mut short_pairs = [atom(&[0])];
        assert_eq!(
            mark_atoms_ap(&mut short_pairs, 0, &[], 1, 0, 6),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short_pairs[0].at_type, 6);
    }

    #[test]
    fn source_port__ichinorm__mark_deriv_agents__line_5743() {
        let pair = |first, second| R2C_ATPAIR {
            at: [first, second],
            atno: first,
        };

        let mut cleared = [atom(&[1]), atom(&[0])];
        cleared[0].at_type = 8;
        cleared[0].cFlags = 9;
        cleared[1].at_type = 10;
        cleared[1].cFlags = 11;
        let mut components = u16::MAX;
        let mut count = i32::MAX;
        assert_eq!(
            mark_deriv_agents(
                &mut cleared,
                &[DERIV_AT::default(), DERIV_AT::default()],
                2,
                &[],
                -1,
                &mut components,
                &mut count,
            ),
            Ok(0)
        );
        assert_eq!((components, count), (0, 0));
        assert!(
            cleared
                .iter()
                .all(|atom| atom.at_type == 0 && atom.cFlags == 0)
        );

        let mut ordinary_atoms = [atom(&[1]), atom(&[0])];
        let mut ordinary_derivatives = [DERIV_AT::default(), DERIV_AT::default()];
        ordinary_derivatives[0].typ[0] = DERIV_BRIDGE_O as i16;
        assert_eq!(
            mark_deriv_agents(
                &mut ordinary_atoms,
                &ordinary_derivatives,
                2,
                &[pair(0, 1)],
                1,
                &mut components,
                &mut count,
            ),
            Ok(0)
        );
        assert_eq!((components, count), (1, 1));
        assert_eq!(ordinary_atoms[0].at_type, 1);
        assert_eq!(ordinary_atoms[1].at_type, 0);

        let mut ring_atoms = [atom(&[1]), atom(&[0])];
        let mut ring_derivatives = [DERIV_AT::default(), DERIV_AT::default()];
        ring_derivatives[0].typ[0] = DERIV_RING_O_OUTSIDE_PRECURSOR as i16;
        assert_eq!(
            mark_deriv_agents(
                &mut ring_atoms,
                &ring_derivatives,
                2,
                &[pair(0, 1)],
                1,
                &mut components,
                &mut count,
            ),
            Ok(0)
        );
        assert_eq!((components, count), (1, 1));
        assert_eq!(ring_atoms[0].at_type, 0);
        assert_eq!(ring_atoms[1].at_type, 1);

        let mut repeated_atoms = [atom(&[1, 2]), atom(&[0]), atom(&[0])];
        let mut repeated_derivatives = [
            DERIV_AT::default(),
            DERIV_AT::default(),
            DERIV_AT::default(),
        ];
        repeated_derivatives[0].typ[0] = DERIV_BRIDGE_O as i16;
        let repeated_pairs = [pair(0, 1), pair(0, 2)];
        assert_eq!(
            mark_deriv_agents(
                &mut repeated_atoms,
                &repeated_derivatives,
                3,
                &repeated_pairs,
                2,
                &mut components,
                &mut count,
            ),
            Ok(0)
        );
        assert_eq!((components, count), (1, 1));
        assert_eq!(repeated_atoms[0].at_type, 1);

        let mut no_mark_atoms = [atom(&[1]), atom(&[0])];
        no_mark_atoms[0].at_type = 55;
        components = 7;
        count = 8;
        assert_eq!(
            mark_deriv_agents(
                &mut no_mark_atoms,
                &[DERIV_AT::default(), DERIV_AT::default()],
                2,
                &[pair(0, 1)],
                1,
                &mut components,
                &mut count,
            ),
            Ok(-3)
        );
        assert_eq!((components, count), (0, 0));
        assert_eq!(no_mark_atoms[0].at_type, 0);

        let mut both_marked = ordinary_derivatives.clone();
        both_marked[1].typ[0] = DERIV_BRIDGE_NH as i16;
        let mut both_marked_atoms = [atom(&[1]), atom(&[0])];
        assert_eq!(
            mark_deriv_agents(
                &mut both_marked_atoms,
                &both_marked,
                2,
                &[pair(0, 1)],
                1,
                &mut components,
                &mut count,
            ),
            Ok(-3)
        );
        assert_eq!((components, count), (0, 0));

        let mut late_error_atoms = [atom(&[1, 2]), atom(&[0]), atom(&[0])];
        let late_error_pairs = [pair(0, 1), pair(1, 2)];
        let late_error_derivatives = [
            ordinary_derivatives[0].clone(),
            DERIV_AT::default(),
            DERIV_AT::default(),
        ];
        components = 9;
        count = 10;
        assert_eq!(
            mark_deriv_agents(
                &mut late_error_atoms,
                &late_error_derivatives,
                3,
                &late_error_pairs,
                2,
                &mut components,
                &mut count,
            ),
            Ok(-3)
        );
        assert_eq!((components, count), (0, 0));
        assert_eq!(late_error_atoms[0].at_type, 1);

        let mut short_atoms = [inp_ATOM::default()];
        short_atoms[0].at_type = 4;
        components = 5;
        count = 6;
        assert_eq!(
            mark_deriv_agents(
                &mut short_atoms,
                &[DERIV_AT::default()],
                2,
                &[],
                0,
                &mut components,
                &mut count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((components, count), (0, 0));
        assert_eq!(short_atoms[0].at_type, 0);

        let mut short_pair_atoms = [inp_ATOM::default()];
        assert_eq!(
            mark_deriv_agents(
                &mut short_pair_atoms,
                &[DERIV_AT::default()],
                1,
                &[],
                1,
                &mut components,
                &mut count,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((components, count), (0, 0));
    }

    #[test]
    fn source_port__ichinorm__replace_arom_bonds__line_5809() {
        fn fixture(replacement: u8) -> ([inp_ATOM; 2], [inp_ATOM; 2]) {
            let mut current = [atom(&[1]), atom(&[0])];
            current[0].orig_at_number = 20;
            current[1].orig_at_number = 10;
            current[0].bond_type[0] = BOND_TYPE_ALTERN as u8;
            current[1].bond_type[0] = BOND_TYPE_ALTERN as u8;

            let mut original = [atom(&[1]), atom(&[0])];
            original[0].orig_at_number = 10;
            original[1].orig_at_number = 20;
            original[0].bond_type[0] = replacement;
            original[1].bond_type[0] = replacement;
            (current, original)
        }

        let (mut current, original) = fixture(BOND_DOUBLE as u8);
        assert_eq!(replace_arom_bonds(&mut current, 2, &original, 2), Ok(0));
        assert_eq!(current[0].bond_type[0], BOND_DOUBLE as u8);
        assert_eq!(current[1].bond_type[0], BOND_DOUBLE as u8);
        assert_eq!(current[0].orig_at_number, 20);
        assert_eq!(current[1].orig_at_number, 10);

        let (mut unchanged, original) = fixture(BOND_SINGLE as u8);
        unchanged[0].bond_type[0] = BOND_TRIPLE as u8;
        unchanged[1].bond_type[0] = BOND_TRIPLE as u8;
        assert_eq!(replace_arom_bonds(&mut unchanged, 2, &original, 2), Ok(0));
        assert_eq!(unchanged[0].bond_type[0], BOND_TRIPLE as u8);
        assert_eq!(unchanged[1].bond_type[0], BOND_TRIPLE as u8);

        let (mut repeated, aromatic_original) = fixture(BOND_TYPE_ALTERN as u8);
        assert_eq!(
            replace_arom_bonds(&mut repeated, 2, &aromatic_original, 2),
            Ok(0)
        );
        assert_eq!(repeated[0].bond_type[0], BOND_TYPE_ALTERN as u8);
        assert_eq!(repeated[1].bond_type[0], BOND_TYPE_ALTERN as u8);

        let (mut missing, original) = fixture(BOND_DOUBLE as u8);
        assert_eq!(
            replace_arom_bonds(&mut missing, 2, &original[..1], 1),
            Ok(2)
        );
        assert_eq!(missing[0].bond_type[0], BOND_TYPE_ALTERN as u8);
        assert_eq!(missing[1].bond_type[0], BOND_TYPE_ALTERN as u8);

        let (mut source_length, original) = fixture(BOND_DOUBLE as u8);
        source_length[1].valence = 0;
        assert_eq!(
            replace_arom_bonds(&mut source_length, 1, &original, 2),
            Ok(1)
        );
        assert_eq!(source_length[0].bond_type[0], BOND_TYPE_ALTERN as u8);

        let mut negative_valence = [inp_ATOM {
            valence: -1,
            ..inp_ATOM::default()
        }];
        assert_eq!(replace_arom_bonds(&mut negative_valence, 1, &[], -1), Ok(0));

        assert_eq!(
            replace_arom_bonds(&mut [], 1, &[], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let (mut bad_neighbor, original) = fixture(BOND_DOUBLE as u8);
        bad_neighbor[0].neighbor[0] = 9;
        assert_eq!(
            replace_arom_bonds(&mut bad_neighbor, 1, &original, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let (mut short_original, original) = fixture(BOND_DOUBLE as u8);
        assert_eq!(
            replace_arom_bonds(&mut short_original, 1, &original[..1], 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__add_explicit_h__line_5868() {
        fn hydrogen(isotope: i8, bond_type: u8) -> inp_ATOM {
            let mut atom = atom(&[0]);
            atom.elname[..3].copy_from_slice(&[b'H' as i8, b'2' as i8, 0]);
            atom.el_number = EL_NUMBER_H;
            atom.iso_atw_diff = isotope;
            atom.bond_type[0] = bond_type;
            atom.bond_stereo[0] = i8::MIN;
            atom.charge = -5;
            atom
        }

        let mut heap = SourceHeap::default();
        let mut heavy = inp_ATOM {
            num_H: 1,
            ..inp_ATOM::default()
        };
        heavy.el_number = EL_NUMBER_C;
        let mut skipped = atom(&[7]);
        skipped.elname[..4].copy_from_slice(&[b'B' as i8, b'A' as i8, b'D' as i8, 0]);
        skipped.el_number = EL_NUMBER_C;
        skipped.bond_type[0] = BOND_TRIPLE as u8;
        skipped.charge = 12;
        let eligible = hydrogen(0, BOND_DOUBLE as u8);
        let atoms = heap
            .allocate_model_storage(vec![heavy, skipped.clone(), eligible.clone()])
            .unwrap();
        let mut input = INP_ATOM_DATA {
            at: atoms,
            num_at: 1,
            num_removed_H: 2,
            ..INP_ATOM_DATA::default()
        };
        heap.trace_source_allocations();
        assert_eq!(add_explicit_H(&mut heap, &mut input), Ok(2));
        assert_eq!(heap.source_allocation_calls(), 1);
        assert_eq!((input.num_at, input.num_removed_H), (2, 0));
        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(values[0].neighbor[0], 1);
        assert_eq!(values[0].bond_type[0], BOND_DOUBLE as u8);
        assert_eq!(values[0].bond_stereo[0], i8::MIN.wrapping_neg());
        assert_eq!(values[0].valence, 1);
        assert_eq!(values[0].chem_bonds_valence, BOND_DOUBLE as i8);
        assert_eq!(values[0].num_H, 0);
        assert_eq!(values[1].elname, eligible.elname);
        assert_eq!(values[1].el_number, EL_NUMBER_H);
        assert_eq!(values[1].neighbor[0], skipped.neighbor[0]);
        assert_eq!(values[1].bond_type[0], skipped.bond_type[0]);
        assert_eq!(values[1].charge, skipped.charge);
        assert_eq!(values[2], eligible);

        for isotope in 1_i8..=3 {
            let mut isotope_heap = SourceHeap::default();
            let mut isotope_heavy = inp_ATOM {
                num_H: 1,
                ..inp_ATOM::default()
            };
            isotope_heavy.num_iso_H[usize::try_from(isotope - 1).unwrap()] = 1;
            let pointer = isotope_heap
                .allocate_model_storage(vec![isotope_heavy, hydrogen(isotope, BOND_SINGLE as u8)])
                .unwrap();
            let mut isotope_input = INP_ATOM_DATA {
                at: pointer,
                num_at: 1,
                num_removed_H: 1,
                ..INP_ATOM_DATA::default()
            };
            assert_eq!(add_explicit_H(&mut isotope_heap, &mut isotope_input), Ok(2));
            let values = isotope_heap.slice(pointer.as_const()).unwrap();
            assert_eq!(values[0].num_H, 0);
            assert_eq!(values[0].num_iso_H, [0, 0, 0]);
        }

        let mut invalid_bond_heap = SourceHeap::default();
        let invalid_bond_pointer = invalid_bond_heap
            .allocate_model_storage(vec![
                inp_ATOM {
                    num_H: 1,
                    ..inp_ATOM::default()
                },
                hydrogen(0, 9),
            ])
            .unwrap();
        let mut invalid_bond_input = INP_ATOM_DATA {
            at: invalid_bond_pointer,
            num_at: 1,
            num_removed_H: 1,
            ..INP_ATOM_DATA::default()
        };
        assert_eq!(
            add_explicit_H(&mut invalid_bond_heap, &mut invalid_bond_input),
            Ok(2)
        );
        let values = invalid_bond_heap
            .slice(invalid_bond_pointer.as_const())
            .unwrap();
        assert_eq!(values[0].bond_type[0], BOND_SINGLE as u8);
        assert_eq!(values[0].chem_bonds_valence, BOND_SINGLE as i8);
        assert_eq!(values[1].bond_type[0], 9);

        let no_add_cases = [
            {
                let mut value = hydrogen(0, BOND_SINGLE as u8);
                value.el_number = EL_NUMBER_C;
                (
                    inp_ATOM {
                        num_H: 1,
                        ..inp_ATOM::default()
                    },
                    value,
                )
            },
            {
                let mut value = hydrogen(0, BOND_SINGLE as u8);
                value.valence = 2;
                (
                    inp_ATOM {
                        num_H: 1,
                        ..inp_ATOM::default()
                    },
                    value,
                )
            },
            {
                let mut value = hydrogen(0, BOND_SINGLE as u8);
                value.neighbor[0] = 1;
                (
                    inp_ATOM {
                        num_H: 1,
                        ..inp_ATOM::default()
                    },
                    value,
                )
            },
            (
                inp_ATOM {
                    num_H: 1,
                    ..inp_ATOM::default()
                },
                hydrogen(-1, BOND_SINGLE as u8),
            ),
            (
                inp_ATOM {
                    num_H: 1,
                    ..inp_ATOM::default()
                },
                hydrogen(4, BOND_SINGLE as u8),
            ),
            (inp_ATOM::default(), hydrogen(0, BOND_SINGLE as u8)),
            (
                inp_ATOM {
                    num_H: 1,
                    num_iso_H: [2, 0, 0],
                    ..inp_ATOM::default()
                },
                hydrogen(0, BOND_SINGLE as u8),
            ),
            (
                inp_ATOM {
                    num_H: 1,
                    num_iso_H: [1, 0, 0],
                    ..inp_ATOM::default()
                },
                hydrogen(0, BOND_SINGLE as u8),
            ),
            (
                inp_ATOM {
                    num_H: 1,
                    ..inp_ATOM::default()
                },
                hydrogen(1, BOND_SINGLE as u8),
            ),
        ];
        for (heavy, removed) in no_add_cases {
            let mut case_heap = SourceHeap::default();
            let pointer = case_heap
                .allocate_model_storage(vec![heavy.clone(), removed])
                .unwrap();
            let mut case_input = INP_ATOM_DATA {
                at: pointer,
                num_at: 1,
                num_removed_H: 1,
                ..INP_ATOM_DATA::default()
            };
            assert_eq!(add_explicit_H(&mut case_heap, &mut case_input), Ok(1));
            assert_eq!((case_input.num_at, case_input.num_removed_H), (1, 1));
            assert_eq!(case_heap.slice(pointer.as_const()).unwrap()[0], heavy);
        }

        let mut failure_heap = SourceHeap::default();
        let failure_pointer = failure_heap
            .allocate_model_storage(vec![
                inp_ATOM {
                    num_H: 1,
                    ..inp_ATOM::default()
                },
                hydrogen(0, BOND_SINGLE as u8),
            ])
            .unwrap();
        let before = failure_heap
            .slice(failure_pointer.as_const())
            .unwrap()
            .to_vec();
        let mut failure_input = INP_ATOM_DATA {
            at: failure_pointer,
            num_at: 1,
            num_removed_H: 1,
            ..INP_ATOM_DATA::default()
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(add_explicit_H(&mut failure_heap, &mut failure_input), Ok(1));
        assert_eq!((failure_input.num_at, failure_input.num_removed_H), (1, 1));
        assert_eq!(
            failure_heap.slice(failure_pointer.as_const()).unwrap(),
            before
        );

        let mut zero_heap = SourceHeap::default();
        let mut zero_input = INP_ATOM_DATA {
            num_at: -7,
            num_removed_H: 0,
            ..INP_ATOM_DATA::default()
        };
        assert_eq!(add_explicit_H(&mut zero_heap, &mut zero_input), Ok(-7));
        assert_eq!(zero_heap.source_allocation_calls(), 0);

        let mut bounds_heap = SourceHeap::default();
        let bounds_pointer = bounds_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut bounds_input = INP_ATOM_DATA {
            at: bounds_pointer,
            num_at: 1,
            num_removed_H: 1,
            ..INP_ATOM_DATA::default()
        };
        assert_eq!(
            add_explicit_H(&mut bounds_heap, &mut bounds_input),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((bounds_input.num_at, bounds_input.num_removed_H), (1, 1));
    }

    #[test]
    fn source_port__ichinorm__fill_out_bond_cuts__line_5628() {
        let sentinel = R2C_ATPAIR {
            at: [91, 92],
            atno: 93,
        };

        let mut empty_output = [sentinel.clone()];
        assert_eq!(
            fill_out_bond_cuts(&[], &[], -1, &mut empty_output, 1),
            Ok(0)
        );
        assert_eq!(empty_output[0], sentinel);

        let mut ring_atom = atom(&[5, 2]);
        ring_atom.neighbor[2] = u16::MAX;
        let mut ring_derivative = DERIV_AT::default();
        ring_derivative.typ[..2].fill(DERIV_RING_O_OUTSIDE_PRECURSOR as i16);
        ring_derivative.ord[..2].copy_from_slice(&[0, 1]);
        ring_derivative.num[..2].copy_from_slice(&[-1, MAX_AT_DERIV as i8]);
        let mut ring_output = [R2C_ATPAIR::default(), R2C_ATPAIR::default()];
        assert_eq!(
            fill_out_bond_cuts(
                &[ring_atom.clone()],
                &[ring_derivative.clone()],
                1,
                &mut ring_output,
                2,
            ),
            Ok(2)
        );
        assert_eq!(ring_output[0].at, [0, 2]);
        assert_eq!(ring_output[1].at, [0, 5]);
        assert!(ring_output.iter().all(|pair| pair.atno == 0));

        let mut no_room = [sentinel.clone(), sentinel.clone()];
        assert_eq!(
            fill_out_bond_cuts(
                &[ring_atom.clone()],
                &[ring_derivative.clone()],
                1,
                &mut no_room,
                1,
            ),
            Ok(-2)
        );
        assert_eq!(no_room, [sentinel.clone(), sentinel.clone()]);

        let mut excessive_size = ring_derivative.clone();
        excessive_size.num[0] = MAX_AT_DERIV as i8 + 1;
        let mut wrong_ring_type = [sentinel.clone(), sentinel.clone()];
        assert_eq!(
            fill_out_bond_cuts(
                &[ring_atom.clone()],
                &[excessive_size],
                1,
                &mut wrong_ring_type,
                2,
            ),
            Ok(-2)
        );
        assert_eq!(wrong_ring_type, [sentinel.clone(), sentinel.clone()]);

        fn dmox_derivatives() -> [DERIV_AT; 2] {
            let mut derivatives = [DERIV_AT::default(), DERIV_AT::default()];
            derivatives[0].typ[0] = DERIV_RING_DMOX_DEOX_N as i16;
            derivatives[0].other_atom = 2;
            derivatives[1].typ[0] = DERIV_RING_DMOX_DEOX_O as i16;
            derivatives[1].other_atom = 1;
            derivatives
        }
        let dmox_atoms = [atom(&[3]), atom(&[0])];
        let dmox = dmox_derivatives();
        let mut dmox_output = [R2C_ATPAIR::default(), R2C_ATPAIR::default()];
        assert_eq!(
            fill_out_bond_cuts(&dmox_atoms, &dmox, 2, &mut dmox_output, 2),
            Ok(2)
        );
        assert_eq!(dmox_output[0].at, [0, 3]);
        assert_eq!(dmox_output[1].at, [0, 1]);
        assert_eq!(dmox_output[0].atno, 0);
        assert_eq!(dmox_output[1].atno, 0);

        let mut dmox_no_room = [sentinel.clone(), sentinel.clone()];
        assert_eq!(
            fill_out_bond_cuts(&dmox_atoms, &dmox, 2, &mut dmox_no_room, 1),
            Ok(-2)
        );
        assert_eq!(dmox_no_room, [sentinel.clone(), sentinel.clone()]);

        let malformed_dmox = [
            {
                let mut value = dmox_derivatives();
                value[0].typ[1] = DERIV_BRIDGE_O as i16;
                value
            },
            {
                let mut value = dmox_derivatives();
                value[0].other_atom = 0;
                value
            },
            {
                let mut value = dmox_derivatives();
                value[0].other_atom = 1;
                value
            },
            {
                let mut value = dmox_derivatives();
                value[1].other_atom = 2;
                value
            },
            {
                let mut value = dmox_derivatives();
                value[1].typ[0] = DERIV_BRIDGE_O as i16;
                value
            },
            {
                let mut value = dmox_derivatives();
                value[1].typ[1] = DERIV_BRIDGE_O as i16;
                value
            },
        ];
        for derivatives in malformed_dmox {
            let mut output = [sentinel.clone(), sentinel.clone()];
            assert_eq!(
                fill_out_bond_cuts(&dmox_atoms, &derivatives, 2, &mut output, 2),
                Ok(-3)
            );
            assert_eq!(output, [sentinel.clone(), sentinel.clone()]);
        }

        let mut ordinary = DERIV_AT::default();
        ordinary.typ[..2].copy_from_slice(&[DERIV_BRIDGE_O as i16, DERIV_BRIDGE_NH as i16]);
        ordinary.ord[..2].copy_from_slice(&[0, 2]);
        let mut ordinary_output = [R2C_ATPAIR::default(), R2C_ATPAIR::default()];
        assert_eq!(
            fill_out_bond_cuts(
                &[ring_atom.clone()],
                &[ordinary.clone()],
                1,
                &mut ordinary_output,
                2,
            ),
            Ok(2)
        );
        assert_eq!(ordinary_output[0].at, [0, 5]);
        assert_eq!(ordinary_output[1].at, [0, u16::MAX]);

        let mut partial_output = [sentinel.clone(), sentinel.clone()];
        assert_eq!(
            fill_out_bond_cuts(
                &[ring_atom.clone()],
                &[ordinary.clone()],
                1,
                &mut partial_output,
                1,
            ),
            Ok(-2)
        );
        assert_eq!(partial_output[0].at, [0, 5]);
        assert_eq!(partial_output[0].atno, 0);
        assert_eq!(partial_output[1], sentinel);

        let mut terminated = ordinary.clone();
        terminated.typ[0] = 0;
        let mut terminated_output = [sentinel.clone()];
        assert_eq!(
            fill_out_bond_cuts(
                &[ring_atom.clone()],
                &[terminated],
                1,
                &mut terminated_output,
                1,
            ),
            Ok(0)
        );
        assert_eq!(terminated_output[0], sentinel);

        assert_eq!(
            fill_out_bond_cuts(&[ring_atom.clone()], &[], 1, &mut [], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut bad_order = ordinary.clone();
        bad_order.ord[0] = -1;
        assert_eq!(
            fill_out_bond_cuts(
                &[ring_atom.clone()],
                &[bad_order],
                1,
                &mut [sentinel.clone()],
                1
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            fill_out_bond_cuts(&[ring_atom], &[ordinary], 1, &mut [], 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut missing_other = dmox_derivatives();
        missing_other[0].other_atom = 3;
        assert_eq!(
            fill_out_bond_cuts(&dmox_atoms, &missing_other, 1, &mut dmox_output, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__unmarkotherindicators__line_681() {
        let mut atoms = [inp_ATOM::default(), inp_ATOM::default()];
        atoms[0].at_type = 17;
        atoms[0].cFlags = 23;
        atoms[0].component = 9;
        atoms[1].at_type = u16::MAX;
        atoms[1].cFlags = i8::MIN;
        atoms[1].component = 11;

        assert_eq!(UnMarkOtherIndicators(&mut atoms, -1), Ok(0));
        assert_eq!(atoms[0].at_type, 17);
        assert_eq!(UnMarkOtherIndicators(&mut atoms, 1), Ok(0));
        assert_eq!(atoms[0].at_type, 0);
        assert_eq!(atoms[0].cFlags, 0);
        assert_eq!(atoms[0].component, 9);
        assert_eq!(atoms[1].at_type, u16::MAX);
        assert_eq!(atoms[1].cFlags, i8::MIN);
        assert_eq!(atoms[1].component, 11);
        assert_eq!(UnMarkOtherIndicators(&mut atoms, 2), Ok(0));
        assert_eq!(atoms[1].at_type, 0);
        assert_eq!(atoms[1].cFlags, 0);
        assert_eq!(atoms[1].component, 11);
        assert_eq!(
            UnMarkOtherIndicators(&mut atoms, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__unmarkdisconnectedcomponents__line_652() {
        let mut heap = SourceHeap::default();
        let mut atoms = [
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.orig_compt_at_numb = u16::try_from(index + 10).unwrap();
            atom.component = u16::try_from(index + 20).unwrap();
            atom.at_type = u16::try_from(index + 30).unwrap();
        }
        let atom_pointer = heap.allocate_model_storage(atoms.to_vec()).unwrap();
        let current_lengths = heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let old_components = heap.allocate_model_storage(vec![7_u16, 8, 9]).unwrap();
        let mut original = ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: 2,
            num_components: 9,
            nCurAtLen: current_lengths,
            nOldCompNumber: old_components,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            UnMarkDisconnectedComponents(&mut heap, &mut original),
            Ok(0)
        );
        assert_eq!(original.num_components, 0);
        assert!(original.nCurAtLen.is_null());
        assert!(original.nOldCompNumber.is_null());
        assert_eq!(
            heap.slice(current_lengths.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(old_components.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        let values = heap.slice(atom_pointer.as_const()).unwrap();
        assert_eq!(values[0].orig_compt_at_numb, 0);
        assert_eq!(values[0].component, 0);
        assert_eq!(values[0].at_type, 30);
        assert_eq!(values[1].orig_compt_at_numb, 0);
        assert_eq!(values[1].component, 0);
        assert_eq!(values[1].at_type, 31);
        assert_eq!(values[2].orig_compt_at_numb, 12);
        assert_eq!(values[2].component, 22);

        let mut negative_heap = SourceHeap::default();
        let negative_current = negative_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let mut negative = ORIG_ATOM_DATA {
            num_inp_atoms: -1,
            num_components: i32::MAX,
            nCurAtLen: negative_current,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            UnMarkDisconnectedComponents(&mut negative_heap, &mut negative),
            Ok(0)
        );
        assert_eq!(negative.num_components, 0);
        assert!(negative.nCurAtLen.is_null());

        let mut bounds_heap = SourceHeap::default();
        let mut bounds_atom = inp_ATOM::default();
        bounds_atom.component = 4;
        bounds_atom.orig_compt_at_numb = 5;
        let bounds_pointer = bounds_heap
            .allocate_model_storage(vec![bounds_atom])
            .unwrap();
        let bounds_current = bounds_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let mut bounds = ORIG_ATOM_DATA {
            at: bounds_pointer,
            num_inp_atoms: 2,
            num_components: 3,
            nCurAtLen: bounds_current,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            UnMarkDisconnectedComponents(&mut bounds_heap, &mut bounds),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let value = &bounds_heap.slice(bounds_pointer.as_const()).unwrap()[0];
        assert_eq!((value.orig_compt_at_numb, value.component), (0, 0));
        assert_eq!(bounds.num_components, 3);
        assert_eq!(bounds.nCurAtLen, bounds_current);

        let mut free_heap = SourceHeap::default();
        let free_atoms = free_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let free_current = free_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let free_old_base = free_heap.allocate_model_storage(vec![2_u16, 3]).unwrap();
        let free_old_interior = free_old_base.offset(1).unwrap();
        let mut free_error = ORIG_ATOM_DATA {
            at: free_atoms,
            num_inp_atoms: 1,
            num_components: 7,
            nCurAtLen: free_current,
            nOldCompNumber: free_old_interior,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            UnMarkDisconnectedComponents(&mut free_heap, &mut free_error),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert!(free_error.nCurAtLen.is_null());
        assert_eq!(free_error.nOldCompNumber, free_old_interior);
        assert_eq!(free_error.num_components, 7);
        assert_eq!(
            free_heap.slice(free_current.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(free_heap.slice(free_old_base.as_const()).unwrap(), &[2, 3]);
    }

    #[test]
    fn source_port__ichinorm__unmarkonecomponent__line_698() {
        let mut atoms = [inp_ATOM::default(), inp_ATOM::default()];
        atoms[0].orig_compt_at_numb = 11;
        atoms[0].component = 12;
        atoms[0].at_type = 13;
        atoms[0].cFlags = 14;
        atoms[1].orig_compt_at_numb = u16::MAX;
        atoms[1].component = u16::MAX;
        atoms[1].at_type = u16::MAX;
        atoms[1].cFlags = i8::MIN;

        assert_eq!(UnMarkOneComponent(&mut atoms, -1), Ok(0));
        assert_eq!(atoms[0].component, 12);
        assert_eq!(UnMarkOneComponent(&mut atoms, 1), Ok(0));
        assert_eq!(atoms[0].orig_compt_at_numb, 0);
        assert_eq!(atoms[0].component, 0);
        assert_eq!(atoms[0].at_type, 13);
        assert_eq!(atoms[0].cFlags, 14);
        assert_eq!(atoms[1].component, u16::MAX);
        assert_eq!(UnMarkOneComponent(&mut atoms, 2), Ok(0));
        assert_eq!(atoms[1].orig_compt_at_numb, 0);
        assert_eq!(atoms[1].component, 0);
        assert_eq!(atoms[1].at_type, u16::MAX);
        assert_eq!(atoms[1].cFlags, i8::MIN);
        assert_eq!(
            UnMarkOneComponent(&mut atoms, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__subtract_dt_from_num_h__line_734() {
        let mut atoms = [
            inp_ATOM {
                num_H: 10,
                num_iso_H: [1, 2, 3],
                at_type: 41,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                num_H: -120,
                num_iso_H: [10, -20, 100],
                at_type: 42,
                ..inp_ATOM::default()
            },
        ];
        let before_negative_at_type = atoms[1].at_type;

        assert_eq!(subtract_DT_from_num_H(&mut atoms, -1), Ok(0));
        assert_eq!(atoms[0].num_H, 10);
        assert_eq!(subtract_DT_from_num_H(&mut atoms, 2), Ok(0));
        assert_eq!(atoms[0].num_H, 4);
        assert_eq!(atoms[1].num_H, 46);
        assert_eq!(atoms[0].num_iso_H, [1, 2, 3]);
        assert_eq!(atoms[0].at_type, 41);
        assert_eq!(atoms[1].at_type, before_negative_at_type);
        assert_eq!(
            subtract_DT_from_num_H(&mut atoms, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__add_inp_atom__line_749() {
        let mut first = inp_ATOM {
            valence: 2,
            at_type: 71,
            num_H: 3,
            ..inp_ATOM::default()
        };
        first.neighbor[0] = 0;
        first.neighbor[1] = u16::MAX;
        let second = inp_ATOM {
            valence: -1,
            at_type: 72,
            num_H: -4,
            neighbor: [9; 20],
            ..inp_ATOM::default()
        };
        let added = [first.clone(), second.clone()];

        let mut unchanged = vec![inp_ATOM::default(); 3];
        unchanged[0].at_type = 99;
        let before = unchanged.clone();
        assert_eq!(add_inp_ATOM(&mut unchanged, 3, -2, &added, 2), Ok(-2));
        assert_eq!(unchanged, before);
        assert_eq!(add_inp_ATOM(&mut unchanged, 3, 0, &added, -3), Ok(-3));
        assert_eq!(unchanged, before);
        assert_eq!(add_inp_ATOM(&mut unchanged, 1, 1, &added, 1), Ok(-1));
        assert_eq!(unchanged, before);

        let mut zero_offset = vec![inp_ATOM::default(); 2];
        assert_eq!(add_inp_ATOM(&mut zero_offset, 2, 0, &added, 2), Ok(2));
        assert_eq!(zero_offset, added);

        let mut destination = vec![inp_ATOM::default(); 4];
        destination[0].at_type = 17;
        assert_eq!(add_inp_ATOM(&mut destination, 4, 1, &added, 2), Ok(3));
        assert_eq!(destination[0].at_type, 17);
        assert_eq!(destination[1].at_type, first.at_type);
        assert_eq!(destination[1].num_H, first.num_H);
        assert_eq!(destination[1].neighbor[0], 1);
        assert_eq!(destination[1].neighbor[1], 0);
        assert_eq!(destination[2], second);
        assert_eq!(destination[3], inp_ATOM::default());

        assert_eq!(
            add_inp_ATOM(&mut destination, 8, 4, &added, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__oad_edit_mergecomponentsandrecreateoad__line_7563() {
        let mut empty_heap = SourceHeap::default();
        let mut empty_original = ORIG_ATOM_DATA::default();
        let mut empty_error = 17;
        assert_eq!(
            OAD_Edit_MergeComponentsAndRecreateOAD(
                &mut empty_heap,
                &mut empty_original,
                &mut [],
                0,
                &mut empty_error,
            ),
            Ok(())
        );
        assert_eq!(empty_error, -999);

        let component_atoms = |base: u16| {
            let mut first = atom(&[1]);
            first.at_type = base + 10;
            first.cFlags = 7;
            first.orig_compt_at_numb = base + 20;
            first.component = base + 30;
            first.bCutVertex = 1;
            first.nRingSystem = base + 40;
            first.nNumAtInRingSystem = 2;
            first.nBlockSystem = base + 50;
            first.num_H = 3;
            first.num_iso_H = [1, 1, 0];
            let mut second = atom(&[0]);
            second.at_type = base + 11;
            second.cFlags = -8;
            second.orig_compt_at_numb = base + 21;
            second.component = base + 31;
            second.bCutVertex = 1;
            second.nRingSystem = base + 40;
            second.nNumAtInRingSystem = 2;
            second.nBlockSystem = base + 50;
            [first, second]
        };

        let mut heap = SourceHeap::default();
        let first = heap
            .allocate_model_storage(component_atoms(0).to_vec())
            .unwrap();
        let second = heap
            .allocate_model_storage(component_atoms(100).to_vec())
            .unwrap();
        let mut current = [
            INP_ATOM_DATA {
                at: first,
                num_at: 2,
                ..INP_ATOM_DATA::default()
            },
            INP_ATOM_DATA {
                at: second,
                num_at: 2,
                ..INP_ATOM_DATA::default()
            },
        ];
        let old_atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let old_coordinates = heap.allocate_model_storage(vec![[5_i8; 32]]).unwrap();
        let old_lengths = heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let old_components = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let mut original = ORIG_ATOM_DATA {
            at: old_atoms,
            num_inp_atoms: 1,
            num_components: 2,
            nCurAtLen: old_lengths,
            nOldCompNumber: old_components,
            szCoord: old_coordinates,
            ..ORIG_ATOM_DATA::default()
        };
        let mut error_code = 73;

        assert_eq!(
            OAD_Edit_MergeComponentsAndRecreateOAD(
                &mut heap,
                &mut original,
                &mut current,
                2,
                &mut error_code,
            ),
            Ok(())
        );
        assert_eq!(error_code, 73);
        assert_eq!(original.num_inp_atoms, 4);
        assert_eq!(original.num_components, 0);
        assert!(original.szCoord.is_null());
        assert!(original.nCurAtLen.is_null());
        assert!(original.nOldCompNumber.is_null());
        assert_eq!(
            heap.slice(old_atoms.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(old_coordinates.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        let merged = heap.slice(original.at.as_const()).unwrap();
        assert_eq!(merged.len(), 4);
        assert_eq!(merged[0].neighbor[0], 1);
        assert_eq!(merged[1].neighbor[0], 0);
        assert_eq!(merged[2].neighbor[0], 3);
        assert_eq!(merged[3].neighbor[0], 2);
        assert_eq!(merged[0].num_H, 1);
        for merged_atom in merged {
            assert_eq!(merged_atom.at_type, 0);
            assert_eq!(merged_atom.cFlags, 0);
            assert_eq!(merged_atom.orig_compt_at_numb, 0);
            assert_eq!(merged_atom.component, 0);
            assert_eq!(merged_atom.bCutVertex, 0);
            assert_eq!(merged_atom.nRingSystem, 0);
            assert_eq!(merged_atom.nNumAtInRingSystem, 0);
            assert_eq!(merged_atom.nBlockSystem, 0);
        }
        assert_eq!(heap.slice(first.as_const()).unwrap()[0].num_H, 1);
        assert_eq!(heap.slice(second.as_const()).unwrap()[0].num_H, 1);
    }

    #[test]
    fn source_port__ichinorm__remove_cut_derivs__line_7666() {
        let mut first = atom(&[1]);
        first.num_H = 2;
        first.num_iso_H = [1, 0, 0];
        first.at_type = 41;
        first.cFlags = 7;
        first.bCutVertex = 1;
        first.nRingSystem = 9;
        let mut second = atom(&[0]);
        second.num_H = 1;
        second.at_type = 42;
        second.cFlags = -7;
        second.bCutVertex = 1;
        second.nRingSystem = 9;
        let mut tritium_derivative = inp_ATOM::default();
        tritium_derivative.num_iso_H[2] = 1;
        tritium_derivative.at_type = 43;
        tritium_derivative.cFlags = 8;

        let mut heap = SourceHeap::default();
        let original_atoms = heap
            .allocate_model_storage(vec![first, second, tritium_derivative])
            .unwrap();
        let mut current = [INP_ATOM_DATA {
            at: original_atoms,
            num_at: 3,
            ..INP_ATOM_DATA::default()
        }];
        let mut error_code = 73;

        assert_eq!(
            remove_cut_derivs(
                &mut heap,
                3,
                original_atoms,
                &mut current,
                0,
                &mut error_code,
            ),
            Ok(())
        );
        assert_eq!(error_code, 73);
        assert_eq!(current[0].num_at, 2);
        assert_ne!(current[0].at, original_atoms);
        assert_eq!(
            heap.slice(original_atoms.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        let merged = heap.slice(current[0].at.as_const()).unwrap();
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].neighbor[0], 1);
        assert_eq!(merged[1].neighbor[0], 0);
        assert_eq!(merged[0].num_H, 2);
        assert_eq!(merged[0].num_iso_H, [1, 0, 0]);
        assert_eq!(merged[0].component, 1);
        assert_eq!(merged[1].component, 1);
        assert_eq!(merged[0].orig_compt_at_numb, 1);
        assert_eq!(merged[1].orig_compt_at_numb, 2);
        for value in merged {
            assert_eq!(value.at_type, 0);
            assert_eq!(value.cFlags, 0);
            assert_eq!(value.bCutVertex, 0);
            assert_eq!(value.nRingSystem, 0);
        }

        let mut no_components_heap = SourceHeap::default();
        let no_components_atoms = no_components_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut no_components_current = [INP_ATOM_DATA {
            at: no_components_atoms,
            num_at: 1,
            ..INP_ATOM_DATA::default()
        }];
        let mut no_components_error = 11;
        assert_eq!(
            remove_cut_derivs(
                &mut no_components_heap,
                0,
                no_components_atoms,
                &mut no_components_current,
                0,
                &mut no_components_error,
            ),
            Ok(())
        );
        assert_eq!(no_components_error, 11);
        assert_eq!(no_components_current[0].num_at, 0);
        assert_eq!(no_components_current[0].at, no_components_atoms);
        assert_eq!(
            no_components_heap
                .slice(no_components_atoms.as_const())
                .map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__mol2atom__freeinpatomdata__line_1058() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default(); 2])
            .unwrap();
        let fixed = heap
            .allocate_model_storage(vec![inp_ATOM::default(); 1])
            .unwrap();
        let mut input = INP_ATOM_DATA {
            at: atoms,
            at_fixed_bonds: fixed,
            num_at: 2,
            num_removed_H: 3,
            num_bonds: 4,
            num_isotopic: 5,
            bExists: 6,
            bDeleted: 7,
            bHasIsotopicLayer: 8,
            bTautomeric: 9,
            bTautPreprocessed: 10,
            nNumRemovedProtons: 11,
            nNumRemovedProtonsIsotopic: [12, 13, 14],
            num_iso_H: [15, 16, 17],
            bTautFlags: 18,
            bTautFlagsDone: 19,
            bNormalizationFlags: 20,
        };
        assert_eq!(FreeInpAtomData(&mut heap, &mut input), Ok(()));
        assert_eq!(input, INP_ATOM_DATA::default());
        assert_eq!(
            heap.slice(atoms.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(fixed.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(FreeInpAtomData(&mut heap, &mut input), Ok(()));
        assert_eq!(input, INP_ATOM_DATA::default());
    }

    #[test]
    fn source_port__ichinorm__free_underiv_temp_data__line_7632() {
        let mut heap = SourceHeap::default();
        let atom_pairs = heap
            .allocate_model_storage(vec![R2C_ATPAIR::default()])
            .unwrap();
        let derivatives = heap
            .allocate_model_storage(vec![DERIV_AT::default()])
            .unwrap();
        let atoms2 = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let first_atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let first_fixed = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let second_atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let second_fixed = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let current = heap
            .allocate_model_storage(vec![
                INP_ATOM_DATA {
                    at: first_atoms,
                    at_fixed_bonds: first_fixed,
                    num_at: 1,
                    ..INP_ATOM_DATA::default()
                },
                INP_ATOM_DATA {
                    at: second_atoms,
                    at_fixed_bonds: second_fixed,
                    num_at: 1,
                    ..INP_ATOM_DATA::default()
                },
            ])
            .unwrap();

        assert_eq!(
            free_underiv_temp_data(&mut heap, atom_pairs, derivatives, atoms2, current, 2,),
            Ok(())
        );
        assert_eq!(
            heap.slice(atom_pairs.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(derivatives.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(atoms2.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(first_atoms.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(first_fixed.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(second_atoms.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(second_fixed.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(current.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );

        assert_eq!(
            free_underiv_temp_data(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                -1,
            ),
            Ok(())
        );
    }

    #[test]
    fn source_port__ichinorm__make_single_cut__line_5495() {
        fn pair() -> (Vec<inp_ATOM>, Vec<DERIV_AT>) {
            (
                vec![atom(&[1]), atom(&[0])],
                vec![DERIV_AT::default(), DERIV_AT::default()],
            )
        }

        let (mut duplicate_atoms, mut duplicate) = pair();
        duplicate[0].typ[0] = (DERIV_BRIDGE_O | DERIV_DUPLIC) as i16;
        duplicate[0].ord[0] = -7;
        assert_eq!(
            make_single_cut(&mut duplicate_atoms, &mut duplicate, 0, 0),
            Ok(0)
        );
        assert_eq!(duplicate_atoms[0].valence, 1);

        let (mut negative_atoms, mut negative) = pair();
        negative[0].typ[0] = DERIV_BRIDGE_O as i16;
        negative[0].ord[0] = -1;
        assert_eq!(
            make_single_cut(&mut negative_atoms, &mut negative, 0, 0),
            Ok(-1)
        );

        let mut missing_atoms = vec![atom(&[1]), atom(&[])];
        let mut missing = vec![DERIV_AT::default(); 2];
        missing[0].typ[0] = DERIV_BRIDGE_O as i16;
        assert_eq!(
            make_single_cut(&mut missing_atoms, &mut missing, 0, 0),
            Ok(-1)
        );

        let mut ordinary_atoms = vec![atom(&[1, 2]), atom(&[0]), atom(&[0])];
        ordinary_atoms[0].num_H = i8::MAX;
        let mut ordinary = vec![DERIV_AT::default(); 3];
        ordinary[0].typ[..2].fill(DERIV_BRIDGE_O as i16);
        ordinary[0].ord[..2].copy_from_slice(&[0, 1]);
        ordinary[1].typ[0] = DERIV_BRIDGE_O as i16;
        ordinary[1].ord[0] = 0;
        assert_eq!(
            make_single_cut(&mut ordinary_atoms, &mut ordinary, 0, 0),
            Ok(1)
        );
        assert_eq!(ordinary_atoms[0].neighbor[0], 2);
        assert_eq!(ordinary_atoms[0].num_H, i8::MIN);
        assert_eq!(ordinary_atoms[1].num_H, 1);
        assert_eq!(ordinary_atoms[1].num_iso_H[2], 1);
        assert_eq!(ordinary[0].ord[..2], [-1, 0]);
        assert_eq!(ordinary[1].ord[0], -1);
        assert_ne!(ordinary[1].typ[0] & DERIV_DUPLIC as i16, 0);

        let (mut ring_atoms, mut ring) = pair();
        ring[0].typ[0] = DERIV_RING_O_OUTSIDE_PRECURSOR as i16;
        assert_eq!(make_single_cut(&mut ring_atoms, &mut ring, 0, 0), Ok(1));
        assert_eq!(ring_atoms[1].num_H, 1);
        assert_eq!(ring_atoms[0].num_H, 1);
        assert_eq!(ring_atoms[0].num_iso_H, [0, 0, 1]);
        assert_eq!(ring_atoms[1].num_iso_H, [0, 0, 0]);

        let (mut replace_o_atoms, mut replace_o) = pair();
        replace_o_atoms[0].elname[0] = b'N' as i8;
        replace_o_atoms[0].el_number = EL_NUMBER_N;
        replace_o[0].typ[0] = DERIV_X_OXIME as i16;
        assert_eq!(
            make_single_cut(&mut replace_o_atoms, &mut replace_o, 0, 0),
            Ok(1)
        );
        assert_eq!(replace_o_atoms[0].elname[0], b'O' as i8);
        assert_eq!(replace_o_atoms[0].el_number, EL_NUMBER_O);
        assert_eq!(replace_o_atoms[0].num_H, 0);
        assert_eq!(replace_o_atoms[1].num_H, 1);
        assert_eq!(replace_o_atoms[1].num_iso_H[2], 1);

        let (mut replace_oh_atoms, mut replace_oh) = pair();
        replace_oh_atoms[0].el_number = EL_NUMBER_N;
        replace_oh[0].typ[0] = DERIV_RING2_PRRLDD_OUTSIDE_PRECUR as i16;
        assert_eq!(
            make_single_cut(&mut replace_oh_atoms, &mut replace_oh, 0, 0),
            Ok(1)
        );
        assert_eq!(replace_oh_atoms[0].el_number, EL_NUMBER_O);
        assert_eq!(replace_oh_atoms[0].num_H, 0);
        assert_eq!(replace_oh_atoms[1].num_H, 1);

        let (mut second_cut_atoms, mut second_cut) = pair();
        second_cut[0].typ[..2].fill(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR as i16);
        second_cut[0].ord[..2].fill(0);
        assert_eq!(
            make_single_cut(&mut second_cut_atoms, &mut second_cut, 0, 1),
            Ok(1)
        );
        assert_eq!(second_cut_atoms[0].el_number, 0);
        assert_eq!(second_cut_atoms[0].num_H, 1);
        assert_eq!(second_cut_atoms[1].num_H, 1);
        assert_eq!(second_cut_atoms[1].num_iso_H[2], 1);

        let mut partial_atoms = vec![inp_ATOM::default(), atom(&[0])];
        partial_atoms[0].neighbor[0] = 1;
        let mut partial = vec![DERIV_AT::default(); 2];
        partial[0].typ[0] = DERIV_BRIDGE_O as i16;
        assert_eq!(
            make_single_cut(&mut partial_atoms, &mut partial, 0, 0),
            Ok(0)
        );
        assert_eq!(partial_atoms[1].valence, 0);
        assert_eq!(partial_atoms[0].num_H, 0);
        assert_eq!(partial[0].ord[0], 0);

        let (mut invalid_atoms, mut invalid) = pair();
        assert_eq!(
            make_single_cut(&mut invalid_atoms, &mut invalid, -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            make_single_cut(&mut invalid_atoms, &mut invalid, 0, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        invalid[0].ord[0] = 20;
        assert_eq!(
            make_single_cut(&mut invalid_atoms, &mut invalid, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__markringsystemsinp__line_59() {
        let mut heap = SourceHeap::default();
        let triangle = heap
            .allocate_model_storage(vec![atom(&[1, 2]), atom(&[0, 2]), atom(&[0, 1])])
            .unwrap();
        assert_eq!(MarkRingSystemsInp(&mut heap, triangle, 3, 0), Ok(1));
        let values = heap.slice(triangle.as_const()).unwrap();
        assert!(values.iter().all(|atom| atom.nRingSystem == 1));
        assert!(values.iter().all(|atom| atom.nNumAtInRingSystem == 3));
        assert!(values.iter().all(|atom| atom.nBlockSystem == 1));
        assert!(values.iter().all(|atom| atom.bCutVertex == 0));

        let chain = heap
            .allocate_model_storage(vec![atom(&[1]), atom(&[0, 2]), atom(&[1])])
            .unwrap();
        assert_eq!(MarkRingSystemsInp(&mut heap, chain, 3, 0), Ok(3));
        let values = heap.slice(chain.as_const()).unwrap();
        assert!(values.iter().all(|atom| atom.nNumAtInRingSystem == 1));
        assert_eq!(values[1].bCutVertex, 1);

        let disconnected = heap
            .allocate_model_storage(vec![
                atom(&[1]),
                atom(&[0]),
                inp_ATOM {
                    nRingSystem: 91,
                    nNumAtInRingSystem: 92,
                    nBlockSystem: 93,
                    bCutVertex: 7,
                    ..inp_ATOM::default()
                },
            ])
            .unwrap();
        assert_eq!(MarkRingSystemsInp(&mut heap, disconnected, 3, 0), Ok(2));
        let untouched = &heap.slice(disconnected.as_const()).unwrap()[2];
        assert_eq!(
            (
                untouched.nRingSystem,
                untouched.nNumAtInRingSystem,
                untouched.nBlockSystem,
                untouched.bCutVertex,
            ),
            (91, 92, 93, 7)
        );

        let empty = heap.allocate_model_storage(Vec::<inp_ATOM>::new()).unwrap();
        assert_eq!(MarkRingSystemsInp(&mut heap, empty, 0, 0), Ok(CT_OVERFLOW));

        for failure in 0..5 {
            let mut failure_heap = SourceHeap::default();
            let atoms = failure_heap
                .allocate_model_storage(vec![atom(&[])])
                .unwrap();
            failure_heap.fail_after_allocations(failure);
            assert_eq!(
                MarkRingSystemsInp(&mut failure_heap, atoms, 1, 0),
                Ok(CT_OUT_OF_RAM)
            );
        }
    }

    #[test]
    fn source_port__ichinorm__mark_atoms_cflags__line_1368() {
        let mut already_marked = vec![inp_ATOM {
            cFlags: 0b0101,
            valence: 1,
            neighbor: {
                let mut value = [0; 20];
                value[0] = 1;
                value
            },
            ..inp_ATOM::default()
        }];
        assert_eq!(
            mark_atoms_cFlags(&mut already_marked, 0, 17, 0b0100),
            Ok(17)
        );
        assert_eq!(already_marked[0].cFlags, 0b0101);

        let mut graph = vec![atom(&[1, 2]), atom(&[0, 2]), atom(&[0, 1]), atom(&[])];
        graph[0].cFlags = 0b0010;
        graph[1].cFlags = 0b1000;
        graph[2].cFlags = 0b0011;
        graph[3].cFlags = 0b0100;
        assert_eq!(mark_atoms_cFlags(&mut graph, 0, -2, 0b0100), Ok(1));
        assert_eq!(graph[0].cFlags, 0b0110);
        assert_eq!(graph[1].cFlags, 0b1100);
        assert_eq!(graph[2].cFlags, 0b0111);
        assert_eq!(graph[3].cFlags, 0b0100);
        assert_eq!(mark_atoms_cFlags(&mut graph, 2, 91, 0b0100), Ok(91));

        let mut signed_flag = vec![atom(&[1]), atom(&[0])];
        signed_flag[0].cFlags = 1;
        signed_flag[1].cFlags = 2;
        assert_eq!(
            mark_atoms_cFlags(&mut signed_flag, 0, i32::MAX, i8::MIN),
            Ok(i32::MIN + 1)
        );
        assert_eq!(signed_flag[0].cFlags, i8::MIN | 1);
        assert_eq!(signed_flag[1].cFlags, i8::MIN | 2);
    }

    #[test]
    fn source_port__ichinorm__unmark_atoms_cflags__line_1388() {
        let mut not_marked = vec![inp_ATOM {
            cFlags: 0b0011,
            valence: 1,
            neighbor: {
                let mut value = [0; 20];
                value[0] = 1;
                value
            },
            ..inp_ATOM::default()
        }];
        assert_eq!(
            unmark_atoms_cFlags(&mut not_marked, 0, 17, 0b0100, !0b0100),
            Ok(17)
        );
        assert_eq!(not_marked[0].cFlags, 0b0011);

        let mut graph = vec![atom(&[1, 2]), atom(&[0, 2]), atom(&[0, 1]), atom(&[])];
        graph[0].cFlags = 0b1111;
        graph[1].cFlags = 0b1100;
        graph[2].cFlags = 0b0111;
        graph[3].cFlags = 0b0100;
        assert_eq!(
            unmark_atoms_cFlags(&mut graph, 0, -2, 0b0100, 0b1010),
            Ok(1)
        );
        assert_eq!(graph[0].cFlags, 0b1010);
        assert_eq!(graph[1].cFlags, 0b1000);
        assert_eq!(graph[2].cFlags, 0b0010);
        assert_eq!(graph[3].cFlags, 0b0100);
        assert_eq!(
            unmark_atoms_cFlags(&mut graph, 1, 91, 0b0100, !0b0100),
            Ok(91)
        );

        let mut signed_flag = vec![atom(&[1]), atom(&[0])];
        signed_flag[0].cFlags = i8::MIN | 1;
        signed_flag[1].cFlags = i8::MIN | 2;
        assert_eq!(
            unmark_atoms_cFlags(&mut signed_flag, 0, i32::MAX, i8::MIN, 0x7f),
            Ok(i32::MIN + 1)
        );
        assert_eq!(signed_flag[0].cFlags, 1);
        assert_eq!(signed_flag[1].cFlags, 2);
    }

    #[test]
    fn source_port__ichinorm__is_c_or_s_db_o__line_1410() {
        fn center(element: u8, charge: i8, radical: i8, neighbors: &[u16]) -> inp_ATOM {
            let mut value = atom(neighbors);
            value.el_number = element;
            value.charge = charge;
            value.radical = radical;
            value
        }
        fn terminal(element: u8, hydrogens: i8, valence: i8, chemical: i8) -> inp_ATOM {
            inp_ATOM {
                el_number: element,
                num_H: hydrogens,
                valence,
                chem_bonds_valence: chemical,
                ..inp_ATOM::default()
            }
        }

        let oxygen = terminal(EL_NUMBER_O as u8, 0, 1, 2);
        for rejected in [
            center(EL_NUMBER_O as u8, 0, 0, &[1]),
            center(EL_NUMBER_C as u8, 1, 0, &[1]),
            center(EL_NUMBER_C as u8, 0, -1, &[1]),
            center(EL_NUMBER_S as u8, -1, 0, &[1]),
            center(EL_NUMBER_S as u8, 0, 2, &[1]),
        ] {
            assert_eq!(is_C_or_S_DB_O(&[rejected, oxygen.clone()], 0), Ok(0));
        }

        for element in [EL_NUMBER_C as u8, EL_NUMBER_S as u8] {
            assert_eq!(
                is_C_or_S_DB_O(&[center(element, 0, 0, &[1]), oxygen.clone()], 0),
                Ok(1)
            );
            assert_eq!(
                is_C_or_S_DB_O(
                    &[
                        center(element, 0, 0, &[1]),
                        terminal(EL_NUMBER_S as u8, 0, 1, 2),
                    ],
                    0,
                ),
                Ok(1)
            );
        }

        for rejected_neighbor in [
            terminal(EL_NUMBER_C as u8, 0, 1, 2),
            terminal(EL_NUMBER_O as u8, 1, 1, 2),
            terminal(EL_NUMBER_O as u8, 0, 2, 2),
            terminal(EL_NUMBER_O as u8, 0, 1, 1),
            terminal(EL_NUMBER_S as u8, -1, 1, 2),
        ] {
            assert_eq!(
                is_C_or_S_DB_O(
                    &[center(EL_NUMBER_C as u8, 0, 0, &[1]), rejected_neighbor,],
                    0,
                ),
                Ok(0)
            );
        }

        let values = [
            center(EL_NUMBER_C as u8, 0, 0, &[1, 2]),
            terminal(EL_NUMBER_C as u8, 0, 1, 1),
            oxygen,
        ];
        assert_eq!(is_C_or_S_DB_O(&values, 0), Ok(1));
    }

    #[test]
    fn source_port__ichinorm__is_c_db_o__line_1434() {
        fn center(element: u8, charge: i8, radical: i8, valence: i8, chemical: i8) -> inp_ATOM {
            let mut value = atom(&[1, 2, 3]);
            value.el_number = element;
            value.charge = charge;
            value.radical = radical;
            value.valence = valence;
            value.chem_bonds_valence = chemical;
            value
        }
        fn terminal(element: u8, hydrogens: i8, valence: i8, chemical: i8) -> inp_ATOM {
            inp_ATOM {
                el_number: element,
                num_H: hydrogens,
                valence,
                chem_bonds_valence: chemical,
                ..inp_ATOM::default()
            }
        }

        let oxygen = terminal(EL_NUMBER_O, 0, 1, 2);
        for rejected in [
            center(EL_NUMBER_S, 0, 0, 3, 4),
            center(EL_NUMBER_C, 1, 0, 3, 4),
            center(EL_NUMBER_C, 0, -1, 3, 4),
            center(EL_NUMBER_C, 0, 0, 2, 4),
            center(EL_NUMBER_C, 0, 0, 3, 3),
        ] {
            assert_eq!(
                is_C_DB_O(
                    &[rejected, oxygen.clone(), oxygen.clone(), oxygen.clone(),],
                    0
                ),
                Ok(0)
            );
        }

        for rejected_neighbor in [
            terminal(EL_NUMBER_C, 0, 1, 2),
            terminal(EL_NUMBER_O, 1, 1, 2),
            terminal(EL_NUMBER_O, 0, 2, 2),
            terminal(EL_NUMBER_O, 0, 1, 1),
        ] {
            assert_eq!(
                is_C_DB_O(
                    &[
                        center(EL_NUMBER_C, 0, 0, 3, 4),
                        rejected_neighbor.clone(),
                        rejected_neighbor.clone(),
                        rejected_neighbor,
                    ],
                    0
                ),
                Ok(0)
            );
        }

        let carbonyl = [
            center(EL_NUMBER_C, 0, 0, 3, 4),
            terminal(EL_NUMBER_C, 0, 1, 1),
            oxygen.clone(),
            oxygen,
        ];
        assert_eq!(is_C_DB_O(&carbonyl, 0), Ok(2));
    }

    #[test]
    fn source_port__ichinorm__is_c_unsat_not_arom__line_1457() {
        fn unsaturated() -> Vec<inp_ATOM> {
            let mut center = atom(&[1, 2, 3]);
            center.el_number = EL_NUMBER_C;
            center.chem_bonds_valence = 4;
            center.bond_type[..3].copy_from_slice(&[
                BOND_DOUBLE as u8,
                BOND_SINGLE as u8,
                BOND_SINGLE as u8,
            ]);
            vec![
                center,
                inp_ATOM::default(),
                inp_ATOM::default(),
                inp_ATOM::default(),
            ]
        }

        assert_eq!(is_C_unsat_not_arom(&unsaturated(), 0), Ok(1));
        for mutation in 0..6 {
            let mut values = unsaturated();
            match mutation {
                0 => values[0].el_number = EL_NUMBER_N,
                1 => values[0].chem_bonds_valence = 3,
                2 => values[0].chem_bonds_valence = 5,
                3 => values[0].num_H = 1,
                4 => values[0].charge = 1,
                5 => values[0].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_C_unsat_not_arom(&values, 0), Ok(0));
        }

        for element in [EL_NUMBER_O, EL_NUMBER_S] {
            let mut values = unsaturated();
            values[1] = inp_ATOM {
                el_number: element,
                valence: 1,
                chem_bonds_valence: 2,
                ..inp_ATOM::default()
            };
            assert_eq!(is_C_unsat_not_arom(&values, 0), Ok(0));
            for mutation in 0..3 {
                let mut not_terminal = values.clone();
                match mutation {
                    0 => not_terminal[1].num_H = 1,
                    1 => not_terminal[1].valence = 2,
                    2 => not_terminal[1].chem_bonds_valence = 1,
                    _ => unreachable!(),
                }
                assert_eq!(is_C_unsat_not_arom(&not_terminal, 0), Ok(1));
            }
        }

        let mut aromatic = unsaturated();
        aromatic[0].bond_type[2] = BOND_TYPE_ALTERN as u8;
        assert_eq!(is_C_unsat_not_arom(&aromatic, 0), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_phenyl__line_1559() {
        fn phenyl() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 7];
            values[0].nRingSystem = 2;
            values[1] = atom(&[0, 2, 6]);
            values[1].el_number = EL_NUMBER_C;
            values[1].nRingSystem = 1;
            values[1].bCutVertex = 1;
            values[1].nNumAtInRingSystem = 6;
            for index in 2..=6 {
                let previous = if index == 2 { 1 } else { index - 1 };
                let next = if index == 6 { 1 } else { index + 1 };
                values[index] = atom(&[previous as u16, next as u16]);
                values[index].el_number = EL_NUMBER_C;
                values[index].num_H = 1;
            }
            values
        }

        assert_eq!(is_Phenyl(&phenyl(), 0, 1), Ok(1));

        for mutate in 0..8 {
            let mut values = phenyl();
            match mutate {
                0 => values[1].el_number = EL_NUMBER_O,
                1 => values[1].valence = 2,
                2 => values[1].num_H = 1,
                3 => values[1].charge = 1,
                4 => values[1].radical = 1,
                5 => values[1].nRingSystem = values[0].nRingSystem,
                6 => values[1].bCutVertex = 0,
                7 => values[1].nNumAtInRingSystem = 5,
                _ => unreachable!(),
            }
            assert_eq!(is_Phenyl(&values, 0, 1), Ok(0));
        }

        let mut no_ring_neighbor = phenyl();
        no_ring_neighbor[1].neighbor[..3].fill(0);
        assert_eq!(is_Phenyl(&no_ring_neighbor, 0, 1), Ok(0));

        for mutate in 0..5 {
            let mut values = phenyl();
            match mutate {
                0 => values[4].el_number = EL_NUMBER_O,
                1 => values[4].valence = 3,
                2 => values[4].num_H = 0,
                3 => values[4].charge = -1,
                4 => values[4].radical = 2,
                _ => unreachable!(),
            }
            assert_eq!(is_Phenyl(&values, 0, 1), Ok(0));
        }

        let mut not_closed = phenyl();
        not_closed[6].neighbor[1] = 0;
        assert_eq!(is_Phenyl(&not_closed, 0, 1), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_pentafluorophenyl__line_1615() {
        fn pentafluorophenyl() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 12];
            values[0].nRingSystem = 2;
            values[1] = atom(&[0, 2, 6]);
            values[1].el_number = EL_NUMBER_C;
            values[1].nRingSystem = 1;
            values[1].bCutVertex = 1;
            values[1].nNumAtInRingSystem = 6;
            for index in 2..=6 {
                let previous = if index == 2 { 1 } else { index - 1 };
                let next = if index == 6 { 1 } else { index + 1 };
                let fluorine = index + 5;
                values[index] = atom(&[previous as u16, next as u16, fluorine as u16]);
                values[index].el_number = EL_NUMBER_C;
                values[fluorine].el_number = EL_NUMBER_F;
                values[fluorine].chem_bonds_valence = 1;
            }
            values
        }

        assert_eq!(is_PentaFluoroPhenyl(&pentafluorophenyl(), 0, 1), Ok(1));

        for mutate in 0..8 {
            let mut values = pentafluorophenyl();
            match mutate {
                0 => values[1].el_number = EL_NUMBER_O,
                1 => values[1].valence = 2,
                2 => values[1].num_H = 1,
                3 => values[1].charge = 1,
                4 => values[1].radical = 1,
                5 => values[1].nRingSystem = values[0].nRingSystem,
                6 => values[1].bCutVertex = 0,
                7 => values[1].nNumAtInRingSystem = 5,
                _ => unreachable!(),
            }
            assert_eq!(is_PentaFluoroPhenyl(&values, 0, 1), Ok(0));
        }

        let mut no_ring_neighbor = pentafluorophenyl();
        no_ring_neighbor[1].neighbor[..3].fill(0);
        assert_eq!(is_PentaFluoroPhenyl(&no_ring_neighbor, 0, 1), Ok(0));

        for mutate in 0..5 {
            let mut values = pentafluorophenyl();
            match mutate {
                0 => values[4].el_number = EL_NUMBER_O,
                1 => values[4].valence = 2,
                2 => values[4].num_H = 1,
                3 => values[4].charge = -1,
                4 => values[4].radical = 2,
                _ => unreachable!(),
            }
            assert_eq!(is_PentaFluoroPhenyl(&values, 0, 1), Ok(0));
        }

        for mutate in 0..5 {
            let mut values = pentafluorophenyl();
            match mutate {
                0 => values[9].el_number = EL_NUMBER_O,
                1 => values[9].chem_bonds_valence = 2,
                2 => values[9].charge = 1,
                3 => values[9].radical = 1,
                4 => values[9].num_H = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_PentaFluoroPhenyl(&values, 0, 1), Ok(0));
        }

        let mut no_continuation = pentafluorophenyl();
        no_continuation[4].neighbor = {
            let mut neighbors = [0; 20];
            neighbors[..3].copy_from_slice(&[3, 9, 3]);
            neighbors
        };
        assert_eq!(is_PentaFluoroPhenyl(&no_continuation, 0, 1), Ok(0));

        let mut not_closed = pentafluorophenyl();
        not_closed[6].neighbor[1] = 0;
        assert_eq!(is_PentaFluoroPhenyl(&not_closed, 0, 1), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_methyl__line_1679() {
        let methyl = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            el_number: EL_NUMBER_C,
            num_H: 3,
            ..inp_ATOM::default()
        };
        assert_eq!(is_Methyl(std::slice::from_ref(&methyl), 0), Ok(1));

        for mutate in 0..6 {
            let mut rejected = methyl.clone();
            match mutate {
                0 => rejected.valence = 2,
                1 => rejected.chem_bonds_valence = 2,
                2 => rejected.el_number = EL_NUMBER_O,
                3 => rejected.num_H = 2,
                4 => rejected.charge = -1,
                5 => rejected.radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_Methyl(std::slice::from_ref(&rejected), 0), Ok(0));
        }
    }

    #[test]
    fn source_port__ichinorm__is_ethyl__line_1694() {
        let mut values = vec![inp_ATOM::default(); 3];
        values[0].el_number = EL_NUMBER_C;
        values[1] = inp_ATOM {
            valence: 2,
            chem_bonds_valence: 2,
            el_number: EL_NUMBER_C,
            num_H: 2,
            neighbor: {
                let mut neighbors = [0; 20];
                neighbors[..2].copy_from_slice(&[0, 2]);
                neighbors
            },
            ..inp_ATOM::default()
        };
        values[2] = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            el_number: EL_NUMBER_C,
            num_H: 3,
            ..inp_ATOM::default()
        };
        assert_eq!(is_Ethyl(&values, 0, 1), Ok(1));
        assert_eq!(is_Ethyl(&values, 2, 1), Ok(0));

        let mut reversed = values.clone();
        reversed[1].neighbor[..2].copy_from_slice(&[2, 0]);
        assert_eq!(is_Ethyl(&reversed, 0, 1), Ok(1));

        for mutate in 0..6 {
            let mut rejected = values.clone();
            match mutate {
                0 => rejected[1].valence = 1,
                1 => rejected[1].chem_bonds_valence = 1,
                2 => rejected[1].el_number = EL_NUMBER_O,
                3 => rejected[1].num_H = 1,
                4 => rejected[1].charge = 1,
                5 => rejected[1].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_Ethyl(&rejected, 0, 1), Ok(0));
        }

        let mut bad_methyl = values.clone();
        bad_methyl[2].num_H = 2;
        assert_eq!(is_Ethyl(&bad_methyl, 0, 1), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_methyl_or_etyl__line_1709() {
        let methyl = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            el_number: EL_NUMBER_C,
            num_H: 3,
            ..inp_ATOM::default()
        };
        assert_eq!(
            is_Methyl_or_Etyl(std::slice::from_ref(&methyl), 0, 0),
            Ok(1)
        );

        let mut ethyl = vec![inp_ATOM::default(); 3];
        ethyl[1] = inp_ATOM {
            valence: 2,
            chem_bonds_valence: 2,
            el_number: EL_NUMBER_C,
            num_H: 2,
            neighbor: {
                let mut neighbors = [0; 20];
                neighbors[..2].copy_from_slice(&[0, 2]);
                neighbors
            },
            ..inp_ATOM::default()
        };
        ethyl[2] = methyl;
        assert_eq!(is_Methyl_or_Etyl(&ethyl, 0, 1), Ok(2));

        ethyl[2].num_H = 2;
        assert_eq!(is_Methyl_or_Etyl(&ethyl, 0, 1), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_si_iv__line_1727() {
        let silicon = inp_ATOM {
            el_number: EL_NUMBER_SI,
            valence: 4,
            chem_bonds_valence: 4,
            ..inp_ATOM::default()
        };
        assert_eq!(is_Si_IV(std::slice::from_ref(&silicon), 0), Ok(1));

        for mutate in 0..5 {
            let mut rejected = silicon.clone();
            match mutate {
                0 => rejected.el_number = EL_NUMBER_C,
                1 => rejected.charge = -1,
                2 => rejected.radical = 1,
                3 => rejected.valence = 3,
                4 => rejected.chem_bonds_valence = 5,
                _ => unreachable!(),
            }
            assert_eq!(is_Si_IV(std::slice::from_ref(&rejected), 0), Ok(0));
        }
    }

    #[test]
    fn source_port__ichinorm__is_deriv_ring_dmox_deox_o__line_1802() {
        fn dmox() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 9];
            let ring = [
                (EL_NUMBER_O, 2, 2, 0, [4, 1, 0, 0]),
                (EL_NUMBER_C, 2, 2, 2, [0, 2, 0, 0]),
                (EL_NUMBER_C, 4, 4, 0, [1, 3, 5, 6]),
                (EL_NUMBER_N, 2, 3, 0, [2, 4, 0, 0]),
                (EL_NUMBER_C, 3, 4, 0, [3, 0, 7, 0]),
            ];
            for (index, (element, valence, bond_valence, hydrogens, neighbors)) in
                ring.into_iter().enumerate()
            {
                values[index].el_number = element;
                values[index].valence = valence;
                values[index].chem_bonds_valence = bond_valence;
                values[index].num_H = hydrogens;
                values[index].neighbor[..usize::try_from(valence).unwrap()]
                    .copy_from_slice(&neighbors[..usize::try_from(valence).unwrap()]);
                values[index].bond_type[1] = if index == 3 {
                    BOND_DOUBLE as u8
                } else {
                    BOND_SINGLE as u8
                };
                values[index].nRingSystem = 1;
            }
            values[0].nNumAtInRingSystem = 5;
            for index in [5, 6] {
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 1,
                    chem_bonds_valence: 1,
                    num_H: 3,
                    neighbor: {
                        let mut neighbors = [0; 20];
                        neighbors[0] = 2;
                        neighbors
                    },
                    ..inp_ATOM::default()
                };
            }
            values[7].el_number = EL_NUMBER_C;
            values
        }

        let mut da = DERIV_AT {
            typ: [91; 4],
            ord: [92; 4],
            num: [93; 4],
            other_atom: 94,
        };
        let original_da = da.clone();
        let mut output = DERIV_AT {
            typ: [11; 4],
            ord: [12; 4],
            num: [13; 4],
            other_atom: 14,
        };
        assert_eq!(
            is_DERIV_RING_DMOX_DEOX_O(&dmox(), 0, 0, Some(&mut da), Some(&mut output)),
            Ok(DERIV_RING_DMOX_DEOX_O as i32)
        );
        assert_eq!(da, original_da);
        assert_eq!(output.typ, [DERIV_RING_DMOX_DEOX_O as i16, 11, 11, 11]);
        assert_eq!(output.ord, [1, 12, 12, 12]);
        assert_eq!(output.num, [4, 13, 13, 13]);
        assert_eq!(output.other_atom, 4);
        assert_eq!(
            is_DERIV_RING_DMOX_DEOX_O(&dmox(), 0, 0, None, None),
            Ok(DERIV_RING_DMOX_DEOX_O as i32)
        );

        let mut no_r = dmox();
        no_r[4].neighbor[2] = 3;
        assert_eq!(
            is_DERIV_RING_DMOX_DEOX_O(&no_r, 0, 0, None, None),
            Ok(DERIV_RING_DMOX_DEOX_O as i32)
        );

        let mut failures = Vec::new();
        let mut value = dmox();
        value[0].el_number = EL_NUMBER_C;
        failures.push(value);
        let mut value = dmox();
        value[0].nNumAtInRingSystem = 4;
        failures.push(value);
        let mut value = dmox();
        value[1].nRingSystem = 2;
        failures.push(value);
        for mutation in 0..7 {
            let mut value = dmox();
            match mutation {
                0 => value[2].charge = 1,
                1 => value[2].radical = 1,
                2 => value[2].bond_type[1] = BOND_DOUBLE as u8,
                3 => value[2].valence = 3,
                4 => value[2].chem_bonds_valence = 3,
                5 => value[2].num_H = 1,
                6 => value[2].el_number = EL_NUMBER_N,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        let mut value = dmox();
        value[8].nRingSystem = 1;
        value[4].neighbor[1] = 8;
        failures.push(value);
        let mut value = dmox();
        value[7].el_number = EL_NUMBER_O;
        failures.push(value);
        let mut value = dmox();
        value[6].num_H = 2;
        failures.push(value);
        let mut value = dmox();
        value[6] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 2,
            chem_bonds_valence: 2,
            num_H: 2,
            neighbor: {
                let mut neighbors = [0; 20];
                neighbors[..2].copy_from_slice(&[2, 8]);
                neighbors
            },
            ..inp_ATOM::default()
        };
        value[8] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 1,
            chem_bonds_valence: 1,
            num_H: 3,
            ..inp_ATOM::default()
        };
        failures.push(value);

        for atoms in failures {
            let mut untouched = DERIV_AT {
                typ: [21; 4],
                ord: [22; 4],
                num: [23; 4],
                other_atom: 24,
            };
            let expected = untouched.clone();
            assert_eq!(
                is_DERIV_RING_DMOX_DEOX_O(&atoms, 0, 0, None, Some(&mut untouched)),
                Ok(0)
            );
            assert_eq!(untouched, expected);
        }
    }

    #[test]
    fn source_port__ichinorm__is_deriv_ring_dmox_deox_n__line_1962() {
        fn dmox_from_nitrogen() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 9];
            let ring = [
                (EL_NUMBER_N, 2, 3, 0, [4, 1, 0, 0]),
                (EL_NUMBER_C, 4, 4, 0, [0, 2, 5, 6]),
                (EL_NUMBER_C, 2, 2, 2, [1, 3, 0, 0]),
                (EL_NUMBER_O, 2, 2, 0, [2, 4, 0, 0]),
                (EL_NUMBER_C, 3, 4, 0, [3, 0, 7, 0]),
            ];
            for (index, (element, valence, bond_valence, hydrogens, neighbors)) in
                ring.into_iter().enumerate()
            {
                values[index].el_number = element;
                values[index].valence = valence;
                values[index].chem_bonds_valence = bond_valence;
                values[index].num_H = hydrogens;
                values[index].neighbor[..usize::try_from(valence).unwrap()]
                    .copy_from_slice(&neighbors[..usize::try_from(valence).unwrap()]);
                values[index].bond_type[1] = if index == 4 {
                    BOND_DOUBLE as u8
                } else {
                    BOND_SINGLE as u8
                };
                values[index].nRingSystem = 1;
            }
            values[0].nNumAtInRingSystem = 5;
            for index in [5, 6] {
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 1,
                    chem_bonds_valence: 1,
                    num_H: 3,
                    neighbor: {
                        let mut neighbors = [0; 20];
                        neighbors[0] = 1;
                        neighbors
                    },
                    ..inp_ATOM::default()
                };
            }
            values[7].el_number = EL_NUMBER_C;
            values
        }

        let mut da = DERIV_AT {
            typ: [31; 4],
            ord: [32; 4],
            num: [33; 4],
            other_atom: 34,
        };
        let expected_da = da.clone();
        let mut output = DERIV_AT {
            typ: [41; 4],
            ord: [42; 4],
            num: [43; 4],
            other_atom: 44,
        };
        assert_eq!(
            is_DERIV_RING_DMOX_DEOX_N(
                &dmox_from_nitrogen(),
                0,
                0,
                Some(&mut da),
                Some(&mut output),
            ),
            Ok(DERIV_RING_DMOX_DEOX_N as i32)
        );
        assert_eq!(da, expected_da);
        assert_eq!(output.typ, [DERIV_RING_DMOX_DEOX_N as i16, 41, 41, 41]);
        assert_eq!(output.ord, [1, 42, 42, 42]);
        assert_eq!(output.num, [4, 43, 43, 43]);
        assert_eq!(output.other_atom, 4);
        assert_eq!(
            is_DERIV_RING_DMOX_DEOX_N(&dmox_from_nitrogen(), 0, 0, None, None),
            Ok(DERIV_RING_DMOX_DEOX_N as i32)
        );

        let mut no_r = dmox_from_nitrogen();
        no_r[4].neighbor[2] = 3;
        assert_eq!(
            is_DERIV_RING_DMOX_DEOX_N(&no_r, 0, 0, None, None),
            Ok(DERIV_RING_DMOX_DEOX_N as i32)
        );

        let mut failures = Vec::new();
        for mutation in 0..4 {
            let mut value = dmox_from_nitrogen();
            match mutation {
                0 => value[0].el_number = EL_NUMBER_O,
                1 => value[0].nNumAtInRingSystem = 4,
                2 => value[0].valence = 1,
                3 => value[0].chem_bonds_valence = 2,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        let mut value = dmox_from_nitrogen();
        value[1].nRingSystem = 2;
        failures.push(value);
        for mutation in 0..7 {
            let mut value = dmox_from_nitrogen();
            match mutation {
                0 => value[2].charge = 1,
                1 => value[2].radical = 1,
                2 => value[2].bond_type[1] = BOND_DOUBLE as u8,
                3 => value[2].valence = 3,
                4 => value[2].chem_bonds_valence = 3,
                5 => value[2].num_H = 1,
                6 => value[2].el_number = EL_NUMBER_N,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        let mut value = dmox_from_nitrogen();
        value[8].nRingSystem = 1;
        value[4].neighbor[1] = 8;
        failures.push(value);
        let mut value = dmox_from_nitrogen();
        value[7].el_number = EL_NUMBER_O;
        failures.push(value);
        let mut value = dmox_from_nitrogen();
        value[6].num_H = 2;
        failures.push(value);
        let mut value = dmox_from_nitrogen();
        value[6] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 2,
            chem_bonds_valence: 2,
            num_H: 2,
            neighbor: {
                let mut neighbors = [0; 20];
                neighbors[..2].copy_from_slice(&[1, 8]);
                neighbors
            },
            ..inp_ATOM::default()
        };
        value[8] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 1,
            chem_bonds_valence: 1,
            num_H: 3,
            ..inp_ATOM::default()
        };
        failures.push(value);

        for atoms in failures {
            let mut untouched = DERIV_AT {
                typ: [51; 4],
                ord: [52; 4],
                num: [53; 4],
                other_atom: 54,
            };
            let expected = untouched.clone();
            assert_eq!(
                is_DERIV_RING_DMOX_DEOX_N(&atoms, 0, 0, None, Some(&mut untouched)),
                Ok(0)
            );
            assert_eq!(untouched, expected);
        }
    }

    #[test]
    fn source_port__ichinorm__is_deriv_ring2_prrldd_pprdn__line_2147() {
        fn ring(ring_size: usize, from_order: usize) -> Vec<inp_ATOM> {
            let substituent = ring_size + 1;
            let oxygen = ring_size + 2;
            let mut values = vec![inp_ATOM::default(); ring_size + 3];
            let mut current_neighbors = [0; 20];
            let ring_start = 2_u16;
            let ring_end = ring_size as u16;
            let ring_orders: Vec<usize> = (0..3).filter(|&order| order != from_order).collect();
            current_neighbors[from_order] = 1;
            current_neighbors[ring_orders[0]] = ring_start;
            current_neighbors[ring_orders[1]] = ring_end;
            values[0] = inp_ATOM {
                el_number: EL_NUMBER_N,
                valence: 3,
                chem_bonds_valence: 3,
                neighbor: current_neighbors,
                nRingSystem: 1,
                nNumAtInRingSystem: ring_size as u16,
                ..inp_ATOM::default()
            };
            values[1] = inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 3,
                chem_bonds_valence: 4,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..3].copy_from_slice(&[0, substituent as u16, oxygen as u16]);
                    neighbors
                },
                bond_type: {
                    let mut bonds = [0; 20];
                    bonds[..3].copy_from_slice(&[
                        BOND_SINGLE as u8,
                        BOND_SINGLE as u8,
                        BOND_DOUBLE as u8,
                    ]);
                    bonds
                },
                nNumAtInRingSystem: 1,
                ..inp_ATOM::default()
            };
            for index in 2..=ring_size {
                let previous = if index == 2 { 0 } else { index - 1 };
                let next = if index == ring_size { 0 } else { index + 1 };
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 2,
                    chem_bonds_valence: 2,
                    num_H: 2,
                    neighbor: {
                        let mut neighbors = [0; 20];
                        neighbors[..2].copy_from_slice(&[previous as u16, next as u16]);
                        neighbors
                    },
                    nRingSystem: 1,
                    ..inp_ATOM::default()
                };
            }
            values[substituent].el_number = EL_NUMBER_C;
            values[oxygen] = inp_ATOM {
                el_number: EL_NUMBER_O,
                valence: 1,
                chem_bonds_valence: 2,
                ..inp_ATOM::default()
            };
            values
        }

        for (size, derivative_type) in [
            (5, DERIV_RING2_PRRLDD_OUTSIDE_PRECUR),
            (6, DERIV_RING2_PPRDN_OUTSIDE_PRECUR),
        ] {
            for from_order in 0..3 {
                let mut output = DERIV_AT {
                    typ: [61; 4],
                    ord: [62; 4],
                    num: [63; 4],
                    other_atom: 64,
                };
                assert_eq!(
                    is_DERIV_RING2_PRRLDD_PPRDN(
                        &ring(size, from_order),
                        0,
                        from_order as i32,
                        None,
                        Some(&mut output),
                    ),
                    Ok(size as i32)
                );
                assert_eq!(
                    output.typ,
                    [derivative_type as i16, derivative_type as i16, 61, 61]
                );
                let mut expected_orders: Vec<i8> = (0..3)
                    .filter(|&order| order != from_order)
                    .map(|v| v as i8)
                    .collect();
                expected_orders.sort_unstable();
                assert_eq!(output.ord[..2], expected_orders);
                assert_eq!(output.ord[2..], [62, 62]);
                assert_eq!(output.num, [(size - 1) as i8, 63, 63, 63]);
                assert_eq!(output.other_atom, 64);
            }
            assert_eq!(
                is_DERIV_RING2_PRRLDD_PPRDN(&ring(size, 0), 0, 0, None, None),
                Ok(size as i32)
            );
        }

        let mut failures = Vec::new();
        for mutation in 0..14 {
            let mut value = ring(5, 0);
            match mutation {
                0 => value[0].el_number = EL_NUMBER_C,
                1 => value[0].nNumAtInRingSystem = 4,
                2 => value[0].nNumAtInRingSystem = 7,
                3 => value[0].valence = 2,
                4 => value[0].chem_bonds_valence = 2,
                5 => value[0].charge = 1,
                6 => value[0].radical = 1,
                7 => value[0].num_H = 1,
                8 => value[1].el_number = EL_NUMBER_N,
                9 => value[1].nNumAtInRingSystem = 2,
                10 => value[1].valence = 2,
                11 => value[1].chem_bonds_valence = 3,
                12 => value[1].charge = 1,
                13 => value[1].num_H = 1,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        for mutation in 0..8 {
            let mut value = ring(5, 0);
            let oxygen = 7;
            match mutation {
                0 => value[6].el_number = EL_NUMBER_O,
                1 => value[1].bond_type[1] = BOND_DOUBLE as u8,
                2 => value[1].bond_type[2] = BOND_SINGLE as u8,
                3 => value[oxygen].el_number = EL_NUMBER_C,
                4 => value[oxygen].valence = 2,
                5 => value[oxygen].chem_bonds_valence = 1,
                6 => value[oxygen].charge = 1,
                7 => value[oxygen].num_H = 1,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        for mutation in 0..7 {
            let mut value = ring(5, 0);
            match mutation {
                0 => value[2].nRingSystem = 2,
                1 => value[3].el_number = EL_NUMBER_N,
                2 => value[3].valence = 1,
                3 => value[3].chem_bonds_valence = 1,
                4 => value[3].num_H = 1,
                5 => value[3].charge = 1,
                6 => value[3].radical = 1,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        let mut open = ring(5, 0);
        open[5].neighbor[1] = 6;
        failures.push(open);

        for atoms in failures {
            let mut output = DERIV_AT {
                typ: [71; 4],
                ord: [72; 4],
                num: [73; 4],
                other_atom: 74,
            };
            let expected = output.clone();
            assert_eq!(
                is_DERIV_RING2_PRRLDD_PPRDN(&atoms, 0, 0, None, Some(&mut output)),
                Ok(0)
            );
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn source_port__ichinorm__check_arom_chain__line_2281() {
        fn aromatic_chain() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 4];
            values[1] = inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 2,
                chem_bonds_valence: 3,
                num_H: 1,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..2].copy_from_slice(&[0, 2]);
                    neighbors
                },
                bond_type: {
                    let mut bonds = [0; 20];
                    bonds[1] = BOND_TYPE_ALTERN as u8;
                    bonds
                },
                ..inp_ATOM::default()
            };
            values[2] = inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 2,
                chem_bonds_valence: 3,
                num_H: 1,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..2].copy_from_slice(&[3, 1]);
                    neighbors
                },
                bond_type: {
                    let mut bonds = [0; 20];
                    bonds[0] = BOND_TYPE_ALTERN as u8;
                    bonds
                },
                ..inp_ATOM::default()
            };
            values[3] = inp_ATOM {
                el_number: EL_NUMBER_O,
                valence: 19,
                chem_bonds_valence: -7,
                num_H: 12,
                ..inp_ATOM::default()
            };
            values
        }

        assert_eq!(check_arom_chain(&aromatic_chain(), 1, 0, 3, 3), Ok(1));
        assert_eq!(check_arom_chain(&aromatic_chain(), 1, 0, 3, 2), Ok(0));
        assert_eq!(check_arom_chain(&aromatic_chain(), 1, 0, 3, 4), Ok(0));

        for mutation in 0..5 {
            let mut values = aromatic_chain();
            match mutation {
                0 => values[2].el_number = EL_NUMBER_N,
                1 => values[2].valence = 3,
                2 => values[2].chem_bonds_valence = 2,
                3 => values[2].num_H = 0,
                4 => values[2].bond_type[0] = BOND_SINGLE as u8,
                _ => unreachable!(),
            }
            assert_eq!(check_arom_chain(&values, 1, 0, 3, 3), Ok(0));
        }

        let mut misses_last = aromatic_chain();
        misses_last[2].neighbor[0] = 0;
        assert_eq!(check_arom_chain(&misses_last, 1, 0, 3, 2), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_dansyl__line_2334() {
        fn dansyl() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 18];
            values[0] = inp_ATOM {
                el_number: EL_NUMBER_O,
                valence: 2,
                chem_bonds_valence: 2,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..2].copy_from_slice(&[17, 1]);
                    neighbors
                },
                nNumAtInRingSystem: 1,
                ..inp_ATOM::default()
            };
            values[1] = inp_ATOM {
                el_number: EL_NUMBER_S,
                valence: 4,
                chem_bonds_valence: 6,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..4].copy_from_slice(&[0, 2, 3, 4]);
                    neighbors
                },
                bond_type: {
                    let mut bonds = [0; 20];
                    bonds[..4].copy_from_slice(&[
                        BOND_SINGLE as u8,
                        BOND_DOUBLE as u8,
                        BOND_DOUBLE as u8,
                        BOND_SINGLE as u8,
                    ]);
                    bonds
                },
                ..inp_ATOM::default()
            };
            for index in [2, 3] {
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_O,
                    valence: 1,
                    chem_bonds_valence: 2,
                    ..inp_ATOM::default()
                };
            }
            for (index, neighbors) in [
                (4, [1, 5, 16]),
                (5, [4, 6, 13]),
                (6, [5, 7, 14]),
                (7, [6, 8, 11]),
            ] {
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 3,
                    chem_bonds_valence: 4,
                    neighbor: {
                        let mut value = [0; 20];
                        value[..3].copy_from_slice(&neighbors);
                        value
                    },
                    bond_type: {
                        let mut value = [0; 20];
                        value[..3].fill(BOND_TYPE_ALTERN as u8);
                        value
                    },
                    ..inp_ATOM::default()
                };
            }
            values[4].bond_type[0] = BOND_SINGLE as u8;
            values[4].nNumAtInRingSystem = 10;
            values[7].bond_type[1] = BOND_SINGLE as u8;
            values[8] = inp_ATOM {
                el_number: EL_NUMBER_N,
                valence: 3,
                chem_bonds_valence: 3,
                neighbor: {
                    let mut value = [0; 20];
                    value[..3].copy_from_slice(&[7, 9, 10]);
                    value
                },
                ..inp_ATOM::default()
            };
            for index in [9, 10] {
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 1,
                    chem_bonds_valence: 1,
                    num_H: 3,
                    ..inp_ATOM::default()
                };
            }
            for (index, neighbors) in [
                (11, [7, 12]),
                (12, [11, 13]),
                (13, [12, 5]),
                (14, [6, 15]),
                (15, [14, 16]),
                (16, [15, 4]),
            ] {
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 2,
                    chem_bonds_valence: 3,
                    num_H: 1,
                    neighbor: {
                        let mut value = [0; 20];
                        value[..2].copy_from_slice(&neighbors);
                        value
                    },
                    bond_type: [BOND_TYPE_ALTERN as u8; 20],
                    ..inp_ATOM::default()
                };
            }
            values
        }

        let mut da = DERIV_AT {
            typ: [81; 4],
            ord: [82; 4],
            num: [83; 4],
            other_atom: 84,
        };
        let expected_da = da.clone();
        let mut output = DERIV_AT {
            typ: [85; 4],
            ord: [86; 4],
            num: [87; 4],
            other_atom: 88,
        };
        assert_eq!(
            is_Dansyl(&dansyl(), 0, 1, Some(&mut da), Some(&mut output)),
            Ok(DERIV_DANSYL as i32)
        );
        assert_eq!(da, expected_da);
        assert_eq!(output.typ, [DERIV_DANSYL as i16, 85, 85, 85]);
        assert_eq!(output.ord, [1, 86, 86, 86]);
        assert_eq!(output.num, [16, 87, 87, 87]);
        assert_eq!(output.other_atom, 88);
        assert_eq!(
            is_Dansyl(&dansyl(), 0, 1, None, None),
            Ok(DERIV_DANSYL as i32)
        );

        for (element, valence, hydrogens, ring_size) in [
            (EL_NUMBER_O, 2, 0, 1),
            (EL_NUMBER_S, 2, 0, 1),
            (EL_NUMBER_N, 2, 1, 1),
            (EL_NUMBER_N, 3, 0, 0),
        ] {
            let mut values = dansyl();
            values[0].el_number = element;
            values[0].valence = valence;
            values[0].chem_bonds_valence = valence;
            values[0].num_H = hydrogens;
            values[0].nNumAtInRingSystem = ring_size;
            assert_eq!(
                is_Dansyl(&values, 0, 1, None, None),
                Ok(DERIV_DANSYL as i32)
            );
        }
        let mut unchecked_fields = dansyl();
        unchecked_fields[0].charge = 3;
        unchecked_fields[0].radical = 2;
        unchecked_fields[8].radical = 4;
        unchecked_fields[8].num_H = 5;
        assert_eq!(
            is_Dansyl(&unchecked_fields, 0, 1, None, None),
            Ok(DERIV_DANSYL as i32)
        );

        let mut failures = Vec::new();
        for mutation in 0..8 {
            let mut value = dansyl();
            match mutation {
                0 => value[0].el_number = EL_NUMBER_C,
                1 => value[0].chem_bonds_valence = 1,
                2 => value[1].el_number = EL_NUMBER_O,
                3 => value[1].valence = 3,
                4 => value[1].chem_bonds_valence = 5,
                5 => value[2].el_number = EL_NUMBER_C,
                6 => value[2].valence = 2,
                7 => value[2].charge = 1,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        for mutation in 0..9 {
            let mut value = dansyl();
            match mutation {
                0 => value[4].el_number = EL_NUMBER_N,
                1 => value[4].nNumAtInRingSystem = 9,
                2 => value[4].valence = 2,
                3 => value[4].chem_bonds_valence = 3,
                4 => value[4].num_H = 1,
                5 => value[4].charge = 1,
                6 => value[4].radical = 1,
                7 => value[4].bond_type[1] = BOND_SINGLE as u8,
                8 => value[1].bond_type[2] = BOND_SINGLE as u8,
                _ => unreachable!(),
            }
            failures.push(value);
        }
        for (index, mutation) in [(5, 0), (5, 1), (6, 0), (6, 1), (7, 0), (7, 1)] {
            let mut value = dansyl();
            if mutation == 0 {
                value[index].el_number = EL_NUMBER_N;
            } else {
                value[index].chem_bonds_valence = 3;
            }
            failures.push(value);
        }
        let mut value = dansyl();
        value[8].el_number = EL_NUMBER_C;
        failures.push(value);
        let mut value = dansyl();
        value[8].charge = 1;
        failures.push(value);
        let mut value = dansyl();
        value[9].num_H = 2;
        failures.push(value);
        let mut value = dansyl();
        value[12].bond_type[1] = BOND_SINGLE as u8;
        failures.push(value);
        let mut value = dansyl();
        value[15].num_H = 0;
        failures.push(value);

        for atoms in failures {
            let mut untouched = DERIV_AT {
                typ: [91; 4],
                ord: [92; 4],
                num: [93; 4],
                other_atom: 94,
            };
            let expected = untouched.clone();
            assert_eq!(is_Dansyl(&atoms, 0, 1, None, Some(&mut untouched)), Ok(0));
            assert_eq!(untouched, expected);
        }
    }

    #[test]
    fn source_port__ichinorm__is_silyl__line_3633() {
        fn methyl(attached_to: u16) -> inp_ATOM {
            inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 1,
                chem_bonds_valence: 1,
                num_H: 3,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[0] = attached_to;
                    neighbors
                },
                ..inp_ATOM::default()
            }
        }
        fn trimethylsilyl() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 5];
            values[0] = inp_ATOM {
                el_number: EL_NUMBER_SI,
                valence: 4,
                chem_bonds_valence: 4,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..4].copy_from_slice(&[1, 2, 3, 4]);
                    neighbors
                },
                ..inp_ATOM::default()
            };
            for index in 2..=4 {
                values[index] = methyl(0);
            }
            values
        }
        fn tert_butyl_dimethylsilyl() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 8];
            values[0] = inp_ATOM {
                el_number: EL_NUMBER_SI,
                valence: 4,
                chem_bonds_valence: 4,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..4].copy_from_slice(&[1, 2, 3, 4]);
                    neighbors
                },
                ..inp_ATOM::default()
            };
            values[2] = methyl(0);
            values[3] = methyl(0);
            values[4] = inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 4,
                chem_bonds_valence: 4,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..4].copy_from_slice(&[0, 5, 6, 7]);
                    neighbors
                },
                ..inp_ATOM::default()
            };
            for index in 5..=7 {
                values[index] = methyl(4);
            }
            values
        }

        assert_eq!(is_silyl(&trimethylsilyl(), 0, 0), Ok(1));
        assert_eq!(is_silyl(&tert_butyl_dimethylsilyl(), 0, 0), Ok(2));

        for mutation in 0..5 {
            let mut values = trimethylsilyl();
            match mutation {
                0 => values[0].el_number = EL_NUMBER_C,
                1 => values[0].valence = 3,
                2 => values[0].chem_bonds_valence = 5,
                3 => values[0].charge = 1,
                4 => values[0].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_silyl(&values, 0, 0), Ok(0));
        }
        for mutation in 0..6 {
            let mut values = trimethylsilyl();
            match mutation {
                0 => values[3].charge = 1,
                1 => values[3].radical = 1,
                2 => values[3].chem_bonds_valence = 2,
                3 => values[3].valence = 2,
                4 => values[3].el_number = EL_NUMBER_O,
                5 => values[3].num_H = 2,
                _ => unreachable!(),
            }
            assert_eq!(is_silyl(&values, 0, 0), Ok(0));
        }

        let mut disilyl = tert_butyl_dimethylsilyl();
        disilyl[4].el_number = EL_NUMBER_SI;
        assert_eq!(is_silyl(&disilyl, 0, 0), Ok(0));
        let mut duplicate_center = tert_butyl_dimethylsilyl();
        duplicate_center[3] = duplicate_center[4].clone();
        assert_eq!(is_silyl(&duplicate_center, 0, 0), Ok(0));

        for mutation in 0..6 {
            let mut values = tert_butyl_dimethylsilyl();
            match mutation {
                0 => values[6].charge = 1,
                1 => values[6].radical = 1,
                2 => values[6].chem_bonds_valence = 2,
                3 => values[6].valence = 2,
                4 => values[6].el_number = EL_NUMBER_O,
                5 => values[6].num_H = 2,
                _ => unreachable!(),
            }
            assert_eq!(is_silyl(&values, 0, 0), Ok(0));
        }
    }

    #[test]
    fn source_port__ichinorm__is_cf3_or_linc3f7__line_4017() {
        fn perfluoro_chain(carbon_count: usize) -> Vec<inp_ATOM> {
            let external = carbon_count;
            let mut values = vec![inp_ATOM::default(); carbon_count + 1];
            let mut next_fluorine = values.len();
            for carbon in 0..carbon_count {
                let previous = if carbon == 0 { external } else { carbon - 1 };
                let mut neighbors = vec![previous as u16];
                if carbon + 1 < carbon_count {
                    neighbors.push((carbon + 1) as u16);
                }
                while neighbors.len() < 4 {
                    neighbors.push(next_fluorine as u16);
                    values.push(inp_ATOM {
                        el_number: EL_NUMBER_F,
                        valence: 1,
                        chem_bonds_valence: 1,
                        ..inp_ATOM::default()
                    });
                    next_fluorine += 1;
                }
                values[carbon] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 4,
                    chem_bonds_valence: 4,
                    neighbor: {
                        let mut value = [0; 20];
                        value[..4].copy_from_slice(&neighbors);
                        value
                    },
                    ..inp_ATOM::default()
                };
            }
            values
        }

        for count in 1..=4 {
            assert_eq!(
                is_CF3_or_linC3F7(&perfluoro_chain(count), 0, 0),
                Ok(if count < 4 { count as i32 } else { 0 })
            );
        }

        let mut missing_backlink = perfluoro_chain(2);
        missing_backlink[1].neighbor[0] = 2;
        assert_eq!(is_CF3_or_linC3F7(&missing_backlink, 0, 0), Ok(-1));

        for mutation in 0..6 {
            let mut values = perfluoro_chain(1);
            match mutation {
                0 => values[0].el_number = EL_NUMBER_N,
                1 => values[0].valence = 3,
                2 => values[0].chem_bonds_valence = 3,
                3 => values[0].num_H = 1,
                4 => values[0].charge = 1,
                5 => values[0].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_CF3_or_linC3F7(&values, 0, 0), Ok(0));
        }

        for mutation in 0..5 {
            let mut values = perfluoro_chain(1);
            match mutation {
                0 => values[2].valence = 2,
                1 => values[2].chem_bonds_valence = 0,
                2 => values[2].num_H = 1,
                3 => values[2].charge = 1,
                4 => values[2].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_CF3_or_linC3F7(&values, 0, 0), Ok(0));
        }

        let mut ignored = perfluoro_chain(1);
        ignored[2] = inp_ATOM {
            el_number: EL_NUMBER_O,
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        assert_eq!(is_CF3_or_linC3F7(&ignored, 0, 0), Ok(0));

        let mut two_carbons = perfluoro_chain(1);
        for index in [2, 3] {
            two_carbons[index] = inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 4,
                chem_bonds_valence: 4,
                ..inp_ATOM::default()
            };
        }
        assert_eq!(is_CF3_or_linC3F7(&two_carbons, 0, 0), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__underiv_buf_clear__line_5147() {
        assert_eq!(underiv_buf_clear(None), Ok(()));
        let mut buffer = [b'a' as i8, b'b' as i8, 0, b'c' as i8];
        assert_eq!(underiv_buf_clear(Some(&mut buffer)), Ok(()));
        assert_eq!(buffer, [0, b'b' as i8, 0, b'c' as i8]);
        assert_eq!(
            underiv_buf_clear(Some(&mut [])),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__underiv_list_add__line_5156() {
        fn c_bytes(bytes: &[u8]) -> Vec<i8> {
            bytes.iter().map(|&byte| byte as i8).collect()
        }

        let added = c_bytes(b"cd\0tail");
        assert_eq!(underiv_list_add(None, 9, Some(&added), b' ' as i8), Ok(0));
        let mut list = c_bytes(b"ab\0xxxxx");
        assert_eq!(
            underiv_list_add(Some(&mut list), 9, None, b' ' as i8),
            Ok(0)
        );
        assert_eq!(list, c_bytes(b"ab\0xxxxx"));
        assert_eq!(
            underiv_list_add(Some(&mut list), 0, Some(&added), b' ' as i8),
            Ok(0)
        );

        let empty = c_bytes(b"\0ignored");
        assert_eq!(
            underiv_list_add(Some(&mut list), 9, Some(&empty), b' ' as i8),
            Ok(0)
        );
        assert_eq!(list, c_bytes(b"ab\0xxxxx"));

        let mut delimited = c_bytes(b"ab\0zzzzz");
        assert_eq!(
            underiv_list_add(Some(&mut delimited), 6, Some(&added), b' ' as i8),
            Ok(5)
        );
        assert_eq!(delimited, c_bytes(b"ab cd\0zz"));

        let mut exact = c_bytes(b"ab\0zzzzz");
        assert_eq!(
            underiv_list_add(Some(&mut exact), 5, Some(&added), b' ' as i8),
            Ok(0)
        );
        assert_eq!(exact, c_bytes(b"ab\0zzzzz"));

        let mut plain = c_bytes(b"ab\0zzzzz");
        assert_eq!(
            underiv_list_add(Some(&mut plain), 5, Some(&added), 0),
            Ok(4)
        );
        assert_eq!(plain, c_bytes(b"abcd\0zzz"));
        let mut plain_exact = c_bytes(b"ab\0zzzzz");
        assert_eq!(
            underiv_list_add(Some(&mut plain_exact), 4, Some(&added), 0),
            Ok(0)
        );

        let mut empty_list = c_bytes(b"\0zzzz");
        assert_eq!(
            underiv_list_add(Some(&mut empty_list), 4, Some(&added), b' ' as i8),
            Ok(2)
        );
        assert_eq!(empty_list, c_bytes(b"cd\0zz"));
        let mut empty_exact = c_bytes(b"\0zzzz");
        assert_eq!(
            underiv_list_add(Some(&mut empty_exact), 3, Some(&added), b' ' as i8),
            Ok(0)
        );

        let mut unterminated = c_bytes(b"abc");
        assert_eq!(
            underiv_list_add(Some(&mut unterminated), 8, Some(&added), 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__underiv_list_get_last__line_5179() {
        fn c_bytes(bytes: &[u8]) -> Vec<i8> {
            bytes.iter().map(|&byte| byte as i8).collect()
        }

        let list = c_bytes(b"one two three\0 ignored later");
        assert_eq!(underiv_list_get_last(Some(&list), b' ' as i8), Ok(Some(7)));
        assert_eq!(underiv_list_get_last(Some(&list), b'/' as i8), Ok(Some(0)));
        assert_eq!(underiv_list_get_last(Some(&list), 0), Ok(None));
        assert_eq!(underiv_list_get_last(None, b' ' as i8), Ok(None));

        let leading = c_bytes(b" alpha\0");
        assert_eq!(
            underiv_list_get_last(Some(&leading), b' ' as i8),
            Ok(Some(0))
        );
        let trailing = c_bytes(b"alpha \0");
        assert_eq!(
            underiv_list_get_last(Some(&trailing), b' ' as i8),
            Ok(Some(5))
        );
        let after_nul = c_bytes(b"alpha\0 beta");
        assert_eq!(
            underiv_list_get_last(Some(&after_nul), b' ' as i8),
            Ok(Some(0))
        );

        let unterminated = c_bytes(b"alpha beta");
        assert_eq!(
            underiv_list_get_last(Some(&unterminated), b' ' as i8),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__is_deriv_chain2__line_4483() {
        fn c_text(bytes: &[i8]) -> String {
            let end = bytes.iter().position(|&byte| byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|&byte| byte as u8 as char)
                .collect()
        }

        fn connected_pair(element: u8, valence: i8, chemical: i8) -> Vec<inp_ATOM> {
            let mut start = atom(&[1]);
            start.el_number = element;
            start.valence = valence;
            start.chem_bonds_valence = chemical;
            start.elname[..2].copy_from_slice(&[element_name(element), 0]);
            let mut neighbor = atom(&[0]);
            neighbor.chem_bonds_valence = 1;
            vec![start, neighbor]
        }

        fn element_name(element: u8) -> i8 {
            match element {
                EL_NUMBER_O => b'O' as i8,
                EL_NUMBER_N => b'N' as i8,
                EL_NUMBER_S => b'S' as i8,
                _ => b'C' as i8,
            }
        }

        let pair = connected_pair(EL_NUMBER_O, 1, 1);
        assert_eq!(
            is_deriv_chain2(&pair, 0, 0, 0, 0, 0, None, 0, None, 0, None),
            Ok(0)
        );
        assert_eq!(
            is_deriv_chain2(
                &pair,
                0,
                DERIV_RING_OUTSIDE_PRECURSOR as i32,
                0,
                0,
                0,
                None,
                0,
                None,
                0,
                None,
            ),
            Ok(0)
        );
        let mut charged = pair.clone();
        charged[0].charge = 1;
        assert_eq!(
            is_deriv_chain2(
                &charged,
                0,
                DERIV_DANSYL as i32,
                0,
                0,
                0,
                None,
                0,
                None,
                0,
                None,
            ),
            Ok(0)
        );
        let mut missing_reverse = pair.clone();
        missing_reverse[1].neighbor[0] = 1;
        assert_eq!(
            is_deriv_chain2(
                &missing_reverse,
                0,
                DERIV_DANSYL as i32,
                0,
                0,
                0,
                None,
                0,
                None,
                0,
                None,
            ),
            Ok(-1)
        );

        let mut silyl = vec![inp_ATOM::default(); 5];
        silyl[0] = atom(&[1]);
        silyl[0].el_number = EL_NUMBER_O;
        silyl[0].chem_bonds_valence = 1;
        silyl[0].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
        silyl[1] = atom(&[0, 2, 3, 4]);
        silyl[1].el_number = EL_NUMBER_SI;
        silyl[1].chem_bonds_valence = 4;
        for methyl in &mut silyl[2..] {
            methyl.el_number = EL_NUMBER_C;
            methyl.valence = 1;
            methyl.chem_bonds_valence = 1;
            methyl.num_H = 3;
        }
        let mut primary = [0_i8; 64];
        let mut secondary = [0_i8; 64];
        let mut bits = 0;
        assert_eq!(
            is_deriv_chain2(
                &silyl,
                0,
                DERIV_BRIDGE_O as i32,
                0,
                0,
                0,
                Some(&mut primary),
                64,
                Some(&mut secondary),
                64,
                Some(&mut bits),
            ),
            Ok(1)
        );
        assert_eq!(c_text(&primary), "RO-TMS");
        assert_eq!(c_text(&secondary), "TMS");
        assert_eq!(bits, DERIV_BIT_TMS as i32);

        let mut amide = vec![inp_ATOM::default(); 7];
        amide[0] = atom(&[1]);
        amide[0].el_number = EL_NUMBER_N;
        amide[0].chem_bonds_valence = 1;
        amide[1] = atom(&[0, 2, 3]);
        amide[1].el_number = EL_NUMBER_C;
        amide[1].chem_bonds_valence = 4;
        amide[2] = atom(&[1]);
        amide[2].el_number = EL_NUMBER_O;
        amide[2].chem_bonds_valence = 2;
        amide[3] = atom(&[1, 4, 5, 6]);
        amide[3].el_number = EL_NUMBER_C;
        amide[3].chem_bonds_valence = 4;
        for fluorine in 4..=6 {
            amide[fluorine] = atom(&[3]);
            amide[fluorine].el_number = EL_NUMBER_F;
            amide[fluorine].chem_bonds_valence = 1;
        }
        primary.fill(0);
        secondary.fill(0);
        bits = 0;
        assert_eq!(
            is_deriv_chain2(
                &amide,
                0,
                DERIV_BRIDGE_NH as i32,
                0,
                0,
                0,
                Some(&mut primary),
                64,
                Some(&mut secondary),
                64,
                Some(&mut bits),
            ),
            Ok(1)
        );
        assert_eq!(c_text(&primary), "RNH-C(O)CF3");
        assert_eq!(c_text(&secondary), "TFA");
        assert_eq!(bits, DERIV_BIT_TFA as i32);

        for (number, expected, name, bit) in [
            (2, "MOX", "MOX", DERIV_BIT_MOX),
            (3, "EtOX", "EtOX", DERIV_BIT_ETOX),
            (5, "OX-TMS", "TMS", DERIV_BIT_TMS),
            (7, "OX-???", "???", DERIV_BIT_UNKNOWN),
        ] {
            let oxime = connected_pair(EL_NUMBER_C, 1, 2);
            primary.fill(0);
            secondary.fill(0);
            bits = 0;
            assert_eq!(
                is_deriv_chain2(
                    &oxime,
                    0,
                    DERIV_X_OXIME as i32,
                    number,
                    0,
                    0,
                    Some(&mut primary),
                    64,
                    Some(&mut secondary),
                    64,
                    Some(&mut bits),
                ),
                Ok(1)
            );
            assert_eq!(c_text(&primary), expected);
            assert_eq!(c_text(&secondary), name);
            assert_eq!(bits, bit as i32);
        }

        for (number, suffix, name, bit) in [
            (3, "Me", "Acetate", DERIV_BIT_ACETATE),
            (6, "CF3", "TFA", DERIV_BIT_TFA),
            (8, "Phe", "Benzoate", DERIV_BIT_BENZOATE),
            (9, "C2F5", "PFP", DERIV_BIT_PFP),
            (12, "C3F7", "HFB", DERIV_BIT_HFB),
        ] {
            let ester = connected_pair(EL_NUMBER_O, 1, 2);
            primary.fill(0);
            secondary.fill(0);
            bits = 0;
            assert_eq!(
                is_deriv_chain2(
                    &ester,
                    0,
                    DERIV_RO_COX as i32,
                    number,
                    0,
                    0,
                    Some(&mut primary),
                    64,
                    Some(&mut secondary),
                    64,
                    Some(&mut bits),
                ),
                Ok(1)
            );
            assert_eq!(c_text(&primary), format!("RO-C(O){suffix}"));
            assert_eq!(c_text(&secondary), name);
            assert_eq!(bits, bit as i32);
        }

        for (type_, number, expected, name, bit) in [
            (DERIV_RING_DMOX_DEOX_N, 4, "DMOX", "DMOX", DERIV_BIT_DMOX),
            (DERIV_RING_DMOX_DEOX_N, 6, "DEOX", "DEOX", DERIV_BIT_DEOX),
            (
                DERIV_RING2_PRRLDD_OUTSIDE_PRECUR,
                4,
                "Pyrrolidide",
                "Pyrrolidide",
                DERIV_BIT_PYRROLIDIDE,
            ),
            (
                DERIV_RING2_PPRDN_OUTSIDE_PRECUR,
                5,
                "Piperidine",
                "Piperidine",
                DERIV_BIT_PIPERIDINE,
            ),
        ] {
            let ring = connected_pair(EL_NUMBER_N, 1, 1);
            primary.fill(0);
            secondary.fill(0);
            bits = 0;
            assert_eq!(
                is_deriv_chain2(
                    &ring,
                    0,
                    type_ as i32,
                    number,
                    0,
                    0,
                    Some(&mut primary),
                    64,
                    Some(&mut secondary),
                    64,
                    Some(&mut bits),
                ),
                Ok(1)
            );
            assert_eq!(c_text(&primary), expected);
            assert_eq!(c_text(&secondary), name);
            assert_eq!(bits, bit as i32);
        }

        let mut dansyl = connected_pair(EL_NUMBER_O, 1, 1);
        dansyl[0].num_H = 1;
        primary.fill(0);
        secondary.fill(0);
        bits = 0;
        assert_eq!(
            is_deriv_chain2(
                &dansyl,
                0,
                DERIV_DANSYL as i32,
                0,
                0,
                0,
                Some(&mut primary),
                64,
                Some(&mut secondary),
                64,
                Some(&mut bits),
            ),
            Ok(1)
        );
        assert_eq!(c_text(&primary), "ROH-Dansyl");
        assert_eq!(c_text(&secondary), "Dansyl");
        assert_eq!(bits, DERIV_BIT_DANSYL as i32);
    }

    #[test]
    fn source_port__ichinorm__is_deriv_chain__line_4980() {
        let mut atoms = vec![atom(&[1]), atom(&[0])];
        atoms[0].el_number = EL_NUMBER_O;
        atoms[0].chem_bonds_valence = 1;
        atoms[0].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
        atoms[1].chem_bonds_valence = 1;

        let mut derivative = DERIV_AT::default();
        assert_eq!(
            is_deriv_chain(
                &atoms,
                0,
                atoms.len() as i32,
                &derivative,
                0,
                None,
                0,
                None,
                0,
                None,
            ),
            Ok(0)
        );

        derivative.typ[0] = DERIV_DANSYL as i16;
        derivative.num[0] = 7;
        derivative.ord[0] = 0;
        let mut primary = [0_i8; 32];
        let mut secondary = [0_i8; 32];
        let mut bits = 0;
        assert_eq!(
            is_deriv_chain(
                &atoms,
                0,
                atoms.len() as i32,
                &derivative,
                0,
                Some(&mut primary),
                32,
                Some(&mut secondary),
                32,
                Some(&mut bits),
            ),
            Ok(1)
        );
        assert_eq!(
            &primary[..10],
            &[
                b'R' as i8, b'O' as i8, b'-' as i8, b'D' as i8, b'a' as i8, b'n' as i8, b's' as i8,
                b'y' as i8, b'l' as i8, 0,
            ]
        );
        assert_eq!(
            &secondary[..7],
            &[
                b'D' as i8, b'a' as i8, b'n' as i8, b's' as i8, b'y' as i8, b'l' as i8, 0,
            ]
        );
        assert_eq!(bits, DERIV_BIT_DANSYL as i32);

        assert_eq!(
            is_deriv_chain(
                &atoms,
                0,
                atoms.len() as i32,
                &derivative,
                -1,
                None,
                0,
                None,
                0,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__is_deriv_chain_or_ring__line_4996() {
        let mut atoms = vec![inp_ATOM::default()];
        atoms[0].charge = 1;

        let empty = DERIV_AT::default();
        let mut index = 0;
        assert_eq!(
            is_deriv_chain_or_ring(&atoms, 0, 1, &empty, &mut index),
            Ok(-1)
        );
        assert_eq!(index, 0);

        let mut chain = DERIV_AT::default();
        chain.typ[0] = DERIV_DANSYL as i16;
        assert_eq!(
            is_deriv_chain_or_ring(&atoms, 0, 1, &chain, &mut index),
            Ok(0)
        );

        let mut malformed = DERIV_AT::default();
        malformed.typ[0] = DERIV_RING_O_OUTSIDE_PRECURSOR as i16;
        assert_eq!(
            is_deriv_chain_or_ring(&atoms, 0, 1, &malformed, &mut index),
            Ok(-1)
        );
        assert_eq!(index, 0);

        let paired = DERIV_AT {
            typ: [DERIV_RING_O_OUTSIDE_PRECURSOR as i16; 4],
            ..DERIV_AT::default()
        };
        index = 1;
        assert_eq!(
            is_deriv_chain_or_ring(&atoms, 0, 1, &paired, &mut index),
            Ok(0)
        );
        assert_eq!(index, 0);
        index = 3;
        assert_eq!(
            is_deriv_chain_or_ring(&atoms, 0, 1, &paired, &mut index),
            Ok(0)
        );
        assert_eq!(index, 2);

        index = -1;
        assert_eq!(
            is_deriv_chain_or_ring(&atoms, 0, 1, &paired, &mut index),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__remove_deriv__line_5046() {
        let mut chain = DERIV_AT {
            typ: [1, 2, 3, 0],
            num: [11, 12, 13, 0],
            ord: [21, 22, 23, 0],
            other_atom: 99,
        };
        assert_eq!(remove_deriv(&mut chain, 1), Ok(0));
        assert_eq!(chain.typ, [1, 3, 0, 0]);
        assert_eq!(chain.num, [11, 13, 0, 0]);
        assert_eq!(chain.ord, [21, 23, 0, 0]);
        assert_eq!(chain.other_atom, 99);

        let ring_type = DERIV_RING_O_OUTSIDE_PRECURSOR as i16;
        let mut ring = DERIV_AT {
            typ: [ring_type, ring_type, 3, 4],
            num: [10, 11, 12, 13],
            ord: [20, 21, 22, 23],
            other_atom: 7,
        };
        assert_eq!(remove_deriv(&mut ring, 1), Ok(0));
        assert_eq!(ring.typ, [3, 4, 0, 0]);
        assert_eq!(ring.num, [12, 13, 0, 0]);
        assert_eq!(ring.ord, [22, 23, 0, 0]);
        assert_eq!(ring.other_atom, 7);

        let mut second_pair = DERIV_AT {
            typ: [3, ring_type, ring_type, 4],
            num: [10, 11, 12, 13],
            ord: [20, 21, 22, 23],
            other_atom: 0,
        };
        assert_eq!(remove_deriv(&mut second_pair, 2), Ok(0));
        assert_eq!(second_pair.typ, [3, 4, 0, 0]);
        assert_eq!(second_pair.num, [10, 13, 0, 0]);
        assert_eq!(second_pair.ord, [20, 23, 0, 0]);

        let mut malformed = DERIV_AT {
            typ: [ring_type, 3, 0, 0],
            num: [1, 2, 0, 0],
            ord: [3, 4, 0, 0],
            other_atom: 5,
        };
        let original = malformed.clone();
        assert_eq!(remove_deriv(&mut malformed, 0), Ok(-1));
        assert_eq!(malformed, original);
        assert_eq!(
            remove_deriv(&mut malformed, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__remove_deriv_mark__line_5108() {
        let mut chain = DERIV_AT {
            typ: [3, 4, 0, 0],
            num: [1, 2, 0, 0],
            ord: [5, 6, 0, 0],
            other_atom: 9,
        };
        assert_eq!(remove_deriv_mark(&mut chain, 1), Ok(0));
        assert_eq!(chain.typ, [3, 4 | DERIV_DUPLIC as i16, 0, 0]);
        assert_eq!(chain.num, [1, 2, 0, 0]);
        assert_eq!(chain.ord, [5, 6, 0, 0]);
        assert_eq!(chain.other_atom, 9);

        let ring_type = DERIV_RING_O_OUTSIDE_PRECURSOR as i16;
        let mut ring = DERIV_AT {
            typ: [3, ring_type, ring_type, 0],
            num: [1, 2, 3, 0],
            ord: [4, 5, 6, 0],
            other_atom: 7,
        };
        assert_eq!(remove_deriv_mark(&mut ring, 2), Ok(0));
        assert_eq!(
            ring.typ,
            [
                3,
                ring_type | DERIV_DUPLIC as i16,
                ring_type | DERIV_DUPLIC as i16,
                0,
            ]
        );
        assert_eq!(ring.num, [1, 2, 3, 0]);
        assert_eq!(ring.ord, [4, 5, 6, 0]);
        assert_eq!(ring.other_atom, 7);

        let mut malformed = DERIV_AT {
            typ: [ring_type, 3, 0, 0],
            ..DERIV_AT::default()
        };
        let original = malformed.clone();
        assert_eq!(remove_deriv_mark(&mut malformed, 0), Ok(-1));
        assert_eq!(malformed, original);
        assert_eq!(
            remove_deriv_mark(&mut malformed, 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__underiv_compare__line_5197() {
        assert_eq!(underiv_compare(&[b'a' as i8, 0], &[b'a' as i8, 0]), Ok(0));
        assert_eq!(underiv_compare(&[b'a' as i8, 0], &[b'b' as i8, 0]), Ok(-1));
        assert_eq!(underiv_compare(&[b'b' as i8, 0], &[b'a' as i8, 0]), Ok(1));
        assert_eq!(
            underiv_compare(&[b'a' as i8, 0, b'z' as i8], &[b'a' as i8, 0]),
            Ok(0)
        );
        assert_eq!(underiv_compare(&[-1, 0], &[1, 0]), Ok(254));
        assert_eq!(
            underiv_compare(&[b'a' as i8], &[b'a' as i8, 0]),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__underiv_list_add_two_cuts__line_5204() {
        fn c_bytes(bytes: &[u8]) -> Vec<i8> {
            bytes.iter().map(|&byte| byte as i8).collect()
        }

        fn output(input: &[u8], delimiter: u8, capacity: usize) -> (Vec<i8>, Vec<i8>) {
            let mut list = vec![b'?' as i8; capacity];
            list[0] = 0;
            let mut underiv = c_bytes(input);
            assert_eq!(
                underiv_list_add_two_cuts(
                    Some(&mut list),
                    capacity as i32,
                    &mut underiv,
                    delimiter as i8,
                ),
                Ok(0)
            );
            (list, underiv)
        }

        let (list, underiv) = output(b"Ac-10 Ac-2\0tail", b' ', 32);
        assert_eq!(&list[..8], &c_bytes(b"Ac-10+2\0"));
        assert_eq!(
            &underiv[..11],
            &[
                b'A' as i8, b'c' as i8, 0, b'1' as i8, b'0' as i8, 0, b'A' as i8, b'c' as i8, 0,
                b'2' as i8, 0,
            ]
        );
        assert_eq!(&underiv[11..], &c_bytes(b"tail"));

        let (list, _) = output(b"R-x R-x\0", b' ', 24);
        assert_eq!(&list[..5], &c_bytes(b"R-2x\0"));
        let (list, _) = output(b"P-z P-a\0", b' ', 24);
        assert_eq!(&list[..6], &c_bytes(b"P-a+z\0"));

        let (list, underiv) = output(b"A-x BB-y\0", b' ', 24);
        assert_eq!(&list[..9], &c_bytes(b"A-x BB-y\0"));
        assert_eq!(&underiv[..9], &c_bytes(b"A-x\0BB-y\0"));
        let (list, _) = output(b"A-x B-y\0", b' ', 24);
        assert_eq!(&list[..8], &c_bytes(b"A-x B-y\0"));
        let (list, _) = output(b"A-x B\0", b' ', 24);
        assert_eq!(&list[..6], &c_bytes(b"A-x B\0"));
        let (list, _) = output(b"solo\0tail", b' ', 24);
        assert_eq!(&list[..5], &c_bytes(b"solo\0"));

        let (list, underiv) = output(b"::A-x::A-y:third\0", b':', 24);
        assert_eq!(&list[..6], &c_bytes(b"A-x+y\0"));
        assert_eq!(underiv[0], b':' as i8);
        assert_eq!(underiv[1], b':' as i8);
        assert_eq!(underiv[5], 0);
        assert_eq!(underiv[6], b':' as i8);
        assert_eq!(underiv[10], 0);
        assert_eq!(&underiv[11..], &c_bytes(b"third\0"));

        let (list, underiv) = output(b"A-x A-y Z-z\0", b' ', 24);
        assert_eq!(&list[..6], &c_bytes(b"A-x+y\0"));
        assert_eq!(&underiv[7..], &c_bytes(b"\0Z-z\0"));

        let (list, underiv) = output(b"whole value\0", 0, 24);
        assert_eq!(&list[..12], &c_bytes(b"whole value\0"));
        assert_eq!(underiv, c_bytes(b"whole value\0"));
        let (list, _) = output(b"   \0", b' ', 8);
        assert_eq!(list[0], 0);

        let mut partial = c_bytes(b"\0????");
        let mut cuts = c_bytes(b"P-a P-b\0");
        assert_eq!(
            underiv_list_add_two_cuts(Some(&mut partial), 5, &mut cuts, b' ' as i8),
            Ok(0)
        );
        assert_eq!(partial, c_bytes(b"P-a+\0"));

        let mut null_target_input = c_bytes(b"Q-b Q-a\0");
        assert_eq!(
            underiv_list_add_two_cuts(None, 0, &mut null_target_input, b' ' as i8),
            Ok(0)
        );
        assert_eq!(
            null_target_input,
            &[b'Q' as i8, 0, b'b' as i8, 0, b'Q' as i8, 0, b'a' as i8, 0,]
        );

        let mut list = c_bytes(b"\0????????");
        let list_capacity = list.len() as i32;
        let mut unterminated = c_bytes(b"first second");
        assert_eq!(
            underiv_list_add_two_cuts(
                Some(&mut list),
                list_capacity,
                &mut unterminated,
                b' ' as i8,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(unterminated[5], 0);
        let mut empty = Vec::<i8>::new();
        assert_eq!(
            underiv_list_add_two_cuts(None, 0, &mut empty, b' ' as i8),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__eliminate_deriv_not_in_list__line_5351() {
        fn report_buffer() -> [i8; 64] {
            [0_i8; 64]
        }

        fn dansyl_fixture() -> Vec<inp_ATOM> {
            let mut atoms = vec![atom(&[1]), atom(&[0])];
            atoms[0].el_number = EL_NUMBER_O;
            atoms[0].chem_bonds_valence = 1;
            atoms[0].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
            atoms[1].chem_bonds_valence = 1;
            atoms
        }

        let atoms = dansyl_fixture();
        let mut empty = vec![DERIV_AT::default()];
        assert_eq!(
            eliminate_deriv_not_in_list(&atoms, &mut empty, -1, None, 0, None, 0, None,),
            Ok(0)
        );
        assert_eq!(
            eliminate_deriv_not_in_list(&atoms, &mut empty, 1, None, 0, None, 0, None,),
            Ok(0)
        );

        let mut too_many = vec![DERIV_AT::default()];
        too_many[0].typ[..3].fill(DERIV_DANSYL as i16);
        assert_eq!(
            eliminate_deriv_not_in_list(&atoms, &mut too_many, 1, None, 0, None, 0, None,),
            Ok(-1)
        );
        assert_eq!(too_many[0].typ[0], DERIV_DANSYL as i16);

        let mut mismatched = vec![DERIV_AT::default()];
        mismatched[0].typ[0] = DERIV_DANSYL as i16;
        mismatched[0].typ[1] = DERIV_BRIDGE_O as i16;
        assert_eq!(
            eliminate_deriv_not_in_list(&atoms, &mut mismatched, 1, None, 0, None, 0, None,),
            Ok(0)
        );
        assert_eq!(mismatched[0].typ[..2], [0, 0]);

        let mut single = vec![DERIV_AT::default()];
        single[0].typ[0] = DERIV_DANSYL as i16;
        single[0].num[0] = 7;
        let mut primary = report_buffer();
        let mut secondary = report_buffer();
        let mut bits = 0;
        assert_eq!(
            eliminate_deriv_not_in_list(
                &atoms,
                &mut single,
                1,
                Some(&mut primary),
                64,
                Some(&mut secondary),
                64,
                Some(&mut bits),
            ),
            Ok(1)
        );
        assert_eq!(&primary[..10], b"RO-Dansyl\0".map(|byte| byte as i8));
        assert_eq!(&secondary[..7], b"Dansyl\0".map(|byte| byte as i8));
        assert_eq!(bits, DERIV_BIT_DANSYL as i32);

        let mut oxime_atoms = vec![atom(&[1]), atom(&[0])];
        oxime_atoms[0].el_number = EL_NUMBER_C;
        oxime_atoms[0].chem_bonds_valence = 2;
        oxime_atoms[1].chem_bonds_valence = 1;
        let mut double = vec![DERIV_AT::default()];
        double[0].typ[..2].fill(DERIV_X_OXIME as i16);
        double[0].num[..2].copy_from_slice(&[2, 3]);
        primary.fill(0);
        secondary.fill(0);
        bits = 0;
        assert_eq!(
            eliminate_deriv_not_in_list(
                &oxime_atoms,
                &mut double,
                1,
                Some(&mut primary),
                64,
                Some(&mut secondary),
                64,
                Some(&mut bits),
            ),
            Ok(2)
        );
        assert_eq!(&primary[..9], b"MOX EtOX\0".map(|byte| byte as i8));
        assert_eq!(&secondary[..9], b"MOX EtOX\0".map(|byte| byte as i8));
        assert_eq!(bits, (DERIV_BIT_MOX | DERIV_BIT_ETOX) as i32);

        let mut no_room = vec![DERIV_AT::default()];
        no_room[0].typ[0] = DERIV_DANSYL as i16;
        let mut short = [0_i8; 4];
        secondary.fill(0);
        bits = 0;
        assert_eq!(
            eliminate_deriv_not_in_list(
                &atoms,
                &mut no_room,
                1,
                Some(&mut short),
                4,
                Some(&mut secondary),
                64,
                Some(&mut bits),
            ),
            Ok(1)
        );
        assert_eq!(short, [0; 4]);
        assert_eq!(bits, DERIV_BIT_DANSYL as i32);

        let mut null_bits = vec![DERIV_AT::default()];
        null_bits[0].typ[0] = DERIV_DANSYL as i16;
        assert_eq!(
            eliminate_deriv_not_in_list(&atoms, &mut null_bits, 1, None, 0, None, 0, None,),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut ring = vec![DERIV_AT::default()];
        ring[0].typ[..2].fill(DERIV_RING_O_OUTSIDE_PRECURSOR as i16);
        let mut charged_atoms = atoms.clone();
        charged_atoms[0].charge = 1;
        assert_eq!(
            eliminate_deriv_not_in_list(&charged_atoms, &mut ring, 1, None, 0, None, 0, None,),
            Ok(0)
        );
        assert_eq!(ring[0].typ[..2], [0, 0]);

        let mut malformed_atoms = vec![atom(&[1]), atom(&[])];
        malformed_atoms[0].chem_bonds_valence = 1;
        for number in [MAX_AT_DERIV as i8, -1_i8] {
            let mut boundary = vec![DERIV_AT::default()];
            boundary[0].typ[0] = DERIV_BRIDGE_O as i16;
            boundary[0].num[0] = number;
            assert_eq!(
                eliminate_deriv_not_in_list(
                    &malformed_atoms,
                    &mut boundary,
                    1,
                    None,
                    0,
                    None,
                    0,
                    None,
                ),
                Ok(-1)
            );
        }
        let mut above = vec![DERIV_AT::default()];
        above[0].typ[0] = DERIV_BRIDGE_O as i16;
        above[0].num[0] = MAX_AT_DERIV as i8 + 1;
        assert_eq!(
            eliminate_deriv_not_in_list(&malformed_atoms, &mut above, 1, None, 0, None, 0, None,),
            Ok(0)
        );
        assert_eq!(above[0].typ[0], 0);

        let mut unexpandable = vec![DERIV_AT::default()];
        unexpandable[0].typ[0] = DERIV_NOT as i16;
        unexpandable[0].num[0] = i8::MAX;
        assert_eq!(
            eliminate_deriv_not_in_list(
                &malformed_atoms,
                &mut unexpandable,
                1,
                None,
                0,
                None,
                0,
                None,
            ),
            Ok(-1)
        );

        let mut too_short = Vec::<DERIV_AT>::new();
        assert_eq!(
            eliminate_deriv_not_in_list(&atoms, &mut too_short, 1, None, 0, None, 0, None,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__sort_merge_underiv__line_5281() {
        fn c_buffer(text: &[u8]) -> [i8; MAX_SDF_VALUE as usize] {
            let mut buffer = [0_i8; MAX_SDF_VALUE as usize];
            for (target, &source) in buffer.iter_mut().zip(text) {
                *target = source as i8;
            }
            buffer
        }

        let mut output = c_buffer(b"");
        let mut list = c_buffer(b"b a a");
        assert_eq!(
            sort_merge_underiv(
                &mut output,
                1,
                &mut list,
                b',' as i8,
                &[b'[' as i8, 0],
                &[b']' as i8, 0],
            ),
            Ok(3)
        );
        assert_eq!(
            &output[..8],
            &[
                b'[' as i8, b'2' as i8, b'a' as i8, b',' as i8, b'b' as i8, b']' as i8, 0, 0,
            ]
        );
        assert_eq!(&list[..6], &[b'b' as i8, 0, b'a' as i8, 0, b'a' as i8, 0]);

        let mut empty_list = c_buffer(b"");
        output.fill(0);
        assert_eq!(
            sort_merge_underiv(
                &mut output,
                1,
                &mut empty_list,
                b',' as i8,
                &[b'(' as i8, 0],
                &[b')' as i8, 0],
            ),
            Ok(0)
        );
        assert_eq!(&output[..2], &[b'(' as i8, 0]);
        output.fill(0);
        assert_eq!(
            sort_merge_underiv(
                &mut output,
                0,
                &mut empty_list,
                b',' as i8,
                &[b'(' as i8, 0],
                &[b')' as i8, 0],
            ),
            Ok(0)
        );
        assert_eq!(&output[..3], &[b'(' as i8, b')' as i8, 0]);

        output = c_buffer(&vec![b'x'; 246]);
        let mut one = c_buffer(b"a");
        assert_eq!(
            sort_merge_underiv(&mut output, 1, &mut one, b',' as i8, &[b'p' as i8, 0], &[0],),
            Ok(-2)
        );
        assert_eq!(output[246], b'p' as i8);
        assert_eq!(output[247], b'!' as i8);
        assert_eq!(output[248], 0);

        let mut too_long = c_buffer(&vec![b'q'; 254]);
        too_long[254] = 0;
        output.fill(0);
        let original = too_long;
        assert_eq!(
            sort_merge_underiv(&mut output, 0, &mut too_long, b',' as i8, &[0], &[0]),
            Ok(0)
        );
        assert_eq!(too_long, original);
        assert_eq!(output[0], 0);
    }

    #[test]
    fn source_port__ichinorm__is_possibly_deriv_neigh__line_2610() {
        let mut ester = vec![inp_ATOM::default(); 4];
        ester[0] = atom(&[1, 2]);
        ester[0].el_number = EL_NUMBER_O;
        ester[0].chem_bonds_valence = 2;
        ester[1] = atom(&[0]);
        ester[1].el_number = EL_NUMBER_C;
        ester[1].chem_bonds_valence = 1;
        ester[2] = atom(&[0, 1, 3]);
        ester[2].el_number = EL_NUMBER_C;
        ester[2].chem_bonds_valence = 4;
        ester[3] = atom(&[2]);
        ester[3].el_number = EL_NUMBER_O;
        ester[3].chem_bonds_valence = 2;
        assert_eq!(
            is_possibly_deriv_neigh(&ester, 0, 1, DERIV_RO_COX as i32, i8::MIN),
            Ok(1)
        );
        assert_eq!(is_possibly_deriv_neigh(&ester, 0, 1, 12345, i8::MAX), Ok(0));

        let mut silyl = vec![inp_ATOM::default(); 6];
        silyl[0] = atom(&[1, 2]);
        silyl[0].el_number = EL_NUMBER_O;
        silyl[0].chem_bonds_valence = 2;
        silyl[1] = atom(&[0]);
        silyl[1].el_number = EL_NUMBER_C;
        silyl[1].chem_bonds_valence = 1;
        silyl[2] = atom(&[0, 3, 4, 5]);
        silyl[2].el_number = EL_NUMBER_SI;
        silyl[2].chem_bonds_valence = 4;
        for methyl in &mut silyl[3..] {
            methyl.el_number = EL_NUMBER_C;
            methyl.valence = 1;
            methyl.chem_bonds_valence = 1;
            methyl.num_H = 3;
        }
        assert_eq!(
            is_possibly_deriv_neigh(&silyl, 0, 1, DERIV_BRIDGE_O as i32, 0),
            Ok(1)
        );

        let mut blocked = ester.clone();
        blocked[1] = atom(&[0, 3, 3, 3]);
        blocked[1].el_number = EL_NUMBER_SI;
        blocked[1].chem_bonds_valence = 4;
        assert_eq!(
            is_possibly_deriv_neigh(&blocked, 0, 1, DERIV_BRIDGE_O as i32, 0),
            Ok(0)
        );

        let mut invalid_order = ester;
        invalid_order[0].neighbor[19] = u16::MAX;
        assert_eq!(
            is_possibly_deriv_neigh(&invalid_order, 0, 20, DERIV_RO_COX as i32, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__is_silyl2__line_3738() {
        fn methyl() -> inp_ATOM {
            inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 1,
                chem_bonds_valence: 1,
                num_H: 3,
                ..inp_ATOM::default()
            }
        }

        let mut tms = vec![inp_ATOM::default(); 5];
        tms[0] = atom(&[1, 2, 3, 4]);
        tms[0].el_number = EL_NUMBER_SI;
        tms[0].chem_bonds_valence = 4;
        tms[2] = methyl();
        tms[3] = methyl();
        tms[4] = methyl();
        assert_eq!(is_silyl2(&tms, 0, 1), Ok(1));

        for mutation in 0..5 {
            let mut rejected = tms.clone();
            match mutation {
                0 => rejected[0].el_number = EL_NUMBER_C,
                1 => rejected[0].valence = 3,
                2 => rejected[0].chem_bonds_valence = 3,
                3 => rejected[0].charge = 1,
                4 => rejected[0].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_silyl2(&rejected, 0, 1), Ok(0));
        }

        let mut tbdms = vec![inp_ATOM::default(); 8];
        tbdms[0] = atom(&[1, 2, 3, 4]);
        tbdms[0].el_number = EL_NUMBER_SI;
        tbdms[0].chem_bonds_valence = 4;
        tbdms[2] = methyl();
        tbdms[3] = methyl();
        tbdms[4] = atom(&[0, 5, 6, 7]);
        tbdms[4].el_number = EL_NUMBER_C;
        tbdms[4].chem_bonds_valence = 4;
        tbdms[5] = methyl();
        tbdms[6] = methyl();
        tbdms[7] = methyl();
        assert_eq!(is_silyl2(&tbdms, 0, 1), Ok(2));

        let mut silicon_chain = tbdms.clone();
        silicon_chain[4].el_number = EL_NUMBER_SI;
        assert_eq!(is_silyl2(&silicon_chain, 0, 1), Ok(0));
        let mut bad_methyl = tbdms.clone();
        bad_methyl[6].num_H = 2;
        assert_eq!(is_silyl2(&bad_methyl, 0, 1), Ok(0));
        let mut duplicate_carbon = tbdms;
        duplicate_carbon[2].valence = 4;
        duplicate_carbon[2].chem_bonds_valence = 4;
        duplicate_carbon[2].num_H = 0;
        assert_eq!(is_silyl2(&duplicate_carbon, 0, 1), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_cf3_or_linc3f7a__line_4000() {
        let mut values = vec![inp_ATOM::default(); 5];
        values[0] = atom(&[1, 2, 3, 4]);
        values[0].el_number = EL_NUMBER_C;
        values[0].chem_bonds_valence = 4;
        values[1] = atom(&[0]);
        values[1].el_number = EL_NUMBER_C;
        values[1].chem_bonds_valence = 1;
        for fluorine in 2..=4 {
            values[fluorine] = atom(&[0]);
            values[fluorine].el_number = EL_NUMBER_F;
            values[fluorine].chem_bonds_valence = 1;
        }
        assert_eq!(is_CF3_or_linC3F7a(&values, 0, 1), Ok(1));
        assert_eq!(is_CF3_or_linC3F7a(&values, 0, 9), Ok(0));
        assert_eq!(
            is_CF3_or_linC3F7a(&values, -1, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut duplicate_previous = values;
        duplicate_previous[0].neighbor[0] = 2;
        duplicate_previous[0].neighbor[1] = 2;
        assert_eq!(is_CF3_or_linC3F7a(&duplicate_previous, 0, 2), Ok(1));
    }

    #[test]
    fn source_port__ichinorm__get_traversed_deriv_type__line_2677() {
        let mut own_mark = vec![inp_ATOM {
            cFlags: 1,
            ..inp_ATOM::default()
        }];
        let mut derivatives = vec![DERIV_AT::default()];
        let mut output = DERIV_AT {
            typ: [7; 4],
            ord: [7; 4],
            num: [7; 4],
            other_atom: 7,
        };
        assert_eq!(
            get_traversed_deriv_type(&mut own_mark, &mut derivatives, 0, &mut output, 1,),
            Ok(0)
        );
        assert_eq!(output, DERIV_AT::default());

        let mut isolated = vec![inp_ATOM::default()];
        assert_eq!(
            get_traversed_deriv_type(&mut isolated, &mut derivatives, 0, &mut output, 1,),
            Ok(-1)
        );

        for element in [EL_NUMBER_O, 26] {
            let mut values = vec![inp_ATOM::default(); 2];
            values[0] = atom(&[1]);
            values[0].el_number = element;
            values[0].num_H = i8::from(element == EL_NUMBER_O);
            values[1].cFlags = 1;
            derivatives.resize(2, DERIV_AT::default());
            assert_eq!(
                get_traversed_deriv_type(&mut values, &mut derivatives, 0, &mut output, 1,),
                Ok(DERIV_NOT as i32)
            );
        }

        let mut unexpandable = vec![atom(&[1]), inp_ATOM::default()];
        unexpandable[1].cFlags = 1;
        derivatives = vec![DERIV_AT::default(); 2];
        derivatives[0].typ[0] = DERIV_X_OXIME as i16;
        assert_eq!(
            get_traversed_deriv_type(&mut unexpandable, &mut derivatives, 0, &mut output, 1,),
            Ok(0)
        );

        let mut oxime = vec![inp_ATOM::default(); 6];
        oxime[0] = atom(&[1, 2]);
        oxime[0].el_number = EL_NUMBER_N;
        oxime[0].chem_bonds_valence = 3;
        oxime[0].nNumAtInRingSystem = 1;
        oxime[0].bond_type[..2].copy_from_slice(&[BOND_DOUBLE as u8, BOND_SINGLE as u8]);
        oxime[1] = atom(&[0, 3, 4]);
        oxime[1].el_number = EL_NUMBER_C;
        oxime[1].chem_bonds_valence = 4;
        oxime[1].cFlags = 1;
        oxime[2] = atom(&[0, 5]);
        oxime[2].el_number = EL_NUMBER_O;
        oxime[2].chem_bonds_valence = 2;
        oxime[3].el_number = EL_NUMBER_C;
        oxime[4].el_number = EL_NUMBER_C;
        oxime[5] = atom(&[2]);
        oxime[5].el_number = EL_NUMBER_C;
        oxime[5].chem_bonds_valence = 1;
        oxime[5].num_H = 3;
        derivatives = vec![DERIV_AT::default(); oxime.len()];
        assert_eq!(
            get_traversed_deriv_type(&mut oxime, &mut derivatives, 0, &mut output, 1,),
            Ok(DERIV_X_OXIME as i32)
        );
        assert_eq!(output.typ[0], DERIV_X_OXIME as i16);
        assert_eq!(output.ord[0], 1);
        assert_eq!(output.num[0], 2);

        let mut ester = vec![inp_ATOM::default(); 7];
        ester[0] = atom(&[1, 2]);
        ester[0].el_number = EL_NUMBER_O;
        ester[0].chem_bonds_valence = 2;
        ester[0].nNumAtInRingSystem = 1;
        ester[1] = atom(&[0, 6]);
        ester[1].el_number = EL_NUMBER_C;
        ester[1].chem_bonds_valence = 2;
        ester[1].cFlags = 1;
        ester[2] = atom(&[0, 3, 4]);
        ester[2].el_number = EL_NUMBER_C;
        ester[2].chem_bonds_valence = 4;
        ester[2].bond_type[..3].copy_from_slice(&[
            BOND_SINGLE as u8,
            BOND_DOUBLE as u8,
            BOND_SINGLE as u8,
        ]);
        ester[3] = atom(&[2]);
        ester[3].el_number = EL_NUMBER_O;
        ester[3].chem_bonds_valence = 2;
        ester[4] = atom(&[2]);
        ester[4].el_number = EL_NUMBER_C;
        ester[4].chem_bonds_valence = 1;
        ester[4].num_H = 3;
        ester[6] = atom(&[1]);
        ester[6].el_number = EL_NUMBER_C;
        ester[6].chem_bonds_valence = 1;
        derivatives = vec![DERIV_AT::default(); ester.len()];
        assert_eq!(
            get_traversed_deriv_type(&mut ester, &mut derivatives, 0, &mut output, 1,),
            Ok(DERIV_RO_COX as i32)
        );
        assert_eq!(output.typ[0], DERIV_RO_COX as i16);
        assert_eq!(output.ord[0], 1);
        assert_eq!(output.num[0], 3);

        let mut amine = vec![inp_ATOM::default(); 3];
        amine[0] = atom(&[1, 2]);
        amine[0].el_number = EL_NUMBER_N;
        amine[0].chem_bonds_valence = 2;
        amine[0].num_H = 1;
        amine[0].bCutVertex = 1;
        amine[1].cFlags = 1;
        derivatives = vec![DERIV_AT::default(); amine.len()];
        assert_eq!(
            get_traversed_deriv_type(&mut amine, &mut derivatives, 0, &mut output, 1,),
            Ok(DERIV_BRIDGE_NH as i32)
        );
        assert_eq!(output.typ[0], DERIV_BRIDGE_NH as i16);
        assert_eq!(output.ord[0], 1);
    }

    #[test]
    fn source_port__ichinorm__mark_atoms_deriv__line_3332() {
        fn bridge_atoms(child: inp_ATOM) -> Vec<inp_ATOM> {
            let mut center = atom(&[1, 2]);
            center.el_number = EL_NUMBER_N;
            center.chem_bonds_valence = 2;
            center.num_H = 1;
            center.bCutVertex = 1;
            let marked = inp_ATOM {
                cFlags: 1,
                ..inp_ATOM::default()
            };
            vec![center, marked, child]
        }

        fn oxime_atoms() -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 6];
            values[0] = atom(&[1, 2]);
            values[0].el_number = EL_NUMBER_N;
            values[0].chem_bonds_valence = 3;
            values[0].nNumAtInRingSystem = 1;
            values[0].bond_type[..2].copy_from_slice(&[BOND_DOUBLE as u8, BOND_SINGLE as u8]);
            values[1] = atom(&[0, 3, 4]);
            values[1].el_number = EL_NUMBER_C;
            values[1].chem_bonds_valence = 4;
            values[1].cFlags = 1;
            values[2] = atom(&[0, 5]);
            values[2].el_number = EL_NUMBER_O;
            values[2].chem_bonds_valence = 2;
            values[3].el_number = EL_NUMBER_C;
            values[4].el_number = EL_NUMBER_C;
            values[5] = atom(&[2]);
            values[5].el_number = EL_NUMBER_C;
            values[5].chem_bonds_valence = 1;
            values[5].num_H = 3;
            values
        }

        fn dmox_atoms(counterpart_is_active: bool) -> Vec<inp_ATOM> {
            let mut values = vec![inp_ATOM::default(); 9];
            let ring = [
                (EL_NUMBER_N, 2, 3, 0, [4, 1, 0, 0]),
                (EL_NUMBER_C, 4, 4, 0, [0, 2, 5, 6]),
                (EL_NUMBER_C, 2, 2, 2, [1, 3, 0, 0]),
                (EL_NUMBER_O, 2, 2, 0, [2, 4, 0, 0]),
                (EL_NUMBER_C, 3, 4, 0, [3, 0, 7, 0]),
            ];
            for (index, (element, valence, bond_valence, hydrogens, neighbors)) in
                ring.into_iter().enumerate()
            {
                values[index].el_number = element;
                values[index].valence = valence;
                values[index].chem_bonds_valence = bond_valence;
                values[index].num_H = hydrogens;
                values[index].neighbor[..usize::try_from(valence).unwrap()]
                    .copy_from_slice(&neighbors[..usize::try_from(valence).unwrap()]);
                values[index].bond_type[1] = if index == 4 {
                    BOND_DOUBLE as u8
                } else {
                    BOND_SINGLE as u8
                };
                values[index].nRingSystem = 1;
            }
            for atom in &mut values[..5] {
                atom.bond_type[0] = BOND_SINGLE as u8;
            }
            values[0].bond_type[0] = BOND_DOUBLE as u8;
            values[0].nNumAtInRingSystem = 5;
            values[3].nNumAtInRingSystem = u16::from(counterpart_is_active) * 5;
            values[4].cFlags = 1;
            for index in [5, 6] {
                values[index] = inp_ATOM {
                    el_number: EL_NUMBER_C,
                    valence: 1,
                    chem_bonds_valence: 1,
                    num_H: 3,
                    neighbor: {
                        let mut neighbors = [0; 20];
                        neighbors[0] = 1;
                        neighbors
                    },
                    ..inp_ATOM::default()
                };
            }
            values[7].el_number = EL_NUMBER_C;
            values
        }

        let mut already_marked = vec![
            inp_ATOM {
                cFlags: 1,
                ..atom(&[1])
            },
            inp_ATOM::default(),
        ];
        let mut derivatives = vec![DERIV_AT::default(); already_marked.len()];
        let mut found = 9;
        assert_eq!(
            mark_atoms_deriv(&mut already_marked, &mut derivatives, 0, 17, 1, &mut found),
            Ok(17)
        );
        assert_eq!(already_marked[1].cFlags, 0);
        assert_eq!(found, 9);

        let mut isolated = vec![inp_ATOM::default()];
        derivatives = vec![DERIV_AT::default()];
        found = 0;
        assert_eq!(
            mark_atoms_deriv(&mut isolated, &mut derivatives, 0, i32::MAX, 1, &mut found,),
            Ok(i32::MIN)
        );
        assert_eq!(isolated[0].cFlags, 1);
        assert_eq!(found, 0);

        let mut prohibited = vec![
            inp_ATOM {
                el_number: EL_NUMBER_O,
                num_H: 1,
                ..atom(&[1])
            },
            inp_ATOM {
                cFlags: 1,
                ..inp_ATOM::default()
            },
        ];
        derivatives = vec![DERIV_AT::default(); prohibited.len()];
        found = i32::MAX;
        assert_eq!(
            mark_atoms_deriv(&mut prohibited, &mut derivatives, 0, 0, 1, &mut found),
            Ok(1)
        );
        assert_eq!(found, i32::MIN);

        let child = inp_ATOM {
            el_number: EL_NUMBER_C,
            ..atom(&[0])
        };
        let mut accepted = bridge_atoms(child);
        derivatives = vec![DERIV_AT::default(); accepted.len()];
        found = 0;
        assert_eq!(
            mark_atoms_deriv(&mut accepted, &mut derivatives, 0, 0, 1, &mut found),
            Ok(2)
        );
        assert_eq!(derivatives[0].typ[0], DERIV_BRIDGE_NH as i16);
        assert_eq!(derivatives[0].ord[0], 1);
        assert_eq!(derivatives[0].num[0], 1);
        assert_eq!(found, 1);

        let child = inp_ATOM {
            el_number: EL_NUMBER_O,
            num_H: 1,
            ..atom(&[0])
        };
        let mut rejected = bridge_atoms(child);
        derivatives = vec![DERIV_AT::default(); rejected.len()];
        found = 0;
        assert_eq!(
            mark_atoms_deriv(&mut rejected, &mut derivatives, 0, 0, 1, &mut found),
            Ok(2)
        );
        assert_eq!(derivatives[0], DERIV_AT::default());
        assert_eq!(found, 1);

        let mut over_limit = bridge_atoms(inp_ATOM::default());
        over_limit.resize(16, inp_ATOM::default());
        for index in 2..16 {
            let mut neighbors = Vec::new();
            neighbors.push(if index == 2 { 0 } else { index - 1 });
            if index + 1 < 16 {
                neighbors.push(index + 1);
            }
            over_limit[index] = atom(
                &neighbors
                    .into_iter()
                    .map(|value| value as u16)
                    .collect::<Vec<_>>(),
            );
            over_limit[index].el_number = EL_NUMBER_C;
            over_limit[index].chem_bonds_valence = over_limit[index].valence;
        }
        derivatives = vec![DERIV_AT::default(); over_limit.len()];
        found = 0;
        assert_eq!(
            mark_atoms_deriv(&mut over_limit, &mut derivatives, 0, 0, 1, &mut found),
            Ok(15)
        );
        assert_eq!(derivatives[0], DERIV_AT::default());
        assert_eq!(found, 0);

        let mut conflict = bridge_atoms(inp_ATOM {
            el_number: EL_NUMBER_C,
            ..atom(&[0])
        });
        derivatives = vec![DERIV_AT::default(); conflict.len()];
        derivatives[0].typ[0] = DERIV_BRIDGE_O as i16;
        derivatives[0].ord[0] = 1;
        derivatives[0].num[0] = 1;
        found = 0;
        assert_eq!(
            mark_atoms_deriv(&mut conflict, &mut derivatives, 0, 0, 1, &mut found),
            Ok(-1)
        );

        let mut oxime = oxime_atoms();
        derivatives = vec![DERIV_AT::default(); oxime.len()];
        found = 0;
        assert_eq!(
            mark_atoms_deriv(&mut oxime, &mut derivatives, 0, 0, 1, &mut found),
            Ok(5)
        );
        assert_eq!(derivatives[0].typ[0], DERIV_X_OXIME as i16);
        assert_eq!(derivatives[0].ord[0], 1);
        assert_eq!(derivatives[0].num[0], 2);
        assert_eq!(found, 1);

        let mut dmox = dmox_atoms(true);
        derivatives = vec![DERIV_AT::default(); dmox.len()];
        let mut nitrogen_output = DERIV_AT::default();
        let mut nitrogen_probe = dmox.clone();
        let mut nitrogen_derivatives = derivatives.clone();
        assert_eq!(
            get_traversed_deriv_type(
                &mut nitrogen_probe,
                &mut nitrogen_derivatives,
                0,
                &mut nitrogen_output,
                1,
            ),
            Ok(DERIV_RING_DMOX_DEOX_N as i32)
        );
        let mut oxygen_output = DERIV_AT::default();
        assert_eq!(
            get_traversed_deriv_type(
                &mut nitrogen_probe,
                &mut nitrogen_derivatives,
                3,
                &mut oxygen_output,
                1,
            ),
            Ok(DERIV_RING_DMOX_DEOX_O as i32)
        );
        found = 0;
        assert!(mark_atoms_deriv(&mut dmox, &mut derivatives, 0, 0, 1, &mut found).is_ok());
        assert_eq!(derivatives[0].typ[0], DERIV_RING_DMOX_DEOX_N as i16);
        assert_eq!(derivatives[3].typ[0], DERIV_RING_DMOX_DEOX_O as i16);
        assert_eq!(found, 1);

        let mut mismatch = dmox_atoms(false);
        derivatives = vec![DERIV_AT::default(); mismatch.len()];
        found = 0;
        assert!(mark_atoms_deriv(&mut mismatch, &mut derivatives, 0, 0, 1, &mut found).is_ok());
        assert_eq!(derivatives[0], DERIV_AT::default());
        assert_eq!(found, 0);
    }

    #[test]
    fn source_port__ichinorm__count_one_bond_atoms__line_3536() {
        let mut atoms = vec![atom(&[1]), atom(&[0])];
        atoms[1].cFlags = 1;
        let mut derivatives = vec![DERIV_AT::default(); atoms.len()];
        let mut found = 7;
        assert_eq!(
            count_one_bond_atoms(&mut atoms, &mut derivatives, 0, 0, 1, &mut found),
            Ok(0)
        );
        assert_eq!(atoms[0].cFlags, 0);
        assert_eq!(found, 7);

        atoms[1].cFlags = 0;
        found = 0;
        assert_eq!(
            count_one_bond_atoms(&mut atoms, &mut derivatives, 0, 0, 1, &mut found),
            Ok(2)
        );
        assert_eq!((atoms[0].cFlags, atoms[1].cFlags), (1, 1));
        assert_eq!(found, 0);

        let mut prohibited = vec![
            inp_ATOM {
                el_number: EL_NUMBER_O,
                num_H: 1,
                ..atom(&[1])
            },
            atom(&[0]),
        ];
        derivatives = vec![DERIV_AT::default(); prohibited.len()];
        found = 0;
        assert_eq!(
            count_one_bond_atoms(&mut prohibited, &mut derivatives, 0, 0, 1, &mut found),
            Ok(2)
        );
        assert_eq!(found, 1);

        assert_eq!(
            count_one_bond_atoms(&mut prohibited, &mut derivatives, 0, -1, 2, &mut found),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            count_one_bond_atoms(&mut prohibited, &mut derivatives, 0, 20, 2, &mut found),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__is_nbutyl__line_3831() {
        fn carbon(neighbors: &[u16], hydrogens: i8) -> inp_ATOM {
            let mut value = atom(neighbors);
            value.el_number = EL_NUMBER_C;
            value.chem_bonds_valence = value.valence;
            value.num_H = hydrogens;
            value
        }

        let chain = vec![
            inp_ATOM::default(),
            carbon(&[0, 2], 2),
            carbon(&[1, 3], 2),
            carbon(&[2, 4], 2),
            carbon(&[3], 3),
        ];
        assert_eq!(is_nButyl(&chain, 1, 0), Ok(1));
        assert_eq!(
            is_nButyl(&chain, -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            is_nButyl(&chain, 1, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut short = chain.clone();
        short[3] = carbon(&[2], 3);
        assert_eq!(is_nButyl(&short, 1, 0), Ok(0));

        let mut long = chain.clone();
        long.push(carbon(&[4], 3));
        long[4] = carbon(&[3, 5], 2);
        assert_eq!(is_nButyl(&long, 1, 0), Ok(0));

        let mut cycle = chain.clone();
        cycle[4] = carbon(&[3, 1], 2);
        assert_eq!(is_nButyl(&cycle, 1, 0), Ok(0));

        for mutation in 0..5 {
            let mut rejected = chain.clone();
            match mutation {
                0 => rejected[2].el_number = EL_NUMBER_N,
                1 => rejected[2].valence = 3,
                2 => rejected[2].chem_bonds_valence = 1,
                3 => rejected[2].charge = 1,
                4 => rejected[2].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_nButyl(&rejected, 1, 0), Ok(0));
        }
    }

    #[test]
    fn source_port__ichinorm__is_me_or_et__line_3862() {
        let methyl = inp_ATOM {
            el_number: EL_NUMBER_C,
            num_H: 4,
            ..inp_ATOM::default()
        };
        assert_eq!(is_Me_or_Et(std::slice::from_ref(&methyl), 0, -1), Ok(1));

        let ethyl = vec![
            inp_ATOM {
                el_number: EL_NUMBER_C,
                chem_bonds_valence: 1,
                num_H: 3,
                ..atom(&[1])
            },
            inp_ATOM {
                el_number: EL_NUMBER_C,
                chem_bonds_valence: 1,
                num_H: 3,
                ..atom(&[0])
            },
        ];
        assert_eq!(is_Me_or_Et(&ethyl, 0, 0), Ok(1));
        assert_eq!(is_Me_or_Et(&ethyl, 0, -1), Ok(2));
        assert_eq!(is_Me_or_Et(&ethyl, 0, 7), Ok(2));

        let mut branched = ethyl.clone();
        branched.push(ethyl[1].clone());
        branched[0] = inp_ATOM {
            el_number: EL_NUMBER_C,
            chem_bonds_valence: 2,
            num_H: 2,
            ..atom(&[1, 2])
        };
        assert_eq!(is_Me_or_Et(&branched, 0, -1), Ok(0));

        for mutation in 0..6 {
            let mut rejected = ethyl.clone();
            match mutation {
                0 => rejected[0].el_number = EL_NUMBER_N,
                1 => rejected[0].valence = 3,
                2 => rejected[0].chem_bonds_valence = 0,
                3 => rejected[1].el_number = EL_NUMBER_O,
                4 => rejected[1].charge = 1,
                5 => rejected[1].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_Me_or_Et(&rejected, 0, -1), Ok(0));
        }
        assert_eq!(
            is_Me_or_Et(&ethyl, -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichinorm__is_phenyl__line_4103() {
        let mut phenyl = vec![inp_ATOM::default(); 7];
        phenyl[0] = inp_ATOM {
            el_number: EL_NUMBER_C,
            valence: 3,
            chem_bonds_valence: 4,
            neighbor: {
                let mut neighbors = [0; 20];
                neighbors[..3].copy_from_slice(&[6, 1, 5]);
                neighbors
            },
            ..inp_ATOM::default()
        };
        for index in 1..=5 {
            let previous = if index == 1 { 0 } else { index - 1 };
            let next = if index == 5 { 0 } else { index + 1 };
            phenyl[index] = inp_ATOM {
                el_number: EL_NUMBER_C,
                valence: 2,
                chem_bonds_valence: 3,
                num_H: 1,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[..2].copy_from_slice(&[previous as u16, next as u16]);
                    neighbors
                },
                ..inp_ATOM::default()
            };
        }
        assert_eq!(is_phenyl(&phenyl, 0, 0), Ok(1));
        assert_eq!(is_phenyl(&phenyl, 0, 2), Ok(0));
        assert_eq!(
            is_phenyl(&phenyl, -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        for mutation in 0..6 {
            let mut rejected = phenyl.clone();
            match mutation {
                0 => rejected[0].el_number = EL_NUMBER_N,
                1 => rejected[0].valence = 2,
                2 => rejected[0].chem_bonds_valence = 3,
                3 => rejected[0].num_H = 2,
                4 => rejected[0].charge = 1,
                5 => rejected[3].radical = 1,
                _ => unreachable!(),
            }
            assert_eq!(is_phenyl(&rejected, 0, 0), Ok(0));
        }
        let mut open = phenyl;
        open[5].neighbor[1] = 6;
        assert_eq!(is_phenyl(&open, 0, 0), Ok(0));
    }

    #[test]
    fn source_port__ichinorm__is_deriv_ring_o_or_nh_outside_precursor__line_4136() {
        fn atom(element: u8, neighbors: &[u16]) -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: element,
                valence: neighbors.len() as i8,
                chem_bonds_valence: neighbors.len() as i8,
                ..inp_ATOM::default()
            };
            atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            atom.bond_type[..neighbors.len()].fill(BOND_SINGLE as u8);
            atom
        }

        let mut atoms = vec![inp_ATOM::default(); 8];
        atoms[0] = atom(EL_NUMBER_B, &[1, 4, 5]);
        atoms[1] = atom(EL_NUMBER_O, &[0, 2]);
        atoms[2] = atom(EL_NUMBER_C, &[1, 3]);
        atoms[3] = atom(EL_NUMBER_C, &[2, 4]);
        atoms[4] = atom(EL_NUMBER_O, &[0, 3]);
        atoms[5] = atom(EL_NUMBER_C, &[0]);
        atoms[1].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
        atoms[4].elname[..2].copy_from_slice(&[b'O' as i8, 0]);

        let derivative = DERIV_AT {
            typ: [
                DERIV_RING_O_OUTSIDE_PRECURSOR as i16,
                DERIV_RING_O_OUTSIDE_PRECURSOR as i16,
                0,
                0,
            ],
            ord: [b'0' as i8, b'1' as i8, 0, 0],
            ..DERIV_AT::default()
        };
        let mut underiv = [0_i8; 16];
        let mut underiv2 = [b'X' as i8, 0, 0, 0];
        let mut bit_underiv = 7;
        assert_eq!(
            is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                &atoms,
                0,
                atoms.len() as i32,
                &derivative,
                0,
                Some(&mut underiv),
                16,
                Some(&mut underiv2),
                4,
                Some(&mut bit_underiv),
            ),
            Ok(-1)
        );
        assert_eq!(
            &underiv[..5],
            &[b'O' as i8, b'O' as i8, b'-' as i8, b'B' as i8, 0]
        );
        assert_eq!(underiv2, [b'X' as i8, 0, 0, 0]);
        assert_eq!(bit_underiv, 7);

        let mut invalid = derivative.clone();
        invalid.typ[0] = 0;
        let mut untouched = [b'Q' as i8, 0, 0, 0];
        assert_eq!(
            is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                &atoms,
                0,
                atoms.len() as i32,
                &invalid,
                0,
                Some(&mut untouched),
                4,
                None,
                0,
                None,
            ),
            Ok(0)
        );
        assert_eq!(untouched, [b'Q' as i8, 0, 0, 0]);

        for mutation in 0..3 {
            let mut rejected = atoms.clone();
            match mutation {
                0 => rejected[0].charge = 1,
                1 => rejected[0].radical = 1,
                2 => rejected[0].chem_bonds_valence = 2,
                _ => unreachable!(),
            }
            assert_eq!(
                is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                    &rejected,
                    0,
                    rejected.len() as i32,
                    &derivative,
                    0,
                    None,
                    0,
                    None,
                    0,
                    None,
                ),
                Ok(0)
            );
        }

        let mut numeric_order = derivative.clone();
        numeric_order.ord[..2].copy_from_slice(&[0, 1]);
        assert_eq!(
            is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                &atoms,
                0,
                atoms.len() as i32,
                &numeric_order,
                0,
                None,
                0,
                None,
                0,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut open_ring = atoms.clone();
        open_ring[4].neighbor[1] = 6;
        assert_eq!(
            is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                &open_ring,
                0,
                open_ring.len() as i32,
                &derivative,
                0,
                None,
                0,
                None,
                0,
                None,
            ),
            Ok(0)
        );

        let mut hetero_rejected = atoms;
        hetero_rejected[2] = atom(EL_NUMBER_C, &[1, 3, 6, 7]);
        hetero_rejected[6] = atom(EL_NUMBER_N, &[2]);
        hetero_rejected[7] = atom(EL_NUMBER_O, &[2]);
        assert_eq!(
            is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(
                &hetero_rejected,
                0,
                hetero_rejected.len() as i32,
                &derivative,
                0,
                None,
                0,
                None,
                0,
                None,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichinorm__add_to_da__line_3245() {
        fn entry(type_: i16, order: i8, number: i8) -> DERIV_AT {
            let mut value = DERIV_AT::default();
            value.typ[0] = type_;
            value.ord[0] = order;
            value.num[0] = number;
            value
        }

        let mut derivative = entry(DERIV_BRIDGE_O as i16, 1, 2);
        let added = entry(DERIV_BRIDGE_NH as i16, 3, 4);
        assert_eq!(add_to_da(&mut derivative, &added), 0);
        assert_eq!(
            (
                &derivative.typ[..2],
                &derivative.ord[..2],
                &derivative.num[..2]
            ),
            (
                &[DERIV_BRIDGE_O as i16, DERIV_BRIDGE_NH as i16][..],
                &[1, 3][..],
                &[2, 4][..],
            )
        );

        let original = entry(DERIV_BRIDGE_O as i16, 1, 2);
        let mut same = original.clone();
        assert_eq!(add_to_da(&mut same, &original), 0);
        let mut different_type = original.clone();
        assert_eq!(
            add_to_da(&mut different_type, &entry(DERIV_BRIDGE_NH as i16, 1, 2),),
            -1
        );
        let mut different_number = original.clone();
        assert_eq!(
            add_to_da(&mut different_number, &entry(DERIV_BRIDGE_O as i16, 1, 3),),
            -2
        );
        let mut different_other = original.clone();
        different_other.other_atom = 4;
        let mut added_other = original.clone();
        added_other.other_atom = 5;
        assert_eq!(add_to_da(&mut different_other, &added_other), -3);

        let mut full = DERIV_AT {
            typ: [DERIV_BRIDGE_O as i16; 4],
            ord: [0, 1, 2, 3],
            num: [1; 4],
            other_atom: 0,
        };
        assert_eq!(
            add_to_da(&mut full, &entry(DERIV_BRIDGE_O as i16, 4, 1)),
            -4
        );

        let mut low_priority = entry(DERIV_BRIDGE_O as i16, 1, 2);
        let high_priority = entry(DERIV_X_OXIME as i16, 2, 5);
        assert_eq!(add_to_da(&mut low_priority, &high_priority), 0);
        assert_eq!(low_priority, high_priority);

        let mut high_kept = high_priority.clone();
        assert_eq!(add_to_da(&mut high_kept, &original), 0);
        assert_eq!(high_kept, high_priority);

        let mut source_condition = original.clone();
        let mut mixed = original.clone();
        mixed.typ[1] = DERIV_X_OXIME as i16;
        mixed.ord[1] = 2;
        mixed.num[1] = 5;
        assert_eq!(add_to_da(&mut source_condition, &mixed), 0);
        assert_eq!(source_condition.typ[0], DERIV_BRIDGE_O as i16);
        assert_eq!(source_condition.typ[1], DERIV_X_OXIME as i16);

        let mut one = original;
        let mut second_with_other = entry(DERIV_BRIDGE_NH as i16, 2, 1);
        second_with_other.other_atom = 9;
        assert_eq!(add_to_da(&mut one, &second_with_other), -3);
    }
}
