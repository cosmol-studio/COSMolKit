use crate::source::base::ichisort::insertions_sort_AT_RANK;
use crate::source::base::ichister::get_opposite_sb_atom_slice;
use crate::source::base::runichi2::bIsSameBond;
use crate::source::base::util::{
    get_el_type, get_el_valence, get_periodic_table_number, has_other_ion_in_sphere_2,
    has_other_ion_neigh, inchi_calloc, inchi_free, ion_el_group, is_in_the_ilist, is_in_the_list,
    n_no_metal_bonds_valence, n_no_metal_neigh_index, n_no_metal_num_bonds,
    n_no_metal_other_neigh_index, n_no_metal_other_neigh_index2, num_of_H,
};
use crate::source_types::{
    AB_MAX_WELL_DEFINED_PARITY, AB_MIN_WELL_DEFINED_PARITY, AB_PARITY_UNDF, AT_NUMB, AT_RANK,
    ATW_H, BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, COMP_ATOM_DATA, CT_OUT_OF_RAM,
    INCHI_MODE, INCHI_T_NUM_MOVABLE, INChI, INChI_Aux, INChI_IsotopicAtom, INChI_IsotopicTGroup,
    INChI_Stereo, MAX_ATOMS, MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS, MAXVAL, MOL_COORD,
    NUM_H_ISOTOPES, ORIG_ATOM_DATA, ORIG_INFO, RADICAL_DOUBLET, RADICAL_SINGLET, REQ_MODE_ISO,
    S_CHAR, SourceHeap, SourceHeapError, SourceMutPointer, TG_FLAG_CHECK_VALENCE_COORD_DONE,
    TG_FLAG_MOVE_CHARGE_COORD_DONE, U_CHAR, inp_ATOM, subgraf, subgraf_edge, subgraf_pathfinder,
};

#[allow(non_snake_case)]
pub(crate) fn add_DT_to_num_H(
    num_atoms: i32,
    atoms: &mut [inp_ATOM],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3705 add_DT_to_num_H
    // INCHI✔️✔️: int add_DT_to_num_H( int num_atoms, inp_ATOM *at )
    // INCHI✔️✔️: /*  assume num_1H, num_D and num_T are not included in num_H */
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j;
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (j = 0; j < NUM_H_ISOTOPES; j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             at[i].num_H += at[i].num_iso_H[j];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: add_DT_to_num_H
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: add_DT_to_num_H
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: add_DT_to_num_H

    let atom_count =
        usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let selected = atoms
        .get_mut(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for atom in selected {
        for isotope_hydrogens in atom.num_iso_H {
            atom.num_H = atom.num_H.wrapping_add(isotope_hydrogens);
        }
    }
    Ok(0)
}

pub(crate) fn cmp_iso_atw_diff_component_no(first: &inp_ATOM, second: &inp_ATOM) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:166 cmp_iso_atw_diff_component_no
    // INCHI✔️✔️: int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int ret = (int) ( (const inp_ATOM*) a1 )->iso_atw_diff - (int) ( (const inp_ATOM*) a2 )->iso_atw_diff;
    // INCHI✔️✔️:     if (!ret) /*  make the sort stable */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ret = (int) ( (const inp_ATOM*) a1 )->component - (int) ( (const inp_ATOM*) a2 )->component;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: cmp_iso_atw_diff_component_no
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: cmp_iso_atw_diff_component_no
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // INCHI✔️✔️: typedef unsigned short AT_NUMB;
    // INCHI✔️✔️:     S_CHAR iso_atw_diff;              /* =0 => natural isotopic abundances                        */
    // INCHI✔️✔️:     AT_NUMB component; /* number of the structure component > 0                    */
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: cmp_iso_atw_diff_component_no

    let isotope_order = i32::from(first.iso_atw_diff) - i32::from(second.iso_atw_diff);
    if isotope_order != 0 {
        isotope_order
    } else {
        i32::from(first.component) - i32::from(second.component)
    }
}

#[allow(non_snake_case)]
pub(crate) fn nFindOneOM(
    atoms: &[inp_ATOM],
    atom_number: i32,
    ord_om: &mut [i32],
    mut num_om: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:932 nFindOneOM
    // INCHI✔️✔️: int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, n_OM, best_value, cur_value, diff; /* djb-rwth: removing redundant variables */
    // INCHI✔️✔️:     int num_best;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (1 == num_OM)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ord_OM[0];
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (1 > num_OM)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* select neighbors with min. number of bonds */
    // INCHI✔️✔️:     num_best = 1;
    // INCHI✔️✔️:     n_OM = (int) at[at_no].neighbor[ord_OM[0]];
    // INCHI✔️✔️:     best_value = (int) at[n_OM].valence;
    // INCHI✔️✔️:     /* compare number of bonds; move indexes of the best neighbors to the first elements of ord_OM[] */
    // INCHI✔️✔️:     for (i = 1; i < num_OM; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         n_OM = at[at_no].neighbor[ord_OM[i]];
    // INCHI✔️✔️:         cur_value = (int) at[n_OM].valence;
    // INCHI✔️✔️:         diff = cur_value - best_value;
    // INCHI✔️✔️:         if (diff < 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:             best_value = cur_value;
    // INCHI✔️✔️:             ord_OM[0] = ord_OM[i];
    // INCHI✔️✔️:             num_best = 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (diff == 0)
    // INCHI✔️✔️:         {    /* was '=', pointed by WDI */
    // INCHI✔️✔️:             ord_OM[num_best++] = ord_OM[i];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     num_OM = num_best;
    // INCHI✔️✔️:     if (1 == num_OM)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ord_OM[0];
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* select neighbors with min. periodic numbers */
    // INCHI✔️✔️:     num_best = 1;
    // INCHI✔️✔️:     n_OM = (int) at[at_no].neighbor[ord_OM[0]];
    // INCHI✔️✔️:     best_value = (int) at[n_OM].el_number;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* compare periodic numbers; move indexes of the best neighbors to the first elements of ord_OM[] */
    // INCHI✔️✔️:     for (i = 1; i < num_OM; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         n_OM = at[at_no].neighbor[ord_OM[i]];
    // INCHI✔️✔️:         cur_value = (int) at[n_OM].el_number;
    // INCHI✔️✔️:         diff = cur_value - best_value;
    // INCHI✔️✔️:         if (diff < 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:             best_value = cur_value;
    // INCHI✔️✔️:             ord_OM[0] = ord_OM[i];
    // INCHI✔️✔️:             num_best = 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (diff == 0)
    // INCHI✔️✔️:         {    /* was '=', pointed by WDI */
    // INCHI✔️✔️:             ord_OM[num_best++] = ord_OM[i];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     num_OM = num_best;
    // INCHI✔️✔️:     if (1 == num_OM)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ord_OM[0];
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* if neighbors are not terminal atoms then reject */
    // INCHI✔️✔️:     if (1 < at[n_OM].valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* if neighbors are terminal atoms then the one without isotope or with lightest isotope */
    // INCHI✔️✔️:     num_best = 1;
    // INCHI✔️✔️:     n_OM = (int) at[at_no].neighbor[ord_OM[0]];
    // INCHI✔️✔️:     best_value = (int) at[n_OM].iso_atw_diff;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* compare periodic numbers; move indexes of the best neighbors to the first elements of ord_OM[] */
    // INCHI✔️✔️:     for (i = 1; i < num_OM; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         n_OM = at[at_no].neighbor[ord_OM[i]];
    // INCHI✔️✔️:         cur_value = (int) at[n_OM].el_number;
    // INCHI✔️✔️:         diff = cur_value - best_value;
    // INCHI✔️✔️:         if (( !cur_value && best_value ) || diff < 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:             best_value = cur_value;
    // INCHI✔️✔️:             ord_OM[0] = ord_OM[i];
    // INCHI✔️✔️:             num_best = 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (diff == 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* was '=', pointed by WDI */
    // INCHI✔️✔️:             ord_OM[num_best++] = ord_OM[i];
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     num_OM = num_best;
    // INCHI✔️✔️:     if (1 == num_OM)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ord_OM[0];
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* return any */
    // INCHI✔️✔️:     return ord_OM[0];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nFindOneOM
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: nFindOneOM
    // INCHI✔️✔️: typedef unsigned short AT_NUMB;
    // INCHI✔️✔️: typedef unsigned char U_CHAR;
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: nFindOneOM

    if num_om == 1 {
        return ord_om
            .first()
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds);
    }
    if num_om < 1 {
        return Ok(-1);
    }

    let candidate_count =
        usize::try_from(num_om).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if candidate_count > ord_om.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let center = atoms
        .get(usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let neighbor = |ordinal: i32| -> Result<usize, SourceHeapError> {
        let ordinal = usize::try_from(ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = center
            .neighbor
            .get(ordinal)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom = usize::from(*atom);
        atoms
            .get(atom)
            .map(|_| atom)
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };

    let mut num_best = 1_usize;
    let mut n_om = neighbor(ord_om[0])?;
    let mut best_value = i32::from(atoms[n_om].valence);
    for index in 1..candidate_count {
        n_om = neighbor(ord_om[index])?;
        let current_value = i32::from(atoms[n_om].valence);
        let difference = current_value - best_value;
        if difference < 0 {
            best_value = current_value;
            ord_om[0] = ord_om[index];
            num_best = 1;
        } else if difference == 0 {
            ord_om[num_best] = ord_om[index];
            num_best += 1;
        }
    }
    num_om = num_best as i32;
    if num_om == 1 {
        return Ok(ord_om[0]);
    }

    num_best = 1;
    n_om = neighbor(ord_om[0])?;
    best_value = i32::from(atoms[n_om].el_number);
    for index in 1..num_om as usize {
        n_om = neighbor(ord_om[index])?;
        let current_value = i32::from(atoms[n_om].el_number);
        let difference = current_value - best_value;
        if difference < 0 {
            best_value = current_value;
            ord_om[0] = ord_om[index];
            num_best = 1;
        } else if difference == 0 {
            ord_om[num_best] = ord_om[index];
            num_best += 1;
        }
    }
    num_om = num_best as i32;
    if num_om == 1 {
        return Ok(ord_om[0]);
    }
    if atoms[n_om].valence > 1 {
        return Ok(-1);
    }

    num_best = 1;
    n_om = neighbor(ord_om[0])?;
    best_value = i32::from(atoms[n_om].iso_atw_diff);
    for index in 1..num_om as usize {
        n_om = neighbor(ord_om[index])?;
        let current_value = i32::from(atoms[n_om].el_number);
        let difference = current_value - best_value;
        if (current_value == 0 && best_value != 0) || difference < 0 {
            best_value = current_value;
            ord_om[0] = ord_om[index];
            num_best = 1;
        } else if difference == 0 {
            ord_om[num_best] = ord_om[index];
            num_best += 1;
        }
    }
    num_om = num_best as i32;
    if num_om == 1 {
        return Ok(ord_om[0]);
    }
    Ok(ord_om[0])
}

pub(crate) fn remove_ion_pairs(
    num_atoms: i32,
    atoms: &mut [inp_ATOM],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:1049 remove_ion_pairs
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: remove_ion_pairs
    // INCHI✔️❌: int remove_ion_pairs( int num_atoms, inp_ATOM *at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int num_changes = 0;
    // INCHI✔️❌: #define MAX_NEIGH 6
    // INCHI✔️❌:
    // INCHI✔️❌:     int i, n, n2, i1, i2, i3, i4, type, chrg;
    // INCHI✔️❌:     int num_C_II = 0, num_C_plus = 0, num_C_minus = 0, num_N_plus = 0, num_N_minus = 0, num_O_plus = 0, num_O_minus = 0, num_All;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef FIX_P_IV_Plus_O_Minus
    // INCHI✔️❌:     int num_P_IV_plus = 0; /* added 2010-03-17 DT */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     inp_ATOM *a;
    // INCHI✔️❌:     /****** count candidates ********/
    // INCHI✔️❌:     for (i = 0, a = at; i < num_atoms; i++, a++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (1 == ( chrg = a->charge ) || -1 == chrg)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             switch (ion_el_group( a->el_number ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case EL_NUMBER_C:
    // INCHI✔️❌:                     if (chrg > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_C_plus++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_C_minus++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case EL_NUMBER_O:
    // INCHI✔️❌:                     if (chrg > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_O_plus++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_O_minus++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case EL_NUMBER_N:
    // INCHI✔️❌:                     if (chrg > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_N_plus++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_N_minus++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #ifdef FIX_P_IV_Plus_O_Minus
    // INCHI✔️❌:                     num_P_IV_plus += a->el_number != EL_NUMBER_N &&
    // INCHI✔️❌:                                      chrg == 1 &&
    // INCHI✔️❌:                                      a->valence == 4 &&
    // INCHI✔️❌:                                      a->chem_bonds_valence == 4; /* added 2010-03-17 DT */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (!chrg && a->chem_bonds_valence + NUMH( a, 0 ) == 2 &&
    // INCHI✔️❌:                   get_el_valence( a->el_number, 0, 0 ) == 4 &&
    // INCHI✔️❌:                   ion_el_group( a->el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_C_II++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     num_All = num_C_II + num_C_plus + num_C_minus + num_N_plus + num_N_minus + num_O_plus + num_O_minus;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* do not add num_P_IV_plus ! -- 2010-03-17 DT */
    // INCHI✔️❌:     if (!num_All)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**************************************************************************/
    // INCHI✔️❌:     /*************************** Terminal ion pairs ***************************/
    // INCHI✔️❌:     /**************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*-------------------------------------------------------------------------
    // INCHI✔️❌:     Pair type 1            N=N,P,As,Sb; O=O,S,Se,Te
    // INCHI✔️❌:     ===========
    // INCHI✔️❌:
    // INCHI✔️❌:     X              X     if X is another -O(-) then neutralize O(-)
    // INCHI✔️❌:     |              |     that has the smallest periodic table number
    // INCHI✔️❌:     O=N(+)-O(-) => O=N=O
    // INCHI✔️❌:     i    n
    // INCHI✔️❌:     --------------------------------------------------------------------------*/
    // INCHI✔️❌:
    // INCHI✔️❌:     for (type = 1; type <= 18; type++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( !type || 1 == type ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_atoms && 0 < num_N_plus && 0 < num_O_minus; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (1 == at[i].charge && 3 == nNoMetalNumBonds( at, i ) &&
    // INCHI✔️❌:                      4 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                      ion_el_group( at[i].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int num_OM = 0, ord_OM[3]; /* -O(-) */
    // INCHI✔️❌:                     int num_O = 0; /* =O    */
    // INCHI✔️❌:                     int num_O_other = 0;
    // INCHI✔️❌:                     for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = at[i].neighbor[i1];
    // INCHI✔️❌:                         if (1 == nNoMetalNumBonds( at, n ) && 0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                             ion_el_group( at[n].el_number) == EL_NUMBER_O) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (BOND_TYPE_SINGLE == at[i].bond_type[i1] &&
    // INCHI✔️❌:                                  -1 == at[n].charge)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ord_OM[num_OM++] = i1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else if (BOND_TYPE_DOUBLE == at[n].bond_type[0] &&
    // INCHI✔️❌:                                       0 == at[n].charge)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 num_O++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 num_O_other++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (num_OM > 0 && num_O > 0 && !num_O_other &&
    // INCHI✔️❌:                          0 <= ( i1 = nFindOneOM( at, i, ord_OM, num_OM ) ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* remove charges and increase bond order */
    // INCHI✔️❌:                         n = at[i].neighbor[i1];
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[i].charge--;
    // INCHI✔️❌:                         at[n].charge++;
    // INCHI✔️❌:                         at[i].radical = 0;
    // INCHI✔️❌:                         at[n].radical = 0;
    // INCHI✔️❌:                         num_changes++;
    // INCHI✔️❌:                         num_N_plus--;
    // INCHI✔️❌:                         num_O_minus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef FIX_P_IV_Plus_O_Minus
    // INCHI✔️❌:
    // INCHI✔️❌:             /*-------------------------------------------------------------------------
    // INCHI✔️❌:             Pair type 1a           P=P,As,Sb; O=O,S,Se,Te  -- added 2010-03-17
    // INCHI✔️❌:             =============
    // INCHI✔️❌:
    // INCHI✔️❌:             X              X     if X, Y, or Z is another -O(-) then neutralize O(-)
    // INCHI✔️❌:             |              |     that has the smallest periodic table number
    // INCHI✔️❌:             Y-P(+)-O(-) => Y-P=O
    // INCHI✔️❌:             |i   n         |
    // INCHI✔️❌:             Z              Z
    // INCHI✔️❌:
    // INCHI✔️❌:             --------------------------------------------------------------------------*/
    // INCHI✔️❌:
    // INCHI✔️❌:             for (i = 0; i < num_atoms && 0 < num_P_IV_plus /*&& 0 < num_N_plus*/ && 0 < num_O_minus; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (1 == at[i].charge && 4 == nNoMetalNumBonds( at, i ) &&
    // INCHI✔️❌:                      4 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                      at[i].el_number != EL_NUMBER_N && ion_el_group( at[i].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int num_OM = 0, ord_OM[4]; /* -O(-) */
    // INCHI✔️❌:                                                /*int num_O  = 0;*/ /* =O    */
    // INCHI✔️❌:                     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:                     for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = at[i].neighbor[i1];
    // INCHI✔️❌:                         if (1 == nNoMetalNumBonds( at, n ) && 0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                             ion_el_group( at[n].el_number) == EL_NUMBER_O) /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (BOND_TYPE_SINGLE == at[i].bond_type[i1] &&
    // INCHI✔️❌:                                  -1 == at[n].charge)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ord_OM[num_OM++] = i1;
    // INCHI✔️❌:                                 /*
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if ( BOND_TYPE_DOUBLE == at[n].bond_type[0] &&
    // INCHI✔️❌:                                 0                == at[n].charge       ) {
    // INCHI✔️❌:                                 num_O ++;
    // INCHI✔️❌:                                 */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (num_OM > 0 /*&& num_O > 0 && !num_O_other*/ &&
    // INCHI✔️❌:                          0 <= ( i1 = nFindOneOM( at, i, ord_OM, num_OM ) ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* remove charges and increase bond order */
    // INCHI✔️❌:                         n = at[i].neighbor[i1];
    // INCHI✔️❌:                         i2 = is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor;
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[i].charge--;
    // INCHI✔️❌:                         at[n].charge++;
    // INCHI✔️❌:                         at[i].radical = 0;
    // INCHI✔️❌:                         at[n].radical = 0;
    // INCHI✔️❌:                         num_changes++;
    // INCHI✔️❌:                         num_N_plus--;
    // INCHI✔️❌:                         num_O_minus--;
    // INCHI✔️❌:                         num_P_IV_plus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif /* FIX_P_IV_Plus_O_Minus */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /*-------------------------------------------------------------------------
    // INCHI✔️❌:         Terminal pair types: 2,3,4,5,6,7,8,9   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
    // INCHI✔️❌:         ====================================
    // INCHI✔️❌:         type #
    // INCHI✔️❌:         2 2:  O=N-C(II)-     => O=N#C-      N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
    // INCHI✔️❌:         3 9:  O=O(+)-C(-)(III) => O=O=C(IV)
    // INCHI✔️❌:         4 3:  O(-)-N(+)(IV)  => O=N(V)  (input structure has at least 1 double bond)
    // INCHI✔️❌:         5 4:  O(-)-O(+)(III) => O=O(IV)
    // INCHI✔️❌:         6 8:  O(-)-O-C(+)(III) => O=O=C(IV)
    // INCHI✔️❌:         7 5:  N(-)=N(+)(IV)  => N#N(V)    allow terminal H on N(-)
    // INCHI✔️❌:         8 6:  N(-)=O(+)(III) => N#O-
    // INCHI✔️❌:         9 7:  N(-)=C(+)(III) => N#C-
    // INCHI✔️❌:         --------------------------------------------------------------------------*/
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!type || (2 <= type && type <= 9)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_atoms && 0 < num_All; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (0 == at[i].charge && 1 == nNoMetalNumBonds( at, i ) && 2 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                      0 == num_of_H( at, i ) &&
    // INCHI✔️❌:                      ion_el_group( at[i].el_number ) == EL_NUMBER_O &&
    // INCHI✔️❌:                      0 <= ( i1 = nNoMetalNeighIndex( at, i ) ) &&
    // INCHI✔️❌:                      at[i].bond_type[i1] <= BOND_TYPE_TRIPLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* terminal O= */
    // INCHI✔️❌:                     n = at[i].neighbor[i1];
    // INCHI✔️❌:                     if (( !type || type == 2 ) && 0 < num_C_II)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* avoid alternating bonds */
    // INCHI✔️❌:                         if (0 == at[n].charge &&
    // INCHI✔️❌:                              2 == nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
    // INCHI✔️❌:                              0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                              ion_el_group( at[n].el_number ) == EL_NUMBER_N &&
    // INCHI✔️❌:                              0 <= ( i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
    // INCHI✔️❌:                              at[n].bond_type[i2] <= BOND_TYPE_TRIPLE)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* i2 = index of opposite to at[i] neighbor of at[n] */
    // INCHI✔️❌:                             /*i2 = (at[n].neighbor[0] == i);*/
    // INCHI✔️❌:                             n2 = at[n].neighbor[i2];
    // INCHI✔️❌:                             if (0 == at[n2].charge &&
    // INCHI✔️❌:                                  2 == at[n2].valence && 2 == at[n2].chem_bonds_valence &&
    // INCHI✔️❌:                                  0 == num_of_H( at, n2 ) &&
    // INCHI✔️❌:                                  ion_el_group( at[n2].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /*       i n n2     */
    // INCHI✔️❌:                                 /* found O=N-C(II)- */
    // INCHI✔️❌:                                 /* convert O=N-C(II)-     => O=N#C- */
    // INCHI✔️❌:
    // INCHI✔️❌:                                 i3 = ( at[n2].neighbor[0] != n ); /* index of at[n] neighbor of n2 */
    // INCHI✔️❌:                                 at[n].chem_bonds_valence = 5; /* N */
    // INCHI✔️❌:                                 at[n2].chem_bonds_valence = 4; /* C */
    // INCHI✔️❌:                                 at[n].bond_type[i2] = BOND_TYPE_TRIPLE;
    // INCHI✔️❌:                                 at[n2].bond_type[i3] = BOND_TYPE_TRIPLE;
    // INCHI✔️❌:                                 at[n2].radical = 0;
    // INCHI✔️❌:                                 num_changes++;
    // INCHI✔️❌:                                 num_C_II--;
    // INCHI✔️❌:                                 num_All--;
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (( !type || type == 3 ) && 0 < num_O_plus && 0 < num_C_minus)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (1 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
    // INCHI✔️❌:                              0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                              ion_el_group( at[n].el_number ) == EL_NUMBER_O &&
    // INCHI✔️❌:                              0 <= ( i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
    // INCHI✔️❌:                              at[n].bond_type[i2] <= BOND_TYPE_TRIPLE)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* found O=O(+)- */
    // INCHI✔️❌:                             /* i2 = index of opposite to at[i] neighbor of at[n] */
    // INCHI✔️❌:                             /*i2 = (at[n].neighbor[0] == i);*/
    // INCHI✔️❌:                             n2 = at[n].neighbor[i2];
    // INCHI✔️❌:                             if (-1 == at[n2].charge && 3 >= nNoMetalNumBonds( at, n2 ) && 3 == nNoMetalBondsValence( at, n2 ) + NUMH( at, n2 ) &&
    // INCHI✔️❌:                                  ion_el_group( at[n2].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /*             i n    n2        */
    // INCHI✔️❌:                                 /* found found O=O(+)-C(-)(III) */
    // INCHI✔️❌:                                 /* convert O=O(+)-C(-)(III)     => O=O=C(IV) */
    // INCHI✔️❌:                                 i3 = ( at[n2].neighbor[0] != n ); /* index of at[n] neighbor of n2 */
    // INCHI✔️❌:                                 at[n].charge--;
    // INCHI✔️❌:                                 at[n2].charge++;
    // INCHI✔️❌:                                 at[n].chem_bonds_valence += 1; /* =O- => =O= */
    // INCHI✔️❌:                                 at[n2].chem_bonds_valence += 1; /* -C  => =C  */
    // INCHI✔️❌:                                 at[n].bond_type[i2] = BOND_TYPE_DOUBLE;
    // INCHI✔️❌:                                 at[n2].bond_type[i3] = BOND_TYPE_DOUBLE;
    // INCHI✔️❌:                                 num_changes++;
    // INCHI✔️❌:                                 num_O_plus--;
    // INCHI✔️❌:                                 num_C_minus--;
    // INCHI✔️❌:                                 num_All -= 2;
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 else if (-1 == at[i].charge &&
    // INCHI✔️❌:                           0 < num_O_minus + num_N_minus &&
    // INCHI✔️❌:                           0 < num_N_plus + num_O_plus + num_C_plus &&
    // INCHI✔️❌:                           1 == nNoMetalNumBonds( at, i ) && 1 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                           0 == num_of_H( at, i ) &&
    // INCHI✔️❌:                           ion_el_group( at[i].el_number ) == EL_NUMBER_O &&
    // INCHI✔️❌:                           0 <= ( i1 = nNoMetalNeighIndex( at, i ) ) &&
    // INCHI✔️❌:                           at[i].bond_type[i1] <= BOND_TYPE_TRIPLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* terminal O(-)- */
    // INCHI✔️❌:                     n = at[i].neighbor[i1];
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (( !type || type == 4 ) && 0 < num_O_minus && 0 < num_N_plus && /* O(-)-N(+)(IV) */
    // INCHI✔️❌:                          1 == at[n].charge && 3 >= nNoMetalNumBonds( at, n ) && 4 == nNoMetalBondsValence( at, n ) &&
    // INCHI✔️❌:                          0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                          ion_el_group( at[n].el_number ) == EL_NUMBER_N /* except >O(+)- */
    // INCHI✔️❌:                          )
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* found O(-)-N(+)(IV) */
    // INCHI✔️❌:                         /* convert O(-)-N(+)(IV)     => O=N(V)  */
    // INCHI✔️❌:
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
    // INCHI✔️❌:                         at[i].charge++;
    // INCHI✔️❌:                         at[n].charge--;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         num_changes++;
    // INCHI✔️❌:                         num_O_minus--;
    // INCHI✔️❌:                         num_N_plus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (( !type || type == 5 ) && 0 < num_O_minus && 0 < num_O_plus &&/* O(-)-O(+)(III) */
    // INCHI✔️❌:                          1 == at[n].charge && 3 >= nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
    // INCHI✔️❌:                          0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                          ion_el_group( at[n].el_number ) == EL_NUMBER_O /* except >O(+)- */
    // INCHI✔️❌:                          )
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* found  O(+)(III) */
    // INCHI✔️❌:                         /* convert O(-)-O(+)(III)    => O=O(IV) */
    // INCHI✔️❌:
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
    // INCHI✔️❌:                         at[i].charge++;
    // INCHI✔️❌:                         at[n].charge--;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         num_changes++;
    // INCHI✔️❌:                         num_O_minus--;
    // INCHI✔️❌:                         num_O_plus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* i    n n2        */
    // INCHI✔️❌:                     if (( !type || type == 6 ) && /* O(-)-O-C(+)(III) */
    // INCHI✔️❌:                          0 < num_O_minus && 0 < num_C_plus &&
    // INCHI✔️❌:                          0 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) && 2 == nNoMetalBondsValence( at, n ) &&
    // INCHI✔️❌:                          0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                          ion_el_group( at[n].el_number ) == EL_NUMBER_O &&
    // INCHI✔️❌:                          0 <= ( i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
    // INCHI✔️❌:                          at[n].bond_type[i2] <= BOND_TYPE_TRIPLE)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* found O(-)-O- */
    // INCHI✔️❌:                         /* i2 = index of opposite to at[i] neighbor of at[n] */
    // INCHI✔️❌:                         /*i2 = (at[n].neighbor[0] == i);*/
    // INCHI✔️❌:                         n2 = at[n].neighbor[i2];
    // INCHI✔️❌:                         if (1 == at[n2].charge && 3 >= nNoMetalNumBonds( at, n2 ) &&
    // INCHI✔️❌:                              3 == nNoMetalBondsValence( at, n2 ) + NUMH( at, n2 ) &&
    // INCHI✔️❌:                              ion_el_group( at[n2].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*       i    n n2  */
    // INCHI✔️❌:                             /* found O(-)-O-C(+)(III) */
    // INCHI✔️❌:                             /* convert O(-)-O-C(+)(III)     => O=O=C(IV) */
    // INCHI✔️❌:                             /*i3 = (at[n2].neighbor[0] != n);*/ /* i3 = index of at[n] neighbor of at[n2] */
    // INCHI✔️❌:                             i3 = (int) ( is_in_the_list( at[n2].neighbor, (AT_NUMB) n, at[n2].valence ) - at[n2].neighbor );
    // INCHI✔️❌:                             /*i4 = index of at[i] in the adjacency list of at[n] */
    // INCHI✔️❌:                             i4 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].charge++;
    // INCHI✔️❌:                             at[n2].charge--;
    // INCHI✔️❌:                             at[i].chem_bonds_valence += 1; /* O-  => O=  */
    // INCHI✔️❌:                             at[n].chem_bonds_valence += 2; /* -O- => =O= */
    // INCHI✔️❌:                             at[n2].chem_bonds_valence += 1; /* -C  => =C  */
    // INCHI✔️❌:                             at[i].bond_type[i1] = BOND_TYPE_DOUBLE;
    // INCHI✔️❌:                             at[n].bond_type[i4] = BOND_TYPE_DOUBLE;
    // INCHI✔️❌:                             at[n].bond_type[i2] = BOND_TYPE_DOUBLE;
    // INCHI✔️❌:                             at[n2].bond_type[i3] = BOND_TYPE_DOUBLE;
    // INCHI✔️❌:                             num_changes++;
    // INCHI✔️❌:                             num_O_minus--;
    // INCHI✔️❌:                             num_C_plus--;
    // INCHI✔️❌:                             num_All -= 2;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (-1 == at[i].charge && 0 < num_N_minus && 0 < num_N_plus + num_O_plus + num_C_plus &&
    // INCHI✔️❌:                           1 == nNoMetalNumBonds( at, i ) && 2 == nNoMetalBondsValence( at, i ) + NUMH( at, i ) &&
    // INCHI✔️❌:                           /*0 == num_of_H( at, i ) &&*/
    // INCHI✔️❌:                           ion_el_group( at[i].el_number ) == EL_NUMBER_N &&
    // INCHI✔️❌:                           0 <= ( i1 = nNoMetalNeighIndex( at, i ) ) &&
    // INCHI✔️❌:                           at[i].bond_type[i1] <= BOND_TYPE_TRIPLE)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* terminal N(-)= */
    // INCHI✔️❌:                     n = at[i].neighbor[i1 = 0];
    // INCHI✔️❌:                     if (( !type || type == 7 ) && 0 < num_N_plus && /* N(-)=N(+)(IV) */
    // INCHI✔️❌:                          1 == at[n].charge && 3 >= nNoMetalNumBonds( at, n ) && 4 == nNoMetalBondsValence( at, n ) &&
    // INCHI✔️❌:                          0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                          ion_el_group( at[n].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* found N(-)-N(+)(IV) */
    // INCHI✔️❌:                         /* convert N(-)=N(+)(IV)     => N#N(V)  */
    // INCHI✔️❌:
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
    // INCHI✔️❌:                         at[i].charge++;
    // INCHI✔️❌:                         at[n].charge--;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         num_changes++;
    // INCHI✔️❌:                         num_N_minus--;
    // INCHI✔️❌:                         num_N_plus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (( !type || type == 8 ) && 0 < num_O_plus && /* N(-)=O(+)(III) */
    // INCHI✔️❌:                          1 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
    // INCHI✔️❌:                          0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                          ion_el_group( at[n].el_number ) == EL_NUMBER_O)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* found N(-)-O(+)(III) */
    // INCHI✔️❌:                         /* convert N(-)=O(+)(III)    => N#O(IV)- */
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
    // INCHI✔️❌:                         at[i].charge++;
    // INCHI✔️❌:                         at[n].charge--;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         num_changes++;
    // INCHI✔️❌:                         num_N_minus--;
    // INCHI✔️❌:                         num_O_plus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (( !type || type == 9 ) && 0 < num_C_plus && /* N(-)=C(+)(III) */
    // INCHI✔️❌:                          1 == at[n].charge && 2 == at[n].valence && 3 == at[n].chem_bonds_valence &&
    // INCHI✔️❌:                          0 == num_of_H( at, n ) &&
    // INCHI✔️❌:                          ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* found N(-)=C(+)(III) */
    // INCHI✔️❌:                         /* convert N(-)=C(+)(III)    => N#C(IV)- */
    // INCHI✔️❌:
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
    // INCHI✔️❌:                         at[i].charge++;
    // INCHI✔️❌:                         at[n].charge--;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         num_changes++;
    // INCHI✔️❌:                         num_N_minus--;
    // INCHI✔️❌:                         num_C_plus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /**************************************************************************/
    // INCHI✔️❌:         /*********************** NON-Terminal ion pairs ***************************/
    // INCHI✔️❌:         /**************************************************************************/
    // INCHI✔️❌:         /*-------------------------------------------------------------------------
    // INCHI✔️❌:         Non-Terminal pair types: 10,11,12,13,14   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
    // INCHI✔️❌:         ========================================
    // INCHI✔️❌:
    // INCHI✔️❌:         10:  N(+)(IV)-C(-)(III)     => N(V)=C(IV)  (N has 3 or 2 bonds)
    // INCHI✔️❌:         11:  N(+)(IV)=C(-)(III)     => N(V)#C(IV)  (N has 3 or 2 bonds)
    // INCHI✔️❌:         12:  N(+)(IV)-N(-)(II)      => N(V)=N(III) (allow terminal H on N(-))
    // INCHI✔️❌:         13: -O(+)-C(-)(III)         => -O=C-
    // INCHI✔️❌:         14: -O(+)=C(-)(III)         => -O#C-
    // INCHI✔️❌:         15:  O(+)(III)-N(-)(II)     => O(IV)=N(III) (allow terminal H on N(-))
    // INCHI✔️❌:         --------------------------------------------------------------------------*/
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!type || (10 <= type && type <= 15)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_atoms && 0 < num_All; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (1 == at[i].charge &&
    // INCHI✔️❌:                      0 < num_N_plus + num_O_plus && 0 < num_C_minus + num_N_minus &&
    // INCHI✔️❌:                      4 >= nNoMetalNumBonds( at, i ) && 4 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                      0 == num_of_H( at, i ) &&
    // INCHI✔️❌:                      ion_el_group( at[i].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* found non-terminal N(+)(IV) */
    // INCHI✔️❌:                     if (( !type || 10 == type ) && 0 < num_N_plus && 0 < num_C_minus)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_neigh = 0, pos_neigh = -1;
    // INCHI✔️❌:                         for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n = at[i].neighbor[i1];
    // INCHI✔️❌:                             if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
    // INCHI✔️❌:                                  /*0 == at[n].num_H &&*/
    // INCHI✔️❌:                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
    // INCHI✔️❌:                                  ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* found N(+)(IV)-C(-)(III); prepare conversion to N(V)=C(IV) */
    // INCHI✔️❌:                                 num_neigh++;
    // INCHI✔️❌:                                 pos_neigh = i1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         i1 = pos_neigh;
    // INCHI✔️❌:                         if (1 == num_neigh &&
    // INCHI✔️❌:                              at[i].bond_type[i1] <= BOND_TYPE_TRIPLE &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, i, n = at[i].neighbor[i1] ) &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, n, i ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*n = at[i].neighbor[i1=pos_neigh];*/
    // INCHI✔️❌:                             i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].charge--;
    // INCHI✔️❌:                             at[n].charge++;
    // INCHI✔️❌:                             at[i].chem_bonds_valence++;
    // INCHI✔️❌:                             at[n].chem_bonds_valence++;
    // INCHI✔️❌:                             at[i].bond_type[i1] ++;
    // INCHI✔️❌:                             at[n].bond_type[i2] ++;
    // INCHI✔️❌:                             num_changes++;
    // INCHI✔️❌:                             num_C_minus--;
    // INCHI✔️❌:                             num_N_plus--;
    // INCHI✔️❌:                             num_All -= 2;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (( !type || 11 == type ) && 0 < num_N_plus && 0 < num_C_minus)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_neigh = 0, pos_neigh = -1;
    // INCHI✔️❌:                         for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n = at[i].neighbor[i1];
    // INCHI✔️❌:                             if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
    // INCHI✔️❌:                                  /*0 == at[n].num_H &&*/
    // INCHI✔️❌:                                  at[i].bond_type[i1] == BOND_TYPE_DOUBLE &&
    // INCHI✔️❌:                                  ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* found N(+)(IV)=C(-)(III); prepare conversion to N(V)#C(IV) */
    // INCHI✔️❌:                                 num_neigh++;
    // INCHI✔️❌:                                 pos_neigh = i1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (1 == num_neigh &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, n, i))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*n = at[i].neighbor[i1=pos_neigh];*/
    // INCHI✔️❌:                             i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].charge--;
    // INCHI✔️❌:                             at[n].charge++;
    // INCHI✔️❌:                             at[i].chem_bonds_valence++;
    // INCHI✔️❌:                             at[n].chem_bonds_valence++;
    // INCHI✔️❌:                             at[i].bond_type[i1] ++;
    // INCHI✔️❌:                             at[n].bond_type[i2] ++;
    // INCHI✔️❌:                             num_changes++;
    // INCHI✔️❌:                             num_C_minus--;
    // INCHI✔️❌:                             num_N_plus--;
    // INCHI✔️❌:                             num_All -= 2;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (!type || (12 == type && 0 < num_N_plus && 0 < num_N_minus)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_neigh = 0, pos_neigh = -1;
    // INCHI✔️❌:                         for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n = at[i].neighbor[i1];
    // INCHI✔️❌:                             if (-1 == at[n].charge && 2 >= nNoMetalNumBonds( at, n ) &&
    // INCHI✔️❌:                                  2 == nNoMetalBondsValence( at, n ) + NUMH( at, n ) &&
    // INCHI✔️❌:                                  /*0 == num_of_H( at, n ) &&*/
    // INCHI✔️❌:                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
    // INCHI✔️❌:                                  ion_el_group( at[n].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* found N(+)(IV)=N(-)(II); prepare conversion to N(V)#N(III) */
    // INCHI✔️❌:                                 num_neigh++;
    // INCHI✔️❌:                                 pos_neigh = i1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (1 == num_neigh &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, n, i))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*n = at[i].neighbor[i1=pos_neigh];*/
    // INCHI✔️❌:                             i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].charge--;
    // INCHI✔️❌:                             at[n].charge++;
    // INCHI✔️❌:                             at[i].chem_bonds_valence++;
    // INCHI✔️❌:                             at[n].chem_bonds_valence++;
    // INCHI✔️❌:                             at[i].bond_type[i1] ++;
    // INCHI✔️❌:                             at[n].bond_type[i2] ++;
    // INCHI✔️❌:                             num_changes++;
    // INCHI✔️❌:                             num_N_minus--;
    // INCHI✔️❌:                             num_N_plus--;
    // INCHI✔️❌:                             num_All -= 2;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (1 == at[i].charge &&
    // INCHI✔️❌:                           0 < num_O_plus && 0 < num_C_minus + num_N_minus &&
    // INCHI✔️❌:                           3 >= nNoMetalNumBonds( at, i ) && 3 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                           0 == num_of_H( at, i ) &&
    // INCHI✔️❌:                           ion_el_group( at[i].el_number ) == EL_NUMBER_O)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* found non-terminal O(+)(III) */
    // INCHI✔️❌:                     if (( !type || 13 == type ) && 0 < num_C_minus)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_neigh = 0, pos_neigh = -1;
    // INCHI✔️❌:                         for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n = at[i].neighbor[i1];
    // INCHI✔️❌:                             if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
    // INCHI✔️❌:                                  /*0 == at[n].num_H &&*/
    // INCHI✔️❌:                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
    // INCHI✔️❌:                                  ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* found O(+)(III)-C(-)(II); prepare conversion to O(IV)=C(IV) */
    // INCHI✔️❌:                                 num_neigh++;
    // INCHI✔️❌:                                 pos_neigh = i1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (1 == num_neigh &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, n, i))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*n = at[i].neighbor[i1=pos_neigh];*/
    // INCHI✔️❌:                             i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].charge--;
    // INCHI✔️❌:                             at[n].charge++;
    // INCHI✔️❌:                             at[i].chem_bonds_valence++;
    // INCHI✔️❌:                             at[n].chem_bonds_valence++;
    // INCHI✔️❌:                             at[i].bond_type[i1] ++;
    // INCHI✔️❌:                             at[n].bond_type[i2] ++;
    // INCHI✔️❌:                             num_changes++;
    // INCHI✔️❌:                             num_C_minus--;
    // INCHI✔️❌:                             num_O_plus--;
    // INCHI✔️❌:                             num_All -= 2;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (( !type || 14 == type ) && 0 < num_C_minus)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_neigh = 0, pos_neigh = -1;
    // INCHI✔️❌:                         for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n = at[i].neighbor[i1];
    // INCHI✔️❌:                             if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
    // INCHI✔️❌:                                  /*0 == at[n].num_H &&*/
    // INCHI✔️❌:                                  at[i].bond_type[i1] == BOND_TYPE_DOUBLE &&
    // INCHI✔️❌:                                  ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* found O(+)(III)=C(-)(III); prepare conversion to O(IV)#C(IV) */
    // INCHI✔️❌:                                 num_neigh++;
    // INCHI✔️❌:                                 pos_neigh = i1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (1 == num_neigh &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, n, i))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*n = at[i].neighbor[i1=pos_neigh];*/
    // INCHI✔️❌:                             i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].charge--;
    // INCHI✔️❌:                             at[n].charge++;
    // INCHI✔️❌:                             at[i].chem_bonds_valence++;
    // INCHI✔️❌:                             at[n].chem_bonds_valence++;
    // INCHI✔️❌:                             at[i].bond_type[i1] ++;
    // INCHI✔️❌:                             at[n].bond_type[i2] ++;
    // INCHI✔️❌:                             num_changes++;
    // INCHI✔️❌:                             num_C_minus--;
    // INCHI✔️❌:                             num_O_plus--;
    // INCHI✔️❌:                             num_All -= 2;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (( !type || 15 == type ) && 0 < num_N_minus)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_neigh = 0, pos_neigh = -1;
    // INCHI✔️❌:                         for (i1 = 0; i1 < at[i].valence; i1++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n = at[i].neighbor[i1];
    // INCHI✔️❌:                             if (-1 == at[n].charge && 2 >= nNoMetalNumBonds( at, n ) &&
    // INCHI✔️❌:                                  2 == nNoMetalBondsValence( at, n ) + NUMH( at, n ) &&
    // INCHI✔️❌:                                  /*0 == num_of_H( at, n ) &&*/
    // INCHI✔️❌:                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
    // INCHI✔️❌:                                  ion_el_group( at[n].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* found O(+)(III)=N(-)(II); prepare conversion to O(IV)#N(III) */
    // INCHI✔️❌:                                 num_neigh++;
    // INCHI✔️❌:                                 pos_neigh = i1;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (1 == num_neigh &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
    // INCHI✔️❌:                              !has_other_ion_neigh( at, n, i))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*n = at[i].neighbor[i1=pos_neigh];*/
    // INCHI✔️❌:                             i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].charge--;
    // INCHI✔️❌:                             at[n].charge++;
    // INCHI✔️❌:                             at[i].chem_bonds_valence++;
    // INCHI✔️❌:                             at[n].chem_bonds_valence++;
    // INCHI✔️❌:                             at[i].bond_type[i1] ++;
    // INCHI✔️❌:                             at[n].bond_type[i2] ++;
    // INCHI✔️❌:                             num_changes++;
    // INCHI✔️❌:                             num_N_minus--;
    // INCHI✔️❌:                             num_O_plus--;
    // INCHI✔️❌:                             num_All -= 2;
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /**************************************************************************/
    // INCHI✔️❌:         /*********************** NON-Terminal ion triples *************************/
    // INCHI✔️❌:         /**************************************************************************/
    // INCHI✔️❌:         /*-------------------------------------------------------------------------
    // INCHI✔️❌:         Non-Terminal triple types: 16, 17, 18   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
    // INCHI✔️❌:         ========================================
    // INCHI✔️❌:         16: C(+)(III)-O-N(-)(II)  => C(IV)=O=N(III)  (allow terminal H on N(-))
    // INCHI✔️❌:         |                     |
    // INCHI✔️❌:         17: C(+)(III)-N-C(-)(III)  => C(IV)=N=C(IV)
    // INCHI✔️❌:
    // INCHI✔️❌:         18: C(-)(III)-N=C(+)(III)  => C(IV)=N#C(IV)   (may have two or no charges)
    // INCHI✔️❌:         C(IV)=N-C(II)          => C(IV)=N#C(IV)
    // INCHI✔️❌:
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (( !type || 16 == type ) && 0 < num_C_plus && 0 < num_N_minus)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int m[2], j[2], k;
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (0 == at[i].charge && 2 == nNoMetalNumBonds( at, i ) && 2 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                      0 == num_of_H( at, i ) &&
    // INCHI✔️❌:                      0 <= ( j[0] = nNoMetalNeighIndex( at, i ) ) &&
    // INCHI✔️❌:                      at[m[0] = at[i].neighbor[j[0]]].charge &&
    // INCHI✔️❌:                      0 <= ( j[1] = nNoMetalOtherNeighIndex( at, i, m[0] ) ) &&
    // INCHI✔️❌:                      0 == at[m[0]].charge + at[m[1] = at[i].neighbor[j[1]]].charge &&
    // INCHI✔️❌:                      5 >= nNoMetalBondsValence( at, m[0] ) + nNoMetalBondsValence( at, m[1] ) &&
    // INCHI✔️❌:                      /*5 >= at[m[0]].chem_bonds_valence + at[m[1]].chem_bonds_valence &&*/
    // INCHI✔️❌:                      ion_el_group( at[i].el_number ) == EL_NUMBER_O)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* found non-terminal A(+)-O-B(-); chem_bond_val of A+B <= 5 */
    // INCHI✔️❌:                     int n_N = -1, n_C = -1, i_C = -1;
    // INCHI✔️❌:                     for (k = 0; k < 2; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = m[k];
    // INCHI✔️❌:                         if (-1 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) + NUMH( at, n ) &&
    // INCHI✔️❌:                              /*0 == num_of_H( at, n ) &&*/
    // INCHI✔️❌:                              ion_el_group( at[n].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n_N = n;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else if (1 == at[n].charge && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
    // INCHI✔️❌:                                   ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n_C = n;
    // INCHI✔️❌:                             i_C = k;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (n_C < 0 || n_N < 0 ||
    // INCHI✔️❌:                          has_other_ion_in_sphere_2( at, n_C, n_N) ||
    // INCHI✔️❌:                          has_other_ion_in_sphere_2( at, n_N, n_C))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* C(+)(III)-O-N(-)(II)  => C(IV)=O=N(III) */
    // INCHI✔️❌:                     for (k = 0; k < 2; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = k ? n_C : n_N;
    // INCHI✔️❌:                         i1 = k ? j[i_C] : j[1 - i_C];
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].charge += ( k ? -1 : 1 );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     num_changes++;
    // INCHI✔️❌:                     num_N_minus--;
    // INCHI✔️❌:                     num_C_plus--;
    // INCHI✔️❌:                     num_All -= 2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (( !type || 17 == type ) && 0 < num_C_plus && 0 < num_C_minus)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int m[3], c[3], j[3], k;
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (0 == at[i].charge && 3 == nNoMetalNumBonds( at, i ) && 3 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                      0 == num_of_H( at, i ) &&
    // INCHI✔️❌:                      0 <= ( j[0] = nNoMetalNeighIndex( at, i ) ) &&
    // INCHI✔️❌:                      0 <= ( j[1] = nNoMetalOtherNeighIndex( at, i, m[0] = at[i].neighbor[j[0]] ) ) &&
    // INCHI✔️❌:                      0 <= ( j[2] = nNoMetalOtherNeighIndex2( at, i, m[0], m[1] = at[i].neighbor[j[1]] ) ) &&
    // INCHI✔️❌:                      1 == !( c[0] = at[m[0]].charge )
    // INCHI✔️❌:                      + !( c[1] = at[m[1]].charge )
    // INCHI✔️❌:                      + !( c[2] = at[m[2] = at[i].neighbor[j[2]]].charge ) &&
    // INCHI✔️❌:                      0 == c[0] + c[1] + c[2] &&
    // INCHI✔️❌:                      2 == ( 3 == ( c[0] ? at[m[0]].chem_bonds_valence + NUMH( at, m[0] ) : 0 ) )
    // INCHI✔️❌:                      + ( 3 == ( c[1] ? at[m[1]].chem_bonds_valence + NUMH( at, m[1] ) : 0 ) )
    // INCHI✔️❌:                      + ( 3 == ( c[2] ? at[m[2]].chem_bonds_valence + NUMH( at, m[2] ) : 0 ) ) &&
    // INCHI✔️❌:                      ion_el_group( at[i].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* found non-terminal A(+)-O-B(-) */
    // INCHI✔️❌:                     int n_Cp = -1, n_Cm = -1, i_Cp = -1, i_Cm = -1; /* p = positive, m = negatice ion C */
    // INCHI✔️❌:                     for (k = 0; k < 3; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (c[k])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n = m[k];
    // INCHI✔️❌:                             if (-1 == at[n].charge &&
    // INCHI✔️❌:                                  ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 n_Cm = n;
    // INCHI✔️❌:                                 i_Cm = k;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else if (1 == at[n].charge &&
    // INCHI✔️❌:                                       ion_el_group( at[n].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 n_Cp = n;
    // INCHI✔️❌:                                 i_Cp = k;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (n_Cp < 0 || n_Cm < 0 ||
    // INCHI✔️❌:                          has_other_ion_in_sphere_2( at, n_Cp, n_Cm) ||
    // INCHI✔️❌:                          has_other_ion_in_sphere_2( at, n_Cm, n_Cp))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /*           |                     |       */
    // INCHI✔️❌:                     /* C(+)(III)-N-C(-)(III)  => C(IV)=N=C(IV) */
    // INCHI✔️❌:                     for (k = 0; k < 2; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = k ? n_Cp : n_Cm;
    // INCHI✔️❌:                         i1 = k ? j[i_Cp] : j[i_Cm];
    // INCHI✔️❌:                         i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                         at[i].bond_type[i1] ++;
    // INCHI✔️❌:                         at[n].bond_type[i2] ++;
    // INCHI✔️❌:                         at[i].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].chem_bonds_valence++;
    // INCHI✔️❌:                         at[n].charge += ( k ? -1 : 1 );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     num_changes++;
    // INCHI✔️❌:                     num_C_minus--;
    // INCHI✔️❌:                     num_C_plus--;
    // INCHI✔️❌:                     num_All -= 2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (( !type || 18 == type ) && ( (0 < num_C_plus && 0 < num_C_minus) || 0 < num_C_II )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int m[2], v[2], j[2], k;
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (0 == at[i].charge && 2 == nNoMetalNumBonds( at, i ) && 3 == nNoMetalBondsValence( at, i ) &&
    // INCHI✔️❌:                      0 == num_of_H( at, i ) &&
    // INCHI✔️❌:                      0 <= ( j[0] = nNoMetalNeighIndex( at, i ) ) &&
    // INCHI✔️❌:                      0 <= ( j[1] = nNoMetalOtherNeighIndex( at, i, m[0] = at[i].neighbor[j[0]] ) ) &&
    // INCHI✔️❌:                      0 == at[m[0]].charge
    // INCHI✔️❌:                      + at[m[1] = at[i].neighbor[j[1]]].charge &&
    // INCHI✔️❌:                      6 == ( v[0] = at[m[0]].chem_bonds_valence + NUMH( at, m[0] ) )
    // INCHI✔️❌:                      + ( v[1] = at[m[1]].chem_bonds_valence + NUMH( at, m[1] ) ) &&
    // INCHI✔️❌:                      2 >= abs( v[0] - v[1] ) &&
    // INCHI✔️❌:                      ion_el_group( at[i].el_number ) == EL_NUMBER_N &&
    // INCHI✔️❌:                      ion_el_group( at[m[0]].el_number ) == EL_NUMBER_C &&
    // INCHI✔️❌:                      ion_el_group( at[m[1]].el_number ) == EL_NUMBER_C)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*                    n_Cm      i n_Cp */
    // INCHI✔️❌:                     /* found non-terminal C(-)(III)-N=C(+)(III) or C(IV)=N-C(II): Cm-N-Cp */
    // INCHI✔️❌:                     /* convert to C(IV)=N#C(IV) */
    // INCHI✔️❌:                     int n_Cp = -1, n_Cm = -1, i_Cp = -1, i_Cm = -1; /* p = positive, m = negatice ion C */
    // INCHI✔️❌:                     for (k = 0; k < 2; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = m[k];
    // INCHI✔️❌:                         if (v[k] == 4 || (v[k] == 3 && at[i].bond_type[j[k]] == BOND_TYPE_SINGLE)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n_Cm = n;
    // INCHI✔️❌:                             i_Cm = k;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else if (v[k] == 2 || (v[k] == 3 && at[i].bond_type[j[k]] == BOND_TYPE_DOUBLE)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             n_Cp = n;
    // INCHI✔️❌:                             i_Cp = k;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (n_Cp < 0 || n_Cm < 0 || at[n_Cp].valence + NUMH( at, n_Cp ) != 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue; /* guarantees at[n_Cp].valence <= 2 */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (v[i_Cp] == 2 || !at[n_Cp].charge)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (at[n_Cp].valence == 2)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* neighbor of at[n_Cp] opposite to at[i] */
    // INCHI✔️❌:                             k = at[n_Cp].neighbor[at[n_Cp].neighbor[0] == i];
    // INCHI✔️❌:                             if (ion_el_group( at[k].el_number ) == EL_NUMBER_N)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (at[n_Cp].charge)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (has_other_ion_in_sphere_2( at, n_Cp, n_Cm) ||
    // INCHI✔️❌:                              has_other_ion_in_sphere_2( at, n_Cm, n_Cp))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         continue; /* unknown case */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /*                                         */
    // INCHI✔️❌:                     /* C(-)(III)-N=C(+)(III)  => C(IV)=N#C(IV) */
    // INCHI✔️❌:                     /* C(IV)=N-C(II)          => C(IV)=N#C(IV) */
    // INCHI✔️❌:                     if (at[n_Cp].charge)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_C_minus--;
    // INCHI✔️❌:                         num_C_plus--;
    // INCHI✔️❌:                         num_All -= 2;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_C_II--;
    // INCHI✔️❌:                         num_All--;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     for (k = 0; k < 2; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n = k ? n_Cp : n_Cm;
    // INCHI✔️❌:                         i3 = k ? i_Cp : i_Cm; /* added to fix the bug */
    // INCHI✔️❌:                                               /*i1 = k? j[i_Cp] : j[i_Cm];*/ /* replaced with next line */
    // INCHI✔️❌:                         i1 = j[i3];
    // INCHI✔️❌:                         if (v[i3 /*was i1*/] < 4)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* WDI found a bug here: bounds violation */
    // INCHI✔️❌:                             int delta = 4 - v[i3 /*was i1*/];
    // INCHI✔️❌:                             i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
    // INCHI✔️❌:                             at[i].bond_type[i1] += delta;
    // INCHI✔️❌:                             at[n].bond_type[i2] += delta;
    // INCHI✔️❌:                             at[i].chem_bonds_valence += delta;
    // INCHI✔️❌:                             at[n].chem_bonds_valence += delta;
    // INCHI✔️❌:                             at[n].charge = 0;
    // INCHI✔️❌:                             at[n].radical = 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     at[i].charge = 0;
    // INCHI✔️❌:                     at[i].radical = 0;
    // INCHI✔️❌:                     num_changes++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_changes;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: remove_ion_pairs
    // END INCHI C FUNCTION: remove_ion_pairs
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: remove_ion_pairs
    // INCHI✔️❌: #define FIX_P_IV_Plus_O_Minus
    // INCHI✔️❌: #define MAX_NEIGH 6
    // INCHI✔️❌: #define EL_NUMBER_C ((U_CHAR)6)
    // INCHI✔️❌: #define EL_NUMBER_N ((U_CHAR)7)
    // INCHI✔️❌: #define EL_NUMBER_O ((U_CHAR)8)
    // INCHI✔️❌: #define BOND_TYPE_SINGLE 1
    // INCHI✔️❌: #define BOND_TYPE_DOUBLE 2
    // INCHI✔️❌: #define BOND_TYPE_TRIPLE 3
    // INCHI✔️❌: #define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
    // INCHI✔️❌: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: remove_ion_pairs

    let atom_count =
        usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let atoms = atoms
        .get_mut(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let numh = |atom: &inp_ATOM| {
        i32::from(atom.num_H)
            + i32::from(atom.num_iso_H[0])
            + i32::from(atom.num_iso_H[1])
            + i32::from(atom.num_iso_H[2])
    };
    let reciprocal =
        |atoms: &[inp_ATOM], atom: usize, neighbor: usize| -> Result<usize, SourceHeapError> {
            let target = atoms.get(atom).ok_or(SourceHeapError::PointerOutOfBounds)?;
            is_in_the_list(
                Some(&target.neighbor),
                neighbor as AT_NUMB,
                i32::from(target.valence),
            )?
            .ok_or(SourceHeapError::PointerOutOfBounds)
        };

    let mut num_changes = 0_i32;
    let mut num_c_ii = 0_i32;
    let mut num_c_plus = 0_i32;
    let mut num_c_minus = 0_i32;
    let mut num_n_plus = 0_i32;
    let mut num_n_minus = 0_i32;
    let mut num_o_plus = 0_i32;
    let mut num_o_minus = 0_i32;
    let mut num_p_iv_plus = 0_i32;

    for atom in atoms.iter() {
        let charge = i32::from(atom.charge);
        if charge == 1 || charge == -1 {
            match ion_el_group(i32::from(atom.el_number)) {
                6 => {
                    if charge > 0 {
                        num_c_plus += 1;
                    } else {
                        num_c_minus += 1;
                    }
                }
                8 => {
                    if charge > 0 {
                        num_o_plus += 1;
                    } else {
                        num_o_minus += 1;
                    }
                }
                7 => {
                    if charge > 0 {
                        num_n_plus += 1;
                    } else {
                        num_n_minus += 1;
                    }
                    num_p_iv_plus += i32::from(
                        atom.el_number != 7
                            && charge == 1
                            && atom.valence == 4
                            && atom.chem_bonds_valence == 4,
                    );
                }
                _ => {}
            }
        } else if charge == 0
            && i32::from(atom.chem_bonds_valence) + numh(atom) == 2
            && get_el_valence(i32::from(atom.el_number), 0, 0)? == 4
            && ion_el_group(i32::from(atom.el_number)) == 6
        {
            num_c_ii += 1;
        }
    }
    let mut num_all =
        num_c_ii + num_c_plus + num_c_minus + num_n_plus + num_n_minus + num_o_plus + num_o_minus;
    if num_all == 0 {
        return Ok(0);
    }

    for pair_type in 1_i32..=18 {
        if pair_type == 1 {
            for i in 0..atom_count {
                if num_n_plus <= 0 || num_o_minus <= 0 {
                    break;
                }
                if atoms[i].charge != 1
                    || n_no_metal_num_bonds(Some(atoms), i as i32)? != 3
                    || n_no_metal_bonds_valence(Some(atoms), i as i32)? != 4
                    || ion_el_group(i32::from(atoms[i].el_number)) != 7
                {
                    continue;
                }
                let mut num_om = 0_usize;
                let mut ord_om = [0_i32; 3];
                let mut num_o = 0_i32;
                let mut num_o_other = 0_i32;
                for i1 in 0..i32::from(atoms[i].valence) {
                    let i1u =
                        usize::try_from(i1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    let n = usize::from(atoms[i].neighbor[i1u]);
                    if n_no_metal_num_bonds(Some(atoms), n as i32)? == 1
                        && num_of_H(atoms, n as i32)? == 0
                        && ion_el_group(i32::from(atoms[n].el_number)) == 8
                    {
                        if atoms[i].bond_type[i1u] == BOND_TYPE_SINGLE as u8
                            && atoms[n].charge == -1
                        {
                            let slot = ord_om
                                .get_mut(num_om)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            *slot = i1;
                            num_om += 1;
                        } else if atoms[n].bond_type[0] == BOND_TYPE_DOUBLE as u8
                            && atoms[n].charge == 0
                        {
                            num_o += 1;
                        } else {
                            num_o_other += 1;
                        }
                    }
                }
                if num_om > 0 && num_o > 0 && num_o_other == 0 {
                    let i1 = nFindOneOM(atoms, i as i32, &mut ord_om, num_om as i32)?;
                    if i1 >= 0 {
                        let i1 = i1 as usize;
                        let n = usize::from(atoms[i].neighbor[i1]);
                        let i2 = reciprocal(atoms, n, i)?;
                        atoms[i].bond_type[i1] = atoms[i].bond_type[i1].wrapping_add(1);
                        atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                        atoms[i].chem_bonds_valence = atoms[i].chem_bonds_valence.wrapping_add(1);
                        atoms[n].chem_bonds_valence = atoms[n].chem_bonds_valence.wrapping_add(1);
                        atoms[i].charge = atoms[i].charge.wrapping_sub(1);
                        atoms[n].charge = atoms[n].charge.wrapping_add(1);
                        atoms[i].radical = 0;
                        atoms[n].radical = 0;
                        num_changes += 1;
                        num_n_plus -= 1;
                        num_o_minus -= 1;
                        num_all -= 2;
                    }
                }
            }

            for i in 0..atom_count {
                if num_p_iv_plus <= 0 || num_o_minus <= 0 {
                    break;
                }
                if atoms[i].charge != 1
                    || n_no_metal_num_bonds(Some(atoms), i as i32)? != 4
                    || n_no_metal_bonds_valence(Some(atoms), i as i32)? != 4
                    || atoms[i].el_number == 7
                    || ion_el_group(i32::from(atoms[i].el_number)) != 7
                {
                    continue;
                }
                let mut num_om = 0_usize;
                let mut ord_om = [0_i32; 4];
                for i1 in 0..i32::from(atoms[i].valence) {
                    let i1u =
                        usize::try_from(i1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    let n = usize::from(atoms[i].neighbor[i1u]);
                    if n_no_metal_num_bonds(Some(atoms), n as i32)? == 1
                        && num_of_H(atoms, n as i32)? == 0
                        && ion_el_group(i32::from(atoms[n].el_number)) == 8
                        && atoms[i].bond_type[i1u] == BOND_TYPE_SINGLE as u8
                        && atoms[n].charge == -1
                    {
                        let slot = ord_om
                            .get_mut(num_om)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        *slot = i1;
                        num_om += 1;
                    }
                }
                if num_om > 0 {
                    let i1 = nFindOneOM(atoms, i as i32, &mut ord_om, num_om as i32)?;
                    if i1 >= 0 {
                        let i1 = i1 as usize;
                        let n = usize::from(atoms[i].neighbor[i1]);
                        let i2 = reciprocal(atoms, n, i)?;
                        atoms[i].bond_type[i1] = atoms[i].bond_type[i1].wrapping_add(1);
                        atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                        atoms[i].chem_bonds_valence = atoms[i].chem_bonds_valence.wrapping_add(1);
                        atoms[n].chem_bonds_valence = atoms[n].chem_bonds_valence.wrapping_add(1);
                        atoms[i].charge = atoms[i].charge.wrapping_sub(1);
                        atoms[n].charge = atoms[n].charge.wrapping_add(1);
                        atoms[i].radical = 0;
                        atoms[n].radical = 0;
                        num_changes += 1;
                        num_n_plus -= 1;
                        num_o_minus -= 1;
                        num_p_iv_plus -= 1;
                        num_all -= 2;
                    }
                }
            }
        }

        if (2..=9).contains(&pair_type) {
            'terminal: for i in 0..atom_count {
                if num_all <= 0 {
                    break;
                }
                let neutral_terminal_o = atoms[i].charge == 0
                    && n_no_metal_num_bonds(Some(atoms), i as i32)? == 1
                    && n_no_metal_bonds_valence(Some(atoms), i as i32)? == 2
                    && num_of_H(atoms, i as i32)? == 0
                    && ion_el_group(i32::from(atoms[i].el_number)) == 8;
                if neutral_terminal_o {
                    let i1 = n_no_metal_neigh_index(Some(atoms), i as i32)?;
                    if i1 >= 0 && atoms[i].bond_type[i1 as usize] <= BOND_TYPE_TRIPLE as u8 {
                        let n = usize::from(atoms[i].neighbor[i1 as usize]);
                        if pair_type == 2
                            && num_c_ii > 0
                            && atoms[n].charge == 0
                            && n_no_metal_num_bonds(Some(atoms), n as i32)? == 2
                            && n_no_metal_bonds_valence(Some(atoms), n as i32)? == 3
                            && num_of_H(atoms, n as i32)? == 0
                            && ion_el_group(i32::from(atoms[n].el_number)) == 7
                        {
                            let i2 = n_no_metal_other_neigh_index(Some(atoms), n as i32, i as i32)?;
                            if i2 >= 0 && atoms[n].bond_type[i2 as usize] <= BOND_TYPE_TRIPLE as u8
                            {
                                let n2 = usize::from(atoms[n].neighbor[i2 as usize]);
                                if atoms[n2].charge == 0
                                    && atoms[n2].valence == 2
                                    && atoms[n2].chem_bonds_valence == 2
                                    && num_of_H(atoms, n2 as i32)? == 0
                                    && ion_el_group(i32::from(atoms[n2].el_number)) == 6
                                {
                                    let i3 = usize::from(atoms[n2].neighbor[0] != n as u16);
                                    atoms[n].chem_bonds_valence = 5;
                                    atoms[n2].chem_bonds_valence = 4;
                                    atoms[n].bond_type[i2 as usize] = BOND_TYPE_TRIPLE as u8;
                                    atoms[n2].bond_type[i3] = BOND_TYPE_TRIPLE as u8;
                                    atoms[n2].radical = 0;
                                    num_changes += 1;
                                    num_c_ii -= 1;
                                    num_all -= 1;
                                    continue 'terminal;
                                }
                            }
                        }
                        if pair_type == 3
                            && num_o_plus > 0
                            && num_c_minus > 0
                            && atoms[n].charge == 1
                            && n_no_metal_num_bonds(Some(atoms), n as i32)? == 2
                            && n_no_metal_bonds_valence(Some(atoms), n as i32)? == 3
                            && num_of_H(atoms, n as i32)? == 0
                            && ion_el_group(i32::from(atoms[n].el_number)) == 8
                        {
                            let i2 = n_no_metal_other_neigh_index(Some(atoms), n as i32, i as i32)?;
                            if i2 >= 0 && atoms[n].bond_type[i2 as usize] <= BOND_TYPE_TRIPLE as u8
                            {
                                let n2 = usize::from(atoms[n].neighbor[i2 as usize]);
                                if atoms[n2].charge == -1
                                    && n_no_metal_num_bonds(Some(atoms), n2 as i32)? <= 3
                                    && n_no_metal_bonds_valence(Some(atoms), n2 as i32)?
                                        + numh(&atoms[n2])
                                        == 3
                                    && ion_el_group(i32::from(atoms[n2].el_number)) == 6
                                {
                                    let i3 = usize::from(atoms[n2].neighbor[0] != n as u16);
                                    atoms[n].charge = atoms[n].charge.wrapping_sub(1);
                                    atoms[n2].charge = atoms[n2].charge.wrapping_add(1);
                                    atoms[n].chem_bonds_valence =
                                        atoms[n].chem_bonds_valence.wrapping_add(1);
                                    atoms[n2].chem_bonds_valence =
                                        atoms[n2].chem_bonds_valence.wrapping_add(1);
                                    atoms[n].bond_type[i2 as usize] = BOND_TYPE_DOUBLE as u8;
                                    atoms[n2].bond_type[i3] = BOND_TYPE_DOUBLE as u8;
                                    num_changes += 1;
                                    num_o_plus -= 1;
                                    num_c_minus -= 1;
                                    num_all -= 2;
                                    continue 'terminal;
                                }
                            }
                        }
                    }
                    continue 'terminal;
                }

                let terminal_o_minus = atoms[i].charge == -1
                    && num_o_minus + num_n_minus > 0
                    && num_n_plus + num_o_plus + num_c_plus > 0
                    && n_no_metal_num_bonds(Some(atoms), i as i32)? == 1
                    && n_no_metal_bonds_valence(Some(atoms), i as i32)? == 1
                    && num_of_H(atoms, i as i32)? == 0
                    && ion_el_group(i32::from(atoms[i].el_number)) == 8;
                if terminal_o_minus {
                    let i1 = n_no_metal_neigh_index(Some(atoms), i as i32)?;
                    if i1 >= 0 && atoms[i].bond_type[i1 as usize] <= BOND_TYPE_TRIPLE as u8 {
                        let i1 = i1 as usize;
                        let n = usize::from(atoms[i].neighbor[i1]);
                        let direct_target = (pair_type == 4
                            && num_o_minus > 0
                            && num_n_plus > 0
                            && atoms[n].charge == 1
                            && n_no_metal_num_bonds(Some(atoms), n as i32)? <= 3
                            && n_no_metal_bonds_valence(Some(atoms), n as i32)? == 4
                            && num_of_H(atoms, n as i32)? == 0
                            && ion_el_group(i32::from(atoms[n].el_number)) == 7)
                            || (pair_type == 5
                                && num_o_minus > 0
                                && num_o_plus > 0
                                && atoms[n].charge == 1
                                && n_no_metal_num_bonds(Some(atoms), n as i32)? <= 3
                                && n_no_metal_bonds_valence(Some(atoms), n as i32)? == 3
                                && num_of_H(atoms, n as i32)? == 0
                                && ion_el_group(i32::from(atoms[n].el_number)) == 8);
                        if direct_target {
                            let i2 = reciprocal(atoms, n, i)?;
                            atoms[i].charge = atoms[i].charge.wrapping_add(1);
                            atoms[n].charge = atoms[n].charge.wrapping_sub(1);
                            atoms[i].chem_bonds_valence =
                                atoms[i].chem_bonds_valence.wrapping_add(1);
                            atoms[n].chem_bonds_valence =
                                atoms[n].chem_bonds_valence.wrapping_add(1);
                            atoms[i].bond_type[i1] = atoms[i].bond_type[i1].wrapping_add(1);
                            atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                            num_changes += 1;
                            num_o_minus -= 1;
                            if pair_type == 4 {
                                num_n_plus -= 1;
                            } else {
                                num_o_plus -= 1;
                            }
                            num_all -= 2;
                            continue 'terminal;
                        }
                        if pair_type == 6
                            && num_o_minus > 0
                            && num_c_plus > 0
                            && atoms[n].charge == 0
                            && n_no_metal_num_bonds(Some(atoms), n as i32)? == 2
                            && n_no_metal_bonds_valence(Some(atoms), n as i32)? == 2
                            && num_of_H(atoms, n as i32)? == 0
                            && ion_el_group(i32::from(atoms[n].el_number)) == 8
                        {
                            let i2 = n_no_metal_other_neigh_index(Some(atoms), n as i32, i as i32)?;
                            if i2 >= 0 && atoms[n].bond_type[i2 as usize] <= BOND_TYPE_TRIPLE as u8
                            {
                                let n2 = usize::from(atoms[n].neighbor[i2 as usize]);
                                if atoms[n2].charge == 1
                                    && n_no_metal_num_bonds(Some(atoms), n2 as i32)? <= 3
                                    && n_no_metal_bonds_valence(Some(atoms), n2 as i32)?
                                        + numh(&atoms[n2])
                                        == 3
                                    && ion_el_group(i32::from(atoms[n2].el_number)) == 6
                                {
                                    let i3 = reciprocal(atoms, n2, n)?;
                                    let i4 = reciprocal(atoms, n, i)?;
                                    atoms[i].charge = atoms[i].charge.wrapping_add(1);
                                    atoms[n2].charge = atoms[n2].charge.wrapping_sub(1);
                                    atoms[i].chem_bonds_valence =
                                        atoms[i].chem_bonds_valence.wrapping_add(1);
                                    atoms[n].chem_bonds_valence =
                                        atoms[n].chem_bonds_valence.wrapping_add(2);
                                    atoms[n2].chem_bonds_valence =
                                        atoms[n2].chem_bonds_valence.wrapping_add(1);
                                    atoms[i].bond_type[i1] = BOND_TYPE_DOUBLE as u8;
                                    atoms[n].bond_type[i4] = BOND_TYPE_DOUBLE as u8;
                                    atoms[n].bond_type[i2 as usize] = BOND_TYPE_DOUBLE as u8;
                                    atoms[n2].bond_type[i3] = BOND_TYPE_DOUBLE as u8;
                                    num_changes += 1;
                                    num_o_minus -= 1;
                                    num_c_plus -= 1;
                                    num_all -= 2;
                                    continue 'terminal;
                                }
                            }
                        }
                    }
                    continue 'terminal;
                }

                let terminal_n_minus = atoms[i].charge == -1
                    && num_n_minus > 0
                    && num_n_plus + num_o_plus + num_c_plus > 0
                    && n_no_metal_num_bonds(Some(atoms), i as i32)? == 1
                    && n_no_metal_bonds_valence(Some(atoms), i as i32)? + numh(&atoms[i]) == 2
                    && ion_el_group(i32::from(atoms[i].el_number)) == 7;
                if terminal_n_minus {
                    let i1_check = n_no_metal_neigh_index(Some(atoms), i as i32)?;
                    if i1_check >= 0
                        && atoms[i].bond_type[i1_check as usize] <= BOND_TYPE_TRIPLE as u8
                    {
                        let i1 = 0_usize;
                        let n = usize::from(atoms[i].neighbor[i1]);
                        let target = (pair_type == 7
                            && num_n_plus > 0
                            && atoms[n].charge == 1
                            && n_no_metal_num_bonds(Some(atoms), n as i32)? <= 3
                            && n_no_metal_bonds_valence(Some(atoms), n as i32)? == 4
                            && num_of_H(atoms, n as i32)? == 0
                            && ion_el_group(i32::from(atoms[n].el_number)) == 7)
                            || (pair_type == 8
                                && num_o_plus > 0
                                && atoms[n].charge == 1
                                && n_no_metal_num_bonds(Some(atoms), n as i32)? == 2
                                && n_no_metal_bonds_valence(Some(atoms), n as i32)? == 3
                                && num_of_H(atoms, n as i32)? == 0
                                && ion_el_group(i32::from(atoms[n].el_number)) == 8)
                            || (pair_type == 9
                                && num_c_plus > 0
                                && atoms[n].charge == 1
                                && atoms[n].valence == 2
                                && atoms[n].chem_bonds_valence == 3
                                && num_of_H(atoms, n as i32)? == 0
                                && ion_el_group(i32::from(atoms[n].el_number)) == 6);
                        if target {
                            let i2 = reciprocal(atoms, n, i)?;
                            atoms[i].charge = atoms[i].charge.wrapping_add(1);
                            atoms[n].charge = atoms[n].charge.wrapping_sub(1);
                            atoms[i].chem_bonds_valence =
                                atoms[i].chem_bonds_valence.wrapping_add(1);
                            atoms[n].chem_bonds_valence =
                                atoms[n].chem_bonds_valence.wrapping_add(1);
                            atoms[i].bond_type[i1] = atoms[i].bond_type[i1].wrapping_add(1);
                            atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                            num_changes += 1;
                            num_n_minus -= 1;
                            if pair_type == 7 {
                                num_n_plus -= 1;
                            } else if pair_type == 8 {
                                num_o_plus -= 1;
                            } else {
                                num_c_plus -= 1;
                            }
                            num_all -= 2;
                            continue 'terminal;
                        }
                    }
                }
            }
        }

        if (10..=15).contains(&pair_type) {
            'nonterminal: for i in 0..atom_count {
                if num_all <= 0 {
                    break;
                }
                if atoms[i].charge == 1
                    && num_n_plus + num_o_plus > 0
                    && num_c_minus + num_n_minus > 0
                    && n_no_metal_num_bonds(Some(atoms), i as i32)? <= 4
                    && n_no_metal_bonds_valence(Some(atoms), i as i32)? == 4
                    && num_of_H(atoms, i as i32)? == 0
                    && ion_el_group(i32::from(atoms[i].el_number)) == 7
                {
                    let expected_bond = if pair_type == 11 { 2_u8 } else { 1_u8 };
                    if (pair_type == 10 || pair_type == 11) && num_n_plus > 0 && num_c_minus > 0 {
                        let mut num_neigh = 0_i32;
                        let mut position = 0_usize;
                        for i1 in 0..i32::from(atoms[i].valence) {
                            let i1 = i1 as usize;
                            let n = usize::from(atoms[i].neighbor[i1]);
                            if atoms[n].charge == -1
                                && atoms[n].valence <= 3
                                && i32::from(atoms[n].chem_bonds_valence) + numh(&atoms[n]) == 3
                                && atoms[i].bond_type[i1] == expected_bond
                                && ion_el_group(i32::from(atoms[n].el_number)) == 6
                            {
                                num_neigh += 1;
                                position = i1;
                            }
                        }
                        if num_neigh == 1 {
                            let n = usize::from(atoms[i].neighbor[position]);
                            if (pair_type != 10
                                || atoms[i].bond_type[position] <= BOND_TYPE_TRIPLE as u8)
                                && has_other_ion_neigh(atoms, i as i32, n as i32)? == 0
                                && has_other_ion_neigh(atoms, n as i32, i as i32)? == 0
                            {
                                let i2 = reciprocal(atoms, n, i)?;
                                atoms[i].charge = atoms[i].charge.wrapping_sub(1);
                                atoms[n].charge = atoms[n].charge.wrapping_add(1);
                                atoms[i].chem_bonds_valence =
                                    atoms[i].chem_bonds_valence.wrapping_add(1);
                                atoms[n].chem_bonds_valence =
                                    atoms[n].chem_bonds_valence.wrapping_add(1);
                                atoms[i].bond_type[position] =
                                    atoms[i].bond_type[position].wrapping_add(1);
                                atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                                num_changes += 1;
                                num_c_minus -= 1;
                                num_n_plus -= 1;
                                num_all -= 2;
                                continue 'nonterminal;
                            }
                        }
                    }
                    if pair_type == 12 && num_n_plus > 0 && num_n_minus > 0 {
                        let mut num_neigh = 0_i32;
                        let mut position = 0_usize;
                        for i1 in 0..i32::from(atoms[i].valence) {
                            let i1 = i1 as usize;
                            let n = usize::from(atoms[i].neighbor[i1]);
                            if atoms[n].charge == -1
                                && n_no_metal_num_bonds(Some(atoms), n as i32)? <= 2
                                && n_no_metal_bonds_valence(Some(atoms), n as i32)?
                                    + numh(&atoms[n])
                                    == 2
                                && atoms[i].bond_type[i1] == BOND_TYPE_SINGLE as u8
                                && ion_el_group(i32::from(atoms[n].el_number)) == 7
                            {
                                num_neigh += 1;
                                position = i1;
                            }
                        }
                        if num_neigh == 1 {
                            let n = usize::from(atoms[i].neighbor[position]);
                            if has_other_ion_neigh(atoms, i as i32, n as i32)? == 0
                                && has_other_ion_neigh(atoms, n as i32, i as i32)? == 0
                            {
                                let i2 = reciprocal(atoms, n, i)?;
                                atoms[i].charge = atoms[i].charge.wrapping_sub(1);
                                atoms[n].charge = atoms[n].charge.wrapping_add(1);
                                atoms[i].chem_bonds_valence =
                                    atoms[i].chem_bonds_valence.wrapping_add(1);
                                atoms[n].chem_bonds_valence =
                                    atoms[n].chem_bonds_valence.wrapping_add(1);
                                atoms[i].bond_type[position] =
                                    atoms[i].bond_type[position].wrapping_add(1);
                                atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                                num_changes += 1;
                                num_n_minus -= 1;
                                num_n_plus -= 1;
                                num_all -= 2;
                                continue 'nonterminal;
                            }
                        }
                    }
                    continue 'nonterminal;
                }

                if atoms[i].charge == 1
                    && num_o_plus > 0
                    && num_c_minus + num_n_minus > 0
                    && n_no_metal_num_bonds(Some(atoms), i as i32)? <= 3
                    && n_no_metal_bonds_valence(Some(atoms), i as i32)? == 3
                    && num_of_H(atoms, i as i32)? == 0
                    && ion_el_group(i32::from(atoms[i].el_number)) == 8
                {
                    let (target_group, max_bonds, target_valence, expected_bond) =
                        if pair_type == 15 {
                            (7_u8, 2_i32, 2_i32, BOND_TYPE_SINGLE as u8)
                        } else {
                            (
                                6_u8,
                                3_i32,
                                3_i32,
                                if pair_type == 14 {
                                    BOND_TYPE_DOUBLE as u8
                                } else {
                                    BOND_TYPE_SINGLE as u8
                                },
                            )
                        };
                    if (pair_type == 13 && num_c_minus > 0)
                        || (pair_type == 14 && num_c_minus > 0)
                        || (pair_type == 15 && num_n_minus > 0)
                    {
                        let mut num_neigh = 0_i32;
                        let mut position = 0_usize;
                        for i1 in 0..i32::from(atoms[i].valence) {
                            let i1 = i1 as usize;
                            let n = usize::from(atoms[i].neighbor[i1]);
                            let bonds = if pair_type == 15 {
                                n_no_metal_num_bonds(Some(atoms), n as i32)?
                            } else {
                                i32::from(atoms[n].valence)
                            };
                            let valence = if pair_type == 15 {
                                n_no_metal_bonds_valence(Some(atoms), n as i32)? + numh(&atoms[n])
                            } else {
                                i32::from(atoms[n].chem_bonds_valence) + numh(&atoms[n])
                            };
                            if atoms[n].charge == -1
                                && bonds <= max_bonds
                                && valence == target_valence
                                && atoms[i].bond_type[i1] == expected_bond
                                && ion_el_group(i32::from(atoms[n].el_number)) == target_group
                            {
                                num_neigh += 1;
                                position = i1;
                            }
                        }
                        if num_neigh == 1 {
                            let n = usize::from(atoms[i].neighbor[position]);
                            if has_other_ion_neigh(atoms, i as i32, n as i32)? == 0
                                && has_other_ion_neigh(atoms, n as i32, i as i32)? == 0
                            {
                                let i2 = reciprocal(atoms, n, i)?;
                                atoms[i].charge = atoms[i].charge.wrapping_sub(1);
                                atoms[n].charge = atoms[n].charge.wrapping_add(1);
                                atoms[i].chem_bonds_valence =
                                    atoms[i].chem_bonds_valence.wrapping_add(1);
                                atoms[n].chem_bonds_valence =
                                    atoms[n].chem_bonds_valence.wrapping_add(1);
                                atoms[i].bond_type[position] =
                                    atoms[i].bond_type[position].wrapping_add(1);
                                atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                                num_changes += 1;
                                if pair_type == 15 {
                                    num_n_minus -= 1;
                                } else {
                                    num_c_minus -= 1;
                                }
                                num_o_plus -= 1;
                                num_all -= 2;
                                continue 'nonterminal;
                            }
                        }
                    }
                }
            }
        }

        if pair_type == 16 && num_c_plus > 0 && num_n_minus > 0 {
            for i in 0..atom_count {
                if atoms[i].charge != 0
                    || n_no_metal_num_bonds(Some(atoms), i as i32)? != 2
                    || n_no_metal_bonds_valence(Some(atoms), i as i32)? != 2
                    || num_of_H(atoms, i as i32)? != 0
                    || ion_el_group(i32::from(atoms[i].el_number)) != 8
                {
                    continue;
                }
                let j0 = n_no_metal_neigh_index(Some(atoms), i as i32)?;
                if j0 < 0 {
                    continue;
                }
                let m0 = usize::from(atoms[i].neighbor[j0 as usize]);
                if atoms[m0].charge == 0 {
                    continue;
                }
                let j1 = n_no_metal_other_neigh_index(Some(atoms), i as i32, m0 as i32)?;
                if j1 < 0 {
                    continue;
                }
                let m1 = usize::from(atoms[i].neighbor[j1 as usize]);
                if i32::from(atoms[m0].charge) + i32::from(atoms[m1].charge) != 0
                    || n_no_metal_bonds_valence(Some(atoms), m0 as i32)?
                        + n_no_metal_bonds_valence(Some(atoms), m1 as i32)?
                        > 5
                {
                    continue;
                }
                let members = [m0, m1];
                let positions = [j0 as usize, j1 as usize];
                let mut n_n = None;
                let mut n_c = None;
                let mut i_c = 0_usize;
                for k in 0..2 {
                    let n = members[k];
                    if atoms[n].charge == -1
                        && n_no_metal_num_bonds(Some(atoms), n as i32)? + numh(&atoms[n]) == 2
                        && ion_el_group(i32::from(atoms[n].el_number)) == 7
                    {
                        n_n = Some(n);
                    } else if atoms[n].charge == 1
                        && i32::from(atoms[n].chem_bonds_valence) + numh(&atoms[n]) == 3
                        && ion_el_group(i32::from(atoms[n].el_number)) == 6
                    {
                        n_c = Some(n);
                        i_c = k;
                    }
                }
                let (Some(n_c), Some(n_n)) = (n_c, n_n) else {
                    continue;
                };
                if has_other_ion_in_sphere_2(atoms, n_c as i32, n_n as i32)? != 0
                    || has_other_ion_in_sphere_2(atoms, n_n as i32, n_c as i32)? != 0
                {
                    continue;
                }
                for k in 0..2 {
                    let n = if k == 1 { n_c } else { n_n };
                    let i1 = if k == 1 {
                        positions[i_c]
                    } else {
                        positions[1 - i_c]
                    };
                    let i2 = reciprocal(atoms, n, i)?;
                    atoms[i].bond_type[i1] = atoms[i].bond_type[i1].wrapping_add(1);
                    atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                    atoms[i].chem_bonds_valence = atoms[i].chem_bonds_valence.wrapping_add(1);
                    atoms[n].chem_bonds_valence = atoms[n].chem_bonds_valence.wrapping_add(1);
                    atoms[n].charge = if k == 1 {
                        atoms[n].charge.wrapping_sub(1)
                    } else {
                        atoms[n].charge.wrapping_add(1)
                    };
                }
                num_changes += 1;
                num_n_minus -= 1;
                num_c_plus -= 1;
                num_all -= 2;
            }
        }

        if pair_type == 17 && num_c_plus > 0 && num_c_minus > 0 {
            for i in 0..atom_count {
                if atoms[i].charge != 0
                    || n_no_metal_num_bonds(Some(atoms), i as i32)? != 3
                    || n_no_metal_bonds_valence(Some(atoms), i as i32)? != 3
                    || num_of_H(atoms, i as i32)? != 0
                    || ion_el_group(i32::from(atoms[i].el_number)) != 7
                {
                    continue;
                }
                let j0 = n_no_metal_neigh_index(Some(atoms), i as i32)?;
                if j0 < 0 {
                    continue;
                }
                let m0 = usize::from(atoms[i].neighbor[j0 as usize]);
                let j1 = n_no_metal_other_neigh_index(Some(atoms), i as i32, m0 as i32)?;
                if j1 < 0 {
                    continue;
                }
                let m1 = usize::from(atoms[i].neighbor[j1 as usize]);
                let j2 =
                    n_no_metal_other_neigh_index2(Some(atoms), i as i32, m0 as i32, m1 as i32)?;
                if j2 < 0 {
                    continue;
                }
                let m2 = usize::from(atoms[i].neighbor[j2 as usize]);
                let members = [m0, m1, m2];
                let positions = [j0 as usize, j1 as usize, j2 as usize];
                let charges = members.map(|n| i32::from(atoms[n].charge));
                if charges.iter().filter(|&&charge| charge == 0).count() != 1
                    || charges.iter().sum::<i32>() != 0
                    || members
                        .iter()
                        .filter(|&&n| {
                            atoms[n].charge != 0
                                && i32::from(atoms[n].chem_bonds_valence) + numh(&atoms[n]) == 3
                        })
                        .count()
                        != 2
                {
                    continue;
                }
                let mut n_cp = None;
                let mut n_cm = None;
                let mut i_cp = 0_usize;
                let mut i_cm = 0_usize;
                for k in 0..3 {
                    let n = members[k];
                    if charges[k] == -1 && ion_el_group(i32::from(atoms[n].el_number)) == 6 {
                        n_cm = Some(n);
                        i_cm = k;
                    } else if charges[k] == 1 && ion_el_group(i32::from(atoms[n].el_number)) == 6 {
                        n_cp = Some(n);
                        i_cp = k;
                    }
                }
                let (Some(n_cp), Some(n_cm)) = (n_cp, n_cm) else {
                    continue;
                };
                if has_other_ion_in_sphere_2(atoms, n_cp as i32, n_cm as i32)? != 0
                    || has_other_ion_in_sphere_2(atoms, n_cm as i32, n_cp as i32)? != 0
                {
                    continue;
                }
                for k in 0..2 {
                    let (n, i1) = if k == 1 {
                        (n_cp, positions[i_cp])
                    } else {
                        (n_cm, positions[i_cm])
                    };
                    let i2 = reciprocal(atoms, n, i)?;
                    atoms[i].bond_type[i1] = atoms[i].bond_type[i1].wrapping_add(1);
                    atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(1);
                    atoms[i].chem_bonds_valence = atoms[i].chem_bonds_valence.wrapping_add(1);
                    atoms[n].chem_bonds_valence = atoms[n].chem_bonds_valence.wrapping_add(1);
                    atoms[n].charge = if k == 1 {
                        atoms[n].charge.wrapping_sub(1)
                    } else {
                        atoms[n].charge.wrapping_add(1)
                    };
                }
                num_changes += 1;
                num_c_minus -= 1;
                num_c_plus -= 1;
                num_all -= 2;
            }
        }

        if pair_type == 18 && ((num_c_plus > 0 && num_c_minus > 0) || num_c_ii > 0) {
            for i in 0..atom_count {
                if atoms[i].charge != 0
                    || n_no_metal_num_bonds(Some(atoms), i as i32)? != 2
                    || n_no_metal_bonds_valence(Some(atoms), i as i32)? != 3
                    || num_of_H(atoms, i as i32)? != 0
                    || ion_el_group(i32::from(atoms[i].el_number)) != 7
                {
                    continue;
                }
                let j0 = n_no_metal_neigh_index(Some(atoms), i as i32)?;
                if j0 < 0 {
                    continue;
                }
                let m0 = usize::from(atoms[i].neighbor[j0 as usize]);
                let j1 = n_no_metal_other_neigh_index(Some(atoms), i as i32, m0 as i32)?;
                if j1 < 0 {
                    continue;
                }
                let m1 = usize::from(atoms[i].neighbor[j1 as usize]);
                if i32::from(atoms[m0].charge) + i32::from(atoms[m1].charge) != 0
                    || ion_el_group(i32::from(atoms[m0].el_number)) != 6
                    || ion_el_group(i32::from(atoms[m1].el_number)) != 6
                {
                    continue;
                }
                let members = [m0, m1];
                let positions = [j0 as usize, j1 as usize];
                let values =
                    members.map(|n| i32::from(atoms[n].chem_bonds_valence) + numh(&atoms[n]));
                if values[0] + values[1] != 6 || (values[0] - values[1]).abs() > 2 {
                    continue;
                }
                let mut n_cp = None;
                let mut n_cm = None;
                let mut i_cp = 0_usize;
                let mut i_cm = 0_usize;
                for k in 0..2 {
                    if values[k] == 4
                        || (values[k] == 3
                            && atoms[i].bond_type[positions[k]] == BOND_TYPE_SINGLE as u8)
                    {
                        n_cm = Some(members[k]);
                        i_cm = k;
                    } else if values[k] == 2
                        || (values[k] == 3
                            && atoms[i].bond_type[positions[k]] == BOND_TYPE_DOUBLE as u8)
                    {
                        n_cp = Some(members[k]);
                        i_cp = k;
                    }
                }
                let (Some(n_cp), Some(n_cm)) = (n_cp, n_cm) else {
                    continue;
                };
                if i32::from(atoms[n_cp].valence) + numh(&atoms[n_cp]) != 2 {
                    continue;
                }
                if values[i_cp] == 2 || atoms[n_cp].charge == 0 {
                    if atoms[n_cp].valence == 2 {
                        let opposite = usize::from(
                            atoms[n_cp].neighbor[usize::from(atoms[n_cp].neighbor[0] == i as u16)],
                        );
                        if ion_el_group(i32::from(atoms[opposite].el_number)) == 7 {
                            continue;
                        }
                    }
                } else if atoms[n_cp].charge != 0 {
                    if has_other_ion_in_sphere_2(atoms, n_cp as i32, n_cm as i32)? != 0
                        || has_other_ion_in_sphere_2(atoms, n_cm as i32, n_cp as i32)? != 0
                    {
                        continue;
                    }
                } else {
                    continue;
                }
                if atoms[n_cp].charge != 0 {
                    num_c_minus -= 1;
                    num_c_plus -= 1;
                    num_all -= 2;
                } else {
                    num_c_ii -= 1;
                    num_all -= 1;
                }
                for k in 0..2 {
                    let (n, member_index) = if k == 1 { (n_cp, i_cp) } else { (n_cm, i_cm) };
                    let i1 = positions[member_index];
                    if values[member_index] < 4 {
                        let delta = 4 - values[member_index];
                        let i2 = reciprocal(atoms, n, i)?;
                        atoms[i].bond_type[i1] = atoms[i].bond_type[i1].wrapping_add(delta as u8);
                        atoms[n].bond_type[i2] = atoms[n].bond_type[i2].wrapping_add(delta as u8);
                        atoms[i].chem_bonds_valence =
                            atoms[i].chem_bonds_valence.wrapping_add(delta as i8);
                        atoms[n].chem_bonds_valence =
                            atoms[n].chem_bonds_valence.wrapping_add(delta as i8);
                        atoms[n].charge = 0;
                        atoms[n].radical = 0;
                    }
                }
                atoms[i].charge = 0;
                atoms[i].radical = 0;
                num_changes += 1;
            }
        }
    }
    Ok(num_changes)
}

#[allow(non_snake_case)]
pub(crate) fn bIsAmmoniumSalt(
    atoms: &[inp_ATOM],
    atom_index: i32,
    oxygen_index: &mut i32,
    neighbor_position: &mut i32,
    num_explicit_h: &mut [S_CHAR; NUM_H_ISOTOPES as usize + 1],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2300 bIsAmmoniumSalt
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: bIsAmmoniumSalt
    // INCHI✔️✔️: int bIsAmmoniumSalt( inp_ATOM *at,
    // INCHI✔️✔️:                      int i,
    // INCHI✔️✔️:                      int *piO,
    // INCHI✔️✔️:                      int *pk,
    // INCHI✔️✔️:                      S_CHAR *num_explicit_H )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* NH4(+charge)-O(-charge)-C -> NH3 + HO-C; any charge including 0, */
    // INCHI✔️✔️:     /* any C except charged or radical F, Cl, Br, I                     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int num_H, num_non_iso_H, num_impl_iso_H, bDisconnect = 1;
    // INCHI✔️✔️:     int j, val, neigh, iO = -1, iC, k = -1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at[i].el_number != EL_NUMBER_N)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* check for NH4-O-C... -> NH3 + HO-C... */
    // INCHI✔️✔️:     val = at[i].valence;
    // INCHI✔️✔️:     num_impl_iso_H = NUM_ISO_H( at, i );
    // INCHI✔️✔️:     num_non_iso_H = at[i].num_H;
    // INCHI✔️✔️:     num_H = num_non_iso_H + num_impl_iso_H;
    // INCHI✔️✔️:     if (val + num_H == 5)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int num_O = 0;
    // INCHI✔️✔️:         memset( num_explicit_H, 0, ( NUM_H_ISOTOPES + 1 ) * sizeof( num_explicit_H[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         for (j = 0; j < val; j++)
    // INCHI✔️✔️:         { /* looking for O: H4N-O-C... */
    // INCHI✔️✔️:             neigh = at[i].neighbor[j];
    // INCHI✔️✔️:             if (at[neigh].num_H ||
    // INCHI✔️✔️:                  (at[neigh].charge && ( at[neigh].el_number != EL_NUMBER_O || at[neigh].charge + at[i].charge )) ||
    // INCHI✔️✔️:                  (at[neigh].radical && at[neigh].radical != RADICAL_SINGLET)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 bDisconnect = 0;
    // INCHI✔️✔️:                 break; /* reject */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (at[neigh].el_number == EL_NUMBER_H && at[neigh].valence == 1 &&
    // INCHI✔️✔️:                  !at[neigh].charge && !at[neigh].radical)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_H++; /* at this point at[].num_H does not include explicit H count */
    // INCHI✔️✔️:                 num_non_iso_H += ( 0 == at[neigh].iso_atw_diff );
    // INCHI✔️✔️:                 num_explicit_H[at[neigh].iso_atw_diff] ++;  /* explicit H on N */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else if (at[neigh].el_number == EL_NUMBER_O && at[neigh].valence == 2 && !num_O)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_O++; /* found O: N-O- */
    // INCHI✔️✔️:                 iO = neigh;
    // INCHI✔️✔️:                 k = j;
    // INCHI✔️✔️:                 iC = at[iO].neighbor[at[iO].neighbor[0] == i];
    // INCHI✔️✔️:                 if (at[iC].el_number != EL_NUMBER_C || /*
    // INCHI✔️✔️:                                                        at[iC].num_H ||
    // INCHI✔️✔️:                                                        at[iC].chem_bonds_valence != 4 || */
    // INCHI✔️✔️:                      at[iC].charge ||
    // INCHI✔️✔️:                      (at[iC].radical && at[iC].radical != RADICAL_SINGLET) /*||
    // INCHI✔️✔️:                                                                          at[iC].valence == at[iC].chem_bonds_valence*/) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     bDisconnect = 0;
    // INCHI✔️✔️:                     break; /* reject */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else if (( at[neigh].el_number == EL_NUMBER_F ||
    // INCHI✔️✔️:                        at[neigh].el_number == EL_NUMBER_CL ||
    // INCHI✔️✔️:                        at[neigh].el_number == EL_NUMBER_BR ||
    // INCHI✔️✔️:                        at[neigh].el_number == EL_NUMBER_I ) &&
    // INCHI✔️✔️:                       at[neigh].valence == 1 && at[neigh].chem_bonds_valence == 1 &&
    // INCHI✔️✔️:                       !at[neigh].charge && !NUMH( at, neigh ) && !num_O)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_O++; /* found O: N-O- */
    // INCHI✔️✔️:                 iO = neigh;
    // INCHI✔️✔️:                 k = j;
    // INCHI✔️✔️:                 /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 bDisconnect = 0;
    // INCHI✔️✔️:                 break;  /* reject */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (bDisconnect && ( num_O != 1 || num_H != 4 ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bDisconnect = 0; /* reject */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         bDisconnect = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (bDisconnect)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *piO = iO;
    // INCHI✔️✔️:         *pk = k;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return bDisconnect;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: bIsAmmoniumSalt
    // END INCHI C FUNCTION: bIsAmmoniumSalt
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsAmmoniumSalt
    // INCHI✔️✔️: #define EL_NUMBER_H  ((U_CHAR)1)
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_N  ((U_CHAR)7)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define EL_NUMBER_F  ((U_CHAR)9)
    // INCHI✔️✔️: #define EL_NUMBER_CL ((U_CHAR)17)
    // INCHI✔️✔️: #define EL_NUMBER_BR ((U_CHAR)35)
    // INCHI✔️✔️: #define EL_NUMBER_I  ((U_CHAR)53)
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // INCHI✔️✔️: #define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsAmmoniumSalt

    let i = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.el_number != 7 {
        return Ok(0);
    }

    let valence = i32::from(atom.valence);
    let mut num_h = i32::from(atom.num_H) + atom.num_iso_H.into_iter().map(i32::from).sum::<i32>();
    if valence + num_h != 5 {
        return Ok(0);
    }

    num_explicit_h.fill(0);
    let mut num_oxygen = 0_i32;
    let mut found_oxygen = -1_i32;
    let mut found_position = -1_i32;
    let mut disconnect = true;
    for j in 0..valence.max(0) as usize {
        let neighbor_index = usize::from(atom.neighbor[j]);
        let neighbor = atoms
            .get(neighbor_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.num_H != 0
            || (neighbor.charge != 0
                && (neighbor.el_number != 8
                    || i32::from(neighbor.charge) + i32::from(atom.charge) != 0))
            || (neighbor.radical != 0 && neighbor.radical != RADICAL_SINGLET as i8)
        {
            disconnect = false;
            break;
        }
        if neighbor.el_number == 1
            && neighbor.valence == 1
            && neighbor.charge == 0
            && neighbor.radical == 0
        {
            num_h += 1;
            let isotope = usize::try_from(neighbor.iso_atw_diff)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let explicit = num_explicit_h
                .get_mut(isotope)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *explicit = explicit.wrapping_add(1);
        } else if neighbor.el_number == 8 && neighbor.valence == 2 && num_oxygen == 0 {
            num_oxygen += 1;
            found_oxygen = neighbor_index as i32;
            found_position = j as i32;
            let opposite_position = usize::from(neighbor.neighbor[0] == i as AT_NUMB);
            let carbon_index = usize::from(neighbor.neighbor[opposite_position]);
            let carbon = atoms
                .get(carbon_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if carbon.el_number != 6
                || carbon.charge != 0
                || (carbon.radical != 0 && carbon.radical != RADICAL_SINGLET as i8)
            {
                disconnect = false;
                break;
            }
        } else if matches!(neighbor.el_number, 9 | 17 | 35 | 53)
            && neighbor.valence == 1
            && neighbor.chem_bonds_valence == 1
            && neighbor.charge == 0
            && neighbor.num_H == 0
            && neighbor.num_iso_H.into_iter().all(|count| count == 0)
            && num_oxygen == 0
        {
            num_oxygen += 1;
            found_oxygen = neighbor_index as i32;
            found_position = j as i32;
        } else {
            disconnect = false;
            break;
        }
    }
    if disconnect && (num_oxygen != 1 || num_h != 4) {
        disconnect = false;
    }
    if disconnect {
        *oxygen_index = found_oxygen;
        *neighbor_position = found_position;
        Ok(1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case)]
pub(crate) fn DisconnectAmmoniumSalt(
    atoms: &mut [inp_ATOM],
    nitrogen_index: i32,
    oxygen_index: i32,
    nitrogen_bond_position: i32,
    num_explicit_h: &[S_CHAR; NUM_H_ISOTOPES as usize + 1],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2398 DisconnectAmmoniumSalt
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: DisconnectAmmoniumSalt
    // INCHI✔️✔️: int DisconnectAmmoniumSalt( inp_ATOM *at,
    // INCHI✔️✔️:                             int iN,
    // INCHI✔️✔️:                             int iO,
    // INCHI✔️✔️:                             int k,
    // INCHI✔️✔️:                             S_CHAR *num_explicit_H )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* disconnect NH4-O from O */
    // INCHI✔️✔️:     /* Note: iO = at[iN].neighbor[k], at[iN] is N, at[iO].neighbor[0] is either N=at[iN] or C=at[iC] */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int nMove_H_iso_diff = -1; /* do not move explicit H */
    // INCHI✔️✔️:     int j, neigh, iso_diff, neigh_pos;
    // INCHI✔️✔️:     int    val = at[iN].valence;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at[iN].charge && !( at[iN].charge + at[iO].charge ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         at[iN].charge = at[iO].charge = 0; /* remove charges */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     neigh_pos = ( at[iO].valence == 2 ) ? ( at[iO].neighbor[1] == iN ) : 0; /* position of at[iN] in the neigh list of iO */
    // INCHI✔️✔️:                                                                             /* disconnect bond O-N */
    // INCHI✔️✔️:     RemoveInpAtBond( at, iO, neigh_pos );
    // INCHI✔️✔️:     RemoveInpAtBond( at, iN, k );
    // INCHI✔️✔️:     val--;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* move 1 H from NH4 to O- or Cl */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* find non-isotopic or the lightest isotopic H to move from N to O */
    // INCHI✔️✔️:     for (iso_diff = 0; iso_diff <= NUM_H_ISOTOPES; iso_diff++) /* djb-rwth: fixing GH PR #72 */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!iso_diff)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* find non-isotopic H */
    // INCHI✔️✔️:             if (at[iN].num_H)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[iN].num_H--;  /* move non-isotopic implicit H */
    // INCHI✔️✔️:                 at[iO].num_H++;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else if (num_explicit_H[0])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nMove_H_iso_diff = 0; /* flag: move explicit non-isotopic H */
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* find isotopic H */
    // INCHI✔️✔️:             /* num_iso_H has length NUM_H_ISOTOPES; do not access out-of-bounds */
    // INCHI✔️✔️:             if ((iso_diff < NUM_H_ISOTOPES) && at[iN].num_iso_H[iso_diff]) /* djb-rwth: fixing GH PR #72 */
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[iN].num_iso_H[iso_diff] --; /* move implicit isotopic H, atw = 1 */
    // INCHI✔️✔️:                 at[iO].num_iso_H[iso_diff] ++;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (num_explicit_H[iso_diff])
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     nMove_H_iso_diff = iso_diff; /* flag: move explicit isotopic H, atw = 1 */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (nMove_H_iso_diff >= 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* move explicit H, it is isotopic if nMove_H_iso_diff > 0 */
    // INCHI✔️✔️:         double dist2_H_O, min_dist2_H_O = -1.0;
    // INCHI✔️✔️:         int    jH = -1, iH = -1;
    // INCHI✔️✔️:         for (j = 0; j < val; j++)
    // INCHI✔️✔️:         { /* looking H in N-H such that H-O is shortest */
    // INCHI✔️✔️:             neigh = at[iN].neighbor[j];
    // INCHI✔️✔️:             if (at[neigh].el_number == EL_NUMBER_H &&
    // INCHI✔️✔️:                  at[neigh].iso_atw_diff == nMove_H_iso_diff)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 dist2_H_O = ( at[neigh].x - at[iO].x ) * ( at[neigh].x - at[iO].x ) +
    // INCHI✔️✔️:                     ( at[neigh].y - at[iO].y ) * ( at[neigh].y - at[iO].y ) +
    // INCHI✔️✔️:                     ( at[neigh].z - at[iO].z ) * ( at[neigh].z - at[iO].z );
    // INCHI✔️✔️:                 if (min_dist2_H_O < 0.0 || min_dist2_H_O > dist2_H_O)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     min_dist2_H_O = dist2_H_O;
    // INCHI✔️✔️:                     iH = neigh;
    // INCHI✔️✔️:                     jH = j;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* reconnect; bonds do not need changes except stereo */
    // INCHI✔️✔️:         neigh_pos = at[iO].valence;
    // INCHI✔️✔️:         at[iO].neighbor[neigh_pos] = iH;
    // INCHI✔️✔️:         at[iO].bond_stereo[neigh_pos] = 0;
    // INCHI✔️✔️:         at[iO].bond_type[neigh_pos] = at[iH].bond_type[0];
    // INCHI✔️✔️:         at[iO].chem_bonds_valence += at[iH].bond_type[0];
    // INCHI✔️✔️:         at[iO].valence++;
    // INCHI✔️✔️:         at[iH].neighbor[0] = iO;
    // INCHI✔️✔️:         at[iH].bond_stereo[0] = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* disconnect H from N */
    // INCHI✔️✔️:         RemoveInpAtBond( at, iN, jH );
    // INCHI✔️✔️:         val--;
    // INCHI✔️✔️:         if (k > jH)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             k--;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: DisconnectAmmoniumSalt
    // END INCHI C FUNCTION: DisconnectAmmoniumSalt
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectAmmoniumSalt
    // INCHI✔️✔️: #define EL_NUMBER_H  ((U_CHAR)1)
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectAmmoniumSalt

    let nitrogen =
        usize::try_from(nitrogen_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let oxygen = usize::try_from(oxygen_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut k = nitrogen_bond_position;
    let mut val = i32::from(
        atoms
            .get(nitrogen)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .valence,
    );
    let oxygen_atom = atoms
        .get(oxygen)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atoms[nitrogen].charge != 0
        && i32::from(atoms[nitrogen].charge) + i32::from(oxygen_atom.charge) == 0
    {
        atoms[nitrogen].charge = 0;
        atoms[oxygen].charge = 0;
    }

    let oxygen_neighbor_position = if atoms[oxygen].valence == 2 {
        i32::from(atoms[oxygen].neighbor[1] == nitrogen as AT_NUMB)
    } else {
        0
    };
    RemoveInpAtBond(atoms, oxygen_index, oxygen_neighbor_position)?;
    RemoveInpAtBond(atoms, nitrogen_index, k)?;
    val -= 1;

    let mut move_explicit_isotope = -1_i32;
    for isotope in 0..=NUM_H_ISOTOPES as usize {
        if isotope == 0 {
            if atoms[nitrogen].num_H != 0 {
                atoms[nitrogen].num_H = atoms[nitrogen].num_H.wrapping_sub(1);
                atoms[oxygen].num_H = atoms[oxygen].num_H.wrapping_add(1);
                break;
            } else if num_explicit_h[0] != 0 {
                move_explicit_isotope = 0;
                break;
            }
        } else if isotope < NUM_H_ISOTOPES as usize && atoms[nitrogen].num_iso_H[isotope] != 0 {
            atoms[nitrogen].num_iso_H[isotope] = atoms[nitrogen].num_iso_H[isotope].wrapping_sub(1);
            atoms[oxygen].num_iso_H[isotope] = atoms[oxygen].num_iso_H[isotope].wrapping_add(1);
            break;
        } else if num_explicit_h[isotope] != 0 {
            move_explicit_isotope = isotope as i32;
            break;
        }
    }

    if move_explicit_isotope >= 0 {
        let mut minimum_distance = -1.0_f64;
        let mut hydrogen = None;
        let mut hydrogen_position = None;
        for j in 0..val.max(0) as usize {
            let neighbor = usize::from(atoms[nitrogen].neighbor[j]);
            let candidate = atoms
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if candidate.el_number == 1
                && i32::from(candidate.iso_atw_diff) == move_explicit_isotope
            {
                let dx = candidate.x - atoms[oxygen].x;
                let dy = candidate.y - atoms[oxygen].y;
                let dz = candidate.z - atoms[oxygen].z;
                let distance = dx * dx + dy * dy + dz * dz;
                if minimum_distance < 0.0 || minimum_distance > distance {
                    minimum_distance = distance;
                    hydrogen = Some(neighbor);
                    hydrogen_position = Some(j);
                }
            }
        }
        let hydrogen = hydrogen.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let hydrogen_position = hydrogen_position.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let oxygen_position = usize::try_from(atoms[oxygen].valence)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if oxygen_position >= MAXVAL as usize {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let hydrogen_bond = atoms[hydrogen].bond_type[0];
        atoms[oxygen].neighbor[oxygen_position] = hydrogen as AT_NUMB;
        atoms[oxygen].bond_stereo[oxygen_position] = 0;
        atoms[oxygen].bond_type[oxygen_position] = hydrogen_bond;
        atoms[oxygen].chem_bonds_valence = atoms[oxygen]
            .chem_bonds_valence
            .wrapping_add(hydrogen_bond as S_CHAR);
        atoms[oxygen].valence = atoms[oxygen].valence.wrapping_add(1);
        atoms[hydrogen].neighbor[0] = oxygen as AT_NUMB;
        atoms[hydrogen].bond_stereo[0] = 0;
        RemoveInpAtBond(atoms, nitrogen_index, hydrogen_position as i32)?;
        val -= 1;
        if k > hydrogen_position as i32 {
            k -= 1;
        }
    }
    let _ = (val, k);
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn bIsMetalSalt(atoms: &[inp_ATOM], atom_index: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2511 bIsMetalSalt
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: bIsMetalSalt
    // INCHI✔️✔️: int bIsMetalSalt( inp_ATOM *at, int i )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int type, val, k, iO, iC, j, neigh;
    // INCHI✔️✔️:     int bDisconnect = 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* check for a metal atom:
    // INCHI✔️✔️:     metal atom should be connected and be a metal */
    // INCHI✔️✔️:     if (!( val = at[i].valence ) ||
    // INCHI✔️✔️:          !( type = get_el_type( at[i].el_number ) ) ||
    // INCHI✔️✔️:          !( type & IS_METAL ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         bDisconnect = 0;  /* reject */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (at[i].num_H)
    // INCHI✔️✔️:         /* metal atom should not have adjacent H or multiple bonds or radical */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         bDisconnect = 0; /* reject */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* check valence */
    // INCHI✔️✔️:         if ((at[i].charge == 0 &&
    // INCHI✔️✔️:             ( (( type & 1 ) && val == get_el_valence( at[i].el_number, 0, 0 )) ||
    // INCHI✔️✔️:              (( type & 2 ) && val == get_el_valence( at[i].el_number, 0, 1 ) ))) ||
    // INCHI✔️✔️:              (at[i].charge > 0 &&
    // INCHI✔️✔️:              ( type & 1 ) && val == get_el_valence( at[i].el_number, at[i].charge, 0 ))) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ; /* accept */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bDisconnect = 0; /* reject */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bDisconnect)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*************************************************************************
    // INCHI✔️✔️:         *                                                                  |    *
    // INCHI✔️✔️:         * check M neighbors. Disconnect if all neighbors are M-O-C# or M-O-C=   *
    // INCHI✔️✔️:         *                                                                  |    *
    // INCHI✔️✔️:         *************************************************************************/
    // INCHI✔️✔️:         for (k = 0; k < at[i].valence; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             iO = at[i].neighbor[k];
    // INCHI✔️✔️:             /* halogenide 2004-07-08 */
    // INCHI✔️✔️:             if (( at[iO].el_number == EL_NUMBER_F ||
    // INCHI✔️✔️:                   at[iO].el_number == EL_NUMBER_CL ||
    // INCHI✔️✔️:                   at[iO].el_number == EL_NUMBER_BR ||
    // INCHI✔️✔️:                   at[iO].el_number == EL_NUMBER_I ) &&
    // INCHI✔️✔️:                  at[iO].valence == 1 && at[iO].chem_bonds_valence == 1 &&
    // INCHI✔️✔️:                  !at[iO].charge && !( at[iO].radical && at[iO].radical != RADICAL_SINGLET ) && !NUMH( at, iO ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ; /* found */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* -O-C= */
    // INCHI✔️✔️:                 if (at[iO].el_number != EL_NUMBER_O ||
    // INCHI✔️✔️:                      NUMH( at, iO ) ||
    // INCHI✔️✔️:                      at[iO].valence != 2 ||
    // INCHI✔️✔️:                      at[iO].charge ||
    // INCHI✔️✔️:                      (at[iO].radical && at[iO].radical != RADICAL_SINGLET) ||
    // INCHI✔️✔️:                      at[iO].valence != at[iO].chem_bonds_valence) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     bDisconnect = 0; /* reject */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 iC = at[iO].neighbor[at[iO].neighbor[0] == i];
    // INCHI✔️✔️:                 if (at[iC].el_number != EL_NUMBER_C ||
    // INCHI✔️✔️:                      at[iC].num_H ||
    // INCHI✔️✔️:                      at[iC].chem_bonds_valence != 4 ||
    // INCHI✔️✔️:                      at[iC].charge ||
    // INCHI✔️✔️:                      (at[iC].radical && at[iC].radical != RADICAL_SINGLET) ||
    // INCHI✔️✔️:                      at[iC].valence == at[iC].chem_bonds_valence) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     bDisconnect = 0; /* reject */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 for (j = 0; j < at[iC].valence; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     neigh = at[iC].neighbor[j];
    // INCHI✔️✔️:                     if (at[neigh].el_number == EL_NUMBER_H)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (j != at[iC].valence)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     bDisconnect = 0; /* reject */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return bDisconnect;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: bIsMetalSalt
    // END INCHI C FUNCTION: bIsMetalSalt
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsMetalSalt
    // INCHI✔️✔️: #define IS_METAL (METAL | METAL2)
    // INCHI✔️✔️: #define METAL 1
    // INCHI✔️✔️: #define METAL2 3
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsMetalSalt

    let i = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let metal = atoms.get(i).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = i32::from(metal.valence);
    let element_type = get_el_type(i32::from(metal.el_number))?;
    if valence == 0 || element_type == 0 || element_type & 3 == 0 || metal.num_H != 0 {
        return Ok(0);
    }

    let valid_valence = if metal.charge == 0 {
        (element_type & 1 != 0 && valence == get_el_valence(i32::from(metal.el_number), 0, 0)?)
            || (element_type & 2 != 0
                && valence == get_el_valence(i32::from(metal.el_number), 0, 1)?)
    } else {
        metal.charge > 0
            && element_type & 1 != 0
            && valence == get_el_valence(i32::from(metal.el_number), i32::from(metal.charge), 0)?
    };
    if !valid_valence {
        return Ok(0);
    }

    for position in 0..valence.max(0) as usize {
        let ligand_index = usize::from(metal.neighbor[position]);
        let ligand = atoms
            .get(ligand_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let ligand_h =
            i32::from(ligand.num_H) + ligand.num_iso_H.into_iter().map(i32::from).sum::<i32>();
        let halogen = matches!(ligand.el_number, 9 | 17 | 35 | 53)
            && ligand.valence == 1
            && ligand.chem_bonds_valence == 1
            && ligand.charge == 0
            && (ligand.radical == 0 || ligand.radical == RADICAL_SINGLET as i8)
            && ligand_h == 0;
        if halogen {
            continue;
        }
        if ligand.el_number != 8
            || ligand_h != 0
            || ligand.valence != 2
            || ligand.charge != 0
            || (ligand.radical != 0 && ligand.radical != RADICAL_SINGLET as i8)
            || ligand.valence != ligand.chem_bonds_valence
        {
            return Ok(0);
        }
        let carbon_position = usize::from(ligand.neighbor[0] == i as AT_NUMB);
        let carbon_index = usize::from(ligand.neighbor[carbon_position]);
        let carbon = atoms
            .get(carbon_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if carbon.el_number != 6
            || carbon.num_H != 0
            || carbon.chem_bonds_valence != 4
            || carbon.charge != 0
            || (carbon.radical != 0 && carbon.radical != RADICAL_SINGLET as i8)
            || carbon.valence == carbon.chem_bonds_valence
        {
            return Ok(0);
        }
        for carbon_position in 0..i32::from(carbon.valence).max(0) as usize {
            let neighbor = usize::from(carbon.neighbor[carbon_position]);
            if atoms
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .el_number
                == 1
            {
                return Ok(0);
            }
        }
    }
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn DisconnectMetalSalt(
    atoms: &mut [inp_ATOM],
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2612 DisconnectMetalSalt
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: DisconnectMetalSalt
    // INCHI✔️✔️: int DisconnectMetalSalt( inp_ATOM *at, int i )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k, iO;
    // INCHI✔️✔️:     /* disconnect metal atom or ion at[i] */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (k = 0; k < at[i].valence; k++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         iO = at[i].neighbor[k];
    // INCHI✔️✔️:         if (at[iO].valence == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at[iO].neighbor[0] == i)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* assuming atom O always has 2 bonds */
    // INCHI✔️✔️:                 /* copy the remaining neighbor to the 0 position */
    // INCHI✔️✔️:                 at[iO].neighbor[0] = at[iO].neighbor[1];
    // INCHI✔️✔️:                 at[iO].bond_stereo[0] = at[iO].bond_stereo[1];
    // INCHI✔️✔️:                 at[iO].bond_type[0] = at[iO].bond_type[1];
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             /* clear neighbor at position 1 */
    // INCHI✔️✔️:             at[iO].neighbor[1] = 0;
    // INCHI✔️✔️:             at[iO].bond_stereo[1] = 0;
    // INCHI✔️✔️:             at[iO].bond_type[1] = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* clear neighbor at position 1 */
    // INCHI✔️✔️:             at[iO].neighbor[0] = 0;
    // INCHI✔️✔️:             at[iO].bond_stereo[0] = 0;
    // INCHI✔️✔️:             at[iO].bond_type[0] = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* make O negatively charged */
    // INCHI✔️✔️:         at[iO].charge = -1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* reduce O valence to account for the removed single bond */
    // INCHI✔️✔️:         at[iO].valence--;
    // INCHI✔️✔️:         at[iO].chem_bonds_valence--;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* clear metal neighbor (O) */
    // INCHI✔️✔️:         at[i].neighbor[k] = 0;
    // INCHI✔️✔️:         at[i].bond_stereo[k] = 0;
    // INCHI✔️✔️:         at[i].bond_type[k] = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* add a positive charge to the metal */
    // INCHI✔️✔️:         at[i].charge++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* set metal valence to zero because it has been disconnected */
    // INCHI✔️✔️:     at[i].valence = 0;
    // INCHI✔️✔️:     at[i].chem_bonds_valence = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return k;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: DisconnectMetalSalt
    // END INCHI C FUNCTION: DisconnectMetalSalt
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectMetalSalt
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectMetalSalt

    let metal_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let metal_valence = i32::from(
        atoms
            .get(metal_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .valence,
    );
    for position in 0..metal_valence.max(0) as usize {
        let ligand_index = usize::from(
            *atoms[metal_index]
                .neighbor
                .get(position)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let ligand_valence = atoms
            .get(ligand_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .valence;
        if ligand_valence == 2 {
            if atoms[ligand_index].neighbor[0] == atom_index as AT_NUMB {
                atoms[ligand_index].neighbor[0] = atoms[ligand_index].neighbor[1];
                atoms[ligand_index].bond_stereo[0] = atoms[ligand_index].bond_stereo[1];
                atoms[ligand_index].bond_type[0] = atoms[ligand_index].bond_type[1];
            }
            atoms[ligand_index].neighbor[1] = 0;
            atoms[ligand_index].bond_stereo[1] = 0;
            atoms[ligand_index].bond_type[1] = 0;
        } else {
            atoms[ligand_index].neighbor[0] = 0;
            atoms[ligand_index].bond_stereo[0] = 0;
            atoms[ligand_index].bond_type[0] = 0;
        }
        atoms[ligand_index].charge = -1;
        atoms[ligand_index].valence = atoms[ligand_index].valence.wrapping_sub(1);
        atoms[ligand_index].chem_bonds_valence =
            atoms[ligand_index].chem_bonds_valence.wrapping_sub(1);
        atoms[metal_index].neighbor[position] = 0;
        atoms[metal_index].bond_stereo[position] = 0;
        atoms[metal_index].bond_type[position] = 0;
        atoms[metal_index].charge = atoms[metal_index].charge.wrapping_add(1);
    }
    atoms[metal_index].valence = 0;
    atoms[metal_index].chem_bonds_valence = 0;
    Ok(metal_valence.max(0))
}

#[allow(non_snake_case)]
pub(crate) fn DisconnectSalts(
    heap: &mut SourceHeap,
    original: &mut ORIG_ATOM_DATA,
    disconnect: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2668 DisconnectSalts
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: DisconnectSalts
    // INCHI✔️✔️: int DisconnectSalts( ORIG_ATOM_DATA *orig_inp_data, int bDisconnect )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, k, iO, num_changes, val;
    // INCHI✔️✔️:     S_CHAR    num_explicit_H[NUM_H_ISOTOPES + 1];
    // INCHI✔️✔️:     inp_ATOM *at = orig_inp_data->at;
    // INCHI✔️✔️:     int num_at = orig_inp_data->num_inp_atoms;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* check each atom */
    // INCHI✔️✔️:     for (i = 0, num_changes = 0; i < num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!( val = at[i].valence ) || /* disconnected atom */
    // INCHI✔️✔️:              val != at[i].chem_bonds_valence || /* a bond has higher multiplicity than 1 */
    // INCHI✔️✔️:              (at[i].radical && at[i].radical != RADICAL_SINGLET) /* radical */) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;   /* reject */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (bDisconnect)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 DisconnectAmmoniumSalt( at, i, iO, k, num_explicit_H );
    // INCHI✔️✔️:                 orig_inp_data->num_inp_bonds--;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:
    // INCHI✔️✔️:             /* count disconnected atoms */
    // INCHI✔️✔️:             num_changes++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (bIsMetalSalt( at, i ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (bDisconnect)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 k = DisconnectMetalSalt( at, i );
    // INCHI✔️✔️:                 orig_inp_data->num_inp_bonds -= k;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             num_changes++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_changes;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: DisconnectSalts
    // END INCHI C FUNCTION: DisconnectSalts
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectSalts
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectSalts

    let atom_pointer = original.at;
    let number_of_atoms = original.num_inp_atoms;
    if number_of_atoms <= 0 {
        return Ok(0);
    }
    let atoms = heap.slice_mut(atom_pointer)?;
    let mut number_of_changes = 0_i32;
    for atom_index in 0..number_of_atoms.max(0) as usize {
        let atom = atoms
            .get(atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = i32::from(atom.valence);
        if valence == 0
            || valence != i32::from(atom.chem_bonds_valence)
            || (atom.radical != 0 && atom.radical != RADICAL_SINGLET as S_CHAR)
        {
            continue;
        }

        let source_index = atom_index as i32;
        let mut ligand_index = 0_i32;
        let mut neighbor_position = 0_i32;
        let mut explicit_hydrogens = [0_i8; NUM_H_ISOTOPES as usize + 1];
        if bIsAmmoniumSalt(
            atoms,
            source_index,
            &mut ligand_index,
            &mut neighbor_position,
            &mut explicit_hydrogens,
        )? != 0
        {
            if disconnect != 0 {
                DisconnectAmmoniumSalt(
                    atoms,
                    source_index,
                    ligand_index,
                    neighbor_position,
                    &explicit_hydrogens,
                )?;
                original.num_inp_bonds = original.num_inp_bonds.wrapping_sub(1);
            }
            number_of_changes = number_of_changes.wrapping_add(1);
        } else if bIsMetalSalt(atoms, source_index)? != 0 {
            if disconnect != 0 {
                let disconnected_bonds = DisconnectMetalSalt(atoms, source_index)?;
                original.num_inp_bonds = original.num_inp_bonds.wrapping_sub(disconnected_bonds);
            }
            number_of_changes = number_of_changes.wrapping_add(1);
        }
    }
    Ok(number_of_changes)
}

#[allow(non_snake_case)]
pub(crate) fn bIsMetalToDisconnect(
    atoms: &[inp_ATOM],
    atom_index: i32,
    check_metal_valence: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2719 bIsMetalToDisconnect
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: bIsMetalToDisconnect
    // INCHI✔️✔️: int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int type, at_valence, num_H;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     if ( !at[i].valence )
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!( type = get_el_type( at[i].el_number ) ) ||
    // INCHI✔️✔️:          !( type & IS_METAL ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     num_H = NUMH( at, i );
    // INCHI✔️✔️:     at_valence = num_H + at[i].chem_bonds_valence;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!at_valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0; /* nothing to disconnect */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bCheckMetalValence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (abs( at[i].charge ) > 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 1; /* multiple charges */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         for (i = 0; i < 2 && ( i & type ); i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at_valence == get_el_valence( at[i].el_number, at[i].charge, i )) /* djb-rwth: fixing coverity ID #499532 -- unresolved issue -- revision required */
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return 2; /* atom has normal valence */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: bIsMetalToDisconnect
    // END INCHI C FUNCTION: bIsMetalToDisconnect
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsMetalToDisconnect
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define IS_METAL (METAL | METAL2)
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // INCHI✔️✔️: int abs(int);
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsMetalToDisconnect

    let index = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let element_type = get_el_type(i32::from(atom.el_number))?;
    if element_type == 0 || element_type & 3 == 0 {
        return Ok(0);
    }
    let number_of_hydrogens =
        i32::from(atom.num_H) + atom.num_iso_H.into_iter().map(i32::from).sum::<i32>();
    let atom_valence = number_of_hydrogens + i32::from(atom.chem_bonds_valence);
    if atom_valence == 0 {
        return Ok(0);
    }
    if check_metal_valence != 0 {
        if i32::from(atom.charge).abs() > 1 {
            return Ok(1);
        }
        let mut source_index = 0_i32;
        while source_index < 2 && source_index & element_type != 0 {
            let source_atom = atoms
                .get(source_index as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom_valence
                == get_el_valence(
                    i32::from(source_atom.el_number),
                    i32::from(source_atom.charge),
                    source_index,
                )?
            {
                return Ok(2);
            }
            source_index += 1;
        }
    }
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn bMayDisconnectMetals(
    heap: &mut SourceHeap,
    original: &mut ORIG_ATOM_DATA,
    check_metal_valence: i32,
    mut taut_flags_done: Option<&mut INCHI_MODE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2762 bMayDisconnectMetals
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: bMayDisconnectMetals
    // INCHI✔️✔️: int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️✔️:                           int bCheckMetalValence,
    // INCHI✔️✔️:                           INCHI_MODE *bTautFlagsDone )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j, k, iO, num_changes, val, bRadOrMultBonds, num_impl_H = 0;
    // INCHI✔️✔️:     S_CHAR    num_explicit_H[NUM_H_ISOTOPES + 1];
    // INCHI✔️✔️:     inp_ATOM *at = orig_inp_data->at;
    // INCHI✔️✔️:     int num_at = orig_inp_data->num_inp_atoms;
    // INCHI✔️✔️:     int *nNumImplH = &orig_inp_data->bDisconnectCoord;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* check each atom */
    // INCHI✔️✔️:     for (i = 0, num_changes = 0; i < num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!( val = at[i].valence ) && !NUMH( at, i ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue; /* disconnected atom */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         bRadOrMultBonds = ( val == 0 ) ||
    // INCHI✔️✔️:             ( val != at[i].chem_bonds_valence ) || /* a bond has higher multiplicity than 1 */
    // INCHI✔️✔️:             ( at[i].radical && at[i].radical != RADICAL_SINGLET ); /* radical */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!bRadOrMultBonds && bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (!bRadOrMultBonds && bIsMetalSalt( at, i ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (1 == ( j = bIsMetalToDisconnect( at, i, bCheckMetalValence ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             num_impl_H += NUMH( at, i );
    // INCHI✔️✔️:             num_changes++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (2 == j && bTautFlagsDone)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *bTautFlagsDone |= TG_FLAG_CHECK_VALENCE_COORD_DONE;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (nNumImplH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *nNumImplH = num_changes ? num_impl_H + 1 : 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_changes;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: bMayDisconnectMetals
    // END INCHI C FUNCTION: bMayDisconnectMetals
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bMayDisconnectMetals
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define TG_FLAG_CHECK_VALENCE_COORD_DONE 0x00000200
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bMayDisconnectMetals

    let number_of_atoms = original.num_inp_atoms;
    if number_of_atoms <= 0 {
        original.bDisconnectCoord = 0;
        return Ok(0);
    }
    let atoms = heap.slice(original.at.as_const())?;
    let mut number_of_changes = 0_i32;
    let mut number_of_implicit_hydrogens = 0_i32;
    for atom_index in 0..number_of_atoms as usize {
        let atom = atoms
            .get(atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = i32::from(atom.valence);
        let number_of_hydrogens =
            i32::from(atom.num_H) + atom.num_iso_H.into_iter().map(i32::from).sum::<i32>();
        if valence == 0 && number_of_hydrogens == 0 {
            continue;
        }
        let radical_or_multiple_bonds = valence == 0
            || valence != i32::from(atom.chem_bonds_valence)
            || (atom.radical != 0 && atom.radical != RADICAL_SINGLET as S_CHAR);
        let source_index = atom_index as i32;
        let mut ligand_index = 0_i32;
        let mut neighbor_position = 0_i32;
        let mut explicit_hydrogens = [0_i8; NUM_H_ISOTOPES as usize + 1];
        if !radical_or_multiple_bonds
            && bIsAmmoniumSalt(
                atoms,
                source_index,
                &mut ligand_index,
                &mut neighbor_position,
                &mut explicit_hydrogens,
            )? != 0
        {
            continue;
        }
        if !radical_or_multiple_bonds && bIsMetalSalt(atoms, source_index)? != 0 {
            continue;
        }
        let metal_status = bIsMetalToDisconnect(atoms, source_index, check_metal_valence)?;
        if metal_status == 1 {
            number_of_implicit_hydrogens =
                number_of_implicit_hydrogens.wrapping_add(number_of_hydrogens);
            number_of_changes = number_of_changes.wrapping_add(1);
        } else if metal_status == 2
            && let Some(flags) = taut_flags_done.as_deref_mut()
        {
            *flags |= TG_FLAG_CHECK_VALENCE_COORD_DONE as INCHI_MODE;
        }
    }
    original.bDisconnectCoord = if number_of_changes != 0 {
        number_of_implicit_hydrogens.wrapping_add(1)
    } else {
        0
    };
    Ok(number_of_changes)
}

#[allow(non_snake_case)]
pub(crate) fn DisconnectMetals(
    heap: &mut SourceHeap,
    original: &mut ORIG_ATOM_DATA,
    check_metal_valence: i32,
    mut taut_flags_done: Option<&mut INCHI_MODE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2865 DisconnectMetals
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: DisconnectMetals
    // INCHI✔️✔️: int DisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️✔️:                       int bCheckMetalValence,
    // INCHI✔️✔️:                       INCHI_MODE *bTautFlagsDone )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j, k, n, iO, num_changes, val, bRadOrMultBonds;
    // INCHI✔️✔️:     int num_impl_H, num_at, err, num_disconnected;
    // INCHI✔️✔️:     S_CHAR num_explicit_H[NUM_H_ISOTOPES + 1];
    // INCHI✔️✔️:     static char elnumber_Heteroat[16] = { '\0', };
    // INCHI✔️✔️:     static int  num_halogens = 0;
    // INCHI✔️✔️:     int num_halogens2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inp_ATOM  *at = NULL;
    // INCHI✔️✔️:     S_CHAR    *bMetal = NULL;
    // INCHI✔️✔️:     inp_ATOM  *atom = orig_inp_data->at;
    // INCHI✔️✔️:     int        num_atoms = orig_inp_data->num_inp_atoms;
    // INCHI✔️✔️:     int        nNumExplH = ( orig_inp_data->bDisconnectCoord > 0 ) ? orig_inp_data->bDisconnectCoord - 1 : 0;
    // INCHI✔️✔️:     AT_NUMB   *nOldCompNumber = orig_inp_data->nOldCompNumber;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     err = 0;
    // INCHI✔️✔️:     num_impl_H = 0;
    // INCHI✔️✔️:     num_at = num_atoms;
    // INCHI✔️✔️:     num_disconnected = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!( at = (inp_ATOM *) inchi_calloc( (long long)num_at + (long long)nNumExplH, sizeof( at[0] ) ) ) || /* djb-rwth: cast operators added */
    // INCHI✔️✔️:          !( bMetal = (S_CHAR    *) inchi_calloc( (long long)num_at + (long long)nNumExplH, sizeof( bMetal[0] ) ) )) /* djb-rwth: cast operators added */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         err = 1;
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!num_halogens) /* if (!elnumber_Heteroat[0] )  */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         i = 0;
    // INCHI✔️✔️:         /* halogens */
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_F; /* 0 */
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_CL;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_BR;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_I;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_AT; /* 4 */
    // INCHI✔️✔️:         num_halogens2 = i;
    // INCHI✔️✔️:         /* other non-metal */
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_N;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_P;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_AS;
    // INCHI✔️✔️:         /*elnumber_Heteroat[i++] = EL_NUMBER_SB;*/ /* metal 10-28-2003 */
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_O;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_S;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_SE;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_TE;
    // INCHI✔️✔️:         /*elnumber_Heteroat[i++] = EL_NUMBER_PO;*/ /* metal 10-28-2003 */
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = (char) EL_NUMBER_B;
    // INCHI✔️✔️:         elnumber_Heteroat[i++] = 0;
    // INCHI✔️✔️:         num_halogens = num_halogens2;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     memcpy(at, atom, num_atoms * sizeof(at[0]));
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* check each atom, mark metals */
    // INCHI✔️✔️:     for (i = 0, k = 0, num_changes = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!( val = at[i].valence ) && !NUMH( at, i ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue; /* disconnected atom */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         bRadOrMultBonds = ( val == 0 ) ||
    // INCHI✔️✔️:             ( val != at[i].chem_bonds_valence ) || /* a bond has higher multiplicity than 1 */
    // INCHI✔️✔️:             ( at[i].radical && at[i].radical != RADICAL_SINGLET ); /* radical */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!bRadOrMultBonds && bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (!bRadOrMultBonds && bIsMetalSalt( at, i ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (1 == ( j = bIsMetalToDisconnect( at, i, bCheckMetalValence ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             num_impl_H += ( k = NUMH( at, i ) );
    // INCHI✔️✔️:             bMetal[i] = 1 + k;
    // INCHI✔️✔️:             num_changes++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (2 == j && bTautFlagsDone)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *bTautFlagsDone |= TG_FLAG_CHECK_VALENCE_COORD_DONE;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_impl_H != nNumExplH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         err = 2;
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* replace implicit H atoms with explicit H atoms */
    // INCHI✔️✔️:     for (i = 0; i < num_atoms && 0 < num_impl_H; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (bMetal[i] <= 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         for (k = 0; k < NUM_H_ISOTOPES + 1; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             n = k ? at[i].num_iso_H[k - 1] : at[i].num_H;
    // INCHI✔️✔️:             for (j = 0; j < n; j++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (num_at >= num_atoms + nNumExplH)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     err = 3;
    // INCHI✔️✔️:                     goto exit_function;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 at[num_at].elname[0] = 'H';
    // INCHI✔️✔️:                 at[num_at].el_number = get_periodic_table_number( at[num_at].elname );
    // INCHI✔️✔️:                 at[num_at].iso_atw_diff = k;
    // INCHI✔️✔️:                 at[num_at].component = at[i].component;
    // INCHI✔️✔️:                 move_explicit_Hcation( at, num_at + 1, i, num_at, 1 );
    // INCHI✔️✔️:                 at[num_at].orig_at_number = num_at + 1;
    // INCHI✔️✔️:                 num_at++;
    // INCHI✔️✔️:                 num_impl_H--;
    // INCHI✔️✔️:                 bMetal[i] --;
    // INCHI✔️✔️:                 if (k)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     at[i].num_iso_H[k - 1] --;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     at[i].num_H--;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (bMetal[i] != 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             err = 4;
    // INCHI✔️✔️:             goto exit_function;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_at != num_atoms + nNumExplH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         err = 5;
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* disconnect metal - ligand bonds */
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!bMetal[i])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* disconnect metal atom M
    // INCHI✔️✔️:
    // INCHI✔️✔️:         Note: Defect in case of bridging ligands:
    // INCHI✔️✔️:
    // INCHI✔️✔️:         M     M                          M     M             M     M(+)
    // INCHI✔️✔️:         \  /   will be transformed to            , not to
    // INCHI✔️✔️:         N(+)                             N(+)                N(-)
    // INCHI✔️✔️:         / \                              / \                 / \
    // INCHI✔️✔️:         R   R                            R   R               R   R
    // INCHI✔️✔️:
    // INCHI✔️✔️:         Non-bridging are OK:
    // INCHI✔️✔️:
    // INCHI✔️✔️:         M     R           M(+)  R
    // INCHI✔️✔️:         \  /                 /
    // INCHI✔️✔️:         N(+)    --->      N
    // INCHI✔️✔️:         / \               / \
    // INCHI✔️✔️:         R   R             R   R
    // INCHI✔️✔️:
    // INCHI✔️✔️:         */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         for (j = at[i].valence - 1; 0 <= j; j--)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (j < at[i].valence && !bMetal[(int) at[i].neighbor[j]])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* do not break metal-metal bond here */
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 num_disconnected += DisconnectOneLigand( at,
    // INCHI✔️✔️:                                                          nOldCompNumber,
    // INCHI✔️✔️:                                                          bMetal,
    // INCHI✔️✔️:                                                          elnumber_Heteroat,
    // INCHI✔️✔️:                                                          num_halogens,
    // INCHI✔️✔️:                                                          num_atoms,
    // INCHI✔️✔️:                                                          i,
    // INCHI✔️✔️:                                                          j,
    // INCHI✔️✔️:                                                          bTautFlagsDone );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* disconnect metal-metal bonds */
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!bMetal[i])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         for (j = at[i].valence - 1; 0 <= j; j--)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (j < at[i].valence && bMetal[(int) at[i].neighbor[j]])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* break metal-metal bond here */
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 num_disconnected += DisconnectOneLigand( at,
    // INCHI✔️✔️:                                                          nOldCompNumber,
    // INCHI✔️✔️:                                                          bMetal,
    // INCHI✔️✔️:                                                          elnumber_Heteroat,
    // INCHI✔️✔️:                                                          num_halogens,
    // INCHI✔️✔️:                                                          num_atoms,
    // INCHI✔️✔️:                                                          i,
    // INCHI✔️✔️:                                                          j,
    // INCHI✔️✔️:                                                          bTautFlagsDone );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_function:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!num_disconnected)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         err = 6;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (at && err)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free( at );
    // INCHI✔️✔️:         at = NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (atom && at)
    // INCHI✔️✔️:     {    /* changed if ( at ) to if ( atom && at ) 2004-04-03 */
    // INCHI✔️✔️:         inchi_free( atom );
    // INCHI✔️✔️:         atom = NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (bMetal)
    // INCHI✔️✔️:         inchi_free( bMetal );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         orig_inp_data->at = at;
    // INCHI✔️✔️:         orig_inp_data->num_inp_atoms = num_at;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return err ? -err : num_disconnected;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: DisconnectMetals
    // END INCHI C FUNCTION: DisconnectMetals
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectMetals
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define NUM_H_ISOTOPES 3
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define TG_FLAG_CHECK_VALENCE_COORD_DONE 0x00000200
    // INCHI✔️✔️: inchi_calloc and inchi_free resolve to the GCC/Linux libc allocation macros.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectMetals
    let number_of_atoms = original.num_inp_atoms;
    let number_of_explicit_hydrogens = if original.bDisconnectCoord > 0 {
        original.bDisconnectCoord - 1
    } else {
        0
    };
    let total_count = i64::from(number_of_atoms) + i64::from(number_of_explicit_hydrogens);
    if total_count < 0 {
        return Ok(-6);
    }
    let total_count = total_count as u64;
    let atom_pointer =
        match inchi_calloc::<inp_ATOM>(heap, total_count, std::mem::size_of::<inp_ATOM>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => return Ok(-6),
            Err(error) => return Err(error),
        };
    let metal_pointer =
        match inchi_calloc::<S_CHAR>(heap, total_count, std::mem::size_of::<S_CHAR>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                inchi_free(heap, atom_pointer)?;
                return Ok(-6);
            }
            Err(error) => {
                inchi_free(heap, atom_pointer)?;
                return Err(error);
            }
        };
    if number_of_atoms < 0 {
        inchi_free(heap, atom_pointer)?;
        inchi_free(heap, metal_pointer)?;
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    let copy_result = heap.with_slice_mut_and_heap(atom_pointer, |destination, heap| {
        let source = heap.slice(original.at.as_const())?;
        let count = number_of_atoms as usize;
        if source.len() < count || destination.len() < count {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        destination[..count].clone_from_slice(&source[..count]);
        Ok(())
    });
    if let Err(error) = copy_result {
        inchi_free(heap, atom_pointer)?;
        inchi_free(heap, metal_pointer)?;
        return Err(error);
    }

    const HETEROATOM_ELEMENTS: [i8; 14] = [9, 17, 35, 53, 85, 7, 15, 33, 8, 16, 34, 52, 5, 0];
    let old_components = (!original.nOldCompNumber.is_null()).then_some(original.nOldCompNumber);
    let processing = heap.with_two_slices_mut_and_optional_third(
        atom_pointer,
        metal_pointer,
        old_components,
        |atoms, metal_marks, mut old_component_numbers| {
            let mut number_of_implicit_hydrogens = 0_i32;
            let mut number_of_changes = 0_i32;
            let mut number_of_disconnections = 0_i32;
            let mut number_of_atoms_after_expansion = number_of_atoms;

            for index in 0..number_of_atoms as usize {
                let valence = i32::from(atoms[index].valence);
                let hydrogen_count = i32::from(atoms[index].num_H)
                    + atoms[index]
                        .num_iso_H
                        .into_iter()
                        .map(i32::from)
                        .sum::<i32>();
                if valence == 0 && hydrogen_count == 0 {
                    continue;
                }
                let radical_or_multiple = valence == 0
                    || valence != i32::from(atoms[index].chem_bonds_valence)
                    || (atoms[index].radical != 0
                        && atoms[index].radical != RADICAL_SINGLET as S_CHAR);
                let mut oxygen = 0_i32;
                let mut ordinal = 0_i32;
                let mut explicit_hydrogens = [0_i8; NUM_H_ISOTOPES as usize + 1];
                if !radical_or_multiple
                    && bIsAmmoniumSalt(
                        atoms,
                        index as i32,
                        &mut oxygen,
                        &mut ordinal,
                        &mut explicit_hydrogens,
                    )? != 0
                {
                } else if !radical_or_multiple && bIsMetalSalt(atoms, index as i32)? != 0 {
                } else {
                    let status = bIsMetalToDisconnect(atoms, index as i32, check_metal_valence)?;
                    if status == 1 {
                        number_of_implicit_hydrogens =
                            number_of_implicit_hydrogens.wrapping_add(hydrogen_count);
                        metal_marks[index] = (1_i32 + hydrogen_count) as S_CHAR;
                        number_of_changes = number_of_changes.wrapping_add(1);
                    } else if status == 2
                        && let Some(flags) = taut_flags_done.as_deref_mut()
                    {
                        *flags |= TG_FLAG_CHECK_VALENCE_COORD_DONE as INCHI_MODE;
                    }
                }
            }
            if number_of_implicit_hydrogens != number_of_explicit_hydrogens {
                return Ok((
                    2_i32,
                    number_of_atoms_after_expansion,
                    number_of_disconnections,
                ));
            }

            for index in 0..number_of_atoms as usize {
                if number_of_implicit_hydrogens <= 0 {
                    break;
                }
                if metal_marks[index] <= 1 {
                    continue;
                }
                for isotope in 0..=NUM_H_ISOTOPES as usize {
                    let count = if isotope == 0 {
                        i32::from(atoms[index].num_H)
                    } else {
                        i32::from(atoms[index].num_iso_H[isotope - 1])
                    };
                    for _ in 0..count.max(0) {
                        if number_of_atoms_after_expansion
                            >= number_of_atoms + number_of_explicit_hydrogens
                        {
                            return Ok((
                                3,
                                number_of_atoms_after_expansion,
                                number_of_disconnections,
                            ));
                        }
                        let hydrogen = number_of_atoms_after_expansion as usize;
                        atoms[hydrogen].elname[0] = b'H' as i8;
                        atoms[hydrogen].el_number =
                            get_periodic_table_number(Some(&atoms[hydrogen].elname))? as U_CHAR;
                        atoms[hydrogen].iso_atw_diff = isotope as S_CHAR;
                        atoms[hydrogen].component = atoms[index].component;
                        let _ = move_explicit_Hcation(
                            atoms,
                            number_of_atoms_after_expansion + 1,
                            index as i32,
                            number_of_atoms_after_expansion,
                            1,
                        )?;
                        atoms[hydrogen].orig_at_number =
                            (number_of_atoms_after_expansion + 1) as AT_NUMB;
                        number_of_atoms_after_expansion += 1;
                        number_of_implicit_hydrogens -= 1;
                        metal_marks[index] = metal_marks[index].wrapping_sub(1);
                        if isotope == 0 {
                            atoms[index].num_H = atoms[index].num_H.wrapping_sub(1);
                        } else {
                            atoms[index].num_iso_H[isotope - 1] =
                                atoms[index].num_iso_H[isotope - 1].wrapping_sub(1);
                        }
                    }
                }
                if metal_marks[index] != 1 {
                    return Ok((4, number_of_atoms_after_expansion, number_of_disconnections));
                }
            }
            if number_of_atoms_after_expansion != number_of_atoms + number_of_explicit_hydrogens {
                return Ok((5, number_of_atoms_after_expansion, number_of_disconnections));
            }

            for disconnect_metal_bonds in [false, true] {
                for index in 0..number_of_atoms as usize {
                    if metal_marks[index] == 0 {
                        continue;
                    }
                    let mut ordinal = i32::from(atoms[index].valence) - 1;
                    while ordinal >= 0 {
                        if ordinal < i32::from(atoms[index].valence) {
                            let neighbor = usize::from(atoms[index].neighbor[ordinal as usize]);
                            let neighbor_is_metal = *metal_marks
                                .get(neighbor)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                != 0;
                            if neighbor_is_metal == disconnect_metal_bonds {
                                number_of_disconnections =
                                    number_of_disconnections.wrapping_add(DisconnectOneLigand(
                                        atoms,
                                        old_component_numbers.as_deref_mut(),
                                        metal_marks,
                                        &HETEROATOM_ELEMENTS,
                                        5,
                                        number_of_atoms,
                                        index as i32,
                                        ordinal,
                                        taut_flags_done.as_deref_mut(),
                                    )?);
                            }
                        }
                        ordinal -= 1;
                    }
                }
            }
            Ok((0, number_of_atoms_after_expansion, number_of_disconnections))
        },
    );

    inchi_free(heap, metal_pointer)?;
    let (mut error, number_of_atoms_after_expansion, number_of_disconnections) = match processing {
        Ok(result) => result,
        Err(error) => {
            inchi_free(heap, atom_pointer)?;
            return Err(error);
        }
    };
    if number_of_disconnections == 0 {
        error = 6;
    }
    if error != 0 {
        inchi_free(heap, atom_pointer)?;
        return Ok(-error);
    }
    inchi_free(heap, original.at)?;
    original.at = atom_pointer;
    original.num_inp_atoms = number_of_atoms_after_expansion;
    Ok(number_of_disconnections)
}

pub(crate) fn get_iat_number(element_number: i32) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4024 get_iat_number
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: get_iat_number
    // INCHI✔️✔️: int get_iat_number( int el_number )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     switch (el_number) {
    // INCHI✔️✔️:         case EL_NUMBER_H:  return IAT_H;
    // INCHI✔️✔️:         case EL_NUMBER_C:  return IAT_C;
    // INCHI✔️✔️:         case EL_NUMBER_N:  return IAT_N;
    // INCHI✔️✔️:         case EL_NUMBER_P:  return IAT_P;
    // INCHI✔️✔️:         case EL_NUMBER_O:  return IAT_O;
    // INCHI✔️✔️:         case EL_NUMBER_S:  return IAT_S;
    // INCHI✔️✔️:         case EL_NUMBER_SE: return IAT_Se;
    // INCHI✔️✔️:         case EL_NUMBER_TE: return IAT_Te;
    // INCHI✔️✔️:         case EL_NUMBER_F:  return IAT_F;
    // INCHI✔️✔️:         case EL_NUMBER_CL: return IAT_Cl;
    // INCHI✔️✔️:         case EL_NUMBER_BR: return IAT_Br;
    // INCHI✔️✔️:         case EL_NUMBER_I:  return IAT_I;
    // INCHI✔️✔️:         default: return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: get_iat_number
    // END INCHI C FUNCTION: get_iat_number
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_iat_number
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: EL_NUMBER_H=1, EL_NUMBER_C=6, EL_NUMBER_N=7, EL_NUMBER_O=8
    // INCHI✔️✔️: EL_NUMBER_F=9, EL_NUMBER_P=15, EL_NUMBER_S=16, EL_NUMBER_CL=17
    // INCHI✔️✔️: EL_NUMBER_SE=34, EL_NUMBER_BR=35, EL_NUMBER_TE=52, EL_NUMBER_I=53
    // INCHI✔️✔️: IAT_H..IAT_I are enum values 0..11; IAT_MAX=12.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_iat_number
    match element_number {
        1 => 0,
        6 => 1,
        7 => 2,
        15 => 3,
        8 => 4,
        16 => 5,
        34 => 6,
        52 => 7,
        9 => 8,
        17 => 9,
        35 => 10,
        53 => 11,
        _ => -1,
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn bHeteroAtomMayHaveXchgIsoH(
    atoms: &[inp_ATOM],
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4047 bHeteroAtomMayHaveXchgIsoH
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat )
{
    inp_ATOM *at = atom + iat, *at2;
    int j, val, is_H = 0, num_H, iat_numb, bAccept; /* djb-rwth: removing redundant variables */

    if (0 > ( iat_numb = get_iat_number( at->el_number ) ))
    {
        return 0;
    }

    if (abs( at->charge ) > 1 || (at->radical && RADICAL_SINGLET != at->radical)) /* djb-rwth: addressing LLVM warning */
    {
        return 0;
    }

    val = -1;
    switch (iat_numb)
    {
        case IAT_N:
        case IAT_P:
            /* djb-rwth: removing redundant code */
            val = 3 + at->charge;
            break;

        case IAT_O:
        case IAT_S:
        case IAT_Se:
        case IAT_Te:
            /* djb-rwth: removing redundant code */
            val = 2 + at->charge;
            break;

        case IAT_F:
        case IAT_Cl:
        case IAT_Br:
        case IAT_I:
            if (at->charge == 0)
            {
                /* djb-rwth: removing redundant code */
                val = 1;
            }
            break;

        case IAT_H:
            if (at->valence == 0 &&
                 at->charge == 1)
            {
                is_H = 1; /* isolated proton */
                val = 0;
            }
    }
    if (val < 0)
    {
        return 0;
    }
    num_H = NUMH( at, 0 );
    if (val != at->chem_bonds_valence + num_H)
    {
        return 0;
    }
    if (is_H)
    {
        return 2; /* H atom */
    }
    else
    {
        /* djb-rwth: removing redundant code */
        for (j = 0, bAccept = 1; j < at->valence && bAccept; j++)
        {
            at2 = atom + (int) at->neighbor[j];
            if ((at2->charge && at->charge) ||
                ( at2->radical && RADICAL_SINGLET != at2->radical )) /* djb-rwth: addressing LLVM warning */
            {
                return 0; /* adjacent charged/radical atoms: do not neutralizate */
            }
        }
    }

    return 1;
}
    */
    // END INCHI C FUNCTION: bHeteroAtomMayHaveXchgIsoH
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bHeteroAtomMayHaveXchgIsoH
    // INCHI✔️✔️: READ_INCHI_STRING=1 includes this function for the selected production target.
    // INCHI✔️✔️: #define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bHeteroAtomMayHaveXchgIsoH

    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ion_atom_type = get_iat_number(i32::from(atom.el_number));
    if ion_atom_type < 0 {
        return Ok(0);
    }
    if i32::from(atom.charge).abs() > 1
        || (atom.radical != 0 && atom.radical != RADICAL_SINGLET as S_CHAR)
    {
        return Ok(0);
    }

    let mut valence = -1_i32;
    let mut is_hydrogen = false;
    match ion_atom_type {
        2 | 3 => valence = 3_i32.wrapping_add(i32::from(atom.charge)),
        4..=7 => valence = 2_i32.wrapping_add(i32::from(atom.charge)),
        8..=11 => {
            if atom.charge == 0 {
                valence = 1;
            }
        }
        0 => {
            if atom.valence == 0 && atom.charge == 1 {
                is_hydrogen = true;
                valence = 0;
            }
        }
        _ => {}
    }
    if valence < 0 {
        return Ok(0);
    }
    let number_of_hydrogens = atom.num_iso_H.into_iter().fold(i32::from(atom.num_H), |sum, value| {
        sum.wrapping_add(i32::from(value))
    });
    if valence != i32::from(atom.chem_bonds_valence).wrapping_add(number_of_hydrogens) {
        return Ok(0);
    }
    if is_hydrogen {
        return Ok(2);
    }

    let neighbor_count = usize::try_from(i32::from(atom.valence).max(0))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for &neighbor in atom
        .neighbor
        .get(..neighbor_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
    {
        let neighbor = atoms
            .get(usize::from(neighbor))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if (neighbor.charge != 0 && atom.charge != 0)
            || (neighbor.radical != 0 && neighbor.radical != RADICAL_SINGLET as S_CHAR)
        {
            return Ok(0);
        }
    }
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn bNumHeterAtomHasIsotopicH(
    atoms: &[inp_ATOM],
    number_of_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4133 bNumHeterAtomHasIsotopicH
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: bNumHeterAtomHasIsotopicH
    // INCHI✔️✔️: int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j, val, is_H = 0, num_H, iat_numb, bAccept, num_iso_H, cur_num_iso_H, num_iso_atoms; /* djb-rwth: removing redundant variables */
    // INCHI✔️✔️:     inp_ATOM *at, *at2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     num_iso_H = 0;
    // INCHI✔️✔️:     num_iso_atoms = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0, at = atom; i < num_atoms; i++, at++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         num_iso_atoms += ( at->iso_atw_diff != 0 || NUM_ISO_H( at, 0 ) );
    // INCHI✔️✔️:         /* isotopic atoms and implicit isotopic H */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (0 >( iat_numb = get_iat_number( at->el_number ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (abs( at->charge ) > 1 || (at->radical && RADICAL_SINGLET != at->radical)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         val = -1;
    // INCHI✔️✔️:         switch (iat_numb)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             case IAT_N:
    // INCHI✔️✔️:             case IAT_P:
    // INCHI✔️✔️:                 /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:                 val = 3 + at->charge;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:
    // INCHI✔️✔️:             case IAT_O:
    // INCHI✔️✔️:             case IAT_S:
    // INCHI✔️✔️:             case IAT_Se:
    // INCHI✔️✔️:             case IAT_Te:
    // INCHI✔️✔️:                 /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:                 val = 2 + at->charge;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:
    // INCHI✔️✔️:             case IAT_F:
    // INCHI✔️✔️:             case IAT_Cl:
    // INCHI✔️✔️:             case IAT_Br:
    // INCHI✔️✔️:             case IAT_I:
    // INCHI✔️✔️:                 if (at->charge == 0)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:                     val = 1;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:
    // INCHI✔️✔️:             case IAT_H:
    // INCHI✔️✔️:                 if (at->valence == 0 &&
    // INCHI✔️✔️:                      at->charge == 1)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     is_H = 1; /* isolated proton */
    // INCHI✔️✔️:                     val = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (val < 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         num_H = NUMH( at, 0 );
    // INCHI✔️✔️:         if (val != at->chem_bonds_valence + num_H)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (is_H)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bAccept = 1;
    // INCHI✔️✔️:             cur_num_iso_H = ( at->iso_atw_diff != 0 );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             cur_num_iso_H = 0;
    // INCHI✔️✔️:             for (j = 0, bAccept = 1; j < at->valence && bAccept; j++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at2 = atom + (int) at->neighbor[j];
    // INCHI✔️✔️:                 if ((at2->charge && at->charge) ||
    // INCHI✔️✔️:                     ( at2->radical && RADICAL_SINGLET != at2->radical )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     bAccept = 0; /* adjacent charged/radical atoms: do not neutralizate */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else if (at2->el_number == EL_NUMBER_H &&
    // INCHI✔️✔️:                           at2->valence == 1 && at2->iso_atw_diff)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     cur_num_iso_H++; /* isotopic explicit H */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:
    // INCHI✔️✔️:             if (bAccept)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_iso_atoms -= cur_num_iso_H;  /* avoid counting explicit H as isotopic atom */
    // INCHI✔️✔️:                 cur_num_iso_H += NUM_ISO_H( at, 0 );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         num_iso_H += ( bAccept && cur_num_iso_H ); /* number of acceptable heteroatoms that have isotopic H */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return
    // INCHI✔️✔️:         ( ( num_iso_H ? 1 : 0 ) | ( num_iso_atoms ? 2 : 0 ) );
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: bNumHeterAtomHasIsotopicH
    // END INCHI C FUNCTION: bNumHeterAtomHasIsotopicH
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bNumHeterAtomHasIsotopicH
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bNumHeterAtomHasIsotopicH
    let atom_count =
        usize::try_from(number_of_atoms.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = atoms
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut is_hydrogen = false;
    let mut number_of_isotopic_hydrogen_atoms = 0_i32;
    let mut number_of_isotopic_atoms = 0_i32;

    for atom in atoms {
        let implicit_isotopic_hydrogens = atom.num_iso_H.into_iter().map(i32::from).sum::<i32>();
        number_of_isotopic_atoms = number_of_isotopic_atoms.wrapping_add(i32::from(
            atom.iso_atw_diff != 0 || implicit_isotopic_hydrogens != 0,
        ));

        let ion_atom_type = get_iat_number(i32::from(atom.el_number));
        if ion_atom_type < 0 {
            continue;
        }
        if i32::from(atom.charge).abs() > 1
            || (atom.radical != 0 && atom.radical != RADICAL_SINGLET as S_CHAR)
        {
            continue;
        }

        let mut valence = -1_i32;
        match ion_atom_type {
            2 | 3 => valence = 3_i32.wrapping_add(i32::from(atom.charge)),
            4..=7 => valence = 2_i32.wrapping_add(i32::from(atom.charge)),
            8..=11 => {
                if atom.charge == 0 {
                    valence = 1;
                }
            }
            0 => {
                if atom.valence == 0 && atom.charge == 1 {
                    is_hydrogen = true;
                    valence = 0;
                }
            }
            _ => {}
        }
        if valence < 0 {
            continue;
        }

        let number_of_hydrogens = i32::from(atom.num_H).wrapping_add(implicit_isotopic_hydrogens);
        if valence != i32::from(atom.chem_bonds_valence).wrapping_add(number_of_hydrogens) {
            continue;
        }

        let mut accept;
        let mut current_isotopic_hydrogens;
        if is_hydrogen {
            accept = true;
            current_isotopic_hydrogens = i32::from(atom.iso_atw_diff != 0);
        } else {
            current_isotopic_hydrogens = 0_i32;
            accept = true;
            let neighbor_count = usize::try_from(i32::from(atom.valence).max(0))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            for &neighbor in atom
                .neighbor
                .get(..neighbor_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                let neighbor = atoms
                    .get(usize::from(neighbor))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if (neighbor.charge != 0 && atom.charge != 0)
                    || (neighbor.radical != 0 && neighbor.radical != RADICAL_SINGLET as S_CHAR)
                {
                    accept = false;
                    break;
                } else if neighbor.el_number == 1
                    && neighbor.valence == 1
                    && neighbor.iso_atw_diff != 0
                {
                    current_isotopic_hydrogens = current_isotopic_hydrogens.wrapping_add(1);
                }
            }
            if accept {
                number_of_isotopic_atoms =
                    number_of_isotopic_atoms.wrapping_sub(current_isotopic_hydrogens);
                current_isotopic_hydrogens =
                    current_isotopic_hydrogens.wrapping_add(implicit_isotopic_hydrogens);
            }
        }
        number_of_isotopic_hydrogen_atoms = number_of_isotopic_hydrogen_atoms
            .wrapping_add(i32::from(accept && current_isotopic_hydrogens != 0));
    }

    Ok(i32::from(number_of_isotopic_hydrogen_atoms != 0)
        | (i32::from(number_of_isotopic_atoms != 0) << 1))
}

pub(crate) fn the_only_doublet_neigh(
    atoms: &[inp_ATOM],
    first_atom: i32,
    first_neighbor_ordinal: &mut i32,
    second_neighbor_ordinal: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:179 the_only_doublet_neigh
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: the_only_doublet_neigh
    // INCHI✔️✔️: int the_only_doublet_neigh( inp_ATOM *at,
    // INCHI✔️✔️:                             int i1,
    // INCHI✔️✔️:                             int *ineigh1,
    // INCHI✔️✔️:                             int *ineigh2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, neigh1, num_rad1 = 0, num_rad2 = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inp_ATOM *a = at + i1, *b;
    // INCHI✔️✔️:     if (RADICAL_DOUBLET != a->radical)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         b = at + ( (int) a->neighbor[i] ); /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:         if (RADICAL_DOUBLET == b->radical)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             num_rad1++;
    // INCHI✔️✔️:             *ineigh1 = i;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (1 == num_rad1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         a = at + ( neigh1 = (int) a->neighbor[*ineigh1] );
    // INCHI✔️✔️:         for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             b = at + (int) a->neighbor[i];
    // INCHI✔️✔️:             if (RADICAL_DOUBLET == b->radical)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 num_rad2++;
    // INCHI✔️✔️:                 *ineigh2 = i;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (1 == num_rad2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return neigh1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return -1;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: the_only_doublet_neigh
    // END INCHI C FUNCTION: the_only_doublet_neigh
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: the_only_doublet_neigh
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define RADICAL_DOUBLET 2
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: the_only_doublet_neigh
    let first_index =
        usize::try_from(first_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let first = atoms
        .get(first_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if first.radical != RADICAL_DOUBLET as S_CHAR {
        return Ok(-1);
    }

    let mut first_doublet_count = 0_i32;
    let first_valence = i32::from(first.valence);
    let mut ordinal = 0_i32;
    while ordinal < first_valence {
        let neighbor_number = *first
            .neighbor
            .get(ordinal as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let neighbor = atoms
            .get(usize::from(neighbor_number))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor.radical == RADICAL_DOUBLET as S_CHAR {
            first_doublet_count = first_doublet_count.wrapping_add(1);
            *first_neighbor_ordinal = ordinal;
        }
        ordinal = ordinal.wrapping_add(1);
    }

    if first_doublet_count == 1 {
        let selected_ordinal = usize::try_from(*first_neighbor_ordinal)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let selected_number = *first
            .neighbor
            .get(selected_ordinal)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let selected = atoms
            .get(usize::from(selected_number))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut second_doublet_count = 0_i32;
        let second_valence = i32::from(selected.valence);
        ordinal = 0;
        while ordinal < second_valence {
            let neighbor_number = *selected
                .neighbor
                .get(ordinal as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor = atoms
                .get(usize::from(neighbor_number))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.radical == RADICAL_DOUBLET as S_CHAR {
                second_doublet_count = second_doublet_count.wrapping_add(1);
                *second_neighbor_ordinal = ordinal;
            }
            ordinal = ordinal.wrapping_add(1);
        }
        if second_doublet_count == 1 {
            return Ok(i32::from(selected_number));
        }
    }
    Ok(-1)
}

#[allow(non_snake_case)]
pub(crate) fn fix_non_uniform_drawn_oxoanions(
    number_of_atoms: i32,
    atoms: &mut [inp_ATOM],
    number_of_changes: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:227 fix_non_uniform_drawn_oxoanions
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: fix_non_uniform_drawn_oxoanions
    // INCHI✔️✔️: int fix_non_uniform_drawn_oxoanions( int num_atoms,
    // INCHI✔️✔️:                                      inp_ATOM *at,
    // INCHI✔️✔️:                                      int *num_changes )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* For central halogen, apply the following
    // INCHI✔️✔️:     correction rules:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     O                             O(-)
    // INCHI✔️✔️:     ||                            |
    // INCHI✔️✔️:     O=Hal(-)=O         ===>        O=Hal=O
    // INCHI✔️✔️:     ||                            ||
    // INCHI✔️✔️:     O                             O
    // INCHI✔️✔️:
    // INCHI✔️✔️:     (perchlorate, etc.)
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     O                             O(-)
    // INCHI✔️✔️:     ||                            |
    // INCHI✔️✔️:     Hal(-)=O         ===>         Hal=O
    // INCHI✔️✔️:     ||                            ||
    // INCHI✔️✔️:     O                             O
    // INCHI✔️✔️:
    // INCHI✔️✔️:     (chlorate, etc.)
    // INCHI✔️✔️:
    // INCHI✔️✔️:     O                             O(-)
    // INCHI✔️✔️:     ||                            |
    // INCHI✔️✔️:     Hal(-)=O         ===>         Hal=O
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     (chlorite, etc.)
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Hal(-)=O         ===>         Hal-O(-)
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     (hypochlorite, etc.)
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     For halcogenes (S, Se, Te)
    // INCHI✔️✔️:     Y                               Y(-)
    // INCHI✔️✔️:     ||                              |
    // INCHI✔️✔️:     RnX(-)            ===>          RnX
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if:
    // INCHI✔️✔️:     1) (X = S, Y = O) || (X = Se, Y = S, O) || (X = Te, Y = O, S, Se)
    // INCHI✔️✔️:     2) valence of X exceeds 6, in initially drawn form
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     So the following is corrected:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     O                               O(-)
    // INCHI✔️✔️:     ||                              |
    // INCHI✔️✔️:     O=S(-)-R            ===>        O=S-R
    // INCHI✔️✔️:     ||                              ||
    // INCHI✔️✔️:     O                               O
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Or
    // INCHI✔️✔️:
    // INCHI✔️✔️:     O                              O-
    // INCHI✔️✔️:     ||                             |
    // INCHI✔️✔️:     F5Te(-)            ===>        F5Te
    // INCHI✔️✔️:
    // INCHI✔️✔️:     but the following remains unchanged:
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     O
    // INCHI✔️✔️:     ||
    // INCHI✔️✔️:     O=S(-)-R
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     The central atom (of IUPAC Group 16-17) is shown as negative but it contains double bond(s)
    // INCHI✔️✔️:     to terminal atom of greater electronegativity (of Group 16).
    // INCHI✔️✔️:     The central atom is halcogen (S,Se,Te) in highest oxidation state or halogen.
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Fix:
    // INCHI✔️✔️:     move negative charge to terminal atom and change double bond to a single one.
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Eligible central atom    Eligible terminal atom at double bond's end
    // INCHI✔️✔️:     Cl              O
    // INCHI✔️✔️:     Br              O
    // INCHI✔️✔️:     I               O
    // INCHI✔️✔️:     [At              S,Se,Te]
    // INCHI✔️✔️:     S               O
    // INCHI✔️✔️:     Se              O,S
    // INCHI✔️✔️:     Te              O, S, Se
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Comments:
    // INCHI✔️✔️:     1. Central atoms of Groups 13-15 are not considered.
    // INCHI✔️✔️:     2. Pauling electronegativities are:
    // INCHI✔️✔️:     F(3.98) > O(3.44) > Cl (3.16) > N (3.04) > Br(2.96) > I(2.66) > S(2.58) > Se(2.55) > At (2.2) > Te(2.1)
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* constants for element numbers */
    // INCHI✔️✔️:     enum elems
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         dNone, dCl = 17, dBr = 35, dI = 53, dAt = 85, dO = 8,
    // INCHI✔️✔️:         dS = 16, dSe = 34, dTe = 52, dP = 15, dC = 6, dN = 7
    // INCHI✔️✔️:     };
    // INCHI✔️✔️:
    // INCHI✔️✔️:     static U_CHAR  allowed_elnums_center_halogen[] = { dCl, dBr, dI, dAt };
    // INCHI✔️✔️:     static U_CHAR  allowed_elnums_center_halcogen[] = { dS, dSe, dTe };
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int en_center;
    // INCHI✔️✔️:     int i, j, k;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Find appropriate central atom. This center should be: ...*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* charged exactly (-1) ... */
    // INCHI✔️✔️:         if (at[i].charge != -1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         en_center = at[i].el_number;
    // INCHI✔️✔️:         /*  from eligible element list ... */
    // INCHI✔️✔️:         if (!memchr( allowed_elnums_center_halogen, en_center, sizeof( allowed_elnums_center_halogen ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* central atom is not not halogen; check if it is halcogen */
    // INCHI✔️✔️:             if (memchr( allowed_elnums_center_halcogen, en_center, sizeof( allowed_elnums_center_halcogen ) ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (at[i].chem_bonds_valence < 7)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* central atom is anionic halcogen, but not in the highest oxidation state */
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* OK, found central halogen or eligible central halcogen. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* non-radical... */
    // INCHI✔️✔️:         if (at[i].radical && ( at[i].radical != RADICAL_SINGLET ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* Center found, now examine the adjacent terminals... */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             int en_term, kk = 0, jj = 0, min_en = 999, iso = 0, min_iso = 999;
    // INCHI✔️✔️:             jj = -1;
    // INCHI✔️✔️:             for (k = 0; k < at[i].valence; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 j = at[i].neighbor[k];
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 /* Terminal should be: ... */
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 /* terminal... */
    // INCHI✔️✔️:                 if (at[j].valence != 1)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /* double-bonded ... */
    // INCHI✔️✔️:                 if (at[i].bond_type[k] != BOND_TYPE_DOUBLE)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /* zero-charged ... */
    // INCHI✔️✔️:                 if (at[j].charge != 0)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /* non-radical */
    // INCHI✔️✔️:                 if (at[j].radical && ( at[j].radical != RADICAL_SINGLET ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /*  of eligible elements list ... */
    // INCHI✔️✔️:                 en_term = at[j].el_number;
    // INCHI✔️✔️:                 switch (en_term)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     case dO:    break;
    // INCHI✔️✔️:                     case dS:    if (( en_center == dSe ) || ( en_center == dAt ) || ( en_center == dTe )) break;  continue;
    // INCHI✔️✔️:                     case dSe:   if (( en_center == dAt ) || ( en_center == dTe )) break;  continue;
    // INCHI✔️✔️:                     case dTe:   if (en_center == dAt) break; continue;
    // INCHI✔️✔️:                     default:    continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 /* From several candidates, select one with less el. number (==more electronegative). */
    // INCHI✔️✔️:                 if (en_term < min_en)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     min_en = en_term; kk = k; jj = j;
    // INCHI✔️✔️:                     min_iso = at[j].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff;
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /* From same-element candidates, select one with less isotopic mass (arbitrary choice). */
    // INCHI✔️✔️:                 else if (en_term == min_en)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     iso = at[j].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff;
    // INCHI✔️✔️:                     if (iso < min_iso)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         min_iso = iso; kk = k; jj = j; continue;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             } /* end of checking nbrs. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:               /* If OK, apply changes. */
    // INCHI✔️✔️:             if (jj >= 0)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[i].charge = 0;
    // INCHI✔️✔️:                 at[jj].charge = -1;
    // INCHI✔️✔️:                 at[i].bond_type[kk] = BOND_TYPE_SINGLE;
    // INCHI✔️✔️:                 at[jj].bond_type[0] = BOND_TYPE_SINGLE;
    // INCHI✔️✔️:                 at[i].bond_stereo[kk] = at[jj].bond_stereo[0] = 0;
    // INCHI✔️✔️:                 at[i].chem_bonds_valence--;
    // INCHI✔️✔️:                 at[jj].chem_bonds_valence--;
    // INCHI✔️✔️:                 ( *num_changes )++;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }  /* end of search for candidate centers. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: fix_non_uniform_drawn_oxoanions
    // END INCHI C FUNCTION: fix_non_uniform_drawn_oxoanions
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: fix_non_uniform_drawn_oxoanions
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define BOND_TYPE_SINGLE 1
    // INCHI✔️✔️: #define BOND_TYPE_DOUBLE 2
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: typedef unsigned char U_CHAR; typedef signed char S_CHAR;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: fix_non_uniform_drawn_oxoanions
    let atom_count =
        usize::try_from(number_of_atoms.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if atoms.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let mut center_index = 0_usize;
    while center_index < atom_count {
        if atoms[center_index].charge != -1 {
            center_index += 1;
            continue;
        }
        let center_element = atoms[center_index].el_number;
        let is_halogen = [17_u8, 35, 53, 85].contains(&center_element);
        if !is_halogen {
            let is_chalcogen = [16_u8, 34, 52].contains(&center_element);
            if !is_chalcogen || atoms[center_index].chem_bonds_valence < 7 {
                center_index += 1;
                continue;
            }
        }
        if atoms[center_index].radical != 0
            && atoms[center_index].radical != RADICAL_SINGLET as S_CHAR
        {
            center_index += 1;
            continue;
        }

        let mut selected_ordinal = 0_usize;
        let mut selected_index: Option<usize> = None;
        let mut minimum_element = 999_i32;
        let mut minimum_isotope = 999_i32;
        let center_valence = i32::from(atoms[center_index].valence);
        let mut ordinal = 0_i32;
        while ordinal < center_valence {
            let ordinal_index = ordinal as usize;
            let neighbor_number = *atoms[center_index]
                .neighbor
                .get(ordinal_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor_index = usize::from(neighbor_number);
            let terminal = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;

            if terminal.valence != 1
                || *atoms[center_index]
                    .bond_type
                    .get(ordinal_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != BOND_TYPE_DOUBLE as u8
                || terminal.charge != 0
                || (terminal.radical != 0 && terminal.radical != RADICAL_SINGLET as S_CHAR)
            {
                ordinal = ordinal.wrapping_add(1);
                continue;
            }

            let terminal_element = terminal.el_number;
            let eligible = match terminal_element {
                8 => true,
                16 => matches!(center_element, 34 | 85 | 52),
                34 => matches!(center_element, 85 | 52),
                52 => center_element == 85,
                _ => false,
            };
            if !eligible {
                ordinal = ordinal.wrapping_add(1);
                continue;
            }

            let terminal_element = i32::from(terminal_element);
            if terminal_element < minimum_element {
                minimum_element = terminal_element;
                selected_ordinal = ordinal_index;
                selected_index = Some(neighbor_index);
                minimum_isotope = if terminal.iso_atw_diff > 0 {
                    i32::from(atoms[center_index].iso_atw_diff).wrapping_sub(1)
                } else {
                    i32::from(atoms[center_index].iso_atw_diff)
                };
            } else if terminal_element == minimum_element {
                let isotope = if terminal.iso_atw_diff > 0 {
                    i32::from(atoms[center_index].iso_atw_diff).wrapping_sub(1)
                } else {
                    i32::from(atoms[center_index].iso_atw_diff)
                };
                if isotope < minimum_isotope {
                    minimum_isotope = isotope;
                    selected_ordinal = ordinal_index;
                    selected_index = Some(neighbor_index);
                }
            }
            ordinal = ordinal.wrapping_add(1);
        }

        if let Some(terminal_index) = selected_index {
            atoms[center_index].charge = 0;
            atoms[terminal_index].charge = -1;
            atoms[center_index].bond_type[selected_ordinal] = BOND_TYPE_SINGLE as u8;
            atoms[terminal_index].bond_type[0] = BOND_TYPE_SINGLE as u8;
            atoms[terminal_index].bond_stereo[0] = 0;
            atoms[center_index].bond_stereo[selected_ordinal] = 0;
            atoms[center_index].chem_bonds_valence =
                atoms[center_index].chem_bonds_valence.wrapping_sub(1);
            atoms[terminal_index].chem_bonds_valence =
                atoms[terminal_index].chem_bonds_valence.wrapping_sub(1);
            *number_of_changes = number_of_changes.wrapping_add(1);
        }
        center_index += 1;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn fix_non_uniform_drawn_amidiniums(
    number_of_atoms: i32,
    atoms: &mut [inp_ATOM],
    number_of_changes: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:452 fix_non_uniform_drawn_amidiniums
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: fix_non_uniform_drawn_amidiniums
    // INCHI✔️✔️: int fix_non_uniform_drawn_amidiniums( int num_atoms,
    // INCHI✔️✔️:                                       inp_ATOM *at,
    // INCHI✔️✔️:                                       int *num_changes )
    // INCHI✔️✔️:
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* Amidines include carboxamidines RC(=NR)NR2,
    // INCHI✔️✔️:     sulfinamidines RS(=NR)NR2 and phosphinamidines, R2P(=NR)NR2.
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     NR                              NR
    // INCHI✔️✔️:     |                               |
    // INCHI✔️✔️:     R"-Y-NHR'         ===>          R"-Y=N(+)HR'
    // INCHI✔️✔️:     (+)
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Y = C, S, P
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Fix:
    // INCHI✔️✔️:     move positive charge to nitrogen and change single bond to a double one.
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     Comment:
    // INCHI✔️✔️:     Fix is applied only if at least one of R's at N is hydrogen
    // INCHI✔️✔️:     (otherwise we have just a '+' delocalization which is already recognized).
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* constants for element numbers */
    // INCHI✔️✔️:     enum elems
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         dNone, dCl = 17, dBr = 35, dI = 53, dAt = 85, dO = 8,
    // INCHI✔️✔️:         dS = 16, dSe = 34, dTe = 52, dP = 15, dC = 6, dN = 7
    // INCHI✔️✔️:     };
    // INCHI✔️✔️:
    // INCHI✔️✔️:     static U_CHAR  allowed_elnums_center[] = { dC, dS, dP };
    // INCHI✔️✔️:     int en_center;
    // INCHI✔️✔️:     int i, j, k, jj, kk, k1;
    // INCHI✔️✔️:     int mismatch = 0, nuH = 0, nuN = 0, nitrogens[MAXVAL];
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Find appropriate central atom. This center should be: ...*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* charged exactly (+1) ... */
    // INCHI✔️✔️:         if (at[i].charge != 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         en_center = at[i].el_number;
    // INCHI✔️✔️:         /*  from eligible element list ... */
    // INCHI✔️✔️:         if (!memchr( allowed_elnums_center, en_center, sizeof( allowed_elnums_center ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* has exactly 3 neighbours connected by single bonds*/
    // INCHI✔️✔️:         if (at[i].valence != 3)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (at[i].chem_bonds_valence != 3)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* non-radical. */
    // INCHI✔️✔️:         if (at[i].radical && ( at[i].radical != RADICAL_SINGLET ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* NB: center must have neutral neighbours, two of them are aliphatic N's of which at least one bears H. */
    // INCHI✔️✔️:         mismatch = nuH = nuN = kk = 0; /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:         memset( nitrogens, 0, sizeof( nitrogens ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         jj = -1;
    // INCHI✔️✔️:         for (k = 0; k < at[i].valence; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             j = at[i].neighbor[k];
    // INCHI✔️✔️:
    // INCHI✔️✔️:             if (at[j].charge != 0)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 mismatch = 1;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (at[j].el_number == dN)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (( at[j].valence > 3 ) || ( at[j].chem_bonds_valence > 3 ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     mismatch = 1;
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 nuH += NUMH( at, j );
    // INCHI✔️✔️:                 nuN++;
    // INCHI✔️✔️:                 if (jj < 0)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     jj = j;
    // INCHI✔️✔️:                     kk = k;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* If OK, apply changes. */
    // INCHI✔️✔️:         if (mismatch)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nuN != 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nuH < 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (jj >= 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             at[i].charge = 0;
    // INCHI✔️✔️:             at[jj].charge = 1;
    // INCHI✔️✔️:             at[i].bond_type[kk] = BOND_TYPE_DOUBLE;
    // INCHI✔️✔️:             for (k1 = 0; k1 < at[jj].valence && i != at[jj].neighbor[k1]; k1++)
    // INCHI✔️✔️:                 ;
    // INCHI✔️✔️:             at[jj].bond_type[k1] = BOND_TYPE_DOUBLE;
    // INCHI✔️✔️:             at[i].chem_bonds_valence++;
    // INCHI✔️✔️:             at[jj].chem_bonds_valence++;
    // INCHI✔️✔️:             /* NB: do nothing with wedge stereo bonds (retain wedge) */
    // INCHI✔️✔️:
    // INCHI✔️✔️:             ( *num_changes )++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }  /* end of search for candidate centers. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: fix_non_uniform_drawn_amidiniums
    // END INCHI C FUNCTION: fix_non_uniform_drawn_amidiniums
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: fix_non_uniform_drawn_amidiniums
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H+NUM_ISO_H(AT,N))
    // INCHI✔️✔️: #define BOND_TYPE_DOUBLE 2
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define MAXVAL 20
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: fix_non_uniform_drawn_amidiniums
    let atom_count =
        usize::try_from(number_of_atoms.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if atoms.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut nitrogens = [0_i32; MAXVAL as usize];

    let mut center_index = 0_usize;
    while center_index < atom_count {
        if atoms[center_index].charge != 1
            || ![6_u8, 16, 15].contains(&atoms[center_index].el_number)
            || atoms[center_index].valence != 3
            || atoms[center_index].chem_bonds_valence != 3
            || (atoms[center_index].radical != 0
                && atoms[center_index].radical != RADICAL_SINGLET as S_CHAR)
        {
            center_index += 1;
            continue;
        }

        let mut mismatch = false;
        let mut number_of_hydrogens = 0_i32;
        let mut number_of_nitrogens = 0_i32;
        let mut first_nitrogen_ordinal = 0_usize;
        nitrogens.fill(0);
        let mut first_nitrogen: Option<usize> = None;
        let center_valence = i32::from(atoms[center_index].valence);
        let mut ordinal = 0_i32;
        while ordinal < center_valence {
            let neighbor_number = *atoms[center_index]
                .neighbor
                .get(ordinal as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor_index = usize::from(neighbor_number);
            let neighbor = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.charge != 0 {
                mismatch = true;
                break;
            }
            if neighbor.el_number == 7 {
                if neighbor.valence > 3 || neighbor.chem_bonds_valence > 3 {
                    mismatch = true;
                    break;
                }
                let neighbor_hydrogens = i32::from(neighbor.num_H)
                    .wrapping_add(neighbor.num_iso_H.into_iter().map(i32::from).sum::<i32>());
                number_of_hydrogens = number_of_hydrogens.wrapping_add(neighbor_hydrogens);
                number_of_nitrogens = number_of_nitrogens.wrapping_add(1);
                if first_nitrogen.is_none() {
                    first_nitrogen = Some(neighbor_index);
                    first_nitrogen_ordinal = ordinal as usize;
                }
            }
            ordinal = ordinal.wrapping_add(1);
        }

        if !mismatch
            && number_of_nitrogens == 2
            && number_of_hydrogens >= 1
            && let Some(nitrogen_index) = first_nitrogen
        {
            atoms[center_index].charge = 0;
            atoms[nitrogen_index].charge = 1;
            atoms[center_index].bond_type[first_nitrogen_ordinal] = BOND_TYPE_DOUBLE as u8;

            let nitrogen_valence = i32::from(atoms[nitrogen_index].valence);
            let mut reciprocal_ordinal = 0_i32;
            while reciprocal_ordinal < nitrogen_valence {
                let neighbor_number = *atoms[nitrogen_index]
                    .neighbor
                    .get(reciprocal_ordinal as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if neighbor_number == center_index as AT_NUMB {
                    break;
                }
                reciprocal_ordinal = reciprocal_ordinal.wrapping_add(1);
            }
            *atoms[nitrogen_index]
                .bond_type
                .get_mut(reciprocal_ordinal as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = BOND_TYPE_DOUBLE as u8;
            atoms[center_index].chem_bonds_valence =
                atoms[center_index].chem_bonds_valence.wrapping_add(1);
            atoms[nitrogen_index].chem_bonds_valence =
                atoms[nitrogen_index].chem_bonds_valence.wrapping_add(1);
            *number_of_changes = number_of_changes.wrapping_add(1);
        }
        center_index += 1;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn fix_odd_things(
    number_of_atoms: i32,
    atoms: &mut [inp_ATOM],
    fix_bug: i32,
    fix_non_uniform_draw: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:603 fix_odd_things
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: fix_odd_things
    // INCHI✔️✔️: int fix_odd_things( int num_atoms,
    // INCHI✔️✔️:                     inp_ATOM *at,
    // INCHI✔️✔️:                     int bFixBug,
    // INCHI✔️✔️:                     int bFixNonUniformDraw )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* N;P;As;Sb;O;S;Se;Te;C;Si */
    // INCHI✔️✔️:     static const U_CHAR en[] = {
    // INCHI✔️✔️:         EL_NUMBER_N,
    // INCHI✔️✔️:         EL_NUMBER_P,
    // INCHI✔️✔️:         EL_NUMBER_AS,
    // INCHI✔️✔️:         EL_NUMBER_SB,
    // INCHI✔️✔️:         EL_NUMBER_O,
    // INCHI✔️✔️:         EL_NUMBER_S,
    // INCHI✔️✔️:         EL_NUMBER_SE,
    // INCHI✔️✔️:         EL_NUMBER_TE
    // INCHI✔️✔️:     };
    // INCHI✔️✔️:     static int ne = sizeof(en)/sizeof(en[0]);
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define FIRST_NEIGHB2  4
    // INCHI✔️✔️: #define FIRST_CENTER2  5
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int i1, i2, k1, k2, c = -1, num_changes = 0;
    // INCHI✔️✔️:     /* djb-rwth: removing redundant variables */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bFixNonUniformDraw)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int ret1; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️✔️:         ret1 = fix_non_uniform_drawn_oxoanions( num_atoms, at, &num_changes ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️✔️:         ret1 = fix_non_uniform_drawn_amidiniums( num_atoms, at, &num_changes ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* H(-)-X  -> H-X(-);  H(+)-X  -> H-X(+) */
    // INCHI✔️✔️:     for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (1 == at[i1].valence &&
    // INCHI✔️✔️:              1 == abs( at[i1].charge ) &&
    // INCHI✔️✔️:              ( 0 == at[i1].radical || RADICAL_SINGLET == at[i1].radical ) &&
    // INCHI✔️✔️:              BOND_TYPE_SINGLE == at[i1].bond_type[0] &&
    // INCHI✔️✔️:              EL_NUMBER_H == at[i1].el_number && EL_NUMBER_H != at[i2 = (int) at[i1].neighbor[0]].el_number &&
    // INCHI✔️✔️:              !NUMH( at, i1 ) && !NUMH( at, i2 ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             at[i2].charge += at[i1].charge;
    // INCHI✔️✔️:             at[i1].charge = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* replace XHm(-)--Y==XHn(+) with XHm==Y--XHn, (n>=0 ,m>=0, X=N,P,As,Sb,O,S,Se,Te) */
    // INCHI✔️✔️:     for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (1 != at[i1].charge ||
    // INCHI✔️✔️:              (at[i1].radical && RADICAL_SINGLET != at[i1].radical) ||
    // INCHI✔️✔️:              at[i1].chem_bonds_valence == at[i1].valence ||
    // INCHI✔️✔️:              !memchr( en, at[i1].el_number, ne ) ||
    // INCHI✔️✔️:              get_el_valence( at[i1].el_number, at[i1].charge, 0 ) != at[i1].chem_bonds_valence + NUMH( at, i1 )) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* found a candidate at[i1] for X in XHn(+) */
    // INCHI✔️✔️:         if (1 == at[i1].valence &&
    // INCHI✔️✔️:              BOND_TYPE_DOUBLE == at[i1].bond_type[0])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             c = (int) at[i1].neighbor[0];
    // INCHI✔️✔️:             for (k2 = 0; k2 < at[c].valence; k2++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 i2 = at[c].neighbor[k2];
    // INCHI✔️✔️:                 if (1 == at[i2].valence &&
    // INCHI✔️✔️:                      -1 == at[i2].charge  &&
    // INCHI✔️✔️:                      at[i2].el_number == at[i1].el_number && /* exact match */
    // INCHI✔️✔️:                      ( 0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical ) &&
    // INCHI✔️✔️:                      BOND_TYPE_SINGLE == at[i2].bond_type[0] &&
    // INCHI✔️✔️:                      /*memchr(en, at[i2].el_number, ne) &&*/
    // INCHI✔️✔️:                      get_el_valence( at[i2].el_number, at[i2].charge, 0 ) == at[i2].chem_bonds_valence + NUMH( at, i2 ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* found both X(-) and X(+); change bonds and remove charges */
    // INCHI✔️✔️:                     for (k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1++)
    // INCHI✔️✔️:                         ;
    // INCHI✔️✔️:                     at[i1].charge = at[i2].charge = 0;
    // INCHI✔️✔️:                     at[i1].bond_type[0] = at[c].bond_type[k1] = BOND_TYPE_SINGLE;
    // INCHI✔️✔️:                     at[i1].chem_bonds_valence--;
    // INCHI✔️✔️:                     at[i2].bond_type[0] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
    // INCHI✔️✔️:                     at[i2].chem_bonds_valence++;
    // INCHI✔️✔️:                     num_changes++;
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* explicit H case: detect H-neighbors and Y */
    // INCHI✔️✔️:             int ineigh, neigh, i1_c, i2_c, num_H_i1, num_H_i2;
    // INCHI✔️✔️:             for (ineigh = 0, num_H_i1 = 0, i1_c = -1; ineigh < at[i1].valence; ineigh++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 neigh = at[i1].neighbor[ineigh];
    // INCHI✔️✔️:                 if (at[neigh].el_number == EL_NUMBER_H)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (at[neigh].chem_bonds_valence == 1 &&
    // INCHI✔️✔️:                         ( 0 == at[neigh].radical || RADICAL_SINGLET == at[neigh].radical ))
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         num_H_i1++; /* found H-neighbor */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         break;  /* wrong neighbor */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else if (at[i1].bond_type[ineigh] == BOND_TYPE_DOUBLE)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* found a candidate for Y; bond must be double */
    // INCHI✔️✔️:                     i1_c = ineigh;
    // INCHI✔️✔️:                     c = neigh;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (i1_c < 0 || num_H_i1 + 1 != at[i1].valence)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             for (k2 = 0; k2 < at[c].valence; k2++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 i2 = at[c].neighbor[k2];
    // INCHI✔️✔️:                 if (-1 == at[i2].charge  &&
    // INCHI✔️✔️:                      at[i2].el_number == at[i1].el_number && /* exact match */
    // INCHI✔️✔️:                      ( 0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical ) &&
    // INCHI✔️✔️:                      get_el_valence( at[i2].el_number, at[i2].charge, 0 ) == at[i2].chem_bonds_valence + NUMH( at, i2 ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     for (ineigh = 0, num_H_i2 = 0, i2_c = -1; ineigh < at[i2].valence; ineigh++)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         neigh = at[i2].neighbor[ineigh];
    // INCHI✔️✔️:                         if (at[neigh].el_number == EL_NUMBER_H)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (at[neigh].chem_bonds_valence == 1 &&
    // INCHI✔️✔️:                                 ( 0 == at[neigh].radical || RADICAL_SINGLET == at[neigh].radical ))
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 num_H_i2++;  /* found H-neighbor */
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             else
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 break; /* wrong neighbor */
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         else
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (c == neigh && at[i2].bond_type[ineigh] == BOND_TYPE_SINGLE)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 i2_c = ineigh; /* position of Y neighbor; bond must be single */
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             else
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     if (num_H_i2 + ( i2_c >= 0 ) != at[i2].valence)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         continue;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     /* found both X(-) and X(+); change bonds and remove charges */
    // INCHI✔️✔️:                     for (k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1++)
    // INCHI✔️✔️:                         ;
    // INCHI✔️✔️:                     if ((i1_c >= 0) && (i2_c >= 0)) /* djb-rwth: fixing coverity ID #499537 */
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         at[i1].charge = at[i2].charge = 0;
    // INCHI✔️✔️:                         at[i1].bond_type[i1_c] = at[c].bond_type[k1] = BOND_TYPE_SINGLE;
    // INCHI✔️✔️:                         at[i1].chem_bonds_valence--;
    // INCHI✔️✔️:                         at[i2].bond_type[i2_c] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
    // INCHI✔️✔️:                         at[i2].chem_bonds_valence++;
    // INCHI✔️✔️:                         num_changes++;
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         continue;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             } /* k2 */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Replace
    // INCHI✔️✔️:
    // INCHI✔️✔️:     X-                X        X=O,S,Se,Te -- terminal atoms  (NEIGHB2)
    // INCHI✔️✔️:     \ |               \ ||
    // INCHI✔️✔️:     >Y++    with      >Y        Y=S,Se,Te   -- central cation  (CENTER2)
    // INCHI✔️✔️:     / |               / ||
    // INCHI✔️✔️:     X-                X        Y valence=4, original Y bond valence = 4
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     --- the following case of P is processed separately in remove_ion_pairs()
    // INCHI✔️✔️:     --- therefire, it has been disabled here, see #ifndef FIX_P_IV_Plus_O_Minus -- 2010-03-17 DT
    // INCHI✔️✔️:
    // INCHI✔️✔️:     X-                X        X=O,S,Se,Te -- terminal atoms  (NEIGHB2)
    // INCHI✔️✔️:     \ |               \ ||
    // INCHI✔️✔️:     >P+     with      >P
    // INCHI✔️✔️:     / |               / |
    // INCHI✔️✔️:     X-                X-       Y valence=4, original Y bond valence = 4
    // INCHI✔️✔️:
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (1 == at[i1].valence &&
    // INCHI✔️✔️:              -1 == at[i1].charge &&
    // INCHI✔️✔️:              ( 0 == at[i1].radical || RADICAL_SINGLET == at[i1].radical ) &&
    // INCHI✔️✔️:              !NUMH( at, i1 ) &&
    // INCHI✔️✔️:              BOND_TYPE_SINGLE == at[i1].bond_type[0] &&
    // INCHI✔️✔️:              memchr( en + FIRST_NEIGHB2, at[i1].el_number, (long long)ne - FIRST_NEIGHB2 )) /* djb-rwth: cast operator added */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             int charge, i;
    // INCHI✔️✔️:             /* found a candidate for X */
    // INCHI✔️✔️:             c = (int) at[i1].neighbor[0]; /* candidate for Y */
    // INCHI✔️✔️:             if (( ( charge = 2 ) == at[c].charge && memchr( en + FIRST_CENTER2, at[c].el_number, (long long)ne - FIRST_CENTER2 ) /* djb-rwth: cast operator added */
    // INCHI✔️✔️:
    // INCHI✔️✔️: #ifndef FIX_P_IV_Plus_O_Minus
    // INCHI✔️✔️:                   || ( charge = 1 ) == at[c].charge && EL_NUMBER_P == at[c].el_number
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:                   ) &&
    // INCHI✔️✔️:                  4 == at[c].valence &&
    // INCHI✔️✔️:                  ( 0 == at[c].radical || RADICAL_SINGLET == at[c].radical ) &&
    // INCHI✔️✔️:                  at[c].valence == at[c].chem_bonds_valence &&
    // INCHI✔️✔️:                  !NUMH( at, c ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ;  /* accept */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 continue; /* ignore at[i1] */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             for (k2 = 0; k2 < at[c].valence; k2++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 i2 = at[c].neighbor[k2];
    // INCHI✔️✔️:                 if (i2 == i1)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (1 == at[i2].valence &&
    // INCHI✔️✔️:                      -1 == at[i2].charge  &&
    // INCHI✔️✔️:                      memchr( en + FIRST_NEIGHB2, at[i2].el_number, (long long)ne - FIRST_NEIGHB2 ) && /* djb-rwth: cast operator added */
    // INCHI✔️✔️:                      /*at[i2].el_number == at[i1].el_number &&*/ /* exact match */
    // INCHI✔️✔️:                      ( 0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical ) &&
    // INCHI✔️✔️:                      !NUMH( at, i2 ) &&
    // INCHI✔️✔️:                      BOND_TYPE_SINGLE == at[i2].bond_type[0])
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* found both X(-) and X(-); change bonds and remove charges */
    // INCHI✔️✔️:                     for (k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1++)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         ;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     for (i = 0; i < charge; i++)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* in case of P it does not matter which X atom is neutralized
    // INCHI✔️✔️:                         because of tautomerism. However, neutral central atom is important
    // INCHI✔️✔️:                         for the neutralization of the components */
    // INCHI✔️✔️:                         switch (i)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 at[i1].charge++; /* = 0; changed 2010-03-17 DT*/
    // INCHI✔️✔️:                                 at[i1].bond_type[0] = at[c].bond_type[k1] = BOND_TYPE_DOUBLE;
    // INCHI✔️✔️:                                 at[i1].bond_stereo[0] = at[c].bond_stereo[k1] = 0;
    // INCHI✔️✔️:                                 at[i1].chem_bonds_valence++;
    // INCHI✔️✔️:                                 at[c].chem_bonds_valence++;
    // INCHI✔️✔️:                                 if (bFixBug) at[c].charge--; /* added 2010-03-17 DT*/
    // INCHI✔️✔️:                                 num_changes++;
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case 1:
    // INCHI✔️✔️:                                 at[i2].charge++; /*= 0; changed 2010-03-17 DT*/
    // INCHI✔️✔️:                                 at[i2].bond_type[0] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
    // INCHI✔️✔️:                                 at[i2].bond_stereo[0] = at[c].bond_stereo[k2] = 0;
    // INCHI✔️✔️:                                 at[i2].chem_bonds_valence++;
    // INCHI✔️✔️:                                 at[c].chem_bonds_valence++;
    // INCHI✔️✔️:                                 if (bFixBug) at[c].charge--; /* added 2010-03-17 DT */
    // INCHI✔️✔️:                                 num_changes++;
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                     /*   -- removed -- 2010-03-17 DT
    // INCHI✔️✔️:                     #if ( FIX_ODD_THINGS_REM_Plus_BUG == 1 )
    // INCHI✔️✔️:                     at[c].charge -= charge;
    // INCHI✔️✔️:                     #else
    // INCHI✔️✔️:                     if ( bFixBug )
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                     at[c].charge -= charge;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     #endif
    // INCHI✔️✔️:                     */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* A(doublet)-B(doublet) -> A=B  (A and B have no other doublet neighbors) */
    // INCHI✔️✔️:     /* A(doublet)=B(doublet) -> A#B  (A and B have no other doublet neighbors) */
    // INCHI✔️✔️:     for (i1 = 0; i1 < num_atoms; i1++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (RADICAL_DOUBLET == at[i1].radical &&
    // INCHI✔️✔️:              0 <= ( i2 = the_only_doublet_neigh( at, i1, &k1, &k2 ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at[i1].bond_type[k1] <= BOND_TYPE_DOUBLE)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 at[i1].bond_type[k1] ++;
    // INCHI✔️✔️:                 at[i1].chem_bonds_valence++;
    // INCHI✔️✔️:                 at[i2].bond_type[k2] ++;
    // INCHI✔️✔️:                 at[i2].chem_bonds_valence++;
    // INCHI✔️✔️:                 at[i1].radical = 0;
    // INCHI✔️✔️:                 at[i2].radical = 0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: #if ( REMOVE_ION_PAIRS_EARLY == 1 )
    // INCHI✔️✔️:     num_changes += remove_ion_pairs( num_atoms, at );
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_changes;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: fix_odd_things
    // END INCHI C FUNCTION: fix_odd_things
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: fix_odd_things
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define FIX_P_IV_Plus_O_Minus
    // INCHI✔️✔️: #define REMOVE_ION_PAIRS_EARLY 1
    // INCHI✔️✔️: #define FIRST_NEIGHB2 4
    // INCHI✔️✔️: #define FIRST_CENTER2 5
    // INCHI✔️✔️: The #ifndef FIX_P_IV_Plus_O_Minus phosphorus branch is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: fix_odd_things
    let atom_count =
        usize::try_from(number_of_atoms.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if atoms.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    const ELIGIBLE_ELEMENTS: [u8; 8] = [7, 15, 33, 51, 8, 16, 34, 52];
    let num_h = |atom: &inp_ATOM| {
        i32::from(atom.num_H).wrapping_add(atom.num_iso_H.into_iter().map(i32::from).sum::<i32>())
    };
    let mut number_of_changes = 0_i32;

    if fix_non_uniform_draw != 0 {
        let _ = fix_non_uniform_drawn_oxoanions(number_of_atoms, atoms, &mut number_of_changes)?;
        let _ = fix_non_uniform_drawn_amidiniums(number_of_atoms, atoms, &mut number_of_changes)?;
    }

    for first_index in 0..atom_count {
        if atoms[first_index].valence == 1
            && i32::from(atoms[first_index].charge).abs() == 1
            && (atoms[first_index].radical == 0
                || atoms[first_index].radical == RADICAL_SINGLET as S_CHAR)
            && atoms[first_index].bond_type[0] == BOND_TYPE_SINGLE as u8
            && atoms[first_index].el_number == 1
        {
            let second_index = usize::from(atoms[first_index].neighbor[0]);
            let second = atoms
                .get(second_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if second.el_number != 1 && num_h(&atoms[first_index]) == 0 && num_h(second) == 0 {
                atoms[second_index].charge = atoms[second_index]
                    .charge
                    .wrapping_add(atoms[first_index].charge);
                atoms[first_index].charge = 0;
            }
        }
    }

    for first_index in 0..atom_count {
        if atoms[first_index].charge != 1
            || (atoms[first_index].radical != 0
                && atoms[first_index].radical != RADICAL_SINGLET as S_CHAR)
            || atoms[first_index].chem_bonds_valence == atoms[first_index].valence
            || !ELIGIBLE_ELEMENTS.contains(&atoms[first_index].el_number)
            || get_el_valence(
                i32::from(atoms[first_index].el_number),
                i32::from(atoms[first_index].charge),
                0,
            )? != i32::from(atoms[first_index].chem_bonds_valence)
                .wrapping_add(num_h(&atoms[first_index]))
        {
            continue;
        }

        if atoms[first_index].valence == 1
            && atoms[first_index].bond_type[0] == BOND_TYPE_DOUBLE as u8
        {
            let center_index = usize::from(atoms[first_index].neighbor[0]);
            let center_valence = i32::from(
                atoms
                    .get(center_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .valence,
            );
            let mut second_ordinal = 0_i32;
            while second_ordinal < center_valence {
                let second_index = usize::from(
                    *atoms[center_index]
                        .neighbor
                        .get(second_ordinal as usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let second = atoms
                    .get(second_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if second.valence == 1
                    && second.charge == -1
                    && second.el_number == atoms[first_index].el_number
                    && (second.radical == 0 || second.radical == RADICAL_SINGLET as S_CHAR)
                    && second.bond_type[0] == BOND_TYPE_SINGLE as u8
                    && get_el_valence(i32::from(second.el_number), i32::from(second.charge), 0)?
                        == i32::from(second.chem_bonds_valence).wrapping_add(num_h(second))
                {
                    let mut first_ordinal = 0_i32;
                    while first_ordinal < center_valence
                        && atoms[center_index].neighbor[first_ordinal as usize]
                            != first_index as AT_NUMB
                    {
                        first_ordinal = first_ordinal.wrapping_add(1);
                    }
                    atoms[second_index].charge = 0;
                    atoms[first_index].charge = 0;
                    atoms[first_index].bond_type[0] = BOND_TYPE_SINGLE as u8;
                    atoms[center_index].bond_type[first_ordinal as usize] = BOND_TYPE_SINGLE as u8;
                    atoms[first_index].chem_bonds_valence =
                        atoms[first_index].chem_bonds_valence.wrapping_sub(1);
                    atoms[second_index].bond_type[0] = BOND_TYPE_DOUBLE as u8;
                    atoms[center_index].bond_type[second_ordinal as usize] = BOND_TYPE_DOUBLE as u8;
                    atoms[second_index].chem_bonds_valence =
                        atoms[second_index].chem_bonds_valence.wrapping_add(1);
                    number_of_changes = number_of_changes.wrapping_add(1);
                    break;
                }
                second_ordinal = second_ordinal.wrapping_add(1);
            }
        } else {
            let first_valence = i32::from(atoms[first_index].valence);
            let mut first_center_ordinal = -1_i32;
            let mut center_index = 0_usize;
            let mut first_hydrogens = 0_i32;
            let mut ordinal = 0_i32;
            while ordinal < first_valence {
                let neighbor_index = usize::from(
                    *atoms[first_index]
                        .neighbor
                        .get(ordinal as usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let neighbor = atoms
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if neighbor.el_number == 1 {
                    if neighbor.chem_bonds_valence == 1
                        && (neighbor.radical == 0 || neighbor.radical == RADICAL_SINGLET as S_CHAR)
                    {
                        first_hydrogens = first_hydrogens.wrapping_add(1);
                    } else {
                        break;
                    }
                } else if atoms[first_index].bond_type[ordinal as usize] == BOND_TYPE_DOUBLE as u8 {
                    first_center_ordinal = ordinal;
                    center_index = neighbor_index;
                }
                ordinal = ordinal.wrapping_add(1);
            }
            if first_center_ordinal < 0 || first_hydrogens.wrapping_add(1) != first_valence {
                continue;
            }

            let center_valence = i32::from(atoms[center_index].valence);
            let mut second_ordinal = 0_i32;
            while second_ordinal < center_valence {
                let second_index = usize::from(
                    *atoms[center_index]
                        .neighbor
                        .get(second_ordinal as usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let second = atoms
                    .get(second_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if second.charge == -1
                    && second.el_number == atoms[first_index].el_number
                    && (second.radical == 0 || second.radical == RADICAL_SINGLET as S_CHAR)
                    && get_el_valence(i32::from(second.el_number), i32::from(second.charge), 0)?
                        == i32::from(second.chem_bonds_valence).wrapping_add(num_h(second))
                {
                    let second_valence = i32::from(second.valence);
                    let mut second_hydrogens = 0_i32;
                    let mut second_center_ordinal = -1_i32;
                    let mut neighbor_ordinal = 0_i32;
                    while neighbor_ordinal < second_valence {
                        let neighbor_index = usize::from(
                            *atoms[second_index]
                                .neighbor
                                .get(neighbor_ordinal as usize)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        let neighbor = atoms
                            .get(neighbor_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if neighbor.el_number == 1 {
                            if neighbor.chem_bonds_valence == 1
                                && (neighbor.radical == 0
                                    || neighbor.radical == RADICAL_SINGLET as S_CHAR)
                            {
                                second_hydrogens = second_hydrogens.wrapping_add(1);
                            } else {
                                break;
                            }
                        } else if center_index == neighbor_index
                            && atoms[second_index].bond_type[neighbor_ordinal as usize]
                                == BOND_TYPE_SINGLE as u8
                        {
                            second_center_ordinal = neighbor_ordinal;
                        } else {
                            break;
                        }
                        neighbor_ordinal = neighbor_ordinal.wrapping_add(1);
                    }
                    if second_hydrogens.wrapping_add(i32::from(second_center_ordinal >= 0))
                        != second_valence
                    {
                        second_ordinal = second_ordinal.wrapping_add(1);
                        continue;
                    }

                    let mut center_first_ordinal = 0_i32;
                    while center_first_ordinal < center_valence
                        && atoms[center_index].neighbor[center_first_ordinal as usize]
                            != first_index as AT_NUMB
                    {
                        center_first_ordinal = center_first_ordinal.wrapping_add(1);
                    }
                    if second_center_ordinal >= 0 {
                        atoms[second_index].charge = 0;
                        atoms[first_index].charge = 0;
                        atoms[first_index].bond_type[first_center_ordinal as usize] =
                            BOND_TYPE_SINGLE as u8;
                        atoms[center_index].bond_type[center_first_ordinal as usize] =
                            BOND_TYPE_SINGLE as u8;
                        atoms[first_index].chem_bonds_valence =
                            atoms[first_index].chem_bonds_valence.wrapping_sub(1);
                        atoms[second_index].bond_type[second_center_ordinal as usize] =
                            BOND_TYPE_DOUBLE as u8;
                        atoms[center_index].bond_type[second_ordinal as usize] =
                            BOND_TYPE_DOUBLE as u8;
                        atoms[second_index].chem_bonds_valence =
                            atoms[second_index].chem_bonds_valence.wrapping_add(1);
                        number_of_changes = number_of_changes.wrapping_add(1);
                        break;
                    }
                }
                second_ordinal = second_ordinal.wrapping_add(1);
            }
        }
    }

    for first_index in 0..atom_count {
        if atoms[first_index].valence != 1
            || atoms[first_index].charge != -1
            || (atoms[first_index].radical != 0
                && atoms[first_index].radical != RADICAL_SINGLET as S_CHAR)
            || num_h(&atoms[first_index]) != 0
            || atoms[first_index].bond_type[0] != BOND_TYPE_SINGLE as u8
            || ![8_u8, 16, 34, 52].contains(&atoms[first_index].el_number)
        {
            continue;
        }
        let center_index = usize::from(atoms[first_index].neighbor[0]);
        let center = atoms
            .get(center_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let charge = 2_i32;
        if center.charge != charge as S_CHAR
            || ![16_u8, 34, 52].contains(&center.el_number)
            || center.valence != 4
            || (center.radical != 0 && center.radical != RADICAL_SINGLET as S_CHAR)
            || center.valence != center.chem_bonds_valence
            || num_h(center) != 0
        {
            continue;
        }
        let center_valence = i32::from(center.valence);
        let mut second_ordinal = 0_i32;
        while second_ordinal < center_valence {
            let second_index = usize::from(
                *atoms[center_index]
                    .neighbor
                    .get(second_ordinal as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            if second_index == first_index {
                second_ordinal = second_ordinal.wrapping_add(1);
                continue;
            }
            let second = atoms
                .get(second_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if second.valence == 1
                && second.charge == -1
                && [8_u8, 16, 34, 52].contains(&second.el_number)
                && (second.radical == 0 || second.radical == RADICAL_SINGLET as S_CHAR)
                && num_h(second) == 0
                && second.bond_type[0] == BOND_TYPE_SINGLE as u8
            {
                let mut first_ordinal = 0_i32;
                while first_ordinal < center_valence
                    && atoms[center_index].neighbor[first_ordinal as usize]
                        != first_index as AT_NUMB
                {
                    first_ordinal = first_ordinal.wrapping_add(1);
                }
                atoms[first_index].charge = atoms[first_index].charge.wrapping_add(1);
                atoms[first_index].bond_type[0] = BOND_TYPE_DOUBLE as u8;
                atoms[center_index].bond_type[first_ordinal as usize] = BOND_TYPE_DOUBLE as u8;
                atoms[first_index].bond_stereo[0] = 0;
                atoms[center_index].bond_stereo[first_ordinal as usize] = 0;
                atoms[first_index].chem_bonds_valence =
                    atoms[first_index].chem_bonds_valence.wrapping_add(1);
                atoms[center_index].chem_bonds_valence =
                    atoms[center_index].chem_bonds_valence.wrapping_add(1);
                if fix_bug != 0 {
                    atoms[center_index].charge = atoms[center_index].charge.wrapping_sub(1);
                }
                number_of_changes = number_of_changes.wrapping_add(1);

                atoms[second_index].charge = atoms[second_index].charge.wrapping_add(1);
                atoms[second_index].bond_type[0] = BOND_TYPE_DOUBLE as u8;
                atoms[center_index].bond_type[second_ordinal as usize] = BOND_TYPE_DOUBLE as u8;
                atoms[second_index].bond_stereo[0] = 0;
                atoms[center_index].bond_stereo[second_ordinal as usize] = 0;
                atoms[second_index].chem_bonds_valence =
                    atoms[second_index].chem_bonds_valence.wrapping_add(1);
                atoms[center_index].chem_bonds_valence =
                    atoms[center_index].chem_bonds_valence.wrapping_add(1);
                if fix_bug != 0 {
                    atoms[center_index].charge = atoms[center_index].charge.wrapping_sub(1);
                }
                number_of_changes = number_of_changes.wrapping_add(1);
                break;
            }
            second_ordinal = second_ordinal.wrapping_add(1);
        }
    }

    for first_index in 0..atom_count {
        if atoms[first_index].radical == RADICAL_DOUBLET as S_CHAR {
            let mut first_ordinal = 0_i32;
            let mut second_ordinal = 0_i32;
            let second_index = the_only_doublet_neigh(
                atoms,
                first_index as i32,
                &mut first_ordinal,
                &mut second_ordinal,
            )?;
            if second_index >= 0 {
                let second_index = second_index as usize;
                if atoms[first_index].bond_type[first_ordinal as usize] <= BOND_TYPE_DOUBLE as u8 {
                    atoms[first_index].bond_type[first_ordinal as usize] =
                        atoms[first_index].bond_type[first_ordinal as usize].wrapping_add(1);
                    atoms[first_index].chem_bonds_valence =
                        atoms[first_index].chem_bonds_valence.wrapping_add(1);
                    atoms[second_index].bond_type[second_ordinal as usize] =
                        atoms[second_index].bond_type[second_ordinal as usize].wrapping_add(1);
                    atoms[second_index].chem_bonds_valence =
                        atoms[second_index].chem_bonds_valence.wrapping_add(1);
                    atoms[first_index].radical = 0;
                    atoms[second_index].radical = 0;
                }
            }
        }
    }

    number_of_changes = number_of_changes.wrapping_add(remove_ion_pairs(number_of_atoms, atoms)?);
    Ok(number_of_changes)
}

pub(crate) fn post_fix_odd_things(_number_of_atoms: i32, _atoms: &mut [inp_ATOM]) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:922 post_fix_odd_things
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: post_fix_odd_things
    // INCHI✔️✔️: int post_fix_odd_things( int num_atoms, inp_ATOM *at )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int num_changes = 0;
    // INCHI✔️✔️:     /* currently does nothing */
    // INCHI✔️✔️:     return
    // INCHI✔️✔️:         num_changes;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: post_fix_odd_things
    // END INCHI C FUNCTION: post_fix_odd_things
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: post_fix_odd_things
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: post_fix_odd_things
    0
}

#[allow(non_snake_case)]
pub(crate) fn DisconnectOneLigand(
    atoms: &mut [inp_ATOM],
    mut old_component_numbers: Option<&mut [AT_NUMB]>,
    metal_marks: &[S_CHAR],
    heteroatom_elements: &[i8],
    number_of_halogens: i32,
    number_of_atoms: i32,
    metal_index: i32,
    ligand_ordinal: i32,
    mut taut_flags_done: Option<&mut INCHI_MODE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3112 DisconnectOneLigand
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: DisconnectOneLigand
    // INCHI✔️✔️: int DisconnectOneLigand( inp_ATOM *at,
    // INCHI✔️✔️:                          AT_NUMB *nOldCompNumber,
    // INCHI✔️✔️:                          S_CHAR *bMetal,
    // INCHI✔️✔️:                          char *elnumber_Heteroat,
    // INCHI✔️✔️:                          int num_halogens,
    // INCHI✔️✔️:                          int num_atoms,
    // INCHI✔️✔️:                          int iMetal,
    // INCHI✔️✔️:                          int jLigand,
    // INCHI✔️✔️:                          INCHI_MODE *bTautFlagsDone )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, j, iLigand, neigh, val;
    // INCHI✔️✔️:     int metal_neigh_ord[MAXVAL], num_neigh_arom_bonds[MAXVAL];
    // INCHI✔️✔️:     int num_metal_neigh, num_disconnections;
    // INCHI✔️✔️:     int num_del_arom_bonds, num_tot_arom_bonds, new_charge;
    // INCHI✔️✔️:     char *p;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     iLigand = at[iMetal].neighbor[jLigand];
    // INCHI✔️✔️:     num_metal_neigh = 0;
    // INCHI✔️✔️:     num_disconnections = 0;
    // INCHI✔️✔️:     num_del_arom_bonds = num_tot_arom_bonds = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* find bonds to disconnect */
    // INCHI✔️✔️:     for (i = 0; i < at[iLigand].valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         num_neigh_arom_bonds[i] = 0;
    // INCHI✔️✔️:         neigh = (int) at[iLigand].neighbor[i];
    // INCHI✔️✔️:         if (neigh < num_atoms && bMetal[neigh])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             metal_neigh_ord[num_metal_neigh++] = i;
    // INCHI✔️✔️:             if (at[iLigand].bond_type[i] > BOND_TYPE_TRIPLE)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* aromatic bond */
    // INCHI✔️✔️:                 for (j = 0; j < at[neigh].valence; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     num_neigh_arom_bonds[i] += ( at[neigh].bond_type[j] > BOND_TYPE_TRIPLE );
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 num_del_arom_bonds++;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         num_tot_arom_bonds += ( at[iLigand].bond_type[i] > BOND_TYPE_TRIPLE );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Disconnect */
    // INCHI✔️✔️:     if (num_del_arom_bonds)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* fix chem_valence of the ligand and its neighbors in case of disconnecting arom. bonds */
    // INCHI✔️✔️:         /* because in this case special care should be taken of updating at[].chem_bonds_valence */
    // INCHI✔️✔️:         for (i = 0; i < num_metal_neigh; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             j = metal_neigh_ord[i];
    // INCHI✔️✔️:             if (num_neigh_arom_bonds[j])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 neigh = at[iLigand].neighbor[j];
    // INCHI✔️✔️:                 at[neigh].chem_bonds_valence -= num_neigh_arom_bonds[j] / 2 - ( num_neigh_arom_bonds[j] - 1 ) / 2;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         at[iLigand].chem_bonds_valence -= num_tot_arom_bonds / 2 - ( num_tot_arom_bonds - num_del_arom_bonds ) / 2;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* disconnect in reverse order, otherwise the metal_neigh_ord[i]
    // INCHI✔️✔️:     becomes invalid after the first disconnection
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:     for (i = num_metal_neigh - 1; 0 <= i; i--)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         num_disconnections += DisconnectInpAtBond( at,
    // INCHI✔️✔️:                                                    nOldCompNumber,
    // INCHI✔️✔️:                                                    iLigand,
    // INCHI✔️✔️:                                                    metal_neigh_ord[i] );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* attempt to change ligand charge to make its valence 'natural' */
    // INCHI✔️✔️:     i = num_tot_arom_bonds - num_del_arom_bonds;
    // INCHI✔️✔️:     if ((i && i != 2 && i != 3) ||
    // INCHI✔️✔️:          (at[iLigand].radical && at[iLigand].radical != RADICAL_SINGLET) ||
    // INCHI✔️✔️:          !( p = strchr( elnumber_Heteroat, at[iLigand].el_number ) )) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto exit_function;  /* non-standard atom */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     val = at[iLigand].chem_bonds_valence + NUMH( at, iLigand );
    // INCHI✔️✔️:     new_charge = MAX_ATOMS; /* impossible value */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!val)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (p - elnumber_Heteroat < num_halogens)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             new_charge = -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = -1; i <= 1; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (val == get_el_valence( at[iLigand].el_number, i, 0 ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 new_charge = i; /* found charge that fits chem. valence */
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (new_charge != MAX_ATOMS)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (( new_charge != at[iLigand].charge ||
    // INCHI✔️✔️:             ( at[iLigand].radical && at[iLigand].radical != RADICAL_SINGLET ) ) &&
    // INCHI✔️✔️:              1 == num_metal_neigh)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (1 == new_charge && 4 == val && 2 == at[iLigand].valence &&
    // INCHI✔️✔️:                  4 == at[iLigand].chem_bonds_valence &&
    // INCHI✔️✔️:                  at[iLigand].bond_type[0] == at[iLigand].bond_type[1])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ; /* do not add +1 charge to disconnected =N=, etc. 2004-10-27 */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (bTautFlagsDone && new_charge != at[iLigand].charge)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     *bTautFlagsDone |= TG_FLAG_MOVE_CHARGE_COORD_DONE;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 at[iMetal].charge -= new_charge - at[iLigand].charge;
    // INCHI✔️✔️:                 at[iLigand].charge = new_charge;
    // INCHI✔️✔️:                 /*at[iLigand].radical = 0;*/
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_function:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_disconnections; /* ret;*/
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: DisconnectOneLigand
    // END INCHI C FUNCTION: DisconnectOneLigand
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectOneLigand
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define MAXVAL 20
    // INCHI✔️✔️: #define BOND_TYPE_TRIPLE 3
    // INCHI✔️✔️: #define NUMH(AT, N) (AT[N].num_H + NUM_ISO_H(AT, N))
    // INCHI✔️✔️: #define RADICAL_SINGLET 1
    // INCHI✔️✔️: #define MAX_ATOMS 32766
    // INCHI✔️✔️: #define TG_FLAG_MOVE_CHARGE_COORD_DONE 0x00000400
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisconnectOneLigand

    let metal = usize::try_from(metal_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let ordinal =
        usize::try_from(ligand_ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let metal_atom = atoms
        .get(metal)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if ligand_ordinal < 0 || ligand_ordinal >= i32::from(metal_atom.valence) {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let ligand = usize::from(
        *metal_atom
            .neighbor
            .get(ordinal)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let ligand_valence = i32::from(
        atoms
            .get(ligand)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .valence,
    );
    if ligand_valence > MAXVAL as i32 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let mut metal_neighbor_ordinals = [0_i32; MAXVAL as usize];
    let mut neighbor_aromatic_bond_counts = [0_i32; MAXVAL as usize];
    let mut number_of_metal_neighbors = 0_usize;
    let mut number_of_disconnections = 0_i32;
    let mut number_of_deleted_aromatic_bonds = 0_i32;
    let mut number_of_total_aromatic_bonds = 0_i32;

    for source_ordinal in 0..ligand_valence.max(0) as usize {
        let neighbor = usize::from(atoms[ligand].neighbor[source_ordinal]);
        let bond_type = atoms[ligand].bond_type[source_ordinal];
        if (neighbor as i32) < number_of_atoms
            && *metal_marks
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                != 0
        {
            metal_neighbor_ordinals[number_of_metal_neighbors] = source_ordinal as i32;
            number_of_metal_neighbors += 1;
            if i32::from(bond_type) > BOND_TYPE_TRIPLE as i32 {
                let neighbor_valence = i32::from(
                    atoms
                        .get(neighbor)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .valence,
                );
                for neighbor_ordinal in 0..neighbor_valence.max(0) as usize {
                    neighbor_aromatic_bond_counts[source_ordinal] += i32::from(
                        atoms[neighbor].bond_type[neighbor_ordinal] > BOND_TYPE_TRIPLE as u8,
                    );
                }
                number_of_deleted_aromatic_bonds += 1;
            }
        }
        number_of_total_aromatic_bonds += i32::from(i32::from(bond_type) > BOND_TYPE_TRIPLE as i32);
    }

    if number_of_deleted_aromatic_bonds != 0 {
        for source_index in 0..number_of_metal_neighbors {
            let source_ordinal = metal_neighbor_ordinals[source_index] as usize;
            let aromatic_count = neighbor_aromatic_bond_counts[source_ordinal];
            if aromatic_count != 0 {
                let neighbor = usize::from(atoms[ligand].neighbor[source_ordinal]);
                let decrement = aromatic_count / 2 - (aromatic_count - 1) / 2;
                atoms
                    .get_mut(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .chem_bonds_valence = atoms[neighbor]
                    .chem_bonds_valence
                    .wrapping_sub(decrement as S_CHAR);
            }
        }
        let decrement = number_of_total_aromatic_bonds / 2
            - (number_of_total_aromatic_bonds - number_of_deleted_aromatic_bonds) / 2;
        atoms[ligand].chem_bonds_valence = atoms[ligand]
            .chem_bonds_valence
            .wrapping_sub(decrement as S_CHAR);
    }

    for source_index in (0..number_of_metal_neighbors).rev() {
        number_of_disconnections = number_of_disconnections.wrapping_add(DisconnectInpAtBond(
            atoms,
            old_component_numbers.as_deref_mut(),
            ligand as i32,
            metal_neighbor_ordinals[source_index],
        )?);
    }

    let remaining_aromatic_bonds =
        number_of_total_aromatic_bonds - number_of_deleted_aromatic_bonds;
    if remaining_aromatic_bonds != 0
        && remaining_aromatic_bonds != 2
        && remaining_aromatic_bonds != 3
    {
        return Ok(number_of_disconnections);
    }
    let ligand_radical = atoms[ligand].radical;
    if ligand_radical != 0 && ligand_radical != RADICAL_SINGLET as S_CHAR {
        return Ok(number_of_disconnections);
    }

    let nul_position = heteroatom_elements
        .iter()
        .position(|&element| element == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let ligand_element = atoms[ligand].el_number;
    let Some(element_position) = heteroatom_elements[..=nul_position]
        .iter()
        .position(|&element| element as U_CHAR == ligand_element)
    else {
        return Ok(number_of_disconnections);
    };

    let atom = &atoms[ligand];
    let valence = i32::from(atom.chem_bonds_valence)
        + i32::from(atom.num_H)
        + atom.num_iso_H.into_iter().map(i32::from).sum::<i32>();
    let mut new_charge = MAX_ATOMS as i32;
    if valence == 0 {
        if (element_position as i32) < number_of_halogens {
            new_charge = -1;
        }
    } else {
        for charge in -1..=1 {
            if valence == get_el_valence(i32::from(ligand_element), charge, 0)? {
                new_charge = charge;
                break;
            }
        }
    }

    if new_charge != MAX_ATOMS as i32
        && (new_charge != i32::from(atoms[ligand].charge)
            || (atoms[ligand].radical != 0 && atoms[ligand].radical != RADICAL_SINGLET as S_CHAR))
        && number_of_metal_neighbors == 1
    {
        let preserve_charge = new_charge == 1
            && valence == 4
            && atoms[ligand].valence == 2
            && atoms[ligand].chem_bonds_valence == 4
            && atoms[ligand].bond_type[0] == atoms[ligand].bond_type[1];
        if !preserve_charge {
            let old_charge = atoms[ligand].charge;
            if let Some(flags) = taut_flags_done.as_deref_mut()
                && new_charge != i32::from(old_charge)
            {
                *flags |= TG_FLAG_MOVE_CHARGE_COORD_DONE as INCHI_MODE;
            }
            atoms[metal].charge = atoms[metal]
                .charge
                .wrapping_sub((new_charge - i32::from(old_charge)) as S_CHAR);
            atoms[ligand].charge = new_charge as S_CHAR;
        }
    }

    Ok(number_of_disconnections)
}

#[allow(non_snake_case)]
pub(crate) fn dist3D(first: &inp_ATOM, second: &inp_ATOM) -> f64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3245 dist3D
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: dist3D
    // INCHI✔️✔️: double dist3D( inp_ATOM *at1, inp_ATOM *at2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     double dx = at1->x - at2->x;
    // INCHI✔️✔️:     double dy = at1->y - at2->y;
    // INCHI✔️✔️:     double dz = at1->z - at2->z;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return sqrt( dx*dx + dy*dy + dz*dz );
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: dist3D
    // END INCHI C FUNCTION: dist3D
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: dist3D
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: sqrt resolves to the C math-library/compiler implementation for double.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: dist3D

    let dx = first.x - second.x;
    let dy = first.y - second.y;
    let dz = first.z - second.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[allow(non_snake_case)]
pub(crate) fn GetMinDistDistribution(
    atoms: &[inp_ATOM],
    number_of_atoms: i32,
    atom_index: i32,
    hydrogen_index: i32,
    in_all_components: i32,
    minimum_distances: &mut [f64],
    number_of_segments: i32,
) -> Result<f64, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3264 GetMinDistDistribution
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: GetMinDistDistribution
    // INCHI✔️✔️: double GetMinDistDistribution( inp_ATOM *at,
    // INCHI✔️✔️:                                int num_at,
    // INCHI✔️✔️:                                int iat,
    // INCHI✔️✔️:                                int iat_H,
    // INCHI✔️✔️:                                int bInAllComponents,
    // INCHI✔️✔️:                                double min_dist[],
    // INCHI✔️✔️:                                int num_segm )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /*    const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
    // INCHI✔️✔️:     const double one_pi = 3.14159265358979323846; /* M_PI */
    // INCHI✔️✔️:     const double two_pi = 2.0*one_pi;
    // INCHI✔️✔️:     const double f_step = two_pi / num_segm;
    // INCHI✔️✔️:     const double h_step = f_step / 2.0;
    // INCHI✔️✔️:     int i, j, k, kk, ki, kn, n, num_bonds;
    // INCHI✔️✔️:     double xi, yi, xn, yn, cross_prod_in, dot_prod_in, xni, yni, rni, tni, rmin;
    // INCHI✔️✔️:     double fi, fk, fn, ft = 0, rt = 0, rk, ri, rn, c, ave_bond_len;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < num_segm; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         min_dist[i] = MAX_BOND_LENGTH; /* more than any distance */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     num_bonds = 0;
    // INCHI✔️✔️:     ave_bond_len = 0.0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (i != iat && i != iat_H &&
    // INCHI✔️✔️:             ( bInAllComponents || at[i].component == at[iat].component ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             for (j = 0; j < at[i].valence; j++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 n = at[i].neighbor[j];
    // INCHI✔️✔️:                 if (( n > i && n != iat ) || n == iat_H)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI✔️✔️:                 if (n == iat)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     int stop = 1;  /* <BRKPT> */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:                 xi = at[i].x - at[iat].x;  /* ri; i != iat */
    // INCHI✔️✔️:                 yi = at[i].y - at[iat].y;
    // INCHI✔️✔️:                 xn = at[n].x - at[iat].x;  /* rn; possibly n == iat */
    // INCHI✔️✔️:                 yn = at[n].y - at[iat].y;
    // INCHI✔️✔️:                 cross_prod_in = xi*yn - xn*yi; /* ((r(i)-r(iat)) x (r(n)-r(iat)) */
    // INCHI✔️✔️:                 if (cross_prod_in < -0.01*MIN_BOND_LENGTH2)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* make sure the r(i)->r(n) vector is counterclockwise around at[iat] */
    // INCHI✔️✔️:                     inchi_swap( (char*) &xi, (char*) &xn, sizeof( xi ) );
    // INCHI✔️✔️:                     inchi_swap( (char*) &yi, (char*) &yn, sizeof( yi ) );
    // INCHI✔️✔️:                     /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 xni = xn - xi; /* r(n)->r(i) */
    // INCHI✔️✔️:                 yni = yn - yi;
    // INCHI✔️✔️:                 rni = xni*xni + yni*yni;
    // INCHI✔️✔️:                 if (rni > 0.01*MIN_BOND_LENGTH2)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* vector length |ri->rn| is not too small */
    // INCHI✔️✔️:                     /* arrowhead of the vector r(t) = ri + (rn-ri)*t; 0 <= t <= 1 points to the bond ri->rn */
    // INCHI✔️✔️:                     /* r(tni) is perpendicular to the bond ri->rn so that min|r(t)| = r(tni) = |tni|*rni */
    // INCHI✔️✔️:                     tni = -( xni*xi + yni*yi ) / rni;
    // INCHI✔️✔️:                     /* find min. distance from n-i bond to at[iat] */
    // INCHI✔️✔️:                     if (tni < 0.0)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         rmin = sqrt( xi*xi + yi*yi );
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else if (tni > 1.0)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         rmin = sqrt( xn*xn + yn*yn );
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         rmin = sqrt( tni*tni*rni );
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     ave_bond_len += sqrt( rni );
    // INCHI✔️✔️:                     num_bonds++;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* zero length i-n bond */
    // INCHI✔️✔️:                     tni = 0.5; /* fake */
    // INCHI✔️✔️:                     rmin = sqrt( xi*xi + yi*yi ); /* arbitrarily choose one */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 if (rmin >= 0.1*MIN_BOND_LENGTH)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* at[iat] does not belong to at[i]-at[n] bond */
    // INCHI✔️✔️:                     int    bCalc_rt = 1;
    // INCHI✔️✔️:                     fi = atan2( yi, xi );
    // INCHI✔️✔️:                     fn = ( n == iat ) ? fi : atan2( yn, xn );
    // INCHI✔️✔️:                     if (fi > fn)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* make sure fn - fi >= 0 */
    // INCHI✔️✔️:                         fn += two_pi;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     if (fi < 0.0)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         fi += two_pi;
    // INCHI✔️✔️:                         fn += two_pi;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     ki = (int) floor( ( fi + h_step ) / f_step );  /* cast does not match function type */
    // INCHI✔️✔️:                     kn = (int) floor( ( fn + h_step ) / f_step );
    // INCHI✔️✔️:
    // INCHI✔️✔️:                     /* the bond may affect several segments */
    // INCHI✔️✔️:                     for (k = ki; k <= kn; k++)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         kk = k % num_segm;
    // INCHI✔️✔️:                         if (min_dist[kk] < rmin)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             continue;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         if (bCalc_rt)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (n == iat)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 ft = fi;
    // INCHI✔️✔️:                                 rt = rmin;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             else
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 double xt, yt;
    // INCHI✔️✔️:                                 xt = xi + xni*tni;
    // INCHI✔️✔️:                                 yt = yi + yni*tni;
    // INCHI✔️✔️:                                 ft = atan2( yt, xt );
    // INCHI✔️✔️:                                 rt = sqrt( xt*xt + yt*yt );
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             bCalc_rt = 0;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         fk = f_step * kk;
    // INCHI✔️✔️:                         c = fabs( cos( fk - ft ) );
    // INCHI✔️✔️:                         if (c < MIN_COS)
    // INCHI✔️✔️:                             c = MIN_COS;
    // INCHI✔️✔️:                         rk = rt / c;
    // INCHI✔️✔️:                         if (min_dist[kk] > rk)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             min_dist[kk] = rk;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* rmin < 0.1*MIN_BOND_LENGTH */
    // INCHI✔️✔️:                     ri = xi*xi + yi*yi;
    // INCHI✔️✔️:                     rn = xn*xn + yn*yn;
    // INCHI✔️✔️:                     if (ri > MIN_BOND_LENGTH2 && rn > MIN_BOND_LENGTH2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         dot_prod_in = xn*xi + yn*yi;
    // INCHI✔️✔️:                         /* a very short bond */
    // INCHI✔️✔️:                         if (dot_prod_in > 0.01*MIN_BOND_LENGTH2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             /* bond does not cross at[iat] */
    // INCHI✔️✔️:                             double fyixi = atan2( yi, xi );
    // INCHI✔️✔️:                             if (fyixi < 0.0) fyixi += two_pi;
    // INCHI✔️✔️:                             kk = (int) floor( ( fyixi + h_step ) / f_step ) % num_segm;
    // INCHI✔️✔️:                             if (min_dist[kk] > rmin)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 min_dist[kk] = rmin;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         else if (dot_prod_in < -0.01*MIN_BOND_LENGTH2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             /* bond does cross at[iat] */
    // INCHI✔️✔️:                             double fyixi = atan2( yi, xi );
    // INCHI✔️✔️:                             if (fyixi < 0.0) fyixi += two_pi;
    // INCHI✔️✔️:                             kk = (int) floor( ( fyixi + h_step ) / f_step ) % num_segm;
    // INCHI✔️✔️:                             if (min_dist[kk] > rmin)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 min_dist[kk] = rmin;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             fyixi += one_pi;
    // INCHI✔️✔️:                             kk = (int) floor( ( fyixi + h_step ) / f_step ) % num_segm;
    // INCHI✔️✔️:                             if (min_dist[kk] > rmin)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 min_dist[kk] = rmin;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         else
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             ; /* error, should not happen */
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else if (ri <= MIN_BOND_LENGTH2 && rn <= MIN_BOND_LENGTH2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* a very short bond coincides with at[iat]; ignore */
    // INCHI✔️✔️:                         ;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* one end of the bond coincides with at[iat] */
    // INCHI✔️✔️:                         fi = ri > rn ? atan2( yi, xi ) : atan2( yn, xn );
    // INCHI✔️✔️:                         if (fi < 0.0) fi += two_pi;
    // INCHI✔️✔️:                         kk = (int) floor( ( fi + h_step ) / f_step ) % num_segm;
    // INCHI✔️✔️:                         if (min_dist[kk] > rmin)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             min_dist[kk] = rmin;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_bonds)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  ave_bond_len / (double) num_bonds;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0.0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: GetMinDistDistribution
    // END INCHI C FUNCTION: GetMinDistDistribution
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: GetMinDistDistribution
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define MIN_BOND_LENGTH (1.0e-6)
    // INCHI✔️✔️: #define MIN_COS (1.0e-6)
    // INCHI✔️✔️: #define MIN_BOND_LENGTH2 (MIN_BOND_LENGTH*MIN_BOND_LENGTH)
    // INCHI✔️✔️: #define MAX_BOND_LENGTH (1.0e30)
    // INCHI✔️✔️: bRELEASE_VERSION == 1, therefore the _DEBUG-only breakpoint block is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: GetMinDistDistribution

    if number_of_segments <= 0 {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    let segment_count = number_of_segments as usize;
    if minimum_distances.len() < segment_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    minimum_distances[..segment_count].fill(1.0e30);
    if number_of_atoms <= 0 {
        return Ok(0.0);
    }
    let atom_count = number_of_atoms as usize;
    if atoms.len() < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let target = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if target >= atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    const ONE_PI: f64 = 3.14159265358979323846;
    const TWO_PI: f64 = 2.0 * ONE_PI;
    const MIN_BOND_LENGTH_SQUARED: f64 = 1.0e-6 * 1.0e-6;
    const MIN_COSINE: f64 = 1.0e-6;
    let full_step = TWO_PI / f64::from(number_of_segments);
    let half_step = full_step / 2.0;
    let mut number_of_bonds = 0_i32;
    let mut average_bond_length = 0.0_f64;

    for source_index in 0..atom_count {
        if source_index as i32 == atom_index || source_index as i32 == hydrogen_index {
            continue;
        }
        if in_all_components == 0 && atoms[source_index].component != atoms[target].component {
            continue;
        }
        let source_valence = i32::from(atoms[source_index].valence);
        for neighbor_ordinal in 0..source_valence.max(0) as usize {
            let neighbor = i32::from(atoms[source_index].neighbor[neighbor_ordinal]);
            if (neighbor > source_index as i32 && neighbor != atom_index)
                || neighbor == hydrogen_index
            {
                continue;
            }
            let neighbor_index =
                usize::try_from(neighbor).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let neighbor_atom = atoms
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut xi = atoms[source_index].x - atoms[target].x;
            let mut yi = atoms[source_index].y - atoms[target].y;
            let mut xn = neighbor_atom.x - atoms[target].x;
            let mut yn = neighbor_atom.y - atoms[target].y;
            let cross_product = xi * yn - xn * yi;
            if cross_product < -0.01 * MIN_BOND_LENGTH_SQUARED {
                std::mem::swap(&mut xi, &mut xn);
                std::mem::swap(&mut yi, &mut yn);
            }

            let xni = xn - xi;
            let yni = yn - yi;
            let squared_bond_length = xni * xni + yni * yni;
            let projection;
            let minimum_distance;
            if squared_bond_length > 0.01 * MIN_BOND_LENGTH_SQUARED {
                projection = -(xni * xi + yni * yi) / squared_bond_length;
                minimum_distance = if projection < 0.0 {
                    (xi * xi + yi * yi).sqrt()
                } else if projection > 1.0 {
                    (xn * xn + yn * yn).sqrt()
                } else {
                    (projection * projection * squared_bond_length).sqrt()
                };
                average_bond_length += squared_bond_length.sqrt();
                number_of_bonds = number_of_bonds.wrapping_add(1);
            } else {
                projection = 0.5;
                minimum_distance = (xi * xi + yi * yi).sqrt();
            }

            if minimum_distance >= 0.1 * 1.0e-6 {
                let mut calculate_projection = true;
                let mut fi = yi.atan2(xi);
                let mut fn_angle = if neighbor == atom_index {
                    fi
                } else {
                    yn.atan2(xn)
                };
                if fi > fn_angle {
                    fn_angle += TWO_PI;
                }
                if fi < 0.0 {
                    fi += TWO_PI;
                    fn_angle += TWO_PI;
                }
                let first_segment = ((fi + half_step) / full_step).floor() as i32;
                let last_segment = ((fn_angle + half_step) / full_step).floor() as i32;
                let mut projection_angle = 0.0_f64;
                let mut projection_distance = 0.0_f64;
                for source_segment in first_segment..=last_segment {
                    let segment = (source_segment % number_of_segments) as usize;
                    if minimum_distances[segment] < minimum_distance {
                        continue;
                    }
                    if calculate_projection {
                        calculate_projection = false;
                        if neighbor == atom_index {
                            projection_angle = fi;
                            projection_distance = minimum_distance;
                        } else {
                            let xt = xi + xni * projection;
                            let yt = yi + yni * projection;
                            projection_angle = yt.atan2(xt);
                            projection_distance = (xt * xt + yt * yt).sqrt();
                        }
                    }
                    let segment_angle = full_step * segment as f64;
                    let mut cosine = (segment_angle - projection_angle).cos().abs();
                    if cosine < MIN_COSINE {
                        cosine = MIN_COSINE;
                    }
                    let radial_distance = projection_distance / cosine;
                    if minimum_distances[segment] > radial_distance {
                        minimum_distances[segment] = radial_distance;
                    }
                }
            } else {
                let source_radius = xi * xi + yi * yi;
                let neighbor_radius = xn * xn + yn * yn;
                if source_radius > MIN_BOND_LENGTH_SQUARED
                    && neighbor_radius > MIN_BOND_LENGTH_SQUARED
                {
                    let dot_product = xn * xi + yn * yi;
                    if dot_product > 0.01 * MIN_BOND_LENGTH_SQUARED {
                        let mut angle = yi.atan2(xi);
                        if angle < 0.0 {
                            angle += TWO_PI;
                        }
                        let segment = (((angle + half_step) / full_step).floor() as i32
                            % number_of_segments) as usize;
                        if minimum_distances[segment] > minimum_distance {
                            minimum_distances[segment] = minimum_distance;
                        }
                    } else if dot_product < -0.01 * MIN_BOND_LENGTH_SQUARED {
                        let mut angle = yi.atan2(xi);
                        if angle < 0.0 {
                            angle += TWO_PI;
                        }
                        let segment = (((angle + half_step) / full_step).floor() as i32
                            % number_of_segments) as usize;
                        if minimum_distances[segment] > minimum_distance {
                            minimum_distances[segment] = minimum_distance;
                        }
                        angle += ONE_PI;
                        let opposite_segment = (((angle + half_step) / full_step).floor() as i32
                            % number_of_segments)
                            as usize;
                        if minimum_distances[opposite_segment] > minimum_distance {
                            minimum_distances[opposite_segment] = minimum_distance;
                        }
                    }
                } else if source_radius <= MIN_BOND_LENGTH_SQUARED
                    && neighbor_radius <= MIN_BOND_LENGTH_SQUARED
                {
                } else {
                    let mut angle = if source_radius > neighbor_radius {
                        yi.atan2(xi)
                    } else {
                        yn.atan2(xn)
                    };
                    if angle < 0.0 {
                        angle += TWO_PI;
                    }
                    let segment = (((angle + half_step) / full_step).floor() as i32
                        % number_of_segments) as usize;
                    if minimum_distances[segment] > minimum_distance {
                        minimum_distances[segment] = minimum_distance;
                    }
                }
            }
        }
    }

    if number_of_bonds != 0 {
        Ok(average_bond_length / f64::from(number_of_bonds))
    } else {
        Ok(0.0)
    }
}

#[allow(non_snake_case)]
pub(crate) fn move_explicit_Hcation(
    atoms: &mut [inp_ATOM],
    number_of_atoms: i32,
    atom_index: i32,
    hydrogen_index: i32,
    in_all_components: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3480 move_explicit_Hcation
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: move_explicit_Hcation
    // INCHI✔️✔️: int move_explicit_Hcation( inp_ATOM *at,
    // INCHI✔️✔️:                            int num_at,
    // INCHI✔️✔️:                            int iat,
    // INCHI✔️✔️:                            int iat_H,
    // INCHI✔️✔️:                            int bInAllComponents )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define NUM_SEGM 20
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*    const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
    // INCHI✔️✔️:     const double one_pi = 3.14159265358979323846; /* M_PI */
    // INCHI✔️✔️:     const double two_pi = 2.0*one_pi;
    // INCHI✔️✔️:     const double f_step = two_pi / NUM_SEGM;
    // INCHI✔️✔️:     const double h_step = f_step / 2.0;
    // INCHI✔️✔️:     double min_dist[NUM_SEGM];
    // INCHI✔️✔️:     int nB, i, k, kk, next, val = 0;
    // INCHI✔️✔️:     double r, r0, xd, yd, zd, xr, yr, zr, ave_bond_len;
    // INCHI✔️✔️:     /*double step = 4.0*atan(1.0)/NUM_SEGM;*/
    // INCHI✔️✔️:     /* find at[iat] neighbors coordinates */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     xd = yd = zd = 0.0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at[iat].valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0, nB = 0, r = 0.0; i < at[iat].valence; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             next = at[iat].neighbor[i];
    // INCHI✔️✔️:             xd += at[next].x;
    // INCHI✔️✔️:             yd += at[next].y;
    // INCHI✔️✔️:             zd += at[next].z;
    // INCHI✔️✔️:             r += dist3D( at + iat, at + next );
    // INCHI✔️✔️:             nB++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         xd /= (double) nB;
    // INCHI✔️✔️:         yd /= (double) nB;
    // INCHI✔️✔️:         zd /= (double) nB;
    // INCHI✔️✔️:         r /= (double) nB;
    // INCHI✔️✔️:         r0 = sqrt( (double) ( xd - at[iat].x )*( xd - at[iat].x )
    // INCHI✔️✔️:                    + (double) ( yd - at[iat].y )*( yd - at[iat].y ) );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (at[iat_H].valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             r = dist3D( at + iat_H, at + (int) at[iat_H].neighbor[0] );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             r = 0.0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         r0 = 0.0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     ave_bond_len = GetMinDistDistribution( at, num_at, iat, iat_H,
    // INCHI✔️✔️:                                            bInAllComponents, min_dist,
    // INCHI✔️✔️:                                            NUM_SEGM );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (r < MIN_BOND_LENGTH && ave_bond_len > MIN_BOND_LENGTH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         r = ave_bond_len; /* ave_bond_len = 0.0 may mean that it is 0D structure */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (r > MIN_BOND_LENGTH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* process non-zero bond lengths */
    // INCHI✔️✔️:         double f;
    // INCHI✔️✔️:         if (10.0*r0 < r)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             xr = -r;     /* arbitrary */
    // INCHI✔️✔️:             yr = 0.0;
    // INCHI✔️✔️:             zr = 0.0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /*
    // INCHI✔️✔️:             if ( r0 < MIN_BOND_LENGTH ) {
    // INCHI✔️✔️:             r0 = 1.0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             */
    // INCHI✔️✔️:             xr = r * ( at[iat].x - xd ) / r0;
    // INCHI✔️✔️:             yr = r * ( at[iat].y - yd ) / r0; /* length = r */
    // INCHI✔️✔️:             zr = r * ( at[iat].z - zd ) / r0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:             /*          -- test: opposire direction --
    // INCHI✔️✔️:             xr =   -r * ( at[iat].x - xd )/r0;
    // INCHI✔️✔️:             yr =   -r * ( at[iat].y - yd )/r0;
    // INCHI✔️✔️:             zr =   -r * ( at[iat].z - zd )/r0;
    // INCHI✔️✔️:             */
    // INCHI✔️✔️:             if (xr*xr + yr*yr < 0.04*r*r)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 xr = -r;
    // INCHI✔️✔️:                 yr = 0.0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         r = sqrt( xr*xr + yr*yr );
    // INCHI✔️✔️:         f = atan2( yr, xr );
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (f < 0.0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             f += two_pi;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         kk = (int) floor( ( f + h_step ) / f_step ) % NUM_SEGM;
    // INCHI✔️✔️:         /* cast does not match function type by design */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (min_dist[kk] < 1.5* r)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             double dist = 1.5*r;
    // INCHI✔️✔️:             int start = -1, len = 0, start_max = -1, len_max = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         again:
    // INCHI✔️✔️:             /* look for longest kk interval with min_dist[kk] >= dist */
    // INCHI✔️✔️:             for (k = 0, start = 0, len = 0, len_max = 0; k < 2 * NUM_SEGM; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 kk = k % NUM_SEGM;
    // INCHI✔️✔️:                 if (min_dist[kk] >= dist)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (!len++)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         start = k;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (len > len_max)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         len_max = len;
    // INCHI✔️✔️:                         start_max = start;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     len = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (!len_max)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (dist > 0.1*r)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     dist *= 0.75;
    // INCHI✔️✔️:                     goto again;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     goto done; /* do it anyway */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* found a good sector */
    // INCHI✔️✔️:                 f = f_step * ( (double)start_max +  ((double)len_max - 1.0 ) / 2.0 ); /* djb-rwth: cast operators added */
    // INCHI✔️✔️:                 r0 = dist / 1.5;
    // INCHI✔️✔️:                 xr = r0 * cos( f );
    // INCHI✔️✔️:                 yr = r0 * sin( f );
    // INCHI✔️✔️:                 zr = zr / r*r0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         xr = yr = zr = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: done:
    // INCHI✔️✔️:     if (at[iat_H].valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* disconnect H */
    // INCHI✔️✔️:         next = at[iat_H].neighbor[0];
    // INCHI✔️✔️:         for (i = 0; i < at[next].valence; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (at[next].neighbor[i] == iat_H)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 RemoveInpAtBond( at, next, i );
    // INCHI✔️✔️:                 i = 0; /* success */
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* isolated H+ cation */
    // INCHI✔️✔️:         next = iat_H;
    // INCHI✔️✔️:         i = 0;
    // INCHI✔️✔️:         at[iat_H].valence = 1;
    // INCHI✔️✔️:         at[iat_H].chem_bonds_valence = 1;
    // INCHI✔️✔️:         at[iat_H].bond_type[0] = BOND_TYPE_SINGLE;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (0 == i /*i < at[next].valence*/)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* move charge */
    // INCHI✔️✔️:         if (at[next].charge > 0 && at[iat].charge < 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             at[next].charge--;
    // INCHI✔️✔️:             at[iat].charge++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* connect H to at[iat] */
    // INCHI✔️✔️:         val = at[iat].valence;
    // INCHI✔️✔️:
    // INCHI✔️✔️: #pragma warning (push)
    // INCHI✔️✔️: #pragma warning (disable: 6386)
    // INCHI✔️✔️:         if (val < MAXVAL)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             at[iat].neighbor[val] = iat_H;
    // INCHI✔️✔️:             at[iat].bond_type[val] = at[iat_H].bond_type[0];
    // INCHI✔️✔️:             at[iat].bond_stereo[val] = 0;
    // INCHI✔️✔️:             at[iat].chem_bonds_valence += at[iat_H].bond_type[0];
    // INCHI✔️✔️:             at[iat].valence = val + 1;
    // INCHI✔️✔️:         };
    // INCHI✔️✔️: #pragma warning (pop)
    // INCHI✔️✔️:
    // INCHI✔️✔️:         at[iat_H].component = at[iat].component;
    // INCHI✔️✔️:         at[iat_H].neighbor[0] = iat;
    // INCHI✔️✔️:         at[iat_H].bond_stereo[0] = 0; /* possible loss of stereo info */
    // INCHI✔️✔️:         at[iat_H].x = at[iat].x + xr;
    // INCHI✔️✔️:         at[iat_H].y = at[iat].y + yr;
    // INCHI✔️✔️:         at[iat_H].z = at[iat].z + zr;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         return 1; /* success */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0; /* failed */
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: move_explicit_Hcation
    // END INCHI C FUNCTION: move_explicit_Hcation
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: move_explicit_Hcation
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // INCHI✔️✔️: #define NUM_SEGM 20
    // INCHI✔️✔️: #define MIN_BOND_LENGTH (1.0e-6)
    // INCHI✔️✔️: #define MAXVAL 20
    // INCHI✔️✔️: #define BOND_TYPE_SINGLE 1
    // INCHI✔️✔️: MSVC #pragma warning directives have no GCC/Linux runtime behavior.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: move_explicit_Hcation

    const NUMBER_OF_SEGMENTS: usize = 20;
    const ONE_PI: f64 = 3.14159265358979323846;
    const TWO_PI: f64 = 2.0 * ONE_PI;
    const MINIMUM_BOND_LENGTH: f64 = 1.0e-6;
    let full_step = TWO_PI / NUMBER_OF_SEGMENTS as f64;
    let half_step = full_step / 2.0;
    let target = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let hydrogen =
        usize::try_from(hydrogen_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if target >= atoms.len() || hydrogen >= atoms.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let mut xd = 0.0_f64;
    let mut yd = 0.0_f64;
    let mut zd = 0.0_f64;
    let mut radius;
    let mut radial_offset;
    if atoms[target].valence != 0 {
        let target_valence = i32::from(atoms[target].valence);
        let mut number_of_bonds = 0_i32;
        radius = 0.0;
        for ordinal in 0..target_valence.max(0) as usize {
            let neighbor = usize::from(atoms[target].neighbor[ordinal]);
            let neighbor_atom = atoms
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            xd += neighbor_atom.x;
            yd += neighbor_atom.y;
            zd += neighbor_atom.z;
            radius += dist3D(&atoms[target], neighbor_atom);
            number_of_bonds = number_of_bonds.wrapping_add(1);
        }
        xd /= f64::from(number_of_bonds);
        yd /= f64::from(number_of_bonds);
        zd /= f64::from(number_of_bonds);
        radius /= f64::from(number_of_bonds);
        radial_offset = ((xd - atoms[target].x) * (xd - atoms[target].x)
            + (yd - atoms[target].y) * (yd - atoms[target].y))
            .sqrt();
    } else {
        radius = if atoms[hydrogen].valence != 0 {
            let neighbor = usize::from(atoms[hydrogen].neighbor[0]);
            dist3D(
                &atoms[hydrogen],
                atoms
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        } else {
            0.0
        };
        radial_offset = 0.0;
    }

    let mut minimum_distances = [0.0_f64; NUMBER_OF_SEGMENTS];
    let average_bond_length = GetMinDistDistribution(
        atoms,
        number_of_atoms,
        atom_index,
        hydrogen_index,
        in_all_components,
        &mut minimum_distances,
        NUMBER_OF_SEGMENTS as i32,
    )?;
    if radius < MINIMUM_BOND_LENGTH && average_bond_length > MINIMUM_BOND_LENGTH {
        radius = average_bond_length;
    }

    let mut xr;
    let mut yr;
    let mut zr;
    if radius > MINIMUM_BOND_LENGTH {
        if 10.0 * radial_offset < radius {
            xr = -radius;
            yr = 0.0;
            zr = 0.0;
        } else {
            xr = radius * (atoms[target].x - xd) / radial_offset;
            yr = radius * (atoms[target].y - yd) / radial_offset;
            zr = radius * (atoms[target].z - zd) / radial_offset;
            if xr * xr + yr * yr < 0.04 * radius * radius {
                xr = -radius;
                yr = 0.0;
            }
        }

        radius = (xr * xr + yr * yr).sqrt();
        let mut angle = yr.atan2(xr);
        if angle < 0.0 {
            angle += TWO_PI;
        }
        let segment =
            (((angle + half_step) / full_step).floor() as i32 % NUMBER_OF_SEGMENTS as i32) as usize;
        if minimum_distances[segment] < 1.5 * radius {
            let mut distance = 1.5 * radius;
            let mut start_max = -1_i32;
            loop {
                let mut start = 0_i32;
                let mut length = 0_i32;
                let mut maximum_length = 0_i32;
                for source_segment in 0..2 * NUMBER_OF_SEGMENTS as i32 {
                    let current_segment = source_segment as usize % NUMBER_OF_SEGMENTS;
                    if minimum_distances[current_segment] >= distance {
                        if length == 0 {
                            start = source_segment;
                        }
                        length += 1;
                    } else {
                        if length > maximum_length {
                            maximum_length = length;
                            start_max = start;
                        }
                        length = 0;
                    }
                }
                if maximum_length == 0 {
                    if distance > 0.1 * radius {
                        distance *= 0.75;
                        continue;
                    }
                    break;
                }
                angle =
                    full_step * (f64::from(start_max) + (f64::from(maximum_length) - 1.0) / 2.0);
                radial_offset = distance / 1.5;
                xr = radial_offset * angle.cos();
                yr = radial_offset * angle.sin();
                zr = zr / radius * radial_offset;
                break;
            }
        }
    } else {
        xr = 0.0;
        yr = 0.0;
        zr = 0.0;
    }

    let next;
    let removal_index;
    if atoms[hydrogen].valence != 0 {
        next = usize::from(atoms[hydrogen].neighbor[0]);
        let next_valence = i32::from(
            atoms
                .get(next)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .valence,
        );
        let mut found = next_valence <= 0;
        for ordinal in 0..next_valence.max(0) as usize {
            if atoms[next].neighbor[ordinal] == hydrogen as AT_NUMB {
                RemoveInpAtBond(atoms, next as i32, ordinal as i32)?;
                found = true;
                break;
            }
        }
        removal_index = if found { 0 } else { next_valence };
    } else {
        next = hydrogen;
        removal_index = 0;
        atoms[hydrogen].valence = 1;
        atoms[hydrogen].chem_bonds_valence = 1;
        atoms[hydrogen].bond_type[0] = BOND_TYPE_SINGLE as U_CHAR;
    }

    if removal_index == 0 {
        if atoms[next].charge > 0 && atoms[target].charge < 0 {
            atoms[next].charge = atoms[next].charge.wrapping_sub(1);
            atoms[target].charge = atoms[target].charge.wrapping_add(1);
        }
        let valence = i32::from(atoms[target].valence);
        if valence < MAXVAL as i32 {
            let ordinal =
                usize::try_from(valence).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            atoms[target].neighbor[ordinal] = hydrogen as AT_NUMB;
            atoms[target].bond_type[ordinal] = atoms[hydrogen].bond_type[0];
            atoms[target].bond_stereo[ordinal] = 0;
            atoms[target].chem_bonds_valence = atoms[target]
                .chem_bonds_valence
                .wrapping_add(atoms[hydrogen].bond_type[0] as S_CHAR);
            atoms[target].valence = (valence + 1) as S_CHAR;
        }
        atoms[hydrogen].component = atoms[target].component;
        atoms[hydrogen].neighbor[0] = target as AT_NUMB;
        atoms[hydrogen].bond_stereo[0] = 0;
        atoms[hydrogen].x = atoms[target].x + xr;
        atoms[hydrogen].y = atoms[target].y + yr;
        atoms[hydrogen].z = atoms[target].z + zr;
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn remove_terminal_HDT(
    heap: &mut SourceHeap,
    num_atoms: i32,
    at: SourceMutPointer<inp_ATOM>,
    b_fix_term_h_charge: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3723 remove_terminal_HDT
    // INCHI✔️❌: int remove_terminal_HDT( int num_atoms, inp_ATOM *at, int bFixTermHChrg )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB   *new_ord;
    // INCHI✔️❌:     inp_ATOM  *new_at;
    // INCHI✔️❌:     char *p;
    // INCHI✔️❌:     static const char szHDT[] = "HDT";
    // INCHI✔️❌:     static const int  kMax = sizeof( szHDT ); /*  = 4 */
    // INCHI✔️❌:     int ret = -1;
    // INCHI✔️❌:     int num_hydrogens = 0, num_H = 0;  /*  number of terminal H, D, T */
    // INCHI✔️❌:     int i, j, k, n, m;
    // INCHI✔️❌:     int val;
    // INCHI✔️❌:     AT_RANK new_HydrogenAt_order[NUM_H_ISOTOPES + 1];
    // INCHI✔️❌:     AT_RANK new_OtherNeigh_order[MAXVAL];
    // INCHI✔️❌:     S_CHAR  old_trans[MAX_NUM_STEREO_BONDS];
    // INCHI✔️❌:
    // INCHI✔️❌:     int  num_OtherNeigh, num_HydrogenAt;
    // INCHI✔️❌:
    // INCHI✔️❌:     new_ord = (AT_NUMB *) inchi_calloc( num_atoms, sizeof( new_ord[0] ) ); /* changed malloc to calloc 9-11-2003 */
    // INCHI✔️❌:     new_at = (inp_ATOM  *) inchi_malloc( sizeof( new_at[0] ) *num_atoms );
    // INCHI✔️❌:     if (!new_ord || !new_at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  move H. D, T to the end of the list of atoms */
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[i].component = i; /*  temporarily save original numbering */
    // INCHI✔️❌:                              /*  get k = temp. hydrogen isotope/non-hydrogen atom type: */
    // INCHI✔️❌:                              /*  k=0:H, k=2:D, k=3:T, k=4=kMax: not a hydrogen */
    // INCHI✔️❌:         k = at[i].elname[1] ? kMax : ( p = (char*) strchr( szHDT, at[i].elname[0] ) ) ? (int) ( p - szHDT ) : kMax;
    // INCHI✔️❌:         /*  set hydrogen isotope atw differences */
    // INCHI✔️❌:         /*  Notes: k-value of isotopic H is incremented to correct iso_atw_diff value later. */
    // INCHI✔️❌:         /*         1H isotope cannot be detected here. */
    // INCHI✔️❌:         if (k == ATW_H || k == ATW_H + 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* D or T, k = 1 or 2 */
    // INCHI✔️❌:             at[i].elname[0] = 'H'; /*  hydrogen isotope */
    // INCHI✔️❌:             at[i].iso_atw_diff = ++k; /*  increment k to make k = iso_atw_diff ( 2 for D, 3 for T ) */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         num_H += ( k != kMax && at[i].valence == 1 && at[i].chem_bonds_valence == 1 && !NUMH( at, i ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* special case: HD, HT, DT, HH: the only non-isotopic H or
    // INCHI✔️❌:     * the lightest isotopic H out of two is removed
    // INCHI✔️❌:     * to become implicit (make the heavier H the "central atom").
    // INCHI✔️❌:     * Note: This must be consistent with MOL_FMT_to_atom()
    // INCHI✔️❌:     * treatment of isotopic Hn aliases.
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (2 == num_H && 2 == num_atoms && !NUMH( at, 0 ) && !NUMH( at, 1 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         if (at[0].iso_atw_diff >= at[1].iso_atw_diff)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             new_ord[0] = 0;
    // INCHI✔️❌:             new_ord[1] = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             new_ord[0] = 1;
    // INCHI✔️❌:             new_ord[1] = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[new_ord[1]].charge)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[new_ord[0]].charge += at[new_ord[1]].charge;
    // INCHI✔️❌:             at[new_ord[1]].charge = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         new_at[new_ord[0]] = at[0];
    // INCHI✔️❌:         new_at[new_ord[1]] = at[1];
    // INCHI✔️❌:         num_hydrogens = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* general case except H-H */
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             k = ( at[i].elname[1] || NUMH( at, i ) ) ? kMax : ( at[i].elname[0] == 'H' ) ? at[i].iso_atw_diff : kMax;
    // INCHI✔️❌:             if (k < kMax && at[i].valence == 1 && at[i].chem_bonds_valence == 1 &&
    // INCHI✔️❌:                  /*  the order of comparison is important */
    // INCHI✔️❌:                 ( ( n = (int) at[i].neighbor[0] ) > i               /* at[n] has not been encountered yet*/ ||
    // INCHI✔️❌:                  (int) new_ord[n] < num_atoms - num_hydrogens ) /* at[n] might have been encountered; it has not been moved */)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*  found an explicit terminal hydrogen */
    // INCHI✔️❌:                 num_hydrogens++;
    // INCHI✔️❌:                 if (k == 0 && ATW_H <= at[i].iso_atw_diff && at[i].iso_atw_diff < ATW_H + NUM_H_ISOTOPES)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     k = at[i].iso_atw_diff; /*  H isotope has already been marked above or elsewhere */ /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[i].charge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  transfer charge from the hydrogen */
    // INCHI✔️❌:                     at[n].charge += at[i].charge;
    // INCHI✔️❌:                     at[i].charge = 0;
    // INCHI✔️❌:                     if (bFixTermHChrg)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*    Fixed bug (July 6, 2008 IPl) :
    // INCHI✔️❌:                         if terminal H was charged (not neutralized before call of remove_terminal_HDT)
    // INCHI✔️❌:                         and had an ordering number > than that of heavy-atom neighbour, then
    // INCHI✔️❌:                         charge on neighbour atom was not adjusted (though charge on H was removed). */
    // INCHI✔️❌:                         if (i > n)
    // INCHI✔️❌:                             /* new_at[new_ord[n]] has been created and filled already */
    // INCHI✔️❌:                             new_at[new_ord[n]].charge = at[n].charge;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 new_ord[i] = num_atoms - num_hydrogens;  /*  move hydrogens to the end of the list */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* atom is not an explicit terminal hydrogen */
    // INCHI✔️❌:                 new_ord[i] = i - num_hydrogens;  /*  adjust non-hydrogens positions */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /*  copy atom to the new position */
    // INCHI✔️❌:             new_at[new_ord[i]] = at[i];
    // INCHI✔️❌:         } /* i */
    // INCHI✔️❌:     } /* general case except H-H */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_hydrogens)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int num_others = num_atoms - num_hydrogens; /*  atoms which are not terminal H, D, T */
    // INCHI✔️❌:         if (num_hydrogens > 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  sort hydrogen isotopes in ascending order, */
    // INCHI✔️❌:             /*  orig, numbers being the secondary sorting key */
    // INCHI✔️❌:             qsort( new_at + num_others, num_hydrogens, sizeof( new_at[0] ), cmp_iso_atw_diff_component_no );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*  save new numbering of hydrogen atoms using temporarily saved orig numbering */
    // INCHI✔️❌:         for (i = num_others; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             new_ord[(int) new_at[i].component] = i;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  renumber neighbors according to new_ord[] and detach terminal hydrogens */
    // INCHI✔️❌:         for (i = 0; i < num_others; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memset( new_HydrogenAt_order, 0, sizeof( new_HydrogenAt_order ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             memset( new_OtherNeigh_order, 0, sizeof( new_OtherNeigh_order ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             num_OtherNeigh = 0;
    // INCHI✔️❌:             num_HydrogenAt = 0;
    // INCHI✔️❌:             num_H = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:             for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 old_trans[m] = 2 - ( new_at[i].sn_ord[m] + new_at[i].sb_ord[m] + ( new_at[i].sn_ord[m] > new_at[i].sb_ord[m] ) ) % 2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             for (k = val = 0; k < new_at[i].valence; k++) /* djb-rwth: removing redundant variables/code */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (num_others <= ( n = new_ord[new_at[i].neighbor[k]] ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  discovered neighbor = disconnected explicit hydrogen
    // INCHI✔️❌:                     *  i = new atom new_at[i] ordering number
    // INCHI✔️❌:                     *  n = new number of the explicit H
    // INCHI✔️❌:                     *  k = ordering number of the explicit H in new_at[i] adjacency list
    // INCHI✔️❌:                     */
    // INCHI✔️❌:                     if (0 < new_at[n].iso_atw_diff && new_at[n].iso_atw_diff < ATW_H + NUM_H_ISOTOPES)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* make explicit isotopic H implicit */
    // INCHI✔️❌:                         new_at[i].num_iso_H[new_at[n].iso_atw_diff - 1] ++; /*  isotopic H */
    // INCHI✔️❌:                         num_HydrogenAt += !new_HydrogenAt_order[new_at[n].iso_atw_diff];
    // INCHI✔️❌:                         new_HydrogenAt_order[new_at[n].iso_atw_diff] = k + 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* make explicit non-isotopic H implicit */
    // INCHI✔️❌:                         new_at[i].num_H++; /*  non-isotopic H */
    // INCHI✔️❌:                         num_HydrogenAt += !num_H;
    // INCHI✔️❌:                         num_H++;
    // INCHI✔️❌:                         new_HydrogenAt_order[0] = k + 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /*  decrement chem. bonds valence because one bond is removed */
    // INCHI✔️❌:                     new_at[i].chem_bonds_valence = inchi_max( 0, new_at[i].chem_bonds_valence - 1 );
    // INCHI✔️❌:                     new_at[n].neighbor[0] = i; /*  update removed hydrogen neighbor number */
    // INCHI✔️❌:                     if (new_at[i].sb_parity[0])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* if the removed H is an SB neighbor then mark it as removed */
    // INCHI✔️❌:                         for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (k == (int) new_at[i].sn_ord[m])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 new_at[i].sn_ord[m] = -( new_at[n].iso_atw_diff + 1 );
    // INCHI✔️❌:                                 /* means the SB neighbor has been removed; (-4)=H, (-3)=1H, (-2)=D, (-1)=T */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* discovered a regular (not an explicit H) neighbor */
    // INCHI✔️❌:                     if (new_at[i].sb_parity[0])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (num_OtherNeigh < MAX_NUM_STEREO_BONDS)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             new_OtherNeigh_order[num_OtherNeigh] = k + 1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         num_OtherNeigh++; /* increment outside of if() to detect overflow */
    // INCHI✔️❌:                         if (val != k)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* store new stereobond and sb-neighbor ordering numbers */
    // INCHI✔️❌:                             for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (k == (int) new_at[i].sb_ord[m])
    // INCHI✔️❌:                                     new_at[i].sb_ord[m] = val;
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                     if (k == (int) new_at[i].sn_ord[m])
    // INCHI✔️❌:                                         new_at[i].sn_ord[m] = val;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     new_at[i].neighbor[val] = new_ord[new_at[i].neighbor[k]];
    // INCHI✔️❌:                     new_at[i].bond_type[val] = new_at[i].bond_type[k];
    // INCHI✔️❌:                     new_at[i].bond_stereo[val] = new_at[i].bond_stereo[k];
    // INCHI✔️❌:                     val++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (new_at[i].valence > val && new_at[i].sb_parity[0])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (num_HydrogenAt == new_at[i].valence - val && num_HydrogenAt + num_OtherNeigh <= MAXVAL)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* recalculate parity so that it would describe neighbor sequence H,1H,D,T,neigh[0],neigh[1]... */
    // INCHI✔️❌:                     memmove(new_OtherNeigh_order + num_HydrogenAt, new_OtherNeigh_order, num_OtherNeigh * sizeof(new_OtherNeigh_order[0]));
    // INCHI✔️❌:                     for (k = 0, j = 1; k <= NUM_H_ISOTOPES; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (new_HydrogenAt_order[k] && (num_HydrogenAt - j < MAXVAL) && (num_HydrogenAt - j >= 0)) /* djb-rwth: fixing buffer overruns */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             new_OtherNeigh_order[num_HydrogenAt - j] = new_HydrogenAt_order[k];
    // INCHI✔️❌:                             for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if ((int) new_at[i].sn_ord[m] == -( k + 1 ))
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     new_at[i].sn_ord[m] = -j;
    // INCHI✔️❌:                                     /* negative means explicit H isotope ord are
    // INCHI✔️❌:                                     (contiguously) in front of the adjacency list */
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             j++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* at this point new_OtherNeigh_order[] contains
    // INCHI✔️❌:                     incremented old ordering numbers in new order */
    // INCHI✔️❌:                     k = insertions_sort_AT_RANK( new_OtherNeigh_order, num_HydrogenAt + num_OtherNeigh ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                                /*if ( k ) {*/
    // INCHI✔️❌:                                /*
    // INCHI✔️❌:                                for ( m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m ++ ) {
    // INCHI✔️❌:                                if ( PARITY_WELL_DEF(new_at[i].sb_parity[m]) ) {
    // INCHI✔️❌:                                if ( old_trans[m] != 2 - (4 + new_at[i].sn_ord[m] + new_at[i].sb_ord[m] + (new_at[i].sn_ord[m] > new_at[i].sb_ord[m]))%2 ) {
    // INCHI✔️❌:                                new_at[i].sb_parity[m] = 3 - new_at[i].sb_parity[m];
    // INCHI✔️❌:                                }
    // INCHI✔️❌:                                }
    // INCHI✔️❌:                                }
    // INCHI✔️❌:                                */
    // INCHI✔️❌:                                /*}*/
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             new_at[i].valence = val;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         memcpy(at, new_at, sizeof(at[0])* num_atoms);
    // INCHI✔️❌:         ret = num_others;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = num_atoms;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (new_ord)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( new_ord );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (new_at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( new_at );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: remove_terminal_HDT
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: remove_terminal_HDT
    // INCHI✔️❌: #define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
    // INCHI✔️❌: #define NUMH(AT,N)     (AT[N].num_H+NUM_ISO_H(AT,N))
    // INCHI✔️❌: #define NUM_H_ISOTOPES 3
    // INCHI✔️❌: #define MAXVAL 20
    // INCHI✔️❌: #define MAX_NUM_STEREO_BONDS 3
    // INCHI✔️❌: #define ATW_H 1
    // INCHI✔️❌: typedef unsigned short AT_NUMB;
    // INCHI✔️❌: typedef unsigned short AT_RANK;
    // INCHI✔️❌: typedef signed char S_CHAR;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: remove_terminal_HDT

    let count = usize::try_from(num_atoms)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let source_atoms = heap
        .slice(at.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();

    let new_ord_pointer =
        match inchi_calloc::<AT_NUMB>(heap, count as u64, std::mem::size_of::<AT_NUMB>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    let new_at_pointer = match heap.allocate(vec![inp_ATOM::default(); count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => {
            inchi_free(heap, new_ord_pointer)?;
            return Err(error);
        }
    };

    let processing = (|| -> Result<i32, SourceHeapError> {
        if new_ord_pointer.is_null() || new_at_pointer.is_null() {
            return Ok(-1);
        }

        let mut atoms = source_atoms;
        let mut new_ord = vec![AT_NUMB::default(); count];
        let mut new_at = vec![inp_ATOM::default(); count];
        let mut num_hydrogens = 0_usize;
        let mut num_h = 0_i32;
        const K_MAX: i32 = 4;

        let total_hydrogens = |atom: &inp_ATOM| {
            i32::from(atom.num_H) + atom.num_iso_H.iter().copied().map(i32::from).sum::<i32>()
        };

        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.component = index as AT_NUMB;
            let mut kind = if atom.elname[1] != 0 {
                K_MAX
            } else {
                [b'H', b'D', b'T']
                    .iter()
                    .position(|symbol| atom.elname[0] == *symbol as i8)
                    .map_or(K_MAX, |value| value as i32)
            };
            if kind == ATW_H as i32 || kind == ATW_H as i32 + 1 {
                atom.elname[0] = b'H' as i8;
                kind += 1;
                atom.iso_atw_diff = kind as S_CHAR;
            }
            num_h += i32::from(
                kind != K_MAX
                    && atom.valence == 1
                    && atom.chem_bonds_valence == 1
                    && total_hydrogens(atom) == 0,
            );
        }

        if num_h == 2
            && count == 2
            && total_hydrogens(&atoms[0]) == 0
            && total_hydrogens(&atoms[1]) == 0
        {
            if atoms[0].iso_atw_diff >= atoms[1].iso_atw_diff {
                new_ord[0] = 0;
                new_ord[1] = 1;
            } else {
                new_ord[0] = 1;
                new_ord[1] = 0;
            }
            let removed = usize::from(new_ord[1]);
            let retained = usize::from(new_ord[0]);
            if atoms[removed].charge != 0 {
                atoms[retained].charge = atoms[retained].charge.wrapping_add(atoms[removed].charge);
                atoms[removed].charge = 0;
            }
            new_at[usize::from(new_ord[0])] = atoms[0].clone();
            new_at[usize::from(new_ord[1])] = atoms[1].clone();
            num_hydrogens = 1;
        } else {
            for index in 0..count {
                let mut kind = if atoms[index].elname[1] != 0 || total_hydrogens(&atoms[index]) != 0
                {
                    K_MAX
                } else if atoms[index].elname[0] == b'H' as i8 {
                    i32::from(atoms[index].iso_atw_diff)
                } else {
                    K_MAX
                };
                let neighbor = usize::from(atoms[index].neighbor[0]);
                let removable = kind < K_MAX
                    && atoms[index].valence == 1
                    && atoms[index].chem_bonds_valence == 1
                    && (neighbor > index
                        || usize::from(
                            *new_ord
                                .get(neighbor)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        ) < count - num_hydrogens);
                if removable {
                    num_hydrogens += 1;
                    if kind == 0
                        && ATW_H as i32 <= i32::from(atoms[index].iso_atw_diff)
                        && i32::from(atoms[index].iso_atw_diff) < (ATW_H + NUM_H_ISOTOPES) as i32
                    {
                        kind = i32::from(atoms[index].iso_atw_diff);
                    }
                    let _ = kind;
                    if atoms[index].charge != 0 {
                        let charge = atoms[index].charge;
                        atoms[neighbor].charge = atoms[neighbor].charge.wrapping_add(charge);
                        atoms[index].charge = 0;
                        if b_fix_term_h_charge != 0 && index > neighbor {
                            new_at[usize::from(new_ord[neighbor])].charge = atoms[neighbor].charge;
                        }
                    }
                    new_ord[index] = (count - num_hydrogens) as AT_NUMB;
                } else {
                    new_ord[index] = (index - num_hydrogens) as AT_NUMB;
                }
                new_at[usize::from(new_ord[index])] = atoms[index].clone();
            }
        }

        let result = if num_hydrogens != 0 {
            let num_others = count - num_hydrogens;
            if num_hydrogens > 1 {
                new_at[num_others..]
                    .sort_by(|first, second| cmp_iso_atw_diff_component_no(first, second).cmp(&0));
            }
            for (index, atom) in new_at.iter().enumerate().skip(num_others) {
                new_ord[usize::from(atom.component)] = index as AT_NUMB;
            }

            for index in 0..num_others {
                let mut hydrogen_order = [AT_RANK::default(); NUM_H_ISOTOPES as usize + 1];
                let mut other_order = [AT_RANK::default(); MAXVAL as usize];
                let mut old_trans = [S_CHAR::default(); MAX_NUM_STEREO_BONDS as usize];
                let mut num_other_neighbors = 0_usize;
                let mut num_hydrogen_atoms = 0_usize;
                let mut ordinary_h_count = 0_i32;

                for (stereo_index, old) in old_trans.iter_mut().enumerate() {
                    if new_at[index].sb_parity[stereo_index] == 0 {
                        break;
                    }
                    let sn = new_at[index].sn_ord[stereo_index];
                    let sb = new_at[index].sb_ord[stereo_index];
                    *old = 2 - (sn + sb + S_CHAR::from(sn > sb)) % 2;
                }

                let old_valence = usize::try_from(new_at[index].valence)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut val = 0_usize;
                for ordinal in 0..old_valence {
                    let original_neighbor = usize::from(new_at[index].neighbor[ordinal]);
                    let neighbor = usize::from(
                        *new_ord
                            .get(original_neighbor)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    if num_others <= neighbor {
                        let isotope = i32::from(new_at[neighbor].iso_atw_diff);
                        if isotope > 0 && isotope < (ATW_H + NUM_H_ISOTOPES) as i32 {
                            let isotope_index = (isotope - 1) as usize;
                            new_at[index].num_iso_H[isotope_index] =
                                new_at[index].num_iso_H[isotope_index].wrapping_add(1);
                            let order_index = isotope as usize;
                            num_hydrogen_atoms += usize::from(hydrogen_order[order_index] == 0);
                            hydrogen_order[order_index] = (ordinal + 1) as AT_RANK;
                        } else {
                            new_at[index].num_H = new_at[index].num_H.wrapping_add(1);
                            num_hydrogen_atoms += usize::from(ordinary_h_count == 0);
                            ordinary_h_count += 1;
                            hydrogen_order[0] = (ordinal + 1) as AT_RANK;
                        }
                        new_at[index].chem_bonds_valence =
                            0.max(new_at[index].chem_bonds_valence - 1);
                        new_at[neighbor].neighbor[0] = index as AT_NUMB;
                        if new_at[index].sb_parity[0] != 0 {
                            for stereo_index in 0..MAX_NUM_STEREO_BONDS as usize {
                                if new_at[index].sb_parity[stereo_index] == 0 {
                                    break;
                                }
                                if ordinal == new_at[index].sn_ord[stereo_index] as usize {
                                    new_at[index].sn_ord[stereo_index] =
                                        -(new_at[neighbor].iso_atw_diff + 1);
                                }
                            }
                        }
                    } else {
                        if new_at[index].sb_parity[0] != 0 {
                            if num_other_neighbors < MAX_NUM_STEREO_BONDS as usize {
                                other_order[num_other_neighbors] = (ordinal + 1) as AT_RANK;
                            }
                            num_other_neighbors += 1;
                            if val != ordinal {
                                for stereo_index in 0..MAX_NUM_STEREO_BONDS as usize {
                                    if new_at[index].sb_parity[stereo_index] == 0 {
                                        break;
                                    }
                                    if ordinal == new_at[index].sb_ord[stereo_index] as usize {
                                        new_at[index].sb_ord[stereo_index] = val as S_CHAR;
                                    } else if ordinal == new_at[index].sn_ord[stereo_index] as usize
                                    {
                                        new_at[index].sn_ord[stereo_index] = val as S_CHAR;
                                    }
                                }
                            }
                        }
                        new_at[index].neighbor[val] = new_ord[original_neighbor];
                        new_at[index].bond_type[val] = new_at[index].bond_type[ordinal];
                        new_at[index].bond_stereo[val] = new_at[index].bond_stereo[ordinal];
                        val += 1;
                    }
                }

                if old_valence > val && new_at[index].sb_parity[0] != 0 {
                    if num_hydrogen_atoms == old_valence - val
                        && num_hydrogen_atoms + num_other_neighbors <= MAXVAL as usize
                    {
                        other_order.copy_within(0..num_other_neighbors, num_hydrogen_atoms);
                        let mut contiguous = 1_i32;
                        for (isotope, source_order) in hydrogen_order.iter().copied().enumerate() {
                            let target = num_hydrogen_atoms as i32 - contiguous;
                            if source_order != 0 && target < MAXVAL as i32 && target >= 0 {
                                other_order[target as usize] = source_order;
                                for stereo_index in 0..MAX_NUM_STEREO_BONDS as usize {
                                    if new_at[index].sb_parity[stereo_index] == 0 {
                                        break;
                                    }
                                    if i32::from(new_at[index].sn_ord[stereo_index])
                                        == -(isotope as i32 + 1)
                                    {
                                        new_at[index].sn_ord[stereo_index] = -contiguous as S_CHAR;
                                    }
                                }
                                contiguous += 1;
                            }
                        }
                        let sort_count = num_hydrogen_atoms + num_other_neighbors;
                        let _ = insertions_sort_AT_RANK(&mut other_order, sort_count as i32)?;
                    }
                }
                new_at[index].valence = val as S_CHAR;
            }
            heap.slice_mut(at)?[..count].clone_from_slice(&new_at);
            num_others as i32
        } else {
            heap.slice_mut(at)?[..count].clone_from_slice(&atoms);
            num_atoms
        };

        heap.slice_mut(new_ord_pointer)?.copy_from_slice(&new_ord);
        heap.slice_mut(new_at_pointer)?.clone_from_slice(&new_at);
        Ok(result)
    })();

    inchi_free(heap, new_ord_pointer)?;
    inchi_free(heap, new_at_pointer)?;
    processing
}

#[allow(non_snake_case)]
pub(crate) fn RemoveInpAtBond(
    atoms: &mut [inp_ATOM],
    atom_index: i32,
    bond_ordinal: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2046 RemoveInpAtBond
    // INCHI✔️❌: int RemoveInpAtBond( inp_ATOM *atom, int iat, int k )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int      i, j, m, m2; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     inp_ATOM *at = atom + iat;
    // INCHI✔️❌:     inp_ATOM *at2 = NULL;
    // INCHI✔️❌:     int      val = at->valence - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (val >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int bond = at->bond_type[k];
    // INCHI✔️❌:         if (bond > BOND_TYPE_TRIPLE)
    // INCHI✔️❌:             bond = BOND_TYPE_SINGLE; /* added 08-06-2003 */
    // INCHI✔️❌:
    // INCHI✔️❌:                                      /* update CML tetrahedral atom parity. */
    // INCHI✔️❌:         if (at->p_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (m = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at->p_orig_at_num[m] == at->orig_at_number)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at->p_parity = 0;
    // INCHI✔️❌:                     break; /* only 3 bonds are present; removing one bond removes stereo */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (at->p_parity /* at->valence == MAX_NUM_STEREO_ATOM_NEIGH*/)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* p_orig_at_num is a fixed size array of MAX_NUM_STEREO_ATOM_NEIGH (4) elements */
    // INCHI✔️❌:                 for (m = 0; m < at->valence && m < MAX_NUM_STEREO_ATOM_NEIGH; m++) /* djb-rwth: fixing GH PR #72 */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (atom[(int) at->neighbor[k]].orig_at_number == at->p_orig_at_num[m])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (m < at->valence && m < MAX_NUM_STEREO_ATOM_NEIGH) /* djb-rwth: fixing GH PR #72 */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at->p_orig_at_num[m] = at->orig_at_number;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at->p_parity = 0; /* wrong neighbors: at->neighbor[k] is not in the list of a stereo neighbors */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* update CML stereogenic bond parities; at this point no removed explicit H exist yet */
    // INCHI✔️❌:         if (at->sb_parity[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (m = 0; m < MAX_NUM_STEREO_BONDS && at->sb_parity[m]; )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (k == at->sb_ord[m] || (k == at->sn_ord[m] && val < 2 && ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* !!! FLAW: does take into account removed H !!! */
    // INCHI✔️❌:                     /* stereogenic bond is being removed OR */
    // INCHI✔️❌:                     /* remove stereogenic bond because its only neighbor is being removed */
    // INCHI✔️❌:                     int pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
    // INCHI✔️❌:                     int len = get_opposite_sb_atom( atom, iat, at->sb_ord[m], &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
    // INCHI✔️❌:                     if (len)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         i = pinxt_sb_parity_ord;
    // INCHI✔️❌:                         at2 = atom + pnxt_atom;
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         i = MAX_NUM_STEREO_BONDS;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /*
    // INCHI✔️❌:                     at2 = atom + at->neighbor[ (int)at->sb_ord[m] ];
    // INCHI✔️❌:                     for ( i = 0; i < MAX_NUM_STEREO_BONDS && at2->sb_parity[i]; i ++ )
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                     if ( iat == at2->neighbor[ (int)at2->sb_ord[i] ] )
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     */
    // INCHI✔️❌:                     if (i < MAX_NUM_STEREO_BONDS && at2->sb_parity[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         m2 = i;
    // INCHI✔️❌:                         /* remove bond parity from at */
    // INCHI✔️❌:                         if (m < MAX_NUM_STEREO_BONDS - 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memmove(at->sb_parity + m, at->sb_parity + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sb_parity[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                             memmove(at->sb_ord + m, at->sb_ord + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sb_ord[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                             memmove(at->sn_ord + m, at->sn_ord + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sn_ord[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                             memmove(at->sn_orig_at_num + m, at->sn_orig_at_num + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sn_orig_at_num[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         at->sb_parity[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         at->sb_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         at->sn_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         at->sn_orig_at_num[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         /* remove bond parity from at2 */
    // INCHI✔️❌:                         if (m2 < MAX_NUM_STEREO_BONDS - 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memmove(at2->sb_parity + m2, at2->sb_parity + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sb_parity[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                             memmove(at2->sb_ord + m2, at2->sb_ord + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sb_ord[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                             memmove(at2->sn_ord + m2, at2->sn_ord + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sn_ord[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                             memmove(at2->sn_orig_at_num + m2, at2->sn_orig_at_num + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sn_orig_at_num[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         at2->sb_parity[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         at2->sb_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         at2->sn_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         at2->sn_orig_at_num[MAX_NUM_STEREO_BONDS - 1] = 0;
    // INCHI✔️❌:                         /* do not increment m here because the array elements have been shifted */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         m++; /* program error: inconsistent stereobond parity */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (k == at->sn_ord[m])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* stereogenic bond neighbor is being removed; another neighbor remains */
    // INCHI✔️❌:                     /* !!! FLAW: does take into account removed H !!! */
    // INCHI✔️❌:                     for (j = 0, i = -1; j < at->valence; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (j != k && j != at->sb_ord[m])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             i = j;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* i is the position of the neighbor that will become a new neighbor */
    // INCHI✔️❌:                     /***************************************************************************
    // INCHI✔️❌:                     *  at->sb_parity[m] is the direction (EVEN=clockwise, ODD=counterclockwise)
    // INCHI✔️❌:                     *  from stereobond to the neighbor. If the neighbor is removed then
    // INCHI✔️❌:                     *  the parity should invert, otherwise it should be unchanged.
    // INCHI✔️❌:                     ***************************************************************************/
    // INCHI✔️❌:                     if (i < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* no alternative neighbor is available */
    // INCHI✔️❌:                         if (ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* parity cannot be not well-defined anymore */
    // INCHI✔️❌:                             int pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
    // INCHI✔️❌:                             int len = get_opposite_sb_atom( atom, iat, at->sb_ord[m], &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
    // INCHI✔️❌:                             if (len > 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 atom[pnxt_atom].sb_parity[pinxt_sb_parity_ord] = at->sb_parity[m] = AB_PARITY_UNDF;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         at->sn_ord[m] = -99; /* sb neighbor has been disconnected */
    // INCHI✔️❌:                         at->sb_ord[m] -= ( at->sb_ord[m] > k ); /* same as above */
    // INCHI✔️❌:                         at->sn_orig_at_num[m] = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (i < at->valence)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* choose another stereogenic bond neighbor, its ord. number is i before bond removal */
    // INCHI✔️❌:                         if (ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* ALL WRONG: 'move' previous stereo bond neighbor to the last position (pos. 2 out of 0,1,2) */
    // INCHI✔️❌:                             /* the parity of the transpositions is (2 - at->sn_ord[m])%2 = at->sn_ord[m] % 2 */
    // INCHI✔️❌:                             /* and replace the neighbor with another; the contribution to the parity is 1 */
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + at->sn_ord[m] + 1 ) % 2;*/
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + k + i +
    // INCHI✔️❌:                             (i > k) + (i > at->sb_ord[m]) ) % 2;*/
    // INCHI✔️❌:                             /*=== parity should be INVERTED ===*/
    // INCHI✔️❌:                             at->sb_parity[m] = 3 - at->sb_parity[m];
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         at->sn_ord[m] = i - ( i > k ); /* ord. number shifted because preceding bond is removed */
    // INCHI✔️❌:                         at->sb_ord[m] -= ( at->sb_ord[m] > k ); /* same as above */
    // INCHI✔️❌:                         at->sn_orig_at_num[m] = atom[(int) at->neighbor[i]].orig_at_number;
    // INCHI✔️❌:                         /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + 1 ) % 2;*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at->sb_parity[m] = 0; /* program error: inconsistent stereobond parity */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     m++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* removing another neighbor, k: first move it to the last position (pos. 2 out of 0,1,2) */
    // INCHI✔️❌:                     if (k < 2 && ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*at->sb_parity[m] =  2 - ( at->sb_parity[m] + k ) % 2;*/
    // INCHI✔️❌:                         /*at->sb_parity[m] =  2 - ( at->sb_parity[m] + (at->sn_ord[m] > k) + (at->sb_ord[m] > k) ) % 2;*/
    // INCHI✔️❌:                         ;/*==== Parity should remain UNCHANGED ===*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (at->sb_ord[m] > k)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at->sb_ord[m] --;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (at->sn_ord[m] > k)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at->sn_ord[m] --;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     m++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (k < val)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memmove(at->neighbor + k, at->neighbor + k + 1, sizeof(at->neighbor[0])* ((long long)val - (long long)k)); /* djb-rwth: cast operators added */
    // INCHI✔️❌:             memmove(at->bond_stereo + k, at->bond_stereo + k + 1, sizeof(at->bond_stereo[0])* ((long long)val - (long long)k)); /* djb-rwth: cast operators added */
    // INCHI✔️❌:             memmove(at->bond_type + k, at->bond_type + k + 1, sizeof(at->bond_type[0])* ((long long)val - (long long)k)); /* djb-rwth: cast operators added */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         at->neighbor[val] = 0;
    // INCHI✔️❌:         at->bond_stereo[val] = 0;
    // INCHI✔️❌:         at->bond_type[val] = 0;
    // INCHI✔️❌:         at->valence = val;
    // INCHI✔️❌:         at->chem_bonds_valence -= bond;
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: RemoveInpAtBond

    fn well_defined(parity: i8) -> bool {
        parity >= AB_MIN_WELL_DEFINED_PARITY as i8 && parity <= AB_MAX_WELL_DEFINED_PARITY as i8
    }

    fn remove_stereo_entry(atom: &mut inp_ATOM, index: usize) {
        let count = MAX_NUM_STEREO_BONDS as usize;
        if index < count - 1 {
            atom.sb_parity.copy_within(index + 1..count, index);
            atom.sb_ord.copy_within(index + 1..count, index);
            atom.sn_ord.copy_within(index + 1..count, index);
            atom.sn_orig_at_num.copy_within(index + 1..count, index);
        }
        atom.sb_parity[count - 1] = 0;
        atom.sb_ord[count - 1] = 0;
        atom.sn_ord[count - 1] = 0;
        atom.sn_orig_at_num[count - 1] = 0;
    }

    let iat = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let val = i32::from(
        atoms
            .get(iat)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .valence,
    ) - 1;
    if val < 0 {
        return Ok(0);
    }
    let k = usize::try_from(bond_ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut bond = *atoms
        .get(iat)
        .and_then(|atom| atom.bond_type.get(k))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if bond > BOND_TYPE_TRIPLE as u8 {
        bond = BOND_TYPE_SINGLE as u8;
    }

    if atoms[iat].p_parity != 0 {
        let original = atoms[iat].orig_at_number;
        if atoms[iat]
            .p_orig_at_num
            .iter()
            .any(|&number| number == original)
        {
            atoms[iat].p_parity = 0;
        }
        if atoms[iat].p_parity != 0 {
            let neighbor = usize::from(
                *atoms[iat]
                    .neighbor
                    .get(k)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let neighbor_original = atoms
                .get(neighbor)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .orig_at_number;
            let limit = usize::try_from(i32::from(atoms[iat].valence))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                .min(MAX_NUM_STEREO_ATOM_NEIGH as usize);
            if let Some(index) = atoms[iat].p_orig_at_num[..limit]
                .iter()
                .position(|&number| number == neighbor_original)
            {
                atoms[iat].p_orig_at_num[index] = original;
            } else {
                atoms[iat].p_parity = 0;
            }
        }
    }

    let mut m = 0_usize;
    while m < MAX_NUM_STEREO_BONDS as usize && atoms[iat].sb_parity[m] != 0 {
        let sb_ord = i32::from(atoms[iat].sb_ord[m]);
        let sn_ord = i32::from(atoms[iat].sn_ord[m]);
        let parity = atoms[iat].sb_parity[m];
        if bond_ordinal == sb_ord || (bond_ordinal == sn_ord && val < 2 && well_defined(parity)) {
            let mut next_atom = 0_i32;
            let mut next_to_current = 0_i32;
            let mut next_parity_ordinal = 0_i32;
            let len = get_opposite_sb_atom_slice(
                atoms,
                atom_index,
                sb_ord,
                Some(&mut next_atom),
                Some(&mut next_to_current),
                Some(&mut next_parity_ordinal),
            )?;
            let opposite_order = if len != 0 {
                usize::try_from(next_parity_ordinal)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?
            } else {
                MAX_NUM_STEREO_BONDS as usize
            };
            let opposite_index = usize::try_from(next_atom).ok();
            if opposite_order < MAX_NUM_STEREO_BONDS as usize
                && opposite_index
                    .and_then(|index| atoms.get(index))
                    .is_some_and(|atom| atom.sb_parity[opposite_order] != 0)
            {
                remove_stereo_entry(&mut atoms[iat], m);
                remove_stereo_entry(
                    atoms
                        .get_mut(opposite_index.unwrap())
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    opposite_order,
                );
            } else {
                m += 1;
            }
        } else if bond_ordinal == sn_ord {
            let valence = i32::from(atoms[iat].valence);
            let mut replacement = -1_i32;
            for candidate in 0..valence {
                if candidate != bond_ordinal && candidate != sb_ord {
                    replacement = candidate;
                    break;
                }
            }
            if replacement < 0 {
                if well_defined(atoms[iat].sb_parity[m]) {
                    let mut next_atom = 0_i32;
                    let mut next_to_current = 0_i32;
                    let mut next_parity_ordinal = 0_i32;
                    let len = get_opposite_sb_atom_slice(
                        atoms,
                        atom_index,
                        sb_ord,
                        Some(&mut next_atom),
                        Some(&mut next_to_current),
                        Some(&mut next_parity_ordinal),
                    )?;
                    if len > 0 {
                        let next = usize::try_from(next_atom)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let order = usize::try_from(next_parity_ordinal)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        atoms
                            .get_mut(next)
                            .and_then(|atom| atom.sb_parity.get_mut(order))
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone_from(&(AB_PARITY_UNDF as i8));
                        atoms[iat].sb_parity[m] = AB_PARITY_UNDF as i8;
                    }
                }
                atoms[iat].sn_ord[m] = -99;
                if i32::from(atoms[iat].sb_ord[m]) > bond_ordinal {
                    atoms[iat].sb_ord[m] = atoms[iat].sb_ord[m].wrapping_sub(1);
                }
                atoms[iat].sn_orig_at_num[m] = 0;
            } else if replacement < valence {
                if well_defined(atoms[iat].sb_parity[m]) {
                    atoms[iat].sb_parity[m] = 3_i8.wrapping_sub(atoms[iat].sb_parity[m]);
                }
                atoms[iat].sn_ord[m] = (replacement - i32::from(replacement > bond_ordinal)) as i8;
                if i32::from(atoms[iat].sb_ord[m]) > bond_ordinal {
                    atoms[iat].sb_ord[m] = atoms[iat].sb_ord[m].wrapping_sub(1);
                }
                let replacement_index = usize::try_from(replacement)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let neighbor = usize::from(atoms[iat].neighbor[replacement_index]);
                atoms[iat].sn_orig_at_num[m] = atoms
                    .get(neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .orig_at_number;
            } else {
                atoms[iat].sb_parity[m] = 0;
            }
            m += 1;
        } else {
            if i32::from(atoms[iat].sb_ord[m]) > bond_ordinal {
                atoms[iat].sb_ord[m] = atoms[iat].sb_ord[m].wrapping_sub(1);
            }
            if i32::from(atoms[iat].sn_ord[m]) > bond_ordinal {
                atoms[iat].sn_ord[m] = atoms[iat].sn_ord[m].wrapping_sub(1);
            }
            m += 1;
        }
    }

    let val_index = usize::try_from(val).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if bond_ordinal < val {
        atoms[iat].neighbor.copy_within(k + 1..=val_index, k);
        atoms[iat].bond_stereo.copy_within(k + 1..=val_index, k);
        atoms[iat].bond_type.copy_within(k + 1..=val_index, k);
    }
    atoms[iat].neighbor[val_index] = 0;
    atoms[iat].bond_stereo[val_index] = 0;
    atoms[iat].bond_type[val_index] = 0;
    atoms[iat].valence = val as i8;
    atoms[iat].chem_bonds_valence = atoms[iat].chem_bonds_valence.wrapping_sub(bond as i8);
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn DisconnectInpAtBond(
    atoms: &mut [inp_ATOM],
    mut old_component_numbers: Option<&mut [u16]>,
    atom_index: i32,
    neighbor_ordinal: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2261 DisconnectInpAtBond
    // INCHI✔️❌: int DisconnectInpAtBond( inp_ATOM *at,
    // INCHI✔️❌:                          AT_NUMB *nOldCompNumber,
    // INCHI✔️❌:                          int iat,
    // INCHI✔️❌:                          int neigh_ord )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int neigh, i, ret = 0;
    // INCHI✔️❌:     int component;
    // INCHI✔️❌:     neigh = at[iat].neighbor[neigh_ord];
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < at[neigh].valence; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (iat == (int) at[neigh].neighbor[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i < at[neigh].valence)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret += RemoveInpAtBond( at, iat, neigh_ord );
    // INCHI✔️❌:         ret += RemoveInpAtBond( at, neigh, i );
    // INCHI✔️❌:         if (nOldCompNumber && ret)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((component = at[iat].component)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nOldCompNumber[component - 1] = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((component = at[neigh].component)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nOldCompNumber[component - 1] = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( ret == 2 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DisconnectInpAtBond

    let iat = usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let ordinal =
        usize::try_from(neighbor_ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let neighbor = usize::from(
        *atoms
            .get(iat)
            .and_then(|atom| atom.neighbor.get(ordinal))
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let opposite = atoms
        .get(neighbor)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let opposite_valence = usize::try_from(i32::from(opposite.valence))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut reverse_ordinal = 0_usize;
    while reverse_ordinal < opposite_valence {
        if i32::from(
            *opposite
                .neighbor
                .get(reverse_ordinal)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ) == atom_index
        {
            break;
        }
        reverse_ordinal += 1;
    }
    if reverse_ordinal >= opposite_valence {
        return Ok(0);
    }

    let mut result = RemoveInpAtBond(atoms, atom_index, neighbor_ordinal)?;
    result += RemoveInpAtBond(atoms, neighbor as i32, reverse_ordinal as i32)?;
    if result != 0 {
        if let Some(numbers) = old_component_numbers.as_deref_mut() {
            let first_component = atoms[iat].component;
            if first_component != 0 {
                *numbers
                    .get_mut(usize::from(first_component - 1))
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            }
            let second_component = atoms[neighbor].component;
            if second_component != 0 {
                *numbers
                    .get_mut(usize::from(second_component - 1))
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            }
        }
    }
    Ok(i32::from(result == 2))
}

#[allow(non_snake_case)]
pub(crate) fn SetConnectedComponentNumber(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
    component_number: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4598 SetConnectedComponentNumber
    // INCHI✔️✔️: int SetConnectedComponentNumber( inp_ATOM *at, int num_at, int component_number )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     for (i = 0; i < num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         at[i].component = (AT_NUMB) component_number;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SetConnectedComponentNumber

    let count =
        usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let selected = atoms
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for atom in selected {
        atom.component = component_number as AT_NUMB;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn Free_INChI_Stereo(
    heap: &mut SourceHeap,
    p_inchi_stereo: SourceMutPointer<INChI_Stereo>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4611 Free_INChI_Stereo
    // INCHI✔❌: int Free_INChI_Stereo( INChI_Stereo *pINChI_Stereo )
    // INCHI✔❌: {
    // INCHI✔❌:     if (pINChI_Stereo)
    // INCHI✔❌:     {
    // INCHI✔❌:         qzfree( pINChI_Stereo->nNumber );
    // INCHI✔❌:         qzfree( pINChI_Stereo->t_parity );
    // INCHI✔❌:         qzfree( pINChI_Stereo->nNumberInv );
    // INCHI✔❌:         qzfree( pINChI_Stereo->t_parityInv );
    // INCHI✔❌:         qzfree( pINChI_Stereo->nBondAtom1 );
    // INCHI✔❌:         qzfree( pINChI_Stereo->nBondAtom2 );
    // INCHI✔❌:         qzfree( pINChI_Stereo->b_parity );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: Free_INChI_Stereo
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:1205 qzfree
    // INCHI✔❌: #define qzfree(X)   do{if(X){inchi_free(X);(X)=NULL;}}while(0)
    // END INCHI ACTIVE MACRO: qzfree

    if p_inchi_stereo.is_null() {
        return Ok(0);
    }

    let pointer = heap.slice(p_inchi_stereo.as_const())?[0].nNumber;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_stereo)?[0].nNumber = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_stereo.as_const())?[0].t_parity;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_stereo)?[0].t_parity = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_stereo.as_const())?[0].nNumberInv;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_stereo)?[0].nNumberInv = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_stereo.as_const())?[0].t_parityInv;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_stereo)?[0].t_parityInv = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_stereo.as_const())?[0].nBondAtom1;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_stereo)?[0].nBondAtom1 = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_stereo.as_const())?[0].nBondAtom2;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_stereo)?[0].nBondAtom2 = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_stereo.as_const())?[0].b_parity;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_stereo)?[0].b_parity = SourceMutPointer::null();
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn Alloc_INChI_Stereo(
    heap: &mut SourceHeap,
    num_atoms: i32,
    num_bonds: i32,
) -> Result<SourceMutPointer<INChI_Stereo>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4629 Alloc_INChI_Stereo
    // INCHI✔️✔️: INChI_Stereo *Alloc_INChI_Stereo( int num_at, int num_bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️:     INChI_Stereo *pINChI_Stereo = (INChI_Stereo *)
    // INCHI✔️✔️:         inchi_calloc( 1, sizeof( INChI_Stereo ) );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (pINChI_Stereo)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (num_at &&
    // INCHI✔️✔️:             ( pINChI_Stereo->nNumber = (AT_NUMB *) inchi_calloc( num_at, sizeof( pINChI_Stereo->nNumber[0] ) ) ) &&
    // INCHI✔️✔️:              ( pINChI_Stereo->t_parity = (S_CHAR  *) inchi_calloc( num_at, sizeof( pINChI_Stereo->t_parity[0] ) ) ) &&
    // INCHI✔️✔️:              ( pINChI_Stereo->nNumberInv = (AT_NUMB *) inchi_calloc( num_at, sizeof( pINChI_Stereo->nNumberInv[0] ) ) ) &&
    // INCHI✔️✔️:              ( pINChI_Stereo->t_parityInv = (S_CHAR  *) inchi_calloc( num_at, sizeof( pINChI_Stereo->t_parityInv[0] ) ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (num_at)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (num_bonds &&
    // INCHI✔️✔️:             ( pINChI_Stereo->nBondAtom1 = (AT_NUMB *) inchi_calloc( num_bonds, sizeof( pINChI_Stereo->nBondAtom1[0] ) ) ) &&
    // INCHI✔️✔️:              ( pINChI_Stereo->nBondAtom2 = (AT_NUMB *) inchi_calloc( num_bonds, sizeof( pINChI_Stereo->nBondAtom2[0] ) ) ) &&
    // INCHI✔️✔️:              ( pINChI_Stereo->b_parity = (S_CHAR  *) inchi_calloc( num_bonds, sizeof( pINChI_Stereo->b_parity[0] ) ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (num_bonds)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         return pINChI_Stereo;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     out_of_RAM:
    // INCHI✔️✔️:
    // INCHI✔️✔️:         Free_INChI_Stereo( pINChI_Stereo );
    // INCHI✔️✔️:         qzfree( pINChI_Stereo );
    // INCHI✔️✔️:     } /* if ( pINChI_Stereo )  */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return NULL;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: Alloc_INChI_Stereo

    fn allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }

    let stereo =
        match inchi_calloc::<INChI_Stereo>(heap, 1, std::mem::size_of::<INChI_Stereo>() as u64) {
            Ok(pointer) => pointer,
            Err(error) if allocation_error(error) => return Ok(SourceMutPointer::null()),
            Err(error) => return Err(error),
        };
    let atom_count = num_atoms as i64 as u64;
    let bond_count = num_bonds as i64 as u64;
    let allocation = (|| {
        if num_atoms != 0 {
            let pointer =
                inchi_calloc::<AT_NUMB>(heap, atom_count, std::mem::size_of::<AT_NUMB>() as u64)?;
            heap.slice_mut(stereo)?[0].nNumber = pointer;
            let pointer =
                inchi_calloc::<S_CHAR>(heap, atom_count, std::mem::size_of::<S_CHAR>() as u64)?;
            heap.slice_mut(stereo)?[0].t_parity = pointer;
            let pointer =
                inchi_calloc::<AT_NUMB>(heap, atom_count, std::mem::size_of::<AT_NUMB>() as u64)?;
            heap.slice_mut(stereo)?[0].nNumberInv = pointer;
            let pointer =
                inchi_calloc::<S_CHAR>(heap, atom_count, std::mem::size_of::<S_CHAR>() as u64)?;
            heap.slice_mut(stereo)?[0].t_parityInv = pointer;
        }
        if num_bonds != 0 {
            let pointer =
                inchi_calloc::<AT_NUMB>(heap, bond_count, std::mem::size_of::<AT_NUMB>() as u64)?;
            heap.slice_mut(stereo)?[0].nBondAtom1 = pointer;
            let pointer =
                inchi_calloc::<AT_NUMB>(heap, bond_count, std::mem::size_of::<AT_NUMB>() as u64)?;
            heap.slice_mut(stereo)?[0].nBondAtom2 = pointer;
            let pointer =
                inchi_calloc::<S_CHAR>(heap, bond_count, std::mem::size_of::<S_CHAR>() as u64)?;
            heap.slice_mut(stereo)?[0].b_parity = pointer;
        }
        Ok(())
    })();
    match allocation {
        Ok(()) => Ok(stereo),
        Err(error) => {
            Free_INChI_Stereo(heap, stereo)?;
            inchi_free(heap, stereo)?;
            if allocation_error(error) {
                Ok(SourceMutPointer::null())
            } else {
                Err(error)
            }
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn Alloc_INChI(
    heap: &mut SourceHeap,
    atoms: &[inp_ATOM],
    num_atoms: i32,
    found_num_bonds: &mut i32,
    found_num_isotopic: &mut i32,
    allocation_mode: i32,
) -> Result<SourceMutPointer<INChI>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4722 Alloc_INChI
    // INCHI✔️✔️: INChI *Alloc_INChI( inp_ATOM *at,
    // INCHI✔️✔️:                     int num_at,
    // INCHI✔️✔️:                     int *found_num_bonds,
    // INCHI✔️✔️:                     int *found_num_isotopic,
    // INCHI✔️✔️:                     int nAllocMode )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int    i, num_bonds, num_isotopic_atoms;
    // INCHI✔️✔️:     INChI  *pINChI;
    // INCHI✔️✔️:     int    bIsotopic = ( nAllocMode & REQ_MODE_ISO );
    // INCHI✔️✔️:     /* int    bTautomeric = (nAllocMode & REQ_MODE_TAUT); */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_at <= 0 ||
    // INCHI✔️✔️:          NULL == ( pINChI = (INChI *) inchi_calloc( 1, sizeof( INChI ) ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0, num_bonds = 0, num_isotopic_atoms = 0; i < num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         num_bonds += at[i].valence;
    // INCHI✔️✔️:         /* if ( bIsotopic ) { */
    // INCHI✔️✔️:         num_isotopic_atoms += ( 0 != at[i].iso_atw_diff ||
    // INCHI✔️✔️:                                 !strcmp( at[i].elname, "D" ) ||
    // INCHI✔️✔️:                                 !strcmp( at[i].elname, "T" ) ||
    // INCHI✔️✔️:                                 at[i].num_iso_H[0] ||
    // INCHI✔️✔️:                                 at[i].num_iso_H[1] ||
    // INCHI✔️✔️:                                 at[i].num_iso_H[2] );
    // INCHI✔️✔️:         /* } */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     num_bonds /= 2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     *found_num_bonds = num_bonds;
    // INCHI✔️✔️:     *found_num_isotopic = num_isotopic_atoms;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (( pINChI->nAtom = (U_CHAR*) inchi_calloc( num_at, sizeof( pINChI->nAtom[0] ) ) ) &&
    // INCHI✔️✔️:         ( pINChI->nConnTable = (AT_NUMB*) inchi_calloc( (long long)num_at + (long long)num_bonds, sizeof( pINChI->nConnTable[0] ) ) ) && /* djb-rwth: cast operator added */
    // INCHI✔️✔️:          ( pINChI->nTautomer = (AT_NUMB*) inchi_calloc( ( ( 3 + INCHI_T_NUM_MOVABLE )*(long long)num_at ) / 2 + 1, sizeof( pINChI->nTautomer[0] ) ) ) && /* djb-rwth: cast operator added */
    // INCHI✔️✔️:          ( pINChI->nNum_H = (S_CHAR*) inchi_calloc( num_at, sizeof( pINChI->nNum_H[0] ) ) ) &&
    // INCHI✔️✔️:          ( pINChI->nNum_H_fixed = (S_CHAR*) inchi_calloc( num_at, sizeof( pINChI->nNum_H_fixed[0] ) ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:         /* nTautomer length: max. number of tautomeric groups is num_at/2
    // INCHI✔️✔️:
    // INCHI✔️✔️:         1 word                     -> number of t-groups
    // INCHI✔️✔️:
    // INCHI✔️✔️:         each group has:
    // INCHI✔️✔️:
    // INCHI✔️✔️:         1 word                     -> number of endpoints+INCHI_T_NUM_MOVABLE
    // INCHI✔️✔️:         INCHI_T_NUM_MOVABLE words   -> number(s) of moveable attachments
    // INCHI✔️✔️:         numbers of endpoints words -> canon. numbers
    // INCHI✔️✔️:
    // INCHI✔️✔️:         max. occurs if each t-group has 2 atoms (num_at/2 t-groups) and all atoms
    // INCHI✔️✔️:         belong to t-groups (num_at endpoints)
    // INCHI✔️✔️:
    // INCHI✔️✔️:         Total: 1 + (number of t-groups)*(1+INCHI_T_NUM_MOVABLE) + (number of endpoints) <=
    // INCHI✔️✔️:         1 + (num_at/2) * (1+INCHI_T_NUM_MOVABLE) + num_at <=
    // INCHI✔️✔️:         1 + (3+INCHI_T_NUM_MOVABLE)*num_at/2 words.
    // INCHI✔️✔️:         */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto out_of_RAM;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pINChI->szHillFormula = NULL; /*  the length is unknown */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bIsotopic)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (num_isotopic_atoms &&
    // INCHI✔️✔️:             ( pINChI->IsotopicAtom = (INChI_IsotopicAtom *) inchi_calloc( num_isotopic_atoms, sizeof( INChI_IsotopicAtom ) ) ) &&
    // INCHI✔️✔️:              ( pINChI->IsotopicTGroup = (INChI_IsotopicTGroup *) inchi_calloc( num_isotopic_atoms, sizeof( INChI_IsotopicTGroup ) ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (num_isotopic_atoms)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (!( pINChI->nPossibleLocationsOfIsotopicH = (AT_NUMB *) inchi_calloc( (long long)num_at + 1, sizeof( pINChI->nPossibleLocationsOfIsotopicH[0] ) ) )) /* djb-rwth: cast operator added */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (( pINChI->Stereo = Alloc_INChI_Stereo( num_at, num_bonds ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto out_of_RAM;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bIsotopic)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (( pINChI->StereoIsotopic = Alloc_INChI_Stereo( num_at, num_bonds ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return pINChI;
    // INCHI✔️✔️:
    // INCHI✔️✔️: out_of_RAM:
    // INCHI✔️✔️:     if (pINChI)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         Free_INChI( &pINChI );
    // INCHI✔️✔️:         /*
    // INCHI✔️✔️:         inchi_free(pINChI);
    // INCHI✔️✔️:         */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return NULL;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: Alloc_INChI

    fn allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }

    if num_atoms <= 0 {
        return Ok(SourceMutPointer::null());
    }
    let mut inchi = match inchi_calloc::<INChI>(heap, 1, std::mem::size_of::<INChI>() as u64) {
        Ok(pointer) => pointer,
        Err(error) if allocation_error(error) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let processing = (|| {
        let count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let selected = atoms
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut num_bonds = 0_i32;
        let mut num_isotopic_atoms = 0_i32;
        for atom in selected {
            num_bonds = num_bonds.wrapping_add(i32::from(atom.valence));
            let is_deuterium = atom.elname[0] == b'D' as i8 && atom.elname[1] == 0;
            let is_tritium = atom.elname[0] == b'T' as i8 && atom.elname[1] == 0;
            if atom.iso_atw_diff != 0
                || is_deuterium
                || is_tritium
                || atom.num_iso_H.iter().any(|value| *value != 0)
            {
                num_isotopic_atoms = num_isotopic_atoms.wrapping_add(1);
            }
        }
        num_bonds /= 2;
        *found_num_bonds = num_bonds;
        *found_num_isotopic = num_isotopic_atoms;

        let atom_count = num_atoms as u64;
        let connection_count = i64::from(num_atoms).wrapping_add(i64::from(num_bonds)) as u64;
        let tautomer_count = ((i64::from(3 + INCHI_T_NUM_MOVABLE as i32) * i64::from(num_atoms))
            / 2)
        .wrapping_add(1) as u64;
        let pointer =
            inchi_calloc::<U_CHAR>(heap, atom_count, std::mem::size_of::<U_CHAR>() as u64)?;
        heap.slice_mut(inchi)?[0].nAtom = pointer;
        let pointer = inchi_calloc::<AT_NUMB>(
            heap,
            connection_count,
            std::mem::size_of::<AT_NUMB>() as u64,
        )?;
        heap.slice_mut(inchi)?[0].nConnTable = pointer;
        let pointer =
            inchi_calloc::<AT_NUMB>(heap, tautomer_count, std::mem::size_of::<AT_NUMB>() as u64)?;
        heap.slice_mut(inchi)?[0].nTautomer = pointer;
        let pointer =
            inchi_calloc::<S_CHAR>(heap, atom_count, std::mem::size_of::<S_CHAR>() as u64)?;
        heap.slice_mut(inchi)?[0].nNum_H = pointer;
        let pointer =
            inchi_calloc::<S_CHAR>(heap, atom_count, std::mem::size_of::<S_CHAR>() as u64)?;
        heap.slice_mut(inchi)?[0].nNum_H_fixed = pointer;
        heap.slice_mut(inchi)?[0].szHillFormula = SourceMutPointer::null();

        let isotopic = allocation_mode & REQ_MODE_ISO as i32;
        if isotopic != 0 {
            if num_isotopic_atoms != 0 {
                let isotope_count = num_isotopic_atoms as i64 as u64;
                let pointer = inchi_calloc::<INChI_IsotopicAtom>(
                    heap,
                    isotope_count,
                    std::mem::size_of::<INChI_IsotopicAtom>() as u64,
                )?;
                heap.slice_mut(inchi)?[0].IsotopicAtom = pointer;
                let pointer = inchi_calloc::<INChI_IsotopicTGroup>(
                    heap,
                    isotope_count,
                    std::mem::size_of::<INChI_IsotopicTGroup>() as u64,
                )?;
                heap.slice_mut(inchi)?[0].IsotopicTGroup = pointer;
            }
            let pointer = inchi_calloc::<AT_NUMB>(
                heap,
                i64::from(num_atoms).wrapping_add(1) as u64,
                std::mem::size_of::<AT_NUMB>() as u64,
            )?;
            heap.slice_mut(inchi)?[0].nPossibleLocationsOfIsotopicH = pointer;
        }

        let stereo = Alloc_INChI_Stereo(heap, num_atoms, num_bonds)?;
        if stereo.is_null() {
            return Err(SourceHeapError::AllocationFailed);
        }
        heap.slice_mut(inchi)?[0].Stereo = stereo;
        if isotopic != 0 {
            let stereo_isotopic = Alloc_INChI_Stereo(heap, num_atoms, num_bonds)?;
            if stereo_isotopic.is_null() {
                return Err(SourceHeapError::AllocationFailed);
            }
            heap.slice_mut(inchi)?[0].StereoIsotopic = stereo_isotopic;
        }
        Ok(())
    })();

    match processing {
        Ok(()) => Ok(inchi),
        Err(error) => {
            Free_INChI(heap, &mut inchi)?;
            if allocation_error(error) {
                Ok(SourceMutPointer::null())
            } else {
                Err(error)
            }
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn Free_INChI(
    heap: &mut SourceHeap,
    pp_inchi: &mut SourceMutPointer<INChI>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4675 Free_INChI
    // INCHI✔️❌: int Free_INChI( INChI **ppINChI )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     INChI *pINChI;
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((pINChI = *ppINChI)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bREUSE_INCHI == 1 )
    // INCHI✔️❌:         if (pINChI->nRefCount-- > 0)
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         Free_INChI_Members( pINChI );
    // INCHI✔️❌:         qzfree( pINChI );
    // INCHI✔️❌:         *ppINChI = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Free_INChI
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:589 bREUSE_INCHI
    // INCHI✔️❌: #define bREUSE_INCHI                1  /* 1=> do not recalulate INChI for components in reconnected
    // INCHI✔️❌:                                         *     structure that are same as in the connected one */
    // END INCHI ACTIVE MACRO: bREUSE_INCHI

    let p_inchi = *pp_inchi;
    if p_inchi.is_null() {
        return Ok(0);
    }
    let old_reference_count = heap.slice(p_inchi.as_const())?[0].nRefCount;
    heap.slice_mut(p_inchi)?[0].nRefCount = old_reference_count.wrapping_sub(1);
    if old_reference_count > 0 {
        return Ok(1);
    }
    Free_INChI_Members(heap, p_inchi)?;
    inchi_free(heap, p_inchi)?;
    *pp_inchi = SourceMutPointer::null();
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn Free_INChI_Members(
    heap: &mut SourceHeap,
    p_inchi: SourceMutPointer<INChI>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4698 Free_INChI_Members
    // INCHI✔❌: int Free_INChI_Members( INChI *pINChI )
    // INCHI✔❌: {
    // INCHI✔❌:     if (pINChI)
    // INCHI✔❌:     {
    // INCHI✔❌:         Free_INChI_Stereo(pINChI->Stereo);
    // INCHI✔❌:         Free_INChI_Stereo(pINChI->StereoIsotopic);
    // INCHI✔❌:         qzfree(pINChI->nAtom);
    // INCHI✔❌:         qzfree(pINChI->nConnTable);
    // INCHI✔❌:         qzfree(pINChI->nTautomer);
    // INCHI✔❌:         qzfree(pINChI->nNum_H);
    // INCHI✔❌:         qzfree(pINChI->nNum_H_fixed);
    // INCHI✔❌:         qzfree(pINChI->IsotopicAtom);
    // INCHI✔❌:         qzfree(pINChI->IsotopicTGroup);
    // INCHI✔❌:         qzfree(pINChI->nPossibleLocationsOfIsotopicH);
    // INCHI✔❌:         qzfree( pINChI->Stereo );
    // INCHI✔❌:         qzfree( pINChI->StereoIsotopic );
    // INCHI✔❌:         qzfree( pINChI->szHillFormula );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: Free_INChI_Members

    if p_inchi.is_null() {
        return Ok(0);
    }
    let stereo = heap.slice(p_inchi.as_const())?[0].Stereo;
    Free_INChI_Stereo(heap, stereo)?;
    let stereo_isotopic = heap.slice(p_inchi.as_const())?[0].StereoIsotopic;
    Free_INChI_Stereo(heap, stereo_isotopic)?;

    let pointer = heap.slice(p_inchi.as_const())?[0].nAtom;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].nAtom = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].nConnTable;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].nConnTable = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].nTautomer;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].nTautomer = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].nNum_H;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].nNum_H = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].nNum_H_fixed;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].nNum_H_fixed = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].IsotopicAtom;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].IsotopicAtom = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].IsotopicTGroup;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].IsotopicTGroup = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].nPossibleLocationsOfIsotopicH;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].nPossibleLocationsOfIsotopicH = SourceMutPointer::null();
    }
    if !stereo.is_null() {
        inchi_free(heap, stereo)?;
        heap.slice_mut(p_inchi)?[0].Stereo = SourceMutPointer::null();
    }
    if !stereo_isotopic.is_null() {
        inchi_free(heap, stereo_isotopic)?;
        heap.slice_mut(p_inchi)?[0].StereoIsotopic = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi.as_const())?[0].szHillFormula;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi)?[0].szHillFormula = SourceMutPointer::null();
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn Free_INChI_Aux(
    heap: &mut SourceHeap,
    pp_inchi_aux: &mut SourceMutPointer<INChI_Aux>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4843 Free_INChI_Aux
    // INCHI✔️❌: int Free_INChI_Aux( INChI_Aux **ppINChI_Aux )
    // INCHI✔️❌: {
    // INCHI✔️❌:     INChI_Aux *pINChI_Aux = *ppINChI_Aux;
    // INCHI✔️❌:     if (pINChI_Aux)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bREUSE_INCHI == 1 )
    // INCHI✔️❌:         if (pINChI_Aux->nRefCount-- > 0)
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         qzfree( pINChI_Aux->nOrigAtNosInCanonOrd );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nOrigAtNosInCanonOrdInv );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv );
    // INCHI✔️❌:         qzfree( pINChI_Aux->szOrigCoord );
    // INCHI✔️❌:         qzfree( pINChI_Aux->OrigInfo );
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         qzfree( pINChI_Aux->nOriginalAtomNumber          );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nCanonicalTGroupNumbers      );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nIsotopicCanonicalTGroupNumbers);
    // INCHI✔️❌:         qzfree( pINChI_Aux->nTautomer                    );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nNontautomericCanonicalNumbers         );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nIsotopicCanonicalNumbers    );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nNontautomericIsotopicCanonicalNumbers );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nNontautomericEquNumbers               );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nNontautomericIsotopicEquNumbers       );
    // INCHI✔️❌:         */
    // INCHI✔️❌:         qzfree( pINChI_Aux->nConstitEquNumbers );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nConstitEquTGroupNumbers );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nConstitEquIsotopicNumbers );
    // INCHI✔️❌:         qzfree( pINChI_Aux->nConstitEquIsotopicTGroupNumbers );
    // INCHI✔️❌:         qzfree( pINChI_Aux );
    // INCHI✔️❌:         *ppINChI_Aux = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Free_INChI_Aux
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:589 bREUSE_INCHI
    // INCHI✔️❌: #define bREUSE_INCHI                1  /* 1=> do not recalulate INChI for components in reconnected
    // INCHI✔️❌:                                         *     structure that are same as in the connected one */
    // END INCHI ACTIVE MACRO: bREUSE_INCHI

    let p_inchi_aux = *pp_inchi_aux;
    if p_inchi_aux.is_null() {
        return Ok(0);
    }
    let old_reference_count = heap.slice(p_inchi_aux.as_const())?[0].nRefCount;
    heap.slice_mut(p_inchi_aux)?[0].nRefCount = old_reference_count.wrapping_sub(1);
    if old_reference_count > 0 {
        return Ok(1);
    }

    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nOrigAtNosInCanonOrd;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nOrigAtNosInCanonOrd = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nIsotopicOrigAtNosInCanonOrd;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nIsotopicOrigAtNosInCanonOrd = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nOrigAtNosInCanonOrdInv;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nOrigAtNosInCanonOrdInv = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nIsotopicOrigAtNosInCanonOrdInv;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nIsotopicOrigAtNosInCanonOrdInv = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].szOrigCoord;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].szOrigCoord = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].OrigInfo;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].OrigInfo = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nConstitEquNumbers;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nConstitEquNumbers = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nConstitEquTGroupNumbers;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nConstitEquTGroupNumbers = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nConstitEquIsotopicNumbers;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nConstitEquIsotopicNumbers = SourceMutPointer::null();
    }
    let pointer = heap.slice(p_inchi_aux.as_const())?[0].nConstitEquIsotopicTGroupNumbers;
    if !pointer.is_null() {
        inchi_free(heap, pointer)?;
        heap.slice_mut(p_inchi_aux)?[0].nConstitEquIsotopicTGroupNumbers = SourceMutPointer::null();
    }
    inchi_free(heap, p_inchi_aux)?;
    *pp_inchi_aux = SourceMutPointer::null();
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn Alloc_INChI_Aux(
    heap: &mut SourceHeap,
    num_atoms: i32,
    num_isotopic_atoms: i32,
    allocation_mode: i32,
    original_coordinates: i32,
) -> Result<SourceMutPointer<INChI_Aux>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4884 Alloc_INChI_Aux
    // INCHI✔️✔️: INChI_Aux *Alloc_INChI_Aux( int num_at,
    // INCHI✔️✔️:                             int num_isotopic_atoms,
    // INCHI✔️✔️:                             int nAllocMode,
    // INCHI✔️✔️:                             int bOrigCoord )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     INChI_Aux     *pINChI_Aux;
    // INCHI✔️✔️:     int    bIsotopic = ( nAllocMode & REQ_MODE_ISO );
    // INCHI✔️✔️:     int    num_at_tg = num_at + num_at / 2;
    // INCHI✔️✔️:     /* int    bTautomeric = (nAllocMode & REQ_MODE_TAUT); */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_at <= 0 ||
    // INCHI✔️✔️:          NULL == ( pINChI_Aux = (INChI_Aux *) inchi_calloc( 1, sizeof( INChI_Aux ) ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (( pINChI_Aux->nOrigAtNosInCanonOrd = (AT_NUMB*)
    // INCHI✔️✔️:           inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nOrigAtNosInCanonOrd[0] ) ) ) &&
    // INCHI✔️✔️:           ( pINChI_Aux->nOrigAtNosInCanonOrdInv = (AT_NUMB*)
    // INCHI✔️✔️:             inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nOrigAtNosInCanonOrd[0] ) ) ) &&
    // INCHI✔️✔️:             ( pINChI_Aux->nConstitEquNumbers = (AT_NUMB*)
    // INCHI✔️✔️:               inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nConstitEquNumbers[0] ) ) ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto out_of_RAM;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_at > 1 &&
    // INCHI✔️✔️:         ( pINChI_Aux->nConstitEquTGroupNumbers = (AT_NUMB*) inchi_calloc( (long long)num_at / 2 + 1, sizeof( pINChI_Aux->nConstitEquTGroupNumbers[0] ) ) )) /* djb-rwth: cast operator added */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (num_at > 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_at > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         pINChI_Aux->OrigInfo = (ORIG_INFO *) inchi_calloc( num_at, sizeof( pINChI_Aux->OrigInfo[0] ) );
    // INCHI✔️✔️:         if (!pINChI_Aux->OrigInfo)
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bOrigCoord && num_at > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         pINChI_Aux->szOrigCoord = (MOL_COORD *) inchi_calloc( num_at, sizeof( pINChI_Aux->szOrigCoord[0] ) );
    // INCHI✔️✔️:         if (!pINChI_Aux->szOrigCoord)
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bIsotopic)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if ( /*num_isotopic_atoms &&*/
    // INCHI✔️✔️:             ( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd = (AT_NUMB*) inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0] ) ) ) &&
    // INCHI✔️✔️:              ( pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv = (AT_NUMB*) inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0] ) ) ) &&
    // INCHI✔️✔️:              ( pINChI_Aux->nConstitEquIsotopicNumbers = (AT_NUMB*) inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nConstitEquIsotopicNumbers[0] ) ) ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (num_isotopic_atoms)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if ( /*num_isotopic_atoms && num_at > 1 &&*/
    // INCHI✔️✔️:             ( pINChI_Aux->nConstitEquIsotopicTGroupNumbers = (AT_NUMB*) inchi_calloc( (long long)num_at / 2 + 1, sizeof( pINChI_Aux->nConstitEquIsotopicTGroupNumbers[0] ) ) )) /* djb-rwth: cast operator added */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else if (num_isotopic_atoms && num_at > 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             goto out_of_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return pINChI_Aux;
    // INCHI✔️✔️:
    // INCHI✔️✔️: out_of_RAM:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (pINChI_Aux)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         Free_INChI_Aux( &pINChI_Aux );
    // INCHI✔️✔️:         /*
    // INCHI✔️✔️:         inchi_free(pINChI_Aux);
    // INCHI✔️✔️:         */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return NULL;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: Alloc_INChI_Aux

    fn allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }
    if num_atoms <= 0 {
        return Ok(SourceMutPointer::null());
    }
    let mut aux = match inchi_calloc::<INChI_Aux>(heap, 1, std::mem::size_of::<INChI_Aux>() as u64)
    {
        Ok(pointer) => pointer,
        Err(error) if allocation_error(error) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let allocation = (|| {
        let atom_count = num_atoms as u64;
        let atom_tgroup_count = num_atoms.wrapping_add(num_atoms / 2) as u64;
        let tgroup_count = (i64::from(num_atoms) / 2).wrapping_add(1) as u64;
        let pointer = inchi_calloc::<AT_NUMB>(heap, atom_tgroup_count, 2)?;
        heap.slice_mut(aux)?[0].nOrigAtNosInCanonOrd = pointer;
        let pointer = inchi_calloc::<AT_NUMB>(heap, atom_tgroup_count, 2)?;
        heap.slice_mut(aux)?[0].nOrigAtNosInCanonOrdInv = pointer;
        let pointer = inchi_calloc::<AT_NUMB>(heap, atom_tgroup_count, 2)?;
        heap.slice_mut(aux)?[0].nConstitEquNumbers = pointer;
        if num_atoms > 1 {
            let pointer = inchi_calloc::<AT_NUMB>(heap, tgroup_count, 2)?;
            heap.slice_mut(aux)?[0].nConstitEquTGroupNumbers = pointer;
        }
        let pointer =
            inchi_calloc::<ORIG_INFO>(heap, atom_count, std::mem::size_of::<ORIG_INFO>() as u64)?;
        heap.slice_mut(aux)?[0].OrigInfo = pointer;
        if original_coordinates != 0 {
            let pointer = inchi_calloc::<MOL_COORD>(
                heap,
                atom_count,
                std::mem::size_of::<MOL_COORD>() as u64,
            )?;
            heap.slice_mut(aux)?[0].szOrigCoord = pointer;
        }
        if allocation_mode & REQ_MODE_ISO as i32 != 0 {
            match inchi_calloc::<AT_NUMB>(heap, atom_tgroup_count, 2) {
                Ok(pointer) => {
                    heap.slice_mut(aux)?[0].nIsotopicOrigAtNosInCanonOrd = pointer;
                    match inchi_calloc::<AT_NUMB>(heap, atom_tgroup_count, 2) {
                        Ok(pointer) => {
                            heap.slice_mut(aux)?[0].nIsotopicOrigAtNosInCanonOrdInv = pointer;
                            match inchi_calloc::<AT_NUMB>(heap, atom_tgroup_count, 2) {
                                Ok(pointer) => {
                                    heap.slice_mut(aux)?[0].nConstitEquIsotopicNumbers = pointer
                                }
                                Err(error)
                                    if num_isotopic_atoms == 0 && allocation_error(error) => {}
                                Err(error) => return Err(error),
                            }
                        }
                        Err(error) if num_isotopic_atoms == 0 && allocation_error(error) => {}
                        Err(error) => return Err(error),
                    }
                }
                Err(error) if num_isotopic_atoms == 0 && allocation_error(error) => {}
                Err(error) => return Err(error),
            }
            match inchi_calloc::<AT_NUMB>(heap, tgroup_count, 2) {
                Ok(pointer) => heap.slice_mut(aux)?[0].nConstitEquIsotopicTGroupNumbers = pointer,
                Err(error)
                    if !(num_isotopic_atoms != 0 && num_atoms > 1) && allocation_error(error) => {}
                Err(error) => return Err(error),
            }
        }
        Ok(())
    })();
    match allocation {
        Ok(()) => Ok(aux),
        Err(error) => {
            Free_INChI_Aux(heap, &mut aux)?;
            if allocation_error(error) {
                Ok(SourceMutPointer::null())
            } else {
                Err(error)
            }
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn imat_free(
    heap: &mut SourceHeap,
    m: i32,
    mut a: SourceMutPointer<SourceMutPointer<i32>>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5036 imat_free
    // INCHI✔❌: void imat_free( int m, int **a )
    // INCHI✔❌: {
    // INCHI✔❌:     int i;
    // INCHI✔❌:     if (NULL != a)
    // INCHI✔❌:     {
    // INCHI✔❌:         for (i = 0; i < m; i++)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (NULL != a[i]) /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free( a[i] );
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         inchi_free( a );
    // INCHI✔❌:         a = NULL;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: imat_free

    if !a.is_null() {
        for i in 0..m {
            let row = heap.slice(a.as_const())?[i as usize];
            if !row.is_null() {
                inchi_free(heap, row)?;
            }
        }
        inchi_free(heap, a)?;
        a = SourceMutPointer::null();
    }
    let _ = a;
    Ok(())
}

pub(crate) fn subgraf_free(
    heap: &mut SourceHeap,
    graph: SourceMutPointer<subgraf>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5155 subgraf_free
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: subgraf_free
    // INCHI✔️❌: void subgraf_free( subgraf *sg )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     if (!sg)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (sg->nodes)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( sg->nodes );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (sg->degrees)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( sg->degrees );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (sg->orig2node)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( sg->orig2node );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (sg->adj)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < sg->nnodes; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (sg->adj[i]) /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( sg->adj[i] );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_free( sg->adj );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_free( sg );
    // INCHI✔️❌:     sg = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: subgraf_free
    // END INCHI C FUNCTION: subgraf_free
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_free
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; inchi_free expands to free.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_free

    if graph.is_null() {
        return Ok(());
    }
    let graph_value = heap
        .slice(graph.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    inchi_free(heap, graph_value.nodes)?;
    inchi_free(heap, graph_value.degrees)?;
    inchi_free(heap, graph_value.orig2node)?;
    if !graph_value.adj.is_null() {
        let row_count = usize::try_from(graph_value.nnodes.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let rows = heap.slice(graph_value.adj.as_const())?;
        if rows.len() < row_count {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let rows = rows[..row_count].to_vec();
        for row in rows {
            inchi_free(heap, row)?;
        }
        inchi_free(heap, graph_value.adj)?;
    }
    inchi_free(heap, graph)
}

pub(crate) fn subgraf_new(
    heap: &mut SourceHeap,
    original_input: &ORIG_ATOM_DATA,
    node_count: i32,
    input_nodes: &[i32],
) -> Result<SourceMutPointer<subgraf>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5064 subgraf_new
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: subgraf_new
    // INCHI✔️❌: subgraf *subgraf_new( ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                       int nnodes,
    // INCHI✔️❌:                       int *nodes )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, iat, nbr, jat, nj, degree, nat, err = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     subgraf *sg = (subgraf *) inchi_calloc( 1, sizeof( subgraf ) );
    // INCHI✔️❌:     if (!sg)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nat = orig_inp_data->num_inp_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* orig2node is mapping of original at numbers --> subgraph node numbers */
    // INCHI✔️❌:     err = 1;
    // INCHI✔️❌:     if (!( sg->orig2node = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ) )) /* djb-rwth: cast operator added */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!( sg->nodes = (int *) inchi_calloc( nnodes, sizeof( int ) ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!( sg->degrees = (int *) inchi_calloc( nnodes, sizeof( int ) ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* NB:    input list of 'nodes' is assumed to be in 'original_atom_numbering domain' which starts from 1.
    // INCHI✔️❌:     Now it is mapped to current atom numbers which starts from 0/connections using at[j].orig_at_number    */
    // INCHI✔️❌:     sg->nnodes = 0;
    // INCHI✔️❌:     for (i = 0; i < nnodes; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sg->nodes[sg->nnodes++] = nodes[i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i <= nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sg->orig2node[i] = -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < nnodes; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sg->orig2node[sg->nodes[i]] = i;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Create and fill subgraph adjacency matrix based on nodes/orig atom numbers
    // INCHI✔️❌:        and connections stored in orig_inp_data */
    // INCHI✔️❌:     sg->adj = (subgraf_edge **) inchi_calloc( nnodes, sizeof( subgraf_edge * ) );
    // INCHI✔️❌:     if (!sg->adj)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < sg->nnodes; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         iat = nodes[i] - 1;    /* current atom  number for this node */
    // INCHI✔️❌:         degree = orig_inp_data->at[iat].valence;
    // INCHI✔️❌:         nj = -1;
    // INCHI✔️❌:         sg->adj[i] = (subgraf_edge *) inchi_calloc( degree, sizeof( subgraf_edge ) );
    // INCHI✔️❌:         if (!sg->adj[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (j = 0; j < degree; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             jat = orig_inp_data->at[iat].neighbor[j]; /* for curr num jat, a (jat+1) would be an orig num */
    // INCHI✔️❌:             nbr = sg->orig2node[jat + 1];
    // INCHI✔️❌:             if (nbr < 0)
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             nj++;
    // INCHI✔️❌:             sg->adj[i][nj].nbr = nbr;
    // INCHI✔️❌:             sg->adj[i][nj].etype = orig_inp_data->at[iat].bond_type[j];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sg->degrees[i] = nj + 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     err = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* subgraf_debug_trace( sg ); */
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         subgraf_free( sg );
    // INCHI✔️❌:         return NULL; /* djb-rwth: avoiding reading from freed memory */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return sg;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: subgraf_new
    // END INCHI C FUNCTION: subgraf_new
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_new
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(int)=4, sizeof(void*)=8, sizeof(subgraf)=40, sizeof(subgraf_edge)=8.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_new

    fn is_allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }

    let graph = match inchi_calloc::<subgraf>(heap, 1, 40) {
        Ok(pointer) => pointer,
        Err(error) if is_allocation_error(error) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let build = (|| -> Result<(), SourceHeapError> {
        let atom_count = original_input.num_inp_atoms;
        let map_count = i64::from(atom_count).wrapping_add(1) as u64;
        let orig2node = inchi_calloc::<i32>(heap, map_count, 4)?;
        heap.slice_mut(graph)?[0].orig2node = orig2node;

        let source_node_count = node_count as u64;
        let nodes = inchi_calloc::<i32>(heap, source_node_count, 4)?;
        heap.slice_mut(graph)?[0].nodes = nodes;
        let degrees = inchi_calloc::<i32>(heap, source_node_count, 4)?;
        heap.slice_mut(graph)?[0].degrees = degrees;

        let count =
            usize::try_from(node_count.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if input_nodes.len() < count {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        for (index, node) in input_nodes[..count].iter().copied().enumerate() {
            heap.slice_mut(nodes)?[index] = node;
            heap.slice_mut(graph)?[0].nnodes = (index as i32).wrapping_add(1);
        }

        let map_length = usize::try_from(atom_count.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let map = heap.slice_mut(orig2node)?;
        if map.len() < map_length {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        map[..map_length].fill(-1);
        for (index, node) in input_nodes[..count].iter().copied().enumerate() {
            let node = usize::try_from(node).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *map.get_mut(node)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = index as i32;
        }

        let adjacency = inchi_calloc::<SourceMutPointer<subgraf_edge>>(heap, source_node_count, 8)?;
        heap.slice_mut(graph)?[0].adj = adjacency;
        let atoms = if original_input.at.is_null() {
            Vec::new()
        } else {
            heap.slice(original_input.at.as_const())?.to_vec()
        };
        for index in 0..count {
            let atom_index = usize::try_from(input_nodes[index].wrapping_sub(1))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = atoms
                .get(atom_index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let source_degree = i32::from(atom.valence);
            let row = inchi_calloc::<subgraf_edge>(heap, source_degree as u64, 8)?;
            heap.slice_mut(adjacency)?[index] = row;
            let degree = usize::try_from(source_degree.max(0))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut included = 0_usize;
            for neighbor_index in 0..degree {
                let original_neighbor = usize::from(atom.neighbor[neighbor_index]) + 1;
                let neighbor = *heap
                    .slice(orig2node.as_const())?
                    .get(original_neighbor)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if neighbor < 0 {
                    continue;
                }
                heap.slice_mut(row)?[included] = subgraf_edge {
                    nbr: neighbor,
                    etype: i32::from(atom.bond_type[neighbor_index]),
                };
                included += 1;
            }
            heap.slice_mut(degrees)?[index] = included as i32;
        }
        Ok(())
    })();
    match build {
        Ok(()) => Ok(graph),
        Err(error) => {
            subgraf_free(heap, graph)?;
            if is_allocation_error(error) {
                Ok(SourceMutPointer::null())
            } else {
                Err(error)
            }
        }
    }
}

pub(crate) fn subgraf_pathfinder_new(
    heap: &mut SourceHeap,
    graph: SourceMutPointer<subgraf>,
    _original_input: &ORIG_ATOM_DATA,
    start: i32,
    end: i32,
) -> Result<SourceMutPointer<subgraf_pathfinder>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5221 subgraf_pathfinder_new
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: subgraf_pathfinder_new
    // INCHI✔️❌: subgraf_pathfinder * subgraf_pathfinder_new( subgraf *sg,
    // INCHI✔️❌:                                              ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                                              int start,
    // INCHI✔️❌:                                              int end )
    // INCHI✔️❌: {
    // INCHI✔️❌:     subgraf_pathfinder *spf = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     spf = (subgraf_pathfinder *) inchi_calloc( 1, sizeof( subgraf_pathfinder ) );
    // INCHI✔️❌:     if (!spf)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     spf->sg = sg;
    // INCHI✔️❌:     spf->start = start;
    // INCHI✔️❌:     spf->end = end;
    // INCHI✔️❌:     spf->nbonds = 0;
    // INCHI✔️❌:     spf->nseen = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     spf->seen = (int *) inchi_calloc( spf->sg->nnodes, sizeof( int ) );
    // INCHI✔️❌:     if (!spf->seen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( spf );
    // INCHI✔️❌:         spf = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     return spf;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: subgraf_pathfinder_new
    // END INCHI C FUNCTION: subgraf_pathfinder_new
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_pathfinder_new
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(subgraf_pathfinder)=40, sizeof(int)=4.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_pathfinder_new

    fn is_allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }
    let pathfinder = match inchi_calloc::<subgraf_pathfinder>(heap, 1, 40) {
        Ok(pointer) => pointer,
        Err(error) if is_allocation_error(error) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let node_count = heap
        .slice(graph.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .nnodes;
    {
        let value = &mut heap.slice_mut(pathfinder)?[0];
        value.sg = graph;
        value.start = start;
        value.end = end;
        value.nbonds = 0;
        value.nseen = 0;
    }
    let seen = match inchi_calloc::<i32>(heap, node_count as u64, 4) {
        Ok(pointer) => pointer,
        Err(error) if is_allocation_error(error) => {
            inchi_free(heap, pathfinder)?;
            return Ok(SourceMutPointer::null());
        }
        Err(error) => return Err(error),
    };
    heap.slice_mut(pathfinder)?[0].seen = seen;
    Ok(pathfinder)
}

pub(crate) fn subgraf_pathfinder_free(
    heap: &mut SourceHeap,
    pathfinder: SourceMutPointer<subgraf_pathfinder>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5252 subgraf_pathfinder_free
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: subgraf_pathfinder_free
    // INCHI✔️✔️: void subgraf_pathfinder_free( subgraf_pathfinder *spf )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (!spf)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (spf->seen)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free( spf->seen );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     inchi_free( spf );
    // INCHI✔️✔️:     spf = NULL;
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: subgraf_pathfinder_free
    // END INCHI C FUNCTION: subgraf_pathfinder_free
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_pathfinder_free
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; inchi_free expands to free.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_pathfinder_free

    if pathfinder.is_null() {
        return Ok(());
    }
    let seen = heap
        .slice(pathfinder.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .seen;
    inchi_free(heap, seen)?;
    inchi_free(heap, pathfinder)
}

pub(crate) fn add_bond_if_unseen(
    heap: &mut SourceHeap,
    pathfinder: SourceMutPointer<subgraf_pathfinder>,
    node0: i32,
    node: i32,
    bond_count: &mut i32,
    bonds: SourceMutPointer<SourceMutPointer<i32>>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5380 add_bond_if_unseen
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: add_bond_if_unseen
    // INCHI✔️✔️: void add_bond_if_unseen( subgraf_pathfinder *spf,
    // INCHI✔️✔️:                          int node0,
    // INCHI✔️✔️:                          int node,
    // INCHI✔️✔️:                          int *nbonds,
    // INCHI✔️✔️:                          int **bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int seen, p, at1, at2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     at1 = spf->sg->nodes[node0];
    // INCHI✔️✔️:     at2 = spf->sg->nodes[node];
    // INCHI✔️✔️: #if 0
    // INCHI✔️✔️:     if (at1 > at2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int tmp = at1;
    // INCHI✔️✔️:         at1 = at2;
    // INCHI✔️✔️:         at2 = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     seen = 0;
    // INCHI✔️✔️:     for (p = 0; p < *nbonds; p++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*if (bonds[p][0] == at1 && bonds[p][1] == at2)*/
    // INCHI✔️✔️:         if (bIsSameBond(at1, at2, bonds[p][0], bonds[p][1]))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             seen = 1;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (!seen)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         bonds[*nbonds][0] = at1;
    // INCHI✔️✔️:         bonds[*nbonds][1] = at2;
    // INCHI✔️✔️:         ( *nbonds )++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: add_bond_if_unseen
    // END INCHI C FUNCTION: add_bond_if_unseen
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: add_bond_if_unseen
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; the #if 0 endpoint-sorting block is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: add_bond_if_unseen

    let graph = heap
        .slice(pathfinder.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .sg;
    let nodes = heap
        .slice(graph.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .nodes;
    let node0 = usize::try_from(node0).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let node = usize::try_from(node).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let source_nodes = heap.slice(nodes.as_const())?;
    let atom1 = *source_nodes
        .get(node0)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atom2 = *source_nodes
        .get(node)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let count = usize::try_from(*bond_count).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for index in 0..count {
        let row = *heap
            .slice(bonds.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let values = heap.slice(row.as_const())?;
        let existing1 = *values.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let existing2 = *values.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if bIsSameBond(atom1, atom2, existing1, existing2) != 0 {
            return Ok(());
        }
    }

    let row = *heap
        .slice(bonds.as_const())?
        .get(count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(row)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = atom1;
    *heap
        .slice_mut(row)?
        .get_mut(1)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = atom2;
    *bond_count = bond_count.wrapping_add(1);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn subgraf_pathfinder_run(
    heap: &mut SourceHeap,
    pathfinder: SourceMutPointer<subgraf_pathfinder>,
    forbidden_count: i32,
    forbidden: SourceMutPointer<i32>,
    bond_count: &mut i32,
    bonds: SourceMutPointer<SourceMutPointer<i32>>,
    atom_count: &mut i32,
    atoms: SourceMutPointer<i32>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5273 subgraf_pathfinder_run
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: subgraf_pathfinder_run
    // INCHI✔️✔️: void subgraf_pathfinder_run( subgraf_pathfinder *spf,
    // INCHI✔️✔️:                              int nforbidden,		/* number of edges forbidden for traversal	*/
    // INCHI✔️✔️:                              int *forbidden,		/* nodes of forbidden edges: [edge1node1,edge1node2, edge2node1,edge2node2, ... ] */
    // INCHI✔️✔️:                              int *nbonds,
    // INCHI✔️✔️:                              int **bonds,			/* collect subgraf bonds here */
    // INCHI✔️✔️:                              int *natoms,
    // INCHI✔️✔️:                              int *atoms				/* if not NULL, collect subgraf atoms here	*/
    // INCHI✔️✔️:                              )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int j, k, node, node0;
    // INCHI✔️✔️:     int f, skip;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (spf->nseen < 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*    Even at very beginning, push start node to seen and set nseen = 1
    // INCHI✔️✔️:         and put end node into subgraf_pathfinder's end                            */
    // INCHI✔️✔️:         return;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     node0 = spf->seen[spf->nseen - 1];
    // INCHI✔️✔️:     for (j = 0; j < spf->sg->degrees[node0]; j++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         node = spf->sg->adj[node0][j].nbr;
    // INCHI✔️✔️:         if (is_in_the_ilist( spf->seen, node, spf->nseen ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nforbidden && forbidden)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             skip = 0;
    // INCHI✔️✔️:             for (f = 0; f < nforbidden; f++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (bIsSameBond(node0, node, forbidden[2 * f], forbidden[2 * f + 1]) )
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     skip = 1;
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (skip)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (node == spf->end)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             spf->seen[spf->nseen++] = node;
    // INCHI✔️✔️:
    // INCHI✔️✔️:             ITRACE_( "\n\tFound path (in orig atom numbers):\t" );
    // INCHI✔️✔️:             for (k = 0; k < spf->nseen; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 int orig_atnum = spf->sg->nodes[spf->seen[k]];
    // INCHI✔️✔️:                 ITRACE_( "%-d ", orig_atnum);
    // INCHI✔️✔️:                 if (atoms && !is_in_the_ilist(atoms, orig_atnum, *natoms))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     atoms[(*natoms)++] = orig_atnum;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             ITRACE_( "\t( In node nums: " );
    // INCHI✔️✔️:             for (k = 1; k < spf->nseen; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 int at1 = spf->seen[k - 1];
    // INCHI✔️✔️:                 int at2 = spf->seen[k];
    // INCHI✔️✔️:                 add_bond_if_unseen( spf, at1, at2, nbonds, bonds );
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 ITRACE_( "%-d ", spf->seen[k] );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             ITRACE_( ")" );
    // INCHI✔️✔️:
    // INCHI✔️✔️:             spf->seen[spf->nseen - 1] = 0;
    // INCHI✔️✔️:             spf->nseen--;    /* pop_back        */
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (j = 0; j < spf->sg->degrees[node0]; j++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         node = spf->sg->adj[node0][j].nbr;
    // INCHI✔️✔️:         if (node == spf->end || is_in_the_ilist( spf->seen, node, spf->nseen ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nforbidden && forbidden)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             skip = 0;
    // INCHI✔️✔️:             for (f = 0; f < nforbidden; f++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (bIsSameBond(node0, node, forbidden[2 * f], forbidden[2 * f + 1]))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     skip = 1;
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (skip)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         spf->seen[spf->nseen++] = node;
    // INCHI✔️✔️:         subgraf_pathfinder_run( spf, 0, NULL, nbonds, bonds, natoms, atoms );
    // INCHI✔️✔️:         spf->seen[spf->nseen - 1] = 0;
    // INCHI✔️✔️:         spf->nseen--;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: subgraf_pathfinder_run
    // END INCHI C FUNCTION: subgraf_pathfinder_run
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_pathfinder_run
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; ITRACE_ is configured diagnostic tracing and does not mutate graph/output state.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: subgraf_pathfinder_run

    fn pathfinder_value(
        heap: &SourceHeap,
        pathfinder: SourceMutPointer<subgraf_pathfinder>,
    ) -> Result<subgraf_pathfinder, SourceHeapError> {
        heap.slice(pathfinder.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn graph_value(
        heap: &SourceHeap,
        graph: SourceMutPointer<subgraf>,
    ) -> Result<subgraf, SourceHeapError> {
        heap.slice(graph.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn neighbor_at(
        heap: &SourceHeap,
        graph: &subgraf,
        node: i32,
        edge: i32,
    ) -> Result<i32, SourceHeapError> {
        let node = usize::try_from(node).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let edge = usize::try_from(edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let row = *heap
            .slice(graph.adj.as_const())?
            .get(node)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        Ok(heap
            .slice(row.as_const())?
            .get(edge)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nbr)
    }
    fn is_forbidden(
        heap: &SourceHeap,
        count: i32,
        forbidden: SourceMutPointer<i32>,
        first: i32,
        second: i32,
    ) -> Result<bool, SourceHeapError> {
        if count == 0 || forbidden.is_null() {
            return Ok(false);
        }
        let values = heap.slice(forbidden.as_const())?;
        for index in 0..count.max(0) {
            let offset = usize::try_from(index)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                .checked_mul(2)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let edge_first = *values
                .get(offset)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let edge_second = *values
                .get(offset + 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if bIsSameBond(first, second, edge_first, edge_second) != 0 {
                return Ok(true);
            }
        }
        Ok(false)
    }

    let initial = pathfinder_value(heap, pathfinder)?;
    if initial.nseen < 1 {
        return Ok(());
    }
    let seen_index = usize::try_from(initial.nseen.wrapping_sub(1))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let node0 = *heap
        .slice(initial.seen.as_const())?
        .get(seen_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let graph = graph_value(heap, initial.sg)?;
    let node0_index = usize::try_from(node0).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let degree = *heap
        .slice(graph.degrees.as_const())?
        .get(node0_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    for edge in 0..degree.max(0) {
        let node = neighbor_at(heap, &graph, node0, edge)?;
        let current = pathfinder_value(heap, pathfinder)?;
        if is_in_the_ilist(
            Some(heap.slice(current.seen.as_const())?),
            node,
            current.nseen,
        )?
        .is_some()
        {
            continue;
        }
        if is_forbidden(heap, forbidden_count, forbidden, node0, node)? {
            continue;
        }
        if node == current.end {
            let append =
                usize::try_from(current.nseen).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(current.seen)?
                .get_mut(append)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = node;
            heap.slice_mut(pathfinder)?[0].nseen = current.nseen.wrapping_add(1);

            let found = pathfinder_value(heap, pathfinder)?;
            for path_index in 0..found.nseen.max(0) {
                let path_index =
                    usize::try_from(path_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let graph_node = *heap
                    .slice(found.seen.as_const())?
                    .get(path_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let graph_node =
                    usize::try_from(graph_node).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let original_atom = *heap
                    .slice(graph.nodes.as_const())?
                    .get(graph_node)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !atoms.is_null()
                    && is_in_the_ilist(
                        Some(heap.slice(atoms.as_const())?),
                        original_atom,
                        *atom_count,
                    )?
                    .is_none()
                {
                    let output_index = usize::try_from(*atom_count)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    *heap
                        .slice_mut(atoms)?
                        .get_mut(output_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = original_atom;
                    *atom_count = atom_count.wrapping_add(1);
                }
            }
            for path_index in 1..found.nseen.max(0) {
                let second =
                    usize::try_from(path_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let first_node = *heap
                    .slice(found.seen.as_const())?
                    .get(second - 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let second_node = *heap
                    .slice(found.seen.as_const())?
                    .get(second)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                add_bond_if_unseen(heap, pathfinder, first_node, second_node, bond_count, bonds)?;
            }
            let after = pathfinder_value(heap, pathfinder)?;
            let pop = usize::try_from(after.nseen.wrapping_sub(1))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(after.seen)?
                .get_mut(pop)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            heap.slice_mut(pathfinder)?[0].nseen = after.nseen.wrapping_sub(1);
            break;
        }
    }

    for edge in 0..degree.max(0) {
        let node = neighbor_at(heap, &graph, node0, edge)?;
        let current = pathfinder_value(heap, pathfinder)?;
        if node == current.end
            || is_in_the_ilist(
                Some(heap.slice(current.seen.as_const())?),
                node,
                current.nseen,
            )?
            .is_some()
        {
            continue;
        }
        if is_forbidden(heap, forbidden_count, forbidden, node0, node)? {
            continue;
        }
        let append =
            usize::try_from(current.nseen).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(current.seen)?
            .get_mut(append)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = node;
        heap.slice_mut(pathfinder)?[0].nseen = current.nseen.wrapping_add(1);
        subgraf_pathfinder_run(
            heap,
            pathfinder,
            0,
            SourceMutPointer::null(),
            bond_count,
            bonds,
            atom_count,
            atoms,
        )?;
        let after = pathfinder_value(heap, pathfinder)?;
        let pop = usize::try_from(after.nseen.wrapping_sub(1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(after.seen)?
            .get_mut(pop)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        heap.slice_mut(pathfinder)?[0].nseen = after.nseen.wrapping_sub(1);
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn subgraf_pathfinder_collect_all(
    heap: &mut SourceHeap,
    pathfinder: SourceMutPointer<subgraf_pathfinder>,
    forbidden_count: i32,
    forbidden: SourceMutPointer<i32>,
    atom_numbers: SourceMutPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5422 subgraf_pathfinder_collect_all
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int subgraf_pathfinder_collect_all(subgraf_pathfinder *spf,
                                   int nforbidden,		/* number of edges forbidden for traversal	*/
                                   int *forbidden,		/* nodes of forbidden edges: [edge1node1,edge1node2, edge2node1, edge2node2, ... ] */
                                   int *atnums          /* 1-based origs# */
                                    )
{
    int j, f, node, next_node, skip;

    node = spf->start;
    spf->seen[spf->nseen] = node;
    atnums[spf->nseen] = spf->sg->nodes[node];
    spf->nseen++;

    for (j = 0; j < spf->sg->degrees[node]; j++)
    {
        next_node = spf->sg->adj[node][j].nbr;
        if (is_in_the_ilist(spf->seen, next_node, spf->nseen))
        {
            continue;
        }
        if (nforbidden && forbidden)
        {
            skip = 0;
            for (f = 0; f < nforbidden; f++)
            {
                if (bIsSameBond(node, next_node, forbidden[2 * f], forbidden[2 * f + 1]))
                {
                    skip = 1;
                    break;
                }
            }
            if (skip)
            {
                continue;
            }
        }
        spf->start = next_node;
        subgraf_pathfinder_collect_all(spf, nforbidden, forbidden, atnums);
    }

    return spf->nseen;
}
    */
    // END INCHI C FUNCTION: subgraf_pathfinder_collect_all
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: subgraf_pathfinder_collect_all
    // INCHI✔️✔️: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this production helper.
    // INCHI✔️✔️: is_in_the_ilist and bIsSameBond use their already-ported source behavior.
    // INCHI✔️✔️: Rust preserves recursive DFS order, permanent seen entries, and the final start value.
    // END INCHI ACTIVE MACRO CONFIGURATION: subgraf_pathfinder_collect_all

    fn pointed_value<T: Clone + 'static>(
        heap: &SourceHeap,
        pointer: SourceMutPointer<T>,
    ) -> Result<T, SourceHeapError> {
        heap.slice(pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    let state = pointed_value(heap, pathfinder)?;
    let graph: subgraf = pointed_value(heap, state.sg)?;
    let node = state.start;
    let node_index = usize::try_from(node).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let seen_index = usize::try_from(state.nseen).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let original_atom_number = *heap
        .slice(graph.nodes.as_const())?
        .get(node_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(state.seen)?
        .get_mut(seen_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = node;
    *heap
        .slice_mut(atom_numbers)?
        .get_mut(seen_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = original_atom_number;
    let next_seen = state
        .nseen
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    heap.slice_mut(pathfinder)?[0].nseen = next_seen;

    let degree = *heap
        .slice(graph.degrees.as_const())?
        .get(node_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for edge_index in 0..degree.max(0) {
        let row = *heap
            .slice(graph.adj.as_const())?
            .get(node_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let next_node = heap
            .slice(row.as_const())?
            .get(usize::try_from(edge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nbr;

        let current = pointed_value(heap, pathfinder)?;
        if is_in_the_ilist(
            Some(heap.slice(current.seen.as_const())?),
            next_node,
            current.nseen,
        )?
        .is_some()
        {
            continue;
        }
        if forbidden_count != 0 && !forbidden.is_null() {
            let forbidden_values = heap.slice(forbidden.as_const())?;
            let mut skip = false;
            for forbidden_index in 0..forbidden_count.max(0) {
                let offset = usize::try_from(forbidden_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                    .checked_mul(2)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let first = *forbidden_values
                    .get(offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let second = *forbidden_values
                    .get(offset + 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if bIsSameBond(node, next_node, first, second) != 0 {
                    skip = true;
                    break;
                }
            }
            if skip {
                continue;
            }
        }
        heap.slice_mut(pathfinder)?[0].start = next_node;
        subgraf_pathfinder_collect_all(
            heap,
            pathfinder,
            forbidden_count,
            forbidden,
            atom_numbers,
        )?;
    }

    Ok(pointed_value(heap, pathfinder)?.nseen)
}

#[allow(non_snake_case)]
pub(crate) fn CompAtomData_GetNumMapping(
    heap: &mut SourceHeap,
    atom_data: &COMP_ATOM_DATA,
    original_numbers: SourceMutPointer<i32>,
    current_numbers: SourceMutPointer<i32>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4986 CompAtomData_GetNumMapping
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: CompAtomData_GetNumMapping
    // INCHI✔️✔️: void CompAtomData_GetNumMapping( COMP_ATOM_DATA *adata, int *orig_num, int *curr_num )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     if (!orig_num || !curr_num)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i = 0; i < adata->num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int orig = adata->at[i].orig_at_number;
    // INCHI✔️✔️:         orig_num[i] = orig;    /* orig's are from 1 */
    // INCHI✔️✔️:         curr_num[orig] = i;    /* curr's are from 0 */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: CompAtomData_GetNumMapping
    // END INCHI C FUNCTION: CompAtomData_GetNumMapping
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CompAtomData_GetNumMapping
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; orig_at_number is AT_NUMB/unsigned short and int is 32-bit.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CompAtomData_GetNumMapping

    if original_numbers.is_null() || current_numbers.is_null() {
        return Ok(());
    }
    let count = usize::try_from(atom_data.num_at.max(0))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if count == 0 {
        return Ok(());
    }
    heap.with_two_slices_mut_and_optional_third(
        original_numbers,
        current_numbers,
        Some(atom_data.at),
        |original_numbers, current_numbers, atoms| {
            let atoms = atoms.ok_or(SourceHeapError::NullPointer)?;
            for index in 0..count {
                let original = i32::from(
                    atoms
                        .get(index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .orig_at_number,
                );
                *original_numbers
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = original;
                let original_index =
                    usize::try_from(original).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                *current_numbers
                    .get_mut(original_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = index as i32;
            }
            Ok(())
        },
    )
}

#[allow(non_snake_case)]
pub(crate) fn imat_new(
    heap: &mut SourceHeap,
    m: i32,
    n: i32,
    a: &mut SourceMutPointer<SourceMutPointer<i32>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5005 imat_new
    // INCHI✔️❌: int imat_new( int m, int n, int ***a )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     if (m == 0 || n == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (*a)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         imat_free(m, *a);
    // INCHI✔️❌:         *a = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *a = (int **) inchi_calloc( m, sizeof( int * ) );
    // INCHI✔️❌:     if (NULL == *a)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < m; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ( *a )[i] = (int *) inchi_calloc( n, sizeof( int ) );
    // INCHI✔️❌:         if (NULL == ( *a )[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: imat_new

    if m == 0 || n == 0 {
        return Ok(0);
    }
    if !a.is_null() {
        imat_free(heap, m, *a)?;
        *a = SourceMutPointer::null();
    }

    let outer = match inchi_calloc::<SourceMutPointer<i32>>(
        heap,
        u64::try_from(m).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
        1,
    ) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(1),
        Err(error) => return Err(error),
    };
    *a = outer;
    let row_count =
        u64::try_from(n).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    for row_index in 0..m {
        let row = match inchi_calloc::<i32>(heap, row_count, 4) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => return Ok(1),
            Err(error) => return Err(error),
        };
        *heap
            .slice_mut(outer.offset(i64::from(row_index))?)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = row;
    }
    Ok(0)
}

pub(crate) fn cmp_components(first: &[AT_NUMB; 3], second: &[AT_NUMB; 3]) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4247 cmp_components
    // INCHI✔️✔️: int cmp_components( const void *a1, const void *a2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int ret;
    // INCHI✔️✔️:     AT_NUMB n1;
    // INCHI✔️✔️:     AT_NUMB n2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     n1 = ( (const AT_NUMB *) a1 )[0];
    // INCHI✔️✔️:     /* number of atoms in the component -- descending order */
    // INCHI✔️✔️:     n2 = ( (const AT_NUMB *) a2 )[0];
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if ((ret = (int) n2 - (int) n1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ret;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* stable sort */
    // INCHI✔️✔️:     n1 = ( (const AT_NUMB *) a1 )[1];
    // INCHI✔️✔️:     /* component ordering number -- ascending order */
    // INCHI✔️✔️:     n2 = ( (const AT_NUMB *) a2 )[1];
    // INCHI✔️✔️:
    // INCHI✔️✔️:     ret = (int) n1 - (int) n2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: cmp_components
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: cmp_components
    // INCHI✔️✔️: typedef unsigned short AT_NUMB;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: cmp_components

    let size_order = i32::from(second[0]) - i32::from(first[0]);
    if size_order != 0 {
        size_order
    } else {
        i32::from(first[1]) - i32::from(second[1])
    }
}

#[allow(non_snake_case)]
pub(crate) fn MarkDisconnectedComponents(
    heap: &mut SourceHeap,
    original: &mut ORIG_ATOM_DATA,
    mut process_old_component_numbers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4277 MarkDisconnectedComponents
    // INCHI✔️❌: int MarkDisconnectedComponents( ORIG_ATOM_DATA *orig_at_data,
    // INCHI✔️❌:                                 int bProcessOldCompNumbers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     typedef AT_NUMB AT_TRIPLE[3];
    // INCHI✔️❌:
    // INCHI✔️❌:     inp_ATOM    *at = orig_at_data->at;
    // INCHI✔️❌:     int         num_at = orig_at_data->num_inp_atoms;
    // INCHI✔️❌:     AT_NUMB     *nCurAtLen = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_NUMB *nNewCompNumber = NULL;
    // INCHI✔️❌:     AT_NUMB *nPrevAtom = NULL;
    // INCHI✔️❌:     S_CHAR  *iNeigh = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_NUMB *nOldCompNumber = NULL;
    // INCHI✔️❌:     int i, j, num_components, ret;
    // INCHI✔️❌:     int new_comp_no;
    // INCHI✔️❌:     AT_NUMB old_comp_no, another_comp_no, no_component;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* component_nbr[i][0] = number of atoms in the component i-1
    // INCHI✔️❌:     * component_nbr[i][1] = original component number (id-1) = i
    // INCHI✔️❌:     * after sorting:
    // INCHI✔️❌:     * component_nbr[j][2] = new number of component #(component_nbr[i][1]+1)
    // INCHI✔️❌:     */
    // INCHI✔️❌:     AT_TRIPLE *component_nbr = NULL;
    // INCHI✔️❌:     int fst_at, nxt_at, cur_at, cur_neq_fst;  /* moved from below 2024-09-01 DT */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* initialize */
    // INCHI✔️❌:     if (bProcessOldCompNumbers && !orig_at_data->nOldCompNumber)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bProcessOldCompNumbers = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     num_components = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     for ( j = 0; j < num_at; j ++ )
    // INCHI✔️❌:     {
    // INCHI✔️❌:     at[j].component = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = -1;
    // INCHI✔️❌:     if (!num_at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nNewCompNumber = (AT_NUMB*)inchi_calloc(num_at, sizeof(nNewCompNumber[0]));
    // INCHI✔️❌:     nPrevAtom = (AT_NUMB*)inchi_calloc(num_at, sizeof(nPrevAtom[0]));
    // INCHI✔️❌:     iNeigh = (S_CHAR*)inchi_calloc(num_at, sizeof(iNeigh[0]));
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!nNewCompNumber || !nPrevAtom || !iNeigh) /* nNewCompNumber: for non-recursive DFS only: */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*printf("\nnum_at = %d\n", num_at);*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Mark and count; avoid deep DFS recursion: it may make verifying software unhappy */
    // INCHI✔️❌:     /* nNewCompNumber[i] will contain new component number for atoms at[i], i=0..num_at-1 */
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0; j < num_at; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!nNewCompNumber[j])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* mark starting with at[j] */
    // INCHI✔️❌:             fst_at = 0;
    // INCHI✔️❌:             nxt_at = 0;
    // INCHI✔️❌:             cur_at = j;
    // INCHI✔️❌:             cur_neq_fst = 1;
    // INCHI✔️❌:             num_components++;
    // INCHI✔️❌:
    // INCHI✔️❌:             /* first time at at[j] */
    // INCHI✔️❌:             fst_at = cur_at;
    // INCHI✔️❌:             nNewCompNumber[fst_at] = (AT_NUMB) num_components;
    // INCHI✔️❌:
    // INCHI✔️❌:             /* find next neighbor */
    // INCHI✔️❌:             do
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (iNeigh[cur_at] < at[cur_at].valence)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int ineigh_incr = (int)iNeigh[cur_at];
    // INCHI✔️❌:                     nxt_at = at[cur_at].neighbor[ineigh_incr];
    // INCHI✔️❌:                     iNeigh[cur_at]++;
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (!nNewCompNumber[nxt_at])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* forward edge: found new atom */
    // INCHI✔️❌:                         nNewCompNumber[nxt_at] = (AT_NUMB) num_components;
    // INCHI✔️❌:                         nPrevAtom[nxt_at] = (AT_NUMB) cur_at;
    // INCHI✔️❌:                         cur_at = nxt_at;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (cur_at == fst_at)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     cur_neq_fst = 0;
    // INCHI✔️❌:                     /* break;  done */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     cur_at = nPrevAtom[cur_at]; /* retract */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             } while (cur_neq_fst);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_free( nPrevAtom );
    // INCHI✔️❌:     nPrevAtom = NULL;
    // INCHI✔️❌:     inchi_free( iNeigh );
    // INCHI✔️❌:     iNeigh = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Allocate more memory */
    // INCHI✔️❌:     i = inchi_max( num_components, orig_at_data->num_components );
    // INCHI✔️❌:
    // INCHI✔️❌:     nCurAtLen = (AT_NUMB*)inchi_calloc((long long)num_components + 1, sizeof(nCurAtLen[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     nOldCompNumber = (AT_NUMB*)inchi_calloc((long long)i + 1, sizeof(nOldCompNumber[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     component_nbr = (AT_TRIPLE*)inchi_calloc((long long)num_components + 1, sizeof(component_nbr[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!nCurAtLen || !nOldCompNumber || !component_nbr)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Count atoms per component and renumber the components */
    // INCHI✔️❌:     for (i = 0; i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         component_nbr[i][0] = 0; /* number of atoms in the component */
    // INCHI✔️❌:         component_nbr[i][1] = i; /* component ordering number */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0; j < num_at; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         component_nbr[(int) nNewCompNumber[j] - 1][0] ++; /* count atoms in each component */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Sort settings
    // INCHI✔️❌:     key: number of atoms
    // INCHI✔️❌:     order: descending
    // INCHI✔️❌:     stable sort
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     qsort( (void*) component_nbr[0], num_components, sizeof( component_nbr[0] ), cmp_components ); /* djb-rwth: fixed buffer overrun */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Invert the transposition */
    // INCHI✔️❌:     for (i = 0; i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nCurAtLen[i] = component_nbr[i][0];
    // INCHI✔️❌:         component_nbr[component_nbr[i][1]][2] = i + 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Renumber the components so that the component with the greatest number of atoms is the first */
    // INCHI✔️❌:     no_component = num_at + 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0; j < num_at; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* new component number for at[j] */
    // INCHI✔️❌:         new_comp_no = component_nbr[(int) nNewCompNumber[j] - 1][2] - 1; /* starts from 0 */
    // INCHI✔️❌:         if (bProcessOldCompNumbers)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* old component number for at[j] */
    // INCHI✔️❌:             old_comp_no = at[j].component;
    // INCHI✔️❌:             /* fill out nOldCompNumber[]; initially it contains zeroes */
    // INCHI✔️❌:             if (!old_comp_no)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nOldCompNumber[new_comp_no] = no_component; /* atom did not have component number */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (nOldCompNumber[new_comp_no] != old_comp_no)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!nOldCompNumber[new_comp_no])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nOldCompNumber[new_comp_no] = old_comp_no;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* at[j] moved from old comp #old_comp_no to old comp #nOldCompNumber[new_comp_no]
    // INCHI✔️❌:                     Both components cannot be equal to any current component */
    // INCHI✔️❌:                     another_comp_no = nOldCompNumber[new_comp_no];
    // INCHI✔️❌:                     for (i = 0; i < num_components; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nOldCompNumber[i] == old_comp_no ||
    // INCHI✔️❌:                              nOldCompNumber[i] == another_comp_no)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nOldCompNumber[i] = no_component;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* nOldCompNumber[new_comp_no] = num_at+1; */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* orig_at_data->nOldCompNumber */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Finally, set the new component number for atom j (NB: starts from 1 ) */
    // INCHI✔️❌:         at[j].component = new_comp_no + 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bProcessOldCompNumbers)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = 0; j < num_components; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nOldCompNumber[j] == no_component)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* the component has atom from another component */
    // INCHI✔️❌:                 nOldCompNumber[j] = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (nOldCompNumber[j] &&
    // INCHI✔️❌:                       !orig_at_data->nOldCompNumber[nOldCompNumber[j] - 1])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* the component has changed in the previous processing  */
    // INCHI✔️❌:                 nOldCompNumber[j] = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = 0; j < num_components; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nOldCompNumber[j] = j + 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = num_components;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNewCompNumber)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nNewCompNumber );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (component_nbr)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( component_nbr );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ret < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (nPrevAtom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( nPrevAtom );
    // INCHI✔️❌:             nPrevAtom = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (iNeigh)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( iNeigh );
    // INCHI✔️❌:             iNeigh = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nCurAtLen)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( nCurAtLen );
    // INCHI✔️❌:             nCurAtLen = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nOldCompNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( nOldCompNumber );
    // INCHI✔️❌:             nOldCompNumber = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         num_components = ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* avoid memory leaks */
    // INCHI✔️❌:     if (orig_at_data->nCurAtLen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( orig_at_data->nCurAtLen );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (orig_at_data->nOldCompNumber)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( orig_at_data->nOldCompNumber );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_at_data->nCurAtLen = nCurAtLen;
    // INCHI✔️❌:     orig_at_data->nOldCompNumber = nOldCompNumber;
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_at_data->num_components = num_components;
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;  /* number of disconnected components;
    // INCHI✔️❌:               1=>single connected structure        */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkDisconnectedComponents

    if process_old_component_numbers != 0 && original.nOldCompNumber.is_null() {
        process_old_component_numbers = 0;
    }
    let num_atoms = original.num_inp_atoms;
    if num_atoms == 0 {
        return Ok(0);
    }
    let allocation_count =
        u64::try_from(num_atoms).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let mut current_lengths = SourceMutPointer::<u16>::null();
    let new_component_numbers =
        inchi_calloc::<u16>(heap, allocation_count, 2).unwrap_or_else(|_| SourceMutPointer::null());
    let mut previous_atoms =
        inchi_calloc::<u16>(heap, allocation_count, 2).unwrap_or_else(|_| SourceMutPointer::null());
    let mut neighbor_indices =
        inchi_calloc::<i8>(heap, allocation_count, 1).unwrap_or_else(|_| SourceMutPointer::null());
    let mut old_component_numbers = SourceMutPointer::<u16>::null();
    let mut component_rows = SourceMutPointer::<[u16; 3]>::null();
    let mut num_components = 0_i32;
    let mut ret = -1_i32;

    let processing = (|| -> Result<bool, SourceHeapError> {
        if new_component_numbers.is_null() || previous_atoms.is_null() || neighbor_indices.is_null()
        {
            return Ok(false);
        }
        let atom_count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        for start in 0..atom_count {
            if heap.slice(new_component_numbers.as_const())?[start] != 0 {
                continue;
            }
            num_components = num_components.wrapping_add(1);
            let first_atom = start;
            let mut current_atom = start;
            heap.slice_mut(new_component_numbers)?[first_atom] = num_components as u16;
            loop {
                let atom = heap
                    .slice(original.at.as_const())?
                    .get(current_atom)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let neighbor_index = heap.slice(neighbor_indices.as_const())?[current_atom];
                if neighbor_index < atom.valence {
                    let neighbor_index = usize::try_from(i32::from(neighbor_index))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let next_atom = usize::from(
                        *atom
                            .neighbor
                            .get(neighbor_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    heap.slice_mut(neighbor_indices)?[current_atom] =
                        heap.slice(neighbor_indices.as_const())?[current_atom].wrapping_add(1);
                    if *heap
                        .slice(new_component_numbers.as_const())?
                        .get(next_atom)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        == 0
                    {
                        *heap
                            .slice_mut(new_component_numbers)?
                            .get_mut(next_atom)
                            .ok_or(SourceHeapError::PointerOutOfBounds)? = num_components as u16;
                        *heap
                            .slice_mut(previous_atoms)?
                            .get_mut(next_atom)
                            .ok_or(SourceHeapError::PointerOutOfBounds)? = current_atom as u16;
                        current_atom = next_atom;
                    }
                } else if current_atom == first_atom {
                    break;
                } else {
                    current_atom =
                        usize::from(heap.slice(previous_atoms.as_const())?[current_atom]);
                }
            }
        }

        inchi_free(heap, previous_atoms)?;
        previous_atoms = SourceMutPointer::null();
        inchi_free(heap, neighbor_indices)?;
        neighbor_indices = SourceMutPointer::null();

        let old_capacity = num_components.max(original.num_components);
        current_lengths = inchi_calloc::<u16>(
            heap,
            u64::try_from(num_components.wrapping_add(1))
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
            2,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
        old_component_numbers = inchi_calloc::<u16>(
            heap,
            u64::try_from(old_capacity.wrapping_add(1))
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
            2,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
        component_rows = inchi_calloc::<[u16; 3]>(
            heap,
            u64::try_from(num_components.wrapping_add(1))
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
            6,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
        if current_lengths.is_null() || old_component_numbers.is_null() || component_rows.is_null()
        {
            return Ok(false);
        }

        let component_count =
            usize::try_from(num_components).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        for index in 0..component_count {
            let row = &mut heap.slice_mut(component_rows)?[index];
            row[0] = 0;
            row[1] = index as u16;
        }
        for atom_index in 0..atom_count {
            let component = usize::from(heap.slice(new_component_numbers.as_const())?[atom_index])
                .checked_sub(1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(component_rows)?[component][0] =
                heap.slice(component_rows.as_const())?[component][0].wrapping_add(1);
        }
        heap.slice_mut(component_rows)?[..component_count]
            .sort_by(|first, second| cmp_components(first, second).cmp(&0));
        for index in 0..component_count {
            let row = heap.slice(component_rows.as_const())?[index];
            heap.slice_mut(current_lengths)?[index] = row[0];
            heap.slice_mut(component_rows)?[usize::from(row[1])][2] = index.wrapping_add(1) as u16;
        }

        let no_component = num_atoms.wrapping_add(1) as u16;
        for atom_index in 0..atom_count {
            let old_index = usize::from(heap.slice(new_component_numbers.as_const())?[atom_index])
                .checked_sub(1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let new_component =
                heap.slice(component_rows.as_const())?[old_index][2].wrapping_sub(1);
            if process_old_component_numbers != 0 {
                let old_component = heap.slice(original.at.as_const())?[atom_index].component;
                let new_index = usize::from(new_component);
                let mapped = heap.slice(old_component_numbers.as_const())?[new_index];
                if old_component == 0 {
                    heap.slice_mut(old_component_numbers)?[new_index] = no_component;
                } else if mapped != old_component {
                    if mapped == 0 {
                        heap.slice_mut(old_component_numbers)?[new_index] = old_component;
                    } else {
                        for component_index in 0..component_count {
                            let value =
                                heap.slice(old_component_numbers.as_const())?[component_index];
                            if value == old_component || value == mapped {
                                heap.slice_mut(old_component_numbers)?[component_index] =
                                    no_component;
                            }
                        }
                    }
                }
            }
            heap.slice_mut(original.at)?[atom_index].component = new_component.wrapping_add(1);
        }
        if process_old_component_numbers != 0 {
            for component_index in 0..component_count {
                let value = heap.slice(old_component_numbers.as_const())?[component_index];
                if value == no_component {
                    heap.slice_mut(old_component_numbers)?[component_index] = 0;
                } else if value != 0 {
                    let old_index = usize::from(value)
                        .checked_sub(1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if *heap
                        .slice(original.nOldCompNumber.as_const())?
                        .get(old_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        == 0
                    {
                        heap.slice_mut(old_component_numbers)?[component_index] = 0;
                    }
                }
            }
        } else {
            for component_index in 0..component_count {
                heap.slice_mut(old_component_numbers)?[component_index] =
                    component_index.wrapping_add(1) as u16;
            }
        }
        ret = num_components;
        Ok(true)
    })();

    let processing_error = match processing {
        Ok(_) => None,
        Err(error) => Some(error),
    };
    if !new_component_numbers.is_null() {
        inchi_free(heap, new_component_numbers)?;
    }
    if !component_rows.is_null() {
        inchi_free(heap, component_rows)?;
    }
    if ret < 0 || processing_error.is_some() {
        if !previous_atoms.is_null() {
            inchi_free(heap, previous_atoms)?;
        }
        if !neighbor_indices.is_null() {
            inchi_free(heap, neighbor_indices)?;
        }
        if !current_lengths.is_null() {
            inchi_free(heap, current_lengths)?;
            current_lengths = SourceMutPointer::null();
        }
        if !old_component_numbers.is_null() {
            inchi_free(heap, old_component_numbers)?;
            old_component_numbers = SourceMutPointer::null();
        }
        num_components = ret;
    }
    if !original.nCurAtLen.is_null() {
        inchi_free(heap, original.nCurAtLen)?;
    }
    if !original.nOldCompNumber.is_null() {
        inchi_free(heap, original.nOldCompNumber)?;
    }
    original.nCurAtLen = current_lengths;
    original.nOldCompNumber = old_component_numbers;
    original.num_components = num_components;
    if let Some(error) = processing_error {
        return Err(error);
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn ExtractConnectedComponent(
    heap: &mut SourceHeap,
    atoms: &[inp_ATOM],
    num_atoms: i32,
    component_number: i32,
    component_atoms: &mut [inp_ATOM],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4558 ExtractConnectedComponent
    // INCHI✔️❌: int ExtractConnectedComponent( inp_ATOM *at,
    // INCHI✔️❌:                                int num_at,
    // INCHI✔️❌:                                int component_number,
    // INCHI✔️❌:                                inp_ATOM *component_at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, num_component_at;
    // INCHI✔️❌:     AT_NUMB *number;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL == ( number = (AT_NUMB*) inchi_calloc( num_at, sizeof( AT_NUMB ) ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return CT_OUT_OF_RAM; /* out of memory */  /*   <BRKPT> */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* copy atoms */
    // INCHI✔️❌:     for (i = 0, num_component_at = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[i].component == component_number)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             number[i] = num_component_at;
    // INCHI✔️❌:             component_at[num_component_at++] = at[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* renumber neighbors */
    // INCHI✔️❌:     for (i = 0; i < num_component_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         component_at[i].orig_compt_at_numb = (AT_NUMB) ( i + 1 );
    // INCHI✔️❌:         for (j = 0; j < component_at[i].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             component_at[i].neighbor[j] = number[(int) component_at[i].neighbor[j]];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_free( number );
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_component_at;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ExtractConnectedComponent

    let Ok(allocation_count) = u64::try_from(num_atoms) else {
        return Ok(CT_OUT_OF_RAM);
    };
    let number = match inchi_calloc::<u16>(heap, allocation_count, 2) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(CT_OUT_OF_RAM),
        Err(error) => return Err(error),
    };
    let processing = (|| -> Result<i32, SourceHeapError> {
        let atom_count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut component_count = 0_usize;
        for (atom_index, atom) in atoms
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .iter()
            .enumerate()
        {
            if i32::from(atom.component) == component_number {
                heap.slice_mut(number)?[atom_index] = component_count as u16;
                *component_atoms
                    .get_mut(component_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = atom.clone();
                component_count += 1;
            }
        }
        for (component_index, atom) in component_atoms
            .get_mut(..component_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .iter_mut()
            .enumerate()
        {
            atom.orig_compt_at_numb = component_index.wrapping_add(1) as u16;
            let valence = usize::try_from(i32::from(atom.valence).max(0))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            for neighbor in atom
                .neighbor
                .get_mut(..valence)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                *neighbor = *heap
                    .slice(number.as_const())?
                    .get(usize::from(*neighbor))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
            }
        }
        i32::try_from(component_count).map_err(|_| SourceHeapError::SourceIntegerOverflow)
    })();
    inchi_free(heap, number)?;
    processing
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{INChI_IsotopicAtom, INChI_IsotopicTGroup, subgraf_edge};

    #[test]
    fn source_port__strutil__cmp_iso_atw_diff_component_no__line_166() {
        let atom = |isotope, component| inp_ATOM {
            iso_atw_diff: isotope,
            component,
            ..inp_ATOM::default()
        };

        let minimum = atom(i8::MIN, u16::MAX);
        let maximum = atom(i8::MAX, 0);
        let minimum_before = minimum.clone();
        let maximum_before = maximum.clone();
        assert_eq!(cmp_iso_atw_diff_component_no(&minimum, &maximum), -255);
        assert_eq!(cmp_iso_atw_diff_component_no(&maximum, &minimum), 255);
        assert_eq!(minimum, minimum_before);
        assert_eq!(maximum, maximum_before);

        assert_eq!(
            cmp_iso_atw_diff_component_no(&atom(7, 0), &atom(7, u16::MAX)),
            -65_535
        );
        assert_eq!(
            cmp_iso_atw_diff_component_no(&atom(7, u16::MAX), &atom(7, 0)),
            65_535
        );
        assert_eq!(cmp_iso_atw_diff_component_no(&atom(7, 23), &atom(7, 23)), 0);

        let mut atoms = [
            atom(1, 9),
            atom(-1, u16::MAX),
            atom(1, 2),
            atom(-1, 0),
            atom(1, 2),
        ];
        atoms.sort_by(|first, second| cmp_iso_atw_diff_component_no(first, second).cmp(&0));
        assert_eq!(
            atoms.map(|entry| (entry.iso_atw_diff, entry.component)),
            [(-1, 0), (-1, u16::MAX), (1, 2), (1, 2), (1, 9)]
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__cmp_iso_atw_diff_component_no__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--cmp-iso-atw-diff-component-no-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "cmp_iso_atw_diff_component_no");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_isotope = |field: &str| {
                i8::try_from(
                    official["input"][field]
                        .as_i64()
                        .unwrap_or_else(|| panic!("{case_id}: {field} must be i8")),
                )
                .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds i8"))
            };
            let parse_component = |field: &str| {
                u16::try_from(
                    official["input"][field]
                        .as_u64()
                        .unwrap_or_else(|| panic!("{case_id}: {field} must be u16")),
                )
                .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds u16"))
            };
            let first = inp_ATOM {
                iso_atw_diff: parse_isotope("first_iso_atw_diff"),
                component: parse_component("first_component"),
                ..inp_ATOM::default()
            };
            let second = inp_ATOM {
                iso_atw_diff: parse_isotope("second_iso_atw_diff"),
                component: parse_component("second_component"),
                ..inp_ATOM::default()
            };
            let first_before = first.clone();
            let second_before = second.clone();
            let rust = cmp_iso_atw_diff_component_no(&first, &second);
            let expected = i32::try_from(
                official["output"]["result"]
                    .as_i64()
                    .expect("result must be an integer"),
            )
            .expect("official result must fit i32");
            assert_eq!(rust, expected, "{case_id}");
            assert_eq!(first, first_before, "{case_id}: Rust first input mutation");
            assert_eq!(
                second, second_before,
                "{case_id}: Rust second input mutation"
            );
            assert_eq!(official["output"]["first_unchanged"], true, "{case_id}");
            assert_eq!(official["output"]["second_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 7);
    }

    #[test]
    fn source_port__strutil__cmp_components__line_4247() {
        assert_eq!(cmp_components(&[9, 4, 0], &[3, 0, 0]), -6);
        assert_eq!(cmp_components(&[3, 0, 0], &[9, 4, 0]), 6);
        assert_eq!(cmp_components(&[u16::MAX, 0, 0], &[0, 0, 0]), -65_535);
        assert_eq!(cmp_components(&[0, 0, 0], &[u16::MAX, 0, 0]), 65_535);

        assert_eq!(cmp_components(&[7, 2, 0], &[7, 8, 0]), -6);
        assert_eq!(cmp_components(&[7, 8, 0], &[7, 2, 0]), 6);
        assert_eq!(cmp_components(&[7, 0, 0], &[7, u16::MAX, 0]), -65_535);
        assert_eq!(cmp_components(&[7, u16::MAX, 0], &[7, 0, 0]), 65_535);

        assert_eq!(cmp_components(&[7, 2, 0], &[7, 2, u16::MAX]), 0);

        let mut rows = [[2, 1, 7], [5, 3, 8], [5, 0, 9], [2, 0, 10]];
        rows.sort_by(|first, second| cmp_components(first, second).cmp(&0));
        assert_eq!(rows, [[5, 0, 9], [5, 3, 8], [2, 0, 10], [2, 1, 7]]);
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__cmp_components__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--cmp-components-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "cmp_components");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_row = |field: &str| {
                let values = official["input"][field]
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be an array"));
                assert_eq!(values.len(), 3, "{case_id}: {field} length");
                let mut row = [0_u16; 3];
                for (index, value) in values.iter().enumerate() {
                    row[index] = u16::try_from(
                        value
                            .as_u64()
                            .unwrap_or_else(|| panic!("{case_id}: {field}[{index}] must be u16")),
                    )
                    .unwrap_or_else(|_| panic!("{case_id}: {field}[{index}] exceeds u16"));
                }
                row
            };
            let first = parse_row("first");
            let second = parse_row("second");
            let first_before = first;
            let second_before = second;
            let rust = cmp_components(&first, &second);
            let expected = i32::try_from(
                official["output"]["result"]
                    .as_i64()
                    .expect("result must be an integer"),
            )
            .expect("official result must fit i32");
            assert_eq!(rust, expected, "{case_id}");
            assert_eq!(first, first_before, "{case_id}: first input mutation");
            assert_eq!(second, second_before, "{case_id}: second input mutation");
            record_count += 1;
        }
        assert_eq!(record_count, 10);
    }

    fn terminal_hdt_atom(symbol: u8, neighbor: u16) -> inp_ATOM {
        let mut atom = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        atom.elname[0] = symbol as i8;
        atom.neighbor[0] = neighbor;
        atom.bond_type[0] = BOND_TYPE_SINGLE as u8;
        atom
    }

    #[test]
    fn source_port__strutil__nfindoneom__line_932() {
        let mut no_candidates = [];
        assert_eq!(nFindOneOM(&[], i32::MIN, &mut no_candidates, 0), Ok(-1));
        assert_eq!(
            nFindOneOM(&[], i32::MAX, &mut no_candidates, i32::MIN),
            Ok(-1)
        );

        let mut one_candidate = [i32::MIN];
        assert_eq!(
            nFindOneOM(&[], i32::MAX, &mut one_candidate, 1),
            Ok(i32::MIN)
        );

        let mut atoms = vec![inp_ATOM::default(); 5];
        atoms[0].neighbor[..4].copy_from_slice(&[1, 2, 3, 4]);
        atoms[1].valence = 4;
        atoms[2].valence = 2;
        atoms[3].valence = 1;
        atoms[4].valence = 3;
        let mut unique_valence = [0, 1, 2, 3];
        assert_eq!(nFindOneOM(&atoms, 0, &mut unique_valence, 4), Ok(2));
        assert_eq!(unique_valence, [2, 1, 2, 3]);

        atoms[1].valence = 1;
        atoms[2].valence = 1;
        atoms[3].valence = 2;
        atoms[4].valence = 1;
        atoms[1].el_number = 16;
        atoms[2].el_number = 8;
        atoms[4].el_number = 15;
        let mut unique_element = [0, 1, 2, 3];
        assert_eq!(nFindOneOM(&atoms, 0, &mut unique_element, 4), Ok(1));
        assert_eq!(unique_element, [1, 1, 3, 3]);

        atoms[1].valence = 2;
        atoms[2].valence = 2;
        atoms[1].el_number = 8;
        atoms[2].el_number = 8;
        let mut nonterminal = [0, 1];
        assert_eq!(nFindOneOM(&atoms, 0, &mut nonterminal, 2), Ok(-1));
        assert_eq!(nonterminal, [0, 1]);

        atoms[1].valence = 1;
        atoms[2].valence = 1;
        atoms[4].valence = 1;
        atoms[1].el_number = 8;
        atoms[2].el_number = 8;
        atoms[4].el_number = 8;
        atoms[1].iso_atw_diff = 127;
        atoms[2].iso_atw_diff = -128;
        atoms[4].iso_atw_diff = 0;
        let mut source_field_asymmetry = [0, 1, 3];
        assert_eq!(nFindOneOM(&atoms, 0, &mut source_field_asymmetry, 3), Ok(1));
        assert_eq!(source_field_asymmetry, [1, 3, 3]);

        atoms[1].el_number = 0;
        atoms[2].el_number = 0;
        atoms[1].iso_atw_diff = 1;
        let mut zero_element = [0, 1];
        assert_eq!(nFindOneOM(&atoms, 0, &mut zero_element, 2), Ok(1));
        assert_eq!(zero_element, [1, 1]);

        let mut too_short = [0];
        assert_eq!(
            nFindOneOM(&atoms, 0, &mut too_short, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut invalid_center = [0, 1];
        assert_eq!(
            nFindOneOM(&atoms, -1, &mut invalid_center, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    fn ion_atom(element: u8, charge: i8) -> inp_ATOM {
        inp_ATOM {
            el_number: element,
            charge,
            ..inp_ATOM::default()
        }
    }

    fn ion_bond(atoms: &mut [inp_ATOM], first: usize, second: usize, bond_type: u8) {
        let first_position = atoms[first].valence as usize;
        let second_position = atoms[second].valence as usize;
        atoms[first].neighbor[first_position] = second as u16;
        atoms[first].bond_type[first_position] = bond_type;
        atoms[first].valence += 1;
        atoms[first].chem_bonds_valence += bond_type as i8;
        atoms[second].neighbor[second_position] = first as u16;
        atoms[second].bond_type[second_position] = bond_type;
        atoms[second].valence += 1;
        atoms[second].chem_bonds_valence += bond_type as i8;
    }

    fn run_one_ion_pair(mut atoms: Vec<inp_ATOM>) -> Vec<inp_ATOM> {
        let count = atoms.len() as i32;
        assert_eq!(remove_ion_pairs(count, &mut atoms), Ok(1));
        assert!(atoms.iter().all(|atom| atom.cFlags == 0));
        atoms
    }

    #[test]
    fn source_port__strutil__remove_ion_pairs__line_1049() {
        assert_eq!(remove_ion_pairs(i32::MIN, &mut []), Ok(0));
        let mut no_candidate = vec![ion_atom(6, 0)];
        assert_eq!(remove_ion_pairs(1, &mut no_candidate), Ok(0));
        assert_eq!(
            remove_ion_pairs(2, &mut no_candidate),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut type_1 = vec![
            ion_atom(7, 1),
            ion_atom(8, 0),
            ion_atom(8, -1),
            ion_atom(6, 0),
        ];
        type_1[0].radical = 2;
        type_1[2].radical = 3;
        ion_bond(&mut type_1, 0, 1, 2);
        ion_bond(&mut type_1, 0, 2, 1);
        ion_bond(&mut type_1, 0, 3, 1);
        let type_1 = run_one_ion_pair(type_1);
        assert_eq!((type_1[0].charge, type_1[2].charge), (0, 0));
        assert_eq!((type_1[0].bond_type[1], type_1[2].bond_type[0]), (2, 2));
        assert_eq!((type_1[0].radical, type_1[2].radical), (0, 0));

        let mut type_1a = vec![
            ion_atom(15, 1),
            ion_atom(8, -1),
            ion_atom(6, 0),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        for neighbor in 1..5 {
            ion_bond(&mut type_1a, 0, neighbor, 1);
        }
        let type_1a = run_one_ion_pair(type_1a);
        assert_eq!((type_1a[0].charge, type_1a[1].charge), (0, 0));
        assert_eq!((type_1a[0].bond_type[0], type_1a[1].bond_type[0]), (2, 2));

        let mut type_2 = vec![
            ion_atom(8, 0),
            ion_atom(7, 0),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        ion_bond(&mut type_2, 0, 1, 2);
        ion_bond(&mut type_2, 1, 2, 1);
        ion_bond(&mut type_2, 2, 3, 1);
        type_2[2].radical = 2;
        let type_2 = run_one_ion_pair(type_2);
        assert_eq!((type_2[1].bond_type[1], type_2[2].bond_type[0]), (3, 3));
        assert_eq!(
            (type_2[1].chem_bonds_valence, type_2[2].chem_bonds_valence),
            (5, 4)
        );
        assert_eq!(type_2[2].radical, 0);

        let mut type_3 = vec![ion_atom(8, 0), ion_atom(8, 1), ion_atom(6, -1)];
        ion_bond(&mut type_3, 0, 1, 2);
        ion_bond(&mut type_3, 1, 2, 1);
        type_3[2].num_H = 2;
        let type_3 = run_one_ion_pair(type_3);
        assert_eq!((type_3[1].charge, type_3[2].charge), (0, 0));
        assert_eq!((type_3[1].bond_type[1], type_3[2].bond_type[0]), (2, 2));

        let mut type_4 = vec![
            ion_atom(8, -1),
            ion_atom(7, 1),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        ion_bond(&mut type_4, 0, 1, 1);
        ion_bond(&mut type_4, 1, 2, 2);
        ion_bond(&mut type_4, 1, 3, 1);
        let type_4 = run_one_ion_pair(type_4);
        assert_eq!((type_4[0].charge, type_4[1].charge), (0, 0));
        assert_eq!((type_4[0].bond_type[0], type_4[1].bond_type[0]), (2, 2));

        let mut type_5 = vec![ion_atom(8, -1), ion_atom(8, 1), ion_atom(6, 0)];
        ion_bond(&mut type_5, 0, 1, 1);
        ion_bond(&mut type_5, 1, 2, 2);
        let type_5 = run_one_ion_pair(type_5);
        assert_eq!((type_5[0].charge, type_5[1].charge), (0, 0));
        assert_eq!((type_5[0].bond_type[0], type_5[1].bond_type[0]), (2, 2));

        let mut type_6 = vec![ion_atom(8, -1), ion_atom(8, 0), ion_atom(6, 1)];
        ion_bond(&mut type_6, 0, 1, 1);
        ion_bond(&mut type_6, 1, 2, 1);
        type_6[2].num_H = 2;
        let type_6 = run_one_ion_pair(type_6);
        assert_eq!((type_6[0].charge, type_6[2].charge), (0, 0));
        assert_eq!((type_6[0].bond_type[0], type_6[1].bond_type[1]), (2, 2));
        assert_eq!(type_6[2].bond_type[0], 2);

        let mut type_7 = vec![
            ion_atom(7, -1),
            ion_atom(7, 1),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        ion_bond(&mut type_7, 0, 1, 2);
        ion_bond(&mut type_7, 1, 2, 1);
        ion_bond(&mut type_7, 1, 3, 1);
        let type_7 = run_one_ion_pair(type_7);
        assert_eq!((type_7[0].charge, type_7[1].charge), (0, 0));
        assert_eq!((type_7[0].bond_type[0], type_7[1].bond_type[0]), (3, 3));

        let mut type_8 = vec![ion_atom(7, -1), ion_atom(8, 1), ion_atom(6, 0)];
        ion_bond(&mut type_8, 0, 1, 2);
        ion_bond(&mut type_8, 1, 2, 1);
        let type_8 = run_one_ion_pair(type_8);
        assert_eq!((type_8[0].charge, type_8[1].charge), (0, 0));
        assert_eq!((type_8[0].bond_type[0], type_8[1].bond_type[0]), (3, 3));

        let mut type_9 = vec![ion_atom(7, -1), ion_atom(6, 1), ion_atom(6, 0)];
        ion_bond(&mut type_9, 0, 1, 2);
        ion_bond(&mut type_9, 1, 2, 1);
        let type_9 = run_one_ion_pair(type_9);
        assert_eq!((type_9[0].charge, type_9[1].charge), (0, 0));
        assert_eq!((type_9[0].bond_type[0], type_9[1].bond_type[0]), (3, 3));

        let mut type_10 = vec![
            ion_atom(7, 1),
            ion_atom(6, -1),
            ion_atom(6, 0),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        for neighbor in 1..5 {
            ion_bond(&mut type_10, 0, neighbor, 1);
        }
        type_10[1].num_H = 2;
        let type_10 = run_one_ion_pair(type_10);
        assert_eq!((type_10[0].charge, type_10[1].charge), (0, 0));
        assert_eq!((type_10[0].bond_type[0], type_10[1].bond_type[0]), (2, 2));

        let mut type_11 = vec![
            ion_atom(7, 1),
            ion_atom(6, -1),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        ion_bond(&mut type_11, 0, 1, 2);
        ion_bond(&mut type_11, 0, 2, 1);
        ion_bond(&mut type_11, 0, 3, 1);
        type_11[1].num_H = 1;
        let type_11 = run_one_ion_pair(type_11);
        assert_eq!((type_11[0].charge, type_11[1].charge), (0, 0));
        assert_eq!((type_11[0].bond_type[0], type_11[1].bond_type[0]), (3, 3));

        let mut type_12 = vec![
            ion_atom(7, 1),
            ion_atom(7, -1),
            ion_atom(6, 0),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        for neighbor in 1..5 {
            ion_bond(&mut type_12, 0, neighbor, 1);
        }
        type_12[1].num_H = 1;
        let type_12 = run_one_ion_pair(type_12);
        assert_eq!((type_12[0].charge, type_12[1].charge), (0, 0));
        assert_eq!((type_12[0].bond_type[0], type_12[1].bond_type[0]), (2, 2));

        let mut type_13 = vec![ion_atom(8, 1), ion_atom(6, -1), ion_atom(6, 0)];
        ion_bond(&mut type_13, 0, 1, 1);
        ion_bond(&mut type_13, 0, 2, 2);
        type_13[1].num_H = 2;
        let type_13 = run_one_ion_pair(type_13);
        assert_eq!((type_13[0].charge, type_13[1].charge), (0, 0));
        assert_eq!((type_13[0].bond_type[0], type_13[1].bond_type[0]), (2, 2));

        let mut type_14 = vec![ion_atom(8, 1), ion_atom(6, -1), ion_atom(6, 0)];
        ion_bond(&mut type_14, 0, 1, 2);
        ion_bond(&mut type_14, 0, 2, 1);
        type_14[1].num_H = 1;
        let type_14 = run_one_ion_pair(type_14);
        assert_eq!((type_14[0].charge, type_14[1].charge), (0, 0));
        assert_eq!((type_14[0].bond_type[0], type_14[1].bond_type[0]), (3, 3));

        let mut type_15 = vec![ion_atom(8, 1), ion_atom(7, -1), ion_atom(6, 0)];
        ion_bond(&mut type_15, 0, 1, 1);
        ion_bond(&mut type_15, 0, 2, 2);
        type_15[1].num_H = 1;
        let type_15 = run_one_ion_pair(type_15);
        assert_eq!((type_15[0].charge, type_15[1].charge), (0, 0));
        assert_eq!((type_15[0].bond_type[0], type_15[1].bond_type[0]), (2, 2));

        let mut type_16 = vec![ion_atom(8, 0), ion_atom(6, 1), ion_atom(7, -1)];
        ion_bond(&mut type_16, 0, 1, 1);
        ion_bond(&mut type_16, 0, 2, 1);
        type_16[1].num_H = 2;
        type_16[2].num_H = 1;
        let type_16 = run_one_ion_pair(type_16);
        assert_eq!((type_16[1].charge, type_16[2].charge), (0, 0));
        assert_eq!((type_16[0].bond_type[0], type_16[1].bond_type[0]), (2, 2));
        assert_eq!((type_16[0].bond_type[1], type_16[2].bond_type[0]), (2, 2));

        let mut type_17 = vec![
            ion_atom(7, 0),
            ion_atom(6, 1),
            ion_atom(6, -1),
            ion_atom(6, 0),
        ];
        for neighbor in 1..4 {
            ion_bond(&mut type_17, 0, neighbor, 1);
        }
        type_17[1].num_H = 2;
        type_17[2].num_H = 2;
        let type_17 = run_one_ion_pair(type_17);
        assert_eq!((type_17[1].charge, type_17[2].charge), (0, 0));
        assert_eq!((type_17[0].bond_type[0], type_17[1].bond_type[0]), (2, 2));
        assert_eq!((type_17[0].bond_type[1], type_17[2].bond_type[0]), (2, 2));

        let mut type_18_charged = vec![ion_atom(7, 0), ion_atom(6, -1), ion_atom(6, 1)];
        ion_bond(&mut type_18_charged, 0, 1, 1);
        ion_bond(&mut type_18_charged, 0, 2, 2);
        type_18_charged[1].num_H = 2;
        type_18_charged[2].num_H = 1;
        type_18_charged[0].radical = 2;
        type_18_charged[1].radical = 2;
        type_18_charged[2].radical = 3;
        let type_18_charged = run_one_ion_pair(type_18_charged);
        assert_eq!(
            (type_18_charged[1].charge, type_18_charged[2].charge),
            (0, 0)
        );
        assert_eq!(
            (
                type_18_charged[0].bond_type[0],
                type_18_charged[1].bond_type[0]
            ),
            (2, 2)
        );
        assert_eq!(
            (
                type_18_charged[0].bond_type[1],
                type_18_charged[2].bond_type[0]
            ),
            (3, 3)
        );
        assert_eq!(
            (
                type_18_charged[0].radical,
                type_18_charged[1].radical,
                type_18_charged[2].radical
            ),
            (0, 0, 0)
        );

        let mut type_18_neutral = vec![ion_atom(7, 0), ion_atom(6, 0), ion_atom(6, 0)];
        ion_bond(&mut type_18_neutral, 0, 1, 2);
        ion_bond(&mut type_18_neutral, 0, 2, 1);
        type_18_neutral[1].num_H = 2;
        type_18_neutral[2].num_H = 1;
        let type_18_neutral = run_one_ion_pair(type_18_neutral);
        assert_eq!(
            (type_18_neutral[1].charge, type_18_neutral[2].charge),
            (0, 0)
        );
        assert_eq!(
            (
                type_18_neutral[0].bond_type[0],
                type_18_neutral[1].bond_type[0]
            ),
            (2, 2)
        );
        assert_eq!(
            (
                type_18_neutral[0].bond_type[1],
                type_18_neutral[2].bond_type[0]
            ),
            (3, 3)
        );
    }

    #[test]
    fn source_port__strutil__bisammoniumsalt__line_2300() {
        let mut oxygen_index = 71;
        let mut neighbor_position = 72;
        let mut explicit_h = [1_i8, 2, 3, 4];
        assert_eq!(
            bIsAmmoniumSalt(
                &[ion_atom(6, 0)],
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(0)
        );
        assert_eq!(
            (oxygen_index, neighbor_position, explicit_h),
            (71, 72, [1, 2, 3, 4])
        );
        assert_eq!(
            bIsAmmoniumSalt(
                &[ion_atom(7, 0)],
                -1,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut valence_mismatch = ion_atom(7, 0);
        valence_mismatch.num_H = 4;
        assert_eq!(
            bIsAmmoniumSalt(
                &[valence_mismatch],
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(0)
        );
        assert_eq!(
            (oxygen_index, neighbor_position, explicit_h),
            (71, 72, [1, 2, 3, 4])
        );

        let mut oxygen_salt = vec![ion_atom(7, 1), ion_atom(8, -1), ion_atom(6, 0)];
        oxygen_salt[0].num_H = 4;
        oxygen_salt[2].radical = RADICAL_SINGLET as i8;
        ion_bond(&mut oxygen_salt, 0, 1, 1);
        ion_bond(&mut oxygen_salt, 1, 2, 1);
        assert_eq!(
            bIsAmmoniumSalt(
                &oxygen_salt,
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(1)
        );
        assert_eq!(
            (oxygen_index, neighbor_position, explicit_h),
            (1, 0, [0; 4])
        );

        for halogen in [9_u8, 17, 35, 53] {
            let mut halogen_salt = vec![ion_atom(7, 0), ion_atom(halogen, 0)];
            halogen_salt[0].num_H = 4;
            ion_bond(&mut halogen_salt, 0, 1, 1);
            oxygen_index = -7;
            neighbor_position = -8;
            explicit_h = [9; 4];
            assert_eq!(
                bIsAmmoniumSalt(
                    &halogen_salt,
                    0,
                    &mut oxygen_index,
                    &mut neighbor_position,
                    &mut explicit_h,
                ),
                Ok(1),
                "halogen {halogen}"
            );
            assert_eq!(
                (oxygen_index, neighbor_position, explicit_h),
                (1, 0, [0; 4])
            );
        }

        let mut explicit_isotopes = vec![
            ion_atom(7, 0),
            ion_atom(8, 0),
            ion_atom(6, 0),
            ion_atom(1, 0),
            ion_atom(1, 0),
            ion_atom(1, 0),
        ];
        explicit_isotopes[0].num_H = 1;
        explicit_isotopes[3].iso_atw_diff = 0;
        explicit_isotopes[4].iso_atw_diff = 1;
        explicit_isotopes[5].iso_atw_diff = 3;
        ion_bond(&mut explicit_isotopes, 0, 1, 1);
        ion_bond(&mut explicit_isotopes, 1, 2, 1);
        ion_bond(&mut explicit_isotopes, 0, 3, 1);
        ion_bond(&mut explicit_isotopes, 0, 4, 1);
        ion_bond(&mut explicit_isotopes, 0, 5, 1);
        explicit_h = [9; 4];
        assert_eq!(
            bIsAmmoniumSalt(
                &explicit_isotopes,
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(1)
        );
        assert_eq!(
            (oxygen_index, neighbor_position, explicit_h),
            (1, 0, [1, 1, 0, 1])
        );

        let mut rejected_neighbor = vec![ion_atom(7, 0), ion_atom(6, 0)];
        rejected_neighbor[0].num_H = 4;
        ion_bond(&mut rejected_neighbor, 0, 1, 1);
        oxygen_index = 81;
        neighbor_position = 82;
        explicit_h = [8; 4];
        assert_eq!(
            bIsAmmoniumSalt(
                &rejected_neighbor,
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(0)
        );
        assert_eq!(
            (oxygen_index, neighbor_position, explicit_h),
            (81, 82, [0; 4])
        );

        let mut bad_charge = oxygen_salt.clone();
        bad_charge[0].charge = i8::MIN;
        bad_charge[1].charge = i8::MIN;
        oxygen_index = 91;
        neighbor_position = 92;
        assert_eq!(
            bIsAmmoniumSalt(
                &bad_charge,
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(0)
        );
        assert_eq!((oxygen_index, neighbor_position), (91, 92));

        let mut bad_carbon = oxygen_salt.clone();
        bad_carbon[2].charge = 1;
        assert_eq!(
            bIsAmmoniumSalt(
                &bad_carbon,
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(0)
        );
        bad_carbon[2].charge = 0;
        bad_carbon[2].radical = 2;
        assert_eq!(
            bIsAmmoniumSalt(
                &bad_carbon,
                0,
                &mut oxygen_index,
                &mut neighbor_position,
                &mut explicit_h,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__strutil__disconnectammoniumsalt__line_2398() {
        let mut implicit = vec![ion_atom(7, 1), ion_atom(8, -1), ion_atom(6, 0)];
        implicit[0].num_H = 4;
        ion_bond(&mut implicit, 0, 1, 1);
        ion_bond(&mut implicit, 1, 2, 1);
        assert_eq!(
            DisconnectAmmoniumSalt(&mut implicit, 0, 1, 0, &[0; 4]),
            Ok(1)
        );
        assert_eq!((implicit[0].charge, implicit[1].charge), (0, 0));
        assert_eq!((implicit[0].num_H, implicit[1].num_H), (3, 1));
        assert_eq!(
            (implicit[0].valence, implicit[0].chem_bonds_valence),
            (0, 0)
        );
        assert_eq!(
            (implicit[1].valence, implicit[1].chem_bonds_valence),
            (1, 1)
        );
        assert_eq!(implicit[1].neighbor[0], 2);

        for isotope in [1_usize, 2] {
            let mut implicit_isotope = vec![ion_atom(7, 0), ion_atom(8, 0), ion_atom(6, 0)];
            implicit_isotope[0].num_iso_H[isotope] = 1;
            ion_bond(&mut implicit_isotope, 0, 1, 1);
            ion_bond(&mut implicit_isotope, 1, 2, 1);
            assert_eq!(
                DisconnectAmmoniumSalt(&mut implicit_isotope, 0, 1, 0, &[0; 4]),
                Ok(1)
            );
            assert_eq!(implicit_isotope[0].num_iso_H[isotope], 0);
            assert_eq!(implicit_isotope[1].num_iso_H[isotope], 1);
        }

        let mut ignored_slot_zero = vec![ion_atom(7, 0), ion_atom(8, 0), ion_atom(6, 0)];
        ignored_slot_zero[0].num_iso_H[0] = 1;
        ion_bond(&mut ignored_slot_zero, 0, 1, 1);
        ion_bond(&mut ignored_slot_zero, 1, 2, 1);
        assert_eq!(
            DisconnectAmmoniumSalt(&mut ignored_slot_zero, 0, 1, 0, &[0; 4]),
            Ok(1)
        );
        assert_eq!(ignored_slot_zero[0].num_iso_H[0], 1);
        assert_eq!(ignored_slot_zero[1].num_iso_H[0], 0);

        for isotope in 0..=NUM_H_ISOTOPES as usize {
            let mut explicit = vec![
                ion_atom(7, 0),
                ion_atom(8, 0),
                ion_atom(6, 0),
                ion_atom(1, 0),
            ];
            explicit[3].iso_atw_diff = isotope as i8;
            explicit[3].bond_stereo[0] = 6;
            ion_bond(&mut explicit, 0, 1, 1);
            ion_bond(&mut explicit, 1, 2, 1);
            ion_bond(&mut explicit, 0, 3, 1);
            let mut explicit_count = [0_i8; 4];
            explicit_count[isotope] = 1;
            assert_eq!(
                DisconnectAmmoniumSalt(&mut explicit, 0, 1, 0, &explicit_count),
                Ok(1),
                "isotope {isotope}"
            );
            assert_eq!((explicit[0].valence, explicit[1].valence), (0, 2));
            assert_eq!((explicit[1].neighbor[0], explicit[1].neighbor[1]), (2, 3));
            assert_eq!(explicit[3].neighbor[0], 1);
            assert_eq!(explicit[3].bond_stereo[0], 0);
        }

        let mut nearest = vec![
            ion_atom(7, 0),
            ion_atom(8, 0),
            ion_atom(6, 0),
            ion_atom(1, 0),
            ion_atom(1, 0),
        ];
        nearest[3].x = 5.0;
        nearest[4].x = 1.0;
        ion_bond(&mut nearest, 0, 1, 1);
        ion_bond(&mut nearest, 1, 2, 1);
        ion_bond(&mut nearest, 0, 3, 1);
        ion_bond(&mut nearest, 0, 4, 1);
        assert_eq!(
            DisconnectAmmoniumSalt(&mut nearest, 0, 1, 0, &[2, 0, 0, 0]),
            Ok(1)
        );
        assert_eq!(nearest[0].neighbor[0], 3);
        assert_eq!((nearest[1].neighbor[0], nearest[1].neighbor[1]), (2, 4));
        assert_eq!(nearest[4].neighbor[0], 1);

        let mut noncancelling_charge = vec![ion_atom(7, 1), ion_atom(8, 1), ion_atom(6, 0)];
        noncancelling_charge[0].num_H = 1;
        ion_bond(&mut noncancelling_charge, 0, 1, 1);
        ion_bond(&mut noncancelling_charge, 1, 2, 1);
        assert_eq!(
            DisconnectAmmoniumSalt(&mut noncancelling_charge, 0, 1, 0, &[0; 4]),
            Ok(1)
        );
        assert_eq!(
            (
                noncancelling_charge[0].charge,
                noncancelling_charge[1].charge
            ),
            (1, 1)
        );
    }

    #[test]
    fn source_port__strutil__bismetalsalt__line_2511() {
        fn halide_salt(metal: u8, charge: i8, halogen: u8) -> Vec<inp_ATOM> {
            let mut atoms = vec![ion_atom(metal, charge), ion_atom(halogen, 0)];
            ion_bond(&mut atoms, 0, 1, BOND_TYPE_SINGLE as u8);
            atoms
        }

        fn oxygen_carbon_salt() -> Vec<inp_ATOM> {
            let mut atoms = vec![
                ion_atom(11, 0),
                ion_atom(8, 0),
                ion_atom(6, 0),
                ion_atom(6, 0),
            ];
            ion_bond(&mut atoms, 0, 1, BOND_TYPE_SINGLE as u8);
            ion_bond(&mut atoms, 1, 2, BOND_TYPE_SINGLE as u8);
            ion_bond(&mut atoms, 2, 3, BOND_TYPE_TRIPLE as u8);
            atoms
        }

        assert_eq!(
            bIsMetalSalt(&[], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            bIsMetalSalt(&[], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(bIsMetalSalt(&[ion_atom(11, 0)], 0), Ok(0));
        assert_eq!(bIsMetalSalt(&halide_salt(6, 0, 9), 0), Ok(0));

        let neutral_metal = halide_salt(11, 0, 9);
        assert_eq!(bIsMetalSalt(&neutral_metal, 0), Ok(1));
        let mut metal_with_h = neutral_metal.clone();
        metal_with_h[0].num_H = 1;
        assert_eq!(bIsMetalSalt(&metal_with_h, 0), Ok(0));
        let mut negative_metal = neutral_metal.clone();
        negative_metal[0].charge = -1;
        assert_eq!(bIsMetalSalt(&negative_metal, 0), Ok(0));

        let positive_metal = halide_salt(12, 1, 9);
        assert_eq!(bIsMetalSalt(&positive_metal, 0), Ok(1));

        let mut second_valence_metal = vec![
            ion_atom(25, 0),
            ion_atom(9, 0),
            ion_atom(17, 0),
            ion_atom(35, 0),
        ];
        for ligand in 1..second_valence_metal.len() {
            ion_bond(&mut second_valence_metal, 0, ligand, BOND_TYPE_SINGLE as u8);
        }
        assert_eq!(bIsMetalSalt(&second_valence_metal, 0), Ok(1));

        for halogen in [9_u8, 17, 35, 53] {
            assert_eq!(
                bIsMetalSalt(&halide_salt(11, 0, halogen), 0),
                Ok(1),
                "halogen {halogen}"
            );
        }

        let mut singlet_halogen = neutral_metal.clone();
        singlet_halogen[1].radical = RADICAL_SINGLET as i8;
        assert_eq!(bIsMetalSalt(&singlet_halogen, 0), Ok(1));

        for (name, mutate) in [
            (
                "valence",
                (|atom: &mut inp_ATOM| atom.valence = 2) as fn(&mut inp_ATOM),
            ),
            ("chemical valence", |atom: &mut inp_ATOM| {
                atom.chem_bonds_valence = 2
            }),
            ("charge", |atom: &mut inp_ATOM| atom.charge = 1),
            ("radical", |atom: &mut inp_ATOM| atom.radical = 2),
            ("ordinary hydrogen", |atom: &mut inp_ATOM| atom.num_H = 1),
            ("isotope hydrogen", |atom: &mut inp_ATOM| {
                atom.num_iso_H[2] = 1
            }),
        ] {
            let mut rejected = neutral_metal.clone();
            mutate(&mut rejected[1]);
            assert_eq!(bIsMetalSalt(&rejected, 0), Ok(0), "{name}");
        }

        let oxygen_salt = oxygen_carbon_salt();
        assert_eq!(bIsMetalSalt(&oxygen_salt, 0), Ok(1));

        for (name, mutate) in [
            (
                "element",
                (|atom: &mut inp_ATOM| atom.el_number = 7) as fn(&mut inp_ATOM),
            ),
            ("ordinary hydrogen", |atom: &mut inp_ATOM| atom.num_H = 1),
            ("isotope hydrogen", |atom: &mut inp_ATOM| {
                atom.num_iso_H[0] = 1
            }),
            ("valence", |atom: &mut inp_ATOM| atom.valence = 1),
            ("charge", |atom: &mut inp_ATOM| atom.charge = -1),
            ("radical", |atom: &mut inp_ATOM| atom.radical = 2),
            ("chemical valence", |atom: &mut inp_ATOM| {
                atom.chem_bonds_valence = 1
            }),
        ] {
            let mut rejected = oxygen_salt.clone();
            mutate(&mut rejected[1]);
            assert_eq!(bIsMetalSalt(&rejected, 0), Ok(0), "oxygen {name}");
        }

        for (name, mutate) in [
            (
                "element",
                (|atom: &mut inp_ATOM| atom.el_number = 7) as fn(&mut inp_ATOM),
            ),
            ("implicit hydrogen", |atom: &mut inp_ATOM| atom.num_H = 1),
            ("chemical valence", |atom: &mut inp_ATOM| {
                atom.chem_bonds_valence = 3
            }),
            ("charge", |atom: &mut inp_ATOM| atom.charge = 1),
            ("radical", |atom: &mut inp_ATOM| atom.radical = 2),
            ("all single bonds", |atom: &mut inp_ATOM| atom.valence = 4),
        ] {
            let mut rejected = oxygen_salt.clone();
            mutate(&mut rejected[2]);
            assert_eq!(bIsMetalSalt(&rejected, 0), Ok(0), "carbon {name}");
        }

        let mut explicit_hydrogen = vec![
            ion_atom(11, 0),
            ion_atom(8, 0),
            ion_atom(6, 0),
            ion_atom(1, 0),
            ion_atom(6, 0),
        ];
        ion_bond(&mut explicit_hydrogen, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut explicit_hydrogen, 1, 2, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut explicit_hydrogen, 2, 3, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut explicit_hydrogen, 2, 4, BOND_TYPE_DOUBLE as u8);
        assert_eq!(bIsMetalSalt(&explicit_hydrogen, 0), Ok(0));

        let mut later_invalid = vec![ion_atom(25, 0), ion_atom(9, 0), ion_atom(6, 0)];
        ion_bond(&mut later_invalid, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut later_invalid, 0, 2, BOND_TYPE_SINGLE as u8);
        assert_eq!(bIsMetalSalt(&later_invalid, 0), Ok(0));

        let mut invalid_neighbor = neutral_metal;
        invalid_neighbor[0].neighbor[0] = u16::MAX;
        assert_eq!(
            bIsMetalSalt(&invalid_neighbor, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__strutil__disconnectmetalsalt__line_2612() {
        assert_eq!(
            DisconnectMetalSalt(&mut [], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            DisconnectMetalSalt(&mut [], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut disconnected = vec![ion_atom(11, i8::MAX)];
        disconnected[0].chem_bonds_valence = 7;
        assert_eq!(DisconnectMetalSalt(&mut disconnected, 0), Ok(0));
        assert_eq!(
            (
                disconnected[0].charge,
                disconnected[0].valence,
                disconnected[0].chem_bonds_valence
            ),
            (i8::MAX, 0, 0)
        );

        let mut terminal = vec![ion_atom(11, 0), ion_atom(8, 0)];
        ion_bond(&mut terminal, 0, 1, BOND_TYPE_SINGLE as u8);
        terminal[0].bond_stereo[0] = 7;
        terminal[1].bond_stereo[0] = -7;
        assert_eq!(DisconnectMetalSalt(&mut terminal, 0), Ok(1));
        assert_eq!(
            (
                terminal[0].charge,
                terminal[0].valence,
                terminal[0].chem_bonds_valence,
                terminal[0].neighbor[0],
                terminal[0].bond_stereo[0],
                terminal[0].bond_type[0]
            ),
            (1, 0, 0, 0, 0, 0)
        );
        assert_eq!(
            (
                terminal[1].charge,
                terminal[1].valence,
                terminal[1].chem_bonds_valence,
                terminal[1].neighbor[0],
                terminal[1].bond_stereo[0],
                terminal[1].bond_type[0]
            ),
            (-1, 0, 0, 0, 0, 0)
        );

        let mut metal_in_slot_zero = vec![ion_atom(11, 0), ion_atom(8, 0), ion_atom(6, 0)];
        ion_bond(&mut metal_in_slot_zero, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut metal_in_slot_zero, 1, 2, BOND_TYPE_DOUBLE as u8);
        metal_in_slot_zero[1].bond_stereo[1] = -8;
        assert_eq!(DisconnectMetalSalt(&mut metal_in_slot_zero, 0), Ok(1));
        assert_eq!(
            (
                metal_in_slot_zero[1].neighbor[0],
                metal_in_slot_zero[1].bond_stereo[0],
                metal_in_slot_zero[1].bond_type[0],
                metal_in_slot_zero[1].neighbor[1],
                metal_in_slot_zero[1].bond_stereo[1],
                metal_in_slot_zero[1].bond_type[1],
                metal_in_slot_zero[1].valence,
                metal_in_slot_zero[1].chem_bonds_valence,
                metal_in_slot_zero[1].charge
            ),
            (2, -8, BOND_TYPE_DOUBLE as u8, 0, 0, 0, 1, 2, -1)
        );

        let mut metal_in_slot_one = vec![ion_atom(11, 0), ion_atom(8, 0), ion_atom(6, 0)];
        ion_bond(&mut metal_in_slot_one, 1, 2, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut metal_in_slot_one, 0, 1, BOND_TYPE_SINGLE as u8);
        metal_in_slot_one[1].bond_stereo[0] = 9;
        metal_in_slot_one[1].bond_stereo[1] = -9;
        assert_eq!(DisconnectMetalSalt(&mut metal_in_slot_one, 0), Ok(1));
        assert_eq!(
            (
                metal_in_slot_one[1].neighbor[0],
                metal_in_slot_one[1].bond_stereo[0],
                metal_in_slot_one[1].bond_type[0],
                metal_in_slot_one[1].neighbor[1],
                metal_in_slot_one[1].bond_stereo[1],
                metal_in_slot_one[1].bond_type[1]
            ),
            (2, 9, BOND_TYPE_DOUBLE as u8, 0, 0, 0)
        );

        let mut multiple = vec![ion_atom(25, i8::MAX), ion_atom(9, 0), ion_atom(17, 0)];
        ion_bond(&mut multiple, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut multiple, 0, 2, BOND_TYPE_SINGLE as u8);
        assert_eq!(DisconnectMetalSalt(&mut multiple, 0), Ok(2));
        assert_eq!(multiple[0].charge, i8::MIN.wrapping_add(1));
        assert_eq!((multiple[1].charge, multiple[2].charge), (-1, -1));
        assert_eq!((multiple[1].valence, multiple[2].valence), (0, 0));

        let mut narrowing = vec![ion_atom(11, 0), ion_atom(8, 0)];
        narrowing[0].valence = 1;
        narrowing[0].neighbor[0] = 1;
        narrowing[1].valence = i8::MIN;
        narrowing[1].chem_bonds_valence = i8::MIN;
        assert_eq!(DisconnectMetalSalt(&mut narrowing, 0), Ok(1));
        assert_eq!(
            (narrowing[1].valence, narrowing[1].chem_bonds_valence),
            (i8::MAX, i8::MAX)
        );

        let mut partial = vec![ion_atom(25, 0), ion_atom(9, 0)];
        ion_bond(&mut partial, 0, 1, BOND_TYPE_SINGLE as u8);
        partial[0].neighbor[1] = u16::MAX;
        partial[0].bond_type[1] = BOND_TYPE_SINGLE as u8;
        partial[0].valence = 2;
        partial[0].chem_bonds_valence = 2;
        assert_eq!(
            DisconnectMetalSalt(&mut partial, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            (
                partial[0].charge,
                partial[0].neighbor[0],
                partial[0].neighbor[1],
                partial[0].valence,
                partial[0].chem_bonds_valence,
                partial[1].charge,
                partial[1].valence
            ),
            (1, 0, u16::MAX, 2, 2, -1, 0)
        );
    }

    #[test]
    fn source_port__strutil__disconnectsalts__line_2668() {
        fn original_with_atoms(
            heap: &mut SourceHeap,
            atoms: Vec<inp_ATOM>,
            number_of_bonds: i32,
        ) -> ORIG_ATOM_DATA {
            let number_of_atoms = atoms.len() as i32;
            ORIG_ATOM_DATA {
                at: heap.allocate_model_storage(atoms).unwrap(),
                num_inp_atoms: number_of_atoms,
                num_inp_bonds: number_of_bonds,
                ..ORIG_ATOM_DATA::default()
            }
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty = ORIG_ATOM_DATA::default();
        assert_eq!(DisconnectSalts(&mut empty_heap, &mut empty, 1), Ok(0));
        empty.num_inp_atoms = -1;
        assert_eq!(DisconnectSalts(&mut empty_heap, &mut empty, 1), Ok(0));
        empty.num_inp_atoms = 1;
        assert_eq!(
            DisconnectSalts(&mut empty_heap, &mut empty, 1),
            Err(SourceHeapError::NullPointer)
        );

        let mut rejected = vec![
            ion_atom(6, 0),
            ion_atom(6, 0),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        rejected[1].valence = 1;
        rejected[1].chem_bonds_valence = 2;
        rejected[2].valence = 1;
        rejected[2].chem_bonds_valence = 1;
        rejected[2].radical = 2;
        rejected[3].valence = 1;
        rejected[3].chem_bonds_valence = 1;
        rejected[3].radical = RADICAL_SINGLET as i8;
        let mut rejected_heap = SourceHeap::default();
        let mut rejected_original = original_with_atoms(&mut rejected_heap, rejected.clone(), 17);
        assert_eq!(
            DisconnectSalts(&mut rejected_heap, &mut rejected_original, 1),
            Ok(0)
        );
        assert_eq!(rejected_original.num_inp_bonds, 17);
        assert_eq!(
            rejected_heap
                .slice(rejected_original.at.as_const())
                .unwrap(),
            rejected
        );

        let mut ammonium = vec![ion_atom(7, 1), ion_atom(8, -1), ion_atom(6, 0)];
        ammonium[0].num_H = 4;
        ion_bond(&mut ammonium, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut ammonium, 1, 2, BOND_TYPE_SINGLE as u8);
        let mut ammonium_dry_heap = SourceHeap::default();
        let mut ammonium_dry = original_with_atoms(&mut ammonium_dry_heap, ammonium.clone(), 2);
        assert_eq!(
            DisconnectSalts(&mut ammonium_dry_heap, &mut ammonium_dry, 0),
            Ok(1)
        );
        assert_eq!(ammonium_dry.num_inp_bonds, 2);
        assert_eq!(
            ammonium_dry_heap.slice(ammonium_dry.at.as_const()).unwrap(),
            ammonium
        );

        let mut ammonium_heap = SourceHeap::default();
        let mut ammonium_original = original_with_atoms(&mut ammonium_heap, ammonium, 2);
        assert_eq!(
            DisconnectSalts(&mut ammonium_heap, &mut ammonium_original, -1),
            Ok(1)
        );
        assert_eq!(ammonium_original.num_inp_bonds, 1);
        let ammonium_output = ammonium_heap
            .slice(ammonium_original.at.as_const())
            .unwrap();
        assert_eq!(
            (
                ammonium_output[0].charge,
                ammonium_output[0].valence,
                ammonium_output[1].charge,
                ammonium_output[1].num_H
            ),
            (0, 0, 0, 1)
        );

        let mut metal = vec![ion_atom(11, 0), ion_atom(9, 0)];
        ion_bond(&mut metal, 0, 1, BOND_TYPE_SINGLE as u8);
        let mut metal_dry_heap = SourceHeap::default();
        let mut metal_dry = original_with_atoms(&mut metal_dry_heap, metal.clone(), 1);
        assert_eq!(
            DisconnectSalts(&mut metal_dry_heap, &mut metal_dry, 0),
            Ok(1)
        );
        assert_eq!(metal_dry.num_inp_bonds, 1);
        assert_eq!(
            metal_dry_heap.slice(metal_dry.at.as_const()).unwrap(),
            metal
        );

        let mut metal_heap = SourceHeap::default();
        let mut metal_original = original_with_atoms(&mut metal_heap, metal, i32::MIN);
        assert_eq!(
            DisconnectSalts(&mut metal_heap, &mut metal_original, 1),
            Ok(1)
        );
        assert_eq!(metal_original.num_inp_bonds, i32::MAX);
        let metal_output = metal_heap.slice(metal_original.at.as_const()).unwrap();
        assert_eq!(
            (
                metal_output[0].charge,
                metal_output[0].valence,
                metal_output[1].charge,
                metal_output[1].valence
            ),
            (1, 0, -1, 0)
        );

        let mut combined = vec![
            ion_atom(7, 1),
            ion_atom(8, -1),
            ion_atom(6, 0),
            ion_atom(11, 0),
            ion_atom(17, 0),
        ];
        combined[0].num_H = 4;
        ion_bond(&mut combined, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut combined, 1, 2, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut combined, 3, 4, BOND_TYPE_SINGLE as u8);
        let mut combined_heap = SourceHeap::default();
        let mut combined_original = original_with_atoms(&mut combined_heap, combined, 3);
        assert_eq!(
            DisconnectSalts(&mut combined_heap, &mut combined_original, 0),
            Ok(2)
        );
        assert_eq!(combined_original.num_inp_bonds, 3);

        let mut partial_atoms = vec![ion_atom(11, 0), ion_atom(9, 0), ion_atom(u8::MAX, 0)];
        ion_bond(&mut partial_atoms, 0, 1, BOND_TYPE_SINGLE as u8);
        partial_atoms[2].valence = 1;
        partial_atoms[2].chem_bonds_valence = 1;
        partial_atoms[2].neighbor[0] = 2;
        let mut partial_heap = SourceHeap::default();
        let mut partial_original = original_with_atoms(&mut partial_heap, partial_atoms, 8);
        assert_eq!(
            DisconnectSalts(&mut partial_heap, &mut partial_original, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(partial_original.num_inp_bonds, 7);
        let partial_output = partial_heap.slice(partial_original.at.as_const()).unwrap();
        assert_eq!(
            (
                partial_output[0].charge,
                partial_output[0].valence,
                partial_output[1].charge,
                partial_output[1].valence
            ),
            (1, 0, -1, 0)
        );
    }

    #[test]
    fn source_port__strutil__bismetaltodisconnect__line_2719() {
        assert_eq!(
            bIsMetalToDisconnect(&[], -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            bIsMetalToDisconnect(&[], 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            bIsMetalToDisconnect(&[ion_atom(u8::MAX, 0)], 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut nonmetal = ion_atom(6, 0);
        nonmetal.chem_bonds_valence = 4;
        assert_eq!(bIsMetalToDisconnect(&[nonmetal], 0, 0), Ok(0));

        let disconnected_metal = ion_atom(11, 0);
        assert_eq!(
            bIsMetalToDisconnect(&[disconnected_metal.clone()], 0, 0),
            Ok(0)
        );
        assert_eq!(
            bIsMetalToDisconnect(&[disconnected_metal], 0, i32::MIN),
            Ok(0)
        );

        let mut cancelled = ion_atom(11, 0);
        cancelled.num_H = 1;
        cancelled.num_iso_H = [1, -1, 0];
        cancelled.chem_bonds_valence = -1;
        assert_eq!(bIsMetalToDisconnect(&[cancelled], 0, 1), Ok(0));

        let mut ordinary_hydrogen = ion_atom(11, 0);
        ordinary_hydrogen.num_H = 1;
        assert_eq!(bIsMetalToDisconnect(&[ordinary_hydrogen], 0, 0), Ok(1));
        let mut isotope_hydrogen = ion_atom(11, 0);
        isotope_hydrogen.num_iso_H = [0, 0, 1];
        assert_eq!(bIsMetalToDisconnect(&[isotope_hydrogen], 0, 1), Ok(1));

        let mut normal_na = ion_atom(11, 0);
        normal_na.chem_bonds_valence = 1;
        assert_eq!(bIsMetalToDisconnect(&[normal_na.clone()], 0, 0), Ok(1));
        assert_eq!(bIsMetalToDisconnect(&[normal_na], 0, 1), Ok(1));

        let mut metal2 = ion_atom(25, 0);
        metal2.chem_bonds_valence = 2;
        assert_eq!(bIsMetalToDisconnect(&[metal2], 0, 1), Ok(1));

        for charge in [i8::MIN, -2, 2, i8::MAX] {
            let mut charged = ion_atom(12, charge);
            charged.chem_bonds_valence = 1;
            assert_eq!(
                bIsMetalToDisconnect(&[charged.clone()], 0, 1),
                Ok(1),
                "charge {charge}"
            );
            assert_eq!(
                bIsMetalToDisconnect(&[charged], 0, i32::MIN),
                Ok(1),
                "truthy check, charge {charge}"
            );
        }

        for element in 0_u8..=120 {
            for charge in [i8::MIN, -2, -1, 0, 1, 2, i8::MAX] {
                for chemical_valence in [-2_i8, -1, 0, 1, 2, 3, 4, 6, i8::MAX] {
                    let mut atom = ion_atom(element, charge);
                    atom.chem_bonds_valence = chemical_valence;
                    assert_ne!(
                        bIsMetalToDisconnect(&[atom], 0, 1).unwrap(),
                        2,
                        "source loop must remain unreachable: element={element} charge={charge} valence={chemical_valence}"
                    );
                }
            }
        }
    }

    #[test]
    fn source_port__strutil__bmaydisconnectmetals__line_2762() {
        fn original_with_atoms(heap: &mut SourceHeap, atoms: Vec<inp_ATOM>) -> ORIG_ATOM_DATA {
            let number_of_atoms = atoms.len() as i32;
            ORIG_ATOM_DATA {
                at: heap.allocate_model_storage(atoms).unwrap(),
                num_inp_atoms: number_of_atoms,
                bDisconnectCoord: 77,
                ..ORIG_ATOM_DATA::default()
            }
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty = ORIG_ATOM_DATA {
            bDisconnectCoord: 77,
            ..ORIG_ATOM_DATA::default()
        };
        let mut empty_flags = 0x40_u64;
        assert_eq!(
            bMayDisconnectMetals(&mut empty_heap, &mut empty, 1, Some(&mut empty_flags)),
            Ok(0)
        );
        assert_eq!((empty.bDisconnectCoord, empty_flags), (0, 0x40));
        empty.num_inp_atoms = -1;
        empty.bDisconnectCoord = 88;
        assert_eq!(
            bMayDisconnectMetals(&mut empty_heap, &mut empty, 0, None),
            Ok(0)
        );
        assert_eq!(empty.bDisconnectCoord, 0);
        empty.num_inp_atoms = 1;
        assert_eq!(
            bMayDisconnectMetals(&mut empty_heap, &mut empty, 0, None),
            Err(SourceHeapError::NullPointer)
        );

        let mut protected = vec![
            ion_atom(7, 1),
            ion_atom(8, -1),
            ion_atom(6, 0),
            ion_atom(11, 0),
            ion_atom(9, 0),
        ];
        protected[0].num_H = 4;
        ion_bond(&mut protected, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut protected, 1, 2, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut protected, 3, 4, BOND_TYPE_SINGLE as u8);
        let mut protected_heap = SourceHeap::default();
        let mut protected_original = original_with_atoms(&mut protected_heap, protected);
        let mut protected_flags = 5_u64;
        assert_eq!(
            bMayDisconnectMetals(
                &mut protected_heap,
                &mut protected_original,
                1,
                Some(&mut protected_flags)
            ),
            Ok(0)
        );
        assert_eq!(
            (protected_original.bDisconnectCoord, protected_flags),
            (0, 5)
        );

        let mut metal_carbon = vec![ion_atom(11, 0), ion_atom(6, 0)];
        ion_bond(&mut metal_carbon, 0, 1, BOND_TYPE_SINGLE as u8);
        let mut carbon_heap = SourceHeap::default();
        let mut carbon_original = original_with_atoms(&mut carbon_heap, metal_carbon);
        let mut carbon_flags = 0x20_u64;
        assert_eq!(
            bMayDisconnectMetals(
                &mut carbon_heap,
                &mut carbon_original,
                1,
                Some(&mut carbon_flags)
            ),
            Ok(1)
        );
        assert_eq!((carbon_original.bDisconnectCoord, carbon_flags), (1, 0x20));

        let mut implicit = ion_atom(12, 0);
        implicit.num_H = 2;
        implicit.num_iso_H = [1, 2, 3];
        let mut implicit_heap = SourceHeap::default();
        let mut implicit_original = original_with_atoms(&mut implicit_heap, vec![implicit]);
        assert_eq!(
            bMayDisconnectMetals(&mut implicit_heap, &mut implicit_original, 0, None),
            Ok(1)
        );
        assert_eq!(implicit_original.bDisconnectCoord, 9);

        let mut radical = vec![ion_atom(25, 0), ion_atom(6, 0)];
        ion_bond(&mut radical, 0, 1, BOND_TYPE_SINGLE as u8);
        radical[0].radical = 2;
        let mut radical_heap = SourceHeap::default();
        let mut radical_original = original_with_atoms(&mut radical_heap, radical);
        assert_eq!(
            bMayDisconnectMetals(&mut radical_heap, &mut radical_original, 1, None),
            Ok(1)
        );
        assert_eq!(radical_original.bDisconnectCoord, 1);

        let mut multiple_bond = vec![ion_atom(25, 0), ion_atom(6, 0)];
        ion_bond(&mut multiple_bond, 0, 1, BOND_TYPE_DOUBLE as u8);
        let mut multiple_heap = SourceHeap::default();
        let mut multiple_original = original_with_atoms(&mut multiple_heap, multiple_bond);
        assert_eq!(
            bMayDisconnectMetals(&mut multiple_heap, &mut multiple_original, 1, None),
            Ok(1)
        );
        assert_eq!(multiple_original.bDisconnectCoord, 1);

        let mut nonmetal = ion_atom(6, 0);
        nonmetal.num_H = 4;
        let mut nonmetal_heap = SourceHeap::default();
        let mut nonmetal_original = original_with_atoms(&mut nonmetal_heap, vec![nonmetal]);
        assert_eq!(
            bMayDisconnectMetals(&mut nonmetal_heap, &mut nonmetal_original, 1, None),
            Ok(0)
        );
        assert_eq!(nonmetal_original.bDisconnectCoord, 0);

        let mut combined = vec![ion_atom(11, 0), ion_atom(12, 0)];
        combined[0].num_H = 1;
        combined[1].num_iso_H = [0, 1, 1];
        let mut combined_heap = SourceHeap::default();
        let mut combined_original = original_with_atoms(&mut combined_heap, combined);
        let mut combined_flags = u64::MAX ^ TG_FLAG_CHECK_VALENCE_COORD_DONE as u64;
        assert_eq!(
            bMayDisconnectMetals(
                &mut combined_heap,
                &mut combined_original,
                i32::MIN,
                Some(&mut combined_flags)
            ),
            Ok(2)
        );
        assert_eq!(combined_original.bDisconnectCoord, 4);
        assert_eq!(combined_flags & TG_FLAG_CHECK_VALENCE_COORD_DONE as u64, 0);

        let mut malformed_heap = SourceHeap::default();
        let mut malformed = original_with_atoms(
            &mut malformed_heap,
            vec![inp_ATOM {
                el_number: u8::MAX,
                num_H: 1,
                ..inp_ATOM::default()
            }],
        );
        assert_eq!(
            bMayDisconnectMetals(&mut malformed_heap, &mut malformed, 0, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(malformed.bDisconnectCoord, 77);
    }

    #[test]
    fn source_port__strutil__dist3d__line_3245() {
        let origin = inp_ATOM::default();
        let point = inp_ATOM {
            x: 3.0,
            y: -4.0,
            z: 12.0,
            ..inp_ATOM::default()
        };
        assert_eq!(dist3D(&point, &origin).to_bits(), 13.0_f64.to_bits());
        assert_eq!(dist3D(&origin, &point).to_bits(), 13.0_f64.to_bits());

        let signed_zero = inp_ATOM {
            x: -0.0,
            y: 0.0,
            z: -0.0,
            ..inp_ATOM::default()
        };
        assert_eq!(dist3D(&signed_zero, &origin).to_bits(), 0.0_f64.to_bits());

        let subnormal = inp_ATOM {
            x: f64::from_bits(1),
            ..inp_ATOM::default()
        };
        assert_eq!(dist3D(&subnormal, &origin).to_bits(), 0.0_f64.to_bits());

        let largest = inp_ATOM {
            x: f64::MAX,
            ..inp_ATOM::default()
        };
        assert_eq!(dist3D(&largest, &origin), f64::INFINITY);

        let positive_infinity = inp_ATOM {
            y: f64::INFINITY,
            ..inp_ATOM::default()
        };
        let negative_infinity = inp_ATOM {
            y: f64::NEG_INFINITY,
            ..inp_ATOM::default()
        };
        assert_eq!(dist3D(&positive_infinity, &origin), f64::INFINITY);
        assert_eq!(
            dist3D(&positive_infinity, &negative_infinity),
            f64::INFINITY
        );
        assert!(dist3D(&positive_infinity, &positive_infinity).is_nan());

        let nan = inp_ATOM {
            z: f64::NAN,
            ..inp_ATOM::default()
        };
        assert!(dist3D(&nan, &origin).is_nan());
    }

    #[test]
    fn source_port__strutil__getmindistdistribution__line_3264() {
        fn distribution(
            first: (f64, f64),
            second: (f64, f64),
            components: (u16, u16, u16),
            all_components: i32,
        ) -> (f64, [f64; 8]) {
            let mut atoms = vec![inp_ATOM::default(); 3];
            atoms[0].component = components.0;
            atoms[1].x = first.0;
            atoms[1].y = first.1;
            atoms[1].component = components.1;
            atoms[2].x = second.0;
            atoms[2].y = second.1;
            atoms[2].component = components.2;
            ion_bond(&mut atoms, 1, 2, BOND_TYPE_SINGLE as u8);
            let mut output = [f64::NAN; 8];
            let average =
                GetMinDistDistribution(&atoms, 3, 0, -1, all_components, &mut output, 8).unwrap();
            (average, output)
        }

        let mut invalid = [];
        assert_eq!(
            GetMinDistDistribution(&[], 0, 0, -1, 0, &mut invalid, 0),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        let mut too_short = [0.0; 1];
        assert_eq!(
            GetMinDistDistribution(&[], 0, 0, -1, 0, &mut too_short, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut empty_output = [7.0; 4];
        assert_eq!(
            GetMinDistDistribution(&[], 0, i32::MIN, i32::MAX, 0, &mut empty_output, 4),
            Ok(0.0)
        );
        assert_eq!(empty_output, [1.0e30; 4]);
        assert_eq!(
            GetMinDistDistribution(&[inp_ATOM::default()], 1, -1, -1, 1, &mut empty_output, 4,),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let (average, ordinary) = distribution((-1.0, 1.0), (1.0, 1.0), (1, 1, 1), 0);
        assert_eq!(average.to_bits(), 2.0_f64.to_bits());
        assert_eq!(ordinary[2].to_bits(), 1.0_f64.to_bits());
        assert!(ordinary.iter().any(|&distance| distance == 1.0e30));
        let (reversed_average, reversed) = distribution((1.0, 1.0), (-1.0, 1.0), (1, 1, 1), 0);
        assert_eq!(reversed_average.to_bits(), average.to_bits());
        assert_eq!(
            reversed.map(f64::to_bits),
            ordinary.map(f64::to_bits),
            "negative cross-product input is normalized by swapping endpoints"
        );

        let (filtered_average, filtered) = distribution((-1.0, 1.0), (1.0, 1.0), (1, 2, 2), 0);
        assert_eq!(filtered_average.to_bits(), 0.0_f64.to_bits());
        assert_eq!(filtered, [1.0e30; 8]);
        let (all_average, all) = distribution((-1.0, 1.0), (1.0, 1.0), (1, 2, 2), -1);
        assert_eq!(all_average.to_bits(), 2.0_f64.to_bits());
        assert_ne!(all, [1.0e30; 8]);

        let mut target_bond = vec![inp_ATOM::default(); 2];
        target_bond[1].x = 3.0;
        target_bond[1].y = 4.0;
        ion_bond(&mut target_bond, 0, 1, BOND_TYPE_SINGLE as u8);
        let mut target_distances = [1.0e30; 4];
        assert_eq!(
            GetMinDistDistribution(&target_bond, 2, 0, -1, 1, &mut target_distances, 4,)
                .unwrap()
                .to_bits(),
            5.0_f64.to_bits()
        );
        assert_eq!(
            target_distances
                .iter()
                .filter(|&&distance| distance != 1.0e30)
                .count(),
            1
        );
        let mut hydrogen_excluded = [7.0; 4];
        assert_eq!(
            GetMinDistDistribution(&target_bond, 2, 0, 1, 1, &mut hydrogen_excluded, 4,),
            Ok(0.0)
        );
        assert_eq!(hydrogen_excluded, [1.0e30; 4]);

        let (_, positive_dot) = distribution((1.0e-6, 10_000.0), (2.0e-6, 0.0), (1, 1, 1), 1);
        assert_eq!(
            positive_dot
                .iter()
                .filter(|&&distance| distance != 1.0e30)
                .count(),
            1
        );
        let (_, negative_dot) = distribution((-1.0e-6, 10_000.0), (2.0e-6, 0.0), (1, 1, 1), 1);
        assert_eq!(
            negative_dot
                .iter()
                .filter(|&&distance| distance != 1.0e30)
                .count(),
            2
        );
        let (_, zero_dot) = distribution((0.0, 10_000.0), (2.0e-6, 0.0), (1, 1, 1), 1);
        assert_eq!(zero_dot, [1.0e30; 8]);

        let (_, both_coincident) = distribution((-5.0e-8, 0.0), (5.0e-8, 0.0), (1, 1, 1), 1);
        assert_eq!(both_coincident, [1.0e30; 8]);
        let (_, one_coincident) = distribution((2.0e-6, 0.0), (0.0, 0.0), (1, 1, 1), 1);
        assert_eq!(
            one_coincident
                .iter()
                .filter(|&&distance| distance != 1.0e30)
                .count(),
            1
        );

        let mut zero_length_atoms = vec![inp_ATOM::default(); 3];
        zero_length_atoms[1].x = 2.0;
        zero_length_atoms[2].x = 2.0;
        ion_bond(&mut zero_length_atoms, 1, 2, BOND_TYPE_SINGLE as u8);
        let mut zero_length_output = [1.0e30; 4];
        assert_eq!(
            GetMinDistDistribution(&zero_length_atoms, 3, 0, -1, 1, &mut zero_length_output, 4,),
            Ok(0.0)
        );
        assert_eq!(
            zero_length_output
                .iter()
                .filter(|&&distance| distance != 1.0e30)
                .count(),
            1
        );

        let mut malformed = vec![inp_ATOM::default(); 2];
        malformed[1].valence = 1;
        malformed[1].neighbor[0] = 0;
        let mut malformed_output = [0.0; 4];
        assert_eq!(
            GetMinDistDistribution(&malformed, 3, 0, -1, 1, &mut malformed_output, 4,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__strutil__move_explicit_hcation__line_3480() {
        assert_eq!(
            move_explicit_Hcation(&mut [], 0, -1, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            move_explicit_Hcation(&mut [inp_ATOM::default()], 1, 0, 1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut isolated = vec![ion_atom(8, -1), ion_atom(1, 1)];
        isolated[0].component = 7;
        isolated[0].x = 3.0;
        isolated[0].y = -2.0;
        isolated[0].z = 5.0;
        isolated[1].bond_stereo[0] = 9;
        assert_eq!(move_explicit_Hcation(&mut isolated, 2, 0, 1, 1), Ok(1));
        assert_eq!(
            (
                isolated[0].charge,
                isolated[0].valence,
                isolated[0].chem_bonds_valence,
                isolated[0].neighbor[0],
                isolated[1].charge,
                isolated[1].valence,
                isolated[1].chem_bonds_valence,
                isolated[1].component,
                isolated[1].neighbor[0],
                isolated[1].bond_stereo[0],
            ),
            (0, 1, 1, 1, 0, 1, 1, 7, 0, 0)
        );
        assert_eq!(
            (isolated[1].x, isolated[1].y, isolated[1].z),
            (3.0, -2.0, 5.0)
        );

        let mut connected = vec![ion_atom(8, -1), ion_atom(7, 1), ion_atom(1, 0)];
        connected[0].component = 3;
        connected[2].x = 1.0;
        ion_bond(&mut connected, 1, 2, BOND_TYPE_SINGLE as u8);
        connected[1].bond_stereo[0] = -5;
        connected[2].bond_stereo[0] = 5;
        assert_eq!(move_explicit_Hcation(&mut connected, 3, 0, 2, 1), Ok(1));
        assert_eq!(
            (
                connected[0].charge,
                connected[0].valence,
                connected[1].charge,
                connected[1].valence,
                connected[1].chem_bonds_valence,
                connected[2].valence,
                connected[2].neighbor[0],
                connected[2].bond_stereo[0],
                connected[2].component,
            ),
            (0, 1, 0, 0, 0, 1, 0, 0, 3)
        );
        assert_eq!(
            (connected[2].x, connected[2].y, connected[2].z),
            (-1.0, 0.0, 0.0)
        );

        let mut missing = vec![
            ion_atom(8, -1),
            ion_atom(7, 1),
            ion_atom(1, 0),
            ion_atom(6, 0),
        ];
        missing[2].valence = 1;
        missing[2].chem_bonds_valence = 1;
        missing[2].neighbor[0] = 1;
        missing[2].bond_type[0] = BOND_TYPE_SINGLE as u8;
        missing[1].valence = 1;
        missing[1].chem_bonds_valence = 1;
        missing[1].neighbor[0] = 3;
        missing[1].bond_type[0] = BOND_TYPE_SINGLE as u8;
        let missing_before = missing.clone();
        assert_eq!(move_explicit_Hcation(&mut missing, 4, 0, 2, 1), Ok(0));
        assert_eq!(missing, missing_before);

        let mut directed = vec![ion_atom(8, 0), ion_atom(6, 0), ion_atom(1, 1)];
        directed[1].x = 1.0;
        ion_bond(&mut directed, 0, 1, BOND_TYPE_SINGLE as u8);
        assert_eq!(move_explicit_Hcation(&mut directed, 3, 0, 2, 1), Ok(1));
        assert_eq!(
            (directed[2].x.to_bits(), directed[2].y.to_bits()),
            ((-1.0_f64).to_bits(), 0.0_f64.to_bits())
        );
        assert_eq!((directed[0].valence, directed[0].neighbor[1]), (2, 2));

        let mut fallback = vec![
            ion_atom(8, 0),
            ion_atom(1, 1),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        fallback[2].x = -1.0;
        fallback[3].x = 1.0;
        ion_bond(&mut fallback, 2, 3, BOND_TYPE_SINGLE as u8);
        assert_eq!(move_explicit_Hcation(&mut fallback, 4, 0, 1, 1), Ok(1));
        let fallback_angle: f64 = (2.0 * 3.14159265358979323846 / 20.0) * 15.0;
        assert_eq!(
            (fallback[1].x.to_bits(), fallback[1].y.to_bits()),
            (
                (2.0_f64 * fallback_angle.cos()).to_bits(),
                (2.0_f64 * fallback_angle.sin()).to_bits(),
            )
        );

        let mut obstructed = vec![
            ion_atom(8, 0),
            ion_atom(6, 0),
            ion_atom(1, 1),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        obstructed[1].x = 1.0;
        ion_bond(&mut obstructed, 0, 1, BOND_TYPE_SINGLE as u8);
        obstructed[3].x = -0.5;
        obstructed[3].y = -0.2;
        obstructed[4].x = -0.5;
        obstructed[4].y = 0.2;
        ion_bond(&mut obstructed, 3, 4, BOND_TYPE_SINGLE as u8);
        assert_eq!(move_explicit_Hcation(&mut obstructed, 5, 0, 2, 1), Ok(1));
        assert_ne!(
            (obstructed[2].x.to_bits(), obstructed[2].y.to_bits()),
            ((-1.0_f64).to_bits(), 0.0_f64.to_bits())
        );
        assert!(obstructed[2].x.is_finite() && obstructed[2].y.is_finite());

        let mut full = vec![ion_atom(6, -1)];
        full.extend((0..MAXVAL as usize).map(|_| ion_atom(6, 0)));
        full.push(ion_atom(1, 1));
        for neighbor in 1..=MAXVAL as usize {
            ion_bond(&mut full, 0, neighbor, BOND_TYPE_SINGLE as u8);
        }
        let hydrogen = MAXVAL as usize + 1;
        full[0].component = 11;
        let target_before = full[0].clone();
        assert_eq!(
            move_explicit_Hcation(&mut full, hydrogen as i32 + 1, 0, hydrogen as i32, 1),
            Ok(1)
        );
        assert_eq!(full[0].valence, MAXVAL as i8);
        assert_eq!(full[0].neighbor, target_before.neighbor);
        assert_eq!(full[0].bond_type, target_before.bond_type);
        assert_eq!(
            (full[hydrogen].component, full[hydrogen].neighbor[0]),
            (11, 0)
        );

        let mut malformed = vec![ion_atom(8, 0), ion_atom(1, 1)];
        malformed[1].valence = 1;
        malformed[1].neighbor[0] = u16::MAX;
        let malformed_before = malformed.clone();
        assert_eq!(
            move_explicit_Hcation(&mut malformed, 2, 0, 1, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(malformed, malformed_before);
    }

    #[test]
    fn source_port__strutil__disconnectmetals__line_2865() {
        fn original(heap: &mut SourceHeap, atoms: Vec<inp_ATOM>) -> ORIG_ATOM_DATA {
            let count = atoms.len() as i32;
            ORIG_ATOM_DATA {
                at: heap.allocate_model_storage(atoms).unwrap(),
                num_inp_atoms: count,
                ..ORIG_ATOM_DATA::default()
            }
        }

        for successful_allocations in [0, 1] {
            let mut heap = SourceHeap::default();
            let mut input = original(&mut heap, vec![ion_atom(11, 0), ion_atom(6, 0)]);
            ion_bond(
                heap.slice_mut(input.at).unwrap(),
                0,
                1,
                BOND_TYPE_SINGLE as u8,
            );
            let original_pointer = input.at;
            heap.fail_after_allocations(successful_allocations);
            assert_eq!(DisconnectMetals(&mut heap, &mut input, 0, None), Ok(-6));
            assert_eq!(input.at, original_pointer);
            assert_eq!(heap.live_allocation_count(), 1);
        }

        let mut no_change_heap = SourceHeap::default();
        let unchanged_atoms = vec![ion_atom(6, 0)];
        let mut no_change = original(&mut no_change_heap, unchanged_atoms.clone());
        let original_pointer = no_change.at;
        assert_eq!(
            DisconnectMetals(&mut no_change_heap, &mut no_change, 1, None),
            Ok(-6)
        );
        assert_eq!(no_change.at, original_pointer);
        assert_eq!(
            no_change_heap.slice(no_change.at.as_const()).unwrap(),
            unchanged_atoms
        );
        assert_eq!(no_change_heap.live_allocation_count(), 1);

        let mut protected_heap = SourceHeap::default();
        let mut protected_atoms = vec![ion_atom(11, 0), ion_atom(9, 0)];
        ion_bond(&mut protected_atoms, 0, 1, BOND_TYPE_SINGLE as u8);
        let mut protected = original(&mut protected_heap, protected_atoms.clone());
        assert_eq!(
            DisconnectMetals(&mut protected_heap, &mut protected, 0, None),
            Ok(-6)
        );
        assert_eq!(
            protected_heap.slice(protected.at.as_const()).unwrap(),
            protected_atoms
        );

        let mut ligand_heap = SourceHeap::default();
        let mut ligand_atoms = vec![ion_atom(11, 0), ion_atom(6, 0)];
        ion_bond(&mut ligand_atoms, 0, 1, BOND_TYPE_SINGLE as u8);
        ligand_atoms[0].component = 1;
        ligand_atoms[1].component = 2;
        let old_components = ligand_heap.allocate_model_storage(vec![9_u16, 10]).unwrap();
        let mut ligand = original(&mut ligand_heap, ligand_atoms);
        ligand.nOldCompNumber = old_components;
        let old_atom_pointer = ligand.at;
        assert_eq!(
            DisconnectMetals(&mut ligand_heap, &mut ligand, 0, None),
            Ok(1)
        );
        assert_ne!(ligand.at, old_atom_pointer);
        assert_eq!(ligand.num_inp_atoms, 2);
        assert_eq!(
            ligand_heap.slice(old_components.as_const()).unwrap(),
            [0, 0]
        );
        let output = ligand_heap.slice(ligand.at.as_const()).unwrap();
        assert_eq!((output[0].valence, output[1].valence), (0, 0));
        assert_eq!(ligand_heap.live_allocation_count(), 2);

        let mut mismatch_heap = SourceHeap::default();
        let mut implicit_atom = ion_atom(11, 0);
        implicit_atom.num_H = 1;
        let mut mismatch = original(&mut mismatch_heap, vec![implicit_atom.clone()]);
        let mismatch_pointer = mismatch.at;
        assert_eq!(
            DisconnectMetals(&mut mismatch_heap, &mut mismatch, 0, None),
            Ok(-6)
        );
        assert_eq!(mismatch.at, mismatch_pointer);
        assert_eq!(
            mismatch_heap.slice(mismatch.at.as_const()).unwrap(),
            [implicit_atom.clone()]
        );

        let mut implicit_heap = SourceHeap::default();
        let mut implicit = original(&mut implicit_heap, vec![implicit_atom]);
        implicit.bDisconnectCoord = 2;
        assert_eq!(
            DisconnectMetals(&mut implicit_heap, &mut implicit, 0, None),
            Ok(1)
        );
        assert_eq!(implicit.num_inp_atoms, 2);
        let expanded = implicit_heap.slice(implicit.at.as_const()).unwrap();
        assert_eq!((expanded[0].num_H, expanded[0].valence), (0, 0));
        assert_eq!(
            (
                expanded[1].el_number,
                expanded[1].iso_atw_diff,
                expanded[1].orig_at_number,
                expanded[1].valence,
            ),
            (1, 0, 2, 0)
        );

        let mut metals_heap = SourceHeap::default();
        let mut metals_atoms = vec![ion_atom(11, 0), ion_atom(12, 0)];
        ion_bond(&mut metals_atoms, 0, 1, BOND_TYPE_SINGLE as u8);
        let mut metals = original(&mut metals_heap, metals_atoms);
        assert_eq!(
            DisconnectMetals(&mut metals_heap, &mut metals, 0, None),
            Ok(1)
        );
        let metals_output = metals_heap.slice(metals.at.as_const()).unwrap();
        assert_eq!((metals_output[0].valence, metals_output[1].valence), (0, 0));
    }

    #[test]
    fn source_port__strutil__get_iat_number__line_4024() {
        let mapped = [
            (1, 0),
            (6, 1),
            (7, 2),
            (15, 3),
            (8, 4),
            (16, 5),
            (34, 6),
            (52, 7),
            (9, 8),
            (17, 9),
            (35, 10),
            (53, 11),
        ];
        for (element_number, ion_atom_type) in mapped {
            assert_eq!(get_iat_number(element_number), ion_atom_type);
        }

        for element_number in [
            i32::MIN,
            -1,
            0,
            2,
            5,
            10,
            14,
            18,
            33,
            36,
            51,
            54,
            118,
            i32::MAX,
        ] {
            assert_eq!(get_iat_number(element_number), -1);
        }
    }

    #[test]
    fn source_port__strutil__bheteroatommayhavexchgisoh__line_4047() {
        assert_eq!(
            bHeteroAtomMayHaveXchgIsoH(&[], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            bHeteroAtomMayHaveXchgIsoH(&[], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let carbon = ion_atom(6, 0);
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[carbon], 0), Ok(0));

        for charge in [-2_i8, 2] {
            let atom = ion_atom(8, charge);
            assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[atom], 0), Ok(0));
        }
        for radical in [RADICAL_DOUBLET as i8, 3, i8::MAX] {
            let mut atom = ion_atom(8, 0);
            atom.chem_bonds_valence = 2;
            atom.radical = radical;
            assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[atom], 0), Ok(0));
        }

        for (element, charge, chemical_valence) in [
            (7, -1, 2),
            (7, 0, 3),
            (7, 1, 4),
            (15, -1, 2),
            (15, 0, 3),
            (15, 1, 4),
            (8, -1, 1),
            (8, 0, 2),
            (8, 1, 3),
            (16, -1, 1),
            (16, 0, 2),
            (16, 1, 3),
            (34, -1, 1),
            (34, 0, 2),
            (34, 1, 3),
            (52, -1, 1),
            (52, 0, 2),
            (52, 1, 3),
        ] {
            let mut atom = ion_atom(element, charge);
            atom.chem_bonds_valence = chemical_valence;
            atom.radical = RADICAL_SINGLET as i8;
            assert_eq!(
                bHeteroAtomMayHaveXchgIsoH(&[atom], 0),
                Ok(1),
                "element={element} charge={charge}"
            );
        }

        for element in [9, 17, 35, 53] {
            let mut neutral = ion_atom(element, 0);
            neutral.num_H = 1;
            assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[neutral], 0), Ok(1));
            for charge in [-1, 1] {
                let charged = ion_atom(element, charge);
                assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[charged], 0), Ok(0));
            }
        }

        let proton = ion_atom(1, 1);
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[proton], 0), Ok(2));
        let mut bonded_proton = ion_atom(1, 1);
        bonded_proton.valence = 1;
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[bonded_proton], 0), Ok(0));
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[ion_atom(1, 0)], 0), Ok(0));

        let mut isotope_sum = ion_atom(8, 0);
        isotope_sum.num_H = 1;
        isotope_sum.num_iso_H = [1, 0, 0];
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[isotope_sum.clone()], 0), Ok(1));
        let mut mismatch = isotope_sum;
        mismatch.num_iso_H[1] = 1;
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[mismatch], 0), Ok(0));

        let mut charged_neighbor = vec![ion_atom(7, 1), ion_atom(6, 1)];
        charged_neighbor[0].num_H = 2;
        ion_bond(&mut charged_neighbor, 0, 1, BOND_TYPE_SINGLE as u8);
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&charged_neighbor, 0), Ok(0));

        let mut neutral_with_charged_neighbor = vec![ion_atom(7, 0), ion_atom(6, 1)];
        neutral_with_charged_neighbor[0].num_H = 2;
        ion_bond(
            &mut neutral_with_charged_neighbor,
            0,
            1,
            BOND_TYPE_SINGLE as u8,
        );
        assert_eq!(
            bHeteroAtomMayHaveXchgIsoH(&neutral_with_charged_neighbor, 0),
            Ok(1)
        );

        for (radical, expected) in [(RADICAL_SINGLET as i8, 1), (RADICAL_DOUBLET as i8, 0)] {
            let mut adjacent_radical = vec![ion_atom(7, 0), ion_atom(6, 0)];
            adjacent_radical[0].num_H = 2;
            adjacent_radical[1].radical = radical;
            ion_bond(&mut adjacent_radical, 0, 1, BOND_TYPE_SINGLE as u8);
            assert_eq!(
                bHeteroAtomMayHaveXchgIsoH(&adjacent_radical, 0),
                Ok(expected)
            );
        }

        let mut negative_valence = ion_atom(7, 0);
        negative_valence.valence = -1;
        negative_valence.chem_bonds_valence = 3;
        assert_eq!(bHeteroAtomMayHaveXchgIsoH(&[negative_valence], 0), Ok(1));

        let malformed = inp_ATOM {
            el_number: 7,
            valence: 1,
            chem_bonds_valence: 3,
            num_H: 0,
            neighbor: [u16::MAX; 20],
            ..inp_ATOM::default()
        };
        assert_eq!(
            bHeteroAtomMayHaveXchgIsoH(&[malformed], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__strutil__bnumheteratomhasisotopich__line_4133() {
        assert_eq!(bNumHeterAtomHasIsotopicH(&[], i32::MIN), Ok(0));
        assert_eq!(bNumHeterAtomHasIsotopicH(&[], 0), Ok(0));
        assert_eq!(
            bNumHeterAtomHasIsotopicH(&[], 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut isotope_only = ion_atom(6, 0);
        isotope_only.iso_atw_diff = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[isotope_only], 1), Ok(2));

        for (element, chemical_valence) in [
            (7, 2),
            (15, 2),
            (8, 1),
            (16, 1),
            (34, 1),
            (52, 1),
            (9, 0),
            (17, 0),
            (35, 0),
            (53, 0),
        ] {
            let mut atom = ion_atom(element, 0);
            atom.chem_bonds_valence = chemical_valence;
            atom.num_iso_H[0] = 1;
            assert_eq!(
                bNumHeterAtomHasIsotopicH(&[atom], 1),
                Ok(3),
                "element {element}"
            );
        }

        let mut positive_nitrogen = ion_atom(7, 1);
        positive_nitrogen.chem_bonds_valence = 3;
        positive_nitrogen.num_iso_H[1] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[positive_nitrogen], 1), Ok(3));
        let mut negative_oxygen = ion_atom(8, -1);
        negative_oxygen.num_iso_H[2] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[negative_oxygen], 1), Ok(3));

        let mut isolated_proton = ion_atom(1, 1);
        isolated_proton.iso_atw_diff = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[isolated_proton], 1), Ok(3));

        let mut charged_halogen = ion_atom(9, 1);
        charged_halogen.num_iso_H[0] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[charged_halogen], 1), Ok(2));
        let mut excessive_charge = ion_atom(8, 2);
        excessive_charge.num_iso_H[0] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[excessive_charge], 1), Ok(2));
        let mut doublet = ion_atom(8, 0);
        doublet.radical = 2;
        doublet.num_iso_H[0] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[doublet], 1), Ok(2));
        let mut valence_mismatch = ion_atom(8, 0);
        valence_mismatch.num_iso_H[0] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[valence_mismatch], 1), Ok(2));

        let mut explicit_isotope = vec![ion_atom(8, 0), ion_atom(1, 0), ion_atom(6, 0)];
        explicit_isotope[1].iso_atw_diff = 1;
        ion_bond(&mut explicit_isotope, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut explicit_isotope, 0, 2, BOND_TYPE_SINGLE as u8);
        assert_eq!(bNumHeterAtomHasIsotopicH(&explicit_isotope, 3), Ok(1));

        let mut adjacent_charge = vec![ion_atom(7, 1), ion_atom(6, 1)];
        ion_bond(&mut adjacent_charge, 0, 1, BOND_TYPE_SINGLE as u8);
        adjacent_charge[0].num_H = 2;
        adjacent_charge[0].num_iso_H[0] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&adjacent_charge, 2), Ok(2));

        let mut adjacent_radical = vec![ion_atom(7, 0), ion_atom(6, 0)];
        ion_bond(&mut adjacent_radical, 0, 1, BOND_TYPE_SINGLE as u8);
        adjacent_radical[0].num_H = 1;
        adjacent_radical[0].num_iso_H[0] = 1;
        adjacent_radical[1].radical = 2;
        assert_eq!(bNumHeterAtomHasIsotopicH(&adjacent_radical, 2), Ok(2));
        adjacent_radical[1].radical = RADICAL_SINGLET as S_CHAR;
        assert_eq!(bNumHeterAtomHasIsotopicH(&adjacent_radical, 2), Ok(3));

        let proton = ion_atom(1, 1);
        let mut after_proton = ion_atom(7, 0);
        after_proton.chem_bonds_valence = 2;
        after_proton.num_iso_H[0] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[proton, after_proton], 2), Ok(2));

        let mut negative_valence = ion_atom(7, 0);
        negative_valence.valence = -1;
        negative_valence.chem_bonds_valence = 2;
        negative_valence.num_iso_H[0] = 1;
        assert_eq!(bNumHeterAtomHasIsotopicH(&[negative_valence], 1), Ok(3));

        let mut malformed_neighbor = ion_atom(8, 0);
        malformed_neighbor.valence = 1;
        malformed_neighbor.chem_bonds_valence = 2;
        malformed_neighbor.neighbor[0] = u16::MAX;
        assert_eq!(
            bNumHeterAtomHasIsotopicH(&[malformed_neighbor], 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__strutil__the_only_doublet_neigh__line_179() {
        let mut first_ordinal = 17;
        let mut second_ordinal = 19;
        assert_eq!(
            the_only_doublet_neigh(&[], -1, &mut first_ordinal, &mut second_ordinal),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((first_ordinal, second_ordinal), (17, 19));

        let ordinary = ion_atom(6, 0);
        assert_eq!(
            the_only_doublet_neigh(&[ordinary], 0, &mut first_ordinal, &mut second_ordinal),
            Ok(-1)
        );
        assert_eq!((first_ordinal, second_ordinal), (17, 19));

        let mut isolated_doublet = ion_atom(6, 0);
        isolated_doublet.radical = RADICAL_DOUBLET as S_CHAR;
        isolated_doublet.valence = -1;
        assert_eq!(
            the_only_doublet_neigh(
                &[isolated_doublet],
                0,
                &mut first_ordinal,
                &mut second_ordinal
            ),
            Ok(-1)
        );
        assert_eq!((first_ordinal, second_ordinal), (17, 19));

        let mut pair = vec![ion_atom(6, 0), ion_atom(7, 0)];
        pair[0].radical = RADICAL_DOUBLET as S_CHAR;
        pair[1].radical = RADICAL_DOUBLET as S_CHAR;
        ion_bond(&mut pair, 0, 1, BOND_TYPE_SINGLE as u8);
        first_ordinal = 17;
        second_ordinal = 19;
        assert_eq!(
            the_only_doublet_neigh(&pair, 0, &mut first_ordinal, &mut second_ordinal),
            Ok(1)
        );
        assert_eq!((first_ordinal, second_ordinal), (0, 0));

        let mut two_first = vec![ion_atom(6, 0), ion_atom(7, 0), ion_atom(8, 0)];
        for atom in &mut two_first {
            atom.radical = RADICAL_DOUBLET as S_CHAR;
        }
        ion_bond(&mut two_first, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut two_first, 0, 2, BOND_TYPE_SINGLE as u8);
        first_ordinal = -7;
        second_ordinal = -9;
        assert_eq!(
            the_only_doublet_neigh(&two_first, 0, &mut first_ordinal, &mut second_ordinal),
            Ok(-1)
        );
        assert_eq!((first_ordinal, second_ordinal), (1, -9));

        let mut two_second = vec![ion_atom(6, 0), ion_atom(7, 0), ion_atom(8, 0)];
        for atom in &mut two_second {
            atom.radical = RADICAL_DOUBLET as S_CHAR;
        }
        ion_bond(&mut two_second, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut two_second, 1, 2, BOND_TYPE_SINGLE as u8);
        first_ordinal = -7;
        second_ordinal = -9;
        assert_eq!(
            the_only_doublet_neigh(&two_second, 0, &mut first_ordinal, &mut second_ordinal),
            Ok(-1)
        );
        assert_eq!((first_ordinal, second_ordinal), (0, 1));

        let mut malformed_first = ion_atom(6, 0);
        malformed_first.radical = RADICAL_DOUBLET as S_CHAR;
        malformed_first.valence = 1;
        malformed_first.neighbor[0] = u16::MAX;
        first_ordinal = 3;
        second_ordinal = 4;
        assert_eq!(
            the_only_doublet_neigh(
                &[malformed_first],
                0,
                &mut first_ordinal,
                &mut second_ordinal
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((first_ordinal, second_ordinal), (3, 4));

        let mut malformed_second = pair;
        malformed_second[1].valence = 2;
        malformed_second[1].neighbor[1] = u16::MAX;
        first_ordinal = 5;
        second_ordinal = 6;
        assert_eq!(
            the_only_doublet_neigh(
                &malformed_second,
                0,
                &mut first_ordinal,
                &mut second_ordinal
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((first_ordinal, second_ordinal), (0, 0));
    }

    #[test]
    fn source_port__strutil__fix_non_uniform_drawn_oxoanions__line_227() {
        fn oxo_pair(center_element: u8, terminal_element: u8) -> Vec<inp_ATOM> {
            let mut atoms = vec![ion_atom(center_element, -1), ion_atom(terminal_element, 0)];
            ion_bond(&mut atoms, 0, 1, BOND_TYPE_DOUBLE as u8);
            atoms[0].bond_stereo[0] = 7;
            atoms[1].bond_stereo[0] = -7;
            if matches!(center_element, 16 | 34 | 52) {
                atoms[0].chem_bonds_valence = 7;
            }
            atoms
        }

        let mut changes = 31;
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(i32::MIN, &mut [], &mut changes),
            Ok(0)
        );
        assert_eq!(changes, 31);
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(1, &mut [], &mut changes),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        for (center, terminal) in [
            (17, 8),
            (35, 8),
            (53, 8),
            (85, 8),
            (16, 8),
            (34, 8),
            (52, 8),
            (34, 16),
            (52, 16),
            (85, 16),
            (52, 34),
            (85, 34),
            (85, 52),
        ] {
            let mut atoms = oxo_pair(center, terminal);
            let center_chemical_valence = atoms[0].chem_bonds_valence;
            let terminal_chemical_valence = atoms[1].chem_bonds_valence;
            let mut local_changes = 0;
            assert_eq!(
                fix_non_uniform_drawn_oxoanions(1, &mut atoms, &mut local_changes),
                Ok(0),
                "center {center}, terminal {terminal}"
            );
            assert_eq!(local_changes, 1);
            assert_eq!((atoms[0].charge, atoms[1].charge), (0, -1));
            assert_eq!(
                (atoms[0].bond_type[0], atoms[1].bond_type[0]),
                (BOND_TYPE_SINGLE as u8, BOND_TYPE_SINGLE as u8)
            );
            assert_eq!((atoms[0].bond_stereo[0], atoms[1].bond_stereo[0]), (0, 0));
            assert_eq!(
                (atoms[0].chem_bonds_valence, atoms[1].chem_bonds_valence),
                (
                    center_chemical_valence.wrapping_sub(1),
                    terminal_chemical_valence.wrapping_sub(1)
                )
            );
        }

        for mut rejected in {
            let mut fixtures = Vec::new();
            let mut wrong_charge = oxo_pair(17, 8);
            wrong_charge[0].charge = 0;
            fixtures.push(wrong_charge);
            fixtures.push(oxo_pair(6, 8));
            let mut low_chalcogen = oxo_pair(16, 8);
            low_chalcogen[0].chem_bonds_valence = 6;
            fixtures.push(low_chalcogen);
            let mut center_radical = oxo_pair(17, 8);
            center_radical[0].radical = RADICAL_DOUBLET as S_CHAR;
            fixtures.push(center_radical);
            let mut nonterminal = oxo_pair(17, 8);
            nonterminal[1].valence = 2;
            fixtures.push(nonterminal);
            let mut single_bond = oxo_pair(17, 8);
            single_bond[0].bond_type[0] = BOND_TYPE_SINGLE as u8;
            fixtures.push(single_bond);
            let mut charged_terminal = oxo_pair(17, 8);
            charged_terminal[1].charge = 1;
            fixtures.push(charged_terminal);
            let mut radical_terminal = oxo_pair(17, 8);
            radical_terminal[1].radical = RADICAL_DOUBLET as S_CHAR;
            fixtures.push(radical_terminal);
            fixtures.push(oxo_pair(17, 16));
            fixtures.push(oxo_pair(16, 16));
            fixtures.push(oxo_pair(34, 34));
            fixtures.push(oxo_pair(52, 52));
            fixtures.push(oxo_pair(85, 7));
            fixtures
        } {
            let before = rejected.clone();
            let mut local_changes = -7;
            assert_eq!(
                fix_non_uniform_drawn_oxoanions(1, &mut rejected, &mut local_changes),
                Ok(0)
            );
            assert_eq!(rejected, before);
            assert_eq!(local_changes, -7);
        }

        let mut singlet = oxo_pair(17, 8);
        singlet[0].radical = RADICAL_SINGLET as S_CHAR;
        singlet[1].radical = RADICAL_SINGLET as S_CHAR;
        changes = 0;
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(1, &mut singlet, &mut changes),
            Ok(0)
        );
        assert_eq!(changes, 1);

        let mut selection = vec![ion_atom(85, -1), ion_atom(16, 0), ion_atom(8, 0)];
        ion_bond(&mut selection, 0, 1, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut selection, 0, 2, BOND_TYPE_DOUBLE as u8);
        changes = 0;
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(1, &mut selection, &mut changes),
            Ok(0)
        );
        assert_eq!((selection[1].charge, selection[2].charge), (0, -1));
        assert_eq!(
            (selection[0].bond_type[0], selection[0].bond_type[1]),
            (BOND_TYPE_DOUBLE as u8, BOND_TYPE_SINGLE as u8)
        );

        let mut isotope_selection = vec![ion_atom(17, -1), ion_atom(8, 0), ion_atom(8, 0)];
        isotope_selection[2].iso_atw_diff = 1;
        ion_bond(&mut isotope_selection, 0, 1, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut isotope_selection, 0, 2, BOND_TYPE_DOUBLE as u8);
        changes = 0;
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(1, &mut isotope_selection, &mut changes),
            Ok(0)
        );
        assert_eq!(
            (isotope_selection[1].charge, isotope_selection[2].charge),
            (0, -1)
        );

        let mut isotope_tie = vec![ion_atom(17, -1), ion_atom(8, 0), ion_atom(8, 0)];
        isotope_tie[1].iso_atw_diff = 1;
        isotope_tie[2].iso_atw_diff = 2;
        ion_bond(&mut isotope_tie, 0, 1, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut isotope_tie, 0, 2, BOND_TYPE_DOUBLE as u8);
        changes = 0;
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(1, &mut isotope_tie, &mut changes),
            Ok(0)
        );
        assert_eq!((isotope_tie[1].charge, isotope_tie[2].charge), (-1, 0));

        let mut wrapping = oxo_pair(17, 8);
        wrapping[0].chem_bonds_valence = i8::MIN;
        wrapping[1].chem_bonds_valence = i8::MIN;
        changes = i32::MAX;
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(1, &mut wrapping, &mut changes),
            Ok(0)
        );
        assert_eq!(
            (
                wrapping[0].chem_bonds_valence,
                wrapping[1].chem_bonds_valence,
                changes
            ),
            (i8::MAX, i8::MAX, i32::MIN)
        );

        let mut malformed = oxo_pair(17, 8);
        malformed[0].valence = 2;
        malformed[0].neighbor[1] = u16::MAX;
        let malformed_before = malformed.clone();
        changes = 4;
        assert_eq!(
            fix_non_uniform_drawn_oxoanions(1, &mut malformed, &mut changes),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(malformed, malformed_before);
        assert_eq!(changes, 4);
    }

    #[test]
    fn source_port__strutil__fix_non_uniform_drawn_amidiniums__line_452() {
        fn amidinium(center_element: u8) -> Vec<inp_ATOM> {
            let mut atoms = vec![
                ion_atom(center_element, 1),
                ion_atom(7, 0),
                ion_atom(7, 0),
                ion_atom(6, 0),
            ];
            ion_bond(&mut atoms, 0, 1, BOND_TYPE_SINGLE as u8);
            ion_bond(&mut atoms, 0, 2, BOND_TYPE_SINGLE as u8);
            ion_bond(&mut atoms, 0, 3, BOND_TYPE_SINGLE as u8);
            atoms[1].num_H = 1;
            atoms[0].bond_stereo[0] = 5;
            atoms[1].bond_stereo[0] = -5;
            atoms
        }

        let mut changes = 41;
        assert_eq!(
            fix_non_uniform_drawn_amidiniums(i32::MIN, &mut [], &mut changes),
            Ok(0)
        );
        assert_eq!(changes, 41);
        assert_eq!(
            fix_non_uniform_drawn_amidiniums(1, &mut [], &mut changes),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        for center_element in [6_u8, 16, 15] {
            let mut atoms = amidinium(center_element);
            let mut local_changes = 0;
            assert_eq!(
                fix_non_uniform_drawn_amidiniums(1, &mut atoms, &mut local_changes),
                Ok(0)
            );
            assert_eq!(local_changes, 1);
            assert_eq!((atoms[0].charge, atoms[1].charge), (0, 1));
            assert_eq!(
                (atoms[0].bond_type[0], atoms[1].bond_type[0]),
                (BOND_TYPE_DOUBLE as u8, BOND_TYPE_DOUBLE as u8)
            );
            assert_eq!(
                (atoms[0].chem_bonds_valence, atoms[1].chem_bonds_valence),
                (4, 2)
            );
            assert_eq!((atoms[0].bond_stereo[0], atoms[1].bond_stereo[0]), (5, -5));
        }

        for mut rejected in {
            let mut fixtures = Vec::new();
            let mut wrong_charge = amidinium(6);
            wrong_charge[0].charge = 0;
            fixtures.push(wrong_charge);
            fixtures.push(amidinium(8));
            let mut wrong_valence = amidinium(6);
            wrong_valence[0].valence = 2;
            fixtures.push(wrong_valence);
            let mut wrong_chemical_valence = amidinium(6);
            wrong_chemical_valence[0].chem_bonds_valence = 4;
            fixtures.push(wrong_chemical_valence);
            let mut center_radical = amidinium(6);
            center_radical[0].radical = RADICAL_DOUBLET as S_CHAR;
            fixtures.push(center_radical);
            let mut charged_neighbor = amidinium(6);
            charged_neighbor[3].charge = -1;
            fixtures.push(charged_neighbor);
            let mut high_nitrogen_valence = amidinium(6);
            high_nitrogen_valence[1].valence = 4;
            fixtures.push(high_nitrogen_valence);
            let mut high_nitrogen_chemical_valence = amidinium(6);
            high_nitrogen_chemical_valence[1].chem_bonds_valence = 4;
            fixtures.push(high_nitrogen_chemical_valence);
            let mut one_nitrogen = amidinium(6);
            one_nitrogen[2].el_number = 8;
            fixtures.push(one_nitrogen);
            let mut three_nitrogens = amidinium(6);
            three_nitrogens[3].el_number = 7;
            fixtures.push(three_nitrogens);
            let mut no_hydrogen = amidinium(6);
            no_hydrogen[1].num_H = 0;
            fixtures.push(no_hydrogen);
            fixtures
        } {
            let before = rejected.clone();
            let mut local_changes = -3;
            assert_eq!(
                fix_non_uniform_drawn_amidiniums(1, &mut rejected, &mut local_changes),
                Ok(0)
            );
            assert_eq!(rejected, before);
            assert_eq!(local_changes, -3);
        }

        let mut singlet_and_isotope = amidinium(6);
        singlet_and_isotope[0].radical = RADICAL_SINGLET as S_CHAR;
        singlet_and_isotope[1].num_H = 0;
        singlet_and_isotope[2].num_iso_H[1] = 1;
        changes = 0;
        assert_eq!(
            fix_non_uniform_drawn_amidiniums(1, &mut singlet_and_isotope, &mut changes),
            Ok(0)
        );
        assert_eq!(changes, 1);
        assert_eq!(
            (singlet_and_isotope[0].charge, singlet_and_isotope[1].charge),
            (0, 1)
        );

        let mut reciprocal_second = amidinium(6);
        reciprocal_second.extend([ion_atom(6, 0), ion_atom(6, 0)]);
        reciprocal_second[1].valence = 3;
        reciprocal_second[1].chem_bonds_valence = 3;
        reciprocal_second[1].neighbor[..3].copy_from_slice(&[4, 5, 0]);
        reciprocal_second[1].bond_type[..3].copy_from_slice(&[1, 1, 1]);
        changes = 0;
        assert_eq!(
            fix_non_uniform_drawn_amidiniums(1, &mut reciprocal_second, &mut changes),
            Ok(0)
        );
        assert_eq!(reciprocal_second[1].bond_type, {
            let mut expected = [0_u8; MAXVAL as usize];
            expected[..3].copy_from_slice(&[1, 1, BOND_TYPE_DOUBLE as u8]);
            expected
        });

        let mut missing_reciprocal = amidinium(6);
        missing_reciprocal.push(ion_atom(6, 0));
        missing_reciprocal[1].neighbor[0] = 4;
        changes = 0;
        assert_eq!(
            fix_non_uniform_drawn_amidiniums(1, &mut missing_reciprocal, &mut changes),
            Ok(0)
        );
        assert_eq!(
            (
                missing_reciprocal[1].bond_type[0],
                missing_reciprocal[1].bond_type[1]
            ),
            (BOND_TYPE_SINGLE as u8, BOND_TYPE_DOUBLE as u8)
        );

        let mut unchecked_initial_bonds = amidinium(6);
        unchecked_initial_bonds[0].bond_type = [BOND_TYPE_TRIPLE as u8; MAXVAL as usize];
        changes = i32::MAX;
        assert_eq!(
            fix_non_uniform_drawn_amidiniums(1, &mut unchecked_initial_bonds, &mut changes),
            Ok(0)
        );
        assert_eq!(changes, i32::MIN);
        assert_eq!(
            unchecked_initial_bonds[0].bond_type[..3],
            [
                BOND_TYPE_DOUBLE as u8,
                BOND_TYPE_TRIPLE as u8,
                BOND_TYPE_TRIPLE as u8
            ]
        );

        let mut malformed = amidinium(6);
        malformed[0].neighbor[2] = u16::MAX;
        let before = malformed.clone();
        changes = 9;
        assert_eq!(
            fix_non_uniform_drawn_amidiniums(1, &mut malformed, &mut changes),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(malformed, before);
        assert_eq!(changes, 9);
    }

    #[test]
    fn source_port__strutil__fix_odd_things__line_603() {
        assert_eq!(fix_odd_things(i32::MIN, &mut [], 0, 0), Ok(0));
        assert_eq!(
            fix_odd_things(1, &mut [], 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut oxo_disabled = {
            let mut atoms = vec![ion_atom(17, -1), ion_atom(8, 0)];
            ion_bond(&mut atoms, 0, 1, BOND_TYPE_DOUBLE as u8);
            atoms
        };
        let oxo_before = oxo_disabled.clone();
        assert_eq!(fix_odd_things(1, &mut oxo_disabled, 0, 0), Ok(0));
        assert_eq!(oxo_disabled, oxo_before);
        assert_eq!(fix_odd_things(1, &mut oxo_disabled, 0, 1), Ok(1));
        assert_eq!((oxo_disabled[0].charge, oxo_disabled[1].charge), (0, -1));

        for source_charge in [-1_i8, 1] {
            let mut charged_hydrogen = vec![ion_atom(1, source_charge), ion_atom(6, 0)];
            ion_bond(&mut charged_hydrogen, 0, 1, BOND_TYPE_SINGLE as u8);
            assert_eq!(fix_odd_things(2, &mut charged_hydrogen, 0, 0), Ok(0));
            assert_eq!(
                (charged_hydrogen[0].charge, charged_hydrogen[1].charge),
                (0, source_charge)
            );
        }
        let mut charge_wrap = vec![ion_atom(1, 1), ion_atom(6, i8::MAX)];
        ion_bond(&mut charge_wrap, 0, 1, BOND_TYPE_SINGLE as u8);
        assert_eq!(fix_odd_things(2, &mut charge_wrap, 0, 0), Ok(0));
        assert_eq!((charge_wrap[0].charge, charge_wrap[1].charge), (0, i8::MIN));

        for mut rejected_hydrogen in {
            let mut fixtures = Vec::new();
            let mut wrong_valence = vec![ion_atom(1, 1), ion_atom(6, 0)];
            ion_bond(&mut wrong_valence, 0, 1, BOND_TYPE_SINGLE as u8);
            wrong_valence[0].valence = 2;
            fixtures.push(wrong_valence);
            let mut neutral = vec![ion_atom(1, 0), ion_atom(6, 0)];
            ion_bond(&mut neutral, 0, 1, BOND_TYPE_SINGLE as u8);
            fixtures.push(neutral);
            let mut radical = vec![ion_atom(1, 1), ion_atom(6, 0)];
            ion_bond(&mut radical, 0, 1, BOND_TYPE_SINGLE as u8);
            radical[0].radical = RADICAL_DOUBLET as S_CHAR;
            fixtures.push(radical);
            let mut double_bond = vec![ion_atom(1, 1), ion_atom(6, 0)];
            ion_bond(&mut double_bond, 0, 1, BOND_TYPE_DOUBLE as u8);
            fixtures.push(double_bond);
            let mut hydrogen_pair = vec![ion_atom(1, 1), ion_atom(1, 0)];
            ion_bond(&mut hydrogen_pair, 0, 1, BOND_TYPE_SINGLE as u8);
            fixtures.push(hydrogen_pair);
            let mut source_h = vec![ion_atom(1, 1), ion_atom(6, 0)];
            ion_bond(&mut source_h, 0, 1, BOND_TYPE_SINGLE as u8);
            source_h[0].num_H = 1;
            fixtures.push(source_h);
            let mut target_h = vec![ion_atom(1, 1), ion_atom(6, 0)];
            ion_bond(&mut target_h, 0, 1, BOND_TYPE_SINGLE as u8);
            target_h[1].num_iso_H[0] = 1;
            fixtures.push(target_h);
            fixtures
        } {
            let before = rejected_hydrogen.clone();
            assert_eq!(
                fix_odd_things(rejected_hydrogen.len() as i32, &mut rejected_hydrogen, 0, 0),
                Ok(0)
            );
            assert_eq!(rejected_hydrogen, before);
        }

        let mut implicit = vec![ion_atom(7, 1), ion_atom(6, 0), ion_atom(7, -1)];
        ion_bond(&mut implicit, 0, 1, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut implicit, 1, 2, BOND_TYPE_SINGLE as u8);
        implicit[0].num_H = 2;
        implicit[2].num_H = 1;
        implicit[0].bond_stereo[0] = 9;
        implicit[2].bond_stereo[0] = -9;
        assert_eq!(fix_odd_things(3, &mut implicit, 0, 0), Ok(1));
        assert_eq!((implicit[0].charge, implicit[2].charge), (0, 0));
        assert_eq!(
            (implicit[0].bond_type[0], implicit[2].bond_type[0]),
            (BOND_TYPE_SINGLE as u8, BOND_TYPE_DOUBLE as u8)
        );
        assert_eq!(
            (
                implicit[0].chem_bonds_valence,
                implicit[2].chem_bonds_valence
            ),
            (1, 2)
        );
        assert_eq!(
            (implicit[0].bond_stereo[0], implicit[2].bond_stereo[0]),
            (9, -9)
        );

        let mut explicit = vec![
            ion_atom(7, 1),
            ion_atom(6, 0),
            ion_atom(7, -1),
            ion_atom(1, 0),
            ion_atom(1, 0),
        ];
        ion_bond(&mut explicit, 0, 1, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut explicit, 1, 2, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut explicit, 0, 3, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut explicit, 2, 4, BOND_TYPE_SINGLE as u8);
        explicit[0].num_H = 1;
        assert_eq!(fix_odd_things(5, &mut explicit, 0, 0), Ok(1));
        assert_eq!((explicit[0].charge, explicit[2].charge), (0, 0));
        assert_eq!(
            (explicit[0].bond_type[0], explicit[2].bond_type[0]),
            (BOND_TYPE_SINGLE as u8, BOND_TYPE_DOUBLE as u8)
        );

        fn divalent_center_fixture() -> Vec<inp_ATOM> {
            let mut atoms = vec![
                ion_atom(16, 2),
                ion_atom(8, -1),
                ion_atom(16, -1),
                ion_atom(6, 0),
                ion_atom(6, 0),
            ];
            for neighbor in 1..5 {
                ion_bond(&mut atoms, 0, neighbor, BOND_TYPE_SINGLE as u8);
            }
            atoms[0].bond_stereo[..2].copy_from_slice(&[7, -7]);
            atoms[1].bond_stereo[0] = 5;
            atoms[2].bond_stereo[0] = -5;
            atoms
        }
        let mut no_bug_fix = divalent_center_fixture();
        assert_eq!(fix_odd_things(5, &mut no_bug_fix, 0, 0), Ok(2));
        assert_eq!(
            (
                no_bug_fix[0].charge,
                no_bug_fix[1].charge,
                no_bug_fix[2].charge
            ),
            (2, 0, 0)
        );
        assert_eq!(
            (
                no_bug_fix[0].chem_bonds_valence,
                no_bug_fix[1].chem_bonds_valence,
                no_bug_fix[2].chem_bonds_valence
            ),
            (6, 2, 2)
        );
        assert_eq!(
            (
                no_bug_fix[0].bond_stereo[0],
                no_bug_fix[0].bond_stereo[1],
                no_bug_fix[1].bond_stereo[0],
                no_bug_fix[2].bond_stereo[0]
            ),
            (0, 0, 0, 0)
        );

        let mut with_bug_fix = divalent_center_fixture();
        assert_eq!(fix_odd_things(5, &mut with_bug_fix, 1, 0), Ok(2));
        assert_eq!(with_bug_fix[0].charge, 0);

        let mut phosphorus_inactive = divalent_center_fixture();
        phosphorus_inactive[0].el_number = 15;
        phosphorus_inactive[0].charge = 1;
        phosphorus_inactive[2].charge = 0;
        assert_eq!(fix_odd_things(5, &mut phosphorus_inactive, 1, 0), Ok(1));
        assert_eq!(
            (
                phosphorus_inactive[0].charge,
                phosphorus_inactive[1].charge,
                phosphorus_inactive[0].bond_type[0],
                phosphorus_inactive[1].bond_type[0],
                phosphorus_inactive[0].chem_bonds_valence,
                phosphorus_inactive[1].chem_bonds_valence,
            ),
            (0, 0, BOND_TYPE_DOUBLE as u8, BOND_TYPE_DOUBLE as u8, 5, 2)
        );

        for initial_bond in [BOND_TYPE_SINGLE as u8, BOND_TYPE_DOUBLE as u8] {
            let mut doublets = vec![ion_atom(6, 0), ion_atom(7, 0)];
            doublets[0].radical = RADICAL_DOUBLET as S_CHAR;
            doublets[1].radical = RADICAL_DOUBLET as S_CHAR;
            ion_bond(&mut doublets, 0, 1, initial_bond);
            assert_eq!(fix_odd_things(2, &mut doublets, 0, 0), Ok(0));
            assert_eq!(doublets[0].bond_type[0], initial_bond + 1);
            assert_eq!((doublets[0].radical, doublets[1].radical), (0, 0));
        }
        let mut triple_doublets = vec![ion_atom(6, 0), ion_atom(7, 0)];
        triple_doublets[0].radical = RADICAL_DOUBLET as S_CHAR;
        triple_doublets[1].radical = RADICAL_DOUBLET as S_CHAR;
        ion_bond(&mut triple_doublets, 0, 1, BOND_TYPE_TRIPLE as u8);
        let triple_before = triple_doublets.clone();
        assert_eq!(fix_odd_things(2, &mut triple_doublets, 0, 0), Ok(0));
        assert_eq!(triple_doublets, triple_before);

        let mut early_ion_pair = vec![
            ion_atom(7, 1),
            ion_atom(8, 0),
            ion_atom(8, -1),
            ion_atom(6, 0),
        ];
        early_ion_pair[0].radical = RADICAL_DOUBLET as S_CHAR;
        early_ion_pair[2].radical = 3;
        ion_bond(&mut early_ion_pair, 0, 1, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut early_ion_pair, 0, 2, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut early_ion_pair, 0, 3, BOND_TYPE_SINGLE as u8);
        assert_eq!(fix_odd_things(4, &mut early_ion_pair, 0, 0), Ok(1));
        assert_eq!((early_ion_pair[0].charge, early_ion_pair[2].charge), (0, 0));

        let mut partial_error = vec![ion_atom(1, 1), ion_atom(6, 0), ion_atom(1, 1)];
        ion_bond(&mut partial_error, 0, 1, BOND_TYPE_SINGLE as u8);
        partial_error[2].valence = 1;
        partial_error[2].bond_type[0] = BOND_TYPE_SINGLE as u8;
        partial_error[2].neighbor[0] = u16::MAX;
        assert_eq!(
            fix_odd_things(3, &mut partial_error, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((partial_error[0].charge, partial_error[1].charge), (0, 1));
        assert_eq!(partial_error[2].charge, 1);
    }

    #[test]
    fn source_port__strutil__post_fix_odd_things__line_922() {
        assert_eq!(post_fix_odd_things(i32::MIN, &mut []), 0);
        assert_eq!(post_fix_odd_things(i32::MAX, &mut []), 0);

        let mut atom = ion_atom(118, i8::MIN);
        atom.radical = i8::MAX;
        atom.valence = i8::MIN;
        atom.chem_bonds_valence = i8::MAX;
        atom.num_H = i8::MIN;
        atom.num_iso_H = [i8::MIN, 0, i8::MAX];
        let before = atom.clone();
        assert_eq!(
            post_fix_odd_things(i32::MIN, std::slice::from_mut(&mut atom)),
            0
        );
        assert_eq!(atom, before);
    }

    #[test]
    fn source_port__strutil__disconnectoneligand__line_3112() {
        assert_eq!(
            DisconnectOneLigand(&mut [], None, &[], &[0], 0, 0, -1, 0, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut no_ligand = vec![ion_atom(11, 0)];
        assert_eq!(
            DisconnectOneLigand(&mut no_ligand, None, &[1], &[0], 0, 1, 0, 0, None),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut halogen = vec![ion_atom(11, 0), ion_atom(9, 0)];
        ion_bond(&mut halogen, 0, 1, BOND_TYPE_SINGLE as u8);
        halogen[0].component = 1;
        halogen[1].component = 2;
        let mut old_components = [7_u16, 8];
        let mut flags = 0x20_u64;
        assert_eq!(
            DisconnectOneLigand(
                &mut halogen,
                Some(&mut old_components),
                &[1, 0],
                &[9, 8, 7, 0],
                1,
                2,
                0,
                0,
                Some(&mut flags),
            ),
            Ok(1)
        );
        assert_eq!(old_components, [0, 0]);
        assert_eq!(flags, 0x20 | TG_FLAG_MOVE_CHARGE_COORD_DONE as u64);
        assert_eq!(
            (
                halogen[0].charge,
                halogen[0].valence,
                halogen[0].chem_bonds_valence,
                halogen[1].charge,
                halogen[1].valence,
                halogen[1].chem_bonds_valence,
            ),
            (1, 0, 0, -1, 0, 0)
        );

        let mut unmarked = vec![ion_atom(11, 0), ion_atom(9, 0)];
        ion_bond(&mut unmarked, 0, 1, BOND_TYPE_SINGLE as u8);
        let unmarked_before = unmarked.clone();
        assert_eq!(
            DisconnectOneLigand(&mut unmarked, None, &[0, 0], &[9, 0], 1, 2, 0, 0, None,),
            Ok(0)
        );
        assert_eq!(unmarked, unmarked_before);

        let mut two_metals = vec![ion_atom(11, 0), ion_atom(9, 0), ion_atom(12, 0)];
        ion_bond(&mut two_metals, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut two_metals, 1, 2, BOND_TYPE_SINGLE as u8);
        assert_eq!(
            DisconnectOneLigand(&mut two_metals, None, &[1, 0, 1], &[9, 0], 1, 3, 0, 0, None,),
            Ok(2)
        );
        assert_eq!(
            two_metals
                .iter()
                .map(|atom| (atom.valence, atom.chem_bonds_valence, atom.charge))
                .collect::<Vec<_>>(),
            vec![(0, 0, 0), (0, 0, 0), (0, 0, 0)]
        );

        let mut outside_count = vec![ion_atom(11, 0), ion_atom(9, 0), ion_atom(12, 0)];
        ion_bond(&mut outside_count, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut outside_count, 1, 2, BOND_TYPE_SINGLE as u8);
        assert_eq!(
            DisconnectOneLigand(
                &mut outside_count,
                None,
                &[1, 0, 1],
                &[9, 0],
                1,
                2,
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!((outside_count[0].valence, outside_count[2].valence), (0, 1));

        for remaining_aromatic_bonds in [0_usize, 1, 2, 3, 4] {
            let mut aromatic = vec![ion_atom(11, 0), ion_atom(7, 0)];
            aromatic.extend((0..remaining_aromatic_bonds).map(|_| ion_atom(6, 0)));
            ion_bond(&mut aromatic, 0, 1, BOND_TYPE_TRIPLE as u8 + 1);
            for neighbor in 2..aromatic.len() {
                ion_bond(&mut aromatic, 1, neighbor, BOND_TYPE_TRIPLE as u8 + 1);
            }
            let mut marks = vec![0_i8; aromatic.len()];
            marks[0] = 1;
            let result = DisconnectOneLigand(
                &mut aromatic,
                None,
                &marks,
                &[7],
                0,
                marks.len() as i32,
                0,
                0,
                None,
            );
            if matches!(remaining_aromatic_bonds, 0 | 2 | 3) {
                assert_eq!(result, Err(SourceHeapError::MissingNulTerminator));
            } else {
                assert_eq!(result, Ok(1));
            }
            assert_eq!(aromatic[0].chem_bonds_valence, 3);
            assert_eq!(aromatic[0].valence, 0);
        }

        let mut invalid_radical = vec![ion_atom(11, 0), ion_atom(8, 0)];
        ion_bond(&mut invalid_radical, 0, 1, BOND_TYPE_SINGLE as u8);
        invalid_radical[1].radical = 2;
        assert_eq!(
            DisconnectOneLigand(&mut invalid_radical, None, &[1, 0], &[8], 0, 2, 0, 0, None,),
            Ok(1)
        );
        assert_eq!(
            (invalid_radical[0].charge, invalid_radical[1].charge),
            (0, 0)
        );

        let mut absent_element = vec![ion_atom(11, 0), ion_atom(7, 0)];
        ion_bond(&mut absent_element, 0, 1, BOND_TYPE_SINGLE as u8);
        assert_eq!(
            DisconnectOneLigand(
                &mut absent_element,
                None,
                &[1, 0],
                &[8, 0, 7],
                0,
                2,
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!((absent_element[0].charge, absent_element[1].charge), (0, 0));

        let mut missing_nul = vec![ion_atom(11, 0), ion_atom(9, 0)];
        ion_bond(&mut missing_nul, 0, 1, BOND_TYPE_SINGLE as u8);
        assert_eq!(
            DisconnectOneLigand(&mut missing_nul, None, &[1, 0], &[9], 1, 2, 0, 0, None,),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!((missing_nul[0].valence, missing_nul[1].valence), (0, 0));

        let mut nonhalogen_zero = vec![ion_atom(11, 0), ion_atom(8, 0)];
        ion_bond(&mut nonhalogen_zero, 0, 1, BOND_TYPE_SINGLE as u8);
        assert_eq!(
            DisconnectOneLigand(
                &mut nonhalogen_zero,
                None,
                &[1, 0],
                &[8, 0],
                0,
                2,
                0,
                0,
                None,
            ),
            Ok(1)
        );
        assert_eq!(
            (nonhalogen_zero[0].charge, nonhalogen_zero[1].charge),
            (0, 0)
        );

        let mut oxygen = vec![ion_atom(11, 0), ion_atom(8, 0), ion_atom(6, 0)];
        ion_bond(&mut oxygen, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut oxygen, 1, 2, BOND_TYPE_SINGLE as u8);
        let mut oxygen_flags = 0_u64;
        assert_eq!(
            DisconnectOneLigand(
                &mut oxygen,
                None,
                &[1, 0, 0],
                &[8, 0],
                0,
                3,
                0,
                0,
                Some(&mut oxygen_flags),
            ),
            Ok(1)
        );
        assert_eq!((oxygen[0].charge, oxygen[1].charge), (1, -1));
        assert_eq!(oxygen_flags, TG_FLAG_MOVE_CHARGE_COORD_DONE as u64);

        let mut isotope_h = vec![ion_atom(11, 0), ion_atom(9, -1)];
        ion_bond(&mut isotope_h, 0, 1, BOND_TYPE_SINGLE as u8);
        isotope_h[1].num_iso_H = [0, 1, 0];
        assert_eq!(
            DisconnectOneLigand(&mut isotope_h, None, &[1, 0], &[9, 0], 1, 2, 0, 0, None,),
            Ok(1)
        );
        assert_eq!((isotope_h[0].charge, isotope_h[1].charge), (-1, 0));

        let mut preserved_imine = vec![
            ion_atom(11, 0),
            ion_atom(7, 0),
            ion_atom(6, 0),
            ion_atom(6, 0),
        ];
        ion_bond(&mut preserved_imine, 0, 1, BOND_TYPE_SINGLE as u8);
        ion_bond(&mut preserved_imine, 1, 2, BOND_TYPE_DOUBLE as u8);
        ion_bond(&mut preserved_imine, 1, 3, BOND_TYPE_DOUBLE as u8);
        let mut preserved_flags = 0x80_u64;
        assert_eq!(
            DisconnectOneLigand(
                &mut preserved_imine,
                None,
                &[1, 0, 0, 0],
                &[7, 0],
                0,
                4,
                0,
                0,
                Some(&mut preserved_flags),
            ),
            Ok(1)
        );
        assert_eq!(
            (preserved_imine[0].charge, preserved_imine[1].charge),
            (0, 0)
        );
        assert_eq!(preserved_flags, 0x80);
        assert_eq!(preserved_imine[1].bond_type[..2], [2, 2]);

        let mut wrapped = vec![ion_atom(11, i8::MIN), ion_atom(9, i8::MIN)];
        ion_bond(&mut wrapped, 0, 1, BOND_TYPE_SINGLE as u8);
        assert_eq!(
            DisconnectOneLigand(&mut wrapped, None, &[1, 0], &[9, 0], 1, 2, 0, 0, None,),
            Ok(1)
        );
        assert_eq!((wrapped[0].charge, wrapped[1].charge), (1, -1));

        let mut partial = vec![ion_atom(11, 0), ion_atom(9, 0), ion_atom(12, 0)];
        partial[0].neighbor[0] = 1;
        partial[0].bond_type[0] = BOND_TYPE_SINGLE as u8;
        partial[0].valence = 1;
        partial[0].chem_bonds_valence = 1;
        partial[1].neighbor[0] = u16::MAX;
        partial[1].bond_type[0] = BOND_TYPE_SINGLE as u8;
        partial[1].neighbor[1] = 2;
        partial[1].bond_type[1] = BOND_TYPE_SINGLE as u8;
        partial[1].valence = 2;
        partial[1].chem_bonds_valence = 2;
        partial[2].neighbor[0] = 1;
        partial[2].bond_type[0] = BOND_TYPE_SINGLE as u8;
        partial[2].valence = 1;
        partial[2].chem_bonds_valence = 1;
        let mut large_marks = vec![0_i8; u16::MAX as usize + 1];
        large_marks[u16::MAX as usize] = 1;
        large_marks[2] = 1;
        assert_eq!(
            DisconnectOneLigand(
                &mut partial,
                None,
                &large_marks,
                &[9, 0],
                1,
                u16::MAX as i32 + 1,
                0,
                0,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((partial[1].valence, partial[2].valence), (1, 0));
    }

    fn source_remove_terminal_hdt(
        atoms: Vec<inp_ATOM>,
        fix_charge: i32,
    ) -> (i32, Vec<inp_ATOM>, SourceHeap) {
        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(atoms).unwrap();
        let count = heap.slice(pointer.as_const()).unwrap().len() as i32;
        let result = remove_terminal_HDT(&mut heap, count, pointer, fix_charge).unwrap();
        let output = heap.slice(pointer.as_const()).unwrap().to_vec();
        (result, output, heap)
    }

    #[test]
    fn source_port__strutil__remove_terminal_hdt__line_3723() {
        let low_isotope = inp_ATOM {
            iso_atw_diff: i8::MIN,
            component: u16::MAX,
            ..inp_ATOM::default()
        };
        let high_isotope = inp_ATOM {
            iso_atw_diff: i8::MAX,
            component: 0,
            ..inp_ATOM::default()
        };
        assert_eq!(
            cmp_iso_atw_diff_component_no(&low_isotope, &high_isotope),
            -255
        );
        assert_eq!(
            cmp_iso_atw_diff_component_no(&high_isotope, &low_isotope),
            255
        );
        let low_component = inp_ATOM {
            iso_atw_diff: 7,
            component: 0,
            ..inp_ATOM::default()
        };
        let high_component = inp_ATOM {
            iso_atw_diff: 7,
            component: u16::MAX,
            ..inp_ATOM::default()
        };
        assert_eq!(
            cmp_iso_atw_diff_component_no(&low_component, &high_component),
            -65_535
        );

        let original = vec![inp_ATOM {
            elname: [b'C' as i8, 0, 0, 0, 0, 0],
            component: 91,
            ..inp_ATOM::default()
        }];
        for successful_allocations in [0, 1] {
            let mut heap = SourceHeap::default();
            let pointer = heap.allocate_model_storage(original.clone()).unwrap();
            let baseline = heap.live_allocation_count();
            heap.fail_after_allocations(successful_allocations);
            assert_eq!(remove_terminal_HDT(&mut heap, 1, pointer, 1), Ok(-1));
            assert_eq!(heap.source_allocation_calls(), 2);
            assert_eq!(heap.live_allocation_count(), baseline);
            assert_eq!(heap.slice(pointer.as_const()).unwrap(), original);
        }

        let mut carbon = inp_ATOM {
            elname: [b'C' as i8, 0, 0, 0, 0, 0],
            valence: 1,
            chem_bonds_valence: 1,
            component: 77,
            ..inp_ATOM::default()
        };
        carbon.neighbor[0] = 1;
        let mut oxygen = inp_ATOM {
            elname: [b'O' as i8, 0, 0, 0, 0, 0],
            valence: 1,
            chem_bonds_valence: 1,
            component: 88,
            ..inp_ATOM::default()
        };
        oxygen.neighbor[0] = 0;
        let (result, unchanged, heap) = source_remove_terminal_hdt(vec![carbon, oxygen], 1);
        assert_eq!(result, 2);
        assert_eq!(heap.live_allocation_count(), 1);
        assert_eq!((unchanged[0].component, unchanged[1].component), (0, 1));
        assert_eq!((unchanged[0].neighbor[0], unchanged[1].neighbor[0]), (1, 0));

        let mut center = inp_ATOM {
            elname: [b'C' as i8, 0, 0, 0, 0, 0],
            valence: 4,
            chem_bonds_valence: 4,
            ..inp_ATOM::default()
        };
        center.neighbor[..4].copy_from_slice(&[1, 2, 3, 4]);
        center.bond_type[..4].fill(BOND_TYPE_SINGLE as u8);
        let hydrogen = terminal_hdt_atom(b'H', 0);
        let deuterium = terminal_hdt_atom(b'D', 0);
        let tritium = terminal_hdt_atom(b'T', 0);
        let mut other = inp_ATOM {
            elname: [b'C' as i8, 0, 0, 0, 0, 0],
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        other.neighbor[0] = 0;
        let (result, reordered, _) =
            source_remove_terminal_hdt(vec![center, hydrogen, deuterium, tritium, other], 1);
        assert_eq!(result, 2);
        assert_eq!(reordered[0].neighbor[0], 1);
        assert_eq!(
            (reordered[0].valence, reordered[0].chem_bonds_valence),
            (1, 1)
        );
        assert_eq!((reordered[0].num_H, reordered[0].num_iso_H), (1, [0, 1, 1]));
        assert_eq!(
            reordered[2..]
                .iter()
                .map(|atom| (atom.iso_atw_diff, atom.component, atom.neighbor[0]))
                .collect::<Vec<_>>(),
            vec![(0, 1, 0), (2, 2, 0), (3, 3, 0)]
        );

        let charge_case = || {
            let mut center = inp_ATOM {
                elname: [b'N' as i8, 0, 0, 0, 0, 0],
                valence: 1,
                chem_bonds_valence: 1,
                charge: 2,
                ..inp_ATOM::default()
            };
            center.neighbor[0] = 1;
            let mut hydrogen = terminal_hdt_atom(b'H', 0);
            hydrogen.charge = -1;
            vec![center, hydrogen]
        };
        let (_, historical, _) = source_remove_terminal_hdt(charge_case(), 0);
        let (_, fixed, _) = source_remove_terminal_hdt(charge_case(), 1);
        assert_eq!((historical[0].charge, historical[1].charge), (2, 0));
        assert_eq!((fixed[0].charge, fixed[1].charge), (1, 0));

        let mut deuterium = terminal_hdt_atom(b'D', 1);
        deuterium.charge = 1;
        let tritium = terminal_hdt_atom(b'T', 0);
        let (result, pair, _) = source_remove_terminal_hdt(vec![deuterium, tritium], 1);
        assert_eq!(result, 1);
        assert_eq!(
            pair.iter()
                .map(|atom| {
                    (
                        atom.elname[0],
                        atom.iso_atw_diff,
                        atom.charge,
                        atom.num_iso_H,
                        atom.valence,
                        atom.neighbor[0],
                        atom.component,
                    )
                })
                .collect::<Vec<_>>(),
            vec![
                (b'H' as i8, 3, 1, [0, 1, 0], 0, 0, 1),
                (b'H' as i8, 2, 0, [0, 0, 0], 1, 0, 0),
            ]
        );

        let mut stereo_center = inp_ATOM {
            elname: [b'C' as i8, 0, 0, 0, 0, 0],
            valence: 3,
            chem_bonds_valence: 3,
            sb_ord: [1, 0, 0],
            sn_ord: [0, 0, 0],
            sb_parity: [2, 0, 0],
            ..inp_ATOM::default()
        };
        stereo_center.neighbor[..3].copy_from_slice(&[1, 2, 3]);
        stereo_center.bond_type[..3].fill(BOND_TYPE_SINGLE as u8);
        stereo_center.bond_stereo[..3].copy_from_slice(&[1, 2, 3]);
        let hydrogen = terminal_hdt_atom(b'H', 0);
        let mut carbon_neighbor = inp_ATOM {
            elname: [b'C' as i8, 0, 0, 0, 0, 0],
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        carbon_neighbor.neighbor[0] = 0;
        let deuterium = terminal_hdt_atom(b'D', 0);
        let (result, stereo, _) = source_remove_terminal_hdt(
            vec![stereo_center, hydrogen, carbon_neighbor, deuterium],
            1,
        );
        assert_eq!(result, 2);
        assert_eq!(
            (
                stereo[0].valence,
                stereo[0].neighbor[0],
                stereo[0].bond_type[0],
                stereo[0].bond_stereo[0],
                stereo[0].num_H,
                stereo[0].num_iso_H,
                stereo[0].sb_ord,
                stereo[0].sn_ord,
                stereo[0].sb_parity,
            ),
            (1, 1, 1, 2, 1, [0, 1, 0], [0, 0, 0], [-1, 0, 0], [2, 0, 0])
        );
        assert_eq!(
            stereo[2..]
                .iter()
                .map(|atom| (atom.iso_atw_diff, atom.neighbor[0]))
                .collect::<Vec<_>>(),
            vec![(0, 0), (2, 0)]
        );
    }

    #[test]
    fn source_port__strutil__add_dt_to_num_h__line_3705() {
        let mut negative = [inp_ATOM {
            num_H: 7,
            num_iso_H: [1, 2, 3],
            ..inp_ATOM::default()
        }];
        assert_eq!(add_DT_to_num_H(-1, &mut negative), Ok(0));
        assert_eq!((negative[0].num_H, negative[0].num_iso_H), (7, [1, 2, 3]));

        let mut atoms = [
            inp_ATOM {
                num_H: 4,
                num_iso_H: [1, 2, 3],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                num_H: i8::MAX,
                num_iso_H: [1, i8::MAX, i8::MIN],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                num_H: 11,
                num_iso_H: [12, 13, 14],
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(add_DT_to_num_H(2, &mut atoms), Ok(0));
        assert_eq!((atoms[0].num_H, atoms[0].num_iso_H), (10, [1, 2, 3]));
        assert_eq!(
            (atoms[1].num_H, atoms[1].num_iso_H),
            (
                i8::MAX
                    .wrapping_add(1)
                    .wrapping_add(i8::MAX)
                    .wrapping_add(i8::MIN),
                [1, i8::MAX, i8::MIN]
            )
        );
        assert_eq!((atoms[2].num_H, atoms[2].num_iso_H), (11, [12, 13, 14]));
        assert_eq!(
            add_DT_to_num_H(4, &mut atoms),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    fn bonded_atom(neighbors: &[u16], bond_types: &[u8]) -> inp_ATOM {
        let mut atom = inp_ATOM {
            valence: neighbors.len() as i8,
            chem_bonds_valence: bond_types.iter().copied().map(i32::from).sum::<i32>() as i8,
            ..inp_ATOM::default()
        };
        atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
        atom.bond_type[..bond_types.len()].copy_from_slice(bond_types);
        atom
    }

    #[test]
    fn source_port__strutil__markdisconnectedcomponents__line_4277() {
        let mut empty_heap = SourceHeap::default();
        let preserved_lengths = empty_heap.allocate(vec![9_u16]).unwrap();
        let preserved_old = empty_heap.allocate(vec![8_u16]).unwrap();
        let mut empty = ORIG_ATOM_DATA {
            num_inp_atoms: 0,
            num_components: 7,
            nCurAtLen: preserved_lengths,
            nOldCompNumber: preserved_old,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            MarkDisconnectedComponents(&mut empty_heap, &mut empty, 1),
            Ok(0)
        );
        assert_eq!(empty.num_components, 7);
        assert_eq!(empty.nCurAtLen, preserved_lengths);
        assert_eq!(empty.nOldCompNumber, preserved_old);

        let mut atoms = vec![
            inp_ATOM::default(),
            bonded_atom(&[2], &[1]),
            bonded_atom(&[1], &[1]),
            bonded_atom(&[4], &[1]),
            bonded_atom(&[3], &[1]),
        ];
        atoms[0].component = 0;
        atoms[1].component = 1;
        atoms[2].component = 1;
        atoms[3].component = 2;
        atoms[4].component = 3;
        let mut heap = SourceHeap::default();
        let atom_pointer = heap.allocate(atoms).unwrap();
        let old_lengths = heap.allocate(vec![5_u16]).unwrap();
        let old_numbers = heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let mut original = ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: 5,
            num_components: 3,
            nCurAtLen: old_lengths,
            nOldCompNumber: old_numbers,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            MarkDisconnectedComponents(&mut heap, &mut original, 1),
            Ok(3)
        );
        assert_eq!(original.num_components, 3);
        assert_eq!(
            heap.slice(atom_pointer.as_const())
                .unwrap()
                .iter()
                .map(|atom| atom.component)
                .collect::<Vec<_>>(),
            vec![3, 1, 1, 2, 2]
        );
        assert_eq!(
            &heap.slice(original.nCurAtLen.as_const()).unwrap()[..3],
            &[2, 2, 1]
        );
        assert_eq!(
            &heap.slice(original.nOldCompNumber.as_const()).unwrap()[..3],
            &[1, 0, 0]
        );
        assert_eq!(
            heap.slice(old_lengths.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(old_numbers.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut failure_heap = SourceHeap::default();
        let failure_atoms = failure_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let failure_lengths = failure_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let failure_old = failure_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let mut failure = ORIG_ATOM_DATA {
            at: failure_atoms,
            num_inp_atoms: 1,
            num_components: 1,
            nCurAtLen: failure_lengths,
            nOldCompNumber: failure_old,
            ..ORIG_ATOM_DATA::default()
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            MarkDisconnectedComponents(&mut failure_heap, &mut failure, 0),
            Ok(-1)
        );
        assert_eq!(failure.num_components, -1);
        assert!(failure.nCurAtLen.is_null());
        assert!(failure.nOldCompNumber.is_null());
        assert_eq!(failure_heap.source_allocation_calls(), 3);
        assert_eq!(
            failure_heap.slice(failure_lengths.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            failure_heap.slice(failure_old.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__strutil__extractconnectedcomponent__line_4558() {
        let mut source = vec![
            bonded_atom(&[2], &[1]),
            bonded_atom(&[3], &[1]),
            bonded_atom(&[0], &[1]),
            bonded_atom(&[1], &[1]),
        ];
        source[0].component = 2;
        source[1].component = 1;
        source[2].component = 2;
        source[3].component = 1;
        source[0].at_type = 41;
        source[2].at_type = 42;
        let sentinel = inp_ATOM {
            at_type: 99,
            orig_compt_at_numb: 88,
            ..inp_ATOM::default()
        };
        let mut result = vec![sentinel.clone(); 4];
        let mut heap = SourceHeap::default();
        assert_eq!(
            ExtractConnectedComponent(&mut heap, &source, 4, 2, &mut result),
            Ok(2)
        );
        assert_eq!(result[0].at_type, 41);
        assert_eq!(result[1].at_type, 42);
        assert_eq!(result[0].orig_compt_at_numb, 1);
        assert_eq!(result[1].orig_compt_at_numb, 2);
        assert_eq!(result[0].neighbor[0], 1);
        assert_eq!(result[1].neighbor[0], 0);
        assert_eq!(result[2], sentinel);
        assert_eq!(result[3], sentinel);

        let before_missing_component = result.clone();
        assert_eq!(
            ExtractConnectedComponent(&mut heap, &source, 4, 99, &mut result),
            Ok(0)
        );
        assert_eq!(result, before_missing_component);

        let mut zero_heap = SourceHeap::default();
        let mut zero_result = vec![inp_ATOM::default()];
        assert_eq!(
            ExtractConnectedComponent(&mut zero_heap, &source, 0, 1, &mut zero_result),
            Ok(0)
        );

        let mut failure_heap = SourceHeap::default();
        failure_heap.fail_after_allocations(0);
        let before_failure = result.clone();
        assert_eq!(
            ExtractConnectedComponent(&mut failure_heap, &source, 4, 2, &mut result),
            Ok(CT_OUT_OF_RAM)
        );
        assert_eq!(result, before_failure);
    }

    #[test]
    fn source_port__strutil__setconnectedcomponentnumber__line_4598() {
        let mut atoms = vec![
            inp_ATOM {
                component: 11,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                component: 22,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                component: 33,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(SetConnectedComponentNumber(&mut atoms, 0, 7), Ok(0));
        assert_eq!(
            atoms.iter().map(|atom| atom.component).collect::<Vec<_>>(),
            [11, 22, 33]
        );
        assert_eq!(SetConnectedComponentNumber(&mut atoms, -1, 8), Ok(0));
        assert_eq!(
            atoms.iter().map(|atom| atom.component).collect::<Vec<_>>(),
            [11, 22, 33]
        );
        assert_eq!(SetConnectedComponentNumber(&mut atoms, 2, i32::MAX), Ok(0));
        assert_eq!(
            atoms.iter().map(|atom| atom.component).collect::<Vec<_>>(),
            [u16::MAX, u16::MAX, 33]
        );
        assert_eq!(SetConnectedComponentNumber(&mut atoms, 3, i32::MIN), Ok(0));
        assert_eq!(
            atoms.iter().map(|atom| atom.component).collect::<Vec<_>>(),
            [0, 0, 0]
        );

        assert_eq!(
            SetConnectedComponentNumber(&mut atoms, 4, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            atoms.iter().map(|atom| atom.component).collect::<Vec<_>>(),
            [0, 0, 0]
        );
    }

    #[test]
    fn source_port__strutil__alloc_inchi_stereo__line_4629() {
        let mut empty_heap = SourceHeap::default();
        empty_heap.trace_source_allocations();
        let empty = Alloc_INChI_Stereo(&mut empty_heap, 0, 0).unwrap();
        assert!(!empty.is_null());
        assert_eq!(empty_heap.source_allocation_calls(), 1);
        assert_eq!(empty_heap.live_allocation_count(), 1);
        assert_eq!(
            empty_heap.slice(empty.as_const()).unwrap()[0],
            INChI_Stereo::default()
        );

        let mut heap = SourceHeap::default();
        heap.trace_source_allocations();
        let stereo = Alloc_INChI_Stereo(&mut heap, 2, 3).unwrap();
        assert!(!stereo.is_null());
        assert_eq!(heap.source_allocation_calls(), 8);
        assert_eq!(heap.live_allocation_count(), 8);
        let value = heap.slice(stereo.as_const()).unwrap()[0].clone();
        assert_eq!(heap.slice(value.nNumber.as_const()).unwrap(), &[0, 0]);
        assert_eq!(heap.slice(value.t_parity.as_const()).unwrap(), &[0, 0]);
        assert_eq!(heap.slice(value.nNumberInv.as_const()).unwrap(), &[0, 0]);
        assert_eq!(heap.slice(value.t_parityInv.as_const()).unwrap(), &[0, 0]);
        assert_eq!(heap.slice(value.nBondAtom1.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(value.nBondAtom2.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(value.b_parity.as_const()).unwrap(), &[0, 0, 0]);

        for successful_allocations in 0..8 {
            let mut failing_heap = SourceHeap::default();
            failing_heap.fail_after_allocations(successful_allocations);
            let result = Alloc_INChI_Stereo(&mut failing_heap, 2, 3).unwrap();
            assert!(result.is_null(), "failure after {successful_allocations}");
            assert_eq!(
                failing_heap.live_allocation_count(),
                0,
                "failure after {successful_allocations}"
            );
            assert_eq!(
                failing_heap.source_allocation_calls(),
                successful_allocations + 1,
                "failure after {successful_allocations}"
            );
        }

        let mut negative_heap = SourceHeap::default();
        let negative = Alloc_INChI_Stereo(&mut negative_heap, -1, 0).unwrap();
        assert!(negative.is_null());
        assert_eq!(negative_heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__strutil__alloc_inchi__line_4722() {
        let mut unchanged_bonds = 71;
        let mut unchanged_isotopes = 72;
        let mut empty_heap = SourceHeap::default();
        assert!(
            Alloc_INChI(
                &mut empty_heap,
                &[],
                0,
                &mut unchanged_bonds,
                &mut unchanged_isotopes,
                REQ_MODE_ISO as i32,
            )
            .unwrap()
            .is_null()
        );
        assert_eq!((unchanged_bonds, unchanged_isotopes), (71, 72));
        assert_eq!(empty_heap.live_allocation_count(), 0);

        let mut atoms = vec![inp_ATOM::default(); 2];
        atoms[0].valence = 1;
        atoms[0].elname[..2].copy_from_slice(&[b'D' as i8, 0]);
        atoms[1].valence = 1;
        atoms[1].num_iso_H[2] = 1;

        let mut bonds = -1;
        let mut isotopes = -1;
        let mut plain_heap = SourceHeap::default();
        let plain = Alloc_INChI(&mut plain_heap, &atoms, 2, &mut bonds, &mut isotopes, 0).unwrap();
        assert!(!plain.is_null());
        assert_eq!((bonds, isotopes), (1, 2));
        let plain_value = plain_heap.slice(plain.as_const()).unwrap()[0].clone();
        assert_eq!(
            plain_heap.slice(plain_value.nAtom.as_const()).unwrap(),
            &[0, 0]
        );
        assert_eq!(
            plain_heap.slice(plain_value.nConnTable.as_const()).unwrap(),
            &[0, 0, 0]
        );
        assert_eq!(
            plain_heap.slice(plain_value.nTautomer.as_const()).unwrap(),
            &[0; 6]
        );
        assert_eq!(
            plain_heap.slice(plain_value.nNum_H.as_const()).unwrap(),
            &[0, 0]
        );
        assert_eq!(
            plain_heap
                .slice(plain_value.nNum_H_fixed.as_const())
                .unwrap(),
            &[0, 0]
        );
        assert!(plain_value.IsotopicAtom.is_null());
        assert!(plain_value.IsotopicTGroup.is_null());
        assert!(plain_value.nPossibleLocationsOfIsotopicH.is_null());
        assert!(!plain_value.Stereo.is_null());
        assert!(plain_value.StereoIsotopic.is_null());

        let mut isotope_heap = SourceHeap::default();
        isotope_heap.trace_source_allocations();
        let isotope = Alloc_INChI(
            &mut isotope_heap,
            &atoms,
            2,
            &mut bonds,
            &mut isotopes,
            REQ_MODE_ISO as i32,
        )
        .unwrap();
        assert!(!isotope.is_null());
        assert_eq!(isotope_heap.source_allocation_calls(), 25);
        assert_eq!(isotope_heap.live_allocation_count(), 25);
        let isotope_value = isotope_heap.slice(isotope.as_const()).unwrap()[0].clone();
        assert_eq!(
            isotope_heap
                .slice(isotope_value.IsotopicAtom.as_const())
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            isotope_heap
                .slice(isotope_value.IsotopicTGroup.as_const())
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            isotope_heap
                .slice(isotope_value.nPossibleLocationsOfIsotopicH.as_const())
                .unwrap(),
            &[0, 0, 0]
        );
        assert!(!isotope_value.Stereo.is_null());
        assert!(!isotope_value.StereoIsotopic.is_null());

        for successful_allocations in 0..25 {
            let mut failing_heap = SourceHeap::default();
            failing_heap.fail_after_allocations(successful_allocations);
            let mut failure_bonds = 91;
            let mut failure_isotopes = 92;
            let result = Alloc_INChI(
                &mut failing_heap,
                &atoms,
                2,
                &mut failure_bonds,
                &mut failure_isotopes,
                REQ_MODE_ISO as i32,
            )
            .unwrap();
            assert!(result.is_null(), "failure after {successful_allocations}");
            assert_eq!(
                failing_heap.live_allocation_count(),
                0,
                "failure after {successful_allocations}"
            );
            if successful_allocations == 0 {
                assert_eq!((failure_bonds, failure_isotopes), (91, 92));
            } else {
                assert_eq!((failure_bonds, failure_isotopes), (1, 2));
            }
        }
    }

    #[test]
    fn source_port__strutil__alloc_inchi_aux__line_4884() {
        let mut empty_heap = SourceHeap::default();
        assert!(
            Alloc_INChI_Aux(&mut empty_heap, 0, 1, REQ_MODE_ISO as i32, 1)
                .unwrap()
                .is_null()
        );
        assert_eq!(empty_heap.live_allocation_count(), 0);

        let mut plain_heap = SourceHeap::default();
        let plain = Alloc_INChI_Aux(&mut plain_heap, 2, 0, 0, 0).unwrap();
        assert!(!plain.is_null());
        let value = plain_heap.slice(plain.as_const()).unwrap()[0].clone();
        assert_eq!(
            plain_heap
                .slice(value.nOrigAtNosInCanonOrd.as_const())
                .unwrap(),
            &[0; 3]
        );
        assert_eq!(
            plain_heap
                .slice(value.nOrigAtNosInCanonOrdInv.as_const())
                .unwrap(),
            &[0; 3]
        );
        assert_eq!(
            plain_heap
                .slice(value.nConstitEquNumbers.as_const())
                .unwrap(),
            &[0; 3]
        );
        assert_eq!(
            plain_heap
                .slice(value.nConstitEquTGroupNumbers.as_const())
                .unwrap(),
            &[0; 2]
        );
        assert_eq!(
            plain_heap.slice(value.OrigInfo.as_const()).unwrap().len(),
            2
        );
        assert!(value.szOrigCoord.is_null());
        assert!(value.nIsotopicOrigAtNosInCanonOrd.is_null());

        let mut full_heap = SourceHeap::default();
        full_heap.trace_source_allocations();
        let full = Alloc_INChI_Aux(&mut full_heap, 2, 1, REQ_MODE_ISO as i32, 1).unwrap();
        assert!(!full.is_null());
        assert_eq!(full_heap.source_allocation_calls(), 11);
        assert_eq!(full_heap.live_allocation_count(), 11);
        let value = full_heap.slice(full.as_const()).unwrap()[0].clone();
        assert_eq!(
            full_heap.slice(value.szOrigCoord.as_const()).unwrap(),
            &[[0_i8; 32]; 2]
        );
        assert_eq!(
            full_heap
                .slice(value.nIsotopicOrigAtNosInCanonOrd.as_const())
                .unwrap(),
            &[0; 3]
        );
        assert_eq!(
            full_heap
                .slice(value.nIsotopicOrigAtNosInCanonOrdInv.as_const())
                .unwrap(),
            &[0; 3]
        );
        assert_eq!(
            full_heap
                .slice(value.nConstitEquIsotopicNumbers.as_const())
                .unwrap(),
            &[0; 3]
        );
        assert_eq!(
            full_heap
                .slice(value.nConstitEquIsotopicTGroupNumbers.as_const())
                .unwrap(),
            &[0; 2]
        );

        for successful_allocations in 0..11 {
            let mut failing_heap = SourceHeap::default();
            failing_heap.fail_after_allocations(successful_allocations);
            let result = Alloc_INChI_Aux(&mut failing_heap, 2, 1, REQ_MODE_ISO as i32, 1).unwrap();
            assert!(result.is_null(), "failure after {successful_allocations}");
            assert_eq!(failing_heap.live_allocation_count(), 0);
        }

        for successful_allocations in 6..10 {
            let mut tolerated_heap = SourceHeap::default();
            tolerated_heap.fail_after_allocations(successful_allocations);
            let result =
                Alloc_INChI_Aux(&mut tolerated_heap, 2, 0, REQ_MODE_ISO as i32, 0).unwrap();
            assert!(!result.is_null(), "failure after {successful_allocations}");
            let value = tolerated_heap.slice(result.as_const()).unwrap()[0].clone();
            match successful_allocations {
                6 => {
                    assert!(value.nIsotopicOrigAtNosInCanonOrd.is_null());
                    assert!(value.nIsotopicOrigAtNosInCanonOrdInv.is_null());
                    assert!(value.nConstitEquIsotopicNumbers.is_null());
                    assert!(!value.nConstitEquIsotopicTGroupNumbers.is_null());
                }
                7 => {
                    assert!(!value.nIsotopicOrigAtNosInCanonOrd.is_null());
                    assert!(value.nIsotopicOrigAtNosInCanonOrdInv.is_null());
                    assert!(value.nConstitEquIsotopicNumbers.is_null());
                    assert!(!value.nConstitEquIsotopicTGroupNumbers.is_null());
                }
                8 => {
                    assert!(!value.nIsotopicOrigAtNosInCanonOrd.is_null());
                    assert!(!value.nIsotopicOrigAtNosInCanonOrdInv.is_null());
                    assert!(value.nConstitEquIsotopicNumbers.is_null());
                    assert!(!value.nConstitEquIsotopicTGroupNumbers.is_null());
                }
                9 => {
                    assert!(!value.nIsotopicOrigAtNosInCanonOrd.is_null());
                    assert!(!value.nIsotopicOrigAtNosInCanonOrdInv.is_null());
                    assert!(!value.nConstitEquIsotopicNumbers.is_null());
                    assert!(value.nConstitEquIsotopicTGroupNumbers.is_null());
                }
                _ => unreachable!(),
            }
        }
    }

    #[test]
    fn source_port__strutil__removeinpatbond__line_2046() {
        let mut isolated = vec![inp_ATOM::default()];
        assert_eq!(RemoveInpAtBond(&mut isolated, 0, i32::MIN), Ok(0));
        assert_eq!(isolated, vec![inp_ATOM::default()]);

        let mut atoms = vec![
            bonded_atom(&[1, 2, 3], &[1, 9, 3]),
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        atoms[0].bond_stereo[..3].copy_from_slice(&[4, 5, 6]);
        assert_eq!(RemoveInpAtBond(&mut atoms, 0, 1), Ok(1));
        assert_eq!(atoms[0].neighbor[..4], [1, 3, 0, 0]);
        assert_eq!(atoms[0].bond_type[..4], [1, 3, 0, 0]);
        assert_eq!(atoms[0].bond_stereo[..4], [4, 6, 0, 0]);
        assert_eq!(atoms[0].valence, 2);
        assert_eq!(atoms[0].chem_bonds_valence, 12);

        let mut tetra = vec![bonded_atom(&[1], &[1]), inp_ATOM::default()];
        tetra[0].orig_at_number = 10;
        tetra[0].p_parity = 1;
        tetra[0].p_orig_at_num = [10, 20, 30, 40];
        assert_eq!(RemoveInpAtBond(&mut tetra, 0, 0), Ok(1));
        assert_eq!(tetra[0].p_parity, 0);

        let mut tetra_replace = vec![
            bonded_atom(&[1, 2], &[1, 1]),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        tetra_replace[0].orig_at_number = 10;
        tetra_replace[0].p_parity = 2;
        tetra_replace[0].p_orig_at_num = [20, 30, 40, 50];
        tetra_replace[1].orig_at_number = 20;
        assert_eq!(RemoveInpAtBond(&mut tetra_replace, 0, 0), Ok(1));
        assert_eq!(tetra_replace[0].p_parity, 2);
        assert_eq!(tetra_replace[0].p_orig_at_num[0], 10);

        let mut tetra_wrong = vec![bonded_atom(&[1], &[1]), inp_ATOM::default()];
        tetra_wrong[0].orig_at_number = 10;
        tetra_wrong[0].p_parity = 2;
        tetra_wrong[0].p_orig_at_num = [20, 30, 40, 50];
        tetra_wrong[1].orig_at_number = 99;
        assert_eq!(RemoveInpAtBond(&mut tetra_wrong, 0, 0), Ok(1));
        assert_eq!(tetra_wrong[0].p_parity, 0);

        let mut paired = vec![
            bonded_atom(&[1, 2], &[2, 1]),
            bonded_atom(&[0], &[2]),
            inp_ATOM::default(),
        ];
        paired[0].sb_parity = [1, 2, 0];
        paired[0].sb_ord = [0, 1, 0];
        paired[0].sn_ord = [1, 0, 0];
        paired[0].sn_orig_at_num = [8, 9, 0];
        paired[1].sb_parity[0] = 2;
        paired[1].sb_ord[0] = 0;
        assert_eq!(RemoveInpAtBond(&mut paired, 0, 0), Ok(1));
        assert_eq!(paired[0].sb_parity, [2, 0, 0]);
        assert_eq!(paired[0].sb_ord, [1, 0, 0]);
        assert_eq!(paired[1].sb_parity, [0, 0, 0]);

        let mut replace_neighbor = vec![
            bonded_atom(&[1, 2, 3], &[2, 1, 1]),
            bonded_atom(&[0], &[2]),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        replace_neighbor[0].sb_parity[0] = 1;
        replace_neighbor[0].sb_ord[0] = 0;
        replace_neighbor[0].sn_ord[0] = 1;
        replace_neighbor[3].orig_at_number = 44;
        assert_eq!(RemoveInpAtBond(&mut replace_neighbor, 0, 1), Ok(1));
        assert_eq!(replace_neighbor[0].sb_parity[0], 2);
        assert_eq!(replace_neighbor[0].sn_ord[0], 1);
        assert_eq!(replace_neighbor[0].sn_orig_at_num[0], 44);

        let mut undefined_neighbor = vec![
            bonded_atom(&[1, 2], &[2, 1]),
            bonded_atom(&[0], &[2]),
            inp_ATOM::default(),
        ];
        undefined_neighbor[0].sb_parity[0] = AB_PARITY_UNDF as i8;
        undefined_neighbor[0].sb_ord[0] = 0;
        undefined_neighbor[0].sn_ord[0] = 1;
        assert_eq!(RemoveInpAtBond(&mut undefined_neighbor, 0, 1), Ok(1));
        assert_eq!(undefined_neighbor[0].sb_parity[0], AB_PARITY_UNDF as i8);
        assert_eq!(undefined_neighbor[0].sn_ord[0], -99);
        assert_eq!(undefined_neighbor[0].sn_orig_at_num[0], 0);

        let mut other_neighbor = vec![
            bonded_atom(&[1, 2, 3, 4], &[1, 1, 2, 1]),
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
            inp_ATOM::default(),
        ];
        other_neighbor[0].sb_parity[0] = 1;
        other_neighbor[0].sb_ord[0] = 2;
        other_neighbor[0].sn_ord[0] = 1;
        assert_eq!(RemoveInpAtBond(&mut other_neighbor, 0, 0), Ok(1));
        assert_eq!(other_neighbor[0].sb_parity[0], 1);
        assert_eq!(other_neighbor[0].sb_ord[0], 1);
        assert_eq!(other_neighbor[0].sn_ord[0], 0);

        let mut invalid = vec![bonded_atom(&[1], &[1]), inp_ATOM::default()];
        assert_eq!(
            RemoveInpAtBond(&mut invalid, -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            RemoveInpAtBond(&mut invalid, 0, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            RemoveInpAtBond(&mut invalid, 0, 20),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__strutil__disconnectinpatbond__line_2261() {
        let mut atoms = vec![bonded_atom(&[1], &[1]), bonded_atom(&[0], &[1])];
        atoms[0].component = 1;
        atoms[1].component = 3;
        let mut components = [11_u16, 22, 33, 44];
        assert_eq!(
            DisconnectInpAtBond(&mut atoms, Some(&mut components), 0, 0),
            Ok(1)
        );
        assert_eq!(atoms[0].valence, 0);
        assert_eq!(atoms[1].valence, 0);
        assert_eq!(atoms[0].chem_bonds_valence, 0);
        assert_eq!(atoms[1].chem_bonds_valence, 0);
        assert_eq!(components, [0, 22, 0, 44]);

        let mut no_components = vec![bonded_atom(&[1], &[2]), bonded_atom(&[0], &[2])];
        assert_eq!(DisconnectInpAtBond(&mut no_components, None, 0, 0), Ok(1));

        let mut missing_reverse = vec![bonded_atom(&[1], &[1]), bonded_atom(&[1], &[1])];
        let before = missing_reverse.clone();
        assert_eq!(DisconnectInpAtBond(&mut missing_reverse, None, 0, 0), Ok(0));
        assert_eq!(missing_reverse, before);

        let mut partial = vec![inp_ATOM::default(), bonded_atom(&[0], &[1])];
        partial[0].neighbor[0] = 1;
        assert_eq!(DisconnectInpAtBond(&mut partial, None, 0, 0), Ok(0));
        assert_eq!(partial[0].valence, 0);
        assert_eq!(partial[1].valence, 0);

        let mut component_overflow = vec![bonded_atom(&[1], &[1]), bonded_atom(&[0], &[1])];
        component_overflow[0].component = 2;
        let mut short_components = [9_u16];
        assert_eq!(
            DisconnectInpAtBond(&mut component_overflow, Some(&mut short_components), 0, 0,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(component_overflow[0].valence, 0);
        assert_eq!(component_overflow[1].valence, 0);

        let mut invalid = vec![inp_ATOM::default()];
        assert_eq!(
            DisconnectInpAtBond(&mut invalid, None, -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            DisconnectInpAtBond(&mut invalid, None, 0, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        invalid[0].neighbor[0] = 7;
        assert_eq!(
            DisconnectInpAtBond(&mut invalid, None, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__strutil__imat_new__line_5005() {
        let mut heap = SourceHeap::default();
        let old_row = heap.allocate(vec![9_i32, 8]).unwrap();
        let old_row_2 = heap.allocate(vec![7_i32, 6]).unwrap();
        let old_outer = heap.allocate(vec![old_row, old_row_2]).unwrap();
        let mut matrix = old_outer;
        assert_eq!(imat_new(&mut heap, 0, 3, &mut matrix), Ok(0));
        assert_eq!(matrix, old_outer);
        assert_eq!(heap.slice(old_row.as_const()).unwrap(), &[9, 8]);
        assert_eq!(heap.slice(old_row_2.as_const()).unwrap(), &[7, 6]);

        assert_eq!(imat_new(&mut heap, 2, 3, &mut matrix), Ok(0));
        assert_ne!(matrix, old_outer);
        let rows = heap.slice(matrix.as_const()).unwrap().to_vec();
        assert_eq!(rows.len(), 2);
        assert_eq!(heap.slice(rows[0].as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(rows[1].as_const()).unwrap(), &[0, 0, 0]);
        heap.slice_mut(rows[1]).unwrap()[2] = 17;
        assert_eq!(heap.slice(rows[1].as_const()).unwrap(), &[0, 0, 17]);
        imat_free(&mut heap, 2, matrix).unwrap();
        assert_eq!(
            heap.slice(old_outer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(old_row.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(old_row_2.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut outer_failure_heap = SourceHeap::default();
        outer_failure_heap.fail_after_allocations(0);
        let mut outer_failure = SourceMutPointer::null();
        assert_eq!(
            imat_new(&mut outer_failure_heap, 2, 2, &mut outer_failure),
            Ok(1)
        );
        assert!(outer_failure.is_null());

        let mut row_failure_heap = SourceHeap::default();
        row_failure_heap.fail_after_allocations(1);
        let mut row_failure = SourceMutPointer::null();
        assert_eq!(
            imat_new(&mut row_failure_heap, 2, 2, &mut row_failure),
            Ok(1)
        );
        assert!(!row_failure.is_null());
        let partial_rows = row_failure_heap.slice(row_failure.as_const()).unwrap();
        assert!(partial_rows[0].is_null());
        imat_free(&mut row_failure_heap, 2, row_failure).unwrap();
    }

    #[test]
    fn source_port__strutil__free_inchi_stereo__line_4611() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            Free_INChI_Stereo(&mut heap, SourceMutPointer::null()),
            Ok(0)
        );

        let n_number = heap.allocate(vec![1_u16, 2]).unwrap();
        let t_parity = heap.allocate(vec![1_i8, 2]).unwrap();
        let n_number_inv = heap.allocate(vec![3_u16, 4]).unwrap();
        let t_parity_inv = heap.allocate(vec![3_i8, 4]).unwrap();
        let n_bond_atom_1 = heap.allocate(vec![5_u16]).unwrap();
        let n_bond_atom_2 = heap.allocate(vec![6_u16]).unwrap();
        let b_parity = heap.allocate(vec![7_i8]).unwrap();
        let stereo = heap
            .allocate(vec![INChI_Stereo {
                nNumberOfStereoCenters: 2,
                nNumber: n_number,
                t_parity,
                nNumberInv: n_number_inv,
                t_parityInv: t_parity_inv,
                nCompInv2Abs: -1,
                bTrivialInv: 1,
                nNumberOfStereoBonds: 1,
                nBondAtom1: n_bond_atom_1,
                nBondAtom2: n_bond_atom_2,
                b_parity,
            }])
            .unwrap();

        assert_eq!(Free_INChI_Stereo(&mut heap, stereo), Ok(0));
        let value = &heap.slice(stereo.as_const()).unwrap()[0];
        assert_eq!(value.nNumberOfStereoCenters, 2);
        assert_eq!(value.nCompInv2Abs, -1);
        assert_eq!(value.bTrivialInv, 1);
        assert_eq!(value.nNumberOfStereoBonds, 1);
        assert!(value.nNumber.is_null());
        assert!(value.t_parity.is_null());
        assert!(value.nNumberInv.is_null());
        assert!(value.t_parityInv.is_null());
        assert!(value.nBondAtom1.is_null());
        assert!(value.nBondAtom2.is_null());
        assert!(value.b_parity.is_null());
        for freed in [
            heap.slice(n_number.as_const()).map(|_| ()),
            heap.slice(t_parity.as_const()).map(|_| ()),
            heap.slice(n_number_inv.as_const()).map(|_| ()),
            heap.slice(t_parity_inv.as_const()).map(|_| ()),
            heap.slice(n_bond_atom_1.as_const()).map(|_| ()),
            heap.slice(n_bond_atom_2.as_const()).map(|_| ()),
            heap.slice(b_parity.as_const()).map(|_| ()),
        ] {
            assert_eq!(freed, Err(SourceHeapError::MissingAllocation));
        }

        assert_eq!(Free_INChI_Stereo(&mut heap, stereo), Ok(0));
        assert!(heap.slice(stereo.as_const()).is_ok());
    }

    #[test]
    fn source_port__strutil__free_inchi__line_4675() {
        let mut heap = SourceHeap::default();
        let mut null = SourceMutPointer::null();
        assert_eq!(Free_INChI(&mut heap, &mut null), Ok(0));
        assert!(null.is_null());

        let retained_member = heap.allocate(vec![6_u8]).unwrap();
        let mut retained = heap
            .allocate(vec![INChI {
                nRefCount: 1,
                nAtom: retained_member,
                ..INChI::default()
            }])
            .unwrap();
        assert_eq!(Free_INChI(&mut heap, &mut retained), Ok(1));
        assert_eq!(heap.slice(retained.as_const()).unwrap()[0].nRefCount, 0);
        assert_eq!(heap.slice(retained_member.as_const()).unwrap(), &[6]);
        assert_eq!(Free_INChI(&mut heap, &mut retained), Ok(0));
        assert!(retained.is_null());
        assert_eq!(
            heap.slice(retained_member.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        for reference_count in [0, -1, i32::MIN] {
            let member = heap.allocate(vec![7_u8]).unwrap();
            let pointer = heap
                .allocate(vec![INChI {
                    nRefCount: reference_count,
                    nAtom: member,
                    ..INChI::default()
                }])
                .unwrap();
            let mut owner = pointer;
            assert_eq!(Free_INChI(&mut heap, &mut owner), Ok(0));
            assert!(owner.is_null());
            assert_eq!(
                heap.slice(pointer.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
            assert_eq!(
                heap.slice(member.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
        }
    }

    #[test]
    fn source_port__strutil__free_inchi_aux__line_4843() {
        let mut heap = SourceHeap::default();
        let mut null = SourceMutPointer::null();
        assert_eq!(Free_INChI_Aux(&mut heap, &mut null), Ok(0));

        let retained_member = heap.allocate(vec![1_u16]).unwrap();
        let mut retained = heap
            .allocate(vec![INChI_Aux {
                nRefCount: 1,
                nOrigAtNosInCanonOrd: retained_member,
                ..INChI_Aux::default()
            }])
            .unwrap();
        assert_eq!(Free_INChI_Aux(&mut heap, &mut retained), Ok(1));
        assert_eq!(heap.slice(retained.as_const()).unwrap()[0].nRefCount, 0);
        assert_eq!(heap.slice(retained_member.as_const()).unwrap(), &[1]);
        assert_eq!(Free_INChI_Aux(&mut heap, &mut retained), Ok(0));
        assert!(retained.is_null());

        let original = heap.allocate(vec![1_u16]).unwrap();
        let isotopic_original = heap.allocate(vec![2_u16]).unwrap();
        let original_inverse = heap.allocate(vec![3_u16]).unwrap();
        let isotopic_inverse = heap.allocate(vec![4_u16]).unwrap();
        let coordinates = heap.allocate(vec![[0_i8; 32]]).unwrap();
        let original_info = heap
            .allocate(vec![crate::source_types::ORIG_INFO::default()])
            .unwrap();
        let equivalent = heap.allocate(vec![5_u16]).unwrap();
        let equivalent_groups = heap.allocate(vec![6_u16]).unwrap();
        let isotopic_equivalent = heap.allocate(vec![7_u16]).unwrap();
        let isotopic_equivalent_groups = heap.allocate(vec![8_u16]).unwrap();
        let pointer = heap
            .allocate(vec![INChI_Aux {
                nRefCount: i32::MIN,
                nOrigAtNosInCanonOrd: original,
                nIsotopicOrigAtNosInCanonOrd: isotopic_original,
                nOrigAtNosInCanonOrdInv: original_inverse,
                nIsotopicOrigAtNosInCanonOrdInv: isotopic_inverse,
                szOrigCoord: coordinates,
                OrigInfo: original_info,
                nConstitEquNumbers: equivalent,
                nConstitEquTGroupNumbers: equivalent_groups,
                nConstitEquIsotopicNumbers: isotopic_equivalent,
                nConstitEquIsotopicTGroupNumbers: isotopic_equivalent_groups,
                ..INChI_Aux::default()
            }])
            .unwrap();
        let mut owner = pointer;
        assert_eq!(Free_INChI_Aux(&mut heap, &mut owner), Ok(0));
        assert!(owner.is_null());
        assert_eq!(
            heap.slice(pointer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(original.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(isotopic_original.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(original_inverse.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(isotopic_inverse.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(coordinates.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(original_info.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(equivalent.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(equivalent_groups.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(isotopic_equivalent.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(isotopic_equivalent_groups.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__strutil__free_inchi_members__line_4698() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            Free_INChI_Members(&mut heap, SourceMutPointer::null()),
            Ok(0)
        );

        let stereo_member = heap.allocate(vec![1_u16]).unwrap();
        let stereo = heap
            .allocate(vec![INChI_Stereo {
                nNumber: stereo_member,
                ..INChI_Stereo::default()
            }])
            .unwrap();
        let isotopic_member = heap.allocate(vec![2_i8]).unwrap();
        let stereo_isotopic = heap
            .allocate(vec![INChI_Stereo {
                b_parity: isotopic_member,
                ..INChI_Stereo::default()
            }])
            .unwrap();
        let hill = heap.allocate(vec![b'C' as i8, 0]).unwrap();
        let atom = heap.allocate(vec![6_u8]).unwrap();
        let connection = heap.allocate(vec![1_u16, 2]).unwrap();
        let tautomer = heap.allocate(vec![3_u16]).unwrap();
        let num_h = heap.allocate(vec![1_i8]).unwrap();
        let num_h_fixed = heap.allocate(vec![2_i8]).unwrap();
        let isotopic_atom = heap.allocate(vec![INChI_IsotopicAtom::default()]).unwrap();
        let isotopic_t_group = heap
            .allocate(vec![INChI_IsotopicTGroup::default()])
            .unwrap();
        let possible_h = heap.allocate(vec![4_u16]).unwrap();
        let inchi = heap
            .allocate(vec![INChI {
                nErrorCode: 17,
                nNumberOfAtoms: 1,
                szHillFormula: hill,
                nAtom: atom,
                nConnTable: connection,
                nTautomer: tautomer,
                nNum_H: num_h,
                nNum_H_fixed: num_h_fixed,
                IsotopicAtom: isotopic_atom,
                IsotopicTGroup: isotopic_t_group,
                Stereo: stereo,
                StereoIsotopic: stereo_isotopic,
                nPossibleLocationsOfIsotopicH: possible_h,
                ..INChI::default()
            }])
            .unwrap();

        assert_eq!(Free_INChI_Members(&mut heap, inchi), Ok(0));
        let value = &heap.slice(inchi.as_const()).unwrap()[0];
        assert_eq!(value.nErrorCode, 17);
        assert_eq!(value.nNumberOfAtoms, 1);
        assert!(value.szHillFormula.is_null());
        assert!(value.nAtom.is_null());
        assert!(value.nConnTable.is_null());
        assert!(value.nTautomer.is_null());
        assert!(value.nNum_H.is_null());
        assert!(value.nNum_H_fixed.is_null());
        assert!(value.IsotopicAtom.is_null());
        assert!(value.IsotopicTGroup.is_null());
        assert!(value.Stereo.is_null());
        assert!(value.StereoIsotopic.is_null());
        assert!(value.nPossibleLocationsOfIsotopicH.is_null());
        for freed in [
            heap.slice(stereo_member.as_const()).map(|_| ()),
            heap.slice(isotopic_member.as_const()).map(|_| ()),
            heap.slice(stereo.as_const()).map(|_| ()),
            heap.slice(stereo_isotopic.as_const()).map(|_| ()),
            heap.slice(hill.as_const()).map(|_| ()),
            heap.slice(atom.as_const()).map(|_| ()),
            heap.slice(connection.as_const()).map(|_| ()),
            heap.slice(tautomer.as_const()).map(|_| ()),
            heap.slice(num_h.as_const()).map(|_| ()),
            heap.slice(num_h_fixed.as_const()).map(|_| ()),
            heap.slice(isotopic_atom.as_const()).map(|_| ()),
            heap.slice(isotopic_t_group.as_const()).map(|_| ()),
            heap.slice(possible_h.as_const()).map(|_| ()),
        ] {
            assert_eq!(freed, Err(SourceHeapError::MissingAllocation));
        }
        assert_eq!(Free_INChI_Members(&mut heap, inchi), Ok(0));
    }

    #[test]
    fn source_port__strutil__imat_free__line_5036() {
        let mut heap = SourceHeap::default();

        assert_eq!(
            imat_free(
                &mut heap,
                i32::MAX,
                SourceMutPointer::<SourceMutPointer<i32>>::null(),
            ),
            Ok(())
        );
        assert_eq!(
            imat_free(
                &mut heap,
                i32::MIN,
                SourceMutPointer::<SourceMutPointer<i32>>::null(),
            ),
            Ok(())
        );

        let zero_count_row = heap.allocate(vec![10_i32]).unwrap();
        let zero_count_outer = heap.allocate(vec![zero_count_row]).unwrap();
        imat_free(&mut heap, 0, zero_count_outer).unwrap();
        assert_eq!(
            heap.slice(zero_count_outer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(zero_count_row.as_const()).unwrap(), &[10]);
        heap.free(zero_count_row).unwrap();

        let negative_count_row = heap.allocate(vec![20_i32]).unwrap();
        let negative_count_outer = heap.allocate(vec![negative_count_row]).unwrap();
        imat_free(&mut heap, i32::MIN, negative_count_outer).unwrap();
        assert_eq!(
            heap.slice(negative_count_outer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(negative_count_row.as_const()).unwrap(), &[20]);
        heap.free(negative_count_row).unwrap();

        let first_row = heap.allocate(vec![1_i32, 2]).unwrap();
        let unvisited_row = heap.allocate(vec![3_i32, 4]).unwrap();
        let partial_outer = heap.allocate(vec![first_row, unvisited_row]).unwrap();
        imat_free(&mut heap, 1, partial_outer).unwrap();
        assert_eq!(
            heap.slice(first_row.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(unvisited_row.as_const()).unwrap(), &[3, 4]);
        assert_eq!(
            heap.slice(partial_outer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        heap.free(unvisited_row).unwrap();

        let row_zero = heap.allocate(vec![5_i32, 6]).unwrap();
        let row_two = heap.allocate(vec![7_i32, 8]).unwrap();
        let complete_outer = heap
            .allocate(vec![row_zero, SourceMutPointer::<i32>::null(), row_two])
            .unwrap();
        imat_free(&mut heap, 3, complete_outer).unwrap();
        assert_eq!(
            heap.slice(row_zero.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(row_two.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(complete_outer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__strutil__subgraf_free__line_5155() {
        let mut heap = SourceHeap::default();
        assert_eq!(subgraf_free(&mut heap, SourceMutPointer::null()), Ok(()));

        let nodes = heap.allocate(vec![1_i32, 2]).unwrap();
        let degrees = heap.allocate(vec![1_i32, 1]).unwrap();
        let orig2node = heap.allocate(vec![-1_i32, 0, 1]).unwrap();
        let first_row = heap
            .allocate(vec![subgraf_edge { nbr: 1, etype: 1 }])
            .unwrap();
        let adjacency = heap
            .allocate(vec![first_row, SourceMutPointer::null()])
            .unwrap();
        let graph = heap
            .allocate(vec![subgraf {
                nnodes: 2,
                nodes,
                degrees,
                orig2node,
                adj: adjacency,
            }])
            .unwrap();
        assert_eq!(subgraf_free(&mut heap, graph), Ok(()));
        for missing in [nodes, degrees, orig2node] {
            assert_eq!(
                heap.slice(missing.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
        }
        assert_eq!(
            heap.slice(first_row.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(adjacency.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(graph.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let retained_row = heap
            .allocate(vec![subgraf_edge { nbr: 7, etype: 3 }])
            .unwrap();
        let negative_adjacency = heap.allocate(vec![retained_row]).unwrap();
        let negative_graph = heap
            .allocate(vec![subgraf {
                nnodes: i32::MIN,
                adj: negative_adjacency,
                ..subgraf::default()
            }])
            .unwrap();
        assert_eq!(subgraf_free(&mut heap, negative_graph), Ok(()));
        assert_eq!(
            heap.slice(negative_adjacency.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(negative_graph.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(retained_row.as_const()).unwrap(),
            &[subgraf_edge { nbr: 7, etype: 3 }]
        );
        heap.free(retained_row).unwrap();

        let partial_graph = heap
            .allocate(vec![subgraf {
                nnodes: i32::MAX,
                ..subgraf::default()
            }])
            .unwrap();
        assert_eq!(subgraf_free(&mut heap, partial_graph), Ok(()));
    }

    #[test]
    fn source_port__strutil__subgraf_new__line_5064() {
        fn graph_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 3];
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[0].bond_type[..2].copy_from_slice(&[1, 2]);
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].bond_type[0] = 1;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].bond_type[0] = 2;
            atoms
        }

        let mut heap = SourceHeap::default();
        let atom_pointer = heap.allocate_model_storage(graph_atoms()).unwrap();
        let original = ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: 3,
            ..ORIG_ATOM_DATA::default()
        };
        let graph = subgraf_new(&mut heap, &original, 2, &[1, 2]).unwrap();
        assert!(!graph.is_null());
        let value = heap.slice(graph.as_const()).unwrap()[0].clone();
        assert_eq!(value.nnodes, 2);
        assert_eq!(heap.slice(value.nodes.as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(value.degrees.as_const()).unwrap(), &[1, 1]);
        assert_eq!(
            heap.slice(value.orig2node.as_const()).unwrap(),
            &[-1, 0, 1, -1]
        );
        let rows = heap.slice(value.adj.as_const()).unwrap();
        assert_eq!(
            heap.slice(rows[0].as_const()).unwrap(),
            &[subgraf_edge { nbr: 1, etype: 1 }, subgraf_edge::default(),]
        );
        assert_eq!(
            heap.slice(rows[1].as_const()).unwrap(),
            &[subgraf_edge { nbr: 0, etype: 1 }]
        );
        subgraf_free(&mut heap, graph).unwrap();

        let empty = subgraf_new(&mut heap, &original, 0, &[]).unwrap();
        assert!(!empty.is_null());
        let empty_value = heap.slice(empty.as_const()).unwrap()[0].clone();
        assert_eq!(empty_value.nnodes, 0);
        assert!(!empty_value.nodes.is_null());
        assert!(!empty_value.degrees.is_null());
        assert!(!empty_value.adj.is_null());
        assert!(heap.slice(empty_value.nodes.as_const()).unwrap().is_empty());
        assert_eq!(
            heap.slice(empty_value.orig2node.as_const()).unwrap(),
            &[-1, -1, -1, -1]
        );
        subgraf_free(&mut heap, empty).unwrap();

        for successful_allocations in 0_u64..7 {
            let mut failure_heap = SourceHeap::default();
            let failure_atoms = failure_heap.allocate_model_storage(graph_atoms()).unwrap();
            let failure_original = ORIG_ATOM_DATA {
                at: failure_atoms,
                num_inp_atoms: 3,
                ..ORIG_ATOM_DATA::default()
            };
            let baseline = failure_heap.live_allocation_count();
            failure_heap.fail_after_allocations(successful_allocations);
            assert!(
                subgraf_new(&mut failure_heap, &failure_original, 2, &[1, 2])
                    .unwrap()
                    .is_null(),
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                failure_heap.live_allocation_count(),
                baseline,
                "allocation ordinal {successful_allocations}"
            );
        }
    }

    #[test]
    fn source_port__strutil__subgraf_pathfinder_new__line_5221() {
        let mut heap = SourceHeap::default();
        let graph = heap
            .allocate_model_storage(vec![subgraf {
                nnodes: 3,
                ..subgraf::default()
            }])
            .unwrap();
        let pathfinder = subgraf_pathfinder_new(
            &mut heap,
            graph,
            &ORIG_ATOM_DATA::default(),
            i32::MIN,
            i32::MAX,
        )
        .unwrap();
        assert!(!pathfinder.is_null());
        let value = heap.slice(pathfinder.as_const()).unwrap()[0].clone();
        assert_eq!(value.sg, graph);
        assert_eq!(value.start, i32::MIN);
        assert_eq!(value.end, i32::MAX);
        assert_eq!(value.maxbonds, 0);
        assert_eq!(value.nbonds, 0);
        assert_eq!(value.nseen, 0);
        assert_eq!(heap.slice(value.seen.as_const()).unwrap(), &[0, 0, 0]);
        inchi_free(&mut heap, value.seen).unwrap();
        inchi_free(&mut heap, pathfinder).unwrap();

        let zero_graph = heap
            .allocate_model_storage(vec![subgraf::default()])
            .unwrap();
        let zero = subgraf_pathfinder_new(&mut heap, zero_graph, &ORIG_ATOM_DATA::default(), 0, 0)
            .unwrap();
        let zero_value = heap.slice(zero.as_const()).unwrap()[0].clone();
        assert!(!zero_value.seen.is_null());
        assert!(heap.slice(zero_value.seen.as_const()).unwrap().is_empty());
        inchi_free(&mut heap, zero_value.seen).unwrap();
        inchi_free(&mut heap, zero).unwrap();

        for successful_allocations in 0_u64..2 {
            let mut failure_heap = SourceHeap::default();
            let failure_graph = failure_heap
                .allocate_model_storage(vec![subgraf {
                    nnodes: 3,
                    ..subgraf::default()
                }])
                .unwrap();
            let baseline = failure_heap.live_allocation_count();
            failure_heap.fail_after_allocations(successful_allocations);
            assert!(
                subgraf_pathfinder_new(
                    &mut failure_heap,
                    failure_graph,
                    &ORIG_ATOM_DATA::default(),
                    1,
                    2,
                )
                .unwrap()
                .is_null()
            );
            assert_eq!(failure_heap.live_allocation_count(), baseline);
        }

        let negative_graph = heap
            .allocate_model_storage(vec![subgraf {
                nnodes: -1,
                ..subgraf::default()
            }])
            .unwrap();
        assert!(
            subgraf_pathfinder_new(&mut heap, negative_graph, &ORIG_ATOM_DATA::default(), 0, 0,)
                .unwrap()
                .is_null()
        );
    }

    #[test]
    fn source_port__strutil__subgraf_pathfinder_free__line_5252() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            subgraf_pathfinder_free(&mut heap, SourceMutPointer::null()),
            Ok(())
        );

        let graph = heap
            .allocate_model_storage(vec![subgraf::default()])
            .unwrap();
        let seen = heap.allocate(vec![1_i32, 2, 3]).unwrap();
        let pathfinder = heap
            .allocate(vec![subgraf_pathfinder {
                sg: graph,
                seen,
                ..subgraf_pathfinder::default()
            }])
            .unwrap();
        assert_eq!(subgraf_pathfinder_free(&mut heap, pathfinder), Ok(()));
        assert_eq!(
            heap.slice(seen.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(pathfinder.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(graph.as_const()).unwrap().len(), 1);

        let without_seen = heap
            .allocate(vec![subgraf_pathfinder {
                sg: graph,
                ..subgraf_pathfinder::default()
            }])
            .unwrap();
        assert_eq!(subgraf_pathfinder_free(&mut heap, without_seen), Ok(()));
        assert_eq!(
            heap.slice(without_seen.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__strutil__add_bond_if_unseen__line_5380() {
        let mut heap = SourceHeap::default();
        let nodes = heap.allocate(vec![i32::MIN, 20_i32, i32::MAX]).unwrap();
        let graph = heap
            .allocate_model_storage(vec![subgraf {
                nnodes: 3,
                nodes,
                ..subgraf::default()
            }])
            .unwrap();
        let pathfinder = heap
            .allocate_model_storage(vec![subgraf_pathfinder {
                sg: graph,
                ..subgraf_pathfinder::default()
            }])
            .unwrap();
        let row0 = heap.allocate(vec![20_i32, i32::MIN]).unwrap();
        let row1 = heap.allocate(vec![111_i32, 222]).unwrap();
        let row2 = heap.allocate(vec![333_i32, 444]).unwrap();
        let bonds = heap.allocate(vec![row0, row1, row2]).unwrap();

        let mut count = 1_i32;
        add_bond_if_unseen(&mut heap, pathfinder, 0, 1, &mut count, bonds).unwrap();
        assert_eq!(count, 1);
        assert_eq!(heap.slice(row1.as_const()).unwrap(), &[111, 222]);

        add_bond_if_unseen(&mut heap, pathfinder, 1, 2, &mut count, bonds).unwrap();
        assert_eq!(count, 2);
        assert_eq!(heap.slice(row1.as_const()).unwrap(), &[20, i32::MAX]);

        add_bond_if_unseen(&mut heap, pathfinder, 2, 0, &mut count, bonds).unwrap();
        assert_eq!(count, 3);
        assert_eq!(heap.slice(row2.as_const()).unwrap(), &[i32::MAX, i32::MIN]);

        let self_row = heap.allocate(vec![9_i32, 8]).unwrap();
        let self_bonds = heap.allocate(vec![self_row]).unwrap();
        let mut self_count = 0_i32;
        add_bond_if_unseen(&mut heap, pathfinder, 1, 1, &mut self_count, self_bonds).unwrap();
        assert_eq!(self_count, 1);
        assert_eq!(heap.slice(self_row.as_const()).unwrap(), &[20, 20]);

        let short_row = heap.allocate(vec![777_i32]).unwrap();
        let short_bonds = heap.allocate(vec![short_row]).unwrap();
        let mut short_count = 0_i32;
        assert_eq!(
            add_bond_if_unseen(&mut heap, pathfinder, 0, 2, &mut short_count, short_bonds,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short_count, 0);
        assert_eq!(heap.slice(short_row.as_const()).unwrap(), &[i32::MIN]);

        let mut negative_count = -1_i32;
        assert_eq!(
            add_bond_if_unseen(&mut heap, pathfinder, 0, 2, &mut negative_count, bonds,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(negative_count, -1);
    }

    #[test]
    fn source_port__strutil__subgraf_pathfinder_run__line_5273() {
        fn diamond(
            heap: &mut SourceHeap,
        ) -> (SourceMutPointer<subgraf_pathfinder>, SourceMutPointer<i32>) {
            let nodes = heap.allocate(vec![101_i32, 102, 103, 104]).unwrap();
            let degrees = heap.allocate(vec![2_i32, 2, 2, 2]).unwrap();
            let rows = [
                vec![
                    subgraf_edge { nbr: 1, etype: 1 },
                    subgraf_edge { nbr: 2, etype: 1 },
                ],
                vec![
                    subgraf_edge { nbr: 0, etype: 1 },
                    subgraf_edge { nbr: 3, etype: 1 },
                ],
                vec![
                    subgraf_edge { nbr: 0, etype: 1 },
                    subgraf_edge { nbr: 3, etype: 1 },
                ],
                vec![
                    subgraf_edge { nbr: 1, etype: 1 },
                    subgraf_edge { nbr: 2, etype: 1 },
                ],
            ];
            let mut row_pointers = Vec::new();
            for row in rows {
                row_pointers.push(heap.allocate(row).unwrap());
            }
            let adjacency = heap.allocate(row_pointers).unwrap();
            let graph = heap
                .allocate_model_storage(vec![subgraf {
                    nnodes: 4,
                    nodes,
                    degrees,
                    adj: adjacency,
                    ..subgraf::default()
                }])
                .unwrap();
            let seen = heap.allocate(vec![0_i32; 4]).unwrap();
            let pathfinder = heap
                .allocate_model_storage(vec![subgraf_pathfinder {
                    sg: graph,
                    start: 0,
                    end: 3,
                    nseen: 1,
                    seen,
                    ..subgraf_pathfinder::default()
                }])
                .unwrap();
            (pathfinder, seen)
        }

        fn outputs(
            heap: &mut SourceHeap,
        ) -> (
            SourceMutPointer<SourceMutPointer<i32>>,
            Vec<SourceMutPointer<i32>>,
            SourceMutPointer<i32>,
        ) {
            let rows = (0..6)
                .map(|_| heap.allocate(vec![-1_i32, -1]).unwrap())
                .collect::<Vec<_>>();
            let bonds = heap.allocate(rows.clone()).unwrap();
            let atoms = heap.allocate(vec![-1_i32; 6]).unwrap();
            (bonds, rows, atoms)
        }

        let mut heap = SourceHeap::default();
        let (pathfinder, seen) = diamond(&mut heap);
        let (bonds, rows, atoms) = outputs(&mut heap);
        let mut bond_count = 0_i32;
        let mut atom_count = 0_i32;
        subgraf_pathfinder_run(
            &mut heap,
            pathfinder,
            0,
            SourceMutPointer::null(),
            &mut bond_count,
            bonds,
            &mut atom_count,
            atoms,
        )
        .unwrap();
        assert_eq!(atom_count, 4);
        assert_eq!(
            heap.slice(atoms.as_const()).unwrap()[..4],
            [101, 102, 104, 103]
        );
        assert_eq!(bond_count, 4);
        let found_bonds = rows[..4]
            .iter()
            .map(|row| heap.slice(row.as_const()).unwrap().to_vec())
            .collect::<Vec<_>>();
        assert_eq!(
            found_bonds,
            vec![
                vec![101, 102],
                vec![102, 104],
                vec![101, 103],
                vec![103, 104]
            ]
        );
        assert_eq!(heap.slice(seen.as_const()).unwrap(), &[0, 0, 0, 0]);
        assert_eq!(heap.slice(pathfinder.as_const()).unwrap()[0].nseen, 1);

        let (filtered_pathfinder, filtered_seen) = diamond(&mut heap);
        let (filtered_bonds, filtered_rows, filtered_atoms) = outputs(&mut heap);
        let forbidden = heap.allocate(vec![1_i32, 0]).unwrap();
        let mut filtered_bond_count = 0_i32;
        let mut filtered_atom_count = 0_i32;
        subgraf_pathfinder_run(
            &mut heap,
            filtered_pathfinder,
            1,
            forbidden,
            &mut filtered_bond_count,
            filtered_bonds,
            &mut filtered_atom_count,
            filtered_atoms,
        )
        .unwrap();
        assert_eq!(filtered_atom_count, 3);
        assert_eq!(
            heap.slice(filtered_atoms.as_const()).unwrap()[..3],
            [101, 103, 104]
        );
        assert_eq!(filtered_bond_count, 2);
        assert_eq!(
            heap.slice(filtered_rows[0].as_const()).unwrap(),
            &[101, 103]
        );
        assert_eq!(
            heap.slice(filtered_rows[1].as_const()).unwrap(),
            &[103, 104]
        );
        assert_eq!(heap.slice(filtered_seen.as_const()).unwrap(), &[0, 0, 0, 0]);

        let (recursive_pathfinder, _) = diamond(&mut heap);
        let (recursive_bonds, _, recursive_atoms) = outputs(&mut heap);
        let recursive_forbidden = heap.allocate(vec![1_i32, 3]).unwrap();
        let mut recursive_bond_count = 0_i32;
        let mut recursive_atom_count = 0_i32;
        subgraf_pathfinder_run(
            &mut heap,
            recursive_pathfinder,
            1,
            recursive_forbidden,
            &mut recursive_bond_count,
            recursive_bonds,
            &mut recursive_atom_count,
            recursive_atoms,
        )
        .unwrap();
        assert_eq!(recursive_atom_count, 4);
        assert_eq!(recursive_bond_count, 4);

        let (no_atoms_pathfinder, no_atoms_seen) = diamond(&mut heap);
        let (no_atoms_bonds, _, _) = outputs(&mut heap);
        let mut no_atoms_bond_count = 0_i32;
        let mut untouched_atom_count = 7_i32;
        subgraf_pathfinder_run(
            &mut heap,
            no_atoms_pathfinder,
            0,
            SourceMutPointer::null(),
            &mut no_atoms_bond_count,
            no_atoms_bonds,
            &mut untouched_atom_count,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_eq!(no_atoms_bond_count, 4);
        assert_eq!(untouched_atom_count, 7);
        assert_eq!(heap.slice(no_atoms_seen.as_const()).unwrap(), &[0, 0, 0, 0]);

        let empty_seen = heap.allocate(vec![9_i32]).unwrap();
        let empty_graph = heap
            .allocate_model_storage(vec![subgraf::default()])
            .unwrap();
        let empty_pathfinder = heap
            .allocate_model_storage(vec![subgraf_pathfinder {
                sg: empty_graph,
                nseen: 0,
                seen: empty_seen,
                ..subgraf_pathfinder::default()
            }])
            .unwrap();
        let mut untouched_bonds = 11_i32;
        let mut untouched_atoms = 12_i32;
        subgraf_pathfinder_run(
            &mut heap,
            empty_pathfinder,
            1,
            SourceMutPointer::null(),
            &mut untouched_bonds,
            SourceMutPointer::null(),
            &mut untouched_atoms,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_eq!((untouched_bonds, untouched_atoms), (11, 12));
        assert_eq!(heap.slice(empty_seen.as_const()).unwrap(), &[9]);
    }

    #[test]
    fn source_port__strutil__subgraf_pathfinder_collect_all__line_5422() {
        fn diamond(
            heap: &mut SourceHeap,
            initial_seen: &[i32],
        ) -> (SourceMutPointer<subgraf_pathfinder>, SourceMutPointer<i32>) {
            let nodes = heap.allocate(vec![101_i32, 102, 103, 104]).unwrap();
            let degrees = heap.allocate(vec![2_i32, 2, 2, 2]).unwrap();
            let rows = [
                vec![
                    subgraf_edge { nbr: 1, etype: 1 },
                    subgraf_edge { nbr: 2, etype: 1 },
                ],
                vec![
                    subgraf_edge { nbr: 0, etype: 1 },
                    subgraf_edge { nbr: 3, etype: 1 },
                ],
                vec![
                    subgraf_edge { nbr: 0, etype: 1 },
                    subgraf_edge { nbr: 3, etype: 1 },
                ],
                vec![
                    subgraf_edge { nbr: 1, etype: 1 },
                    subgraf_edge { nbr: 2, etype: 1 },
                ],
            ];
            let row_pointers = rows
                .into_iter()
                .map(|row| heap.allocate(row).unwrap())
                .collect::<Vec<_>>();
            let adjacency = heap.allocate(row_pointers).unwrap();
            let graph = heap
                .allocate_model_storage(vec![subgraf {
                    nnodes: 4,
                    nodes,
                    degrees,
                    adj: adjacency,
                    ..subgraf::default()
                }])
                .unwrap();
            let mut seen_values = vec![-9_i32; 6];
            seen_values[..initial_seen.len()].copy_from_slice(initial_seen);
            let seen = heap.allocate(seen_values).unwrap();
            let pathfinder = heap
                .allocate_model_storage(vec![subgraf_pathfinder {
                    sg: graph,
                    start: 0,
                    nseen: initial_seen.len() as i32,
                    seen,
                    ..subgraf_pathfinder::default()
                }])
                .unwrap();
            (pathfinder, seen)
        }

        let mut heap = SourceHeap::default();
        let (pathfinder, seen) = diamond(&mut heap, &[]);
        let atom_numbers = heap.allocate(vec![-7_i32; 6]).unwrap();
        let allocation_count = heap.live_allocation_count();
        assert_eq!(
            subgraf_pathfinder_collect_all(
                &mut heap,
                pathfinder,
                0,
                SourceMutPointer::null(),
                atom_numbers,
            ),
            Ok(4)
        );
        assert_eq!(heap.live_allocation_count(), allocation_count);
        assert_eq!(heap.slice(seen.as_const()).unwrap(), &[0, 1, 3, 2, -9, -9]);
        assert_eq!(
            heap.slice(atom_numbers.as_const()).unwrap(),
            &[101, 102, 104, 103, -7, -7]
        );
        let state = &heap.slice(pathfinder.as_const()).unwrap()[0];
        assert_eq!((state.start, state.nseen), (2, 4));

        let (filtered_pathfinder, filtered_seen) = diamond(&mut heap, &[]);
        let filtered_atoms = heap.allocate(vec![55_i32; 6]).unwrap();
        let forbidden = heap.allocate(vec![77_i32, 88, 3, 1, 0, 2]).unwrap();
        assert_eq!(
            subgraf_pathfinder_collect_all(
                &mut heap,
                filtered_pathfinder,
                3,
                forbidden,
                filtered_atoms,
            ),
            Ok(2)
        );
        assert_eq!(
            heap.slice(filtered_seen.as_const()).unwrap(),
            &[0, 1, -9, -9, -9, -9]
        );
        assert_eq!(
            heap.slice(filtered_atoms.as_const()).unwrap(),
            &[101, 102, 55, 55, 55, 55]
        );
        let state = &heap.slice(filtered_pathfinder.as_const()).unwrap()[0];
        assert_eq!((state.start, state.nseen), (1, 2));

        for (forbidden_count, forbidden_pointer) in [(1, SourceMutPointer::null()), (-1, forbidden)]
        {
            let (unfiltered_pathfinder, unfiltered_seen) = diamond(&mut heap, &[]);
            let unfiltered_atoms = heap.allocate(vec![-1_i32; 6]).unwrap();
            assert_eq!(
                subgraf_pathfinder_collect_all(
                    &mut heap,
                    unfiltered_pathfinder,
                    forbidden_count,
                    forbidden_pointer,
                    unfiltered_atoms,
                ),
                Ok(4)
            );
            assert_eq!(
                heap.slice(unfiltered_seen.as_const()).unwrap()[..4],
                [0, 1, 3, 2]
            );
        }

        let (preseen_pathfinder, preseen) = diamond(&mut heap, &[3]);
        let preseen_atoms = heap.allocate(vec![999_i32; 6]).unwrap();
        assert_eq!(
            subgraf_pathfinder_collect_all(
                &mut heap,
                preseen_pathfinder,
                0,
                SourceMutPointer::null(),
                preseen_atoms,
            ),
            Ok(4)
        );
        assert_eq!(
            heap.slice(preseen.as_const()).unwrap(),
            &[3, 0, 1, 2, -9, -9]
        );
        assert_eq!(
            heap.slice(preseen_atoms.as_const()).unwrap(),
            &[999, 101, 102, 103, 999, 999]
        );
        let state = &heap.slice(preseen_pathfinder.as_const()).unwrap()[0];
        assert_eq!((state.start, state.nseen), (2, 4));

        let isolated_nodes = heap.allocate(vec![i32::MIN]).unwrap();
        let isolated_degrees = heap.allocate(vec![-1_i32]).unwrap();
        let isolated_graph = heap
            .allocate_model_storage(vec![subgraf {
                nnodes: 1,
                nodes: isolated_nodes,
                degrees: isolated_degrees,
                ..subgraf::default()
            }])
            .unwrap();
        let isolated_seen = heap.allocate(vec![71_i32, 72]).unwrap();
        let isolated_pathfinder = heap
            .allocate_model_storage(vec![subgraf_pathfinder {
                sg: isolated_graph,
                start: 0,
                seen: isolated_seen,
                ..subgraf_pathfinder::default()
            }])
            .unwrap();
        let isolated_atoms = heap.allocate(vec![81_i32, 82]).unwrap();
        assert_eq!(
            subgraf_pathfinder_collect_all(
                &mut heap,
                isolated_pathfinder,
                0,
                SourceMutPointer::null(),
                isolated_atoms,
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(isolated_seen.as_const()).unwrap(), &[0, 72]);
        assert_eq!(
            heap.slice(isolated_atoms.as_const()).unwrap(),
            &[i32::MIN, 82]
        );

        let invalid_seen = heap.allocate(vec![17_i32]).unwrap();
        let invalid_pathfinder = heap
            .allocate_model_storage(vec![subgraf_pathfinder {
                sg: isolated_graph,
                start: -1,
                seen: invalid_seen,
                ..subgraf_pathfinder::default()
            }])
            .unwrap();
        assert_eq!(
            subgraf_pathfinder_collect_all(
                &mut heap,
                invalid_pathfinder,
                0,
                SourceMutPointer::null(),
                isolated_atoms,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(invalid_seen.as_const()).unwrap(), &[17]);
        assert_eq!(
            heap.slice(invalid_pathfinder.as_const()).unwrap()[0].nseen,
            0
        );
    }

    #[test]
    fn source_port__strutil__compatomdata_getnummapping__line_4986() {
        let mut heap = SourceHeap::default();
        let atoms = [3_u16, 1, 3, 0, u16::MAX]
            .into_iter()
            .map(|orig_at_number| inp_ATOM {
                orig_at_number,
                ..inp_ATOM::default()
            })
            .collect::<Vec<_>>();
        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let data = COMP_ATOM_DATA {
            at: atom_pointer,
            num_at: 5,
            ..COMP_ATOM_DATA::default()
        };
        let original = heap.allocate(vec![-7_i32; 6]).unwrap();
        let current = heap
            .allocate(vec![-9_i32; usize::from(u16::MAX) + 1])
            .unwrap();
        CompAtomData_GetNumMapping(&mut heap, &data, original, current).unwrap();
        assert_eq!(
            heap.slice(original.as_const()).unwrap(),
            &[3, 1, 3, 0, i32::from(u16::MAX), -7]
        );
        let current_values = heap.slice(current.as_const()).unwrap();
        assert_eq!(current_values[0], 3);
        assert_eq!(current_values[1], 1);
        assert_eq!(current_values[3], 2);
        assert_eq!(current_values[usize::from(u16::MAX)], 4);
        assert_eq!(current_values[2], -9);

        let inaccessible = COMP_ATOM_DATA {
            at: SourceMutPointer::null(),
            num_at: i32::MAX,
            ..COMP_ATOM_DATA::default()
        };
        let untouched = heap.allocate(vec![11_i32, 12]).unwrap();
        CompAtomData_GetNumMapping(
            &mut heap,
            &inaccessible,
            SourceMutPointer::null(),
            untouched,
        )
        .unwrap();
        assert_eq!(heap.slice(untouched.as_const()).unwrap(), &[11, 12]);
        CompAtomData_GetNumMapping(
            &mut heap,
            &inaccessible,
            untouched,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_eq!(heap.slice(untouched.as_const()).unwrap(), &[11, 12]);

        let negative = COMP_ATOM_DATA {
            at: SourceMutPointer::null(),
            num_at: i32::MIN,
            ..COMP_ATOM_DATA::default()
        };
        let empty_original = heap.allocate(Vec::<i32>::new()).unwrap();
        let empty_current = heap.allocate(Vec::<i32>::new()).unwrap();
        CompAtomData_GetNumMapping(&mut heap, &negative, empty_original, empty_current).unwrap();

        let one_atom = heap
            .allocate_model_storage(vec![inp_ATOM {
                orig_at_number: 2,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let one = COMP_ATOM_DATA {
            at: one_atom,
            num_at: 1,
            ..COMP_ATOM_DATA::default()
        };
        let partial_original = heap.allocate(vec![-1_i32]).unwrap();
        let short_current = heap.allocate(vec![-2_i32; 2]).unwrap();
        assert_eq!(
            CompAtomData_GetNumMapping(&mut heap, &one, partial_original, short_current),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial_original.as_const()).unwrap(), &[2]);
        assert_eq!(heap.slice(short_current.as_const()).unwrap(), &[-2, -2]);
    }
}
