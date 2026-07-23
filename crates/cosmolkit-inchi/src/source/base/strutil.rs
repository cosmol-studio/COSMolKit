use crate::source::base::ichisort::insertions_sort_AT_RANK;
use crate::source::base::ichister::get_opposite_sb_atom_slice;
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    AB_MAX_WELL_DEFINED_PARITY, AB_MIN_WELL_DEFINED_PARITY, AB_PARITY_UNDF, AT_NUMB, AT_RANK,
    ATW_H, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, CT_OUT_OF_RAM, INChI, INChI_Aux, INChI_Stereo,
    MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS, MAXVAL, NUM_H_ISOTOPES, ORIG_ATOM_DATA,
    S_CHAR, SourceHeap, SourceHeapError, SourceMutPointer, inp_ATOM,
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

    let isotope_order = i32::from(first.iso_atw_diff) - i32::from(second.iso_atw_diff);
    if isotope_order != 0 {
        isotope_order
    } else {
        i32::from(first.component) - i32::from(second.component)
    }
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
        heap.slice_mut(component_rows)?[..component_count].sort_by(|first, second| {
            let size_order = i32::from(second[0]) - i32::from(first[0]);
            if size_order != 0 {
                size_order.cmp(&0)
            } else {
                (i32::from(first[1]) - i32::from(second[1])).cmp(&0)
            }
        });
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
    use crate::source_types::{INChI_IsotopicAtom, INChI_IsotopicTGroup};

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
}
