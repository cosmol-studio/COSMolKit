use crate::source::base::mol2atom::{CreateInpAtom, FreeCompAtomData};
use crate::source::base::util::inchi_calloc;
use crate::source_types::{
    AT_NUMB, COMP_ATOM_DATA, INP_ATOM_DATA2, SourceHeap, SourceHeapError, SourceMutPointer,
    TAUT_INI, TAUT_NUM, TAUT_YES,
};

const SOURCE_SIZEOF_AT_NUMB: u64 = 2;

#[allow(non_snake_case)]
pub(crate) fn CreateCompAtomData(
    heap: &mut SourceHeap,
    input_atom_data: &mut COMP_ATOM_DATA,
    num_atoms: i32,
    num_components: i32,
    intermediate_tautomer: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2232 CreateCompAtomData
    // INCHI✔️❌: int CreateCompAtomData( COMP_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                         int num_atoms,
    // INCHI✔️❌:                         int num_components,
    // INCHI✔️❌:                         int bIntermediateTaut )
    // INCHI✔️❌: {
    // INCHI✔️❌:     FreeCompAtomData( inp_at_data );
    // INCHI✔️❌:     if (( inp_at_data->at = CreateInpAtom( num_atoms ) ) &&
    // INCHI✔️❌:         ( num_components <= 1 || bIntermediateTaut ||
    // INCHI✔️❌:         ( inp_at_data->nOffsetAtAndH = (AT_NUMB*) inchi_calloc( 2 * ( (long long)num_components + 1 ), sizeof( inp_at_data->nOffsetAtAndH[0] ) ) ) )) /* djb-rwth: cast operator added */
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         inp_at_data->num_at = num_atoms;
    // INCHI✔️❌:         inp_at_data->num_components = ( num_components > 1 ) ? num_components : 0;
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     FreeCompAtomData( inp_at_data );
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateCompAtomData
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompAtomData
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(long long) == 8; sizeof(AT_NUMB) == 2
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompAtomData

    FreeCompAtomData(heap, input_atom_data)?;
    input_atom_data.at = CreateInpAtom(heap, num_atoms)?;
    if !input_atom_data.at.is_null() {
        if num_components <= 1 || intermediate_tautomer != 0 {
            input_atom_data.num_at = num_atoms;
            input_atom_data.num_components = if num_components > 1 {
                num_components
            } else {
                0
            };
            return Ok(1);
        }

        let offset_count = 2_u64
            .checked_mul(
                u64::try_from(i64::from(num_components) + 1)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
            )
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        input_atom_data.nOffsetAtAndH =
            match inchi_calloc::<AT_NUMB>(heap, offset_count, SOURCE_SIZEOF_AT_NUMB) {
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                result => result?,
            };
        if !input_atom_data.nOffsetAtAndH.is_null() {
            input_atom_data.num_at = num_atoms;
            input_atom_data.num_components = num_components;
            return Ok(1);
        }
    }

    FreeCompAtomData(heap, input_atom_data)?;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CreateCompositeNormAtom(
    heap: &mut SourceHeap,
    composite_norm_data: &mut [COMP_ATOM_DATA; TAUT_NUM as usize + 1],
    all_inp_norm_data: &[INP_ATOM_DATA2],
    num_components: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1973 CreateCompositeNormAtom
    // INCHI✔️❌: int CreateCompositeNormAtom( COMP_ATOM_DATA *composite_norm_data,
    // INCHI✔️❌:                              INP_ATOM_DATA2 *all_inp_norm_data,
    // INCHI✔️❌:                              int num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, jj, k, n, m, tot_num_at, tot_num_H, cur_num_at, cur_num_H; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int num_comp[TAUT_NUM + 1], num_taut[TAUT_NUM + 1], num_del[TAUT_NUM + 1], num_at[TAUT_NUM + 1], num_inp_at[TAUT_NUM + 1];
    // INCHI✔️❌:     int ret = 0, indicator = 1;
    // INCHI✔️❌:     inp_ATOM *at, *at_from;
    // INCHI✔️❌:     memset( num_comp, 0, sizeof( num_comp ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( num_taut, 0, sizeof( num_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( num_del, 0, sizeof( num_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     /* count taut and non-taut components */
    // INCHI✔️❌:     for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_comp[j] = num_taut[j] = 0;
    // INCHI✔️❌:         for (i = 0; i < num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_del[j] += ( 0 != all_inp_norm_data[i][j].bDeleted );
    // INCHI✔️❌:                 num_comp[j] ++;
    // INCHI✔️❌:                 num_taut[j] += ( 0 != all_inp_norm_data[i][j].bTautomeric );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count intermediate taut structure components */
    // INCHI✔️❌:     if (num_comp[TAUT_YES] > num_del[TAUT_YES] && num_taut[TAUT_YES])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         num_comp[TAUT_INI] = num_comp[TAUT_YES] - num_del[TAUT_YES];
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0, j = TAUT_YES; i < num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (all_inp_norm_data[i][j].bExists &&
    // INCHI✔️❌:                 ( all_inp_norm_data[i][j].bDeleted ||
    // INCHI✔️❌:                     (all_inp_norm_data[i][j].bTautomeric &&
    // INCHI✔️❌:                     all_inp_norm_data[i][j].at_fixed_bonds &&
    // INCHI✔️❌:                     all_inp_norm_data[i][j].bTautPreprocessed) )) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_comp[TAUT_INI] ++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count atoms and allocate composite atom data */
    // INCHI✔️❌:     for (jj = 0; jj <= TAUT_INI; jj++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_at[jj] = num_inp_at[jj] = 0;
    // INCHI✔️❌:         j = inchi_min( jj, TAUT_YES );
    // INCHI✔️❌:         if (num_comp[jj])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_components; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (all_inp_norm_data[i][j].bDeleted)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* find k = the normaized structure index */
    // INCHI✔️❌:                 if (jj == TAUT_INI)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists &&
    // INCHI✔️❌:                          all_inp_norm_data[i][j].at_fixed_bonds)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted &&
    // INCHI✔️❌:                                 !all_inp_norm_data[i][j].bDeleted)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 k = j;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 num_inp_at[jj] += all_inp_norm_data[i][k].num_at; /* all atoms including terminal H */
    // INCHI✔️❌:                 num_at[jj] += all_inp_norm_data[i][k].num_at - all_inp_norm_data[i][k].num_removed_H;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num_inp_at[jj])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!CreateCompAtomData( composite_norm_data + jj, num_inp_at[jj], num_components, jj == TAUT_INI ))
    // INCHI✔️❌:                     goto exit_error;
    // INCHI✔️❌:                 composite_norm_data[jj].num_removed_H = num_inp_at[jj] - num_at[jj];
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* fill out composite atom */
    // INCHI✔️❌:     for (jj = 0; jj <= TAUT_INI; jj++, indicator <<= 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = inchi_min( jj, TAUT_YES );
    // INCHI✔️❌:         if (num_comp[jj])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             tot_num_at = 0;
    // INCHI✔️❌:             tot_num_H = 0;
    // INCHI✔️❌:             for (i = 0; i < num_components; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (all_inp_norm_data[i][j].bDeleted)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][j].nNumRemovedProtons;
    // INCHI✔️❌:                     for (n = 0; n < NUM_H_ISOTOPES; n++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][j].nNumRemovedProtonsIsotopic[n];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                 k = TAUT_NUM; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔️❌:                 /* find k = the normaized structure index */
    // INCHI✔️❌:                 if (jj == TAUT_INI)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists && all_inp_norm_data[i][j].at_fixed_bonds)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (all_inp_norm_data[i][j].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 k = j;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* copy main atoms */
    // INCHI✔️❌:                 cur_num_H = all_inp_norm_data[i][k].num_removed_H;       /* number of terminal H atoms */
    // INCHI✔️❌:                 cur_num_at = all_inp_norm_data[i][k].num_at - cur_num_H;  /* number of all but explicit terminal H atoms */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (( tot_num_at + cur_num_at ) > num_at[jj] ||
    // INCHI✔️❌:                     ( num_at[jj] + tot_num_H + cur_num_H ) > num_inp_at[jj])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_error; /* miscount */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 at = composite_norm_data[jj].at + tot_num_at; /* points to the 1st destination atom */
    // INCHI✔️❌:                 at_from = ( jj == TAUT_INI && k == TAUT_YES && all_inp_norm_data[i][k].at_fixed_bonds ) ?
    // INCHI✔️❌:                     all_inp_norm_data[i][k].at_fixed_bonds : all_inp_norm_data[i][k].at;
    // INCHI✔️❌:                 memcpy( at, at_from, sizeof( composite_norm_data[0].at[0] ) * cur_num_at ); /* copy atoms except terminal H */
    // INCHI✔️❌:                 /* shift neighbors of main atoms */
    // INCHI✔️❌:                 for (n = 0; n < cur_num_at; n++, at++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (m = 0; m < at->valence; m++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at->neighbor[m] += tot_num_at;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* copy explicit H */
    // INCHI✔️❌:                 if (cur_num_H)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at = composite_norm_data[jj].at + num_at[jj] + tot_num_H; /* points to the 1st destination atom */
    // INCHI✔️❌:                     memcpy( at, at_from + cur_num_at, sizeof( composite_norm_data[0].at[0] ) * cur_num_H );
    // INCHI✔️❌:                     /* shift neighbors of explicit H atoms */
    // INCHI✔️❌:                     for (n = 0; n < cur_num_H; n++, at++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         for (m = 0; m < at->valence; m++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             at->neighbor[m] += tot_num_at;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* composite counts */
    // INCHI✔️❌:                 composite_norm_data[jj].bHasIsotopicLayer |= all_inp_norm_data[i][k].bHasIsotopicLayer;
    // INCHI✔️❌:                 composite_norm_data[jj].num_isotopic += all_inp_norm_data[i][k].num_isotopic;
    // INCHI✔️❌:                 composite_norm_data[jj].num_bonds += all_inp_norm_data[i][k].num_bonds;
    // INCHI✔️❌:                 composite_norm_data[jj].bTautomeric += ( j == jj ) && all_inp_norm_data[i][k].bTautomeric;
    // INCHI✔️❌:                 composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][k].nNumRemovedProtons;
    // INCHI✔️❌:                 for (n = 0; n < NUM_H_ISOTOPES; n++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][k].nNumRemovedProtonsIsotopic[n];
    // INCHI✔️❌:                     composite_norm_data[jj].num_iso_H[n] += all_inp_norm_data[i][k].num_iso_H[n];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                 composite_norm_data[j].num_at            += cur_num_at + cur_num_H;
    // INCHI✔️❌:                 composite_norm_data[j].num_removed_H     += cur_num_H;
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 /* total count */
    // INCHI✔️❌:                 tot_num_at += cur_num_at;
    // INCHI✔️❌:                 tot_num_H += cur_num_H;
    // INCHI✔️❌:                 /* offset for the next component */
    // INCHI✔️❌:                 if (composite_norm_data[jj].nOffsetAtAndH)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     composite_norm_data[jj].nOffsetAtAndH[2 * i] = tot_num_at;
    // INCHI✔️❌:                     composite_norm_data[jj].nOffsetAtAndH[2 * i + 1] = num_at[jj] + tot_num_H;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (tot_num_at != num_at[jj] ||
    // INCHI✔️❌:                  num_at[jj] + tot_num_H != num_inp_at[jj])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error; /* miscount */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             composite_norm_data[jj].bExists = ( tot_num_at > 0 );
    // INCHI✔️❌:             ret |= indicator;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateCompositeNormAtom
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompositeNormAtom
    // INCHI✔️❌: #define ALT_TAUT(X) ((X)>TAUT_YES? TAUT_YES : 1-(X))
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; TAUT_NUM == 2; TAUT_YES == 1; TAUT_INI == 2
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompositeNormAtom

    let component_count = if num_components > 0 {
        usize::try_from(num_components).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    if component_count > all_inp_norm_data.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let taut_num = TAUT_NUM as usize;
    let taut_yes = TAUT_YES as usize;
    let taut_ini = TAUT_INI as usize;
    let mut num_comp = [0_i32; TAUT_NUM as usize + 1];
    let mut num_taut = [0_i32; TAUT_NUM as usize + 1];
    let mut num_del = [0_i32; TAUT_NUM as usize + 1];
    let mut num_at = [0_i32; TAUT_NUM as usize + 1];
    let mut num_inp_at = [0_i32; TAUT_NUM as usize + 1];
    let mut ret = 0_i32;

    for j in 0..taut_num {
        num_comp[j] = 0;
        num_taut[j] = 0;
        for component in &all_inp_norm_data[..component_count] {
            if component[j].bExists != 0 {
                num_del[j] = num_del[j].wrapping_add(i32::from(component[j].bDeleted != 0));
                num_comp[j] = num_comp[j].wrapping_add(1);
                num_taut[j] = num_taut[j].wrapping_add(i32::from(component[j].bTautomeric != 0));
            }
        }
    }

    if num_comp[taut_yes] > num_del[taut_yes] && num_taut[taut_yes] != 0 {
        let j = taut_yes;
        for component in &all_inp_norm_data[..component_count] {
            if component[j].bExists != 0
                && (component[j].bDeleted != 0
                    || (component[j].bTautomeric != 0
                        && !component[j].at_fixed_bonds.is_null()
                        && component[j].bTautPreprocessed != 0))
            {
                num_comp[taut_ini] = num_comp[taut_ini].wrapping_add(1);
            }
        }
    }

    for jj in 0..=taut_ini {
        num_at[jj] = 0;
        num_inp_at[jj] = 0;
        let j = jj.min(taut_yes);
        let alternate = if j > taut_yes { taut_yes } else { 1 - j };
        if num_comp[jj] != 0 {
            for component in &all_inp_norm_data[..component_count] {
                if component[j].bDeleted != 0 {
                    continue;
                }
                let k;
                if jj == taut_ini {
                    if component[j].bExists != 0 && !component[j].at_fixed_bonds.is_null() {
                        k = j;
                    } else if component[alternate].bExists != 0
                        && component[alternate].bDeleted == 0
                        && component[j].bDeleted == 0
                    {
                        k = alternate;
                    } else if component[j].bExists != 0 {
                        k = j;
                    } else {
                        continue;
                    }
                } else if component[j].bExists != 0 {
                    k = j;
                } else if component[alternate].bExists != 0 && component[alternate].bDeleted == 0 {
                    k = alternate;
                } else {
                    continue;
                }
                num_inp_at[jj] = num_inp_at[jj].wrapping_add(component[k].num_at);
                num_at[jj] = num_at[jj]
                    .wrapping_add(component[k].num_at.wrapping_sub(component[k].num_removed_H));
            }
            if num_inp_at[jj] != 0 {
                if CreateCompAtomData(
                    heap,
                    &mut composite_norm_data[jj],
                    num_inp_at[jj],
                    num_components,
                    i32::from(jj == taut_ini),
                )? == 0
                {
                    return Ok(ret);
                }
                composite_norm_data[jj].num_removed_H = num_inp_at[jj].wrapping_sub(num_at[jj]);
            }
        }
    }

    let mut indicator = 1_i32;
    for jj in 0..=taut_ini {
        let j = jj.min(taut_yes);
        let alternate = if j > taut_yes { taut_yes } else { 1 - j };
        if num_comp[jj] != 0 {
            let mut total_atoms = 0_i32;
            let mut total_hydrogens = 0_i32;
            for (i, component) in all_inp_norm_data[..component_count].iter().enumerate() {
                if component[j].bDeleted != 0 {
                    composite_norm_data[jj].nNumRemovedProtons = composite_norm_data[jj]
                        .nNumRemovedProtons
                        .wrapping_add(component[j].nNumRemovedProtons);
                    for n in 0..composite_norm_data[jj].nNumRemovedProtonsIsotopic.len() {
                        composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] =
                            composite_norm_data[jj].nNumRemovedProtonsIsotopic[n]
                                .wrapping_add(component[j].nNumRemovedProtonsIsotopic[n]);
                    }
                    continue;
                }

                let k;
                if jj == taut_ini {
                    if component[j].bExists != 0 && !component[j].at_fixed_bonds.is_null() {
                        k = j;
                    } else if component[alternate].bExists != 0 {
                        k = alternate;
                    } else if component[j].bExists != 0 && component[alternate].bDeleted == 0 {
                        k = j;
                    } else {
                        continue;
                    }
                } else if component[j].bExists != 0 {
                    k = j;
                } else if component[alternate].bExists != 0 && component[alternate].bDeleted == 0 {
                    k = alternate;
                } else {
                    continue;
                }

                let current_hydrogens = component[k].num_removed_H;
                let current_atoms = component[k].num_at.wrapping_sub(current_hydrogens);
                if total_atoms.wrapping_add(current_atoms) > num_at[jj]
                    || num_at[jj]
                        .wrapping_add(total_hydrogens)
                        .wrapping_add(current_hydrogens)
                        > num_inp_at[jj]
                {
                    return Ok(ret);
                }

                let source_pointer =
                    if jj == taut_ini && k == taut_yes && !component[k].at_fixed_bonds.is_null() {
                        component[k].at_fixed_bonds
                    } else {
                        component[k].at
                    };
                let source_count = usize::try_from(component[k].num_at)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let source_atoms = heap
                    .slice(source_pointer.as_const())?
                    .get(..source_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let current_atom_count = usize::try_from(current_atoms)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let current_hydrogen_count = usize::try_from(current_hydrogens)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let main_start = usize::try_from(total_atoms)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let hydrogen_start = usize::try_from(num_at[jj].wrapping_add(total_hydrogens))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;

                {
                    let destination = heap.slice_mut(composite_norm_data[jj].at)?;
                    let main_end = main_start
                        .checked_add(current_atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let main = destination
                        .get_mut(main_start..main_end)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    main.clone_from_slice(&source_atoms[..current_atom_count]);
                    for atom in main {
                        let valence = usize::try_from(atom.valence.max(0))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if valence > atom.neighbor.len() {
                            return Err(SourceHeapError::PointerOutOfBounds);
                        }
                        for neighbor in &mut atom.neighbor[..valence] {
                            *neighbor = i32::from(*neighbor).wrapping_add(total_atoms) as AT_NUMB;
                        }
                    }

                    if current_hydrogens != 0 {
                        let hydrogen_end = hydrogen_start
                            .checked_add(current_hydrogen_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let hydrogens = destination
                            .get_mut(hydrogen_start..hydrogen_end)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        hydrogens.clone_from_slice(
                            &source_atoms
                                [current_atom_count..current_atom_count + current_hydrogen_count],
                        );
                        for atom in hydrogens {
                            let valence = usize::try_from(atom.valence.max(0))
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            if valence > atom.neighbor.len() {
                                return Err(SourceHeapError::PointerOutOfBounds);
                            }
                            for neighbor in &mut atom.neighbor[..valence] {
                                *neighbor =
                                    i32::from(*neighbor).wrapping_add(total_atoms) as AT_NUMB;
                            }
                        }
                    }
                }

                composite_norm_data[jj].bHasIsotopicLayer |= component[k].bHasIsotopicLayer;
                composite_norm_data[jj].num_isotopic = composite_norm_data[jj]
                    .num_isotopic
                    .wrapping_add(component[k].num_isotopic);
                composite_norm_data[jj].num_bonds = composite_norm_data[jj]
                    .num_bonds
                    .wrapping_add(component[k].num_bonds);
                composite_norm_data[jj].bTautomeric = composite_norm_data[jj]
                    .bTautomeric
                    .wrapping_add(i32::from(j == jj && component[k].bTautomeric != 0));
                composite_norm_data[jj].nNumRemovedProtons = composite_norm_data[jj]
                    .nNumRemovedProtons
                    .wrapping_add(component[k].nNumRemovedProtons);
                for n in 0..composite_norm_data[jj].nNumRemovedProtonsIsotopic.len() {
                    composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] = composite_norm_data[jj]
                        .nNumRemovedProtonsIsotopic[n]
                        .wrapping_add(component[k].nNumRemovedProtonsIsotopic[n]);
                    composite_norm_data[jj].num_iso_H[n] = composite_norm_data[jj].num_iso_H[n]
                        .wrapping_add(component[k].num_iso_H[n]);
                }
                total_atoms = total_atoms.wrapping_add(current_atoms);
                total_hydrogens = total_hydrogens.wrapping_add(current_hydrogens);
                if !composite_norm_data[jj].nOffsetAtAndH.is_null() {
                    let offsets = heap.slice_mut(composite_norm_data[jj].nOffsetAtAndH)?;
                    offsets[2 * i] = total_atoms as AT_NUMB;
                    offsets[2 * i + 1] = num_at[jj].wrapping_add(total_hydrogens) as AT_NUMB;
                }
            }
            if total_atoms != num_at[jj]
                || num_at[jj].wrapping_add(total_hydrogens) != num_inp_at[jj]
            {
                return Ok(ret);
            }
            composite_norm_data[jj].bExists = i32::from(total_atoms > 0);
            ret |= indicator;
        }
        indicator <<= 1;
    }
    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::mol2atom::FreeCompAtomData;
    use crate::source_types::inp_ATOM;

    fn atom(number: AT_NUMB, neighbor: AT_NUMB) -> inp_ATOM {
        inp_ATOM {
            orig_at_number: number,
            neighbor: {
                let mut neighbors = [0; 20];
                neighbors[0] = neighbor;
                neighbors
            },
            valence: 1,
            ..inp_ATOM::default()
        }
    }

    fn input_data(
        heap: &mut SourceHeap,
        atoms: Vec<inp_ATOM>,
        removed_hydrogens: i32,
    ) -> crate::source_types::INP_ATOM_DATA {
        let atom_count = i32::try_from(atoms.len()).unwrap();
        crate::source_types::INP_ATOM_DATA {
            at: heap.allocate_model_storage(atoms).unwrap(),
            num_at: atom_count,
            num_removed_H: removed_hydrogens,
            bExists: 1,
            ..crate::source_types::INP_ATOM_DATA::default()
        }
    }

    fn free_composites(
        heap: &mut SourceHeap,
        composites: &mut [COMP_ATOM_DATA; TAUT_NUM as usize + 1],
    ) {
        for composite in composites {
            FreeCompAtomData(heap, composite).unwrap();
        }
    }

    #[test]
    fn source_port__inchi_dll_a2__createcompatomdata__line_2232() {
        let mut simple_heap = SourceHeap::default();
        let old_atoms = simple_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let old_offsets = simple_heap
            .allocate_model_storage(vec![91_u16, 92_u16])
            .unwrap();
        let mut simple = COMP_ATOM_DATA {
            at: old_atoms,
            num_at: 77,
            num_removed_H: 76,
            nOffsetAtAndH: old_offsets,
            num_components: 75,
            ..COMP_ATOM_DATA::default()
        };
        assert_eq!(
            CreateCompAtomData(&mut simple_heap, &mut simple, 2, 1, 0),
            Ok(1)
        );
        assert_eq!(simple.num_at, 2);
        assert_eq!(simple.num_components, 0);
        assert_eq!(simple.num_removed_H, 0);
        assert!(!simple.at.is_null());
        assert!(simple.nOffsetAtAndH.is_null());
        assert_eq!(simple_heap.slice(simple.at.as_const()).unwrap().len(), 2);
        assert_eq!(simple_heap.live_allocation_count(), 1);
        FreeCompAtomData(&mut simple_heap, &mut simple).unwrap();

        for (components, intermediate, expected_components) in [(-7, 0, 0), (2, 1, 2), (2, -1, 2)] {
            let mut heap = SourceHeap::default();
            let mut data = COMP_ATOM_DATA::default();
            assert_eq!(
                CreateCompAtomData(&mut heap, &mut data, 1, components, intermediate),
                Ok(1)
            );
            assert_eq!(data.num_at, 1);
            assert_eq!(data.num_components, expected_components);
            assert!(!data.at.is_null());
            assert!(data.nOffsetAtAndH.is_null());
            assert_eq!(heap.live_allocation_count(), 1);
            FreeCompAtomData(&mut heap, &mut data).unwrap();
        }

        let mut composite_heap = SourceHeap::default();
        let mut composite = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut composite_heap, &mut composite, 3, 4, 0),
            Ok(1)
        );
        assert_eq!(composite.num_at, 3);
        assert_eq!(composite.num_components, 4);
        assert_eq!(
            composite_heap
                .slice(composite.nOffsetAtAndH.as_const())
                .unwrap(),
            &[0_u16; 10]
        );
        assert_eq!(composite_heap.live_allocation_count(), 2);
        FreeCompAtomData(&mut composite_heap, &mut composite).unwrap();

        let mut zero_heap = SourceHeap::default();
        let mut zero = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut zero_heap, &mut zero, 0, 0, 0),
            Ok(1)
        );
        assert_eq!(zero.num_at, 0);
        assert!(!zero.at.is_null());
        assert_eq!(zero_heap.slice(zero.at.as_const()).unwrap().len(), 0);
        FreeCompAtomData(&mut zero_heap, &mut zero).unwrap();

        let mut negative_heap = SourceHeap::default();
        let mut negative = COMP_ATOM_DATA {
            num_at: 19,
            ..COMP_ATOM_DATA::default()
        };
        assert_eq!(
            CreateCompAtomData(&mut negative_heap, &mut negative, -1, 2, 0),
            Ok(0)
        );
        assert_eq!(negative, COMP_ATOM_DATA::default());
        assert_eq!(negative_heap.live_allocation_count(), 0);

        let mut first_failure_heap = SourceHeap::default();
        first_failure_heap.fail_after_allocations(0);
        let mut first_failure = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut first_failure_heap, &mut first_failure, 1, 2, 0),
            Ok(0)
        );
        assert_eq!(first_failure, COMP_ATOM_DATA::default());
        assert_eq!(first_failure_heap.source_allocation_calls(), 1);
        assert_eq!(first_failure_heap.live_allocation_count(), 0);

        let mut second_failure_heap = SourceHeap::default();
        second_failure_heap.fail_after_allocations(1);
        let mut second_failure = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut second_failure_heap, &mut second_failure, 1, 2, 0),
            Ok(0)
        );
        assert_eq!(second_failure, COMP_ATOM_DATA::default());
        assert_eq!(second_failure_heap.source_allocation_calls(), 2);
        assert_eq!(second_failure_heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__inchi_dll_a2__createcompositenormatom__line_1973() {
        let mut heap = SourceHeap::default();
        let mut component0_non = input_data(&mut heap, vec![atom(10, 0), atom(11, 0)], 1);
        component0_non.num_bonds = 1;
        component0_non.num_isotopic = 2;
        component0_non.bHasIsotopicLayer = 1;
        component0_non.nNumRemovedProtons = 3;
        component0_non.nNumRemovedProtonsIsotopic = [1, 2, 3];
        component0_non.num_iso_H = [4, 5, 6];

        let mut component0_taut = input_data(&mut heap, vec![atom(20, 0), atom(21, 0)], 1);
        component0_taut.at_fixed_bonds = heap
            .allocate_model_storage(vec![atom(30, 0), atom(31, 0)])
            .unwrap();
        component0_taut.bTautomeric = 1;
        component0_taut.bTautPreprocessed = 1;
        component0_taut.num_bonds = 2;
        component0_taut.num_isotopic = 3;
        component0_taut.bHasIsotopicLayer = 2;
        component0_taut.nNumRemovedProtons = 4;
        component0_taut.nNumRemovedProtonsIsotopic = [2, 3, 4];
        component0_taut.num_iso_H = [5, 6, 7];

        let mut component1_non = input_data(&mut heap, vec![atom(40, 0)], 0);
        component1_non.num_bonds = 3;
        let component1_taut = crate::source_types::INP_ATOM_DATA {
            bExists: 1,
            bDeleted: 1,
            nNumRemovedProtons: 7,
            nNumRemovedProtonsIsotopic: [1, 1, 1],
            ..crate::source_types::INP_ATOM_DATA::default()
        };

        let mut component2_non = input_data(&mut heap, vec![atom(50, 0)], 0);
        component2_non.num_bonds = 4;
        component2_non.num_isotopic = 5;
        component2_non.bHasIsotopicLayer = 4;
        component2_non.nNumRemovedProtons = 6;
        component2_non.nNumRemovedProtonsIsotopic = [3, 4, 5];
        component2_non.num_iso_H = [6, 7, 8];

        let inputs = [
            [component0_non, component0_taut],
            [component1_non, component1_taut],
            [
                component2_non,
                crate::source_types::INP_ATOM_DATA::default(),
            ],
        ];
        let mut composites = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(&mut heap, &mut composites, &inputs, 3),
            Ok(7)
        );

        let non_atoms = heap.slice(composites[0].at.as_const()).unwrap();
        assert_eq!(
            non_atoms
                .iter()
                .map(|entry| entry.orig_at_number)
                .collect::<Vec<_>>(),
            vec![10, 40, 50, 11]
        );
        assert_eq!(non_atoms[1].neighbor[0], 1);
        assert_eq!(non_atoms[2].neighbor[0], 2);
        assert_eq!(non_atoms[3].neighbor[0], 0);
        assert_eq!(composites[0].num_at, 4);
        assert_eq!(composites[0].num_removed_H, 1);
        assert_eq!(composites[0].num_bonds, 8);
        assert_eq!(composites[0].bHasIsotopicLayer, 5);
        assert_eq!(
            heap.slice(composites[0].nOffsetAtAndH.as_const()).unwrap(),
            &[1, 4, 2, 4, 3, 4, 0, 0]
        );

        let taut_atoms = heap.slice(composites[1].at.as_const()).unwrap();
        assert_eq!(
            taut_atoms
                .iter()
                .map(|entry| entry.orig_at_number)
                .collect::<Vec<_>>(),
            vec![20, 50, 21]
        );
        assert_eq!(taut_atoms[1].neighbor[0], 1);
        assert_eq!(composites[1].bTautomeric, 1);
        assert_eq!(composites[1].nNumRemovedProtons, 17);
        assert_eq!(composites[1].nNumRemovedProtonsIsotopic, [6, 8, 10]);
        assert_eq!(composites[1].num_iso_H, [11, 13, 15]);
        assert_eq!(
            heap.slice(composites[1].nOffsetAtAndH.as_const()).unwrap(),
            &[1, 3, 0, 0, 2, 3, 0, 0]
        );

        let intermediate_atoms = heap.slice(composites[2].at.as_const()).unwrap();
        assert_eq!(
            intermediate_atoms
                .iter()
                .map(|entry| entry.orig_at_number)
                .collect::<Vec<_>>(),
            vec![30, 50, 31]
        );
        assert_eq!(intermediate_atoms[1].neighbor[0], 1);
        assert_eq!(composites[2].bTautomeric, 0);
        assert_eq!(composites[2].nNumRemovedProtons, 17);
        assert!(composites[2].nOffsetAtAndH.is_null());
        free_composites(&mut heap, &mut composites);

        let mut no_component_heap = SourceHeap::default();
        let mut no_component = std::array::from_fn(|index| COMP_ATOM_DATA {
            num_at: index as i32 + 70,
            ..COMP_ATOM_DATA::default()
        });
        let before = no_component.clone();
        assert_eq!(
            CreateCompositeNormAtom(&mut no_component_heap, &mut no_component, &[], -1),
            Ok(0)
        );
        assert_eq!(no_component, before);

        let mut first_failure_heap = SourceHeap::default();
        let first_failure_input = [[
            input_data(&mut first_failure_heap, vec![atom(1, 0)], 0),
            crate::source_types::INP_ATOM_DATA::default(),
        ]];
        first_failure_heap.fail_after_allocations(0);
        let mut first_failure = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(
                &mut first_failure_heap,
                &mut first_failure,
                &first_failure_input,
                1,
            ),
            Ok(0)
        );
        assert_eq!(
            first_failure,
            std::array::from_fn(|_| COMP_ATOM_DATA::default())
        );

        let mut later_failure_heap = SourceHeap::default();
        let later_failure_input = [[input_data(&mut later_failure_heap, vec![atom(1, 0)], 0), {
            let mut data = input_data(&mut later_failure_heap, vec![atom(2, 0)], 0);
            data.bTautomeric = 0;
            data
        }]];
        later_failure_heap.fail_after_allocations(1);
        let mut later_failure = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(
                &mut later_failure_heap,
                &mut later_failure,
                &later_failure_input,
                1,
            ),
            Ok(0)
        );
        assert!(!later_failure[0].at.is_null());
        assert!(later_failure[1].at.is_null());
        free_composites(&mut later_failure_heap, &mut later_failure);

        let mut miscount_heap = SourceHeap::default();
        let mut miscount_mobile = input_data(&mut miscount_heap, vec![atom(10, 0), atom(11, 0)], 0);
        miscount_mobile.bTautomeric = 1;
        let mut miscount_deleted_alt = input_data(
            &mut miscount_heap,
            vec![atom(20, 0), atom(21, 0), atom(22, 0), atom(23, 0)],
            0,
        );
        miscount_deleted_alt.bDeleted = 1;
        let trigger_non = crate::source_types::INP_ATOM_DATA::default();
        let mut trigger_mobile = input_data(&mut miscount_heap, vec![atom(30, 0)], 0);
        trigger_mobile.bTautomeric = 1;
        trigger_mobile.bTautPreprocessed = 1;
        trigger_mobile.at_fixed_bonds = miscount_heap
            .allocate_model_storage(vec![atom(31, 0)])
            .unwrap();
        let miscount_inputs = [
            [miscount_deleted_alt, miscount_mobile],
            [trigger_non, trigger_mobile],
        ];
        let mut miscount = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(&mut miscount_heap, &mut miscount, &miscount_inputs, 2),
            Ok(3)
        );
        assert_eq!(miscount[0].bExists, 1);
        assert_eq!(
            miscount_heap.slice(miscount[0].at.as_const()).unwrap()[0].orig_at_number,
            30
        );
        assert_eq!(miscount[1].bExists, 1);
        assert_eq!(miscount[2].bExists, 0);
        free_composites(&mut miscount_heap, &mut miscount);
    }
}
