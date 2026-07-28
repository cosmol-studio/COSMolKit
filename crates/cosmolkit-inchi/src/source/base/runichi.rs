use crate::source::api::inchi_dll_a2::CreateCompositeNormAtom;
use crate::source::base::ichicano::{InchiTimeAddMsec, InchiTimeElapsed, InchiTimeGet};
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichimake::Create_INChI;
use crate::source::base::ichinorm::FreeInpAtomData;
use crate::source::base::ichiprt1::{EditINCHI_HidePolymerZz, OrigStruct_FillOut, OrigStruct_Free};
use crate::source::base::mol_fmt2::MolfileSaveCopy;
use crate::source::base::mol_fmt4::{IntArray_AppendIfAbsent, OrigAtData_WriteToSDfile};
use crate::source::base::mol2atom::{CreateInpAtomData, FreeCompAtomData};
use crate::source::base::runichi2::{
    DoOneStructureEarlyPreprocessing, GetOneComponent, POSEContext_Free, POSEContext_Init,
    TreatErrorsInReadTheStructure, bIsSameBond, eprint_bytes, sdf_label_value,
};
use crate::source::base::runichi3::{
    OAD_Polymer_CyclizeCloseableUnits, OAD_Polymer_DebugTrace, OAD_Polymer_FindBackbones,
    OAD_Polymer_GetRepresentation, OAD_Polymer_SmartReopenCyclizedUnits,
    OAD_ValidatePolymerAndPseudoElementData, OrigAtData_AddBond, OrigAtData_DebugTrace,
    OrigAtData_RemoveBond, PreprocessOneStructure,
};
use crate::source::base::runichi4::{
    GetProcessingWarningsOneComponentInChI, SortAndPrintINChI, TreatCreateINChIWarning,
    TreatErrorsInCreateOneComponentINChI, bIsStructChiral,
};
use crate::source::base::strutil::{
    Alloc_INChI, Alloc_INChI_Aux, SetConnectedComponentNumber, subgraf_free, subgraf_new,
    subgraf_pathfinder_collect_all, subgraf_pathfinder_free, subgraf_pathfinder_new,
};
use crate::source::base::util::{
    extract_auxinfo_substring, extract_inchi_substring, inchi_calloc, inchi_free, is_in_the_ilist,
};
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, _IS_OKAY, CANON_GLOBALS, CMP_COMPONENTS, COMP_ATOM_DATA, CT_OUT_OF_RAM,
    CT_USER_QUIT_ERR, FLAG_INP_AT_CHIRAL, INCHI_BAS, INCHI_CLOCK, INCHI_IOS_STRING,
    INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM, INCHI_MODE, INCHI_NUM, INCHI_OUT_NO_AUX_INFO,
    INCHI_OUT_SAVEOPT, INCHI_OUT_SDFILE_ATOMS_DT, INCHI_OUT_SDFILE_ONLY, INCHI_OUT_SHORT_AUX_INFO,
    INCHI_REC, INP_ATOM_DATA, INP_ATOM_DATA2, INPUT_PARMS, LOG_MASK_ALL, NORM_CANON_FLAGS,
    OAD_StructureEdits, ORIG_ATOM_DATA, ORIG_STRUCT, PINChI_Aux2, PINChI2,
    POLYMER_REPRESENTATION_MIXED, POLYMER_REPRESENTATION_STRUCTURE_BASED, POLYMERS_LEGACY,
    POLYMERS_LEGACY_PLUS, POLYMERS_MODERN, REQ_MODE_BASIC, REQ_MODE_DIFF_UU_STEREO, REQ_MODE_ISO,
    REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO, REQ_MODE_SB_IGN_ALL_UU,
    REQ_MODE_SC_IGN_ALL_UU, REQ_MODE_STEREO, REQ_MODE_TAUT, SAVE_OPT_15T, SAVE_OPT_FIXEDH,
    SAVE_OPT_KET, SAVE_OPT_PT_06_00, SAVE_OPT_PT_13_00, SAVE_OPT_PT_16_00, SAVE_OPT_PT_18_00,
    SAVE_OPT_PT_22_00, SAVE_OPT_PT_39_00, SAVE_OPT_RECMET, SAVE_OPT_SLUUD, SAVE_OPT_SUU,
    STRUCT_DATA, SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, TAUT_NON,
    TAUT_NUM, TAUT_YES, TG_FLAG_1_5_TAUT, TG_FLAG_DISCONNECT_COORD_DONE,
    TG_FLAG_DISCONNECT_SALTS_DONE, TG_FLAG_FIX_ODD_THINGS_DONE, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE,
    TG_FLAG_FOUND_ISOTOPIC_H_DONE, TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_MOVE_CHARGE_COORD_DONE,
    TG_FLAG_MOVE_HPLUS2NEUTR_DONE, TG_FLAG_MOVE_POS_CHARGES_DONE, TG_FLAG_PT_06_00,
    TG_FLAG_PT_13_00, TG_FLAG_PT_16_00, TG_FLAG_PT_18_00, TG_FLAG_PT_22_00, TG_FLAG_PT_39_00,
    TG_FLAG_RECONNECT_COORD, clock_t, inchiTime, inp_ATOM, tagFrameShifScheme_FSS_NONE,
    tagFrameShifScheme_FSS_STARS_CYCLED, tagFrameShifScheme_FSS_STARS_CYCLED_SORTED,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE, tagINCHIBondType_INCHI_BOND_TYPE_NONE,
    tagInputType_INPUT_INCHI, tagInputType_INPUT_MOLFILE, tagInputType_INPUT_SDFILE,
};

fn runichi_c_text(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<Vec<u8>, SourceHeapError> {
    if pointer.is_null() {
        return Ok(Vec::new());
    }
    let bytes = heap.slice(pointer)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..length].iter().map(|byte| *byte as u8).collect())
}

// Two 96-byte `INP_ATOM_DATA` records under the selected GCC/LP64 ABI.
const INP_ATOM_DATA2_GCC_LP64_SIZE: u64 = 192;
const PINCHI2_GCC_LP64_SIZE: u64 = 16;

fn overwrite_source_c_buffer(
    destination: Option<&mut [i8]>,
    bytes: &[u8],
) -> Result<(), SourceHeapError> {
    let destination = destination.ok_or(SourceHeapError::PointerOutOfBounds)?;
    if destination.len() <= bytes.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    for (target, source) in destination.iter_mut().zip(bytes.iter()) {
        *target = *source as i8;
    }
    destination[bytes.len()] = 0;
    Ok(())
}

#[rustfmt::skip]
pub(crate) fn swap_atoms_xyz(
    heap: &mut SourceHeap,
    original_atom_data: Option<&ORIG_ATOM_DATA>,
    first_atom: i32,
    second_atom: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2401 swap_atoms_xyz
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void swap_atoms_xyz( ORIG_ATOM_DATA *orig_at_data, int ia1, int ia2 )
{
    double x, y, z;
    
    if (ia1 != ia2)
    {
        x = orig_at_data->at[ia1].x; 
        y = orig_at_data->at[ia1].y; 
        z = orig_at_data->at[ia1].z;

        orig_at_data->at[ia1].x = orig_at_data->at[ia2].x;
        orig_at_data->at[ia1].y = orig_at_data->at[ia2].y;
        orig_at_data->at[ia1].z = orig_at_data->at[ia2].z;
        
        orig_at_data->at[ia2].x = x;
        orig_at_data->at[ia2].y = y;
        orig_at_data->at[ia2].z = z;
    }

    return;
}
    */
    // END INCHI C FUNCTION: swap_atoms_xyz

    if first_atom == second_atom {
        return Ok(());
    }
    let original_atom_data = original_atom_data.ok_or(SourceHeapError::NullPointer)?;
    let first = usize::try_from(first_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let second = usize::try_from(second_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = heap.slice_mut(original_atom_data.at)?;
    if first >= atoms.len() || second >= atoms.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let x = atoms[first].x;
    let y = atoms[first].y;
    let z = atoms[first].z;

    atoms[first].x = atoms[second].x;
    atoms[first].y = atoms[second].y;
    atoms[first].z = atoms[second].z;

    atoms[second].x = x;
    atoms[second].y = y;
    atoms[second].z = z;
    Ok(())
}

fn set_renumbered_or_delete_from_copy(
    number: &mut [i32],
    copied: &[i32],
    renumbering: &[i32],
    base: i32,
) -> Result<i32, SourceHeapError> {
    let deletion = base
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let mut new_number_of_elements = 0_usize;
    for old_number in copied.iter().copied() {
        let renumbering_index = old_number
            .checked_sub(base)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let renumbering_index =
            usize::try_from(renumbering_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mapped = *renumbering
            .get(renumbering_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let new_number = mapped
            .checked_add(base)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        if new_number == deletion {
            continue;
        }
        *number
            .get_mut(new_number_of_elements)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = new_number;
        new_number_of_elements += 1;
    }
    i32::try_from(new_number_of_elements).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_StructureEdits_Apply(
    heap: &mut SourceHeap,
    _structure_data: &mut STRUCT_DATA,
    _input_parameters: &INPUT_PARMS,
    original_atom_data: &mut ORIG_ATOM_DATA,
    edits: &OAD_StructureEdits,
    return_status: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2427 OAD_StructureEdits_Apply
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_StructureEdits_Apply( STRUCT_DATA *sd, 
                              INPUT_PARMS *ip, 
                              ORIG_ATOM_DATA *orig_at_data, 
                              OAD_StructureEdits *ed, 
                              int *ret)
{
    int ok = 0, fail;
    int i, j, old_a1, old_a2, new_a1, new_a2, n_edits = 0;
    int n_del_atom, n_del_bond, n_new_bond, n_mod_bond, n_mod_coord;
    int a1, a2;
    int bond_type = INCHI_BOND_TYPE_NONE, bond_stereo = INCHI_BOND_STEREO_NONE;
    inp_ATOM *at = orig_at_data->at;
    OAD_Polymer *p = orig_at_data->polymer;
    int *at_renum = NULL;
    int *ibuf = NULL;
    int n_max_stored=-1;

    *ret = _IS_OKAY;

    n_del_atom = ed->del_atom->used;
    n_del_bond = ed->del_bond->used/2;
    n_new_bond = ed->new_bond->used/3;
    n_mod_bond = ed->mod_bond->used/4;
    n_mod_coord = ed->mod_coord->used/2;
    if (n_del_atom + n_del_bond + n_new_bond + n_mod_bond + n_mod_coord < 1)
    {
        return 0;
    }

    n_max_stored = inchi_max(2 * (orig_at_data->num_inp_atoms + 1), 2 * (orig_at_data->num_inp_bonds + 1)); /* set all-purpose buffer */

    ITRACE_("\n***************************\nOrig_at_data BEFORE EDITS:\n");
    OrigAtData_DebugTrace(orig_at_data);
    OAD_Polymer_DebugTrace(orig_at_data->polymer);

    /* Delete bonds */
    if (n_del_bond)
    {
        for (i = 0; i < 2 * n_del_bond; i += 2)
        {
            a1 = ed->del_bond->item[i] - 1;
            a2 = ed->del_bond->item[i + 1] - 1;
            ok = OrigAtData_RemoveBond(a1, a2, at, &bond_type, &bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            n_edits++;
        }
    }
    
    /* Add bonds */
    if (n_new_bond)
    {
        for (i = 0; i < 2 * n_new_bond; i += 2)
        {
            a1 = ed->new_bond->item[i] - 1;
            a2 = ed->new_bond->item[i + 1] - 1;
            /* TODO: consider real bond_type, bond_stereo */
            /* OrigAtData_AddSingleStereolessBond( a1, a2, at, &dummy ); */
            ok = OrigAtData_AddBond(a1, a2, at, bond_type, bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            n_edits++;
        }
    }
    
    /* Modify bonds */
    if (n_mod_bond)
    {
        for (j = 0; j < 4 * n_mod_bond; j += 4)
        {
            old_a1 = ed->mod_bond->item[j];
            old_a2 = ed->mod_bond->item[j + 1];
            new_a1 = ed->mod_bond->item[j + 2];
            new_a2 = ed->mod_bond->item[j + 3];
            if ((old_a1 == new_a1&&old_a2 == new_a2) || (old_a2 == new_a1&&old_a1 == new_a2))
            {
                continue;
            }
            ok = OrigAtData_RemoveBond(old_a1 - 1, old_a2 - 1, at, &bond_type, &bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            ok = OrigAtData_AddBond(new_a1 - 1, new_a2 - 1, at, bond_type, bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            /* Correct CRU blist lists */
            for (i = 0; i < p->n; i++)
            {
                OAD_PolymerUnit *u = p->units[i];
                if (!u->blist)
                {
                    /* No crossing bonds in the unit */
                    continue;
                }
                if ( bIsSameBond(u->blist[0], u->blist[1], old_a1, old_a2) )
                {
                    u->blist[0] = new_a1;
                    u->blist[1] = new_a2;
                }
                else if ( bIsSameBond(u->blist[2], u->blist[3], old_a1, old_a2) )
                {
                    u->blist[2] = new_a1;
                    u->blist[3] = new_a2;
                }
            }

            n_edits++;
        }
    }

    /* Modify coordinates */
    if (n_mod_coord)
    {
        for (j = 0; j < 2 * n_mod_coord; j += 2)
        {
            old_a1 = ed->mod_coord->item[j];
            new_a1 = ed->mod_bond->item[j + 1];
            swap_atoms_xyz(orig_at_data, old_a1 - 1, new_a1 - 1);
            n_edits++;
        }
    }

    
    /* Delete atoms */
    if (n_del_atom)
    {
        int nat0, nat, nacc;
        inp_ATOM *new_at = NULL, *new_at0=NULL;

        at_renum = (int *)inchi_calloc(n_max_stored, sizeof(int));
        if (!at_renum)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }
        /* all-purpose buffer */
        ibuf = (int *)inchi_calloc(n_max_stored, sizeof(int));
        if (!ibuf)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }

        fail = mark_atoms_to_delete_or_renumber(orig_at_data, ed, at_renum); 
        if (fail)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }
        /* Now remove atom by atom */

        nat0 = orig_at_data->num_inp_atoms;
        nat = nat0 - ed->del_atom->used;
        new_at = (inp_ATOM  *)inchi_calloc(nat, sizeof(new_at[0]));
        if (!new_at)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }

        for (i = 0, nacc = 0; i < nat0; i++)
        {
            AT_NUMB nbr0[MAXVAL];
            U_CHAR btype0[MAXVAL];
            int m, macc, valen;
            int new_num = at_renum[i];				
            if (new_num == -1)
            {
                /* Skip removed atom */
                continue;
            }

            /* Atom to keep; copy it */ 
            new_at0 = new_at + nacc;
            ++nacc;
            memcpy(new_at0, orig_at_data->at + i, sizeof(new_at[0]));
            /* Correct its own number(s) */
            new_at0->orig_at_number = new_num + 1;
            
            /* Correct its nbr number(s) */
            valen = new_at0->valence;
            memcpy(nbr0, new_at0->neighbor, valen * sizeof(AT_NUMB));
            memcpy(btype0, new_at0->bond_type, valen);
            memset(new_at0->neighbor, 0, valen); /* djb-rwth: memset_s C11/Annex K variant? */
            for (m = 0, macc=0; m < valen; m++)
            {
                int num2 = nbr0[m];
                int renum2 = at_renum[num2];
                if (renum2 == num2)
                {
                    /* keep old */
                    new_at0->neighbor[macc++] = num2; 
                }
                else if (renum2 == -1)
                {
                    /* skip and decrement valences */
                    new_at0->chem_bonds_valence -= btype0[m];
                    new_at0->valence--;
                }
                else 
                {
                    /* set renumbered */
                    new_at0->neighbor[macc++] = renum2;
                }
            }
        }

        if (new_at)
        {
            inchi_free(orig_at_data->at);
            orig_at_data->at = new_at;
        }

        orig_at_data->num_inp_atoms = nacc;
        orig_at_data->num_inp_bonds = 0;
        for (i = 0; i < nacc; i++)
        {
            orig_at_data->num_inp_bonds += new_at[i].valence;
        }
        orig_at_data->num_inp_bonds /= 2;

        /* Correct other data */

        /* Correct polymer data */
        if (p)
        {
            int iu;

            for (iu = 0; iu < p->n; iu++)
            {
                OAD_PolymerUnit* u = p->units[iu];
                int new_na, new_nb, new_bb;
                memset(ibuf, 0, n_max_stored * sizeof(ibuf[0])); /* djb-rwth: memset_s C11/Annex K variant? */
                if (u)
                {
                    if (u->alist)
                    {
                        memcpy(ibuf, u->alist, u->na * sizeof(ibuf[0]));
                        new_na = set_renumbered_or_delete(u->alist, ibuf, u->na, at_renum, 1);
                        if (new_na == -1)
                        {
                            *ret = _IS_ERROR;
                            goto exit_function;
                        }
                        u->na = new_na;
                    }
                    if (u->blist)
                    {
                        memcpy(ibuf, u->blist, 2 * (long long)u->nb * sizeof(int)); /* djb-rwth: cast operator added */
                        new_nb = set_renumbered_or_delete(u->blist, ibuf, 2*u->nb, at_renum, 1);
                        new_nb /= 2;
                        if (new_nb == -1)
                        {
                            *ret = _IS_ERROR;
                            goto exit_function;
                        }
                        u->nb = new_nb;
                    }


                    if (u->bkbonds)
                    {
                        int b, new_nbkbonds;
                        for (b = 0, new_nbkbonds = 0; b < u->nbkbonds; b++)
                        {
                            int bnd[2];
                            bnd[0] = u->bkbonds[b][0];
                            bnd[1] = u->bkbonds[b][1];
                            memcpy(ibuf, bnd, 2 * sizeof(ibuf[0]));
                            memset(u->bkbonds[b], 0, 2 * sizeof(u->bkbonds[0][0])); /* djb-rwth: memset_s C11/Annex K variant? */
                            new_bb = set_renumbered_or_delete(bnd, ibuf, 2, at_renum, 1);
                            if (new_bb == -1)
                            {
                                *ret = _IS_ERROR;
                                goto exit_function;
                            }
                            else if (new_bb == 2)
                            {
                                u->bkbonds[new_nbkbonds][0] = bnd[0];
                                u->bkbonds[new_nbkbonds][1] = bnd[1];
                                new_nbkbonds++;
                                /* OK, settled new nums or kept old */
                            }
                            else
                            {
                                continue;
                            }
                        }
                        u->nbkbonds = new_nbkbonds;
                    }

                    if (u->blist) /* djb-rwth: fixing a NULL pointer dereference */
                    {
                        u->cap1 = u->blist[0];
                        u->end_atom1 = u->blist[1];
                        u->cap2 = u->blist[2];
                        u->end_atom2 = u->blist[3];
                    }
                    if (u->cap1 < 1 || u->cap2 < 1 || u->end_atom1 < 1 || u->end_atom2 < 1)
                    {
                        *ret = _IS_ERROR;
                        goto exit_function;
                    }


                }
            }
        } /* if (p) */

        /* Correct V300 data */
        if (orig_at_data->v3000)
        {
            ;
        }


    } /* if (n_del_atom) */
    

exit_function:
    if (ibuf)
    {
        inchi_free(ibuf);
    }
    if (at_renum)
    {
        inchi_free(at_renum);
    }
    return n_edits;
}
    */
    // END INCHI C FUNCTION: OAD_StructureEdits_Apply
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_Apply
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64.
    // INCHI✔️❌: ITRACE_ expands to 0 && _inchi_trace; called debug helpers retain their active control flow.
    // INCHI✔️❌: Rust preserves the source's two-item new-bond stride, mod_bond coordinate lookup,
    // INCHI✔️❌: byte-count neighbor memset, partial mutation, and single exit cleanup order.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_Apply

    fn pointed_value<T: Clone + 'static>(
        heap: &SourceHeap,
        pointer: SourceMutPointer<T>,
    ) -> Result<T, SourceHeapError> {
        heap.slice(pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn allocation_failed(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }
    fn checked_double(value: i32) -> Result<i32, SourceHeapError> {
        value
            .checked_mul(2)
            .ok_or(SourceHeapError::SourceIntegerOverflow)
    }
    fn clear_neighbor_prefix(neighbors: &mut [u16; 20], byte_count: usize) {
        let mut bytes = [0_u8; 40];
        for (index, neighbor) in neighbors.iter().enumerate() {
            bytes[index * 2..index * 2 + 2].copy_from_slice(&neighbor.to_ne_bytes());
        }
        bytes[..byte_count].fill(0);
        for (index, neighbor) in neighbors.iter_mut().enumerate() {
            *neighbor = u16::from_ne_bytes([bytes[index * 2], bytes[index * 2 + 1]]);
        }
    }

    *return_status = _IS_OKAY as i32;
    let delete_atoms = pointed_value(heap, edits.del_atom)?;
    let delete_bonds = pointed_value(heap, edits.del_bond)?;
    let new_bonds = pointed_value(heap, edits.new_bond)?;
    let modified_bonds = pointed_value(heap, edits.mod_bond)?;
    let modified_coordinates = pointed_value(heap, edits.mod_coord)?;
    let number_of_deleted_atoms = delete_atoms.used;
    let number_of_deleted_bonds = delete_bonds.used / 2;
    let number_of_new_bonds = new_bonds.used / 3;
    let number_of_modified_bonds = modified_bonds.used / 4;
    let number_of_modified_coordinates = modified_coordinates.used / 2;
    let edit_count = number_of_deleted_atoms
        .checked_add(number_of_deleted_bonds)
        .and_then(|value| value.checked_add(number_of_new_bonds))
        .and_then(|value| value.checked_add(number_of_modified_bonds))
        .and_then(|value| value.checked_add(number_of_modified_coordinates))
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if edit_count < 1 {
        return Ok(0);
    }

    let atom_buffer_size = original_atom_data
        .num_inp_atoms
        .checked_add(1)
        .and_then(|value| value.checked_mul(2))
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let bond_buffer_size = original_atom_data
        .num_inp_bonds
        .checked_add(1)
        .and_then(|value| value.checked_mul(2))
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let maximum_stored = atom_buffer_size.max(bond_buffer_size);

    OrigAtData_DebugTrace(heap, original_atom_data)?;
    OAD_Polymer_DebugTrace(heap, original_atom_data.polymer)?;

    let atoms = original_atom_data.at;
    let polymer = original_atom_data.polymer;
    let mut atom_renumbering = SourceMutPointer::<i32>::null();
    let mut integer_buffer = SourceMutPointer::<i32>::null();
    let mut number_of_edits = 0_i32;
    let mut bond_type = tagINCHIBondType_INCHI_BOND_TYPE_NONE as i32;
    let mut bond_stereo = tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE as i32;

    let operation = (|| -> Result<i32, SourceHeapError> {
        if number_of_deleted_bonds != 0 {
            let end = checked_double(number_of_deleted_bonds)?;
            let values = heap.slice(delete_bonds.item.as_const())?.to_vec();
            let mut index = 0_i32;
            while index < end {
                let offset = usize::try_from(index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom1 = values
                    .get(offset)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let atom2 = values
                    .get(offset + 1)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let removed = OrigAtData_RemoveBond(
                    heap,
                    atom1,
                    atom2,
                    atoms,
                    &mut bond_type,
                    &mut bond_stereo,
                    &mut original_atom_data.num_inp_bonds,
                )?;
                if removed == 0 {
                    *return_status = _IS_ERROR as i32;
                    return Ok(number_of_edits);
                }
                number_of_edits = number_of_edits.wrapping_add(1);
                index = index.wrapping_add(2);
            }
        }

        if number_of_new_bonds != 0 {
            let end = checked_double(number_of_new_bonds)?;
            let values = heap.slice(new_bonds.item.as_const())?.to_vec();
            let mut index = 0_i32;
            while index < end {
                let offset = usize::try_from(index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom1 = values
                    .get(offset)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let atom2 = values
                    .get(offset + 1)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let added = OrigAtData_AddBond(
                    heap,
                    atom1,
                    atom2,
                    atoms,
                    bond_type,
                    bond_stereo,
                    &mut original_atom_data.num_inp_bonds,
                )?;
                if added == 0 {
                    *return_status = _IS_ERROR as i32;
                    return Ok(number_of_edits);
                }
                number_of_edits = number_of_edits.wrapping_add(1);
                index = index.wrapping_add(2);
            }
        }

        if number_of_modified_bonds != 0 {
            let end = number_of_modified_bonds
                .checked_mul(4)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let values = heap.slice(modified_bonds.item.as_const())?.to_vec();
            let mut index = 0_i32;
            while index < end {
                let offset = usize::try_from(index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let old_atom1 = *values
                    .get(offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let old_atom2 = *values
                    .get(offset + 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let new_atom1 = *values
                    .get(offset + 2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let new_atom2 = *values
                    .get(offset + 3)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if (old_atom1 == new_atom1 && old_atom2 == new_atom2)
                    || (old_atom2 == new_atom1 && old_atom1 == new_atom2)
                {
                    index = index.wrapping_add(4);
                    continue;
                }
                let removed = OrigAtData_RemoveBond(
                    heap,
                    old_atom1
                        .checked_sub(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    old_atom2
                        .checked_sub(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    atoms,
                    &mut bond_type,
                    &mut bond_stereo,
                    &mut original_atom_data.num_inp_bonds,
                )?;
                if removed == 0 {
                    *return_status = _IS_ERROR as i32;
                    return Ok(number_of_edits);
                }
                let added = OrigAtData_AddBond(
                    heap,
                    new_atom1
                        .checked_sub(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    new_atom2
                        .checked_sub(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    atoms,
                    bond_type,
                    bond_stereo,
                    &mut original_atom_data.num_inp_bonds,
                )?;
                if added == 0 {
                    *return_status = _IS_ERROR as i32;
                    return Ok(number_of_edits);
                }

                if polymer.is_null() {
                    return Err(SourceHeapError::UnsupportedSourceBehavior);
                }
                let polymer_value = pointed_value(heap, polymer)?;
                let unit_count = usize::try_from(polymer_value.n.max(0))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let unit_pointers = heap.slice(polymer_value.units.as_const())?.to_vec();
                for unit_index in 0..unit_count {
                    let unit_pointer = *unit_pointers
                        .get(unit_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let unit = pointed_value(heap, unit_pointer)?;
                    if unit.blist.is_null() {
                        continue;
                    }
                    let mut list = heap.slice(unit.blist.as_const())?.to_vec();
                    if list.len() < 4 {
                        return Err(SourceHeapError::PointerOutOfBounds);
                    }
                    if bIsSameBond(list[0], list[1], old_atom1, old_atom2) != 0 {
                        list[0] = new_atom1;
                        list[1] = new_atom2;
                        heap.slice_mut(unit.blist)?[..2].copy_from_slice(&list[..2]);
                    } else if bIsSameBond(list[2], list[3], old_atom1, old_atom2) != 0 {
                        list[2] = new_atom1;
                        list[3] = new_atom2;
                        heap.slice_mut(unit.blist)?[2..4].copy_from_slice(&list[2..4]);
                    }
                }
                number_of_edits = number_of_edits.wrapping_add(1);
                index = index.wrapping_add(4);
            }
        }

        if number_of_modified_coordinates != 0 {
            let end = checked_double(number_of_modified_coordinates)?;
            let coordinate_values = heap.slice(modified_coordinates.item.as_const())?.to_vec();
            let bond_values = heap.slice(modified_bonds.item.as_const())?.to_vec();
            let mut index = 0_i32;
            while index < end {
                let offset = usize::try_from(index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let old_atom = coordinate_values
                    .get(offset)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let new_atom = bond_values
                    .get(offset + 1)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                swap_atoms_xyz(
                    heap,
                    Some(original_atom_data),
                    old_atom
                        .checked_sub(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    new_atom
                        .checked_sub(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                )?;
                number_of_edits = number_of_edits.wrapping_add(1);
                index = index.wrapping_add(2);
            }
        }

        if number_of_deleted_atoms != 0 {
            atom_renumbering = match inchi_calloc::<i32>(heap, maximum_stored as u64, 4) {
                Ok(pointer) => pointer,
                Err(error) if allocation_failed(error) => {
                    *return_status = _IS_ERROR as i32;
                    return Ok(number_of_edits);
                }
                Err(error) => return Err(error),
            };
            integer_buffer = match inchi_calloc::<i32>(heap, maximum_stored as u64, 4) {
                Ok(pointer) => pointer,
                Err(error) if allocation_failed(error) => {
                    *return_status = _IS_ERROR as i32;
                    return Ok(number_of_edits);
                }
                Err(error) => return Err(error),
            };
            let failure = mark_atoms_to_delete_or_renumber(
                heap,
                original_atom_data,
                edits,
                atom_renumbering,
            )?;
            if failure != 0 {
                *return_status = _IS_ERROR as i32;
                return Ok(number_of_edits);
            }

            let original_atom_count = original_atom_data.num_inp_atoms;
            let current_deletions = pointed_value(heap, edits.del_atom)?;
            let new_atom_count = original_atom_count
                .checked_sub(current_deletions.used)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let new_atoms = match inchi_calloc::<inp_ATOM>(heap, new_atom_count as u64, 176) {
                Ok(pointer) => pointer,
                Err(error) if allocation_failed(error) => {
                    *return_status = _IS_ERROR as i32;
                    return Ok(number_of_edits);
                }
                Err(error) => return Err(error),
            };

            let source_atoms = heap.slice(original_atom_data.at.as_const())?.to_vec();
            let source_count = usize::try_from(original_atom_count)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mut accepted = 0_usize;
            for atom_index in 0..source_count {
                let new_number = *heap
                    .slice(atom_renumbering.as_const())?
                    .get(atom_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if new_number == -1 {
                    continue;
                }
                let mut atom = source_atoms
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                atom.orig_at_number = new_number.wrapping_add(1) as u16;
                let valence = usize::try_from(atom.valence)
                    .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                if valence > atom.neighbor.len() {
                    return Err(SourceHeapError::UnsupportedSourceBehavior);
                }
                let neighbors = atom.neighbor[..valence].to_vec();
                let bond_types = atom.bond_type[..valence].to_vec();
                clear_neighbor_prefix(&mut atom.neighbor, valence);
                let mut retained = 0_usize;
                for bond_index in 0..valence {
                    let old_neighbor = usize::from(neighbors[bond_index]);
                    let renumbered_neighbor = *heap
                        .slice(atom_renumbering.as_const())?
                        .get(old_neighbor)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if renumbered_neighbor == old_neighbor as i32 {
                        atom.neighbor[retained] = old_neighbor as u16;
                        retained += 1;
                    } else if renumbered_neighbor == -1 {
                        atom.chem_bonds_valence = atom
                            .chem_bonds_valence
                            .wrapping_sub(bond_types[bond_index] as i8);
                        atom.valence = atom.valence.wrapping_sub(1);
                    } else {
                        atom.neighbor[retained] = renumbered_neighbor as u16;
                        retained += 1;
                    }
                }
                *heap
                    .slice_mut(new_atoms)?
                    .get_mut(accepted)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = atom;
                accepted += 1;
            }

            inchi_free(heap, original_atom_data.at)?;
            original_atom_data.at = new_atoms;
            original_atom_data.num_inp_atoms =
                i32::try_from(accepted).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            original_atom_data.num_inp_bonds = 0;
            for atom in heap
                .slice(new_atoms.as_const())?
                .get(..accepted)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                original_atom_data.num_inp_bonds = original_atom_data
                    .num_inp_bonds
                    .wrapping_add(i32::from(atom.valence));
            }
            original_atom_data.num_inp_bonds /= 2;

            if !polymer.is_null() {
                let polymer_value = pointed_value(heap, polymer)?;
                let unit_count = usize::try_from(polymer_value.n.max(0))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let unit_pointers = heap.slice(polymer_value.units.as_const())?.to_vec();
                for unit_index in 0..unit_count {
                    heap.slice_mut(integer_buffer)?.fill(0);
                    let unit_pointer = *unit_pointers
                        .get(unit_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if unit_pointer.is_null() {
                        continue;
                    }
                    let mut unit = pointed_value(heap, unit_pointer)?;
                    if !unit.alist.is_null() {
                        let count = usize::try_from(unit.na)
                            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                        let values = heap
                            .slice(unit.alist.as_const())?
                            .get(..count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .to_vec();
                        heap.slice_mut(integer_buffer)?[..count].copy_from_slice(&values);
                        unit.na = set_renumbered_or_delete(
                            heap,
                            unit.alist,
                            integer_buffer,
                            unit.na,
                            atom_renumbering.as_const(),
                            1,
                        )?;
                    }
                    if !unit.blist.is_null() {
                        let count = unit
                            .nb
                            .checked_mul(2)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                        let count = usize::try_from(count)
                            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                        let values = heap
                            .slice(unit.blist.as_const())?
                            .get(..count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .to_vec();
                        heap.slice_mut(integer_buffer)?[..count].copy_from_slice(&values);
                        let mut new_count = set_renumbered_or_delete(
                            heap,
                            unit.blist,
                            integer_buffer,
                            i32::try_from(count)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                            atom_renumbering.as_const(),
                            1,
                        )?;
                        new_count /= 2;
                        if new_count == -1 {
                            *return_status = _IS_ERROR as i32;
                            return Ok(number_of_edits);
                        }
                        unit.nb = new_count;
                    }
                    if !unit.bkbonds.is_null() {
                        let backbone_count = usize::try_from(unit.nbkbonds)
                            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                        let rows = heap.slice(unit.bkbonds.as_const())?.to_vec();
                        let renumbering = heap.slice(atom_renumbering.as_const())?.to_vec();
                        let mut new_backbone_count = 0_usize;
                        for backbone_index in 0..backbone_count {
                            let row = *rows
                                .get(backbone_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut bond = [
                                *heap
                                    .slice(row.as_const())?
                                    .first()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                *heap
                                    .slice(row.as_const())?
                                    .get(1)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            ];
                            heap.slice_mut(integer_buffer)?[..2].copy_from_slice(&bond);
                            heap.slice_mut(row)?
                                .get_mut(..2)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .fill(0);
                            let copied = bond;
                            bond.fill(0);
                            heap.slice_mut(integer_buffer)?[..2].copy_from_slice(&copied);
                            let new_count = set_renumbered_or_delete_from_copy(
                                &mut bond,
                                &copied,
                                &renumbering,
                                1,
                            )?;
                            if new_count == -1 {
                                *return_status = _IS_ERROR as i32;
                                return Ok(number_of_edits);
                            } else if new_count == 2 {
                                let target = *rows
                                    .get(new_backbone_count)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                heap.slice_mut(target)?
                                    .get_mut(..2)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                    .copy_from_slice(&bond);
                                new_backbone_count += 1;
                            }
                        }
                        unit.nbkbonds = i32::try_from(new_backbone_count)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    }
                    if !unit.blist.is_null() {
                        let list = heap.slice(unit.blist.as_const())?;
                        if list.len() < 4 {
                            return Err(SourceHeapError::PointerOutOfBounds);
                        }
                        unit.cap1 = list[0];
                        unit.end_atom1 = list[1];
                        unit.cap2 = list[2];
                        unit.end_atom2 = list[3];
                    }
                    heap.slice_mut(unit_pointer)?[0] = unit.clone();
                    if unit.cap1 < 1
                        || unit.cap2 < 1
                        || unit.end_atom1 < 1
                        || unit.end_atom2 < 1
                    {
                        *return_status = _IS_ERROR as i32;
                        return Ok(number_of_edits);
                    }
                }
            }

            if !original_atom_data.v3000.is_null() {
                // The selected source branch intentionally has an empty body.
            }
        }

        Ok(number_of_edits)
    })();

    let buffer_cleanup = inchi_free(heap, integer_buffer);
    let renumbering_cleanup = inchi_free(heap, atom_renumbering);
    match (operation, buffer_cleanup, renumbering_cleanup) {
        (Err(error), _, _) => Err(error),
        (Ok(_), Err(error), _) | (Ok(_), Ok(()), Err(error)) => Err(error),
        (Ok(result), Ok(()), Ok(())) => Ok(result),
    }
}

#[rustfmt::skip]
pub(crate) fn set_renumbered_or_delete(
    heap: &mut SourceHeap,
    number: SourceMutPointer<i32>,
    buffer: SourceMutPointer<i32>,
    number_of_elements: i32,
    renumbering: SourceConstPointer<i32>,
    base: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2775 set_renumbered_or_delete
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int set_renumbered_or_delete( int *number,	/* numbers to correct */
                              int *buf,		/* must be enough size to hold nelems elements */
                              int nelems,	/* initial size of numbers */
                              int *renum,	/* new numbers in order of old numbers; always 0-based */
                              int base)
{
    int i, new_nelems;
    memcpy(buf, number, nelems * sizeof(int));
    memset(number, 0, nelems * sizeof(int)); /* djb-rwth: memset_s C11/Annex K variant? */
    for (i = 0, new_nelems = 0; i < nelems; i++)
    {
        int new_num = renum[ buf[i]-base ] + base;
        if (new_num == (base-1))
        {
            continue;
        }
        else
        {
            number[new_nelems++] = new_num;
        }
    }
    return new_nelems;
}
    */
    // END INCHI C FUNCTION: set_renumbered_or_delete

    let count = usize::try_from(number_of_elements)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    if count == 0 {
        return Ok(0);
    }

    let copied = heap
        .slice(number.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(buffer)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(&copied);
    heap.slice_mut(number)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(0);

    let copied = heap.slice(buffer.as_const())?[..count].to_vec();
    let renumbering = heap.slice(renumbering)?.to_vec();
    set_renumbered_or_delete_from_copy(
        &mut heap.slice_mut(number)?[..count],
        &copied,
        &renumbering,
        base,
    )
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ProcessOneStructureExCore(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &mut INPUT_PARMS,
    mut title: Option<&mut [i8]>,
    inchi_components: &mut [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    aux_components: &mut [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    mut input_file: Option<&mut INCHI_IOSTREAM>,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    mut output_file: Option<&mut INCHI_IOSTREAM>,
    mut problem_file: Option<&mut INCHI_IOSTREAM>,
    original_input_pointer: SourceMutPointer<ORIG_ATOM_DATA>,
    prepared_input_pointer: SourceMutPointer<ORIG_ATOM_DATA>,
    input_number: i64,
    mut string_buffer: Option<&mut INCHI_IOS_STRING>,
    save_option_bits: u8,
    stdout: SourceMutPointer<crate::source_types::FILE>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2803 ProcessOneStructureExCore
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int ProcessOneStructureExCore( struct tagINCHI_CLOCK *ic,
                               struct tagCANON_GLOBALS  *CG,
                               STRUCT_DATA *sd,
                               INPUT_PARMS *ip,
                               char *szTitle,
                               PINChI2 *pINChI2[INCHI_NUM],
                               PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                               INCHI_IOSTREAM *inp_file,
                               INCHI_IOSTREAM *log_file,
                               INCHI_IOSTREAM *out_file,
                               INCHI_IOSTREAM *prb_file,
                               ORIG_ATOM_DATA *orig_inp_data,
                               ORIG_ATOM_DATA *prep_inp_data,
                               long num_inp,
                               INCHI_IOS_STRING *strbuf,
                               unsigned char save_opt_bits )
{
    int res = _IS_OKAY;
    int mind_polymers;

#ifdef TARGET_LIB_FOR_WINCHI
    inchi_ios_free_str( out_file );
    inchi_ios_print(out_file, "Structure: %d\n", num_inp);
#endif

    /* Polymer and pseudoelement specific */
    res = ValidateAndPreparePolymerAndPseudoatoms( ic, CG, sd, ip, szTitle, pINChI2, pINChI_Aux2,
                                                    inp_file, log_file, out_file, prb_file,
                                                    orig_inp_data, prep_inp_data, num_inp, strbuf,
                                                    save_opt_bits, &mind_polymers);
    if (res== _IS_ERROR || res == _IS_FATAL )
    {
        return res;
    }

    /* Call the very actual worker placed under this (ProcessOneStructureCore) wrapper */
    res = ProcessOneStructure( ic, CG, sd, ip, szTitle, pINChI2, pINChI_Aux2,
                               inp_file, log_file, out_file, prb_file,
                               orig_inp_data, prep_inp_data, num_inp, strbuf,
                               save_opt_bits);

    if ( (res == _IS_OKAY  || res == _IS_WARNING ) && mind_polymers )
    {
        /* Post-edit the polymer layer at older polymer treatment modes (1.05, 1.05+) */
        if (ip->bPolymers == POLYMERS_LEGACY || ip->bPolymers == POLYMERS_LEGACY_PLUS)
        {
            /* Cut and hide "Zz" and related things in InChI (AuxInfo has a specifics). */
            int n_pzz = 0, n_zy = orig_inp_data->n_zy;
            if (orig_inp_data->polymer)
            {
                n_pzz = orig_inp_data->polymer->n_pzz;
            }

            EditINCHI_HidePolymerZz(out_file, n_pzz, n_zy);

        }
    }

#ifdef TARGET_LIB_FOR_WINCHI
/*	if ( res == _IS_ERROR || res == _IS_FATAL )
    {
        inchi_ios_print(out_file, "Error %d (%s)\n", sd->nErrorCode, sd->pStrErrStruct);
    }
*/
    /*push_to_winchi_text_window(out_file); */
    /*inchi_ios_free_str(out_file);*/
    /*inchi_ios_flush(out_file);*/
#endif

    return res;
}
    */
    // END INCHI C FUNCTION: ProcessOneStructureExCore
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ProcessOneStructureExCore
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; TARGET_LIB_FOR_WINCHI is inactive.
    // INCHI✔️❌: The modeled worker receives explicit stdout and clock values for its active libc globals.
    // INCHI✔️❌: SourceHeap temporarily removes and reinserts the original-data allocation to preserve aliasing
    // INCHI✔️❌: while obtaining the mutable source object, which is materially slower than a C pointer dereference.
    // END INCHI ACTIVE MACRO CONFIGURATION: ProcessOneStructureExCore

    let mut mind_polymers = 0_i32;
    let validation_result = if original_input_pointer.is_null() {
        ValidateAndPreparePolymerAndPseudoatoms(
            heap,
            clock,
            canonical_globals,
            structure_data,
            input_parameters,
            title.as_deref_mut(),
            inchi_components,
            aux_components,
            input_file.as_deref_mut(),
            log_file.as_deref_mut(),
            output_file.as_deref_mut(),
            problem_file.as_deref_mut(),
            None,
            None,
            input_number,
            string_buffer.as_deref_mut(),
            save_option_bits,
            &mut mind_polymers,
        )?
    } else {
        heap.with_slice_mut_and_heap_mut(original_input_pointer, |original, heap| {
            let original = original
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            ValidateAndPreparePolymerAndPseudoatoms(
                heap,
                clock,
                canonical_globals,
                structure_data,
                input_parameters,
                title.as_deref_mut(),
                inchi_components,
                aux_components,
                input_file.as_deref_mut(),
                log_file.as_deref_mut(),
                output_file.as_deref_mut(),
                problem_file.as_deref_mut(),
                Some(original),
                None,
                input_number,
                string_buffer.as_deref_mut(),
                save_option_bits,
                &mut mind_polymers,
            )
        })?
    };
    if validation_result == _IS_ERROR as i32 || validation_result == _IS_FATAL as i32 {
        return Ok(validation_result);
    }

    let result = ProcessOneStructure(
        heap,
        clock,
        canonical_globals,
        structure_data,
        input_parameters,
        title.as_deref_mut(),
        inchi_components,
        aux_components,
        input_file.as_deref_mut(),
        log_file.as_deref_mut(),
        output_file.as_deref_mut(),
        problem_file.as_deref_mut(),
        original_input_pointer,
        prepared_input_pointer,
        input_number,
        string_buffer.as_deref_mut(),
        save_option_bits,
        stdout,
        clock_result,
    )?;

    if (result == _IS_OKAY as i32 || result == crate::source_types::_IS_WARNING as i32)
        && mind_polymers != 0
        && (input_parameters.bPolymers == POLYMERS_LEGACY as i32
            || input_parameters.bPolymers == POLYMERS_LEGACY_PLUS as i32)
    {
        let original = heap
            .slice(original_input_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let number_of_placeholder_zy = original.n_zy;
        let number_of_polymer_zz = if original.polymer.is_null() {
            0
        } else {
            heap.slice(original.polymer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .n_pzz
        };
        EditINCHI_HidePolymerZz(
            heap,
            output_file
                .as_deref_mut()
                .ok_or(SourceHeapError::NullPointer)?,
            number_of_polymer_zz,
            number_of_placeholder_zy,
        )?;
    }

    Ok(result)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ValidateAndPreparePolymerAndPseudoatoms(
    heap: &mut SourceHeap,
    _clock: SourceMutPointer<INCHI_CLOCK>,
    _canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    _title: Option<&mut [i8]>,
    _inchi_components: &mut [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    _aux_components: &mut [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    _input_file: Option<&mut INCHI_IOSTREAM>,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    _output_file: Option<&mut INCHI_IOSTREAM>,
    _problem_file: Option<&mut INCHI_IOSTREAM>,
    original_atom_data: Option<&mut ORIG_ATOM_DATA>,
    _prepared_atom_data: Option<&mut ORIG_ATOM_DATA>,
    input_number: i64,
    _string_buffer: Option<&mut INCHI_IOS_STRING>,
    _save_option_bits: u8,
    mind_polymers: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2879 ValidateAndPreparePolymerAndPseudoatoms
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int ValidateAndPreparePolymerAndPseudoatoms( struct tagINCHI_CLOCK *ic,
                                             struct tagCANON_GLOBALS *CG,
                                             STRUCT_DATA *sd,
                                             INPUT_PARMS *ip,
                                             char *szTitle,
                                             PINChI2 *pINChI2[INCHI_NUM],
                                             PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                                             INCHI_IOSTREAM *inp_file,
                                             INCHI_IOSTREAM *log_file,
                                             INCHI_IOSTREAM *out_file,
                                             INCHI_IOSTREAM *prb_file,
                                             ORIG_ATOM_DATA *orig_inp_data,
                                             ORIG_ATOM_DATA *prep_inp_data,
                                             long num_inp,
                                             INCHI_IOS_STRING *strbuf,
                                             unsigned char save_opt_bits,
                                             int *mind_polymers )
{
    int res = _IS_OKAY;

    int mind_pseudoelements = 0;
    
    /* djb-rwth: fixing coverity ID #499512 */
    if (!orig_inp_data)
    {
        goto exit_function;
    }

    *mind_polymers = orig_inp_data->polymer && orig_inp_data->polymer->n > 0;
    *mind_polymers = *mind_polymers && orig_inp_data->valid_polymer &&
        (ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE);
    mind_pseudoelements = (ip->bNPZz == 1) || (ip->bPolymers != POLYMERS_NO);


    /* Validate the data */
    res = OAD_ValidatePolymerAndPseudoElementData(orig_inp_data,
        ip->bPolymers,
        ip->bNPZz,
        sd->pStrErrStruct,
        ip->bNoWarnings);
    if (res)
    {
        sd->nErrorCode = res;
        inchi_ios_eprint( log_file, "Error %d (%s) structure #%ld.%s%s%s%s\n",
                          sd->nErrorCode, sd->pStrErrStruct, num_inp,
                          SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        res = _IS_ERROR;
        orig_inp_data->num_inp_atoms = -1; /* djb-rwth: fixing coverity ID #499522 */
        goto exit_function;
    }

    if (*mind_polymers || mind_pseudoelements)
    {
        /*OrigAtData_DebugTrace(orig_inp_data);*/
        if (*mind_polymers &&
            ip->bPolymers != POLYMERS_MODERN &&
            (ip->bFrameShiftScheme == FSS_STARS_CYCLED || ip->bFrameShiftScheme == FSS_STARS_CYCLED_SORTED))
        {
            /*  Analyze and cyclize frame-shift eligible CRUs using InChI canonical numbers 
                (do this only at older polymer treatment modes 1.05, 1.05+)					
            */
            res = OAD_Polymer_CyclizeCloseableUnits( orig_inp_data,
                                                     ip->bPolymers,
                                                     sd->pStrErrStruct,
                                                     ip->bNoWarnings );
            if (res)
            {
                sd->nErrorCode = res;
                AddErrorMessage(sd->pStrErrStruct, "Error while processing polymer-related input");
                res = _IS_ERROR;
                orig_inp_data->num_inp_atoms = -1;
                goto exit_function;
            }
            /*OrigAtData_DebugTrace(orig_inp_data);*/
        }
    }

exit_function:

    return res;
}
    */
    // END INCHI C FUNCTION: ValidateAndPreparePolymerAndPseudoatoms
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ValidateAndPreparePolymerAndPseudoatoms
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; long is modeled as i64.
    // INCHI✔️❌: SDF_LBL_VAL and WarningMessage/AddErrorMessage use the already-ported active header behavior.
    // INCHI✔️❌: Parameters not referenced by the active C body remain unevaluated in Rust.
    // INCHI✔️❌: Performance is materially worse because formatted log emission builds temporary owned byte buffers.
    // END INCHI ACTIVE MACRO CONFIGURATION: ValidateAndPreparePolymerAndPseudoatoms

    let Some(original_atom_data) = original_atom_data else {
        return Ok(_IS_OKAY as i32);
    };

    let has_polymer_units = if original_atom_data.polymer.is_null() {
        false
    } else {
        heap.slice(original_atom_data.polymer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .n
            > 0
    };
    *mind_polymers = i32::from(
        has_polymer_units
            && original_atom_data.valid_polymer != 0
            && (input_parameters.nInputType == tagInputType_INPUT_MOLFILE
                || input_parameters.nInputType == tagInputType_INPUT_SDFILE),
    );
    let mind_pseudoelements = input_parameters.bNPZz == 1
        || input_parameters.bPolymers != crate::source_types::POLYMERS_NO as i32;

    let mut result = OAD_ValidatePolymerAndPseudoElementData(
        heap,
        original_atom_data,
        input_parameters.bPolymers,
        input_parameters.bNPZz,
        Some(&mut structure_data.pStrErrStruct),
        input_parameters.bNoWarnings,
    )?;
    if result != 0 {
        structure_data.nErrorCode = result;
        let error_end = structure_data
            .pStrErrStruct
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let mut output = format!("Error {} (", structure_data.nErrorCode).into_bytes();
        output.extend(
            structure_data.pStrErrStruct[..error_end]
                .iter()
                .map(|byte| *byte as u8),
        );
        output.extend_from_slice(b") structure #");
        output.extend_from_slice(input_number.to_string().as_bytes());
        output.push(b'.');
        output.extend_from_slice(&sdf_label_value(heap, input_parameters)?);
        output.push(b'\n');
        let _ = eprint_bytes(heap, log_file.as_deref_mut(), &output)?;
        result = _IS_ERROR as i32;
        original_atom_data.num_inp_atoms = -1;
        return Ok(result);
    }

    if (*mind_polymers != 0 || mind_pseudoelements)
        && *mind_polymers != 0
        && input_parameters.bPolymers != POLYMERS_MODERN as i32
        && (input_parameters.bFrameShiftScheme
            == tagFrameShifScheme_FSS_STARS_CYCLED as i32
            || input_parameters.bFrameShiftScheme
                == tagFrameShifScheme_FSS_STARS_CYCLED_SORTED as i32)
    {
        result = OAD_Polymer_CyclizeCloseableUnits(
            heap,
            original_atom_data,
            input_parameters.bPolymers,
            Some(&mut structure_data.pStrErrStruct),
            input_parameters.bNoWarnings,
        )?;
        if result != 0 {
            structure_data.nErrorCode = result;
            let message = b"Error while processing polymer-related input\0"
                .map(|byte| byte as i8);
            let _ = AddErrorMessage(
                Some(&mut structure_data.pStrErrStruct),
                Some(&message),
            )?;
            result = _IS_ERROR as i32;
            original_atom_data.num_inp_atoms = -1;
            return Ok(result);
        }
    }

    Ok(result)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_ProcessOneStructureNoEdits(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    structure_data: &STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    title: &[i8],
    inchi_components: Option<&[SourceMutPointer<PINChI2>; INCHI_NUM as usize]>,
    aux_components: Option<&[SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize]>,
    input_file: SourceMutPointer<INCHI_IOSTREAM>,
    log_file: SourceMutPointer<INCHI_IOSTREAM>,
    output_file: SourceMutPointer<INCHI_IOSTREAM>,
    problem_file: SourceMutPointer<INCHI_IOSTREAM>,
    original_input: &ORIG_ATOM_DATA,
    prepared_input: &[ORIG_ATOM_DATA],
    input_number: i64,
    string_buffer: Option<&INCHI_IOS_STRING>,
    save_option_bits: u8,
    number_of_polymer_zz: &mut i32,
    inchi_output: &mut SourceMutPointer<i8>,
    aux_output: &mut SourceMutPointer<i8>,
    stdout: SourceMutPointer<crate::source_types::FILE>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2967 OAD_ProcessOneStructureNoEdits
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_ProcessOneStructureNoEdits( struct tagINCHI_CLOCK    *ic,
                                    struct tagCANON_GLOBALS  *CG,
                                    STRUCT_DATA              *sd,
                                    INPUT_PARMS              *ip,
                                    char                     *szTitle,
                                    PINChI2                  *pINChI2[INCHI_NUM],
                                    PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                                    INCHI_IOSTREAM           *inp_file,
                                    INCHI_IOSTREAM           *log_file,
                                    INCHI_IOSTREAM           *out_file,
                                    INCHI_IOSTREAM           *prb_file,
                                    ORIG_ATOM_DATA           *orig_inp_data,
                                    ORIG_ATOM_DATA           *prep_inp_data,
                                    long                     num_inp,
                                    INCHI_IOS_STRING         *strbuf,
                                    unsigned char            save_opt_bits,
                                    int                      *n_pzz,
                                    char					 **sinchi, 
                                    char					 **saux)
{
    size_t slen;
    int ret = _IS_OKAY, dup_fail = 0;
    POSEContext dup_context, *dup = &dup_context;

    *n_pzz = 0;
    *sinchi = NULL;
    *saux = NULL;

    /* Make a working copy of all the native input */
    dup_fail = POSEContext_Init(dup, sd, ip, szTitle, pINChI2, pINChI_Aux2,
        inp_file, log_file, out_file, prb_file,
        orig_inp_data, prep_inp_data,
        num_inp, strbuf, save_opt_bits);
    if (dup_fail)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    /* Set necessary for this specific case options */
    dup->orig_inp_data->polymer->treat = dup->ip.bPolymers = POLYMERS_MODERN;
    dup->ip.bFoldPolymerSRU = 0;
    dup->ip.bFrameShiftScheme = FSS_NONE;
    dup->ip.bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO);
    dup->ip.bDisplay = dup->ip.bDisplayCompositeResults = dup->ip.bDisplayEachComponentINChI = 0;
    dup->ip.bFoldPolymerSRU = 0;
    /*dup->ip.bTautFlags |= TG_FLAG_RECONNECT_COORD;*/

    /* Calculate */
    ret = ProcessOneStructureExCore( ic, CG, &dup->sd, &dup->ip, dup->szTitle,
                                     dup->pINChI2, dup->pINChI_Aux2,
                                     dup->inp_file, dup->log_file,
                                     dup->out_file, dup->prb_file,
                                     dup->orig_inp_data, dup->prep_inp_data,
                                     dup->num_inp, dup->strbuf, dup->save_opt_bits );

    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }
    
    *n_pzz = dup->orig_inp_data->polymer->n_pzz;
    /* Extract InChI */
    slen = dup->out_file->s.nUsedLength;
    extract_inchi_substring(sinchi, dup->out_file->s.pStr, slen);
    if (!*sinchi)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }
    /* Extract AuxInfo */
    slen = dup->out_file->s.nUsedLength;
    extract_auxinfo_substring(saux, dup->out_file->s.pStr, slen);
    if (!*saux)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }

exit_function:
    POSEContext_Free(dup);

    return ret;
}
    */
    // END INCHI C FUNCTION: OAD_ProcessOneStructureNoEdits
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_ProcessOneStructureNoEdits
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; POLYMERS_MODERN=1 and FSS_NONE=1.
    // INCHI✔️❌: Explicit stdout and clock_result model active globals passed through the worker.
    // INCHI✔️❌: Cloned stream and strbuf records add copying and allocation-map lookups absent in C.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_ProcessOneStructureNoEdits

    *number_of_polymer_zz = 0;
    *inchi_output = SourceMutPointer::null();
    *aux_output = SourceMutPointer::null();

    let mut duplicate = crate::source_types::POSEContext::default();
    let execution = (|| -> Result<i32, SourceHeapError> {
        let duplicate_failure = POSEContext_Init(
            heap, &mut duplicate, Some(structure_data), Some(input_parameters), title,
            inchi_components, aux_components, input_file, log_file, output_file, problem_file,
            Some(original_input), Some(prepared_input), input_number, string_buffer,
            save_option_bits,
        )?;
        if duplicate_failure != 0 {
            return Ok(_IS_ERROR as i32);
        }

        let polymer = heap
            .slice(duplicate.orig_inp_data.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .polymer;
        heap.slice_mut(polymer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .treat = POLYMERS_MODERN as i32;
        duplicate.ip.bPolymers = POLYMERS_MODERN as i32;
        duplicate.ip.bFoldPolymerSRU = 0;
        duplicate.ip.bFrameShiftScheme = tagFrameShifScheme_FSS_NONE as i32;
        duplicate.ip.bINChIOutputOptions &=
            !((INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO) as i32);
        duplicate.ip.bDisplay = 0;
        duplicate.ip.bDisplayCompositeResults = 0;
        duplicate.ip.bDisplayEachComponentINChI = 0;
        duplicate.ip.bFoldPolymerSRU = 0;

        let mut input_stream = if duplicate.inp_file.is_null() {
            None
        } else {
            Some(heap.slice(duplicate.inp_file.as_const())?
                .first().ok_or(SourceHeapError::PointerOutOfBounds)?.clone())
        };
        let stored_streams = heap.slice(duplicate.out_file.as_const())?;
        let mut streams: [INCHI_IOSTREAM; 3] = std::array::from_fn(|index| {
            stored_streams[index].clone()
        });
        let mut string = heap.slice(duplicate.strbuf.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?.clone();
        let (output_streams, remainder) = streams.split_at_mut(1);
        let (log_streams, problem_streams) = remainder.split_at_mut(1);
        let worker_result = ProcessOneStructureExCore(
            heap, clock, canonical_globals, &mut duplicate.sd, &mut duplicate.ip,
            Some(duplicate.szTitle.as_mut_slice()), &mut duplicate.pINChI2,
            &mut duplicate.pINChI_Aux2, input_stream.as_mut(), log_streams.first_mut(),
            output_streams.first_mut(), problem_streams.first_mut(), duplicate.orig_inp_data,
            duplicate.prep_inp_data, duplicate.num_inp, Some(&mut string),
            duplicate.save_opt_bits, stdout, clock_result,
        );

        if let Some(input_stream) = input_stream {
            *heap.slice_mut(duplicate.inp_file)?.first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = input_stream;
        }
        heap.slice_mut(duplicate.out_file)?[..3].clone_from_slice(&streams);
        duplicate.inchi_file.clone_from_slice(&streams);
        *heap.slice_mut(duplicate.strbuf)?.first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = string.clone();
        duplicate.temp_string_container = string;

        let ret = worker_result?;
        if ret == _IS_FATAL as i32 || ret == _IS_ERROR as i32 {
            return Ok(ret);
        }

        let original = heap.slice(duplicate.orig_inp_data.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        *number_of_polymer_zz = heap.slice(original.polymer.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?.n_pzz;

        let output_stream = heap.slice(duplicate.out_file.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let string_length = u64::try_from(output_stream.s.nUsedLength)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        extract_inchi_substring(heap, inchi_output, output_stream.s.pStr.as_const(), string_length)?;
        if inchi_output.is_null() {
            return Ok(_IS_ERROR as i32);
        }

        let output_stream = heap.slice(duplicate.out_file.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let string_length = u64::try_from(output_stream.s.nUsedLength)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        extract_auxinfo_substring(heap, aux_output, output_stream.s.pStr.as_const(), string_length)?;
        if aux_output.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        Ok(ret)
    })();

    let cleanup = POSEContext_Free(heap, &mut duplicate);
    match (execution, cleanup) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(ret), Ok(())) => Ok(ret),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_ProcessOneStructure105Plus(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    structure_data: &STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    title: &[i8],
    inchi_components: Option<&[SourceMutPointer<PINChI2>; INCHI_NUM as usize]>,
    aux_components: Option<&[SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize]>,
    input_file: SourceMutPointer<INCHI_IOSTREAM>,
    log_file: SourceMutPointer<INCHI_IOSTREAM>,
    output_file: SourceMutPointer<INCHI_IOSTREAM>,
    problem_file: SourceMutPointer<INCHI_IOSTREAM>,
    original_input: &ORIG_ATOM_DATA,
    prepared_input: &[ORIG_ATOM_DATA],
    input_number: i64,
    string_buffer: Option<&INCHI_IOS_STRING>,
    save_option_bits: u8,
    inchi_output: &mut SourceMutPointer<i8>,
    aux_output: &mut SourceMutPointer<i8>,
    stdout: SourceMutPointer<crate::source_types::FILE>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3062 OAD_ProcessOneStructure105Plus
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_ProcessOneStructure105Plus( struct tagINCHI_CLOCK    *ic,
                                    struct tagCANON_GLOBALS  *CG,
                                    STRUCT_DATA              *sd,
                                    INPUT_PARMS              *ip,
                                    char                     *szTitle,
                                    PINChI2                  *pINChI2[INCHI_NUM],
                                    PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                                    INCHI_IOSTREAM           *inp_file,
                                    INCHI_IOSTREAM           *log_file,
                                    INCHI_IOSTREAM           *out_file,
                                    INCHI_IOSTREAM           *prb_file,
                                    ORIG_ATOM_DATA           *orig_inp_data,
                                    ORIG_ATOM_DATA           *prep_inp_data,
                                    long                     num_inp,
                                    INCHI_IOS_STRING         *strbuf,
                                    unsigned char            save_opt_bits,
                                    char					 **sinchi,
                                    char					 **saux)
{
    size_t slen;
    int ret = _IS_OKAY, dup_fail = 0;
    POSEContext dup_context, *dup = &dup_context;

    *sinchi = NULL;
    *saux = NULL;

    /* Make a working copy of all the native input */
    dup_fail = POSEContext_Init(dup, sd, ip, szTitle, pINChI2, pINChI_Aux2,
                                inp_file, log_file, out_file, prb_file,
                                orig_inp_data, prep_inp_data,
                                num_inp, strbuf, save_opt_bits);
    if (dup_fail)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    /* Set necessary for this specific case options */
    /* 1.05+ polymer treatment mode: hidden Zz, senior bkbond comes the first */
    dup->orig_inp_data->polymer->treat = dup->ip.bPolymers = POLYMERS_LEGACY_PLUS;
    /* include full aux info */
    dup->ip.bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO);
    /* request no /D display in inchi-1 executable */
    dup->ip.bDisplay = dup->ip.bDisplayCompositeResults = dup->ip.bDisplayEachComponentINChI = 0;

    /* Calculate */
    ret = ProcessOneStructureExCore(ic, CG, &dup->sd, &dup->ip, dup->szTitle,
                                    dup->pINChI2, dup->pINChI_Aux2,
                                    dup->inp_file, dup->log_file,
                                    dup->out_file, dup->prb_file,
                                    dup->orig_inp_data, dup->prep_inp_data,
                                    dup->num_inp, dup->strbuf, dup->save_opt_bits);

    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }

    /* Extract InChI */
    slen = dup->out_file->s.nUsedLength;
    extract_inchi_substring(sinchi, dup->out_file->s.pStr, slen);
    if (!*sinchi)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }
    /* Extract AuxInfo */
    slen = dup->out_file->s.nUsedLength;
    extract_auxinfo_substring(saux, dup->out_file->s.pStr, slen);
    if (!*saux)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }


exit_function:
    POSEContext_Free(dup);

    return ret;
}
    */
    // END INCHI C FUNCTION: OAD_ProcessOneStructure105Plus
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_ProcessOneStructure105Plus
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; POLYMERS_LEGACY_PLUS=3.
    // INCHI✔️❌: Explicit stdout and clock_result model active globals passed through the worker.
    // INCHI✔️❌: Cloned stream and strbuf records add copying and allocation-map lookups absent in C.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_ProcessOneStructure105Plus

    *inchi_output = SourceMutPointer::null();
    *aux_output = SourceMutPointer::null();

    let mut duplicate = crate::source_types::POSEContext::default();
    let execution = (|| -> Result<i32, SourceHeapError> {
        let duplicate_failure = POSEContext_Init(
            heap, &mut duplicate, Some(structure_data), Some(input_parameters), title,
            inchi_components, aux_components, input_file, log_file, output_file, problem_file,
            Some(original_input), Some(prepared_input), input_number, string_buffer,
            save_option_bits,
        )?;
        if duplicate_failure != 0 {
            return Ok(_IS_ERROR as i32);
        }

        let polymer = heap.slice(duplicate.orig_inp_data.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?.polymer;
        heap.slice_mut(polymer)?.first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?.treat = POLYMERS_LEGACY_PLUS as i32;
        duplicate.ip.bPolymers = POLYMERS_LEGACY_PLUS as i32;
        duplicate.ip.bINChIOutputOptions &=
            !((INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO) as i32);
        duplicate.ip.bDisplay = 0;
        duplicate.ip.bDisplayCompositeResults = 0;
        duplicate.ip.bDisplayEachComponentINChI = 0;

        let mut input_stream = if duplicate.inp_file.is_null() {
            None
        } else {
            Some(heap.slice(duplicate.inp_file.as_const())?
                .first().ok_or(SourceHeapError::PointerOutOfBounds)?.clone())
        };
        let stored_streams = heap.slice(duplicate.out_file.as_const())?;
        let mut streams: [INCHI_IOSTREAM; 3] =
            std::array::from_fn(|index| stored_streams[index].clone());
        let mut string = heap.slice(duplicate.strbuf.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?.clone();
        let (output_streams, remainder) = streams.split_at_mut(1);
        let (log_streams, problem_streams) = remainder.split_at_mut(1);
        let worker_result = ProcessOneStructureExCore(
            heap, clock, canonical_globals, &mut duplicate.sd, &mut duplicate.ip,
            Some(duplicate.szTitle.as_mut_slice()), &mut duplicate.pINChI2,
            &mut duplicate.pINChI_Aux2, input_stream.as_mut(), log_streams.first_mut(),
            output_streams.first_mut(), problem_streams.first_mut(), duplicate.orig_inp_data,
            duplicate.prep_inp_data, duplicate.num_inp, Some(&mut string),
            duplicate.save_opt_bits, stdout, clock_result,
        );

        if let Some(input_stream) = input_stream {
            *heap.slice_mut(duplicate.inp_file)?.first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = input_stream;
        }
        heap.slice_mut(duplicate.out_file)?[..3].clone_from_slice(&streams);
        duplicate.inchi_file.clone_from_slice(&streams);
        *heap.slice_mut(duplicate.strbuf)?.first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = string.clone();
        duplicate.temp_string_container = string;

        let ret = worker_result?;
        if ret == _IS_FATAL as i32 || ret == _IS_ERROR as i32 {
            return Ok(ret);
        }

        let output_stream = heap.slice(duplicate.out_file.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let string_length = u64::try_from(output_stream.s.nUsedLength)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        extract_inchi_substring(heap, inchi_output, output_stream.s.pStr.as_const(), string_length)?;
        if inchi_output.is_null() {
            return Ok(_IS_ERROR as i32);
        }

        let output_stream = heap.slice(duplicate.out_file.as_const())?
            .first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let string_length = u64::try_from(output_stream.s.nUsedLength)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        extract_auxinfo_substring(heap, aux_output, output_stream.s.pStr.as_const(), string_length)?;
        if aux_output.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        Ok(ret)
    })();

    let cleanup = POSEContext_Free(heap, &mut duplicate);
    match (execution, cleanup) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(ret), Ok(())) => Ok(ret),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn mark_atoms_to_delete_or_renumber(
    heap: &mut SourceHeap,
    original_atom_data: &ORIG_ATOM_DATA,
    edits: &OAD_StructureEdits,
    atom_renumbering: SourceMutPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3155 mark_atoms_to_delete_or_renumber
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int mark_atoms_to_delete_or_renumber( ORIG_ATOM_DATA *orig_at_data,
                                      OAD_StructureEdits *ed,
                                      int *at_renum)
{
    int i, j;
    int fail=0, ret = 0;
    size_t *atnums = NULL; /* djb-rwth: needs to be size_t type */
    size_t max_atoms = orig_at_data->num_inp_atoms;

    /* NB:	new/old ORIG_ATOM_DATA atom numbers are 0-based (==orig_number-1) 
            while those in ed->... are just 1-based orig_numbers */
    
    for (i = 0; (size_t)i < max_atoms; i++)
    {
        at_renum[i] = i;
    }

    if (ed->del_side_chains)
    {
        /* Extend list of atoms to be deleted with those connected to original ones
        (i.e., delete a whole connected component(s) comprising original atoms)
        */
        int natnums = 0;
        atnums = (size_t *)inchi_calloc(max_atoms, sizeof(size_t)); /* djb-rwth: size_t type used for max_atoms to fit the definition of inchi_calloc  */
        if (!atnums)
        {
            return _IS_ERROR;
        }
        for (i = 0; (size_t)i < max_atoms; i++)
        {
            int iatom = ed->del_atom->item[i] - 1;
            subgraf *sg = NULL;
            subgraf_pathfinder *spf = NULL;
            if (i >= ed->del_atom->used) /* NB: ed->del_atom->used may be increased within this loop */
            {
                break;
            }
            for (j = 0; (size_t)j < max_atoms; j++)
            {
                atnums[j] = orig_at_data->at[j].orig_at_number; /*j+1;*/
            }
            sg = subgraf_new(orig_at_data, max_atoms, (int*)atnums);
            memset(atnums, 0, max_atoms * sizeof(int)); /* djb-rwth: memset_s C11/Annex K variant? */
            if (!sg)
            {
                ret = _IS_ERROR;
                goto exit_function;
            }
            spf = subgraf_pathfinder_new(sg, orig_at_data, iatom, iatom);
            if (!spf)
            {
                ret = _IS_ERROR;
                goto exit_function;
            }
            spf->start = iatom;
            spf->nseen = 0;
            natnums = subgraf_pathfinder_collect_all(spf, 0, NULL, (int*)atnums);
            if (natnums)
            {
                for (j = 0; j < natnums && j < max_atoms; j++) /* djb-rwth: fixing buffer overruns */
                {
                    fail = IntArray_AppendIfAbsent(ed->del_atom, atnums[j]);
                    if (fail)
                    {
                        ret = _IS_ERROR;
                        goto exit_function;
                    }
                }
            }
            subgraf_free(sg);
            subgraf_pathfinder_free(spf);
        }

    } /* if (ed->del_side_chains) */

    for (i = max_atoms - 1; i >= 0; i--)
    {
        int orig_num = i + 1;	/* NB: ed->del_atom->item contains orig# which are (OAD# + 1) */
        if (is_in_the_ilist(ed->del_atom->item, orig_num, ed->del_atom->used)) 
        {
            /* mark as deleted atnum */
            at_renum[i] = -1;
            /* shift other atnums */
            for (j = max_atoms - 1; j > i; j--)
            {
                if (at_renum[j] != -1)
                {
                    at_renum[j]--;
                }
            }
        }
    }


exit_function:
    if (atnums)
    {
        inchi_free(atnums);
    }
    return ret;
}
    */
    // END INCHI C FUNCTION: mark_atoms_to_delete_or_renumber
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: mark_atoms_to_delete_or_renumber
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4 and sizeof(size_t)=8.
    // INCHI✔️❌: With del_side_chains and max_atoms > 1, the active (int *)atnums cast makes nodes[1]
    // INCHI✔️❌: the zero high half of atnums[0], and subgraf_new then evaluates at[-1]; that C domain is undefined.
    // INCHI✔️❌: Rust returns structured UnsupportedSourceBehavior for that undefined domain instead of inspecting heap-prefix bytes.
    // END INCHI ACTIVE MACRO CONFIGURATION: mark_atoms_to_delete_or_renumber

    fn pointed_value<T: Clone + 'static>(
        heap: &SourceHeap,
        pointer: SourceMutPointer<T>,
    ) -> Result<T, SourceHeapError> {
        heap.slice(pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn is_allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }

    let max_atoms = usize::try_from(original_atom_data.num_inp_atoms)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    if max_atoms != 0 {
        let renumbering = heap.slice_mut(atom_renumbering)?;
        let initialized = renumbering
            .get_mut(..max_atoms)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for (index, value) in initialized.iter_mut().enumerate() {
            *value = i32::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        }
    }

    if edits.del_side_chains != 0 {
        let atom_numbers = match inchi_calloc::<i32>(heap, max_atoms as u64, 8) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };
        if max_atoms > 1 {
            inchi_free(heap, atom_numbers)?;
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }

        let expanded = (|| -> Result<i32, SourceHeapError> {
            let mut index = 0_usize;
            while index < max_atoms {
                let deletion = pointed_value(heap, edits.del_atom)?;
                let atom_to_delete = heap
                    .slice(deletion.item.as_const())?
                    .get(index)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                if i32::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    >= deletion.used
                {
                    break;
                }

                let atoms = heap.slice(original_atom_data.at.as_const())?.to_vec();
                let output = heap.slice_mut(atom_numbers)?;
                for atom_index in 0..max_atoms {
                    output[atom_index] = i32::from(
                        atoms
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .orig_at_number,
                    );
                }

                let graph_nodes = heap.slice(atom_numbers.as_const())?.to_vec();
                let graph = subgraf_new(
                    heap,
                    original_atom_data,
                    original_atom_data.num_inp_atoms,
                    &graph_nodes,
                )?;
                heap.slice_mut(atom_numbers)?.fill(0);
                if graph.is_null() {
                    return Ok(_IS_ERROR as i32);
                }
                let pathfinder = subgraf_pathfinder_new(
                    heap,
                    graph,
                    original_atom_data,
                    atom_to_delete,
                    atom_to_delete,
                )?;
                if pathfinder.is_null() {
                    return Ok(_IS_ERROR as i32);
                }
                {
                    let state = &mut heap.slice_mut(pathfinder)?[0];
                    state.start = atom_to_delete;
                    state.nseen = 0;
                }
                let collected = subgraf_pathfinder_collect_all(
                    heap,
                    pathfinder,
                    0,
                    SourceMutPointer::null(),
                    atom_numbers,
                )?;
                if collected != 0 {
                    let append_count = usize::try_from(collected)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                        .min(max_atoms);
                    for append_index in 0..append_count {
                        let new_item = heap.slice(atom_numbers.as_const())?[append_index];
                        let mut deletion = pointed_value(heap, edits.del_atom)?;
                        let append_result =
                            IntArray_AppendIfAbsent(heap, &mut deletion, new_item);
                        heap.slice_mut(edits.del_atom)?[0] = deletion;
                        if append_result? != 0 {
                            return Ok(_IS_ERROR as i32);
                        }
                    }
                }
                subgraf_free(heap, graph)?;
                subgraf_pathfinder_free(heap, pathfinder)?;
                index += 1;
            }
            Ok(0)
        })();
        inchi_free(heap, atom_numbers)?;
        let expanded = expanded?;
        if expanded != 0 {
            return Ok(expanded);
        }
    }

    let deletion = pointed_value(heap, edits.del_atom)?;
    let deletion_values = if deletion.used == 0 {
        None
    } else {
        Some(heap.slice(deletion.item.as_const())?.to_vec())
    };
    for index in (0..max_atoms).rev() {
        let original_number = i32::try_from(index)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        if is_in_the_ilist(
            deletion_values.as_deref(),
            original_number,
            deletion.used,
        )?
        .is_some()
        {
            heap.slice_mut(atom_renumbering)?[index] = -1;
            for shifted_index in ((index + 1)..max_atoms).rev() {
                let value = &mut heap.slice_mut(atom_renumbering)?[shifted_index];
                if *value != -1 {
                    *value = value
                        .checked_sub(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ProcessOneStructure(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &mut INPUT_PARMS,
    mut title: Option<&mut [i8]>,
    inchi_components: &mut [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    aux_components: &mut [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    mut input_file: Option<&mut INCHI_IOSTREAM>,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    mut output_file: Option<&mut INCHI_IOSTREAM>,
    mut problem_file: Option<&mut INCHI_IOSTREAM>,
    original_input_pointer: SourceMutPointer<ORIG_ATOM_DATA>,
    prepared_input_pointer: SourceMutPointer<ORIG_ATOM_DATA>,
    input_number: i64,
    mut string_buffer: Option<&mut INCHI_IOS_STRING>,
    mut save_option_bits: u8,
    stdout: SourceMutPointer<crate::source_types::FILE>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:218 ProcessOneStructure
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: ProcessOneStructure
    // INCHI✔️❌: int ProcessOneStructure( INCHI_CLOCK            *ic,
    // INCHI✔️❌:                          CANON_GLOBALS          *pCG,
    // INCHI✔️❌:                          STRUCT_DATA            *sd,
    // INCHI✔️❌:                          INPUT_PARMS            *ip,
    // INCHI✔️❌:                          char                   *szTitle,
    // INCHI✔️❌:                          PINChI2                *pINChI[INCHI_NUM],
    // INCHI✔️❌:                          PINChI_Aux2            *pINChI_Aux[INCHI_NUM],
    // INCHI✔️❌:                          INCHI_IOSTREAM         *inp_file,
    // INCHI✔️❌:                          INCHI_IOSTREAM         *log_file,
    // INCHI✔️❌:                          INCHI_IOSTREAM         *out_file,
    // INCHI✔️❌:                          INCHI_IOSTREAM         *prb_file,
    // INCHI✔️❌:                          ORIG_ATOM_DATA         *orig_inp_data,
    // INCHI✔️❌:                          ORIG_ATOM_DATA         *prep_inp_data,
    // INCHI✔️❌:                          long                   num_inp,
    // INCHI✔️❌:                          INCHI_IOS_STRING       *strbuf,
    // INCHI✔️❌:                          unsigned char          save_opt_bits )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet = 0,
    // INCHI✔️❌:         nRet1, i, k,
    // INCHI✔️❌:         maxINChI = 0,
    // INCHI✔️❌:         bSortPrintINChIFlags = 0;
    // INCHI✔️❌:     COMP_ATOM_DATA
    // INCHI✔️❌:         composite_norm_data[INCHI_NUM][TAUT_NUM + 1];    /*    [0]:non-taut,
    // INCHI✔️❌:                                                         [1]:taut,
    // INCHI✔️❌:                                                         [2]:intermediate taut struct */
    // INCHI✔️❌:     NORM_CANON_FLAGS ncFlags;
    // INCHI✔️❌:     NORM_CANON_FLAGS *pncFlags = &ncFlags;
    // INCHI✔️❌:     ORIG_STRUCT OrigStruct;
    // INCHI✔️❌:     ORIG_STRUCT *pOrigStruct = NULL;
    // INCHI✔️❌:     int err, ret1 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #ifdef GHI100_FIX
    // INCHI✔️❌: #if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
    // INCHI✔️❌:     setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    1. Preliminary work */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: fixing coverity ID #499508 */
    // INCHI✔️❌:     if (!orig_inp_data)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     int is_polymer = orig_inp_data->valid_polymer
    // INCHI✔️❌:                      && orig_inp_data->polymer
    // INCHI✔️❌:                      && orig_inp_data->polymer->n ;
    // INCHI✔️❌:
    // INCHI✔️❌:     int is_polymer2inchi = is_polymer && ( ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE );
    // INCHI✔️❌:
    // INCHI✔️❌:     sd->bUserQuitComponent = 0;
    // INCHI✔️❌:     sd->bUserQuitComponentDisplay = 0;
    // INCHI✔️❌:     memset( composite_norm_data, 0, sizeof( composite_norm_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( pncFlags, 0, sizeof( *pncFlags ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For experimental purposes only */
    // INCHI✔️❌:     /*ret1 = DoOneStructureEarlyPreprocessing( num_inp, sd, ip, inp_file,
    // INCHI✔️❌:                                              log_file, out_file, prb_file,
    // INCHI✔️❌:                                              orig_inp_data, prep_inp_data );
    // INCHI✔️❌:     */
    // INCHI✔️❌:     ret1 = DoOneStructureEarlyPreprocessing( ic, pCG, num_inp, sd, ip,
    // INCHI✔️❌:                                              inp_file, log_file, out_file, prb_file,
    // INCHI✔️❌:                                              orig_inp_data, prep_inp_data );
    // INCHI✔️❌:     switch (ret1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         case _IS_SKIP:
    // INCHI✔️❌:         case _IS_ERROR:
    // INCHI✔️❌:         case _IS_FATAL:
    // INCHI✔️❌:             nRet = ret1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ret1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (is_polymer)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Polymer house-keeping related to possible CRU frame shift(s) */
    // INCHI✔️❌:
    // INCHI✔️❌:         orig_inp_data->polymer->frame_shift_scheme = ip->bFrameShiftScheme;
    // INCHI✔️❌:         orig_inp_data->polymer->treat = ip->bPolymers;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!is_polymer2inchi)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Polymer structure is being restored from InChI string            */
    // INCHI✔️❌:             /* If CRUs were pre-cyclized, re-open them in preferred forms here  */
    // INCHI✔️❌:             if (orig_inp_data->polymer->frame_shift_scheme == FSS_STARS_CYCLED)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 OAD_Polymer_SmartReopenCyclizedUnits( orig_inp_data->polymer,
    // INCHI✔️❌:                                                       orig_inp_data->at,
    // INCHI✔️❌:                                                       orig_inp_data->num_inp_atoms,
    // INCHI✔️❌:                                                       &orig_inp_data->num_inp_bonds );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     ret1 = OrigAtData_SaveMolfile( orig_inp_data, sd, ip, num_inp, out_file );
    // INCHI✔️❌:     if (ret1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     pOrigStruct = &OrigStruct;
    // INCHI✔️❌:     memset( pOrigStruct, 0, sizeof( *pOrigStruct ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     OrigAtData_StoreNativeInput( pCG, &nRet, sd,  ip,  orig_inp_data, pOrigStruct );
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    2. Create INChI for the whole disconnected or original structure */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nRet1 = CreateOneStructureINChI( pCG, ic, sd, ip, szTitle,
    // INCHI✔️❌:                                          pINChI, pINChI_Aux, INCHI_BAS,
    // INCHI✔️❌:                                          inp_file, log_file, out_file, prb_file,
    // INCHI✔️❌:                                          orig_inp_data, prep_inp_data,
    // INCHI✔️❌:                                          composite_norm_data,
    // INCHI✔️❌:                                          num_inp, strbuf, pncFlags );
    // INCHI✔️❌:         nRet = inchi_max( nRet, nRet1 );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* If we create InChI from polymer-containing structure */
    // INCHI✔️❌:         if (is_polymer2inchi)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int polymer_repr_type = OAD_Polymer_GetRepresentation( orig_inp_data->polymer );
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef ALLOW_MIXED_SRU_AND_MON
    // INCHI✔️❌:             if (polymer_repr_type == POLYMER_REPRESENTATION_STRUCTURE_BASED ||
    // INCHI✔️❌:                  polymer_repr_type == POLYMER_REPRESENTATION_MIXED)
    // INCHI✔️❌: #else
    // INCHI✔️❌:             if (polymer_repr_type == POLYMER_REPRESENTATION_STRUCTURE_BASED)
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Temporarily copy ptr to polymer data to prep_inp_data */
    // INCHI✔️❌:                 OAD_Polymer *prep_polymer = prep_inp_data->polymer; /* may be NULL */
    // INCHI✔️❌:                 prep_inp_data->polymer = orig_inp_data->polymer;
    // INCHI✔️❌:
    // INCHI✔️❌:                 OAD_Polymer_FindBackbones( prep_inp_data, /* NB: not orig_inp_data! */
    // INCHI✔️❌:                                            &( composite_norm_data[INCHI_BAS][TAUT_YES] ),
    // INCHI✔️❌:                                            &err, sd->pStrErrStruct );
    // INCHI✔️❌:                 if (err)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret1 = _IS_ERROR;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 nRet = inchi_max( nRet, ret1 );
    // INCHI✔️❌:                 prep_inp_data->polymer = prep_polymer;    /* restore temp copied*/
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         maxINChI = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* 3. Create INChI for the whole metal-reconnected structure */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nRet != _IS_FATAL                                              &&
    // INCHI✔️❌:          nRet != _IS_ERROR &&
    // INCHI✔️❌:          ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
    // INCHI✔️❌:          ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         nRet1 = CreateOneStructureINChI( pCG, ic, sd, ip, szTitle,
    // INCHI✔️❌:                                          pINChI, pINChI_Aux, INCHI_REC, inp_file,
    // INCHI✔️❌:                                          log_file, out_file, prb_file,
    // INCHI✔️❌:                                          orig_inp_data, prep_inp_data,
    // INCHI✔️❌:                                          composite_norm_data,
    // INCHI✔️❌:                                          num_inp, strbuf, pncFlags );
    // INCHI✔️❌:         nRet = inchi_max( nRet, nRet1 );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (is_polymer2inchi)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret1 = 0;
    // INCHI✔️❌:             /* temporarily copy ptr to polymer data to prep_inp_data */
    // INCHI✔️❌:             prep_inp_data->polymer = orig_inp_data->polymer;
    // INCHI✔️❌:
    // INCHI✔️❌:             OAD_Polymer_FindBackbones( prep_inp_data, /* NB: not orig_inp_data! */
    // INCHI✔️❌:                                        &( composite_norm_data[INCHI_REC][TAUT_YES] ),
    // INCHI✔️❌:                                        &err, sd->pStrErrStruct );
    // INCHI✔️❌:             if (err)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ret1 = _IS_ERROR;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             nRet = inchi_max( nRet, ret1 );
    // INCHI✔️❌:             prep_inp_data->polymer = NULL;    /* remove temp copied */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             maxINChI = 2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( sd->bChiralFlag & FLAG_INP_AT_CHIRAL ) &&
    // INCHI✔️❌:             ( ip->nMode & REQ_MODE_STEREO ) &&
    // INCHI✔️❌:              !( ip->nMode & ( REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO ) ) &&
    // INCHI✔️❌:              !bIsStructChiral( pINChI, sd->num_components ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!ip->bNoWarnings)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 WarningMessage( sd->pStrErrStruct, "Not chiral" );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!sd->bUserQuitComponent && !sd->bUserQuit)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nRet1 = TreatCreateINChIWarning( sd, ip, prep_inp_data, num_inp,
    // INCHI✔️❌:                                              inp_file, log_file, out_file, prb_file );
    // INCHI✔️❌:             nRet = inchi_max( nRet, nRet1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    4. Sort and print INChI for the whole structure */
    // INCHI✔️❌:
    // INCHI✔️❌:     PrepareSaveOptBits( &save_opt_bits, ip );
    // INCHI✔️❌:     if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nRet1 = SortAndPrintINChI( pCG, out_file, strbuf, log_file, ip,
    // INCHI✔️❌:                                    orig_inp_data, prep_inp_data,
    // INCHI✔️❌:                                    composite_norm_data,
    // INCHI✔️❌:                                    pOrigStruct, sd->num_components,
    // INCHI✔️❌:                                    sd->num_non_taut, sd->num_taut,
    // INCHI✔️❌:                                    sd->bTautFlags, sd->bTautFlagsDone,
    // INCHI✔️❌:                                    pncFlags, num_inp,
    // INCHI✔️❌:                                    pINChI, pINChI_Aux,
    // INCHI✔️❌:                                    &bSortPrintINChIFlags, save_opt_bits ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    5. Post-process */
    // INCHI✔️❌:
    // INCHI✔️❌:     DisplayOrigAndResultStructuresAndComponents( nRet, ic, pCG, sd, ip, szTitle,
    // INCHI✔️❌:                                                  pINChI, pINChI_Aux,
    // INCHI✔️❌:                                                  inp_file, log_file, out_file,
    // INCHI✔️❌:                                                  orig_inp_data, prep_inp_data,
    // INCHI✔️❌:                                                  num_inp, maxINChI,
    // INCHI✔️❌:                                                  composite_norm_data );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     SaveOkProcessedMolfile( nRet, sd, ip, prb_file, inp_file );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Cleanup */
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < INCHI_NUM; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (k = 0; k < TAUT_NUM + 1; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             FreeCompAtomData( &composite_norm_data[i][k] );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     OrigStruct_Free( pOrigStruct );
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: ProcessOneStructure
    // END INCHI C FUNCTION: ProcessOneStructure
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: ProcessOneStructure
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; GHI100_FIX is undefined.
    // INCHI✔️❌: ALLOW_MIXED_SRU_AND_MON=1; INCHI_NUM=2; TAUT_NUM=2.
    // INCHI✔️❌: SourceHeap borrowing requires temporary struct snapshots and model-storage stack slots.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: ProcessOneStructure

    if original_input_pointer.is_null() {
        return Ok(0);
    }

    let mut canonical = heap
        .slice(canonical_globals.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let mut original_input = heap
        .slice(original_input_pointer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let mut prepared_input = heap.slice(prepared_input_pointer.as_const())?.to_vec();

    let result = (|| -> Result<i32, SourceHeapError> {
        let mut result = 0_i32;
        let mut max_inchi = 0_i32;
        let mut composite_normalized_data: [[COMP_ATOM_DATA; TAUT_NUM as usize + 1];
            INCHI_NUM as usize] =
            std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
        let mut normalization_flags = NORM_CANON_FLAGS::default();
        let mut original_structure = ORIG_STRUCT::default();

        let is_polymer = if original_input.valid_polymer != 0 && !original_input.polymer.is_null() {
            heap.slice(original_input.polymer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .n
                != 0
        } else {
            false
        };
        let is_polymer_to_inchi = is_polymer
            && (input_parameters.nInputType == tagInputType_INPUT_MOLFILE
                || input_parameters.nInputType == tagInputType_INPUT_SDFILE);

        structure_data.bUserQuitComponent = 0;
        structure_data.bUserQuitComponentDisplay = 0;
        let preprocessing_result = DoOneStructureEarlyPreprocessing(
            heap,
            clock,
            &mut canonical,
            input_number,
            structure_data,
            input_parameters,
            input_file.as_deref_mut(),
            log_file.as_deref_mut(),
            output_file.as_deref_mut(),
            problem_file.as_deref_mut(),
            &mut original_input,
            prepared_input
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
            clock_result,
        )?;
        if preprocessing_result == crate::source_types::_IS_SKIP as i32
            || preprocessing_result == _IS_ERROR as i32
            || preprocessing_result == _IS_FATAL as i32
        {
            result = preprocessing_result;
        }
        if preprocessing_result != 0 {
            return Ok(result);
        }

        if is_polymer {
            {
                let polymer = heap
                    .slice_mut(original_input.polymer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                polymer.frame_shift_scheme = input_parameters.bFrameShiftScheme;
                polymer.treat = input_parameters.bPolymers;
            }
            if !is_polymer_to_inchi
                && heap.slice(original_input.polymer.as_const())?[0].frame_shift_scheme
                    == tagFrameShifScheme_FSS_STARS_CYCLED as i32
            {
                OAD_Polymer_SmartReopenCyclizedUnits(
                    heap,
                    original_input.polymer,
                    original_input.at,
                    original_input.num_inp_atoms,
                    &mut original_input.num_inp_bonds,
                )?;
            }
        }

        let mut fallback_output = INCHI_IOSTREAM::default();
        let save_result = OrigAtData_SaveMolfile(
            heap,
            &original_input,
            structure_data,
            input_parameters,
            input_number,
            output_file.as_deref_mut().unwrap_or(&mut fallback_output),
        )?;
        if save_result != 0 {
            return Ok(result);
        }

        OrigAtData_StoreNativeInput(
            heap,
            &mut canonical,
            &mut result,
            structure_data,
            input_parameters,
            &mut original_input,
            &mut original_structure,
        )?;

        if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
            let created = CreateOneStructureINChI(
                heap,
                &mut canonical,
                clock,
                structure_data,
                input_parameters,
                title.as_deref_mut(),
                inchi_components,
                aux_components,
                INCHI_BAS as i32,
                input_file.as_deref_mut(),
                log_file.as_deref_mut(),
                output_file.as_deref_mut(),
                problem_file.as_deref_mut(),
                &mut original_input,
                &mut prepared_input,
                &mut composite_normalized_data,
                input_number,
                string_buffer.as_deref_mut(),
                &mut normalization_flags,
                clock_result,
            )?;
            result = result.max(created);
            if is_polymer_to_inchi {
                let representation = OAD_Polymer_GetRepresentation(heap, original_input.polymer)?;
                if representation == POLYMER_REPRESENTATION_STRUCTURE_BASED as i32
                    || representation == POLYMER_REPRESENTATION_MIXED as i32
                {
                    let prepared = prepared_input
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let saved_polymer = prepared.polymer;
                    prepared.polymer = original_input.polymer;
                    let mut error = 0_i32;
                    OAD_Polymer_FindBackbones(
                        heap,
                        prepared,
                        Some(&composite_normalized_data[INCHI_BAS as usize][TAUT_YES as usize]),
                        &mut error,
                        Some(&mut structure_data.pStrErrStruct),
                    )?;
                    if error != 0 {
                        result = result.max(_IS_ERROR as i32);
                    }
                    prepared.polymer = saved_polymer;
                }
            }
        }
        if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
            max_inchi = 1;
        }

        if result != _IS_FATAL as i32
            && result != _IS_ERROR as i32
            && structure_data.bTautFlagsDone[INCHI_BAS as usize]
                & u64::from(TG_FLAG_DISCONNECT_COORD_DONE)
                != 0
            && input_parameters.bTautFlags & u64::from(TG_FLAG_RECONNECT_COORD) != 0
        {
            let created = CreateOneStructureINChI(
                heap,
                &mut canonical,
                clock,
                structure_data,
                input_parameters,
                title.as_deref_mut(),
                inchi_components,
                aux_components,
                INCHI_REC as i32,
                input_file.as_deref_mut(),
                log_file.as_deref_mut(),
                output_file.as_deref_mut(),
                problem_file.as_deref_mut(),
                &mut original_input,
                &mut prepared_input,
                &mut composite_normalized_data,
                input_number,
                string_buffer.as_deref_mut(),
                &mut normalization_flags,
                clock_result,
            )?;
            result = result.max(created);
            if is_polymer_to_inchi {
                let prepared = prepared_input
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                prepared.polymer = original_input.polymer;
                let mut error = 0_i32;
                OAD_Polymer_FindBackbones(
                    heap,
                    prepared,
                    Some(&composite_normalized_data[INCHI_REC as usize][TAUT_YES as usize]),
                    &mut error,
                    Some(&mut structure_data.pStrErrStruct),
                )?;
                if error != 0 {
                    result = result.max(_IS_ERROR as i32);
                }
                prepared.polymer = SourceMutPointer::null();
            }
            if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
                max_inchi = 2;
            }
        }

        if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
            if structure_data.bChiralFlag & FLAG_INP_AT_CHIRAL as i32 != 0
                && input_parameters.nMode & u64::from(REQ_MODE_STEREO) != 0
                && input_parameters.nMode
                    & u64::from(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO)
                    == 0
                && bIsStructChiral(heap, *inchi_components, structure_data.num_components)? == 0
                && input_parameters.bNoWarnings == 0
            {
                let message = b"Not chiral\0".map(|byte| byte as i8);
                AddErrorMessage(Some(&mut structure_data.pStrErrStruct), Some(&message))?;
            }
            if structure_data.bUserQuitComponent == 0 && structure_data.bUserQuit == 0 {
                let warning = TreatCreateINChIWarning(
                    heap,
                    structure_data,
                    input_parameters,
                    prepared_input
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    input_number,
                    input_file.as_deref_mut(),
                    log_file.as_deref_mut(),
                    output_file.as_deref_mut(),
                    problem_file.as_deref_mut(),
                )?;
                result = result.max(warning);
            }
        }

        PrepareSaveOptBitsForStructure(&mut save_option_bits, input_parameters);
        if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
            heap.slice_mut(canonical_globals)?[0] = canonical.clone();
            heap.slice_mut(original_input_pointer)?[0] = original_input.clone();
            heap.slice_mut(prepared_input_pointer)?[..prepared_input.len()]
                .clone_from_slice(&prepared_input);
            let original_structure_slot =
                heap.allocate_model_storage(vec![original_structure.clone()])?;
            let sort_flags = heap.allocate_model_storage(vec![0_i32])?;
            let mut fallback_output = INCHI_IOSTREAM::default();
            let mut fallback_log = INCHI_IOSTREAM::default();
            let sort_result = SortAndPrintINChI(
                heap,
                canonical_globals,
                output_file.as_deref_mut().unwrap_or(&mut fallback_output),
                string_buffer.as_deref_mut(),
                log_file.as_deref_mut().unwrap_or(&mut fallback_log),
                input_parameters,
                original_input_pointer.as_const(),
                prepared_input_pointer.as_const(),
                Some(&mut composite_normalized_data),
                original_structure_slot.as_const(),
                &structure_data.num_components,
                &structure_data.num_non_taut,
                &structure_data.num_taut,
                &mut structure_data.bTautFlags,
                &mut structure_data.bTautFlagsDone,
                &normalization_flags,
                input_number,
                *inchi_components,
                *aux_components,
                sort_flags,
                save_option_bits,
                stdout,
            );
            heap.free(sort_flags)?;
            heap.free(original_structure_slot)?;
            let _ = sort_result?;
            canonical = heap.slice(canonical_globals.as_const())?[0].clone();
        }

        let _ = max_inchi;
        DisplayOrigAndResultStructuresAndComponents(input_parameters);
        SaveOkProcessedMolfile(
            heap,
            result,
            structure_data,
            input_parameters,
            problem_file.as_deref(),
            input_file.as_deref_mut(),
        )?;
        for domain in &mut composite_normalized_data {
            for representation in domain {
                FreeCompAtomData(heap, representation)?;
            }
        }
        OrigStruct_Free(heap, Some(&mut original_structure))?;
        Ok(result)
    })();

    heap.slice_mut(canonical_globals)?[0] = canonical;
    heap.slice_mut(original_input_pointer)?[0] = original_input;
    let stored_prepared = heap.slice_mut(prepared_input_pointer)?;
    if stored_prepared.len() < prepared_input.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    stored_prepared[..prepared_input.len()].clone_from_slice(&prepared_input);
    result
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_SaveMolfile(
    heap: &mut SourceHeap,
    original_input: &ORIG_ATOM_DATA,
    structure_data: &STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    input_number: i64,
    output: &mut INCHI_IOSTREAM,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:564 OrigAtData_SaveMolfile
    // INCHI✔️❌: int OrigAtData_SaveMolfile( ORIG_ATOM_DATA  *orig_inp_data,
    // INCHI✔️❌:                             STRUCT_DATA     *sd,
    // INCHI✔️❌:                             INPUT_PARMS     *ip,
    // INCHI✔️❌:                             long            num_inp,
    // INCHI✔️❌:                             INCHI_IOSTREAM  *out_file )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_OKAY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char szNumber[256];
    // INCHI✔️❌:         sprintf(szNumber, "Structure #%ld. %s%s%s%s", num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));
    // INCHI✔️❌:         ret = OrigAtData_WriteToSDfile( orig_inp_data, out_file, szNumber, NULL,
    // INCHI✔️❌:                                         ( sd->bChiralFlag&FLAG_INP_AT_CHIRAL ) ? 1 : 0,
    // INCHI✔️❌:                                         ( ip->bINChIOutputOptions&INCHI_OUT_SDFILE_ATOMS_DT ) ? 1 : 0,
    // INCHI✔️❌:                                         ip->pSdfLabel, ip->pSdfValue );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_SaveMolfile
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_SaveMolfile
    // INCHI✔️❌: #define SDF_LBL_VAL(L, V) ((L) && (L)[0]) ? gsSpace : gsEmpty, ((L) && (L)[0]) ? L : gsEmpty, ((L) && (L)[0]) ? (((V) && (V)[0]) ? gsEqual : gsSpace) : gsEmpty, ((V) && (V)[0]) ? V : ((L) && (L)[0]) ? gsMissing : gsEmpty
    // INCHI✔️❌: const char gsMissing[] = "is missing"; const char gsEmpty[] = "";
    // INCHI✔️❌: const char gsSpace[] = " "; const char gsEqual[] = "=";
    // INCHI✔️❌: #define INCHI_OUT_SDFILE_ONLY 0x0010
    // INCHI✔️❌: #define INCHI_OUT_SDFILE_ATOMS_DT 0x0800
    // INCHI✔️❌: #define FLAG_INP_AT_CHIRAL 1
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(long) == 8
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_SaveMolfile

    if input_parameters.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 == 0 {
        return Ok(_IS_OKAY as i32);
    }

    let label = runichi_c_text(heap, input_parameters.pSdfLabel.as_const())?;
    let value = runichi_c_text(heap, input_parameters.pSdfValue.as_const())?;
    let mut title = format!("Structure #{input_number}. ").into_bytes();
    if !label.is_empty() {
        title.push(b' ');
        title.extend_from_slice(&label);
        if value.is_empty() {
            title.push(b' ');
            title.extend_from_slice(b"is missing");
        } else {
            title.push(b'=');
            title.extend_from_slice(&value);
        }
    } else if !value.is_empty() {
        title.extend_from_slice(&value);
    }
    if title.len() >= 256 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let title = heap.allocate_model_storage(
        title
            .iter()
            .copied()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let result = OrigAtData_WriteToSDfile(
        heap,
        original_input,
        Some(output),
        SourceMutPointer::null(),
        title.as_const(),
        SourceConstPointer::null(),
        i32::from(structure_data.bChiralFlag & FLAG_INP_AT_CHIRAL as i32 != 0),
        i32::from(input_parameters.bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT as i32 != 0),
        input_parameters.pSdfLabel.as_const(),
        input_parameters.pSdfValue.as_const(),
    );
    heap.free(title)?;
    result
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_StoreNativeInput<'a>(
    heap: &mut SourceHeap,
    canon_globals: &mut CANON_GLOBALS,
    return_code: &mut i32,
    structure_data: &mut STRUCT_DATA,
    _input_parameters: &INPUT_PARMS,
    original_atom_data: &mut ORIG_ATOM_DATA,
    original_structure: &'a mut ORIG_STRUCT,
) -> Result<&'a mut ORIG_STRUCT, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:593 OrigAtData_StoreNativeInput
    // INCHI✔️✔️: ORIG_STRUCT * OrigAtData_StoreNativeInput( CANON_GLOBALS    *pCG,
    // INCHI✔️✔️:                                            int              *nRet,
    // INCHI✔️✔️:                                            STRUCT_DATA      *sd,
    // INCHI✔️✔️:                                            INPUT_PARMS      *ip,
    // INCHI✔️✔️:                                            ORIG_ATOM_DATA   *orig_inp_data,
    // INCHI✔️✔️:                                            ORIG_STRUCT      *pOrigStruct )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*    v. 1.05 always create and fill OrigStruc as it may be used to store e.g. polymer info    */
    // INCHI✔️✔️:     /*    If normal AuxInfo is requested, create full reversibility information from native inp data
    // INCHI✔️✔️:     if ( ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO))
    // INCHI✔️✔️:         return NULL; */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (OrigStruct_FillOut( pCG, orig_inp_data, pOrigStruct, sd ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         AddErrorMessage( sd->pStrErrStruct, "Cannot interpret reversibility information" );
    // INCHI✔️✔️:         sd->nStructReadError = 99;
    // INCHI✔️✔️:         sd->nErrorType = _IS_ERROR;
    // INCHI✔️✔️:         *nRet = _IS_ERROR;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return pOrigStruct;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_StoreNativeInput
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_StoreNativeInput
    // INCHI✔️✔️: #define _IS_ERROR 2
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_StoreNativeInput

    if OrigStruct_FillOut(
        heap,
        canon_globals,
        original_atom_data,
        original_structure,
        structure_data,
    )? != 0
    {
        let message = b"Cannot interpret reversibility information\0".map(|byte| byte as i8);
        AddErrorMessage(Some(&mut structure_data.pStrErrStruct), Some(&message))?;
        structure_data.nStructReadError = 99;
        structure_data.nErrorType = _IS_ERROR as i32;
        *return_code = _IS_ERROR as i32;
    }

    Ok(original_structure)
}

#[allow(non_snake_case)]
pub(crate) fn PrepareSaveOptBitsForStructure(
    save_option_bits: &mut u8,
    input_parameters: &INPUT_PARMS,
) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:621 PrepareSaveOptBits
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: PrepareSaveOptBits
    // INCHI✔️✔️: void PrepareSaveOptBits( unsigned char *save_opt_bits, INPUT_PARMS *ip )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (ip->nInputType != INPUT_INCHI)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *save_opt_bits = 0;
    // INCHI✔️✔️:         if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ( *save_opt_bits ) |= SAVE_OPT_RECMET;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (0 != ( ip->nMode & REQ_MODE_BASIC ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ( *save_opt_bits ) |= SAVE_OPT_FIXEDH;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ( *save_opt_bits ) |= SAVE_OPT_SLUUD;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (0 == ( ip->nMode & ( REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU ) ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ( *save_opt_bits ) |= SAVE_OPT_SUU;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (0 != ( ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ( *save_opt_bits ) |= SAVE_OPT_KET;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (0 != ( ip->bTautFlags & TG_FLAG_1_5_TAUT ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ( *save_opt_bits ) |= SAVE_OPT_15T;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             /* djb-rwth: addressing coverity ID #499536 -- despite different bit-sizes, works properly */
    // INCHI✔️✔️:             if (0 != (ip->bTautFlags & TG_FLAG_PT_22_00))
    // INCHI✔️✔️:                 (*save_opt_bits) |= SAVE_OPT_PT_22_00;
    // INCHI✔️✔️:             if (0 != (ip->bTautFlags & TG_FLAG_PT_16_00))
    // INCHI✔️✔️:                 (*save_opt_bits) |= SAVE_OPT_PT_16_00;
    // INCHI✔️✔️:             if (0 != (ip->bTautFlags & TG_FLAG_PT_06_00))
    // INCHI✔️✔️:                 (*save_opt_bits) |= SAVE_OPT_PT_06_00;
    // INCHI✔️✔️:             if (0 != (ip->bTautFlags & TG_FLAG_PT_39_00))
    // INCHI✔️✔️:                 (*save_opt_bits) |= SAVE_OPT_PT_39_00;
    // INCHI✔️✔️:             if (0 != (ip->bTautFlags & TG_FLAG_PT_13_00))
    // INCHI✔️✔️:                 (*save_opt_bits) |= SAVE_OPT_PT_13_00;
    // INCHI✔️✔️:             if (0 != (ip->bTautFlags & TG_FLAG_PT_18_00))
    // INCHI✔️✔️:                 (*save_opt_bits) |= SAVE_OPT_PT_18_00;
    // INCHI✔️✔️:             /* Check if /SNon requested and turn OFF stereo bits if so */
    // INCHI✔️✔️:             if (!( ip->nMode & REQ_MODE_STEREO ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 ( *save_opt_bits ) &= ~SAVE_OPT_SUU;
    // INCHI✔️✔️:                 ( *save_opt_bits ) &= ~SAVE_OPT_SLUUD;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: PrepareSaveOptBits
    // END INCHI C FUNCTION: PrepareSaveOptBits
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: PrepareSaveOptBits
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; unsigned char is 8-bit.
    // INCHI✔️✔️: SAVE_OPT values 0x100 and above are truncated on assignment to unsigned char.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: PrepareSaveOptBits

    if input_parameters.nInputType != tagInputType_INPUT_INCHI {
        *save_option_bits = 0;
        if input_parameters.bINChIOutputOptions & INCHI_OUT_SAVEOPT as i32 != 0 {
            if input_parameters.bTautFlags & u64::from(TG_FLAG_RECONNECT_COORD) != 0 {
                *save_option_bits |= SAVE_OPT_RECMET as u8;
            }
            if input_parameters.nMode & u64::from(REQ_MODE_BASIC) != 0 {
                *save_option_bits |= SAVE_OPT_FIXEDH as u8;
            }
            if input_parameters.nMode & u64::from(REQ_MODE_DIFF_UU_STEREO) != 0 {
                *save_option_bits |= SAVE_OPT_SLUUD as u8;
            }
            if input_parameters.nMode & u64::from(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)
                == 0
            {
                *save_option_bits |= SAVE_OPT_SUU as u8;
            }
            if input_parameters.bTautFlags & u64::from(TG_FLAG_KETO_ENOL_TAUT) != 0 {
                *save_option_bits |= SAVE_OPT_KET as u8;
            }
            if input_parameters.bTautFlags & u64::from(TG_FLAG_1_5_TAUT) != 0 {
                *save_option_bits |= SAVE_OPT_15T as u8;
            }
            for (taut_flag, save_flag) in [
                (TG_FLAG_PT_22_00, SAVE_OPT_PT_22_00),
                (TG_FLAG_PT_16_00, SAVE_OPT_PT_16_00),
                (TG_FLAG_PT_06_00, SAVE_OPT_PT_06_00),
                (TG_FLAG_PT_39_00, SAVE_OPT_PT_39_00),
                (TG_FLAG_PT_13_00, SAVE_OPT_PT_13_00),
                (TG_FLAG_PT_18_00, SAVE_OPT_PT_18_00),
            ] {
                if input_parameters.bTautFlags & u64::from(taut_flag) != 0 {
                    *save_option_bits |= save_flag as u8;
                }
            }
            if input_parameters.nMode & u64::from(REQ_MODE_STEREO) == 0 {
                *save_option_bits &= !(SAVE_OPT_SUU as u8);
                *save_option_bits &= !(SAVE_OPT_SLUUD as u8);
            }
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn DisplayOrigAndResultStructuresAndComponents(input_parameters: &mut INPUT_PARMS) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:679 DisplayOrigAndResultStructuresAndComponents
    // INCHI✔️✔️: void DisplayOrigAndResultStructuresAndComponents( int               nRet,
    // INCHI✔️✔️:                                                   INCHI_CLOCK       *ic,
    // INCHI✔️✔️:                                                   CANON_GLOBALS     *pCG,
    // INCHI✔️✔️:                                                   STRUCT_DATA       *sd,
    // INCHI✔️✔️:                                                   INPUT_PARMS       *ip,
    // INCHI✔️✔️:                                                   char              *szTitle,
    // INCHI✔️✔️:                                                   PINChI2           *pINChI[INCHI_NUM],
    // INCHI✔️✔️:                                                   PINChI_Aux2       *pINChI_Aux[INCHI_NUM],
    // INCHI✔️✔️:                                                   INCHI_IOSTREAM    *inp_file,
    // INCHI✔️✔️:                                                   INCHI_IOSTREAM    *log_file,
    // INCHI✔️✔️:                                                   INCHI_IOSTREAM    *out_file,
    // INCHI✔️✔️:                                                   ORIG_ATOM_DATA    *orig_inp_data,
    // INCHI✔️✔️:                                                   ORIG_ATOM_DATA    *prep_inp_data,
    // INCHI✔️✔️:                                                   long              num_inp,
    // INCHI✔️✔️:                                                   int               maxINChI,
    // INCHI✔️✔️:                                                   COMP_ATOM_DATA    composite_norm_data[INCHI_NUM][TAUT_NUM + 1] )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (ip->bDisplay)    ip->bDisplayCompositeResults = 1;    /* v. 1.05 */
    // INCHI✔️✔️:
    // INCHI❌❌: #ifndef COMPILE_ANSI_ONLY /* { */
    // INCHI❌❌:
    // INCHI❌❌:     /* Display equivalent components on original or preprocessed structure(s) */
    // INCHI❌❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI❌❌:     if (nRet != _IS_FATAL && nRet != _IS_ERROR && /*ip->bDisplay &&*/
    // INCHI❌❌:         ( ip->bCompareComponents & CMP_COMPONENTS ) && !sd->bUserQuit && !sd->bUserQuitComponent)
    // INCHI❌❌:     {
    // INCHI❌❌:         int j, ret, ord;
    // INCHI❌❌:         int bDisplaySaved = ip->bDisplay;
    // INCHI❌❌:         ORIG_ATOM_DATA *inp_data;
    // INCHI❌❌:         AT_NUMB         nEquSet;
    // INCHI❌❌:         for (ord = -1; ord < INCHI_NUM; ord++)
    // INCHI❌❌:         {
    // INCHI❌❌:             switch (ord)
    // INCHI❌❌:             {
    // INCHI❌❌:                 case -1:
    // INCHI❌❌:                     j = INCHI_BAS;  /* preprocessed non-tautomeric */
    // INCHI❌❌:                     break;
    // INCHI❌❌:                 case 0:
    // INCHI❌❌:                     j = INCHI_REC;  /* preprocessed tautomeric */
    // INCHI❌❌:                     break;
    // INCHI❌❌:                 case 1:
    // INCHI❌❌:                     j = -1;        /* original input */
    // INCHI❌❌:                     break;
    // INCHI❌❌:                 default:
    // INCHI❌❌:                     continue;
    // INCHI❌❌:             }
    // INCHI❌❌:             inp_data = j < 0 ? orig_inp_data : prep_inp_data + j;
    // INCHI❌❌:             if (inp_data && inp_data->num_inp_atoms && inp_data->at &&
    // INCHI❌❌:                  inp_data->nEquLabels &&
    // INCHI❌❌:                  inp_data->nNumEquSets)
    // INCHI❌❌:             {
    // INCHI❌❌:                 for (nEquSet = 1; nEquSet <= inp_data->nNumEquSets; nEquSet++)
    // INCHI❌❌:                 {
    // INCHI❌❌:                     ip->dp.nEquLabels = inp_data->nEquLabels;
    // INCHI❌❌:                     ip->dp.nCurEquLabel = nEquSet;
    // INCHI❌❌:                     ip->dp.nNumEquSets = inp_data->nNumEquSets;
    // INCHI❌❌:                     ip->bDisplay = 1; /* force display if it was not requested */
    // INCHI❌❌:                     ret = DisplayTheWholeStructure( pCG, ic, sd, ip, szTitle,
    // INCHI❌❌:                                                     inp_file, log_file, inp_data,
    // INCHI❌❌:                                                     num_inp, j, 1 /*bShowStructure*/, 0 );
    // INCHI❌❌:                     ip->dp.nEquLabels = NULL;
    // INCHI❌❌:                     ip->dp.nCurEquLabel = 0;
    // INCHI❌❌:                     ip->dp.nNumEquSets = 0;
    // INCHI❌❌:                     ip->bDisplay = bDisplaySaved; /* restore display option */
    // INCHI❌❌:                     if (ret)
    // INCHI❌❌:                     {
    // INCHI❌❌:             /* user pressed Esc */
    // INCHI❌❌:                         goto exit_loop;
    // INCHI❌❌:                     }
    // INCHI❌❌:                 }
    // INCHI❌❌:             }
    // INCHI❌❌:         }
    // INCHI❌❌:     exit_loop:;
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI❌❌:
    // INCHI❌❌:     /* Display composite results and equivalent components on composite results */
    // INCHI❌❌:     if (nRet != _IS_FATAL && nRet != _IS_ERROR && /*ip->bDisplay &&*/ ip->bDisplayCompositeResults)
    // INCHI❌❌:     {
    // INCHI❌❌:         int iINChI;
    // INCHI❌❌:         for (iINChI = 0; iINChI < maxINChI && !sd->bUserQuitComponentDisplay; iINChI++)
    // INCHI❌❌:         {
    // INCHI❌❌:             DisplayTheWholeCompositeStructure( pCG, ic, ip, sd, num_inp,
    // INCHI❌❌:                                                iINChI, pINChI[iINChI], pINChI_Aux[iINChI],
    // INCHI❌❌:                                                orig_inp_data, prep_inp_data, composite_norm_data[iINChI] );
    // INCHI❌❌:         }
    // INCHI❌❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI❌❌:         if (!ip->bDisplay && sd->bUserQuitComponentDisplay)
    // INCHI❌❌:         {
    // INCHI❌❌:             sd->bUserQuit = 1;
    // INCHI❌❌:         }
    // INCHI❌❌: #endif
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI✔️✔️: #endif /* } COMPILE_ANSI_ONLY */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: DisplayOrigAndResultStructuresAndComponents
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisplayOrigAndResultStructuresAndComponents
    // INCHI✔️✔️: #define COMPILE_ANSI_ONLY
    // INCHI✔️✔️: TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisplayOrigAndResultStructuresAndComponents

    if input_parameters.bDisplay != 0 {
        input_parameters.bDisplayCompositeResults = 1;
    }
}

#[allow(non_snake_case)]
pub(crate) fn SaveOkProcessedMolfile(
    heap: &mut SourceHeap,
    return_code: i32,
    structure_data: &STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    problem_file: Option<&INCHI_IOSTREAM>,
    input_file: Option<&mut INCHI_IOSTREAM>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:785 SaveOkProcessedMolfile
    // INCHI✔️❌: void SaveOkProcessedMolfile( int            nRet,
    // INCHI✔️❌:                              STRUCT_DATA    *sd,
    // INCHI✔️❌:                              INPUT_PARMS    *ip,
    // INCHI✔️❌:                              INCHI_IOSTREAM *prb_file,
    // INCHI✔️❌:                              INCHI_IOSTREAM *inp_file )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (ip->bSaveAllGoodStructsAsProblem &&
    // INCHI✔️❌:          nRet != _IS_FATAL                &&
    // INCHI✔️❌:          nRet != _IS_ERROR                &&
    // INCHI✔️❌:          prb_file                         &&
    // INCHI✔️❌:          prb_file->f &&
    // INCHI✔️❌:          0L <= sd->fPtrStart              &&
    // INCHI✔️❌:          sd->fPtrStart < sd->fPtrEnd)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, 0 ); /* djb-rwth: addressing coverity ID #499510 -- return values handled properly */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SaveOkProcessedMolfile
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: SaveOkProcessedMolfile
    // INCHI✔️❌: #define _IS_FATAL 3
    // INCHI✔️❌: #define _IS_ERROR 2
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(long) == 8
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: SaveOkProcessedMolfile

    if input_parameters.bSaveAllGoodStructsAsProblem != 0
        && return_code != _IS_FATAL as i32
        && return_code != _IS_ERROR as i32
        && problem_file.is_some()
        && !problem_file.expect("checked above").f.is_null()
        && structure_data.fPtrStart >= 0
        && structure_data.fPtrStart < structure_data.fPtrEnd
    {
        let _ = MolfileSaveCopy(
            heap,
            input_file,
            structure_data.fPtrStart,
            structure_data.fPtrEnd,
            problem_file.expect("checked above").f,
            0,
        )?;
    }
    Ok(())
}

fn component_elapsed(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    start: &inchiTime,
    clock_result: clock_t,
) -> Result<i64, SourceHeapError> {
    let mut value = heap
        .slice(clock.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let elapsed = InchiTimeElapsed(&mut value, Some(start), clock_result);
    *heap
        .slice_mut(clock)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
    Ok(elapsed)
}

fn free_structure_normalized_data(
    heap: &mut SourceHeap,
    data: SourceMutPointer<INP_ATOM_DATA2>,
    component_count: usize,
) -> Result<(), SourceHeapError> {
    if data.is_null() {
        return Ok(());
    }
    for component in 0..component_count {
        for representation in 0..TAUT_NUM as usize {
            let mut normalized = {
                let rows = heap.slice_mut(data)?;
                std::mem::take(&mut rows[component][representation])
            };
            FreeInpAtomData(heap, &mut normalized)?;
        }
    }
    inchi_free(heap, data)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CreateOneStructureINChI(
    heap: &mut SourceHeap,
    canonical_globals: &mut CANON_GLOBALS,
    clock: SourceMutPointer<INCHI_CLOCK>,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &mut INPUT_PARMS,
    mut title: Option<&mut [i8]>,
    inchi_components: &mut [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    aux_components: &mut [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    inchi_kind: i32,
    mut input_file: Option<&mut INCHI_IOSTREAM>,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    mut output_file: Option<&mut INCHI_IOSTREAM>,
    mut problem_file: Option<&mut INCHI_IOSTREAM>,
    original_input: &mut ORIG_ATOM_DATA,
    prepared_input: &mut [ORIG_ATOM_DATA],
    composite_normalized_data: &mut [[COMP_ATOM_DATA; TAUT_NUM as usize + 1]; INCHI_NUM as usize],
    input_number: i64,
    _string_buffer: Option<&mut INCHI_IOS_STRING>,
    normalization_flags: &mut NORM_CANON_FLAGS,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:809 CreateOneStructureINChI
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: CreateOneStructureINChI
    // INCHI✔️❌: int CreateOneStructureINChI( CANON_GLOBALS          *pCG,
    // INCHI✔️❌:                              INCHI_CLOCK            *ic,
    // INCHI✔️❌:                              STRUCT_DATA            *sd,
    // INCHI✔️❌:                              INPUT_PARMS            *ip,
    // INCHI✔️❌:                              char                   *szTitle,
    // INCHI✔️❌:                              PINChI2                *pINChI2[INCHI_NUM],
    // INCHI✔️❌:                              PINChI_Aux2            *pINChI_Aux2[INCHI_NUM],
    // INCHI✔️❌:                              int                    iINChI,
    // INCHI✔️❌:                              INCHI_IOSTREAM         *inp_file,
    // INCHI✔️❌:                              INCHI_IOSTREAM         *log_file,
    // INCHI✔️❌:                              INCHI_IOSTREAM         *out_file,
    // INCHI✔️❌:                              INCHI_IOSTREAM         *prb_file,
    // INCHI✔️❌:                              ORIG_ATOM_DATA         *orig_inp_data,
    // INCHI✔️❌:                              ORIG_ATOM_DATA         *prep_inp_data,
    // INCHI✔️❌:                              COMP_ATOM_DATA         composite_norm_data2[][TAUT_NUM + 1],
    // INCHI✔️❌:                              long                   num_inp,
    // INCHI✔️❌:                              INCHI_IOS_STRING		*strbuf,
    // INCHI✔️❌:                              NORM_CANON_FLAGS       *pncFlags )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, nRet = 0, n = 0l;
    // INCHI✔️❌: #if defined (TARGET_EXE_STANDALONE) && defined(_WIN32)
    // INCHI✔️❌:     int err_display;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef GHI100_FIX
    // INCHI✔️❌: #if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
    // INCHI✔️❌:     setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     PINChI2     *pINChI = NULL;
    // INCHI✔️❌:     PINChI_Aux2 *pINChI_Aux = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     INP_ATOM_DATA InpCurAtData;
    // INCHI✔️❌:     INP_ATOM_DATA *inp_cur_data;
    // INCHI✔️❌:
    // INCHI✔️❌:     INP_ATOM_DATA InpNormAtData, InpNormTautData;
    // INCHI✔️❌:     INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
    // INCHI✔️❌:     ORIG_ATOM_DATA *cur_prep_inp_data = prep_inp_data + iINChI;
    // INCHI✔️❌:     inchiTime      ulTStart;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Always create info data structures (but do not display them always )
    // INCHI✔️❌:     #ifndef COMPILE_ANSI_ONLY
    // INCHI✔️❌:     */
    // INCHI✔️❌:     int            bShowStructure = 0;
    // INCHI✔️❌:     int            bStructurePreprocessed = 0; /* All changes except disconnection */
    // INCHI✔️❌:     int            bStructureDisconnected = 0;
    // INCHI✔️❌:     int            bAlsoOutputReconnected = 0, bINCHI_LIB_Flag = 0;
    // INCHI✔️❌:     COMP_ATOM_DATA *composite_norm_data = composite_norm_data2[iINChI];
    // INCHI✔️❌:     INP_ATOM_DATA2 *all_inp_norm_data = NULL;
    // INCHI✔️❌:     /*#endif*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*        Order of actions:
    // INCHI✔️❌:
    // INCHI✔️❌:         if ( orig_inp_data is NOT empty AND
    // INCHI✔️❌:              prep_inp_data[0] IS empty ) then do
    // INCHI✔️❌:              in PreprocessOneStructure()        :
    // INCHI✔️❌:
    // INCHI✔️❌:             1. copy orig_inp_data --> prep_inp_data[0]
    // INCHI✔️❌:             2. fix odd things in prep_inp_data[0]
    // INCHI✔️❌:             3. if( orig_inp_data->bDisconnectSalts ) then
    // INCHI✔️❌:                   -- disconnect salts in prep_inp_data[0]
    // INCHI✔️❌:             4. move protons to neutralize charges on heteroatoms
    // INCHI✔️❌:             5. if( orig_inp_data->bDisconnectCoord ) then
    // INCHI✔️❌:                   -- copy prep_inp_data[0] --> prep_inp_data[1]
    // INCHI✔️❌:                   -- disconnect metals in prep_inp_data[0]
    // INCHI✔️❌:
    // INCHI✔️❌:         iINChI = 0
    // INCHI✔️❌:         =========
    // INCHI✔️❌:         (normal/disconnected layer)
    // INCHI✔️❌:
    // INCHI✔️❌:             1. normalize prep_inp_data[0] in inp_norm_data[0,1]
    // INCHI✔️❌:             2. create INChI[ iINChI ] out of inp_norm_data[0,1]
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         iINChI = 1 AND orig_inp_data->bDisconnectCoord > 0
    // INCHI✔️❌:         =================================================
    // INCHI✔️❌:         (reconnected layer)
    // INCHI✔️❌:
    // INCHI✔️❌:             1. normalize prep_inp_data[1] in inp_norm_data[0,1]
    // INCHI✔️❌:             2. create INChI[ iINChI ] out of inp_norm_data[0,1]
    // INCHI✔️❌:
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     ip->msec_LeftTime = ip->msec_MaxTime; /* start timeout countdown for each component */
    // INCHI✔️❌:
    // INCHI✔️❌:     inp_cur_data = &InpCurAtData;
    // INCHI✔️❌:     inp_norm_data[TAUT_NON] = &InpNormAtData;
    // INCHI✔️❌:     inp_norm_data[TAUT_YES] = &InpNormTautData;
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( inp_cur_data, 0, sizeof( *inp_cur_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( inp_norm_data[TAUT_NON], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( inp_norm_data[TAUT_YES], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*#ifndef COMPILE_ANSI_ONLY*/
    // INCHI✔️❌:         memset( composite_norm_data + TAUT_NON, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         memset( composite_norm_data + TAUT_YES, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         memset( composite_norm_data + TAUT_INI, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }   /*#endif*/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bAllowEmptyStructure && !orig_inp_data->at && !orig_inp_data->num_inp_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (!orig_inp_data->at || !orig_inp_data->num_inp_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* nothing to do */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (iINChI == 1 && orig_inp_data->bDisconnectCoord <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* m = iINChI; */ /* orig_inp_data index */
    // INCHI✔️❌:     if (iINChI != INCHI_BAS && iINChI != INCHI_REC)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AddErrorMessage( sd->pStrErrStruct, "Fatal undetermined program error" );
    // INCHI✔️❌:         sd->nStructReadError = 97;
    // INCHI✔️❌:         nRet = sd->nErrorType = _IS_FATAL;
    // INCHI✔️❌:         inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*******************************************************************
    // INCHI✔️❌:      *                                                                 *
    // INCHI✔️❌:      *                                                                 *
    // INCHI✔️❌:      *  Whole structure preprocessing: 1st step of the normalization   *
    // INCHI✔️❌:      *                                                                 *
    // INCHI✔️❌:      *  Happen only on the first call to CreateOneStructureINChI()      *
    // INCHI✔️❌:      *                                                                 *
    // INCHI✔️❌:      *                                                                 *
    // INCHI✔️❌:      *******************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (( !prep_inp_data->at || !prep_inp_data->num_inp_atoms ) &&
    // INCHI✔️❌:          orig_inp_data->num_inp_atoms > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* The structure has not been preprocessed */
    // INCHI✔️❌:         if (ip->msec_MaxTime)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             InchiTimeGet( &ulTStart );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         PreprocessOneStructure( ic, sd, ip, orig_inp_data, prep_inp_data );
    // INCHI✔️❌:
    // INCHI✔️❌:         pncFlags->bTautFlags[iINChI][TAUT_YES] =
    // INCHI✔️❌:                 pncFlags->bTautFlags[iINChI][TAUT_NON] =
    // INCHI✔️❌:                     sd->bTautFlags[INCHI_BAS] | ip->bTautFlags;
    // INCHI✔️❌:
    // INCHI✔️❌:         pncFlags->bTautFlagsDone[iINChI][TAUT_YES] =
    // INCHI✔️❌:             pncFlags->bTautFlagsDone[iINChI][TAUT_NON] =
    // INCHI✔️❌:             sd->bTautFlagsDone[INCHI_BAS] | ip->bTautFlagsDone;
    // INCHI✔️❌:
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*#ifndef COMPILE_ANSI_ONLY*/
    // INCHI✔️❌:             /* in this location the call happens once for each input structure, before preprocessing */
    // INCHI✔️❌:             bStructurePreprocessed = ( 0 != ( sd->bTautFlagsDone[INCHI_BAS] & (
    // INCHI✔️❌:                 TG_FLAG_MOVE_HPLUS2NEUTR_DONE |
    // INCHI✔️❌:                 TG_FLAG_DISCONNECT_SALTS_DONE |
    // INCHI✔️❌:                 TG_FLAG_MOVE_POS_CHARGES_DONE |
    // INCHI✔️❌:                 TG_FLAG_FIX_ODD_THINGS_DONE ) ) );
    // INCHI✔️❌:
    // INCHI✔️❌:             bStructureDisconnected = ( 0 != ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) );
    // INCHI✔️❌:
    // INCHI✔️❌:             bShowStructure = ( bStructurePreprocessed ||
    // INCHI✔️❌:                                bStructureDisconnected ||
    // INCHI✔️❌:                                prep_inp_data[0].num_components > 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:             /* sd->bTautFlags[] contains output flags
    // INCHI✔️❌:                ip->bTautFlags   contains input flags
    // INCHI✔️❌:             */
    // INCHI✔️❌:             bAlsoOutputReconnected = ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
    // INCHI✔️❌:                 ( ip->bTautFlags               & TG_FLAG_RECONNECT_COORD );
    // INCHI✔️❌:             bINCHI_LIB_Flag = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:             /*************** output structures to TARGET_LIB_FOR_WINCHI conditions *********************
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  Send to TARGET_LIB_FOR_WINCHI:
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  type                      component  conditions
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  COMPONENT_ORIGINAL              #0:  (num_components > 1)
    // INCHI✔️❌:              *  COMPONENT_ORIGINAL_PREPROCESSED #0:  (num_components > 1) && (preprocessed)
    // INCHI✔️❌:              *  COMPONENT_ORIGINAL              #1:  (num_components = 1) && (preprocessed)
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  Flags explanation:
    // INCHI✔️❌:              *        MAIN => iINChI=0,  RECN => iINChI=1 (Reconnected)
    // INCHI✔️❌:              *        ORIG => Original, PREP => Preprocessed
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  Possible flags:           k
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  COMP_ORIG_0_MAIN  0x0001  0  COMPONENT_ORIGINAL, bMain, component #0
    // INCHI✔️❌:              *  COMP_ORIG_0_RECN  0x0002  1  COMPONENT_ORIGINAL, bRecn, component #0
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  COMP_PREP_0_MAIN  0x0004  2  COMPONENT_ORIGINAL_PREPROCESSED, bMain, component #0
    // INCHI✔️❌:              *  COMP_PREP_0_RECN  0x0008  3  COMPONENT_ORIGINAL_PREPROCESSED, bRecn, component #0
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  COMP_ORIG_1_MAIN  0x0010  4  COMPONENT_ORIGINAL, bMain, component #1
    // INCHI✔️❌:              *  COMP_ORIG_1_RECN  0x0020  5  COMPONENT_ORIGINAL, bRecn, component #1
    // INCHI✔️❌:              *
    // INCHI✔️❌:              *  bReconnected  = k%2     (0 or 1)
    // INCHI✔️❌:              *  nComponent    = k/4     (0 or 1)
    // INCHI✔️❌:              *  bPreprocessed = (k/2)%2 (0 or 1)
    // INCHI✔️❌:              *
    // INCHI✔️❌:              ******************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:             /* Original -> Main, component #0, Original */
    // INCHI✔️❌:             if (prep_inp_data[INCHI_BAS].num_components > 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bINCHI_LIB_Flag |= COMP_ORIG_0_MAIN;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Original -> Main, component #1, Original */
    // INCHI✔️❌:                 if (prep_inp_data[INCHI_BAS].num_components == 1 && bStructurePreprocessed)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bINCHI_LIB_Flag |= COMP_ORIG_1_MAIN;
    // INCHI✔️❌:                     /* preprocessed will be added when output canonicalization results */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (bAlsoOutputReconnected)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Original -> Reconnected, component #0, Original */
    // INCHI✔️❌:                 if (prep_inp_data[INCHI_REC].num_components > 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bINCHI_LIB_Flag |= COMP_ORIG_0_RECN;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (prep_inp_data[INCHI_BAS].num_components == 1 && bStructurePreprocessed)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* Original -> Reconnected, component #1, Original */
    // INCHI✔️❌:                     bINCHI_LIB_Flag |= COMP_ORIG_1_RECN;
    // INCHI✔️❌:                     /* preprocessed will be added when output canonicalization results */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (ip->msec_MaxTime)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ip->msec_LeftTime -= InchiTimeElapsed( ic, &ulTStart );
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* display the ORIGINAL, UN-PREPROCESSED structure */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (ip->bDisplay)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (DisplayTheWholeStructure( pCG, ic, sd, ip, szTitle,
    // INCHI✔️❌:                     inp_file, log_file, orig_inp_data, num_inp,
    // INCHI✔️❌:                     -1, bShowStructure, bINCHI_LIB_Flag ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:                     goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         } /*#endif */
    // INCHI✔️❌:
    // INCHI✔️❌:         switch (sd->nErrorType)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             case _IS_ERROR:
    // INCHI✔️❌:             case _IS_FATAL:
    // INCHI✔️❌:                 /* error message */
    // INCHI✔️❌:                 nRet = TreatErrorsInReadTheStructure( sd, ip,
    // INCHI✔️❌:                                                       LOG_MASK_ALL,
    // INCHI✔️❌:                                                       inp_file, log_file, out_file, prb_file,
    // INCHI✔️❌:                                                       prep_inp_data, &num_inp );
    // INCHI✔️❌:                 goto exit_cycle;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* tranfer flags from INChI_Aux to sd */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*#ifndef COMPILE_ANSI_ONLY */ /* { */
    // INCHI✔️❌:
    // INCHI✔️❌:         /******************************************/
    // INCHI✔️❌:         /*      Displaying the structures         */
    // INCHI✔️❌:         /*          Only under WIN32              */
    // INCHI✔️❌:         /******************************************/
    // INCHI✔️❌:         if ( /* ip->bDisplayCompositeResults && !sd->bUserQuitComponentDisplay && */
    // INCHI✔️❌:              prep_inp_data[iINChI].num_components > 1
    // INCHI✔️❌:            )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             all_inp_norm_data = (INP_ATOM_DATA2 *) inchi_calloc( prep_inp_data[iINChI].num_components, sizeof( all_inp_norm_data[0] ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Display the input structure AFTER PREPROCESSING */
    // INCHI✔️❌:         switch (iINChI)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             case INCHI_BAS:
    // INCHI✔️❌:                 /*------------ Possibly disconnected structure -------------------*/
    // INCHI✔️❌:                 bStructurePreprocessed = 0 != ( sd->bTautFlagsDone[iINChI] & (
    // INCHI✔️❌:                     TG_FLAG_MOVE_HPLUS2NEUTR_DONE |
    // INCHI✔️❌:                     TG_FLAG_DISCONNECT_SALTS_DONE |
    // INCHI✔️❌:                     TG_FLAG_MOVE_POS_CHARGES_DONE |
    // INCHI✔️❌:                     TG_FLAG_MOVE_CHARGE_COORD_DONE |
    // INCHI✔️❌:                     TG_FLAG_DISCONNECT_COORD_DONE |
    // INCHI✔️❌:                     TG_FLAG_FIX_ODD_THINGS_DONE ) );
    // INCHI✔️❌:                 bINCHI_LIB_Flag = 0;
    // INCHI✔️❌:                 /* Preprocessed/Main -> Main, component #0, Preprocessed */
    // INCHI✔️❌:                 if (prep_inp_data[iINChI].num_components > 1 &&
    // INCHI✔️❌:                      bStructurePreprocessed)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bINCHI_LIB_Flag |= COMP_PREP_0_MAIN;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 bShowStructure = ( bStructurePreprocessed &&
    // INCHI✔️❌:                                    prep_inp_data[iINChI].num_components > 1 );
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             case INCHI_REC:
    // INCHI✔️❌:                 /*------------ Reconnected structure ------------------------------*/
    // INCHI✔️❌:                 bAlsoOutputReconnected =
    // INCHI✔️❌:                     ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
    // INCHI✔️❌:                     ( ip->bTautFlags               & TG_FLAG_RECONNECT_COORD );
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (!bAlsoOutputReconnected)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 bStructurePreprocessed = 0 != ( sd->bTautFlagsDone[iINChI] & (
    // INCHI✔️❌:                     TG_FLAG_MOVE_HPLUS2NEUTR_DONE |
    // INCHI✔️❌:                     TG_FLAG_DISCONNECT_SALTS_DONE |
    // INCHI✔️❌:                     TG_FLAG_MOVE_POS_CHARGES_DONE |
    // INCHI✔️❌:                     TG_FLAG_FIX_ODD_THINGS_DONE ) );
    // INCHI✔️❌:                 bINCHI_LIB_Flag = 0;
    // INCHI✔️❌:                 /* Preprocessed/Reconnected -> Reconnected, component #0, Preprocessed */
    // INCHI✔️❌:                 if (prep_inp_data[iINChI].num_components > 1 && bStructurePreprocessed)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bINCHI_LIB_Flag |= COMP_PREP_0_RECN;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 bShowStructure = ( bStructurePreprocessed &&
    // INCHI✔️❌:                                    prep_inp_data[iINChI].num_components > 1 );
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 bShowStructure = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ip->bDisplay && prep_inp_data[iINChI].num_inp_atoms > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (DisplayTheWholeStructure( pCG, ic, sd, ip, szTitle,
    // INCHI✔️❌:                 inp_file, log_file,
    // INCHI✔️❌:                 prep_inp_data + iINChI,
    // INCHI✔️❌:                 num_inp,
    // INCHI✔️❌:                 iINChI,
    // INCHI✔️❌:                 bShowStructure,
    // INCHI✔️❌:                 bINCHI_LIB_Flag ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } /* #endif */ /* } ifndef COMPILE_ANSI_ONLY */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* allocate pINChI[iINChI] and pINChI_Aux2[iINChI] -- arrays of pointers to INChI and INChI_Aux */
    // INCHI✔️❌:     /* assign values to sd->num_components[]                                                  */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: MYREALLOC2 has been replaced and the whole block rewritten to address memory leaks and reading from freed memory locations */
    // INCHI✔️❌:     do {
    // INCHI✔️❌:         if( (sd->num_components[iINChI]) <= ((long long)cur_prep_inp_data->num_components) ) {
    // INCHI✔️❌:             PINChI2* newPTR1 = (PINChI2*)inchi_calloc( ((long long)cur_prep_inp_data->num_components)+1, sizeof(PINChI2) );
    // INCHI✔️❌:             PINChI_Aux2* newPTR2 = (PINChI_Aux2*)inchi_calloc( ((long long)cur_prep_inp_data->num_components)+1, sizeof(PINChI_Aux2) );
    // INCHI✔️❌:             if ( newPTR1 && newPTR2 ) {
    // INCHI✔️❌:                 if ( (pINChI2[iINChI]) && (sd->num_components[iINChI]) > 0 )
    // INCHI✔️❌:                     memcpy(newPTR1, pINChI2[iINChI], (sd->num_components[iINChI]) * sizeof(PINChI2));
    // INCHI✔️❌:                 if ((pINChI_Aux2[iINChI]) && (sd->num_components[iINChI]) > 0)
    // INCHI✔️❌:                     memcpy(newPTR2, pINChI_Aux2[iINChI], (sd->num_components[iINChI]) * sizeof(PINChI_Aux2));
    // INCHI✔️❌:                 if (pINChI2[iINChI])
    // INCHI✔️❌:                     inchi_free(pINChI2[iINChI]);
    // INCHI✔️❌:                 if (pINChI_Aux2[iINChI])
    // INCHI✔️❌:                     inchi_free(pINChI_Aux2[iINChI]);
    // INCHI✔️❌:                 pINChI2[iINChI] = newPTR1;
    // INCHI✔️❌:                 pINChI_Aux2[iINChI] = newPTR2;
    // INCHI✔️❌:                 sd->num_components[iINChI] = cur_prep_inp_data->num_components;
    // INCHI✔️❌:                 k = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else {
    // INCHI✔️❌:                 inchi_free(newPTR1);
    // INCHI✔️❌:                 inchi_free(newPTR2);
    // INCHI✔️❌:                 k = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else { k = 0; }
    // INCHI✔️❌:     } while (0);
    // INCHI✔️❌:
    // INCHI✔️❌:     if (k)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AddErrorMessage( sd->pStrErrStruct, "Cannot allocate output data. Terminating" );
    // INCHI✔️❌:         sd->nStructReadError = 99;
    // INCHI✔️❌:         sd->nErrorType = _IS_FATAL;
    // INCHI✔️❌:         inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     pINChI = pINChI2[iINChI];
    // INCHI✔️❌:     pINChI_Aux = pINChI_Aux2[iINChI];
    // INCHI✔️❌:
    // INCHI✔️❌:     /**************************************************************************/
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*   M A I N   C Y C L E:   P R O C E S S    C O M P O N E N T S          */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*                     O N E   B Y   O N E                                */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /**************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0, nRet = 0;
    // INCHI✔️❌:             !sd->bUserQuitComponent && i < cur_prep_inp_data->num_components;
    // INCHI✔️❌:                 i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ip->msec_MaxTime)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             InchiTimeGet( &ulTStart );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #ifndef TARGET_LIB_FOR_WINCHI  /* { */
    // INCHI✔️❌: #if ( bREUSE_INCHI == 1 )
    // INCHI✔️❌:
    // INCHI✔️❌:         if ((iINChI == INCHI_REC &&
    // INCHI✔️❌:              /*( !ip->bDisplay &&
    // INCHI✔️❌:                !ip->bDisplayCompositeResults && */
    // INCHI✔️❌:                !( ip->bCompareComponents & CMP_COMPONENTS )) ||
    // INCHI✔️❌:                sd->bUserQuitComponentDisplay) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Reconnected structure (06-20-2005: added "&& !ip->bDisplayCompositeResults" to display composite structure) */
    // INCHI✔️❌:             int m = iINChI - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:             /* Find whether we have already calculated this INChI in basic (disconnected) layer */
    // INCHI✔️❌:             for (j = n = 0; j < prep_inp_data[m].num_components; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i + 1 == prep_inp_data[m].nOldCompNumber[j] &&
    // INCHI✔️❌:                     ( pINChI2[m][j][TAUT_NON] || pINChI2[m][j][TAUT_YES] ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* Yes, we have already done this */
    // INCHI✔️❌:                     if (!n++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         memcpy(pINChI + i, pINChI2[m] + j, sizeof(pINChI[0]));
    // INCHI✔️❌:                         memcpy(pINChI_Aux + i, pINChI_Aux2[m] + j, sizeof(pINChI_Aux[0]));
    // INCHI✔️❌:                         for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (pINChI[i][k])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 pINChI[i][k]->nRefCount++;
    // INCHI✔️❌:                                 if (pINChI[i][k]->nNumberOfAtoms > 0)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     switch (k)
    // INCHI✔️❌:                                     {
    // INCHI✔️❌:                                         case TAUT_NON:
    // INCHI✔️❌:                                             sd->num_non_taut[iINChI] ++;
    // INCHI✔️❌:                                             break;
    // INCHI✔️❌:                                         case TAUT_YES:
    // INCHI✔️❌:                                             if (pINChI[i][k]->lenTautomer > 0)
    // INCHI✔️❌:                                             {
    // INCHI✔️❌:                                                 sd->num_taut[iINChI] ++;
    // INCHI✔️❌:                                             }
    // INCHI✔️❌:                                             else
    // INCHI✔️❌:                                                 if (!pINChI[i][TAUT_NON] ||
    // INCHI✔️❌:                                                      !pINChI[i][TAUT_NON]->nNumberOfAtoms)
    // INCHI✔️❌:                                                 {
    // INCHI✔️❌:                                                     sd->num_non_taut[iINChI] ++;
    // INCHI✔️❌:                                                 }
    // INCHI✔️❌:                                             break;
    // INCHI✔️❌:                                     }
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (pINChI_Aux[i][k])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 pINChI_Aux[i][k]->nRefCount++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (n == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (n > 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* ith component is equivalent to more than one another component */
    // INCHI✔️❌:                 AddErrorMessage( sd->pStrErrStruct, "Cannot distinguish components" );
    // INCHI✔️❌:                 sd->nStructReadError = 99;
    // INCHI✔️❌:                 sd->nErrorType = _IS_ERROR;
    // INCHI✔️❌:                 inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #endif /* } TARGET_LIB_FOR_WINCHI */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*****************************************************/
    // INCHI✔️❌:         /*  a) Allocate memory and extract current component */
    // INCHI✔️❌:         /*****************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:         nRet = GetOneComponent( ic, sd, ip,
    // INCHI✔️❌:                                 log_file, out_file,
    // INCHI✔️❌:                                 inp_cur_data, cur_prep_inp_data,
    // INCHI✔️❌:                                 i, num_inp );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ip->msec_MaxTime)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ip->msec_LeftTime -= InchiTimeElapsed( ic, &ulTStart );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         switch (nRet)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             case _IS_ERROR:
    // INCHI✔️❌:             case _IS_FATAL:
    // INCHI✔️❌:                 goto exit_cycle;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if !defined(TARGET_API_LIB) && !defined(COMPILE_ANSI_ONLY)
    // INCHI✔️❌:         /*  console request: Display the component? */
    // INCHI✔️❌:         if (ip->bDisplay && inp_file->f != stdin)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (user_quit( ic, "Enter=Display Component, Esc=Stop ?", ip->ulDisplTime ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sd->bUserQuitComponent = 1;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /*#ifndef COMPILE_ANSI_ONLY
    // INCHI✔️❌:         { */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  b) Display the extracted original component structure */
    // INCHI✔️❌:         if (ip->bDisplay && inp_cur_data->at && !sd->bUserQuitComponentDisplay)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_prep_inp_data->num_components == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sprintf(szTitle, "%sInput Structure #%ld.%s%s%s%s%s",
    // INCHI✔️❌:                     bStructurePreprocessed ? "Preprocessed " : "",
    // INCHI✔️❌:                     num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sprintf(szTitle, "Component #%d of %d, Input Structure #%ld.%s%s%s%s%s",
    // INCHI✔️❌:                     i + 1, cur_prep_inp_data->num_components,
    // INCHI✔️❌:                     num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if defined (TARGET_EXE_STANDALONE) && defined(_WIN32)
    // INCHI✔️❌:             err_display = DisplayStructure( pCG,
    // INCHI✔️❌:                                     inp_cur_data->at,
    // INCHI✔️❌:                                     inp_cur_data->num_at,
    // INCHI✔️❌:                                     NULL, /* OAD_Polymer *polymer, */
    // INCHI✔️❌:                                     0,
    // INCHI✔️❌:                                     1,
    // INCHI✔️❌:                                     0,
    // INCHI✔️❌:                                     NULL,
    // INCHI✔️❌:                                     1           /*isotopic*/,
    // INCHI✔️❌:                                     0           /*taut*/,
    // INCHI✔️❌:                                     NULL,
    // INCHI✔️❌:                                     NULL,
    // INCHI✔️❌:                                     ip->bAbcNumbers,
    // INCHI✔️❌:                                     &ip->dp,
    // INCHI✔️❌:                                     ip->nMode,
    // INCHI✔️❌:                                     szTitle );
    // INCHI✔️❌:
    // INCHI✔️❌:             sd->bUserQuitComponentDisplay = (err_display == ESC_KEY );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!err_display)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_fprintf( stderr, "Cannot display the structure\n" );
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:             if (DRAWDATA && DRAWDATA_EXISTS)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 struct DrawData vDrawData;
    // INCHI✔️❌:                 int    nType = COMPONENT_ORIGINAL;
    // INCHI✔️❌:                 vDrawData.pWindowData = CreateWinData_( pCG,
    // INCHI✔️❌:                                                         inp_cur_data->at,
    // INCHI✔️❌:                                                         inp_cur_data->num_at,
    // INCHI✔️❌:                                                         0,
    // INCHI✔️❌:                                                         1 /* bAdd_DT_to_num_H */,
    // INCHI✔️❌:                                                         0,
    // INCHI✔️❌:                                                         NULL,
    // INCHI✔️❌:                                                         1 /* display isotopic if present */,
    // INCHI✔️❌:                                                         0,
    // INCHI✔️❌:                                                         NULL,
    // INCHI✔️❌:                                                         NULL,
    // INCHI✔️❌:                                                         ip->bAbcNumbers,
    // INCHI✔️❌:                                                         &ip->dp,
    // INCHI✔️❌:                                                         ip->nMode );
    // INCHI✔️❌:                 if (vDrawData.pWindowData != NULL)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (DRAWDATA_EXISTS( i + 1, nType, iINChI ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* i = component number */
    // INCHI✔️❌:                         nType = COMPONENT_ORIGINAL_PREPROCESSED;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     vDrawData.nComponent = i + 1;
    // INCHI✔️❌:                     vDrawData.nType = nType;
    // INCHI✔️❌:                     vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
    // INCHI✔️❌:                     vDrawData.szTitle = inchi__strdup( szTitle );
    // INCHI✔️❌:                     vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
    // INCHI✔️❌:                     DRAWDATA( &vDrawData );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*#endif */  /* } COMPILE_ANSI_ONLY */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /*******************************************************************************/
    // INCHI✔️❌:         /*                                                                             */
    // INCHI✔️❌:         /*  N O R M A L I Z A T I O N    a n d     C A N O N I C A L I Z A T I O N     */
    // INCHI✔️❌:         /*                                                                             */
    // INCHI✔️❌:         /*         (both tautomeric and non-tautomeric if requested)                   */
    // INCHI✔️❌:         /*                                                                             */
    // INCHI✔️❌:         /*******************************************************************************/
    // INCHI✔️❌:         /*  c) Create the component's INChI ( copies ip->bTautFlags into sd->bTautFlags)*/
    // INCHI✔️❌:         /*******************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:         nRet = CreateOneComponentINChI( pCG, ic, sd, ip,
    // INCHI✔️❌:                                         inp_cur_data, orig_inp_data,
    // INCHI✔️❌:                                         pINChI/*2[iINChI]*/,
    // INCHI✔️❌:                                         pINChI_Aux/*2[iINChI]*/,
    // INCHI✔️❌:                                         iINChI, i, num_inp,
    // INCHI✔️❌:                                         inp_norm_data, pncFlags, log_file );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  d) Display one component structure and/or INChI results only if there was no error */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* #ifndef COMPILE_ANSI_ONLY */ /* { */
    // INCHI✔️❌:         if (!nRet)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  output one component INChI to the stdout if requested */
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             if ( ip->bDisplayEachComponentINChI ) {
    // INCHI✔️❌:                 int cur_num_non_taut = (pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms>0);
    // INCHI✔️❌:                 int cur_num_taut     = (pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms>0);
    // INCHI✔️❌:                 if ( ip->bDisplayEachComponentINChI && cur_num_non_taut + cur_num_taut ) {
    // INCHI✔️❌:                     SortAndPrintINChI( pCG, stdout, pStr, nStrLen, NULL,
    // INCHI✔️❌:                                        ip, 1, cur_num_non_taut, cur_num_taut,
    // INCHI✔️❌:                                        num_inp, pINChI+i, pINChI_Aux+i,
    // INCHI✔️❌:                                        save_opt_bits);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             */
    // INCHI✔️❌:             /**************************************************************************
    // INCHI✔️❌:              * display from one up to 4 structure pictures-results for each component *
    // INCHI✔️❌:              * Enable buttons:                                                        *
    // INCHI✔️❌:              * BN (non-tautomeric non-isotopic): inp_norm_data[0]->bExists            *
    // INCHI✔️❌:              * TN (tautomeric non-isotopic):     inp_norm_data[1]->bExists            *
    // INCHI✔️❌:              * BI (non-tautomeric isotopic):     inp_norm_data[0]->bExists &&         *
    // INCHI✔️❌:              *                                   inp_norm_data[0]->bHasIsotopicLayer  *
    // INCHI✔️❌:              * TI (tautomeric isotopic):         inp_norm_data[1]->bExists &&         *
    // INCHI✔️❌:              *                                   inp_norm_data[1]->bHasIsotopicLayer  *
    // INCHI✔️❌:              **************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:             int bIsotopic, bTautomeric, bDisplayTaut, bHasIsotopicLayer, bFixedBondsTaut, m_max, m, nNumDisplayedFixedBondTaut = 0; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:
    // INCHI✔️❌:             for ( j = 0;
    // INCHI✔️❌:                   ip->bDisplay && !sd->bUserQuitComponentDisplay && j < TAUT_NUM;
    // INCHI✔️❌:                   j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (inp_norm_data[j]->bExists && !inp_norm_data[j]->bDeleted)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bTautomeric = ( pINChI[i][j]->lenTautomer > 0 );
    // INCHI✔️❌:                      /* same as (inp_norm_data[j]->bTautomeric > 0) */
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* If requested tautomeric and no tautmerism found then do not say mobile or fixed H. 2004-10-27 */
    // INCHI✔️❌:                     bDisplayTaut = ( !( ip->nMode & REQ_MODE_BASIC ) && !bTautomeric ) ? -1 : bTautomeric;
    // INCHI✔️❌:                     bHasIsotopicLayer = ( inp_norm_data[j]->bHasIsotopicLayer > 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:                     for (k = 0; k <= bHasIsotopicLayer; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bIsotopic = ( k > 0 );
    // INCHI✔️❌:                         m_max = inp_norm_data[j]->at_fixed_bonds && inp_norm_data[j]->bTautPreprocessed ? 1 : 0;
    // INCHI✔️❌:                         for (m = m_max; 0 <= m; m--)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bFixedBondsTaut = ( m > 0 );
    // INCHI✔️❌:                             nNumDisplayedFixedBondTaut += bFixedBondsTaut;
    // INCHI✔️❌:                                 /* display only one time */
    // INCHI✔️❌:
    // INCHI✔️❌:                             /*  Added number of components, added another format for a single component case - DCh */
    // INCHI✔️❌:                             if (cur_prep_inp_data->num_components > 1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 sprintf(szTitle, "%s Component #%d of %d, Structure #%ld%s%s.%s%s%s%s%s",
    // INCHI✔️❌:                                     bFixedBondsTaut ? "Preprocessed" : "Result for",
    // INCHI✔️❌:                                     i + 1, cur_prep_inp_data->num_components, num_inp,
    // INCHI✔️❌:                                     bDisplayTaut == 1 ? ", mobile H" : bDisplayTaut == 0 ? ", fixed H" : "",
    // INCHI✔️❌:                                     bIsotopic ? ", isotopic" : "",
    // INCHI✔️❌:                                     SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 sprintf(szTitle, "%s Structure #%ld%s%s.%s%s%s%s%s",
    // INCHI✔️❌:                                     bFixedBondsTaut ? "Preprocessed" : "Result for",
    // INCHI✔️❌:                                     num_inp,
    // INCHI✔️❌:                                     bDisplayTaut == 1 ? ", mobile H" : bDisplayTaut == 0 ? ", fixed H" : "",
    // INCHI✔️❌:                                     bIsotopic ? ", isotopic" : "",
    // INCHI✔️❌:                                     SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if defined (TARGET_EXE_STANDALONE) && defined(_WIN32)
    // INCHI✔️❌:                             if (bFixedBondsTaut && nNumDisplayedFixedBondTaut != 1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             if (ip->bDisplay)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (bFixedBondsTaut)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     err_display = DisplayStructure( pCG,
    // INCHI✔️❌:                                                             inp_norm_data[j]->at_fixed_bonds,
    // INCHI✔️❌:                                                             inp_norm_data[j]->num_at,
    // INCHI✔️❌:                                                             NULL, /* OAD_Polymer *polymer, */
    // INCHI✔️❌:                                                             inp_norm_data[j]->num_removed_H,
    // INCHI✔️❌:                                                             0 /*bAdd_DT_to_num_H*/,
    // INCHI✔️❌:                                                             inp_norm_data[j]->nNumRemovedProtons,
    // INCHI✔️❌:                                                             inp_norm_data[j]->nNumRemovedProtonsIsotopic,
    // INCHI✔️❌:                                                             bHasIsotopicLayer,
    // INCHI✔️❌:                                                             j,
    // INCHI✔️❌:                                                             NULL,
    // INCHI✔️❌:                                                             NULL,
    // INCHI✔️❌:                                                             ip->bAbcNumbers,
    // INCHI✔️❌:                                                             &ip->dp,
    // INCHI✔️❌:                                                             ip->nMode,
    // INCHI✔️❌:                                                             szTitle );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     err_display = DisplayStructure( pCG,
    // INCHI✔️❌:                                                             inp_norm_data[j]->at,
    // INCHI✔️❌:                                                             inp_norm_data[j]->num_at,
    // INCHI✔️❌:                                                             NULL, /* OAD_Polymer *polymer, */
    // INCHI✔️❌:                                                             0,
    // INCHI✔️❌:                                                             0 /*bAdd_DT_to_num_H*/,
    // INCHI✔️❌:                                                             0,
    // INCHI✔️❌:                                                             NULL,
    // INCHI✔️❌:                                                             k,
    // INCHI✔️❌:                                                             j,
    // INCHI✔️❌:                                                             pINChI[i],
    // INCHI✔️❌:                                                             pINChI_Aux[i],
    // INCHI✔️❌:                                                             ip->bAbcNumbers,
    // INCHI✔️❌:                                                             &ip->dp,
    // INCHI✔️❌:                                                             ip->nMode, szTitle );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if ((sd->bUserQuitComponentDisplay = (err_display == ESC_KEY ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     break;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (DRAWDATA && !bFixedBondsTaut)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 struct DrawData vDrawData;
    // INCHI✔️❌:                                 vDrawData.pWindowData =
    // INCHI✔️❌:                                     CreateWinData_( pCG,
    // INCHI✔️❌:                                                     inp_norm_data[j]->at,
    // INCHI✔️❌:                                                     inp_norm_data[j]->num_at,
    // INCHI✔️❌:                                                     0,
    // INCHI✔️❌:                                                     0 /* bAdd_DT_to_num_H */,
    // INCHI✔️❌:                                                     0,
    // INCHI✔️❌:                                                     NULL,
    // INCHI✔️❌:                                                     k,
    // INCHI✔️❌:                                                     j,
    // INCHI✔️❌:                                                     pINChI[i],
    // INCHI✔️❌:                                                     pINChI_Aux[i],
    // INCHI✔️❌:                                                     ip->bAbcNumbers,
    // INCHI✔️❌:                                                     &ip->dp,
    // INCHI✔️❌:                                                     ip->nMode );
    // INCHI✔️❌:
    // INCHI✔️❌:                                 if (vDrawData.pWindowData != NULL)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     int nType;
    // INCHI✔️❌:                                     vDrawData.nComponent = i + 1;
    // INCHI✔️❌:                                     if (bTautomeric == 0)
    // INCHI✔️❌:                                         nType = ( bIsotopic == 0 ) ? COMPONENT_BN
    // INCHI✔️❌:                                         : COMPONENT_BI;
    // INCHI✔️❌:                                     else
    // INCHI✔️❌:                                         nType = ( bIsotopic == 0 ) ? COMPONENT_TN
    // INCHI✔️❌:                                         : COMPONENT_TI;
    // INCHI✔️❌:                                     vDrawData.nType = nType;
    // INCHI✔️❌:
    // INCHI✔️❌:                                     vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
    // INCHI✔️❌:                                     vDrawData.szTitle = inchi__strdup( szTitle );
    // INCHI✔️❌:                                     vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
    // INCHI✔️❌:                                     DRAWDATA( &vDrawData );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else if (DRAWDATA && bFixedBondsTaut)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 struct DrawData vDrawData;
    // INCHI✔️❌:                                 if (( ip->bCompareComponents & CMP_COMPONENTS ) &&
    // INCHI✔️❌:                                      !( ip->bCompareComponents & CMP_COMPONENTS_NONTAUT ) &&
    // INCHI✔️❌:                                      !bIsotopic == !inp_norm_data[j]->bHasIsotopicLayer)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:
    // INCHI✔️❌:                                     vDrawData.pWindowData =
    // INCHI✔️❌:                                         CreateWinData_( pCG,
    // INCHI✔️❌:                                                         inp_norm_data[j]->at_fixed_bonds,
    // INCHI✔️❌:                                                         inp_norm_data[j]->num_at,
    // INCHI✔️❌:                                                         inp_norm_data[j]->num_removed_H,
    // INCHI✔️❌:                                                         0 /* bAdd_DT_to_num_H */,
    // INCHI✔️❌:                                                         inp_norm_data[j]->nNumRemovedProtons,
    // INCHI✔️❌:                                                         inp_norm_data[j]->nNumRemovedProtonsIsotopic,
    // INCHI✔️❌:                                                         k,
    // INCHI✔️❌:                                                         j,
    // INCHI✔️❌:                                                         NULL,
    // INCHI✔️❌:                                                         NULL,
    // INCHI✔️❌:                                                         ip->bAbcNumbers,
    // INCHI✔️❌:                                                         &ip->dp,
    // INCHI✔️❌:                                                         ip->nMode );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     continue;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 if (vDrawData.pWindowData != NULL)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     vDrawData.nComponent = i + 1;
    // INCHI✔️❌:                                     vDrawData.nType = COMPONENT_ORIGINAL_PREPROCESSED;
    // INCHI✔️❌:                                     vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
    // INCHI✔️❌:                                     vDrawData.szTitle = inchi__strdup( szTitle );
    // INCHI✔️❌:                                     vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
    // INCHI✔️❌:                                     DRAWDATA( &vDrawData );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* Save normalized components for composite display */
    // INCHI✔️❌:             if ( /*ip->bDisplayCompositeResults && */
    // INCHI✔️❌:                  all_inp_norm_data
    // INCHI✔️❌:                 )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (inp_norm_data[j]->bExists)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         all_inp_norm_data[i][j] = *inp_norm_data[j];
    // INCHI✔️❌:                         memset( inp_norm_data[j], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* #endif */ /* } COMPILE_ANSI_ONLY */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nRet)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nRet = TreatErrorsInCreateOneComponentINChI( sd, ip,
    // INCHI✔️❌:                                                          cur_prep_inp_data,
    // INCHI✔️❌:                                                          i, num_inp, inp_file,
    // INCHI✔️❌:                                                          log_file, out_file, prb_file );
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /**************************************************************************/
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*   E N D   O F   T H E    M A I N   C Y C L E   P R O C E S S I N G     */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*          C O M P O N E N T S    O N E   B Y   O N E                    */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /*                                                                        */
    // INCHI✔️❌:     /**************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌: exit_cycle:
    // INCHI✔️❌:     switch (nRet)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         case _IS_FATAL:
    // INCHI✔️❌:         case _IS_ERROR:
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         default:
    // INCHI✔️❌:
    // INCHI✔️❌:         /* #ifndef COMPILE_ANSI_ONLY *//* { */
    // INCHI✔️❌:             /* composite results picture(s) */
    // INCHI✔️❌:             if (all_inp_norm_data)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 CreateCompositeNormAtom( composite_norm_data,
    // INCHI✔️❌:                                          all_inp_norm_data,
    // INCHI✔️❌:                                          prep_inp_data[iINChI].num_components );
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                 for ( i = 0; i < prep_inp_data[iINChI].num_components; i ++ ) {
    // INCHI✔️❌:                     for ( k = 0; k < TAUT_NUM; k ++ ) {
    // INCHI✔️❌:                        FreeInpAtomData( &all_inp_norm_data[i][k] );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_free( all_inp_norm_data );
    // INCHI✔️❌:                 all_inp_norm_data = NULL;
    // INCHI✔️❌:                 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* #endif */ /* } COMPILE_ANSI_ONLY */
    // INCHI✔️❌:
    // INCHI✔️❌:             break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*#ifndef COMPILE_ANSI_ONLY*/ /* { */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* avoid memory leaks in case of error */
    // INCHI✔️❌:     if (all_inp_norm_data)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < prep_inp_data[iINChI].num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 FreeInpAtomData( &all_inp_norm_data[i][k] );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_free( all_inp_norm_data );
    // INCHI✔️❌:         all_inp_norm_data = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌: /*#endif */ /* } COMPILE_ANSI_ONLY */
    // INCHI✔️❌:
    // INCHI✔️❌:     FreeInpAtomData( inp_cur_data );
    // INCHI✔️❌:     for (i = 0; i < TAUT_NUM; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         FreeInpAtomData( inp_norm_data[i] );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: CreateOneStructureINChI
    // END INCHI C FUNCTION: CreateOneStructureINChI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateOneStructureINChI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; bREUSE_INCHI=1.
    // INCHI✔️❌: TARGET_LIB_FOR_WINCHI and standalone Windows display blocks are inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateOneStructureINChI
    input_parameters.msec_LeftTime = input_parameters.msec_MaxTime;

    let mut current_input = INP_ATOM_DATA::default();
    let mut normalized_data: INP_ATOM_DATA2 = std::array::from_fn(|_| INP_ATOM_DATA::default());
    let mut structure_preprocessed = false;

    if original_input.at.is_null() && original_input.num_inp_atoms == 0 {
        if input_parameters.bAllowEmptyStructure == 0 {
            return Ok(0);
        }
    } else if original_input.at.is_null() || original_input.num_inp_atoms == 0 {
        return Ok(0);
    }
    if inchi_kind == INCHI_REC as i32 && original_input.bDisconnectCoord <= 0 {
        return Ok(0);
    }
    let kind = match inchi_kind {
        value if value == INCHI_BAS as i32 => INCHI_BAS as usize,
        value if value == INCHI_REC as i32 => INCHI_REC as usize,
        _ => {
            let message = b"Fatal undetermined program error\0".map(|byte| byte as i8);
            AddErrorMessage(Some(&mut structure_data.pStrErrStruct), Some(&message))?;
            structure_data.nStructReadError = 97;
            structure_data.nErrorType = _IS_FATAL as i32;
            return Ok(_IS_FATAL as i32);
        }
    };
    prepared_input
        .get(kind)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    composite_normalized_data[kind] = std::array::from_fn(|_| COMP_ATOM_DATA::default());

    let mut result = 0_i32;
    if (prepared_input[0].at.is_null() || prepared_input[0].num_inp_atoms == 0)
        && original_input.num_inp_atoms > 0
    {
        let mut start = inchiTime::default();
        if input_parameters.msec_MaxTime != 0 {
            InchiTimeGet(&mut start, clock_result);
        }
        let mut clock_value = heap
            .slice(clock.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let _ = PreprocessOneStructure(
            heap,
            Some(&mut clock_value),
            structure_data,
            input_parameters,
            original_input,
            prepared_input,
        )?;
        heap.slice_mut(clock)?[0] = clock_value;
        structure_preprocessed = structure_data.bTautFlagsDone[INCHI_BAS as usize]
            & u64::from(
                TG_FLAG_MOVE_HPLUS2NEUTR_DONE
                    | TG_FLAG_DISCONNECT_SALTS_DONE
                    | TG_FLAG_MOVE_POS_CHARGES_DONE
                    | TG_FLAG_FIX_ODD_THINGS_DONE,
            )
            != 0;
        let base_flags =
            structure_data.bTautFlags[INCHI_BAS as usize] | input_parameters.bTautFlags;
        normalization_flags.bTautFlags[kind] = [base_flags; TAUT_NUM as usize];
        let base_flags_done =
            structure_data.bTautFlagsDone[INCHI_BAS as usize] | input_parameters.bTautFlagsDone;
        normalization_flags.bTautFlagsDone[kind] = [base_flags_done; TAUT_NUM as usize];
        if input_parameters.msec_MaxTime != 0 {
            let elapsed = component_elapsed(heap, clock, &start, clock_result)?;
            input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
        }
        if matches!(structure_data.nErrorType, value if value == _IS_ERROR as i32 || value == _IS_FATAL as i32)
        {
            let mut number = input_number;
            result = TreatErrorsInReadTheStructure(
                heap,
                structure_data,
                input_parameters,
                LOG_MASK_ALL as i32,
                input_file.as_deref_mut(),
                log_file.as_deref_mut(),
                output_file.as_deref_mut(),
                problem_file.as_deref_mut(),
                &prepared_input[0],
                &mut number,
            )?;
            FreeInpAtomData(heap, &mut current_input)?;
            for normalized in &mut normalized_data {
                FreeInpAtomData(heap, normalized)?;
            }
            return Ok(result);
        }
    }

    let component_count_i32 = prepared_input[kind].num_components;
    let component_count = usize::try_from(component_count_i32.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut all_normalized = SourceMutPointer::<INP_ATOM_DATA2>::null();
    if component_count_i32 > 1 {
        all_normalized = match inchi_calloc::<INP_ATOM_DATA2>(
            heap,
            component_count as u64,
            INP_ATOM_DATA2_GCC_LP64_SIZE,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }

    if structure_data.num_components[kind] <= component_count_i32 {
        let row_count = i64::from(component_count_i32)
            .checked_add(1)
            .and_then(|value| u64::try_from(value).ok())
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let new_inchi = match inchi_calloc::<PINChI2>(heap, row_count, PINCHI2_GCC_LP64_SIZE) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        let new_aux = match inchi_calloc::<PINChI_Aux2>(heap, row_count, PINCHI2_GCC_LP64_SIZE) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if new_inchi.is_null() || new_aux.is_null() {
            inchi_free(heap, new_inchi)?;
            inchi_free(heap, new_aux)?;
            let message = b"Cannot allocate output data. Terminating\0".map(|byte| byte as i8);
            AddErrorMessage(Some(&mut structure_data.pStrErrStruct), Some(&message))?;
            structure_data.nStructReadError = 99;
            structure_data.nErrorType = _IS_FATAL as i32;
            free_structure_normalized_data(heap, all_normalized, component_count)?;
            return Ok(result);
        }
        let old_count = usize::try_from(structure_data.num_components[kind].max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if !inchi_components[kind].is_null() && old_count > 0 {
            let old = heap.slice(inchi_components[kind].as_const())?[..old_count].to_vec();
            heap.slice_mut(new_inchi)?[..old_count].copy_from_slice(&old);
        }
        if !aux_components[kind].is_null() && old_count > 0 {
            let old = heap.slice(aux_components[kind].as_const())?[..old_count].to_vec();
            heap.slice_mut(new_aux)?[..old_count].copy_from_slice(&old);
        }
        inchi_free(heap, inchi_components[kind])?;
        inchi_free(heap, aux_components[kind])?;
        inchi_components[kind] = new_inchi;
        aux_components[kind] = new_aux;
        structure_data.num_components[kind] = component_count_i32;
    }

    match kind {
        value if value == INCHI_BAS as usize => {
            structure_preprocessed = structure_data.bTautFlagsDone[kind]
                & u64::from(
                    TG_FLAG_MOVE_HPLUS2NEUTR_DONE
                        | TG_FLAG_DISCONNECT_SALTS_DONE
                        | TG_FLAG_MOVE_POS_CHARGES_DONE
                        | TG_FLAG_MOVE_CHARGE_COORD_DONE
                        | TG_FLAG_DISCONNECT_COORD_DONE
                        | TG_FLAG_FIX_ODD_THINGS_DONE,
                )
                != 0;
        }
        value if value == INCHI_REC as usize => {
            let also_output_reconnected = structure_data.bTautFlagsDone[INCHI_BAS as usize]
                & u64::from(TG_FLAG_DISCONNECT_COORD_DONE)
                != 0
                && input_parameters.bTautFlags & u64::from(TG_FLAG_RECONNECT_COORD) != 0;
            if also_output_reconnected {
                structure_preprocessed = structure_data.bTautFlagsDone[kind]
                    & u64::from(
                        TG_FLAG_MOVE_HPLUS2NEUTR_DONE
                            | TG_FLAG_DISCONNECT_SALTS_DONE
                            | TG_FLAG_MOVE_POS_CHARGES_DONE
                            | TG_FLAG_FIX_ODD_THINGS_DONE,
                    )
                    != 0;
            }
        }
        _ => unreachable!(),
    }

    result = 0;
    for component in 0..component_count {
        if structure_data.bUserQuitComponent != 0 {
            break;
        }

        if (kind == INCHI_REC as usize
            && (input_parameters.bCompareComponents & CMP_COMPONENTS as i32) == 0)
            || structure_data.bUserQuitComponentDisplay != 0
        {
            let mut matches = 0_i32;
            let base_count =
                usize::try_from(prepared_input[INCHI_BAS as usize].num_components.max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            if !prepared_input[INCHI_BAS as usize].nOldCompNumber.is_null()
                && !inchi_components[INCHI_BAS as usize].is_null()
            {
                for base_component in 0..base_count {
                    let old_number = heap
                        .slice(prepared_input[INCHI_BAS as usize].nOldCompNumber.as_const())?
                        [base_component];
                    let base_row = heap.slice(inchi_components[INCHI_BAS as usize].as_const())?
                        [base_component];
                    if usize::from(old_number) == component + 1
                        && base_row.iter().any(|pointer| !pointer.is_null())
                    {
                        matches = matches.wrapping_add(1);
                        if matches == 1 {
                            heap.slice_mut(inchi_components[kind])?[component] = base_row;
                            let base_aux = heap
                                .slice(aux_components[INCHI_BAS as usize].as_const())?
                                [base_component];
                            heap.slice_mut(aux_components[kind])?[component] = base_aux;
                            for representation in 0..TAUT_NUM as usize {
                                let inchi_pointer = base_row[representation];
                                if !inchi_pointer.is_null() {
                                    let inchi = &mut heap.slice_mut(inchi_pointer)?[0];
                                    inchi.nRefCount = inchi.nRefCount.wrapping_add(1);
                                    if inchi.nNumberOfAtoms > 0 {
                                        if representation == TAUT_NON as usize {
                                            structure_data.num_non_taut[kind] =
                                                structure_data.num_non_taut[kind].wrapping_add(1);
                                        } else if inchi.lenTautomer > 0 {
                                            structure_data.num_taut[kind] =
                                                structure_data.num_taut[kind].wrapping_add(1);
                                        } else if base_row[TAUT_NON as usize].is_null()
                                            || heap.slice(base_row[TAUT_NON as usize].as_const())?
                                                [0]
                                            .nNumberOfAtoms
                                                == 0
                                        {
                                            structure_data.num_non_taut[kind] =
                                                structure_data.num_non_taut[kind].wrapping_add(1);
                                        }
                                    }
                                }
                                let aux_pointer = base_aux[representation];
                                if !aux_pointer.is_null() {
                                    let aux = &mut heap.slice_mut(aux_pointer)?[0];
                                    aux.nRefCount = aux.nRefCount.wrapping_add(1);
                                }
                            }
                        }
                    }
                }
            }
            if matches == 1 {
                continue;
            }
            if matches > 1 {
                let message = b"Cannot distinguish components\0".map(|byte| byte as i8);
                AddErrorMessage(Some(&mut structure_data.pStrErrStruct), Some(&message))?;
                structure_data.nStructReadError = 99;
                structure_data.nErrorType = _IS_ERROR as i32;
                inchi_free(heap, all_normalized)?;
                return Ok(result);
            }
        }

        let mut start = inchiTime::default();
        if input_parameters.msec_MaxTime != 0 {
            InchiTimeGet(&mut start, clock_result);
        }
        let mut clock_value = heap
            .slice(clock.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        result = GetOneComponent(
            heap,
            &mut clock_value,
            structure_data,
            input_parameters,
            log_file.as_deref_mut(),
            output_file.as_deref_mut(),
            &mut current_input,
            &prepared_input[kind],
            component as i32,
            input_number,
            clock_result,
            clock_result,
        )?;
        heap.slice_mut(clock)?[0] = clock_value;
        if input_parameters.msec_MaxTime != 0 {
            let elapsed = component_elapsed(heap, clock, &start, clock_result)?;
            input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
        }
        if result == _IS_ERROR as i32 || result == _IS_FATAL as i32 {
            break;
        }

        if input_parameters.bDisplay != 0
            && !current_input.at.is_null()
            && structure_data.bUserQuitComponentDisplay == 0
        {
            let suffix = sdf_label_value(heap, input_parameters)?;
            let mut bytes = if component_count_i32 == 1 {
                let prefix = if structure_preprocessed {
                    "Preprocessed "
                } else {
                    ""
                };
                format!("{prefix}Input Structure #{input_number}.").into_bytes()
            } else {
                format!(
                    "Component #{} of {}, Input Structure #{input_number}.",
                    component + 1,
                    component_count_i32
                )
                .into_bytes()
            };
            bytes.extend_from_slice(&suffix);
            if kind == INCHI_REC as usize {
                bytes.extend_from_slice(b" (Reconnected)");
            }
            overwrite_source_c_buffer(title.as_deref_mut(), &bytes)?;
        }

        let mut fallback_log = INCHI_IOSTREAM::default();
        let component_log = log_file.as_deref_mut().unwrap_or(&mut fallback_log);
        let mut inchi_rows = heap.slice(inchi_components[kind].as_const())?.to_vec();
        let mut aux_rows = heap.slice(aux_components[kind].as_const())?.to_vec();
        result = CreateOneComponentINChI(
            heap,
            canonical_globals,
            clock,
            structure_data,
            input_parameters,
            &mut current_input,
            Some(original_input),
            &mut inchi_rows,
            &mut aux_rows,
            inchi_kind,
            component as i32,
            input_number,
            &mut normalized_data,
            normalization_flags,
            component_log,
            clock_result,
        )?;
        heap.slice_mut(inchi_components[kind])?[..inchi_rows.len()].copy_from_slice(&inchi_rows);
        heap.slice_mut(aux_components[kind])?[..aux_rows.len()].copy_from_slice(&aux_rows);

        if result == 0 && input_parameters.bDisplay != 0 {
            let suffix = sdf_label_value(heap, input_parameters)?;
            for representation in 0..TAUT_NUM as usize {
                if normalized_data[representation].bExists == 0
                    || normalized_data[representation].bDeleted != 0
                {
                    continue;
                }
                let inchi_pointer = inchi_rows[component][representation];
                let tautomeric =
                    i32::from(heap.slice(inchi_pointer.as_const())?[0].lenTautomer > 0);
                let display_taut =
                    if input_parameters.nMode & u64::from(REQ_MODE_BASIC) == 0 && tautomeric == 0 {
                        -1
                    } else {
                        tautomeric
                    };
                let has_isotopic_layer =
                    i32::from(normalized_data[representation].bHasIsotopicLayer > 0);
                for isotopic in 0..=has_isotopic_layer {
                    let fixed_max = i32::from(
                        !normalized_data[representation].at_fixed_bonds.is_null()
                            && normalized_data[representation].bTautPreprocessed != 0,
                    );
                    for fixed_bonds in (0..=fixed_max).rev() {
                        let prefix = if fixed_bonds != 0 {
                            "Preprocessed"
                        } else {
                            "Result for"
                        };
                        let mut bytes = if component_count_i32 > 1 {
                            format!(
                                "{prefix} Component #{} of {}, Structure #{input_number}",
                                component + 1,
                                component_count_i32
                            )
                            .into_bytes()
                        } else {
                            format!("{prefix} Structure #{input_number}").into_bytes()
                        };
                        if display_taut == 1 {
                            bytes.extend_from_slice(b", mobile H");
                        } else if display_taut == 0 {
                            bytes.extend_from_slice(b", fixed H");
                        }
                        if isotopic != 0 {
                            bytes.extend_from_slice(b", isotopic");
                        }
                        bytes.push(b'.');
                        bytes.extend_from_slice(&suffix);
                        if kind == INCHI_REC as usize {
                            bytes.extend_from_slice(b" (Reconnected)");
                        }
                        overwrite_source_c_buffer(title.as_deref_mut(), &bytes)?;
                    }
                }
            }
        }

        if result == 0 && !all_normalized.is_null() {
            for representation in 0..TAUT_NUM as usize {
                if normalized_data[representation].bExists != 0 {
                    heap.slice_mut(all_normalized)?[component][representation] =
                        std::mem::take(&mut normalized_data[representation]);
                }
            }
        }
        if result != 0 {
            result = TreatErrorsInCreateOneComponentINChI(
                heap,
                structure_data,
                input_parameters,
                &prepared_input[kind],
                component as i32,
                input_number,
                input_file.as_deref_mut(),
                log_file.as_deref_mut(),
                output_file.as_deref_mut(),
                problem_file.as_deref_mut(),
            )?;
            break;
        }
    }

    if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 && !all_normalized.is_null() {
        let all = heap.slice(all_normalized.as_const())?.to_vec();
        let _ = CreateCompositeNormAtom(
            heap,
            &mut composite_normalized_data[kind],
            &all,
            component_count_i32,
        )?;
    }
    free_structure_normalized_data(heap, all_normalized, component_count)?;
    FreeInpAtomData(heap, &mut current_input)?;
    for normalized in &mut normalized_data {
        FreeInpAtomData(heap, normalized)?;
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CreateOneComponentINChI(
    heap: &mut SourceHeap,
    canonical_globals: &mut CANON_GLOBALS,
    clock: SourceMutPointer<INCHI_CLOCK>,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &mut INPUT_PARMS,
    current_input: &mut INP_ATOM_DATA,
    original_input: Option<&ORIG_ATOM_DATA>,
    inchi_rows: &mut [PINChI2],
    aux_rows: &mut [PINChI_Aux2],
    inchi_kind: i32,
    component_index: i32,
    _input_number: i64,
    normalized_data: &mut [INP_ATOM_DATA; TAUT_NUM as usize],
    normalization_flags: &mut NORM_CANON_FLAGS,
    _log_file: &mut INCHI_IOSTREAM,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:1747 CreateOneComponentINChI
    // INCHI✔️❌: complete source frame follows verbatim; SourceHeap pointer validation and modeled allocations add overhead.
    /*
    int CreateOneComponentINChI( CANON_GLOBALS      *pCG,
                                 INCHI_CLOCK        *ic,
                                 STRUCT_DATA        *sd,
                                 INPUT_PARMS        *ip,
                                 INP_ATOM_DATA      *inp_cur_data,
                                 ORIG_ATOM_DATA     *orig_inp_data,
                                 PINChI2            *pINChI,
                                 PINChI_Aux2        *pINChI_Aux,
                                 int                iINChI,
                                 int                i,
                                 long               num_inp,
                                 INP_ATOM_DATA      **inp_norm_data,
                                 NORM_CANON_FLAGS   *pncFlags,
                                 INCHI_IOSTREAM     *log_file )
    {
        inchiTime     ulTStart, ulTEnd, *pulTEnd = NULL;
        int           k, num_at, ret = 0;
        int           bOrigCoord;
        INCHI_MODE     bTautFlags = ip->bTautFlags;
        INCHI_MODE     bTautFlagsDone = ( ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS] );
        INChI       *cur_INChI[TAUT_NUM];
        INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
        long          lElapsedTime;

        int nAllocMode = 0;  /* moved from below 2024-09-01 DT */

    #ifdef GHI100_FIX
    #if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
        setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
    #endif
    #endif

        InchiTimeGet( &ulTStart );
        bOrigCoord =
            !( ip->bINChIOutputOptions & ( INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO ) );

        for (k = 0; k < TAUT_NUM; k++)
        {
            cur_INChI[k] = NULL;
            cur_INChI_Aux[k] = NULL;
        }

        /*  Allocate memory for non-tautomeric (k=0) and tautomeric (k=1) results */
        for (k = 0; k < TAUT_NUM; k++)
        {
            /* djb-rwth: introducing variables for correct nAllocMode expression */
            int nAM1 = 0, nAM2 = 0;

            if (k == TAUT_YES)
                nAM1 = REQ_MODE_TAUT;

            if (bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))
                nAM2 = ip->nMode & REQ_MODE_ISO;

            nAllocMode = nAM1 | nAM2; /* djb-rwth: original sequence of bit-wise operations had to be rewritten */


            if ((k == TAUT_NON && ( ip->nMode & REQ_MODE_BASIC )) ||
                 (k == TAUT_YES && ( ip->nMode & REQ_MODE_TAUT ))) /* djb-rwth: addressing LLVM warning */
            {
                /*  alloc INChI and INChI_Aux */
                cur_INChI[k] = Alloc_INChI( inp_cur_data->at,
                                                inp_cur_data->num_at,
                                                &inp_cur_data->num_bonds,
                                                &inp_cur_data->num_isotopic,
                                                nAllocMode );

                cur_INChI_Aux[k] = Alloc_INChI_Aux( inp_cur_data->num_at,
                                                    inp_cur_data->num_isotopic,
                                                    nAllocMode,
                                                    bOrigCoord );
                if (cur_INChI_Aux[k])
                {
                    cur_INChI_Aux[k]->bIsIsotopic = inp_cur_data->num_isotopic;
                }
                /*  alloc memory for the output structure: non-tautomeric and tautomeric (for displaying) */

                CreateInpAtomData( inp_norm_data[k], inp_cur_data->num_at, k );
            }
            else
            {
                FreeInpAtomData( inp_norm_data[k] );
            }
        }

        lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
        if (ip->msec_MaxTime)
        {
            ip->msec_LeftTime -= lElapsedTime;
        }
        sd->ulStructTime += lElapsedTime;

        /*^^^#if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) ) */
        #if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) && !defined(TARGET_EXE_STANDALONE) )
        #endif

        /******************************************************
         *
         *  Get one component canonical numberings, etc.
         *
         ******************************************************/

        /* Create_INChI() return value:
         * num_at <= 0: error code
         * num_at >  0: number of atoms (excluding terminal hydrogen atoms)
         * inp_norm_data[0] => non-tautomeric, inp_norm_data[1] => tautomeric */

        InchiTimeGet( &ulTStart );

        if (ip->msec_MaxTime)
        {
            ulTEnd = ulTStart;
            pulTEnd = &ulTEnd;
            if (ip->msec_LeftTime > 0)
            {
                InchiTimeAddMsec( ic, pulTEnd, ip->msec_LeftTime );
            }
        }

        num_at = Create_INChI( pCG, ic, ip,
                               cur_INChI, cur_INChI_Aux,
                               orig_inp_data/* not used */,
                               inp_cur_data->at, inp_norm_data, inp_cur_data->num_at,
                               ip->nMode,
                               &bTautFlags, &bTautFlagsDone,
                               pulTEnd, NULL, sd->pStrErrStruct );

        SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, i + 1 );
                            /*  NB: normalization alters structure component number */

        for (k = 0; k < TAUT_NUM; k++)
        {
            if (cur_INChI_Aux[k] && cur_INChI_Aux[k]->nNumberOfAtoms > 0)
            {
                pncFlags->bNormalizationFlags[iINChI][k] |=
                    cur_INChI_Aux[k]->bNormalizationFlags;
                pncFlags->bTautFlags[iINChI][k] |=
                    cur_INChI_Aux[k]->bTautFlags;
                pncFlags->bTautFlagsDone[iINChI][k] |=
                    cur_INChI_Aux[k]->bTautFlagsDone;
                pncFlags->nCanonFlags[iINChI][k] |=
                    cur_INChI_Aux[k]->nCanonFlags;
            }
        }

        /*  Detect errors */
        if (num_at < 0)
        {
            sd->nErrorCode = num_at;
        }
        else if (num_at == 0)
        {
            sd->nErrorCode = -1;
        }
        else if (cur_INChI[TAUT_NON] && cur_INChI[TAUT_NON]->nErrorCode)
        {   /*  non-tautomeric error */
            sd->nErrorCode = cur_INChI[TAUT_NON]->nErrorCode;
        }
        else if (cur_INChI[TAUT_YES] && cur_INChI[TAUT_YES]->nErrorCode)
        {   /*  tautomeric error */
            sd->nErrorCode = cur_INChI[TAUT_YES]->nErrorCode;
        }


    #if ( bRELEASE_VERSION == 0 )
        if (cur_INChI[TAUT_NON]) sd->bExtract |= cur_INChI[TAUT_NON]->bExtract;
        if (cur_INChI[TAUT_YES]) sd->bExtract |= cur_INChI[TAUT_YES]->bExtract;
        if (( TG_FLAG_TEST_TAUT3_SALTS_DONE & bTautFlagsDone ))
        {
            sd->bExtract |= EXTR_TEST_TAUT3_SALTS_DONE;
        }
    #endif

        /*  Detect and store stereo warnings */
        if (!sd->nErrorCode)
            GetProcessingWarningsOneComponentInChI( cur_INChI, inp_norm_data, sd, ip->bNoWarnings );

        lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
        if (ip->msec_MaxTime)
        {
            ip->msec_LeftTime -= lElapsedTime;
        }
        sd->ulStructTime += lElapsedTime;

    #if !defined(TARGET_API_LIB) && !defined(COMPILE_ANSI_ONLY)
        /*  Display the results */
        if (ip->bDisplay)
        {
            eat_keyboard_input( );
        }
    #endif

        /*  a) No matter what happened save the allocated INChI pointers */
        /*  save the INChI of the current component */

        InchiTimeGet( &ulTStart );
        for (k = 0; k < TAUT_NUM; k++)
        {
            pINChI[i][k] = cur_INChI[k];
            pINChI_Aux[i][k] = cur_INChI_Aux[k];
            cur_INChI[k] = NULL;
            cur_INChI_Aux[k] = NULL;
        }

        /*  b) Count one component structure and/or INChI results only
               if there was no error
               Set inp_norm_data[j]->num_removed_H = number of removed explicit H
        */

        if (!sd->nErrorCode)
        {
            /*  find where the current processed structure is located */
            int cur_is_in_non_taut = ( pINChI[i][TAUT_NON] &&
                                       pINChI[i][TAUT_NON]->nNumberOfAtoms > 0 );
            int cur_is_in_taut = ( pINChI[i][TAUT_YES] &&
                                       pINChI[i][TAUT_YES]->nNumberOfAtoms > 0 );

            int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[i][TAUT_NON]->lenTautomer) ||
                (cur_is_in_taut && 0 == pINChI[i][TAUT_YES]->lenTautomer); /* djb-rwth: addressing LLVM warnings */
            int cur_is_taut = cur_is_in_taut && 0 < pINChI[i][TAUT_YES]->lenTautomer;

            /*
            sd->bTautFlags[iINChI]     |= bTautFlags;
            sd->bTautFlagsDone[iINChI] |= bTautFlagsDone;
            */

            if (cur_is_non_taut + cur_is_taut)
            {
                /*  count tautomeric and non-tautomeric components of the structures */
                int j1 = cur_is_in_non_taut ? TAUT_NON : TAUT_YES;
                int j2 = cur_is_in_taut ? TAUT_YES : TAUT_NON;
                int j;
                sd->num_non_taut[iINChI] += cur_is_non_taut;
                sd->num_taut[iINChI] += cur_is_taut;
                for (j = j1; j <= j2; j++)
                {
                    int bIsotopic = ( pINChI[i][j]->nNumberOfIsotopicAtoms ||
                                      pINChI[i][j]->nNumberOfIsotopicTGroups ||
                                      (pINChI[i][j]->nPossibleLocationsOfIsotopicH && pINChI[i][j]->nPossibleLocationsOfIsotopicH[0] > 1) ); /* djb-rwth: addressing LLVM warning */
                    if (pINChI_Aux[i][j] && (j == TAUT_YES)) /* djb-rwth: fixing a NULL pointer dereference */
                    {
                        bIsotopic |= ( 0 < pINChI_Aux[i][j]->nNumRemovedIsotopicH[0] +
                                          pINChI_Aux[i][j]->nNumRemovedIsotopicH[1] +
                                          pINChI_Aux[i][j]->nNumRemovedIsotopicH[2] );
                    }

                    inp_norm_data[j]->bExists = 1; /*  j=0: non-taut exists, j=1: taut exists */
                    inp_norm_data[j]->bHasIsotopicLayer = bIsotopic;
                    /*inp_norm_data[j]->num_removed_H = inp_norm_data[j]->num_at - num_at;*/
                }
            }
        }

        /* return (sd->nErrorCode==CT_OUT_OF_RAM || sd->nErrorCode==CT_USER_QUIT_ERR)? _IS_FATAL :
                sd->nErrorCode? _IS_ERROR : 0; */

        if (sd->nErrorCode == CT_OUT_OF_RAM || sd->nErrorCode == CT_USER_QUIT_ERR)
        {
            ret = _IS_FATAL;
        }
        else if (sd->nErrorCode)
        {
            ret = _IS_ERROR;
        }

        lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
        if (ip->msec_MaxTime)
        {
            ip->msec_LeftTime -= lElapsedTime;
        }
        sd->ulStructTime += lElapsedTime;

        return ret;
    }
    */
    // END INCHI C FUNCTION: CreateOneComponentINChI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateOneComponentINChI
    // INCHI✔️❌: #define COMPILE_ANSI_ONLY
    // INCHI✔️❌: #define TARGET_API_LIB
    // INCHI✔️❌: GCC/Linux; bRELEASE_VERSION != 0; GHI100_FIX inactive
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateOneComponentINChI
    let mut start = inchiTime::default();
    InchiTimeGet(&mut start, clock_result);
    let original_coordinates = i32::from(
        input_parameters.bINChIOutputOptions
            & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO) as i32
            == 0,
    );
    let mut taut_flags: INCHI_MODE = input_parameters.bTautFlags;
    let kind = usize::try_from(inchi_kind)
        .ok()
        .filter(|index| *index < INCHI_NUM as usize)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut taut_flags_done =
        input_parameters.bTautFlagsDone | structure_data.bTautFlagsDone[INCHI_BAS as usize];
    let mut current_inchi = [SourceMutPointer::null(); TAUT_NUM as usize];
    let mut current_aux = [SourceMutPointer::null(); TAUT_NUM as usize];
    let atom_count = usize::try_from(current_input.num_at.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let atoms = if current_input.at.is_null() {
        Vec::new()
    } else {
        heap.slice(current_input.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };

    for representation in 0..TAUT_NUM as usize {
        let taut_mode = if representation == TAUT_YES as usize {
            REQ_MODE_TAUT as i32
        } else {
            0
        };
        let isotope_mode = if taut_flags_done
            & u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
            != 0
        {
            input_parameters.nMode as i32 & REQ_MODE_ISO as i32
        } else {
            0
        };
        let allocation_mode = taut_mode | isotope_mode;
        let enabled = (representation == TAUT_NON as usize
            && input_parameters.nMode & u64::from(REQ_MODE_BASIC) != 0)
            || (representation == TAUT_YES as usize
                && input_parameters.nMode & u64::from(REQ_MODE_TAUT) != 0);
        if enabled {
            current_inchi[representation] = Alloc_INChI(
                heap,
                &atoms,
                current_input.num_at,
                &mut current_input.num_bonds,
                &mut current_input.num_isotopic,
                allocation_mode,
            )?;
            current_aux[representation] = Alloc_INChI_Aux(
                heap,
                current_input.num_at,
                current_input.num_isotopic,
                allocation_mode,
                original_coordinates,
            )?;
            if !current_aux[representation].is_null() {
                heap.slice_mut(current_aux[representation])?[0].bIsIsotopic =
                    current_input.num_isotopic;
            }
            let _ = CreateInpAtomData(
                heap,
                &mut normalized_data[representation],
                current_input.num_at,
                representation as i32,
            )?;
        } else {
            FreeInpAtomData(heap, &mut normalized_data[representation])?;
        }
    }

    let mut elapsed = component_elapsed(heap, clock, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    InchiTimeGet(&mut start, clock_result);
    let maximum_time = if input_parameters.msec_MaxTime != 0 {
        let pointer = heap.allocate_model_storage(vec![start.clone()])?;
        if input_parameters.msec_LeftTime > 0 {
            let mut clock_value = heap.slice(clock.as_const())?[0].clone();
            InchiTimeAddMsec(
                &mut clock_value,
                heap.slice_mut(pointer)?.first_mut(),
                input_parameters.msec_LeftTime as u64,
            );
            heap.slice_mut(clock)?[0] = clock_value;
        }
        pointer
    } else {
        SourceMutPointer::null()
    };
    let create_result = Create_INChI(
        heap,
        canonical_globals,
        clock,
        input_parameters,
        current_inchi,
        current_aux,
        original_input,
        current_input.at,
        normalized_data,
        current_input.num_at,
        input_parameters.nMode,
        &mut taut_flags,
        &mut taut_flags_done,
        maximum_time,
        None,
        Some(&mut structure_data.pStrErrStruct),
        clock_result,
    );
    if !maximum_time.is_null() {
        heap.free(maximum_time)?;
    }
    let num_atoms = create_result?;

    if !current_input.at.is_null() {
        SetConnectedComponentNumber(
            heap.slice_mut(current_input.at)?,
            current_input.num_at,
            component_index.wrapping_add(1),
        )?;
    }
    for representation in 0..TAUT_NUM as usize {
        if !current_aux[representation].is_null() {
            let aux = &heap.slice(current_aux[representation].as_const())?[0];
            if aux.nNumberOfAtoms > 0 {
                normalization_flags.bNormalizationFlags[kind][representation] |=
                    aux.bNormalizationFlags;
                normalization_flags.bTautFlags[kind][representation] |= aux.bTautFlags;
                normalization_flags.bTautFlagsDone[kind][representation] |= aux.bTautFlagsDone;
                normalization_flags.nCanonFlags[kind][representation] |= aux.nCanonFlags;
            }
        }
    }
    if num_atoms < 0 {
        structure_data.nErrorCode = num_atoms;
    } else if num_atoms == 0 {
        structure_data.nErrorCode = -1;
    } else if !current_inchi[TAUT_NON as usize].is_null()
        && heap.slice(current_inchi[TAUT_NON as usize].as_const())?[0].nErrorCode != 0
    {
        structure_data.nErrorCode =
            heap.slice(current_inchi[TAUT_NON as usize].as_const())?[0].nErrorCode;
    } else if !current_inchi[TAUT_YES as usize].is_null()
        && heap.slice(current_inchi[TAUT_YES as usize].as_const())?[0].nErrorCode != 0
    {
        structure_data.nErrorCode =
            heap.slice(current_inchi[TAUT_YES as usize].as_const())?[0].nErrorCode;
    }
    if structure_data.nErrorCode == 0 {
        let _ = GetProcessingWarningsOneComponentInChI(
            heap,
            &current_inchi,
            normalized_data,
            structure_data,
            input_parameters.bNoWarnings,
        )?;
    }
    elapsed = component_elapsed(heap, clock, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    InchiTimeGet(&mut start, clock_result);
    let row = usize::try_from(component_index)
        .ok()
        .filter(|index| *index < inchi_rows.len() && *index < aux_rows.len())
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    inchi_rows[row] = current_inchi;
    aux_rows[row] = current_aux;

    if structure_data.nErrorCode == 0 {
        let non_taut_in = !inchi_rows[row][TAUT_NON as usize].is_null()
            && heap.slice(inchi_rows[row][TAUT_NON as usize].as_const())?[0].nNumberOfAtoms > 0;
        let taut_in = !inchi_rows[row][TAUT_YES as usize].is_null()
            && heap.slice(inchi_rows[row][TAUT_YES as usize].as_const())?[0].nNumberOfAtoms > 0;
        let non_taut = (non_taut_in
            && heap.slice(inchi_rows[row][TAUT_NON as usize].as_const())?[0].lenTautomer == 0)
            || (taut_in
                && heap.slice(inchi_rows[row][TAUT_YES as usize].as_const())?[0].lenTautomer == 0);
        let taut = taut_in
            && heap.slice(inchi_rows[row][TAUT_YES as usize].as_const())?[0].lenTautomer > 0;
        if non_taut || taut {
            structure_data.num_non_taut[kind] =
                structure_data.num_non_taut[kind].wrapping_add(i32::from(non_taut));
            structure_data.num_taut[kind] =
                structure_data.num_taut[kind].wrapping_add(i32::from(taut));
            let first = if non_taut_in { TAUT_NON } else { TAUT_YES } as usize;
            let last = if taut_in { TAUT_YES } else { TAUT_NON } as usize;
            for representation in first..=last {
                let inchi = &heap.slice(inchi_rows[row][representation].as_const())?[0];
                let possible_isotope_h = !inchi.nPossibleLocationsOfIsotopicH.is_null()
                    && heap.slice(inchi.nPossibleLocationsOfIsotopicH.as_const())?[0] > 1;
                let mut isotopic = inchi.nNumberOfIsotopicAtoms != 0
                    || inchi.nNumberOfIsotopicTGroups != 0
                    || possible_isotope_h;
                if representation == TAUT_YES as usize && !aux_rows[row][representation].is_null() {
                    let aux = &heap.slice(aux_rows[row][representation].as_const())?[0];
                    isotopic |= aux
                        .nNumRemovedIsotopicH
                        .iter()
                        .fold(0_i32, |sum, value| sum.wrapping_add(i32::from(*value)))
                        > 0;
                }
                normalized_data[representation].bExists = 1;
                normalized_data[representation].bHasIsotopicLayer = i32::from(isotopic);
            }
        }
    }
    let result = if structure_data.nErrorCode == CT_OUT_OF_RAM
        || structure_data.nErrorCode == CT_USER_QUIT_ERR
    {
        _IS_FATAL as i32
    } else if structure_data.nErrorCode != 0 {
        _IS_ERROR as i32
    } else {
        0
    };
    elapsed = component_elapsed(heap, clock, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        INCHI_IOS_TYPE_FILE, INChI, INT_ARRAY, OAD_Polymer, OAD_PolymerUnit, POLYMER_CONN_HT,
        POLYMER_SST_NON, POLYMER_STY_SRU, POLYMERS_LEGACY, SourceFile, inp_ATOM,
    };

    fn string_stream() -> INCHI_IOSTREAM {
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING::default(),
            ..INCHI_IOSTREAM::default()
        }
    }

    fn stream_bytes(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> Vec<u8> {
        if stream.s.pStr.is_null() {
            Vec::new()
        } else {
            heap.slice(stream.s.pStr.as_const()).unwrap()[..stream.s.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }
    }

    fn c_text(heap: &mut SourceHeap, bytes: &[u8]) -> SourceMutPointer<i8> {
        heap.allocate_model_storage(
            bytes
                .iter()
                .copied()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect(),
        )
        .unwrap()
    }

    fn array_c_bytes(bytes: &[i8]) -> Vec<u8> {
        bytes[..bytes.iter().position(|byte| *byte == 0).unwrap()]
            .iter()
            .map(|byte| *byte as u8)
            .collect()
    }

    fn structure_carbon(component: u16, original_number: u16) -> inp_ATOM {
        let mut atom = inp_ATOM {
            el_number: 6,
            num_H: 4,
            component,
            orig_at_number: original_number,
            ..inp_ATOM::default()
        };
        atom.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
        atom
    }

    #[test]
    fn source_port__runichi__swap_atoms_xyz__line_2401() {
        fn coordinate_bits(atom: &inp_ATOM) -> [u64; 3] {
            [atom.x.to_bits(), atom.y.to_bits(), atom.z.to_bits()]
        }

        fn without_coordinates(mut atom: inp_ATOM) -> inp_ATOM {
            atom.x = 0.0;
            atom.y = 0.0;
            atom.z = 0.0;
            atom
        }

        let mut heap = SourceHeap::default();
        assert_eq!(swap_atoms_xyz(&mut heap, None, i32::MIN, i32::MIN), Ok(()));
        assert_eq!(
            swap_atoms_xyz(&mut heap, None, 0, 1),
            Err(SourceHeapError::NullPointer)
        );

        let null_atoms = ORIG_ATOM_DATA::default();
        assert_eq!(swap_atoms_xyz(&mut heap, Some(&null_atoms), 7, 7), Ok(()));
        assert_eq!(
            swap_atoms_xyz(&mut heap, Some(&null_atoms), 0, 1),
            Err(SourceHeapError::NullPointer)
        );

        let mut first = inp_ATOM {
            x: f64::from_bits(0x7ff8_0000_0000_0042),
            y: -0.0,
            z: f64::INFINITY,
            charge: -3,
            radical: 2,
            orig_at_number: 17,
            ..inp_ATOM::default()
        };
        first.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
        first.neighbor[0] = 1;
        first.valence = 1;
        let mut second = inp_ATOM {
            x: f64::from_bits(0xfff8_0000_0000_0084),
            y: 0.0,
            z: f64::NEG_INFINITY,
            charge: 4,
            radical: 3,
            orig_at_number: 29,
            ..inp_ATOM::default()
        };
        second.elname[..2].copy_from_slice(&[b'N' as i8, 0]);
        second.neighbor[0] = 0;
        second.valence = 1;
        let first_before = first.clone();
        let second_before = second.clone();
        let atoms = heap.allocate_model_storage(vec![first, second]).unwrap();
        let coordinates = heap
            .allocate_model_storage(vec![[b'1' as i8; 32], [b'2' as i8; 32]])
            .unwrap();
        let data = ORIG_ATOM_DATA {
            at: atoms,
            szCoord: coordinates,
            num_inp_atoms: 2,
            ..ORIG_ATOM_DATA::default()
        };

        assert_eq!(swap_atoms_xyz(&mut heap, Some(&data), 0, 0), Ok(()));
        assert_eq!(
            coordinate_bits(&heap.slice(atoms.as_const()).unwrap()[0]),
            coordinate_bits(&first_before)
        );

        assert_eq!(swap_atoms_xyz(&mut heap, Some(&data), 0, 1), Ok(()));
        let swapped = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            coordinate_bits(&swapped[0]),
            coordinate_bits(&second_before)
        );
        assert_eq!(coordinate_bits(&swapped[1]), coordinate_bits(&first_before));
        assert_eq!(
            without_coordinates(swapped[0].clone()),
            without_coordinates(first_before.clone())
        );
        assert_eq!(
            without_coordinates(swapped[1].clone()),
            without_coordinates(second_before.clone())
        );
        assert_eq!(
            heap.slice(coordinates.as_const()).unwrap(),
            &[[b'1' as i8; 32], [b'2' as i8; 32]]
        );

        assert_eq!(swap_atoms_xyz(&mut heap, Some(&data), 1, 0), Ok(()));
        let restored = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            coordinate_bits(&restored[0]),
            coordinate_bits(&first_before)
        );
        assert_eq!(
            coordinate_bits(&restored[1]),
            coordinate_bits(&second_before)
        );

        for (first_index, second_index) in [(-1, 0), (0, -1), (0, 2), (2, 0)] {
            let before = heap.slice(atoms.as_const()).unwrap().to_vec();
            assert_eq!(
                swap_atoms_xyz(&mut heap, Some(&data), first_index, second_index),
                Err(SourceHeapError::PointerOutOfBounds)
            );
            let after = heap.slice(atoms.as_const()).unwrap();
            for (actual, expected) in after.iter().zip(&before) {
                assert_eq!(coordinate_bits(actual), coordinate_bits(expected));
                assert_eq!(
                    without_coordinates(actual.clone()),
                    without_coordinates(expected.clone())
                );
            }
        }
    }

    #[test]
    fn source_port__runichi__set_renumbered_or_delete__line_2775() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                SourceConstPointer::null(),
                i32::MIN,
            ),
            Ok(0)
        );
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                -1,
                SourceConstPointer::null(),
                0,
            ),
            Err(SourceHeapError::AllocationElementCountOutOfRange)
        );

        let numbers0 = heap.allocate_model_storage(vec![2, 0, 1, 99]).unwrap();
        let buffer0 = heap.allocate_model_storage(vec![7, 7, 7, 88]).unwrap();
        let renumbering0 = heap.allocate_model_storage(vec![1, -1, 0]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(&mut heap, numbers0, buffer0, 3, renumbering0.as_const(), 0,),
            Ok(2)
        );
        assert_eq!(heap.slice(numbers0.as_const()).unwrap(), &[0, 1, 0, 99]);
        assert_eq!(heap.slice(buffer0.as_const()).unwrap(), &[2, 0, 1, 88]);
        assert_eq!(heap.slice(renumbering0.as_const()).unwrap(), &[1, -1, 0]);

        let numbers1 = heap.allocate_model_storage(vec![3, 1, 2, 77]).unwrap();
        let buffer1 = heap.allocate_model_storage(vec![0; 4]).unwrap();
        let renumbering1 = heap.allocate_model_storage(vec![1, -1, 0]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(&mut heap, numbers1, buffer1, 3, renumbering1.as_const(), 1,),
            Ok(2)
        );
        assert_eq!(heap.slice(numbers1.as_const()).unwrap(), &[1, 2, 0, 77]);
        assert_eq!(heap.slice(buffer1.as_const()).unwrap(), &[3, 1, 2, 0]);

        let unusual_numbers = heap.allocate_model_storage(vec![-3, -2, 41]).unwrap();
        let unusual_buffer = heap.allocate_model_storage(vec![0; 3]).unwrap();
        let unusual_renumbering = heap.allocate_model_storage(vec![10, -1]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                unusual_numbers,
                unusual_buffer,
                2,
                unusual_renumbering.as_const(),
                -3,
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(unusual_numbers.as_const()).unwrap(), &[7, 0, 41]);
        assert_eq!(heap.slice(unusual_buffer.as_const()).unwrap(), &[-3, -2, 0]);

        let invalid_numbers = heap.allocate_model_storage(vec![3, 2, 1]).unwrap();
        let invalid_buffer = heap.allocate_model_storage(vec![9, 9, 9]).unwrap();
        let short_renumbering = heap.allocate_model_storage(vec![0]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                invalid_numbers,
                invalid_buffer,
                3,
                short_renumbering.as_const(),
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(invalid_numbers.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(invalid_buffer.as_const()).unwrap(), &[3, 2, 1]);

        let overflow_numbers = heap.allocate_model_storage(vec![0]).unwrap();
        let overflow_buffer = heap.allocate_model_storage(vec![7]).unwrap();
        let overflow_renumbering = heap.allocate_model_storage(vec![0]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                overflow_numbers,
                overflow_buffer,
                1,
                overflow_renumbering.as_const(),
                i32::MIN,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(heap.slice(overflow_numbers.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(overflow_buffer.as_const()).unwrap(), &[0]);

        let add_overflow_numbers = heap.allocate_model_storage(vec![1]).unwrap();
        let add_overflow_buffer = heap.allocate_model_storage(vec![0]).unwrap();
        let add_overflow_renumbering = heap.allocate_model_storage(vec![i32::MAX]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                add_overflow_numbers,
                add_overflow_buffer,
                1,
                add_overflow_renumbering.as_const(),
                1,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(heap.slice(add_overflow_numbers.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(add_overflow_buffer.as_const()).unwrap(), &[1]);

        let short_numbers = heap.allocate_model_storage(vec![1]).unwrap();
        let adequate_buffer = heap.allocate_model_storage(vec![5, 6]).unwrap();
        let adequate_renumbering = heap.allocate_model_storage(vec![0, 1]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                short_numbers,
                adequate_buffer,
                2,
                adequate_renumbering.as_const(),
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(adequate_buffer.as_const()).unwrap(), &[5, 6]);

        let adequate_numbers = heap.allocate_model_storage(vec![1, 2]).unwrap();
        let short_buffer = heap.allocate_model_storage(vec![5]).unwrap();
        assert_eq!(
            set_renumbered_or_delete(
                &mut heap,
                adequate_numbers,
                short_buffer,
                2,
                adequate_renumbering.as_const(),
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(adequate_numbers.as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(short_buffer.as_const()).unwrap(), &[5]);
    }

    #[test]
    fn source_port__runichi__validateandpreparepolymerandpseudoatoms__line_2879() {
        fn invoke(
            heap: &mut SourceHeap,
            structure: &mut STRUCT_DATA,
            parameters: &INPUT_PARMS,
            log: Option<&mut INCHI_IOSTREAM>,
            original: Option<&mut ORIG_ATOM_DATA>,
            input_number: i64,
            mind_polymers: &mut i32,
        ) -> Result<i32, SourceHeapError> {
            let mut inchi_components = [SourceMutPointer::<PINChI2>::null(); INCHI_NUM as usize];
            let mut aux_components = [SourceMutPointer::<PINChI_Aux2>::null(); INCHI_NUM as usize];
            ValidateAndPreparePolymerAndPseudoatoms(
                heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                structure,
                parameters,
                None,
                &mut inchi_components,
                &mut aux_components,
                None,
                log,
                None,
                None,
                original,
                None,
                input_number,
                None,
                u8::MAX,
                mind_polymers,
            )
        }

        fn named_atom(name: &[u8], element: u8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: element,
                ..inp_ATOM::default()
            };
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn polymer_case(
            heap: &mut SourceHeap,
        ) -> (ORIG_ATOM_DATA, SourceMutPointer<OAD_PolymerUnit>) {
            let atoms = heap
                .allocate_model_storage(vec![
                    named_atom(b"Zz\0", 0),
                    named_atom(b"C\0", 6),
                    named_atom(b"N\0", 7),
                    named_atom(b"Zz\0", 0),
                ])
                .unwrap();
            let mut bond_count = 0_i32;
            OrigAtData_AddBond(heap, 0, 1, atoms, 1, 0, &mut bond_count).unwrap();
            OrigAtData_AddBond(heap, 3, 2, atoms, 1, 0, &mut bond_count).unwrap();
            let alist = heap.allocate_model_storage(vec![2_i32, 3]).unwrap();
            let blist = heap.allocate_model_storage(vec![1_i32, 2, 4, 3]).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    type_: POLYMER_STY_SRU as i32,
                    subtype: POLYMER_SST_NON as i32,
                    conn: POLYMER_CONN_HT as i32,
                    na: 2,
                    nb: 2,
                    alist,
                    blist,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            let units = heap.allocate_model_storage(vec![unit]).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units,
                    n: 1,
                    treat: POLYMERS_LEGACY as i32,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            (
                ORIG_ATOM_DATA {
                    at: atoms,
                    num_inp_atoms: 4,
                    num_inp_bonds: bond_count,
                    polymer,
                    valid_polymer: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                unit,
            )
        }

        let mut heap = SourceHeap::default();
        let mut structure = STRUCT_DATA {
            nErrorCode: 71,
            ..STRUCT_DATA::default()
        };
        let parameters = INPUT_PARMS::default();
        let mut mind_polymers = 93_i32;
        assert_eq!(
            invoke(
                &mut heap,
                &mut structure,
                &parameters,
                None,
                None,
                i64::MIN,
                &mut mind_polymers,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(mind_polymers, 93);
        assert_eq!(structure.nErrorCode, 71);

        let mut ordinary = ORIG_ATOM_DATA {
            valid_polymer: 81,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            invoke(
                &mut heap,
                &mut structure,
                &parameters,
                None,
                Some(&mut ordinary),
                0,
                &mut mind_polymers,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(mind_polymers, 0);
        assert_eq!(ordinary.valid_polymer, 0);

        let pseudo_atoms = heap
            .allocate_model_storage(vec![named_atom(b"Zz\0", 0)])
            .unwrap();
        let mut invalid_pseudo = ORIG_ATOM_DATA {
            at: pseudo_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let label = c_text(&mut heap, b"ID");
        let value = c_text(&mut heap, b"sample-7");
        let error_parameters = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            nInputType: tagInputType_INPUT_SDFILE,
            ..INPUT_PARMS::default()
        };
        let mut log = string_stream();
        structure = STRUCT_DATA::default();
        assert_eq!(
            invoke(
                &mut heap,
                &mut structure,
                &error_parameters,
                Some(&mut log),
                Some(&mut invalid_pseudo),
                -17,
                &mut mind_polymers,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(mind_polymers, 0);
        assert_eq!(structure.nErrorCode, 75);
        assert_eq!(invalid_pseudo.num_inp_atoms, -1);
        assert_eq!(
            array_c_bytes(&structure.pStrErrStruct),
            b"Invalid element(s): Zz; Polymer-unrelated pseudoatoms are not allowed"
        );
        assert_eq!(
            stream_bytes(&heap, &log),
            b"Error 75 (Invalid element(s): Zz; Polymer-unrelated pseudoatoms are not allowed) structure #-17. ID=sample-7\n"
        );

        let allowed_atoms = heap
            .allocate_model_storage(vec![named_atom(b"Zz\0", 0)])
            .unwrap();
        let mut allowed_pseudo = ORIG_ATOM_DATA {
            at: allowed_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let allowed_parameters = INPUT_PARMS {
            bNPZz: 1,
            nInputType: tagInputType_INPUT_MOLFILE,
            ..INPUT_PARMS::default()
        };
        structure = STRUCT_DATA::default();
        assert_eq!(
            invoke(
                &mut heap,
                &mut structure,
                &allowed_parameters,
                None,
                Some(&mut allowed_pseudo),
                i64::MAX,
                &mut mind_polymers,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(mind_polymers, 0);
        assert_eq!(allowed_pseudo.n_zy, 1);
        assert_eq!(
            &heap.slice(allowed_pseudo.at.as_const()).unwrap()[0].elname[..3],
            &[b'Z' as i8, b'y' as i8, 0]
        );

        let legacy_parameters = INPUT_PARMS {
            nInputType: tagInputType_INPUT_MOLFILE,
            bPolymers: POLYMERS_LEGACY as i32,
            bNPZz: 1,
            bFrameShiftScheme: tagFrameShifScheme_FSS_STARS_CYCLED as i32,
            ..INPUT_PARMS::default()
        };
        let (mut legacy, legacy_unit) = polymer_case(&mut heap);
        structure = STRUCT_DATA::default();
        assert_eq!(
            invoke(
                &mut heap,
                &mut structure,
                &legacy_parameters,
                None,
                Some(&mut legacy),
                1,
                &mut mind_polymers,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(mind_polymers, 1);
        assert_eq!(legacy.num_inp_bonds, 1);
        assert_eq!(heap.slice(legacy_unit.as_const()).unwrap()[0].cyclized, 1);

        let (mut modern, modern_unit) = polymer_case(&mut heap);
        let modern_parameters = INPUT_PARMS {
            bPolymers: POLYMERS_MODERN as i32,
            ..legacy_parameters.clone()
        };
        assert_eq!(
            invoke(
                &mut heap,
                &mut structure,
                &modern_parameters,
                None,
                Some(&mut modern),
                2,
                &mut mind_polymers,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(mind_polymers, 1);
        assert_eq!(modern.num_inp_bonds, 2);
        assert_eq!(heap.slice(modern_unit.as_const()).unwrap()[0].cyclized, 0);

        let (mut inchi_input, inchi_unit) = polymer_case(&mut heap);
        let inchi_parameters = INPUT_PARMS {
            nInputType: tagInputType_INPUT_INCHI,
            ..legacy_parameters
        };
        assert_eq!(
            invoke(
                &mut heap,
                &mut structure,
                &inchi_parameters,
                None,
                Some(&mut inchi_input),
                3,
                &mut mind_polymers,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(mind_polymers, 0);
        assert_eq!(inchi_input.num_inp_bonds, 2);
        assert_eq!(heap.slice(inchi_unit.as_const()).unwrap()[0].cyclized, 0);
    }

    #[test]
    fn source_port__runichi__processonestructureexcore__line_2803() {
        fn invoke(
            heap: &mut SourceHeap,
            clock: SourceMutPointer<INCHI_CLOCK>,
            canonical: SourceMutPointer<CANON_GLOBALS>,
            structure: &mut STRUCT_DATA,
            parameters: &mut INPUT_PARMS,
            output: Option<&mut INCHI_IOSTREAM>,
            log: Option<&mut INCHI_IOSTREAM>,
            original: SourceMutPointer<ORIG_ATOM_DATA>,
            prepared: SourceMutPointer<ORIG_ATOM_DATA>,
            stdout: SourceMutPointer<SourceFile>,
        ) -> Result<i32, SourceHeapError> {
            let mut inchi_components = [SourceMutPointer::<PINChI2>::null(); INCHI_NUM as usize];
            let mut aux_components = [SourceMutPointer::<PINChI_Aux2>::null(); INCHI_NUM as usize];
            ProcessOneStructureExCore(
                heap,
                clock,
                canonical,
                structure,
                parameters,
                None,
                &mut inchi_components,
                &mut aux_components,
                None,
                log,
                output,
                None,
                original,
                prepared,
                -37,
                None,
                0xa5,
                stdout,
                0,
            )
        }

        fn source_polymer(heap: &mut SourceHeap) -> ORIG_ATOM_DATA {
            let atoms = heap
                .allocate_model_storage(vec![structure_carbon(0, 1)])
                .unwrap();
            let atom_list = heap.allocate_model_storage(vec![1_i32]).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    type_: crate::source_types::POLYMER_STY_MON as i32,
                    subtype: POLYMER_SST_NON as i32,
                    conn: crate::source_types::POLYMER_CONN_NON as i32,
                    na: 1,
                    nb: 0,
                    alist: atom_list,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            let units = heap.allocate_model_storage(vec![unit]).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units,
                    n: 1,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            ORIG_ATOM_DATA {
                at: atoms,
                num_inp_atoms: 1,
                polymer,
                valid_polymer: 1,
                ..ORIG_ATOM_DATA::default()
            }
        }

        let mut null_heap = SourceHeap::default();
        let mut null_structure = STRUCT_DATA {
            bUserQuitComponent: 73,
            ..STRUCT_DATA::default()
        };
        assert_eq!(
            invoke(
                &mut null_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut null_structure,
                &mut INPUT_PARMS::default(),
                None,
                None,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(null_structure.bUserQuitComponent, 73);

        let mut error_heap = SourceHeap::default();
        let mut invalid_atom = inp_ATOM::default();
        invalid_atom.elname[..3].copy_from_slice(&[b'Z' as i8, b'z' as i8, 0]);
        let invalid_atoms = error_heap
            .allocate_model_storage(vec![invalid_atom])
            .unwrap();
        let invalid_original = error_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA {
                at: invalid_atoms,
                num_inp_atoms: 1,
                ..ORIG_ATOM_DATA::default()
            }])
            .unwrap();
        let mut error_structure = STRUCT_DATA {
            bUserQuitComponent: 79,
            ..STRUCT_DATA::default()
        };
        let mut error_parameters = INPUT_PARMS {
            nInputType: tagInputType_INPUT_SDFILE,
            ..INPUT_PARMS::default()
        };
        let mut error_log = string_stream();
        assert_eq!(
            invoke(
                &mut error_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut error_structure,
                &mut error_parameters,
                None,
                Some(&mut error_log),
                invalid_original,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(error_structure.nErrorCode, 75);
        assert_eq!(error_structure.bUserQuitComponent, 79);
        assert_eq!(
            error_heap.slice(invalid_original.as_const()).unwrap()[0].num_inp_atoms,
            -1
        );
        assert_eq!(
            stream_bytes(&error_heap, &error_log),
            b"Error 75 (Invalid element(s): Zz; Polymer-unrelated pseudoatoms are not allowed) structure #-37.\n"
        );

        let mut skip_heap = SourceHeap::default();
        let skip_canonical = skip_heap
            .allocate_model_storage(vec![CANON_GLOBALS::default()])
            .unwrap();
        let skip_original_value = source_polymer(&mut skip_heap);
        let skip_original = skip_heap
            .allocate_model_storage(vec![skip_original_value])
            .unwrap();
        let skip_prepared = skip_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA::default(); INCHI_NUM as usize])
            .unwrap();
        let mut skip_parameters = INPUT_PARMS {
            nInputType: tagInputType_INPUT_MOLFILE,
            bPolymers: POLYMERS_LEGACY as i32,
            bNPZz: 1,
            bIgnoreUnchanged: 1,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            invoke(
                &mut skip_heap,
                SourceMutPointer::null(),
                skip_canonical,
                &mut STRUCT_DATA::default(),
                &mut skip_parameters,
                None,
                None,
                skip_original,
                skip_prepared,
                SourceMutPointer::null(),
            ),
            Ok(crate::source_types::_IS_SKIP as i32)
        );

        let mut ordinary_heap = SourceHeap::default();
        let ordinary_clock = ordinary_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let ordinary_canonical = ordinary_heap
            .allocate_model_storage(vec![CANON_GLOBALS::default()])
            .unwrap();
        let ordinary_original = ordinary_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA::default()])
            .unwrap();
        let ordinary_prepared = ordinary_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA::default(); INCHI_NUM as usize])
            .unwrap();
        let ordinary_stdout = ordinary_heap
            .allocate_model_storage(vec![SourceFile::default()])
            .unwrap();
        let mut ordinary_output = string_stream();
        let mut ordinary_structure = STRUCT_DATA {
            bChiralFlag: FLAG_INP_AT_CHIRAL as i32,
            bUserQuitComponent: 83,
            ..STRUCT_DATA::default()
        };
        let mut ordinary_parameters = INPUT_PARMS {
            bAllowEmptyStructure: 1,
            nMode: u64::from(REQ_MODE_STEREO),
            bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            invoke(
                &mut ordinary_heap,
                ordinary_clock,
                ordinary_canonical,
                &mut ordinary_structure,
                &mut ordinary_parameters,
                Some(&mut ordinary_output),
                None,
                ordinary_original,
                ordinary_prepared,
                ordinary_stdout,
            ),
            Ok(crate::source_types::_IS_WARNING as i32)
        );
        assert_eq!(ordinary_structure.bUserQuitComponent, 0);
        assert_eq!(
            array_c_bytes(&ordinary_structure.pStrErrStruct),
            b"Not chiral"
        );

        for (mode, expects_post_edit) in [
            (POLYMERS_MODERN as i32, false),
            (POLYMERS_LEGACY as i32, true),
            (POLYMERS_LEGACY_PLUS as i32, true),
        ] {
            let mut heap = SourceHeap::default();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let canonical = heap
                .allocate_model_storage(vec![CANON_GLOBALS::default()])
                .unwrap();
            let original_value = source_polymer(&mut heap);
            let original = heap.allocate_model_storage(vec![original_value]).unwrap();
            let prepared = heap
                .allocate_model_storage(vec![ORIG_ATOM_DATA::default(); INCHI_NUM as usize])
                .unwrap();
            let stdout = heap
                .allocate_model_storage(vec![SourceFile::default()])
                .unwrap();
            let mut structure = STRUCT_DATA::default();
            let mut parameters = INPUT_PARMS {
                nInputType: tagInputType_INPUT_MOLFILE,
                bPolymers: mode,
                bNPZz: 1,
                nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
                bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
                ..INPUT_PARMS::default()
            };
            let result = invoke(
                &mut heap,
                clock,
                canonical,
                &mut structure,
                &mut parameters,
                None,
                None,
                original,
                prepared,
                stdout,
            );
            if expects_post_edit {
                assert_eq!(result, Err(SourceHeapError::NullPointer));
            } else {
                assert!(
                    matches!(
                        result,
                        Ok(value)
                            if value == _IS_OKAY as i32
                                || value == crate::source_types::_IS_WARNING as i32
                    ),
                    "mode={mode} result={result:?} error_code={} error_type={} message={:?}",
                    structure.nErrorCode,
                    structure.nErrorType,
                    array_c_bytes(&structure.pStrErrStruct),
                );
            }
        }
    }

    #[test]
    fn source_port__runichi__oad_processonestructurenoedits__line_2967() {
        fn invoke(
            heap: &mut SourceHeap,
            original: &ORIG_ATOM_DATA,
            parameters: &INPUT_PARMS,
            supplied_inchi: Option<&[SourceMutPointer<PINChI2>; INCHI_NUM as usize]>,
            n_pzz: &mut i32,
            inchi: &mut SourceMutPointer<i8>,
            aux: &mut SourceMutPointer<i8>,
        ) -> Result<i32, SourceHeapError> {
            let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
            let canonical = heap.allocate_model_storage(vec![CANON_GLOBALS::default()])?;
            let stdout = heap.allocate_model_storage(vec![SourceFile {
                is_standard_stream: true,
                ..SourceFile::default()
            }])?;
            let prepared = [original.clone(), ORIG_ATOM_DATA::default()];
            OAD_ProcessOneStructureNoEdits(
                heap,
                clock,
                canonical,
                &STRUCT_DATA::default(),
                parameters,
                &[0],
                supplied_inchi,
                None,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                original,
                &prepared,
                i64::MIN,
                None,
                0xa5,
                n_pzz,
                inchi,
                aux,
                stdout,
                0,
            )
        }

        fn pointer_bytes(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Vec<u8> {
            let bytes = heap.slice(pointer.as_const()).unwrap();
            bytes[..bytes.iter().position(|byte| *byte == 0).unwrap()]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        fn complete_polymer(heap: &mut SourceHeap, n_pzz: i32) -> SourceMutPointer<OAD_Polymer> {
            if n_pzz == 0 {
                return heap
                    .allocate_model_storage(vec![OAD_Polymer {
                        treat: POLYMERS_LEGACY as i32,
                        ..OAD_Polymer::default()
                    }])
                    .unwrap();
            }
            let atom_list = heap.allocate_model_storage(vec![1_i32]).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    type_: crate::source_types::POLYMER_STY_MON as i32,
                    subtype: POLYMER_SST_NON as i32,
                    conn: crate::source_types::POLYMER_CONN_NON as i32,
                    na: 1,
                    nb: 0,
                    alist: atom_list,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            let units = heap.allocate_model_storage(vec![unit]).unwrap();
            let pzz = if n_pzz > 0 {
                heap.allocate_model_storage((1..=n_pzz).collect()).unwrap()
            } else {
                SourceMutPointer::null()
            };
            heap.allocate_model_storage(vec![OAD_Polymer {
                units,
                n: 1,
                n_pzz,
                pzz,
                treat: POLYMERS_LEGACY as i32,
                ..OAD_Polymer::default()
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let source_polymer = complete_polymer(&mut heap, 0);
        let original = ORIG_ATOM_DATA {
            polymer: source_polymer,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let parameters = INPUT_PARMS {
            nInputType: tagInputType_INPUT_MOLFILE,
            bPolymers: POLYMERS_LEGACY as i32,
            bFrameShiftScheme: tagFrameShifScheme_FSS_STARS_CYCLED as i32,
            bFoldPolymerSRU: 9,
            bDisplay: 8,
            bDisplayCompositeResults: 7,
            bDisplayEachComponentINChI: 6,
            bAllowEmptyStructure: 1,
            nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
            bINChIOutputOptions: (INCHI_OUT_NO_AUX_INFO
                | INCHI_OUT_SHORT_AUX_INFO
                | crate::source_types::INCHI_OUT_PLAIN_TEXT)
                as i32,
            ..INPUT_PARMS::default()
        };
        let old_inchi = c_text(&mut heap, b"old-inchi");
        let old_aux = c_text(&mut heap, b"old-aux");
        let mut n_pzz = -1;
        let mut inchi = old_inchi;
        let mut aux = old_aux;
        let result = invoke(
            &mut heap,
            &original,
            &parameters,
            None,
            &mut n_pzz,
            &mut inchi,
            &mut aux,
        );
        assert!(
            matches!(result, Ok(value) if value == _IS_OKAY as i32 || value == crate::source_types::_IS_WARNING as i32),
            "result={result:?} n_pzz={n_pzz} inchi_null={} aux_null={} inchi={:?} aux={:?}",
            inchi.is_null(),
            aux.is_null(),
            (!inchi.is_null()).then(|| pointer_bytes(&heap, inchi)),
            (!aux.is_null()).then(|| pointer_bytes(&heap, aux)),
        );
        assert_eq!(n_pzz, 0);
        assert_eq!(pointer_bytes(&heap, inchi), b"InChI=1//");
        assert_eq!(pointer_bytes(&heap, aux), b"AuxInfo=1//");
        assert_eq!(pointer_bytes(&heap, old_inchi), b"old-inchi");
        assert_eq!(pointer_bytes(&heap, old_aux), b"old-aux");
        assert_eq!(
            heap.slice(source_polymer.as_const()).unwrap()[0].treat,
            POLYMERS_LEGACY as i32
        );
        assert_eq!(parameters.bPolymers, POLYMERS_LEGACY as i32);
        assert_eq!(parameters.bFoldPolymerSRU, 9);
        assert_eq!(parameters.bDisplay, 8);
        assert_eq!(parameters.bINChIOutputOptions, 67);

        let mut reject_heap = SourceHeap::default();
        let reject_polymer = reject_heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let reject_original = ORIG_ATOM_DATA {
            polymer: reject_polymer,
            ..ORIG_ATOM_DATA::default()
        };
        let nonnull = reject_heap
            .allocate_model_storage(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        let supplied = [nonnull, SourceMutPointer::null()];
        n_pzz = 91;
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut reject_heap,
                &reject_original,
                &parameters,
                Some(&supplied),
                &mut n_pzz,
                &mut inchi,
                &mut aux,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(n_pzz, 0);
        assert!(inchi.is_null());
        assert!(aux.is_null());

        let mut failure_heap = SourceHeap::default();
        let failure_atoms = failure_heap
            .allocate_model_storage(vec![structure_carbon(0, 1)])
            .unwrap();
        let failure_polymer = complete_polymer(&mut failure_heap, 0);
        let failure_original = ORIG_ATOM_DATA {
            at: failure_atoms,
            num_inp_atoms: 1,
            polymer: failure_polymer,
            ..ORIG_ATOM_DATA::default()
        };
        failure_heap.fail_after_allocations(0);
        n_pzz = 92;
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut failure_heap,
                &failure_original,
                &parameters,
                None,
                &mut n_pzz,
                &mut inchi,
                &mut aux,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(n_pzz, 0);
        assert!(inchi.is_null());
        assert!(aux.is_null());
        assert!(failure_heap.slice(failure_atoms.as_const()).is_ok());

        let mut invalid_heap = SourceHeap::default();
        let mut invalid_atom = inp_ATOM::default();
        invalid_atom.elname[..3].copy_from_slice(&[b'Z' as i8, b'z' as i8, 0]);
        let invalid_atoms = invalid_heap
            .allocate_model_storage(vec![invalid_atom.clone()])
            .unwrap();
        let invalid_polymer = complete_polymer(&mut invalid_heap, 0);
        let invalid_original = ORIG_ATOM_DATA {
            at: invalid_atoms,
            num_inp_atoms: 1,
            polymer: invalid_polymer,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        };
        n_pzz = 93;
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut invalid_heap,
                &invalid_original,
                &INPUT_PARMS {
                    nInputType: tagInputType_INPUT_SDFILE,
                    ..parameters.clone()
                },
                None,
                &mut n_pzz,
                &mut inchi,
                &mut aux,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(n_pzz, 0);
        assert!(inchi.is_null());
        assert!(aux.is_null());
        assert_eq!(
            invalid_heap.slice(invalid_atoms.as_const()).unwrap()[0],
            invalid_atom
        );

        let mut null_polymer_heap = SourceHeap::default();
        n_pzz = 94;
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut null_polymer_heap,
                &ORIG_ATOM_DATA::default(),
                &parameters,
                None,
                &mut n_pzz,
                &mut inchi,
                &mut aux,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(n_pzz, 0);
        assert!(inchi.is_null());
        assert!(aux.is_null());
    }

    #[test]
    fn source_port__runichi__oad_processonestructure105plus__line_3062() {
        fn invoke(
            heap: &mut SourceHeap,
            original: &ORIG_ATOM_DATA,
            parameters: &INPUT_PARMS,
            supplied_inchi: Option<&[SourceMutPointer<PINChI2>; INCHI_NUM as usize]>,
            inchi: &mut SourceMutPointer<i8>,
            aux: &mut SourceMutPointer<i8>,
        ) -> Result<i32, SourceHeapError> {
            let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
            let canonical = heap.allocate_model_storage(vec![CANON_GLOBALS::default()])?;
            let stdout = heap.allocate_model_storage(vec![SourceFile {
                is_standard_stream: true,
                ..SourceFile::default()
            }])?;
            let prepared = [original.clone(), ORIG_ATOM_DATA::default()];
            OAD_ProcessOneStructure105Plus(
                heap,
                clock,
                canonical,
                &STRUCT_DATA::default(),
                parameters,
                &[0],
                supplied_inchi,
                None,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                original,
                &prepared,
                i64::MIN,
                None,
                0xa5,
                inchi,
                aux,
                stdout,
                0,
            )
        }

        fn pointer_bytes(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Vec<u8> {
            let bytes = heap.slice(pointer.as_const()).unwrap();
            bytes[..bytes.iter().position(|byte| *byte == 0).unwrap()]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        fn empty_polymer(heap: &mut SourceHeap) -> SourceMutPointer<OAD_Polymer> {
            heap.allocate_model_storage(vec![OAD_Polymer {
                treat: POLYMERS_MODERN as i32,
                ..OAD_Polymer::default()
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let source_polymer = empty_polymer(&mut heap);
        let original = ORIG_ATOM_DATA {
            polymer: source_polymer,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let parameters = INPUT_PARMS {
            nInputType: tagInputType_INPUT_MOLFILE,
            bPolymers: POLYMERS_MODERN as i32,
            bFrameShiftScheme: tagFrameShifScheme_FSS_STARS_CYCLED as i32,
            bFoldPolymerSRU: 9,
            bDisplay: 8,
            bDisplayCompositeResults: 7,
            bDisplayEachComponentINChI: 6,
            bAllowEmptyStructure: 1,
            nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
            bINChIOutputOptions: (INCHI_OUT_NO_AUX_INFO
                | INCHI_OUT_SHORT_AUX_INFO
                | crate::source_types::INCHI_OUT_PLAIN_TEXT)
                as i32,
            ..INPUT_PARMS::default()
        };
        let old_inchi = c_text(&mut heap, b"old-inchi");
        let old_aux = c_text(&mut heap, b"old-aux");
        let mut inchi = old_inchi;
        let mut aux = old_aux;
        let result = invoke(
            &mut heap,
            &original,
            &parameters,
            None,
            &mut inchi,
            &mut aux,
        );
        assert!(
            matches!(result, Ok(value) if value == _IS_OKAY as i32 || value == crate::source_types::_IS_WARNING as i32),
            "result={result:?} inchi_null={} aux_null={} inchi={:?} aux={:?}",
            inchi.is_null(),
            aux.is_null(),
            (!inchi.is_null()).then(|| pointer_bytes(&heap, inchi)),
            (!aux.is_null()).then(|| pointer_bytes(&heap, aux)),
        );
        assert_eq!(pointer_bytes(&heap, inchi), b"InChI=1//");
        assert_eq!(pointer_bytes(&heap, aux), b"AuxInfo=1//");
        assert_eq!(pointer_bytes(&heap, old_inchi), b"old-inchi");
        assert_eq!(pointer_bytes(&heap, old_aux), b"old-aux");
        assert_eq!(
            heap.slice(source_polymer.as_const()).unwrap()[0].treat,
            POLYMERS_MODERN as i32
        );
        assert_eq!(parameters.bPolymers, POLYMERS_MODERN as i32);
        assert_eq!(
            parameters.bFrameShiftScheme,
            tagFrameShifScheme_FSS_STARS_CYCLED as i32
        );
        assert_eq!(parameters.bFoldPolymerSRU, 9);
        assert_eq!(parameters.bDisplay, 8);
        assert_eq!(parameters.bDisplayCompositeResults, 7);
        assert_eq!(parameters.bDisplayEachComponentINChI, 6);
        assert_eq!(parameters.bINChIOutputOptions, 67);

        let mut reject_heap = SourceHeap::default();
        let reject_polymer = empty_polymer(&mut reject_heap);
        let reject_original = ORIG_ATOM_DATA {
            polymer: reject_polymer,
            ..ORIG_ATOM_DATA::default()
        };
        let nonnull = reject_heap
            .allocate_model_storage(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        let supplied = [nonnull, SourceMutPointer::null()];
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut reject_heap,
                &reject_original,
                &parameters,
                Some(&supplied),
                &mut inchi,
                &mut aux,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert!(inchi.is_null());
        assert!(aux.is_null());

        let mut failure_heap = SourceHeap::default();
        let failure_atoms = failure_heap
            .allocate_model_storage(vec![structure_carbon(0, 1)])
            .unwrap();
        let failure_polymer = empty_polymer(&mut failure_heap);
        let failure_original = ORIG_ATOM_DATA {
            at: failure_atoms,
            num_inp_atoms: 1,
            polymer: failure_polymer,
            ..ORIG_ATOM_DATA::default()
        };
        failure_heap.fail_after_allocations(0);
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut failure_heap,
                &failure_original,
                &parameters,
                None,
                &mut inchi,
                &mut aux,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert!(inchi.is_null());
        assert!(aux.is_null());
        assert!(failure_heap.slice(failure_atoms.as_const()).is_ok());

        let mut invalid_heap = SourceHeap::default();
        let mut invalid_atom = inp_ATOM::default();
        invalid_atom.elname[..3].copy_from_slice(&[b'Z' as i8, b'z' as i8, 0]);
        let invalid_atoms = invalid_heap
            .allocate_model_storage(vec![invalid_atom.clone()])
            .unwrap();
        let invalid_polymer = empty_polymer(&mut invalid_heap);
        let invalid_original = ORIG_ATOM_DATA {
            at: invalid_atoms,
            num_inp_atoms: 1,
            polymer: invalid_polymer,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        };
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut invalid_heap,
                &invalid_original,
                &INPUT_PARMS {
                    nInputType: tagInputType_INPUT_SDFILE,
                    ..parameters.clone()
                },
                None,
                &mut inchi,
                &mut aux,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert!(inchi.is_null());
        assert!(aux.is_null());
        assert_eq!(
            invalid_heap.slice(invalid_atoms.as_const()).unwrap()[0],
            invalid_atom
        );

        let mut null_polymer_heap = SourceHeap::default();
        inchi = old_inchi;
        aux = old_aux;
        assert_eq!(
            invoke(
                &mut null_polymer_heap,
                &ORIG_ATOM_DATA::default(),
                &parameters,
                None,
                &mut inchi,
                &mut aux,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert!(inchi.is_null());
        assert!(aux.is_null());
    }

    #[test]
    fn source_port__runichi__mark_atoms_to_delete_or_renumber__line_3155() {
        fn edits(
            heap: &mut SourceHeap,
            items: Vec<i32>,
            used: i32,
            increment: i32,
            side_chains: i32,
        ) -> OAD_StructureEdits {
            let allocated = items.len() as i32;
            let item = heap.allocate_model_storage(items).unwrap();
            let deletion = heap
                .allocate_model_storage(vec![INT_ARRAY {
                    item,
                    allocated,
                    used,
                    increment,
                }])
                .unwrap();
            OAD_StructureEdits {
                del_atom: deletion,
                del_side_chains: side_chains,
                ..OAD_StructureEdits::default()
            }
        }

        let mut heap = SourceHeap::default();
        let deletion_edits = edits(&mut heap, vec![2, 4, 99], 2, 3, 0);
        let renumbering = heap.allocate_model_storage(vec![77_i32; 6]).unwrap();
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: 5,
                    ..ORIG_ATOM_DATA::default()
                },
                &deletion_edits,
                renumbering,
            ),
            Ok(0)
        );
        assert_eq!(
            heap.slice(renumbering.as_const()).unwrap(),
            &[0, -1, 1, -1, 2, 77]
        );

        let ignored_edits = edits(&mut heap, vec![0, 6, 3], 3, 1, 0);
        let ignored_renumbering = heap.allocate_model_storage(vec![-8_i32; 5]).unwrap();
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: 4,
                    ..ORIG_ATOM_DATA::default()
                },
                &ignored_edits,
                ignored_renumbering,
            ),
            Ok(0)
        );
        assert_eq!(
            heap.slice(ignored_renumbering.as_const()).unwrap(),
            &[0, 1, -1, 2, -8]
        );

        let empty_edits = edits(&mut heap, Vec::new(), 0, 0, 0);
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA::default(),
                &empty_edits,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        let negative_count = heap.allocate_model_storage(vec![31_i32]).unwrap();
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: -1,
                    ..ORIG_ATOM_DATA::default()
                },
                &empty_edits,
                negative_count,
            ),
            Err(SourceHeapError::AllocationElementCountOutOfRange)
        );
        assert_eq!(heap.slice(negative_count.as_const()).unwrap(), &[31]);

        let one_atom = heap
            .allocate_model_storage(vec![inp_ATOM {
                orig_at_number: 1,
                valence: 0,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let side_edits = edits(&mut heap, vec![1, 91], 1, 2, 1);
        let side_renumbering = heap.allocate_model_storage(vec![44_i32, 45]).unwrap();
        let baseline = heap.live_allocation_count();
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA {
                    at: one_atom,
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &side_edits,
                side_renumbering,
            ),
            Ok(0)
        );
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(heap.slice(side_renumbering.as_const()).unwrap(), &[-1, 45]);
        let side_state = &heap.slice(side_edits.del_atom.as_const()).unwrap()[0];
        assert_eq!((side_state.used, side_state.allocated), (1, 2));
        assert_eq!(heap.slice(side_state.item.as_const()).unwrap(), &[1, 91]);

        let no_seed_edits = edits(&mut heap, vec![73], 0, 1, 1);
        let no_seed_renumbering = heap.allocate_model_storage(vec![88_i32]).unwrap();
        let baseline = heap.live_allocation_count();
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA {
                    at: one_atom,
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &no_seed_edits,
                no_seed_renumbering,
            ),
            Ok(0)
        );
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(heap.slice(no_seed_renumbering.as_const()).unwrap(), &[0]);

        let mut allocation_failure_heap = SourceHeap::default();
        let allocation_failure_edits = edits(&mut allocation_failure_heap, vec![1], 1, 1, 1);
        let allocation_failure_atoms = allocation_failure_heap
            .allocate_model_storage(vec![inp_ATOM {
                orig_at_number: 1,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let allocation_failure_renumbering = allocation_failure_heap
            .allocate_model_storage(vec![99_i32])
            .unwrap();
        let baseline = allocation_failure_heap.live_allocation_count();
        allocation_failure_heap.fail_after_allocations(0);
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut allocation_failure_heap,
                &ORIG_ATOM_DATA {
                    at: allocation_failure_atoms,
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &allocation_failure_edits,
                allocation_failure_renumbering,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(allocation_failure_heap.live_allocation_count(), baseline);
        assert_eq!(
            allocation_failure_heap
                .slice(allocation_failure_renumbering.as_const())
                .unwrap(),
            &[0]
        );

        let mut graph_failure_heap = SourceHeap::default();
        let graph_failure_edits = edits(&mut graph_failure_heap, vec![1], 1, 1, 1);
        let graph_failure_atoms = graph_failure_heap
            .allocate_model_storage(vec![inp_ATOM {
                orig_at_number: 1,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let graph_failure_renumbering = graph_failure_heap
            .allocate_model_storage(vec![99_i32])
            .unwrap();
        let baseline = graph_failure_heap.live_allocation_count();
        graph_failure_heap.fail_after_allocations(1);
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut graph_failure_heap,
                &ORIG_ATOM_DATA {
                    at: graph_failure_atoms,
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &graph_failure_edits,
                graph_failure_renumbering,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(graph_failure_heap.live_allocation_count(), baseline);

        let unsupported_edits = edits(&mut heap, vec![1, 2], 1, 2, 1);
        let unsupported_renumbering = heap.allocate_model_storage(vec![61_i32, 62]).unwrap();
        let baseline = heap.live_allocation_count();
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: 2,
                    ..ORIG_ATOM_DATA::default()
                },
                &unsupported_edits,
                unsupported_renumbering,
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(
            heap.slice(unsupported_renumbering.as_const()).unwrap(),
            &[0, 1]
        );

        let negative_used_edits = edits(&mut heap, vec![1], -1, 1, 0);
        let negative_used_renumbering = heap.allocate_model_storage(vec![90_i32]).unwrap();
        assert_eq!(
            mark_atoms_to_delete_or_renumber(
                &mut heap,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &negative_used_edits,
                negative_used_renumbering,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            heap.slice(negative_used_renumbering.as_const()).unwrap(),
            &[0]
        );
    }

    #[test]
    fn source_port__runichi__oad_structureedits_apply__line_2427() {
        fn int_array(
            heap: &mut SourceHeap,
            values: Vec<i32>,
            used: i32,
        ) -> SourceMutPointer<INT_ARRAY> {
            let allocated = values.len() as i32;
            let item = heap.allocate_model_storage(values).unwrap();
            heap.allocate_model_storage(vec![INT_ARRAY {
                item,
                allocated,
                used,
                increment: 1,
            }])
            .unwrap()
        }

        fn edits(
            heap: &mut SourceHeap,
            del_atom: (Vec<i32>, i32),
            del_bond: (Vec<i32>, i32),
            new_bond: (Vec<i32>, i32),
            mod_bond: (Vec<i32>, i32),
            mod_coord: (Vec<i32>, i32),
        ) -> OAD_StructureEdits {
            let del_atom = int_array(heap, del_atom.0, del_atom.1);
            let del_bond = int_array(heap, del_bond.0, del_bond.1);
            let new_bond = int_array(heap, new_bond.0, new_bond.1);
            let mod_bond = int_array(heap, mod_bond.0, mod_bond.1);
            let mod_coord = int_array(heap, mod_coord.0, mod_coord.1);
            OAD_StructureEdits {
                del_atom,
                del_bond,
                new_bond,
                mod_bond,
                mod_coord,
                del_side_chains: 0,
            }
        }

        fn original(
            heap: &mut SourceHeap,
            atom_count: usize,
            bonds: &[(i32, i32, i32, i32)],
        ) -> ORIG_ATOM_DATA {
            let atoms = (0..atom_count)
                .map(|index| inp_ATOM {
                    orig_at_number: (index + 1) as u16,
                    x: index as f64 + 0.25,
                    y: index as f64 + 0.5,
                    z: index as f64 + 0.75,
                    ..inp_ATOM::default()
                })
                .collect();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let mut number_of_bonds = 0;
            for &(atom1, atom2, bond_type, bond_stereo) in bonds {
                assert_eq!(
                    OrigAtData_AddBond(
                        heap,
                        atom1,
                        atom2,
                        atoms,
                        bond_type,
                        bond_stereo,
                        &mut number_of_bonds,
                    ),
                    Ok(1)
                );
            }
            ORIG_ATOM_DATA {
                at: atoms,
                num_inp_atoms: atom_count as i32,
                num_inp_bonds: number_of_bonds,
                ..ORIG_ATOM_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let empty = edits(
            &mut heap,
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
        );
        let mut empty_original = ORIG_ATOM_DATA::default();
        let mut status = 91;
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut empty_original,
                &empty,
                &mut status,
            ),
            Ok(0)
        );
        assert_eq!(status, _IS_OKAY as i32);

        let mut bond_original = original(&mut heap, 4, &[(0, 1, 2, 3), (2, 3, 1, 0)]);
        let crossing_bonds = heap.allocate_model_storage(vec![3_i32, 4, 9, 10]).unwrap();
        let unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                blist: crossing_bonds,
                nb: 2,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate_model_storage(vec![unit]).unwrap();
        bond_original.polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let bond_edits = edits(
            &mut heap,
            (Vec::new(), 0),
            (vec![1, 2], 2),
            (vec![1, 3, 77], 3),
            (vec![3, 4, 2, 4], 4),
            (Vec::new(), 0),
        );
        status = 91;
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut bond_original,
                &bond_edits,
                &mut status,
            ),
            Ok(3)
        );
        assert_eq!(status, _IS_OKAY as i32);
        assert_eq!(bond_original.num_inp_bonds, 2);
        let atoms = heap.slice(bond_original.at.as_const()).unwrap();
        assert_eq!((atoms[0].valence, atoms[0].neighbor[0]), (1, 2));
        assert_eq!((atoms[0].bond_type[0], atoms[0].bond_stereo[0]), (2, 3));
        assert_eq!((atoms[1].valence, atoms[1].neighbor[0]), (1, 3));
        assert_eq!((atoms[1].bond_type[0], atoms[1].bond_stereo[0]), (1, 0));
        assert_eq!((atoms[2].valence, atoms[2].neighbor[0]), (1, 0));
        assert_eq!((atoms[3].valence, atoms[3].neighbor[0]), (1, 1));
        assert_eq!(
            heap.slice(crossing_bonds.as_const()).unwrap(),
            &[2, 4, 9, 10]
        );

        let mut stride_original = original(&mut heap, 4, &[]);
        let stride_edits = edits(
            &mut heap,
            (Vec::new(), 0),
            (Vec::new(), 0),
            (vec![1, 2, 3, 4, 4, 77], 6),
            (Vec::new(), 0),
            (Vec::new(), 0),
        );
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut stride_original,
                &stride_edits,
                &mut status,
            ),
            Ok(2)
        );
        assert_eq!(stride_original.num_inp_bonds, 2);
        let atoms = heap.slice(stride_original.at.as_const()).unwrap();
        assert_eq!((atoms[0].neighbor[0], atoms[1].neighbor[0]), (1, 0));
        assert_eq!((atoms[2].neighbor[0], atoms[3].neighbor[0]), (3, 2));

        let mut skip_original = original(&mut heap, 2, &[(0, 1, 1, 0)]);
        let skip_edits = edits(
            &mut heap,
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (vec![1, 2, 1, 2, 1, 2, 2, 1], 8),
            (Vec::new(), 0),
        );
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut skip_original,
                &skip_edits,
                &mut status,
            ),
            Ok(0)
        );
        assert_eq!(skip_original.num_inp_bonds, 1);

        let mut null_polymer_original = original(&mut heap, 3, &[(0, 1, 1, 0)]);
        let null_polymer_edits = edits(
            &mut heap,
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (vec![1, 2, 1, 3], 4),
            (Vec::new(), 0),
        );
        status = 91;
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut null_polymer_original,
                &null_polymer_edits,
                &mut status,
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(status, _IS_OKAY as i32);
        assert_eq!(null_polymer_original.num_inp_bonds, 1);
        let atoms = heap.slice(null_polymer_original.at.as_const()).unwrap();
        assert_eq!((atoms[0].valence, atoms[0].neighbor[0]), (1, 2));
        assert_eq!(atoms[1].valence, 0);
        assert_eq!((atoms[2].valence, atoms[2].neighbor[0]), (1, 0));

        let mut failed_modify_original = original(&mut heap, 3, &[(0, 1, 1, 0)]);
        heap.slice_mut(failed_modify_original.at).unwrap()[2].valence = 20;
        let failed_modify_edits = edits(
            &mut heap,
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (vec![1, 2, 1, 3], 4),
            (Vec::new(), 0),
        );
        status = 91;
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut failed_modify_original,
                &failed_modify_edits,
                &mut status,
            ),
            Ok(0)
        );
        assert_eq!(status, _IS_ERROR as i32);
        assert_eq!(failed_modify_original.num_inp_bonds, 0);
        let atoms = heap.slice(failed_modify_original.at.as_const()).unwrap();
        assert_eq!((atoms[0].valence, atoms[1].valence), (0, 0));
        assert_eq!(atoms[2].valence, 20);

        let mut coordinate_original = original(&mut heap, 2, &[]);
        let coordinate_edits = edits(
            &mut heap,
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (vec![91, 2], 2),
            (vec![1, 99], 2),
        );
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut coordinate_original,
                &coordinate_edits,
                &mut status,
            ),
            Ok(1)
        );
        let atoms = heap.slice(coordinate_original.at.as_const()).unwrap();
        assert_eq!([atoms[0].x, atoms[0].y, atoms[0].z], [1.25, 1.5, 1.75]);
        assert_eq!([atoms[1].x, atoms[1].y, atoms[1].z], [0.25, 0.5, 0.75]);

        let mut deletion_original = original(&mut heap, 3, &[(0, 1, 1, 0), (0, 2, 2, 0)]);
        let old_atoms = deletion_original.at;
        let deletion_edits = edits(
            &mut heap,
            (vec![3], 1),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
        );
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut deletion_original,
                &deletion_edits,
                &mut status,
            ),
            Ok(0)
        );
        assert_eq!(status, _IS_OKAY as i32);
        assert_eq!(
            heap.slice(old_atoms.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            (
                deletion_original.num_inp_atoms,
                deletion_original.num_inp_bonds
            ),
            (2, 1)
        );
        let atoms = heap.slice(deletion_original.at.as_const()).unwrap();
        assert_eq!((atoms[0].orig_at_number, atoms[0].valence), (1, 1));
        assert_eq!((atoms[0].neighbor[0], atoms[0].neighbor[1]), (1, 2));
        assert_eq!((atoms[0].bond_type[0], atoms[0].bond_type[1]), (1, 2));
        assert_eq!((atoms[1].orig_at_number, atoms[1].neighbor[0]), (2, 0));

        let mut polymer_original = original(&mut heap, 4, &[]);
        let v3000 = heap
            .allocate_model_storage(vec![crate::source_types::OAD_V3000::default()])
            .unwrap();
        polymer_original.v3000 = v3000;
        let alist = heap.allocate_model_storage(vec![1_i32, 2, 3]).unwrap();
        let blist = heap.allocate_model_storage(vec![1_i32, 2, 2, 3]).unwrap();
        let first_backbone = heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let second_backbone = heap.allocate_model_storage(vec![3_i32, 4]).unwrap();
        let backbone = heap
            .allocate_model_storage(vec![first_backbone, second_backbone])
            .unwrap();
        let polymer_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 3,
                nb: 2,
                alist,
                blist,
                maxbkbonds: 2,
                nbkbonds: 2,
                bkbonds: backbone,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let polymer_units = heap.allocate_model_storage(vec![polymer_unit]).unwrap();
        polymer_original.polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units: polymer_units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let polymer_edits = edits(
            &mut heap,
            (vec![4], 1),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
        );
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut polymer_original,
                &polymer_edits,
                &mut status,
            ),
            Ok(0)
        );
        assert_eq!(polymer_original.num_inp_atoms, 3);
        assert_eq!(polymer_original.v3000, v3000);
        let unit_state = &heap.slice(polymer_unit.as_const()).unwrap()[0];
        assert_eq!(
            (unit_state.na, unit_state.nb, unit_state.nbkbonds),
            (3, 2, 1)
        );
        assert_eq!(heap.slice(alist.as_const()).unwrap(), &[1, 2, 3]);
        assert_eq!(heap.slice(blist.as_const()).unwrap(), &[1, 2, 2, 3]);
        assert_eq!(heap.slice(first_backbone.as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(second_backbone.as_const()).unwrap(), &[0, 0]);
        assert_eq!(
            (
                unit_state.cap1,
                unit_state.end_atom1,
                unit_state.cap2,
                unit_state.end_atom2,
            ),
            (1, 2, 2, 3)
        );

        let mut renumbered_neighbor_original = original(&mut heap, 3, &[(0, 2, 1, 0)]);
        let renumbered_neighbor_edits = edits(
            &mut heap,
            (vec![2], 1),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
        );
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut renumbered_neighbor_original,
                &renumbered_neighbor_edits,
                &mut status,
            ),
            Ok(0)
        );
        let atoms = heap
            .slice(renumbered_neighbor_original.at.as_const())
            .unwrap();
        assert_eq!(renumbered_neighbor_original.num_inp_atoms, 2);
        assert_eq!((atoms[0].neighbor[0], atoms[1].neighbor[0]), (1, 0));

        for failure_ordinal in 0..=2 {
            let mut failure_heap = SourceHeap::default();
            let mut failure_original = original(&mut failure_heap, 1, &[]);
            let original_pointer = failure_original.at;
            let failure_edits = edits(
                &mut failure_heap,
                (vec![1], 1),
                (Vec::new(), 0),
                (Vec::new(), 0),
                (Vec::new(), 0),
                (Vec::new(), 0),
            );
            let baseline = failure_heap.live_allocation_count();
            failure_heap.fail_after_allocations(failure_ordinal);
            status = 91;
            assert_eq!(
                OAD_StructureEdits_Apply(
                    &mut failure_heap,
                    &mut STRUCT_DATA::default(),
                    &INPUT_PARMS::default(),
                    &mut failure_original,
                    &failure_edits,
                    &mut status,
                ),
                Ok(0)
            );
            assert_eq!(status, _IS_ERROR as i32);
            assert_eq!(failure_heap.live_allocation_count(), baseline);
            assert_eq!(failure_original.at, original_pointer);
            assert_eq!(failure_original.num_inp_atoms, 1);
        }

        let mut missing_bond_original = original(&mut heap, 2, &[]);
        let missing_bond_edits = edits(
            &mut heap,
            (Vec::new(), 0),
            (vec![1, 2], 2),
            (Vec::new(), 0),
            (Vec::new(), 0),
            (Vec::new(), 0),
        );
        status = 91;
        assert_eq!(
            OAD_StructureEdits_Apply(
                &mut heap,
                &mut STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut missing_bond_original,
                &missing_bond_edits,
                &mut status,
            ),
            Ok(0)
        );
        assert_eq!(status, _IS_ERROR as i32);
    }

    fn prepared_structure(
        heap: &mut SourceHeap,
        atoms: Vec<inp_ATOM>,
        component_lengths: Vec<u16>,
        old_component_numbers: Vec<u16>,
    ) -> ORIG_ATOM_DATA {
        ORIG_ATOM_DATA {
            num_inp_atoms: atoms.len() as i32,
            at: heap.allocate_model_storage(atoms).unwrap(),
            num_components: component_lengths.len() as i32,
            nCurAtLen: heap.allocate_model_storage(component_lengths).unwrap(),
            nOldCompNumber: heap.allocate_model_storage(old_component_numbers).unwrap(),
            ..ORIG_ATOM_DATA::default()
        }
    }

    #[test]
    fn source_port__runichi__processonestructure__line_218() {
        assert_eq!(
            ProcessOneStructure(
                &mut SourceHeap::default(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut STRUCT_DATA {
                    bUserQuitComponent: 7,
                    ..STRUCT_DATA::default()
                },
                &mut INPUT_PARMS::default(),
                None,
                &mut [SourceMutPointer::null(); INCHI_NUM as usize],
                &mut [SourceMutPointer::null(); INCHI_NUM as usize],
                None,
                None,
                None,
                None,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                i64::MIN,
                None,
                0xff,
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );

        let mut skip_heap = SourceHeap::default();
        let skip_clock = skip_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let skip_canonical = skip_heap
            .allocate_model_storage(vec![CANON_GLOBALS::default()])
            .unwrap();
        let skip_original = skip_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA::default()])
            .unwrap();
        let skip_prepared = skip_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA::default(); INCHI_NUM as usize])
            .unwrap();
        let mut skip_structure = STRUCT_DATA {
            bUserQuitComponent: 11,
            bUserQuitComponentDisplay: 13,
            ..STRUCT_DATA::default()
        };
        assert_eq!(
            ProcessOneStructure(
                &mut skip_heap,
                skip_clock,
                skip_canonical,
                &mut skip_structure,
                &mut INPUT_PARMS {
                    bIgnoreUnchanged: 1,
                    ..INPUT_PARMS::default()
                },
                None,
                &mut [SourceMutPointer::null(); INCHI_NUM as usize],
                &mut [SourceMutPointer::null(); INCHI_NUM as usize],
                None,
                None,
                None,
                None,
                skip_original,
                skip_prepared,
                -9,
                None,
                0xa5,
                SourceMutPointer::null(),
                0,
            ),
            Ok(crate::source_types::_IS_SKIP as i32)
        );
        assert_eq!(skip_structure.bUserQuitComponent, 0);
        assert_eq!(skip_structure.bUserQuitComponentDisplay, 0);

        let mut full_heap = SourceHeap::default();
        let full_clock = full_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let full_canonical = full_heap
            .allocate_model_storage(vec![CANON_GLOBALS::default()])
            .unwrap();
        let full_original = full_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA::default()])
            .unwrap();
        let full_prepared = full_heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA::default(); INCHI_NUM as usize])
            .unwrap();
        let stdout = full_heap
            .allocate_model_storage(vec![SourceFile::default()])
            .unwrap();
        let mut output = string_stream();
        let mut log = string_stream();
        let mut title = [0_i8; 256];
        let mut inchi = [SourceMutPointer::null(); INCHI_NUM as usize];
        let mut aux = [SourceMutPointer::null(); INCHI_NUM as usize];
        let mut structure = STRUCT_DATA {
            bChiralFlag: FLAG_INP_AT_CHIRAL as i32,
            bUserQuitComponent: 17,
            bUserQuitComponentDisplay: 19,
            ..STRUCT_DATA::default()
        };
        let mut parameters = INPUT_PARMS {
            bAllowEmptyStructure: 1,
            nMode: u64::from(REQ_MODE_STEREO),
            bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            ProcessOneStructure(
                &mut full_heap,
                full_clock,
                full_canonical,
                &mut structure,
                &mut parameters,
                Some(&mut title),
                &mut inchi,
                &mut aux,
                None,
                Some(&mut log),
                Some(&mut output),
                None,
                full_original,
                full_prepared,
                5,
                None,
                0,
                stdout,
                0,
            ),
            Ok(crate::source_types::_IS_WARNING as i32)
        );
        assert_eq!(structure.bUserQuitComponent, 0);
        assert_eq!(structure.bUserQuitComponentDisplay, 0);
        assert_eq!(array_c_bytes(&structure.pStrErrStruct), b"Not chiral");
        assert_eq!(
            stream_bytes(&full_heap, &log),
            b"Warning (Not chiral) structure #5.\n"
        );
        let full_output = String::from_utf8(stream_bytes(&full_heap, &output)).unwrap();
        assert_eq!(full_output, "");
        assert!(!inchi[INCHI_BAS as usize].is_null());
        assert!(!aux[INCHI_BAS as usize].is_null());
    }

    #[test]
    fn source_port__runichi__createonestructureinchi__line_809() {
        let mut disabled_heap = SourceHeap::default();
        let disabled_clock = disabled_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut disabled_outputs = [SourceMutPointer::null(); INCHI_NUM as usize];
        let mut disabled_aux = [SourceMutPointer::null(); INCHI_NUM as usize];
        let mut disabled_parameters = INPUT_PARMS {
            msec_MaxTime: 17,
            msec_LeftTime: -1,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            CreateOneStructureINChI(
                &mut disabled_heap,
                &mut CANON_GLOBALS::default(),
                disabled_clock,
                &mut STRUCT_DATA::default(),
                &mut disabled_parameters,
                None,
                &mut disabled_outputs,
                &mut disabled_aux,
                INCHI_BAS as i32,
                None,
                None,
                None,
                None,
                &mut ORIG_ATOM_DATA::default(),
                &mut [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()],
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                1,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(disabled_outputs, [SourceMutPointer::null(); 2]);
        assert_eq!(disabled_parameters.msec_LeftTime, 17);

        let mut empty_heap = SourceHeap::default();
        let empty_clock = empty_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut empty_outputs = [SourceMutPointer::null(); INCHI_NUM as usize];
        let mut empty_aux = [SourceMutPointer::null(); INCHI_NUM as usize];
        let mut empty_structure = STRUCT_DATA::default();
        let mut empty_parameters = INPUT_PARMS {
            bAllowEmptyStructure: 1,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            CreateOneStructureINChI(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                empty_clock,
                &mut empty_structure,
                &mut empty_parameters,
                None,
                &mut empty_outputs,
                &mut empty_aux,
                INCHI_BAS as i32,
                None,
                None,
                None,
                None,
                &mut ORIG_ATOM_DATA::default(),
                &mut [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()],
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                2,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(0)
        );
        assert!(!empty_outputs[INCHI_BAS as usize].is_null());
        assert!(!empty_aux[INCHI_BAS as usize].is_null());
        assert_eq!(empty_structure.num_components[INCHI_BAS as usize], 0);

        let mut invalid_heap = SourceHeap::default();
        let invalid_clock = invalid_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let invalid_atoms = invalid_heap
            .allocate_model_storage(vec![structure_carbon(0, 1)])
            .unwrap();
        let mut invalid_original = ORIG_ATOM_DATA {
            at: invalid_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut invalid_structure = STRUCT_DATA::default();
        assert_eq!(
            CreateOneStructureINChI(
                &mut invalid_heap,
                &mut CANON_GLOBALS::default(),
                invalid_clock,
                &mut invalid_structure,
                &mut INPUT_PARMS::default(),
                None,
                &mut [SourceMutPointer::null(); 2],
                &mut [SourceMutPointer::null(); 2],
                -1,
                None,
                None,
                None,
                None,
                &mut invalid_original,
                &mut [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()],
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                3,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(invalid_structure.nStructReadError, 97);

        for successful_allocations in [0_u64, 1] {
            let mut failure_heap = SourceHeap::default();
            let failure_clock = failure_heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            failure_heap.fail_after_allocations(successful_allocations);
            let mut failure_outputs = [SourceMutPointer::null(); 2];
            let mut failure_aux = [SourceMutPointer::null(); 2];
            let mut failure_structure = STRUCT_DATA::default();
            assert_eq!(
                CreateOneStructureINChI(
                    &mut failure_heap,
                    &mut CANON_GLOBALS::default(),
                    failure_clock,
                    &mut failure_structure,
                    &mut INPUT_PARMS {
                        bAllowEmptyStructure: 1,
                        ..INPUT_PARMS::default()
                    },
                    None,
                    &mut failure_outputs,
                    &mut failure_aux,
                    INCHI_BAS as i32,
                    None,
                    None,
                    None,
                    None,
                    &mut ORIG_ATOM_DATA::default(),
                    &mut [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()],
                    &mut std::array::from_fn(|_| {
                        std::array::from_fn(|_| COMP_ATOM_DATA::default())
                    }),
                    4,
                    None,
                    &mut NORM_CANON_FLAGS::default(),
                    0,
                ),
                Ok(0),
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(failure_structure.nStructReadError, 99);
            assert_eq!(failure_structure.nErrorType, _IS_FATAL as i32);
            assert_eq!(failure_outputs, [SourceMutPointer::null(); 2]);
            assert_eq!(failure_aux, [SourceMutPointer::null(); 2]);
        }

        let mut single_heap = SourceHeap::default();
        let single_clock = single_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let single_atoms = single_heap
            .allocate_model_storage(vec![structure_carbon(0, 1)])
            .unwrap();
        let mut single_original = ORIG_ATOM_DATA {
            at: single_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut single_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut single_outputs = [SourceMutPointer::null(); 2];
        let mut single_aux = [SourceMutPointer::null(); 2];
        let mut single_structure = STRUCT_DATA::default();
        let mut single_flags = NORM_CANON_FLAGS::default();
        let mut single_composite =
            std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
        assert_eq!(
            CreateOneStructureINChI(
                &mut single_heap,
                &mut CANON_GLOBALS::default(),
                single_clock,
                &mut single_structure,
                &mut INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
                    bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
                    ..INPUT_PARMS::default()
                },
                None,
                &mut single_outputs,
                &mut single_aux,
                INCHI_BAS as i32,
                None,
                None,
                None,
                None,
                &mut single_original,
                &mut single_prepared,
                &mut single_composite,
                5,
                None,
                &mut single_flags,
                0,
            ),
            Ok(0)
        );
        assert_eq!(single_prepared[0].num_components, 1);
        assert_eq!(single_structure.num_components[0], 1);
        let single_rows = single_heap
            .slice(single_outputs[INCHI_BAS as usize].as_const())
            .unwrap();
        assert!(!single_rows[0][TAUT_YES as usize].is_null());
        assert_eq!(single_structure.num_non_taut[0], 1);

        let mut reuse_heap = SourceHeap::default();
        let reuse_clock = reuse_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut reuse_original = ORIG_ATOM_DATA {
            at: reuse_heap
                .allocate_model_storage(vec![structure_carbon(1, 1)])
                .unwrap(),
            num_inp_atoms: 1,
            bDisconnectCoord: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut reuse_prepared = [
            prepared_structure(
                &mut reuse_heap,
                vec![structure_carbon(1, 1)],
                vec![1],
                vec![1],
            ),
            prepared_structure(
                &mut reuse_heap,
                vec![structure_carbon(1, 1)],
                vec![1],
                vec![1],
            ),
        ];
        let mut reuse_outputs = [SourceMutPointer::null(); 2];
        let mut reuse_aux = [SourceMutPointer::null(); 2];
        let mut reuse_structure = STRUCT_DATA::default();
        let mut reuse_parameters = INPUT_PARMS {
            nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
            bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
            bTautFlags: TG_FLAG_RECONNECT_COORD as INCHI_MODE,
            ..INPUT_PARMS::default()
        };
        let mut reuse_composite =
            std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
        let mut reuse_flags = NORM_CANON_FLAGS::default();
        for kind in [INCHI_BAS as i32, INCHI_REC as i32] {
            assert_eq!(
                CreateOneStructureINChI(
                    &mut reuse_heap,
                    &mut CANON_GLOBALS::default(),
                    reuse_clock,
                    &mut reuse_structure,
                    &mut reuse_parameters,
                    None,
                    &mut reuse_outputs,
                    &mut reuse_aux,
                    kind,
                    None,
                    None,
                    None,
                    None,
                    &mut reuse_original,
                    &mut reuse_prepared,
                    &mut reuse_composite,
                    6,
                    None,
                    &mut reuse_flags,
                    0,
                ),
                Ok(0)
            );
        }
        let base_row = reuse_heap.slice(reuse_outputs[0].as_const()).unwrap()[0];
        let reconnected_row = reuse_heap.slice(reuse_outputs[1].as_const()).unwrap()[0];
        assert_eq!(base_row, reconnected_row);
        let shared = base_row[TAUT_YES as usize];
        assert_eq!(reuse_heap.slice(shared.as_const()).unwrap()[0].nRefCount, 1);
        assert_eq!(reuse_structure.num_non_taut, [1, 1]);

        let mut preprocess_failure_heap = SourceHeap::default();
        let preprocess_failure_clock = preprocess_failure_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut preprocess_failure_original = ORIG_ATOM_DATA {
            at: preprocess_failure_heap
                .allocate_model_storage(vec![structure_carbon(0, 1)])
                .unwrap(),
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut preprocess_failure_prepared =
            [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut preprocess_failure_structure = STRUCT_DATA::default();
        preprocess_failure_heap.fail_after_allocations(0);
        assert_eq!(
            CreateOneStructureINChI(
                &mut preprocess_failure_heap,
                &mut CANON_GLOBALS::default(),
                preprocess_failure_clock,
                &mut preprocess_failure_structure,
                &mut INPUT_PARMS::default(),
                None,
                &mut [SourceMutPointer::null(); 2],
                &mut [SourceMutPointer::null(); 2],
                INCHI_BAS as i32,
                None,
                None,
                None,
                None,
                &mut preprocess_failure_original,
                &mut preprocess_failure_prepared,
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                7,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(preprocess_failure_structure.nStructReadError, 99);
        assert_eq!(
            array_c_bytes(&preprocess_failure_structure.pStrErrStruct),
            b"Out of RAM"
        );

        for fail_composite_storage in [false, true] {
            let mut multi_heap = SourceHeap::default();
            let multi_clock = multi_heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let mut multi_original = ORIG_ATOM_DATA {
                at: multi_heap
                    .allocate_model_storage(vec![structure_carbon(1, 1), structure_carbon(2, 2)])
                    .unwrap(),
                num_inp_atoms: 2,
                ..ORIG_ATOM_DATA::default()
            };
            let mut multi_prepared = [
                prepared_structure(
                    &mut multi_heap,
                    vec![structure_carbon(1, 1), structure_carbon(2, 2)],
                    vec![1, 1],
                    vec![1, 2],
                ),
                ORIG_ATOM_DATA::default(),
            ];
            let mut multi_outputs = [SourceMutPointer::null(); 2];
            let mut multi_aux = [SourceMutPointer::null(); 2];
            let mut multi_structure = STRUCT_DATA::default();
            let mut multi_composite =
                std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
            if fail_composite_storage {
                multi_heap.fail_after_allocations(0);
            }
            assert_eq!(
                CreateOneStructureINChI(
                    &mut multi_heap,
                    &mut CANON_GLOBALS::default(),
                    multi_clock,
                    &mut multi_structure,
                    &mut INPUT_PARMS {
                        nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
                        bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
                        ..INPUT_PARMS::default()
                    },
                    None,
                    &mut multi_outputs,
                    &mut multi_aux,
                    INCHI_BAS as i32,
                    None,
                    None,
                    None,
                    None,
                    &mut multi_original,
                    &mut multi_prepared,
                    &mut multi_composite,
                    8,
                    None,
                    &mut NORM_CANON_FLAGS::default(),
                    0,
                ),
                Ok(0),
                "all normalized allocation failure={fail_composite_storage}"
            );
            assert_eq!(multi_structure.num_components[0], 2);
            let rows = multi_heap.slice(multi_outputs[0].as_const()).unwrap();
            assert!(!rows[0][TAUT_YES as usize].is_null());
            assert!(!rows[1][TAUT_YES as usize].is_null());
            assert_eq!(
                multi_composite[0][TAUT_YES as usize].bExists,
                i32::from(!fail_composite_storage)
            );
            if !fail_composite_storage {
                assert_eq!(multi_composite[0][TAUT_YES as usize].num_at, 2);
            }
        }

        let mut title_heap = SourceHeap::default();
        let title_clock = title_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let title_atom = title_heap
            .allocate_model_storage(vec![structure_carbon(1, 1)])
            .unwrap();
        let mut title_original = ORIG_ATOM_DATA {
            at: title_atom,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut title_prepared = [
            prepared_structure(
                &mut title_heap,
                vec![structure_carbon(1, 1)],
                vec![1],
                vec![1],
            ),
            ORIG_ATOM_DATA::default(),
        ];
        let title_label = c_text(&mut title_heap, b"ID");
        let title_value = c_text(&mut title_heap, b"abc");
        let mut title = [0x5a_i8; 256];
        assert_eq!(
            CreateOneStructureINChI(
                &mut title_heap,
                &mut CANON_GLOBALS::default(),
                title_clock,
                &mut STRUCT_DATA::default(),
                &mut INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
                    bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
                    bDisplay: 1,
                    pSdfLabel: title_label,
                    pSdfValue: title_value,
                    ..INPUT_PARMS::default()
                },
                Some(&mut title),
                &mut [SourceMutPointer::null(); 2],
                &mut [SourceMutPointer::null(); 2],
                INCHI_BAS as i32,
                None,
                None,
                None,
                None,
                &mut title_original,
                &mut title_prepared,
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                9,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(0)
        );
        let expected_title = b"Result for Structure #9. ID=abc";
        assert_eq!(&array_c_bytes(&title), expected_title);
        assert_eq!(title[expected_title.len()], 0);
        assert_eq!(title[expected_title.len() + 1], 0x5a);

        let mut quit_heap = SourceHeap::default();
        let quit_clock = quit_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut quit_original = ORIG_ATOM_DATA {
            at: quit_heap
                .allocate_model_storage(vec![structure_carbon(1, 1)])
                .unwrap(),
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut quit_prepared = [
            prepared_structure(
                &mut quit_heap,
                vec![structure_carbon(1, 1)],
                vec![1],
                vec![1],
            ),
            ORIG_ATOM_DATA::default(),
        ];
        let mut quit_outputs = [SourceMutPointer::null(); 2];
        let mut quit_aux = [SourceMutPointer::null(); 2];
        let mut quit_structure = STRUCT_DATA {
            bUserQuitComponent: 1,
            ..STRUCT_DATA::default()
        };
        assert_eq!(
            CreateOneStructureINChI(
                &mut quit_heap,
                &mut CANON_GLOBALS::default(),
                quit_clock,
                &mut quit_structure,
                &mut INPUT_PARMS::default(),
                None,
                &mut quit_outputs,
                &mut quit_aux,
                INCHI_BAS as i32,
                None,
                None,
                None,
                None,
                &mut quit_original,
                &mut quit_prepared,
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                10,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(0)
        );
        assert!(
            quit_heap.slice(quit_outputs[0].as_const()).unwrap()[0]
                .iter()
                .all(|pointer| pointer.is_null())
        );

        let mut capacity_heap = SourceHeap::default();
        let capacity_clock = capacity_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let retained = capacity_heap
            .allocate_model_storage(vec![INChI::default()])
            .unwrap();
        let capacity_rows = capacity_heap
            .allocate_model_storage(vec![
                [SourceMutPointer::null(); 2],
                [SourceMutPointer::null(); 2],
                [retained, SourceMutPointer::null()],
            ])
            .unwrap();
        let capacity_aux_rows = capacity_heap
            .allocate_model_storage(vec![[SourceMutPointer::null(); 2]; 3])
            .unwrap();
        let mut capacity_outputs = [capacity_rows, SourceMutPointer::null()];
        let mut capacity_aux = [capacity_aux_rows, SourceMutPointer::null()];
        let mut capacity_original = ORIG_ATOM_DATA {
            at: capacity_heap
                .allocate_model_storage(vec![structure_carbon(1, 1)])
                .unwrap(),
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut capacity_prepared = [
            prepared_structure(
                &mut capacity_heap,
                vec![structure_carbon(1, 1)],
                vec![1],
                vec![1],
            ),
            ORIG_ATOM_DATA::default(),
        ];
        let mut capacity_structure = STRUCT_DATA::default();
        capacity_structure.num_components[0] = 2;
        assert_eq!(
            CreateOneStructureINChI(
                &mut capacity_heap,
                &mut CANON_GLOBALS::default(),
                capacity_clock,
                &mut capacity_structure,
                &mut INPUT_PARMS {
                    nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
                    bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
                    ..INPUT_PARMS::default()
                },
                None,
                &mut capacity_outputs,
                &mut capacity_aux,
                INCHI_BAS as i32,
                None,
                None,
                None,
                None,
                &mut capacity_original,
                &mut capacity_prepared,
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                12,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(capacity_outputs[0], capacity_rows);
        assert_eq!(capacity_aux[0], capacity_aux_rows);
        assert_eq!(capacity_structure.num_components[0], 2);
        assert_eq!(
            capacity_heap.slice(capacity_rows.as_const()).unwrap()[2][0],
            retained
        );

        let mut duplicate_heap = SourceHeap::default();
        let duplicate_clock = duplicate_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let first_inchi = duplicate_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let second_inchi = duplicate_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let base_rows = duplicate_heap
            .allocate_model_storage(vec![
                [SourceMutPointer::null(), first_inchi],
                [SourceMutPointer::null(), second_inchi],
            ])
            .unwrap();
        let base_aux = duplicate_heap
            .allocate_model_storage(vec![[SourceMutPointer::null(); 2]; 2])
            .unwrap();
        let mut duplicate_outputs = [base_rows, SourceMutPointer::null()];
        let mut duplicate_aux = [base_aux, SourceMutPointer::null()];
        let mut duplicate_original = ORIG_ATOM_DATA {
            at: duplicate_heap
                .allocate_model_storage(vec![structure_carbon(1, 1)])
                .unwrap(),
            num_inp_atoms: 1,
            bDisconnectCoord: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut duplicate_prepared = [
            prepared_structure(
                &mut duplicate_heap,
                vec![structure_carbon(1, 1), structure_carbon(2, 2)],
                vec![1, 1],
                vec![1, 1],
            ),
            prepared_structure(
                &mut duplicate_heap,
                vec![structure_carbon(1, 1)],
                vec![1],
                vec![1],
            ),
        ];
        let mut duplicate_structure = STRUCT_DATA::default();
        duplicate_structure.num_components[0] = 2;
        assert_eq!(
            CreateOneStructureINChI(
                &mut duplicate_heap,
                &mut CANON_GLOBALS::default(),
                duplicate_clock,
                &mut duplicate_structure,
                &mut INPUT_PARMS::default(),
                None,
                &mut duplicate_outputs,
                &mut duplicate_aux,
                INCHI_REC as i32,
                None,
                None,
                None,
                None,
                &mut duplicate_original,
                &mut duplicate_prepared,
                &mut std::array::from_fn(|_| {
                    std::array::from_fn(|_| COMP_ATOM_DATA::default())
                }),
                11,
                None,
                &mut NORM_CANON_FLAGS::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(duplicate_structure.nErrorType, _IS_ERROR as i32);
        assert_eq!(duplicate_structure.nStructReadError, 99);
        assert_eq!(
            array_c_bytes(&duplicate_structure.pStrErrStruct),
            b"Cannot distinguish components"
        );
        assert_eq!(
            duplicate_heap.slice(first_inchi.as_const()).unwrap()[0].nRefCount,
            1
        );
    }

    #[test]
    fn source_port__runichi__createonecomponentinchi__line_1747() {
        let mut empty_heap = SourceHeap::default();
        let empty_clock = empty_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut empty_structure = STRUCT_DATA::default();
        let mut empty_parameters = INPUT_PARMS::default();
        let mut empty_input = INP_ATOM_DATA::default();
        let mut empty_inchi = vec![[SourceMutPointer::null(); TAUT_NUM as usize]];
        let mut empty_aux = vec![[SourceMutPointer::null(); TAUT_NUM as usize]];
        let mut empty_normalized = std::array::from_fn(|_| INP_ATOM_DATA::default());
        assert_eq!(
            CreateOneComponentINChI(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                empty_clock,
                &mut empty_structure,
                &mut empty_parameters,
                &mut empty_input,
                None,
                &mut empty_inchi,
                &mut empty_aux,
                INCHI_BAS as i32,
                0,
                -1,
                &mut empty_normalized,
                &mut NORM_CANON_FLAGS::default(),
                &mut INCHI_IOSTREAM::default(),
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(empty_structure.nErrorCode, -1);
        assert_eq!(
            empty_inchi[0],
            [SourceMutPointer::null(); TAUT_NUM as usize]
        );

        let mut heap = SourceHeap::default();
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut carbon = inp_ATOM {
            el_number: 6,
            num_H: 4,
            orig_at_number: 1,
            component: 99,
            ..inp_ATOM::default()
        };
        carbon.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
        let atoms = heap.allocate_model_storage(vec![carbon]).unwrap();
        let mut input = INP_ATOM_DATA {
            at: atoms,
            num_at: 1,
            ..INP_ATOM_DATA::default()
        };
        let mut parameters = INPUT_PARMS {
            nMode: u64::from(REQ_MODE_TAUT | crate::source_types::REQ_MODE_NON_ISO),
            msec_MaxTime: 10,
            msec_LeftTime: 10,
            bINChIOutputOptions: INCHI_OUT_NO_AUX_INFO as i32,
            ..INPUT_PARMS::default()
        };
        let mut structure = STRUCT_DATA::default();
        let mut inchi = vec![[SourceMutPointer::null(); TAUT_NUM as usize]];
        let mut aux = vec![[SourceMutPointer::null(); TAUT_NUM as usize]];
        let mut normalized = std::array::from_fn(|_| INP_ATOM_DATA::default());
        let mut flags = NORM_CANON_FLAGS::default();
        assert_eq!(
            CreateOneComponentINChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &mut structure,
                &mut parameters,
                &mut input,
                None,
                &mut inchi,
                &mut aux,
                INCHI_BAS as i32,
                0,
                i64::MAX,
                &mut normalized,
                &mut flags,
                &mut INCHI_IOSTREAM::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(structure.nErrorCode, 0);
        assert_eq!(structure.num_non_taut[INCHI_BAS as usize], 1);
        assert_eq!(structure.num_taut[INCHI_BAS as usize], 0);
        assert_eq!(parameters.msec_LeftTime, 10);
        assert_eq!(heap.slice(atoms.as_const()).unwrap()[0].component, 1);
        assert!(inchi[0][TAUT_NON as usize].is_null());
        assert!(!inchi[0][TAUT_YES as usize].is_null());
        assert!(!aux[0][TAUT_YES as usize].is_null());
        let generated = &heap.slice(inchi[0][TAUT_YES as usize].as_const()).unwrap()[0];
        assert_eq!(generated.nNumberOfAtoms, 1);
        assert_eq!(heap.slice(generated.nAtom.as_const()).unwrap(), &[6]);
        assert_eq!(heap.slice(generated.nNum_H.as_const()).unwrap(), &[4]);
        assert_eq!(normalized[TAUT_YES as usize].bExists, 1);
        assert_eq!(normalized[TAUT_YES as usize].bHasIsotopicLayer, 0);
    }

    #[test]
    fn source_port__runichi__origatdata_savemolfile__line_564() {
        let mut disabled_heap = SourceHeap::default();
        let disabled = INPUT_PARMS {
            pSdfLabel: SourceMutPointer::null(),
            pSdfValue: SourceMutPointer::null(),
            ..INPUT_PARMS::default()
        };
        let mut disabled_stream = string_stream();
        assert_eq!(
            OrigAtData_SaveMolfile(
                &mut disabled_heap,
                &ORIG_ATOM_DATA::default(),
                &STRUCT_DATA::default(),
                &disabled,
                i64::MIN,
                &mut disabled_stream,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(stream_bytes(&disabled_heap, &disabled_stream), b"");

        for (label_text, value_text, expected_title) in [
            (b"".as_slice(), b"".as_slice(), "Structure #-7. "),
            (b"".as_slice(), b"V".as_slice(), "Structure #-7. V"),
            (
                b"L".as_slice(),
                b"".as_slice(),
                "Structure #-7.  L is missing",
            ),
            (b"L".as_slice(), b"V".as_slice(), "Structure #-7.  L=V"),
        ] {
            let mut heap = SourceHeap::default();
            let label = c_text(&mut heap, label_text);
            let value = c_text(&mut heap, value_text);
            let parameters = INPUT_PARMS {
                pSdfLabel: label,
                pSdfValue: value,
                bINChIOutputOptions: INCHI_OUT_SDFILE_ONLY as i32,
                ..INPUT_PARMS::default()
            };
            let mut stream = string_stream();
            assert_eq!(
                OrigAtData_SaveMolfile(
                    &mut heap,
                    &ORIG_ATOM_DATA::default(),
                    &STRUCT_DATA::default(),
                    &parameters,
                    -7,
                    &mut stream,
                ),
                Ok(0)
            );
            let output = String::from_utf8(stream_bytes(&heap, &stream)).unwrap();
            assert_eq!(output.lines().next(), Some(expected_title));
            assert!(output.ends_with("$$$$\n"));
            if value_text.is_empty() {
                assert!(!output.contains("> <"));
            } else {
                let field = if label_text.is_empty() { "ID" } else { "L" };
                assert!(output.contains(&format!("> <{field}>\n V\n\n")));
            }
        }

        let mut heap = SourceHeap::default();
        let label = c_text(&mut heap, b"ISO");
        let value = c_text(&mut heap, b"D");
        let mut atom = inp_ATOM::default();
        atom.elname[..2].copy_from_slice(&[b'H' as i8, 0]);
        atom.el_number = 1;
        atom.iso_atw_diff = 2;
        atom.valence = 1;
        atom.chem_bonds_valence = 1;
        atom.bond_type[0] = 1;
        let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
        let input = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let parameters = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            bINChIOutputOptions: (INCHI_OUT_SDFILE_ONLY | INCHI_OUT_SDFILE_ATOMS_DT) as i32,
            ..INPUT_PARMS::default()
        };
        let structure = STRUCT_DATA {
            bChiralFlag: FLAG_INP_AT_CHIRAL as i32,
            ..STRUCT_DATA::default()
        };
        let mut stream = string_stream();
        assert_eq!(
            OrigAtData_SaveMolfile(&mut heap, &input, &structure, &parameters, 9, &mut stream),
            Ok(0)
        );
        let output = String::from_utf8(stream_bytes(&heap, &stream)).unwrap();
        assert!(output.starts_with("Structure #9.  ISO=D\n"));
        assert!(output.contains("  1  0  0  0  1"));
        assert!(output.contains(" D   0  0"));
    }

    #[test]
    fn source_port__runichi__origatdata_storenativeinput__line_593() {
        fn error_text(structure_data: &STRUCT_DATA) -> String {
            let end = structure_data
                .pStrErrStruct
                .iter()
                .position(|byte| *byte == 0)
                .unwrap();
            String::from_utf8(
                structure_data.pStrErrStruct[..end]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap()
        }

        let parameters = INPUT_PARMS {
            bINChIOutputOptions: i32::MAX,
            ..INPUT_PARMS::default()
        };

        let mut success_heap = SourceHeap::default();
        let mut success_return = 71;
        let mut success_data = STRUCT_DATA::default();
        let mut success_input = ORIG_ATOM_DATA {
            n_zy: 12,
            ..ORIG_ATOM_DATA::default()
        };
        let mut success_output = ORIG_STRUCT {
            num_atoms: 88,
            ..ORIG_STRUCT::default()
        };
        let success_address = std::ptr::addr_of_mut!(success_output);
        let returned = OrigAtData_StoreNativeInput(
            &mut success_heap,
            &mut CANON_GLOBALS::default(),
            &mut success_return,
            &mut success_data,
            &parameters,
            &mut success_input,
            &mut success_output,
        )
        .unwrap() as *mut ORIG_STRUCT;
        assert_eq!(returned, success_address);
        assert_eq!(success_return, 71);
        assert_eq!(success_data.nStructReadError, 0);
        assert_eq!(success_data.nErrorType, 0);
        assert_eq!(error_text(&success_data), "");
        assert_eq!(success_output.num_atoms, 0);
        assert_eq!(success_output.n_zy, 12);
        assert!(!success_output.szAtoms.is_null());
        assert!(!success_output.szBonds.is_null());
        crate::source::base::ichiprt1::OrigStruct_Free(
            &mut success_heap,
            Some(&mut success_output),
        )
        .unwrap();

        let mut failure_heap = SourceHeap::default();
        failure_heap.fail_after_allocations(0);
        let mut failure_return = 73;
        let mut failure_data = STRUCT_DATA::default();
        failure_data.pStrErrStruct[..6].copy_from_slice(&[
            b'p' as i8, b'r' as i8, b'i' as i8, b'o' as i8, b'r' as i8, 0,
        ]);
        let mut failure_input = ORIG_ATOM_DATA {
            n_zy: -7,
            ..ORIG_ATOM_DATA::default()
        };
        let mut failure_output = ORIG_STRUCT {
            num_atoms: 91,
            ..ORIG_STRUCT::default()
        };
        let failure_address = std::ptr::addr_of_mut!(failure_output);
        let returned = OrigAtData_StoreNativeInput(
            &mut failure_heap,
            &mut CANON_GLOBALS::default(),
            &mut failure_return,
            &mut failure_data,
            &parameters,
            &mut failure_input,
            &mut failure_output,
        )
        .unwrap() as *mut ORIG_STRUCT;
        assert_eq!(returned, failure_address);
        assert_eq!(failure_return, _IS_ERROR as i32);
        assert_eq!(failure_data.nStructReadError, 99);
        assert_eq!(failure_data.nErrorType, _IS_ERROR as i32);
        assert_eq!(
            error_text(&failure_data),
            "prior; Cannot interpret reversibility information"
        );
        assert_eq!(failure_output.num_atoms, 91);
        assert_eq!(failure_output.n_zy, -7);
        assert!(failure_output.szAtoms.is_null());
        assert!(failure_output.szBonds.is_null());
    }

    #[test]
    fn source_port__runichi__preparesaveoptbits__line_621() {
        let mut preserved = 0xa5;
        PrepareSaveOptBitsForStructure(
            &mut preserved,
            &INPUT_PARMS {
                nInputType: tagInputType_INPUT_INCHI,
                bINChIOutputOptions: INCHI_OUT_SAVEOPT as i32,
                nMode: u64::MAX,
                bTautFlags: u64::MAX,
                ..INPUT_PARMS::default()
            },
        );
        assert_eq!(preserved, 0xa5);

        let mut cleared = 0xff;
        PrepareSaveOptBitsForStructure(&mut cleared, &INPUT_PARMS::default());
        assert_eq!(cleared, 0);

        let mut ordinary_stereo = 0xff;
        PrepareSaveOptBitsForStructure(
            &mut ordinary_stereo,
            &INPUT_PARMS {
                bINChIOutputOptions: INCHI_OUT_SAVEOPT as i32,
                nMode: u64::from(REQ_MODE_STEREO),
                ..INPUT_PARMS::default()
            },
        );
        assert_eq!(ordinary_stereo, SAVE_OPT_SUU as u8);

        let mut stereo_and_fixed = 0;
        PrepareSaveOptBitsForStructure(
            &mut stereo_and_fixed,
            &INPUT_PARMS {
                bINChIOutputOptions: INCHI_OUT_SAVEOPT as i32,
                nMode: u64::from(REQ_MODE_STEREO | REQ_MODE_DIFF_UU_STEREO | REQ_MODE_BASIC),
                ..INPUT_PARMS::default()
            },
        );
        assert_eq!(
            stereo_and_fixed,
            (SAVE_OPT_SUU | SAVE_OPT_SLUUD | SAVE_OPT_FIXEDH) as u8
        );

        for ignored_mode in [REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU] {
            let mut ignored = 0xff;
            PrepareSaveOptBitsForStructure(
                &mut ignored,
                &INPUT_PARMS {
                    bINChIOutputOptions: INCHI_OUT_SAVEOPT as i32,
                    nMode: u64::from(REQ_MODE_STEREO | ignored_mode),
                    ..INPUT_PARMS::default()
                },
            );
            assert_eq!(ignored, 0);
        }

        let mut tautomer_options = 0;
        PrepareSaveOptBitsForStructure(
            &mut tautomer_options,
            &INPUT_PARMS {
                bINChIOutputOptions: INCHI_OUT_SAVEOPT as i32,
                nMode: u64::from(REQ_MODE_STEREO | REQ_MODE_SB_IGN_ALL_UU),
                bTautFlags: u64::from(
                    TG_FLAG_RECONNECT_COORD
                        | TG_FLAG_KETO_ENOL_TAUT
                        | TG_FLAG_1_5_TAUT
                        | TG_FLAG_PT_22_00
                        | TG_FLAG_PT_16_00
                        | TG_FLAG_PT_06_00
                        | TG_FLAG_PT_39_00
                        | TG_FLAG_PT_13_00
                        | TG_FLAG_PT_18_00,
                ),
                ..INPUT_PARMS::default()
            },
        );
        assert_eq!(
            tautomer_options,
            (SAVE_OPT_RECMET | SAVE_OPT_KET | SAVE_OPT_15T | SAVE_OPT_PT_22_00 | SAVE_OPT_PT_16_00)
                as u8
        );

        let mut stereo_disabled = 0xff;
        PrepareSaveOptBitsForStructure(
            &mut stereo_disabled,
            &INPUT_PARMS {
                bINChIOutputOptions: INCHI_OUT_SAVEOPT as i32,
                nMode: u64::from(REQ_MODE_DIFF_UU_STEREO),
                ..INPUT_PARMS::default()
            },
        );
        assert_eq!(stereo_disabled, 0);
    }

    #[test]
    fn source_port__runichi__displayorigandresultstructuresandcomponents__line_679() {
        for (display, initial_composite, expected_composite) in
            [(0, 17, 17), (1, -9, 1), (-1, 23, 1)]
        {
            let mut parameters = INPUT_PARMS {
                bDisplay: display,
                bDisplayCompositeResults: initial_composite,
                ..INPUT_PARMS::default()
            };
            DisplayOrigAndResultStructuresAndComponents(&mut parameters);
            assert_eq!(parameters.bDisplay, display);
            assert_eq!(parameters.bDisplayCompositeResults, expected_composite);
        }
    }

    #[test]
    fn source_port__runichi__saveokprocessedmolfile__line_785() {
        fn run_case(
            enabled: i32,
            return_code: i32,
            start: i64,
            end: i64,
            problem_present: bool,
            output_present: bool,
        ) -> (Result<(), SourceHeapError>, Vec<u8>, u64) {
            let mut heap = SourceHeap::default();
            let input = heap
                .allocate_model_storage(vec![SourceFile {
                    bytes: b"first\nsecond\n".to_vec(),
                    ..SourceFile::default()
                }])
                .unwrap();
            let output = if output_present {
                heap.allocate_model_storage(vec![SourceFile::default()])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let mut input_stream = INCHI_IOSTREAM {
                f: input,
                type_: INCHI_IOS_TYPE_FILE as i32,
                ..INCHI_IOSTREAM::default()
            };
            let problem_stream = INCHI_IOSTREAM {
                f: output,
                type_: INCHI_IOS_TYPE_FILE as i32,
                ..INCHI_IOSTREAM::default()
            };
            let result = SaveOkProcessedMolfile(
                &mut heap,
                return_code,
                &STRUCT_DATA {
                    fPtrStart: start,
                    fPtrEnd: end,
                    ..STRUCT_DATA::default()
                },
                &INPUT_PARMS {
                    bSaveAllGoodStructsAsProblem: enabled,
                    ..INPUT_PARMS::default()
                },
                problem_present.then_some(&problem_stream),
                Some(&mut input_stream),
            );
            let output_bytes = if output.is_null() {
                Vec::new()
            } else {
                heap.slice(output.as_const()).unwrap()[0].bytes.clone()
            };
            let input_position = heap.slice(input.as_const()).unwrap()[0].position;
            (result, output_bytes, input_position)
        }

        assert_eq!(
            run_case(1, _IS_OKAY as i32, 0, 6, true, true),
            (Ok(()), b"first\n".to_vec(), 6)
        );
        for case in [
            (0, _IS_OKAY as i32, 0, 6, true, true),
            (1, _IS_FATAL as i32, 0, 6, true, true),
            (1, _IS_ERROR as i32, 0, 6, true, true),
            (1, _IS_OKAY as i32, 0, 6, false, true),
            (1, _IS_OKAY as i32, 0, 6, true, false),
            (1, _IS_OKAY as i32, -1, 6, true, true),
            (1, _IS_OKAY as i32, 6, 6, true, true),
            (1, _IS_OKAY as i32, 7, 6, true, true),
        ] {
            assert_eq!(
                run_case(case.0, case.1, case.2, case.3, case.4, case.5),
                (Ok(()), Vec::new(), 0)
            );
        }
    }
}
