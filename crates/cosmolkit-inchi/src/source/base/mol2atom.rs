use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichinorm::FreeInpAtomData;
use crate::source::base::mol_fmt2::{MolfileGetXYZDimAndNormFactors, MolfileHasNoChemStruc};
use crate::source::base::runichi3::OAD_Polymer_Free;
use crate::source::base::util::{
    detect_unusual_el_valence, extract_h_atoms, get_num_H, get_periodic_table_number, inchi_calloc,
    inchi_free, is_el_a_metal, is_in_the_list, mystrncpy_slice, n_bonds_val_to_metal,
};
use crate::source_types::{
    BOND_TYPE_ALTERN, BOND_TYPE_SINGLE, COMP_ATOM_DATA, INP_ATOM_DATA, INPUT_STEREO_DBLE_EITHER,
    INPUT_STEREO_SNGL_DOWN, INPUT_STEREO_SNGL_EITHER, INPUT_STEREO_SNGL_UP, MAX_INPUT_BOND_TYPE,
    MAXVAL, MIN_INPUT_BOND_TYPE, MOL_FMT_DATA, NUM_H_ISOTOPES, OAD_Polymer, OAD_V3000,
    ORIG_ATOM_DATA, POLYMERS_MODERN, RADICAL_SINGLET, RADICAL_TRIPLET, STEREO_DBLE_EITHER,
    STEREO_SNGL_DOWN, STEREO_SNGL_EITHER, STEREO_SNGL_UP, SourceHeap, SourceHeapError,
    SourceMutPointer, ZERO_ATW_DIFF, inp_ATOM, local_util::ERR_ELEM,
    tagFrameShifScheme_FSS_STARS_CYCLED,
};

const SOURCE_SIZEOF_INP_ATOM: u64 = 176;
const SOURCE_SIZEOF_INT: u64 = 4;
const SOURCE_SIZEOF_POINTER: u64 = 8;
const SOURCE_SIZEOF_OAD_POLYMER: u64 = 48;
const SOURCE_SIZEOF_OAD_POLYMER_UNIT: u64 = 240;
const SOURCE_SIZEOF_OAD_V3000: u64 = 104;

#[allow(non_snake_case)]
pub(crate) fn FreeInpAtom(
    heap: &mut SourceHeap,
    atom_slot: Option<&mut SourceMutPointer<inp_ATOM>>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1046 FreeInpAtom
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
    // END INCHI C FUNCTION: FreeInpAtom

    if let Some(atom_slot) = atom_slot
        && !atom_slot.is_null()
    {
        inchi_free(heap, *atom_slot)?;
        *atom_slot = SourceMutPointer::null();
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn CreateInpAtom(
    heap: &mut SourceHeap,
    num_atoms: i32,
) -> Result<SourceMutPointer<inp_ATOM>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1033 CreateInpAtom
    // INCHI✔️❌: inp_ATOM *CreateInpAtom(int num_atoms)
    // INCHI✔️❌: {
    // INCHI✔️❌:     void *p = NULL;
    // INCHI✔️❌:     if (num_atoms >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         p = inchi_calloc(num_atoms, sizeof(inp_ATOM));
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     // void *p = inchi_calloc(num_atoms, sizeof(inp_ATOM));
    // INCHI✔️❌:     return (inp_ATOM *)p;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateInpAtom

    if num_atoms < 0 {
        return Ok(SourceMutPointer::null());
    }
    match inchi_calloc(
        heap,
        u64::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
        SOURCE_SIZEOF_INP_ATOM,
    ) {
        Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
        result => result,
    }
}

#[allow(non_snake_case)]
pub(crate) fn CreateInpAtomData(
    heap: &mut SourceHeap,
    input: &mut INP_ATOM_DATA,
    num_atoms: i32,
    create_at_fixed_bonds: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1071 CreateInpAtomData
    // INCHI✔️❌: int CreateInpAtomData(INP_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                       int num_atoms,
    // INCHI✔️❌:                       int create_at_fixed_bonds)
    // INCHI✔️❌: {
    // INCHI✔️❌:     FreeInpAtomData(inp_at_data);
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((inp_at_data->at = CreateInpAtom(num_atoms)) &&
    // INCHI✔️❌:         (!create_at_fixed_bonds || (inp_at_data->at_fixed_bonds = CreateInpAtom(num_atoms))))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inp_at_data->num_at = num_atoms;
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     FreeInpAtomData(inp_at_data);
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateInpAtomData

    FreeInpAtomData(heap, input)?;
    input.at = CreateInpAtom(heap, num_atoms)?;
    if !input.at.is_null() {
        if create_at_fixed_bonds == 0 {
            input.num_at = num_atoms;
            return Ok(1);
        }
        input.at_fixed_bonds = CreateInpAtom(heap, num_atoms)?;
        if !input.at_fixed_bonds.is_null() {
            input.num_at = num_atoms;
            return Ok(1);
        }
    }
    FreeInpAtomData(heap, input)?;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn FreeCompAtomData(
    heap: &mut SourceHeap,
    input_atom_data: &mut COMP_ATOM_DATA,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1090 FreeCompAtomData
    // INCHI✔️❌: void FreeCompAtomData(COMP_ATOM_DATA *inp_at_data)
    // INCHI✔️❌: {
    // INCHI✔️❌:     FreeInpAtom(&inp_at_data->at);
    // INCHI✔️❌:
    // INCHI✔️❌:     if (inp_at_data->nOffsetAtAndH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(inp_at_data->nOffsetAtAndH);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     memset(inp_at_data, 0, sizeof(*inp_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeCompAtomData

    FreeInpAtom(heap, Some(&mut input_atom_data.at))?;
    if !input_atom_data.nOffsetAtAndH.is_null() {
        inchi_free(heap, input_atom_data.nOffsetAtAndH)?;
    }
    *input_atom_data = COMP_ATOM_DATA::default();
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn FreeOrigAtData(
    heap: &mut SourceHeap,
    orig_at_data: Option<&mut ORIG_ATOM_DATA>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1224 FreeOrigAtData
    // INCHI✔️❌: void FreeOrigAtData(ORIG_ATOM_DATA *orig_at_data)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (!orig_at_data)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     FreeInpAtom(&orig_at_data->at);
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_at_data->nCurAtLen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(orig_at_data->nCurAtLen);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_at_data->nOldCompNumber)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(orig_at_data->nOldCompNumber);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_at_data->szCoord)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(orig_at_data->szCoord);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_at_data->nEquLabels)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(orig_at_data->nEquLabels);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_at_data->nSortedOrder)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(orig_at_data->nSortedOrder);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* v 1.05 */
    // INCHI✔️❌:     FreeExtOrigAtData(orig_at_data->polymer, orig_at_data->v3000);
    // INCHI✔️❌:
    // INCHI✔️❌:     memset(orig_at_data, 0, sizeof(*orig_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeOrigAtData

    let Some(orig_at_data) = orig_at_data else {
        return Ok(());
    };

    FreeInpAtom(heap, Some(&mut orig_at_data.at))?;
    if !orig_at_data.nCurAtLen.is_null() {
        inchi_free(heap, orig_at_data.nCurAtLen)?;
    }
    if !orig_at_data.nOldCompNumber.is_null() {
        inchi_free(heap, orig_at_data.nOldCompNumber)?;
    }
    if !orig_at_data.szCoord.is_null() {
        inchi_free(heap, orig_at_data.szCoord)?;
    }
    if !orig_at_data.nEquLabels.is_null() {
        inchi_free(heap, orig_at_data.nEquLabels)?;
    }
    if !orig_at_data.nSortedOrder.is_null() {
        inchi_free(heap, orig_at_data.nSortedOrder)?;
    }
    FreeExtOrigAtData(heap, orig_at_data.polymer, orig_at_data.v3000)?;
    *orig_at_data = ORIG_ATOM_DATA::default();
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn FreeExtOrigAtData(
    heap: &mut SourceHeap,
    mut pd: SourceMutPointer<OAD_Polymer>,
    v3k: SourceMutPointer<OAD_V3000>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1269 FreeExtOrigAtData
    // INCHI✔❌: void FreeExtOrigAtData(OAD_Polymer *pd, OAD_V3000 *v3k)
    // INCHI✔❌: {
    // INCHI✔❌:     int k;
    // INCHI✔❌:
    // INCHI✔❌:     OAD_Polymer_Free(pd);
    // INCHI✔❌:     pd = NULL;
    // INCHI✔❌:
    // INCHI✔❌:     if (v3k)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (v3k->atom_index_orig)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free(v3k->atom_index_orig);
    // INCHI✔❌:             v3k->atom_index_orig = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (v3k->atom_index_fin)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free(v3k->atom_index_fin);
    // INCHI✔❌:             v3k->atom_index_fin = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (v3k->n_haptic_bonds && v3k->lists_haptic_bonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             for (k = 0; k < v3k->n_haptic_bonds; k++)
    // INCHI✔❌:                 if (v3k->lists_haptic_bonds[k])
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free(v3k->lists_haptic_bonds[k]);
    // INCHI✔❌:                     v3k->lists_haptic_bonds[k] = NULL;
    // INCHI✔❌:                 }
    // INCHI✔❌:             inchi_free(v3k->lists_haptic_bonds);
    // INCHI✔❌:             v3k->lists_haptic_bonds = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (v3k->n_steabs && v3k->lists_steabs)
    // INCHI✔❌:         {
    // INCHI✔❌:             for (k = 0; k < v3k->n_steabs; k++)
    // INCHI✔❌:                 if (v3k->lists_steabs[k])
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free(v3k->lists_steabs[k]);
    // INCHI✔❌:                     v3k->lists_steabs[k] = NULL;
    // INCHI✔❌:                 }
    // INCHI✔❌:             inchi_free(v3k->lists_steabs);
    // INCHI✔❌:             v3k->lists_steabs = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (v3k->n_sterel && v3k->lists_sterel)
    // INCHI✔❌:         {
    // INCHI✔❌:             for (k = 0; k < v3k->n_sterel; k++)
    // INCHI✔❌:                 if (v3k->lists_sterel[k])
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free(v3k->lists_sterel[k]);
    // INCHI✔❌:                     v3k->lists_sterel[k] = NULL;
    // INCHI✔❌:                 }
    // INCHI✔❌:             inchi_free(v3k->lists_sterel);
    // INCHI✔❌:             v3k->lists_sterel = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (v3k->n_sterac && v3k->lists_sterac)
    // INCHI✔❌:         {
    // INCHI✔❌:             for (k = 0; k < v3k->n_sterac; k++)
    // INCHI✔❌:                 if (v3k->lists_sterac[k])
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free(v3k->lists_sterac[k]);
    // INCHI✔❌:                     v3k->lists_sterac[k] = NULL;
    // INCHI✔❌:                 }
    // INCHI✔❌:             inchi_free(v3k->lists_sterac);
    // INCHI✔❌:             v3k->lists_sterac = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         memset(v3k, 0, sizeof(*v3k)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:         inchi_free(v3k);
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: FreeExtOrigAtData

    OAD_Polymer_Free(heap, pd)?;
    pd = SourceMutPointer::null();
    let _ = pd;

    if !v3k.is_null() {
        let (atom_index_orig, atom_index_fin) = {
            let value = heap
                .slice(v3k.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (value.atom_index_orig, value.atom_index_fin)
        };
        if !atom_index_orig.is_null() {
            inchi_free(heap, atom_index_orig)?;
            heap.slice_mut(v3k)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .atom_index_orig = SourceMutPointer::null();
        }
        if !atom_index_fin.is_null() {
            inchi_free(heap, atom_index_fin)?;
            heap.slice_mut(v3k)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .atom_index_fin = SourceMutPointer::null();
        }

        let (count, lists) = {
            let value = heap
                .slice(v3k.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (value.n_haptic_bonds, value.lists_haptic_bonds)
        };
        if count != 0 && !lists.is_null() {
            for k in 0..count {
                let item = heap.slice(lists.as_const())?[k as usize];
                if !item.is_null() {
                    inchi_free(heap, item)?;
                    heap.slice_mut(lists)?[k as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, lists)?;
            heap.slice_mut(v3k)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .lists_haptic_bonds = SourceMutPointer::null();
        }

        let (count, lists) = {
            let value = heap
                .slice(v3k.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (value.n_steabs, value.lists_steabs)
        };
        if count != 0 && !lists.is_null() {
            for k in 0..count {
                let item = heap.slice(lists.as_const())?[k as usize];
                if !item.is_null() {
                    inchi_free(heap, item)?;
                    heap.slice_mut(lists)?[k as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, lists)?;
            heap.slice_mut(v3k)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .lists_steabs = SourceMutPointer::null();
        }

        let (count, lists) = {
            let value = heap
                .slice(v3k.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (value.n_sterel, value.lists_sterel)
        };
        if count != 0 && !lists.is_null() {
            for k in 0..count {
                let item = heap.slice(lists.as_const())?[k as usize];
                if !item.is_null() {
                    inchi_free(heap, item)?;
                    heap.slice_mut(lists)?[k as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, lists)?;
            heap.slice_mut(v3k)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .lists_sterel = SourceMutPointer::null();
        }

        let (count, lists) = {
            let value = heap
                .slice(v3k.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (value.n_sterac, value.lists_sterac)
        };
        if count != 0 && !lists.is_null() {
            for k in 0..count {
                let item = heap.slice(lists.as_const())?[k as usize];
                if !item.is_null() {
                    inchi_free(heap, item)?;
                    heap.slice_mut(lists)?[k as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, lists)?;
            heap.slice_mut(v3k)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .lists_sterac = SourceMutPointer::null();
        }

        *heap
            .slice_mut(v3k)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = OAD_V3000::default();
        inchi_free(heap, v3k)?;
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn SetExtOrigAtDataByMolfileExtInput(
    heap: &mut SourceHeap,
    mfdata: &MOL_FMT_DATA,
    pp_polymer: &mut SourceMutPointer<OAD_Polymer>,
    pp_v3000: &mut SourceMutPointer<OAD_V3000>,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1340 SetExtOrigAtDataByMolfileExtInput
    // INCHI✔️❌: complete source frame follows verbatim; Rust ownership requires checked heap accesses and cloned source records.
    /*
int SetExtOrigAtDataByMolfileExtInput(MOL_FMT_DATA *mfdata,
                                      OAD_Polymer **ppPolymer,
                                      OAD_V3000 **ppV3000,
                                      char *pStrErr)
{
    int k, m, err = 0;
    OAD_V3000 *pv = NULL;
    int nsgroups = mfdata->ctab.sgroups.used;

    /* djb-rwth: addressing coverity ID #499476 -- TREAT_ERR properly used in all cases */
    /* Polymers */
    if (nsgroups > 0)
    {
        /* Prepare OAD_Polymer container */
        *ppPolymer = (OAD_Polymer *)inchi_calloc(1, sizeof(OAD_Polymer));
        if (!(*ppPolymer))
        {
            TREAT_ERR(err, 9001, "Out of RAM");
            goto exit_function;
        }

        /* Convert Molfile's Sgroup's to OAD_PolymerUnit's */
        (*ppPolymer)->units = (OAD_PolymerUnit **)inchi_calloc(nsgroups, sizeof((*ppPolymer)->units[0]));
        if (!(*ppPolymer)->units)
        {
            TREAT_ERR(err, 9001, "Out of RAM");
            goto exit_function;
        }
        memset((*ppPolymer)->units, 0, sizeof(*(*ppPolymer)->units)); /* djb-rwth: memset_s C11/Annex K variant? */

        (*ppPolymer)->n = nsgroups;
        (*ppPolymer)->is_in_reconn = 0;
        (*ppPolymer)->pzz = NULL;
        (*ppPolymer)->edit_repeats = 0;
        (*ppPolymer)->really_do_frame_shift = 0;
        (*ppPolymer)->frame_shift_scheme = FSS_STARS_CYCLED;
        (*ppPolymer)->treat = POLYMERS_MODERN;

        for (k = 0; k < nsgroups; k++)
        {
            int q = 0;
            MOL_FMT_SGROUP *groupk = mfdata->ctab.sgroups.group[k];

            OAD_PolymerUnit *unitk = (*ppPolymer)->units[k] = (OAD_PolymerUnit *)inchi_calloc(1, sizeof(OAD_PolymerUnit));

            if (!unitk)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }

            memset(unitk, 0, sizeof(*unitk)); /* djb-rwth: memset_s C11/Annex K variant? */
            unitk->id = groupk->id;
            unitk->type = groupk->type;
            unitk->subtype = groupk->subtype;
            unitk->conn = groupk->conn;
            unitk->label = groupk->label;

            for (q = 0; q < 4; q++)
            {
                unitk->xbr1[q] = groupk->xbr1[q];
                unitk->xbr2[q] = groupk->xbr2[q];
            }
            strcpy(unitk->smt, groupk->smt);
            unitk->na = groupk->alist.used;
            unitk->alist = (int *)inchi_calloc(unitk->na, sizeof(int));
            if (!unitk->alist)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < unitk->na; m++)
            {
                unitk->alist[m] = groupk->alist.item[m];
            }
            unitk->nb = groupk->blist.used;
            if (unitk->nb > 0)
            {
                unitk->blist = (int *)inchi_calloc(2 * (long long)unitk->nb, sizeof(int)); /* djb-rwth: cast operator added */
                if (!unitk->blist)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (m = 0; m < groupk->blist.used; m++)
                {
                    int ib, ia1, ia2;
                    ib = groupk->blist.item[m];
                    if (ib < 1 || ib > mfdata->ctab.n_bonds)
                    {
                        TREAT_ERR(err, 9004, "Polymer unit in Molfile refers to invalid bond");
                        goto exit_function;
                    }
                    ia1 = mfdata->ctab.bonds[ib - 1].atnum1;
                    ia2 = mfdata->ctab.bonds[ib - 1].atnum2;
                    unitk->blist[2 * m] = ia1;
                    unitk->blist[2 * m + 1] = ia2;
                    if (!strcmp(mfdata->ctab.atoms[ia1 - 1].symbol, "H") ||
                        !strcmp(mfdata->ctab.atoms[ia2 - 1].symbol, "H"))
                    {
                        TREAT_ERR(err, 9002, "Hydrogen as polymer end group is not supported");
                        goto exit_function;
                    }
                }
            }
            else
            {
                unitk->blist = NULL;
            }
        }
    }

    /* V3000 Extensions */
    if (mfdata->ctab.v3000)
    {
        int nn;
        MOL_FMT_v3000 *mpv = mfdata->ctab.v3000;

        *ppV3000 = (OAD_V3000 *)inchi_calloc(1, sizeof(OAD_V3000));
        pv = *ppV3000;
        if (!pv)
        {
            TREAT_ERR(err, 9001, "Out of RAM");
            goto exit_function;
        }
        memset(pv, 0, sizeof(*pv)); /* djb-rwth: memset_s C11/Annex K variant? */

        pv->n_collections = mpv->n_collections;
        pv->n_haptic_bonds = mpv->n_haptic_bonds;
        pv->n_non_haptic_bonds = mpv->n_non_haptic_bonds;
        pv->n_sgroups = mpv->n_sgroups;
        pv->n_non_star_atoms = mpv->n_non_star_atoms;
        pv->n_star_atoms = mpv->n_star_atoms;
        pv->n_steabs = mpv->n_steabs;
        pv->n_sterac = mpv->n_sterac;
        pv->n_sterel = mpv->n_sterel;
        pv->n_3d_constraints = mpv->n_3d_constraints;

        if (mpv->atom_index_orig)
        {
            pv->atom_index_orig = (int *)inchi_calloc(mfdata->ctab.n_atoms, sizeof(int));
            if (NULL == pv->atom_index_orig)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            memcpy(pv->atom_index_orig, mpv->atom_index_orig, mfdata->ctab.n_atoms);
        }
        if (mpv->atom_index_fin)
        {
            pv->atom_index_fin = (int *)inchi_calloc(mfdata->ctab.n_atoms, sizeof(int));
            if (NULL == pv->atom_index_fin)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            memcpy(pv->atom_index_fin, mpv->atom_index_fin, mfdata->ctab.n_atoms);
        }
        if (mpv->n_haptic_bonds && mpv->haptic_bonds)
        {
            pv->lists_haptic_bonds = (int **)inchi_calloc(mpv->n_haptic_bonds, sizeof(int *));
            if (NULL == pv->lists_haptic_bonds)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_haptic_bonds; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->haptic_bonds->lists[m];
                nn = mol_lst[2] + 3;
                lst = pv->lists_haptic_bonds[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (mpv->n_steabs && mpv->steabs)
        {
            pv->lists_steabs = (int **)inchi_calloc(mpv->n_steabs, sizeof(int *));
            if (NULL == pv->lists_steabs)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_steabs; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->steabs->lists[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_steabs[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (mpv->n_sterac && mpv->sterac)
        {
            pv->lists_sterac = (int **)inchi_calloc(mpv->n_sterac, sizeof(int *));
            if (NULL == pv->lists_sterac)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_sterac; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->sterac->lists[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_sterac[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (mpv->n_sterel && mpv->sterel)
        {
            pv->lists_sterel = (int **)inchi_calloc(mpv->n_sterel, sizeof(int *));
            if (NULL == pv->lists_sterel)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_sterel; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->sterel->lists[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_sterel[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
    }

exit_function:
    if (err)
    {
        FreeExtOrigAtData((*ppPolymer), pv);
        *ppPolymer = NULL;
    }

    return err;
}
    */
    // END INCHI C FUNCTION: SetExtOrigAtDataByMolfileExtInput

    let mut err = 0;
    let mut pv = SourceMutPointer::null();
    let nsgroups = mfdata.ctab.sgroups.used;

    macro_rules! treat_err {
        ($new_err:expr, $message:literal) => {{
            if err == 0 && $new_err != 0 {
                err = $new_err;
            }
            let message = $message
                .bytes()
                .map(|byte| byte as i8)
                .chain([0])
                .collect::<Vec<_>>();
            AddErrorMessage(p_str_err.as_deref_mut(), Some(&message))?;
        }};
    }

    'conversion: {
        if nsgroups > 0 {
            *pp_polymer = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_POLYMER) {
                Ok(pointer) => pointer,
                Err(_) => {
                    treat_err!(9001, "Out of RAM");
                    SourceMutPointer::null()
                }
            };
            if pp_polymer.is_null() {
                break 'conversion;
            }

            let units = match inchi_calloc(
                heap,
                u64::try_from(nsgroups).unwrap_or(u64::MAX),
                SOURCE_SIZEOF_POINTER,
            ) {
                Ok(pointer) => pointer,
                Err(_) => {
                    treat_err!(9001, "Out of RAM");
                    break 'conversion;
                }
            };
            {
                let polymer = &mut heap.slice_mut(*pp_polymer)?[0];
                polymer.units = units;
                polymer.n = nsgroups;
                polymer.is_in_reconn = 0;
                polymer.pzz = SourceMutPointer::null();
                polymer.edit_repeats = 0;
                polymer.really_do_frame_shift = 0;
                polymer.frame_shift_scheme = tagFrameShifScheme_FSS_STARS_CYCLED as i32;
                polymer.treat = POLYMERS_MODERN as i32;
            }

            for k in 0..nsgroups {
                let group_pointer = heap.slice(mfdata.ctab.sgroups.group.as_const())?[k as usize];
                let group = heap
                    .slice(group_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let unit = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_POLYMER_UNIT) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        treat_err!(9001, "Out of RAM");
                        break 'conversion;
                    }
                };
                heap.slice_mut(units)?[k as usize] = unit;
                {
                    let output = &mut heap.slice_mut(unit)?[0];
                    output.id = group.id;
                    output.type_ = group.type_;
                    output.subtype = group.subtype;
                    output.conn = group.conn;
                    output.label = group.label;
                    output.xbr1 = group.xbr1;
                    output.xbr2 = group.xbr2;
                    let nul = group
                        .smt
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or(SourceHeapError::MissingNulTerminator)?;
                    output.smt[..=nul].copy_from_slice(&group.smt[..=nul]);
                    output.na = group.alist.used;
                }

                let alist = match inchi_calloc(
                    heap,
                    u64::try_from(group.alist.used).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_INT,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        treat_err!(9001, "Out of RAM");
                        break 'conversion;
                    }
                };
                heap.slice_mut(unit)?[0].alist = alist;
                for m in 0..group.alist.used {
                    heap.slice_mut(alist)?[m as usize] =
                        heap.slice(group.alist.item.as_const())?[m as usize];
                }

                heap.slice_mut(unit)?[0].nb = group.blist.used;
                if group.blist.used > 0 {
                    let blist_count = u64::try_from(i64::from(group.blist.used) * 2)
                        .unwrap_or(u64::MAX);
                    let blist = match inchi_calloc(heap, blist_count, SOURCE_SIZEOF_INT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            treat_err!(9001, "Out of RAM");
                            break 'conversion;
                        }
                    };
                    heap.slice_mut(unit)?[0].blist = blist;
                    for m in 0..group.blist.used {
                        let ib = heap.slice(group.blist.item.as_const())?[m as usize];
                        if ib < 1 || ib > mfdata.ctab.n_bonds {
                            treat_err!(9004, "Polymer unit in Molfile refers to invalid bond");
                            break 'conversion;
                        }
                        let bond = heap.slice(mfdata.ctab.bonds.as_const())?[(ib - 1) as usize].clone();
                        let ia1 = i32::from(bond.atnum1);
                        let ia2 = i32::from(bond.atnum2);
                        heap.slice_mut(blist)?[(2 * m) as usize] = ia1;
                        heap.slice_mut(blist)?[(2 * m + 1) as usize] = ia2;
                        let atoms = heap.slice(mfdata.ctab.atoms.as_const())?;
                        let first_is_h = atoms[(ia1 - 1) as usize].symbol[0] == b'H' as i8
                            && atoms[(ia1 - 1) as usize].symbol[1] == 0;
                        let second_is_h = atoms[(ia2 - 1) as usize].symbol[0] == b'H' as i8
                            && atoms[(ia2 - 1) as usize].symbol[1] == 0;
                        if first_is_h || second_is_h {
                            treat_err!(9002, "Hydrogen as polymer end group is not supported");
                            break 'conversion;
                        }
                    }
                } else {
                    heap.slice_mut(unit)?[0].blist = SourceMutPointer::null();
                }
            }
        }

        if !mfdata.ctab.v3000.is_null() {
            let input = heap
                .slice(mfdata.ctab.v3000.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            *pp_v3000 = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_V3000) {
                Ok(pointer) => pointer,
                Err(_) => {
                    treat_err!(9001, "Out of RAM");
                    SourceMutPointer::null()
                }
            };
            pv = *pp_v3000;
            if pv.is_null() {
                break 'conversion;
            }
            {
                let output = &mut heap.slice_mut(pv)?[0];
                output.n_collections = input.n_collections;
                output.n_haptic_bonds = input.n_haptic_bonds;
                output.n_non_haptic_bonds = input.n_non_haptic_bonds;
                output.n_sgroups = input.n_sgroups;
                output.n_non_star_atoms = input.n_non_star_atoms;
                output.n_star_atoms = input.n_star_atoms;
                output.n_steabs = input.n_steabs;
                output.n_sterac = input.n_sterac;
                output.n_sterel = input.n_sterel;
                output.n_3d_constraints = input.n_3d_constraints;
            }

            let atom_count = u64::try_from(mfdata.ctab.n_atoms).unwrap_or(u64::MAX);
            if !input.atom_index_orig.is_null() {
                let destination = match inchi_calloc(heap, atom_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        treat_err!(9001, "Out of RAM");
                        break 'conversion;
                    }
                };
                heap.slice_mut(pv)?[0].atom_index_orig = destination;
                for byte_index in 0..usize::try_from(mfdata.ctab.n_atoms).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_orig.as_const())?[word] as u32 >> shift) & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }
            if !input.atom_index_fin.is_null() {
                let destination = match inchi_calloc(heap, atom_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        treat_err!(9001, "Out of RAM");
                        break 'conversion;
                    }
                };
                heap.slice_mut(pv)?[0].atom_index_fin = destination;
                for byte_index in 0..usize::try_from(mfdata.ctab.n_atoms).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_fin.as_const())?[word] as u32 >> shift) & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }

            macro_rules! copy_lists {
                ($count:expr, $source_lists:expr, $destination_field:ident, $length_index:expr, $extra:expr) => {
                    if $count != 0 && !$source_lists.is_null() {
                        let outer = match inchi_calloc(
                            heap,
                            u64::try_from($count).unwrap_or(u64::MAX),
                            SOURCE_SIZEOF_POINTER,
                        ) {
                            Ok(pointer) => pointer,
                            Err(_) => {
                                treat_err!(9001, "Out of RAM");
                                break 'conversion;
                            }
                        };
                        heap.slice_mut(pv)?[0].$destination_field = outer;
                        let source_outer = heap
                            .slice($source_lists.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .lists;
                        for m in 0..$count {
                            let source = heap.slice(source_outer.as_const())?[m as usize];
                            let nn = heap.slice(source.as_const())?[$length_index].wrapping_add($extra);
                            let row = match inchi_calloc(
                                heap,
                                u64::try_from(nn).unwrap_or(u64::MAX),
                                SOURCE_SIZEOF_INT,
                            ) {
                                Ok(pointer) => pointer,
                                Err(_) => {
                                    treat_err!(9001, "Out of RAM");
                                    break 'conversion;
                                }
                            };
                            heap.slice_mut(outer)?[m as usize] = row;
                            for k in 0..nn {
                                heap.slice_mut(row)?[k as usize] =
                                    heap.slice(source.as_const())?[k as usize];
                            }
                        }
                    }
                };
            }

            copy_lists!(input.n_haptic_bonds, input.haptic_bonds, lists_haptic_bonds, 2, 3);
            copy_lists!(input.n_steabs, input.steabs, lists_steabs, 1, 2);
            copy_lists!(input.n_sterac, input.sterac, lists_sterac, 1, 2);
            copy_lists!(input.n_sterel, input.sterel, lists_sterel, 1, 2);
        }
    }

    if err != 0 {
        FreeExtOrigAtData(heap, *pp_polymer, pv)?;
        *pp_polymer = SourceMutPointer::null();
    }
    Ok(err)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MakeInpAtomsFromMolfileData(
    heap: &mut SourceHeap,
    mfdata: &MOL_FMT_DATA,
    num_atoms: &mut i32,
    num_bonds: &mut i32,
    at_inp: SourceMutPointer<inp_ATOM>,
    b_do_not_add_h: i32,
    err: &mut i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<SourceMutPointer<inp_ATOM>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:508 MakeInpAtomsFromMolfileData
    // INCHI✔️❌: complete source frame follows verbatim; checked heap slices and the temporary atom copy add overhead.
    /*
inp_ATOM *MakeInpAtomsFromMolfileData(MOL_FMT_DATA *mfdata,
                                      int *num_atoms,
                                      int *num_bonds,
                                      inp_ATOM *at_inp,
                                      int bDoNotAddH,
                                      int *err,
                                      char *pStrErr)
{
    inp_ATOM *at = NULL;
    /* char      *bond_stereo = NULL; */
    AT_NUMB *p1, *p2;
    int i, a1, a2, n1, n2, bonds, iso_atw_diff;
    char bond_stereo, bond_type;

    *err = 0;

    *num_atoms = mfdata->ctab.n_atoms;
    *num_bonds = 0;

    /*(@nnuk : Nauman Ullah Khan) */
    LOG_NO_ARGS("\n############### (L579:mol2atom.c) ################\n");
    LOG_MULT_ARGS("Number of atoms : %d\n", *num_atoms);
    LOG_NO_ARGS("####################################################\n");

    if (MolfileHasNoChemStruc(mfdata))
    {
        goto exit_function;
    }

    /* Allocate memory if necessary */
    if (at_inp)
    {
        at = at_inp;
    }
    else
    {
        at = CreateInpAtom(*num_atoms);
        if (!at)
        {
            *err = -1;
            TREAT_ERR_AND_FIN(*err, -1, exit_function, "Out of RAM");
        }
    }

    /* Copy atoms info */

    /*(@nnuk : Nauman Ullah Khan) */
    LOG_NO_ARGS("\n##################### Atoms Data #########################\n");

    for (i = 0; i < *num_atoms; i++)
    {

        mystrncpy(at[i].elname, mfdata->ctab.atoms[i].symbol, sizeof(at->elname));
        /* at[i].chem_bonds_valence = mfdata->ctab.atoms[i].valence; */
        /*  MOLfile valence; will change */

        at[i].orig_at_number = (AT_NUMB)(i + 1);
        at[i].iso_atw_diff = mfdata->ctab.atoms[i].mass_difference;
        at[i].charge = mfdata->ctab.atoms[i].charge;
        at[i].radical = mfdata->ctab.atoms[i].radical;

        /* see SetInpAtomXYZ()
        at[i].x                  = mfdata->ctab.atoms[i].fx;
        at[i].y                  = mfdata->ctab.atoms[i].fy;
        at[i].z                  = mfdata->ctab.atoms[i].fz;
        */

        iso_atw_diff = mfdata->ctab.atoms[i].mass_difference;

        at[i].iso_atw_diff = iso_atw_diff == ZERO_ATW_DIFF
                                 ? 1
                             : iso_atw_diff > 0 ? iso_atw_diff + 1
                                                : iso_atw_diff;

#if (SINGLET_IS_TRIPLET == 1)
        if (at[i].radical == RADICAL_SINGLET)
        {
            at[i].radical = RADICAL_TRIPLET;
        }
#endif
#if (bRELEASE_VERSION != 1)
        if (isdigit(at[i].elname[0]))
        { /*  for testing */
            mystrncpy(at[i].elname, "C", sizeof(at->elname));
        }
#endif

        if (ERR_ELEM == (n1 = get_periodic_table_number(at[i].elname)))
        {
            /*  Case when elname contains more than 1 element: extract number of H if possible */
            at[i].num_H = extract_H_atoms(at[i].elname, at[i].num_iso_H);

            if (!at[i].elname[0] && NUMH(at, i))
            {
                /* alias contains only H. Added 2004-07-21, fixed 2004-07-22
                 * move the heaviest isotope to the "central atom"
                 * Note: this must be consistent with H-H treatment in remove_terminal_HDT()
                 */
                strcpy(at[i].elname, "H");
                if (NUM_ISO_H(at, i))
                {
                    int j;
                    for (j = NUM_H_ISOTOPES - 1; 0 <= j; j--)
                    {
                        if (at[i].num_iso_H[j])
                        {
                            at[i].num_iso_H[j]--;
                            at[i].iso_atw_diff = 1 + j;
                            break;
                        }
                    }
                }
                else
                {
                    at[i].num_H--;
                }
            }
            if (ERR_ELEM == (n1 = get_periodic_table_number(at[i].elname)))
            {
                n1 = 0;
            }
        } /* if ( ERR_ELEM == */

        at[i].el_number = (U_CHAR)n1;
        if (!n1)
        {
            *err = -2;
            TREAT_ERR(*err, -2, "Unknown element(s):");
            TREAT_ERR_AND_FIN(*err, -2, exit_function, at[i].elname);
        }
        else
        {
            /* replace explicit D or T with isotopic H (added 2003-06-02) */
            if ((int)EL_NUMBER_H == n1 && !at[i].iso_atw_diff)
            {
                switch (at[i].elname[0])
                {
                case 'D':
                    at[i].iso_atw_diff = 2;
                    mystrncpy(at[i].elname, "H", sizeof(at->elname));
                    break;
                case 'T':
                    at[i].iso_atw_diff = 3;
                    mystrncpy(at[i].elname, "H", sizeof(at->elname));
                    break;
                }
            }
        }

        /*(@nnuk : Nauman Ullah Khan) */
        LOG_NO_ARGS("\n############################## (L714:mol2atom.c) #####################################\n");
        LOG_MULT_ARGS("Atom %d: element=%s, x=%f, y=%f, z=%f, chrg=%d, rad=%d, iso=%d\n", i, at[i].elname, at[i].x, at[i].y, at[i].z, at[i].charge, at[i].radical, at[i].iso_atw_diff);
        LOG_NO_ARGS("########################################################################################\n");

    } /* eof copy atom info */

    /*---------------- stereo information notes. ------------------------

    Currently:  1. stereo sign
    =========   --------------
    MOLfile     (atom number = MOLfile atom number - 1, no stdata as an intermediate)
    |        if mfdata->ctab.bonds[i].atnum1 < mfdata->ctab.bonds[i].atnum2
    v        then
    inp_ATOM        stereo > 0
    else
    stereo < 0

    2. neighbor z-coordinate
    ------------------------
    neighbor z-coord > 0 for Up if sign(stdata_bond_no) = sign(at[i].neighbor[j]-i)

    --------------------------------------------------------------------*/

    /* Copy bond info */

    LOG_NO_ARGS("\n######################### Bonds Data ###############################\n");

    for (i = 0, bonds = 0; i < mfdata->ctab.n_bonds; i++)
    {

        bond_stereo = mfdata->ctab.bonds[i].bond_stereo;
        bond_type = mfdata->ctab.bonds[i].bond_type;

        a1 = mfdata->ctab.bonds[i].atnum1 - 1;
        a2 = mfdata->ctab.bonds[i].atnum2 - 1;

        if (a1 < 0 || a1 >= *num_atoms ||
            a2 < 0 || a2 >= *num_atoms ||
            a1 == a2)
        {
            *err |= 1; /*  bond for impossible atom number(s); ignored */
            TREAT_ERR(*err, 0, "Bond to nonexistent atom");
            continue;
        }

        /*  check for multiple bonds between same atoms */
        p1 = is_in_the_list(at[a1].neighbor, (AT_NUMB)a2, at[a1].valence);
        p2 = is_in_the_list(at[a2].neighbor, (AT_NUMB)a1, at[a2].valence);

        /*(@nnuk : Nauman Ullah Khan) */
        LOG_NO_ARGS("\n################## (L771:mol2atom.c) ##################\n");
        LOG_MULT_ARGS("Valence = %d\n", at[i].valence);
        LOG_NO_ARGS("#########################################################\n");

        if ((p1 || p2) && (p1 || at[a1].valence < MAXVAL) && (p2 || at[a2].valence < MAXVAL))
        {
            n1 = p1 ? (p1 - at[a1].neighbor) : at[a1].valence++;
            n2 = p2 ? (p2 - at[a2].neighbor) : at[a2].valence++;
            TREAT_ERR(*err, 0, "Multiple bonds between two atoms");
            *err |= 2; /*  multiple bonds between atoms */
        }
        else if (!p1 && !p2 && at[a1].valence < MAXVAL && at[a2].valence < MAXVAL)
        {
            n1 = at[a1].valence++;
            n2 = at[a2].valence++;
            bonds++;
        }
        else
        {
            char szMsg[64];
            *err |= 4; /*  too large number of bonds. Some bonds ignored. */
            sprintf(szMsg, "Atom '%s' has more than %d bonds",
                    at[a1].valence >= MAXVAL ? at[a1].elname : at[a2].elname, MAXVAL);
            TREAT_ERR(*err, 0, szMsg);
            continue;
        }

        if (bond_type < MIN_INPUT_BOND_TYPE || bond_type > MAX_INPUT_BOND_TYPE)
        {
            char szBondType[16];
            sprintf(szBondType, "%d", bond_type);
            bond_type = 1;
            TREAT_ERR(*err, 0, "Unrecognized bond type:");
            TREAT_ERR(*err, 0, szBondType);
            *err |= 8; /*  Unrecognized Bond type replaced with single bond */
        }

        /* bond type */
        if (n1 < MAXVAL && n2 < MAXVAL) /* djb-rwth: buffer overrun prevention */
        {
            at[a1].bond_type[n1] = at[a2].bond_type[n2] = bond_type;

            /* connection */
            at[a1].neighbor[n1] = (AT_NUMB)a2;
            at[a2].neighbor[n2] = (AT_NUMB)a1;

            /* stereo */
            if (bond_stereo == INPUT_STEREO_DBLE_EITHER /* 3 */)
            {
                at[a1].bond_stereo[n1] =
                    at[a2].bond_stereo[n2] =
                        STEREO_DBLE_EITHER;
            }
            else if (bond_stereo == INPUT_STEREO_SNGL_UP ||     /* 1 */
                     bond_stereo == INPUT_STEREO_SNGL_EITHER || /* 4 */
                     bond_stereo == INPUT_STEREO_SNGL_DOWN /* 6 */)
            {
                char cStereo;
                switch (bond_stereo)
                {
                case INPUT_STEREO_SNGL_UP:
                    cStereo = STEREO_SNGL_UP;
                    break;
                case INPUT_STEREO_SNGL_EITHER:
                    cStereo = STEREO_SNGL_EITHER;
                    break;
                case INPUT_STEREO_SNGL_DOWN:
                    cStereo = STEREO_SNGL_DOWN;
                    break;
                }
                at[a1].bond_stereo[n1] = cStereo;  /*  >0: the wedge (pointed) end is at this atom, a1 */
                at[a2].bond_stereo[n2] = -cStereo; /*  <0: the wedge (pointed) end is at the opposite atom, a1 */
            }
            else if (bond_stereo)
            {
                *err |= 16; /*  Ignored unrecognized Bond stereo */
                TREAT_ERR(*err, 0, "Unrecognized bond stereo");
                continue;
            }
        }

        /*(@nnuk : Nauman Ullah Khan) */
        LOG_NO_ARGS("\n################ (L862:mol2atom.c) ##################\n");
        LOG_MULT_ARGS("Bond %d: atom1=%d, atom2=%d, type=%d, stereo=%d\n", i, a1, a2, bond_type, bond_stereo);
        LOG_NO_ARGS("#######################################################\n");

    } /* eof copy bond info */

    *num_bonds = bonds;

    /* special valences */
    calculate_valences(mfdata, at, num_atoms, bDoNotAddH, err, pStrErr);

exit_function:;

    return at;
}
    */
    // END INCHI C FUNCTION: MakeInpAtomsFromMolfileData
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MakeInpAtomsFromMolfileData
    // INCHI✔️❌: #define SINGLET_IS_TRIPLET 1; the singlet-to-triplet branch is active.
    // INCHI✔️❌: #define bRELEASE_VERSION 1; the digit-prefixed debug-element branch is inactive.
    // INCHI✔️❌: NUM_ISO_H and NUMH are active header macros and are reproduced inline below.
    // INCHI✔️❌: LOG_NO_ARGS and LOG_MULT_ARGS compile to production logging macros with no chemical-state effect.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MakeInpAtomsFromMolfileData

    macro_rules! add_error {
        ($message:expr) => {{
            let text = format!("{}", $message);
            let message = text
                .bytes()
                .map(|byte| byte as i8)
                .chain([0])
                .collect::<Vec<_>>();
            AddErrorMessage(p_str_err.as_deref_mut(), Some(&message))?;
        }};
    }

    *err = 0;
    *num_atoms = mfdata.ctab.n_atoms;
    *num_bonds = 0;

    if MolfileHasNoChemStruc(Some(mfdata)) != 0 {
        return Ok(SourceMutPointer::null());
    }

    let atom_count = usize::try_from(*num_atoms)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let source_atoms = heap
        .slice(mfdata.ctab.atoms.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let bond_count = if mfdata.ctab.n_bonds > 0 {
        usize::try_from(mfdata.ctab.n_bonds)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    let source_bonds = if bond_count == 0 {
        Vec::new()
    } else {
        heap.slice(mfdata.ctab.bonds.as_const())?
            .get(..bond_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };

    let at = if at_inp.is_null() {
        let allocated = CreateInpAtom(heap, *num_atoms)?;
        if allocated.is_null() {
            *err = -1;
            add_error!("Out of RAM");
            return Ok(SourceMutPointer::null());
        }
        allocated
    } else {
        at_inp
    };

    let mut atoms = heap
        .slice(at.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let mut fatal_atom_error = false;

    for (index, source) in source_atoms.iter().enumerate() {
        mystrncpy_slice(Some(&mut atoms[index].elname), Some(&source.symbol), 6)?;
        atoms[index].orig_at_number = index.wrapping_add(1) as u16;
        atoms[index].iso_atw_diff = source.mass_difference;
        atoms[index].charge = source.charge;
        atoms[index].radical = source.radical;

        let isotope_difference = i32::from(source.mass_difference);
        atoms[index].iso_atw_diff = if isotope_difference == ZERO_ATW_DIFF as i32 {
            1
        } else if isotope_difference > 0 {
            isotope_difference.wrapping_add(1) as i8
        } else {
            isotope_difference as i8
        };

        if atoms[index].radical == RADICAL_SINGLET as i8 {
            atoms[index].radical = RADICAL_TRIPLET as i8;
        }

        let mut element_number = get_periodic_table_number(Some(&atoms[index].elname))?;
        if element_number == ERR_ELEM {
            let atom = &mut atoms[index];
            atom.num_H =
                extract_h_atoms(Some(&mut atom.elname), Some(&mut atom.num_iso_H))? as i8;
            let isotope_hydrogens = atoms[index]
                .num_iso_H
                .iter()
                .fold(0_i32, |sum, value| sum.wrapping_add(i32::from(*value)));
            if atoms[index].elname[0] == 0
                && i32::from(atoms[index].num_H).wrapping_add(isotope_hydrogens) != 0
            {
                atoms[index].elname = [b'H' as i8, 0, 0, 0, 0, 0];
                if isotope_hydrogens != 0 {
                    for isotope in (0..NUM_H_ISOTOPES as usize).rev() {
                        if atoms[index].num_iso_H[isotope] != 0 {
                            atoms[index].num_iso_H[isotope] =
                                atoms[index].num_iso_H[isotope].wrapping_sub(1);
                            atoms[index].iso_atw_diff = isotope.wrapping_add(1) as i8;
                            break;
                        }
                    }
                } else {
                    atoms[index].num_H = atoms[index].num_H.wrapping_sub(1);
                }
            }
            element_number = get_periodic_table_number(Some(&atoms[index].elname))?;
            if element_number == ERR_ELEM {
                element_number = 0;
            }
        }

        atoms[index].el_number = element_number as u8;
        if element_number == 0 {
            *err = -2;
            add_error!("Unknown element(s):");
            let element = atoms[index]
                .elname
                .iter()
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8 as char)
                .collect::<String>();
            add_error!(element);
            fatal_atom_error = true;
            break;
        }

        if element_number == 1 && atoms[index].iso_atw_diff == 0 {
            match atoms[index].elname[0] as u8 {
                b'D' => {
                    atoms[index].iso_atw_diff = 2;
                    mystrncpy_slice(
                        Some(&mut atoms[index].elname),
                        Some(&[b'H' as i8, 0]),
                        6,
                    )?;
                }
                b'T' => {
                    atoms[index].iso_atw_diff = 3;
                    mystrncpy_slice(
                        Some(&mut atoms[index].elname),
                        Some(&[b'H' as i8, 0]),
                        6,
                    )?;
                }
                _ => {}
            }
        }
    }

    if !fatal_atom_error {
        let mut bonds = 0_i32;
        for source in &source_bonds {
            let bond_stereo = source.bond_stereo;
            let mut bond_type = source.bond_type;
            let a1 = i32::from(source.atnum1).wrapping_sub(1);
            let a2 = i32::from(source.atnum2).wrapping_sub(1);

            if a1 < 0 || a1 >= *num_atoms || a2 < 0 || a2 >= *num_atoms || a1 == a2 {
                *err |= 1;
                add_error!("Bond to nonexistent atom");
                continue;
            }
            let a1 = a1 as usize;
            let a2 = a2 as usize;
            let position1 = is_in_the_list(
                Some(&atoms[a1].neighbor),
                a2 as u16,
                i32::from(atoms[a1].valence),
            )?;
            let position2 = is_in_the_list(
                Some(&atoms[a2].neighbor),
                a1 as u16,
                i32::from(atoms[a2].valence),
            )?;

            let (n1, n2) = if (position1.is_some() || position2.is_some())
                && (position1.is_some() || i32::from(atoms[a1].valence) < MAXVAL as i32)
                && (position2.is_some() || i32::from(atoms[a2].valence) < MAXVAL as i32)
            {
                let n1 = position1.unwrap_or_else(|| {
                    let value = atoms[a1].valence as usize;
                    atoms[a1].valence = atoms[a1].valence.wrapping_add(1);
                    value
                });
                let n2 = position2.unwrap_or_else(|| {
                    let value = atoms[a2].valence as usize;
                    atoms[a2].valence = atoms[a2].valence.wrapping_add(1);
                    value
                });
                add_error!("Multiple bonds between two atoms");
                *err |= 2;
                (n1, n2)
            } else if position1.is_none()
                && position2.is_none()
                && i32::from(atoms[a1].valence) < MAXVAL as i32
                && i32::from(atoms[a2].valence) < MAXVAL as i32
            {
                let n1 = atoms[a1].valence as usize;
                let n2 = atoms[a2].valence as usize;
                atoms[a1].valence = atoms[a1].valence.wrapping_add(1);
                atoms[a2].valence = atoms[a2].valence.wrapping_add(1);
                bonds = bonds.wrapping_add(1);
                (n1, n2)
            } else {
                *err |= 4;
                let element = if i32::from(atoms[a1].valence) >= MAXVAL as i32 {
                    &atoms[a1].elname
                } else {
                    &atoms[a2].elname
                };
                let element = element
                    .iter()
                    .take_while(|byte| **byte != 0)
                    .map(|byte| *byte as u8 as char)
                    .collect::<String>();
                add_error!(format!("Atom '{element}' has more than {MAXVAL} bonds"));
                continue;
            };

            if bond_type < MIN_INPUT_BOND_TYPE as i8 || bond_type > MAX_INPUT_BOND_TYPE as i8 {
                let original_type = bond_type;
                bond_type = BOND_TYPE_SINGLE as i8;
                add_error!("Unrecognized bond type:");
                add_error!(original_type);
                *err |= 8;
            }

            if n1 < MAXVAL as usize && n2 < MAXVAL as usize {
                atoms[a1].bond_type[n1] = bond_type as u8;
                atoms[a2].bond_type[n2] = bond_type as u8;
                atoms[a1].neighbor[n1] = a2 as u16;
                atoms[a2].neighbor[n2] = a1 as u16;

                if bond_stereo == INPUT_STEREO_DBLE_EITHER as i8 {
                    atoms[a1].bond_stereo[n1] = STEREO_DBLE_EITHER as i8;
                    atoms[a2].bond_stereo[n2] = STEREO_DBLE_EITHER as i8;
                } else if bond_stereo == INPUT_STEREO_SNGL_UP as i8
                    || bond_stereo == INPUT_STEREO_SNGL_EITHER as i8
                    || bond_stereo == INPUT_STEREO_SNGL_DOWN as i8
                {
                    let stereo = match bond_stereo {
                        value if value == INPUT_STEREO_SNGL_UP as i8 => STEREO_SNGL_UP as i8,
                        value if value == INPUT_STEREO_SNGL_EITHER as i8 => {
                            STEREO_SNGL_EITHER as i8
                        }
                        _ => STEREO_SNGL_DOWN as i8,
                    };
                    atoms[a1].bond_stereo[n1] = stereo;
                    atoms[a2].bond_stereo[n2] = stereo.wrapping_neg();
                } else if bond_stereo != 0 {
                    *err |= 16;
                    add_error!("Unrecognized bond stereo");
                    continue;
                }
            }
        }
        *num_bonds = bonds;
        calculate_valences(
            heap,
            Some(mfdata),
            &mut atoms,
            *num_atoms,
            b_do_not_add_h,
            err,
            p_str_err.as_deref_mut(),
        )?;
    }

    heap.slice_mut(at)?[..atom_count].clone_from_slice(&atoms);
    Ok(at)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn calculate_valences(
    heap: &SourceHeap,
    mfdata: Option<&MOL_FMT_DATA>,
    at: &mut [inp_ATOM],
    num_atoms: i32,
    b_do_not_add_h: i32,
    err: &mut i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:807 calculate_valences
    // INCHI✔️❌: complete source frame follows verbatim; checked Rust indexing and message materialization add overhead.
    /*
void calculate_valences(MOL_FMT_DATA *mfdata,
                        inp_ATOM *at,
                        int *num_atoms,
                        int bDoNotAddH,
                        int *err,
                        char *pStrErr)
{
    int bNonMetal;
    int a1, a2, n1, n2, valence;
    AT_NUMB *p1;

    /* special valences */

    for (bNonMetal = 0; bNonMetal < 2; bNonMetal++)
    {
        for (a1 = 0; a1 < *num_atoms; a1++)
        {
            int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1],
                bond_type,
                bHasMetalNeighbor;
            /* should the "!=" be replaced with "==" ??? */
            if (bNonMetal == is_el_a_metal(at[a1].el_number))
            {
                /* first process all metals, after that all non-metals */
                continue;
            }

            memset(num_bond_type, 0, sizeof(num_bond_type)); /* djb-rwth: memset_s C11/Annex K variant? */

            /* valence = at[a1].chem_bonds_valence; */ /*  save atom valence if available */
                                                       /* 2006-08-31: fix for uncharged >N(IV)- in an aromatic ring */

            valence =
                (mfdata && mfdata->ctab.atoms) ? mfdata->ctab.atoms[a1].valence
                                               : at[a1].chem_bonds_valence;

            at[a1].chem_bonds_valence = 0;
            bHasMetalNeighbor = 0;
            for (n1 = 0; n1 < at[a1].valence; n1++)
            {
                bond_type = at[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
                if (bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE)
                {
                    bond_type = 0;
                    TREAT_ERR(*err, 0, "Unknown bond type in MOLfile assigned as a single bond");
                }
                num_bond_type[bond_type]++;
                /* -- too a radical solution -- removed from next to ver 1.12B --- */
            }

            for (n1 = 0;
                 MIN_INPUT_BOND_TYPE + n1 <= 3 &&
                 MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE;
                 n1++)
            {
                /* add all bond orders except for "aromatic" bonds */
                at[a1].chem_bonds_valence += (MIN_INPUT_BOND_TYPE + n1) * num_bond_type[n1];
            }

            n2 = 0; /* djb-rwth: ignoring LLVM warning: value used */
            if (MIN_INPUT_BOND_TYPE <= BOND_TYPE_ALTERN &&
                BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
                (n2 = num_bond_type[BOND_TYPE_ALTERN - MIN_INPUT_BOND_TYPE]))
            {
                /* accept input aromatic bonds for now */
                switch (n2)
                {
                case 2:
                    at[a1].chem_bonds_valence += 3; /* =A- */
                    break;
                case 3:
                    at[a1].chem_bonds_valence += 4; /* =A< */
                    break;
                default:
                    /*  if 1 or >= 4 aromatic bonds then replace   */
                    /* such bonds with single bonds                */
                    /*  and detect an error in the input structure */
                    for (n1 = 0; n1 < at[a1].valence; n1++)
                    {
                        if (at[a1].bond_type[n1] == BOND_TYPE_ALTERN)
                        {
                            a2 = at[a1].neighbor[n1];
                            p1 = is_in_the_list(at[a2].neighbor, (AT_NUMB)a1,
                                                at[a2].valence);
                            if (p1)
                            {
                                at[a1].bond_type[n1] =
                                    at[a2].bond_type[p1 - at[a2].neighbor] = BOND_TYPE_SINGLE;
                            }
                            else
                            {
                                *err = -2; /*  Program error */
                                TREAT_ERR(*err, 0, "Program error interpreting MOLfile");
                                return; /*  no structure */
                            }
                        }
                    }

                    at[a1].chem_bonds_valence += n2;
                    *err |= 32;
                    TREAT_ERR(*err, 0, "Atom has 1 or more than 3 aromatic bonds");
                    n2 = 0;
                    break;
                }
            }

            if (n2 && !valence)
            {
                /* atom has aromatic bonds AND the chemical valence is not known */

                int num_H = NUMH(at, a1);
                /* bug fix 2006-08-25: aliased H result in num_H > 0 */
                /* => wrong call to detect_unusual_el_valence() */
                int chem_valence = at[a1].chem_bonds_valence /*+ num_H*/;

                int bUnusualValenceArom =
                    detect_unusual_el_valence((int)at[a1].el_number, at[a1].charge,
                                              at[a1].radical, chem_valence,
                                              num_H, at[a1].valence);
                int bUnusualValenceNoArom =
                    detect_unusual_el_valence((int)at[a1].el_number, at[a1].charge,
                                              at[a1].radical, chem_valence - 1,
                                              num_H, at[a1].valence);

#if (CHECK_AROMBOND2ALT == 1)
                if (bUnusualValenceArom &&
                    !bUnusualValenceNoArom &&
                    0 == nBondsValToMetal(at, a1))
#else
                if (bUnusualValenceArom && !bUnusualValenceNoArom)
#endif
                {
                    /* typically NH in 5-member aromatic ring */
                    at[a1].chem_bonds_valence--;
                }
            }
            else if (n2 && valence)
            {
                /* atom has aromatic bonds AND the chemical valence is known */
                int num_H = NUMH(at, a1);
                int chem_valence = at[a1].chem_bonds_valence + num_H;
                if (valence == chem_valence - 1)
                {
                    /* typically NH in 5-member aromatic ring */
                    at[a1].chem_bonds_valence--;
                }
            }

            /*  Set number of hydrogen atoms */
            if (mfdata)
            {
                at[a1].num_H = get_num_H(at[a1].elname,
                                         at[a1].num_H,
                                         at[a1].num_iso_H,
                                         at[a1].charge, at[a1].radical,
                                         at[a1].chem_bonds_valence,
                                         mfdata->ctab.atoms[a1].valence, /* instead of valence */
                                         mfdata->ctab.atoms[a1].atom_aliased_flag,
                                         bDoNotAddH,
                                         bHasMetalNeighbor);
            }
        }
    } /* for ( bNonMetal = ... */

    return;
}
    */
    // END INCHI C FUNCTION: calculate_valences
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: calculate_valences
    // INCHI✔️❌: #define CHECK_AROMBOND2ALT          1
    // INCHI✔️❌: #define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
    // INCHI✔️❌: #define NUMH(AT, N) (AT[N].num_H + NUM_ISO_H(AT, N))
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: calculate_valences

    macro_rules! add_error {
        ($message:literal) => {{
            let message = $message
                .bytes()
                .map(|byte| byte as i8)
                .chain([0])
                .collect::<Vec<_>>();
            AddErrorMessage(p_str_err.as_deref_mut(), Some(&message))?;
        }};
    }

    let atom_count = if num_atoms <= 0 {
        0
    } else {
        num_atoms as usize
    };
    if atom_count > at.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    for b_non_metal in 0..2 {
        for a1 in 0..atom_count {
            if b_non_metal == is_el_a_metal(i32::from(at[a1].el_number))? {
                continue;
            }

            let mut num_bond_type = [0_i32; (MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1) as usize];
            let source_valence = if let Some(data) = mfdata {
                if data.ctab.atoms.is_null() {
                    i32::from(at[a1].chem_bonds_valence)
                } else {
                    i32::from(heap.slice(data.ctab.atoms.as_const())?[a1].valence)
                }
            } else {
                i32::from(at[a1].chem_bonds_valence)
            };

            at[a1].chem_bonds_valence = 0;
            let valence = i32::from(at[a1].valence);
            for n1 in 0..valence {
                let mut bond_type = i32::from(at[a1].bond_type[n1 as usize])
                    .wrapping_sub(MIN_INPUT_BOND_TYPE as i32);
                if bond_type < 0
                    || bond_type > (MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE) as i32
                {
                    bond_type = 0;
                    add_error!("Unknown bond type in MOLfile assigned as a single bond");
                }
                num_bond_type[bond_type as usize] =
                    num_bond_type[bond_type as usize].wrapping_add(1);
            }

            let mut n1 = 0_i32;
            while MIN_INPUT_BOND_TYPE as i32 + n1 <= 3
                && MIN_INPUT_BOND_TYPE as i32 + n1 <= MAX_INPUT_BOND_TYPE as i32
            {
                let contribution = (MIN_INPUT_BOND_TYPE as i32 + n1)
                    .wrapping_mul(num_bond_type[n1 as usize]);
                at[a1].chem_bonds_valence =
                    at[a1].chem_bonds_valence.wrapping_add(contribution as i8);
                n1 = n1.wrapping_add(1);
            }

            let mut n2 = 0_i32;
            if MIN_INPUT_BOND_TYPE <= BOND_TYPE_ALTERN
                && BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE
            {
                n2 = num_bond_type[(BOND_TYPE_ALTERN - MIN_INPUT_BOND_TYPE) as usize];
            }
            if n2 != 0 {
                match n2 {
                    2 => {
                        at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_add(3);
                    }
                    3 => {
                        at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_add(4);
                    }
                    _ => {
                        for n1 in 0..valence {
                            if u32::from(at[a1].bond_type[n1 as usize]) == BOND_TYPE_ALTERN {
                                let a2 = usize::from(at[a1].neighbor[n1 as usize]);
                                let reciprocal = {
                                    let neighbor = at
                                        .get(a2)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    is_in_the_list(
                                        Some(&neighbor.neighbor),
                                        a1 as u16,
                                        i32::from(neighbor.valence),
                                    )?
                                };
                                let Some(reciprocal) = reciprocal else {
                                    *err = -2;
                                    add_error!("Program error interpreting MOLfile");
                                    return Ok(());
                                };
                                at[a1].bond_type[n1 as usize] = BOND_TYPE_SINGLE as u8;
                                at[a2].bond_type[reciprocal] = BOND_TYPE_SINGLE as u8;
                            }
                        }
                        at[a1].chem_bonds_valence =
                            at[a1].chem_bonds_valence.wrapping_add(n2 as i8);
                        *err |= 32;
                        add_error!("Atom has 1 or more than 3 aromatic bonds");
                        n2 = 0;
                    }
                }
            }

            if n2 != 0 && source_valence == 0 {
                let atom = &at[a1];
                let num_h = i32::from(atom.num_H)
                    .wrapping_add(atom.num_iso_H.iter().copied().map(i32::from).sum::<i32>());
                let chemical_valence = i32::from(atom.chem_bonds_valence);
                let unusual_aromatic = detect_unusual_el_valence(
                    i32::from(atom.el_number),
                    i32::from(atom.charge),
                    i32::from(atom.radical),
                    chemical_valence,
                    num_h,
                    i32::from(atom.valence),
                )?;
                let unusual_non_aromatic = detect_unusual_el_valence(
                    i32::from(atom.el_number),
                    i32::from(atom.charge),
                    i32::from(atom.radical),
                    chemical_valence.wrapping_sub(1),
                    num_h,
                    i32::from(atom.valence),
                )?;
                if unusual_aromatic != 0
                    && unusual_non_aromatic == 0
                    && n_bonds_val_to_metal(Some(at), a1 as i32)? == 0
                {
                    at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_sub(1);
                }
            } else if n2 != 0 && source_valence != 0 {
                let atom = &at[a1];
                let num_h = i32::from(atom.num_H)
                    .wrapping_add(atom.num_iso_H.iter().copied().map(i32::from).sum::<i32>());
                let chemical_valence = i32::from(atom.chem_bonds_valence).wrapping_add(num_h);
                if source_valence == chemical_valence.wrapping_sub(1) {
                    at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_sub(1);
                }
            }

            if let Some(data) = mfdata {
                let mol_atom = heap.slice(data.ctab.atoms.as_const())?[a1].clone();
                let atom = &at[a1];
                let hydrogen_count = get_num_H(
                    Some(&atom.elname),
                    i32::from(atom.num_H),
                    Some(&atom.num_iso_H),
                    i32::from(atom.charge),
                    i32::from(atom.radical),
                    i32::from(atom.chem_bonds_valence),
                    i32::from(mol_atom.valence),
                    i32::from(mol_atom.atom_aliased_flag),
                    b_do_not_add_h,
                    0,
                )?;
                at[a1].num_H = hydrogen_count as i8;
            }
        }
    }

    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn SetInpAtomsXYZ(
    heap: &mut SourceHeap,
    mfdata: &MOL_FMT_DATA,
    num_atoms: i32,
    at: SourceMutPointer<inp_ATOM>,
    err: &mut i32,
    p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:975 SetInpAtomsXYZ
    // INCHI✔️❌: complete source frame follows verbatim; checked heap access and source-coordinate cloning add overhead.
    /*
int SetInpAtomsXYZ(MOL_FMT_DATA *mfdata,
                   int num_atoms,
                   inp_ATOM *at,
                   int *err,
                   char *pStrErr)
{
    int i, num_dimensions = 0;

#if (NORMALIZE_INP_COORD == 1)
    int do_scale_xyz = 1;
#else
    int do_scale_xyz = 0;
#endif
    double x0, y0, z0, xmin, ymin, zmin, scaler;

    num_dimensions = MolfileGetXYZDimAndNormFactors(mfdata, do_scale_xyz,
                                                    &x0, &y0, &z0,
                                                    &xmin, &ymin, &zmin,
                                                    &scaler, err, pStrErr);

    if (num_dimensions == 0)
    {
        goto exit_function;
    }

    for (i = 0; i < num_atoms; i++)
    {

        double x = mfdata->ctab.atoms[i].fx;
        double y = mfdata->ctab.atoms[i].fy;
        double z = mfdata->ctab.atoms[i].fz;

        if (!do_scale_xyz)
        {
            at[i].x = x;
            at[i].y = y;
            at[i].z = z;
        }
        else
        {
            x = (x - xmin) * scaler + x0;
            y = (y - ymin) * scaler + y0;
            z = (z - zmin) * scaler + z0;
            /* floor() behavior is not well defined for negative arguments.
             * Use positive arguments only to get nearest integer.
             */
            at[i].x = (x >= 0.0) ? (int)floor(x + 0.5) : -(int)floor(-x + 0.5);
            at[i].y = (y >= 0.0) ? (int)floor(y + 0.5) : -(int)floor(-y + 0.5);
            at[i].z = (z >= 0.0) ? (int)floor(z + 0.5) : -(int)floor(-z + 0.5);
        }
    }

exit_function:;

    return num_dimensions;
}
    */
    // END INCHI C FUNCTION: SetInpAtomsXYZ
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: SetInpAtomsXYZ
    // INCHI✔️❌: #define NORMALIZE_INP_COORD 0; do_scale_xyz is 0 and the integer normalization branch is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: SetInpAtomsXYZ

    let (mut x0, mut y0, mut z0) = (0.0, 0.0, 0.0);
    let (mut xmin, mut ymin, mut zmin) = (0.0, 0.0, 0.0);
    let mut scaler = 0.0;
    let num_dimensions = MolfileGetXYZDimAndNormFactors(
        heap,
        Some(mfdata),
        0,
        &mut x0,
        &mut y0,
        &mut z0,
        &mut xmin,
        &mut ymin,
        &mut zmin,
        &mut scaler,
        err,
        p_str_err,
    )?;
    if num_dimensions == 0 {
        return Ok(0);
    }

    let count = if num_atoms > 0 {
        usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    let coordinates = heap
        .slice(mfdata.ctab.atoms.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .iter()
        .map(|atom| (atom.fx, atom.fy, atom.fz))
        .collect::<Vec<_>>();
    let output = heap
        .slice_mut(at)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for (atom, &(x, y, z)) in output.iter_mut().zip(&coordinates) {
        atom.x = x;
        atom.y = y;
        atom.z = z;
    }
    Ok(num_dimensions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        INT_ARRAY, MOL_FMT_ATOM, MOL_FMT_BOND, MOL_FMT_SGROUP, MOL_FMT_SGROUPS, MOL_FMT_v3000,
        NUM_LISTS,
    };

    fn assert_missing<T: 'static>(heap: &SourceHeap, pointer: SourceMutPointer<T>) {
        assert_eq!(
            heap.slice(pointer.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__mol2atom__freeinpatom__line_1046() {
        let mut heap = SourceHeap::default();

        FreeInpAtom(&mut heap, None).unwrap();

        let mut null_slot = SourceMutPointer::null();
        FreeInpAtom(&mut heap, Some(&mut null_slot)).unwrap();
        assert!(null_slot.is_null());

        let atom = heap.allocate(vec![inp_ATOM::default(); 2]).unwrap();
        let stale_alias = atom;
        let mut atom_slot = atom;
        FreeInpAtom(&mut heap, Some(&mut atom_slot)).unwrap();
        assert!(atom_slot.is_null());
        assert_missing(&heap, stale_alias);

        let mut stale_slot = stale_alias;
        assert_eq!(
            FreeInpAtom(&mut heap, Some(&mut stale_slot)),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(stale_slot, stale_alias);
    }

    #[test]
    fn source_port__mol2atom__createinpatom__line_1033() {
        let mut heap = SourceHeap::default();
        assert!(CreateInpAtom(&mut heap, -1).unwrap().is_null());
        assert!(CreateInpAtom(&mut heap, i32::MIN).unwrap().is_null());

        let empty = CreateInpAtom(&mut heap, 0).unwrap();
        assert_eq!(heap.slice(empty.as_const()).unwrap(), &[]);

        let atoms = CreateInpAtom(&mut heap, 3).unwrap();
        assert_eq!(
            heap.slice(atoms.as_const()).unwrap(),
            &[
                inp_ATOM::default(),
                inp_ATOM::default(),
                inp_ATOM::default()
            ]
        );
        let mut atom_slot = atoms;
        FreeInpAtom(&mut heap, Some(&mut atom_slot)).unwrap();
        assert!(atom_slot.is_null());
        assert_missing(&heap, atoms);
    }

    #[test]
    fn source_port__mol2atom__createinpatomdata__line_1071() {
        let mut heap = SourceHeap::default();
        let old_atoms = heap.allocate(vec![inp_ATOM::default()]).unwrap();
        let old_fixed = heap.allocate(vec![inp_ATOM::default()]).unwrap();
        let mut input = INP_ATOM_DATA {
            at: old_atoms,
            at_fixed_bonds: old_fixed,
            num_at: 91,
            num_bonds: 92,
            bExists: 1,
            ..INP_ATOM_DATA::default()
        };
        assert_eq!(CreateInpAtomData(&mut heap, &mut input, -1, 1), Ok(0));
        assert_eq!(input, INP_ATOM_DATA::default());
        assert_missing(&heap, old_atoms);
        assert_missing(&heap, old_fixed);

        assert_eq!(CreateInpAtomData(&mut heap, &mut input, 2, 0), Ok(1));
        assert_eq!(input.num_at, 2);
        assert!(!input.at.is_null());
        assert!(input.at_fixed_bonds.is_null());
        assert_eq!(
            heap.slice(input.at.as_const()).unwrap(),
            &[inp_ATOM::default(), inp_ATOM::default()]
        );

        assert_eq!(CreateInpAtomData(&mut heap, &mut input, 3, -7), Ok(1));
        assert_eq!(input.num_at, 3);
        assert_eq!(heap.slice(input.at.as_const()).unwrap().len(), 3);
        assert_eq!(
            heap.slice(input.at_fixed_bonds.as_const()).unwrap().len(),
            3
        );

        let mut first_failure_heap = SourceHeap::default();
        first_failure_heap.fail_after_allocations(0);
        let mut first_failure = INP_ATOM_DATA {
            num_at: 8,
            ..INP_ATOM_DATA::default()
        };
        assert_eq!(
            CreateInpAtomData(&mut first_failure_heap, &mut first_failure, 2, 1),
            Ok(0)
        );
        assert_eq!(first_failure, INP_ATOM_DATA::default());
        assert_eq!(first_failure_heap.source_allocation_calls(), 1);

        let mut second_failure_heap = SourceHeap::default();
        second_failure_heap.fail_after_allocations(1);
        let mut second_failure = INP_ATOM_DATA::default();
        assert_eq!(
            CreateInpAtomData(&mut second_failure_heap, &mut second_failure, 2, 1),
            Ok(0)
        );
        assert_eq!(second_failure, INP_ATOM_DATA::default());
        assert_eq!(second_failure_heap.source_allocation_calls(), 2);
    }

    #[test]
    fn source_port__mol2atom__freecompatomdata__line_1090() {
        let mut heap = SourceHeap::default();
        let mut empty = COMP_ATOM_DATA {
            num_at: 7,
            num_removed_H: 8,
            num_bonds: 9,
            num_isotopic: 10,
            bExists: 1,
            bDeleted: 2,
            bHasIsotopicLayer: 3,
            bTautomeric: 4,
            nNumRemovedProtons: 5,
            nNumRemovedProtonsIsotopic: [6, 7, 8],
            num_iso_H: [9, 10, 11],
            num_components: 12,
            ..COMP_ATOM_DATA::default()
        };
        assert_eq!(FreeCompAtomData(&mut heap, &mut empty), Ok(()));
        assert_eq!(empty, COMP_ATOM_DATA::default());

        let atoms = heap.allocate(vec![inp_ATOM::default(); 3]).unwrap();
        let offsets = heap.allocate(vec![1_u16, 3, 5, 7]).unwrap();
        let mut full = COMP_ATOM_DATA {
            at: atoms,
            num_at: 3,
            num_removed_H: 2,
            num_bonds: 1,
            num_isotopic: 4,
            bExists: 1,
            bDeleted: 1,
            bHasIsotopicLayer: 1,
            bTautomeric: 1,
            nNumRemovedProtons: -2,
            nNumRemovedProtonsIsotopic: [-1, 2, 3],
            num_iso_H: [4, 5, 6],
            nOffsetAtAndH: offsets,
            num_components: 2,
        };
        assert_eq!(FreeCompAtomData(&mut heap, &mut full), Ok(()));
        assert_eq!(full, COMP_ATOM_DATA::default());
        assert_missing(&heap, atoms);
        assert_missing(&heap, offsets);

        assert_eq!(FreeCompAtomData(&mut heap, &mut full), Ok(()));
        assert_eq!(full, COMP_ATOM_DATA::default());
    }

    #[test]
    fn source_port__mol2atom__freeorigatdata__line_1224() {
        let mut heap = SourceHeap::default();
        FreeOrigAtData(&mut heap, None).unwrap();

        let mut empty = ORIG_ATOM_DATA {
            num_dimensions: 3,
            num_inp_atoms: 17,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        };
        FreeOrigAtData(&mut heap, Some(&mut empty)).unwrap();
        assert_eq!(empty, ORIG_ATOM_DATA::default());

        let atom = heap.allocate(vec![inp_ATOM::default(); 2]).unwrap();
        let cur = heap.allocate(vec![1_u16, 2]).unwrap();
        let old = heap.allocate(vec![3_u16, 4]).unwrap();
        let coord = heap.allocate(vec![Default::default(); 2]).unwrap();
        let equ = heap.allocate(vec![5_u16, 6]).unwrap();
        let sorted = heap.allocate(vec![7_u16, 8]).unwrap();
        let polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        let v3000 = heap.allocate(vec![OAD_V3000::default()]).unwrap();
        let mut full = ORIG_ATOM_DATA {
            at: atom,
            num_dimensions: 2,
            num_inp_bonds: 1,
            num_inp_atoms: 2,
            num_components: 1,
            nCurAtLen: cur,
            nOldCompNumber: old,
            nEquLabels: equ,
            nSortedOrder: sorted,
            szCoord: coord,
            polymer,
            v3000,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        };
        FreeOrigAtData(&mut heap, Some(&mut full)).unwrap();
        assert_eq!(full, ORIG_ATOM_DATA::default());
        assert_missing(&heap, atom);
        assert_missing(&heap, cur);
        assert_missing(&heap, old);
        assert_missing(&heap, coord);
        assert_missing(&heap, equ);
        assert_missing(&heap, sorted);
        assert_missing(&heap, polymer);
        assert_missing(&heap, v3000);
    }

    #[test]
    fn source_port__mol2atom__freeextorigatdata__line_1269() {
        let mut heap = SourceHeap::default();
        FreeExtOrigAtData(
            &mut heap,
            SourceMutPointer::null(),
            SourceMutPointer::null(),
        )
        .unwrap();

        let polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        FreeExtOrigAtData(&mut heap, polymer, SourceMutPointer::null()).unwrap();
        assert_missing(&heap, polymer);

        let orig = heap.allocate(vec![1_i32, 2]).unwrap();
        let fin = heap.allocate(vec![2_i32, 1]).unwrap();
        let mut rows = Vec::new();
        let mut lists = Vec::new();
        for base in [10_i32, 20, 30, 40] {
            let first = heap.allocate(vec![base, base + 1]).unwrap();
            let outer = heap
                .allocate(vec![first, SourceMutPointer::null()])
                .unwrap();
            rows.push(first);
            lists.push(outer);
        }
        let full = heap
            .allocate(vec![OAD_V3000 {
                atom_index_orig: orig,
                atom_index_fin: fin,
                n_haptic_bonds: 2,
                lists_haptic_bonds: lists[0],
                n_steabs: 2,
                lists_steabs: lists[1],
                n_sterel: 2,
                lists_sterel: lists[2],
                n_sterac: 2,
                lists_sterac: lists[3],
                ..OAD_V3000::default()
            }])
            .unwrap();
        FreeExtOrigAtData(&mut heap, SourceMutPointer::null(), full).unwrap();
        assert_missing(&heap, orig);
        assert_missing(&heap, fin);
        for pointer in rows {
            assert_missing(&heap, pointer);
        }
        for pointer in lists {
            assert_missing(&heap, pointer);
        }
        assert_missing(&heap, full);

        let mut zero_rows = Vec::new();
        let mut zero_lists = Vec::new();
        for base in [50_i32, 60, 70, 80] {
            let row = heap.allocate(vec![base]).unwrap();
            let outer = heap.allocate(vec![row]).unwrap();
            zero_rows.push(row);
            zero_lists.push(outer);
        }
        let zero = heap
            .allocate(vec![OAD_V3000 {
                lists_haptic_bonds: zero_lists[0],
                lists_steabs: zero_lists[1],
                lists_sterel: zero_lists[2],
                lists_sterac: zero_lists[3],
                ..OAD_V3000::default()
            }])
            .unwrap();
        FreeExtOrigAtData(&mut heap, SourceMutPointer::null(), zero).unwrap();
        for pointer in &zero_rows {
            assert!(heap.slice(pointer.as_const()).is_ok());
        }
        for pointer in &zero_lists {
            assert!(heap.slice(pointer.as_const()).is_ok());
        }
        for pointer in zero_rows {
            heap.free(pointer).unwrap();
        }
        for pointer in zero_lists {
            heap.free(pointer).unwrap();
        }

        let mut negative_rows = Vec::new();
        let mut negative_lists = Vec::new();
        for base in [90_i32, 100, 110, 120] {
            let row = heap.allocate(vec![base]).unwrap();
            let outer = heap.allocate(vec![row]).unwrap();
            negative_rows.push(row);
            negative_lists.push(outer);
        }
        let negative = heap
            .allocate(vec![OAD_V3000 {
                n_haptic_bonds: i32::MIN,
                lists_haptic_bonds: negative_lists[0],
                n_steabs: i32::MIN,
                lists_steabs: negative_lists[1],
                n_sterel: i32::MIN,
                lists_sterel: negative_lists[2],
                n_sterac: i32::MIN,
                lists_sterac: negative_lists[3],
                ..OAD_V3000::default()
            }])
            .unwrap();
        FreeExtOrigAtData(&mut heap, SourceMutPointer::null(), negative).unwrap();
        for pointer in negative_lists {
            assert_missing(&heap, pointer);
        }
        for pointer in negative_rows {
            assert!(heap.slice(pointer.as_const()).is_ok());
            heap.free(pointer).unwrap();
        }

        let null_lists = heap
            .allocate(vec![OAD_V3000 {
                n_haptic_bonds: i32::MAX,
                n_steabs: i32::MAX,
                n_sterel: i32::MAX,
                n_sterac: i32::MAX,
                ..OAD_V3000::default()
            }])
            .unwrap();
        FreeExtOrigAtData(&mut heap, SourceMutPointer::null(), null_lists).unwrap();
        assert_missing(&heap, null_lists);
    }

    #[test]
    fn source_port__mol2atom__setextorigatdatabymolfileextinput__line_1340() {
        fn symbol(value: &[u8]) -> [i8; 6] {
            let mut result = [0_i8; 6];
            for (destination, source) in result.iter_mut().zip(value.iter().copied()) {
                *destination = source as i8;
            }
            result
        }

        fn error_text(buffer: &[i8]) -> Vec<u8> {
            buffer
                .iter()
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8)
                .collect()
        }

        fn fixture(
            heap: &mut SourceHeap,
            bond_index: i32,
            first_symbol: &[u8],
            second_symbol: &[u8],
            include_v3000: bool,
        ) -> MOL_FMT_DATA {
            let alist = heap.allocate(vec![1_i32, 2]).unwrap();
            let blist = heap.allocate(vec![bond_index]).unwrap();
            let mut smt = [0_i8; 80];
            smt[..7].copy_from_slice(&[b'A' as i8, b'B' as i8, 0, 9, 8, 7, 6]);
            let group = heap
                .allocate(vec![MOL_FMT_SGROUP {
                    id: 11,
                    type_: 12,
                    subtype: 13,
                    conn: 14,
                    label: 15,
                    xbr1: [1.0, 2.0, 3.0, 4.0],
                    xbr2: [5.0, 6.0, 7.0, 8.0],
                    smt,
                    alist: INT_ARRAY {
                        item: alist,
                        allocated: 2,
                        used: 2,
                        increment: 1,
                    },
                    blist: INT_ARRAY {
                        item: blist,
                        allocated: 1,
                        used: 1,
                        increment: 1,
                    },
                }])
                .unwrap();
            let groups = heap.allocate(vec![group]).unwrap();
            let atoms = heap
                .allocate(vec![
                    MOL_FMT_ATOM {
                        symbol: symbol(first_symbol),
                        ..MOL_FMT_ATOM::default()
                    },
                    MOL_FMT_ATOM {
                        symbol: symbol(second_symbol),
                        ..MOL_FMT_ATOM::default()
                    },
                ])
                .unwrap();
            let bonds = heap
                .allocate(vec![MOL_FMT_BOND {
                    atnum1: 1,
                    atnum2: 2,
                    bond_type: 1,
                    bond_stereo: 0,
                }])
                .unwrap();

            let v3000 = if include_v3000 {
                let atom_index_orig = heap
                    .allocate(vec![0x1122_3344_i32, 0x5566_7788_i32])
                    .unwrap();
                let atom_index_fin = heap
                    .allocate(vec![0x1020_3040_i32, 0x5060_7080_i32])
                    .unwrap();

                let haptic_row = heap.allocate(vec![7_i32, 8, 2, 11, 12]).unwrap();
                let steabs_row = heap.allocate(vec![0_i32, 2, 21, 22]).unwrap();
                let sterac_row = heap.allocate(vec![31_i32, 2, 32, 33]).unwrap();
                let sterel_row = heap.allocate(vec![41_i32, 2, 42, 43]).unwrap();
                let haptic_lists = heap.allocate(vec![haptic_row]).unwrap();
                let steabs_lists = heap.allocate(vec![steabs_row]).unwrap();
                let sterac_lists = heap.allocate(vec![sterac_row]).unwrap();
                let sterel_lists = heap.allocate(vec![sterel_row]).unwrap();
                let haptic = heap
                    .allocate(vec![NUM_LISTS {
                        lists: haptic_lists,
                        allocated: 1,
                        used: 1,
                        increment: 1,
                    }])
                    .unwrap();
                let steabs = heap
                    .allocate(vec![NUM_LISTS {
                        lists: steabs_lists,
                        allocated: 1,
                        used: 1,
                        increment: 1,
                    }])
                    .unwrap();
                let sterac = heap
                    .allocate(vec![NUM_LISTS {
                        lists: sterac_lists,
                        allocated: 1,
                        used: 1,
                        increment: 1,
                    }])
                    .unwrap();
                let sterel = heap
                    .allocate(vec![NUM_LISTS {
                        lists: sterel_lists,
                        allocated: 1,
                        used: 1,
                        increment: 1,
                    }])
                    .unwrap();
                heap.allocate(vec![MOL_FMT_v3000 {
                    n_non_star_atoms: 51,
                    n_star_atoms: 52,
                    atom_index_orig,
                    atom_index_fin,
                    n_sgroups: 53,
                    n_3d_constraints: 54,
                    n_collections: 55,
                    n_non_haptic_bonds: 56,
                    n_haptic_bonds: 1,
                    haptic_bonds: haptic,
                    n_steabs: 1,
                    steabs,
                    n_sterel: 1,
                    sterel,
                    n_sterac: 1,
                    sterac,
                }])
                .unwrap()
            } else {
                SourceMutPointer::null()
            };

            MOL_FMT_DATA {
                ctab: crate::source_types::MOL_FMT_CTAB {
                    n_atoms: 5,
                    n_bonds: 1,
                    atoms,
                    bonds,
                    sgroups: MOL_FMT_SGROUPS {
                        group: groups,
                        allocated: 1,
                        used: 1,
                        increment: 1,
                    },
                    v3000,
                    ..crate::source_types::MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            }
        }

        let mut untouched_heap = SourceHeap::default();
        let old_polymer = untouched_heap
            .allocate(vec![OAD_Polymer::default()])
            .unwrap();
        let old_v3000 = untouched_heap.allocate(vec![OAD_V3000::default()]).unwrap();
        let mut polymer_slot = old_polymer;
        let mut v3000_slot = old_v3000;
        assert_eq!(
            SetExtOrigAtDataByMolfileExtInput(
                &mut untouched_heap,
                &MOL_FMT_DATA::default(),
                &mut polymer_slot,
                &mut v3000_slot,
                None,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, old_polymer);
        assert_eq!(v3000_slot, old_v3000);
        FreeExtOrigAtData(&mut untouched_heap, old_polymer, old_v3000).unwrap();

        let mut heap = SourceHeap::default();
        let input = fixture(&mut heap, 1, b"C", b"O", true);
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        let mut errors = [0_i8; 256];
        assert_eq!(
            SetExtOrigAtDataByMolfileExtInput(
                &mut heap,
                &input,
                &mut polymer_slot,
                &mut v3000_slot,
                Some(&mut errors),
            ),
            Ok(0)
        );
        assert!(error_text(&errors).is_empty());
        let polymer = &heap.slice(polymer_slot.as_const()).unwrap()[0];
        assert_eq!(
            (
                polymer.n,
                polymer.is_in_reconn,
                polymer.edit_repeats,
                polymer.really_do_frame_shift,
                polymer.frame_shift_scheme,
                polymer.treat,
            ),
            (
                1,
                0,
                0,
                0,
                tagFrameShifScheme_FSS_STARS_CYCLED as i32,
                POLYMERS_MODERN as i32,
            )
        );
        assert!(polymer.pzz.is_null());
        let unit_pointer = heap.slice(polymer.units.as_const()).unwrap()[0];
        let unit = &heap.slice(unit_pointer.as_const()).unwrap()[0];
        assert_eq!(
            (unit.id, unit.type_, unit.subtype, unit.conn, unit.label),
            (11, 12, 13, 14, 15)
        );
        assert_eq!(unit.xbr1, [1.0, 2.0, 3.0, 4.0]);
        assert_eq!(unit.xbr2, [5.0, 6.0, 7.0, 8.0]);
        assert_eq!(&unit.smt[..7], &[b'A' as i8, b'B' as i8, 0, 0, 0, 0, 0]);
        assert_eq!(heap.slice(unit.alist.as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(unit.blist.as_const()).unwrap(), &[1, 2]);

        let output = &heap.slice(v3000_slot.as_const()).unwrap()[0];
        assert_eq!(
            (
                output.n_non_star_atoms,
                output.n_star_atoms,
                output.n_sgroups,
                output.n_3d_constraints,
                output.n_collections,
                output.n_non_haptic_bonds,
            ),
            (51, 52, 53, 54, 55, 56)
        );
        assert_eq!(
            heap.slice(output.atom_index_orig.as_const()).unwrap(),
            &[0x1122_3344, 0x88, 0, 0, 0]
        );
        assert_eq!(
            heap.slice(output.atom_index_fin.as_const()).unwrap(),
            &[0x1020_3040, 0x80, 0, 0, 0]
        );
        let haptic = heap.slice(output.lists_haptic_bonds.as_const()).unwrap()[0];
        let steabs = heap.slice(output.lists_steabs.as_const()).unwrap()[0];
        let sterac = heap.slice(output.lists_sterac.as_const()).unwrap()[0];
        let sterel = heap.slice(output.lists_sterel.as_const()).unwrap()[0];
        assert_eq!(heap.slice(haptic.as_const()).unwrap(), &[7, 8, 2, 11, 12]);
        assert_eq!(heap.slice(steabs.as_const()).unwrap(), &[0, 2, 21, 22]);
        assert_eq!(heap.slice(sterac.as_const()).unwrap(), &[31, 2, 32, 33]);
        assert_eq!(heap.slice(sterel.as_const()).unwrap(), &[41, 2, 42, 43]);
        FreeExtOrigAtData(&mut heap, polymer_slot, v3000_slot).unwrap();

        for invalid_bond in [0, 2] {
            let mut heap = SourceHeap::default();
            let input = fixture(&mut heap, invalid_bond, b"C", b"O", false);
            let old_v3000 = heap.allocate(vec![OAD_V3000::default()]).unwrap();
            let mut polymer_slot = SourceMutPointer::null();
            let mut v3000_slot = old_v3000;
            let mut errors = [0_i8; 256];
            errors[..6].copy_from_slice(&[
                b'p' as i8, b'r' as i8, b'i' as i8, b'o' as i8, b'r' as i8, 0,
            ]);
            assert_eq!(
                SetExtOrigAtDataByMolfileExtInput(
                    &mut heap,
                    &input,
                    &mut polymer_slot,
                    &mut v3000_slot,
                    Some(&mut errors),
                ),
                Ok(9004)
            );
            assert!(polymer_slot.is_null());
            assert_eq!(v3000_slot, old_v3000);
            assert_eq!(
                error_text(&errors),
                b"prior; Polymer unit in Molfile refers to invalid bond"
            );
            inchi_free(&mut heap, old_v3000).unwrap();
        }

        for (first, second) in [(b"H".as_slice(), b"C".as_slice()), (b"C", b"H")] {
            let mut heap = SourceHeap::default();
            let input = fixture(&mut heap, 1, first, second, false);
            let mut polymer_slot = SourceMutPointer::null();
            let mut v3000_slot = SourceMutPointer::null();
            let mut errors = [0_i8; 256];
            assert_eq!(
                SetExtOrigAtDataByMolfileExtInput(
                    &mut heap,
                    &input,
                    &mut polymer_slot,
                    &mut v3000_slot,
                    Some(&mut errors),
                ),
                Ok(9002)
            );
            assert!(polymer_slot.is_null());
            assert_eq!(
                error_text(&errors),
                b"Hydrogen as polymer end group is not supported"
            );
        }

        let mut no_bond_heap = SourceHeap::default();
        let no_bond_input = fixture(&mut no_bond_heap, 1, b"C", b"O", false);
        let group_pointer = no_bond_heap
            .slice(no_bond_input.ctab.sgroups.group.as_const())
            .unwrap()[0];
        no_bond_heap.slice_mut(group_pointer).unwrap()[0].blist.used = 0;
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetExtOrigAtDataByMolfileExtInput(
                &mut no_bond_heap,
                &no_bond_input,
                &mut polymer_slot,
                &mut v3000_slot,
                None,
            ),
            Ok(0)
        );
        let polymer = &no_bond_heap.slice(polymer_slot.as_const()).unwrap()[0];
        let unit = no_bond_heap.slice(polymer.units.as_const()).unwrap()[0];
        assert!(
            no_bond_heap.slice(unit.as_const()).unwrap()[0]
                .blist
                .is_null()
        );
        FreeExtOrigAtData(&mut no_bond_heap, polymer_slot, SourceMutPointer::null()).unwrap();

        for failure_after in 0..16_u64 {
            let mut heap = SourceHeap::default();
            let input = fixture(&mut heap, 1, b"C", b"O", true);
            heap.fail_after_allocations(failure_after);
            let old_polymer = heap
                .allocate_model_storage(vec![OAD_Polymer::default()])
                .unwrap();
            let old_v3000 = heap
                .allocate_model_storage(vec![OAD_V3000::default()])
                .unwrap();
            let mut polymer_slot = old_polymer;
            let mut v3000_slot = old_v3000;
            let mut errors = [0_i8; 256];
            assert_eq!(
                SetExtOrigAtDataByMolfileExtInput(
                    &mut heap,
                    &input,
                    &mut polymer_slot,
                    &mut v3000_slot,
                    Some(&mut errors),
                ),
                Ok(9001),
                "allocation ordinal {failure_after}"
            );
            assert!(polymer_slot.is_null(), "allocation ordinal {failure_after}");
            assert_eq!(error_text(&errors), b"Out of RAM");
            if failure_after < 5 {
                assert_eq!(v3000_slot, old_v3000, "allocation ordinal {failure_after}");
            } else if failure_after == 5 {
                assert!(v3000_slot.is_null(), "allocation ordinal {failure_after}");
            } else {
                assert_eq!(
                    heap.slice(v3000_slot.as_const()),
                    Err(SourceHeapError::MissingAllocation),
                    "allocation ordinal {failure_after}"
                );
            }
        }
    }

    #[test]
    fn source_port__mol2atom__setinpatomsxyz__line_975() {
        fn fixture(
            heap: &mut SourceHeap,
            coordinates: &[(f64, f64, f64)],
            bonds: Vec<MOL_FMT_BOND>,
        ) -> MOL_FMT_DATA {
            let atoms = heap
                .allocate(
                    coordinates
                        .iter()
                        .map(|&(fx, fy, fz)| MOL_FMT_ATOM {
                            fx,
                            fy,
                            fz,
                            ..MOL_FMT_ATOM::default()
                        })
                        .collect(),
                )
                .unwrap();
            let bond_count = bonds.len() as i32;
            let bond_pointer = if bonds.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate(bonds).unwrap()
            };
            MOL_FMT_DATA {
                ctab: crate::source_types::MOL_FMT_CTAB {
                    n_atoms: coordinates.len() as i32,
                    n_bonds: bond_count,
                    atoms,
                    bonds: bond_pointer,
                    ..crate::source_types::MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let point = fixture(&mut heap, &[(7.0, 8.0, 9.0)], vec![]);
        let point_output = heap
            .allocate(vec![inp_ATOM {
                x: 1.0,
                y: 2.0,
                z: 3.0,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let mut err = 41;
        assert_eq!(
            SetInpAtomsXYZ(&mut heap, &point, 1, point_output, &mut err, None),
            Ok(0)
        );
        assert_eq!(err, 41);
        assert_eq!(
            (
                heap.slice(point_output.as_const()).unwrap()[0].x,
                heap.slice(point_output.as_const()).unwrap()[0].y,
                heap.slice(point_output.as_const()).unwrap()[0].z,
            ),
            (1.0, 2.0, 3.0)
        );

        let planar = fixture(&mut heap, &[(-1.25, 2.5, -0.0), (4.75, -6.5, 0.0)], vec![]);
        let planar_output = heap
            .allocate(vec![inp_ATOM::default(), inp_ATOM::default()])
            .unwrap();
        assert_eq!(
            SetInpAtomsXYZ(&mut heap, &planar, 1, planar_output, &mut err, None),
            Ok(2)
        );
        let atoms = heap.slice(planar_output.as_const()).unwrap();
        assert_eq!((atoms[0].x, atoms[0].y), (-1.25, 2.5));
        assert_eq!(atoms[0].z.to_bits(), (-0.0_f64).to_bits());
        assert_eq!((atoms[1].x, atoms[1].y, atoms[1].z), (0.0, 0.0, 0.0));

        heap.slice_mut(planar_output).unwrap()[0].x = 19.0;
        assert_eq!(
            SetInpAtomsXYZ(&mut heap, &planar, -1, planar_output, &mut err, None),
            Ok(2)
        );
        assert_eq!(heap.slice(planar_output.as_const()).unwrap()[0].x, 19.0);

        let spatial = fixture(
            &mut heap,
            &[(1.0, 2.0, 3.0), (-4.0, 5.0, -6.0)],
            vec![MOL_FMT_BOND {
                atnum1: 0,
                atnum2: 1,
                ..MOL_FMT_BOND::default()
            }],
        );
        let spatial_output = heap
            .allocate(vec![inp_ATOM::default(), inp_ATOM::default()])
            .unwrap();
        let mut errors = [0_i8; 256];
        err = 0;
        assert_eq!(
            SetInpAtomsXYZ(
                &mut heap,
                &spatial,
                2,
                spatial_output,
                &mut err,
                Some(&mut errors),
            ),
            Ok(3)
        );
        let atoms = heap.slice(spatial_output.as_const()).unwrap();
        assert_eq!((atoms[0].x, atoms[0].y, atoms[0].z), (1.0, 2.0, 3.0));
        assert_eq!((atoms[1].x, atoms[1].y, atoms[1].z), (-4.0, 5.0, -6.0));
        assert_eq!(err, 1);
        let error_length = errors.iter().position(|byte| *byte == 0).unwrap();
        assert_eq!(
            errors[..error_length]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            b"Bond to nonexistent atom"
        );
    }

    #[test]
    fn source_port__mol2atom__makeinpatomsfrommolfiledata__line_508() {
        fn symbol(value: &[u8]) -> [i8; 6] {
            assert!(value.len() < 6);
            let mut result = [0_i8; 6];
            for (target, source) in result.iter_mut().zip(value.iter().copied()) {
                *target = source as i8;
            }
            result
        }

        fn mol_atom(value: &[u8]) -> MOL_FMT_ATOM {
            MOL_FMT_ATOM {
                symbol: symbol(value),
                ..MOL_FMT_ATOM::default()
            }
        }

        fn fixture(
            heap: &mut SourceHeap,
            atoms: Vec<MOL_FMT_ATOM>,
            bonds: Vec<MOL_FMT_BOND>,
        ) -> MOL_FMT_DATA {
            let atom_count = atoms.len() as i32;
            let bond_count = bonds.len() as i32;
            let atom_pointer = heap.allocate(atoms).unwrap();
            let bond_pointer = if bonds.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate(bonds).unwrap()
            };
            MOL_FMT_DATA {
                ctab: crate::source_types::MOL_FMT_CTAB {
                    n_atoms: atom_count,
                    n_bonds: bond_count,
                    atoms: atom_pointer,
                    bonds: bond_pointer,
                    ..crate::source_types::MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            }
        }

        fn error_text(buffer: &[i8]) -> String {
            String::from_utf8(
                buffer
                    .iter()
                    .take_while(|byte| **byte != 0)
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap()
        }

        let mut no_structure_heap = SourceHeap::default();
        let caller_atoms = no_structure_heap
            .allocate(vec![inp_ATOM {
                charge: 17,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let mut num_atoms = 99;
        let mut num_bonds = 88;
        let mut err = 77;
        assert!(
            MakeInpAtomsFromMolfileData(
                &mut no_structure_heap,
                &MOL_FMT_DATA::default(),
                &mut num_atoms,
                &mut num_bonds,
                caller_atoms,
                0,
                &mut err,
                None,
            )
            .unwrap()
            .is_null()
        );
        assert_eq!((num_atoms, num_bonds, err), (0, 0, 0));
        assert_eq!(
            no_structure_heap.slice(caller_atoms.as_const()).unwrap()[0].charge,
            17
        );

        let mut allocation_heap = SourceHeap::default();
        let allocation_input = fixture(&mut allocation_heap, vec![mol_atom(b"C")], vec![]);
        allocation_heap.fail_after_allocations(0);
        let mut errors = [0_i8; 256];
        let mut num_atoms = 0;
        let mut num_bonds = 0;
        let mut err = 0;
        assert!(
            MakeInpAtomsFromMolfileData(
                &mut allocation_heap,
                &allocation_input,
                &mut num_atoms,
                &mut num_bonds,
                SourceMutPointer::null(),
                0,
                &mut err,
                Some(&mut errors),
            )
            .unwrap()
            .is_null()
        );
        assert_eq!((num_atoms, num_bonds, err), (1, 0, -1));
        assert_eq!(error_text(&errors), "Out of RAM");

        let mut atom_heap = SourceHeap::default();
        let mut deuterium = mol_atom(b"D");
        deuterium.radical = RADICAL_SINGLET as i8;
        let mut zero_difference = mol_atom(b"C");
        zero_difference.mass_difference = ZERO_ATW_DIFF as i8;
        let mut positive_difference = mol_atom(b"O");
        positive_difference.mass_difference = 1;
        let atom_input = fixture(
            &mut atom_heap,
            vec![deuterium, zero_difference, positive_difference],
            vec![],
        );
        let caller_output = atom_heap
            .allocate(vec![
                inp_ATOM {
                    x: 11.0,
                    ..inp_ATOM::default()
                };
                3
            ])
            .unwrap();
        let mut num_atoms = 0;
        let mut num_bonds = 0;
        let mut err = 0;
        let output = MakeInpAtomsFromMolfileData(
            &mut atom_heap,
            &atom_input,
            &mut num_atoms,
            &mut num_bonds,
            caller_output,
            1,
            &mut err,
            None,
        )
        .unwrap();
        assert_eq!(output, caller_output);
        assert_eq!((num_atoms, num_bonds, err), (3, 0, 0));
        let atoms = atom_heap.slice(output.as_const()).unwrap();
        assert_eq!(&atoms[0].elname[..2], &[b'H' as i8, 0]);
        assert_eq!(atoms[0].el_number, 1);
        assert_eq!(atoms[0].iso_atw_diff, 2);
        assert_eq!(atoms[0].radical, RADICAL_TRIPLET as i8);
        assert_eq!(atoms[0].orig_at_number, 1);
        assert_eq!(atoms[0].x, 11.0);
        assert_eq!(atoms[1].iso_atw_diff, 1);
        assert_eq!(atoms[2].iso_atw_diff, 2);

        let mut unknown_heap = SourceHeap::default();
        let unknown_input = fixture(
            &mut unknown_heap,
            vec![mol_atom(b"C"), mol_atom(b"Q"), mol_atom(b"O")],
            vec![],
        );
        let unknown_output = unknown_heap
            .allocate(vec![
                inp_ATOM {
                    charge: 91,
                    ..inp_ATOM::default()
                };
                3
            ])
            .unwrap();
        let mut errors = [0_i8; 256];
        let mut num_atoms = 0;
        let mut num_bonds = 0;
        let mut err = 0;
        assert_eq!(
            MakeInpAtomsFromMolfileData(
                &mut unknown_heap,
                &unknown_input,
                &mut num_atoms,
                &mut num_bonds,
                unknown_output,
                0,
                &mut err,
                Some(&mut errors),
            )
            .unwrap(),
            unknown_output
        );
        assert_eq!((num_atoms, num_bonds, err), (3, 0, -2));
        let atoms = unknown_heap.slice(unknown_output.as_const()).unwrap();
        assert_eq!(atoms[0].el_number, 6);
        assert_eq!(atoms[1].el_number, 0);
        assert_eq!(atoms[2].charge, 91);
        assert_eq!(error_text(&errors), "Unknown element(s): Q");

        for (input_stereo, expected_first, expected_second) in [
            (0_i8, 0_i8, 0_i8),
            (
                INPUT_STEREO_DBLE_EITHER as i8,
                STEREO_DBLE_EITHER as i8,
                STEREO_DBLE_EITHER as i8,
            ),
            (
                INPUT_STEREO_SNGL_UP as i8,
                STEREO_SNGL_UP as i8,
                -(STEREO_SNGL_UP as i8),
            ),
            (
                INPUT_STEREO_SNGL_EITHER as i8,
                STEREO_SNGL_EITHER as i8,
                -(STEREO_SNGL_EITHER as i8),
            ),
            (
                INPUT_STEREO_SNGL_DOWN as i8,
                STEREO_SNGL_DOWN as i8,
                -(STEREO_SNGL_DOWN as i8),
            ),
        ] {
            let mut heap = SourceHeap::default();
            let input = fixture(
                &mut heap,
                vec![mol_atom(b"C"), mol_atom(b"O")],
                vec![MOL_FMT_BOND {
                    atnum1: 1,
                    atnum2: 2,
                    bond_type: 2,
                    bond_stereo: input_stereo,
                }],
            );
            let mut num_atoms = 0;
            let mut num_bonds = 0;
            let mut err = 0;
            let output = MakeInpAtomsFromMolfileData(
                &mut heap,
                &input,
                &mut num_atoms,
                &mut num_bonds,
                SourceMutPointer::null(),
                1,
                &mut err,
                None,
            )
            .unwrap();
            assert_eq!((num_atoms, num_bonds, err), (2, 1, 0));
            let atoms = heap.slice(output.as_const()).unwrap();
            assert_eq!(atoms[0].neighbor[0], 1);
            assert_eq!(atoms[1].neighbor[0], 0);
            assert_eq!(atoms[0].bond_type[0], 2);
            assert_eq!(atoms[1].bond_type[0], 2);
            assert_eq!(atoms[0].bond_stereo[0], expected_first);
            assert_eq!(atoms[1].bond_stereo[0], expected_second);
            assert_eq!(atoms[0].chem_bonds_valence, 2);
            assert_eq!(atoms[1].chem_bonds_valence, 2);
        }

        let mut warning_heap = SourceHeap::default();
        let warning_input = fixture(
            &mut warning_heap,
            vec![mol_atom(b"C"), mol_atom(b"O"), mol_atom(b"N")],
            vec![
                MOL_FMT_BOND {
                    atnum1: 1,
                    atnum2: 2,
                    bond_type: 9,
                    bond_stereo: 0,
                },
                MOL_FMT_BOND {
                    atnum1: 1,
                    atnum2: 2,
                    bond_type: 2,
                    bond_stereo: INPUT_STEREO_SNGL_EITHER as i8,
                },
                MOL_FMT_BOND {
                    atnum1: 2,
                    atnum2: 3,
                    bond_type: 1,
                    bond_stereo: 9,
                },
                MOL_FMT_BOND {
                    atnum1: 3,
                    atnum2: 3,
                    bond_type: 1,
                    bond_stereo: 0,
                },
            ],
        );
        let mut errors = [0_i8; 512];
        let mut num_atoms = 0;
        let mut num_bonds = 0;
        let mut err = 0;
        let output = MakeInpAtomsFromMolfileData(
            &mut warning_heap,
            &warning_input,
            &mut num_atoms,
            &mut num_bonds,
            SourceMutPointer::null(),
            1,
            &mut err,
            Some(&mut errors),
        )
        .unwrap();
        assert_eq!((num_atoms, num_bonds, err), (3, 2, 1 | 2 | 8 | 16));
        let atoms = warning_heap.slice(output.as_const()).unwrap();
        assert_eq!(atoms[0].valence, 1);
        assert_eq!(atoms[1].valence, 2);
        assert_eq!(atoms[2].valence, 1);
        assert_eq!(atoms[0].bond_type[0], 2);
        assert_eq!(atoms[0].bond_stereo[0], STEREO_SNGL_EITHER as i8);
        assert_eq!(atoms[1].bond_stereo[0], -(STEREO_SNGL_EITHER as i8));
        assert_eq!(atoms[1].neighbor[1], 2);
        assert_eq!(atoms[2].neighbor[0], 1);
        assert_eq!(atoms[1].bond_stereo[1], 0);
        assert_eq!(atoms[2].bond_stereo[0], 0);
        assert_eq!(
            error_text(&errors),
            "Unrecognized bond type: 9; Multiple bonds between two atoms; Unrecognized bond stereo; Bond to nonexistent atom"
        );
    }

    #[test]
    fn source_port__mol2atom__calculate_valences__line_807() {
        fn name(value: &[u8]) -> [i8; 6] {
            let mut result = [0_i8; 6];
            for (destination, source) in result.iter_mut().zip(value.iter().copied()) {
                *destination = source as i8;
            }
            result
        }

        fn error_text(buffer: &[i8]) -> Vec<u8> {
            buffer
                .iter()
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8)
                .collect()
        }

        fn atom(
            element: &[u8],
            element_number: u8,
            neighbors: &[u16],
            bond_types: &[u8],
        ) -> inp_ATOM {
            let mut result = inp_ATOM {
                elname: name(element),
                el_number: element_number,
                valence: neighbors.len() as i8,
                ..inp_ATOM::default()
            };
            result.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            result.bond_type[..bond_types.len()].copy_from_slice(bond_types);
            result
        }

        let heap = SourceHeap::default();
        let mut untouched = vec![inp_ATOM {
            chem_bonds_valence: 77,
            ..inp_ATOM::default()
        }];
        let mut err = 19;
        calculate_valences(&heap, None, &mut untouched, -1, 0, &mut err, None).unwrap();
        assert_eq!(untouched[0].chem_bonds_valence, 77);
        assert_eq!(err, 19);

        let mut ordinary = vec![atom(b"C", 6, &[1, 2, 3], &[1, 2, 3]); 4];
        ordinary[0].chem_bonds_valence = 99;
        let mut err = 0;
        calculate_valences(&heap, None, &mut ordinary, 1, 0, &mut err, None).unwrap();
        assert_eq!(ordinary[0].chem_bonds_valence, 6);
        assert_eq!(err, 0);

        for unknown in [0_u8, 5] {
            let mut atoms = vec![atom(b"C", 6, &[0], &[unknown])];
            let mut err = 7;
            let mut errors = [0_i8; 256];
            calculate_valences(&heap, None, &mut atoms, 1, 0, &mut err, Some(&mut errors)).unwrap();
            assert_eq!(atoms[0].chem_bonds_valence, 1);
            assert_eq!(err, 7);
            assert_eq!(
                error_text(&errors),
                b"Unknown bond type in MOLfile assigned as a single bond"
            );
        }

        let mut one_aromatic = vec![
            atom(b"C", 6, &[1], &[BOND_TYPE_ALTERN as u8]),
            atom(b"C", 6, &[0], &[BOND_TYPE_ALTERN as u8]),
        ];
        let mut err = 0;
        let mut errors = [0_i8; 256];
        calculate_valences(
            &heap,
            None,
            &mut one_aromatic,
            2,
            0,
            &mut err,
            Some(&mut errors),
        )
        .unwrap();
        assert_eq!(err, 32);
        assert_eq!(one_aromatic[0].bond_type[0], BOND_TYPE_SINGLE as u8);
        assert_eq!(one_aromatic[1].bond_type[0], BOND_TYPE_SINGLE as u8);
        assert_eq!(one_aromatic[0].chem_bonds_valence, 1);
        assert_eq!(one_aromatic[1].chem_bonds_valence, 1);
        assert_eq!(
            error_text(&errors),
            b"Atom has 1 or more than 3 aromatic bonds"
        );

        let mut broken = vec![
            atom(b"C", 6, &[1], &[BOND_TYPE_ALTERN as u8]),
            atom(b"C", 6, &[], &[]),
        ];
        let mut err = 32;
        let mut errors = [0_i8; 256];
        calculate_valences(&heap, None, &mut broken, 2, 0, &mut err, Some(&mut errors)).unwrap();
        assert_eq!(err, -2);
        assert_eq!(broken[0].bond_type[0], BOND_TYPE_ALTERN as u8);
        assert_eq!(error_text(&errors), b"Program error interpreting MOLfile");

        let mut two_aromatic = vec![
            atom(b"C", 6, &[1, 2], &[4, 4]),
            atom(b"C", 6, &[0], &[4]),
            atom(b"C", 6, &[0], &[4]),
        ];
        let mut err = 0;
        calculate_valences(&heap, None, &mut two_aromatic, 1, 0, &mut err, None).unwrap();
        assert_eq!(two_aromatic[0].chem_bonds_valence, 3);

        let mut three_aromatic = vec![
            atom(b"C", 6, &[1, 2, 3], &[4, 4, 4]),
            atom(b"C", 6, &[0], &[4]),
            atom(b"C", 6, &[0], &[4]),
            atom(b"C", 6, &[0], &[4]),
        ];
        let mut err = 0;
        calculate_valences(&heap, None, &mut three_aromatic, 1, 0, &mut err, None).unwrap();
        assert_eq!(three_aromatic[0].chem_bonds_valence, 4);

        let mut aromatic_n = vec![
            atom(b"N", 7, &[1, 2], &[4, 4]),
            atom(b"C", 6, &[0], &[4]),
            atom(b"C", 6, &[0], &[4]),
        ];
        aromatic_n[0].num_H = 1;
        let mut err = 0;
        calculate_valences(&heap, None, &mut aromatic_n, 1, 0, &mut err, None).unwrap();
        assert_eq!(aromatic_n[0].chem_bonds_valence, 2);

        let mut metal_neighbor = vec![
            atom(b"N", 7, &[1, 2], &[4, 4]),
            atom(b"Fe", 26, &[0], &[4]),
            atom(b"C", 6, &[0], &[4]),
        ];
        metal_neighbor[0].num_H = 1;
        let mut err = 0;
        calculate_valences(&heap, None, &mut metal_neighbor, 1, 0, &mut err, None).unwrap();
        assert_eq!(metal_neighbor[0].chem_bonds_valence, 3);

        let mut mf_heap = SourceHeap::default();
        let mol_atoms = mf_heap
            .allocate(vec![
                MOL_FMT_ATOM {
                    valence: 2,
                    ..MOL_FMT_ATOM::default()
                },
                MOL_FMT_ATOM::default(),
                MOL_FMT_ATOM::default(),
            ])
            .unwrap();
        let mfdata = MOL_FMT_DATA {
            ctab: crate::source_types::MOL_FMT_CTAB {
                atoms: mol_atoms,
                ..crate::source_types::MOL_FMT_CTAB::default()
            },
            ..MOL_FMT_DATA::default()
        };
        let mut known_aromatic = vec![
            atom(b"C", 6, &[1, 2], &[4, 4]),
            atom(b"C", 6, &[0], &[4]),
            atom(b"C", 6, &[0], &[4]),
        ];
        let mut err = 0;
        calculate_valences(
            &mf_heap,
            Some(&mfdata),
            &mut known_aromatic,
            1,
            0,
            &mut err,
            None,
        )
        .unwrap();
        assert_eq!(known_aromatic[0].chem_bonds_valence, 2);
        assert_eq!(known_aromatic[0].num_H, 0);

        let mut hydrogen_heap = SourceHeap::default();
        let mol_atoms = hydrogen_heap
            .allocate(vec![
                MOL_FMT_ATOM::default(),
                MOL_FMT_ATOM::default(),
                MOL_FMT_ATOM {
                    atom_aliased_flag: 1,
                    ..MOL_FMT_ATOM::default()
                },
            ])
            .unwrap();
        let mfdata = MOL_FMT_DATA {
            ctab: crate::source_types::MOL_FMT_CTAB {
                atoms: mol_atoms,
                ..crate::source_types::MOL_FMT_CTAB::default()
            },
            ..MOL_FMT_DATA::default()
        };
        let mut hydrogen_atoms = vec![
            atom(b"C", 6, &[], &[]),
            atom(b"C", 6, &[], &[]),
            atom(b"C", 6, &[], &[]),
        ];
        hydrogen_atoms[2].num_H = 2;
        let mut err = 0;
        calculate_valences(
            &hydrogen_heap,
            Some(&mfdata),
            &mut hydrogen_atoms,
            3,
            0,
            &mut err,
            None,
        )
        .unwrap();
        assert_eq!(hydrogen_atoms[0].num_H, 4);
        assert_eq!(hydrogen_atoms[2].num_H, 2);

        let mut no_h_atoms = vec![atom(b"C", 6, &[], &[])];
        let single_mol_atom = hydrogen_heap
            .allocate(vec![MOL_FMT_ATOM::default()])
            .unwrap();
        let single_data = MOL_FMT_DATA {
            ctab: crate::source_types::MOL_FMT_CTAB {
                atoms: single_mol_atom,
                ..crate::source_types::MOL_FMT_CTAB::default()
            },
            ..MOL_FMT_DATA::default()
        };
        calculate_valences(
            &hydrogen_heap,
            Some(&single_data),
            &mut no_h_atoms,
            1,
            1,
            &mut err,
            None,
        )
        .unwrap();
        assert_eq!(no_h_atoms[0].num_H, 0);
    }
}
