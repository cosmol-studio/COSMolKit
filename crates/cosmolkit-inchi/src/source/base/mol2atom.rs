use crate::source::base::ichinorm::FreeInpAtomData;
use crate::source::base::runichi3::OAD_Polymer_Free;
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    COMP_ATOM_DATA, INP_ATOM_DATA, OAD_Polymer, OAD_V3000, ORIG_ATOM_DATA, SourceHeap,
    SourceHeapError, SourceMutPointer, inp_ATOM,
};

const SOURCE_SIZEOF_INP_ATOM: u64 = 176;

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

#[cfg(test)]
mod tests {
    use super::*;

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
}
