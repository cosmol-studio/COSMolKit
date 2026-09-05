use crate::source::api::inchi_dll::FreeInChIExtInput;
use crate::source::base::mol2atom::FreeExtOrigAtData;
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    OAD_Polymer, OAD_V3000, SourceHeap, SourceHeapError, SourceMutPointer, inchi_Input_Polymer,
    inchi_Input_PolymerUnit, inchi_Input_V3000,
};

const SOURCE_SIZEOF_OAD_POLYMER: u64 = 48;
const SOURCE_SIZEOF_OAD_POLYMER_UNIT: u64 = 240;
const SOURCE_SIZEOF_OAD_V3000: u64 = 104;
const SOURCE_SIZEOF_POINTER: u64 = 8;
const SOURCE_SIZEOF_INT: u64 = 4;
const SOURCE_SIZEOF_INCHI_INPUT_POLYMER: u64 = 16;
const SOURCE_SIZEOF_INCHI_INPUT_POLYMER_UNIT: u64 = 192;
const SOURCE_SIZEOF_INCHI_INPUT_V3000: u64 = 104;

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn SetInChIExtInputByExtOrigAtData(
    heap: &mut SourceHeap,
    orp: SourceMutPointer<OAD_Polymer>,
    orv: SourceMutPointer<OAD_V3000>,
    iip: &mut SourceMutPointer<inchi_Input_Polymer>,
    iiv: &mut SourceMutPointer<inchi_Input_V3000>,
    nat: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3258 SetInChIExtInputByExtOrigAtData
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int SetInChIExtInputByExtOrigAtData( OAD_Polymer     *orp,
                                     OAD_V3000     *orv,
                                     inchi_Input_Polymer **iip,
                                     inchi_Input_V3000     **iiv,
                                     int nat )
{
    int    k, m, err = 0;

        /* Polymers */
    if (orp && orp->n > 0)
    {
        /* djb-rwth: fixing oss-fuzz issue #67695, #66748 */
        inchi_Input_Polymer* iip_tmp = (inchi_Input_Polymer*) inchi_calloc( 1, sizeof( inchi_Input_Polymer ) );
        inchi_Input_PolymerUnit** units_tmp = (inchi_Input_PolymerUnit**) inchi_calloc( orp->n, sizeof( ( *iip )->units[0] ) );
        int** uk_al_tmp = (int**)inchi_malloc((orp->n) * sizeof(int*));
        inchi_Input_PolymerUnit** unitk = (inchi_Input_PolymerUnit**)inchi_malloc((orp->n) * sizeof(inchi_Input_PolymerUnit*));

        if (!iip_tmp || !units_tmp || !uk_al_tmp || !unitk)
        {
            err = 9001;
            goto preexitf;
        }
        /* *iip = iip_tmp; */
        iip_tmp->n = orp->n;
        iip_tmp->units = units_tmp;
        memset(units_tmp, 0, sizeof( *units_tmp) ); /* djb-rwth: memset_s C11/Annex K variant? */
        for (k = 0; k < orp->n; k++)
        {
            int q = 0;
            unitk[k] = (inchi_Input_PolymerUnit*)inchi_calloc(1, sizeof(inchi_Input_PolymerUnit));
            OAD_PolymerUnit    *groupk = orp->units[k];
            /* unitk = ( *iip )->units[k]; */
            if (!unitk[k])
            {
                err = 9001; 
                goto preexitf;
            }
            iip_tmp->units[k] = unitk[k];
            memset( unitk[k], 0, sizeof(*unitk[k])); /* djb-rwth: memset_s C11/Annex K variant? */
            unitk[k]->id = groupk->id;
            unitk[k]->type = groupk->type;
            unitk[k]->subtype = groupk->subtype;
            unitk[k]->conn = groupk->conn;
            unitk[k]->label = groupk->label;
            for (q = 0; q < 4; q++)
            {
                unitk[k]->xbr1[q] = groupk->xbr1[q];
                unitk[k]->xbr2[q] = groupk->xbr2[q];
            }
            strcpy( unitk[k]->smt, groupk->smt);
            unitk[k]->na = groupk->na;
            uk_al_tmp[k] = (int*)inchi_calloc(unitk[k]->na, sizeof(int));
            if (!uk_al_tmp[k])
            {
                err = 9001; 
                goto preexitf;
            }
            unitk[k]->alist = uk_al_tmp[k];
            for (m = 0; m < unitk[k]->na; m++)
            {
                uk_al_tmp[k][m] = groupk->alist[m];
            }
            unitk[k]->nb = groupk->nb;
            if (unitk[k]->nb > 0)
            {
                unitk[k]->blist = (int*)inchi_calloc(2 * (long long)unitk[k]->nb, sizeof(int)); /* djb-rwth: cast operator added */
                if (!unitk[k]->blist)
                {
                    err = 9001;
                    goto preexitf;
                }
                for (m = 0; m < 2 * groupk->nb; m++)
                {
                    unitk[k]->blist[m] = groupk->blist[m];
                }
            }
            else
            {
                unitk[k]->blist = NULL;
            }

            inchi_free(unitk[k]);
            inchi_free(uk_al_tmp[k]);
        }
    /* djb-rwth: avoiding memory leak */
    preexitf:
        /* djb-rwth: fixing GHI #165 */
        if (iip_tmp && *iip) /* djb-rwth: fixing oss-fuzz issue #455987437 */
        {
            memcpy(*iip, iip_tmp, sizeof(inchi_Input_Polymer));
        }
        inchi_free(iip_tmp);
        inchi_free(units_tmp);
        inchi_free(uk_al_tmp);
        inchi_free(unitk);
        if (err)
        {
            goto exitf;
        }
    }

    if (orv)
    {
        int nn;
        *iiv = (inchi_Input_V3000 *) inchi_calloc( 1, sizeof(inchi_Input_V3000) ); /* djb-rwth: fixing the incorrect type of variable */
        if (!*iiv)
        {
            err = 9001;
            goto exitf;
        }
        memset( *iiv, 0, sizeof( **iiv ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        ( *iiv )->n_collections = orv->n_collections;
        ( *iiv )->n_haptic_bonds = orv->n_haptic_bonds;
        ( *iiv )->n_non_haptic_bonds = orv->n_non_haptic_bonds;
        ( *iiv )->n_sgroups = orv->n_sgroups;
        ( *iiv )->n_non_star_atoms = orv->n_non_star_atoms;
        ( *iiv )->n_star_atoms = orv->n_star_atoms;
        ( *iiv )->n_steabs = orv->n_steabs;
        ( *iiv )->n_sterac = orv->n_sterac;
        ( *iiv )->n_sterel = orv->n_sterel;
        ( *iiv )->n_3d_constraints = orv->n_3d_constraints;

        if (orv->atom_index_orig)
        {
            ( *iiv )->atom_index_orig = (int *) inchi_calloc( nat, sizeof( int ) );
            if (NULL == ( *iiv )->atom_index_orig)
            {
                err = 9001;
                goto exitf;
            }
            memcpy( ( *iiv )->atom_index_orig, orv->atom_index_orig, nat );
        }
        if (orv->atom_index_fin)
        {
            ( *iiv )->atom_index_fin = (int *) inchi_calloc( nat, sizeof( int ) );
            if (NULL == ( *iiv )->atom_index_fin)
            {
                err = 9001;
                goto exitf;
            }
            memcpy( ( *iiv )->atom_index_fin, orv->atom_index_fin, nat );
        }
        if (orv->n_haptic_bonds && orv->lists_haptic_bonds)
        {
            ( *iiv )->lists_haptic_bonds = (int **) inchi_calloc( orv->n_haptic_bonds, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_haptic_bonds)
            {
                err = 9001;
                goto exitf;
            }
            for (m = 0; m < orv->n_haptic_bonds; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_haptic_bonds[m];
                nn = mol_lst[2] + 3;
                lst = ( *iiv )->lists_haptic_bonds[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (orv->n_steabs && orv->lists_steabs)
        {
            ( *iiv )->lists_steabs = (int **) inchi_calloc( orv->n_steabs, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_steabs) { err = 9001; goto exitf; }
            for (m = 0; m < orv->n_steabs; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_steabs[m];
                nn = mol_lst[1] + 2;
                lst = ( *iiv )->lists_steabs[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (orv->n_sterac && orv->lists_sterac)
        {
            ( *iiv )->lists_sterac = (int **) inchi_calloc( orv->n_sterac, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_sterac) { err = 9001; goto exitf; }
            for (m = 0; m < orv->n_sterac; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_sterac[m];
                nn = mol_lst[1] + 2;
                lst = ( *iiv )->lists_sterac[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (orv->n_sterel && orv->lists_sterel)
        {
            ( *iiv )->lists_sterel = (int **) inchi_calloc( orv->n_sterel, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_sterel) { err = 9001; goto exitf; }
            for (m = 0; m < orv->n_sterel; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_sterel[m];
                nn = mol_lst[1] + 2;
                lst = ( *iiv )->lists_sterel[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
    }

exitf:
    if (err)
    {
        FreeInChIExtInput( *iip, *iiv );
    }

    return err;
}
    */
    // END INCHI C FUNCTION: SetInChIExtInputByExtOrigAtData
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: SetInChIExtInputByExtOrigAtData
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_malloc/calloc/free are active libc macros.
    // INCHI✔️❌: The source frees each temporary polymer unit and atom list before copying only its header into an already non-NULL output.
    // INCHI✔️❌: The source memcpy calls for atom indices copy nat bytes into nat-int allocations; Rust reproduces that byte count exactly.
    // INCHI✔️❌: SourceHeap validation, map lookup, cloning, and bytewise memcpy modeling are materially slower than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: SetInChIExtInputByExtOrigAtData

    let mut err = 0_i32;

    'conversion: {
        if !orp.is_null() {
            let polymer = heap
                .slice(orp.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if polymer.n > 0 {
                let iip_tmp: SourceMutPointer<inchi_Input_Polymer> = match inchi_calloc(
                    heap,
                    1,
                    SOURCE_SIZEOF_INCHI_INPUT_POLYMER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                let units_tmp: SourceMutPointer<SourceMutPointer<inchi_Input_PolymerUnit>> = match inchi_calloc(
                    heap,
                    u64::try_from(polymer.n).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                let uk_al_tmp: SourceMutPointer<SourceMutPointer<i32>> = match inchi_calloc(
                    heap,
                    u64::try_from(polymer.n).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                let unitk: SourceMutPointer<SourceMutPointer<inchi_Input_PolymerUnit>> =
                    match inchi_calloc(
                        heap,
                        u64::try_from(polymer.n).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_POINTER,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => SourceMutPointer::null(),
                    };

                if iip_tmp.is_null()
                    || units_tmp.is_null()
                    || uk_al_tmp.is_null()
                    || unitk.is_null()
                {
                    err = 9001;
                } else {
                    {
                        let temporary = &mut heap.slice_mut(iip_tmp)?[0];
                        temporary.n = polymer.n;
                        temporary.units = units_tmp;
                    }

                    for k in 0..polymer.n {
                        let allocated_unit = match inchi_calloc(
                            heap,
                            1,
                            SOURCE_SIZEOF_INCHI_INPUT_POLYMER_UNIT,
                        ) {
                            Ok(pointer) => pointer,
                            Err(_) => SourceMutPointer::null(),
                        };
                        heap.slice_mut(unitk)?[k as usize] = allocated_unit;
                        let group_pointer = heap.slice(polymer.units.as_const())?[k as usize];
                        if allocated_unit.is_null() {
                            err = 9001;
                            break;
                        }
                        let group = heap
                            .slice(group_pointer.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone();
                        heap.slice_mut(units_tmp)?[k as usize] = allocated_unit;
                        {
                            let output = &mut heap.slice_mut(allocated_unit)?[0];
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
                            output.na = group.na;
                        }
                        let atom_list = match inchi_calloc(
                            heap,
                            u64::try_from(group.na).unwrap_or(u64::MAX),
                            SOURCE_SIZEOF_INT,
                        ) {
                            Ok(pointer) => pointer,
                            Err(_) => SourceMutPointer::null(),
                        };
                        heap.slice_mut(uk_al_tmp)?[k as usize] = atom_list;
                        if atom_list.is_null() {
                            err = 9001;
                            break;
                        }
                        heap.slice_mut(allocated_unit)?[0].alist = atom_list;
                        for m in 0..group.na {
                            heap.slice_mut(atom_list)?[m as usize] =
                                heap.slice(group.alist.as_const())?[m as usize];
                        }
                        heap.slice_mut(allocated_unit)?[0].nb = group.nb;
                        if group.nb > 0 {
                            let bond_list = match inchi_calloc(
                                heap,
                                u64::try_from(i64::from(group.nb) * 2).unwrap_or(u64::MAX),
                                SOURCE_SIZEOF_INT,
                            ) {
                                Ok(pointer) => pointer,
                                Err(_) => SourceMutPointer::null(),
                            };
                            heap.slice_mut(allocated_unit)?[0].blist = bond_list;
                            if bond_list.is_null() {
                                err = 9001;
                                break;
                            }
                            for m in 0..group.nb.wrapping_mul(2) {
                                heap.slice_mut(bond_list)?[m as usize] =
                                    heap.slice(group.blist.as_const())?[m as usize];
                            }
                        } else {
                            heap.slice_mut(allocated_unit)?[0].blist =
                                SourceMutPointer::null();
                        }

                        inchi_free(heap, allocated_unit)?;
                        inchi_free(heap, atom_list)?;
                    }
                }

                if !iip_tmp.is_null() && !iip.is_null() {
                    let copied = heap
                        .slice(iip_tmp.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    *heap
                        .slice_mut(*iip)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = copied;
                }
                inchi_free(heap, iip_tmp)?;
                inchi_free(heap, units_tmp)?;
                inchi_free(heap, uk_al_tmp)?;
                inchi_free(heap, unitk)?;
                if err != 0 {
                    break 'conversion;
                }
            }
        }

        if !orv.is_null() {
            let input = heap
                .slice(orv.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            *iiv = match inchi_calloc(heap, 1, SOURCE_SIZEOF_INCHI_INPUT_V3000) {
                Ok(pointer) => pointer,
                Err(_) => SourceMutPointer::null(),
            };
            if iiv.is_null() {
                err = 9001;
                break 'conversion;
            }
            {
                let output = &mut heap.slice_mut(*iiv)?[0];
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

            let nat_count = u64::try_from(nat).unwrap_or(u64::MAX);
            if !input.atom_index_orig.is_null() {
                let destination = match inchi_calloc(heap, nat_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                heap.slice_mut(*iiv)?[0].atom_index_orig = destination;
                if destination.is_null() {
                    err = 9001;
                    break 'conversion;
                }
                for byte_index in 0..usize::try_from(nat).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_orig.as_const())?[word] as u32
                        >> shift)
                        & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }
            if !input.atom_index_fin.is_null() {
                let destination = match inchi_calloc(heap, nat_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                heap.slice_mut(*iiv)?[0].atom_index_fin = destination;
                if destination.is_null() {
                    err = 9001;
                    break 'conversion;
                }
                for byte_index in 0..usize::try_from(nat).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_fin.as_const())?[word] as u32
                        >> shift)
                        & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }

            if input.n_haptic_bonds != 0 && !input.lists_haptic_bonds.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_haptic_bonds).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                heap.slice_mut(*iiv)?[0].lists_haptic_bonds = outer;
                if outer.is_null() {
                    err = 9001;
                    break 'conversion;
                }
                for m in 0..input.n_haptic_bonds {
                    let source = heap.slice(input.lists_haptic_bonds.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[2].wrapping_add(3);
                    let list = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => SourceMutPointer::null(),
                    };
                    heap.slice_mut(outer)?[m as usize] = list;
                    if list.is_null() {
                        err = 9001;
                        break 'conversion;
                    }
                    for k in 0..nn {
                        heap.slice_mut(list)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
                if err != 0 {
                    break 'conversion;
                }
            }
            if input.n_steabs != 0 && !input.lists_steabs.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_steabs).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                heap.slice_mut(*iiv)?[0].lists_steabs = outer;
                if outer.is_null() {
                    err = 9001;
                    break 'conversion;
                }
                for m in 0..input.n_steabs {
                    let source = heap.slice(input.lists_steabs.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let list = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => SourceMutPointer::null(),
                    };
                    heap.slice_mut(outer)?[m as usize] = list;
                    if list.is_null() {
                        err = 9001;
                        break 'conversion;
                    }
                    for k in 0..nn {
                        heap.slice_mut(list)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
                if err != 0 {
                    break 'conversion;
                }
            }
            if input.n_sterac != 0 && !input.lists_sterac.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_sterac).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                heap.slice_mut(*iiv)?[0].lists_sterac = outer;
                if outer.is_null() {
                    err = 9001;
                    break 'conversion;
                }
                for m in 0..input.n_sterac {
                    let source = heap.slice(input.lists_sterac.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let list = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => SourceMutPointer::null(),
                    };
                    heap.slice_mut(outer)?[m as usize] = list;
                    if list.is_null() {
                        err = 9001;
                        break 'conversion;
                    }
                    for k in 0..nn {
                        heap.slice_mut(list)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
                if err != 0 {
                    break 'conversion;
                }
            }
            if input.n_sterel != 0 && !input.lists_sterel.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_sterel).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => SourceMutPointer::null(),
                };
                heap.slice_mut(*iiv)?[0].lists_sterel = outer;
                if outer.is_null() {
                    err = 9001;
                    break 'conversion;
                }
                for m in 0..input.n_sterel {
                    let source = heap.slice(input.lists_sterel.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let list = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => SourceMutPointer::null(),
                    };
                    heap.slice_mut(outer)?[m as usize] = list;
                    if list.is_null() {
                        err = 9001;
                        break 'conversion;
                    }
                    for k in 0..nn {
                        heap.slice_mut(list)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
                if err != 0 {
                    break 'conversion;
                }
            }
        }
    }

    if err != 0 {
        FreeInChIExtInput(heap, *iip, *iiv)?;
    }
    Ok(err)
}

#[allow(non_snake_case)]
pub(crate) fn SetExtOrigAtDataByInChIExtInput(
    heap: &mut SourceHeap,
    pp_polymer: &mut SourceMutPointer<OAD_Polymer>,
    pp_v3000: &mut SourceMutPointer<OAD_V3000>,
    iep: SourceMutPointer<inchi_Input_Polymer>,
    iev: SourceMutPointer<inchi_Input_V3000>,
    nat: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3021 SetExtOrigAtDataByInChIExtInput
    // INCHI✔❌: int SetExtOrigAtDataByInChIExtInput( OAD_Polymer **ppPolymer,
    // INCHI✔❌:                                      OAD_V3000 **ppV3000,
    // INCHI✔❌:                                      inchi_Input_Polymer *iep,
    // INCHI✔❌:                                      inchi_Input_V3000 *iev,
    // INCHI✔❌:                                      int nat )
    // INCHI✔❌: {
    // INCHI✔❌:     int    k, m, err = 0;
    // INCHI✔❌:     OAD_V3000 *pv = NULL;
    // INCHI✔❌:
    // INCHI✔❌:     /* Polymers */
    // INCHI✔❌:     if (iep && iep->n)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Prepare OAD_Polymer container */
    // INCHI✔❌:         *ppPolymer = (OAD_Polymer *) inchi_calloc( 1, sizeof( OAD_Polymer ) );
    // INCHI✔❌:         if (!*ppPolymer)
    // INCHI✔❌:         {
    // INCHI✔❌:             err = 9001;
    // INCHI✔❌:             goto exitf;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         /* Convert Molfile's Sgroup's to OAD_PolymerUnit's */
    // INCHI✔❌:         ( *ppPolymer )->units = (OAD_PolymerUnit**) inchi_calloc( iep->n, sizeof( ( *ppPolymer )->units[0] ) );
    // INCHI✔❌:         if (!( *ppPolymer )->units)
    // INCHI✔❌:         {
    // INCHI✔❌:             err = 9001;
    // INCHI✔❌:             goto exitf;
    // INCHI✔❌:         }
    // INCHI✔❌:         memset( ( *ppPolymer )->units, 0, sizeof( *( *ppPolymer )->units ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:         ( *ppPolymer )->n = iep->n;
    // INCHI✔❌:         /*( *ppPolymer )->valid = -1;*/
    // INCHI✔❌:         ( *ppPolymer )->really_do_frame_shift = 0;
    // INCHI✔❌:
    // INCHI✔❌:         for (k = 0; k < iep->n; k++)
    // INCHI✔❌:         {
    // INCHI✔❌:             int q = 0;
    // INCHI✔❌:             OAD_PolymerUnit *unitk;
    // INCHI✔❌:
    // INCHI✔❌:             inchi_Input_PolymerUnit *groupk = iep->units[k];
    // INCHI✔❌:             ( *ppPolymer )->units[k] = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
    // INCHI✔❌:             unitk = ( *ppPolymer )->units[k];
    // INCHI✔❌:             if (!unitk)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             memset( unitk, 0, sizeof( *unitk ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:             unitk->id = groupk->id;
    // INCHI✔❌:             unitk->type = groupk->type;
    // INCHI✔❌:             unitk->subtype = groupk->subtype;
    // INCHI✔❌:             unitk->conn = groupk->conn;
    // INCHI✔❌:             unitk->label = groupk->label;
    // INCHI✔❌:
    // INCHI✔❌:             for (q = 0; q < 4; q++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->xbr1[q] = groupk->xbr1[q];
    // INCHI✔❌:                 unitk->xbr2[q] = groupk->xbr2[q];
    // INCHI✔❌:             }
    // INCHI✔❌:             strcpy( unitk->smt, groupk->smt );
    // INCHI✔❌:             unitk->na = groupk->na;
    // INCHI✔❌:             unitk->alist = (int *) inchi_calloc( unitk->na, sizeof( int ) );
    // INCHI✔❌:             if (!unitk->alist )
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             for (m = 0; m < unitk->na; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->alist[m] = groupk->alist[m];
    // INCHI✔❌:             }
    // INCHI✔❌:             unitk->nb = groupk->nb;
    // INCHI✔❌:             if (unitk->nb > 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->blist = (int *) inchi_calloc( 2 * (long long)unitk->nb, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔❌:                 if (!unitk->blist )
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (m = 0; m < 2 * groupk->nb; m++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     unitk->blist[m] = groupk->blist[m];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->blist = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* V3000 Extensions */
    // INCHI✔❌:     if (iev)
    // INCHI✔❌:     {
    // INCHI✔❌:         int nn;
    // INCHI✔❌:         *ppV3000 = (OAD_V3000 *) inchi_calloc( 1, sizeof( OAD_V3000 ) );
    // INCHI✔❌:         pv = *ppV3000;
    // INCHI✔❌:         if (!pv)
    // INCHI✔❌:         {
    // INCHI✔❌:             err = 9001;
    // INCHI✔❌:             goto exitf;
    // INCHI✔❌:         }
    // INCHI✔❌:         memset( pv, 0, sizeof( *pv ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:         pv->n_collections = iev->n_collections;
    // INCHI✔❌:         pv->n_haptic_bonds = iev->n_haptic_bonds;
    // INCHI✔❌:         pv->n_non_haptic_bonds = iev->n_non_haptic_bonds;
    // INCHI✔❌:         pv->n_sgroups = iev->n_sgroups;
    // INCHI✔❌:         pv->n_non_star_atoms = iev->n_non_star_atoms;
    // INCHI✔❌:         pv->n_star_atoms = iev->n_star_atoms;
    // INCHI✔❌:         pv->n_steabs = iev->n_steabs;
    // INCHI✔❌:         pv->n_sterac = iev->n_sterac;
    // INCHI✔❌:         pv->n_sterel = iev->n_sterel;
    // INCHI✔❌:         pv->n_3d_constraints = iev->n_3d_constraints;
    // INCHI✔❌:
    // INCHI✔❌:         if (iev->atom_index_orig)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->atom_index_orig = (int *) inchi_calloc( nat, sizeof( int ) );
    // INCHI✔❌:             if (NULL == pv->atom_index_orig)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             memcpy( pv->atom_index_orig, iev->atom_index_orig, nat );
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->atom_index_fin)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->atom_index_fin = (int *) inchi_calloc( nat, sizeof( int ) );
    // INCHI✔❌:             if (NULL == pv->atom_index_fin)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             memcpy( pv->atom_index_fin, iev->atom_index_fin, nat );
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_haptic_bonds && iev->lists_haptic_bonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_haptic_bonds = (int **) inchi_calloc( iev->n_haptic_bonds, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_haptic_bonds)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             for (m = 0; m < iev->n_haptic_bonds; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_haptic_bonds[m];
    // INCHI✔❌:                 nn = mol_lst[2] + 3;
    // INCHI✔❌:                 lst = pv->lists_haptic_bonds[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_steabs && iev->lists_steabs)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_steabs = (int **) inchi_calloc( iev->n_steabs, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_steabs) { err = 9001; goto exitf; }
    // INCHI✔❌:             for (m = 0; m < iev->n_steabs; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_steabs[m];
    // INCHI✔❌:                 nn = mol_lst[1] + 2;
    // INCHI✔❌:                 lst = pv->lists_steabs[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_sterac && iev->lists_sterac)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_sterac = (int **) inchi_calloc( iev->n_sterac, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_sterac) { err = 9001; goto exitf; }
    // INCHI✔❌:             for (m = 0; m < iev->n_sterac; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_sterac[m];
    // INCHI✔❌:                 nn = mol_lst[1] + 2;
    // INCHI✔❌:                 lst = pv->lists_sterac[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_sterel && iev->lists_sterel)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_sterel = (int **) inchi_calloc( iev->n_sterel, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_sterel) { err = 9001; goto exitf; }
    // INCHI✔❌:             for (m = 0; m < iev->n_sterel; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_sterel[m];
    // INCHI✔❌:                 nn = mol_lst[1] + 2;
    // INCHI✔❌:                 lst = pv->lists_sterel[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌: exitf:
    // INCHI✔❌:     if (err)
    // INCHI✔❌:     {
    // INCHI✔❌:         FreeExtOrigAtData( *ppPolymer, pv );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return err;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: SetExtOrigAtDataByInChIExtInput

    let mut err = 0;
    let mut pv = SourceMutPointer::null();

    'port: {
        if !iep.is_null() {
            let input_polymer = heap
                .slice(iep.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if input_polymer.n != 0 {
                *pp_polymer = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_POLYMER) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        SourceMutPointer::null()
                    }
                };
                if pp_polymer.is_null() {
                    break 'port;
                }
                let unit_count = u64::try_from(input_polymer.n).unwrap_or(u64::MAX);
                let units = match inchi_calloc(heap, unit_count, SOURCE_SIZEOF_POINTER) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                {
                    let polymer = &mut heap.slice_mut(*pp_polymer)?[0];
                    polymer.units = units;
                    polymer.n = input_polymer.n;
                    polymer.really_do_frame_shift = 0;
                }
                for k in 0..input_polymer.n {
                    let group_pointer = heap.slice(input_polymer.units.as_const())?[k as usize];
                    let group = heap
                        .slice(group_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let unit = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_POLYMER_UNIT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
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
                        output.na = group.na;
                    }
                    let atom_count = u64::try_from(group.na).unwrap_or(u64::MAX);
                    let alist = match inchi_calloc(heap, atom_count, SOURCE_SIZEOF_INT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(unit)?[0].alist = alist;
                    for m in 0..group.na {
                        heap.slice_mut(alist)?[m as usize] = heap.slice(group.alist.as_const())?[m as usize];
                    }
                    heap.slice_mut(unit)?[0].nb = group.nb;
                    if group.nb > 0 {
                        let bond_count = u64::try_from(i64::from(group.nb) * 2).unwrap_or(u64::MAX);
                        let blist = match inchi_calloc(heap, bond_count, SOURCE_SIZEOF_INT) {
                            Ok(pointer) => pointer,
                            Err(_) => {
                                err = 9001;
                                break 'port;
                            }
                        };
                        heap.slice_mut(unit)?[0].blist = blist;
                        for m in 0..group.nb * 2 {
                            heap.slice_mut(blist)?[m as usize] = heap.slice(group.blist.as_const())?[m as usize];
                        }
                    } else {
                        heap.slice_mut(unit)?[0].blist = SourceMutPointer::null();
                    }
                }
            }
        }

        if !iev.is_null() {
            let input = heap
                .slice(iev.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            *pp_v3000 = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_V3000) {
                Ok(pointer) => pointer,
                Err(_) => {
                    err = 9001;
                    SourceMutPointer::null()
                }
            };
            pv = *pp_v3000;
            if pv.is_null() {
                break 'port;
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
            let nat_count = u64::try_from(nat).unwrap_or(u64::MAX);
            if !input.atom_index_orig.is_null() {
                let destination = match inchi_calloc(heap, nat_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].atom_index_orig = destination;
                for byte_index in 0..usize::try_from(nat).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_orig.as_const())?[word] as u32 >> shift) & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }
            if !input.atom_index_fin.is_null() {
                let destination = match inchi_calloc(heap, nat_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].atom_index_fin = destination;
                for byte_index in 0..usize::try_from(nat).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_fin.as_const())?[word] as u32 >> shift) & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }

            if input.n_haptic_bonds != 0 && !input.lists_haptic_bonds.is_null() {
                let count = u64::try_from(input.n_haptic_bonds).unwrap_or(u64::MAX);
                let outer = match inchi_calloc(heap, count, SOURCE_SIZEOF_POINTER) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_haptic_bonds = outer;
                for m in 0..input.n_haptic_bonds {
                    let source = heap.slice(input.lists_haptic_bonds.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[2].wrapping_add(3);
                    let row = match inchi_calloc(heap, u64::try_from(nn).unwrap_or(u64::MAX), SOURCE_SIZEOF_INT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] = heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
            if input.n_steabs != 0 && !input.lists_steabs.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_steabs).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_steabs = outer;
                for m in 0..input.n_steabs {
                    let source = heap.slice(input.lists_steabs.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let row = match inchi_calloc(heap, u64::try_from(nn).unwrap_or(u64::MAX), SOURCE_SIZEOF_INT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] = heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
            if input.n_sterac != 0 && !input.lists_sterac.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_sterac).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_sterac = outer;
                for m in 0..input.n_sterac {
                    let source = heap.slice(input.lists_sterac.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let row = match inchi_calloc(heap, u64::try_from(nn).unwrap_or(u64::MAX), SOURCE_SIZEOF_INT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] = heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
            if input.n_sterel != 0 && !input.lists_sterel.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_sterel).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_sterel = outer;
                for m in 0..input.n_sterel {
                    let source = heap.slice(input.lists_sterel.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let row = match inchi_calloc(heap, u64::try_from(nn).unwrap_or(u64::MAX), SOURCE_SIZEOF_INT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] = heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
        }
    }

    if err != 0 {
        FreeExtOrigAtData(heap, *pp_polymer, pv)?;
    }
    Ok(err)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{OAD_PolymerUnit, inchi_Input_PolymerUnit};

    struct Fixture {
        polymer: SourceMutPointer<inchi_Input_Polymer>,
        v3000: SourceMutPointer<inchi_Input_V3000>,
    }

    struct ReverseFixture {
        polymer: SourceMutPointer<OAD_Polymer>,
        v3000: SourceMutPointer<OAD_V3000>,
    }

    fn allocate_reverse_fixture(heap: &mut SourceHeap, bonds: i32) -> ReverseFixture {
        let alist = heap.allocate(vec![3_i32, 4]).unwrap();
        let blist = if bonds > 0 {
            heap.allocate(vec![3_i32, 4]).unwrap()
        } else {
            SourceMutPointer::null()
        };
        let mut smt = [0_i8; 80];
        smt[..7].copy_from_slice(&[b'A' as i8, b'B' as i8, 0, 9, 8, 7, 6]);
        let unit = heap
            .allocate(vec![OAD_PolymerUnit {
                id: 11,
                type_: 12,
                subtype: 13,
                conn: 14,
                label: 15,
                na: 2,
                nb: bonds,
                xbr1: [1.0, 2.0, 3.0, 4.0],
                xbr2: [5.0, 6.0, 7.0, 8.0],
                smt,
                alist,
                blist,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate(vec![unit]).unwrap();
        let polymer = heap
            .allocate(vec![OAD_Polymer {
                units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();

        let atom_index_orig = heap.allocate(vec![0x1122_3344_i32, 0x5566_7788_i32]).unwrap();
        let atom_index_fin = heap.allocate(vec![0x1020_3040_i32, 0x5060_7080_i32]).unwrap();
        let haptic_row = heap.allocate(vec![7_i32, 8, 2, 11, 12]).unwrap();
        let steabs_row = heap.allocate(vec![0_i32, 2, 21, 22]).unwrap();
        let sterac_row = heap.allocate(vec![31_i32, 2, 32, 33]).unwrap();
        let sterel_row = heap.allocate(vec![41_i32, 2, 42, 43]).unwrap();
        let lists_haptic_bonds = heap.allocate(vec![haptic_row]).unwrap();
        let lists_steabs = heap.allocate(vec![steabs_row]).unwrap();
        let lists_sterac = heap.allocate(vec![sterac_row]).unwrap();
        let lists_sterel = heap.allocate(vec![sterel_row]).unwrap();
        let v3000 = heap
            .allocate(vec![OAD_V3000 {
                n_non_star_atoms: 51,
                n_star_atoms: 52,
                atom_index_orig,
                atom_index_fin,
                n_sgroups: 53,
                n_3d_constraints: 54,
                n_collections: 55,
                n_non_haptic_bonds: 56,
                n_haptic_bonds: 1,
                lists_haptic_bonds,
                n_steabs: 1,
                lists_steabs,
                n_sterel: 1,
                lists_sterel,
                n_sterac: 1,
                lists_sterac,
            }])
            .unwrap();
        ReverseFixture { polymer, v3000 }
    }

    fn allocate_fixture(heap: &mut SourceHeap, bonds: i32) -> Fixture {
        let alist = heap.allocate(vec![3_i32, 4]).unwrap();
        let blist = if bonds > 0 {
            heap.allocate(vec![3_i32, 4]).unwrap()
        } else {
            SourceMutPointer::null()
        };
        let mut smt = [0_i8; 80];
        smt[..7].copy_from_slice(&[b'A' as i8, b'B' as i8, 0, 9, 8, 7, 6]);
        let unit = heap
            .allocate(vec![inchi_Input_PolymerUnit {
                id: 11,
                type_: 12,
                subtype: 13,
                conn: 14,
                label: 15,
                na: 2,
                nb: bonds,
                xbr1: [1.0, 2.0, 3.0, 4.0],
                xbr2: [5.0, 6.0, 7.0, 8.0],
                smt,
                alist,
                blist,
            }])
            .unwrap();
        let units = heap.allocate(vec![unit]).unwrap();
        let polymer = heap.allocate(vec![inchi_Input_Polymer { units, n: 1 }]).unwrap();

        let atom_index_orig = heap.allocate(vec![0x1122_3344_i32, 0x5566_7788_i32]).unwrap();
        let atom_index_fin = heap.allocate(vec![0x1020_3040_i32, 0x5060_7080_i32]).unwrap();
        let haptic_row = heap.allocate(vec![7_i32, 8, 2, 11, 12]).unwrap();
        let steabs_row = heap.allocate(vec![0_i32, 2, 21, 22]).unwrap();
        let sterac_row = heap.allocate(vec![31_i32, 2, 32, 33]).unwrap();
        let sterel_row = heap.allocate(vec![41_i32, 2, 42, 43]).unwrap();
        let lists_haptic_bonds = heap.allocate(vec![haptic_row]).unwrap();
        let lists_steabs = heap.allocate(vec![steabs_row]).unwrap();
        let lists_sterac = heap.allocate(vec![sterac_row]).unwrap();
        let lists_sterel = heap.allocate(vec![sterel_row]).unwrap();
        let v3000 = heap
            .allocate(vec![inchi_Input_V3000 {
                n_non_star_atoms: 51,
                n_star_atoms: 52,
                atom_index_orig,
                atom_index_fin,
                n_sgroups: 53,
                n_3d_constraints: 54,
                n_collections: 55,
                n_non_haptic_bonds: 56,
                n_haptic_bonds: 1,
                lists_haptic_bonds,
                n_steabs: 1,
                lists_steabs,
                n_sterel: 1,
                lists_sterel,
                n_sterac: 1,
                lists_sterac,
            }])
            .unwrap();
        Fixture { polymer, v3000 }
    }

    #[test]
    fn source_port__inchi_dll__setinchiextinputbyextorigatdata__line_3258() {
        let mut null_heap = SourceHeap::default();
        let old_polymer = null_heap.allocate(vec![inchi_Input_Polymer::default()]).unwrap();
        let old_v3000 = null_heap.allocate(vec![inchi_Input_V3000::default()]).unwrap();
        let mut polymer_slot = old_polymer;
        let mut v3000_slot = old_v3000;
        assert_eq!(
            SetInChIExtInputByExtOrigAtData(
                &mut null_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut polymer_slot,
                &mut v3000_slot,
                i32::MIN,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, old_polymer);
        assert_eq!(v3000_slot, old_v3000);
        inchi_free(&mut null_heap, old_polymer).unwrap();
        inchi_free(&mut null_heap, old_v3000).unwrap();

        let mut zero_heap = SourceHeap::default();
        let zero_polymer = zero_heap
            .allocate(vec![OAD_Polymer {
                n: 0,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let old_output = zero_heap.allocate(vec![inchi_Input_Polymer::default()]).unwrap();
        let mut polymer_slot = old_output;
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetInChIExtInputByExtOrigAtData(
                &mut zero_heap,
                zero_polymer,
                SourceMutPointer::null(),
                &mut polymer_slot,
                &mut v3000_slot,
                0,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, old_output);
        assert_eq!(zero_heap.source_allocation_calls(), 2);
        inchi_free(&mut zero_heap, zero_polymer).unwrap();
        inchi_free(&mut zero_heap, old_output).unwrap();

        let mut no_output_heap = SourceHeap::default();
        let fixture = allocate_reverse_fixture(&mut no_output_heap, 1);
        no_output_heap.trace_source_allocations();
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetInChIExtInputByExtOrigAtData(
                &mut no_output_heap,
                fixture.polymer,
                fixture.v3000,
                &mut polymer_slot,
                &mut v3000_slot,
                5,
            ),
            Ok(0)
        );
        assert!(polymer_slot.is_null());
        assert_eq!(no_output_heap.source_allocation_calls(), 18);
        let output = &no_output_heap.slice(v3000_slot.as_const()).unwrap()[0];
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
            no_output_heap.slice(output.atom_index_orig.as_const()).unwrap(),
            &[0x1122_3344, 0x88, 0, 0, 0]
        );
        assert_eq!(
            no_output_heap.slice(output.atom_index_fin.as_const()).unwrap(),
            &[0x1020_3040, 0x80, 0, 0, 0]
        );
        let haptic = no_output_heap.slice(output.lists_haptic_bonds.as_const()).unwrap()[0];
        let steabs = no_output_heap.slice(output.lists_steabs.as_const()).unwrap()[0];
        let sterac = no_output_heap.slice(output.lists_sterac.as_const()).unwrap()[0];
        let sterel = no_output_heap.slice(output.lists_sterel.as_const()).unwrap()[0];
        assert_eq!(no_output_heap.slice(haptic.as_const()).unwrap(), &[7, 8, 2, 11, 12]);
        assert_eq!(no_output_heap.slice(steabs.as_const()).unwrap(), &[0, 2, 21, 22]);
        assert_eq!(no_output_heap.slice(sterac.as_const()).unwrap(), &[31, 2, 32, 33]);
        assert_eq!(no_output_heap.slice(sterel.as_const()).unwrap(), &[41, 2, 42, 43]);
        FreeInChIExtInput(&mut no_output_heap, SourceMutPointer::null(), v3000_slot).unwrap();
        FreeExtOrigAtData(&mut no_output_heap, fixture.polymer, fixture.v3000).unwrap();
        assert_eq!(no_output_heap.live_source_allocation_count(), 1);
        assert_eq!(no_output_heap.live_source_allocations_of::<i32>(), 1);

        let mut no_bond_heap = SourceHeap::default();
        let fixture = allocate_reverse_fixture(&mut no_bond_heap, 0);
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetInChIExtInputByExtOrigAtData(
                &mut no_bond_heap,
                fixture.polymer,
                SourceMutPointer::null(),
                &mut polymer_slot,
                &mut v3000_slot,
                0,
            ),
            Ok(0)
        );
        assert!(polymer_slot.is_null());
        FreeExtOrigAtData(&mut no_bond_heap, fixture.polymer, fixture.v3000).unwrap();
        assert_eq!(no_bond_heap.live_source_allocation_count(), 0);

        let mut existing_heap = SourceHeap::default();
        let fixture = allocate_reverse_fixture(&mut existing_heap, 0);
        let output_polymer = existing_heap.allocate(vec![inchi_Input_Polymer::default()]).unwrap();
        let mut polymer_slot = output_polymer;
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetInChIExtInputByExtOrigAtData(
                &mut existing_heap,
                fixture.polymer,
                SourceMutPointer::null(),
                &mut polymer_slot,
                &mut v3000_slot,
                0,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, output_polymer);
        let copied = &existing_heap.slice(output_polymer.as_const()).unwrap()[0];
        assert_eq!(copied.n, 1);
        assert_eq!(
            existing_heap.slice(copied.units.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        FreeExtOrigAtData(&mut existing_heap, fixture.polymer, fixture.v3000).unwrap();
        inchi_free(&mut existing_heap, output_polymer).unwrap();
        assert_eq!(existing_heap.live_source_allocation_count(), 0);

        let mut zero_lists_heap = SourceHeap::default();
        let source_row = zero_lists_heap.allocate(vec![1_i32, 0]).unwrap();
        let source_list = zero_lists_heap.allocate(vec![source_row]).unwrap();
        let zero_v3000 = zero_lists_heap
            .allocate(vec![OAD_V3000 {
                lists_haptic_bonds: source_list,
                lists_steabs: source_list,
                lists_sterac: source_list,
                lists_sterel: source_list,
                ..OAD_V3000::default()
            }])
            .unwrap();
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetInChIExtInputByExtOrigAtData(
                &mut zero_lists_heap,
                SourceMutPointer::null(),
                zero_v3000,
                &mut polymer_slot,
                &mut v3000_slot,
                0,
            ),
            Ok(0)
        );
        let output = &zero_lists_heap.slice(v3000_slot.as_const()).unwrap()[0];
        assert!(output.lists_haptic_bonds.is_null());
        assert!(output.lists_steabs.is_null());
        assert!(output.lists_sterac.is_null());
        assert!(output.lists_sterel.is_null());
        FreeInChIExtInput(&mut zero_lists_heap, SourceMutPointer::null(), v3000_slot).unwrap();
        inchi_free(&mut zero_lists_heap, source_row).unwrap();
        inchi_free(&mut zero_lists_heap, source_list).unwrap();
        inchi_free(&mut zero_lists_heap, zero_v3000).unwrap();

        let mut negative_heap = SourceHeap::default();
        let negative_list = negative_heap.allocate(vec![SourceMutPointer::<i32>::null()]).unwrap();
        let negative_source = negative_heap
            .allocate(vec![OAD_V3000 {
                n_haptic_bonds: -1,
                lists_haptic_bonds: negative_list,
                ..OAD_V3000::default()
            }])
            .unwrap();
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetInChIExtInputByExtOrigAtData(
                &mut negative_heap,
                SourceMutPointer::null(),
                negative_source,
                &mut polymer_slot,
                &mut v3000_slot,
                0,
            ),
            Ok(9001)
        );
        assert_eq!(
            negative_heap.slice(v3000_slot.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        inchi_free(&mut negative_heap, negative_list).unwrap();
        inchi_free(&mut negative_heap, negative_source).unwrap();

        for failure_after in 0..18_u64 {
            let mut heap = SourceHeap::default();
            let fixture = allocate_reverse_fixture(&mut heap, 1);
            heap.fail_after_allocations(failure_after);
            let mut polymer_slot = SourceMutPointer::null();
            let mut v3000_slot = SourceMutPointer::null();
            assert_eq!(
                SetInChIExtInputByExtOrigAtData(
                    &mut heap,
                    fixture.polymer,
                    fixture.v3000,
                    &mut polymer_slot,
                    &mut v3000_slot,
                    5,
                ),
                Ok(9001),
                "allocation ordinal {failure_after}"
            );
            assert!(polymer_slot.is_null(), "allocation ordinal {failure_after}");
            assert_eq!(
                heap.source_allocation_calls(),
                (failure_after + 1).max(4),
                "allocation ordinal {failure_after}"
            );
            if !v3000_slot.is_null() {
                assert_eq!(
                    heap.slice(v3000_slot.as_const()),
                    Err(SourceHeapError::MissingAllocation),
                    "allocation ordinal {failure_after}"
                );
            }
        }
    }

    #[test]
    fn source_port__inchi_dll__setextorigatdatabyinchiextinput__line_3021() {
        let mut heap = SourceHeap::default();
        let old_polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        let old_v3000 = heap.allocate(vec![OAD_V3000::default()]).unwrap();
        let mut polymer_slot = old_polymer;
        let mut v3000_slot = old_v3000;
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                i32::MIN,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, old_polymer);
        assert_eq!(v3000_slot, old_v3000);
        FreeExtOrigAtData(&mut heap, old_polymer, old_v3000).unwrap();

        let zero_input = heap
            .allocate(vec![inchi_Input_Polymer {
                units: SourceMutPointer::null(),
                n: 0,
            }])
            .unwrap();
        let untouched = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        let mut polymer_slot = untouched;
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                zero_input,
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, untouched);
        FreeExtOrigAtData(&mut heap, untouched, SourceMutPointer::null()).unwrap();

        let fixture = allocate_fixture(&mut heap, 1);
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                fixture.polymer,
                fixture.v3000,
                5,
            ),
            Ok(0)
        );
        let polymer = &heap.slice(polymer_slot.as_const()).unwrap()[0];
        assert_eq!(polymer.n, 1);
        assert_eq!(polymer.really_do_frame_shift, 0);
        let unit_pointer = heap.slice(polymer.units.as_const()).unwrap()[0];
        let unit = &heap.slice(unit_pointer.as_const()).unwrap()[0];
        assert_eq!(
            (unit.id, unit.type_, unit.subtype, unit.conn, unit.label),
            (11, 12, 13, 14, 15)
        );
        assert_eq!(unit.xbr1, [1.0, 2.0, 3.0, 4.0]);
        assert_eq!(unit.xbr2, [5.0, 6.0, 7.0, 8.0]);
        assert_eq!(&unit.smt[..7], &[b'A' as i8, b'B' as i8, 0, 0, 0, 0, 0]);
        assert_eq!(heap.slice(unit.alist.as_const()).unwrap(), &[3, 4]);
        assert_eq!(heap.slice(unit.blist.as_const()).unwrap(), &[3, 4]);

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

        let no_bond_fixture = allocate_fixture(&mut heap, 0);
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                no_bond_fixture.polymer,
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );
        let polymer = &heap.slice(polymer_slot.as_const()).unwrap()[0];
        let unit_pointer = heap.slice(polymer.units.as_const()).unwrap()[0];
        assert!(heap.slice(unit_pointer.as_const()).unwrap()[0].blist.is_null());
        FreeExtOrigAtData(&mut heap, polymer_slot, SourceMutPointer::null()).unwrap();

        for successful_allocations in 0..16_u64 {
            let mut heap = SourceHeap::default();
            let fixture = allocate_fixture(&mut heap, 1);
            heap.fail_after_allocations(successful_allocations);
            let mut polymer_slot = SourceMutPointer::null();
            let mut v3000_slot = SourceMutPointer::null();
            assert_eq!(
                SetExtOrigAtDataByInChIExtInput(
                    &mut heap,
                    &mut polymer_slot,
                    &mut v3000_slot,
                    fixture.polymer,
                    fixture.v3000,
                    5,
                ),
                Ok(9001),
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                heap.source_allocation_calls(),
                successful_allocations + 1,
                "allocation ordinal {successful_allocations}"
            );
            if !polymer_slot.is_null() {
                assert_eq!(
                    heap.slice(polymer_slot.as_const()),
                    Err(SourceHeapError::MissingAllocation),
                    "allocation ordinal {successful_allocations}"
                );
            }
            if !v3000_slot.is_null() {
                assert_eq!(
                    heap.slice(v3000_slot.as_const()),
                    Err(SourceHeapError::MissingAllocation),
                    "allocation ordinal {successful_allocations}"
                );
            }
        }
    }
}
