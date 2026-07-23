use crate::source::base::mol2atom::FreeExtOrigAtData;
use crate::source::base::strutil::Free_INChI_Members;
use crate::source::base::util::inchi_free;
use crate::source_types::{InpInChI, SourceHeap, SourceHeapError, SourceMutPointer};

#[allow(non_snake_case)]
pub(crate) fn FreeInpInChI(
    heap: &mut SourceHeap,
    p_one_input: &mut InpInChI,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1438 FreeInpInChI
    // INCHI✔❌: void FreeInpInChI(InpInChI* pOneInput)
    // INCHI✔❌: {
    // INCHI✔❌:     int iINChI, k, j;
    // INCHI✔❌:     for (iINChI = 0; iINChI < INCHI_NUM; iINChI++)
    // INCHI✔❌:     {
    // INCHI✔❌:         for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (pOneInput->pInpInChI[iINChI][j])
    // INCHI✔❌:             {
    // INCHI✔❌:                 for (k = 0; k < pOneInput->nNumComponents[iINChI][j]; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌: #if (FIX_OSS_FUZZ_25734_28139 == 1)
    // INCHI✔❌:                     U_CHAR* k_nAtom = (&pOneInput->pInpInChI[iINChI][j][k])->nAtom;
    // INCHI✔❌:                     AT_NUMB* k_nConnTable = (&pOneInput->pInpInChI[iINChI][j][k])->nConnTable;
    // INCHI✔❌:                     AT_NUMB* k_nTautomer = (&pOneInput->pInpInChI[iINChI][j][k])->nTautomer;
    // INCHI✔❌:                     S_CHAR* k_nNum_H = (&pOneInput->pInpInChI[iINChI][j][k])->nNum_H;
    // INCHI✔❌:                     S_CHAR* k_nNum_H_fixed = (&pOneInput->pInpInChI[iINChI][j][k])->nNum_H_fixed;
    // INCHI✔❌:                     char* k_szHillFormula = (&pOneInput->pInpInChI[iINChI][j][k])->szHillFormula;
    // INCHI✔❌:                     AT_NUMB* k_nPossibleLocationsOfIsotopicH = (&pOneInput->pInpInChI[iINChI][j][k])->nPossibleLocationsOfIsotopicH;
    // INCHI✔❌:                     INChI_IsotopicAtom* k_IsotopicAtom = (&pOneInput->pInpInChI[iINChI][j][k])->IsotopicAtom;
    // INCHI✔❌:                     INChI_IsotopicTGroup* k_IsotopicTGroup = (&pOneInput->pInpInChI[iINChI][j][k])->IsotopicTGroup;
    // INCHI✔❌:                     INChI_Stereo* k_Stereo = (&pOneInput->pInpInChI[iINChI][j][k])->Stereo;
    // INCHI✔❌:                     INChI_Stereo* k_StereoIsotopic = (&pOneInput->pInpInChI[iINChI][j][k])->StereoIsotopic;
    // INCHI✔❌:
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:                     Free_INChI_Members(&pOneInput->pInpInChI[iINChI][j][k]);
    // INCHI✔❌:
    // INCHI✔❌: #if (FIX_OSS_FUZZ_25734_28139 == 1)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* prevent erroneous repeated freeing in copied pInpInChIp[][][kk] */
    // INCHI✔❌:                         int kk;
    // INCHI✔❌:                         for (kk = k + 1; kk < pOneInput->nNumComponents[iINChI][j]; kk++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (k_nAtom == (&pOneInput->pInpInChI[iINChI][j][kk])->nAtom)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nAtom = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nConnTable == (&pOneInput->pInpInChI[iINChI][j][kk])->nConnTable)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nConnTable = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nTautomer == (&pOneInput->pInpInChI[iINChI][j][kk])->nTautomer)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nTautomer = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nNum_H == (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nNum_H_fixed == (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H_fixed)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H_fixed = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_szHillFormula == (&pOneInput->pInpInChI[iINChI][j][kk])->szHillFormula)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->szHillFormula = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nPossibleLocationsOfIsotopicH == (&pOneInput->pInpInChI[iINChI][j][kk])->nPossibleLocationsOfIsotopicH)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nPossibleLocationsOfIsotopicH = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_IsotopicAtom == (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicAtom)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicAtom = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_IsotopicTGroup == (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicTGroup)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicTGroup = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_Stereo == (&pOneInput->pInpInChI[iINChI][j][kk])->Stereo)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->Stereo = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_StereoIsotopic == (&pOneInput->pInpInChI[iINChI][j][kk])->StereoIsotopic)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->StereoIsotopic = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:                 }
    // INCHI✔❌:                 inchi_free(pOneInput->pInpInChI[iINChI][j]);
    // INCHI✔❌:                 pOneInput->pInpInChI[iINChI][j] = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:             if (pOneInput->nNumProtons[iINChI][j].pNumProtons)
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free(pOneInput->nNumProtons[iINChI][j].pNumProtons);
    // INCHI✔❌:                 pOneInput->nNumProtons[iINChI][j].pNumProtons = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     if (pOneInput->atom)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free(pOneInput->atom);
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     FreeExtOrigAtData(pOneInput->polymer, pOneInput->v3000);
    // INCHI✔❌:
    // INCHI✔❌:     memset(pOneInput, 0, sizeof(*pOneInput)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌: }
    // END INCHI C FUNCTION: FreeInpInChI

    for i_inchi in 0..2_usize {
        for j in 0..2_usize {
            let components = p_one_input.pInpInChI[i_inchi][j];
            let component_count = p_one_input.nNumComponents[i_inchi][j];
            if !components.is_null() {
                for k in 0..component_count {
                    let k = i64::from(k);
                    let component = components.offset(k)?;
                    let snapshot = heap
                        .slice(component.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    Free_INChI_Members(heap, component)?;
                    for kk in (k + 1)..i64::from(component_count) {
                        let later = components.offset(kk)?;
                        let value = heap
                            .slice_mut(later)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if snapshot.nAtom == value.nAtom {
                            value.nAtom = SourceMutPointer::null();
                        }
                        if snapshot.nConnTable == value.nConnTable {
                            value.nConnTable = SourceMutPointer::null();
                        }
                        if snapshot.nTautomer == value.nTautomer {
                            value.nTautomer = SourceMutPointer::null();
                        }
                        if snapshot.nNum_H == value.nNum_H {
                            value.nNum_H = SourceMutPointer::null();
                        }
                        if snapshot.nNum_H_fixed == value.nNum_H_fixed {
                            value.nNum_H_fixed = SourceMutPointer::null();
                        }
                        if snapshot.szHillFormula == value.szHillFormula {
                            value.szHillFormula = SourceMutPointer::null();
                        }
                        if snapshot.nPossibleLocationsOfIsotopicH
                            == value.nPossibleLocationsOfIsotopicH
                        {
                            value.nPossibleLocationsOfIsotopicH = SourceMutPointer::null();
                        }
                        if snapshot.IsotopicAtom == value.IsotopicAtom {
                            value.IsotopicAtom = SourceMutPointer::null();
                        }
                        if snapshot.IsotopicTGroup == value.IsotopicTGroup {
                            value.IsotopicTGroup = SourceMutPointer::null();
                        }
                        if snapshot.Stereo == value.Stereo {
                            value.Stereo = SourceMutPointer::null();
                        }
                        if snapshot.StereoIsotopic == value.StereoIsotopic {
                            value.StereoIsotopic = SourceMutPointer::null();
                        }
                    }
                }
                inchi_free(heap, components)?;
                p_one_input.pInpInChI[i_inchi][j] = SourceMutPointer::null();
            }
            let protons = p_one_input.nNumProtons[i_inchi][j].pNumProtons;
            if !protons.is_null() {
                inchi_free(heap, protons)?;
                p_one_input.nNumProtons[i_inchi][j].pNumProtons = SourceMutPointer::null();
            }
        }
    }
    if !p_one_input.atom.is_null() {
        inchi_free(heap, p_one_input.atom)?;
    }
    FreeExtOrigAtData(heap, p_one_input.polymer, p_one_input.v3000)?;
    *p_one_input = InpInChI::default();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        COMPONENT_REM_PROTONS, INChI, INChI_IsotopicAtom, INChI_IsotopicTGroup, INChI_Stereo,
        OAD_Polymer, OAD_V3000, inp_ATOM,
    };

    #[test]
    fn source_port__ichirvr7__freeinpinchi__line_1438() {
        let mut heap = SourceHeap::default();
        let mut empty = InpInChI::default();
        FreeInpInChI(&mut heap, &mut empty).unwrap();
        assert_eq!(empty, InpInChI::default());

        let hill = heap.allocate(vec![b'C' as i8, 0]).unwrap();
        let atom_member = heap.allocate(vec![6_u8]).unwrap();
        let connection = heap.allocate(vec![1_u16]).unwrap();
        let tautomer = heap.allocate(vec![2_u16]).unwrap();
        let num_h = heap.allocate(vec![1_i8]).unwrap();
        let num_h_fixed = heap.allocate(vec![2_i8]).unwrap();
        let possible_h = heap.allocate(vec![3_u16]).unwrap();
        let isotopic_atom = heap.allocate(vec![INChI_IsotopicAtom::default()]).unwrap();
        let isotopic_t_group = heap
            .allocate(vec![INChI_IsotopicTGroup::default()])
            .unwrap();
        let stereo = heap.allocate(vec![INChI_Stereo::default()]).unwrap();
        let stereo_isotopic = heap.allocate(vec![INChI_Stereo::default()]).unwrap();
        let shared = INChI {
            szHillFormula: hill,
            nAtom: atom_member,
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
        };
        let components = heap.allocate(vec![shared.clone(), shared]).unwrap();
        let second_components = heap.allocate(vec![INChI::default()]).unwrap();
        let proton_00 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let proton_01 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let proton_10 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let proton_11 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let atoms = heap.allocate(vec![inp_ATOM::default()]).unwrap();
        let polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        let v3000 = heap.allocate(vec![OAD_V3000::default()]).unwrap();

        let mut input = InpInChI::default();
        input.pInpInChI[0][0] = components;
        input.nNumComponents[0][0] = 2;
        input.pInpInChI[1][1] = second_components;
        input.nNumComponents[1][1] = 1;
        input.nNumProtons[0][0].pNumProtons = proton_00;
        input.nNumProtons[0][1].pNumProtons = proton_01;
        input.nNumProtons[1][0].pNumProtons = proton_10;
        input.nNumProtons[1][1].pNumProtons = proton_11;
        input.atom = atoms;
        input.num_atoms = 1;
        input.num_inp = 77;
        input.polymer = polymer;
        input.v3000 = v3000;

        FreeInpInChI(&mut heap, &mut input).unwrap();
        assert_eq!(input, InpInChI::default());
        for freed in [
            heap.slice(hill.as_const()).map(|_| ()),
            heap.slice(atom_member.as_const()).map(|_| ()),
            heap.slice(connection.as_const()).map(|_| ()),
            heap.slice(tautomer.as_const()).map(|_| ()),
            heap.slice(num_h.as_const()).map(|_| ()),
            heap.slice(num_h_fixed.as_const()).map(|_| ()),
            heap.slice(possible_h.as_const()).map(|_| ()),
            heap.slice(isotopic_atom.as_const()).map(|_| ()),
            heap.slice(isotopic_t_group.as_const()).map(|_| ()),
            heap.slice(stereo.as_const()).map(|_| ()),
            heap.slice(stereo_isotopic.as_const()).map(|_| ()),
            heap.slice(components.as_const()).map(|_| ()),
            heap.slice(second_components.as_const()).map(|_| ()),
            heap.slice(proton_00.as_const()).map(|_| ()),
            heap.slice(proton_01.as_const()).map(|_| ()),
            heap.slice(proton_10.as_const()).map(|_| ()),
            heap.slice(proton_11.as_const()).map(|_| ()),
            heap.slice(atoms.as_const()).map(|_| ()),
            heap.slice(polymer.as_const()).map(|_| ()),
            heap.slice(v3000.as_const()).map(|_| ()),
        ] {
            assert_eq!(freed, Err(SourceHeapError::MissingAllocation));
        }
    }
}
