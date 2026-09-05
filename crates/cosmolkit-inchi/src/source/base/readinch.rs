use crate::source::base::ichi_io::inchi_ios_getsTab1;
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichister::{FixUnkn0DStereoBonds, ReconcileAllCmlBondParities};
use crate::source::base::util::{inchi_calloc, inchi_free, is_in_the_list};
use crate::source_types::{
    AB_PARITY_EVEN, AB_PARITY_NONE, AB_PARITY_ODD, AB_PARITY_UNDF, BOND_TYPE_DOUBLE,
    INCHI_IOSTREAM, MAX_NUM_STEREO_BONDS, NO_ATOM, SB_PARITY_FLAG, SB_PARITY_MASK, SB_PARITY_SHFT,
    SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, inchi_Stereo0D, inp_ATOM,
    tagINCHIStereoParity0D_INCHI_PARITY_EVEN, tagINCHIStereoParity0D_INCHI_PARITY_NONE,
    tagINCHIStereoParity0D_INCHI_PARITY_ODD, tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED,
    tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN, tagINCHIStereoType0D_INCHI_StereoType_Allene,
    tagINCHIStereoType0D_INCHI_StereoType_DoubleBond, tagINCHIStereoType0D_INCHI_StereoType_None,
    tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral,
};

// `inchi_api.h` under the pinned LP64 ABI: five 2-byte `AT_NUM` values and
// two 1-byte `S_CHAR` values total 12 bytes.
const SOURCE_SIZEOF_INCHI_STEREO0D: u64 = 5 * 2 + 2;

pub(crate) fn CreateInchi_Stereo0D(
    heap: &mut SourceHeap,
    num_stereo0d: i32,
) -> Result<SourceMutPointer<inchi_Stereo0D>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:118 CreateInchi_Stereo0D
    // INCHI✔️❌: inchi_Stereo0D * CreateInchi_Stereo0D( int num_stereo0D )
    // INCHI✔️❌: {
    // INCHI✔️❌:     return (inchi_Stereo0D*) inchi_calloc( num_stereo0D, sizeof( inchi_Stereo0D ) );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateInchi_Stereo0D

    inchi_calloc(heap, num_stereo0d as u64, SOURCE_SIZEOF_INCHI_STEREO0D)
}

pub(crate) fn FreeInchi_Stereo0D(
    heap: &mut SourceHeap,
    stereo_slot: Option<&mut SourceMutPointer<inchi_Stereo0D>>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:125 FreeInchi_Stereo0D
    // INCHI✔️❌: void FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (stereo0D && *stereo0D)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( *stereo0D );
    // INCHI✔️❌:         *stereo0D = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeInchi_Stereo0D

    let Some(stereo_slot) = stereo_slot else {
        return Ok(());
    };
    if stereo_slot.is_null() {
        return Ok(());
    }
    inchi_free(heap, *stereo_slot)?;
    *stereo_slot = SourceMutPointer::null();
    Ok(())
}

fn add_extract_0d_error(
    heap: &mut SourceHeap,
    error_buffer: Option<SourceMutPointer<i8>>,
    message: &str,
) -> Result<(), SourceHeapError> {
    let Some(error_buffer) = error_buffer else {
        return Ok(());
    };
    let mut message_bytes: Vec<i8> = message.bytes().map(|byte| byte as i8).collect();
    message_bytes.push(0);
    AddErrorMessage(Some(heap.slice_mut(error_buffer)?), Some(&message_bytes))?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
#[allow(non_snake_case)]
pub(crate) fn Extract0DParities(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    stereo0d: SourceMutPointer<inchi_Stereo0D>,
    num_stereo0d: i32,
    error_buffer: Option<SourceMutPointer<i8>>,
    mut error: Option<&mut i32>,
    unknown_parity: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:136 Extract0DParities
    // INCHI✔️❌: int Extract0DParities( inp_ATOM *at,
    // INCHI✔️❌:                        int nNumAtoms,
    // INCHI✔️❌:                        inchi_Stereo0D *stereo0D,
    // INCHI✔️❌:                        int num_stereo0D,
    // INCHI✔️❌:                        char *pStrErr,
    // INCHI✔️❌:                        int *err,
    // INCHI✔️❌:                        int vABParityUnknown )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         vABParityUnknown holds actual value of an internal constant signifying
    // INCHI✔️❌:         unknown parity: either the same as for undefined parity (default==standard)
    // INCHI✔️❌:         or a specific one (non-std; requested by SLUUD switch).
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (stereo0D && num_stereo0D > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i0D, a2, k, k_prev, type, j, j1, j2, len, parity, parityNM;
    // INCHI✔️❌:         int sb_ord_from_i1, sb_ord_from_i2, sn_ord_from_i1, sn_ord_from_i2;
    // INCHI✔️❌:         AT_NUMB i1n, i2n, i1, i2;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i0D = 0; i0D < num_stereo0D; i0D++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             parity = ( stereo0D[i0D].parity & SB_PARITY_MASK );
    // INCHI✔️❌:             parityNM = ( stereo0D[i0D].parity & SB_PARITY_FLAG ) >> SB_PARITY_SHFT;
    // INCHI✔️❌:
    // INCHI✔️❌:             if (parity == INCHI_PARITY_NONE ||
    // INCHI✔️❌:                  (parity != INCHI_PARITY_ODD && parity != INCHI_PARITY_EVEN &&
    // INCHI✔️❌:                  parity != INCHI_PARITY_UNKNOWN && parity != INCHI_PARITY_UNDEFINED)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 char szTemp[16];
    // INCHI✔️❌:                 sprintf(szTemp, "#%d", i0D + 1);
    // INCHI✔️❌:                 TREAT_ERR( *err, 0, "Wrong 0D stereo descriptor(s):" );
    // INCHI✔️❌:                 TREAT_ERR( *err, 0, szTemp );
    // INCHI✔️❌:                 continue; /* warning */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             type = stereo0D[i0D].type;
    // INCHI✔️❌:             a2 = stereo0D[i0D].central_atom; /* central atom or -1 */
    // INCHI✔️❌:             j = -1;
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             sb_ord_from_i1 = sb_ord_from_i2 = sn_ord_from_i1 = sn_ord_from_i2 = -1;
    // INCHI✔️❌:             i1n = i2n = i1 = i2 = MAX_ATOMS + 1;
    // INCHI✔️❌:
    // INCHI✔️❌:             if (( type == INCHI_StereoType_Tetrahedral ||
    // INCHI✔️❌:                 ((type == INCHI_StereoType_Allene ) &&
    // INCHI✔️❌:                   0 <= a2 && a2 < nNumAtoms)) ||
    // INCHI✔️❌:                   (type == INCHI_StereoType_DoubleBond &&
    // INCHI✔️❌:                   a2 == NO_ATOM)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* test the quadruplet */
    // INCHI✔️❌:                 for (j = 0, k_prev = -1; j < 4; j++, k_prev = k)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     k = stereo0D[i0D].neighbor[j];
    // INCHI✔️❌:                     if (k < 0 || k >= nNumAtoms || k_prev == k)
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     /* tetrahedral atom connectivity test */
    // INCHI✔️❌:                     if (type == INCHI_StereoType_Tetrahedral &&
    // INCHI✔️❌:                          k != a2 &&
    // INCHI✔️❌:                          !is_in_the_list( at[a2].neighbor, (AT_NUMB) k, at[a2].valence ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* Double bond, Cumulene and allene are tested in the next if() */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* Find in the adjacency lists the double bond neighbor that leads to the opposite atom */
    // INCHI✔️❌:             if (j == 4 && ( type == INCHI_StereoType_Allene ||
    // INCHI✔️❌:                 type == INCHI_StereoType_DoubleBond ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 AT_NUMB *p1 = NULL, *p2 = NULL, *q1 = NULL, *q2 = NULL;
    // INCHI✔️❌:                 i1n = (AT_NUMB) stereo0D[i0D].neighbor[0];
    // INCHI✔️❌:                 i1 = (AT_NUMB) stereo0D[i0D].neighbor[1];
    // INCHI✔️❌:                 i2 = (AT_NUMB) stereo0D[i0D].neighbor[2];
    // INCHI✔️❌:                 i2n = (AT_NUMB) stereo0D[i0D].neighbor[3];
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* find q1 and q2 */
    // INCHI✔️❌:                 if (!( q1 = is_in_the_list( at[i1].neighbor, i1n, at[i1].valence ) ) ||
    // INCHI✔️❌:                      !( q2 = is_in_the_list( at[i2].neighbor, i2n, at[i2].valence ) ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j = -2; /* error flag */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* allene or cumulene; follow double bonds from i1 to i2 */
    // INCHI✔️❌:                     if (!( p1 = is_in_the_list( at[i1].neighbor, i2, at[i1].valence ) ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* at[i1] and at[i2] are not connected: can be only allene or cumulene */
    // INCHI✔️❌:
    // INCHI✔️❌:                         AT_NUMB prev, cur, next;
    // INCHI✔️❌:                         int     num_dbond, i, next_ord, half_len;
    // INCHI✔️❌:
    // INCHI✔️❌:                         cur = next = i1;
    // INCHI✔️❌:                         len = half_len = 0;
    // INCHI✔️❌:                         while (len < 20)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* arbitrary very high upper limit to prevent infinite loop */
    // INCHI✔️❌:                             prev = cur;
    // INCHI✔️❌:                             cur = next;
    // INCHI✔️❌:
    // INCHI✔️❌:                             for (i = 0, num_dbond = 0; i < at[cur].valence; i++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* follow double bond path && avoid going back */
    // INCHI✔️❌:                                 if (at[cur].bond_type[i] == BOND_TYPE_DOUBLE &&
    // INCHI✔️❌:                                      prev != at[cur].neighbor[i])
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     next = at[cur].neighbor[i];
    // INCHI✔️❌:                                     next_ord = i;
    // INCHI✔️❌:                                     num_dbond++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (num_dbond == 1 && next != i1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 len++;
    // INCHI✔️❌:                                 if (len == 1)
    // INCHI✔️❌:                                     sb_ord_from_i1 = next_ord;
    // INCHI✔️❌:
    // INCHI✔️❌:                                 if (type == INCHI_StereoType_Allene && next == (AT_NUMB) a2)
    // INCHI✔️❌:                                     half_len = len;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (cur == i2 && prev != cur && 0 == num_dbond && len > 1 &&
    // INCHI✔️❌:                             ( p2 = is_in_the_list( at[i2].neighbor, prev, at[i2].valence ) ) &&
    // INCHI✔️❌:                             ( type != INCHI_StereoType_Allene || len == 2 * half_len ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             sb_ord_from_i2 = p2 - at[i2].neighbor;
    // INCHI✔️❌:                             sn_ord_from_i1 = q1 - at[i1].neighbor;
    // INCHI✔️❌:                             sn_ord_from_i2 = q2 - at[i2].neighbor;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             j = -5; /* error flag */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* allene must have been already processed, otherwise error */
    // INCHI✔️❌:                         if (type == INCHI_StereoType_Allene)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* error: atoms #1 and #2 of allene are connected */
    // INCHI✔️❌:                             j = -3; /* error flag */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* double bond only; the bond type is not checked because at the end
    // INCHI✔️❌:                                of the normalization it may happen to be alternating */
    // INCHI✔️❌:                             if (type == INCHI_StereoType_DoubleBond &&
    // INCHI✔️❌:                                 ( p2 = is_in_the_list( at[i2].neighbor, i1, at[i2].valence ) ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 sb_ord_from_i1 = p1 - at[i1].neighbor;
    // INCHI✔️❌:                                 sb_ord_from_i2 = p2 - at[i2].neighbor;
    // INCHI✔️❌:                                 sn_ord_from_i1 = q1 - at[i1].neighbor;
    // INCHI✔️❌:                                 sn_ord_from_i2 = q2 - at[i2].neighbor;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 j = -4; /* error flag */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (j != 4)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 char szTemp[16];
    // INCHI✔️❌:                 sprintf(szTemp, "#%d", i0D + 1);
    // INCHI✔️❌:                 TREAT_ERR( *err, 0, "Wrong 0D stereo descriptor(s):" );
    // INCHI✔️❌:                 TREAT_ERR( *err, 0, szTemp );
    // INCHI✔️❌:                 continue; /* error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             switch (type)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case INCHI_StereoType_None:
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 case INCHI_StereoType_DoubleBond:
    // INCHI✔️❌:                 case INCHI_StereoType_Allene:
    // INCHI✔️❌:                     for (j1 = 0; j1 < MAX_NUM_STEREO_BONDS && at[i1].sb_parity[j1]; j1++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     for (j2 = 0; j2 < MAX_NUM_STEREO_BONDS && at[i2].sb_parity[j2]; j2++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (j1 < MAX_NUM_STEREO_BONDS && j2 < MAX_NUM_STEREO_BONDS &&
    // INCHI✔️❌:                          sb_ord_from_i1 >= 0 && sb_ord_from_i2 >= 0 &&
    // INCHI✔️❌:                          sn_ord_from_i1 >= 0 && sn_ord_from_i2 >= 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         switch (parity)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case INCHI_PARITY_ODD:
    // INCHI✔️❌:                                 at[i1].sb_parity[j1] = AB_PARITY_ODD;
    // INCHI✔️❌:                                 at[i2].sb_parity[j2] = AB_PARITY_EVEN;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case INCHI_PARITY_EVEN:
    // INCHI✔️❌:                                 at[i1].sb_parity[j1] = AB_PARITY_ODD;
    // INCHI✔️❌:                                 at[i2].sb_parity[j2] = AB_PARITY_ODD;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case INCHI_PARITY_UNDEFINED:
    // INCHI✔️❌:                                 at[i1].sb_parity[j1] = AB_PARITY_UNDF;
    // INCHI✔️❌:                                 at[i2].sb_parity[j2] = AB_PARITY_UNDF;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             default:
    // INCHI✔️❌:                                 if (parity == INCHI_PARITY_UNKNOWN)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     at[i1].sb_parity[j1] = vABParityUnknown;
    // INCHI✔️❌:                                     at[i2].sb_parity[j2] = vABParityUnknown;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     at[i1].sb_parity[j1] = AB_PARITY_NONE;
    // INCHI✔️❌:                                     at[i2].sb_parity[j2] = AB_PARITY_NONE;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         switch (parityNM)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case INCHI_PARITY_ODD:
    // INCHI✔️❌:                                 at[i1].sb_parity[j1] |= AB_PARITY_ODD << SB_PARITY_SHFT;
    // INCHI✔️❌:                                 at[i2].sb_parity[j2] |= AB_PARITY_EVEN << SB_PARITY_SHFT;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case INCHI_PARITY_EVEN:
    // INCHI✔️❌:                                 at[i1].sb_parity[j1] |= AB_PARITY_ODD << SB_PARITY_SHFT;
    // INCHI✔️❌:                                 at[i2].sb_parity[j2] |= AB_PARITY_ODD << SB_PARITY_SHFT;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case INCHI_PARITY_UNDEFINED:
    // INCHI✔️❌:                                 at[i1].sb_parity[j1] |= AB_PARITY_UNDF << SB_PARITY_SHFT;
    // INCHI✔️❌:                                 at[i2].sb_parity[j2] |= AB_PARITY_UNDF << SB_PARITY_SHFT;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             default:
    // INCHI✔️❌:                                 if (parityNM == INCHI_PARITY_UNKNOWN)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     at[i1].sb_parity[j1] |= vABParityUnknown << SB_PARITY_SHFT;
    // INCHI✔️❌:                                     at[i2].sb_parity[j2] |= vABParityUnknown << SB_PARITY_SHFT;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         at[i1].sb_ord[j1] = sb_ord_from_i1;
    // INCHI✔️❌:                         at[i1].sn_ord[j1] = sn_ord_from_i1;
    // INCHI✔️❌:                         at[i1].sn_orig_at_num[j1] = at[i1n].orig_at_number;
    // INCHI✔️❌:
    // INCHI✔️❌:                         at[i2].sb_ord[j2] = sb_ord_from_i2;
    // INCHI✔️❌:                         at[i2].sn_ord[j2] = sn_ord_from_i2;
    // INCHI✔️❌:                         at[i2].sn_orig_at_num[j2] = at[i2n].orig_at_number;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case INCHI_StereoType_Tetrahedral:
    // INCHI✔️❌:                     switch (parity)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         case INCHI_PARITY_ODD:
    // INCHI✔️❌:                             at[a2].p_parity = AB_PARITY_ODD;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case INCHI_PARITY_EVEN:
    // INCHI✔️❌:                             at[a2].p_parity = AB_PARITY_EVEN;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         case INCHI_PARITY_UNDEFINED:
    // INCHI✔️❌:                             at[a2].p_parity = AB_PARITY_UNDF;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         default:
    // INCHI✔️❌:                             if (parity == INCHI_PARITY_UNKNOWN)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 at[a2].p_parity = vABParityUnknown;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     for (j = 0; j < 4; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = stereo0D[i0D].neighbor[j];
    // INCHI✔️❌:                         at[a2].p_orig_at_num[j] = at[k].orig_at_number;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:
    // INCHI✔️❌:                 default:
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* take care of Unknown stereobonds:                                     */
    // INCHI✔️❌:         /* copy their Unknown stereo descriptors to at->bond_stereo (2005-03-01) */
    // INCHI✔️❌:         /* Note: to this stage, unk/undef set to what was requested              */
    // INCHI✔️❌:         /*( through vABParityUnknown )  (2009-12-12)                             */
    // INCHI✔️❌:         FixUnkn0DStereoBonds( at, nNumAtoms );
    // INCHI✔️❌:
    // INCHI✔️❌:         if ((k = ReconcileAllCmlBondParities( at, nNumAtoms, 0 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char szErrCode[16];
    // INCHI✔️❌:             sprintf( szErrCode, "%d", k );
    // INCHI✔️❌:             AddErrorMessage( pStrErr, "0D Parities Reconciliation failed:" );
    // INCHI✔️❌:             AddErrorMessage( pStrErr, szErrCode );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Extract0DParities

    if stereo0d.is_null() || num_stereo0d <= 0 {
        return Ok(0);
    }

    for descriptor_index in 0..num_stereo0d {
        let descriptor_index_usize =
            usize::try_from(descriptor_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let descriptor = heap
            .slice(stereo0d.as_const())?
            .get(descriptor_index_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let parity = i32::from(descriptor.parity) & SB_PARITY_MASK as i32;
        let parity_nm = (i32::from(descriptor.parity) & SB_PARITY_FLAG as i32) >> SB_PARITY_SHFT;
        if parity == tagINCHIStereoParity0D_INCHI_PARITY_NONE as i32
            || (parity != tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32
                && parity != tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32
                && parity != tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN as i32
                && parity != tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i32)
        {
            let _ = error.as_deref_mut().ok_or(SourceHeapError::NullPointer)?;
            add_extract_0d_error(heap, error_buffer, "Wrong 0D stereo descriptor(s):")?;
            add_extract_0d_error(heap, error_buffer, &format!("#{}", descriptor_index + 1))?;
            continue;
        }

        let stereo_type = i32::from(descriptor.type_);
        let central_atom = i32::from(descriptor.central_atom);
        let mut j = -1_i32;
        let (mut sb_ord_from_i1, mut sb_ord_from_i2) = (-1_i32, -1_i32);
        let (mut sn_ord_from_i1, mut sn_ord_from_i2) = (-1_i32, -1_i32);
        let (mut i1n, mut i2n, mut i1, mut i2) = (32767_u16, 32767_u16, 32767_u16, 32767_u16);

        if stereo_type == tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i32
            || (stereo_type == tagINCHIStereoType0D_INCHI_StereoType_Allene as i32
                && 0 <= central_atom
                && central_atom < num_atoms)
            || (stereo_type == tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i32
                && central_atom == NO_ATOM)
        {
            let mut previous = -1_i32;
            j = 0;
            while j < 4 {
                let k = i32::from(descriptor.neighbor[j as usize]);
                if k < 0 || k >= num_atoms || previous == k {
                    break;
                }
                if stereo_type == tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i32
                    && k != central_atom
                {
                    let center_index = usize::try_from(central_atom)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let center = heap
                        .slice(atoms.as_const())?
                        .get(center_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if is_in_the_list(Some(&center.neighbor), k as u16, i32::from(center.valence))?
                        .is_none()
                    {
                        break;
                    }
                }
                previous = k;
                j += 1;
            }
        }

        if j == 4
            && (stereo_type == tagINCHIStereoType0D_INCHI_StereoType_Allene as i32
                || stereo_type == tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i32)
        {
            i1n = descriptor.neighbor[0] as u16;
            i1 = descriptor.neighbor[1] as u16;
            i2 = descriptor.neighbor[2] as u16;
            i2n = descriptor.neighbor[3] as u16;
            let i1_index = usize::from(i1);
            let i2_index = usize::from(i2);
            let atom_values = heap.slice(atoms.as_const())?;
            let atom_i1 = atom_values
                .get(i1_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom_i2 = atom_values
                .get(i2_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let q1 = is_in_the_list(Some(&atom_i1.neighbor), i1n, i32::from(atom_i1.valence))?;
            let q2 = is_in_the_list(Some(&atom_i2.neighbor), i2n, i32::from(atom_i2.valence))?;
            if q1.is_none() || q2.is_none() {
                j = -2;
            } else {
                let p1 = is_in_the_list(Some(&atom_i1.neighbor), i2, i32::from(atom_i1.valence))?;
                if p1.is_none() {
                    let mut current = i1;
                    let mut next = i1;
                    let mut previous = i1;
                    let mut num_double_bonds = 0_i32;
                    let mut next_ordinal = 0_i32;
                    let mut length = 0_i32;
                    let mut half_length = 0_i32;
                    while length < 20 {
                        previous = current;
                        current = next;
                        num_double_bonds = 0;
                        let current_atom = heap
                            .slice(atoms.as_const())?
                            .get(usize::from(current))
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        for ordinal in 0..i32::from(current_atom.valence) {
                            let ordinal_usize = usize::try_from(ordinal)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            if current_atom.bond_type[ordinal_usize] == BOND_TYPE_DOUBLE as u8
                                && previous != current_atom.neighbor[ordinal_usize]
                            {
                                next = current_atom.neighbor[ordinal_usize];
                                next_ordinal = ordinal;
                                num_double_bonds += 1;
                            }
                        }
                        if num_double_bonds == 1 && next != i1 {
                            length += 1;
                            if length == 1 {
                                sb_ord_from_i1 = next_ordinal;
                            }
                            if stereo_type == tagINCHIStereoType0D_INCHI_StereoType_Allene as i32
                                && next == central_atom as u16
                            {
                                half_length = length;
                            }
                        } else {
                            break;
                        }
                    }
                    let final_atom = heap
                        .slice(atoms.as_const())?
                        .get(i2_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let p2 = is_in_the_list(
                        Some(&final_atom.neighbor),
                        previous,
                        i32::from(final_atom.valence),
                    )?;
                    if current == i2
                        && previous != current
                        && num_double_bonds == 0
                        && length > 1
                        && p2.is_some()
                        && (stereo_type != tagINCHIStereoType0D_INCHI_StereoType_Allene as i32
                            || length == 2 * half_length)
                    {
                        sb_ord_from_i2 = p2.unwrap() as i32;
                        sn_ord_from_i1 = q1.unwrap() as i32;
                        sn_ord_from_i2 = q2.unwrap() as i32;
                    } else {
                        j = -5;
                    }
                } else if stereo_type == tagINCHIStereoType0D_INCHI_StereoType_Allene as i32 {
                    j = -3;
                } else {
                    let p2 =
                        is_in_the_list(Some(&atom_i2.neighbor), i1, i32::from(atom_i2.valence))?;
                    if stereo_type == tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i32
                        && p2.is_some()
                    {
                        sb_ord_from_i1 = p1.unwrap() as i32;
                        sb_ord_from_i2 = p2.unwrap() as i32;
                        sn_ord_from_i1 = q1.unwrap() as i32;
                        sn_ord_from_i2 = q2.unwrap() as i32;
                    } else {
                        j = -4;
                    }
                }
            }
        }

        if j != 4 {
            let _ = error.as_deref_mut().ok_or(SourceHeapError::NullPointer)?;
            add_extract_0d_error(heap, error_buffer, "Wrong 0D stereo descriptor(s):")?;
            add_extract_0d_error(heap, error_buffer, &format!("#{}", descriptor_index + 1))?;
            continue;
        }

        match stereo_type as u32 {
            tagINCHIStereoType0D_INCHI_StereoType_None => continue,
            tagINCHIStereoType0D_INCHI_StereoType_DoubleBond
            | tagINCHIStereoType0D_INCHI_StereoType_Allene => {
                let i1_index = usize::from(i1);
                let i2_index = usize::from(i2);
                let j1 = heap.slice(atoms.as_const())?[i1_index]
                    .sb_parity
                    .iter()
                    .take(MAX_NUM_STEREO_BONDS as usize)
                    .position(|value| *value == 0)
                    .unwrap_or(MAX_NUM_STEREO_BONDS as usize);
                let j2 = heap.slice(atoms.as_const())?[i2_index]
                    .sb_parity
                    .iter()
                    .take(MAX_NUM_STEREO_BONDS as usize)
                    .position(|value| *value == 0)
                    .unwrap_or(MAX_NUM_STEREO_BONDS as usize);
                if j1 < MAX_NUM_STEREO_BONDS as usize
                    && j2 < MAX_NUM_STEREO_BONDS as usize
                    && sb_ord_from_i1 >= 0
                    && sb_ord_from_i2 >= 0
                    && sn_ord_from_i1 >= 0
                    && sn_ord_from_i2 >= 0
                {
                    let (left, right) = match parity as u32 {
                        tagINCHIStereoParity0D_INCHI_PARITY_ODD => {
                            (AB_PARITY_ODD as i8, AB_PARITY_EVEN as i8)
                        }
                        tagINCHIStereoParity0D_INCHI_PARITY_EVEN => {
                            (AB_PARITY_ODD as i8, AB_PARITY_ODD as i8)
                        }
                        tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED => {
                            (AB_PARITY_UNDF as i8, AB_PARITY_UNDF as i8)
                        }
                        tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN => {
                            (unknown_parity as i8, unknown_parity as i8)
                        }
                        _ => (AB_PARITY_NONE as i8, AB_PARITY_NONE as i8),
                    };
                    heap.slice_mut(atoms)?[i1_index].sb_parity[j1] = left;
                    heap.slice_mut(atoms)?[i2_index].sb_parity[j2] = right;
                    let (left_nm, right_nm) = match parity_nm as u32 {
                        tagINCHIStereoParity0D_INCHI_PARITY_ODD => (
                            (AB_PARITY_ODD << SB_PARITY_SHFT) as i8,
                            (AB_PARITY_EVEN << SB_PARITY_SHFT) as i8,
                        ),
                        tagINCHIStereoParity0D_INCHI_PARITY_EVEN => (
                            (AB_PARITY_ODD << SB_PARITY_SHFT) as i8,
                            (AB_PARITY_ODD << SB_PARITY_SHFT) as i8,
                        ),
                        tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED => (
                            (AB_PARITY_UNDF << SB_PARITY_SHFT) as i8,
                            (AB_PARITY_UNDF << SB_PARITY_SHFT) as i8,
                        ),
                        tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN => (
                            (unknown_parity << SB_PARITY_SHFT) as i8,
                            (unknown_parity << SB_PARITY_SHFT) as i8,
                        ),
                        _ => (0, 0),
                    };
                    heap.slice_mut(atoms)?[i1_index].sb_parity[j1] |= left_nm;
                    heap.slice_mut(atoms)?[i2_index].sb_parity[j2] |= right_nm;
                    heap.slice_mut(atoms)?[i1_index].sb_ord[j1] = sb_ord_from_i1 as i8;
                    heap.slice_mut(atoms)?[i1_index].sn_ord[j1] = sn_ord_from_i1 as i8;
                    let left_original =
                        heap.slice(atoms.as_const())?[usize::from(i1n)].orig_at_number;
                    heap.slice_mut(atoms)?[i1_index].sn_orig_at_num[j1] = left_original;
                    heap.slice_mut(atoms)?[i2_index].sb_ord[j2] = sb_ord_from_i2 as i8;
                    heap.slice_mut(atoms)?[i2_index].sn_ord[j2] = sn_ord_from_i2 as i8;
                    let right_original =
                        heap.slice(atoms.as_const())?[usize::from(i2n)].orig_at_number;
                    heap.slice_mut(atoms)?[i2_index].sn_orig_at_num[j2] = right_original;
                }
            }
            tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral => {
                let center_index = usize::try_from(central_atom)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_parity = match parity as u32 {
                    tagINCHIStereoParity0D_INCHI_PARITY_ODD => AB_PARITY_ODD as i8,
                    tagINCHIStereoParity0D_INCHI_PARITY_EVEN => AB_PARITY_EVEN as i8,
                    tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED => AB_PARITY_UNDF as i8,
                    tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN => unknown_parity as i8,
                    _ => continue,
                };
                heap.slice_mut(atoms)?[center_index].p_parity = atom_parity;
                for ordinal in 0..4 {
                    let neighbor_index = usize::try_from(descriptor.neighbor[ordinal])
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let original = heap.slice(atoms.as_const())?[neighbor_index].orig_at_number;
                    heap.slice_mut(atoms)?[center_index].p_orig_at_num[ordinal] = original;
                }
            }
            _ => {}
        }
    }

    FixUnkn0DStereoBonds(heap, atoms, num_atoms)?;
    let reconciliation = ReconcileAllCmlBondParities(heap, atoms, num_atoms, 0)?;
    if reconciliation != 0 {
        add_extract_0d_error(heap, error_buffer, "0D Parities Reconciliation failed:")?;
        add_extract_0d_error(heap, error_buffer, &reconciliation.to_string())?;
    }
    Ok(0)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn FindToken(
    heap: &mut SourceHeap,
    input: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
    token: SourceConstPointer<i8>,
    token_length: i32,
    line: SourceMutPointer<i8>,
    line_capacity: i32,
    mut cursor: SourceMutPointer<i8>,
    result_length: &mut i32,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:445 FindToken
    // INCHI✔️❌: char* FindToken( INCHI_IOSTREAM *inp_file,
    // INCHI✔️❌:                  int *bTooLongLine,
    // INCHI✔️❌:                  const char *sToken,
    // INCHI✔️❌:                  int lToken,
    // INCHI✔️❌:                  char *szLine,
    // INCHI✔️❌:                  int nLenLine,
    // INCHI✔️❌:                  char *p,
    // INCHI✔️❌:                  int *res )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char *q;
    // INCHI✔️❌:     int   res2;
    // INCHI✔️❌:
    // INCHI✔️❌:     while (!( q = strstr( p, sToken ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( q = strrchr( p, '/' ) ) && ( q + lToken > szLine + *res ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *res -= q - szLine; /* res = the length of the szLine to be left in */
    // INCHI✔️❌:             memmove(szLine, q, (long long)*res + 1); /* djb-rwth: cast operator added */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *res = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         res2 = inchi_ios_getsTab1( szLine + *res, nLenLine - *res - 1,
    // INCHI✔️❌:                                    inp_file, bTooLongLine );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!*bTooLongLine || 0 > res2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* the line is over or end of file */
    // INCHI✔️❌:             return NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *res += res2;
    // INCHI✔️❌:             p = szLine;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return q + lToken;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FindToken

    loop {
        let (token_offset, last_slash_offset) = {
            let haystack = heap.slice(cursor.as_const())?;
            let haystack_length = haystack
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            let haystack = &haystack[..haystack_length];
            let needle = heap.slice(token)?;
            let needle_length = needle
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            let needle = &needle[..needle_length];
            let token_offset = if needle.is_empty() {
                Some(0)
            } else {
                haystack
                    .windows(needle.len())
                    .position(|window| window == needle)
            };
            (
                token_offset,
                haystack.iter().rposition(|byte| *byte == b'/' as i8),
            )
        };

        if let Some(token_offset) = token_offset {
            let token_position = cursor.offset(
                i64::try_from(token_offset)
                    .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?,
            )?;
            return token_position.offset(i64::from(token_length));
        }

        if let Some(last_slash_offset) = last_slash_offset {
            let slash = cursor.offset(
                i64::try_from(last_slash_offset)
                    .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?,
            )?;
            let slash_from_line = slash.difference(line)?;
            let slash_plus_token = slash_from_line
                .checked_add(i64::from(token_length))
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            if slash_plus_token > i64::from(*result_length) {
                let slash_from_line = i32::try_from(slash_from_line)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                *result_length = result_length
                    .checked_sub(slash_from_line)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let source_start = usize::try_from(slash_from_line)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let count = usize::try_from(*result_length)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let source_end = source_start
                    .checked_add(count)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let bytes = heap.slice_mut(line)?;
                if source_end > bytes.len() || count > bytes.len() {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                bytes.copy_within(source_start..source_end, 0);
            } else {
                *result_length = 0;
            }
        } else {
            *result_length = 0;
        }

        let destination = line.offset(i64::from(*result_length))?;
        let read_capacity = line_capacity
            .checked_sub(*result_length)
            .and_then(|value| value.checked_sub(1))
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let result = inchi_ios_getsTab1(heap, destination, read_capacity, input, too_long_line)?;
        if *too_long_line == 0 || result < 0 {
            return Ok(SourceMutPointer::null());
        }
        *result_length = result_length
            .checked_add(result)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        cursor = line;
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn LoadLine(
    heap: &mut SourceHeap,
    input: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
    item_is_over: &mut i32,
    slash_slot: &mut SourceMutPointer<i8>,
    line: SourceMutPointer<i8>,
    line_capacity: i32,
    minimum_load_capacity: i32,
    mut cursor: SourceMutPointer<i8>,
    result_length: &mut i32,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:489 LoadLine
    // INCHI✔️❌: char *LoadLine( INCHI_IOSTREAM *inp_file,
    // INCHI✔️❌:                 int *bTooLongLine,
    // INCHI✔️❌:                 int *bItemIsOver,
    // INCHI✔️❌:                 char **s,
    // INCHI✔️❌:                 char *szLine,
    // INCHI✔️❌:                 int nLenLine,
    // INCHI✔️❌:                 int nMinLen2Load,
    // INCHI✔️❌:                 char *p,
    // INCHI✔️❌:                 int *res )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     int pos = p - szLine, res2;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!*bItemIsOver && nLenLine - ( *res - pos ) > nMinLen2Load)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* load the next portion if possible */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (pos)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *res -= pos;
    // INCHI✔️❌:             memmove(szLine, p, (long long)*res + 1); /* djb-rwth: cast operator added */
    // INCHI✔️❌:             p = szLine;
    // INCHI✔️❌:             if (*s)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *s -= pos;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         res2 = inchi_ios_getsTab1( szLine + *res,
    // INCHI✔️❌:                                    nLenLine - *res - 1,
    // INCHI✔️❌:                                    inp_file, bTooLongLine );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (res2 > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *bItemIsOver = ( ( *s = strchr( p + *res, '/' ) ) || !*bTooLongLine );
    // INCHI✔️❌:             *res += res2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *bItemIsOver = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return p;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: LoadLine

    let position = i32::try_from(cursor.difference(line)?)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let retained = result_length
        .checked_sub(position)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let available = line_capacity
        .checked_sub(retained)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if *item_is_over == 0 && available > minimum_load_capacity {
        if position != 0 {
            *result_length = retained;
            let source_start =
                usize::try_from(position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let count = usize::try_from(*result_length)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let source_end = source_start
                .checked_add(count)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let bytes = heap.slice_mut(line)?;
            if source_end > bytes.len() || count > bytes.len() {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            bytes.copy_within(source_start..source_end, 0);
            cursor = line;
            if !slash_slot.is_null() {
                *slash_slot = slash_slot.offset(-i64::from(position))?;
            }
        }

        let destination = line.offset(i64::from(*result_length))?;
        let read_capacity = line_capacity
            .checked_sub(*result_length)
            .and_then(|value| value.checked_sub(1))
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let result = inchi_ios_getsTab1(heap, destination, read_capacity, input, too_long_line)?;
        if result > 0 {
            let appended = cursor.offset(i64::from(*result_length))?;
            let slash_offset = {
                let bytes = heap.slice(appended.as_const())?;
                let length = bytes
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                bytes[..length].iter().position(|byte| *byte == b'/' as i8)
            };
            *slash_slot = if let Some(slash_offset) = slash_offset {
                appended.offset(
                    i64::try_from(slash_offset)
                        .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?,
                )?
            } else {
                SourceMutPointer::null()
            };
            *item_is_over = i32::from(!slash_slot.is_null() || *too_long_line == 0);
            *result_length = result_length
                .checked_add(result)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        } else {
            *item_is_over = 1;
        }
    }

    Ok(cursor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{AB_PARITY_UNKN, INCHI_IOS_STRING, INCHI_IOS_TYPE_STRING};
    use crate::test_support::allocate_source_fixture;

    #[test]
    fn source_port__readinch__findtoken__line_445() {
        let mut heap = SourceHeap::default();
        let token = allocate_source_fixture(
            &mut heap,
            b"InChI="
                .iter()
                .map(|byte| *byte as i8)
                .chain([0])
                .collect(),
        );

        let line = allocate_source_fixture(
            &mut heap,
            b"xx/InChI=abc"
                .iter()
                .map(|byte| *byte as i8)
                .chain([0])
                .collect(),
        );
        let empty_input = allocate_source_fixture(&mut heap, Vec::<i8>::new());
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: empty_input,
                nAllocatedLength: 0,
                nUsedLength: 0,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let mut too_long = 0;
        let mut length = 12;
        let found = FindToken(
            &mut heap,
            &mut ios,
            &mut too_long,
            token.as_const(),
            6,
            line,
            13,
            line,
            &mut length,
        )
        .unwrap();
        assert_eq!(found.difference(line), Ok(9));
        assert_eq!(
            &heap.slice(found.as_const()).unwrap()[..4],
            &[b'a' as i8, b'b' as i8, b'c' as i8, 0]
        );
        assert_eq!(ios.s.nPtr, 0);
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
        assert_eq!(inchi_free(&mut heap, empty_input), Ok(()));

        let line = allocate_source_fixture(&mut heap, {
            let mut bytes = vec![0_i8; 16];
            bytes[..5]
                .copy_from_slice(&[b'/' as i8, b'I' as i8, b'n' as i8, b'C' as i8, b'h' as i8]);
            bytes
        });
        let continuation = allocate_source_fixture(
            &mut heap,
            b"I=restXX".iter().map(|byte| *byte as i8).collect(),
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: continuation,
                nAllocatedLength: 8,
                nUsedLength: 8,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let mut too_long = 1;
        let mut length = 5;
        let found = FindToken(
            &mut heap,
            &mut ios,
            &mut too_long,
            token.as_const(),
            6,
            line,
            16,
            line,
            &mut length,
        )
        .unwrap();
        assert_eq!(found.difference(line), Ok(7));
        assert_eq!(length, 13);
        assert_eq!(too_long, 1);
        assert_eq!(ios.s.nPtr, 8);
        assert_eq!(
            &heap.slice(found.as_const()).unwrap()[..7],
            &[
                b'r' as i8, b'e' as i8, b's' as i8, b't' as i8, b'X' as i8, b'X' as i8, 0,
            ]
        );
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
        assert_eq!(inchi_free(&mut heap, continuation), Ok(()));
        assert_eq!(inchi_free(&mut heap, token), Ok(()));
    }

    #[test]
    fn source_port__readinch__loadline__line_489() {
        let mut heap = SourceHeap::default();

        let line = allocate_source_fixture(&mut heap, {
            let mut bytes = vec![0_i8; 20];
            bytes[..5]
                .copy_from_slice(&[b'x' as i8, b'x' as i8, b'a' as i8, b'b' as i8, b'c' as i8]);
            bytes
        });
        let empty_input = allocate_source_fixture(&mut heap, Vec::<i8>::new());
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: empty_input,
                nAllocatedLength: 0,
                nUsedLength: 0,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let mut too_long = 1;
        let mut item_is_over = 0;
        let mut slash = line.offset(4).unwrap();
        let mut length = 5;
        let cursor = LoadLine(
            &mut heap,
            &mut ios,
            &mut too_long,
            &mut item_is_over,
            &mut slash,
            line,
            20,
            2,
            line.offset(2).unwrap(),
            &mut length,
        )
        .unwrap();
        assert_eq!(cursor, line);
        assert_eq!(length, 3);
        assert_eq!(slash.difference(line), Ok(2));
        assert_eq!(item_is_over, 1);
        assert_eq!(too_long, 0);
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..4],
            &[b'a' as i8, b'b' as i8, b'c' as i8, 0]
        );
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
        assert_eq!(inchi_free(&mut heap, empty_input), Ok(()));

        let line = allocate_source_fixture(&mut heap, {
            let mut bytes = vec![0_i8; 20];
            bytes[..5]
                .copy_from_slice(&[b'x' as i8, b'x' as i8, b'a' as i8, b'b' as i8, b'c' as i8]);
            bytes
        });
        let continuation = allocate_source_fixture(
            &mut heap,
            b"de/f\t".iter().map(|byte| *byte as i8).collect(),
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: continuation,
                nAllocatedLength: 5,
                nUsedLength: 5,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let mut too_long = 1;
        let mut item_is_over = 0;
        let mut slash = SourceMutPointer::null();
        let mut length = 5;
        let cursor = LoadLine(
            &mut heap,
            &mut ios,
            &mut too_long,
            &mut item_is_over,
            &mut slash,
            line,
            20,
            2,
            line.offset(2).unwrap(),
            &mut length,
        )
        .unwrap();
        assert_eq!(cursor, line);
        assert_eq!(length, 7);
        assert_eq!(slash.difference(line), Ok(5));
        assert_eq!(item_is_over, 1);
        assert_eq!(too_long, 0);
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..8],
            &[
                b'a' as i8, b'b' as i8, b'c' as i8, b'd' as i8, b'e' as i8, b'/' as i8, b'f' as i8,
                0,
            ]
        );
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
        assert_eq!(inchi_free(&mut heap, continuation), Ok(()));
    }

    fn extract_error_text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> String {
        let bytes = heap.slice(pointer.as_const()).unwrap();
        let length = bytes.iter().position(|byte| *byte == 0).unwrap();
        String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect()).unwrap()
    }

    fn tetrahedral_atoms() -> Vec<inp_ATOM> {
        let mut atoms = vec![inp_ATOM::default(); 4];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.orig_at_number = 10 + index as u16;
        }
        atoms[1].neighbor[..3].copy_from_slice(&[0, 2, 3]);
        atoms[1].valence = 3;
        atoms
    }

    fn direct_double_bond_atoms() -> Vec<inp_ATOM> {
        let mut atoms = vec![inp_ATOM::default(); 4];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.orig_at_number = 20 + index as u16;
        }
        atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        atoms[1].bond_type[..2].copy_from_slice(&[1, BOND_TYPE_DOUBLE as u8]);
        atoms[1].valence = 2;
        atoms[2].neighbor[..2].copy_from_slice(&[1, 3]);
        atoms[2].bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, 1]);
        atoms[2].valence = 2;
        atoms
    }

    #[test]
    fn source_port__readinch__extract0dparities__line_136() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            Extract0DParities(
                &mut heap,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                None,
                None,
                AB_PARITY_UNKN as i32,
            )
            .unwrap(),
            0
        );

        for (parity, expected) in [
            (tagINCHIStereoParity0D_INCHI_PARITY_ODD, AB_PARITY_ODD),
            (tagINCHIStereoParity0D_INCHI_PARITY_EVEN, AB_PARITY_EVEN),
            (
                tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED,
                AB_PARITY_UNDF,
            ),
            (tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN, AB_PARITY_UNKN),
        ] {
            let atoms = heap.allocate(tetrahedral_atoms()).unwrap();
            let descriptor = heap
                .allocate(vec![inchi_Stereo0D {
                    neighbor: [0, 1, 2, 3],
                    central_atom: 1,
                    type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
                    parity: parity as i8,
                }])
                .unwrap();
            let mut error = 7;
            Extract0DParities(
                &mut heap,
                atoms,
                4,
                descriptor,
                1,
                None,
                Some(&mut error),
                AB_PARITY_UNKN as i32,
            )
            .unwrap();
            let center = &heap.slice(atoms.as_const()).unwrap()[1];
            assert_eq!(center.p_parity, expected as i8);
            assert_eq!(center.p_orig_at_num, [10, 11, 12, 13]);
            assert_eq!(error, 7);
        }

        let invalid_atoms = heap.allocate(tetrahedral_atoms()).unwrap();
        let invalid = heap
            .allocate(vec![inchi_Stereo0D {
                neighbor: [0, 1, 2, 3],
                central_atom: 1,
                type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
                parity: tagINCHIStereoParity0D_INCHI_PARITY_NONE as i8,
            }])
            .unwrap();
        let error_buffer = heap.allocate(vec![0_i8; 256]).unwrap();
        let mut error = 9;
        Extract0DParities(
            &mut heap,
            invalid_atoms,
            4,
            invalid,
            1,
            Some(error_buffer),
            Some(&mut error),
            AB_PARITY_UNKN as i32,
        )
        .unwrap();
        assert_eq!(error, 9);
        assert_eq!(
            extract_error_text(&heap, error_buffer),
            "Wrong 0D stereo descriptor(s): #1"
        );

        let disconnected_atoms = heap.allocate(tetrahedral_atoms()).unwrap();
        let disconnected_descriptor = heap
            .allocate(vec![inchi_Stereo0D {
                neighbor: [0, 1, 2, 2],
                central_atom: 1,
                type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
                parity: tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8,
            }])
            .unwrap();
        let disconnected_error = heap.allocate(vec![0_i8; 256]).unwrap();
        Extract0DParities(
            &mut heap,
            disconnected_atoms,
            4,
            disconnected_descriptor,
            1,
            Some(disconnected_error),
            Some(&mut error),
            AB_PARITY_UNKN as i32,
        )
        .unwrap();
        assert_eq!(
            extract_error_text(&heap, disconnected_error),
            "Wrong 0D stereo descriptor(s): #1"
        );

        let double_atoms = heap.allocate(direct_double_bond_atoms()).unwrap();
        let double_descriptor = heap
            .allocate(vec![inchi_Stereo0D {
                neighbor: [0, 1, 2, 3],
                central_atom: NO_ATOM as i16,
                type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                parity: (tagINCHIStereoParity0D_INCHI_PARITY_ODD
                    | (tagINCHIStereoParity0D_INCHI_PARITY_EVEN << SB_PARITY_SHFT))
                    as i8,
            }])
            .unwrap();
        Extract0DParities(
            &mut heap,
            double_atoms,
            4,
            double_descriptor,
            1,
            None,
            Some(&mut error),
            AB_PARITY_UNKN as i32,
        )
        .unwrap();
        let double_values = heap.slice(double_atoms.as_const()).unwrap();
        assert_eq!(double_values[1].sb_parity[0], 9);
        assert_eq!(double_values[2].sb_parity[0], 10);
        assert_eq!(double_values[1].sb_ord[0], 1);
        assert_eq!(double_values[2].sb_ord[0], 0);
        assert_eq!(double_values[1].sn_ord[0], 0);
        assert_eq!(double_values[2].sn_ord[0], 1);
        assert_eq!(double_values[1].sn_orig_at_num[0], 20);
        assert_eq!(double_values[2].sn_orig_at_num[0], 23);

        let unknown_atoms = heap.allocate(direct_double_bond_atoms()).unwrap();
        let unknown_descriptor = heap
            .allocate(vec![inchi_Stereo0D {
                neighbor: [0, 1, 2, 3],
                central_atom: NO_ATOM as i16,
                type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                parity: tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN as i8,
            }])
            .unwrap();
        Extract0DParities(
            &mut heap,
            unknown_atoms,
            4,
            unknown_descriptor,
            1,
            None,
            Some(&mut error),
            AB_PARITY_UNKN as i32,
        )
        .unwrap();
        let unknown_values = heap.slice(unknown_atoms.as_const()).unwrap();
        assert_eq!(unknown_values[1].sb_parity[0], AB_PARITY_UNKN as i8);
        assert_eq!(unknown_values[2].sb_parity[0], AB_PARITY_UNKN as i8);
        assert_eq!(
            unknown_values[1].bond_stereo[1],
            crate::source_types::STEREO_DBLE_EITHER as i8
        );
        assert_eq!(
            unknown_values[2].bond_stereo[0],
            crate::source_types::STEREO_DBLE_EITHER as i8
        );

        let mut cumulene = vec![inp_ATOM::default(); 5];
        for (index, atom) in cumulene.iter_mut().enumerate() {
            atom.orig_at_number = 30 + index as u16;
        }
        cumulene[1].neighbor[..2].copy_from_slice(&[0, 2]);
        cumulene[1].bond_type[..2].copy_from_slice(&[1, BOND_TYPE_DOUBLE as u8]);
        cumulene[1].valence = 2;
        cumulene[2].neighbor[..2].copy_from_slice(&[1, 3]);
        cumulene[2].bond_type[..2]
            .copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_TYPE_DOUBLE as u8]);
        cumulene[2].valence = 2;
        cumulene[3].neighbor[..2].copy_from_slice(&[2, 4]);
        cumulene[3].bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, 1]);
        cumulene[3].valence = 2;
        let cumulene = heap.allocate(cumulene).unwrap();
        let allene = heap
            .allocate(vec![inchi_Stereo0D {
                neighbor: [0, 1, 3, 4],
                central_atom: 2,
                type_: tagINCHIStereoType0D_INCHI_StereoType_Allene as i8,
                parity: tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8,
            }])
            .unwrap();
        Extract0DParities(
            &mut heap,
            cumulene,
            5,
            allene,
            1,
            None,
            Some(&mut error),
            AB_PARITY_UNKN as i32,
        )
        .unwrap();
        let cumulene_values = heap.slice(cumulene.as_const()).unwrap();
        assert_eq!(cumulene_values[1].sb_ord[0], 1);
        assert_eq!(cumulene_values[3].sb_ord[0], 0);
        assert_eq!(cumulene_values[1].sn_orig_at_num[0], 30);
        assert_eq!(cumulene_values[3].sn_orig_at_num[0], 34);

        let allocation_atoms = heap.allocate(tetrahedral_atoms()).unwrap();
        let allocation_descriptor = heap
            .allocate(vec![inchi_Stereo0D {
                neighbor: [0, 1, 2, 3],
                central_atom: 1,
                type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
                parity: tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8,
            }])
            .unwrap();
        let allocation_error = heap.allocate(vec![0_i8; 256]).unwrap();
        heap.fail_after_allocations(0);
        Extract0DParities(
            &mut heap,
            allocation_atoms,
            4,
            allocation_descriptor,
            1,
            Some(allocation_error),
            Some(&mut error),
            AB_PARITY_UNKN as i32,
        )
        .unwrap();
        assert_eq!(
            extract_error_text(&heap, allocation_error),
            "0D Parities Reconciliation failed: -1"
        );
    }

    #[test]
    fn source_port__readinch__createinchi_stereo0d__line_118() {
        let mut heap = SourceHeap::default();

        let stereo = CreateInchi_Stereo0D(&mut heap, 2).unwrap();
        assert_eq!(
            heap.slice(stereo.as_const()).unwrap(),
            &[inchi_Stereo0D::default(), inchi_Stereo0D::default()]
        );
        assert_eq!(inchi_free(&mut heap, stereo), Ok(()));

        assert_eq!(
            CreateInchi_Stereo0D(&mut heap, -1),
            Err(SourceHeapError::AllocationSizeOverflow)
        );
    }

    #[test]
    fn source_port__readinch__freeinchi_stereo0d__line_125() {
        let mut heap = SourceHeap::default();
        assert_eq!(FreeInchi_Stereo0D(&mut heap, None), Ok(()));

        let mut null_slot = SourceMutPointer::null();
        assert_eq!(FreeInchi_Stereo0D(&mut heap, Some(&mut null_slot)), Ok(()));
        assert!(null_slot.is_null());

        let mut stereo_slot = allocate_source_fixture(
            &mut heap,
            vec![inchi_Stereo0D::default(), inchi_Stereo0D::default()],
        );
        let alias = stereo_slot.as_const();
        assert_eq!(
            FreeInchi_Stereo0D(&mut heap, Some(&mut stereo_slot)),
            Ok(())
        );
        assert!(stereo_slot.is_null());
        assert_eq!(heap.slice(alias), Err(SourceHeapError::MissingAllocation));

        let allocation = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default(); 2]);
        let mut interior_slot = allocation.offset(1).unwrap();
        assert_eq!(
            FreeInchi_Stereo0D(&mut heap, Some(&mut interior_slot)),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert_eq!(interior_slot, allocation.offset(1).unwrap());
        assert_eq!(heap.slice(allocation.as_const()).unwrap().len(), 2);
        assert_eq!(inchi_free(&mut heap, allocation), Ok(()));
    }
}
