#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "mode.h"
#include "inchi_api.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "ichierr.h"
#include "ichicant.h"
#include "ichirvrs.h"
#include "ichitaut.h"
#include "inpdef.h"
#include "ichi.h"
#include "ichi_io.h"
#include "ichitime.h"
#include "inchi_dll_b.h"
#include "readinch.h"
#include "sha2.h"
#include "util.h"

void *__real_malloc(size_t size);
void *__real_calloc(size_t count, size_t size);
void *__real_realloc(void *pointer, size_t size);
void __real_free(void *pointer);

static int ORACLE_ALLOCATION_FAILURE_ENABLED = 0;
static int ORACLE_ALLOCATION_ORDINAL = 0;
static int ORACLE_ALLOCATION_CALLS = 0;
static int ORACLE_DEFER_FREES = 0;
static void *ORACLE_DEFERRED_FREES[4096];
static size_t ORACLE_DEFERRED_FREE_COUNT = 0;

static int ORACLE_NORMALIZE_ACTIVE = 0;
static int ORACLE_NORMALIZE_FORCED_REBUILD_RETURN = 0;
static int ORACLE_NORMALIZE_FORCE_REBUILD = 0;
static int ORACLE_NORMALIZE_FORCE_COMPARE = 0;
static int ORACLE_NORMALIZE_PREFREE_EXACT = 1;
static char ORACLE_NORMALIZE_EVENTS[8192];
static size_t ORACLE_NORMALIZE_EVENTS_LENGTH = 0;
static INChI **ORACLE_NORMALIZE_INCHI_HOLDERS[TAUT_NUM];
static INChI_Aux **ORACLE_NORMALIZE_AUX_HOLDERS[TAUT_NUM];
static INP_ATOM_DATA *ORACLE_NORMALIZE_NORM_HOLDERS[TAUT_NUM];
static INChI *ORACLE_NORMALIZE_REVERSED[TAUT_NUM];
static INChI *ORACLE_NORMALIZE_ORIGINAL[TAUT_NUM];
static int ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[TAUT_NUM];
static int ORACLE_NORMALIZE_ZZ_ACTIVE = 0;
static int ORACLE_NORMALIZE_ZZ_FAILURE_STAGE = 0;
static int ORACLE_NORMALIZE_ZZ_CALLOC_ORDINAL = 0;
static int ORACLE_NORMALIZE_ZZ_FORCE_GROWTH = 0;
static int ORACLE_NORMALIZE_ZZ_SUCCESSFUL_ALLOCATIONS = 0;
static int ORACLE_NORMALIZE_ZZ_FREES = 0;
static int ORACLE_NORMALIZE_ZZ_REALLOC_CALLS = 0;
static char ORACLE_NORMALIZE_ZZ_FORMULA[256];
static int ORACLE_NORMALIZE_FIRST_LESS_ACTIVE = 0;
static int ORACLE_NORMALIZE_FINAL_ACTIVE = 0;
static int ORACLE_NORMALIZE_REBUILD_RETURNS[8];
static int ORACLE_NORMALIZE_REBUILD_PREP[8];
static int ORACLE_NORMALIZE_REBUILD_NORM[8];
static int ORACLE_NORMALIZE_REBUILD_TINFO[8];
static int ORACLE_NORMALIZE_REBUILD_LENGTH = 0;
static int ORACLE_NORMALIZE_REBUILD_POSITION = 0;
static INCHI_MODE ORACLE_NORMALIZE_COMPARE_FLAGS[8];
static int ORACLE_NORMALIZE_COMPARE_ERRORS[8];
static int ORACLE_NORMALIZE_COMPARE_H1[8];
static int ORACLE_NORMALIZE_COMPARE_H2[8];
static int ORACLE_NORMALIZE_COMPARE_ENDPOINTS[8];
static int ORACLE_NORMALIZE_COMPARE_LENGTH = 0;
static int ORACLE_NORMALIZE_COMPARE_POSITION = 0;
static int ORACLE_NORMALIZE_FILL_RETURNS[8];
static int ORACLE_NORMALIZE_FILL_LENGTH = 0;
static int ORACLE_NORMALIZE_FILL_POSITION = 0;
static int ORACLE_NORMALIZE_FIX_LESS_RETURNS[8];
static int ORACLE_NORMALIZE_FIX_LESS_LENGTH = 0;
static int ORACLE_NORMALIZE_FIX_LESS_POSITION = 0;
static int ORACLE_NORMALIZE_FIX_MORE_RETURNS[8];
static int ORACLE_NORMALIZE_FIX_MORE_LENGTH = 0;
static int ORACLE_NORMALIZE_FIX_MORE_POSITION = 0;
static int ORACLE_NORMALIZE_FIX_EXTRA_RETURNS[8];
static int ORACLE_NORMALIZE_FIX_EXTRA_LENGTH = 0;
static int ORACLE_NORMALIZE_FIX_EXTRA_POSITION = 0;
static int ORACLE_NORMALIZE_FIX_FIXED_RETURNS[8];
static int ORACLE_NORMALIZE_FIX_FIXED_LENGTH = 0;
static int ORACLE_NORMALIZE_FIX_FIXED_POSITION = 0;
static int ORACLE_NORMALIZE_FIX_MOBILE_RETURNS[8];
static int ORACLE_NORMALIZE_FIX_MOBILE_LENGTH = 0;
static int ORACLE_NORMALIZE_FIX_MOBILE_POSITION = 0;
static int ORACLE_NORMALIZE_FIX_STEREO_RETURNS[8];
static int ORACLE_NORMALIZE_FIX_STEREO_LENGTH = 0;
static int ORACLE_NORMALIZE_FIX_STEREO_POSITION = 0;
static inp_ATOM *ORACLE_NORMALIZE_EXPECTED_NORM_AT[TAUT_NUM];
static int ORACLE_NORMALIZE_EXPECTED_TINFO_PRESENT = 0;
static int ORACLE_NORMALIZE_EXPECTED_TINFO_LENGTH = 0;
static AT_RANK ORACLE_NORMALIZE_EXPECTED_TINFO_H[2];
static AT_NUMB ORACLE_NORMALIZE_EXPECTED_TINFO_ENDPOINTS[2];

typedef struct tagOracleNormalizeAllocation
{
    void *pointer;
    const char *label;
} ORACLE_NORMALIZE_ALLOCATION;

static ORACLE_NORMALIZE_ALLOCATION ORACLE_NORMALIZE_ALLOCATIONS[16];
static size_t ORACLE_NORMALIZE_ALLOCATION_COUNT = 0;

static void oracle_normalize_event(const char *event)
{
    int written;
    size_t remaining;
    if (!ORACLE_NORMALIZE_ACTIVE)
    {
        return;
    }
    remaining = sizeof(ORACLE_NORMALIZE_EVENTS) -
                ORACLE_NORMALIZE_EVENTS_LENGTH;
    written = snprintf(ORACLE_NORMALIZE_EVENTS + ORACLE_NORMALIZE_EVENTS_LENGTH,
                       remaining, "%s\"%s\"",
                       ORACLE_NORMALIZE_EVENTS_LENGTH ? "," : "", event);
    if (written < 0 || (size_t) written >= remaining)
    {
        fputs("NormalizeAndCompare event buffer exceeded\n", stderr);
        abort();
    }
    ORACLE_NORMALIZE_EVENTS_LENGTH += (size_t) written;
}

static void oracle_normalize_register_allocation(void *pointer,
                                                 const char *label)
{
    if (!pointer || ORACLE_NORMALIZE_ALLOCATION_COUNT >=
                        sizeof(ORACLE_NORMALIZE_ALLOCATIONS) /
                            sizeof(ORACLE_NORMALIZE_ALLOCATIONS[0]))
    {
        fputs("NormalizeAndCompare allocation registry exceeded\n", stderr);
        abort();
    }
    ORACLE_NORMALIZE_ALLOCATIONS[ORACLE_NORMALIZE_ALLOCATION_COUNT].pointer =
        pointer;
    ORACLE_NORMALIZE_ALLOCATIONS[ORACLE_NORMALIZE_ALLOCATION_COUNT].label =
        label;
    ORACLE_NORMALIZE_ALLOCATION_COUNT++;
}

static void oracle_normalize_record_free(void *pointer)
{
    size_t index;
    char event[96];
    if (!ORACLE_NORMALIZE_ACTIVE)
    {
        return;
    }
    for (index = 0; index < ORACLE_NORMALIZE_ALLOCATION_COUNT; index++)
    {
        if (ORACLE_NORMALIZE_ALLOCATIONS[index].pointer == pointer)
        {
            snprintf(event, sizeof(event), "free:%s",
                     ORACLE_NORMALIZE_ALLOCATIONS[index].label);
            oracle_normalize_event(event);
            return;
        }
    }
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        return;
    }
    ORACLE_NORMALIZE_PREFREE_EXACT = 0;
    oracle_normalize_event("free:unregistered");
}

static int oracle_allocation_should_fail(void)
{
    if (!ORACLE_ALLOCATION_FAILURE_ENABLED)
    {
        return 0;
    }
    ORACLE_ALLOCATION_CALLS++;
    return ORACLE_ALLOCATION_ORDINAL > 0 &&
           ORACLE_ALLOCATION_CALLS == ORACLE_ALLOCATION_ORDINAL;
}

void *__wrap_malloc(size_t size)
{
    if (oracle_allocation_should_fail())
    {
        return NULL;
    }
    return __real_malloc(size);
}

void *__wrap_calloc(size_t count, size_t size)
{
    void *pointer;
    int should_fail = oracle_allocation_should_fail();
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        ORACLE_NORMALIZE_ZZ_CALLOC_ORDINAL++;
        should_fail = should_fail ||
                      ORACLE_NORMALIZE_ZZ_FAILURE_STAGE ==
                          ORACLE_NORMALIZE_ZZ_CALLOC_ORDINAL;
    }
    if (should_fail)
    {
        return NULL;
    }
    pointer = __real_calloc(count, size);
    if (pointer && ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        ORACLE_NORMALIZE_ZZ_SUCCESSFUL_ALLOCATIONS++;
    }
    return pointer;
}

void *__wrap_realloc(void *pointer, size_t size)
{
    void *replacement;
    size_t index;
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        ORACLE_NORMALIZE_ZZ_REALLOC_CALLS++;
        oracle_normalize_event("inchi_realloc:begin");
        if (ORACLE_NORMALIZE_ZZ_FAILURE_STAGE == 4)
        {
            oracle_normalize_event("inchi_realloc:null");
            return NULL;
        }
    }
    replacement = __real_realloc(pointer, size);
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        oracle_normalize_event(replacement ? "inchi_realloc:non-null"
                                           : "inchi_realloc:null");
        if (replacement)
        {
            for (index = 0; index < ORACLE_NORMALIZE_ALLOCATION_COUNT; index++)
            {
                if (ORACLE_NORMALIZE_ALLOCATIONS[index].pointer == pointer)
                {
                    ORACLE_NORMALIZE_ALLOCATIONS[index].pointer = replacement;
                    break;
                }
            }
        }
    }
    return replacement;
}

void __wrap_free(void *pointer)
{
    size_t i;
    if (!pointer)
    {
        return;
    }
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        ORACLE_NORMALIZE_ZZ_FREES++;
    }
    oracle_normalize_record_free(pointer);
    for (i = 0; i < ORACLE_DEFERRED_FREE_COUNT; i++)
    {
        if (ORACLE_DEFERRED_FREES[i] == pointer)
        {
            if (ORACLE_DEFER_FREES)
            {
                return;
            }
            ORACLE_DEFERRED_FREES[i] =
                ORACLE_DEFERRED_FREES[--ORACLE_DEFERRED_FREE_COUNT];
            __real_free(pointer);
            return;
        }
    }
    if (ORACLE_DEFER_FREES)
    {
        if (ORACLE_DEFERRED_FREE_COUNT >=
            sizeof(ORACLE_DEFERRED_FREES) / sizeof(ORACLE_DEFERRED_FREES[0]))
        {
            fputs("official oracle deferred-free capacity exceeded\n", stderr);
            abort();
        }
        ORACLE_DEFERRED_FREES[ORACLE_DEFERRED_FREE_COUNT++] = pointer;
        return;
    }
    __real_free(pointer);
}

int __real_MakeOneInChIOutOfStrFromINChI2(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, StrFromINChI *pStruct,
    inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA,
    ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **t_group_info,
    inp_ATOM **at_norm, inp_ATOM **at_prep);

int __wrap_MakeOneInChIOutOfStrFromINChI2(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, StrFromINChI *pStruct,
    inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA,
    ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **t_group_info,
    inp_ATOM **at_norm, inp_ATOM **at_prep)
{
    if ((ORACLE_NORMALIZE_FIRST_LESS_ACTIVE ||
         ORACLE_NORMALIZE_FINAL_ACTIVE) &&
        ORACLE_NORMALIZE_REBUILD_POSITION < ORACLE_NORMALIZE_REBUILD_LENGTH)
    {
        int position = ORACLE_NORMALIZE_REBUILD_POSITION++;
        oracle_normalize_event("MakeOneInChIOutOfStrFromINChI2");
        *at_norm = ORACLE_NORMALIZE_REBUILD_NORM[position] ? at : NULL;
        *at_prep = ORACLE_NORMALIZE_REBUILD_PREP[position] ? at : NULL;
        *t_group_info = ORACLE_NORMALIZE_REBUILD_TINFO[position]
                            ? &pStruct->One_ti
                            : NULL;
        return ORACLE_NORMALIZE_REBUILD_RETURNS[position];
    }
    if (ORACLE_NORMALIZE_ACTIVE && ORACLE_NORMALIZE_FORCE_REBUILD)
    {
        oracle_normalize_event("MakeOneInChIOutOfStrFromINChI2");
        return ORACLE_NORMALIZE_FORCED_REBUILD_RETURN;
    }
    return __real_MakeOneInChIOutOfStrFromINChI2(
        pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
        t_group_info, at_norm, at_prep);
}

INCHI_MODE __real_CompareReversedINChI2(INChI *i1, INChI *i2,
                                       INChI_Aux *a1, INChI_Aux *a2,
                                       ICR *picr, int *err);

INCHI_MODE __wrap_CompareReversedINChI2(INChI *i1, INChI *i2,
                                       INChI_Aux *a1, INChI_Aux *a2,
                                       ICR *picr, int *err)
{
    if (ORACLE_NORMALIZE_ACTIVE)
    {
        int reverse_index = i1 == ORACLE_NORMALIZE_REVERSED[1] ? 1 : 0;
        int original_index = i2 == ORACLE_NORMALIZE_ORIGINAL[1] ? 1 : 0;
        char event[72];
        snprintf(event, sizeof(event),
                 "CompareReversedINChI2:r%d:o%d", reverse_index,
                 original_index);
        oracle_normalize_event(event);
    }
    if ((ORACLE_NORMALIZE_FIRST_LESS_ACTIVE ||
         ORACLE_NORMALIZE_FINAL_ACTIVE) &&
        ORACLE_NORMALIZE_COMPARE_POSITION < ORACLE_NORMALIZE_COMPARE_LENGTH)
    {
        int position = ORACLE_NORMALIZE_COMPARE_POSITION++;
        memset(picr, 0, sizeof(*picr));
        picr->tot_num_H1 = ORACLE_NORMALIZE_COMPARE_H1[position];
        picr->tot_num_H2 = ORACLE_NORMALIZE_COMPARE_H2[position];
        picr->num_endp_in1_only =
            ORACLE_NORMALIZE_COMPARE_ENDPOINTS[position];
        *err = ORACLE_NORMALIZE_COMPARE_ERRORS[position];
        (void) a1;
        (void) a2;
        return ORACLE_NORMALIZE_COMPARE_FLAGS[position];
    }
    if (ORACLE_NORMALIZE_ACTIVE && ORACLE_NORMALIZE_FORCE_COMPARE)
    {
        memset(picr, 0, sizeof(*picr));
        *err = 0;
        (void) a1;
        (void) a2;
        return IDIF_PROBLEM;
    }
    return __real_CompareReversedINChI2(i1, i2, a1, a2, picr, err);
}

int __real_inchi_strbuf_init(INCHI_IOS_STRING *buf, int start_size,
                             int incr_size);
void __real_inchi_strbuf_close(INCHI_IOS_STRING *buf);
int __real_MergeZzInHillFormula(INCHI_IOS_STRING *strbuf);

int __wrap_inchi_strbuf_init(INCHI_IOS_STRING *buf, int start_size,
                             int incr_size)
{
    int result;
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        oracle_normalize_event("inchi_strbuf_init:begin");
    }
    result = __real_inchi_strbuf_init(buf, start_size, incr_size);
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        oracle_normalize_event(result > 0 ? "inchi_strbuf_init:positive"
                                          : "inchi_strbuf_init:non-positive");
    }
    return result;
}

int __wrap_MergeZzInHillFormula(INCHI_IOS_STRING *strbuf)
{
    int result;
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        oracle_normalize_event("MergeZzInHillFormula:begin");
    }
    result = __real_MergeZzInHillFormula(strbuf);
    if (ORACLE_NORMALIZE_ZZ_ACTIVE && ORACLE_NORMALIZE_ZZ_FORCE_GROWTH &&
        result == 0)
    {
        strbuf->nUsedLength = strbuf->nAllocatedLength + 1;
        oracle_normalize_event("MergeZzInHillFormula:forced-growth");
    }
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        oracle_normalize_event(result == 0 ? "MergeZzInHillFormula:zero"
                                           : "MergeZzInHillFormula:negative");
    }
    return result;
}

void __wrap_inchi_strbuf_close(INCHI_IOS_STRING *buf)
{
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        oracle_normalize_event("inchi_strbuf_close:begin");
    }
    __real_inchi_strbuf_close(buf);
    if (ORACLE_NORMALIZE_ZZ_ACTIVE)
    {
        oracle_normalize_event("inchi_strbuf_close:end");
    }
}

int oracle_test_FillOutExtraFixedHDataRestr(StrFromINChI *pStruct)
{
    if ((ORACLE_NORMALIZE_FIRST_LESS_ACTIVE ||
         ORACLE_NORMALIZE_FINAL_ACTIVE) &&
        ORACLE_NORMALIZE_FILL_POSITION < ORACLE_NORMALIZE_FILL_LENGTH)
    {
        oracle_normalize_event("FillOutExtraFixedHDataRestr");
        return ORACLE_NORMALIZE_FILL_RETURNS[ORACLE_NORMALIZE_FILL_POSITION++];
    }
    return FillOutExtraFixedHDataRestr(pStruct);
}

int oracle_test_FixLessHydrogenInFormula(
    BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
    inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA,
    ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta,
    int forbidden_edge_mask)
{
    if (ORACLE_NORMALIZE_FIRST_LESS_ACTIVE &&
        ORACLE_NORMALIZE_FIX_LESS_POSITION <
            ORACLE_NORMALIZE_FIX_LESS_LENGTH)
    {
        oracle_normalize_event("FixLessHydrogenInFormula");
        return ORACLE_NORMALIZE_FIX_LESS_RETURNS
            [ORACLE_NORMALIZE_FIX_LESS_POSITION++];
    }
    return FixLessHydrogenInFormula(
        pBNS, pBD, pStruct, at, at2, atf, pVA, pTCGroups, pnNumRunBNS,
        pnTotalDelta, forbidden_edge_mask);
}

int oracle_test_FixMoreHydrogenInFormula(
    BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
    inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA,
    ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta,
    int forbidden_edge_mask)
{
    if (ORACLE_NORMALIZE_FIRST_LESS_ACTIVE &&
        ORACLE_NORMALIZE_FIX_MORE_POSITION <
            ORACLE_NORMALIZE_FIX_MORE_LENGTH)
    {
        oracle_normalize_event("FixMoreHydrogenInFormula");
        return ORACLE_NORMALIZE_FIX_MORE_RETURNS
            [ORACLE_NORMALIZE_FIX_MORE_POSITION++];
    }
    return FixMoreHydrogenInFormula(
        pBNS, pBD, pStruct, at, at2, atf, pVA, pTCGroups, pnNumRunBNS,
        pnTotalDelta, forbidden_edge_mask);
}

int oracle_test_FixRemoveExtraTautEndpoints(
    BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
    inp_ATOM *at2, inp_ATOM *atf, inp_ATOM *atn, VAL_AT *pVA,
    ALL_TC_GROUPS *pTCGroups, ICR *picr, int *pnNumRunBNS,
    int *pnTotalDelta, int forbidden_edge_mask)
{
    if (ORACLE_NORMALIZE_FIRST_LESS_ACTIVE &&
        ORACLE_NORMALIZE_FIX_EXTRA_POSITION <
            ORACLE_NORMALIZE_FIX_EXTRA_LENGTH)
    {
        oracle_normalize_event("FixRemoveExtraTautEndpoints");
        return ORACLE_NORMALIZE_FIX_EXTRA_RETURNS
            [ORACLE_NORMALIZE_FIX_EXTRA_POSITION++];
    }
    return FixRemoveExtraTautEndpoints(
        pBNS, pBD, pStruct, at, at2, atf, atn, pVA, pTCGroups, picr,
        pnNumRunBNS, pnTotalDelta, forbidden_edge_mask);
}

int __real_FixFixedHRestoredStructure(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
    StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ti,
    inp_ATOM **at_norm, inp_ATOM **at_prep, INChI *pInChI[], long num_inp,
    int bHasSomeFixedH, int *pnNumRunBNS, int *pnTotalDelta,
    int forbidden_edge_mask, int forbidden_stereo_edge_mask);

int __wrap_FixFixedHRestoredStructure(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
    StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ti,
    inp_ATOM **at_norm, inp_ATOM **at_prep, INChI *pInChI[], long num_inp,
    int bHasSomeFixedH, int *pnNumRunBNS, int *pnTotalDelta,
    int forbidden_edge_mask, int forbidden_stereo_edge_mask)
{
    if (ORACLE_NORMALIZE_FINAL_ACTIVE &&
        ORACLE_NORMALIZE_FIX_FIXED_POSITION <
            ORACLE_NORMALIZE_FIX_FIXED_LENGTH)
    {
        char event[96];
        snprintf(event, sizeof(event),
                 "FixFixedHRestoredStructure:num_inp=%ld", num_inp);
        oracle_normalize_event(event);
        return ORACLE_NORMALIZE_FIX_FIXED_RETURNS
            [ORACLE_NORMALIZE_FIX_FIXED_POSITION++];
    }
    return __real_FixFixedHRestoredStructure(
        pCG, ic, ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA,
        pTCGroups, ti, at_norm, at_prep, pInChI, num_inp, bHasSomeFixedH,
        pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
        forbidden_stereo_edge_mask);
}

int __real_FixMobileHRestoredStructure(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
    StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ti,
    inp_ATOM **at_norm, inp_ATOM **at_prep, INChI *pInChI[], long num_inp,
    int bHasSomeFixedH, int *pnNumRunBNS, int *pnTotalDelta,
    int forbidden_edge_mask, int forbidden_stereo_edge_mask);

int __wrap_FixMobileHRestoredStructure(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
    StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ti,
    inp_ATOM **at_norm, inp_ATOM **at_prep, INChI *pInChI[], long num_inp,
    int bHasSomeFixedH, int *pnNumRunBNS, int *pnTotalDelta,
    int forbidden_edge_mask, int forbidden_stereo_edge_mask)
{
    if (ORACLE_NORMALIZE_FINAL_ACTIVE &&
        ORACLE_NORMALIZE_FIX_MOBILE_POSITION <
            ORACLE_NORMALIZE_FIX_MOBILE_LENGTH)
    {
        char event[96];
        snprintf(event, sizeof(event),
                 "FixMobileHRestoredStructure:num_inp=%ld", num_inp);
        oracle_normalize_event(event);
        return ORACLE_NORMALIZE_FIX_MOBILE_RETURNS
            [ORACLE_NORMALIZE_FIX_MOBILE_POSITION++];
    }
    return __real_FixMobileHRestoredStructure(
        pCG, ic, ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA,
        pTCGroups, ti, at_norm, at_prep, pInChI, num_inp, bHasSomeFixedH,
        pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
        forbidden_stereo_edge_mask);
}

int __real_FixRestoredStructureStereo(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INCHI_MODE cmpInChI, ICR *icr,
    INCHI_MODE cmpInChI2, ICR *icr2, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
    StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ti,
    inp_ATOM **at_norm, inp_ATOM **at_prep, INChI *pInChI[], long num_inp,
    int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask,
    int forbidden_stereo_edge_mask);

int __wrap_FixRestoredStructureStereo(
    CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INCHI_MODE cmpInChI, ICR *icr,
    INCHI_MODE cmpInChI2, ICR *icr2, const INPUT_PARMS *ip,
    STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
    StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ti,
    inp_ATOM **at_norm, inp_ATOM **at_prep, INChI *pInChI[], long num_inp,
    int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask,
    int forbidden_stereo_edge_mask)
{
    if (ORACLE_NORMALIZE_FINAL_ACTIVE &&
        ORACLE_NORMALIZE_FIX_STEREO_POSITION <
            ORACLE_NORMALIZE_FIX_STEREO_LENGTH)
    {
        oracle_normalize_event("FixRestoredStructureStereo");
        return ORACLE_NORMALIZE_FIX_STEREO_RETURNS
            [ORACLE_NORMALIZE_FIX_STEREO_POSITION++];
    }
    return __real_FixRestoredStructureStereo(
        pCG, ic, cmpInChI, icr, cmpInChI2, icr2, ip, sd, pBNS, pBD,
        pStruct, at, at2, at3, pVA, pTCGroups, ti, at_norm, at_prep,
        pInChI, num_inp, pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
        forbidden_stereo_edge_mask);
}

int __real_Free_INChI(INChI **inchi);
int __real_Free_INChI_Aux(INChI_Aux **aux);
void __real_FreeInpAtomData(INP_ATOM_DATA *data);
int __real_free_t_group_info(T_GROUP_INFO *info);

int __wrap_Free_INChI(INChI **inchi)
{
    int index;
    if (ORACLE_NORMALIZE_ACTIVE)
    {
        char event[64];
        for (index = 0; index < TAUT_NUM; index++)
        {
            if (ORACLE_NORMALIZE_INCHI_HOLDERS[index] == inchi)
            {
                snprintf(event, sizeof(event), "Free_INChI[%d]:%s", index,
                         *inchi ? "present" : "null");
                oracle_normalize_event(event);
                if (*inchi && ((*inchi)->nErrorCode != 100 + index ||
                               (*inchi)->nNumberOfAtoms !=
                                   ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[index] ||
                               (*inchi)->nRefCount != 0))
                {
                    ORACLE_NORMALIZE_PREFREE_EXACT = 0;
                }
                break;
            }
        }
        if (*inchi && (*inchi)->szHillFormula)
        {
            snprintf(ORACLE_NORMALIZE_ZZ_FORMULA,
                     sizeof(ORACLE_NORMALIZE_ZZ_FORMULA), "%s",
                     (*inchi)->szHillFormula);
        }
    }
    return __real_Free_INChI(inchi);
}

int __wrap_Free_INChI_Aux(INChI_Aux **aux)
{
    int index;
    if (ORACLE_NORMALIZE_ACTIVE)
    {
        char event[72];
        for (index = 0; index < TAUT_NUM; index++)
        {
            if (ORACLE_NORMALIZE_AUX_HOLDERS[index] == aux)
            {
                snprintf(event, sizeof(event), "Free_INChI_Aux[%d]:%s", index,
                         *aux ? "present" : "null");
                oracle_normalize_event(event);
                if (*aux && ((*aux)->nErrorCode != 200 + index ||
                             (*aux)->nNumberOfAtoms != 20 + index ||
                             (*aux)->nRefCount != 0))
                {
                    ORACLE_NORMALIZE_PREFREE_EXACT = 0;
                }
                break;
            }
        }
    }
    return __real_Free_INChI_Aux(aux);
}

void __wrap_FreeInpAtomData(INP_ATOM_DATA *data)
{
    int index;
    if (ORACLE_NORMALIZE_ACTIVE)
    {
        char event[72];
        for (index = 0; index < TAUT_NUM; index++)
        {
            if (ORACLE_NORMALIZE_NORM_HOLDERS[index] == data ||
                (!data && !ORACLE_NORMALIZE_NORM_HOLDERS[index]))
            {
                snprintf(event, sizeof(event), "FreeInpAtomData[%d]:%s", index,
                         data ? "present" : "null");
                oracle_normalize_event(event);
                if (data &&
                    (data->num_at != 30 + index ||
                     data->num_bonds != 40 + index ||
                     data->at != ORACLE_NORMALIZE_EXPECTED_NORM_AT[index] ||
                     data->at_fixed_bonds))
                {
                    ORACLE_NORMALIZE_PREFREE_EXACT = 0;
                }
                ORACLE_NORMALIZE_NORM_HOLDERS[index] =
                    (INP_ATOM_DATA *) (uintptr_t) 1;
                break;
            }
        }
    }
    __real_FreeInpAtomData(data);
}

int __wrap_free_t_group_info(T_GROUP_INFO *info)
{
    if (ORACLE_NORMALIZE_ACTIVE)
    {
        int index;
        oracle_normalize_event(info && info->t_group
                                   ? "free_t_group_info:present"
                                   : "free_t_group_info:null");
        if (ORACLE_NORMALIZE_FINAL_ACTIVE)
        {
            if (!info || (!!info->t_group !=
                          ORACLE_NORMALIZE_EXPECTED_TINFO_PRESENT) ||
                info->num_t_groups != ORACLE_NORMALIZE_EXPECTED_TINFO_LENGTH)
            {
                ORACLE_NORMALIZE_PREFREE_EXACT = 0;
            }
            else
            {
                for (index = 0;
                     index < ORACLE_NORMALIZE_EXPECTED_TINFO_LENGTH;
                     index++)
                {
                    if (info->t_group[index].num[0] !=
                            ORACLE_NORMALIZE_EXPECTED_TINFO_H[index] ||
                        info->t_group[index].nNumEndpoints !=
                            ORACLE_NORMALIZE_EXPECTED_TINFO_ENDPOINTS[index])
                    {
                        ORACLE_NORMALIZE_PREFREE_EXACT = 0;
                    }
                }
            }
        }
        else if (!info || !info->t_group || info->num_t_groups != 1 ||
                 info->t_group[0].nNumEndpoints != 51 ||
                 info->t_group[0].num[0] != 52 ||
                 info->t_group[0].num[1] != 53)
        {
            ORACLE_NORMALIZE_PREFREE_EXACT = 0;
        }
    }
    return __real_free_t_group_info(info);
}

static void oracle_flush_deferred_frees(void)
{
    while (ORACLE_DEFERRED_FREE_COUNT)
    {
        __real_free(ORACLE_DEFERRED_FREES[--ORACLE_DEFERRED_FREE_COUNT]);
    }
}

int SetAtomProperties(inp_ATOM *at,
                      MOL_COORD *szCoord,
                      inchi_Atom *ati,
                      int a1,
                      int *nDim,
                      char *pStrErr,
                      int *err);
int SetBondProperties(inp_ATOM *at,
                      inchi_Atom *ati,
                      int a1,
                      int j,
                      int nNumAtoms,
                      int *nNumBonds,
                      char *pStrErr,
                      int *err);
int InchiToInchi_Input(INCHI_IOSTREAM *inp_molfile,
                       inchi_Input *orig_at_data,
                       int bMergeAllInputStructures,
                       int bDoNotAddH,
                       int vABParityUnknown,
                       INPUT_TYPE nInputType,
                       char *pSdfLabel,
                       char *pSdfValue,
                       long *lSdfId,
                       INCHI_MODE *pInpAtomFlags,
                       int *err,
                       char *pStrErr);
int cmp_components(const void *a1, const void *a2);
int cmp_rad_endpoints(const void *a1, const void *a2);
int cmp_iso_atw_diff_component_no(const void *a1, const void *a2);
int CompTGroupNumber(const void *tg1, const void *tg2, void *p);
int CompCGroupNumber(const void *cg1, const void *cg2, void *p);
int CompNeighLists(const void *a1, const void *a2, void *p);
int CompNeighListsUpToMaxRank(const void *a1, const void *a2, void *p);
int cmp_charge_val(const void *a1, const void *a2, void *p);
int comp_cc_cand(const void *a1, const void *a2);
const char *base26_triplet_1(const unsigned char *a);
const char *base26_triplet_2(const unsigned char *a);
const char *base26_triplet_3(const unsigned char *a);
const char *base26_triplet_4(const unsigned char *a);
const char *base26_dublet_for_bits_28_to_36(unsigned char *a);
const char *base26_dublet_for_bits_56_to_64(unsigned char *a);
void get_xtra_hash_major_hex(const unsigned char *a, char *szXtra);
void get_xtra_hash_minor_hex(const unsigned char *a, char *szXtra);
void oracle_sha2_starts(sha2_context *ctx);
void oracle_sha2_process(sha2_context *ctx, unsigned char data[64]);
void oracle_sha2_update(sha2_context *ctx, unsigned char *input, int ilen);
void oracle_sha2_finish(sha2_context *ctx, unsigned char output[32]);
void oracle_sha2_csum(unsigned char *input, int ilen,
                      unsigned char output[32]);

static const char *EXPECTED_VERSION = "1.07.5";
static const char *EXPECTED_DESCRIPTION =
    "InChI version 1, Software 1.07.5 (API Library)";

typedef struct tagOracleCase
{
    const char *case_id;
    const char *text;
    INPUT_TYPE input_type;
    int file_mode;
    int generated_text;
    int do_not_add_h;
    int caller_atom_capacity;
    int caller_stereo_capacity;
    int allocation_failure_ordinal;
    int omit_atom_output;
    int omit_stereo_output;
    int omit_label_output;
    int omit_value_output;
    int omit_id_output;
    int omit_flags_output;
    int pass_zero_atom_capacity;
    int pass_zero_stereo_count;
    int initial_error_code;
} ORACLE_CASE;

static const ORACLE_CASE INCHI_TO_ATOM_CASES[] = {
    {
        "plain-labeled-string",
        "Structure: 42. LABEL =value\n"
        "AuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:0,0,0;1,1,1;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-labeled-file",
        "Structure: 42. LABEL =value\n"
        "AuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:0,0,0;1,1,1;\n",
        INPUT_INCHI_PLAIN,
        1,
    },
    {
        "plain-integer-boundaries",
        "AuxInfo=1/rA:1C4+255.128i32768oh255d2t3/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-metal-stereo",
        "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d-+3;s4;/rC:;;;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-aromatic-warning",
        "AuxInfo=1/rA:2CO/rB:a1;/rC:0,0,0;1,0,0;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "xml-labeled",
        "<structure number=\"7\" id.name=\"xml-label\" id.value=\"xml-value\">\n"
        "<identifier.auxiliary-info>\n"
        "<reversibility>\n"
        "<atoms>\n2nCO\n</atoms>\n"
        "<bonds>\ns1;\n</bonds>\n"
        "<xyz>\n0,0,0;1,1,1;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-error",
        "<structure number=\"8\">\n"
        "<message type=\"error (no InChI)\" value=\"source error\">\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "plain-malformed-coordinates",
        "AuxInfo=1/rA:2CO/rB:s1;/rC:0,0x;1,0,0;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    { "plain-eof", "", INPUT_INCHI_PLAIN, 0 },
    {
        "xml-fatal",
        "<structure number=\"9\">\n"
        "<message type=\"fatal (aborted)\" value=\"source fatal\">\n",
        INPUT_INCHI_XML,
        0,
    },
    { "plain-missing-atom", "AuxInfo=1/no-reversibility\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-wrong-bonds", "AuxInfo=1/rA:2CO/rB:x/rC:;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-undefined-coordinates", "AuxInfo=1/rA:2CO/rB:s1;/rC:;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-2d-coordinates", "AuxInfo=1/rA:2CO/rB:s1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-v", "AuxInfo=1/rA:2CO/rB:v1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-V", "AuxInfo=1/rA:2CO/rB:V1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-w", "AuxInfo=1/rA:2CO/rB:w1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-t", "AuxInfo=1/rA:2CO/rB:t1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-p", "AuxInfo=1/rA:2CO/rB:p1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-P", "AuxInfo=1/rA:2CO/rB:P1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-n", "AuxInfo=1/rA:2CO/rB:n1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bond-N", "AuxInfo=1/rA:2CO/rB:N1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-tetrahedral-odd", "AuxInfo=1/rA:4CCCClo/rB:;;s1s2s3;/rC:;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-tetrahedral-even", "AuxInfo=1/rA:4CCCCle/rB:;;s1s2s3;/rC:;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-tetrahedral-unknown", "AuxInfo=1/rA:4CCCClu/rB:;;s1s2s3;/rC:;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-tetrahedral-undefined", "AuxInfo=1/rA:4CCCCl?/rB:;;s1s2s3;/rC:;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-odd-odd", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d--3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-odd-even", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d-+3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-odd-unknown", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d-u3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-odd-undefined", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d-?3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-even-odd", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d+-3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-even-even", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d++3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-even-unknown", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d+u3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-even-undefined", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d+?3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-unknown-odd", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;du-3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-unknown-even", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;du+3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-unknown-unknown", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;duu3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-unknown-undefined", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;du?3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-undefined-odd", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d?-3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-undefined-even", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d?+3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-undefined-unknown", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d?u3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-parity-undefined-undefined", "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d??3;s4;/rC:;;;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-do-not-add-h", "AuxInfo=1/rA:1C/rB:/rC:;\n", INPUT_INCHI_PLAIN, 0, 0, 1 },
    { "plain-caller-capacity", "AuxInfo=1/rA:2CCle/rB:s1;/rC:0,0,0;1,0,0;\n", INPUT_INCHI_PLAIN, 0, 0, 0, 3, 3 },
    { "plain-allocation-failure-atom", "AuxInfo=1/rA:1C/rB:/rC:;\n", INPUT_INCHI_PLAIN, 0, 0, 0, 0, 0, 1 },
    { "plain-allocation-failure-stereo", "AuxInfo=1/rA:1Cle/rB:/rC:;\n", INPUT_INCHI_PLAIN, 0, 0, 0, 0, 0, 2 },
    { "plain-cross-buffer", NULL, INPUT_INCHI_PLAIN, 0, 1 },
    { "plain-token-final-segment", NULL, INPUT_INCHI_PLAIN, 0, 2 },
    { "plain-empty-reversibility", "AuxInfo=1//\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-nonempty-double-slash", "AuxInfo=1//x\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-missing-bonds", "AuxInfo=1/rA:1C/rC:;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-missing-coordinates", "AuxInfo=1/rA:1C/rB:/\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-bonds-no-section-terminator", "AuxInfo=1/rA:1C/rB:\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-single-atom-trailing-bond-data", "AuxInfo=1/rA:1C/rB:x/rC:;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-wrong-atom-count", "AuxInfo=1/rA:2C/rB:;/rC:;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-wrong-bond-type", "AuxInfo=1/rA:2CO/rB:q1;/rC:;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-nonexistent-bond", "AuxInfo=1/rA:2CO/rB:s3;/rC:;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-wrong-bond-count", "AuxInfo=1/rA:3CCO/rB:s1/rC:;;;\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-wrong-coordinate-count", "AuxInfo=1/rA:2CO/rB:s1;/rC:;;;/tail\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-too-few-coordinates-section-terminator", "AuxInfo=1/rA:2CO/rB:s1;/rC:;/tail\n", INPUT_INCHI_PLAIN, 0 },
    { "plain-trailing-coordinate-data", "AuxInfo=1/rA:2CO/rB:s1;/rC:;;;\n", INPUT_INCHI_PLAIN, 0 },
    {
        "plain-long-integer-narrowing",
        "Structure: 9223372036854775807. LIMIT =value\n"
        "AuxInfo=1/rA:1C4+9223372036854775808.9223372036854775808i9223372036854775808oh9223372036854775808d9223372036854775808t9223372036854775808/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "xml-missing-reversibility",
        "<structure number=\"10\">\n"
        "<identifier.auxiliary-info>\n"
        "</identifier.auxiliary-info>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-reversibility-without-auxinfo",
        "<structure number=\"105\">\n<reversibility>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "plain-find-next-only",
        .text = "Structure: 12. FIND =next\n"
                "AuxInfo=1/rA:2CO/rB:s1;/rC:0,0,0;1,0,0;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_atom_output = 1,
    },
    {
        .case_id = "plain-no-stereo-output",
        .text = "AuxInfo=1/rA:4CCCCle/rB:;;s1s2s3;/rC:;;;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_stereo_output = 1,
    },
    {
        .case_id = "plain-error-no-stereo-output",
        .text = "AuxInfo=1/rA:2CO/rB:x/rC:;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_stereo_output = 1,
    },
    {
        .case_id = "plain-no-metadata-outputs",
        .text = "Structure: 13. OMIT =metadata\n"
                "AuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_label_output = 1,
        .omit_value_output = 1,
        .omit_id_output = 1,
        .omit_flags_output = 1,
    },
    {
        "plain-caller-atom-too-small",
        "AuxInfo=1/rA:2CO/rB:s1;/rC:0,0,0;1,0,0;\n",
        INPUT_INCHI_PLAIN,
        0, 0, 0, 1,
    },
    {
        "plain-caller-stereo-too-small",
        "AuxInfo=1/rA:5FeOCCCle/rB:;s1s2;d-+3;s4;/rC:;;;;;\n",
        INPUT_INCHI_PLAIN,
        0, 0, 0, 0, 1,
    },
    {
        "plain-header-no-dot-space",
        "Structure: 14 LABEL =value\n"
        "AuxInfo=1/rA:1C/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-header-is-missing",
        "Structure: 15. LABEL is missing\n"
        "AuxInfo=1/rA:1C/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-nonchiral-isolated",
        "AuxInfo=1/rA:1nC0/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-empty-atom-quantities",
        "AuxInfo=1/rA:1C+.iohdt/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-aromatic-two",
        "AuxInfo=1/rA:3CCC/rB:a1;a1a2;/rC:;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-aromatic-three",
        "AuxInfo=1/rA:4CCCC/rB:a1;a1a2;a1a2a3;/rC:;;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-aromatic-four",
        "AuxInfo=1/rA:5CCCCC/rB:a1;a1a2;a1a2a3;a1a2a3a4;/rC:;;;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        .case_id = "plain-allene",
        .text = "AuxInfo=1/rA:5CCCCC/rB:s1;d2;d-+3;s4;/rC:;;;;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .do_not_add_h = 1,
    },
    {
        .case_id = "plain-cumulene",
        .text = "AuxInfo=1/rA:6CCCCCC/rB:s1;d2;d3;d-+4;s5;/rC:;;;;;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .do_not_add_h = 1,
    },
    {
        "xml-eof-after-reversibility",
        "<structure number=\"16\">\n"
        "<identifier.auxiliary-info>\n"
        "<reversibility>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-wrong-atoms-tag",
        "<structure number=\"17\">\n"
        "<identifier.auxiliary-info>\n"
        "<reversibility>\n"
        "<not-atoms>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-eof-after-atoms-tag",
        "<structure number=\"18\">\n"
        "<identifier.auxiliary-info>\n"
        "<reversibility>\n"
        "<atoms>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-no-structure-number",
        "<structure id.name=\"only-label\">\n"
        "<identifier.auxiliary-info>\n"
        "<reversibility>\n"
        "<atoms>\n1cC\n</atoms>\n"
        "<bonds>\n\n</bonds>\n"
        "<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-caller-capacities-too-small",
        "<structure number=\"19\">\n"
        "<identifier.auxiliary-info>\n"
        "<reversibility>\n"
        "<atoms>\n2nCO\n</atoms>\n"
        "<bonds>\ns1;\n</bonds>\n"
        "<xyz>\n0,0,0;1,0,0;\n</xyz>\n",
        INPUT_INCHI_XML,
        0, 0, 0, 1, 1,
    },
    {
        "plain-zero-atom-count",
        "AuxInfo=1/rA:0/rB:/rC:\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-bond-nonalphabetic",
        "AuxInfo=1/rA:2CO/rB:1;/rC:;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        .case_id = "plain-stereo-reallocation",
        .text = "AuxInfo=1/rA:3CloCloClo/rB:d--1;d--1d--2;/rC:;;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .do_not_add_h = 1,
    },
    {
        .case_id = "plain-allocation-failure-stereo-reallocation",
        .text = "AuxInfo=1/rA:3CloCloClo/rB:d--1;d--1d--2;/rC:;;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .do_not_add_h = 1,
        .allocation_failure_ordinal = 3,
    },
    {
        .case_id = "plain-allocation-failure-coordinate",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .allocation_failure_ordinal = 3,
    },
    {
        "plain-partial-coordinate-components",
        "AuxInfo=1/rA:1C/rB:/rC:1,;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-empty-coordinate-components",
        "AuxInfo=1/rA:1C/rB:/rC:,1,;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "xml-rich-atom",
        "<structure number=\"20\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
        "1nC0+2.3i4oh2d3t4\n</atoms>\n"
        "<bonds>\n\n</bonds>\n<xyz>\n1,;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-all-bond-codes",
        "<structure number=\"21\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
        "12CCCCCCCCCCCC\n</atoms>\n<bonds>\n"
        "v1;V1;w1;s1;d1;t1;a1;p1;P1;n1;N1;\n</bonds>\n<xyz>\n"
        ";;;;;;;;;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-tetrahedral",
        "<structure number=\"22\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
        "4CCCClo\n</atoms>\n<bonds>\n;s1s2s3;\n</bonds>\n"
        "<xyz>\n;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-allene",
        .text = "<structure number=\"23\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "5CCCCC\n</atoms>\n<bonds>\n"
                "s1;d2;d-+3;s4;\n</bonds>\n<xyz>\n;;;;;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .do_not_add_h = 1,
    },
    {
        "xml-aromatic-two",
        "<structure number=\"24\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
        "3CCC\n</atoms>\n<bonds>\na1;a1a2;\n</bonds>\n"
        "<xyz>\n;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-aromatic-three",
        "<structure number=\"25\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
        "4CCCC\n</atoms>\n<bonds>\na1;a1a2;a1a2a3;\n</bonds>\n"
        "<xyz>\n;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-aromatic-four",
        "<structure number=\"26\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
        "5CCCCC\n</atoms>\n<bonds>\na1;a1a2;a1a2a3;a1a2a3a4;\n</bonds>\n"
        "<xyz>\n;;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-caller-capacities-reused",
        .text = "<structure number=\"27\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "4CCCClo\n</atoms>\n<bonds>\n;;s1s2s3;\n</bonds>\n"
                "<xyz>\n;;;;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .caller_atom_capacity = 5,
        .caller_stereo_capacity = 5,
    },
    {
        .case_id = "xml-caller-capacities-reused-error",
        .text = "<structure number=\"107\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "4CCCClo\n</atoms>\n<bonds>\n;s1s2s3;\n</bonds>\n"
                "<xyz>\n;;;;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .caller_atom_capacity = 5,
        .caller_stereo_capacity = 5,
    },
    {
        .case_id = "xml-no-stereo-output",
        .text = "<structure number=\"28\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "4CCCClo\n</atoms>\n<bonds>\n;;s1s2s3;\n</bonds>\n"
                "<xyz>\n;;;;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .omit_stereo_output = 1,
    },
    {
        "xml-metal-stereo",
        "<structure number=\"106\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n5FeOCCC\n</atoms>\n<bonds>\n;s1s2;d-+3;s4;\n</bonds>\n"
        "<xyz>\n;;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-find-next-only",
        .text = "<structure number=\"29\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "2CO\n</atoms>\n<bonds>\ns1;\n</bonds>\n"
                "<xyz>\n0,0,0;1,0,0;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .omit_atom_output = 1,
    },
    {
        .case_id = "xml-allocation-failure-atom",
        .text = "<structure number=\"30\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .allocation_failure_ordinal = 1,
    },
    {
        .case_id = "xml-allocation-failure-stereo",
        .text = "<structure number=\"31\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .allocation_failure_ordinal = 2,
    },
    {
        .case_id = "xml-allocation-failure-coordinate",
        .text = "<structure number=\"32\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n"
                "1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .allocation_failure_ordinal = 3,
    },
    {
        "xml-zero-atom-count",
        "<structure number=\"33\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n0\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-wrong-atom-count",
        "<structure number=\"34\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n2C\n</atoms>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-bond-nonalphabetic",
        "<structure number=\"35\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n2CO\n</atoms>\n"
        "<bonds>\n1;\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-missing-bond-neighbor",
        "<structure number=\"36\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n2CO\n</atoms>\n"
        "<bonds>\ns;\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-nonexistent-bond",
        "<structure number=\"37\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n2CO\n</atoms>\n"
        "<bonds>\ns3;\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-wrong-bond-count",
        "<structure number=\"38\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n3CCO\n</atoms>\n"
        "<bonds>\ns1\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-missing-coordinate-tag",
        "<structure number=\"39\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C\n</atoms>\n"
        "<bonds>\n\n</bonds>\n<not-xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-malformed-coordinate",
        "<structure number=\"40\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C\n</atoms>\n"
        "<bonds>\n\n</bonds>\n<xyz>\n0x;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-nested-error-cleanup",
        "<structure number=\"41\">\n"
        "<message type=\"error (no InChI)\" value=\"nested error\">\n"
        "<structure number=\"42\">\n</structure>\n</structure>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-success-structure-cleanup",
        "<structure number=\"43\">\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C\n</atoms>\n"
        "<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n</structure>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "plain-header-empty-label-negative-charge",
        "Structure: 44. =value\nAuxInfo=1/rA:1C-2/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-header-no-label",
        "Structure: 45. \nAuxInfo=1/rA:1C/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "xml-atom-parity-even",
        "<structure number=\"46\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1Cle\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-atom-parity-unknown",
        "<structure number=\"47\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1Clu\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-atom-parity-undefined",
        "<structure number=\"48\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1Cl?\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-tetrahedral-three-neighbors-valid",
        "<structure number=\"49\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n4CCCClo\n</atoms>\n<bonds>\n;;s1s2s3;\n</bonds>\n"
        "<xyz>\n;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-tetrahedral-four-neighbors-valid",
        "<structure number=\"50\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n5CCCCClo\n</atoms>\n<bonds>\n;;;s1s2s3s4;\n</bonds>\n"
        "<xyz>\n;;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-bond-parity-even-odd",
        "<structure number=\"51\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\nd+-1;\n</bonds>\n<xyz>\n;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-bond-parity-unknown-undefined",
        "<structure number=\"52\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\ndu?1;\n</bonds>\n<xyz>\n;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-wrong-bond-type",
        "<structure number=\"53\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\nx1;\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-empty-bonds-direct-end",
        "<structure number=\"54\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-empty-coordinates-direct-end",
        "<structure number=\"55\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-coordinate-leading-commas",
        "<structure number=\"56\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n,,;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-stereo-reallocation",
        .text = "<structure number=\"57\">\n<identifier.auxiliary-info>\n<reversibility>\n"
                "<atoms>\n3CloCloClo\n</atoms>\n<bonds>\nd--1;d--1d--2;\n</bonds>\n"
                "<xyz>\n;;;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .do_not_add_h = 1,
    },
    {
        .case_id = "xml-allocation-failure-stereo-reallocation",
        .text = "<structure number=\"58\">\n<identifier.auxiliary-info>\n<reversibility>\n"
                "<atoms>\n3CloCloClo\n</atoms>\n<bonds>\nd--1;d--1d--2;\n</bonds>\n"
                "<xyz>\n;;;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .do_not_add_h = 1,
        .allocation_failure_ordinal = 3,
    },
    {
        "xml-eof-before-bonds-header",
        "<structure number=\"59\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-eof-after-bonds-header",
        "<structure number=\"60\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-atom-trailing-data",
        "<structure number=\"61\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1CC\n</atoms>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-bond-trailing-data",
        "<structure number=\"62\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\ns1;;\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-coordinate-trailing-data",
        "<structure number=\"63\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-cross-line-sections",
        "<structure number=\"64\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n4C\nC\nC\nC\n</atoms>\n"
        "<bonds>\ns1;\ns1;\ns1;\n</bonds>\n"
        "<xyz>\n;\n;\n;\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-bond-parity-undefined-unknown",
        "<structure number=\"65\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\nd?u1;\n</bonds>\n<xyz>\n;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-explicit-valence",
        "<structure number=\"66\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C4\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "plain-caller-error-cleanup",
        .text = "AuxInfo=1/rA:2CloC/rB:x1;/rC:;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .caller_atom_capacity = 2,
        .caller_stereo_capacity = 2,
    },
    {
        .case_id = "plain-reverse-allene",
        .text = "AuxInfo=1/rA:3CCC/rB:d--1;d2;/rC:;;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .do_not_add_h = 1,
    },
    {
        .case_id = "xml-even-cumulene",
        .text = "<structure number=\"67\">\n<identifier.auxiliary-info>\n<reversibility>\n"
                "<atoms>\n4CCCC\n</atoms>\n<bonds>\nd--1;d2;d3;\n</bonds>\n"
                "<xyz>\n;;;;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .do_not_add_h = 1,
    },
    {
        "xml-wrong-bonds-header-eof",
        "<structure number=\"68\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<not-bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-wrong-bonds-header-next-structure",
        "<structure number=\"70\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<not-bonds>\n</structure>\n"
        "<structure number=\"71\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1O\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-structure-end-as-atom-data",
        "<structure number=\"72\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n</structure>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-nested-structure-as-atom-data",
        "<structure number=\"73\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n<structure number=\"74\">\n</structure>\n</structure>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "plain-caller-zero-passed-capacities",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .caller_atom_capacity = 1,
        .caller_stereo_capacity = 1,
        .pass_zero_atom_capacity = 1,
        .pass_zero_stereo_count = 1,
    },
    {
        .case_id = "plain-no-header-no-metadata-outputs",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_label_output = 1,
        .omit_value_output = 1,
        .omit_id_output = 1,
        .omit_flags_output = 1,
    },
    { "plain-empty-structure-header", "Structure:", INPUT_INCHI_PLAIN, 0 },
    { "plain-empty-label-value", "Structure: 75. LABEL =", INPUT_INCHI_PLAIN, 0 },
    {
        .case_id = "plain-missing-label-no-output",
        .text = "Structure: 76. SAMPLE is missing\nAuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_label_output = 1,
    },
    {
        "plain-skip-unrecognized-line",
        "not a recognized header\nAuxInfo=1/rA:1C/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    { "plain-empty-reversibility-no-newline", "AuxInfo=1//", INPUT_INCHI_PLAIN, 0 },
    { "plain-positive-atom-count-no-payload", "AuxInfo=1/rA:1", INPUT_INCHI_PLAIN, 0 },
    {
        "plain-extra-atom",
        "AuxInfo=1/rA:1CC/rB:/rC:;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-secondary-parity-default",
        "AuxInfo=1/rA:2CC/rB:d-x1;/rC:;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-tetrahedral-four-neighbors-valid",
        "AuxInfo=1/rA:5CCCCCle/rB:;;;s1s2s3s4;/rC:;;;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        .case_id = "xml-no-metadata-outputs",
        .text = "<structure number=\"77\" id.name=\"name\" id.value=\"value\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C\n</atoms>\n"
                "<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .omit_label_output = 1,
        .omit_value_output = 1,
        .omit_id_output = 1,
        .omit_flags_output = 1,
    },
    {
        "xml-number-not-terminated-at-endptr",
        "<structure number=\"78x\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-id-name-unterminated",
        "<structure number=\"79\" id.name=\"unterminated>\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C\n</atoms>\n"
        "<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-id-value-unterminated",
        "<structure number=\"80\" id.value=\"unterminated>\n"
        "<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C\n</atoms>\n"
        "<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-error-find-next-only",
        .text = "<structure number=\"81\">\n"
                "<message type=\"error (no InChI)\" value=\"hidden\">\n",
        .input_type = INPUT_INCHI_XML,
        .omit_atom_output = 1,
    },
    {
        "xml-error-unterminated",
        "<structure number=\"82\">\n"
        "<message type=\"error (no InChI)\" value=\"unterminated>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-atom-optional-numeric-suffixes",
        "<structure number=\"83\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C+.ihdt\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-secondary-parity-default",
        "<structure number=\"84\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\nd-x1;\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-wrong-bonds-header-initial-error",
        .text = "<structure number=\"85\">\n<identifier.auxiliary-info>\n<reversibility>\n"
                "<atoms>\n1C\n</atoms>\n<not-bonds>\n",
        .input_type = INPUT_INCHI_XML,
        .initial_error_code = INCHI_INP_ERROR_ERR,
    },
    { "plain-auxinfo-no-slash", "AuxInfo=1", INPUT_INCHI_PLAIN, 0 },
    {
        "plain-mixed-aromatic-single",
        "AuxInfo=1/rA:3CCC/rB:a1;s1;/rC:;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-cumulene-max-chain",
        "AuxInfo=1/rA:22CCCCCCCCCCCCCCCCCCCCCC/"
        "rB:d--1;d2;d3;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18;d19;d20;d21;/"
        "rC:;;;;;;;;;;;;;;;;;;;;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-nonzero-identical-bond-coordinates",
        "AuxInfo=1/rA:2CC/rB:s1;/rC:1,1,1;1,1,1;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "xml-negative-charge",
        "<structure number=\"86\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C-\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-cumulene-max-chain",
        "<structure number=\"87\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n22CCCCCCCCCCCCCCCCCCCCCC\n</atoms>\n<bonds>\n"
        "d--1;d2;d3;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18;d19;d20;d21;\n"
        "</bonds>\n<xyz>\n;;;;;;;;;;;;;;;;;;;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-eof-before-coordinate-header",
        "<structure number=\"88\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-eof-after-coordinate-header",
        "<structure number=\"89\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-caller-capacity-reused-valid",
        .text = "<structure number=\"90\">\n<identifier.auxiliary-info>\n<reversibility>\n"
                "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .caller_atom_capacity = 2,
        .caller_stereo_capacity = 2,
    },
    {
        "xml-nonzero-identical-bond-coordinates",
        "<structure number=\"91\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\ns1;\n</bonds>\n"
        "<xyz>\n1,1,1;1,1,1;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-cross-buffer-all-sections",
        .input_type = INPUT_INCHI_XML,
        .generated_text = 3,
        .omit_atom_output = 1,
        .omit_stereo_output = 1,
    },
    {
        "xml-eof-mid-atoms",
        "<structure number=\"93\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2C",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-eof-mid-bonds",
        "<structure number=\"94\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n3CCC\n</atoms>\n<bonds>\ns1;",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-eof-mid-coordinates",
        "<structure number=\"95\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\ns1;\n</bonds>\n<xyz>\n;",
        INPUT_INCHI_XML,
        0,
    },
    { .case_id = "xml-long-open-atoms", .input_type = INPUT_INCHI_XML, .generated_text = 4 },
    { .case_id = "xml-long-open-bonds", .input_type = INPUT_INCHI_XML, .generated_text = 5 },
    { .case_id = "xml-long-open-coordinates", .input_type = INPUT_INCHI_XML, .generated_text = 6 },
    { .case_id = "plain-long-open-atoms", .input_type = INPUT_INCHI_PLAIN, .generated_text = 7 },
    { .case_id = "plain-long-open-bonds", .input_type = INPUT_INCHI_PLAIN, .generated_text = 8 },
    { .case_id = "plain-long-open-coordinates", .input_type = INPUT_INCHI_PLAIN, .generated_text = 9 },
    {
        .case_id = "plain-eof-initial-error",
        .text = "",
        .input_type = INPUT_INCHI_PLAIN,
        .initial_error_code = INCHI_INP_ERROR_ERR,
    },
    {
        "plain-stereo-simple-endpoints",
        "AuxInfo=1/rA:2CC/rB:d--1;/rC:;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-stereo-degree-two-single-double",
        "AuxInfo=1/rA:3CCC/rB:s1;d--1;/rC:;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-stereo-middle-isotope-h",
        "AuxInfo=1/rA:3CC4hC/rB:d--1;d2;/rC:;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-stereo-middle-isotope-h-reverse",
        "AuxInfo=1/rA:3CC4hC/rB:d1;d--2;/rC:;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-stereo-later-interior-isotope-h",
        "AuxInfo=1/rA:4CCC4hC/rB:d--1;d2;d3;/rC:;;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "plain-stereo-central-cumulene",
        "AuxInfo=1/rA:4CCCC/rB:d1;d--2;d3;/rC:;;;;\n",
        INPUT_INCHI_PLAIN,
        0,
    },
    {
        "xml-stereo-simple-endpoints",
        "<structure number=\"99\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n2CC\n</atoms>\n<bonds>\nd--1;\n</bonds>\n<xyz>\n;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-stereo-degree-two-single-double",
        "<structure number=\"100\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n3CCC\n</atoms>\n<bonds>\ns1;d--1;\n</bonds>\n<xyz>\n;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-stereo-middle-isotope-h",
        "<structure number=\"101\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n3CC4hC\n</atoms>\n<bonds>\nd--1;d2;\n</bonds>\n<xyz>\n;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-stereo-middle-isotope-h-reverse",
        "<structure number=\"103\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n3CC4hC\n</atoms>\n<bonds>\nd1;d--2;\n</bonds>\n<xyz>\n;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-stereo-later-interior-isotope-h",
        "<structure number=\"104\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n4CCC4hC\n</atoms>\n<bonds>\nd--1;d2;d3;\n</bonds>\n<xyz>\n;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        "xml-stereo-central-cumulene",
        "<structure number=\"102\">\n<identifier.auxiliary-info>\n<reversibility>\n"
        "<atoms>\n4CCCC\n</atoms>\n<bonds>\nd1;d--2;d3;\n</bonds>\n<xyz>\n;;;;\n</xyz>\n",
        INPUT_INCHI_XML,
        0,
    },
    {
        .case_id = "xml-success-file-mode",
        .text = "<structure number=\"69\">\n<identifier.auxiliary-info>\n<reversibility>\n"
                "<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
        .file_mode = 1,
    },
    {
        "invalid-input-type",
        "not an InChI input\n",
        INPUT_NONE,
        0,
    },
};

static void print_json_string(const char *text)
{
    const unsigned char *p = (const unsigned char *) text;
    putchar('"');
    while (*p)
    {
        switch (*p)
        {
            case '"':
                fputs("\\\"", stdout);
                break;
            case '\\':
                fputs("\\\\", stdout);
                break;
            case '\b':
                fputs("\\b", stdout);
                break;
            case '\f':
                fputs("\\f", stdout);
                break;
            case '\n':
                fputs("\\n", stdout);
                break;
            case '\r':
                fputs("\\r", stdout);
                break;
            case '\t':
                fputs("\\t", stdout);
                break;
            default:
                if (*p < 0x20)
                {
                    printf("\\u%04x", (unsigned int) *p);
                }
                else
                {
                    putchar(*p);
                }
                break;
        }
        p++;
    }
    putchar('"');
}

static uint64_t double_bits(double value)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static double double_from_bits(uint64_t bits)
{
    double value;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

static void print_u8_array(const unsigned char *values, int length)
{
    int i;
    putchar('[');
    for (i = 0; i < length; i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("%u", (unsigned int) values[i]);
    }
    putchar(']');
}

static void print_ulong_array(const unsigned long *values, int length)
{
    int index;
    putchar('[');
    for (index = 0; index < length; index++)
    {
        if (index)
        {
            putchar(',');
        }
        printf("%" PRIu64, (uint64_t) values[index]);
    }
    putchar(']');
}

static void print_write_coord_record(const char *case_id, double value)
{
    unsigned char output[64];
    memset(output, 0xa5, sizeof(output));
    WriteCoord((char *) output, value);

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
          "\"api_version\":\"1.07.5\"},\"case_id\":", stdout);
    print_json_string(case_id);
    printf(",\"operation\":\"write_coord\",\"input\":{\"bits\":%" PRIu64
           ",\"locale\":\"C\"},\"output\":{\"bytes\":",
           double_bits(value));
    print_u8_array(output, (int) sizeof(output));
    fputs("}}\n", stdout);
}

static int print_write_coord_records(void)
{
    static const struct
    {
        const char *case_id;
        double value;
    } thresholds[] = {
        { "negative-scientific", -9999999.9 },
        { "negative-fixed-2", -999999.99 },
        { "negative-fixed-3", -99999.999 },
        { "positive-fixed-4", 99999.9999 },
        { "positive-fixed-3", 999999.999 },
        { "positive-fixed-2", 9999999.99 },
        { "positive-fixed-1", 99999999.9 },
    };
    static const struct
    {
        const char *case_id;
        uint64_t bits;
    } special[] = {
        { "positive-zero", UINT64_C(0x0000000000000000) },
        { "negative-zero", UINT64_C(0x8000000000000000) },
        { "positive-min-subnormal", UINT64_C(0x0000000000000001) },
        { "negative-min-subnormal", UINT64_C(0x8000000000000001) },
        { "positive-max-subnormal", UINT64_C(0x000fffffffffffff) },
        { "negative-max-subnormal", UINT64_C(0x800fffffffffffff) },
        { "positive-min-normal", UINT64_C(0x0010000000000000) },
        { "negative-min-normal", UINT64_C(0x8010000000000000) },
        { "positive-max-finite", UINT64_C(0x7fefffffffffffff) },
        { "negative-max-finite", UINT64_C(0xffefffffffffffff) },
        { "positive-infinity", UINT64_C(0x7ff0000000000000) },
        { "negative-infinity", UINT64_C(0xfff0000000000000) },
        { "positive-quiet-nan", UINT64_C(0x7ff8000000000000) },
        { "negative-quiet-nan", UINT64_C(0xfff8000000000000) },
        { "positive-signaling-nan", UINT64_C(0x7ff0000000000001) },
        { "negative-signaling-nan", UINT64_C(0xfff0000000000001) },
    };
    static const struct
    {
        const char *case_id;
        double value;
    } rounding[] = {
        { "fixed-4-nearest-below-even", 1.23445 },
        { "fixed-4-nearest-above-even", 1.23455 },
        { "fixed-4-negative-nearest-below-even", -1.23445 },
        { "fixed-4-negative-nearest-above-even", -1.23455 },
        { "negative-fixed-3-nearest", -100000.0005 },
        { "positive-fixed-3-nearest", 100000.0005 },
        { "negative-fixed-2-nearest", -1000000.005 },
        { "positive-fixed-2-nearest", 1000000.005 },
        { "positive-fixed-1-nearest", 10000000.05 },
        { "negative-scientific-2-exact-tie-even", -10050000.0 },
        { "negative-scientific-2-exact-tie-odd", -10150000.0 },
        { "positive-scientific-3-exact-tie-even", 100050000.0 },
        { "positive-scientific-3-exact-tie-odd", 100150000.0 },
    };
    size_t i;
    char case_id[96];
    uint64_t random_bits = UINT64_C(0x6a09e667f3bcc909);

    if (!setlocale(LC_NUMERIC, "C") || strcmp(setlocale(LC_NUMERIC, NULL), "C") != 0)
    {
        fprintf(stderr, "failed to select the C numeric locale\n");
        return 67;
    }

    for (i = 0; i < sizeof(thresholds) / sizeof(thresholds[0]); i++)
    {
        snprintf(case_id, sizeof(case_id), "%s-predecessor", thresholds[i].case_id);
        print_write_coord_record(case_id, nextafter(thresholds[i].value, -INFINITY));
        snprintf(case_id, sizeof(case_id), "%s-equality", thresholds[i].case_id);
        print_write_coord_record(case_id, thresholds[i].value);
        snprintf(case_id, sizeof(case_id), "%s-successor", thresholds[i].case_id);
        print_write_coord_record(case_id, nextafter(thresholds[i].value, INFINITY));
    }
    for (i = 0; i < sizeof(special) / sizeof(special[0]); i++)
    {
        print_write_coord_record(special[i].case_id, double_from_bits(special[i].bits));
    }
    for (i = 0; i < sizeof(rounding) / sizeof(rounding[0]); i++)
    {
        double value = rounding[i].value;
        snprintf(case_id, sizeof(case_id), "%s-predecessor", rounding[i].case_id);
        print_write_coord_record(case_id, nextafter(value, -INFINITY));
        snprintf(case_id, sizeof(case_id), "%s-value", rounding[i].case_id);
        print_write_coord_record(case_id, value);
        snprintf(case_id, sizeof(case_id), "%s-successor", rounding[i].case_id);
        print_write_coord_record(case_id, nextafter(value, INFINITY));
    }
    for (i = 0; i < 4096; i++)
    {
        random_bits ^= random_bits << 13;
        random_bits ^= random_bits >> 7;
        random_bits ^= random_bits << 17;
        snprintf(case_id, sizeof(case_id), "deterministic-bits-%04zu", i);
        print_write_coord_record(case_id, double_from_bits(random_bits));
    }
    return 0;
}

static void print_parse_options_record(const char *case_id,
                                       const unsigned char *input,
                                       size_t input_length,
                                       int maxargs)
{
    unsigned char command[512];
    const char *argv[40];
    int count;
    int i;

    if (input_length + 16 > sizeof(command) || maxargs > (int) (sizeof(argv) / sizeof(argv[0])))
    {
        fprintf(stderr, "parse_options_string oracle fixture is too large\n");
        exit(68);
    }
    memset(command, 0xa5, sizeof(command));
    memcpy(command, input, input_length);
    for (i = 0; i < (int) (sizeof(argv) / sizeof(argv[0])); i++)
    {
        argv[i] = (const char *) (uintptr_t) 1;
    }

    count = parse_options_string((char *) command, argv, maxargs);

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
          "\"api_version\":\"1.07.5\"},\"case_id\":", stdout);
    print_json_string(case_id);
    printf(",\"operation\":\"parse_options_string\",\"input\":{\"maxargs\":%d,"
           "\"bytes\":",
           maxargs);
    print_u8_array(input, (int) input_length);
    printf("},\"output\":{\"count\":%d,\"argv_offsets\":[", count);
    for (i = 0; i <= count; i++)
    {
        intptr_t offset;
        if (i)
        {
            putchar(',');
        }
        if (argv[i] == NULL)
        {
            offset = -1;
        }
        else if (i == 0)
        {
            offset = -2;
        }
        else
        {
            offset = (intptr_t) ((uintptr_t) argv[i] - (uintptr_t) command);
        }
        printf("%" PRIdPTR, offset);
    }
    fputs("],\"argv_bytes\":[", stdout);
    for (i = 0; i < count; i++)
    {
        if (i)
        {
            putchar(',');
        }
        print_u8_array((const unsigned char *) argv[i], (int) strlen(argv[i]));
    }
    fputs("],\"buffer_bytes\":", stdout);
    print_u8_array(command, (int) input_length + 16);
    printf(",\"terminal_null\":%s}}\n", argv[count] == NULL ? "true" : "false");
}

static int print_parse_options_records(void)
{
    static const struct
    {
        const char *case_id;
        const char *input;
        int maxargs;
    } cases[] = {
        { "empty", "", 8 },
        { "whitespace", "  \t \t", 8 },
        { "unquoted", "alpha beta\tgamma", 8 },
        { "quoted", "\"alpha beta\" \"gamma\tdelta\"", 8 },
        { "empty-quoted", "\"\" tail", 8 },
        { "quoted-unquoted-concatenation", "ab\"cd ef\"gh", 8 },
        { "consecutive-quotes", "\"a\"\"b\" \"\"\"\"", 8 },
        { "one-backslash-before-quote", "a\\\"b c", 8 },
        { "two-backslashes-before-quote", "a\\\\\"b c\"", 8 },
        { "three-backslashes-before-quote", "a\\\\\\\"b c", 8 },
        { "four-backslashes-before-quote", "a\\\\\\\\\"b c\"", 8 },
        { "trailing-backslash", "alpha\\", 8 },
        { "unterminated-quote", "\"alpha beta", 8 },
        { "leading-and-trailing-space", " \talpha beta  ", 8 },
        { "maxargs-two", "one two", 2 },
        { "maxargs-three-saturation", "one two three", 3 },
        { "maxargs-eight-saturation", "0 1 2 3 4 5 6 7 8", 8 },
        { "maxargs-thirty-two-saturation",
          "00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 "
          "20 21 22 23 24 25 26 27 28 29 30 31 32 33 34",
          32 },
    };
    size_t i;
    for (i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
    {
        print_parse_options_record(cases[i].case_id,
                                   (const unsigned char *) cases[i].input,
                                   strlen(cases[i].input) + 1,
                                   cases[i].maxargs);
    }
    return 0;
}

static void print_element_lookup_record(const char *case_id,
                                        const unsigned char *symbol,
                                        size_t symbol_length)
{
    int result = el_number_in_internal_ref_table((const char *) symbol);
    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
          "\"api_version\":\"1.07.5\"},\"case_id\":", stdout);
    print_json_string(case_id);
    fputs(",\"operation\":\"el_number_in_internal_ref_table\",\"input\":{\"symbol_bytes\":",
          stdout);
    print_u8_array(symbol, (int) symbol_length);
    printf("},\"output\":{\"result\":%d,\"n_el_data_len\":%d,\"err_elem\":%d}}\n",
           result, nElDataLen, ERR_ELEM);
}

static int print_element_lookup_records(void)
{
    unsigned char symbol[4];
    char case_id[48];
    int first;
    int second;
    static const struct
    {
        const char *case_id;
        unsigned char symbol[8];
        size_t length;
    } boundary_cases[] = {
        { "empty-sentinel", { 0 }, 1 },
        { "hydrogen", { 'H', 0 }, 2 },
        { "deuterium", { 'D', 0 }, 2 },
        { "tritium", { 'T', 0 }, 2 },
        { "pseudo-zy", { 'Z', 'y', 0 }, 3 },
        { "pseudo-zz", { 'Z', 'z', 0 }, 3 },
        { "inactive-pseudo-zu", { 'Z', 'u', 0 }, 3 },
        { "inactive-pseudo-zv", { 'Z', 'v', 0 }, 3 },
        { "inactive-pseudo-zw", { 'Z', 'w', 0 }, 3 },
        { "inactive-pseudo-zx", { 'Z', 'x', 0 }, 3 },
        { "lowercase", { 'h', 0 }, 2 },
        { "long-name", { 'C', 'a', 'r', 'b', 'o', 'n', 0 }, 7 },
        { "digit-suffix", { 'C', '1', 0 }, 3 },
        { "question-marks", { '?', '?', 0 }, 3 },
        { "high-byte", { 0xff, 0 }, 2 },
        { "embedded-nul-tail", { 'C', 0, 'l', 0 }, 4 },
    };
    size_t i;

    symbol[1] = 0;
    for (first = 'A'; first <= 'Z'; first++)
    {
        symbol[0] = (unsigned char) first;
        snprintf(case_id, sizeof(case_id), "candidate-%c", first);
        print_element_lookup_record(case_id, symbol, 2);
    }
    symbol[2] = 0;
    for (first = 'A'; first <= 'Z'; first++)
    {
        symbol[0] = (unsigned char) first;
        for (second = 'a'; second <= 'z'; second++)
        {
            symbol[1] = (unsigned char) second;
            snprintf(case_id, sizeof(case_id), "candidate-%c%c", first, second);
            print_element_lookup_record(case_id, symbol, 3);
        }
    }
    for (i = 0; i < sizeof(boundary_cases) / sizeof(boundary_cases[0]); i++)
    {
        print_element_lookup_record(boundary_cases[i].case_id,
                                    boundary_cases[i].symbol,
                                    boundary_cases[i].length);
    }
    return 0;
}

static void print_periodic_lookup_record(const char *case_id,
                                         const unsigned char *symbol,
                                         size_t symbol_length)
{
    int result = get_periodic_table_number((const char *) symbol);
    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
          "\"api_version\":\"1.07.5\"},\"case_id\":", stdout);
    print_json_string(case_id);
    fputs(",\"operation\":\"get_periodic_table_number\",\"input\":{\"symbol_bytes\":",
          stdout);
    print_u8_array(symbol, (int) symbol_length);
    printf(",\"is_null\":false},\"output\":{\"result\":%d,\"err_elem\":%d}}\n",
           result, ERR_ELEM);
}

static int print_periodic_lookup_records(void)
{
    unsigned char symbol[4];
    char case_id[48];
    int first;
    int second;
    static const struct
    {
        const char *case_id;
        unsigned char symbol[8];
        size_t length;
    } boundary_cases[] = {
        { "empty-sentinel", { 0 }, 1 },
        { "hydrogen", { 'H', 0 }, 2 },
        { "deuterium", { 'D', 0 }, 2 },
        { "tritium", { 'T', 0 }, 2 },
        { "pseudo-zy", { 'Z', 'y', 0 }, 3 },
        { "pseudo-zz", { 'Z', 'z', 0 }, 3 },
        { "inactive-pseudo-zu", { 'Z', 'u', 0 }, 3 },
        { "inactive-pseudo-zv", { 'Z', 'v', 0 }, 3 },
        { "inactive-pseudo-zw", { 'Z', 'w', 0 }, 3 },
        { "inactive-pseudo-zx", { 'Z', 'x', 0 }, 3 },
        { "lowercase", { 'h', 0 }, 2 },
        { "long-name", { 'C', 'a', 'r', 'b', 'o', 'n', 0 }, 7 },
        { "digit-suffix", { 'C', '1', 0 }, 3 },
        { "question-marks", { '?', '?', 0 }, 3 },
        { "high-byte", { 0xff, 0 }, 2 },
        { "embedded-nul-tail", { 'C', 0, 'l', 0 }, 4 },
    };
    size_t i;

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
          "\"api_version\":\"1.07.5\"},\"case_id\":\"null\"," 
          "\"operation\":\"get_periodic_table_number\",\"input\":{\"symbol_bytes\":[],"
          "\"is_null\":true},\"output\":{\"result\":255,\"err_elem\":255}}\n", stdout);

    symbol[1] = 0;
    for (first = 'A'; first <= 'Z'; first++)
    {
        symbol[0] = (unsigned char) first;
        snprintf(case_id, sizeof(case_id), "candidate-%c", first);
        print_periodic_lookup_record(case_id, symbol, 2);
    }
    symbol[2] = 0;
    for (first = 'A'; first <= 'Z'; first++)
    {
        symbol[0] = (unsigned char) first;
        for (second = 'a'; second <= 'z'; second++)
        {
            symbol[1] = (unsigned char) second;
            snprintf(case_id, sizeof(case_id), "candidate-%c%c", first, second);
            print_periodic_lookup_record(case_id, symbol, 3);
        }
    }
    for (i = 0; i < sizeof(boundary_cases) / sizeof(boundary_cases[0]); i++)
    {
        print_periodic_lookup_record(boundary_cases[i].case_id,
                                     boundary_cases[i].symbol,
                                     boundary_cases[i].length);
    }
    return 0;
}

static int print_el_valence_records(void)
{
    int periodic_values[124];
    int periodic_count = 0;
    int periodic;
    int charge;
    int slot;
    int i;

    periodic_values[periodic_count++] = INT32_MIN;
    periodic_values[periodic_count++] = -1;
    periodic_values[periodic_count++] = 0;
    periodic_values[periodic_count++] = 1;
    for (periodic = 2; periodic <= 121; periodic++)
    {
        periodic_values[periodic_count++] = periodic;
    }
    for (i = 0; i < periodic_count; i++)
    {
        periodic = periodic_values[i];
        for (charge = -3; charge <= 3; charge++)
        {
            for (slot = 0; slot <= 5; slot++)
            {
                int result = get_el_valence(periodic, charge, slot);
                printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
                       "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
                       "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
                       "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
                       "\"api_version\":\"1.07.5\"},"
                       "\"case_id\":\"periodic-%d-charge-%d-slot-%d\"," 
                       "\"operation\":\"get_el_valence\"," 
                       "\"input\":{\"periodic_number\":%d,\"charge\":%d,\"slot\":%d},"
                       "\"output\":{\"result\":%d}}\n",
                       periodic, charge, slot, periodic, charge, slot, result);
            }
        }
    }
    return 0;
}

static int print_metal_records(void)
{
    int periodic_number;
    for (periodic_number = -1; periodic_number <= 121; periodic_number++)
    {
        int result = is_el_a_metal(periodic_number);
        printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
               "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
               "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
               "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
               "\"api_version\":\"1.07.5\"},"
               "\"case_id\":\"periodic-%d\",\"operation\":\"is_el_a_metal\"," 
               "\"input\":{\"periodic_number\":%d},\"output\":{\"result\":%d}}\n",
               periodic_number,
               periodic_number,
               result);
    }
    return 0;
}

static void print_atomic_mass_record(int atomic_number)
{
    int result = get_atomic_mass_from_elnum(atomic_number);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
           "\"api_version\":\"1.07.5\"},"
           "\"case_id\":\"atomic-number-%d\",\"operation\":\"get_atomic_mass_from_elnum\"," 
           "\"input\":{\"atomic_number\":%d},\"output\":{\"result\":%d}}\n",
           atomic_number,
           atomic_number,
           result);
}

static int print_atomic_mass_records(void)
{
    int atomic_number;
    print_atomic_mass_record(INT32_MIN + 1);
    for (atomic_number = -10; atomic_number <= 130; atomic_number++)
    {
        print_atomic_mass_record(atomic_number);
    }
    print_atomic_mass_record(INT32_MAX - 1);
    return 0;
}

typedef struct tagDetectUnusualValenceCase
{
    const char *case_id;
    int periodic_number;
    int charge;
    int radical;
    int bonds_valence;
    int hydrogen_count;
    int bond_count;
} DETECT_UNUSUAL_VALENCE_CASE;

static void print_detect_unusual_valence_record(const char *case_id,
                                                int periodic_number,
                                                int charge,
                                                int radical,
                                                int bonds_valence,
                                                int hydrogen_count,
                                                int bond_count)
{
    int result = detect_unusual_el_valence(periodic_number,
                                           charge,
                                           radical,
                                           bonds_valence,
                                           hydrogen_count,
                                           bond_count);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
           "\"api_version\":\"1.07.5\"},\"case_id\":");
    print_json_string(case_id);
    printf(",\"operation\":\"detect_unusual_el_valence\"," 
           "\"input\":{\"periodic_number\":%d,\"charge\":%d,"
           "\"radical\":%d,\"bonds_valence\":%d,"
           "\"hydrogen_count\":%d,\"bond_count\":%d},"
           "\"output\":{\"result\":%d}}\n",
           periodic_number,
           charge,
           radical,
           bonds_valence,
           hydrogen_count,
           bond_count,
           result);
}

static int print_detect_unusual_valence_records(void)
{
    static const int radicals[] = {
        INT32_MIN, -1, 0, 1, 2, 3, 4, INT32_MAX,
    };
    static const DETECT_UNUSUAL_VALENCE_CASE branch_cases[] = {
        {"early-valid", INT32_MAX, 0, INT32_MAX, INT32_MAX, 0, 0},
        {"invalid-charge-equal-min", INT32_MAX, INT32_MIN, INT32_MIN, INT32_MIN, 0, INT32_MIN},
        {"invalid-charge-different-max", INT32_MAX, INT32_MAX, INT32_MAX, INT32_MAX, INT32_MIN, INT32_MIN},
        {"zero-first-all-single", 2, 0, 0, 1, 0, 1},
        {"nonzero-first-all-single", 1, 0, 0, 1, 0, 1},
        {"later-slot-match", 7, 0, 0, 0, 5, 1},
        {"no-match", 6, 0, 0, 0, 2, 1},
        {"radical-default-negative", 6, 0, -1, 0, 4, 1},
        {"radical-default-zero", 6, 0, 0, 0, 4, 1},
        {"radical-default-positive", 6, 0, 4, 0, 4, 1},
        {"radical-doublet", 6, 0, 2, 0, 3, 1},
        {"radical-triplet", 6, 0, 3, 0, 2, 1},
        {"radical-singlet", 6, 0, 1, 0, 2, 1},
        {"known-nonpositive", 2, 0, 3, 0, 1, 1},
        {"sum-min", 6, 0, 0, INT32_MIN, 0, 1},
        {"sum-max", 6, 0, 0, INT32_MAX, 0, 1},
        {"sum-min-plus-max", 6, 0, 0, INT32_MIN, INT32_MAX, 1},
        {"sum-max-plus-min", 6, 0, 0, INT32_MAX, INT32_MIN, 1},
        {"periodic-sentinel", 121, 0, 0, 0, 7, 1},
        {"periodic-negative-h-row", INT32_MIN, 0, 0, 0, 1, 1},
    };
    int periodic_values[124];
    int periodic_count = 0;
    int periodic;
    int charge;
    int radical_index;
    int chemical_valence;
    int i;
    char case_id[160];

    periodic_values[periodic_count++] = INT32_MIN;
    periodic_values[periodic_count++] = -1;
    periodic_values[periodic_count++] = 0;
    periodic_values[periodic_count++] = 1;
    for (periodic = 2; periodic <= 121; periodic++)
    {
        periodic_values[periodic_count++] = periodic;
    }
    for (i = 0; i < periodic_count; i++)
    {
        periodic = periodic_values[i];
        for (charge = -3; charge <= 3; charge++)
        {
            for (radical_index = 0;
                 radical_index < (int) (sizeof(radicals) / sizeof(radicals[0]));
                 radical_index++)
            {
                for (chemical_valence = -2; chemical_valence <= 16; chemical_valence++)
                {
                    snprintf(case_id,
                             sizeof(case_id),
                             "table-periodic-%d-charge-%d-radical-%d-valence-%d",
                             periodic,
                             charge,
                             radicals[radical_index],
                             chemical_valence);
                    print_detect_unusual_valence_record(case_id,
                                                        periodic,
                                                        charge,
                                                        radicals[radical_index],
                                                        0,
                                                        chemical_valence,
                                                        1);
                }
            }
        }
    }
    for (i = 0; i < (int) (sizeof(branch_cases) / sizeof(branch_cases[0])); i++)
    {
        const DETECT_UNUSUAL_VALENCE_CASE *test_case = branch_cases + i;
        print_detect_unusual_valence_record(test_case->case_id,
                                            test_case->periodic_number,
                                            test_case->charge,
                                            test_case->radical,
                                            test_case->bonds_valence,
                                            test_case->hydrogen_count,
                                            test_case->bond_count);
    }
    return 0;
}

static void print_i16_array(const AT_NUM *values, int length)
{
    int i;
    putchar('[');
    for (i = 0; i < length; i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("%d", (int) values[i]);
    }
    putchar(']');
}

static void print_i8_array(const S_CHAR *values, int length)
{
    int i;
    putchar('[');
    for (i = 0; i < length; i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("%d", (int) values[i]);
    }
    putchar(']');
}

static void print_char_array(const char *values, int length)
{
    int i;
    putchar('[');
    for (i = 0; i < length; i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("%d", (int) values[i]);
    }
    putchar(']');
}

typedef struct tagExtractChargeCase
{
    const char *case_id;
    const char *text;
} EXTRACT_CHARGE_CASE;

static void print_extract_charge_record(const char *case_id, const char *text)
{
    char buffer[128];
    int radical = -1901;
    int charge = -1902;
    int status;
    size_t length = strlen(text);

    if (length + 1 > sizeof(buffer))
    {
        fprintf(stderr, "extract charge oracle input is too long\n");
        exit(67);
    }
    memset(buffer, 0x55, sizeof(buffer));
    memcpy(buffer, text, length + 1);
    status = extract_charges_and_radicals(buffer, &radical, &charge);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
           "\"api_version\":\"1.07.5\"},\"case_id\":");
    print_json_string(case_id);
    fputs(",\"operation\":\"extract_charges_and_radicals\",\"input\":{\"text\":", stdout);
    print_json_string(text);
    printf("},\"output\":{\"status\":%d,\"radical\":%d,\"charge\":%d,"
           "\"nul_offset\":%d,\"buffer\":",
           status,
           radical,
           charge,
           (int) strlen(buffer));
    print_char_array(buffer, (int) sizeof(buffer));
    fputs("}}\n", stdout);
}

static int print_extract_charge_records(void)
{
    static const EXTRACT_CHARGE_CASE cases[] = {
        {"empty", ""},
        {"plain", "C"},
        {"only-plus", "+"},
        {"only-minus", "-"},
        {"only-caret", "^"},
        {"only-colon", ":"},
        {"only-dot", "."},
        {"only-two-dots", ".."},
        {"only-three-dots", "..."},
        {"only-four-dots", "...."},
        {"single-plus", "C+"},
        {"single-minus", "C-"},
        {"two-plus", "C++"},
        {"three-minus", "C---"},
        {"plus-minus", "C+-"},
        {"minus-plus", "C-+"},
        {"four-plus", "C++++"},
        {"four-minus", "C----"},
        {"plus-zero", "C+0"},
        {"minus-zero", "C-0"},
        {"plus-one", "C+1"},
        {"minus-one", "C-1"},
        {"plus-two", "C+2"},
        {"minus-two", "C-2"},
        {"plus-int-max", "C+2147483647"},
        {"minus-int-max", "C-2147483647"},
        {"narrow-minus-int-max", "C+2147483649"},
        {"narrow-minus-one", "C+4294967295"},
        {"narrow-zero", "C+4294967296"},
        {"narrow-one", "C+4294967297"},
        {"long-max", "C+9223372036854775807"},
        {"long-positive-overflow", "C+9223372036854775808"},
        {"long-min", "C+ -9223372036854775808"},
        {"long-negative-overflow", "C+ -9223372036854775809"},
        {"decimal-stop", "C+12x"},
        {"space-decimal-stop", "C+ 12x"},
        {"tab-negative", "C+\t-2"},
        {"malformed-alpha", "C+abc"},
        {"malformed-space-alpha", "C+  abc"},
        {"malformed-sign", "C+ +x"},
        {"repeated-fields", "C+2H-3"},
        {"repeated-same-sign", "C+2+3"},
        {"caret-separated", "C^x^^"},
        {"caret-dot", "C^."},
        {"two-caret-dot", "C^^."},
        {"two-caret-colon", "C^^:"},
        {"colon-middle", "C:O"},
        {"dot-middle", "C.O"},
        {"caret-then-plus", "C^+2"},
        {"plus-then-caret", "C+2^"},
        {"charge-caret-dots", "C-2^.."},
        {"signs-after-letters", "Fe+2H-3^^:"},
    };
    static const char *suffixes[] = {
        "", "0", "1", "2", "12", "x", " 2", "\t2",
    };
    char text[128];
    char case_id[128];
    int case_index;
    int length;
    int pattern;
    int suffix_index;
    int i;

    for (i = 0; i < (int) (sizeof(cases) / sizeof(cases[0])); i++)
    {
        print_extract_charge_record(cases[i].case_id, cases[i].text);
    }
    case_index = 0;
    for (length = 1; length <= 6; length++)
    {
        for (pattern = 0; pattern < (1 << length); pattern++)
        {
            for (suffix_index = 0;
                 suffix_index < (int) (sizeof(suffixes) / sizeof(suffixes[0]));
                 suffix_index++)
            {
                int cursor = 0;
                text[cursor++] = 'C';
                for (i = 0; i < length; i++)
                {
                    text[cursor++] = (pattern & (1 << i)) ? '-' : '+';
                }
                strcpy(text + cursor, suffixes[suffix_index]);
                snprintf(case_id,
                         sizeof(case_id),
                         "generated-%04d-length-%d-pattern-%d-suffix-%d",
                         case_index++,
                         length,
                         pattern,
                         suffix_index);
                print_extract_charge_record(case_id, text);
            }
        }
    }
    return 0;
}

typedef struct tagExtractHydrogenCase
{
    const char *case_id;
    const char *text;
    S_CHAR initial_isotopes[3];
} EXTRACT_HYDROGEN_CASE;

static int print_extract_hydrogen_records(void)
{
    static const EXTRACT_HYDROGEN_CASE cases[] = {
        {"empty", "", {7, 0, 0}},
        {"plain", "C", {7, 0, 0}},
        {"h", "H", {7, 0, 0}},
        {"d", "D", {7, 0, 0}},
        {"t", "T", {7, 0, 0}},
        {"helium", "He", {7, 0, 0}},
        {"d-lower", "Dh", {7, 0, 0}},
        {"t-lower", "Th", {7, 0, 0}},
        {"h-zero", "H0", {7, 0, 0}},
        {"h-one", "H1", {7, 0, 0}},
        {"h-two", "H2", {7, 0, 0}},
        {"d-two", "D2", {7, 0, 0}},
        {"t-two", "T2", {7, 0, 0}},
        {"h-repeat", "H2H3", {7, 0, 0}},
        {"d-repeat", "D2D3", {7, 0, 0}},
        {"t-repeat", "T2T3", {7, 0, 0}},
        {"bare-mixed", "HDT", {7, 0, 0}},
        {"counted-mixed", "H2D3T4", {7, 0, 0}},
        {"suffix-h", "CH2", {7, 0, 0}},
        {"prefix-h", "H2C", {7, 0, 0}},
        {"digit-before-h", "C2H", {7, 0, 0}},
        {"h-stop", "H12x", {7, 0, 0}},
        {"d-stop", "D12x", {7, 0, 0}},
        {"t-stop", "T12x", {7, 0, 0}},
        {"h-plus", "H+", {7, 0, 0}},
        {"h-minus", "H-", {7, 0, 0}},
        {"h-caret", "H^", {7, 0, 0}},
        {"h-dot", "H.", {7, 0, 0}},
        {"h-colon", "H:", {7, 0, 0}},
        {"alias-ph4d", "pH4d", {7, 0, 0}},
        {"alias-ah2b", "aH2b", {7, 0, 0}},
        {"two-h", "HH", {7, 0, 0}},
        {"two-d", "DD", {7, 0, 0}},
        {"two-t", "TT", {7, 0, 0}},
        {"leading-zero", "H00", {7, 0, 0}},
        {"leading-zero-two", "H0002", {7, 0, 0}},
        {"int-max", "H2147483647", {7, 0, 0}},
        {"narrow-int-min", "H2147483648", {7, 0, 0}},
        {"narrow-int-min-plus-one", "H2147483649", {7, 0, 0}},
        {"narrow-minus-one", "H4294967295", {7, 0, 0}},
        {"narrow-zero", "H4294967296", {7, 0, 0}},
        {"long-max", "H9223372036854775807", {7, 0, 0}},
        {"long-overflow", "H9223372036854775808", {7, 0, 0}},
        {"d-int-min", "D2147483648", {7, 0, 0}},
        {"t-minus-one", "T4294967295", {7, 0, 0}},
        {"mixed-boundary", "H4294967296D2147483649T4294967295", {7, 0, 0}},
        {"lower-h", "h", {7, 0, 0}},
        {"lower-d", "d", {7, 0, 0}},
        {"lower-t", "t", {7, 0, 0}},
        {"element-lower-h", "Nh", {7, 0, 0}},
        {"embedded-counts", "H2xD3yT4", {7, 0, 0}},
        {"unchanged-isotope-zero", "D0T0", {7, 11, 13}},
        {"d-wrap-positive", "D1", {7, 127, 0}},
        {"t-wrap-positive", "T255", {7, 0, -128}},
        {"all-initial-fields", "H2", {-7, -11, 13}},
        {"no-token-initial-fields", "Xe", {-7, -11, 13}},
    };
    char buffer[128];
    int i;

    for (i = 0; i < (int) (sizeof(cases) / sizeof(cases[0])); i++)
    {
        const EXTRACT_HYDROGEN_CASE *test_case = cases + i;
        S_CHAR isotopes[3];
        int hydrogen_count;
        size_t length = strlen(test_case->text);
        if (length + 1 > sizeof(buffer))
        {
            return 68;
        }
        memset(buffer, 0x55, sizeof(buffer));
        memcpy(buffer, test_case->text, length + 1);
        memcpy(isotopes, test_case->initial_isotopes, sizeof(isotopes));
        hydrogen_count = extract_H_atoms(buffer, isotopes);
        printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
               "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
               "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
               "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
               "\"api_version\":\"1.07.5\"},\"case_id\":");
        print_json_string(test_case->case_id);
        fputs(",\"operation\":\"extract_H_atoms\",\"input\":{\"text\":", stdout);
        print_json_string(test_case->text);
        fputs(",\"initial_isotopes\":", stdout);
        print_i8_array(test_case->initial_isotopes, 3);
        printf("},\"output\":{\"hydrogen_count\":%d,\"isotopes\":",
               hydrogen_count);
        print_i8_array(isotopes, 3);
        printf(",\"nul_offset\":%d,\"buffer\":", (int) strlen(buffer));
        print_char_array(buffer, (int) sizeof(buffer));
        fputs("}}\n", stdout);
    }
    return 0;
}

static void print_list_record(const char *case_id,
                              AT_NUMB *path,
                              int array_length,
                              AT_NUMB target,
                              int path_length)
{
    AT_NUMB *result = is_in_the_list(path, target, path_length);
    int offset = result ? (int) (result - path) : -1;
    int i;
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
           "\"api_version\":\"1.07.5\"},\"case_id\":");
    print_json_string(case_id);
    fputs(",\"operation\":\"is_in_the_list\",\"input\":{\"path\":[", stdout);
    for (i = 0; i < array_length; i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("%u", (unsigned int) path[i]);
    }
    printf("],\"target\":%u,\"path_length\":%d},\"output\":{\"offset\":%d}}\n",
           (unsigned int) target,
           path_length,
           offset);
}

static int print_list_records(void)
{
    static const AT_NUMB alphabet[] = {0, 1, 65535};
    static const AT_NUMB targets[] = {0, 1, 2, 65535};
    AT_NUMB path[5];
    char case_id[128];
    int length;
    int combination;
    int combinations;
    int target_index;
    int path_length;
    int i;

    print_list_record("null-zero", NULL, 0, 7, 0);
    for (length = 0; length <= 5; length++)
    {
        combinations = 1;
        for (i = 0; i < length; i++)
        {
            combinations *= 3;
        }
        for (combination = 0; combination < combinations; combination++)
        {
            int encoded = combination;
            for (i = 0; i < length; i++)
            {
                path[i] = alphabet[encoded % 3];
                encoded /= 3;
            }
            for (target_index = 0;
                 target_index < (int) (sizeof(targets) / sizeof(targets[0]));
                 target_index++)
            {
                for (path_length = 0; path_length <= length; path_length++)
                {
                    snprintf(case_id,
                             sizeof(case_id),
                             "length-%d-combination-%d-target-%u-prefix-%d",
                             length,
                             combination,
                             (unsigned int) targets[target_index],
                             path_length);
                    print_list_record(case_id,
                                      path,
                                      length,
                                      targets[target_index],
                                      path_length);
                }
            }
        }
    }
    return 0;
}

static void print_bonds_to_metal_record(const char *case_id,
                                        int valence,
                                        const unsigned char *elements,
                                        const unsigned char *bond_types)
{
    inp_ATOM atoms[21];
    int result;
    int i;
    memset(atoms, 0, sizeof(atoms));
    atoms[0].valence = (S_CHAR) valence;
    for (i = 0; i < valence; i++)
    {
        atoms[0].neighbor[i] = (AT_NUMB) (i + 1);
        atoms[0].bond_type[i] = bond_types[i];
        atoms[i + 1].el_number = elements[i];
    }
    result = nBondsValToMetal(atoms, 0);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
           "\"api_version\":\"1.07.5\"},\"case_id\":");
    print_json_string(case_id);
    printf(",\"operation\":\"nBondsValToMetal\",\"input\":{\"valence\":%d,\"elements\":[",
           valence);
    for (i = 0; i < (valence > 0 ? valence : 0); i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("%u", (unsigned int) elements[i]);
    }
    fputs("],\"bond_types\":[", stdout);
    for (i = 0; i < (valence > 0 ? valence : 0); i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("%u", (unsigned int) bond_types[i]);
    }
    printf("]},\"output\":{\"result\":%d}}\n", result);
}

static int print_bonds_to_metal_records(void)
{
    unsigned char elements[20];
    unsigned char bond_types[20];
    char case_id[96];
    int element;
    int bond_type;
    int i;

    for (element = 0; element <= 121; element++)
    {
        for (bond_type = 0; bond_type <= 255; bond_type++)
        {
            elements[0] = (unsigned char) element;
            bond_types[0] = (unsigned char) bond_type;
            snprintf(case_id,
                     sizeof(case_id),
                     "single-element-%d-bond-%d",
                     element,
                     bond_type);
            print_bonds_to_metal_record(case_id, 1, elements, bond_types);
        }
    }
    print_bonds_to_metal_record("negative-valence", -1, elements, bond_types);
    print_bonds_to_metal_record("zero-valence", 0, elements, bond_types);
    {
        static const unsigned char mixed_elements[] = {3, 6, 26, 8, 118};
        static const unsigned char mixed_bonds[] = {1, 255, 2, 255, 3};
        print_bonds_to_metal_record("mixed-metal-nonmetal", 5, mixed_elements, mixed_bonds);
    }
    {
        static const unsigned char metal_elements[] = {3, 26};
        static const unsigned char early_bonds[] = {4, 1};
        static const unsigned char late_bonds[] = {1, 4};
        print_bonds_to_metal_record("invalid-early", 2, metal_elements, early_bonds);
        print_bonds_to_metal_record("invalid-late", 2, metal_elements, late_bonds);
    }
    for (i = 0; i < 20; i++)
    {
        elements[i] = 3;
        bond_types[i] = 3;
    }
    print_bonds_to_metal_record("maxval-triples", 20, elements, bond_types);
    elements[0] = 121;
    bond_types[0] = 255;
    print_bonds_to_metal_record("sentinel-invalid-ignored", 1, elements, bond_types);
    return 0;
}

static void print_atom(const inchi_Atom *atom)
{
    printf("{\"x_bits\":%" PRIu64 ",\"y_bits\":%" PRIu64
           ",\"z_bits\":%" PRIu64 ",\"neighbor\":",
           double_bits(atom->x), double_bits(atom->y), double_bits(atom->z));
    print_i16_array(atom->neighbor, MAXVAL);
    fputs(",\"bond_type\":", stdout);
    print_i8_array(atom->bond_type, MAXVAL);
    fputs(",\"bond_stereo\":", stdout);
    print_i8_array(atom->bond_stereo, MAXVAL);
    fputs(",\"elname_bytes\":", stdout);
    print_char_array(atom->elname, ATOM_EL_LEN);
    printf(",\"num_bonds\":%d,\"num_iso_h\":", (int) atom->num_bonds);
    print_i8_array(atom->num_iso_H, NUM_H_ISOTOPES + 1);
    printf(",\"isotopic_mass\":%d,\"radical\":%d,\"charge\":%d}",
           (int) atom->isotopic_mass, (int) atom->radical, (int) atom->charge);
}

static void print_inp_atom_logical(const inp_ATOM *atom)
{
    int first = 1;
    int i;
#define EMIT_ATOM_INT(value) \
    do \
    { \
        if (!first) \
        { \
            putchar(','); \
        } \
        printf("%lld", (long long) (value)); \
        first = 0; \
    } while (0)
    fputs("{\"integer_fields\":[", stdout);
    for (i = 0; i < 6; i++) EMIT_ATOM_INT(atom->elname[i]);
    EMIT_ATOM_INT(atom->el_number);
    for (i = 0; i < MAXVAL; i++) EMIT_ATOM_INT(atom->neighbor[i]);
    EMIT_ATOM_INT(atom->orig_at_number);
    EMIT_ATOM_INT(atom->orig_compt_at_numb);
    for (i = 0; i < MAXVAL; i++) EMIT_ATOM_INT(atom->bond_stereo[i]);
    for (i = 0; i < MAXVAL; i++) EMIT_ATOM_INT(atom->bond_type[i]);
    EMIT_ATOM_INT(atom->valence);
    EMIT_ATOM_INT(atom->chem_bonds_valence);
    EMIT_ATOM_INT(atom->num_H);
    for (i = 0; i < 3; i++) EMIT_ATOM_INT(atom->num_iso_H[i]);
    EMIT_ATOM_INT(atom->iso_atw_diff);
    EMIT_ATOM_INT(atom->charge);
    EMIT_ATOM_INT(atom->radical);
    EMIT_ATOM_INT(atom->bAmbiguousStereo);
    EMIT_ATOM_INT(atom->cFlags);
    EMIT_ATOM_INT(atom->at_type);
    EMIT_ATOM_INT(atom->component);
    EMIT_ATOM_INT(atom->endpoint);
    EMIT_ATOM_INT(atom->c_point);
    EMIT_ATOM_INT(atom->bUsed0DParity);
    EMIT_ATOM_INT(atom->p_parity);
    for (i = 0; i < 4; i++) EMIT_ATOM_INT(atom->p_orig_at_num[i]);
    for (i = 0; i < 3; i++) EMIT_ATOM_INT(atom->sb_ord[i]);
    for (i = 0; i < 3; i++) EMIT_ATOM_INT(atom->sn_ord[i]);
    for (i = 0; i < 3; i++) EMIT_ATOM_INT(atom->sb_parity[i]);
    for (i = 0; i < 3; i++) EMIT_ATOM_INT(atom->sn_orig_at_num[i]);
    EMIT_ATOM_INT(atom->bCutVertex);
    EMIT_ATOM_INT(atom->nRingSystem);
    EMIT_ATOM_INT(atom->nNumAtInRingSystem);
    EMIT_ATOM_INT(atom->nBlockSystem);
    printf("],\"coordinate_bits\":[%" PRIu64 ",%" PRIu64 ",%" PRIu64 "]}",
           double_bits(atom->x),
           double_bits(atom->y),
           double_bits(atom->z));
#undef EMIT_ATOM_INT
}

static void print_set_atom_record(const char *case_id,
                                  int radical,
                                  int charge,
                                  const char *element,
                                  double x,
                                  double y,
                                  double z,
                                  int use_coordinates,
                                  int atom_index,
                                  int initial_dimensions,
                                  int initial_error,
                                  const char *initial_error_text)
{
    inp_ATOM targets[2];
    inchi_Atom inputs[2];
    MOL_COORD coordinates[2];
    char error_text[STR_ERR_LEN];
    int dimensions = initial_dimensions;
    int error = initial_error;
    int status;

    memset(targets, 0, sizeof(targets));
    memset(inputs, 0, sizeof(inputs));
    memset(coordinates, '#', sizeof(coordinates));
    memset(error_text, 0, sizeof(error_text));
    memset(targets[atom_index].elname, 'Z', sizeof(targets[atom_index].elname));
    strncpy(error_text, initial_error_text, sizeof(error_text) - 1);
    strncpy(inputs[atom_index].elname, element, sizeof(inputs[atom_index].elname) - 1);
    inputs[atom_index].radical = (S_CHAR) radical;
    inputs[atom_index].charge = (S_CHAR) charge;
    inputs[atom_index].x = x;
    inputs[atom_index].y = y;
    inputs[atom_index].z = z;
    status = SetAtomProperties(targets,
                               use_coordinates ? coordinates : NULL,
                               inputs,
                               atom_index,
                               &dimensions,
                               error_text,
                               &error);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
           "\"api_version\":\"1.07.5\"},\"case_id\":");
    print_json_string(case_id);
    printf(",\"operation\":\"SetAtomProperties\",\"input\":{\"radical\":%d,"
           "\"charge\":%d,\"element\":",
           radical,
           charge);
    print_json_string(element);
    printf(",\"coordinate_bits\":[%" PRIu64 ",%" PRIu64 ",%" PRIu64 "],"
           "\"use_coordinates\":%s,\"atom_index\":%d,"
           "\"initial_dimensions\":%d,\"initial_error\":%d,"
           "\"initial_error_text\":",
           double_bits(x),
           double_bits(y),
           double_bits(z),
           use_coordinates ? "true" : "false",
           atom_index,
           initial_dimensions,
           initial_error);
    print_json_string(initial_error_text);
    printf("},\"output\":{\"status\":%d,\"atom\":", status);
    print_inp_atom_logical(targets + atom_index);
    fputs(",\"coordinate_buffer\":", stdout);
    print_char_array(coordinates[atom_index], (int) sizeof(MOL_COORD));
    printf(",\"dimensions\":%d,\"error\":%d,\"error_text\":",
           dimensions,
           error);
    print_json_string(error_text);
    fputs("}}\n", stdout);
}

static int print_set_atom_records(void)
{
    int radical;
    char case_id[64];
    for (radical = -128; radical <= 127; radical++)
    {
        snprintf(case_id, sizeof(case_id), "radical-%d", radical);
        print_set_atom_record(case_id,
                              radical,
                              radical,
                              "C",
                              0.0,
                              -0.0,
                              0.0,
                              radical & 1,
                              0,
                              4,
                              4,
                              "");
    }
    print_set_atom_record("threshold-x-below", 0, -128, "H", nextafter(1.0e-6, 0.0), 0.0, 0.0, 1, 0, 4, 0, "");
    print_set_atom_record("threshold-x-equal", 1, 127, "He", 1.0e-6, 0.0, 0.0, 1, 0, 4, 0, "");
    print_set_atom_record("threshold-x-above", 2, 0, "Li", nextafter(1.0e-6, INFINITY), 0.0, 0.0, 1, 0, 4, 0, "");
    print_set_atom_record("threshold-y-negative", 3, 0, "Be", 0.0, -nextafter(1.0e-6, INFINITY), 0.0, 1, 0, 1, 0, "");
    print_set_atom_record("threshold-z-above", 4, 0, "B", 0.0, 0.0, nextafter(1.0e-6, INFINITY), 1, 0, 4, 0, "");
    print_set_atom_record("nan", 5, 0, "N", NAN, NAN, NAN, 1, 0, 0, 0, "");
    print_set_atom_record("infinity", 6, 0, "O", INFINITY, -INFINITY, 0.0, 1, 0, 0, 0, "");
    print_set_atom_record("index-one", 7, -7, "Clxyz", 1.25, -2.5, 3.0e-6, 1, 1, 4, 16, "seed");
    return 0;
}

static void setup_set_bond_atoms(inp_ATOM *atoms, int mode)
{
    int i;
    memset(atoms, 0, 2 * sizeof(*atoms));
    strcpy(atoms[0].elname, "C");
    strcpy(atoms[1].elname, "O");
    if (mode == 4 || mode == 5 || mode == 15)
    {
        atoms[0].valence = atoms[1].valence = 1;
        atoms[0].neighbor[0] = 1;
        atoms[1].neighbor[0] = 0;
        atoms[0].bond_type[0] = atoms[1].bond_type[0] = mode == 5 ? 2 : 1;
        atoms[0].bond_stereo[0] = mode == 15 ? 1 : 0;
        atoms[1].bond_stereo[0] = mode == 15 ? -1 : 0;
    }
    else if (mode == 6 || mode == 13)
    {
        atoms[0].valence = 2;
        atoms[0].neighbor[0] = atoms[0].neighbor[1] = 1;
        atoms[0].bond_type[0] = atoms[0].bond_type[1] = 1;
        atoms[1].valence = 1;
        atoms[1].neighbor[0] = 0;
        atoms[1].bond_type[0] = 1;
    }
    else if (mode == 7 || mode == 8 || mode == 12 || mode == 14)
    {
        atoms[0].valence = 1;
        atoms[0].neighbor[0] = 1;
        atoms[0].bond_type[0] = (mode == 8 || mode == 14) ? 2 : 1;
        if (mode == 12)
        {
            atoms[1].valence = MAXVAL;
            for (i = 0; i < MAXVAL; i++) atoms[1].neighbor[i] = 99;
        }
    }
    else if (mode == 9)
    {
        atoms[1].valence = 1;
        atoms[1].neighbor[0] = 0;
        atoms[1].bond_type[0] = 1;
    }
    else if (mode == 10 || mode == 11)
    {
        inp_ATOM *full = atoms + (mode == 10 ? 0 : 1);
        full->valence = MAXVAL;
        for (i = 0; i < MAXVAL; i++) full->neighbor[i] = 99;
    }
}

static void print_set_bond_record(const char *case_id,
                                  int mode,
                                  int bond_type,
                                  int bond_stereo)
{
    inp_ATOM atoms[2];
    inchi_Atom inputs[2];
    char error_text[STR_ERR_LEN];
    int num_bonds = 5;
    int error = 16;
    int status;
    int neighbor = mode == 1 ? -1 : mode == 2 ? 2 : mode == 3 ? 0 : 1;

    setup_set_bond_atoms(atoms, mode);
    memset(inputs, 0, sizeof(inputs));
    memset(error_text, 0, sizeof(error_text));
    inputs[0].neighbor[0] = (AT_NUM) neighbor;
    inputs[0].bond_type[0] = (S_CHAR) bond_type;
    inputs[0].bond_stereo[0] = (S_CHAR) bond_stereo;
    status = SetBondProperties(atoms,
                               inputs,
                               0,
                               0,
                               2,
                               &num_bonds,
                               error_text,
                               &error);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\"," 
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\"," 
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\"," 
           "\"api_version\":\"1.07.5\"},\"case_id\":");
    print_json_string(case_id);
    printf(",\"operation\":\"SetBondProperties\","
           "\"input\":{\"mode\":%d,\"bond_type\":%d,\"bond_stereo\":%d},"
           "\"output\":{\"status\":%d,\"atoms\":[",
           mode,
           bond_type,
           bond_stereo,
           status);
    print_inp_atom_logical(atoms);
    putchar(',');
    print_inp_atom_logical(atoms + 1);
    printf("],\"num_bonds\":%d,\"error\":%d,\"error_text\":",
           num_bonds,
           error);
    print_json_string(error_text);
    fputs("}}\n", stdout);
}

static int print_set_bond_records(void)
{
    char case_id[64];
    int value;
    int mode;
    for (value = -128; value <= 127; value++)
    {
        snprintf(case_id, sizeof(case_id), "bond-type-%d", value);
        print_set_bond_record(case_id, 0, value, 0);
    }
    for (value = -128; value <= 127; value++)
    {
        snprintf(case_id, sizeof(case_id), "bond-stereo-%d", value);
        print_set_bond_record(case_id, 0, 1, value);
    }
    for (mode = 1; mode <= 15; mode++)
    {
        snprintf(case_id, sizeof(case_id), "state-mode-%d", mode);
        print_set_bond_record(case_id, mode, 1, 0);
    }
    return 0;
}

static void print_atoms(const inchi_Atom *atoms, int atom_count)
{
    int i;
    putchar('[');
    for (i = 0; i < atom_count; i++)
    {
        if (i)
        {
            putchar(',');
        }
        print_atom(atoms + i);
    }
    putchar(']');
}

static void print_bond_fields(const inchi_Atom *atoms, int atom_count)
{
    int atom_index;
    int slot;
    int first = 1;
    putchar('[');
    for (atom_index = 0; atom_index < atom_count; atom_index++)
    {
        for (slot = 0; slot < atoms[atom_index].num_bonds; slot++)
        {
            if (!first)
            {
                putchar(',');
            }
            first = 0;
            printf("{\"atom_index\":%d,\"slot\":%d,\"neighbor\":%d,"
                   "\"bond_type\":%d,\"bond_stereo\":%d}",
                   atom_index, slot, (int) atoms[atom_index].neighbor[slot],
                   (int) atoms[atom_index].bond_type[slot],
                   (int) atoms[atom_index].bond_stereo[slot]);
        }
    }
    putchar(']');
}

static void print_stereo(const inchi_Stereo0D *stereo, int stereo_count)
{
    int i;
    putchar('[');
    for (i = 0; i < stereo_count; i++)
    {
        if (i)
        {
            putchar(',');
        }
        printf("{\"central_atom\":%d,\"neighbor\":",
               (int) stereo[i].central_atom);
        print_i16_array(stereo[i].neighbor, 4);
        printf(",\"type\":%d,\"parity\":%d}",
               (int) stereo[i].type, (int) stereo[i].parity);
    }
    putchar(']');
}

static const char *pointer_state(const void *pointer, const void *original);

typedef struct tagInchiToInputOracleCase
{
    const char *case_id;
    const char *text;
    INPUT_TYPE input_type;
    int file_mode;
    int merge_all;
    int do_not_add_h;
    int ab_parity_unknown;
    int initial_atom_count;
    int initial_atom_storage;
    int initial_stereo_count;
    int initial_stereo_storage;
    int allocation_failure_ordinal;
    int initial_error;
    int omit_output;
    int omit_optional_outputs;
} INCHI_TO_INPUT_ORACLE_CASE;

static int oracle_pointer_was_freed(const void *pointer)
{
    size_t i;
    for (i = 0; i < ORACLE_DEFERRED_FREE_COUNT; i++)
    {
        if (ORACLE_DEFERRED_FREES[i] == pointer)
        {
            return 1;
        }
    }
    return 0;
}

static const INCHI_TO_INPUT_ORACLE_CASE INCHI_TO_INPUT_CASES[] = {
    {
        .case_id = "single-labeled-string",
        .text = "Structure: 42. LABEL =value\n"
                "AuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:0,0,0;1,1,1;\n",
        .input_type = INPUT_INCHI_PLAIN,
    },
    {
        .case_id = "single-labeled-file",
        .text = "Structure: 42. LABEL =value\n"
                "AuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:0,0,0;1,1,1;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .file_mode = 1,
    },
    {
        .case_id = "single-stops-before-second",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n"
                "AuxInfo=1/rA:1N/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
    },
    {
        .case_id = "merge-three-eof",
        .text = "Structure: 1. FIRST =one\n"
                "AuxInfo=1/rA:1C/rB:/rC:;\n"
                "Structure: 2. SECOND =two\n"
                "AuxInfo=1/rA:2NO/rB:s1;/rC:;;\n"
                "Structure: 3. THIRD =three\n"
                "AuxInfo=1/rA:1Coe/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .merge_all = 1,
    },
    {
        .case_id = "merge-preexisting-remap",
        .text = "AuxInfo=1/rA:2NO0o/rB:s1;/rC:0,0,0;1,0,0;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .merge_all = 1,
        .initial_atom_count = 1,
        .initial_atom_storage = 1,
        .initial_stereo_count = 1,
        .initial_stereo_storage = 1,
    },
    {
        .case_id = "count-only",
        .text = "AuxInfo=1/rA:2CO/rB:s1;/rC:;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_output = 1,
    },
    {
        .case_id = "empty-structure",
        .text = "AuxInfo=1//\n",
        .input_type = INPUT_INCHI_PLAIN,
    },
    {
        .case_id = "malformed-cleans-old",
        .text = "AuxInfo=1/rA:2CO/rB:x/rC:;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .initial_atom_count = 1,
        .initial_atom_storage = 1,
        .initial_stereo_count = 1,
        .initial_stereo_storage = 1,
    },
    {
        .case_id = "max-atoms",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .initial_atom_count = MAX_ATOMS - 1,
        .initial_atom_storage = 1,
    },
    {
        .case_id = "seeded-error-unknown-text",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .initial_error = 77,
    },
    {
        .case_id = "xml-single",
        .text = "<structure number=\"7\" id.name=\"xml-label\" id.value=\"xml-value\">\n"
                "<identifier.auxiliary-info>\n<reversibility>\n"
                "<atoms>\n2nCO\n</atoms>\n"
                "<bonds>\ns1;\n</bonds>\n"
                "<xyz>\n0,0,0;1,1,1;\n</xyz>\n",
        .input_type = INPUT_INCHI_XML,
    },
    {
        .case_id = "standalone-eof",
        .text = "",
        .input_type = INPUT_INCHI_PLAIN,
    },
    {
        .case_id = "warning-text-success",
        .text = "AuxInfo=1/rA:2CO/rB:a1;/rC:;;\n",
        .input_type = INPUT_INCHI_PLAIN,
    },
    {
        .case_id = "xml-fatal-cleanup",
        .text = "<structure number=\"9\">\n"
                "<message type=\"fatal (aborted)\" value=\"source fatal\">\n",
        .input_type = INPUT_INCHI_XML,
        .initial_atom_count = 1,
        .initial_atom_storage = 1,
        .initial_stereo_count = 1,
        .initial_stereo_storage = 1,
    },
    {
        .case_id = "merge-double-bond-stereo",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n"
                "AuxInfo=1/rA:2CC/rB:d--1;/rC:;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .merge_all = 1,
    },
    {
        .case_id = "non-null-zero-count-old",
        .text = "AuxInfo=1/rA:1C/rB:/rC:;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .initial_atom_storage = 1,
    },
    {
        .case_id = "null-optional-eof",
        .text = "",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_optional_outputs = 1,
    },
    {
        .case_id = "null-optional-malformed",
        .text = "AuxInfo=1/rA:2CO/rB:x/rC:;;\n",
        .input_type = INPUT_INCHI_PLAIN,
        .omit_optional_outputs = 1,
    },
};

static int run_inchi_to_input_case(const INCHI_TO_INPUT_ORACLE_CASE *test_case)
{
    INCHI_IOSTREAM input;
    inchi_Input output;
    inchi_Input *output_pointer = test_case->omit_output ? NULL : &output;
    inchi_Atom *original_atoms = NULL;
    inchi_Stereo0D *original_stereo = NULL;
    char *original_options = NULL;
    char *input_copy = NULL;
    FILE *input_file = NULL;
    char sdf_label[MAX_SDF_HEADER] = { 0 };
    char sdf_value[MAX_SDF_VALUE] = { 0 };
    char error_text[STR_ERR_LEN] = { 0 };
    long sdf_id = -1;
    INCHI_MODE atom_flags = 0;
    int error = test_case->initial_error;
    int status;
    int allocation_calls;
    int original_atom_freed;
    int original_stereo_freed;
    long input_position;
    size_t input_length = strlen(test_case->text);

    memset(&input, 0, sizeof(input));
    memset(&output, 0, sizeof(output));
    if (test_case->file_mode)
    {
        input_file = tmpfile();
        if (!input_file ||
            fwrite(test_case->text, 1, input_length, input_file) != input_length)
        {
            fputs("failed to create InchiToInchi_Input oracle file\n", stderr);
            return 70;
        }
        rewind(input_file);
        input.f = input_file;
        input.type = INCHI_IOS_TYPE_FILE;
    }
    else
    {
        input_copy = (char *) malloc(input_length + 1);
        if (!input_copy)
        {
            fputs("failed to allocate InchiToInchi_Input oracle input\n", stderr);
            return 71;
        }
        memcpy(input_copy, test_case->text, input_length + 1);
        input.s.pStr = input_copy;
        input.s.nAllocatedLength = (int) input_length + 1;
        input.s.nUsedLength = (int) input_length;
        input.s.nPtr = 0;
        input.type = INCHI_IOS_TYPE_STRING;
    }

    if (output_pointer)
    {
        if (test_case->initial_atom_storage > 0)
        {
            original_atoms = (inchi_Atom *) calloc(
                (size_t) test_case->initial_atom_storage, sizeof(*original_atoms));
            if (!original_atoms)
            {
                return 72;
            }
            original_atoms[0].x = -0.0;
            original_atoms[0].y = 2.5;
            original_atoms[0].z = -3.5;
            memcpy(original_atoms[0].elname, "F", 2);
            original_atoms[0].num_iso_H[0] = 4;
            original_atoms[0].charge = -1;
            output.atom = original_atoms;
        }
        if (test_case->initial_stereo_storage > 0)
        {
            original_stereo = (inchi_Stereo0D *) calloc(
                (size_t) test_case->initial_stereo_storage, sizeof(*original_stereo));
            if (!original_stereo)
            {
                return 73;
            }
            original_stereo[0].central_atom = 0;
            original_stereo[0].neighbor[0] = 0;
            original_stereo[0].neighbor[1] = 0;
            original_stereo[0].neighbor[2] = 0;
            original_stereo[0].neighbor[3] = 0;
            original_stereo[0].type = INCHI_StereoType_Tetrahedral;
            original_stereo[0].parity = INCHI_PARITY_ODD;
            output.stereo0D = original_stereo;
        }
        original_options = (char *) malloc(4);
        if (!original_options)
        {
            return 74;
        }
        memcpy(original_options, "-X\0", 3);
        original_options[3] = (char) 0x5a;
        output.szOptions = original_options;
        output.num_atoms = (AT_NUM) test_case->initial_atom_count;
        output.num_stereo0D = (AT_NUM) test_case->initial_stereo_count;
    }

    ORACLE_ALLOCATION_CALLS = 0;
    ORACLE_ALLOCATION_ORDINAL = test_case->allocation_failure_ordinal;
    ORACLE_ALLOCATION_FAILURE_ENABLED = 1;
    ORACLE_DEFER_FREES = 1;
    status = InchiToInchi_Input(&input, output_pointer, test_case->merge_all,
                                test_case->do_not_add_h,
                                test_case->ab_parity_unknown,
                                test_case->input_type,
                                test_case->omit_optional_outputs ? NULL : sdf_label,
                                test_case->omit_optional_outputs ? NULL : sdf_value,
                                test_case->omit_optional_outputs ? NULL : &sdf_id,
                                test_case->omit_optional_outputs ? NULL : &atom_flags,
                                &error,
                                test_case->omit_optional_outputs ? NULL : error_text);
    ORACLE_ALLOCATION_FAILURE_ENABLED = 0;
    allocation_calls = ORACLE_ALLOCATION_CALLS;
    input_position = test_case->file_mode ? ftell(input_file) : input.s.nPtr;
    original_atom_freed = oracle_pointer_was_freed(original_atoms);
    original_stereo_freed = oracle_pointer_was_freed(original_stereo);

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\","
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\","
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\","
          "\"api_version\":\"1.07.5\"},\"case_id\":", stdout);
    print_json_string(test_case->case_id);
    fputs(",\"operation\":\"inchi_to_inchi_input\",\"input\":{\"mode\":", stdout);
    print_json_string(test_case->file_mode ? "file" : "string");
    fputs(",\"text\":", stdout);
    print_json_string(test_case->text);
    printf(",\"input_type\":%d,\"merge_all\":%d,\"do_not_add_h\":%d,"
           "\"ab_parity_unknown\":%d,\"initial_atom_count\":%d,"
           "\"initial_atom_storage\":%d,\"initial_stereo_count\":%d,"
           "\"initial_stereo_storage\":%d,\"allocation_failure_ordinal\":%d,"
           "\"initial_error\":%d,\"omit_output\":%d,"
           "\"omit_optional_outputs\":%d},\"output\":{"
           "\"status\":%d,\"error_code\":%d,\"error_text\":",
           (int) test_case->input_type, test_case->merge_all,
           test_case->do_not_add_h, test_case->ab_parity_unknown,
           test_case->initial_atom_count, test_case->initial_atom_storage,
           test_case->initial_stereo_count, test_case->initial_stereo_storage,
           test_case->allocation_failure_ordinal, test_case->initial_error,
           test_case->omit_output, test_case->omit_optional_outputs, status, error);
    if (test_case->omit_optional_outputs) fputs("null", stdout); else print_json_string(error_text);
    printf(",\"input_position\":%ld,\"sdf_label\":", input_position);
    if (test_case->omit_optional_outputs) fputs("null", stdout); else print_json_string(sdf_label);
    fputs(",\"sdf_value\":", stdout);
    if (test_case->omit_optional_outputs) fputs("null", stdout); else print_json_string(sdf_value);
    if (test_case->omit_optional_outputs)
    {
        printf(",\"sdf_id\":null,\"atom_flags\":null,\"allocation_calls\":%d,",
               allocation_calls);
    }
    else
    {
        printf(",\"sdf_id\":%ld,\"atom_flags\":%llu,\"allocation_calls\":%d,",
               sdf_id, (unsigned long long) atom_flags, allocation_calls);
    }
    if (!output_pointer)
    {
        fputs("\"num_atoms\":null,\"num_stereo0d\":null,\"atoms\":[],"
              "\"bond_fields\":[],\"stereo0d\":[],\"atom_pointer_state\":null,"
              "\"stereo_pointer_state\":null,\"options_pointer_state\":null,",
              stdout);
    }
    else
    {
        printf("\"num_atoms\":%d,\"num_stereo0d\":%d,\"atoms\":",
               (int) output.num_atoms, (int) output.num_stereo0D);
        print_atoms(output.atom, output.atom && output.num_atoms > 0 ? output.num_atoms : 0);
        fputs(",\"bond_fields\":", stdout);
        print_bond_fields(output.atom, output.atom && output.num_atoms > 0 ? output.num_atoms : 0);
        fputs(",\"stereo0d\":", stdout);
        print_stereo(output.stereo0D,
                     output.stereo0D && output.num_stereo0D > 0 ? output.num_stereo0D : 0);
        fputs(",\"atom_pointer_state\":", stdout);
        print_json_string(pointer_state(output.atom, original_atoms));
        fputs(",\"stereo_pointer_state\":", stdout);
        print_json_string(pointer_state(output.stereo0D, original_stereo));
        fputs(",\"options_pointer_state\":", stdout);
        print_json_string(pointer_state(output.szOptions, original_options));
        putchar(',');
    }
    printf("\"original_atom_freed\":%s,\"original_stereo_freed\":%s}}\n",
           original_atom_freed ? "true" : "false",
           original_stereo_freed ? "true" : "false");

    if (output_pointer)
    {
        FreeInchi_Input(output_pointer);
    }
    ORACLE_DEFER_FREES = 0;
    oracle_flush_deferred_frees();
    free(original_options);
    free(input_copy);
    if (input_file)
    {
        fclose(input_file);
    }
    return 0;
}

static int print_inchi_to_input_records(void)
{
    size_t i;
    int ordinal;
    for (i = 0; i < sizeof(INCHI_TO_INPUT_CASES) / sizeof(INCHI_TO_INPUT_CASES[0]); i++)
    {
        int result = run_inchi_to_input_case(INCHI_TO_INPUT_CASES + i);
        if (result)
        {
            return result;
        }
    }
    for (ordinal = 1; ordinal <= 16; ordinal++)
    {
        char case_id[64];
        INCHI_TO_INPUT_ORACLE_CASE test_case;
        memset(&test_case, 0, sizeof(test_case));
        snprintf(case_id, sizeof(case_id), "allocation-atom-%d", ordinal);
        test_case.case_id = case_id;
        test_case.text = "AuxInfo=1/rA:1C/rB:/rC:;\n"
                         "AuxInfo=1/rA:1N/rB:/rC:;\n";
        test_case.input_type = INPUT_INCHI_PLAIN;
        test_case.merge_all = 1;
        test_case.allocation_failure_ordinal = ordinal;
        if (run_inchi_to_input_case(&test_case))
        {
            return 75;
        }
    }
    for (ordinal = 1; ordinal <= 16; ordinal++)
    {
        char case_id[64];
        INCHI_TO_INPUT_ORACLE_CASE test_case;
        memset(&test_case, 0, sizeof(test_case));
        snprintf(case_id, sizeof(case_id), "allocation-stereo-%d", ordinal);
        test_case.case_id = case_id;
        test_case.text = "AuxInfo=1/rA:1Coe/rB:/rC:;\n"
                         "AuxInfo=1/rA:1C0o/rB:/rC:;\n";
        test_case.input_type = INPUT_INCHI_PLAIN;
        test_case.merge_all = 1;
        test_case.allocation_failure_ordinal = ordinal;
        if (run_inchi_to_input_case(&test_case))
        {
            return 76;
        }
    }
    return 0;
}

typedef struct tagAuxInputOracleCase
{
    const char *case_id;
    const char *text;
    int generated_text;
    int do_not_add_h;
    int distinguish_unknown_undefined_stereo;
    int allocation_failure_ordinal;
    int output_mode;
    int seed_old_allocations;
} AUX_INPUT_ORACLE_CASE;

static const AUX_INPUT_ORACLE_CASE AUX_INPUT_CASES[] = {
    { "null-outer", "AuxInfo=1/rA:1C/rB:/rC:;\n", 0, 0, 0, 0, 1, 0 },
    { "null-inner", "AuxInfo=1/rA:1C/rB:/rC:;\n", 0, 0, 0, 0, 2, 0 },
    { "single-default-h", "AuxInfo=1/rA:1C/rB:/rC:;\n", 0, 0, 0, 0, 0, 0 },
    { "single-no-add-h", "AuxInfo=1/rA:1C/rB:/rC:;\n", 0, 1, 0, 0, 0, 0 },
    {
        "single-chiral",
        "Structure: 42. LABEL =value\n"
        "AuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:0,0,0;1,1,1;\n",
        0, 0, 0, 0, 0, 0
    },
    {
        "merge-remap-chiral",
        "Structure: 1. FIRST =one\nAuxInfo=1/0/N:1/rA:1C/rB:/rC:;\n"
        "Structure: 2. SECOND =two\nAuxInfo=1/rA:2NO0o/rB:s1;/rC:;;\n",
        0, 0, 0, 0, 0, 0
    },
    {
        "tetrahedral-unknown-default",
        "AuxInfo=1/rA:4CCCClu/rB:;;s1s2s3;/rC:;;;;\n",
        0, 0, 0, 0, 0, 0
    },
    {
        "tetrahedral-unknown-distinct",
        "AuxInfo=1/rA:4CCCClu/rB:;;s1s2s3;/rC:;;;;\n",
        0, 0, 1, 0, 0, 0
    },
    { "empty-structure", "AuxInfo=1//\n", 0, 0, 0, 0, 0, 0 },
    { "eof", "", 0, 0, 0, 0, 0, 0 },
    { "missing-atom", "AuxInfo=1/no-reversibility\n", 0, 0, 0, 0, 0, 0 },
    { "wrong-bonds", "AuxInfo=1/rA:2CO/rB:x/rC:;;\n", 0, 0, 0, 0, 0, 0 },
    {
        "aromatic-error",
        "AuxInfo=1/rA:2CO/rB:a1;/rC:;;\n",
        0, 0, 0, 0, 0, 0
    },
    { "old-reset-success", "AuxInfo=1/rA:1N/rB:/rC:;\n", 0, 0, 0, 0, 0, 1 },
    {
        "old-reset-error",
        "AuxInfo=1/rA:2CO/rB:x/rC:;;\n",
        0, 0, 0, 0, 0, 1
    },
    { "max-atoms", NULL, 1, 0, 0, 0, 0, 0 },
};

static char *generated_aux_input_text(int kind)
{
    const size_t atom_count = MAX_ATOMS;
    const char prefix[] = "AuxInfo=1/rA:32766";
    const char bonds[] = "/rB:";
    const char coordinates[] = "/rC:";
    size_t length;
    char *text;
    char *cursor;
    if (kind != 1)
    {
        return NULL;
    }
    length = sizeof(prefix) - 1 + atom_count + sizeof(bonds) - 1 + atom_count - 1 +
             sizeof(coordinates) - 1 + atom_count + 1;
    text = (char *) __real_malloc(length + 1);
    if (!text)
    {
        return NULL;
    }
    cursor = text;
    memcpy(cursor, prefix, sizeof(prefix) - 1);
    cursor += sizeof(prefix) - 1;
    memset(cursor, 'C', atom_count);
    cursor += atom_count;
    memcpy(cursor, bonds, sizeof(bonds) - 1);
    cursor += sizeof(bonds) - 1;
    memset(cursor, ';', atom_count - 1);
    cursor += atom_count - 1;
    memcpy(cursor, coordinates, sizeof(coordinates) - 1);
    cursor += sizeof(coordinates) - 1;
    memset(cursor, ';', atom_count);
    cursor += atom_count;
    *cursor++ = '\n';
    *cursor = '\0';
    return text;
}

static int run_aux_input_case(const AUX_INPUT_ORACLE_CASE *test_case, int standard)
{
    InchiInpData output;
    InchiInpData *output_pointer = test_case->output_mode == 1 ? NULL : &output;
    inchi_Input input;
    inchi_Atom original_atom_value;
    inchi_Stereo0D original_stereo_value;
    inchi_Atom *original_atoms = NULL;
    inchi_Stereo0D *original_stereo = NULL;
    char *original_options = NULL;
    char *generated_text = NULL;
    char *input_copy;
    const char *case_text;
    size_t input_length;
    int status;
    int allocation_calls;
    int original_atom_freed;
    int original_stereo_freed;
    int original_atom_unchanged;
    int original_stereo_unchanged;

    generated_text = test_case->generated_text ?
        generated_aux_input_text(test_case->generated_text) : NULL;
    case_text = generated_text ? generated_text : test_case->text;
    if (!case_text)
    {
        return 77;
    }
    input_length = strlen(case_text);
    input_copy = (char *) __real_malloc(input_length + 1);
    if (!input_copy)
    {
        __real_free(generated_text);
        return 78;
    }
    memcpy(input_copy, case_text, input_length + 1);

    memset(&output, 0x59, sizeof(output));
    memset(&input, 0, sizeof(input));
    memset(&original_atom_value, 0, sizeof(original_atom_value));
    memset(&original_stereo_value, 0, sizeof(original_stereo_value));
    if (test_case->output_mode == 2)
    {
        output.pInp = NULL;
    }
    else
    {
        output.pInp = &input;
        original_options = (char *) __real_malloc(4);
        if (!original_options)
        {
            __real_free(input_copy);
            __real_free(generated_text);
            return 79;
        }
        memcpy(original_options, "-X\0Z", 4);
        input.szOptions = original_options;
        if (test_case->seed_old_allocations)
        {
            original_atoms = (inchi_Atom *) __real_calloc(1, sizeof(*original_atoms));
            original_stereo = (inchi_Stereo0D *) __real_calloc(1, sizeof(*original_stereo));
            if (!original_atoms || !original_stereo)
            {
                __real_free(original_atoms);
                __real_free(original_stereo);
                __real_free(original_options);
                __real_free(input_copy);
                __real_free(generated_text);
                return 80;
            }
            original_atoms[0].x = -0.0;
            original_atoms[0].y = 2.5;
            original_atoms[0].z = -3.5;
            memcpy(original_atoms[0].elname, "F", 2);
            original_atoms[0].num_iso_H[0] = 4;
            original_atoms[0].charge = -1;
            original_stereo[0].central_atom = 0;
            original_stereo[0].neighbor[0] = 0;
            original_stereo[0].neighbor[1] = 0;
            original_stereo[0].neighbor[2] = 0;
            original_stereo[0].neighbor[3] = 0;
            original_stereo[0].type = INCHI_StereoType_Tetrahedral;
            original_stereo[0].parity = INCHI_PARITY_ODD;
            original_atom_value = original_atoms[0];
            original_stereo_value = original_stereo[0];
            input.atom = original_atoms;
            input.stereo0D = original_stereo;
            input.num_atoms = 1;
            input.num_stereo0D = 1;
        }
    }

    ORACLE_ALLOCATION_CALLS = 0;
    ORACLE_ALLOCATION_ORDINAL = test_case->allocation_failure_ordinal;
    ORACLE_ALLOCATION_FAILURE_ENABLED = 1;
    ORACLE_DEFER_FREES = 1;
    status = standard ?
        Get_std_inchi_Input_FromAuxInfo(input_copy, test_case->do_not_add_h, output_pointer) :
        Get_inchi_Input_FromAuxInfo(input_copy, test_case->do_not_add_h,
                                   test_case->distinguish_unknown_undefined_stereo,
                                   output_pointer);
    ORACLE_ALLOCATION_FAILURE_ENABLED = 0;
    allocation_calls = ORACLE_ALLOCATION_CALLS;
    original_atom_freed = oracle_pointer_was_freed(original_atoms);
    original_stereo_freed = oracle_pointer_was_freed(original_stereo);
    original_atom_unchanged = !original_atoms ||
        memcmp(original_atoms, &original_atom_value, sizeof(original_atom_value)) == 0;
    original_stereo_unchanged = !original_stereo ||
        memcmp(original_stereo, &original_stereo_value, sizeof(original_stereo_value)) == 0;

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\","
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\","
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\","
          "\"api_version\":\"1.07.5\"},\"case_id\":", stdout);
    print_json_string(test_case->case_id);
    fputs(",\"operation\":", stdout);
    print_json_string(standard ? "get_std_inchi_input_from_aux_info" :
                                 "get_inchi_input_from_aux_info");
    fputs(",\"input\":{\"text\":", stdout);
    if (test_case->generated_text) fputs("null", stdout); else print_json_string(case_text);
    printf(",\"generated_text\":%d,\"do_not_add_h\":%d,"
           "\"distinguish_unknown_undefined_stereo\":%d,"
           "\"allocation_failure_ordinal\":%d,\"output_mode\":%d,"
           "\"seed_old_allocations\":%d},\"output\":{\"status\":%d,"
           "\"allocation_calls\":%d,\"input_bytes\":",
           test_case->generated_text, test_case->do_not_add_h,
           test_case->distinguish_unknown_undefined_stereo,
           test_case->allocation_failure_ordinal, test_case->output_mode,
           test_case->seed_old_allocations, status, allocation_calls);
    print_u8_array((const unsigned char *) input_copy, (int) input_length + 1);
    if (!output_pointer)
    {
        fputs(",\"input_pointer_state\":null,\"b_chiral\":null,"
              "\"error_text\":null,\"error_text_terminated\":null,"
              "\"error_bytes\":null,\"num_atoms\":null,"
              "\"num_stereo0d\":null,\"atoms\":[],\"bond_fields\":[],"
              "\"stereo0d\":[],\"atom_pointer_state\":null,"
              "\"stereo_pointer_state\":null,\"options_pointer_state\":null,",
              stdout);
    }
    else
    {
        fputs(",\"input_pointer_state\":", stdout);
        print_json_string(output.pInp == &input ? "reused" :
                          output.pInp ? "other" : "null");
        printf(",\"b_chiral\":%d,\"error_text\":", output.bChiral);
        if (memchr(output.szErrMsg, '\0', STR_ERR_LEN))
        {
            print_json_string(output.szErrMsg);
        }
        else
        {
            fputs("null", stdout);
        }
        printf(",\"error_text_terminated\":%s",
               memchr(output.szErrMsg, '\0', STR_ERR_LEN) ? "true" : "false");
        fputs(",\"error_bytes\":", stdout);
        print_u8_array((const unsigned char *) output.szErrMsg, STR_ERR_LEN);
        if (!output.pInp)
        {
            fputs(",\"num_atoms\":null,\"num_stereo0d\":null,\"atoms\":[],"
                  "\"bond_fields\":[],\"stereo0d\":[],"
                  "\"atom_pointer_state\":null,\"stereo_pointer_state\":null,"
                  "\"options_pointer_state\":null,", stdout);
        }
        else
        {
            printf(",\"num_atoms\":%d,\"num_stereo0d\":%d,\"atoms\":",
                   (int) input.num_atoms, (int) input.num_stereo0D);
            print_atoms(input.atom, input.atom && input.num_atoms > 0 ? input.num_atoms : 0);
            fputs(",\"bond_fields\":", stdout);
            print_bond_fields(input.atom, input.atom && input.num_atoms > 0 ? input.num_atoms : 0);
            fputs(",\"stereo0d\":", stdout);
            print_stereo(input.stereo0D,
                         input.stereo0D && input.num_stereo0D > 0 ? input.num_stereo0D : 0);
            fputs(",\"atom_pointer_state\":", stdout);
            print_json_string(pointer_state(input.atom, original_atoms));
            fputs(",\"stereo_pointer_state\":", stdout);
            print_json_string(pointer_state(input.stereo0D, original_stereo));
            fputs(",\"options_pointer_state\":", stdout);
            print_json_string(pointer_state(input.szOptions, original_options));
            putchar(',');
        }
    }
    printf("\"original_atom_freed\":%s,\"original_stereo_freed\":%s,"
           "\"original_atom_unchanged\":%s,\"original_stereo_unchanged\":%s}}\n",
           original_atom_freed ? "true" : "false",
           original_stereo_freed ? "true" : "false",
           original_atom_unchanged ? "true" : "false",
           original_stereo_unchanged ? "true" : "false");

    if (output_pointer && output.pInp)
    {
        Free_inchi_Input(output.pInp);
    }
    ORACLE_DEFER_FREES = 0;
    oracle_flush_deferred_frees();
    __real_free(original_atoms);
    __real_free(original_stereo);
    __real_free(original_options);
    __real_free(input_copy);
    __real_free(generated_text);
    return 0;
}

static int print_aux_input_records(int standard)
{
    size_t i;
    int ordinal;
    for (i = 0; i < sizeof(AUX_INPUT_CASES) / sizeof(AUX_INPUT_CASES[0]); i++)
    {
        AUX_INPUT_ORACLE_CASE test_case = AUX_INPUT_CASES[i];
        if (standard && test_case.distinguish_unknown_undefined_stereo)
        {
            continue;
        }
        if (run_aux_input_case(&test_case, standard))
        {
            return 81;
        }
    }
    for (ordinal = 1; ordinal <= 12; ordinal++)
    {
        char case_id[64];
        AUX_INPUT_ORACLE_CASE test_case;
        memset(&test_case, 0, sizeof(test_case));
        snprintf(case_id, sizeof(case_id), "allocation-merge-%d", ordinal);
        test_case.case_id = case_id;
        test_case.text = "AuxInfo=1/rA:1C/rB:/rC:;\n"
                         "AuxInfo=1/rA:2NO0o/rB:s1;/rC:;;\n";
        test_case.allocation_failure_ordinal = ordinal;
        if (run_aux_input_case(&test_case, standard))
        {
            return 82;
        }
    }
    return 0;
}

static char *generated_inchi_to_atom_text(int kind)
{
    const char prefix[] = "AuxInfo=";
    const char middle[] = "/rA:1C/rB:/rC:;";
    size_t first_run = (size_t) INCHI_LINE_LEN + 32;
    size_t second_run = kind == 1 ? (size_t) INCHI_LINE_LEN : 0;
    size_t length = sizeof(prefix) - 1 + first_run + sizeof(middle) - 1 +
                    (second_run ? 1 : 0) + second_run + 1;
    char *text;
    char *cursor;
    if (kind == 3)
    {
        const char xml_prefix[] =
            "<structure number=\"92\">\n<identifier.auxiliary-info>\n"
            "<reversibility>\n<atoms>\n262200";
        const char atoms_end[] = "\n</atoms>\n<bonds>\n";
        const char bonds_end[] = "\n</bonds>\n<xyz>\n";
        const char xml_end[] = "\n</xyz>\n";
        const size_t atom_count = 262200;
        length = sizeof(xml_prefix) - 1 + atom_count + sizeof(atoms_end) - 1 +
                 atom_count - 1 + sizeof(bonds_end) - 1 + atom_count +
                 sizeof(xml_end) - 1;
        text = (char *) malloc(length + 1);
        if (!text)
        {
            return NULL;
        }
        cursor = text;
        memcpy(cursor, xml_prefix, sizeof(xml_prefix) - 1);
        cursor += sizeof(xml_prefix) - 1;
        memset(cursor, 'C', atom_count);
        cursor += atom_count;
        memcpy(cursor, atoms_end, sizeof(atoms_end) - 1);
        cursor += sizeof(atoms_end) - 1;
        memset(cursor, ';', atom_count - 1);
        cursor += atom_count - 1;
        memcpy(cursor, bonds_end, sizeof(bonds_end) - 1);
        cursor += sizeof(bonds_end) - 1;
        memset(cursor, ';', atom_count);
        cursor += atom_count;
        memcpy(cursor, xml_end, sizeof(xml_end) - 1);
        cursor += sizeof(xml_end) - 1;
        *cursor = '\0';
        return text;
    }
    if (kind >= 4 && kind <= 9)
    {
        static const char *const malformed_prefixes[] = {
            "<structure number=\"96\">\n<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C",
            "<structure number=\"97\">\n<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n2CC\n</atoms>\n<bonds>\ns1;",
            "<structure number=\"98\">\n<identifier.auxiliary-info>\n<reversibility>\n<atoms>\n1C\n</atoms>\n<bonds>\n\n</bonds>\n<xyz>\n;",
            "AuxInfo=1/rA:1C",
            "AuxInfo=1/rA:2CC/rB:s1;",
            "AuxInfo=1/rA:1C/rB:/rC:;",
        };
        static const char *const malformed_suffixes[] = {
            "",
            "",
            "",
            "/rB:/rC:;\n",
            "/rC:;;\n",
            "\n",
        };
        const size_t index = (size_t) (kind - 4);
        const size_t prefix_length = strlen(malformed_prefixes[index]);
        const size_t suffix_length = strlen(malformed_suffixes[index]);
        const size_t fill_length = (size_t) INCHI_LINE_LEN;
        const int atom_fill = kind == 4 || kind == 7;
        length = prefix_length + fill_length + suffix_length;
        text = (char *) malloc(length + 1);
        if (!text)
        {
            return NULL;
        }
        cursor = text;
        memcpy(cursor, malformed_prefixes[index], prefix_length);
        cursor += prefix_length;
        memset(cursor, atom_fill ? 'C' : 'X', fill_length);
        cursor += fill_length;
        memcpy(cursor, malformed_suffixes[index], suffix_length);
        cursor += suffix_length;
        *cursor = '\0';
        return text;
    }
    text = (char *) malloc(length + 1);
    if (!text)
    {
        return NULL;
    }
    cursor = text;
    memcpy(cursor, prefix, sizeof(prefix) - 1);
    cursor += sizeof(prefix) - 1;
    memset(cursor, 'X', first_run);
    cursor += first_run;
    memcpy(cursor, middle, sizeof(middle) - 1);
    cursor += sizeof(middle) - 1;
    if (second_run)
    {
        *cursor++ = '/';
        memset(cursor, 'X', second_run);
        cursor += second_run;
    }
    *cursor++ = '\n';
    *cursor = '\0';
    return text;
}

static const char *pointer_state(const void *pointer, const void *original)
{
    if (!pointer)
    {
        return "null";
    }
    return original && pointer == original ? "reused" : "allocated";
}

static int run_inchi_to_atom_case(const ORACLE_CASE *test_case)
{
    INCHI_IOSTREAM input;
    inchi_Atom *atoms = NULL;
    inchi_Stereo0D *stereo = NULL;
    char sdf_label[MAX_SDF_HEADER] = { 0 };
    char sdf_value[MAX_SDF_VALUE] = { 0 };
    char error_text[STR_ERR_LEN] = { 0 };
    char *input_copy = NULL;
    FILE *input_file = NULL;
    int stereo_count = 0;
    int dimensions = -1;
    int bond_count = -1;
    int error_code = 0;
    int status;
    int atom_count;
    int atom_storage_count;
    int stereo_storage_count;
    int passed_atom_capacity;
    long id = -1;
    INCHI_MODE atom_flags = 0;
    char *generated_text = NULL;
    const char *case_text = test_case->text;
    inchi_Atom *original_atoms = NULL;
    inchi_Stereo0D *original_stereo = NULL;
    long input_position;
    size_t input_length;

    if (test_case->generated_text)
    {
        generated_text = generated_inchi_to_atom_text(test_case->generated_text);
        if (!generated_text)
        {
            fprintf(stderr, "failed to allocate generated oracle input\n");
            return 69;
        }
        case_text = generated_text;
    }
    input_length = strlen(case_text);

    memset(&input, 0, sizeof(input));
    if (test_case->file_mode)
    {
        input_file = tmpfile();
        if (!input_file || fwrite(case_text, 1, input_length, input_file) != input_length)
        {
            fprintf(stderr, "failed to create oracle input file\n");
            return 70;
        }
        rewind(input_file);
        input.f = input_file;
        input.type = INCHI_IOS_TYPE_FILE;
    }
    else
    {
        input_copy = (char *) malloc(input_length + 1);
        if (!input_copy)
        {
            fprintf(stderr, "failed to allocate oracle input string\n");
            return 71;
        }
        memcpy(input_copy, case_text, input_length + 1);
        input.s.pStr = input_copy;
        input.s.nAllocatedLength = (int) input_length + 1;
        input.s.nUsedLength = (int) input_length;
        input.s.nPtr = 0;
        input.type = INCHI_IOS_TYPE_STRING;
    }

    if (!test_case->omit_atom_output && test_case->caller_atom_capacity > 0)
    {
        atoms = (inchi_Atom *) calloc((size_t) test_case->caller_atom_capacity,
                                     sizeof(*atoms));
        if (!atoms)
        {
            fprintf(stderr, "failed to allocate caller atom storage\n");
            return 72;
        }
        memset(atoms, 0x5a,
               (size_t) test_case->caller_atom_capacity * sizeof(*atoms));
        original_atoms = atoms;
    }
    if (!test_case->omit_stereo_output && test_case->caller_stereo_capacity > 0)
    {
        stereo = (inchi_Stereo0D *) calloc((size_t) test_case->caller_stereo_capacity,
                                           sizeof(*stereo));
        if (!stereo)
        {
            fprintf(stderr, "failed to allocate caller stereo storage\n");
            return 73;
        }
        memset(stereo, 0x5a,
               (size_t) test_case->caller_stereo_capacity * sizeof(*stereo));
        original_stereo = stereo;
        stereo_count = test_case->caller_stereo_capacity;
    }
    if (test_case->pass_zero_stereo_count)
    {
        stereo_count = 0;
    }
    passed_atom_capacity = test_case->pass_zero_atom_capacity
                               ? 0
                               : test_case->caller_atom_capacity;
    error_code = test_case->initial_error_code;

    ORACLE_ALLOCATION_CALLS = 0;
    ORACLE_ALLOCATION_ORDINAL = test_case->allocation_failure_ordinal;
    ORACLE_ALLOCATION_FAILURE_ENABLED = test_case->allocation_failure_ordinal > 0;
    ORACLE_DEFER_FREES = 1;

    status = InchiToInchiAtom(&input,
                              test_case->omit_stereo_output ? NULL : &stereo,
                              &stereo_count,
                              test_case->do_not_add_h, 0,
                              test_case->input_type,
                              test_case->omit_atom_output ? NULL : &atoms,
                              passed_atom_capacity, &dimensions,
                              &bond_count,
                              test_case->omit_label_output ? NULL : sdf_label,
                              test_case->omit_value_output ? NULL : sdf_value,
                              test_case->omit_id_output ? NULL : &id,
                              test_case->omit_flags_output ? NULL : &atom_flags,
                              &error_code, error_text);
    ORACLE_ALLOCATION_FAILURE_ENABLED = 0;
    ORACLE_DEFER_FREES = 0;
    atom_count = status > 0 ? status : 0;
    atom_storage_count = atom_count;
    if (!atoms)
    {
        atom_storage_count = 0;
    }
    if (atoms == original_atoms && test_case->caller_atom_capacity > atom_storage_count)
    {
        atom_storage_count = test_case->caller_atom_capacity;
    }
    stereo_storage_count = stereo_count > 0 ? stereo_count : 0;
    if (!stereo)
    {
        stereo_storage_count = 0;
    }
    if (stereo == original_stereo &&
        test_case->caller_stereo_capacity > stereo_storage_count)
    {
        stereo_storage_count = test_case->caller_stereo_capacity;
    }
    input_position = test_case->file_mode ? ftell(input_file) : input.s.nPtr;

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
          "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\","
          "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\","
          "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\","
          "\"api_version\":\"1.07.5\"},\"case_id\":", stdout);
    print_json_string(test_case->case_id);
    fputs(",\"operation\":\"inchi_to_inchi_atom\",\"input\":{\"mode\":", stdout);
    print_json_string(test_case->file_mode ? "file" : "string");
    printf(",\"input_type\":%d,\"text\":", (int) test_case->input_type);
    print_json_string(case_text);
    printf(",\"do_not_add_h\":%d,\"caller_atom_capacity\":%d,"
           "\"caller_stereo_capacity\":%d,\"allocation_failure_ordinal\":%d,"
           "\"omit_atom_output\":%d,\"omit_stereo_output\":%d,"
           "\"omit_label_output\":%d,\"omit_value_output\":%d,"
           "\"omit_id_output\":%d,\"omit_flags_output\":%d,"
           "\"pass_zero_atom_capacity\":%d,\"pass_zero_stereo_count\":%d,"
           "\"initial_error_code\":%d",
           test_case->do_not_add_h, test_case->caller_atom_capacity,
           test_case->caller_stereo_capacity,
           test_case->allocation_failure_ordinal,
           test_case->omit_atom_output, test_case->omit_stereo_output,
           test_case->omit_label_output, test_case->omit_value_output,
           test_case->omit_id_output, test_case->omit_flags_output,
           test_case->pass_zero_atom_capacity,
           test_case->pass_zero_stereo_count,
           test_case->initial_error_code);
    printf("},\"output\":{\"status\":%d,\"error_code\":%d,\"error_text\":",
           status, error_code);
    print_json_string(error_text);
    printf(",\"atom_count\":%d,\"bond_count\":%d,"
           "\"coordinate_dimensions\":%d,\"atoms\":",
           atom_count, bond_count, dimensions);
    print_atoms(atoms, atom_storage_count);
    fputs(",\"bond_fields\":", stdout);
    print_bond_fields(atoms, atoms ? atom_count : 0);
    fputs(",\"stereo0d\":", stdout);
    print_stereo(stereo, stereo_storage_count);
    fputs(",\"sdf_label\":", stdout);
    if (test_case->omit_label_output)
    {
        fputs("null", stdout);
    }
    else
    {
        print_json_string(sdf_label);
    }
    fputs(",\"sdf_value\":", stdout);
    if (test_case->omit_value_output)
    {
        fputs("null", stdout);
    }
    else
    {
        print_json_string(sdf_value);
    }
    fputs(",\"id\":", stdout);
    if (test_case->omit_id_output)
    {
        fputs("null", stdout);
    }
    else
    {
        printf("%ld", id);
    }
    fputs(",\"atom_flags\":", stdout);
    if (test_case->omit_flags_output)
    {
        fputs("null", stdout);
    }
    else
    {
        printf("%llu", (unsigned long long) atom_flags);
    }
    printf(",\"stereo_count\":%d,\"atom_pointer_state\":", stereo_count);
    print_json_string(pointer_state(atoms, original_atoms));
    fputs(",\"stereo_pointer_state\":", stdout);
    print_json_string(pointer_state(stereo, original_stereo));
    printf(",\"input_position\":%ld}}\n", input_position);

    FreeInchi_Atom(&atoms);
    FreeInchi_Stereo0D(&stereo);
    oracle_flush_deferred_frees();
    free(input_copy);
    if (input_file)
    {
        fclose(input_file);
    }
    free(generated_text);
    return 0;
}

static int print_inchi_to_atom_records(void)
{
    size_t i;
    int result;
    for (i = 0; i < sizeof(INCHI_TO_ATOM_CASES) / sizeof(INCHI_TO_ATOM_CASES[0]); i++)
    {
        result = run_inchi_to_atom_case(INCHI_TO_ATOM_CASES + i);
        if (result)
        {
            return result;
        }
    }
    return 0;
}

typedef struct tagOracleComponentCase
{
    const char *case_id;
    AT_NUMB first[3];
    AT_NUMB second[3];
} ORACLE_COMPONENT_CASE;

static const ORACLE_COMPONENT_CASE COMPONENT_CASES[] = {
    { "size-descending", { 9, 4, 0 }, { 3, 0, 0 } },
    { "size-ascending", { 3, 0, 0 }, { 9, 4, 0 } },
    { "size-max-min", { 65535, 0, 0 }, { 0, 0, 0 } },
    { "size-min-max", { 0, 0, 0 }, { 65535, 0, 0 } },
    { "order-ascending", { 7, 2, 0 }, { 7, 8, 0 } },
    { "order-descending", { 7, 8, 0 }, { 7, 2, 0 } },
    { "order-min-max", { 7, 0, 0 }, { 7, 65535, 0 } },
    { "order-max-min", { 7, 65535, 0 }, { 7, 0, 0 } },
    { "equal-keys-third-diff", { 7, 2, 0 }, { 7, 2, 65535 } },
    { "all-equal", { 65535, 65535, 65535 }, { 65535, 65535, 65535 } },
};

static int print_cmp_components_records(void)
{
    size_t i;
    for (i = 0; i < sizeof(COMPONENT_CASES) / sizeof(COMPONENT_CASES[0]); i++)
    {
        const ORACLE_COMPONENT_CASE *test_case = COMPONENT_CASES + i;
        int result = cmp_components(test_case->first, test_case->second);
        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"cmp_components\","
               "\"input\":{\"first\":[%u,%u,%u],"
               "\"second\":[%u,%u,%u]},"
               "\"output\":{\"result\":%d}}\n",
               (unsigned int) test_case->first[0],
               (unsigned int) test_case->first[1],
               (unsigned int) test_case->first[2],
               (unsigned int) test_case->second[0],
               (unsigned int) test_case->second[1],
               (unsigned int) test_case->second[2],
               result);
    }
    return 0;
}

typedef struct tagOracleRadEndpointCase
{
    const char *case_id;
    int first[2];
    int second[2];
} ORACLE_RAD_ENDPOINT_CASE;

static const ORACLE_RAD_ENDPOINT_CASE RAD_ENDPOINT_CASES[] = {
    { "first-field-less", { INT_MIN, INT_MAX }, { 0, 0 } },
    { "first-field-greater", { INT_MAX, INT_MIN }, { 0, 0 } },
    { "second-field-less", { 7, INT_MIN }, { 7, 0 } },
    { "second-field-greater", { 7, INT_MAX }, { 7, 0 } },
    { "both-fields-equal", { 7, 11 }, { 7, 11 } },
};

static int print_cmp_rad_endpoints_records(void)
{
    size_t i;
    for (i = 0; i < sizeof(RAD_ENDPOINT_CASES) / sizeof(RAD_ENDPOINT_CASES[0]); i++)
    {
        const ORACLE_RAD_ENDPOINT_CASE *test_case = RAD_ENDPOINT_CASES + i;
        int result = cmp_rad_endpoints(test_case->first, test_case->second);
        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"cmp_rad_endpoints\","
               "\"input\":{\"first\":[%d,%d],"
               "\"second\":[%d,%d]},"
               "\"output\":{\"result\":%d}}\n",
               test_case->first[0],
               test_case->first[1],
               test_case->second[0],
               test_case->second[1],
               result);
    }
    return 0;
}

typedef struct tagOracleTGroupCase
{
    const char *case_id;
    AT_NUMB first_group_number;
    AT_NUMB second_group_number;
    int context_nonnull;
} ORACLE_T_GROUP_CASE;

static const ORACLE_T_GROUP_CASE T_GROUP_CASES[] = {
    { "minimum-minus-maximum-null-context", 0, 65535, 0 },
    { "maximum-minus-minimum-nonnull-context", 65535, 0, 1 },
    { "equal-null-context", 17, 17, 0 },
    { "ascending-nonnull-context", 11, 29, 1 },
    { "descending-null-context", 29, 11, 0 },
};

static int print_comp_t_group_number_records(void)
{
    size_t i;
    for (i = 0; i < sizeof(T_GROUP_CASES) / sizeof(T_GROUP_CASES[0]); i++)
    {
        const ORACLE_T_GROUP_CASE *test_case = T_GROUP_CASES + i;
        T_GROUP first;
        T_GROUP second;
        T_GROUP first_before;
        T_GROUP second_before;
        int context = 0x13579bdf;
        int context_before = context;
        int result;

        memset(&first, 0xa5, sizeof(first));
        memset(&second, 0x5a, sizeof(second));
        first.nGroupNumber = test_case->first_group_number;
        second.nGroupNumber = test_case->second_group_number;
        memcpy(&first_before, &first, sizeof(first));
        memcpy(&second_before, &second, sizeof(second));
        result = CompTGroupNumber(
            &first,
            &second,
            test_case->context_nonnull ? &context : NULL);

        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"CompTGroupNumber\","
               "\"input\":{\"first_group_number\":%u,"
               "\"second_group_number\":%u,"
               "\"context_nonnull\":%s},"
               "\"output\":{\"result\":%d,"
               "\"first_unchanged\":%s,"
               "\"second_unchanged\":%s,"
               "\"context_unchanged\":%s}}\n",
               (unsigned int) test_case->first_group_number,
               (unsigned int) test_case->second_group_number,
               test_case->context_nonnull ? "true" : "false",
               result,
               memcmp(&first, &first_before, sizeof(first)) == 0 ? "true" : "false",
               memcmp(&second, &second_before, sizeof(second)) == 0 ? "true" : "false",
               context == context_before ? "true" : "false");
    }
    return 0;
}

typedef struct tagOracleCGroupCase
{
    const char *case_id;
    AT_NUMB first_group_number;
    AT_NUMB second_group_number;
    int context_nonnull;
} ORACLE_C_GROUP_CASE;

static const ORACLE_C_GROUP_CASE C_GROUP_CASES[] = {
    { "minimum-minus-maximum-null-context", 0, 65535, 0 },
    { "maximum-minus-minimum-nonnull-context", 65535, 0, 1 },
    { "equal-null-context", 23, 23, 0 },
    { "ascending-nonnull-context", 11, 29, 1 },
    { "descending-null-context", 29, 11, 0 },
};

static int print_comp_c_group_number_records(void)
{
    size_t i;
    for (i = 0; i < sizeof(C_GROUP_CASES) / sizeof(C_GROUP_CASES[0]); i++)
    {
        const ORACLE_C_GROUP_CASE *test_case = C_GROUP_CASES + i;
        C_GROUP first;
        C_GROUP second;
        C_GROUP first_before;
        C_GROUP second_before;
        int context = 0x2468ace;
        int context_before = context;
        int result;

        memset(&first, 0xa5, sizeof(first));
        memset(&second, 0x5a, sizeof(second));
        first.nGroupNumber = test_case->first_group_number;
        second.nGroupNumber = test_case->second_group_number;
        memcpy(&first_before, &first, sizeof(first));
        memcpy(&second_before, &second, sizeof(second));
        result = CompCGroupNumber(
            &first,
            &second,
            test_case->context_nonnull ? &context : NULL);

        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"CompCGroupNumber\","
               "\"input\":{\"first_group_number\":%u,"
               "\"second_group_number\":%u,"
               "\"context_nonnull\":%s},"
               "\"output\":{\"result\":%d,"
               "\"first_unchanged\":%s,"
               "\"second_unchanged\":%s,"
               "\"context_unchanged\":%s}}\n",
               (unsigned int) test_case->first_group_number,
               (unsigned int) test_case->second_group_number,
               test_case->context_nonnull ? "true" : "false",
               result,
               memcmp(&first, &first_before, sizeof(first)) == 0 ? "true" : "false",
               memcmp(&second, &second_before, sizeof(second)) == 0 ? "true" : "false",
               context == context_before ? "true" : "false");
    }
    return 0;
}

typedef struct tagOracleIsoComponentCase
{
    const char *case_id;
    S_CHAR first_iso_atw_diff;
    AT_NUMB first_component;
    S_CHAR second_iso_atw_diff;
    AT_NUMB second_component;
} ORACLE_ISO_COMPONENT_CASE;

static const ORACLE_ISO_COMPONENT_CASE ISO_COMPONENT_CASES[] = {
    { "isotope-minimum-minus-maximum", -128, 65535, 127, 0 },
    { "isotope-maximum-minus-minimum", 127, 0, -128, 65535 },
    { "isotope-priority-over-component", -1, 65535, 0, 0 },
    { "component-minimum-minus-maximum", 7, 0, 7, 65535 },
    { "component-maximum-minus-minimum", 7, 65535, 7, 0 },
    { "component-ascending", 7, 11, 7, 29 },
    { "both-fields-equal", 7, 23, 7, 23 },
};

static int print_cmp_iso_atw_diff_component_no_records(void)
{
    size_t i;
    for (i = 0; i < sizeof(ISO_COMPONENT_CASES) / sizeof(ISO_COMPONENT_CASES[0]); i++)
    {
        const ORACLE_ISO_COMPONENT_CASE *test_case = ISO_COMPONENT_CASES + i;
        inp_ATOM first;
        inp_ATOM second;
        inp_ATOM first_before;
        inp_ATOM second_before;
        int result;

        memset(&first, 0xa5, sizeof(first));
        memset(&second, 0x5a, sizeof(second));
        first.iso_atw_diff = test_case->first_iso_atw_diff;
        first.component = test_case->first_component;
        second.iso_atw_diff = test_case->second_iso_atw_diff;
        second.component = test_case->second_component;
        memcpy(&first_before, &first, sizeof(first));
        memcpy(&second_before, &second, sizeof(second));
        result = cmp_iso_atw_diff_component_no(&first, &second);

        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"cmp_iso_atw_diff_component_no\","
               "\"input\":{\"first_iso_atw_diff\":%d,"
               "\"first_component\":%u,"
               "\"second_iso_atw_diff\":%d,"
               "\"second_component\":%u},"
               "\"output\":{\"result\":%d,"
               "\"first_unchanged\":%s,"
               "\"second_unchanged\":%s}}\n",
               (int) test_case->first_iso_atw_diff,
               (unsigned int) test_case->first_component,
               (int) test_case->second_iso_atw_diff,
               (unsigned int) test_case->second_component,
               result,
               memcmp(&first, &first_before, sizeof(first)) == 0 ? "true" : "false",
               memcmp(&second, &second_before, sizeof(second)) == 0 ? "true" : "false");
    }
    return 0;
}

typedef struct tagOracleNeighListsCase
{
    const char *case_id;
    AT_RANK first_list[3];
    AT_RANK second_list[3];
} ORACLE_NEIGH_LISTS_CASE;

static const ORACLE_NEIGH_LISTS_CASE NEIGH_LISTS_CASES[] = {
    { "first-rank-less", { 1, 1, 0 }, { 1, 2, 0 } },
    { "first-rank-greater", { 1, 2, 0 }, { 1, 1, 0 } },
    { "later-rank-less", { 2, 1, 2 }, { 2, 1, 3 } },
    { "equal", { 2, 1, 2 }, { 2, 1, 2 } },
    { "first-shorter", { 1, 1, 0 }, { 2, 1, 2 } },
    { "first-longer", { 2, 1, 2 }, { 1, 1, 0 } },
    { "rank-boundary", { 1, 0, 0 }, { 1, 4, 0 } },
};

static int print_comp_neigh_lists_records(void)
{
    static const AT_RANK RANK_TEMPLATE[5] = { 0, 10, 20, 30, 65535 };
    size_t i;
    for (i = 0; i < sizeof(NEIGH_LISTS_CASES) / sizeof(NEIGH_LISTS_CASES[0]); i++)
    {
        const ORACLE_NEIGH_LISTS_CASE *test_case = NEIGH_LISTS_CASES + i;
        AT_RANK first_list[3];
        AT_RANK second_list[3];
        AT_RANK first_before[3];
        AT_RANK second_before[3];
        AT_RANK ranks[5];
        AT_RANK ranks_before[5];
        NEIGH_LIST lists[2];
        NEIGH_LIST lists_before[2];
        AT_RANK first_index = 0;
        AT_RANK second_index = 1;
        AT_RANK first_index_before = first_index;
        AT_RANK second_index_before = second_index;
        CANON_GLOBALS globals;
        CANON_GLOBALS globals_before;
        int result;

        memcpy(first_list, test_case->first_list, sizeof(first_list));
        memcpy(second_list, test_case->second_list, sizeof(second_list));
        memcpy(ranks, RANK_TEMPLATE, sizeof(ranks));
        lists[0] = first_list;
        lists[1] = second_list;
        memset(&globals, 0xa5, sizeof(globals));
        globals.m_pNeighList_RankForSort = lists;
        globals.m_pn_RankForSort = ranks;
        memcpy(first_before, first_list, sizeof(first_list));
        memcpy(second_before, second_list, sizeof(second_list));
        memcpy(ranks_before, ranks, sizeof(ranks));
        memcpy(lists_before, lists, sizeof(lists));
        memcpy(&globals_before, &globals, sizeof(globals));

        result = CompNeighLists(&first_index, &second_index, &globals);

        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"CompNeighLists\","
               "\"input\":{\"first_list\":[%u,%u,%u],"
               "\"second_list\":[%u,%u,%u],"
               "\"ranks\":[%u,%u,%u,%u,%u],"
               "\"first_index\":%u,\"second_index\":%u},"
               "\"output\":{\"result\":%d,"
               "\"first_list_unchanged\":%s,"
               "\"second_list_unchanged\":%s,"
               "\"ranks_unchanged\":%s,"
               "\"list_pointers_unchanged\":%s,"
               "\"indices_unchanged\":%s,"
               "\"globals_unchanged\":%s}}\n",
               (unsigned int) first_list[0],
               (unsigned int) first_list[1],
               (unsigned int) first_list[2],
               (unsigned int) second_list[0],
               (unsigned int) second_list[1],
               (unsigned int) second_list[2],
               (unsigned int) ranks[0],
               (unsigned int) ranks[1],
               (unsigned int) ranks[2],
               (unsigned int) ranks[3],
               (unsigned int) ranks[4],
               (unsigned int) first_index,
               (unsigned int) second_index,
               result,
               memcmp(first_list, first_before, sizeof(first_list)) == 0 ? "true" : "false",
               memcmp(second_list, second_before, sizeof(second_list)) == 0 ? "true" : "false",
               memcmp(ranks, ranks_before, sizeof(ranks)) == 0 ? "true" : "false",
               memcmp(lists, lists_before, sizeof(lists)) == 0 ? "true" : "false",
               first_index == first_index_before && second_index == second_index_before ? "true" : "false",
               memcmp(&globals, &globals_before, sizeof(globals)) == 0 ? "true" : "false");
    }
    return 0;
}

typedef struct tagOracleNeighListsUpToMaxRankCase
{
    const char *case_id;
    AT_RANK first_list[4];
    AT_RANK second_list[4];
    AT_RANK max_rank;
} ORACLE_NEIGH_LISTS_UP_TO_MAX_RANK_CASE;

static const ORACLE_NEIGH_LISTS_UP_TO_MAX_RANK_CASE NEIGH_LISTS_UP_TO_MAX_RANK_CASES[] = {
    { "first-rank-less", { 3, 0, 1, 3 }, { 3, 0, 2, 4 }, 10 },
    { "first-rank-greater", { 3, 0, 2, 4 }, { 3, 0, 1, 3 }, 10 },
    { "trimmed-equal", { 3, 0, 1, 3 }, { 3, 0, 2, 4 }, 1 },
    { "first-shorter", { 1, 0, 0, 0 }, { 2, 0, 1, 0 }, 2 },
    { "first-longer", { 2, 0, 1, 0 }, { 1, 0, 0, 0 }, 2 },
    { "all-trimmed", { 1, 3, 0, 0 }, { 1, 4, 0, 0 }, 0 },
    { "rank-maximum-included", { 1, 5, 0, 0 }, { 1, 4, 0, 0 }, 65535 },
    { "rank-maximum-excluded", { 1, 5, 0, 0 }, { 1, 4, 0, 0 }, 65534 },
};

static int print_comp_neigh_lists_up_to_max_rank_records(void)
{
    static const AT_RANK RANK_TEMPLATE[6] = { 1, 2, 3, 9, 10, 65535 };
    size_t i;
    for (i = 0;
         i < sizeof(NEIGH_LISTS_UP_TO_MAX_RANK_CASES) /
                 sizeof(NEIGH_LISTS_UP_TO_MAX_RANK_CASES[0]);
         i++)
    {
        const ORACLE_NEIGH_LISTS_UP_TO_MAX_RANK_CASE *test_case =
            NEIGH_LISTS_UP_TO_MAX_RANK_CASES + i;
        AT_RANK first_list[4];
        AT_RANK second_list[4];
        AT_RANK first_before[4];
        AT_RANK second_before[4];
        AT_RANK ranks[6];
        AT_RANK ranks_before[6];
        NEIGH_LIST lists[2];
        NEIGH_LIST lists_before[2];
        AT_RANK first_index = 0;
        AT_RANK second_index = 1;
        AT_RANK first_index_before = first_index;
        AT_RANK second_index_before = second_index;
        CANON_GLOBALS globals;
        CANON_GLOBALS globals_before;
        int result;

        memcpy(first_list, test_case->first_list, sizeof(first_list));
        memcpy(second_list, test_case->second_list, sizeof(second_list));
        memcpy(ranks, RANK_TEMPLATE, sizeof(ranks));
        lists[0] = first_list;
        lists[1] = second_list;
        memset(&globals, 0xa5, sizeof(globals));
        globals.m_pNeighList_RankForSort = lists;
        globals.m_pn_RankForSort = ranks;
        globals.m_nMaxAtNeighRankForSort = test_case->max_rank;
        memcpy(first_before, first_list, sizeof(first_list));
        memcpy(second_before, second_list, sizeof(second_list));
        memcpy(ranks_before, ranks, sizeof(ranks));
        memcpy(lists_before, lists, sizeof(lists));
        memcpy(&globals_before, &globals, sizeof(globals));

        result = CompNeighListsUpToMaxRank(&first_index, &second_index, &globals);

        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"CompNeighListsUpToMaxRank\","
               "\"input\":{\"first_list\":[%u,%u,%u,%u],"
               "\"second_list\":[%u,%u,%u,%u],"
               "\"ranks\":[%u,%u,%u,%u,%u,%u],"
               "\"max_rank\":%u,"
               "\"first_index\":%u,\"second_index\":%u},"
               "\"output\":{\"result\":%d,"
               "\"first_list_unchanged\":%s,"
               "\"second_list_unchanged\":%s,"
               "\"ranks_unchanged\":%s,"
               "\"list_pointers_unchanged\":%s,"
               "\"indices_unchanged\":%s,"
               "\"globals_unchanged\":%s}}\n",
               (unsigned int) first_list[0],
               (unsigned int) first_list[1],
               (unsigned int) first_list[2],
               (unsigned int) first_list[3],
               (unsigned int) second_list[0],
               (unsigned int) second_list[1],
               (unsigned int) second_list[2],
               (unsigned int) second_list[3],
               (unsigned int) ranks[0],
               (unsigned int) ranks[1],
               (unsigned int) ranks[2],
               (unsigned int) ranks[3],
               (unsigned int) ranks[4],
               (unsigned int) ranks[5],
               (unsigned int) globals.m_nMaxAtNeighRankForSort,
               (unsigned int) first_index,
               (unsigned int) second_index,
               result,
               memcmp(first_list, first_before, sizeof(first_list)) == 0 ? "true" : "false",
               memcmp(second_list, second_before, sizeof(second_list)) == 0 ? "true" : "false",
               memcmp(ranks, ranks_before, sizeof(ranks)) == 0 ? "true" : "false",
               memcmp(lists, lists_before, sizeof(lists)) == 0 ? "true" : "false",
               first_index == first_index_before && second_index == second_index_before ? "true" : "false",
               memcmp(&globals, &globals_before, sizeof(globals)) == 0 ? "true" : "false");
    }
    return 0;
}

typedef struct tagOracleChargeValCase
{
    const char *case_id;
    int first_valence;
    int first_charge;
    int first_order;
    int second_valence;
    int second_charge;
    int second_order;
} ORACLE_CHARGE_VAL_CASE;

static const ORACLE_CHARGE_VAL_CASE CHARGE_VAL_CASES[] = {
    { "valence-less", 1, 0, 0, 2, 0, 0 },
    { "valence-greater", 2, 0, 0, 1, 0, 0 },
    { "absolute-charge-less", 2, 0, 0, 2, 1, 0 },
    { "absolute-charge-greater", 2, 2, 0, 2, 1, 0 },
    { "positive-before-negative", 2, 1, 0, 2, -1, 0 },
    { "negative-after-positive", 2, -1, 0, 2, 1, 0 },
    { "ordering-less", 2, 1, 3, 2, 1, 5 },
    { "ordering-greater", 2, 1, 5, 2, 1, 3 },
    { "equal", 2, 1, 3, 2, 1, 3 },
    { "valence-minimum", INT_MIN, 0, 0, 0, 0, 0 },
    { "valence-maximum", INT_MAX, 0, 0, 0, 0, 0 },
    { "absolute-charge-maximum", 0, INT_MAX, 0, 0, 0, 0 },
    { "ordering-minimum", 0, 0, INT_MIN, 0, 0, 0 },
    { "ordering-maximum", 0, 0, INT_MAX, 0, 0, 0 },
};

static int print_cmp_charge_val_records(void)
{
    size_t i;
    for (i = 0; i < sizeof(CHARGE_VAL_CASES) / sizeof(CHARGE_VAL_CASES[0]); i++)
    {
        const ORACLE_CHARGE_VAL_CASE *test_case = CHARGE_VAL_CASES + i;
        CHARGE_VAL first;
        CHARGE_VAL second;
        CHARGE_VAL first_before;
        CHARGE_VAL second_before;
        int result;

        first.nValence = test_case->first_valence;
        first.nCharge = test_case->first_charge;
        first.nValenceOrderingNumber = test_case->first_order;
        second.nValence = test_case->second_valence;
        second.nCharge = test_case->second_charge;
        second.nValenceOrderingNumber = test_case->second_order;
        memcpy(&first_before, &first, sizeof(first));
        memcpy(&second_before, &second, sizeof(second));
        result = cmp_charge_val(&first, &second, NULL);

        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        printf(",\"operation\":\"cmp_charge_val\","
               "\"input\":{\"first\":{\"valence\":%d,\"charge\":%d,\"order\":%d},"
               "\"second\":{\"valence\":%d,\"charge\":%d,\"order\":%d}},"
               "\"output\":{\"result\":%d,"
               "\"first_unchanged\":%s,"
               "\"second_unchanged\":%s}}\n",
               first.nValence,
               first.nCharge,
               first.nValenceOrderingNumber,
               second.nValence,
               second.nCharge,
               second.nValenceOrderingNumber,
               result,
               memcmp(&first, &first_before, sizeof(first)) == 0 ? "true" : "false",
               memcmp(&second, &second_before, sizeof(second)) == 0 ? "true" : "false");
    }
    return 0;
}

typedef struct tagOracleCcCandCase
{
    const char *case_id;
    CC_CAND first;
    CC_CAND second;
} ORACLE_CC_CAND_CASE;

static const ORACLE_CC_CAND_CASE CC_CAND_CASES[] = {
    { "metal-min-max",
      { 7, 2, 3, SCHAR_MIN, 0, 4, 2, 1, 6 },
      { 7, 2, 3, SCHAR_MAX, 0, 4, 2, 1, 6 } },
    { "metal-max-min",
      { 7, 2, 3, SCHAR_MAX, 0, 4, 2, 1, 6 },
      { 7, 2, 3, SCHAR_MIN, 0, 4, 2, 1, 6 } },
    { "bonds-to-metal",
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 },
      { 7, 2, 3, 0, 5, 4, 2, 1, 6 } },
    { "periodic-row",
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 },
      { 7, 2, 3, 0, 0, 4, 5, 1, 6 } },
    { "bond-count",
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 },
      { 7, 5, 3, 0, 0, 4, 2, 1, 6 } },
    { "chemical-valence",
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 },
      { 7, 2, 5, 0, 0, 4, 2, 1, 6 } },
    { "first-no-valence-electrons",
      { 7, 2, 3, 0, 0, 0, 2, 1, 6 },
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 } },
    { "second-no-valence-electrons",
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 },
      { 7, 2, 3, 0, 0, 0, 2, 1, 6 } },
    { "different-nonzero-valence-electrons",
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 },
      { 7, 2, 3, 0, 0, -4, 2, 1, 6 } },
    { "vertex-minimum-neighbor",
      { INT_MIN, 2, 3, 0, 0, 4, 2, 1, 6 },
      { INT_MIN + 1, 2, 3, 0, 0, 4, 2, 1, 6 } },
    { "vertex-maximum-neighbor",
      { INT_MAX, 2, 3, 0, 0, 4, 2, 1, 6 },
      { INT_MAX - 1, 2, 3, 0, 0, 4, 2, 1, 6 } },
    { "ignored-fields-and-equal",
      { 7, 2, 3, 0, 0, 4, 2, 1, 6 },
      { 7, 2, 3, 0, 0, 4, 2, SCHAR_MAX, UCHAR_MAX } },
};

static void print_cc_cand_json(const CC_CAND *candidate)
{
    printf("{\"iat\":%d,\"num_bonds\":%d,\"chem_valence\":%d,"
           "\"metal\":%d,\"bonds_to_metal\":%d,"
           "\"valence_electrons\":%d,\"periodic_row\":%d,"
           "\"charge_states\":%d,\"element\":%u}",
           (int) candidate->iat,
           (int) candidate->num_bonds,
           (int) candidate->chem_valence,
           (int) candidate->cMetal,
           (int) candidate->cNumBondsToMetal,
           (int) candidate->cNumValenceElectrons,
           (int) candidate->cPeriodicRowNumber,
           (int) candidate->cNumChargeStates,
           (unsigned int) candidate->el_number);
}

static int print_comp_cc_cand_records(void)
{
    size_t i;
    for (i = 0; i < sizeof(CC_CAND_CASES) / sizeof(CC_CAND_CASES[0]); i++)
    {
        const ORACLE_CC_CAND_CASE *test_case = CC_CAND_CASES + i;
        CC_CAND first = test_case->first;
        CC_CAND second = test_case->second;
        CC_CAND first_before;
        CC_CAND second_before;
        int result;

        memcpy(&first_before, &first, sizeof(first));
        memcpy(&second_before, &second, sizeof(second));
        result = comp_cc_cand(&first, &second);

        fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
              "\"case_id\":", stdout);
        print_json_string(test_case->case_id);
        fputs(",\"operation\":\"comp_cc_cand\",\"input\":{\"first\":", stdout);
        print_cc_cand_json(&first);
        fputs(",\"second\":", stdout);
        print_cc_cand_json(&second);
        printf("},\"output\":{\"result\":%d,"
               "\"first_unchanged\":%s,"
               "\"second_unchanged\":%s}}\n",
               result,
               memcmp(&first, &first_before, sizeof(first)) == 0 ? "true" : "false",
               memcmp(&second, &second_before, sizeof(second)) == 0 ? "true" : "false");
    }
    return 0;
}

static void print_base26_triplet_1_record(const char *case_prefix,
                                          unsigned int ordinal,
                                          unsigned char first,
                                          unsigned char second)
{
    unsigned char input[2] = { first, second };
    unsigned char input_before[2];
    const char *result;

    memcpy(input_before, input, sizeof(input));
    result = base26_triplet_1(input);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s-%u\",\"operation\":\"base26_triplet_1\","
           "\"input\":{\"bytes\":[%u,%u]},"
           "\"output\":{\"bytes\":[%u,%u,%u],\"nul\":%u,"
           "\"input_unchanged\":%s}}\n",
           case_prefix,
           ordinal,
           (unsigned int) input[0],
           (unsigned int) input[1],
           (unsigned int) (unsigned char) result[0],
           (unsigned int) (unsigned char) result[1],
           (unsigned int) (unsigned char) result[2],
           (unsigned int) (unsigned char) result[3],
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_base26_triplet_1_records(void)
{
    unsigned int index;
    unsigned int high_bits;
    for (index = 0; index < 16384; index++)
    {
        print_base26_triplet_1_record("index",
                                      index,
                                      (unsigned char) index,
                                      (unsigned char) (index >> 8));
    }
    for (high_bits = 0; high_bits < 4; high_bits++)
    {
        print_base26_triplet_1_record("high-mask",
                                      high_bits,
                                      0x5a,
                                      (unsigned char) (0x15 | (high_bits << 6)));
    }
    return 0;
}

static void print_base26_triplet_2_record(const char *case_prefix,
                                          unsigned int ordinal,
                                          unsigned char first,
                                          unsigned char second,
                                          unsigned char third,
                                          unsigned char fourth)
{
    unsigned char input[4] = { first, second, third, fourth };
    unsigned char input_before[4];
    const char *result;

    memcpy(input_before, input, sizeof(input));
    result = base26_triplet_2(input);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s-%u\",\"operation\":\"base26_triplet_2\","
           "\"input\":{\"bytes\":[%u,%u,%u,%u]},"
           "\"output\":{\"bytes\":[%u,%u,%u],\"nul\":%u,"
           "\"input_unchanged\":%s}}\n",
           case_prefix,
           ordinal,
           (unsigned int) input[0],
           (unsigned int) input[1],
           (unsigned int) input[2],
           (unsigned int) input[3],
           (unsigned int) (unsigned char) result[0],
           (unsigned int) (unsigned char) result[1],
           (unsigned int) (unsigned char) result[2],
           (unsigned int) (unsigned char) result[3],
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_base26_triplet_2_records(void)
{
    static const unsigned char ignored_bits[4][3] = {
        { 0x00, 0x00, 0x00 },
        { 0xff, 0x3f, 0xf0 },
        { 0x5a, 0x15, 0xa0 },
        { 0xa5, 0x2a, 0x50 }
    };
    unsigned int index;
    unsigned int variant;
    unsigned int fixed_index = 0x155a;
    unsigned char second = (unsigned char) ((fixed_index & 0x03) << 6);
    unsigned char third = (unsigned char) ((fixed_index >> 2) & 0xff);
    unsigned char fourth = (unsigned char) ((fixed_index >> 10) & 0x0f);

    for (index = 0; index < 16384; index++)
    {
        print_base26_triplet_2_record(
            "index",
            index,
            0,
            (unsigned char) ((index & 0x03) << 6),
            (unsigned char) ((index >> 2) & 0xff),
            (unsigned char) ((index >> 10) & 0x0f));
    }
    for (variant = 0; variant < 4; variant++)
    {
        print_base26_triplet_2_record(
            "ignored-bits",
            variant,
            ignored_bits[variant][0],
            (unsigned char) (second | ignored_bits[variant][1]),
            third,
            (unsigned char) (fourth | ignored_bits[variant][2]));
    }
    return 0;
}

static void print_base26_triplet_3_record(const char *case_prefix,
                                          unsigned int ordinal,
                                          const unsigned char input_value[6])
{
    unsigned char input[6];
    unsigned char input_before[6];
    const char *result;

    memcpy(input, input_value, sizeof(input));
    memcpy(input_before, input, sizeof(input));
    result = base26_triplet_3(input);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s-%u\",\"operation\":\"base26_triplet_3\","
           "\"input\":{\"bytes\":[%u,%u,%u,%u,%u,%u]},"
           "\"output\":{\"bytes\":[%u,%u,%u],\"nul\":%u,"
           "\"input_unchanged\":%s}}\n",
           case_prefix,
           ordinal,
           (unsigned int) input[0],
           (unsigned int) input[1],
           (unsigned int) input[2],
           (unsigned int) input[3],
           (unsigned int) input[4],
           (unsigned int) input[5],
           (unsigned int) (unsigned char) result[0],
           (unsigned int) (unsigned char) result[1],
           (unsigned int) (unsigned char) result[2],
           (unsigned int) (unsigned char) result[3],
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_base26_triplet_3_records(void)
{
    static const unsigned char ignored_bits[4][5] = {
        { 0x00, 0x00, 0x00, 0x00, 0x00 },
        { 0xff, 0xff, 0xff, 0x0f, 0xfc },
        { 0x5a, 0xa5, 0x3c, 0x05, 0xa8 },
        { 0xa5, 0x5a, 0xc3, 0x0a, 0x54 }
    };
    unsigned int index;
    unsigned int variant;
    unsigned int fixed_index = 0x155a;
    unsigned char canonical[6] = {
        0,
        0,
        0,
        (unsigned char) ((fixed_index & 0x0f) << 4),
        (unsigned char) ((fixed_index >> 4) & 0xff),
        (unsigned char) ((fixed_index >> 12) & 0x03)
    };

    for (index = 0; index < 16384; index++)
    {
        unsigned char input[6] = {
            0,
            0,
            0,
            (unsigned char) ((index & 0x0f) << 4),
            (unsigned char) ((index >> 4) & 0xff),
            (unsigned char) ((index >> 12) & 0x03)
        };
        print_base26_triplet_3_record("index", index, input);
    }
    for (variant = 0; variant < 4; variant++)
    {
        unsigned char input[6] = {
            ignored_bits[variant][0],
            ignored_bits[variant][1],
            ignored_bits[variant][2],
            (unsigned char) (canonical[3] | ignored_bits[variant][3]),
            canonical[4],
            (unsigned char) (canonical[5] | ignored_bits[variant][4])
        };
        print_base26_triplet_3_record("ignored-bits", variant, input);
    }
    return 0;
}

static void print_base26_triplet_4_record(const char *case_prefix,
                                          unsigned int ordinal,
                                          const unsigned char input_value[7])
{
    unsigned char input[7];
    unsigned char input_before[7];
    const char *result;

    memcpy(input, input_value, sizeof(input));
    memcpy(input_before, input, sizeof(input));
    result = base26_triplet_4(input);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s-%u\",\"operation\":\"base26_triplet_4\","
           "\"input\":{\"bytes\":[%u,%u,%u,%u,%u,%u,%u]},"
           "\"output\":{\"bytes\":[%u,%u,%u],\"nul\":%u,"
           "\"input_unchanged\":%s}}\n",
           case_prefix,
           ordinal,
           (unsigned int) input[0],
           (unsigned int) input[1],
           (unsigned int) input[2],
           (unsigned int) input[3],
           (unsigned int) input[4],
           (unsigned int) input[5],
           (unsigned int) input[6],
           (unsigned int) (unsigned char) result[0],
           (unsigned int) (unsigned char) result[1],
           (unsigned int) (unsigned char) result[2],
           (unsigned int) (unsigned char) result[3],
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_base26_triplet_4_records(void)
{
    static const unsigned char ignored_bits[4][6] = {
        { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
        { 0xff, 0xff, 0xff, 0xff, 0xff, 0x03 },
        { 0x5a, 0xa5, 0x3c, 0xc3, 0x69, 0x01 },
        { 0xa5, 0x5a, 0xc3, 0x3c, 0x96, 0x02 }
    };
    unsigned int index;
    unsigned int variant;
    unsigned int fixed_index = 0x155a;
    unsigned char canonical[7] = {
        0,
        0,
        0,
        0,
        0,
        (unsigned char) ((fixed_index & 0x3f) << 2),
        (unsigned char) ((fixed_index >> 6) & 0xff)
    };

    for (index = 0; index < 16384; index++)
    {
        unsigned char input[7] = {
            0,
            0,
            0,
            0,
            0,
            (unsigned char) ((index & 0x3f) << 2),
            (unsigned char) ((index >> 6) & 0xff)
        };
        print_base26_triplet_4_record("index", index, input);
    }
    for (variant = 0; variant < 4; variant++)
    {
        unsigned char input[7] = {
            ignored_bits[variant][0],
            ignored_bits[variant][1],
            ignored_bits[variant][2],
            ignored_bits[variant][3],
            ignored_bits[variant][4],
            (unsigned char) (canonical[5] | ignored_bits[variant][5]),
            canonical[6]
        };
        print_base26_triplet_4_record("ignored-bits", variant, input);
    }
    return 0;
}

static void print_base26_dublet_for_bits_28_to_36_record(
    const char *case_prefix,
    unsigned int ordinal,
    const unsigned char input_value[5])
{
    unsigned char input[5];
    unsigned char input_before[5];
    const char *result;

    memcpy(input, input_value, sizeof(input));
    memcpy(input_before, input, sizeof(input));
    result = base26_dublet_for_bits_28_to_36(input);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s-%u\","
           "\"operation\":\"base26_dublet_for_bits_28_to_36\","
           "\"input\":{\"bytes\":[%u,%u,%u,%u,%u]},"
           "\"output\":{\"bytes\":[%u,%u],\"nul\":%u,"
           "\"input_unchanged\":%s}}\n",
           case_prefix,
           ordinal,
           (unsigned int) input[0],
           (unsigned int) input[1],
           (unsigned int) input[2],
           (unsigned int) input[3],
           (unsigned int) input[4],
           (unsigned int) (unsigned char) result[0],
           (unsigned int) (unsigned char) result[1],
           (unsigned int) (unsigned char) result[2],
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_base26_dublet_for_bits_28_to_36_records(void)
{
    static const unsigned char ignored_bits[4][5] = {
        { 0x00, 0x00, 0x00, 0x00, 0x00 },
        { 0xff, 0xff, 0xff, 0x0f, 0xe0 },
        { 0x5a, 0xa5, 0x3c, 0x05, 0xa0 },
        { 0xa5, 0x5a, 0xc3, 0x0a, 0x40 }
    };
    unsigned int index;
    unsigned int variant;
    unsigned int fixed_index = 0x15a;
    unsigned char canonical[5] = {
        0,
        0,
        0,
        (unsigned char) ((fixed_index & 0x0f) << 4),
        (unsigned char) ((fixed_index >> 4) & 0x1f)
    };

    for (index = 0; index < 512; index++)
    {
        unsigned char input[5] = {
            0,
            0,
            0,
            (unsigned char) ((index & 0x0f) << 4),
            (unsigned char) ((index >> 4) & 0x1f)
        };
        print_base26_dublet_for_bits_28_to_36_record("index", index, input);
    }
    for (variant = 0; variant < 4; variant++)
    {
        unsigned char input[5] = {
            ignored_bits[variant][0],
            ignored_bits[variant][1],
            ignored_bits[variant][2],
            (unsigned char) (canonical[3] | ignored_bits[variant][3]),
            (unsigned char) (canonical[4] | ignored_bits[variant][4])
        };
        print_base26_dublet_for_bits_28_to_36_record(
            "ignored-bits", variant, input);
    }
    return 0;
}

static void print_base26_dublet_for_bits_56_to_64_record(
    const char *case_prefix,
    unsigned int ordinal,
    const unsigned char input_value[9])
{
    unsigned char input[9];
    unsigned char input_before[9];
    const char *result;

    memcpy(input, input_value, sizeof(input));
    memcpy(input_before, input, sizeof(input));
    result = base26_dublet_for_bits_56_to_64(input);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s-%u\","
           "\"operation\":\"base26_dublet_for_bits_56_to_64\","
           "\"input\":{\"bytes\":[%u,%u,%u,%u,%u,%u,%u,%u,%u]},"
           "\"output\":{\"bytes\":[%u,%u],\"nul\":%u,"
           "\"input_unchanged\":%s}}\n",
           case_prefix,
           ordinal,
           (unsigned int) input[0],
           (unsigned int) input[1],
           (unsigned int) input[2],
           (unsigned int) input[3],
           (unsigned int) input[4],
           (unsigned int) input[5],
           (unsigned int) input[6],
           (unsigned int) input[7],
           (unsigned int) input[8],
           (unsigned int) (unsigned char) result[0],
           (unsigned int) (unsigned char) result[1],
           (unsigned int) (unsigned char) result[2],
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_base26_dublet_for_bits_56_to_64_records(void)
{
    static const unsigned char ignored_bits[4][8] = {
        { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
        { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe },
        { 0x5a, 0xa5, 0x3c, 0xc3, 0x69, 0x96, 0x0f, 0xaa },
        { 0xa5, 0x5a, 0xc3, 0x3c, 0x96, 0x69, 0xf0, 0x54 }
    };
    unsigned int index;
    unsigned int variant;
    unsigned int fixed_index = 0x15a;
    unsigned char canonical[9] = {
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (unsigned char) (fixed_index & 0xff),
        (unsigned char) ((fixed_index >> 8) & 0x01)
    };

    for (index = 0; index < 512; index++)
    {
        unsigned char input[9] = {
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            (unsigned char) (index & 0xff),
            (unsigned char) ((index >> 8) & 0x01)
        };
        print_base26_dublet_for_bits_56_to_64_record("index", index, input);
    }
    for (variant = 0; variant < 4; variant++)
    {
        unsigned char input[9] = {
            ignored_bits[variant][0],
            ignored_bits[variant][1],
            ignored_bits[variant][2],
            ignored_bits[variant][3],
            ignored_bits[variant][4],
            ignored_bits[variant][5],
            ignored_bits[variant][6],
            canonical[7],
            (unsigned char) (canonical[8] | ignored_bits[variant][7])
        };
        print_base26_dublet_for_bits_56_to_64_record(
            "ignored-bits", variant, input);
    }
    return 0;
}

static void print_get_xtra_hash_major_hex_record(unsigned int seed)
{
    unsigned char input[32];
    unsigned char input_before[32];
    unsigned char output[64];
    unsigned int index;

    for (index = 0; index < 32; index++)
    {
        input[index] = (unsigned char) (seed + index * 17);
    }
    input[8] = (unsigned char) seed;
    memcpy(input_before, input, sizeof(input));
    memset(output, 0xa5, sizeof(output));
    get_xtra_hash_major_hex(input, (char *) output);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"seed-%u\","
           "\"operation\":\"get_xtra_hash_major_hex\","
           "\"input\":{\"bytes\":",
           seed);
    print_u8_array(input, (int) sizeof(input));
    fputs("},\"output\":{\"bytes\":", stdout);
    print_u8_array(output, (int) sizeof(output));
    printf(",\"input_unchanged\":%s}}\n",
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_get_xtra_hash_major_hex_records(void)
{
    unsigned int seed;
    for (seed = 0; seed < 256; seed++)
    {
        print_get_xtra_hash_major_hex_record(seed);
    }
    return 0;
}

static void print_get_xtra_hash_minor_hex_record(unsigned int seed)
{
    unsigned char input[32];
    unsigned char input_before[32];
    unsigned char output[64];
    unsigned int index;

    for (index = 0; index < 32; index++)
    {
        input[index] = (unsigned char) (seed + index * 17);
    }
    input[4] = (unsigned char) seed;
    memcpy(input_before, input, sizeof(input));
    memset(output, 0xa5, sizeof(output));
    get_xtra_hash_minor_hex(input, (char *) output);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"seed-%u\","
           "\"operation\":\"get_xtra_hash_minor_hex\","
           "\"input\":{\"bytes\":",
           seed);
    print_u8_array(input, (int) sizeof(input));
    fputs("},\"output\":{\"bytes\":", stdout);
    print_u8_array(output, (int) sizeof(output));
    printf(",\"input_unchanged\":%s}}\n",
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_get_xtra_hash_minor_hex_records(void)
{
    unsigned int seed;
    for (seed = 0; seed < 256; seed++)
    {
        print_get_xtra_hash_minor_hex_record(seed);
    }
    return 0;
}

static void print_sha2_starts_record(unsigned int seed)
{
    sha2_context ctx;
    sha2_context before;
    unsigned int index;

    ctx.total[0] =
        (unsigned long) (UINT64_C(0x0123456789abcdef) ^
                         ((uint64_t) seed << 32) ^ seed);
    ctx.total[1] =
        (unsigned long) (UINT64_C(0xfedcba9876543210) ^
                         ((uint64_t) seed << 40) ^ ((uint64_t) seed << 8));
    for (index = 0; index < 8; index++)
    {
        ctx.state[index] =
            (unsigned long) (UINT64_C(0x9e3779b97f4a7c15) *
                             (uint64_t) (index + 1) ^
                             ((uint64_t) seed << ((index * 7) % 56)) ^
                             (uint64_t) (seed + index * 29));
    }
    for (index = 0; index < 64; index++)
    {
        ctx.buffer[index] =
            (unsigned char) (0xa5U ^ seed ^ (index * 37U));
    }
    memcpy(&before, &ctx, sizeof(ctx));

    oracle_sha2_starts(&ctx);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"seed-%u\",\"operation\":\"sha2_starts\","
           "\"input\":{\"total\":",
           seed);
    print_ulong_array(before.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(before.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(before.buffer, 64);
    fputs("},\"output\":{\"total\":", stdout);
    print_ulong_array(ctx.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(ctx.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(ctx.buffer, 64);
    printf(",\"buffer_unchanged\":%s}}\n",
           memcmp(ctx.buffer, before.buffer, sizeof(ctx.buffer)) == 0
               ? "true"
               : "false");
}

static int print_sha2_starts_records(void)
{
    unsigned int seed;
    for (seed = 0; seed < 32; seed++)
    {
        print_sha2_starts_record(seed);
    }
    return 0;
}

typedef struct tagSha2CsumCase
{
    const char *case_id;
    int ilen;
    int null_input;
} Sha2CsumCase;

static const Sha2CsumCase SHA2_CSUM_CASES[] = {
    {"negative-min-null", INT_MIN, 1},
    {"negative-one-null", -1, 1},
    {"zero-null", 0, 1},
    {"zero-nonnull", 0, 0},
    {"length-1", 1, 0},
    {"length-54", 54, 0},
    {"length-55", 55, 0},
    {"length-56", 56, 0},
    {"length-57", 57, 0},
    {"length-63", 63, 0},
    {"length-64", 64, 0},
    {"length-65", 65, 0},
    {"length-119", 119, 0},
    {"length-120", 120, 0},
    {"length-121", 121, 0},
    {"length-127", 127, 0},
    {"length-128", 128, 0},
    {"length-129", 129, 0},
    {"length-255", 255, 0},
    {"length-256", 256, 0},
    {"length-257", 257, 0},
    {"length-511", 511, 0},
    {"length-512", 512, 0},
    {"length-513", 513, 0},
    {"length-1023", 1023, 0},
    {"length-1024", 1024, 0},
};

static void print_sha2_csum_record(const Sha2CsumCase *test_case,
                                   unsigned int case_index)
{
    unsigned char input[1024];
    unsigned char input_before[1024];
    unsigned char output[32];
    unsigned char output_before[32];
    unsigned char *input_pointer;
    unsigned int index;

    for (index = 0; index < sizeof(input); index++)
    {
        input[index] =
            (unsigned char) (case_index * 17U + index * 29U + (index >> 3));
    }
    for (index = 0; index < sizeof(output); index++)
    {
        output[index] = (unsigned char) (0xa5U ^ case_index ^ index);
    }
    memcpy(input_before, input, sizeof(input));
    memcpy(output_before, output, sizeof(output));
    input_pointer = test_case->null_input ? NULL : input;

    oracle_sha2_csum(input_pointer, test_case->ilen, output);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s\",\"operation\":\"sha2_csum\","
           "\"input\":{\"ilen\":%d,\"input_pointer_null\":%s,\"bytes\":",
           test_case->case_id, test_case->ilen,
           test_case->null_input ? "true" : "false");
    print_u8_array(input_before, (int) sizeof(input_before));
    fputs(",\"output\":", stdout);
    print_u8_array(output_before, (int) sizeof(output_before));
    fputs("},\"output\":{\"bytes\":", stdout);
    print_u8_array(input, (int) sizeof(input));
    fputs(",\"digest\":", stdout);
    print_u8_array(output, (int) sizeof(output));
    printf(",\"input_unchanged\":%s}}\n",
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_sha2_csum_records(void)
{
    unsigned int index;
    for (index = 0;
         index < sizeof(SHA2_CSUM_CASES) / sizeof(SHA2_CSUM_CASES[0]);
         index++)
    {
        print_sha2_csum_record(&SHA2_CSUM_CASES[index], index);
    }
    return 0;
}

typedef struct tagInchiKeyCase
{
    const char *case_id;
    const char *text;
    int xtra1;
    int xtra2;
    int null_input;
    int null_xtra1;
    int null_xtra2;
    int allocation_failure_ordinal;
    int generated_minor_length;
} InchiKeyCase;

static const InchiKeyCase INCHI_KEY_CASES[] = {
    {"null-input", "", 1, 1, 1, 0, 0, 0, 0},
    {"short-prefix", "InChI=", 1, 1, 0, 0, 0, 0, 0},
    {"wrong-prefix", "inchi=1S/CH4", 1, 1, 0, 0, 0, 0, 0},
    {"wrong-version", "InChI=2S/CH4", 1, 1, 0, 0, 0, 0, 0},
    {"wrong-slash", "InChI=1S:CH4", 1, 1, 0, 0, 0, 0, 0},
    {"invalid-first-body-byte", "InChI=1S/!bad", 1, 1, 0, 0, 0, 0, 0},
    {"standard-empty", "InChI=1S//", 1, 1, 0, 0, 0, 0, 0},
    {"standard-error", "InChI=1S/?", 1, 1, 0, 0, 0, 0, 0},
    {"standard-methane-xtra", "InChI=1S/CH4/h1H4", 1, 1, 0, 0, 0, 0, 0},
    {"nonstandard-methane-no-xtra", "InChI=1/CH4/h1H4", 0, 0, 0, 0, 0, 0, 0},
    {"beta-methane-null-xtra", "InChI=1B/CH4/h1H4", 1, 1, 0, 1, 1, 0, 0},
    {"major-c-h-q", "InChI=1/C2H6/c1-2/h1-2H3/q+1", 1, 1, 0, 0, 0, 0, 0},
    {"default-minor", "InChI=1/CH4/t1-/m0/s1", 1, 1, 0, 0, 0, 0, 0},
    {"nonstandard-fixed-h", "InChI=1/CH4/fh1H4", 1, 1, 0, 0, 0, 0, 0},
    {"nonstandard-reconnected", "InChI=1/CH4/rC/h1H4", 1, 1, 0, 0, 0, 0, 0},
    {"standard-fixed-h-error", "InChI=1S/CH4/fh1H4", 1, 1, 0, 0, 0, 0, 0},
    {"standard-reconnected-error", "InChI=1S/CH4/rC", 1, 1, 0, 0, 0, 0, 0},
    {"proton-plus-1", "InChI=1S/CH4/h1H4/p+1", 1, 1, 0, 0, 0, 0, 0},
    {"proton-plus-12", "InChI=1S/CH4/h1H4/p+12", 1, 1, 0, 0, 0, 0, 0},
    {"proton-plus-13", "InChI=1S/CH4/h1H4/p+13", 1, 1, 0, 0, 0, 0, 0},
    {"proton-minus-1", "InChI=1S/CH4/h1H4/p-1", 1, 1, 0, 0, 0, 0, 0},
    {"proton-minus-12", "InChI=1S/CH4/h1H4/p-12", 1, 1, 0, 0, 0, 0, 0},
    {"proton-minus-13", "InChI=1S/CH4/h1H4/p-13", 1, 1, 0, 0, 0, 0, 0},
    {"proton-zero", "InChI=1S/CH4/p0", 1, 1, 0, 0, 0, 0, 0},
    {"proton-empty-at-end", "InChI=1S/CH4/p", 1, 1, 0, 0, 0, 0, 0},
    {"proton-short-before-minor", "InChI=1S/CH4/p/t1-", 1, 1, 0, 0, 0, 0, 0},
    {"proton-stops-at-nondigit", "InChI=1S/CH4/p+1junk/t1-", 1, 1, 0, 0, 0, 0, 0},
    {"proton-positive-long-overflow", "InChI=1S/CH4/p9223372036854775808",
     1, 1, 0, 0, 0, 0, 0},
    {"proton-negative-long-overflow", "InChI=1S/CH4/p-9223372036854775809",
     1, 1, 0, 0, 0, 0, 0},
    {"trailing-invalid-substring-byte", "InChI=1S/CH4/h1H4 trailing bytes",
     1, 1, 0, 0, 0, 0, 0},
    {"minor-length-254", NULL, 1, 1, 0, 0, 0, 0, 254},
    {"minor-length-255", NULL, 1, 1, 0, 0, 0, 0, 255},
    {"null-major-xtra-output", "InChI=1S/CH4/h1H4", 1, 1, 0, 1, 0, 0, 0},
    {"null-minor-xtra-output", "InChI=1S/CH4/h1H4", 1, 1, 0, 0, 1, 0, 0},
    {"allocation-failure-smajor", "InChI=1S/CH4/h1H4", 1, 1, 0, 0, 0, 2, 0},
    {"allocation-failure-sminor", "InChI=1S/CH4/h1H4", 1, 1, 0, 0, 0, 3, 0},
    {"allocation-failure-stmp", "InChI=1S/CH4/h1H4", 1, 1, 0, 0, 0, 4, 0},
    {"allocation-failure-sproto", "InChI=1S/CH4/h1H4", 1, 1, 0, 0, 0, 5, 0},
};

static void print_inchi_key_record(const InchiKeyCase *test_case,
                                   unsigned int case_index)
{
    char input[640];
    char input_before[640];
    char key[64];
    char key_before[64];
    char xtra1[64];
    char xtra1_before[64];
    char xtra2[64];
    char xtra2_before[64];
    const char *source;
    int status;
    int allocation_calls;
    int deferred_frees;
    int allocation_failed;
    int successful_allocations;
    int live_allocations;
    unsigned int index;

    for (index = 0; index < sizeof(input); index++)
    {
        input[index] = (char) (0x31U + ((case_index + index * 17U) & 0x3fU));
    }
    for (index = 0; index < sizeof(key); index++)
    {
        key[index] = (char) (0xa5U ^ case_index ^ index);
        xtra1[index] = (char) (0x5aU ^ case_index ^ index);
        xtra2[index] = (char) (0x3cU ^ case_index ^ index);
    }
    if (test_case->generated_minor_length)
    {
        const char *prefix = "InChI=1/C/";
        size_t prefix_length = strlen(prefix);
        size_t minor_body_length =
            (size_t) test_case->generated_minor_length - 1U;
        memcpy(input, prefix, prefix_length);
        memset(input + prefix_length, 'x', minor_body_length);
        input[prefix_length + minor_body_length] = '\0';
    }
    else
    {
        size_t text_length = strlen(test_case->text);
        memcpy(input, test_case->text, text_length + 1);
    }
    memcpy(input_before, input, sizeof(input));
    memcpy(key_before, key, sizeof(key));
    memcpy(xtra1_before, xtra1, sizeof(xtra1));
    memcpy(xtra2_before, xtra2, sizeof(xtra2));

    source = test_case->null_input ? NULL : input;
    ORACLE_ALLOCATION_CALLS = 0;
    ORACLE_ALLOCATION_ORDINAL = test_case->allocation_failure_ordinal;
    ORACLE_ALLOCATION_FAILURE_ENABLED = 1;
    ORACLE_DEFER_FREES = 1;
    status = GetINCHIKeyFromINCHI(
        source, test_case->xtra1, test_case->xtra2, key,
        test_case->null_xtra1 ? NULL : xtra1,
        test_case->null_xtra2 ? NULL : xtra2);
    ORACLE_ALLOCATION_FAILURE_ENABLED = 0;
    ORACLE_DEFER_FREES = 0;

    allocation_calls = ORACLE_ALLOCATION_CALLS;
    deferred_frees = (int) ORACLE_DEFERRED_FREE_COUNT;
    allocation_failed =
        test_case->allocation_failure_ordinal > 0 &&
        allocation_calls >= test_case->allocation_failure_ordinal;
    successful_allocations = allocation_calls - allocation_failed;
    live_allocations = successful_allocations - deferred_frees;

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
          "\"case_id\":", stdout);
    print_json_string(test_case->case_id);
    printf(",\"operation\":\"get_inchi_key_from_inchi\","
           "\"input\":{\"xtra1\":%d,\"xtra2\":%d,"
           "\"input_pointer_null\":%s,\"xtra1_pointer_null\":%s,"
           "\"xtra2_pointer_null\":%s,\"allocation_failure_ordinal\":%d,"
           "\"bytes\":",
           test_case->xtra1, test_case->xtra2,
           test_case->null_input ? "true" : "false",
           test_case->null_xtra1 ? "true" : "false",
           test_case->null_xtra2 ? "true" : "false",
           test_case->allocation_failure_ordinal);
    print_u8_array((const unsigned char *) input_before, (int) sizeof(input_before));
    fputs(",\"key\":", stdout);
    print_u8_array((const unsigned char *) key_before, (int) sizeof(key_before));
    fputs(",\"xtra1_bytes\":", stdout);
    print_u8_array((const unsigned char *) xtra1_before, (int) sizeof(xtra1_before));
    fputs(",\"xtra2_bytes\":", stdout);
    print_u8_array((const unsigned char *) xtra2_before, (int) sizeof(xtra2_before));
    printf("},\"output\":{\"status\":%d,\"bytes\":", status);
    print_u8_array((const unsigned char *) input, (int) sizeof(input));
    fputs(",\"key\":", stdout);
    print_u8_array((const unsigned char *) key, (int) sizeof(key));
    fputs(",\"xtra1_bytes\":", stdout);
    print_u8_array((const unsigned char *) xtra1, (int) sizeof(xtra1));
    fputs(",\"xtra2_bytes\":", stdout);
    print_u8_array((const unsigned char *) xtra2, (int) sizeof(xtra2));
    printf(",\"input_unchanged\":%s,\"allocation_calls\":%d,"
           "\"deferred_frees\":%d,\"live_allocations\":%d}}\n",
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false",
           allocation_calls, deferred_frees, live_allocations);

    oracle_flush_deferred_frees();
}

static int print_inchi_key_records(void)
{
    unsigned int index;
    for (index = 0;
         index < sizeof(INCHI_KEY_CASES) / sizeof(INCHI_KEY_CASES[0]);
         index++)
    {
        print_inchi_key_record(&INCHI_KEY_CASES[index], index);
    }
    return 0;
}

static void print_sha2_process_record(unsigned int seed)
{
    sha2_context ctx;
    sha2_context before;
    unsigned char data[64];
    unsigned char data_before[64];
    unsigned int index;

    ctx.total[0] =
        (unsigned long) (UINT64_C(0x0123456789abcdef) ^
                         ((uint64_t) seed << 32) ^ seed);
    ctx.total[1] =
        (unsigned long) (UINT64_C(0xfedcba9876543210) ^
                         ((uint64_t) seed << 40) ^ ((uint64_t) seed << 8));
    for (index = 0; index < 8; index++)
    {
        ctx.state[index] =
            (unsigned long) (UINT64_C(0x9e3779b97f4a7c15) *
                             (uint64_t) (index + 1) ^
                             ((uint64_t) seed << ((index * 7) % 56)) ^
                             (uint64_t) (seed + index * 29));
    }
    for (index = 0; index < 64; index++)
    {
        ctx.buffer[index] =
            (unsigned char) (0xa5U ^ seed ^ (index * 37U));
        data[index] =
            (unsigned char) (seed + index * 17U + (index >> 1));
    }
    memcpy(&before, &ctx, sizeof(ctx));
    memcpy(data_before, data, sizeof(data));

    oracle_sha2_process(&ctx, data);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"seed-%u\",\"operation\":\"sha2_process\","
           "\"input\":{\"total\":",
           seed);
    print_ulong_array(before.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(before.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(before.buffer, 64);
    fputs(",\"data\":", stdout);
    print_u8_array(data_before, 64);
    fputs("},\"output\":{\"total\":", stdout);
    print_ulong_array(ctx.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(ctx.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(ctx.buffer, 64);
    fputs(",\"data\":", stdout);
    print_u8_array(data, 64);
    printf(",\"total_unchanged\":%s,\"buffer_unchanged\":%s,"
           "\"data_unchanged\":%s}}\n",
           memcmp(ctx.total, before.total, sizeof(ctx.total)) == 0
               ? "true"
               : "false",
           memcmp(ctx.buffer, before.buffer, sizeof(ctx.buffer)) == 0
               ? "true"
               : "false",
           memcmp(data, data_before, sizeof(data)) == 0 ? "true" : "false");
}

static int print_sha2_process_records(void)
{
    unsigned int seed;
    for (seed = 0; seed < 256; seed++)
    {
        print_sha2_process_record(seed);
    }
    return 0;
}

typedef struct tagSha2UpdateCase
{
    const char *case_id;
    unsigned long total0;
    unsigned long total1;
    int ilen;
} Sha2UpdateCase;

static const Sha2UpdateCase SHA2_UPDATE_CASES[] = {
    {"negative-min", 17UL, 23UL, INT_MIN},
    {"negative-one", 17UL, 23UL, -1},
    {"zero", 17UL, 23UL, 0},
    {"empty-left-one", 0UL, 7UL, 1},
    {"empty-left-63", 0UL, 7UL, 63},
    {"empty-left-64", 0UL, 7UL, 64},
    {"empty-left-65", 0UL, 7UL, 65},
    {"empty-left-127", 0UL, 7UL, 127},
    {"empty-left-128", 0UL, 7UL, 128},
    {"empty-left-129", 0UL, 7UL, 129},
    {"left-1-fill-minus-one", 1UL, 11UL, 62},
    {"left-1-fill", 1UL, 11UL, 63},
    {"left-1-fill-plus-one", 1UL, 11UL, 64},
    {"left-1-fill-plus-block", 1UL, 11UL, 127},
    {"left-1-fill-plus-block-tail", 1UL, 11UL, 128},
    {"left-13-fill-minus-one", 13UL, 19UL, 50},
    {"left-13-fill", 13UL, 19UL, 51},
    {"left-13-fill-plus-one", 13UL, 19UL, 52},
    {"left-13-fill-plus-block", 13UL, 19UL, 115},
    {"left-13-fill-plus-block-tail", 13UL, 19UL, 116},
    {"left-48-fill-minus-one", 48UL, 29UL, 15},
    {"left-48-fill", 48UL, 29UL, 16},
    {"left-48-fill-plus-one", 48UL, 29UL, 17},
    {"left-48-fill-plus-block", 48UL, 29UL, 80},
    {"left-48-fill-plus-block-tail", 48UL, 29UL, 81},
    {"left-48-two-blocks-tail", 48UL, 29UL, 145},
    {"left-63-fill", 63UL, 31UL, 1},
    {"left-63-fill-plus-one", 63UL, 31UL, 2},
    {"left-63-fill-plus-block", 63UL, 31UL, 65},
    {"left-63-fill-plus-block-tail", 63UL, 31UL, 66},
    {"low-counter-no-carry", 0xfffffff0UL, 0x1234UL, 15},
    {"low-counter-carry-at-boundary", 0xfffffff0UL, 0x1234UL, 16},
    {"low-counter-carry-plus-one", 0xfffffff0UL, 0x1234UL, 17},
    {"low-counter-carry-multiple-blocks", 0xfffffff0UL, 0x1234UL, 145},
    {"high-counter-wrap", 0xffffffffUL, ULONG_MAX, 1},
    {"lp64-high-total-bits", (unsigned long) UINT64_C(0x123456780000003f),
     (unsigned long) UINT64_C(0xfedcba9876543210), 1},
};

static void print_sha2_update_record(const Sha2UpdateCase *test_case,
                                     unsigned int case_index)
{
    sha2_context ctx;
    sha2_context before;
    unsigned char input[256];
    unsigned char input_before[256];
    unsigned int index;

    ctx.total[0] = test_case->total0;
    ctx.total[1] = test_case->total1;
    for (index = 0; index < 8; index++)
    {
        ctx.state[index] =
            (unsigned long) (UINT64_C(0x6a09e667f3bcc908) ^
                             (UINT64_C(0x9e3779b97f4a7c15) *
                              (uint64_t) (index + case_index + 1)));
    }
    for (index = 0; index < 64; index++)
    {
        ctx.buffer[index] =
            (unsigned char) (0xa5U ^ (case_index * 11U) ^ (index * 37U));
    }
    for (index = 0; index < sizeof(input); index++)
    {
        input[index] =
            (unsigned char) (case_index * 17U + index * 29U + (index >> 2));
    }
    memcpy(&before, &ctx, sizeof(ctx));
    memcpy(input_before, input, sizeof(input));

    oracle_sha2_update(&ctx, input, test_case->ilen);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s\",\"operation\":\"sha2_update\","
           "\"input\":{\"ilen\":%d,\"total\":",
           test_case->case_id, test_case->ilen);
    print_ulong_array(before.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(before.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(before.buffer, 64);
    fputs(",\"bytes\":", stdout);
    print_u8_array(input_before, (int) sizeof(input_before));
    fputs("},\"output\":{\"total\":", stdout);
    print_ulong_array(ctx.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(ctx.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(ctx.buffer, 64);
    fputs(",\"bytes\":", stdout);
    print_u8_array(input, (int) sizeof(input));
    printf(",\"input_unchanged\":%s}}\n",
           memcmp(input, input_before, sizeof(input)) == 0 ? "true" : "false");
}

static int print_sha2_update_records(void)
{
    unsigned int index;
    for (index = 0;
         index < sizeof(SHA2_UPDATE_CASES) / sizeof(SHA2_UPDATE_CASES[0]);
         index++)
    {
        print_sha2_update_record(&SHA2_UPDATE_CASES[index], index);
    }
    return 0;
}

typedef struct tagSha2FinishCase
{
    const char *case_id;
    unsigned long total0;
    unsigned long total1;
} Sha2FinishCase;

static const Sha2FinishCase SHA2_FINISH_CASES[] = {
    {"last-0", 0UL, 0UL},
    {"last-1", 1UL, 3UL},
    {"last-54", 54UL, 5UL},
    {"last-55", 55UL, 7UL},
    {"last-56", 56UL, 11UL},
    {"last-57", 57UL, 13UL},
    {"last-63", 63UL, 17UL},
    {"carry-in-message-length", 0xfffffff7UL, 0x1234UL},
    {"carry-in-padding", 0xffffffffUL, 0x5678UL},
    {"carry-wraps-high-counter", 0xffffffffUL, ULONG_MAX},
    {"lp64-high-total-last-55",
     (unsigned long) UINT64_C(0x1234567800000037),
     (unsigned long) UINT64_C(0xfedcba9876543210)},
    {"lp64-high-total-last-56",
     (unsigned long) UINT64_C(0xabcdef0100000038),
     (unsigned long) UINT64_C(0x0123456789abcdef)},
};

static void print_sha2_finish_record(const Sha2FinishCase *test_case,
                                     unsigned int case_index)
{
    sha2_context ctx;
    sha2_context before;
    unsigned char output[32];
    unsigned char output_before[32];
    unsigned int index;

    ctx.total[0] = test_case->total0;
    ctx.total[1] = test_case->total1;
    for (index = 0; index < 8; index++)
    {
        ctx.state[index] =
            (unsigned long) (UINT64_C(0x6a09e667f3bcc908) ^
                             (UINT64_C(0x9e3779b97f4a7c15) *
                              (uint64_t) (index + case_index + 1)));
    }
    for (index = 0; index < 64; index++)
    {
        ctx.buffer[index] =
            (unsigned char) (0x5aU ^ (case_index * 19U) ^ (index * 41U));
    }
    for (index = 0; index < sizeof(output); index++)
    {
        output[index] = (unsigned char) (0xa5U ^ index ^ case_index);
    }
    memcpy(&before, &ctx, sizeof(ctx));
    memcpy(output_before, output, sizeof(output));

    oracle_sha2_finish(&ctx, output);

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s\",\"operation\":\"sha2_finish\","
           "\"input\":{\"total\":",
           test_case->case_id);
    print_ulong_array(before.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(before.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(before.buffer, 64);
    fputs(",\"output\":", stdout);
    print_u8_array(output_before, 32);
    fputs("},\"output\":{\"total\":", stdout);
    print_ulong_array(ctx.total, 2);
    fputs(",\"state\":", stdout);
    print_ulong_array(ctx.state, 8);
    fputs(",\"buffer\":", stdout);
    print_u8_array(ctx.buffer, 64);
    fputs(",\"digest\":", stdout);
    print_u8_array(output, 32);
    fputs("}}\n", stdout);
}

static int print_sha2_finish_records(void)
{
    unsigned int index;
    for (index = 0;
         index < sizeof(SHA2_FINISH_CASES) / sizeof(SHA2_FINISH_CASES[0]);
         index++)
    {
        print_sha2_finish_record(&SHA2_FINISH_CASES[index], index);
    }
    return 0;
}

static void print_nullable_json_string(const char *text)
{
    if (text)
    {
        print_json_string(text);
    }
    else
    {
        fputs("null", stdout);
    }
}

static void print_root_atom(const inchi_Atom *atom)
{
    printf("{\"coordinate_bits\":[%" PRIu64 ",%" PRIu64 ",%" PRIu64 "],"
           "\"neighbor\":",
           double_bits(atom->x), double_bits(atom->y), double_bits(atom->z));
    print_i16_array(atom->neighbor, MAXVAL);
    fputs(",\"bond_type\":", stdout);
    print_i8_array(atom->bond_type, MAXVAL);
    fputs(",\"bond_stereo\":", stdout);
    print_i8_array(atom->bond_stereo, MAXVAL);
    fputs(",\"elname\":", stdout);
    print_char_array(atom->elname, ATOM_EL_LEN);
    printf(",\"num_bonds\":%d,\"num_iso_h\":", (int) atom->num_bonds);
    print_i8_array(atom->num_iso_H, NUM_H_ISOTOPES + 1);
    printf(",\"isotopic_mass\":%d,\"radical\":%d,\"charge\":%d}",
           (int) atom->isotopic_mass, (int) atom->radical, (int) atom->charge);
}

static void print_root_stereo(const inchi_Stereo0D *stereo)
{
    fputs("{\"neighbor\":", stdout);
    print_i16_array(stereo->neighbor, 4);
    printf(",\"central_atom\":%d,\"type\":%d,\"parity\":%d}",
           (int) stereo->central_atom, (int) stereo->type, (int) stereo->parity);
}

static void set_root_element(inchi_Atom *atom, const char *element)
{
    size_t length = strlen(element);
    if (length >= ATOM_EL_LEN)
    {
        abort();
    }
    memcpy(atom->elname, element, length + 1);
}

typedef struct tagRootGetInchiCase
{
    const char *case_id;
    int kind;
    const char *options;
} RootGetInchiCase;

static const RootGetInchiCase ROOT_GET_INCHI_CASES[] = {
    {"generate-methane-standard", 0, NULL},
    {"generate-methane-fixedh-options", 0, "-FixedH -RecMet -SUU -SLUUD"},
    {"generate-tetrahedral-0d", 1, NULL},
    {"generate-tetrahedral-0d-relative", 1, "-AuxNone -SRel"},
    {"generate-tetrahedral-0d-racemic", 1, "-AuxNone -SRac"},
    {"generate-carbon-2d", 2, NULL},
    {"generate-carbon-3d-isotope-charge-radical", 3, "-FixedH"},
    {"generate-carbon-iron-reconnected", 4, "-RecMet -FixedH"},
    {"generate-pseudoatom-error", 5, NULL},
};

static void print_root_get_inchi_record(const RootGetInchiCase *test_case)
{
    inchi_Atom atoms[5];
    inchi_Atom atoms_before[5];
    inchi_Stereo0D stereo[1];
    inchi_Stereo0D stereo_before[1];
    char options[96];
    char options_before[96];
    inchi_Input input;
    inchi_Output output;
    int status;
    int atom_count = 1;
    int stereo_count = 0;
    int input_unchanged;

    memset(atoms, 0, sizeof(atoms));
    memset(stereo, 0, sizeof(stereo));
    memset(options, 0xa5, sizeof(options));
    if (test_case->options)
    {
        size_t length = strlen(test_case->options);
        memcpy(options, test_case->options, length + 1);
    }

    if (test_case->kind == 0)
    {
        set_root_element(&atoms[0], "C");
        atoms[0].num_iso_H[0] = -1;
    }
    else if (test_case->kind == 1)
    {
        int index;
        static const char *elements[] = {"C", "F", "Cl", "Br", "I"};
        atom_count = 5;
        stereo_count = 1;
        for (index = 0; index < atom_count; index++)
        {
            set_root_element(&atoms[index], elements[index]);
        }
        for (index = 0; index < 4; index++)
        {
            atoms[0].neighbor[index] = (AT_NUM) (index + 1);
            atoms[0].bond_type[index] = INCHI_BOND_TYPE_SINGLE;
            atoms[index + 1].neighbor[0] = 0;
            atoms[index + 1].bond_type[0] = INCHI_BOND_TYPE_SINGLE;
            atoms[index + 1].num_bonds = 1;
            stereo[0].neighbor[index] = (AT_NUM) (index + 1);
        }
        atoms[0].num_bonds = 4;
        stereo[0].central_atom = 0;
        stereo[0].type = INCHI_StereoType_Tetrahedral;
        stereo[0].parity = INCHI_PARITY_ODD;
    }
    else if (test_case->kind == 2)
    {
        set_root_element(&atoms[0], "C");
        atoms[0].num_iso_H[0] = -1;
        atoms[0].x = 1.25;
        atoms[0].y = -2.5;
    }
    else if (test_case->kind == 3)
    {
        set_root_element(&atoms[0], "C");
        atoms[0].num_iso_H[0] = -1;
        atoms[0].x = 1.25;
        atoms[0].y = -2.5;
        atoms[0].z = 3.75;
        atoms[0].isotopic_mass = 13;
        atoms[0].charge = 1;
        atoms[0].radical = INCHI_RADICAL_DOUBLET;
    }
    else if (test_case->kind == 4)
    {
        atom_count = 2;
        set_root_element(&atoms[0], "C");
        set_root_element(&atoms[1], "Fe");
        atoms[0].num_iso_H[0] = -1;
        atoms[0].neighbor[0] = 1;
        atoms[0].bond_type[0] = INCHI_BOND_TYPE_SINGLE;
        atoms[0].num_bonds = 1;
        atoms[1].neighbor[0] = 0;
        atoms[1].bond_type[0] = INCHI_BOND_TYPE_SINGLE;
        atoms[1].num_bonds = 1;
    }
    else
    {
        set_root_element(&atoms[0], "*");
    }

    memcpy(atoms_before, atoms, sizeof(atoms));
    memcpy(stereo_before, stereo, sizeof(stereo));
    memcpy(options_before, options, sizeof(options));
    memset(&input, 0, sizeof(input));
    input.atom = atoms;
    input.num_atoms = (AT_NUM) atom_count;
    input.stereo0D = stereo_count ? stereo : NULL;
    input.num_stereo0D = (AT_NUM) stereo_count;
    input.szOptions = test_case->options ? options : NULL;
    memset(&output, 0, sizeof(output));

    status = GetINCHI(&input, &output);
    input_unchanged =
        memcmp(atoms, atoms_before, sizeof(atoms)) == 0 &&
        memcmp(stereo, stereo_before, sizeof(stereo)) == 0 &&
        memcmp(options, options_before, sizeof(options)) == 0;

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\",\"case_id\":",
          stdout);
    print_json_string(test_case->case_id);
    printf(",\"operation\":\"get_inchi_root\",\"input\":{\"kind\":%d,"
           "\"options\":",
           test_case->kind);
    print_nullable_json_string(test_case->options);
    printf("},\"output\":{\"status\":%d,\"inchi\":", status);
    print_nullable_json_string(output.szInChI);
    fputs(",\"aux\":", stdout);
    print_nullable_json_string(output.szAuxInfo);
    fputs(",\"message\":", stdout);
    print_nullable_json_string(output.szMessage);
    fputs(",\"log\":", stdout);
    print_nullable_json_string(output.szLog);
    printf(",\"input_unchanged\":%s", input_unchanged ? "true" : "false");

    FreeINCHI(&output);
    printf(",\"after_free\":{\"inchi_null\":%s,\"aux_null\":%s,"
           "\"message_null\":%s,\"log_null\":%s}}}\n",
           output.szInChI ? "false" : "true",
           output.szAuxInfo ? "false" : "true",
           output.szMessage ? "false" : "true",
           output.szLog ? "false" : "true");
}

typedef struct tagRootGetStructCase
{
    const char *case_id;
    const char *inchi;
    const char *options;
} RootGetStructCase;

static const RootGetStructCase ROOT_GET_STRUCT_CASES[] = {
    {"parse-methane", "InChI=1S/CH4/h1H4", NULL},
    {"parse-warning",
     "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1", NULL},
    {"parse-isotope-charge", "InChI=1S/CH5N/c1-2/h2H2,1H3/p+1/i1+1", NULL},
    {"parse-tetrahedral", "InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m1/s1", NULL},
    {"parse-fixed-h",
     "InChI=1/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/"
     "h1H,(H4,6,7,8,9,10,11)/f/h8,10H,6H2",
     NULL},
    {"parse-reconnected", "InChI=1/CH4/h1H4/rCH4/h1H4", NULL},
    {"parse-phosphoserine",
     "InChI=1S/C3H8NO6P/c4-2(3(5)6)1-10-11(7,8)9/h2H,1,4H2,"
     "(H,5,6)(H2,7,8,9)/t2-/m0/s1",
     NULL},
    {"parse-n-methylpyridinium",
     "InChI=1S/C6H8N/c1-7-5-3-2-4-6-7/h2-6H,1H3/q+1", NULL},
    {"parse-polymer",
     "InChI=1B/C4H4N4.2Zz/c1-5-2-7-4-8-3-6-1;;/h1-4H;;/"
     "z101-1-8(9,10-8,3,1,6,2,5,2,7,3,6,1,5,4,7,4,8)/"
     "b5-1-,5-2+,6-1+,6-3-,7-2+,7-4+,8-3+,8-4+;;",
     NULL},
    {"parse-malformed", "bad", NULL},
};

static void print_root_get_struct_record(const RootGetStructCase *test_case)
{
    char input_text[512];
    char input_before[512];
    char options[96];
    char options_before[96];
    inchi_InputINCHI input;
    inchi_OutputStruct output;
    int status;
    int index;
    int input_unchanged;

    memset(input_text, 0xa5, sizeof(input_text));
    memcpy(input_text, test_case->inchi, strlen(test_case->inchi) + 1);
    memcpy(input_before, input_text, sizeof(input_text));
    memset(options, 0xa5, sizeof(options));
    if (test_case->options)
    {
        memcpy(options, test_case->options, strlen(test_case->options) + 1);
    }
    memcpy(options_before, options, sizeof(options));
    input.szInChI = input_text;
    input.szOptions = test_case->options ? options : NULL;
    memset(&output, 0xa5, sizeof(output));

    status = GetStructFromINCHI(&input, &output);
    input_unchanged =
        memcmp(input_text, input_before, sizeof(input_text)) == 0 &&
        memcmp(options, options_before, sizeof(options)) == 0;

    fputs("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\",\"case_id\":",
          stdout);
    print_json_string(test_case->case_id);
    fputs(",\"operation\":\"get_struct_from_inchi_root\",\"input\":{\"inchi\":",
          stdout);
    print_json_string(test_case->inchi);
    fputs(",\"options\":", stdout);
    print_nullable_json_string(test_case->options);
    printf("},\"output\":{\"status\":%d,\"num_atoms\":%d,"
           "\"num_stereo\":%d,\"message\":",
           status, (int) output.num_atoms, (int) output.num_stereo0D);
    print_nullable_json_string(output.szMessage);
    fputs(",\"log\":", stdout);
    print_nullable_json_string(output.szLog);
    fputs(",\"warning_flags\":[", stdout);
    printf("%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 "],\"atoms\":[",
           (uint64_t) output.WarningFlags[0][0],
           (uint64_t) output.WarningFlags[0][1],
           (uint64_t) output.WarningFlags[1][0],
           (uint64_t) output.WarningFlags[1][1]);
    for (index = 0; index < output.num_atoms; index++)
    {
        if (index)
        {
            putchar(',');
        }
        print_root_atom(output.atom + index);
    }
    fputs("],\"stereo\":[", stdout);
    for (index = 0; index < output.num_stereo0D; index++)
    {
        if (index)
        {
            putchar(',');
        }
        print_root_stereo(output.stereo0D + index);
    }
    printf("],\"input_unchanged\":%s", input_unchanged ? "true" : "false");

    FreeStructFromINCHI(&output);
    printf(",\"after_free\":{\"atom_null\":%s,\"stereo_null\":%s,"
           "\"message_null\":%s,\"log_null\":%s,\"num_atoms\":%d,"
           "\"num_stereo\":%d,\"warning_flags\":[%" PRIu64 ",%" PRIu64
           ",%" PRIu64 ",%" PRIu64 "]}}}\n",
           output.atom ? "false" : "true",
           output.stereo0D ? "false" : "true",
           output.szMessage ? "false" : "true",
           output.szLog ? "false" : "true",
           (int) output.num_atoms, (int) output.num_stereo0D,
           (uint64_t) output.WarningFlags[0][0],
           (uint64_t) output.WarningFlags[0][1],
           (uint64_t) output.WarningFlags[1][0],
           (uint64_t) output.WarningFlags[1][1]);
}

static int print_rdkit_core_root_records(void)
{
    unsigned int index;
    for (index = 0;
         index < sizeof(ROOT_GET_INCHI_CASES) / sizeof(ROOT_GET_INCHI_CASES[0]);
         index++)
    {
        print_root_get_inchi_record(ROOT_GET_INCHI_CASES + index);
    }
    for (index = 0;
         index < sizeof(ROOT_GET_STRUCT_CASES) / sizeof(ROOT_GET_STRUCT_CASES[0]);
         index++)
    {
        print_root_get_struct_record(ROOT_GET_STRUCT_CASES + index);
    }
    return 0;
}

static void cn_list_hash_i32(uint64_t *hash, int value)
{
    uint32_t bits = (uint32_t) value;
    unsigned int byte_index;
    for (byte_index = 0; byte_index < 4; byte_index++)
    {
        *hash ^= (bits >> (byte_index * 8)) & UINT32_C(0xff);
        *hash *= UINT64_C(1099511628211);
    }
}

static int print_cn_list_record(void)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    int list_index;
    for (list_index = 0; list_index < cnListNumEl; list_index++)
    {
        int node_index;
        cn_list_hash_i32(&hash, cnList[list_index].bits);
        for (node_index = 0; node_index < cnList[list_index].len; node_index++)
        {
            int edge_index;
            const C_NODE *node = cnList[list_index].pCN + node_index;
            cn_list_hash_i32(&hash, node->v.type);
            cn_list_hash_i32(&hash, node->v.cap);
            cn_list_hash_i32(&hash, node->v.flow);
            for (edge_index = 0; edge_index < MAX_CN_VAL; edge_index++)
            {
                const ECF *edge = node->e + edge_index;
                cn_list_hash_i32(&hash, edge->neigh);
                cn_list_hash_i32(&hash, edge->cap);
                cn_list_hash_i32(&hash, edge->bForbiddenEdge);
                cn_list_hash_i32(&hash, edge->flow);
            }
        }
    }
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"case_id\":\"cn-list-active-table\",\"operation\":\"cn_list\"," 
           "\"output\":{\"count\":%d,\"hash\":\"%" PRIu64 "\"," 
           "\"cn_p_edge\":{\"neighbor\":%d,\"cap\":%d," 
           "\"forbidden\":%d,\"flow\":%d}}}\n",
           cnListNumEl, hash,
           (int) cnList[14].pCN[0].e[0].neigh,
           (int) cnList[14].pCN[0].e[0].cap,
           (int) cnList[14].pCN[0].e[0].bForbiddenEdge,
           (int) cnList[14].pCN[0].e[0].flow);
    return 0;
}

static int print_normalize_and_compare_negative_case(unsigned int holder_mask,
                                                     int forced_return)
{
    CANON_GLOBALS canonical_globals;
    CANON_GLOBALS canonical_globals_before;
    struct tagINCHI_CLOCK clock;
    struct tagINCHI_CLOCK clock_before;
    INPUT_PARMS parameters;
    INPUT_PARMS parameters_before;
    STRUCT_DATA data;
    STRUCT_DATA data_before;
    BN_STRUCT bns;
    BN_STRUCT bns_before;
    BN_DATA bns_data;
    BN_DATA bns_data_before;
    BNS_VERTEX vertex;
    BNS_VERTEX vertex_before;
    StrFromINChI structure;
    StrFromINChI expected_structure;
    inp_ATOM atoms[1];
    inp_ATOM atoms_before[1];
    inp_ATOM atoms2[1];
    inp_ATOM atoms2_before[1];
    inp_ATOM atoms3[1];
    inp_ATOM atoms3_before[1];
    VAL_AT valence[1];
    VAL_AT valence_before[1];
    ALL_TC_GROUPS groups;
    ALL_TC_GROUPS groups_before;
    SRM restore_mode;
    INChI *input[TAUT_NUM] = {NULL, NULL};
    static const char *inchi_labels[TAUT_NUM] = {
        "pOneINChI[0]", "pOneINChI[1]"};
    static const char *aux_labels[TAUT_NUM] = {
        "pOneINChI_Aux[0]", "pOneINChI_Aux[1]"};
    static const char *norm_labels[TAUT_NUM] = {
        "pOne_norm_data[0]", "pOne_norm_data[1]"};
    char case_id[96];
    int index;
    int result;
    int runs = 17;
    int delta = -19;
    int complete_caller_state_exact;
    int allocation_free_exact;
    int holders_null;
    int one_ti_cleared;
    size_t source_allocations_before;
    size_t source_allocations_after;

    memset(&canonical_globals, 0, sizeof(canonical_globals));
    memset(&clock, 0, sizeof(clock));
    memset(&parameters, 0, sizeof(parameters));
    memset(&data, 0, sizeof(data));
    memset(&bns, 0, sizeof(bns));
    memset(&bns_data, 0, sizeof(bns_data));
    memset(&vertex, 0, sizeof(vertex));
    memset(&structure, 0, sizeof(structure));
    memset(atoms, 0, sizeof(atoms));
    memset(atoms2, 0, sizeof(atoms2));
    memset(atoms3, 0, sizeof(atoms3));
    memset(valence, 0, sizeof(valence));
    memset(&groups, 0, sizeof(groups));
    memset(&restore_mode, 0, sizeof(restore_mode));

    parameters.nMode = REQ_MODE_TAUT | REQ_MODE_NON_ISO;
    parameters.bTautFlags = TG_FLAG_FIX_ISO_FIXEDH_BUG |
                            TG_FLAG_FIX_TERM_H_CHRG_BUG;
    bns.num_atoms = 1;
    bns.num_vertices = 1;
    bns.vert = &vertex;
    atoms[0].el_number = 6;
    atoms[0].num_H = 4;
    atoms[0].orig_at_number = 1;
    atoms[0].component = 1;
    strcpy(atoms[0].elname, "C");
    structure.at = atoms;
    structure.num_atoms = 1;
    structure.bMobileH = TAUT_YES;
    structure.iMobileH = TAUT_YES;
    structure.pSrm = &restore_mode;

    ORACLE_NORMALIZE_ALLOCATION_COUNT = 0;
    for (index = 0; index < TAUT_NUM; index++)
    {
        unsigned int shift = (unsigned int) index * 3U;
        if (holder_mask & (1U << shift))
        {
            structure.pOneINChI[index] = calloc(1, sizeof(INChI));
            if (!structure.pOneINChI[index])
            {
                return 70;
            }
            structure.pOneINChI[index]->nErrorCode = 100 + index;
            structure.pOneINChI[index]->nNumberOfAtoms = 10 + index;
            oracle_normalize_register_allocation(structure.pOneINChI[index],
                                                 inchi_labels[index]);
        }
        if (holder_mask & (1U << (shift + 1U)))
        {
            structure.pOneINChI_Aux[index] = calloc(1, sizeof(INChI_Aux));
            if (!structure.pOneINChI_Aux[index])
            {
                return 70;
            }
            structure.pOneINChI_Aux[index]->nErrorCode = 200 + index;
            structure.pOneINChI_Aux[index]->nNumberOfAtoms = 20 + index;
            oracle_normalize_register_allocation(
                structure.pOneINChI_Aux[index], aux_labels[index]);
        }
        if (holder_mask & (1U << (shift + 2U)))
        {
            structure.pOne_norm_data[index] =
                calloc(1, sizeof(INP_ATOM_DATA));
            if (!structure.pOne_norm_data[index])
            {
                return 70;
            }
            structure.pOne_norm_data[index]->num_at = 30 + index;
            structure.pOne_norm_data[index]->num_bonds = 40 + index;
            oracle_normalize_register_allocation(
                structure.pOne_norm_data[index], norm_labels[index]);
        }
    }
    structure.One_ti.t_group = calloc(1, sizeof(T_GROUP));
    if (!structure.One_ti.t_group)
    {
        return 70;
    }
    structure.One_ti.num_t_groups = 1;
    structure.One_ti.t_group[0].nNumEndpoints = 51;
    structure.One_ti.t_group[0].num[0] = 52;
    structure.One_ti.t_group[0].num[1] = 53;
    oracle_normalize_register_allocation(structure.One_ti.t_group,
                                         "One_ti.t_group");

    canonical_globals_before = canonical_globals;
    clock_before = clock;
    parameters_before = parameters;
    data_before = data;
    bns_before = bns;
    bns_data_before = bns_data;
    vertex_before = vertex;
    memcpy(atoms_before, atoms, sizeof(atoms));
    memcpy(atoms2_before, atoms2, sizeof(atoms2));
    memcpy(atoms3_before, atoms3, sizeof(atoms3));
    memcpy(valence_before, valence, sizeof(valence));
    groups_before = groups;
    expected_structure = structure;
    for (index = 0; index < TAUT_NUM; index++)
    {
        expected_structure.pOneINChI[index] = NULL;
        expected_structure.pOneINChI_Aux[index] = NULL;
        expected_structure.pOne_norm_data[index] = NULL;
    }
    memset(&expected_structure.One_ti, 0, sizeof(expected_structure.One_ti));

    ORACLE_NORMALIZE_EVENTS[0] = '\0';
    ORACLE_NORMALIZE_EVENTS_LENGTH = 0;
    ORACLE_NORMALIZE_PREFREE_EXACT = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[0] = 10;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[1] = 11;
    ORACLE_NORMALIZE_FORCED_REBUILD_RETURN = forced_return;
    ORACLE_NORMALIZE_FORCE_REBUILD = 1;
    for (index = 0; index < TAUT_NUM; index++)
    {
        ORACLE_NORMALIZE_INCHI_HOLDERS[index] =
            &structure.pOneINChI[index];
        ORACLE_NORMALIZE_AUX_HOLDERS[index] =
            &structure.pOneINChI_Aux[index];
        ORACLE_NORMALIZE_NORM_HOLDERS[index] =
            structure.pOne_norm_data[index];
    }
    ORACLE_DEFERRED_FREE_COUNT = 0;
    ORACLE_DEFER_FREES = 1;
    ORACLE_NORMALIZE_ACTIVE = 1;
    result = NormalizeAndCompare(
        &canonical_globals, &clock, &parameters, &data, &bns, &bns_data,
        &structure, atoms, atoms2, atoms3, valence, &groups, input, LONG_MAX,
        0, &runs, &delta, 4, 0);
    ORACLE_NORMALIZE_ACTIVE = 0;

    holders_null = structure.pOneINChI[0] == NULL &&
                   structure.pOneINChI[1] == NULL &&
                   structure.pOneINChI_Aux[0] == NULL &&
                   structure.pOneINChI_Aux[1] == NULL &&
                   structure.pOne_norm_data[0] == NULL &&
                   structure.pOne_norm_data[1] == NULL;
    one_ti_cleared = structure.One_ti.t_group == NULL &&
                     structure.One_ti.num_t_groups == 0;
    complete_caller_state_exact =
        memcmp(&canonical_globals, &canonical_globals_before,
               sizeof(canonical_globals)) == 0 &&
        memcmp(&clock, &clock_before, sizeof(clock)) == 0 &&
        memcmp(&parameters, &parameters_before, sizeof(parameters)) == 0 &&
        memcmp(&data, &data_before, sizeof(data)) == 0 &&
        memcmp(&bns, &bns_before, sizeof(bns)) == 0 &&
        memcmp(&bns_data, &bns_data_before, sizeof(bns_data)) == 0 &&
        memcmp(&vertex, &vertex_before, sizeof(vertex)) == 0 &&
        memcmp(atoms, atoms_before, sizeof(atoms)) == 0 &&
        memcmp(atoms2, atoms2_before, sizeof(atoms2)) == 0 &&
        memcmp(atoms3, atoms3_before, sizeof(atoms3)) == 0 &&
        memcmp(valence, valence_before, sizeof(valence)) == 0 &&
        memcmp(&groups, &groups_before, sizeof(groups)) == 0 &&
        memcmp(&structure, &expected_structure, sizeof(structure)) == 0 &&
        input[0] == NULL && input[1] == NULL;
    source_allocations_before = ORACLE_NORMALIZE_ALLOCATION_COUNT;
    source_allocations_after = source_allocations_before >=
                                       ORACLE_DEFERRED_FREE_COUNT
                                   ? source_allocations_before -
                                         ORACLE_DEFERRED_FREE_COUNT
                                   : source_allocations_before;
    allocation_free_exact =
        ORACLE_DEFERRED_FREE_COUNT == source_allocations_before;
    snprintf(case_id, sizeof(case_id), "normalize-initial-%d-mask-%02x",
             forced_return, holder_mask);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s\",\"operation\":\"NormalizeAndCompare\","
           "\"input\":{\"holder_mask\":%u,\"forced_return\":%d},"
           "\"output\":{\"result\":%d,\"runs\":%d,\"delta\":%d,"
           "\"complete_caller_state_exact\":%s,"
           "\"prefree_state_exact\":%s,\"allocation_free_exact\":%s,"
           "\"source_allocations_before\":%zu,"
           "\"source_allocations_after\":%zu,\"holders_null\":%s,"
           "\"one_ti_cleared\":%s,\"input_pointers_preserved\":%s,"
           "\"cleanup_events\":[%s]}}\n",
           case_id, holder_mask, forced_return, result, runs, delta,
           complete_caller_state_exact ? "true" : "false",
           ORACLE_NORMALIZE_PREFREE_EXACT ? "true" : "false",
           allocation_free_exact ? "true" : "false",
           source_allocations_before, source_allocations_after,
           holders_null ? "true" : "false",
           one_ti_cleared ? "true" : "false",
           input[0] == NULL && input[1] == NULL ? "true" : "false",
           ORACLE_NORMALIZE_EVENTS);
    fflush(stdout);

    ORACLE_DEFER_FREES = 0;
    oracle_flush_deferred_frees();
    return 0;
}

static int print_normalize_and_compare_layer_case(int mobile_h,
                                                  int original_state,
                                                  int reversed_state,
                                                  int force_compare)
{
    CANON_GLOBALS canonical_globals;
    struct tagINCHI_CLOCK clock;
    INPUT_PARMS parameters;
    STRUCT_DATA data;
    BN_STRUCT bns;
    BN_DATA bns_data;
    BNS_VERTEX vertex;
    StrFromINChI structure;
    StrFromINChI expected_structure;
    inp_ATOM atoms[1];
    inp_ATOM atoms2[1];
    inp_ATOM atoms3[1];
    VAL_AT valence[1];
    ALL_TC_GROUPS groups;
    SRM restore_mode;
    INChI original[TAUT_NUM];
    INChI *input[TAUT_NUM];
    int result;
    int runs = 17;
    int delta = -19;
    int expected_original_index;
    int expected_reversed_index;
    int holders_null;
    int complete_caller_state_exact;
    char case_id[128];

    memset(&canonical_globals, 0, sizeof(canonical_globals));
    memset(&clock, 0, sizeof(clock));
    memset(&parameters, 0, sizeof(parameters));
    memset(&data, 0, sizeof(data));
    memset(&bns, 0, sizeof(bns));
    memset(&bns_data, 0, sizeof(bns_data));
    memset(&vertex, 0, sizeof(vertex));
    memset(&structure, 0, sizeof(structure));
    memset(atoms, 0, sizeof(atoms));
    memset(atoms2, 0, sizeof(atoms2));
    memset(atoms3, 0, sizeof(atoms3));
    memset(valence, 0, sizeof(valence));
    memset(&groups, 0, sizeof(groups));
    memset(&restore_mode, 0, sizeof(restore_mode));
    memset(original, 0, sizeof(original));

    parameters.nMode = REQ_MODE_TAUT | REQ_MODE_NON_ISO;
    bns.num_atoms = 1;
    bns.num_vertices = 1;
    bns.vert = &vertex;
    atoms[0].el_number = 6;
    atoms[0].num_H = 4;
    atoms[0].orig_at_number = 1;
    atoms[0].component = 1;
    strcpy(atoms[0].elname, "C");
    structure.at = atoms;
    structure.num_atoms = 1;
    structure.bMobileH = (char) mobile_h;
    structure.iMobileH = TAUT_YES;
    structure.pSrm = &restore_mode;
    original[0].nNumberOfAtoms = 1;
    original[0].nErrorCode = 70;
    original[1].nErrorCode = 71;
    input[0] = original;
    input[1] = NULL;
    if (original_state != 0)
    {
        original[1].nNumberOfAtoms = original_state == 1 ? 0 : 1;
        original[1].bDeleted = original_state == 2;
        input[1] = original + 1;
    }

    ORACLE_NORMALIZE_ALLOCATION_COUNT = 0;
    structure.pOneINChI[0] = calloc(1, sizeof(INChI));
    if (!structure.pOneINChI[0])
    {
        return 70;
    }
    structure.pOneINChI[0]->nNumberOfAtoms = 1;
    structure.pOneINChI[0]->nErrorCode = 100;
    oracle_normalize_register_allocation(structure.pOneINChI[0],
                                         "pOneINChI[0]");
    if (reversed_state != 0)
    {
        structure.pOneINChI[1] = calloc(1, sizeof(INChI));
        if (!structure.pOneINChI[1])
        {
            return 70;
        }
        structure.pOneINChI[1]->nNumberOfAtoms =
            reversed_state == 1 ? 0 : 1;
        structure.pOneINChI[1]->bDeleted = reversed_state == 2;
        structure.pOneINChI[1]->nErrorCode = 101;
        oracle_normalize_register_allocation(structure.pOneINChI[1],
                                             "pOneINChI[1]");
    }
    if (!force_compare)
    {
        structure.pOne_norm_data[0] = calloc(1, sizeof(INP_ATOM_DATA));
        if (!structure.pOne_norm_data[0])
        {
            return 70;
        }
        structure.pOne_norm_data[0]->num_at = 30;
        structure.pOne_norm_data[0]->num_bonds = 40;
        oracle_normalize_register_allocation(structure.pOne_norm_data[0],
                                             "pOne_norm_data[0]");
    }

    expected_original_index = mobile_h == TAUT_NON && original_state == 3;
    expected_reversed_index = mobile_h == TAUT_NON && reversed_state == 3;
    expected_structure = structure;
    expected_structure.pOneINChI[0] = NULL;
    expected_structure.pOneINChI[1] = NULL;
    expected_structure.pOne_norm_data[0] = NULL;

    ORACLE_NORMALIZE_EVENTS[0] = '\0';
    ORACLE_NORMALIZE_EVENTS_LENGTH = 0;
    ORACLE_NORMALIZE_PREFREE_EXACT = 1;
    ORACLE_NORMALIZE_FORCED_REBUILD_RETURN = 0;
    ORACLE_NORMALIZE_FORCE_REBUILD = 1;
    ORACLE_NORMALIZE_FORCE_COMPARE = force_compare;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[0] = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[1] =
        reversed_state == 1 ? 0 : 1;
    ORACLE_NORMALIZE_REVERSED[0] = structure.pOneINChI[0];
    ORACLE_NORMALIZE_REVERSED[1] = structure.pOneINChI[1];
    ORACLE_NORMALIZE_ORIGINAL[0] = input[0];
    ORACLE_NORMALIZE_ORIGINAL[1] = input[1];
    ORACLE_NORMALIZE_INCHI_HOLDERS[0] = &structure.pOneINChI[0];
    ORACLE_NORMALIZE_INCHI_HOLDERS[1] = &structure.pOneINChI[1];
    ORACLE_NORMALIZE_AUX_HOLDERS[0] = &structure.pOneINChI_Aux[0];
    ORACLE_NORMALIZE_AUX_HOLDERS[1] = &structure.pOneINChI_Aux[1];
    ORACLE_NORMALIZE_NORM_HOLDERS[0] = structure.pOne_norm_data[0];
    ORACLE_NORMALIZE_NORM_HOLDERS[1] = NULL;
    ORACLE_DEFERRED_FREE_COUNT = 0;
    ORACLE_DEFER_FREES = 1;
    ORACLE_NORMALIZE_ACTIVE = 1;
    result = NormalizeAndCompare(
        &canonical_globals, &clock, &parameters, &data, &bns, &bns_data,
        &structure, atoms, atoms2, atoms3, valence, &groups, input, LONG_MAX,
        0, &runs, &delta, 4, 0);
    ORACLE_NORMALIZE_ACTIVE = 0;
    ORACLE_NORMALIZE_FORCE_COMPARE = 0;

    holders_null = !structure.pOneINChI[0] && !structure.pOneINChI[1] &&
                   !structure.pOneINChI_Aux[0] &&
                   !structure.pOneINChI_Aux[1] &&
                   !structure.pOne_norm_data[0] &&
                   !structure.pOne_norm_data[1];
    complete_caller_state_exact =
        memcmp(&structure, &expected_structure, sizeof(structure)) == 0 &&
        input[0] == original &&
        input[1] == (original_state ? original + 1 : NULL) &&
        original[0].nErrorCode == 70 && original[1].nErrorCode == 71;
    snprintf(case_id, sizeof(case_id),
             force_compare ? "normalize-layer-m%d-o%d-r%d"
                           : "normalize-common-success",
             mobile_h, original_state, reversed_state);
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"case_id\":\"%s\",\"operation\":\"NormalizeAndCompare\","
           "\"family\":\"%s\","
           "\"input\":{\"mobile_h\":%d,\"original_state\":%d,"
           "\"reversed_state\":%d},"
           "\"output\":{\"result\":%d,\"runs\":%d,\"delta\":%d,"
           "\"selected_original_index\":%d,"
           "\"selected_reversed_index\":%d,"
           "\"complete_caller_state_exact\":%s,"
           "\"prefree_state_exact\":%s,\"allocation_free_exact\":%s,"
           "\"source_allocations_before\":%zu,"
           "\"source_allocations_after\":0,\"holders_null\":%s,"
           "\"one_ti_cleared\":true,\"input_pointers_preserved\":true,"
           "\"cleanup_events\":[%s]}}\n",
           case_id, force_compare ? "layer-selection" : "common-success",
           mobile_h, original_state, reversed_state, result, runs,
           delta, expected_original_index, expected_reversed_index,
           complete_caller_state_exact ? "true" : "false",
           ORACLE_NORMALIZE_PREFREE_EXACT ? "true" : "false",
           ORACLE_DEFERRED_FREE_COUNT == ORACLE_NORMALIZE_ALLOCATION_COUNT
               ? "true"
               : "false",
           ORACLE_NORMALIZE_ALLOCATION_COUNT,
           holders_null ? "true" : "false", ORACLE_NORMALIZE_EVENTS);
    fflush(stdout);
    ORACLE_DEFER_FREES = 0;
    oracle_flush_deferred_frees();
    return 0;
}

static int run_normalize_and_compare_zz_case(int n_zy, int n_pzz,
                                             const char *formula,
                                             int failure_stage,
                                             int force_growth,
                                             int emit_record)
{
    CANON_GLOBALS canonical_globals;
    CANON_GLOBALS canonical_globals_before;
    struct tagINCHI_CLOCK clock;
    struct tagINCHI_CLOCK clock_before;
    INPUT_PARMS parameters;
    INPUT_PARMS parameters_before;
    STRUCT_DATA data;
    STRUCT_DATA data_before;
    BN_STRUCT bns;
    BN_STRUCT bns_before;
    BN_DATA bns_data;
    BN_DATA bns_data_before;
    BNS_VERTEX vertex;
    BNS_VERTEX vertex_before;
    StrFromINChI structure;
    StrFromINChI expected_structure;
    inp_ATOM atoms[1];
    inp_ATOM atoms_before[1];
    inp_ATOM atoms2[1];
    inp_ATOM atoms2_before[1];
    inp_ATOM atoms3[1];
    inp_ATOM atoms3_before[1];
    VAL_AT valence[1];
    VAL_AT valence_before[1];
    ALL_TC_GROUPS groups;
    ALL_TC_GROUPS groups_before;
    SRM restore_mode;
    INChI original;
    INChI *input[TAUT_NUM] = {&original, NULL};
    int result;
    int runs = 17;
    int delta = -19;
    int complete_caller_state_exact;
    int holders_null;
    char case_id[160];

    memset(&canonical_globals, 0, sizeof(canonical_globals));
    memset(&clock, 0, sizeof(clock));
    memset(&parameters, 0, sizeof(parameters));
    memset(&data, 0, sizeof(data));
    memset(&bns, 0, sizeof(bns));
    memset(&bns_data, 0, sizeof(bns_data));
    memset(&vertex, 0, sizeof(vertex));
    memset(&structure, 0, sizeof(structure));
    memset(atoms, 0, sizeof(atoms));
    memset(atoms2, 0, sizeof(atoms2));
    memset(atoms3, 0, sizeof(atoms3));
    memset(valence, 0, sizeof(valence));
    memset(&groups, 0, sizeof(groups));
    memset(&restore_mode, 0, sizeof(restore_mode));
    memset(&original, 0, sizeof(original));

    parameters.nMode = REQ_MODE_TAUT | REQ_MODE_NON_ISO;
    bns.num_atoms = 1;
    bns.num_vertices = 1;
    bns.vert = &vertex;
    atoms[0].el_number = 6;
    atoms[0].num_H = 4;
    atoms[0].orig_at_number = 1;
    atoms[0].component = 1;
    strcpy(atoms[0].elname, "C");
    structure.at = atoms;
    structure.num_atoms = 1;
    structure.bMobileH = TAUT_YES;
    structure.iMobileH = TAUT_YES;
    structure.pSrm = &restore_mode;
    structure.n_zy = n_zy;
    structure.n_pzz = n_pzz;
    original.nNumberOfAtoms = 1;
    original.nErrorCode = 70;

    ORACLE_NORMALIZE_ALLOCATION_COUNT = 0;
    structure.pOneINChI[0] = calloc(1, sizeof(INChI));
    if (!structure.pOneINChI[0])
    {
        return 70;
    }
    structure.pOneINChI[0]->nNumberOfAtoms = 1;
    structure.pOneINChI[0]->nErrorCode = 100;
    oracle_normalize_register_allocation(structure.pOneINChI[0],
                                         "pOneINChI[0]");
    if (formula)
    {
        size_t formula_size = strlen(formula) + 1;
        structure.pOneINChI[0]->szHillFormula = calloc(formula_size, 1);
        if (!structure.pOneINChI[0]->szHillFormula)
        {
            return 70;
        }
        memcpy(structure.pOneINChI[0]->szHillFormula, formula, formula_size);
        oracle_normalize_register_allocation(
            structure.pOneINChI[0]->szHillFormula, "formula");
    }

    canonical_globals_before = canonical_globals;
    clock_before = clock;
    parameters_before = parameters;
    data_before = data;
    bns_before = bns;
    bns_data_before = bns_data;
    vertex_before = vertex;
    memcpy(atoms_before, atoms, sizeof(atoms));
    memcpy(atoms2_before, atoms2, sizeof(atoms2));
    memcpy(atoms3_before, atoms3, sizeof(atoms3));
    memcpy(valence_before, valence, sizeof(valence));
    groups_before = groups;
    expected_structure = structure;
    expected_structure.pOneINChI[0] = NULL;

    ORACLE_NORMALIZE_EVENTS[0] = '\0';
    ORACLE_NORMALIZE_EVENTS_LENGTH = 0;
    ORACLE_NORMALIZE_PREFREE_EXACT = 1;
    ORACLE_NORMALIZE_ZZ_FORMULA[0] = '\0';
    ORACLE_NORMALIZE_FORCED_REBUILD_RETURN = 0;
    ORACLE_NORMALIZE_FORCE_REBUILD = 1;
    ORACLE_NORMALIZE_FORCE_COMPARE = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[0] = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[1] = 0;
    ORACLE_NORMALIZE_REVERSED[0] = structure.pOneINChI[0];
    ORACLE_NORMALIZE_REVERSED[1] = NULL;
    ORACLE_NORMALIZE_ORIGINAL[0] = input[0];
    ORACLE_NORMALIZE_ORIGINAL[1] = NULL;
    ORACLE_NORMALIZE_INCHI_HOLDERS[0] = &structure.pOneINChI[0];
    ORACLE_NORMALIZE_INCHI_HOLDERS[1] = &structure.pOneINChI[1];
    ORACLE_NORMALIZE_AUX_HOLDERS[0] = &structure.pOneINChI_Aux[0];
    ORACLE_NORMALIZE_AUX_HOLDERS[1] = &structure.pOneINChI_Aux[1];
    ORACLE_NORMALIZE_NORM_HOLDERS[0] = NULL;
    ORACLE_NORMALIZE_NORM_HOLDERS[1] = NULL;
    ORACLE_NORMALIZE_ZZ_FAILURE_STAGE = failure_stage;
    ORACLE_NORMALIZE_ZZ_CALLOC_ORDINAL = 0;
    ORACLE_NORMALIZE_ZZ_FORCE_GROWTH = force_growth;
    ORACLE_NORMALIZE_ZZ_SUCCESSFUL_ALLOCATIONS = formula ? 2 : 1;
    ORACLE_NORMALIZE_ZZ_FREES = 0;
    ORACLE_NORMALIZE_ZZ_REALLOC_CALLS = 0;
    ORACLE_NORMALIZE_ZZ_ACTIVE = 1;
    ORACLE_NORMALIZE_ACTIVE = 1;
    result = NormalizeAndCompare(
        &canonical_globals, &clock, &parameters, &data, &bns, &bns_data,
        &structure, atoms, atoms2, atoms3, valence, &groups, input, LONG_MAX,
        0, &runs, &delta, 4, 0);
    ORACLE_NORMALIZE_ACTIVE = 0;
    ORACLE_NORMALIZE_ZZ_ACTIVE = 0;
    ORACLE_NORMALIZE_FORCE_COMPARE = 0;

    holders_null = !structure.pOneINChI[0] && !structure.pOneINChI[1] &&
                   !structure.pOneINChI_Aux[0] &&
                   !structure.pOneINChI_Aux[1] &&
                   !structure.pOne_norm_data[0] &&
                   !structure.pOne_norm_data[1];
    complete_caller_state_exact =
        memcmp(&canonical_globals, &canonical_globals_before,
               sizeof(canonical_globals)) == 0 &&
        memcmp(&clock, &clock_before, sizeof(clock)) == 0 &&
        memcmp(&parameters, &parameters_before, sizeof(parameters)) == 0 &&
        memcmp(&data, &data_before, sizeof(data)) == 0 &&
        memcmp(&bns, &bns_before, sizeof(bns)) == 0 &&
        memcmp(&bns_data, &bns_data_before, sizeof(bns_data)) == 0 &&
        memcmp(&vertex, &vertex_before, sizeof(vertex)) == 0 &&
        memcmp(atoms, atoms_before, sizeof(atoms)) == 0 &&
        memcmp(atoms2, atoms2_before, sizeof(atoms2)) == 0 &&
        memcmp(atoms3, atoms3_before, sizeof(atoms3)) == 0 &&
        memcmp(valence, valence_before, sizeof(valence)) == 0 &&
        memcmp(&groups, &groups_before, sizeof(groups)) == 0 &&
        memcmp(&structure, &expected_structure, sizeof(structure)) == 0 &&
        input[0] == &original && !input[1] && original.nErrorCode == 70;

    if (emit_record)
    {
        snprintf(case_id, sizeof(case_id),
                 "normalize-zz-nzy%d-npzz%d-f%s-fail%d-grow%d", n_zy,
                 n_pzz, formula ? formula : "null", failure_stage,
                 force_growth);
        printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
               "\"case_id\":\"%s\",\"operation\":\"NormalizeAndCompare\"," 
               "\"family\":\"zz-zy\",\"input\":{\"n_zy\":%d," 
               "\"n_pzz\":%d,\"formula_present\":%s," 
               "\"formula\":\"%s\"," 
               "\"failure_stage\":%d,\"force_growth\":%s}," 
               "\"output\":{\"result\":%d,\"runs\":%d,\"delta\":%d," 
               "\"formula_before_cleanup\":\"%s\"," 
               "\"complete_caller_state_exact\":%s," 
               "\"prefree_state_exact\":%s,\"holders_null\":%s," 
               "\"source_allocation_calls\":%d," 
               "\"successful_allocations\":%d,\"frees\":%d," 
               "\"allocation_free_exact\":%s,\"cleanup_events\":[%s]}}\n",
               case_id, n_zy, n_pzz, formula ? "true" : "false",
               formula ? formula : "", failure_stage,
               force_growth ? "true" : "false", result, runs, delta,
               ORACLE_NORMALIZE_ZZ_FORMULA,
               complete_caller_state_exact ? "true" : "false",
               ORACLE_NORMALIZE_PREFREE_EXACT ? "true" : "false",
               holders_null ? "true" : "false",
               ORACLE_NORMALIZE_ZZ_CALLOC_ORDINAL +
                   ORACLE_NORMALIZE_ZZ_REALLOC_CALLS,
               ORACLE_NORMALIZE_ZZ_SUCCESSFUL_ALLOCATIONS,
               ORACLE_NORMALIZE_ZZ_FREES,
               ORACLE_NORMALIZE_ZZ_SUCCESSFUL_ALLOCATIONS ==
                       ORACLE_NORMALIZE_ZZ_FREES
                   ? "true"
                   : "false",
               ORACLE_NORMALIZE_EVENTS);
        fflush(stdout);
    }
    return 0;
}

static int print_normalize_and_compare_zz_undefined_case(void)
{
    pid_t child;
    int status;
    fflush(stdout);
    child = fork();
    if (child < 0)
    {
        return 71;
    }
    if (child == 0)
    {
        (void) run_normalize_and_compare_zz_case(1, 1, "CZzZz2", 1, 0,
                                                  0);
        _exit(0);
    }
    if (waitpid(child, &status, 0) != child)
    {
        return 72;
    }
    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"case_id\":\"normalize-zz-initial-buffer-allocation-undefined\"," 
           "\"operation\":\"NormalizeAndCompare\"," 
           "\"family\":\"zz-zy-undefined\"," 
           "\"input\":{\"n_zy\":1,\"n_pzz\":1," 
           "\"formula\":\"CZzZz2\",\"failure_stage\":1}," 
           "\"output\":{\"source_defined\":false," 
           "\"termination_kind\":\"%s\",\"termination_value\":%d}}\n",
           WIFSIGNALED(status) ? "signal" : "exit",
           WIFSIGNALED(status) ? WTERMSIG(status) : WEXITSTATUS(status));
    fflush(stdout);
    return 0;
}

typedef struct tagOracleNormalizeFirstLessCase
{
    const char *name;
    int mobile_h;
    int rebuild_len;
    int rebuild_ret[3];
    int rebuild_prep[3];
    int compare_len;
    INCHI_MODE compare_flags[4];
    int compare_error[4];
    int compare_h1[4];
    int compare_h2[4];
    int fill_len;
    int fill_ret[2];
    int fix_len;
    int fix_ret[2];
    int repair_kind;
    int rebuild_norm[3];
    int compare_endpoints[4];
} ORACLE_NORMALIZE_FIRST_LESS_CASE;

static const ORACLE_NORMALIZE_FIRST_LESS_CASE
    ORACLE_NORMALIZE_FIRST_LESS_CASES[] = {
        {"fill-skip", TAUT_YES, 1, {0}, {0}, 1, {IDIF_PROBLEM}},
        {"fill-negative", TAUT_NON, 1, {0}, {0}, 0, {0}, {0}, {0}, {0},
         1, {-7}},
        {"fill-positive", TAUT_NON, 1, {0}, {0}, 0, {0}, {0}, {0}, {0},
         1, {7}},
        {"compare-problem-precedence", TAUT_NON, 1, {0}, {0}, 1,
         {IDIF_PROBLEM}, {9}, {0}, {0}, 1, {0}},
        {"compare-error", TAUT_NON, 1, {0}, {0}, 1, {0}, {9}, {0}, {0},
         1, {0}},
        {"less-flag-false", TAUT_YES, 1, {0}, {1}, 2,
         {0, IDIF_PROBLEM}},
        {"less-prep-null", TAUT_YES, 1, {0}, {0}, 2,
         {IDIF_LESS_H, IDIF_PROBLEM}, {0}, {0}, {2}},
        {"less-delta-zero", TAUT_YES, 1, {0}, {1}, 2,
         {IDIF_LESS_H, IDIF_PROBLEM}, {0}, {2}, {2}},
        {"less-delta-negative", TAUT_YES, 1, {0}, {1}, 2,
         {IDIF_LESS_H, IDIF_PROBLEM}, {0}, {3}, {2}},
        {"less-fixer-negative", TAUT_YES, 1, {0}, {1}, 1,
         {IDIF_LESS_H}, {0}, {0}, {3}, 0, {0}, 1, {-7}},
        {"less-fixer-zero", TAUT_YES, 1, {0}, {1}, 2,
         {IDIF_LESS_H, IDIF_PROBLEM}, {0}, {0}, {3}, 0, {0}, 1, {0}},
        {"less-rebuild-negative", TAUT_YES, 2, {0, -8}, {1, 1}, 1,
         {IDIF_LESS_H}, {0}, {0}, {3}, 0, {0}, 1, {1}},
        {"less-repeat-fill-positive", TAUT_NON, 2, {0, 0}, {1, 1}, 1,
         {IDIF_LESS_H}, {0}, {0}, {3}, 2, {0, 7}, 1, {1}},
        {"less-exit-cleared", TAUT_YES, 2, {0, 0}, {1, 1}, 3,
         {IDIF_LESS_H, 0, IDIF_PROBLEM}, {0}, {0, 0, 0}, {3, 2, 0},
         0, {0}, 1, {1}},
        {"less-exit-prep-null", TAUT_YES, 2, {0, 0}, {1, 0}, 3,
         {IDIF_LESS_H, IDIF_LESS_H, IDIF_PROBLEM}, {0}, {0, 0, 0},
         {3, 2, 0}, 0, {0}, 1, {1}},
        {"less-exit-zero-delta", TAUT_YES, 2, {0, 0}, {1, 1}, 3,
         {IDIF_LESS_H, IDIF_LESS_H, IDIF_PROBLEM}, {0}, {0, 2, 0},
         {3, 2, 0}, 0, {0}, 1, {1}},
        {"less-exit-nondecreasing", TAUT_YES, 2, {0, 0}, {1, 1}, 3,
         {IDIF_LESS_H, IDIF_LESS_H, IDIF_PROBLEM}, {0}, {0, 0, 0},
         {3, 3, 0}, 0, {0}, 1, {1}},
        {"less-continue-decreasing", TAUT_YES, 2, {0, 0}, {1, 1}, 3,
        {IDIF_LESS_H, IDIF_LESS_H, IDIF_PROBLEM}, {0}, {0, 0, 0},
         {3, 2, 0}, 0, {0}, 2, {1, 0}}};

enum
{
    ORACLE_NORMALIZE_REPAIR_LESS = 0,
    ORACLE_NORMALIZE_REPAIR_MORE = 1,
    ORACLE_NORMALIZE_REPAIR_EXTRA = 2
};

static const ORACLE_NORMALIZE_FIRST_LESS_CASE
    ORACLE_NORMALIZE_MORE_EXTRA_CASES[] = {
        {.name = "more-flag-false", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .compare_len = 2, .compare_flags = {0, IDIF_PROBLEM},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-prep-null", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {0},
         .compare_len = 2,
         .compare_flags = {IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {2}, .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-delta-zero", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .compare_len = 2,
         .compare_flags = {IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {2}, .compare_h2 = {2},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-delta-negative", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .compare_len = 2,
         .compare_flags = {IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {2}, .compare_h2 = {3},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-fixer-negative", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .compare_len = 1, .compare_flags = {IDIF_MORE_H},
         .compare_h1 = {3}, .fix_len = 1, .fix_ret = {-7},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-fixer-zero", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .compare_len = 2,
         .compare_flags = {IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {3}, .fix_len = 1, .fix_ret = {0},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-rebuild-negative", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, -8},
         .rebuild_prep = {1, 1}, .compare_len = 1,
         .compare_flags = {IDIF_MORE_H}, .compare_h1 = {3},
         .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-repeat-fill-positive", .mobile_h = TAUT_NON,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .compare_len = 1,
         .compare_flags = {IDIF_MORE_H}, .compare_h1 = {3},
         .fill_len = 2, .fill_ret = {0, 7},
         .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-exit-cleared", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .compare_len = 3,
         .compare_flags = {IDIF_MORE_H, 0, IDIF_PROBLEM},
         .compare_h1 = {3, 2}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-exit-prep-null", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 0}, .compare_len = 3,
         .compare_flags = {IDIF_MORE_H, IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {3, 2}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-exit-zero-delta", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .compare_len = 3,
         .compare_flags = {IDIF_MORE_H, IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {3, 2}, .compare_h2 = {0, 2},
         .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-exit-nondecreasing", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .compare_len = 3,
         .compare_flags = {IDIF_MORE_H, IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {3, 3}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "more-continue-decreasing", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .compare_len = 3,
         .compare_flags = {IDIF_MORE_H, IDIF_MORE_H, IDIF_PROBLEM},
         .compare_h1 = {3, 2}, .fix_len = 2, .fix_ret = {1, 0},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_MORE},
        {.name = "extra-flag-false", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .rebuild_norm = {1}, .compare_len = 2,
         .compare_flags = {0, IDIF_PROBLEM},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-norm-null", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .compare_len = 2,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_PROBLEM},
         .compare_endpoints = {2},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-delta-zero", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .rebuild_norm = {1}, .compare_len = 2,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_PROBLEM},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-delta-negative", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .rebuild_norm = {1}, .compare_len = 2,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_PROBLEM},
         .compare_endpoints = {-1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-fixer-negative", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .rebuild_norm = {1}, .compare_len = 1,
         .compare_flags = {IDIF_EXTRA_TG_ENDP},
         .compare_endpoints = {3}, .fix_len = 1, .fix_ret = {-7},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-fixer-zero", .mobile_h = TAUT_YES,
         .rebuild_len = 1, .rebuild_ret = {0}, .rebuild_prep = {1},
         .rebuild_norm = {1}, .compare_len = 2,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_PROBLEM},
         .compare_endpoints = {3}, .fix_len = 1, .fix_ret = {0},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-rebuild-negative", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, -8},
         .rebuild_prep = {1, 1}, .rebuild_norm = {1, 1},
         .compare_len = 1, .compare_flags = {IDIF_EXTRA_TG_ENDP},
         .compare_endpoints = {3}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-repeat-fill-positive", .mobile_h = TAUT_NON,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .rebuild_norm = {1, 1},
         .compare_len = 1, .compare_flags = {IDIF_EXTRA_TG_ENDP},
         .compare_endpoints = {3}, .fill_len = 2, .fill_ret = {0, 7},
         .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-exit-cleared", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .rebuild_norm = {1, 1},
         .compare_len = 3,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, 0, IDIF_PROBLEM},
         .compare_endpoints = {3, 2}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-exit-norm-null", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .rebuild_norm = {1, 0},
         .compare_len = 3,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_EXTRA_TG_ENDP,
                           IDIF_PROBLEM},
         .compare_endpoints = {3, 2}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-exit-zero-delta", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .rebuild_norm = {1, 1},
         .compare_len = 3,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_EXTRA_TG_ENDP,
                           IDIF_PROBLEM},
         .compare_endpoints = {3, 0}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-exit-nondecreasing", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .rebuild_norm = {1, 1},
         .compare_len = 3,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_EXTRA_TG_ENDP,
                           IDIF_PROBLEM},
         .compare_endpoints = {3, 3}, .fix_len = 1, .fix_ret = {1},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA},
        {.name = "extra-continue-decreasing", .mobile_h = TAUT_YES,
         .rebuild_len = 2, .rebuild_ret = {0, 0},
         .rebuild_prep = {1, 1}, .rebuild_norm = {1, 1},
         .compare_len = 3,
         .compare_flags = {IDIF_EXTRA_TG_ENDP, IDIF_EXTRA_TG_ENDP,
                           IDIF_PROBLEM},
         .compare_endpoints = {3, 2}, .fix_len = 2, .fix_ret = {1, 0},
         .repair_kind = ORACLE_NORMALIZE_REPAIR_EXTRA}};

static void print_normalize_int_array(const int *values, int length)
{
    int index;
    putchar('[');
    for (index = 0; index < length; index++)
    {
        if (index)
        {
            putchar(',');
        }
        printf("%d", values[index]);
    }
    putchar(']');
}

static void print_normalize_mode_array(const INCHI_MODE *values, int length)
{
    int index;
    putchar('[');
    for (index = 0; index < length; index++)
    {
        if (index)
        {
            putchar(',');
        }
        printf("%" PRIu64, (uint64_t) values[index]);
    }
    putchar(']');
}

static int print_normalize_and_compare_first_less_case(
    const ORACLE_NORMALIZE_FIRST_LESS_CASE *test_case)
{
    CANON_GLOBALS canonical_globals;
    struct tagINCHI_CLOCK clock;
    INPUT_PARMS parameters;
    STRUCT_DATA data;
    BN_STRUCT bns;
    BN_DATA bns_data;
    BNS_VERTEX vertex;
    StrFromINChI structure;
    StrFromINChI expected_structure;
    inp_ATOM atoms[1], atoms2[1], atoms3[1];
    VAL_AT valence[1];
    ALL_TC_GROUPS groups;
    SRM restore_mode;
    INChI original;
    INChI *input[TAUT_NUM] = {&original, NULL};
    int result, runs = 17, delta = -19, index;
    int holders_null, complete_caller_state_exact;

    memset(&canonical_globals, 0, sizeof(canonical_globals));
    memset(&clock, 0, sizeof(clock));
    memset(&parameters, 0, sizeof(parameters));
    memset(&data, 0, sizeof(data));
    memset(&bns, 0, sizeof(bns));
    memset(&bns_data, 0, sizeof(bns_data));
    memset(&vertex, 0, sizeof(vertex));
    memset(&structure, 0, sizeof(structure));
    memset(atoms, 0, sizeof(atoms));
    memset(atoms2, 0, sizeof(atoms2));
    memset(atoms3, 0, sizeof(atoms3));
    memset(valence, 0, sizeof(valence));
    memset(&groups, 0, sizeof(groups));
    memset(&restore_mode, 0, sizeof(restore_mode));
    memset(&original, 0, sizeof(original));
    parameters.nMode = REQ_MODE_TAUT | REQ_MODE_NON_ISO;
    bns.num_atoms = bns.num_vertices = 1;
    bns.vert = &vertex;
    atoms[0].el_number = 6;
    atoms[0].num_H = 4;
    atoms[0].orig_at_number = atoms[0].component = 1;
    strcpy(atoms[0].elname, "C");
    structure.at = atoms;
    structure.num_atoms = 1;
    structure.bMobileH = TAUT_INI;
    structure.iMobileH = (char) test_case->mobile_h;
    structure.pSrm = &restore_mode;
    original.nNumberOfAtoms = 1;
    original.nErrorCode = 70;
    ORACLE_NORMALIZE_ALLOCATION_COUNT = 0;
    structure.pOneINChI[0] = calloc(1, sizeof(INChI));
    if (!structure.pOneINChI[0])
    {
        return 70;
    }
    structure.pOneINChI[0]->nNumberOfAtoms = 1;
    structure.pOneINChI[0]->nErrorCode = 100;
    oracle_normalize_register_allocation(structure.pOneINChI[0],
                                         "pOneINChI[0]");
    expected_structure = structure;
    expected_structure.pOneINChI[0] = NULL;

    ORACLE_NORMALIZE_REBUILD_LENGTH = test_case->rebuild_len;
    ORACLE_NORMALIZE_COMPARE_LENGTH = test_case->compare_len;
    ORACLE_NORMALIZE_FILL_LENGTH = test_case->fill_len;
    ORACLE_NORMALIZE_FIX_LESS_LENGTH =
        test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_LESS
            ? test_case->fix_len : 0;
    ORACLE_NORMALIZE_FIX_MORE_LENGTH =
        test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_MORE
            ? test_case->fix_len : 0;
    ORACLE_NORMALIZE_FIX_EXTRA_LENGTH =
        test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_EXTRA
            ? test_case->fix_len : 0;
    for (index = 0; index < test_case->rebuild_len; index++)
    {
        ORACLE_NORMALIZE_REBUILD_RETURNS[index] = test_case->rebuild_ret[index];
        ORACLE_NORMALIZE_REBUILD_PREP[index] = test_case->rebuild_prep[index];
        ORACLE_NORMALIZE_REBUILD_NORM[index] = test_case->rebuild_norm[index];
    }
    for (index = 0; index < test_case->compare_len; index++)
    {
        ORACLE_NORMALIZE_COMPARE_FLAGS[index] = test_case->compare_flags[index];
        ORACLE_NORMALIZE_COMPARE_ERRORS[index] = test_case->compare_error[index];
        ORACLE_NORMALIZE_COMPARE_H1[index] = test_case->compare_h1[index];
        ORACLE_NORMALIZE_COMPARE_H2[index] = test_case->compare_h2[index];
        ORACLE_NORMALIZE_COMPARE_ENDPOINTS[index] =
            test_case->compare_endpoints[index];
    }
    memcpy(ORACLE_NORMALIZE_FILL_RETURNS, test_case->fill_ret,
           sizeof(test_case->fill_ret));
    memcpy(ORACLE_NORMALIZE_FIX_LESS_RETURNS, test_case->fix_ret,
           sizeof(test_case->fix_ret));
    memcpy(ORACLE_NORMALIZE_FIX_MORE_RETURNS, test_case->fix_ret,
           sizeof(test_case->fix_ret));
    memcpy(ORACLE_NORMALIZE_FIX_EXTRA_RETURNS, test_case->fix_ret,
           sizeof(test_case->fix_ret));
    ORACLE_NORMALIZE_REBUILD_POSITION = ORACLE_NORMALIZE_COMPARE_POSITION = 0;
    ORACLE_NORMALIZE_FILL_POSITION = ORACLE_NORMALIZE_FIX_LESS_POSITION = 0;
    ORACLE_NORMALIZE_FIX_MORE_POSITION = 0;
    ORACLE_NORMALIZE_FIX_EXTRA_POSITION = 0;
    ORACLE_NORMALIZE_EVENTS[0] = '\0';
    ORACLE_NORMALIZE_EVENTS_LENGTH = 0;
    ORACLE_NORMALIZE_PREFREE_EXACT = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[0] = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[1] = 0;
    ORACLE_NORMALIZE_REVERSED[0] = structure.pOneINChI[0];
    ORACLE_NORMALIZE_REVERSED[1] = NULL;
    ORACLE_NORMALIZE_ORIGINAL[0] = input[0];
    ORACLE_NORMALIZE_ORIGINAL[1] = NULL;
    ORACLE_NORMALIZE_INCHI_HOLDERS[0] = &structure.pOneINChI[0];
    ORACLE_NORMALIZE_INCHI_HOLDERS[1] = &structure.pOneINChI[1];
    ORACLE_NORMALIZE_AUX_HOLDERS[0] = &structure.pOneINChI_Aux[0];
    ORACLE_NORMALIZE_AUX_HOLDERS[1] = &structure.pOneINChI_Aux[1];
    ORACLE_NORMALIZE_NORM_HOLDERS[0] = ORACLE_NORMALIZE_NORM_HOLDERS[1] = NULL;
    ORACLE_DEFERRED_FREE_COUNT = 0;
    ORACLE_DEFER_FREES = 1;
    ORACLE_NORMALIZE_FIRST_LESS_ACTIVE = ORACLE_NORMALIZE_ACTIVE = 1;
    result = NormalizeAndCompare(
        &canonical_globals, &clock, &parameters, &data, &bns, &bns_data,
        &structure, atoms, atoms2, atoms3, valence, &groups, input, LONG_MAX,
        0, &runs, &delta, 4, 0);
    ORACLE_NORMALIZE_ACTIVE = ORACLE_NORMALIZE_FIRST_LESS_ACTIVE = 0;
    holders_null = !structure.pOneINChI[0] && !structure.pOneINChI[1] &&
                   !structure.pOneINChI_Aux[0] && !structure.pOneINChI_Aux[1] &&
                   !structure.pOne_norm_data[0] && !structure.pOne_norm_data[1];
    complete_caller_state_exact =
        memcmp(&structure, &expected_structure, sizeof(structure)) == 0 &&
        input[0] == &original && !input[1] && original.nErrorCode == 70;

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"case_id\":\"normalize-%s-%s\"," 
           "\"operation\":\"NormalizeAndCompare\"," 
           "\"family\":\"%s\",\"input\":{\"repair_kind\":\"%s\"," 
           "\"mobile_h\":%d,\"rebuild_ret\":",
           test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_LESS
               ? "first-less" : "more-extra",
           test_case->name,
           test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_LESS
               ? "first-less" : "more-extra",
           test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_LESS
               ? "less" :
               (test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_MORE
                    ? "more" : "extra"),
           test_case->mobile_h);
    print_normalize_int_array(test_case->rebuild_ret, test_case->rebuild_len);
    printf(",\"rebuild_prep\":");
    print_normalize_int_array(test_case->rebuild_prep, test_case->rebuild_len);
    printf(",\"rebuild_norm\":");
    print_normalize_int_array(test_case->rebuild_norm, test_case->rebuild_len);
    printf(",\"compare_flags\":");
    print_normalize_mode_array(test_case->compare_flags,
                               test_case->compare_len);
    printf(",\"compare_error\":");
    print_normalize_int_array(test_case->compare_error, test_case->compare_len);
    printf(",\"compare_h1\":");
    print_normalize_int_array(test_case->compare_h1, test_case->compare_len);
    printf(",\"compare_h2\":");
    print_normalize_int_array(test_case->compare_h2, test_case->compare_len);
    printf(",\"compare_endpoints\":");
    print_normalize_int_array(test_case->compare_endpoints,
                              test_case->compare_len);
    printf(",\"fill_ret\":");
    print_normalize_int_array(test_case->fill_ret, test_case->fill_len);
    printf(",\"fix_ret\":");
    print_normalize_int_array(test_case->fix_ret, test_case->fix_len);
    printf("},\"output\":{\"result\":%d,\"runs\":%d,\"delta\":%d," 
           "\"rebuild_used\":%d,\"compare_used\":%d," 
           "\"fill_used\":%d,\"fix_used\":%d," 
           "\"complete_caller_state_exact\":%s," 
           "\"prefree_state_exact\":%s,\"holders_null\":%s," 
           "\"allocation_free_exact\":%s,\"cleanup_events\":[%s]}}\n",
           result, runs, delta, ORACLE_NORMALIZE_REBUILD_POSITION,
           ORACLE_NORMALIZE_COMPARE_POSITION, ORACLE_NORMALIZE_FILL_POSITION,
           test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_LESS
               ? ORACLE_NORMALIZE_FIX_LESS_POSITION
               : (test_case->repair_kind == ORACLE_NORMALIZE_REPAIR_MORE
                      ? ORACLE_NORMALIZE_FIX_MORE_POSITION
                      : ORACLE_NORMALIZE_FIX_EXTRA_POSITION),
           complete_caller_state_exact ? "true" : "false",
           ORACLE_NORMALIZE_PREFREE_EXACT ? "true" : "false",
           holders_null ? "true" : "false",
           ORACLE_DEFERRED_FREE_COUNT == ORACLE_NORMALIZE_ALLOCATION_COUNT
               ? "true" : "false",
           ORACLE_NORMALIZE_EVENTS);
    fflush(stdout);
    ORACLE_DEFER_FREES = 0;
    oracle_flush_deferred_frees();
    return 0;
}

typedef struct tagOracleNormalizeFinalCase
{
    const char *name;
    int b_mobile_h;
    int i_mobile_h;
    int original_layer1;
    int reversed_layer1;
    long num_inp;
    int fixed_len;
    int fixed_ret[4];
    int mobile_len;
    int mobile_ret[1];
    INCHI_MODE primary_flags;
    int primary_error;
    INCHI_MODE optional_flags;
    int optional_error;
    int stereo_len;
    int stereo_ret[1];
    int tc_len;
    int tc_edges[2];
    int tinfo_present;
    int tinfo_len;
    int tinfo_h[2];
    int tinfo_endpoints[2];
} ORACLE_NORMALIZE_FINAL_CASE;

static const ORACLE_NORMALIZE_FINAL_CASE ORACLE_NORMALIZE_FINAL_CASES[] = {
    {.name = "fixed-negative-long-min", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .num_inp = LONG_MIN,
     .fixed_len = 1, .fixed_ret = {-7}},
    {.name = "fixed-zero", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .fixed_len = 1, .fixed_ret = {0},
     .primary_flags = IDIF_PROBLEM},
    {.name = "fixed-one-positive", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .fixed_len = 2, .fixed_ret = {1, 0},
     .primary_flags = IDIF_PROBLEM},
    {.name = "fixed-two-positive", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .fixed_len = 3, .fixed_ret = {1, 1, 0},
     .primary_flags = IDIF_PROBLEM},
    {.name = "fixed-three-positive", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .fixed_len = 3, .fixed_ret = {1, 1, 1},
     .primary_flags = IDIF_PROBLEM},
    {.name = "mobile-negative-long-max", .b_mobile_h = TAUT_YES,
     .i_mobile_h = TAUT_YES, .num_inp = LONG_MAX,
     .mobile_len = 1, .mobile_ret = {-7}},
    {.name = "mobile-zero", .b_mobile_h = TAUT_YES,
     .i_mobile_h = TAUT_YES, .mobile_len = 1, .mobile_ret = {0},
     .primary_flags = IDIF_PROBLEM},
    {.name = "mobile-positive", .b_mobile_h = TAUT_YES,
     .i_mobile_h = TAUT_YES, .mobile_len = 1, .mobile_ret = {7},
     .primary_flags = IDIF_PROBLEM},
    {.name = "non-enum-mobile-h", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .primary_flags = IDIF_PROBLEM},
    {.name = "primary-problem-precedence", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .primary_flags = IDIF_PROBLEM,
     .primary_error = 9},
    {.name = "primary-error", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .primary_error = 9},
    {.name = "primary-ordinary-flags", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .primary_flags = IDIF_MORE_H,
     .stereo_len = 1, .stereo_ret = {0}},
    {.name = "optional-skipped", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .fixed_len = 1, .fixed_ret = {0},
     .stereo_len = 1, .stereo_ret = {0}},
    {.name = "optional-original-only", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .original_layer1 = 1,
     .fixed_len = 1, .fixed_ret = {0}, .stereo_len = 1,
     .stereo_ret = {0}},
    {.name = "optional-reversed-only", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .reversed_layer1 = 1,
     .fixed_len = 1, .fixed_ret = {0}, .stereo_len = 1,
     .stereo_ret = {0}},
    {.name = "optional-both", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .original_layer1 = 1,
     .reversed_layer1 = 1, .fixed_len = 1, .fixed_ret = {0},
     .stereo_len = 1, .stereo_ret = {0}},
    {.name = "optional-problem-primary-check", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .original_layer1 = 1,
     .reversed_layer1 = 1, .fixed_len = 1, .fixed_ret = {0},
     .optional_flags = IDIF_PROBLEM, .stereo_len = 1,
     .stereo_ret = {7}},
    {.name = "optional-error", .b_mobile_h = TAUT_NON,
     .i_mobile_h = TAUT_NON, .original_layer1 = 1,
     .reversed_layer1 = 1, .fixed_len = 1, .fixed_ret = {0},
     .optional_error = 9},
    {.name = "stereo-negative", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {-7}},
    {.name = "stereo-zero", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {0}},
    {.name = "stereo-positive", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {7}},
    {.name = "endpoint-zero-null-tinfo", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {0}},
    {.name = "endpoint-mixed", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {0},
     .tc_len = 2, .tc_edges = {2, 3}, .tinfo_present = 1,
     .tinfo_len = 2, .tinfo_h = {0, 1},
     .tinfo_endpoints = {5, 7}},
    {.name = "endpoint-all-zero-h", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {0},
     .tinfo_present = 1, .tinfo_len = 2,
     .tinfo_endpoints = {USHRT_MAX, USHRT_MAX}},
    {.name = "endpoint-signed-max", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {0},
     .tc_len = 1, .tc_edges = {INT_MAX}, .tinfo_present = 1,
     .tinfo_len = 1, .tinfo_h = {USHRT_MAX},
     .tinfo_endpoints = {USHRT_MAX}},
    {.name = "endpoint-signed-min", .b_mobile_h = TAUT_INI,
     .i_mobile_h = TAUT_INI, .stereo_len = 1, .stereo_ret = {0},
     .tc_len = 1, .tc_edges = {INT_MIN}}
};

static int print_normalize_and_compare_final_case(
    const ORACLE_NORMALIZE_FINAL_CASE *test_case)
{
    CANON_GLOBALS canonical_globals;
    struct tagINCHI_CLOCK clock;
    INPUT_PARMS parameters;
    STRUCT_DATA data;
    BN_STRUCT bns;
    BN_DATA bns_data;
    BNS_VERTEX vertex;
    StrFromINChI structure, expected_structure;
    inp_ATOM atoms[1], atoms2[1], atoms3[1];
    VAL_AT valence[1];
    ALL_TC_GROUPS groups;
    TC_GROUP tc_groups[2];
    SRM restore_mode;
    INChI original[2];
    INChI *input[TAUT_NUM];
    int result, runs = 17, delta = -19, index;
    int optional_active, holders_null, complete_caller_state_exact;

    memset(&canonical_globals, 0, sizeof(canonical_globals));
    memset(&clock, 0, sizeof(clock));
    memset(&parameters, 0, sizeof(parameters));
    memset(&data, 0, sizeof(data));
    memset(&bns, 0, sizeof(bns));
    memset(&bns_data, 0, sizeof(bns_data));
    memset(&vertex, 0, sizeof(vertex));
    memset(&structure, 0, sizeof(structure));
    memset(atoms, 0, sizeof(atoms));
    memset(atoms2, 0, sizeof(atoms2));
    memset(atoms3, 0, sizeof(atoms3));
    memset(valence, 0, sizeof(valence));
    memset(&groups, 0, sizeof(groups));
    memset(tc_groups, 0, sizeof(tc_groups));
    memset(&restore_mode, 0, sizeof(restore_mode));
    memset(original, 0, sizeof(original));

    parameters.nMode = REQ_MODE_TAUT | REQ_MODE_NON_ISO;
    bns.num_atoms = bns.num_vertices = 1;
    bns.vert = &vertex;
    atoms[0].el_number = 6;
    atoms[0].num_H = 4;
    atoms[0].orig_at_number = atoms[0].component = 1;
    strcpy(atoms[0].elname, "C");
    structure.at = atoms;
    structure.num_atoms = 1;
    structure.bMobileH = (char) test_case->b_mobile_h;
    structure.iMobileH = (char) test_case->i_mobile_h;
    structure.pSrm = &restore_mode;
    original[0].nNumberOfAtoms = 1;
    original[0].nErrorCode = 70;
    original[1].nNumberOfAtoms = test_case->original_layer1 ? 1 : 0;
    original[1].nErrorCode = 71;
    input[0] = &original[0];
    input[1] = test_case->original_layer1 ? &original[1] : NULL;

    ORACLE_NORMALIZE_ALLOCATION_COUNT = 0;
    structure.pOneINChI[0] = calloc(1, sizeof(INChI));
    structure.pOne_norm_data[0] = calloc(1, sizeof(INP_ATOM_DATA));
    if (!structure.pOneINChI[0] || !structure.pOne_norm_data[0])
    {
        return 70;
    }
    structure.pOneINChI[0]->nNumberOfAtoms = 1;
    structure.pOneINChI[0]->nErrorCode = 100;
    structure.pOne_norm_data[0]->num_at = 30;
    structure.pOne_norm_data[0]->num_bonds = 40;
    structure.pOne_norm_data[0]->at = calloc(1, sizeof(inp_ATOM));
    if (!structure.pOne_norm_data[0]->at)
    {
        return 70;
    }
    oracle_normalize_register_allocation(structure.pOneINChI[0],
                                         "pOneINChI[0]");
    oracle_normalize_register_allocation(structure.pOne_norm_data[0],
                                         "pOne_norm_data[0]");
    oracle_normalize_register_allocation(structure.pOne_norm_data[0]->at,
                                         "pOne_norm_data[0].at");
    if (test_case->reversed_layer1)
    {
        structure.pOneINChI[1] = calloc(1, sizeof(INChI));
        if (!structure.pOneINChI[1])
        {
            return 70;
        }
        structure.pOneINChI[1]->nNumberOfAtoms = 1;
        structure.pOneINChI[1]->nErrorCode = 101;
        oracle_normalize_register_allocation(structure.pOneINChI[1],
                                             "pOneINChI[1]");
    }
    groups.num_tgroups = test_case->tc_len;
    groups.pTCG = test_case->tc_len ? tc_groups : NULL;
    for (index = 0; index < test_case->tc_len; index++)
    {
        tc_groups[index].num_edges = test_case->tc_edges[index];
    }
    if (test_case->tinfo_present)
    {
        structure.One_ti.t_group =
            calloc((size_t) test_case->tinfo_len, sizeof(T_GROUP));
        if (!structure.One_ti.t_group)
        {
            return 70;
        }
        structure.One_ti.num_t_groups = test_case->tinfo_len;
        oracle_normalize_register_allocation(structure.One_ti.t_group,
                                             "One_ti.t_group");
        for (index = 0; index < test_case->tinfo_len; index++)
        {
            structure.One_ti.t_group[index].num[0] =
                (AT_RANK) test_case->tinfo_h[index];
            structure.One_ti.t_group[index].nNumEndpoints =
                (AT_NUMB) test_case->tinfo_endpoints[index];
        }
    }

    expected_structure = structure;
    expected_structure.pOneINChI[0] = expected_structure.pOneINChI[1] = NULL;
    expected_structure.pOne_norm_data[0] = NULL;
    memset(&expected_structure.One_ti, 0, sizeof(expected_structure.One_ti));

    ORACLE_NORMALIZE_REBUILD_LENGTH = 1;
    ORACLE_NORMALIZE_REBUILD_RETURNS[0] = 0;
    ORACLE_NORMALIZE_REBUILD_PREP[0] = 1;
    ORACLE_NORMALIZE_REBUILD_NORM[0] = 1;
    ORACLE_NORMALIZE_REBUILD_TINFO[0] = test_case->tinfo_present;
    optional_active = test_case->b_mobile_h == TAUT_NON &&
                      (test_case->original_layer1 ||
                       test_case->reversed_layer1);
    ORACLE_NORMALIZE_COMPARE_LENGTH = 2 + optional_active;
    memset(ORACLE_NORMALIZE_COMPARE_FLAGS, 0,
           sizeof(ORACLE_NORMALIZE_COMPARE_FLAGS));
    memset(ORACLE_NORMALIZE_COMPARE_ERRORS, 0,
           sizeof(ORACLE_NORMALIZE_COMPARE_ERRORS));
    ORACLE_NORMALIZE_COMPARE_FLAGS[1] = test_case->primary_flags;
    ORACLE_NORMALIZE_COMPARE_ERRORS[1] = test_case->primary_error;
    if (optional_active)
    {
        ORACLE_NORMALIZE_COMPARE_FLAGS[2] = test_case->optional_flags;
        ORACLE_NORMALIZE_COMPARE_ERRORS[2] = test_case->optional_error;
    }
    ORACLE_NORMALIZE_FILL_LENGTH =
        test_case->i_mobile_h == TAUT_NON ? 1 : 0;
    ORACLE_NORMALIZE_FILL_RETURNS[0] = 0;
    ORACLE_NORMALIZE_FIX_FIXED_LENGTH = test_case->fixed_len;
    ORACLE_NORMALIZE_FIX_MOBILE_LENGTH = test_case->mobile_len;
    ORACLE_NORMALIZE_FIX_STEREO_LENGTH = test_case->stereo_len;
    memcpy(ORACLE_NORMALIZE_FIX_FIXED_RETURNS, test_case->fixed_ret,
           sizeof(test_case->fixed_ret));
    memcpy(ORACLE_NORMALIZE_FIX_MOBILE_RETURNS, test_case->mobile_ret,
           sizeof(test_case->mobile_ret));
    memcpy(ORACLE_NORMALIZE_FIX_STEREO_RETURNS, test_case->stereo_ret,
           sizeof(test_case->stereo_ret));
    ORACLE_NORMALIZE_REBUILD_POSITION = ORACLE_NORMALIZE_COMPARE_POSITION = 0;
    ORACLE_NORMALIZE_FILL_POSITION = ORACLE_NORMALIZE_FIX_FIXED_POSITION = 0;
    ORACLE_NORMALIZE_FIX_MOBILE_POSITION = 0;
    ORACLE_NORMALIZE_FIX_STEREO_POSITION = 0;
    ORACLE_NORMALIZE_EVENTS[0] = '\0';
    ORACLE_NORMALIZE_EVENTS_LENGTH = 0;
    ORACLE_NORMALIZE_PREFREE_EXACT = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[0] = 1;
    ORACLE_NORMALIZE_EXPECTED_INCHI_ATOMS[1] =
        test_case->reversed_layer1 ? 1 : 0;
    ORACLE_NORMALIZE_REVERSED[0] = structure.pOneINChI[0];
    ORACLE_NORMALIZE_REVERSED[1] = structure.pOneINChI[1];
    ORACLE_NORMALIZE_ORIGINAL[0] = input[0];
    ORACLE_NORMALIZE_ORIGINAL[1] = input[1];
    ORACLE_NORMALIZE_INCHI_HOLDERS[0] = &structure.pOneINChI[0];
    ORACLE_NORMALIZE_INCHI_HOLDERS[1] = &structure.pOneINChI[1];
    ORACLE_NORMALIZE_AUX_HOLDERS[0] = &structure.pOneINChI_Aux[0];
    ORACLE_NORMALIZE_AUX_HOLDERS[1] = &structure.pOneINChI_Aux[1];
    ORACLE_NORMALIZE_NORM_HOLDERS[0] = structure.pOne_norm_data[0];
    ORACLE_NORMALIZE_NORM_HOLDERS[1] = NULL;
    ORACLE_NORMALIZE_EXPECTED_NORM_AT[0] =
        structure.pOne_norm_data[0]->at;
    ORACLE_NORMALIZE_EXPECTED_NORM_AT[1] = NULL;
    ORACLE_NORMALIZE_EXPECTED_TINFO_PRESENT = test_case->tinfo_present;
    ORACLE_NORMALIZE_EXPECTED_TINFO_LENGTH = test_case->tinfo_len;
    for (index = 0; index < 2; index++)
    {
        ORACLE_NORMALIZE_EXPECTED_TINFO_H[index] =
            (AT_RANK) test_case->tinfo_h[index];
        ORACLE_NORMALIZE_EXPECTED_TINFO_ENDPOINTS[index] =
            (AT_NUMB) test_case->tinfo_endpoints[index];
    }
    ORACLE_DEFERRED_FREE_COUNT = 0;
    ORACLE_DEFER_FREES = 1;
    ORACLE_NORMALIZE_FINAL_ACTIVE = ORACLE_NORMALIZE_ACTIVE = 1;
    result = NormalizeAndCompare(
        &canonical_globals, &clock, &parameters, &data, &bns, &bns_data,
        &structure, atoms, atoms2, atoms3, valence, &groups, input,
        test_case->num_inp, 0, &runs, &delta, 4, 0);
    ORACLE_NORMALIZE_ACTIVE = ORACLE_NORMALIZE_FINAL_ACTIVE = 0;

    holders_null = !structure.pOneINChI[0] && !structure.pOneINChI[1] &&
                   !structure.pOneINChI_Aux[0] &&
                   !structure.pOneINChI_Aux[1] &&
                   !structure.pOne_norm_data[0] &&
                   !structure.pOne_norm_data[1];
    complete_caller_state_exact =
        memcmp(&structure, &expected_structure, sizeof(structure)) == 0 &&
        input[0] == &original[0] &&
        input[1] == (test_case->original_layer1 ? &original[1] : NULL) &&
        original[0].nErrorCode == 70 && original[1].nErrorCode == 71;

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\"," 
           "\"case_id\":\"normalize-final-%s\"," 
           "\"operation\":\"NormalizeAndCompare\"," 
           "\"family\":\"final\",\"input\":{" 
           "\"b_mobile_h\":%d,\"i_mobile_h\":%d," 
           "\"original_layer1\":%s,\"reversed_layer1\":%s," 
           "\"num_inp\":%ld,\"fixed_ret\":",
           test_case->name, test_case->b_mobile_h, test_case->i_mobile_h,
           test_case->original_layer1 ? "true" : "false",
           test_case->reversed_layer1 ? "true" : "false",
           test_case->num_inp);
    print_normalize_int_array(test_case->fixed_ret, test_case->fixed_len);
    printf(",\"mobile_ret\":");
    print_normalize_int_array(test_case->mobile_ret, test_case->mobile_len);
    printf(",\"primary_flags\":%" PRIu64 ",\"primary_error\":%d," 
           "\"optional_flags\":%" PRIu64 ",\"optional_error\":%d," 
           "\"stereo_ret\":",
           (uint64_t) test_case->primary_flags, test_case->primary_error,
           (uint64_t) test_case->optional_flags, test_case->optional_error);
    print_normalize_int_array(test_case->stereo_ret, test_case->stereo_len);
    printf(",\"tc_edges\":");
    print_normalize_int_array(test_case->tc_edges, test_case->tc_len);
    printf(",\"tinfo_present\":%s,\"tinfo_h\":",
           test_case->tinfo_present ? "true" : "false");
    print_normalize_int_array(test_case->tinfo_h, test_case->tinfo_len);
    printf(",\"tinfo_endpoints\":");
    print_normalize_int_array(test_case->tinfo_endpoints,
                              test_case->tinfo_len);
    printf("},\"output\":{\"result\":%d,\"runs\":%d,\"delta\":%d," 
           "\"rebuild_used\":%d,\"compare_used\":%d," 
           "\"fill_used\":%d,\"fixed_used\":%d," 
           "\"mobile_used\":%d,\"stereo_used\":%d," 
           "\"complete_caller_state_exact\":%s," 
           "\"prefree_state_exact\":%s,\"holders_null\":%s," 
           "\"one_ti_cleared\":%s,\"input_pointers_preserved\":%s," 
           "\"allocation_free_exact\":%s,\"source_allocations_before\":%zu," 
           "\"source_allocations_after\":%d,\"cleanup_events\":[%s]}}\n",
           result, runs, delta, ORACLE_NORMALIZE_REBUILD_POSITION,
           ORACLE_NORMALIZE_COMPARE_POSITION, ORACLE_NORMALIZE_FILL_POSITION,
           ORACLE_NORMALIZE_FIX_FIXED_POSITION,
           ORACLE_NORMALIZE_FIX_MOBILE_POSITION,
           ORACLE_NORMALIZE_FIX_STEREO_POSITION,
           complete_caller_state_exact ? "true" : "false",
           ORACLE_NORMALIZE_PREFREE_EXACT ? "true" : "false",
           holders_null ? "true" : "false",
           !structure.One_ti.t_group && !structure.One_ti.num_t_groups
               ? "true" : "false",
           input[0] == &original[0] &&
                   input[1] ==
                       (test_case->original_layer1 ? &original[1] : NULL)
               ? "true" : "false",
           ORACLE_DEFERRED_FREE_COUNT == ORACLE_NORMALIZE_ALLOCATION_COUNT
               ? "true" : "false",
           ORACLE_NORMALIZE_ALLOCATION_COUNT, 0,
           ORACLE_NORMALIZE_EVENTS);
    fflush(stdout);
    ORACLE_DEFER_FREES = 0;
    oracle_flush_deferred_frees();
    return 0;
}

static int print_normalize_and_compare_records(void)
{
    static const int forced_returns[] = {RI_ERR_ALLOC, RI_ERR_PROGR};
    unsigned int return_index;
    unsigned int holder_mask;
    for (return_index = 0;
         return_index < sizeof(forced_returns) / sizeof(forced_returns[0]);
         return_index++)
    {
        for (holder_mask = 0; holder_mask < 64; holder_mask++)
        {
            int status = print_normalize_and_compare_negative_case(
                holder_mask, forced_returns[return_index]);
            if (status != 0)
            {
                return status;
            }
        }
    }
    for (return_index = 0; return_index < 2; return_index++)
    {
        int mobile_h = return_index == 0 ? TAUT_NON : TAUT_YES;
        int original_state;
        for (original_state = 0; original_state < 4; original_state++)
        {
            int reversed_state;
            for (reversed_state = 0; reversed_state < 4; reversed_state++)
            {
                int status = print_normalize_and_compare_layer_case(
                    mobile_h, original_state, reversed_state, 1);
                if (status != 0)
                {
                    return status;
                }
            }
        }
    }
    if (print_normalize_and_compare_layer_case(TAUT_INI, 0, 0, 0) != 0)
    {
        return 70;
    }
    if (run_normalize_and_compare_zz_case(0, 1, "CZzZz2", 0, 0, 1) ||
        run_normalize_and_compare_zz_case(1, 0, "CZzZz2", 0, 0, 1) ||
        run_normalize_and_compare_zz_case(1, 1, NULL, 0, 0, 1) ||
        run_normalize_and_compare_zz_case(1, 1, "C2H6", 0, 0, 1) ||
        run_normalize_and_compare_zz_case(1, 1, "CZzZz2", 0, 0, 1) ||
        run_normalize_and_compare_zz_case(1, 1, "CZzZz2", 2, 0, 1) ||
        run_normalize_and_compare_zz_case(1, 1, "CZzZz2", 3, 0, 1) ||
        run_normalize_and_compare_zz_case(1, 1, "CZzZz2", 0, 1, 1) ||
        run_normalize_and_compare_zz_case(1, 1, "CZzZz2", 4, 1, 1))
    {
        return 70;
    }
    if (print_normalize_and_compare_zz_undefined_case() != 0)
    {
        return 70;
    }
    for (return_index = 0;
         return_index < sizeof(ORACLE_NORMALIZE_FIRST_LESS_CASES) /
                            sizeof(ORACLE_NORMALIZE_FIRST_LESS_CASES[0]);
         return_index++)
    {
        if (print_normalize_and_compare_first_less_case(
                ORACLE_NORMALIZE_FIRST_LESS_CASES + return_index) != 0)
        {
            return 70;
        }
    }
    for (return_index = 0;
         return_index < sizeof(ORACLE_NORMALIZE_MORE_EXTRA_CASES) /
                            sizeof(ORACLE_NORMALIZE_MORE_EXTRA_CASES[0]);
         return_index++)
    {
        if (print_normalize_and_compare_first_less_case(
                ORACLE_NORMALIZE_MORE_EXTRA_CASES + return_index) != 0)
        {
            return 70;
        }
    }
    for (return_index = 0;
         return_index < sizeof(ORACLE_NORMALIZE_FINAL_CASES) /
                            sizeof(ORACLE_NORMALIZE_FINAL_CASES[0]);
         return_index++)
    {
        if (print_normalize_and_compare_final_case(
                ORACLE_NORMALIZE_FINAL_CASES + return_index) != 0)
        {
            return 70;
        }
    }
    return 0;
}

static int print_version_record(void)
{
    const char *version = IXA_INCHIBUILDER_GetInChIVersion(IXA_FALSE);
    const char *description = IXA_INCHIBUILDER_GetInChIVersion(IXA_TRUE);
    if (!version || strcmp(version, EXPECTED_VERSION) != 0)
    {
        fprintf(stderr, "unexpected official InChI version\n");
        return 65;
    }
    if (!description || strcmp(description, EXPECTED_DESCRIPTION) != 0)
    {
        fprintf(stderr, "unexpected official InChI description\n");
        return 66;
    }

    printf("{\"schema_version\":\"cosmolkit-inchi-official-c-v1\","
           "\"oracle\":{\"project\":\"IUPAC InChI\",\"tag\":\"v1.07.5\","
           "\"commit\":\"11a87982bb518f57ac013f0b258c283655e1ea1d\","
           "\"tree_sha256\":\"4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd\","
           "\"api_version\":\"%s\"},\"case_id\":\"oracle-version\","
           "\"operation\":\"version\",\"input\":{\"full_description\":false},"
           "\"output\":{\"status\":0,\"value\":\"%s\"}}\n",
           version, version);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        fprintf(stderr,
                "usage: official-c-oracle --version-record|--inchi-to-inchi-atom-records|--inchi-to-inchi-input-records|--get-inchi-input-from-aux-info-records|--get-std-inchi-input-from-aux-info-records|--write-coord-records|--parse-options-records|--element-lookup-records|--periodic-lookup-records|--el-valence-records|--metal-records|--atomic-mass-records|--detect-unusual-valence-records|--extract-charge-records|--extract-hydrogen-records|--list-records|--bonds-to-metal-records|--set-atom-records|--set-bond-records|--cmp-components-records|--cmp-rad-endpoints-records|--comp-t-group-number-records|--comp-c-group-number-records|--cmp-iso-atw-diff-component-no-records|--comp-neigh-lists-records|--comp-neigh-lists-up-to-max-rank-records|--cmp-charge-val-records|--comp-cc-cand-records|--base26-triplet-1-records|--base26-triplet-2-records|--base26-triplet-3-records|--base26-triplet-4-records|--base26-dublet-for-bits-28-to-36-records|--base26-dublet-for-bits-56-to-64-records|--get-xtra-hash-major-hex-records|--get-xtra-hash-minor-hex-records|--sha2-starts-records|--sha2-process-records|--sha2-update-records|--sha2-finish-records|--sha2-csum-records|--inchi-key-records|--rdkit-core-root-records|--normalize-and-compare-records|--cn-list-record\n");
        return 64;
    }
    if (strcmp(argv[1], "--version-record") == 0)
    {
        return print_version_record();
    }
    if (strcmp(argv[1], "--inchi-to-inchi-atom-records") == 0)
    {
        return print_inchi_to_atom_records();
    }
    if (strcmp(argv[1], "--inchi-to-inchi-input-records") == 0)
    {
        return print_inchi_to_input_records();
    }
    if (strcmp(argv[1], "--get-inchi-input-from-aux-info-records") == 0)
    {
        return print_aux_input_records(0);
    }
    if (strcmp(argv[1], "--get-std-inchi-input-from-aux-info-records") == 0)
    {
        return print_aux_input_records(1);
    }
    if (strcmp(argv[1], "--write-coord-records") == 0)
    {
        return print_write_coord_records();
    }
    if (strcmp(argv[1], "--parse-options-records") == 0)
    {
        return print_parse_options_records();
    }
    if (strcmp(argv[1], "--element-lookup-records") == 0)
    {
        return print_element_lookup_records();
    }
    if (strcmp(argv[1], "--periodic-lookup-records") == 0)
    {
        return print_periodic_lookup_records();
    }
    if (strcmp(argv[1], "--el-valence-records") == 0)
    {
        return print_el_valence_records();
    }
    if (strcmp(argv[1], "--metal-records") == 0)
    {
        return print_metal_records();
    }
    if (strcmp(argv[1], "--atomic-mass-records") == 0)
    {
        return print_atomic_mass_records();
    }
    if (strcmp(argv[1], "--detect-unusual-valence-records") == 0)
    {
        return print_detect_unusual_valence_records();
    }
    if (strcmp(argv[1], "--extract-charge-records") == 0)
    {
        return print_extract_charge_records();
    }
    if (strcmp(argv[1], "--extract-hydrogen-records") == 0)
    {
        return print_extract_hydrogen_records();
    }
    if (strcmp(argv[1], "--list-records") == 0)
    {
        return print_list_records();
    }
    if (strcmp(argv[1], "--bonds-to-metal-records") == 0)
    {
        return print_bonds_to_metal_records();
    }
    if (strcmp(argv[1], "--set-atom-records") == 0)
    {
        return print_set_atom_records();
    }
    if (strcmp(argv[1], "--set-bond-records") == 0)
    {
        return print_set_bond_records();
    }
    if (strcmp(argv[1], "--cmp-components-records") == 0)
    {
        return print_cmp_components_records();
    }
    if (strcmp(argv[1], "--cmp-rad-endpoints-records") == 0)
    {
        return print_cmp_rad_endpoints_records();
    }
    if (strcmp(argv[1], "--comp-t-group-number-records") == 0)
    {
        return print_comp_t_group_number_records();
    }
    if (strcmp(argv[1], "--comp-c-group-number-records") == 0)
    {
        return print_comp_c_group_number_records();
    }
    if (strcmp(argv[1], "--cmp-iso-atw-diff-component-no-records") == 0)
    {
        return print_cmp_iso_atw_diff_component_no_records();
    }
    if (strcmp(argv[1], "--comp-neigh-lists-records") == 0)
    {
        return print_comp_neigh_lists_records();
    }
    if (strcmp(argv[1], "--comp-neigh-lists-up-to-max-rank-records") == 0)
    {
        return print_comp_neigh_lists_up_to_max_rank_records();
    }
    if (strcmp(argv[1], "--cmp-charge-val-records") == 0)
    {
        return print_cmp_charge_val_records();
    }
    if (strcmp(argv[1], "--comp-cc-cand-records") == 0)
    {
        return print_comp_cc_cand_records();
    }
    if (strcmp(argv[1], "--base26-triplet-1-records") == 0)
    {
        return print_base26_triplet_1_records();
    }
    if (strcmp(argv[1], "--base26-triplet-2-records") == 0)
    {
        return print_base26_triplet_2_records();
    }
    if (strcmp(argv[1], "--base26-triplet-3-records") == 0)
    {
        return print_base26_triplet_3_records();
    }
    if (strcmp(argv[1], "--base26-triplet-4-records") == 0)
    {
        return print_base26_triplet_4_records();
    }
    if (strcmp(argv[1], "--base26-dublet-for-bits-28-to-36-records") == 0)
    {
        return print_base26_dublet_for_bits_28_to_36_records();
    }
    if (strcmp(argv[1], "--base26-dublet-for-bits-56-to-64-records") == 0)
    {
        return print_base26_dublet_for_bits_56_to_64_records();
    }
    if (strcmp(argv[1], "--get-xtra-hash-major-hex-records") == 0)
    {
        return print_get_xtra_hash_major_hex_records();
    }
    if (strcmp(argv[1], "--get-xtra-hash-minor-hex-records") == 0)
    {
        return print_get_xtra_hash_minor_hex_records();
    }
    if (strcmp(argv[1], "--sha2-starts-records") == 0)
    {
        return print_sha2_starts_records();
    }
    if (strcmp(argv[1], "--sha2-process-records") == 0)
    {
        return print_sha2_process_records();
    }
    if (strcmp(argv[1], "--sha2-update-records") == 0)
    {
        return print_sha2_update_records();
    }
    if (strcmp(argv[1], "--sha2-finish-records") == 0)
    {
        return print_sha2_finish_records();
    }
    if (strcmp(argv[1], "--sha2-csum-records") == 0)
    {
        return print_sha2_csum_records();
    }
    if (strcmp(argv[1], "--inchi-key-records") == 0)
    {
        return print_inchi_key_records();
    }
    if (strcmp(argv[1], "--rdkit-core-root-records") == 0)
    {
        return print_rdkit_core_root_records();
    }
    if (strcmp(argv[1], "--normalize-and-compare-records") == 0)
    {
        return print_normalize_and_compare_records();
    }
    if (strcmp(argv[1], "--cn-list-record") == 0)
    {
        return print_cn_list_record();
    }
    fprintf(stderr, "unknown official C oracle operation\n");
    return 64;
}
