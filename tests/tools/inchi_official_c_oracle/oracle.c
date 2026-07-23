#include <float.h>
#include <inttypes.h>
#include <locale.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mode.h"
#include "inchi_api.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "ichierr.h"
#include "inpdef.h"
#include "ichi.h"
#include "ichi_io.h"
#include "inchi_dll_b.h"
#include "readinch.h"
#include "util.h"

void *__real_malloc(size_t size);
void *__real_calloc(size_t count, size_t size);
void __real_free(void *pointer);

static int ORACLE_ALLOCATION_FAILURE_ENABLED = 0;
static int ORACLE_ALLOCATION_ORDINAL = 0;
static int ORACLE_ALLOCATION_CALLS = 0;
static int ORACLE_DEFER_FREES = 0;
static void *ORACLE_DEFERRED_FREES[4096];
static size_t ORACLE_DEFERRED_FREE_COUNT = 0;

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
    if (oracle_allocation_should_fail())
    {
        return NULL;
    }
    return __real_calloc(count, size);
}

void __wrap_free(void *pointer)
{
    size_t i;
    if (!pointer)
    {
        return;
    }
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
                "usage: official-c-oracle --version-record|--inchi-to-inchi-atom-records|--inchi-to-inchi-input-records|--get-inchi-input-from-aux-info-records|--get-std-inchi-input-from-aux-info-records|--write-coord-records|--parse-options-records|--element-lookup-records|--periodic-lookup-records|--el-valence-records|--metal-records|--atomic-mass-records|--detect-unusual-valence-records|--extract-charge-records|--extract-hydrogen-records|--list-records|--bonds-to-metal-records|--set-atom-records|--set-bond-records\n");
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
    fprintf(stderr, "unknown official C oracle operation\n");
    return 64;
}
