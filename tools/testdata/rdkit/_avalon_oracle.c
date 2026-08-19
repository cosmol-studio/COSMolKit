/*
 * Test-only Avalon oracle harness.
 *
 * The harness calls the pinned RDKit Avalon C ABI directly. It intentionally
 * mirrors the adapter's public nBits/8 conversion and internal four-byte
 * rounding so the generated records describe the selected GetAvalonFP path.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reaccs.h"
#include "ssmatch.h"
#include "smi2mol.h"
#include "utilities.h"

static void emit(const char *line) {
  char *end;
  unsigned long nbits;
  unsigned long flags;
  long query;
  char *smiles = strdup(line);
  char *field;
  struct reaccs_molecule_t *molecule;
  unsigned nbytes;
  unsigned internal_nbytes;
  unsigned char *bits;
  int *counts;
  unsigned i;

  if (smiles == NULL) {
    puts("error\tallocation");
    return;
  }
  field = strchr(smiles, '\t');
  if (field == NULL) {
    puts("error\tmalformed_input");
    free(smiles);
    return;
  }
  *field++ = '\0';
  nbits = strtoul(field, &end, 10);
  if (*end != '\t') {
    puts("error\tmalformed_nbits");
    free(smiles);
    return;
  }
  query = strtol(end + 1, &end, 10);
  if (*end != '\t') {
    puts("error\tmalformed_query");
    free(smiles);
    return;
  }
  flags = strtoul(end + 1, &end, 16);
  if (*end != '\0') {
    puts("error\tmalformed_flags");
    free(smiles);
    return;
  }

  nbytes = (unsigned)(nbits / 8UL);
  internal_nbytes = (nbytes + 3U) & ~3U;
  bits = (unsigned char *)calloc(internal_nbytes ? internal_nbytes : 1U, 1U);
  counts = (int *)calloc((internal_nbytes ? internal_nbytes : 1U) * 8U,
                         sizeof(int));
  molecule = SMIToMOL(smiles, DY_AROMATICITY);
  if (bits == NULL || counts == NULL || molecule == NULL) {
    puts(molecule == NULL ? "error\tconversion" : "error\tallocation");
    FreeMolecule(molecule);
    free(counts);
    free(bits);
    free(smiles);
    return;
  }

  SetFingerprintBits(molecule, (char *)bits, internal_nbytes, (int)flags,
                     (int)query, 0);
  if (!query) {
    SetFingerprintBits(molecule, (char *)bits, internal_nbytes, (int)flags, 0,
                       ACCUMULATE_BITS | USE_DY_AROMATICITY);
  }
  SetFingerprintCountsWithFocus(molecule, counts, internal_nbytes * 8U,
                                (int)flags, (int)query, 0, 0);

  printf("ok\t%s\t%lu\t%ld\t%06lx\t", smiles, nbits, query, flags);
  for (i = 0; i < nbytes; ++i) printf("%02x", bits[i]);
  putchar('\t');
  for (i = 0; i < nbytes * 8U; ++i) {
    if (counts[i] != 0) printf("%u:%d,", i, counts[i]);
  }
  putchar('\n');

  FreeMolecule(molecule);
  free(counts);
  free(bits);
  free(smiles);
}

int main(void) {
  char line[65536];
  while (fgets(line, sizeof(line), stdin) != NULL) {
    size_t length = strlen(line);
    if (length > 0 && line[length - 1] == '\n') line[length - 1] = '\0';
    if (line[0] != '\0') emit(line);
  }
  return 0;
}
