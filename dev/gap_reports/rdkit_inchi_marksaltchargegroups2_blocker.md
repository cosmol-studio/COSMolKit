# `MarkSaltChargeGroups2` Candidate-Storage Audit

## Scope

- Current planned Port step: 1963
- C source: `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2961-3453`
- Selected build: `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`, GCC/Linux
- Active macro branch: `IGNORE_TGROUP_WITHOUT_H == 1`

## Audited Read

The active source at `ichitaut.c:3067` reads the first slot after the newly
collected candidates when capacity remains:

```c
if ((nNumCandidates < nMaxNumCandidates) &&
    (2 == s_candidate[nNumCandidates].type))
```

This read is not an uninitialized allocation read in the selected production
call path.

## Allocation And Lifetime Proof

1. `ichi_bns.c:5283` clears the owning `S_GROUP_INFO` before setup.
2. `ichi_bns.c:5344-5346` is the production allocation site for
   `s_group_info->s_candidate`; it uses
   `inchi_calloc(num_atoms, sizeof(s_group_info->s_candidate[0]))`.
3. Successful allocation sets `max_num_candidates = num_atoms` at
   `ichi_bns.c:5347-5349`.
4. The same allocation is reused across the salt-processing loop. The call at
   `ichi_bns.c:5755` runs `MarkSaltChargeGroups`, and the call at
   `ichi_bns.c:5791` subsequently runs `MarkSaltChargeGroups2`; neither call
   reallocates or clears the candidate array.
5. Candidate-producing functions overwrite only their current prefix. They
   reset counters such as `num_candidates`, but do not zero the unused tail.
   A later invocation may therefore observe a `type` value left by an earlier
   invocation in slot `nNumCandidates`.
6. The allocation remains live until `ichi_bns.c:5949-5951`, where it is freed
   on function exit.

The first invocation starts with a zeroed spare slot because of `inchi_calloc`.
Later invocations have deterministic stateful behavior: the spare slot is the
previously initialized array element, possibly carrying a prior candidate's
`type`. This is observable source behavior because a value of `2` disables
every current candidate in the loop at `ichitaut.c:3064-3085`.

## Rust Reproduction Requirement

The Rust port must:

- allocate production candidate storage with zero initialization;
- reuse the same `SourceHeap` allocation across calls;
- overwrite only the candidate prefix written by the C loop;
- read the preserved slot at index `nNumCandidates` when
  `nNumCandidates < nMaxNumCandidates`;
- avoid clearing the unused tail on entry or exit;
- test both an initially zero spare slot and a preserved spare slot whose
  `type` is `2`.

## Disposition

The previous blocker conclusion was incorrect. The selected production path
defines the allocation state and lifetime needed to reproduce the read.
`MarkSaltChargeGroups2` is therefore not blocked on allocation provenance and
must be ported with the stateful spare-slot behavior above.
