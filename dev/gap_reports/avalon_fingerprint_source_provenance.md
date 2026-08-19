# Avalon Fingerprint Source Provenance

Status: Step 51 audit complete. The Avalon fingerprint remains fail-closed;
this report records the source and licensing boundary required before any Rust
port is presented as supported.

## Pinned sources

The RDKit adapter is from the vendored RDKit source at the pinned project
revision (`2026.03.1`, revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`). The adapter files are:

- `third_party/rdkit/External/AvalonTools/AvalonTools.h`
- `third_party/rdkit/External/AvalonTools/AvalonTools.cpp`
- `third_party/rdkit/External/AvalonTools/Wrap/pyAvalonTools.cpp`
- `third_party/rdkit/External/AvalonTools/Wrap/testAvalonTools.py`
- `third_party/rdkit/External/AvalonTools/CMakeLists.txt`

The engine source is `rdkit/ava-formake`, tag
`AvalonToolkit_2.0.5-pre.3`, downloaded from:

`https://github.com/rdkit/ava-formake/archive/refs/tags/AvalonToolkit_2.0.5-pre.3.tar.gz`

The archive was read-only audited in a system temporary directory. The
verified checksums are:

```text
MD5:    7a20c25a7e79f3344e0f9f49afa03351
SHA256: a553cc9e04c2485891a2e215e309b71db62161b4cb2bff8efbb074c72441c487
```

The archive has not yet been copied into the repository's ignored source
staging directory. That is the separate Step 53 materialization task.

## Adapter defaults and call boundary

`AvalonTools.h` declares the C++ molecule overload with:

```text
nBits=512, isQuery=false, resetVect=true, bitFlags=avalonSSSBits
avalonSSSBits=0x007FFF
avalonSimilarityBits=0xF07FFF
```

The Python wrapper in `Wrap/pyAvalonTools.cpp` intentionally exposes a
different public default for `GetAvalonFP`:

```text
nBits=512, isQuery=false, resetVect=false, bitFlags=avalonSimilarityBits
```

The selected COSMolKit API must follow the pinned Python wrapper for Python
defaults and must model the source fields `nBits`, `isQuery`, and `bitFlags`.
`resetVect` is an adapter-internal concern for a newly allocated result and is
not a useful public option on a value-returning COSMolKit method. Count
fingerprints, string-input overloads, canonical SMILES, structure checking,
and coordinate generation are outside this selected bit-vector port.

The wrapper's molecule path is:

```text
GetAvalonFP(ROMol, ...)
  -> AvalonTools::getAvalonFP(ROMol, ExplicitBitVect&, ...)
  -> molToReaccs(ROMol)
  -> MolToMolBlock(ROMol, true)
  -> LocaleSwitcher
  -> MolStr2Mol()
  -> getFp()
  -> SetFingerprintBits()
  -> SetFingerprintCountsWithFocus()
```

For non-query molecules, the adapter performs a second accumulation pass with
`ACCUMULATE_BITS | USE_DY_AROMATICITY`. The adapter rounds its internal byte
buffer to a multiple of four before word packing, while the public `nBits`
conversion uses `nBits / 8` and warns for a non-byte-aligned size. These are
source behavior, not opportunities for normalization or fallback.

## Engine license and redistribution

The Avalon archive root `license.txt` grants source and binary redistribution
under a BSD-like license attributed to Copyright 2001-2011 Novartis Pharma AG.
It requires retaining the copyright, conditions, and disclaimer in source and
binary distributions and prohibits endorsement by the named holders.

The vendored RDKit adapter is covered by the RDKit BSD-3-Clause license in
`third_party/rdkit/license.txt`. The adapter source files also carry RDKit
copyright notices. The active Avalon C implementation consists of the 26 C
files selected in RDKit's `External/AvalonTools/CMakeLists.txt`:

```text
common/layout.c             common/symboltable.c       common/patclean.c
common/utilities.c          common/symbol_lists.c      common/stereo.c
common/set.c                common/perceive.c          common/local.c
common/graph.c              common/geometry.c          common/forio.c
common/depictutil.c         common/denormal.c          common/casutils.c
common/ssmatch.c            common/rtutils.c           common/smi2mol.c
common/didepict.c           common/pattern.c           common/canonizer.c
common/aacheck.c             common/fixcharges.c        programs/struchk.c
common/reaccsio.c            common/hashcode.c
```

The reachable headers are under `src/main/C/include`. The source tree carries
the same Novartis terms on the normal C and header files; `casutils.c` has an
older historical file header without a repeated license block, and
`programs/struchk.c` places the Novartis block after its historical header.
Both remain covered by the archive root license and must not be omitted from
the notice audit. The port must preserve both the Avalon notice and the RDKit
notice in the repository and in any binary redistribution metadata.

## Current COSMolKit mismatch

`crates/cosmolkit-core/src/properties/avalon_fingerprint.rs` currently exposes
invented path-fingerprint fields (`min_path`, `max_path`,
`n_bits_per_hash`, `use_bond_order`, `use_hs`, `tautomeric_fingerprint`, and
`from_atoms`) and fails closed. `python/src/lib.rs` mirrors that wrong shape.
Those fields are not part of `GetAvalonFP(ROMol, ...)`; they must be replaced
by a typed bit-flag mask plus `nBits` and `isQuery` before the implementation
is enabled. No approximate path hash, compatibility shim, or molecule-specific
fallback is permitted.

## Porting consequence

The next ordered tasks are to materialize this exact archive under ignored
`tmp/parity-audit/sources/avalon/`, inventory the complete active call graph
and line ranges, then correct the public parameter types while retaining the
structured unsupported result. Only after those artifacts and focused tests
exist may the C source bodies be ported line by line with two-axis markers.
