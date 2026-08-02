# `cosmolkit-inchi` Crate-Boundary Audit

## Audit Scope

This report audits the existing `crates/cosmolkit-inchi` skeleton against the
approved reusable-crate boundary in `dev/archive/plans/rdkit_inchi_full_port_plan.md`. It
records the state before the module, type, allocation, engine, adapter, and
public API implementation steps. It does not treat planned behavior as already
implemented.

Audit date: `2026-07-19`.

## Approved Boundary

`cosmolkit-inchi` is the toolkit-neutral, pure-Rust source-port boundary for
the official IUPAC InChI engine and the narrow RDKit InChI adapter behavior.
The completed crate owns toolkit-neutral atom, bond, coordinate, isotope,
radical, hydrogen, stereo, status, diagnostics, option, AuxInfo, InChIKey,
engine, and adapter behavior.

The crate must not own or depend on:

- COSMolKit `Molecule` or `cosmolkit-core`;
- the molecule operation registry, molecule cache, Python, or PyO3 types;
- RDKit binaries or runtime libraries;
- Open Babel;
- a native or system InChI library in production;
- subprocesses, external executables, or network access;
- SMILES or MolBlock regeneration as a substitute for InChI behavior;
- `crates/cosmolkit-core-old`.

Conversion between COSMolKit `Molecule` and the toolkit-neutral input belongs
in `cosmolkit-core`, outside this crate. A native official-C oracle may exist
only as non-production test tooling and may not become the implementation
under test.

## Existing Skeleton

The crate currently contains only:

- `crates/cosmolkit-inchi/Cargo.toml`;
- `crates/cosmolkit-inchi/README.md`;
- `crates/cosmolkit-inchi/src/lib.rs`.

The manifest uses the workspace Rust 2024 edition and version, sets
`publish = false`, declares `MIT OR Apache-2.0`, and has no dependencies. The
crate is registered as a workspace member and appears in `Cargo.lock`.
`cargo check -p cosmolkit-inchi` succeeds, and
`cargo tree -p cosmolkit-inchi --edges normal` reports no normal dependency.

No other current core, facade, Python, or test source depends on
`cosmolkit-inchi`. The library root contains documentation only. It declares
no modules, types, constants, traits, callbacks, functions, or feature flags
and exposes no InChI generation, parsing, conversion, AuxInfo, or InChIKey
API.

## Conformance

The skeleton currently satisfies these structural requirements:

- it is a separate reusable workspace crate;
- it has no dependency on `cosmolkit-core`, RDKit, Python, PyO3, Open Babel,
  native libraries, subprocess tooling, network clients, SMILES, or MolBlock
  conversion;
- it does not expose an approximate, placeholder, FFI-backed, or
  regeneration-backed production API;
- its library documentation states that public generation and parsing remain
  unavailable until the selected source call graph is ported and tested;
- it distinguishes the official engine port from COSMolKit integration.

These are boundary properties only. They do not establish any implemented
InChI behavior or source parity.

## Gaps

The approved crate boundary is not yet implemented beyond the empty shell:

- there is no module map for source types, memory/allocation behavior, IO,
  normalization, canonicalization, tautomerism, stereochemistry, parsing,
  serialization, AuxInfo, InChIKey, engine entry points, or RDKit adapter
  behavior;
- there are no owned Rust equivalents for production `libinchi` structs,
  unions, enums, constants, callback shapes, or fixed-width integer
  semantics;
- there is no explicit allocation and ownership policy mapping C allocation,
  reallocation, freeing, pointer aliasing, nullability, and buffer lifetime to
  Rust;
- there is no implemented `BEGIN INCHI C FUNCTION` / `END INCHI C FUNCTION`
  framing policy in Rust functions and no `INCHI` two-axis marker-bearing
  source frames;
- there is no test harness, official-C oracle tooling, provenance-checked
  golden schema, source-derived focused test, or parity suite in this crate;
- there is no engine or adapter implementation and intentionally no public
  production API;
- there is no `cosmolkit-core` conversion layer yet; this remains external to
  the reusable crate and must not be added here.

The README's source-alignment section is stale. It says the official release,
commit, checksum, and license approval are unpinned and that
`third_party/InChI` is absent. The approved and vendored source is now IUPAC
InChI `v1.07.5`, annotated tag object
`fa28fcb6fbba554c952ad321b5d7797a5fb001a5`, peeled commit
`11a87982bb518f57ac013f0b258c283655e1ea1d`, archive SHA-256
`88532b3f599d125940e91af5d3135b31b0392b4c5a6e3f25e6418d9b56c5d5e3`,
and local path `third_party/InChI`. MIT vendoring was approved, and the
upstream license is present at `third_party/InChI/LICENSE`.

The README also leaves the RDKit adapter commit unspecified. The source
inventory records the local adapter baseline as the vendored RDKit
`2026.03.1` tree. Adapter source remains independently traceable and must use
`RDKit` markers rather than `INCHI` markers.

## Required Follow-On Ownership

The next architecture step must update the crate README with the pinned
provenance, complete module map, owned type boundary, allocation policy, and
source-frame policy. Later function-level steps must implement each inventoried
production function with its complete source frame and focused test before a
production API is exposed.

Integration must preserve the following ownership split:

| Concern | Owner |
|---|---|
| Official InChI source-derived types and algorithms | `cosmolkit-inchi` |
| Toolkit-neutral input/output, diagnostics, options, AuxInfo, InChIKey | `cosmolkit-inchi` |
| Source-derived RDKit adapter behavior over toolkit-neutral data | `cosmolkit-inchi` |
| COSMolKit `Molecule` conversion and operation integration | `cosmolkit-core` |
| Rust facade exposure | `cosmolkit` |
| Python/PyO3 exposure | `python` / `cosmolkit-py` |
| Official C execution used to create or verify goldens | non-production test tooling only |

## Audit Result

The existing crate is correctly isolated and appropriately refuses to expose
unimplemented behavior. It is nevertheless only a documentation skeleton.
All functional parts of the approved boundary remain explicit downstream
work; none may be replaced by a heuristic, fallback, native production link,
or external process.
