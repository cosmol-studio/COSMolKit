# cosmolkit-inchi

The primary user-facing API is documented in [`cosmolkit`](../cosmolkit/README.md).

`cosmolkit-inchi` is a toolkit-neutral, pure Rust source port of the official
IUPAC InChI engine for the four scalar operations used by COSMolKit:

- `mol_to_inchi`
- `mol_to_inchi_key`
- `inchi_to_inchi_key`
- `mol_from_inchi`

The crate owns neutral atom, bond, coordinate, isotope, charge, hydrogen, and
stereo representations. It does not depend on `cosmolkit-core`; conversion
between a COSMolKit `Molecule` and the neutral InChI graph belongs to the core
crate.

## Examples

The direct InChIKey operation needs no toolkit adapter:

```bash
cargo run -p cosmolkit-inchi --example inchi_to_key
cargo run -p cosmolkit-inchi --example inchi_to_key -- 'InChI=1S/H2O/h1H2'
```

The neutral graph operations use downstream toolkit traits for element,
valence, sanitization, hydrogen, and stereochemistry behavior. These examples
include an explicitly methane-only adapter that rejects every operation outside
that narrow boundary:

```bash
cargo run -p cosmolkit-inchi --example neutral_graph_to_inchi
cargo run -p cosmolkit-inchi --example inchi_to_neutral_graph
```

## Behavior Boundary

The four public operations reproduce pinned official InChI v1.07.5 and RDKit
2026.03.1 behavior exactly where the official C source defines behavior. The
audited `NormalizeAndCompare` initial-buffer allocation-failure path executes
undefined C behavior; the Rust API returns a deterministic structured
allocation error for that path.

Generation and parsing return the source-defined status, diagnostics, message,
log, AuxInfo, graph, and key fields represented by their public output types.
Nonfatal diagnostics remain explicit and source failures are not converted into
empty or chemically plausible fallback results.

## Production Constraints

Production code uses only owned Rust data and the Rust allocator. It does not:

- link or dynamically load an official/system InChI library;
- invoke external executables or subprocesses;
- route through SMILES, MolBlock, or SDF as an implementation substitute;
- use RDKit, Open Babel, Python, or PyO3 at runtime;
- silently fall back when source behavior is unavailable.

Official C is compiled only by independent test-oracle fixtures and never enters
the production dependency graph.

The live official-InChI and RDKit C/C++ oracle tests are explicitly ignored by
ordinary `cargo test` runs because their pinned upstream source trees are
development inputs, not repository test dependencies. When those source trees
are present, run a named oracle test with `--ignored --exact`; the test still
performs the complete byte/field comparison and fails on any mismatch. For
example:

```bash
cargo test -p cosmolkit-inchi --release \
  official_c_oracle__getinchikeyfrominchi__exact -- --ignored --exact

cargo test -p cosmolkit-inchi --release \
  rdkit_cpp_oracle__inchitoinchikey__exact -- --ignored --exact
```

## Validation Status

Focused source-port tests cover the active official-engine call graph, and
independent official-C oracle tests compare observable output fields at the
function boundaries. The public COSMolKit adapter is also compared exactly
with pinned RDKit output on the committed small corpus and the 5000-row InChI
corpus. Tests do not ignore fields, use compatibility thresholds, or route
production calls through the C oracle.

## Source Alignment

The official engine source is pinned to IUPAC InChI `v1.07.5`:

- upstream repository: <https://github.com/IUPAC-InChI/InChI>
- peeled commit: `11a87982bb518f57ac013f0b258c283655e1ea1d`
- commit tree: `214d60b92a94b03508a6b262590262bb4f89e275`
- source archive SHA-256:
  `88532b3f599d125940e91af5d3135b31b0392b4c5a6e3f25e6418d9b56c5d5e3`
- local source path: `third_party/InChI`
- approved production target:
  `INCHI-1-SRC/INCHI_API/libinchi/src/CMakeLists.txt`

The COSMolKit adapter is aligned with RDKit `Release_2026.03.1`, independently
traceable to `third_party/rdkit/External/INCHI-API/inchi.cpp` and `inchi.h`.
Official engine and RDKit adapter source frames remain separately attributed.

Every source-backed Rust function retains its complete C/C++ source frame and
two-axis reproduction marker according to
[`../../dev/source_reproduction_protocol.md`](../../dev/source_reproduction_protocol.md).
