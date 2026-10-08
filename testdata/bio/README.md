# Structural-Biology Test Data

Committed structural inputs live in `fixtures/`. The Gemmi mmCIF writer
profile is `gemmi_mmcif_writer_profile.json`; it binds the selected fixtures,
input formats, writer arguments, and expected output names to Gemmi commit
`5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e` (Gemmi 0.7.5).

This profile records the original fixture provenance, not the 0.5.0
preparation procedure. Current BIO corpus and special-regression preparation
use the [standard test runner](../../parity-tests_fixed/README.md), with pinned
Gemmi Python bindings and package-owned inputs and generated references.
Do not compile a separate C++ oracle for that workflow.

The profile covers a PDB-to-mmCIF path and a represented mmCIF rewrite path.
Together they exercise category ordering, CIF quoting and serialization, atom
and anisotropic output, secondary structural categories, connections,
assemblies, crystallographic transforms, and metadata categories.
