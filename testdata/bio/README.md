# Structural-Biology Test Data

Committed structural inputs live in `fixtures/`. The Gemmi mmCIF writer
profile is `gemmi_mmcif_writer_profile.json`; it binds the selected fixtures,
input formats, writer arguments, and expected output names to Gemmi commit
`5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e` (Gemmi 0.7.5).

Prepare exact writer output with:

```bash
.venv/bin/python tools/testdata/gemmi/generate_all.py \
  --profile bio_mmcif_writer \
  --suite mmcif_writer \
  --jobs 4
```

The generator verifies hashes of the pinned Gemmi reader, writer, serializer,
build, and version sources before compiling a narrow source-owned oracle under
`target/`. The oracle directly calls `read_structure_file`,
`make_mmcif_document`, and the default `write_cif_to_stream`; it does not add
CLI-specific serialization options. The generator does not infer source
identity from an installed package and it does not invoke Git. Generated
output and `manifest.json` are written under
`testdata/bio/expected/gemmi/bio_mmcif_writer/`; that directory is ignored and
must be regenerated when any declared identity changes.

The profile covers a PDB-to-mmCIF path and a represented mmCIF rewrite path.
Together they exercise category ordering, CIF quoting and serialization, atom
and anisotropic output, secondary structural categories, connections,
assemblies, crystallographic transforms, and metadata categories.
