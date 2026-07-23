# Official InChI Production Type Inventory

## Configuration

- Upstream release: `v1.07.5`
- Upstream commit: `11a87982bb518f57ac013f0b258c283655e1ea1d`
- Production target: official `libinchi`
- Production C translation units: `60`
- Audited production headers in the target wrapper: `41`
- Preprocessor configuration: `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`
- Audited ABI arithmetic profile: GCC/Clang Linux LP64
- Declaration baseline generator: `bindgen-cli 0.72.1`
- Rust declaration transformer: `dev/inchi-source-type-generator`
- Reproduction driver: `dev/generate_inchi_source_types.py`
- Owned Rust output:
  `crates/cosmolkit-inchi/src/source_types/generated.rs`

The driver processes the shared production headers and then every production
C translation unit separately. Header declarations occupy the shared source
type module. C-file-local declarations occupy a module named for their
translation unit. This is required because macro values can differ by
translation unit; for example, `RESET_EDGE_FORBIDDEN_MASK` is configured as
both `0` and `1` in different source scopes. Conflicting local values are not
flattened or selected heuristically.

## Generated Coverage

The configured output contains:

- `183` unique generated structure names;
- `2` source unions represented manually by owned bit-semantic Rust types;
- `1` function-local/header-defined `stbsp__context` represented manually;
- `262` generated type-alias declarations across shared and local scopes;
- `40` configured source enum representations, with every enumerator retained
  as its exact integer constant;
- `1,837` generated constant declarations across shared and translation-unit
  scopes;
- `37` non-empty translation-unit-local declaration modules;
- the configured callback typedefs, nullable callback variables, and callback
  fields as safe Rust function shapes.

The other 23 translation units define functions and data but no additional
local type or macro declarations after the shared-header declarations are
excluded.

## Ownership Translation

The generated output is layout-independent and contains no `repr(C)`, FFI
block, native C function call, native raw pointer, or `unsafe` block.

- C signed/unsigned integer spelling is fixed to the audited LP64 width.
- C pointer fields become nullable `SourceMutPointer<T>` or
  `SourceConstPointer<T>` values containing an allocation identity and an
  element offset.
- `SourceHeap` owns typed allocations and preserves null, alias, interior,
  one-past, dangling-after-free, and type-mismatch states without exposing a
  native address.
- Fixed arrays retain the exact configured source length.
- Every generated structure has a safe explicit zero/empty `Default`; arrays
  of any source length are initialized recursively without `zeroed()`.
- C enum storage remains a raw signed or unsigned integer alias so unknown
  values and exact error paths are preserved. Validated public enums will be a
  separate facade boundary.
- `union tagSplitLong` and `union BnsAltPath` retain their configured
  little-endian bit interpretation through explicit accessors.

## Cross-Checks

The independently compiled official `libinchi.so` debug information was
queried for configured structs and enums. Every source-owned structure visible
in DWARF is present in the Rust inventory. The three deliberately absent DWARF
structures are:

- `_IO_FILE`, replaced by the owned byte-oriented `SourceFile`;
- `__va_list_tag`, replaced by `SourceVaList` and typed format arguments;
- `stbsp__context`, present as a manual owned Rust structure because its C
  declaration is local to a header-defined function section and was not
  emitted by the header declaration baseline.

The per-translation-unit extraction additionally found source declarations
not retained in optimized DWARF, including `tagNodeValues`; these are present
in the Rust inventory. GCC syntax-checks the declaration wrapper under the same
defines, and `cargo check -p cosmolkit-inchi` compiles the owned Rust result.

## Regeneration

Install or provide `bindgen-cli 0.72.1`, then run from the repository root:

```bash
BINDGEN=/path/to/bindgen \
LIBCLANG_PATH=/usr/lib/llvm-18/lib \
.venv/bin/python dev/generate_inchi_source_types.py
```

The script rejects a different bindgen version and rejects a production C-file
count other than 60. The AST transformer rejects conflicting declarations
within one configured scope and asserts that no raw pointer or C ABI remains
in its output.
