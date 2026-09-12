# RDKit XYZ Runtime Adapter Migration

## Outcome

`cosmolkit-io::read_xyz_detached()` is now the only production XYZ parser.
The compatibility entry point in `cosmolkit::io::xyz` contains only:

1. delegation to the detached reader;
2. installation of the returned topology, coordinate, and property blocks;
3. assignment of `TopologyTrust::CoordinateOnly`; and
4. structured translation of parse and installation errors.

The former duplicate runtime tokenization, unsigned-count conversion, element
normalization, coordinate parsing, and detached graph assembly were removed.

## Contract Preservation

The runtime adapter retains the existing public
`read_xyz_from_str() -> Result<Molecule, XyzReadError>` boundary. It does not
perform chemistry finalization, infer bonds, mutate caches, or reinterpret
parser errors. The detached owner retains the verbatim RDKit source anchors and
two-axis status markers for `MolFromXYZBlock`, `MolFromXYZDataStream`,
`ParseXYZFileAtomLine`, and `FileParserUtils` numeric conversion.

## Focused Validation

Both debug and release strict-contract runs passed all seven retained runtime
tests:

```text
cargo test -p cosmolkit --features op-contracts-strict xyz_reader
7 passed, 0 failed

cargo test -p cosmolkit --release --features op-contracts-strict xyz_reader
7 passed, 0 failed
```

The cases cover normal graph/coordinate/property installation, uppercase
two-letter element normalization, extra-field rejection, scientific-notation
rejection, RDKit unsigned-count whitespace and prefix behavior, unsigned
overflow-to-zero behavior, and bounded allocation for a truncated maximum
count record.

## Remaining IO Work

Typed V2000/V3000 SGroup and enhanced-stereo collection lowering now belongs to
`cosmolkit-io`; see `rdkit_sdf_typed_sgroup_migration.md`. The query-bearing
MolBlock result contract is also detached; see
`rdkit_sdf_query_graph_contract.md`. The legacy runtime MolBlock/SDF
implementation must not be deleted until the remaining parameter,
property-list, finalization, and supplier behavior has detached ownership.
