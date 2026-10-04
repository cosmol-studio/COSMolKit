# Parity task table

One task tests one chemistry operation. A parameter profile is an explicit
argument combination, not a corpus size, language, repetition count or old
audit phase. The reference is RDKit 2026.03.1.

The Rust definitions live under `registry::molecule_plan`. They are a design
table, **not passing evidence**. The parent registry activates the two fuzzy
tasks and nine molecular tasks through typed input/result adapters, reference
preparation and full-selection preflight.
Activation must consume these definitions, not create another matrix in Python.

## Category map

Categories describe behavior, not implementation readiness. Concrete molecules,
queries, fingerprint values and biological structures remain distinct inputs.
The initial 12 molecular tasks below are the first frozen matrices, not the
boundary of the project.

| Category | Initial matrix already defined | Further task scopes |
|---|---|---|
| Notation | SmilesRead | `smiles_write`, `cxsmiles_read`, `cxsmiles_write`, `smarts_read`, `smarts_write` |
| MolecularIo | — | `molblock_read`, `molblock_write`, `sdf_read`, `sdf_write`, `mol2_read`, `binary_roundtrip` |
| Chemistry | Sanitize, Kekulize, AddHydrogens, RemoveHydrogens, Valence | `sanitize_stages`, `aromaticity`, `radicals`, `conjugation`, `hybridization`, `cleanup`, `rings_fast`, `rings_sssr`, `rings_symmetrized`, `ring_families`, `fragments`, `distance_matrix`, `hydrogen_options` |
| Stereo | CipLabels, PotentialStereo | `structure_stereo`, `stereo_enumeration` |
| Descriptors | MolecularWeight, ExactMolecularWeight, MolecularFormula | `connectivity`, `lipinski_counts`, `mqn`, `crippen`, `surface_area` |
| FingerprintGeneration | — | `morgan`, `maccs`, `rdk_fingerprint`, `avalon`, `atom_pair`, `topological_torsion`, `pattern`, `layered` |
| FingerprintValues | fuzzy_and, fuzzy_or (executable) | `fingerprint_arithmetic`, `fingerprint_similarity`, `fingerprint_serialization` |
| Query | — | `substructure_match`, `compiled_query`, `generic_groups`, `mcs` |
| Depiction | Coordinates2d | `depiction_options`, `depiction_normalize`, `svg`, `png` |
| Conformers | — | `bounds`, `embedding`, `conformer_pruning` |
| ForceFields | — | `uff_parameters`, `uff_energy_gradient`, `uff_optimize`, `mmff_parameters`, `mmff_energy_gradient`, `mmff_optimize` |
| Alignment | — | `alignment`, `best_rmsd` |
| Tautomers | — | `tautomer_enumeration`, `tautomer_canonicalization` |
| Identifiers | — | `inchi_generate`, `inchi_read`, `inchi_key` |
| Bio | — | `pdb_read`, `mmcif_read`, `bio_write`, `bio_selection`, `bio_sequence`, `bio_geometry` |
| Composition | — | `batch`, `operation_sequences`, `concurrent_reads` |
| Native | — | `confseq` |

Readiness is separate from category:

- **Executable:** two fuzzy tasks and nine molecular tasks (20 profiles):
  SmilesRead, Sanitize, Kekulize, MolecularWeight, ExactMolecularWeight,
  MolecularFormula, AddHydrogens, RemoveHydrogens and Coordinates2d.
- **Matrix defined, not wired:** CipLabels (4), PotentialStereo (8), Valence (2).
  CIP needs one preparation-time resolution of selected indices shared by both
  executors. Potential stereo needs a reference adapter for the complete output,
  including ranks/relations, not just Python's stereo-info subset. Valence needs
  a public result readout. None is silently reduced to a weaker comparison.
- Molecular error rows currently remain explicit non-passes with stage and
  diagnostics. Matching rejection/text does not establish error-category parity.
- **Future scope, matrix not frozen:** the 76 named rows below. Some owners
  already implement parts; others need implementation or a public projection.
  This label is a test-integration state, not a claim that chemistry is absent.

The future scopes are also declared in Rust as `FUTURE_TASKS`. Their count is
`None`, not zero: enumerating parameter names is not a completed matrix.
Promotion requires exact typed values and combinations, a pinned reference,
a public callable/result boundary, error mapping and executable adapters.
Do not multiply every named option blindly or adopt old profile counts without
recovering the actual corresponding options.

## Future task scopes by category

These scopes are included now so that missing public APIs or unfinished ports
do not disappear from the testing plan. They do not authorize production work
outside the sole architecture execution plan.

### Notation

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `smiles_write` | canonical, isomeric, kekule, explicit bonds/H, root and random traversal | Text and atom/bond output order | Rdkit |
| `cxsmiles_read` | CX fields, strictness, names, coordinate records | Typed topology, query predicates, properties and coordinates | Rdkit |
| `cxsmiles_write` | CX field mask, coordinate selector, canonical/isomeric options | Text and output mappings; approved coordinate-selection difference explicit | Rdkit |
| `smarts_read` | merge H, replacements, CX/name policy | Ordered query graph and predicate trees | Rdkit |
| `smarts_write` | root, isomeric and CX fields | Exact text and query semantics | Rdkit |

### MolecularIo

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `molblock_read` | V2000/V3000, sanitize, remove H, strict parsing, coordinate mode | Concrete/query result, graph, stereo, SGroups, coordinates and errors | Rdkit |
| `molblock_write` | V2000/V3000, stereo, kekule, coordinate selector | Text and preserved graph fields | Rdkit |
| `sdf_read` | record boundaries, invalid records, duplicate data fields, reader options | Ordered records, errors, data fields and graph category | Rdkit |
| `sdf_write` | record/data-field order, molblock options | Serialized records and fields | Rdkit |
| `mol2_read` | sanitize, remove H, substructure cleanup | Topology, atom types, charge and outcome | Rdkit |
| `binary_roundtrip` | serialization version, property/coordinate selection | Before/after CK public state; no RDKit byte-compatibility claim | Roundtrip |

### Chemistry

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `sanitize_stages` | named single stages, prerequisite sequences and approved stage combinations | Stage result and structured failure; separate from initial ALL profile | Rdkit |
| `aromaticity` | each supported model, prepared kekule input | Atom/bond aromatic flags and bond orders | Rdkit |
| `radicals` | raw versus prepared input | Per-atom radical counts and errors | Rdkit |
| `conjugation` | bond types and prepared chemistry state | Per-bond conjugation | Rdkit |
| `hybridization` | prepared valence/aromaticity state | Per-atom hybridization | Rdkit |
| `cleanup` | ordinary, organometallic and atropisomer cleanup as separate profiles | Topology changes and outcome | Rdkit |
| `rings_fast` | uninitialized and existing ring state | Ring rows, membership and find type | Rdkit |
| `rings_sssr` | dative/hydrogen bond inclusion | Ordered rings and membership | Rdkit |
| `rings_symmetrized` | dative/hydrogen bond inclusion | Symmetrized rings and membership | Rdkit |
| `ring_families` | bond inclusion, fused/spiro/bridged inputs | Family rows, counts and relevant cycles | Rdkit |
| `fragments` | index groups versus molecule outputs, sanitize fragments | Components, source mappings and output graphs | Rdkit |
| `distance_matrix` | bond-order weighting, atom weights, topological versus supplied 3D | Dimensions and every matrix entry | Rdkit |
| `hydrogen_options` | each remaining AddHs/RemoveHs option and source-required interactions | Topology, isotopes, mapping, coordinates and stereo | Rdkit |

### Stereo

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `structure_stereo` | 2D/3D selection, replace existing tags | Atom tags and bond stereo | Rdkit |
| `stereo_enumeration` | unassigned-only, unique, enhanced groups, embedding, limits and fixed seed | Enumerated structures, ordering and completion state | Rdkit |

### Descriptors

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `connectivity` | Chi families/orders, Hall-Kier, Kappa, Phi | Scalar descriptor values | Rdkit |
| `lipinski_counts` | HBD/HBA, rotatable-bond definitions, atom/ring classifications | Exact integer counts | Rdkit |
| `mqn` | complete component vector | Every vector element | Rdkit |
| `crippen` | hydrogen handling, per-atom contributions | LogP/MR totals and contribution arrays | Rdkit |
| `surface_area` | Labute hydrogen policy; VSA family and default/custom bins | Scalar area, bins and contributions | Rdkit |

### FingerprintGeneration

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `morgan` | radius, bit/count/sparse output, chirality, features, invariants, roots, provenance | Vector and complete additional output | Rdkit |
| `maccs` | standard key definition | Exact bit vector | Rdkit |
| `rdk_fingerprint` | path limits, branched paths, H, bond order, size, roots, provenance | Vector and atom/bond provenance | Rdkit |
| `avalon` | size, query mode and supported flag profiles | Exact bit vector | Rdkit |
| `atom_pair` | distance range, 2D/3D, chirality, roots/ignored atoms, bit/count/sparse output | Vector and additional output | Rdkit |
| `topological_torsion` | torsion length, chirality, roots, count simulation, output representation | Vector and additional output | Rdkit |
| `pattern` | size, tautomer mode, atom counts and set-only mask | Vector and updated counts | Rdkit |
| `layered` | layer mask, paths, roots, count seeds, set-only mask | Vector and updated counts | Rdkit |

### FingerprintValues

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `fingerprint_arithmetic` | representation/index width, binary/scalar operation, signed/zero counts | Length, exact entries and error | Rdkit |
| `fingerprint_similarity` | named metric, count/bit semantics, empty vectors, metric parameters | Scalar result and error | Rdkit |
| `fingerprint_serialization` | supported text/binary form and index width | Bytes/text and reconstructed entries | Rdkit |

### Query

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `substructure_match` | chirality, query-query, recursion, uniqueness, max matches, properties, callbacks | Boolean, ordered mappings and errors | Rdkit |
| `compiled_query` | same matching profiles as ordinary matching | Same result as reference and ordinary CK path | Rdkit |
| `generic_groups` | each supported generic label and generic-matcher option | Match results | Rdkit |
| `mcs` | atom/bond comparators, ring/stereo constraints, objective, seed, threshold, limits, ties | Counts, query, SMARTS, degenerate results and completion | Rdkit |

### Depiction

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `depiction_options` | orientation, templates, constrained maps, sampling/seed, mimic distances | Coordinates and source topology | Rdkit |
| `depiction_normalize` | normalization, straightening, conformer selection and scaling | Coordinates and returned transform/scale | Rdkit |
| `svg` | labels, highlights, stereo, SGroups, annotations, dimensions | Declared SVG structure/text boundary; approved CK branding explicit | Rdkit |
| `png` | same drawing profiles, image size and renderer | Declared decoded pixel/rendering boundary, not unexplained byte equality | Rdkit |

### Conformers

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `bounds` | bond policy, smoothing, macrocycle and fragment settings | Every bound entry and outcome | Rdkit |
| `embedding` | ETKDG version, count, fixed seed, chirality, random coordinates, constraints | Outcome, conformer count, coordinates and diagnostics | Rdkit |
| `conformer_pruning` | RMS threshold, symmetry, heavy atoms and terminal groups | Retained conformers and ordering | Rdkit |

### ForceFields

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `uff_parameters` | element/bond environments and interfragment options | Parameter availability and values | Rdkit |
| `uff_energy_gradient` | fixed input geometry and constraints | Energy and gradient vector | Rdkit |
| `uff_optimize` | iteration limit, tolerances, constraints and conformer selection | Status, final energy and coordinates | Rdkit |
| `mmff_parameters` | MMFF94/MMFF94s and parameter families | Types, charges, parameter availability and values | Rdkit |
| `mmff_energy_gradient` | variant, fixed geometry and constraints | Energy and gradient vector | Rdkit |
| `mmff_optimize` | variant, iteration limit, tolerances and constraints | Status, final energy and coordinates | Rdkit |

### Alignment

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `alignment` | mapping, weights, reflection, conformer IDs and iteration limit | RMSD, transform and coordinates | Rdkit |
| `best_rmsd` | symmetry, terminal groups, mapping limits and H policy | Best mapping and RMSD | Rdkit |

### Tautomers

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `tautomer_enumeration` | default/v1 transforms, enumeration limits, stereo/isotope retention and reassignment | Ordered results, modified atoms/bonds and completion | Rdkit |
| `tautomer_canonicalization` | scoring profiles and enumeration policy | Scores and canonical selected structure | Rdkit |

### Identifiers

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `inchi_generate` | standard/nonstandard, stereo/isotope/fixed-H options | InChI, AuxInfo, return status and messages | Rdkit |
| `inchi_read` | sanitize, remove H and supported parse options | Molecule state and errors | Rdkit |
| `inchi_key` | valid/invalid InChI and supported key modes | Exact key and outcome | Rdkit |

### Bio

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `pdb_read` | models, altloc, connectivity, element inference and metadata | Typed hierarchy, source IDs, coordinates and metadata | Gemmi |
| `mmcif_read` | models, altloc, author/label IDs, entities, assemblies and missing values | Typed hierarchy, identifiers, coordinates and metadata | Gemmi |
| `bio_write` | PDB/mmCIF, model/altloc selection and formatting policy | Source-backed preserved fields and roundtrip | Gemmi |
| `bio_selection` | protein/DNA/RNA/water/ligand, chains, residues, models and altloc | Selected original IDs and topology | Gemmi |
| `bio_sequence` | protein/DNA/RNA, modified residues, missing/unknown residues | Sequence, residue classification and source mapping | Gemmi |
| `bio_geometry` | distances, neighbors, crystal/assembly transforms and selection | Coordinates, contacts and mappings | Gemmi |

### Composition

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `batch` | same chemistry profiles, invalid-record policy, ordering and worker count | Every indexed result/error versus scalar Rust | RustScalar |
| `operation_sequences` | explicit source-backed operation sequences | State after each stage, not only the last result | Rdkit |
| `concurrent_reads` | declared read operations and fixed schedules | Same outputs as scalar Rust, unchanged input | RustScalar |

### Native

| Task | Parameter axes to freeze | Compare | Reference |
|---|---|---|---|
| `confseq` | native sequence/geometry parameters | Local correctness and binding consistency; not RDKit parity | NoExternalReference |

### Reference boundaries

- `Rdkit`: the pinned RDKit chemistry interface, including explicitly
  registered CK differences. For direct InChI engine ABI-level comparison,
  freeze the official InChI reference separately; do not mix its fields with
  RDKit wrapper behavior.
- `Gemmi`: a source-backed Gemmi behavior with an exact version/commit to be
  frozen per task. Bio behavior sourced from another library needs its own
  reference profile; this label does not claim all CK Bio APIs exist in Gemmi.
- `RustScalar`: batch/concurrency consistency against the same Rust scalar
  operation, not independent proof of chemistry parity.
- `Roundtrip`: preservation of declared CK state, not RDKit binary format
  compatibility.
- `NoExternalReference`: native behavior such as ConfSeq. Keep local
  correctness regressions and binding consistency checks, never label it
  RDKit parity. It appears here to make that boundary explicit, not as a new
  corpus-parity chemistry task.

Performance and language bindings are cross-cutting validation dimensions,
not chemistry categories. A performance task needs a named reference baseline,
matching inputs/options and a declared metric. Python/JS checks consume the
same Rust task IDs and profile identities rather than define duplicate tasks.
Corpus scale and thread schedules do not create new chemistry features.

## Initial frozen parameter matrices

`false/true` means both values; `x` means the full Cartesian product.

| Task | Chemistry parameters | Profiles per input | Result to compare |
|---|---|---:|---|
| `fuzzy_and` | index width: u32/u64 | 2 | Exact length and ordered sparse entries |
| `fuzzy_or` | index width: u32/u64 | 2 | Exact length and ordered sparse entries |
| `SmilesRead` | sanitize false/true x remove_hydrogens false/true | 4 | Parse outcome and topology |
| `Sanitize` | operations = ALL | 1 | Outcome, failing stage and resulting topology |
| `Kekulize` | clear_aromatic_flags false/true | 2 | Outcome, bond orders and aromatic flags |
| `MolecularWeight` | only_heavy false/true | 2 | f64 bits |
| `ExactMolecularWeight` | only_heavy false/true | 2 | f64 bits |
| `MolecularFormula` | separate_isotopes false/true x abbreviate_h_isotopes false/true | 4 | Exact formula text |
| `AddHydrogens` | explicit_only false/true | 2 | Resulting atom/bond rows and stereo state |
| `RemoveHydrogens` | sanitize false/true | 2 | Outcome and resulting atom/bond rows and stereo state |
| `CipLabels` | Four explicit profiles below; not a product | 4 | Atom/bond CIP labels, absence and outcome |
| `PotentialStereo` | clean false/true x flag_possible false/true x allow_nontetrahedral false/true | 8 | Ordered stereo records and cleaned topology |
| `Coordinates2d` | One explicit baseline below | 1 | Coordinate shape, coordinates and preserved topology |
| `Valence` | strict false/true; model = RdkitLike | 2 | Explicit valence and implicit-H rows, violation and outcome |

There are 12 planned molecular tasks and 34 chemistry profiles per eligible
input. These are initial named boundaries, not exhaustive coverage of every
upstream overload. In particular, one ALL sanitize profile does not cover
every operation-bit subset, and two hydrogen profiles do not cover every
RemoveHs option. Add further complete named profiles explicitly.

### CIP selections

| Profile | atoms | bonds | max_recursive_iterations |
|---|---|---|---:|
| All | None | None | 0 |
| FirstTaggedAtom | first atom with a non-unspecified chiral tag, or empty | None | 0 |
| FirstStereoBond | None | first bond with non-NONE stereo, or empty | 0 |
| Empty | empty | empty | 1 |

Choose the lowest original index satisfying each predicate. Preparation must
record the explicit selected indices once; pass the same indices to CK and
RDKit. Do not select independently from potentially divergent parsed results.
None means all; empty means none. Missing eligible centers do not skip a case.
The empty/limit-one profile reproduces the old matrix; it does NOT prove
iteration-limit exhaustion. A nonempty exhaustion profile is separate work.

## Input preparation and fixed options

- `SmilesRead` consumes original text. `Sanitize` and `Valence` start from
  sanitize=false, remove_hydrogens=false. Other molecular tasks start from
  sanitize=true, remove_hydrogens=true, with no generated coordinates.
- For those parser preparations, allow_cxsmiles=true, strict_cxsmiles=true,
  parse_name=true, skip_cleanup=false, debug_parse=false, replacements=empty.
  Every preparation outcome is recorded; a CK-only failure is never filtered.
- `RemoveHydrogens` additionally applies ordinary AddHs with explicit_only=false
  before the measured removal. This is a declared preparation dependency, not
  an extra hidden profile or a claim about all hydrogen-bearing input forms.
- AddHs fixes add_coords=false, add_residue_info=false, skip_queries=false and
  only_on_atoms=None. Spatial/residue/query-selection profiles are not included.
- RemoveHs fixes remove_with_wedged_bond, remove_mapped, remove_in_sgroups,
  show_warnings and remove_nonimplicit to true; every other boolean except the
  varied sanitize is false. No stale valence cache is compared as molecule state.
- The initial 2D profile explicitly uses canonical_orientation=false,
  clear_existing_2d=true, an empty coordinate map, flips_per_sample=0, samples=0,
  sample_seed=0, permute_degree_four=false, force_rdkit=false and
  use_ring_templates=false. Pin the reference's preferCoordGen=false too.
  These are the current CK default arguments, not an assumption that all
  RDKit defaults match. The old corpus comparison used absolute error <=1e-8;
  retain that named numerical boundary, shape and non-finite checks, never
  describe it as bit equality. Mixed 2D/3D storage is not this profile's input.

Preparation differences must be reported distinctly from operation differences.
No denominator called "mutually parseable" may conceal CK-only parser rejects.
The eventual runner must serialize resolved options in reference identities;
calling whatever Default means in a later release is not a frozen profile.

## Comparison boundaries

Topology means ordered atom IDs, element/atomic number, isotope, formal charge,
explicit H, noImplicit, radicals, aromaticity, hybridization and chiral state;
ordered bond IDs/endpoints, order, aromaticity, conjugation, direction, stereo
and stereo references. Keep property absence distinct from an explicit value.
Define the typed schema before activating any molecular task. Do not substitute
SMILES serialization or an aggregate hash for inspecting these fields.

Cache internals are not topology. `Valence` additionally needs canonical public
read-only access to its result rows: that interface must be registered before
implementation. Until then, status-only checks must not pass as full valence
parity. Ring/ring-family tasks likewise await a complete public result boundary.

CIP compares label presence and values for every atom/bond. Potential stereo
compares ordered center/type/specified/descriptor/permutation/controlling-atom
records and any cleaned topology. CK-only rank/ring-relation fields require a
source-backed reference projection before claiming those additional fields.

Every task compares success versus error. Structured error categories, atom
indices and sanitize stage must have an explicit source-backed mapping;
arbitrary cross-language exception text is not the contract. RDKit partial
mutation on failure is not permission to violate CK transaction semantics.

Value-style input preservation and available in-place projections are separate
API checks against the same chemistry expectation, not extra chemistry profiles.
Do not invent missing in-place APIs (for example sanitize_) for an old script.

## Corpora, activation and exclusions

Corpus selection is independent: the existing 152-row and 5,000-row SMILES
sets, then an explicitly selected large corpus. Fingerprint tasks consume pairs,
not SMILES; no task silently changes input kinds. Counts are input rows x the
declared profiles, with preparation failures accounted for separately.

All selected references must be prepared and globally preflighted before CK
execution. Large corpora require bounded-memory sharding; the current pilot's
in-memory vectors are not a million-row implementation. Binding checks reuse
Rust profiles later; no Python CK execution path is added here.

The future scopes above are included in the catalog but not in the initial
34-profile total. Each needs its own complete typed matrix before activation;
unavailable public APIs must not be bypassed through domain-crate dependencies.

Old `audit_core.py`, `audit_surfaces.py` and `audit_stereo.py` supply reference
branches and evidence, not task categories or inherited passing status. This
table extends some matrices explicitly (for example both non-tetrahedral modes)
and does not claim those additions were covered by the old audit.
