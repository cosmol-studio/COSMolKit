# RDKit SMARTS Full Call Path and Single-Core Audit

## Purpose

This audit fixes the source boundary, complete upstream call path, current
COSMolKit implementation inventory, and consolidation decision for the full
SMARTS port.  It is the source map for
`dev/archive/plans/rdkit_smarts_full_port_plan.md`; it is not a parity claim.

The governing rule is that existing SMARTS behavior is input to the same
source-ordered port as missing behavior.  Existing code is not a parallel
completed branch, and no historical parser, matcher, writer, recursive-query
representation, or consumer-local SMARTS decoder may survive as a second core.

## Pinned Source and Provenance

- Chemistry reference: RDKit `2026.03.1`.
- Source revision: `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- Restored local source root: `third_party/rdkit/`.
- Initial planning fallback: official RDKit GitHub files at the exact revision
  above were used only before the local tree was restored.
- Verified manifest:
  `dev/source_provenance/rdkit_smarts_2026_03_1.md` records the official archive
  hash and the hashes of all selected SMARTS source and test files.

The source prerequisite is now satisfied. All line ranges below refer to the
restored, hash-pinned local files. Future port steps must copy their verbatim
anchors from this local tree; the earlier network reads remain planning
evidence only and are not acceptable in-function source anchors.

## Full-Port Boundary

The in-scope closure is ordinary molecule SMARTS as exposed by GraphMol:

1. molecule, atom, and bond SMARTS parsing, including parser parameters,
   replacements, CXSMARTS/name handling, query-H merging, bond-direction
   stereo finalization, and cleanup;
2. the typed atom, bond, composite, range, set, property, recursive, and
   generic-group query model reached by parsing and writing;
3. ordinary `ROMol`-against-`ROMol` substructure matching, recursive matching,
   query-query matching, chirality, enhanced stereo, generic matchers,
   properties, callbacks, uniqueness, match limits, ordering, and VF2;
4. atom, bond, molecule, fragment, CXSMARTS, and fragment-CXSMARTS writing;
5. Rust and Python public SMARTS/query/match/write APIs and exact parity tests;
6. every existing descriptor, fingerprint, force-field, coordinate-template,
   MolBlock/SDF, and Python consumer routed through the one compiled query core.

Reaction SMARTS/SMIRKS is a separate ChemicalReaction source closure and is
explicitly outside this plan.  `MolBundle`, `ResonanceMolSupplier`, and
SubstructLibrary overloads are also outside this ordinary-molecule closure
because COSMolKit does not expose corresponding container models; the plan must
record structured unsupported boundaries rather than imitate them.  Generic
groups are in scope because `useGenericMatchers` is part of the ordinary
`SubstructMatchParameters` path.

## Upstream Entry-to-Core Call Path

### Python and public C++ entry points

`Code/GraphMol/Wrap/rdmolfiles.cpp` exposes `MolFromSmarts`,
`MolFromSmartsHelper`, `SmartsParserParams`, `MolToSmarts`, and
`MolToCXSmarts`.  The wrappers normalize Python strings and replacement maps,
then call the C++ V1 parse/write entry points.  `Wrap/substructmethods.h`
provides `pyFinalMatchFunctor`, `pyMatchFunctor`, conversion helpers,
`HasSubstructMatch`, `GetSubstructMatch`, `GetSubstructMatches`, and their help
dispatchers.

`Code/GraphMol/SmilesParse/SmilesParse.h` defines `SmartsParserParams` and the
public `MolFromSmarts`, `AtomFromSmarts`, and `BondFromSmarts` adapters.
Defaults are `allowCXSMILES=true`, `strictCXSMILES=true`, `parseName=true`,
`mergeHs=false`, `skipCleanup=false`, and `debugParse=false`, with a replacement
map.

### Parser dispatch and finalization

`Code/GraphMol/SmilesParse/SmilesParse.cpp` calls the following functions in
source order:

| Function | Role in the SMARTS path |
| --- | --- |
| `generic_parse_helper` | Create scanner state, buffer the input, invoke the selected grammar, and destroy scanner resources. |
| `smarts_parse_helper` | Invoke molecule SMARTS grammar mode. |
| `smarts_bond_parse` | Invoke bond-only SMARTS grammar mode. |
| `smarts_atom_parse` | Invoke atom-only SMARTS grammar mode. |
| `smarts_parse` | Dispatch molecule, atom, or bond start tokens. |
| `labelRecursivePatterns` | Assign recursive-query serial numbers across the typed query graph. |
| `toMol` | Run grammar, close rings, check chirality, set unspecified bonds, adjust chirality flags, and clean failures. |
| `toAtom` | Run the same scanner/parser core in atom mode. |
| `toBond` | Run the same scanner/parser core in bond mode. |
| `preprocessSmiles<T>` | Apply ordered replacement expansion without creating a second parser. |
| `handleCXPartAndName<T>` | Split and parse CX data and names under the parser flags. |
| `AtomFromSmarts` | Public atom-query adapter. |
| `BondFromSmarts` | Public bond-query adapter. |
| `MolFromSmarts` | Preprocess, parse, label recursion, handle CX/name, merge query Hs, set bond stereo, and clean properties. |

`MolFromSmarts` reaches `MolOps::mergeQueryHs` and `MolOps::hasQueryHs` in
`Code/GraphMol/AddHs.cpp`, including recursion into compiled recursive-query
molecules.  It reaches `MolOps::setBondStereoFromDirections` in
`Code/GraphMol/Chirality.cpp`.  These are shared molecular/query operations and
must be reused, not copied into a SMARTS-only substitute.

### Lexer and grammar

`Code/GraphMol/SmilesParse/smarts.ll` defines the `IN_ATOM`, `IN_BRANCH`, and
`IN_RECURSION` scanner states and tokens for elements, aromatic atoms,
wildcards, `a/A`, isotope and charge syntax, map numbers, `D/d/X/x/v/z/Z/h/H`,
`R/r/k`, hybridization `^0` through `^5`, tetrahedral and non-tetrahedral
chirality, every SMARTS bond operator, recursive `$(`, ranges, ring labels,
branches, dots, negation, and boolean operators.

`Code/GraphMol/SmilesParse/smarts.yy` defines these source grammar units:

- `meta_start`, `bad_atom_def`, `mol`, `atomd`, and `hydrogen_atom`;
- `atom_expr`, `point_query`, `recursive_query`, `atom_query`,
  `possible_range_query`, and `simple_atom`;
- `bond_expr`, `bond_query`, and `bondd`;
- `charge_spec`, `ring_number`, `number`, `nonzero_number`, `digit`, and
  `branch_open_token`.

The grammar directly constructs typed query nodes.  Important source semantics
include implicit high-precedence AND by adjacency, explicit high-precedence
`&`, low-precedence `;`, OR `,`, negation, range queries, hydrogen ambiguity,
map numbers as atom properties, recursive subqueries as compiled molecules,
dative-bond direction, ring bookmarks, fragment dots, and CX bond indices.

### Shared parse operations

`Code/GraphMol/SmilesParse/SmilesParseOps.cpp` supplies the shared functions
`ClearAtomChemicalProps`, `CheckRingClosureBranchStatus`, `ReportParseError`,
`CleanupAfterParseError`, `couldBeRingClosure`, `AddFragToMol`,
`GetBondOrdering`, `AdjustAtomChiralityFlags`, `GetUnspecifiedBondType`,
`SetUnspecifiedBondTypes`, `swapBondDirIfNeeded`, `checkChiralPermutation`,
`CheckChiralitySpecifications`, `CloseMolRings`, `CleanupAfterParsing`,
`getUnspecifiedQueryBond`, and `detail::printSyntaxErrorMessage`.

These operations are shared with SMILES.  A SMARTS port must extend or reuse
their single Rust equivalents and may not clone ring closure, fragment merge,
chirality correction, unspecified-bond, or error-cleanup logic into another
module.

## Typed Query Model Call Path

`Code/GraphMol/QueryOps.h/.cpp`, `QueryAtom.h/.cpp`, and `QueryBond.h/.cpp`
provide the nodes and matching functions reached by the grammar and matcher:

- null, equality, greater, greater-or-equal, less, less-or-equal, range, set,
  AND, OR, and XOR query nodes;
- atom data functions for element/type, isotope/mass, charge, aromaticity,
  degree variants, H count variants, valence, hybridization, unsaturation,
  radical/chiral state, ring membership/count/size, and hetero-neighbor counts;
- bond data functions for order, direction, ring state/count, and stereo;
- simple/range atom query factories, atom type/aliphatic/aromatic factories,
  formal-charge and isotope factories, every degree/H/valence/ring/hybridization
  factory, null query factories, bond-order/ring/single-or-aromatic factories;
- `AtomRingQuery`, `RecursiveStructureQuery`, property-presence queries,
  property-value queries, and complex A/AH/Q/QH/X/XH/M/MH queries;
- `QueryAtom::expandQuery`, `QueryAtom::Match`, `QueryAtom::QueryMatch`,
  `QueryBond::expandQuery`, `QueryBond::Match`, and `QueryBond::QueryMatch`;
- `completeMolQueries`, `finalizeQueryFromDescription`,
  `replaceAtomWithQueryAtom`, `hasBondTypeQuery`, complex-query classification,
  and query-description conversion helpers.

This typed graph is the only valid parser output and writer input.  Atom map
numbers are properties, not always-true query leaves.  Recursive queries own a
compiled typed query molecule and serial, not a SMARTS string that is reparsed
during matching.

## Matcher Call Path

`Code/GraphMol/Substruct/SubstructMatch.h` defines the ordinary match contract:
`useChirality`, `useEnhancedStereo`, `aromaticMatchesConjugated`,
`useQueryQueryMatches`, `useGenericMatchers`, `recursionPossible`, `uniquify`,
`maxMatches`, `numThreads`, atom and bond property lists, final/atom/bond
callbacks with override modes, `maxRecursiveMatches`,
`specifiedStereoQueryMatchesUnspecified`, and
`aromaticMatchesSingleOrDouble`.

`Code/GraphMol/Substruct/SubstructUtils.cpp` provides `propertyCompat`,
`atomCompat`, `chiralAtomCompat`, `bondCompat`, duplicate removal, and match
ordering helpers.  `Code/GraphMol/Substruct/SubstructMatch.cpp` then calls:

1. `hasChiralLabel` and `enhancedStereoIsOK`;
2. `MatchSubqueries` and `RecursiveMatcher` under `RecursiveLocker`;
3. `AtomLabelFunctor`, `BondLabelFunctor`, and `MolMatchFinalCheckFunctor`;
4. `boost::vf2_all` through the ordinary `SubstructMatch(ROMol, ROMol, params)`
   entry point;
5. `insertIfNeeded`/`tryToInsert` to enforce final acceptance and uniqueness;
6. `AtomCoordsMatchFunctor` where coordinate-based final checks are requested.

`Code/GraphMol/Substruct/vf2.hpp` supplies `NodeInfo`, `Pair`, `getOtherIdx`,
node sorting, every `VF2SubState` method, recursive `Match`/`MatchAll`, free
`match` overloads, `vf2`, and `vf2_all`.  The current Rust backtracker must be
source-compared function by function; passing existing feature patterns is not
evidence that VF2 ordering, pruning, uniqueness, and limit behavior are exact.

`Code/GraphMol/GenericGroups/GenericGroups.cpp` is reached when
`useGenericMatchers=true`.  Its closure includes `IsHydrogen`, `AllAtomsMatch`,
all Group/GroupH/Group-star, Alkyl, Alkenyl, Alkynyl, cyclic/acyclic,
carbo/hetero, aryl, alkoxy, D/T/H+/Pol/R matchers, ring helper functions,
`genericAtomMatcher`, `adjustQueryPropertiesWithGenericGroups`,
`convertGenericQueriesToSubstanceGroups`, and
`setGenericQueriesFromProperties`.

## Writer Call Path

`Code/GraphMol/SmilesParse/SmartsWrite.cpp` contains `_combineChildSmarts`,
`describeQuery`, `getAtomSmartsSimple`,
`getRecursiveStructureQuerySmarts`, `getBasicBondRepr`,
`getBondSmartsSimple`, `_recurseGetSmarts`, `_recurseBondSmarts`,
`FragmentSmartsConstruct`, `getNonQueryAtomSmarts`,
`getNonQueryBondSmarts`, the private `molToSmarts`,
`SmartsWrite::GetAtomSmarts`, `SmartsWrite::GetBondSmarts`, `MolToSmarts`,
`MolFragmentToSmarts`, `MolToCXSmarts`, and `MolFragmentToCXSmarts`.

The writer consumes the same typed query graph and reuses SMILES
canonicalization, traversal, and CX-extension machinery.  COSMolKit currently
has no SMARTS writer, so the port must add writer behavior to the one query
model rather than serialize the existing `SmartsMolecule` metadata tree.

## Current COSMolKit Inventory

### Duplicate parser cores

1. `crates/cosmolkit-core/src/search/query.rs::parse_smarts` contains an older
   token/parser implementation.  Repository search finds no production caller;
   its behavior is retained only by local tests.
2. `crates/cosmolkit-core/src/search/smarts_parse.rs::parse_smarts` is the
   production parser.  It first returns a separate `SmartsMolecule` metadata
   graph and `build_query_molecule` later converts that graph to `Molecule`.

The second parser is the migration base, not proof of parity.  Its
`SmartsParseParams` argument is substantially ignored, several source markers
remain partial or unsupported, and the intermediate metadata graph has become
a public Python value even though it cannot be passed directly to Python
substructure functions.

### Query and matcher divergences

- `AtomQueryPredicate::RecursiveSmarts(String)` stores source text and
  `search/substruct.rs` reparses and caches it during matching; RDKit stores a
  compiled `RecursiveStructureQuery` with a serial number.
- `AtomMapNumber` is represented as an always-true predicate; RDKit stores it
  as an atom property.
- unsupported predicate leaves can return `false` during matching; complete
  parsing and typed match preflight must instead reject unsupported state with
  a structured error.
- the Rust public match parameters expose only a small subset of RDKit's
  parameter object, and enhanced stereo/query-query/generic/callback/property
  branches remain absent or marker-open.
- existing Python match APIs accept a molecule query, while Python
  `parse_smarts` returns a different, non-matchable `SmartsMolecule` object.
- no atom, bond, molecule, fragment, or CXSMARTS writer is implemented.

### Consumer-local SMARTS behavior to remove

| Location | Historical branch | Required consolidation |
| --- | --- | --- |
| `chemistry/coordinates.rs` | `parse_rdkit_template_graph_model`, `expand_template_smarts_bonds`, and template-query construction manually decode SMARTS topology. | Read topology and maps from the one compiled query molecule. |
| `chemistry/forcefield/crystalff/torsion_preferences.rs` | Patterns are parsed twice; `map_pattern_atom_indices`, `unspecified_smarts_bond_query`, and `expand_crystalff_smarts_bonds` implement another decoder. | Compile once, retain the compiled query, and derive maps/bonds from it. |
| `chemistry/forcefield/torsion_query.rs` | Calls the conversion helper. | Call the canonical compiler directly. |
| `properties/fingerprint.rs` | Builds query molecules per pattern through the conversion helper. | Use cached canonical compiled queries. |
| `properties/descriptors.rs` | `rdkit_count_smarts_matches` reparses static patterns. | Cache canonical compiled queries without changing descriptor semantics. |
| `io/sdf.rs` and `io/sdf/postprocess.rs` | SMARTSQ creates string recursive predicates. | Parse once into the canonical compiled recursive-query type. |
| `search/substruct.rs` | Recursive matching reparses stored strings. | Match the owned compiled recursive query directly. |

Existing descriptor, fingerprint, torsion, coordinate-template, and SDF tests
are mandatory regression consumers of the core port.  They do not own SMARTS
semantics and must not preserve local fallback behavior after migration.

## Single-Core Architecture Decision

The implementation plan must converge on these invariants:

1. There is one scanner/parser engine with molecule, atom, and bond start
   modes; atom/bond public functions do not contain separate parsers.
2. Parsing returns the canonical query-bearing `Molecule`/`Atom`/`Bond` model
   directly; `SmartsMolecule` is removed rather than retained as a compatibility
   representation.
3. There is one typed atom/bond query algebra used by SMARTS, MolBlock/SDF
   query input, matching, writing, descriptors, fingerprints, force fields,
   and coordinates.
4. Recursive query nodes own compiled query molecules and serials; no match
   path reparses SMARTS text.
5. There is one ordinary-molecule matcher and one VF2 state engine; convenience
   APIs only delegate with parameters.
6. There is one SMARTS writer consuming the canonical typed query graph and
   shared SMILES/CX traversal machinery.
7. Parser replacements, generic labels, atom maps, and CX/name annotations are
   typed properties or parser metadata, never fake always-true predicates.
8. Unsupported out-of-scope source surfaces return structured unsupported
   errors and cannot silently no-op or evaluate false.
9. The older `query.rs::parse_smarts`, the intermediate `SmartsMolecule`, every
   consumer-local decoder, and recursive string reparse/cache are deleted after
   their tests and callers are migrated; no compatibility shim remains.

## Existing Work Integration Rule

Every upstream function or grammar unit in the execution plan has a source
comparison step regardless of current marker status.  If the current Rust
behavior is exact, that step retains or moves it into the canonical core,
copies the exact source anchors into the corresponding function, records the
local complexity review, and adds/retains parity evidence.  If it is partial,
the same step finishes it.  Nothing is skipped merely because it was needed by
a past fingerprint, descriptor, force-field, or coordinate task.

This produces a source-ordered main line: current code is reconciled in place,
historical duplicates are removed at explicit migration gates, and a function
cannot be considered complete while any alternate implementation of the same
SMARTS behavior remains reachable.

## Pinned Local Source Line-Range Index

All ranges are inclusive and relative to `third_party/rdkit/`. Multiple ranges denote source overloads or template/concrete definitions. This index was generated from the restored hash-pinned files; it corrects the earlier conceptual `NullQuery::Match` and `RecursiveStructureQuery::Match` labels to the actual source functions.

### `Code/GraphMol/AddHs.cpp`

- `isQueryH`: `1175-1218`
- `needsHs`: `1340-1348`
- `mergeQueryHs(RWMol)`: `1235-1332`, `1333-1338`
- `mergeQueryHs(ROMol)`: `1235-1332`, `1333-1338`
- `hasQueryHs`: `1350-1399`

### `Code/GraphMol/Chirality.cpp`

- `MolOps::setBondStereoFromDirections`: `3707-3746`

### `Code/GraphMol/GenericGroups/GenericGroups.cpp`

- `IsHydrogen`: `25-33`
- `AllAtomsMatch`: `35-76`
- `GroupAtomMatcher`: `78-91`
- `GroupHAtomMatcher`: `93-100`
- `GroupStarAtomMatcher`: `102-117`
- `GroupStarHAtomMatcher`: `119-126`
- `AlkylAtomMatcher`: `128-146`
- `AlkylHAtomMatcher`: `148-163`
- `AlkenylAtomMatcher`: `314-317`
- `AlkenylHAtomMatcher`: `319-326`
- `AlkynylAtomMatcher`: `328-331`
- `AlkynylHAtomMatcher`: `333-340`
- `AcyclicAtomMatcher`: `165-179`
- `AcyclicHAtomMatcher`: `181-190`
- `CarboacyclicAtomMatcher`: `192-205`
- `CarboacyclicHAtomMatcher`: `207-217`
- `HeteroacyclicAtomMatcher`: `219-233`
- `HeteroacyclicHAtomMatcher`: `235-242`
- `AlkoxyacyclicAtomMatcher`: `244-275`
- `AlkoxyacyclicHAtomMatcher`: `277-284`
- `checkAtomRing`: `343-358`
- `checkBondRing`: `359-371`
- `FusedRingMatch`: `373-447`
- `CarbocycloalkylAtomMatcher`: `450-459`
- `CarbocycloalkylHAtomMatcher`: `461-468`
- `CarbocycloalkenylAtomMatcher`: `470-483`
- `CarbocycloalkenylHAtomMatcher`: `485-492`
- `CarboarylAtomMatcher`: `494-503`
- `CarboarylHAtomMatcher`: `505-512`
- `CarbocyclicAtomMatcher`: `514-520`
- `CarbocyclicHAtomMatcher`: `522-529`
- `NoCarbonRingAtomMatcher`: `531-537`
- `NoCarbonRingHAtomMatcher`: `539-546`
- `HeterocyclicAtomMatcher`: `548-559`
- `HeterocyclicHAtomMatcher`: `561-568`
- `HeteroarylAtomMatcher`: `570-583`
- `HeteroarylHAtomMatcher`: `585-592`
- `CyclicAtomMatcher`: `594-603`
- `CyclicHAtomMatcher`: `605-612`
- `DAtomMatcher`: `614-621`
- `TAtomMatcher`: `623-630`
- `HplusAtomMatcher`: `632-639`
- `PolAtomMatcher`: `641-650`
- `RAtomMatcher`: `652-655`
- `genericAtomMatcher`: `659-682`
- `adjustQueryPropertiesWithGenericGroups`: `684-695`
- `convertGenericQueriesToSubstanceGroups`: `697-709`
- `setGenericQueriesFromProperties`: `711-758`

### `Code/GraphMol/QueryAtom.cpp`

- `QueryAtom::copy`: `22-25`
- `QueryAtom::expandQuery`: `27-66`
- `localMatch`: `69-75`
- `queriesMatch`: `76-178`
- `QueryAtom::Match`: `181-185`
- `QueryAtom::QueryMatch`: `186-194`

### `Code/GraphMol/QueryBond.cpp`

- `QueryBond::QueryBond`: `15-21`, `23-26`
- `QueryBond::copy`: `42-45`
- `QueryBond::setBondType`: `47-54`
- `QueryBond::setBondDir`: `56-70`
- `QueryBond::expandQuery`: `72-110`
- `localMatch`: `113-119`
- `queriesMatch`: `121-203`
- `QueryBond::Match`: `206-210`
- `QueryBond::QueryMatch`: `211-219`
- `QueryBond::getValenceContrib`: `221-226`

### `Code/GraphMol/QueryOps.cpp`

- `queryIsAtomBridgehead`: `22-87`
- `queryAtomBondProduct`: `198-208`
- `queryAtomAllBondProduct`: `209-223`
- `isMetal`: `1180-1183`
- `makeAtomNumQuery`: `245-248`
- `makeAtomTypeQuery`: `250-253`
- `makeAtomImplicitValenceQuery`: `225-230`
- `makeAtomExplicitValenceQuery`: `231-236`
- `makeAtomTotalValenceQuery`: `238-243`
- `makeAtomExplicitDegreeQuery`: `254-259`
- `makeAtomTotalDegreeQuery`: `261-266`
- `makeAtomHeavyAtomDegreeQuery`: `268-273`
- `makeAtomHCountQuery`: `275-279`
- `makeAtomHasImplicitHQuery`: `286-291`
- `makeAtomImplicitHCountQuery`: `280-285`
- `makeAtomAromaticQuery`: `293-297`
- `makeAtomAliphaticQuery`: `299-303`
- `makeAtomMassQuery`: `312-317`
- `makeAtomIsotopeQuery`: `319-323`
- `makeAtomFormalChargeQuery`: `325-330`
- `makeAtomNegativeFormalChargeQuery`: `332-337`
- `makeAtomHybridizationQuery`: `339-344`
- `makeAtomNumRadicalElectronsQuery`: `346-351`
- `makeAtomHasChiralTagQuery`: `353-358`
- `makeAtomMissingChiralTagQuery`: `360-365`
- `makeAtomUnsaturatedQuery`: `305-310`
- `makeAtomInRingQuery`: `367-371`
- `makeAtomInNRingsQuery`: `504-509`
- `makeAtomInRingOfSizeQuery`: `97-104`, `106-117`
- `makeAtomMinRingSizeQuery`: `183-189`
- `makeAtomRingBondCountQuery`: `90-95`
- `makeAtomHasRingBondQuery`: `511-516`
- `makeAtomNumHeteroatomNbrsQuery`: `518-523`
- `makeAtomHasHeteroatomNbrsQuery`: `525-530`
- `makeAtomNumAliphaticHeteroatomNbrsQuery`: `531-536`
- `makeAtomHasAliphaticHeteroatomNbrsQuery`: `538-543`
- `makeAtomNonHydrogenDegreeQuery`: `545-550`
- `makeAtomIsBridgeheadQuery`: `373-378`
- `makeQAtomQuery`: `380-390`
- `makeQHAtomQuery`: `391-396`
- `makeAAtomQuery`: `397-402`
- `makeAHAtomQuery`: `403-407`
- `makeXAtomQuery`: `409-425`
- `makeXHAtomQuery`: `426-432`
- `makeMAtomQuery`: `434-446`
- `makeMHAtomQuery`: `447-502`
- `makeBondOrderEqualsQuery`: `552-559`
- `makeSingleOrAromaticBondQuery`: `561-568`
- `makeDoubleOrAromaticBondQuery`: `570-577`
- `makeSingleOrDoubleBondQuery`: `579-586`
- `makeSingleOrDoubleOrAromaticBondQuery`: `589-596`
- `makeBondDirEqualsQuery`: `656-662`
- `makeBondHasStereoQuery`: `664-670`
- `makeBondIsInRingQuery`: `672-678`
- `makeBondInNRingsQuery`: `680-686`
- `makeBondInRingOfSizeQuery`: `119-181`
- `makeBondMinRingSizeQuery`: `190-196`
- `makeBondNullQuery`: `688-694`
- `makeAtomNullQuery`: `696-702`
- `convertComplexNameToQuery`: `704-726`
- `hasBondTypeQuery`: `605-623`
- `hasComplexBondTypeQuery`: `650-653`
- `isComplexQuery(Bond)`: `728-766`, `892-919`
- `isComplexQuery(Atom)`: `728-766`, `892-919`
- `isAtomListQuery`: `819-846`
- `getAtomListQueryVals`: `848-890`
- `completeQueryAndChildren`: `963-978`
- `completeMolQueries`: `980-987`
- `replaceAtomWithQueryAtom`: `989-1004`
- `finalizeAtomRingSizeQuery`: `1012-1046`
- `finalizeQueryFromDescription(Atom)`: `1048-1138`, `1140-1178`
- `finalizeQueryFromDescription(Bond)`: `1048-1138`, `1140-1178`

### `Code/GraphMol/QueryOps.h`

- `queryAtomAromatic`: `77-79`
- `queryAtomAliphatic`: `80-82`
- `queryAtomExplicitDegree`: `83-85`
- `queryAtomTotalDegree`: `86-88`
- `queryAtomNonHydrogenDegree`: `90-101`
- `queryAtomHeavyAtomDegree`: `103-114`
- `queryAtomHCount`: `115-117`
- `queryAtomImplicitHCount`: `118-120`
- `queryAtomHasImplicitH`: `121-123`
- `queryAtomImplicitValence`: `124-126`
- `queryAtomExplicitValence`: `127-129`
- `queryAtomTotalValence`: `130-132`
- `queryAtomUnsaturated`: `133-135`
- `queryAtomNum`: `136-136`
- `makeAtomType`: `137-139`
- `parseAtomType`: `140-148`
- `getAtomTypeIsAromatic`: `149-149`
- `getAtomTypeAtomicNum`: `150-155`
- `queryAtomType`: `157-159`
- `queryAtomMass`: `161-164`
- `queryAtomIsotope`: `165-167`
- `queryAtomFormalCharge`: `168-170`
- `queryAtomNegativeFormalCharge`: `171-173`
- `queryAtomHybridization`: `174-176`
- `queryAtomNumRadicalElectrons`: `177-179`
- `queryAtomHasChiralTag`: `180-182`
- `queryAtomMissingChiralTag`: `183-186`
- `queryAtomHasHeteroatomNbrs`: `188-199`
- `queryAtomNumHeteroatomNbrs`: `201-213`
- `queryAtomHasAliphaticHeteroatomNbrs`: `215-227`
- `queryAtomNumAliphaticHeteroatomNbrs`: `229-242`
- `queryBondOrder`: `250-252`
- `queryBondIsSingleOrAromatic`: `253-256`
- `queryBondIsDoubleOrAromatic`: `257-260`
- `queryBondIsSingleOrDouble`: `261-264`
- `queryBondIsSingleOrDoubleOrAromatic`: `265-269`
- `queryBondDir`: `270-272`
- `queryIsBondInNRings`: `273-275`
- `queryBondHasStereo`: `276-278`
- `queryIsAtomInNRings`: `283-285`
- `queryIsAtomInRing`: `286-288`
- `queryAtomHasRingBond`: `289-300`
- `queryIsBondInRing`: `303-305`
- `queryAtomMinRingSize`: `306-308`
- `queryBondMinRingSize`: `309-311`
- `queryAtomRingBondCount`: `313-326`
- `queryAtomIsInRingOfSize`: `328-334`, `340-359`, `361-367`
- `queryBondIsInRingOfSize`: `369-376`
- `queryAtomRingMembership`: `741-744`
- `nullDataFun`: `852-854`
- `nullQueryFun`: `856-858`
- `isAtomDummy`: `1168-1172`
- `makeAtomSimpleQuery`: `379-386`
- `makeAtomRangeQuery`: `388-397`
- `makeAtomNumQuery`: `401-403`
- `makeAtomTypeQuery`: `409-412`
- `makeAtomImplicitValenceQuery`: `419-421`
- `makeAtomExplicitValenceQuery`: `427-429`
- `makeAtomTotalValenceQuery`: `435-437`
- `makeAtomExplicitDegreeQuery`: `443-445`
- `makeAtomTotalDegreeQuery`: `451-453`
- `makeAtomHeavyAtomDegreeQuery`: `459-461`
- `makeAtomHCountQuery`: `467-469`
- `makeAtomHasImplicitHQuery`: `475-477`
- `makeAtomImplicitHCountQuery`: `483-485`
- `makeAtomAromaticQuery`: `491-493`
- `makeAtomAliphaticQuery`: `499-501`
- `makeAtomMassQuery`: `507-510`
- `makeAtomIsotopeQuery`: `516-518`
- `makeAtomFormalChargeQuery`: `524-526`
- `makeAtomNegativeFormalChargeQuery`: `533-535`
- `makeAtomHybridizationQuery`: `542-544`
- `makeAtomNumRadicalElectronsQuery`: `550-552`
- `makeAtomHasChiralTagQuery`: `560-562`
- `makeAtomMissingChiralTagQuery`: `569-571`
- `makeAtomUnsaturatedQuery`: `577-579`
- `makeAtomInRingQuery`: `585-587`
- `makeAtomInNRingsQuery`: `593-595`
- `makeAtomMinRingSizeQuery`: `610-612`
- `makeAtomRingBondCountQuery`: `618-620`
- `makeAtomHasRingBondQuery`: `650-652`
- `makeAtomNumHeteroatomNbrsQuery`: `658-660`
- `makeAtomHasHeteroatomNbrsQuery`: `667-669`
- `makeAtomNumAliphaticHeteroatomNbrsQuery`: `675-678`
- `makeAtomHasAliphaticHeteroatomNbrsQuery`: `685-687`
- `makeAtomNonHydrogenDegreeQuery`: `694-696`
- `makeAtomIsBridgeheadQuery`: `703-705`
- `makeHasPropQuery`: `902-905`
- `makePropQuery`: `1142-1145`, `1148-1153`
- `AtomRingQuery::Match`: `767-781`
- `RecursiveStructureQuery::RecursiveStructureQuery()`: `796-799`
- `RecursiveStructureQuery::RecursiveStructureQuery(ROMol,unsigned int)`: `805-811`
- `RecursiveStructureQuery::getAtIdx`: `813-816`
- `RecursiveStructureQuery::setQueryMol`: `823-823`
- `RecursiveStructureQuery::getQueryMol`: `825-825`
- `RecursiveStructureQuery::copy`: `828-840`
- `RecursiveStructureQuery::getSerialNumber`: `841-841`

### `Code/GraphMol/SmilesParse/SmartsWrite.cpp`

- `_combineChildSmarts`: `45-91`
- `describeQuery`: `94-100`
- `getAtomSmartsSimple`: `103-341`
- `getRecursiveStructureQuerySmarts`: `343-356`
- `getBasicBondRepr`: `358-411`
- `getBondSmartsSimple`: `413-469`
- `_recurseGetSmarts`: `471-567`
- `_recurseBondSmarts`: `569-664`
- `FragmentSmartsConstruct`: `666-732`
- `getNonQueryAtomSmarts`: `736-802`
- `getNonQueryBondSmarts`: `806-830`
- `molToSmarts`: `832-896`
- `SmartsWrite::GetAtomSmarts`: `901-963`
- `SmartsWrite::GetBondSmarts`: `965-1007`
- `MolToSmarts`: `1011-1021`
- `MolFragmentToSmarts`: `1023-1058`
- `MolToCXSmarts`: `1060-1071`
- `MolFragmentToCXSmarts`: `1073-1085`

### `Code/GraphMol/SmilesParse/SmilesParse.cpp`

- `generic_parse_helper`: `70-114`
- `smarts_parse_helper`: `116-128`
- `smarts_bond_parse`: `129-134`
- `smarts_atom_parse`: `136-141`
- `smarts_parse`: `143-148`
- `labelRecursivePatterns`: `186-239`
- `toMol`: `242-289`
- `toAtom`: `291-310`
- `toBond`: `312-331`
- `preprocessSmiles<SmartsParserParams>`: `336-372`
- `handleCXPartAndName<SmartsParserParams>`: `389-417`
- `AtomFromSmarts`: `536-540`
- `BondFromSmarts`: `542-546`
- `MolFromSmarts`: `548-576`

### `Code/GraphMol/SmilesParse/SmilesParseOps.cpp`

- `ClearAtomChemicalProps`: `25-30`
- `CheckRingClosureBranchStatus`: `32-64`
- `ReportParseError`: `66-73`
- `CleanupAfterParseError`: `75-84`
- `couldBeRingClosure`: `87-87`
- `AddFragToMol`: `93-199`
- `GetBondOrdering`: `214-273`
- `AdjustAtomChiralityFlags`: `275-347`
- `GetUnspecifiedBondType`: `349-361`
- `SetUnspecifiedBondTypes`: `362-375`
- `swapBondDirIfNeeded`: `378-396`
- `checkChiralPermutation`: `406-413`
- `CheckChiralitySpecifications`: `415-446`
- `CloseMolRings`: `448-633`
- `CleanupAfterParsing`: `635-678`
- `getUnspecifiedQueryBond`: `680-693`
- `detail::printSyntaxErrorMessage`: `697-729`

### `Code/GraphMol/Substruct/SubstructMatch.cpp`

- `hasChiralLabel`: `38-42`
- `enhancedStereoIsOK`: `44-112`
- `insertIfNeeded`: `129-154`
- `tryToInsert`: `156-167`
- `MolMatchFinalCheckFunctor constructor`: `179-192`
- `MolMatchFinalCheckFunctor::operator()`: `194-382`, `392-408`, `420-436`, `720-733`
- `AtomLabelFunctor::operator()`: `194-382`, `392-408`, `420-436`, `720-733`
- `BondLabelFunctor::operator()`: `194-382`, `392-408`, `420-436`, `720-733`
- `RecursiveLocker constructor/destructor`: `461-465`, `467-474`
- `SubstructMatch(ROMol,ROMol,params)`: `481-525`, `527-535`, `537-545`, `547-557`, `564-606`
- `RecursiveMatcher`: `609-663`
- `MatchSubqueries(recursive execution)`: `665-716`
- `AtomCoordsMatchFunctor::operator()`: `194-382`, `392-408`, `420-436`, `720-733`

### `Code/GraphMol/Substruct/SubstructUtils.cpp`

- `propertyCompat`: `104-124`
- `atomCompat`: `126-156`
- `chiralAtomCompat`: `158-175`
- `bondCompat`: `177-241`
- `removeDuplicates`: `243-272`
- `getMostSubstitutedCoreMatch`: `51-58`, `274-280`
- `sortMatchesByDegreeOfCoreSubstitution`: `59-70`, `282-288`
- `isAtomTerminalRGroupOrQueryHydrogen`: `290-295`
- `updateSubstructMatchParamsFromJSON`: `300-320`
- `substructMatchParamsToJSON`: `322-340`

### `Code/GraphMol/Substruct/vf2.hpp`

- `nodeInfoComp1`: `45-59`
- `nodeInfoComp2`: `65-85`
- `getOtherIdx`: `88-95`
- `SortNodesByFrequency`: `110-145`
- `VF2SubState::VF2SubState`: `175-214`, `216-239`, `241-252`
- `VF2SubState copy constructor`: `175-214`, `216-239`, `241-252`
- `VF2SubState::IsGoal`: `254-254`
- `VF2SubState::MatchChecks`: `255-257`
- `VF2SubState::IsDead`: `258-258`
- `VF2SubState::CoreLen`: `259-259`
- `VF2SubState::NextPair`: `263-357`
- `VF2SubState::IsFeasiblePair`: `358-435`
- `VF2SubState::AddPair`: `436-478`
- `VF2SubState::GetCoreSet`: `479-488`
- `VF2SubState::Clone`: `489-489`
- `VF2SubState::BackTrack`: `490-525`
- `VF2SubState::Match`: `526-549`
- `VF2SubState::MatchAll`: `552-582`
- `match(single)`: `595-602`, `614-618`
- `match(all)`: `595-602`, `614-618`
- `vf2`: `632-653`
- `vf2_all`: `663-677`

### `Code/GraphMol/Wrap/rdmolfiles.cpp`

- `MolFromSmarts`: `83-101`
- `MolFromSmartsHelper`: `538-547`

### `Code/GraphMol/Wrap/substructmethods.h`

- `pyFinalMatchFunctor`: `26-26`
- `pyMatchFunctor`: `43-43`
- `convertMatches`: `55-61`
- `convertMatchesToTupleOfPairs`: `63-73`
- `HasSubstructMatch`: `76-83`
- `GetSubstructMatch`: `86-95`
- `GetSubstructMatches`: `98-114`
- `pySubstructHelper`: `117-124`
- `helpHasSubstructMatch`: `126-133`
- `helpGetSubstructMatch`: `136-147`
- `helpGetSubstructMatches`: `150-159`

### `Code/Query/AndQuery.h`

- `AndQuery::Match`: `27-41`

### `Code/Query/EqualityQuery.h`

- `EqualityQuery::Match`: `51-59`

### `Code/Query/GreaterEqualQuery.h`

- `GreaterEqualQuery::Match`: `38-46`

### `Code/Query/GreaterQuery.h`

- `GreaterQuery::Match`: `38-46`

### `Code/Query/LessEqualQuery.h`

- `LessEqualQuery::Match`: `38-46`

### `Code/Query/LessQuery.h`

- `LessQuery::Match`: `38-46`

### `Code/Query/NullQueryAlgebra.h`

- `mergeBothNullQ`: `19-41`
- `mergeNullQFirst`: `44-62`
- `mergeNullQueries`: `67-83`

### `Code/Query/OrQuery.h`

- `OrQuery::Match`: `26-40`

### `Code/Query/Query.h`

- `Query::addChild`: `108-108`
- `Query::setNegation`: `62-62`

### `Code/Query/RangeQuery.h`

- `RangeQuery::Match`: `61-84`

### `Code/Query/SetQuery.h`

- `SetQuery::Match`: `42-46`

### `Code/Query/XOrQuery.h`

- `XOrQuery::Match`: `27-45`

### `Code/GraphMol/SmilesParse/smarts.yy`

- `meta_start`: `189-234`
- `bad_atom_def`: `235-246`
- `mol`: `247-395`
- `atomd`: `396-413`
- `hydrogen_atom`: `425-487`
- `atom_expr`: `488-522`
- `point_query`: `523-533`
- `recursive_query`: `534-574`
- `atom_query`: `575-763`
- `possible_range_query`: `764-797`
- `simple_atom`: `798-819`
- `bond_expr`: `820-837`
- `bond_query`: `838-846`
- `bondd`: `847-877`
- `charge_spec`: `878-886`
- `ring_number`: `887-897`
- `number`: `898-902`
- `nonzero_number`: `903-912`
- `digit`: `913-917`
- `branch_open_token`: `918-920`

### `Code/GraphMol/SmilesParse/smarts.ll`

- `INITIAL and punctuation tokens`: `98-116, 288-468`
- `IN_ATOM state`: `94-94, 117-286, 324-330, 395-415`
- `IN_BRANCH state`: `95-95, 397-398`
- `IN_RECURSION state`: `96-96, 395-399`
- `element and aromatic-element tokens`: `117-222, 291-330`
- `atom primitive tokens`: `224-350`
- `bond operator tokens`: `354-390`
- `chirality and hybridization tokens`: `108-114, 432-464`
- `range and ring-number tokens`: `401-430`
- `scanner error and end-of-input rules`: `465-468`

### `Code/GraphMol/Wrap/rdmolfiles.cpp` registration blocks

- `BOOST_PYTHON_MODULE(rdmolfiles)`: `884-2914`
- SMARTS parser parameters and `MolFromSmarts` registrations: `1641-1767`
- `MolToSmarts` and `MolToCXSmarts` registrations: `2083-2140`

Every function-labelled port item resolves to at least one local pinned definition range.

## Validation Evidence Required for Closure

- pinned upstream parser tests from `smatest.cpp` and
  `smarts_catch_tests.cpp` covering success, rejection, typed query trees,
  properties, maps, recursion, ranges, chirality, CX/name, and merge-H behavior;
- pinned upstream matcher tests from `testSubstructMatch.cpp`, `catch_tests.cpp`,
  and `Wrap/testSubstructureMatch.py` covering exact match mappings, order,
  counts, limits, parameters, query-query behavior, stereo, recursion, generic
  groups, properties, and callbacks;
- writer tests covering atom/bond/molecule/fragment output and parse-write-parse
  semantic equivalence, including CXSMARTS;
- regression tests for every migrated internal consumer and for MolBlock/SDF
  query interoperability;
- strict core release checks, strict workspace release tests, Python binding
  rebuild/tests/type checks/stubs/docs, and a final no-duplicate-symbol audit.

No aggregate pass upgrades a marker whose corresponding function lacks exact
source anchors and local performance review.  No passing consumer test proves
general SMARTS parity outside the behavior it observes.
