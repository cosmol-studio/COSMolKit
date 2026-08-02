# RDKit MMFF Symbol Inventory

Step: `dev/archive/plans/rdkit_forcefield_full_port_plan.md` Step 675.

Scope audited:
- `third_party/rdkit/Code/ForceField/MMFF/Params.h`
- `third_party/rdkit/Code/ForceField/MMFF/Params.cpp`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/MMFF.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/Builder.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/Builder.cpp`

This is an inventory only. It does not declare any port complete. Every later
port step must still copy the corresponding RDKit source lines into the Rust
implementation and keep two-axis source-reproduction markers.

## Parameter Data Types

From `ForceField/MMFF/Params.h`:

- `DEG2RAD`, `RAD2DEG`, `MDYNE_A_TO_KCAL_MOL`.
- `isDoubleZero(double)`.
- `clipToOne(double &)`.
- `MMFFDef`: `eqLevel[4]`.
- `MMFFProp`: `atno`, `crd`, `val`, `pilp`, `mltb`, `arom`, `linh`, `sbmb`.
- `MMFFPBCI`: `pbci`, `fcadj`.
- `MMFFChg`: `bci`.
- `MMFFBond`: `kb`, `r0`.
- `MMFFHerschbachLaurie`: `a_ij`, `d_ij`, `dp_ij`.
- `MMFFCovRadPauEle`: `r0`, `chi`.
- `MMFFAngle`: `ka`, `theta0`.
- `MMFFStbn`: `kbaIJK`, `kbaKJI`.
- `MMFFOop`: `koop`.
- `MMFFTor`: `V1`, `V2`, `V3`.
- `MMFFVdW`: `alpha_i`, `N_i`, `A_i`, `G_i`, `R_star`, `DA`.
- `MMFFVdWRijstarEps`: `R_ij_starUnscaled`, `epsilonUnscaled`,
  `R_ij_star`, `epsilon`.

## Parameter Collections And Queries

From `ForceField/MMFF/Params.h` and `Params.cpp`:

- `MMFFAromCollection`
  - Constructor takes optional `std::vector<std::uint8_t> *`; empty argument
    uses `defaultMMFFArom`.
  - Query: `isMMFFAromatic(atomType)`.
  - Storage: `d_params`.
- `MMFFDefCollection`
  - Constructor parses `defaultMMFFDef` or caller string.
  - Query: `operator()(atomType)`.
  - Branches: skip lines beginning `*`; skip duplicate atom types via
    `oldAtomType`; map lookup under `RDKIT_MMFF_PARAMS_USE_STD_MAP`; vector
    lookup otherwise; missing or zero atom type returns null in vector mode.
- `MMFFPropCollection`
  - Constructor parses `defaultMMFFProp` or caller string.
  - Query: `operator()(atomType)`.
  - Branches: comment-line skip; map lookup or vector `equal_range` over
    `d_iAtomType`; missing atom type returns null.
- `MMFFPBCICollection`
  - Constructor parses `defaultMMFFPBCI` or caller string.
  - Query: `operator()(atomType)`.
  - Branches: comment-line skip; map mode skips first field before atom type;
    vector mode skips first two fields and relies on table order; null for
    out-of-range atom type.
- `MMFFChgCollection`
  - Constructor parses `defaultMMFFChg` or caller string.
  - Query: `getMMFFChgParams(bondType, iAtomType, jAtomType)`.
  - Branches: canonicalize atom-type order and return sign `-1` or `1`;
    map or vector `equal_range` lookup over atom type and bond type; missing
    parameter returns null with the canonical sign.
- `MMFFBondCollection`
  - Constructor parses `defaultMMFFBond` or caller string.
  - Query: `operator()(bondType, atomType, nbrAtomType)`.
  - Branches: canonicalize atom-type order; map or vector 3-key lookup;
    missing parameter returns null.
- `MMFFBndkCollection`
  - Constructor parses `defaultMMFFBndk` or caller string.
  - Query: `operator()(atomicNum, nbrAtomicNum)`.
  - Branches: canonicalize atomic-number order; map or vector 2-key lookup;
    missing parameter returns null.
- `MMFFHerschbachLaurieCollection`
  - Constructor parses `defaultMMFFHerschbachLaurie` or caller string.
  - Query: `operator()(iRow, jRow)`.
  - Branches: canonicalize row order; map or vector 2-key lookup; missing
    parameter returns null.
- `MMFFCovRadPauEleCollection`
  - Constructor parses `defaultMMFFCovRadPauEle` or caller string.
  - Query: `operator()(atomicNum)`.
  - Branches: map lookup or vector `equal_range`; missing parameter returns
    null.
- `MMFFAngleCollection`
  - Constructor parses caller string or concatenates `defaultMMFFAngleData[]`
    until `"EOS"`.
  - Query: `operator()(mmffDef, angleType, iAtomType, jAtomType, kAtomType)`.
  - Branches: four-stage equivalence lookup using MMFF definition levels;
    terminal atoms canonicalized by swapping `i/k`; map or vector lookup;
    missing parameter returns null.
- `MMFFStbnCollection`
  - Constructor parses `defaultMMFFStbn` or caller string.
  - Query: `getMMFFStbnParams(stretchBendType, bondType1, bondType2,
    iAtomType, jAtomType, kAtomType)`.
  - Branches: canonicalize `i/k`; return swap flag; equal atom types use
    `bondType1 < bondType2` to decide swap; map or vector lookup; missing
    parameter returns null.
- `MMFFDfsbCollection`
  - Constructor parses `defaultMMFFDfsb` or caller string.
  - Query: `getMMFFDfsbParams(periodicTableRow1, periodicTableRow2,
    periodicTableRow3)`.
  - Branches: canonicalize row 1/3 and return swap flag; map-only lookup;
    missing parameter returns null.
- `MMFFOopCollection`
  - Constructor takes `isMMFFs` and parses `defaultMMFFOop` or
    `defaultMMFFsOop` unless caller data is supplied.
  - Query: `operator()(mmffDef, iAtomType, jAtomType, kAtomType, lAtomType)`.
  - Branches: four-stage equivalence lookup; sort non-central atom types;
    map or vector lookup; MMFF94/MMFF94s data selection; missing parameter
    returns null.
- `MMFFTorCollection`
  - Constructor takes `isMMFFs` and parses `defaultMMFFTor` or
    `defaultMMFFsTor` unless caller data is supplied.
  - Query: `getMMFFTorParams(mmffDef, torType, iAtomType, jAtomType,
    kAtomType, lAtomType)`.
  - Branches: five-stage torsion equivalence lookup; half-wild-card stages;
    central-pair and terminal-pair canonicalization; special retry using
    secondary torsion type after type 5 fallback; map or vector lookup;
    break behavior when `maxIter == 4`; missing parameter returns null with
    selected torsion type.
- `MMFFVdWCollection`
  - Constructor parses `defaultMMFFVdW` or caller string.
  - Query: `operator()(atomType)`.
  - Branches: first non-comment row stores global `power`, `B`, `Beta`,
    `DARAD`, `DAEPS`; remaining rows store atom data and computed `R_star`;
    map lookup or vector `equal_range`; missing parameter returns null.

## Default Parameter Tables

From `ForceField/MMFF/Params.cpp`:

- `defaultMMFFArom`.
- `defaultMMFFDef`.
- `defaultMMFFProp`.
- `defaultMMFFPBCI`.
- `defaultMMFFChg`.
- `defaultMMFFBond`.
- `defaultMMFFBndk`.
- `defaultMMFFHerschbachLaurie`.
- `defaultMMFFCovRadPauEle`.
- `defaultMMFFAngleData[]` with `"EOS"` sentinel.
- `defaultMMFFStbn`.
- `defaultMMFFDfsb`.
- `defaultMMFFOop`.
- `defaultMMFFsOop`.
- `defaultMMFFTor`.
- `defaultMMFFsTor`.
- `defaultMMFFVdW`.

## Default Parameter Singletons

From `GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h` and `.cpp` namespace
`RDKit::MMFF::DefaultParameters`:

- `getMMFFProp()`.
- `getMMFFArom()`.
- `getMMFFBndk()` in `.cpp` and required by empirical bond stretching.
- `getMMFFBond()`.
- `getMMFFChg()`.
- `getMMFFDef()`.
- `getMMFFHerschbachLaurie()`.
- `getMMFFPBCI()`.
- `getMMFFAngle()`.
- `getMMFFStbn()`.
- `getMMFFDfsb()`.
- `getMMFFCovRadPauEle()` in `.cpp` and required by empirical bond stretching.
- `getMMFFTor(isMMFFs)` with separate static `MMFF94` and `MMFF94s` instances.
- `getMMFFOop(isMMFFs)` with separate static `MMFF94` and `MMFF94s` instances.
- `getMMFFVdW()`.

## Atom-Typer Data Types And Properties

From `GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h` and `.cpp`:

- `RingMembership`
  - `getIsInAromaticRing()`.
  - `setIsInAromaticRing(bool)`.
  - `getRingIdxSet() const`.
  - `getRingIdxSet()`.
  - Data: `d_isInAromaticRing`, `d_ringIdxSet`.
- `RingMembershipSize`
  - Type aliases: `RingMembershipMap`, `RingSizeMembershipMap`.
  - Constant: `IS_AROMATIC_BIT`.
  - Constructor indexes molecule ring membership by ring size and atom index.
  - Methods: `isAtomInAromaticRingOfSize`, `areAtomsInSameAromaticRing`,
    `areAtomsInSameRingOfSize`.
  - Branches: maximum ring count precondition; aromatic-ring flagging; insert
    per-size and per-atom maps; set intersection for common ring checks;
    varargs ring-size query.
- `MMFFAtomProperties`
  - Data: `mmffAtomType`, `mmffFormalCharge`, `mmffPartialCharge`.
  - Pointer alias: `MMFFAtomPropertiesPtr`.
- Dielectric model enum values: `CONSTANT`, `DISTANCE`.
- Verbosity enum values: `MMFF_VERBOSITY_NONE`, `MMFF_VERBOSITY_LOW`,
  `MMFF_VERBOSITY_HIGH`.
- `MMFFMolProperties`
  - Constructor: sanitizes/kekulizes as needed, sets MMFF aromaticity, assigns
    heavy atom types before hydrogen types, computes charges if valid, and
    emits high-verbosity atom table.
  - Atom-property getters: `getMMFFAtomType`, `getMMFFFormalCharge`,
    `getMMFFPartialCharge`.
  - Term toggles: `set/getMMFFBondTerm`, `set/getMMFFAngleTerm`,
    `set/getMMFFStretchBendTerm`, `set/getMMFFOopTerm`,
    `set/getMMFFTorsionTerm`, `set/getMMFFVdWTerm`, `set/getMMFFEleTerm`.
  - Variant: `setMMFFVariant`, `getMMFFVariant`.
  - Dielectric: `setMMFFDielectricConstant`, `getMMFFDielectricConstant`,
    `setMMFFDielectricModel`, `getMMFFDielectricModel`.
  - Verbosity/output: `setMMFFVerbosity`, `getMMFFVerbosity`,
    `setMMFFOStream`, `getMMFFOStream`.
  - Validity: `isValid`.
  - Parameter methods: `getMMFFBondType`, `getMMFFAngleType`,
    `getMMFFTorsionType`, `computeMMFFCharges`,
    `getMMFFTorsionEmpiricalRuleParams`,
    `getMMFFBondStretchEmpiricalRuleParams`,
    `getMMFFBondStretchParams`, `getMMFFAngleBendParams`,
    `getMMFFStretchBendParams`, `getMMFFTorsionParams`,
    `getMMFFOopBendParams`, `getMMFFVdWParams`.
  - Private helpers: `setMMFFHeavyAtomType`, `setMMFFHydrogenType`,
    `setMMFFFormalCharge`, `setMMFFPartialCharge`.

## Atom-Typer Free Functions

From `GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h` and `.cpp`:

- `getPeriodicTableRow(atomicNum)`.
- `getPeriodicTableRowHL(atomicNum)` from `.cpp`, used by empirical bond
  stretching.
- `isAromaticAtomType(atomType)` from `.cpp`.
- `isRingAromatic(mol, ringIndxVect)`.
- `isAtomInAromaticRingOfSize(atom, ringSize)`.
- `isAtomNOxide(atom)`.
- `isAngleInRingOfSize3or4(mol, idx1, idx2, idx3)`.
- `isTorsionInRingOfSize4or5(mol, idx1, idx2, idx3, idx4)`.
- `areAtomsInSameRingOfSize(mol, ringSize, numAtoms, ...)`.
- `areAtomsInSameAromaticRing(mol, idx1, idx2)`.
- `sanitizeMMFFMol(RWMol &)`.
- `getMMFFStretchBendType(angleType, bondType1, bondType2)`.
- `getMMFFAngleBendEmpiricalRuleParams(...)`.

## Atom-Typer Branch Families To Port

- Heavy atom typing:
  - Aromatic 5-ring and 6-ring branches.
  - Ring alpha/beta heteroatom analysis.
  - Aromatic carbon, nitrogen, oxygen, and sulfur assignments.
  - Aliphatic carbon by degree, ring size, double/triple bonds, guanidinium,
    carboxylate, carbonate, carbonyl, isonitrile.
  - Nitrogen by degree, N-oxide, nitro/nitrate/nitroso, imine/azo, amide,
    sulfonamide, azide, isonitrile, guanidinium, pyridinium, anions.
  - Oxygen by degree/valence, water, oxonium/oxenium, hydroxyl categories,
    oxide, carboxylate/nitro/phosphate/sulfate/perchlorate terminal oxygen,
    carbonyl/sulfoxide/nitroso.
  - Halogens, silicon, phosphorus, sulfur, chlorine perchlorate, alkali/alkaline
    metals, iron, copper, zinc, bromide/iodide branches.
  - Invalid atom-type branch sets `d_valid = false`.
- Hydrogen typing:
  - Hydrogen attached to carbon/silicon/phosphorus/sulfur.
  - Hydrogen attached to nitrogen by heavy atom type.
  - Hydrogen attached to oxygen by heavy atom type and neighbor context.
  - Invalid hydrogen-type branch sets `d_valid = false`.
- Charge computation:
  - Up-front formal charge assignment for anionic terminal O/S, five-ring
    aromatic anions, conjugated positive nitrogens, diazonium secondary
    nitrogens, simple cations, simple dications/trications, and simple anions.
  - Partial charge equation using `MMFFPBCI` and `MMFFChg`.
  - Negative-neighbor formal charge influence when `fcadj` is zero.
  - Special anionic divalent nitrogen adjustment next to positive neighbors.
  - Signed bond-charge increments and fallback to PBCI difference.
- Parameter selection:
  - Bond type from single `sbmb`/aromatic property match.
  - Angle type from bond-type sum and ring size 3/4.
  - Stretch-bend type switch for angle types 1 through 8.
  - Torsion type empirical correction for single central bond and adjacent
    class-1 bonds.
  - Ring torsion type 4/5 with fused-three-ring guard and hydrogen-containing
    five-membered ring branch.
  - Empirical bond stretch rule, including Bndk present/missing branches and
    Herschbach-Laurie fallback.
  - Empirical angle bend rule, including missing/zero tabulated parameters,
    central-atom coordination/valence/linearity/ring-size branches, and atom
    constant switch.
  - Empirical torsion rule branches `(a)` through `(h)`.
  - OOP missing-parameter branch excludes the term.
  - VdW branch requires both atom parameter records and applies MMFF scaling.

## Convenience API

From `GraphMol/ForceFieldHelpers/MMFF/MMFF.h`:

- `MMFFOptimizeMolecule(mol, maxIters, mmffVariant, nonBondedThresh, confId,
  ignoreInterfragInteractions)`.
  - Branches: construct `MMFFMolProperties`; if invalid, return `(-1, -1)`;
    otherwise construct force field and call `OptimizeMolecule`.
- `MMFFOptimizeMoleculeConfs(mol, res, numThreads, maxIters, mmffVariant,
  nonBondedThresh, ignoreInterfragInteractions)`.
  - Branches: construct `MMFFMolProperties`; if valid, construct shared
    conformer force field with `confId = -1` and call `OptimizeMoleculeConfs`;
    if invalid, resize results to conformer count and fill `(-1, -1)`.
  - Thread branch is delegated to `ForceFieldsHelper::OptimizeMoleculeConfs`.

## Builder Data Types And Helpers

From `GraphMol/ForceFieldHelpers/MMFF/Builder.h` and `.cpp`:

- `constructForceField(mol, nonBondedThresh, confId,
  ignoreInterfragInteractions)`.
- `constructForceField(mol, mmffMolProperties, nonBondedThresh, confId,
  ignoreInterfragInteractions)`.
- `Tools::DefaultTorsionBondSmarts`
  - `string()`.
  - `query()`.
  - Private `create()`.
  - Static data: `ds_string`, `ds_instance`, and `ds_flag` under
    `RDK_BUILD_THREADSAFE_SSS`.
- Neighbor relation enum:
  - `RELATION_1_2`.
  - `RELATION_1_3`.
  - `RELATION_1_4`.
  - `RELATION_1_X`.
- `Tools::twoBitCellPos(nAtoms, i, j)`.
- `Tools::setTwoBitCell(res, pos, value)`.
- `Tools::getTwoBitCell(res, pos)`.
- `Tools::buildNeighborMatrix(mol)`.
- `Tools::addBonds(mol, mmffMolProperties, field)`.
- `Tools::addAngles(mol, mmffMolProperties, field)`.
- `Tools::addStretchBend(mol, mmffMolProperties, field)`.
- `Tools::addOop(mol, mmffMolProperties, field)`.
- `Tools::addTorsions(mol, mmffMolProperties, field, torsionBondSmarts)`.
- `Tools::addVdW(mol, confId, mmffMolProperties, field, neighborMatrix,
  nonBondedThresh, ignoreInterfragInteractions)`.
- `Tools::addEle(mol, confId, mmffMolProperties, field, neighborMatrix,
  nonBondedThresh, ignoreInterfragInteractions)`.
- `Tools::addNonbonded(mol, confId, mmffMolProperties, field, neighborMatrix,
  nonBondedThresh, ignoreInterfragInteractions)`.

## Builder Branch Families To Port

- Precondition branches for force field pointer, MMFF properties pointer, and
  valid atom types.
- Verbosity branches:
  - No output, low output with totals, and high output with per-term tables.
  - Separate stream accumulation in `addNonbonded`.
- `addBonds`:
  - Iterate all bonds.
  - Add term only when `getMMFFBondStretchParams` succeeds.
  - Push contrib only when at least one term exists.
- `buildNeighborMatrix`:
  - Initialize packed relation matrix with `RELATION_1_X`.
  - Set diagonal to `RELATION_1_X`.
  - Branch graph distances 1, 2, 3, connected greater than 3, and disconnected.
- `addAngles`:
  - Skip degree-1 central atoms.
  - Enumerate unique neighbor pairs by iterator ordering.
  - Add term only when `getMMFFAngleBendParams` succeeds.
- `addStretchBend`:
  - Skip degree-1 and linear central atoms.
  - Enumerate unique neighbor pairs by local counters.
  - Add term only when stretch-bend, bond, and angle parameters succeed.
  - Emit one or two high-verbosity rows depending on nonzero force constants.
- `addOop`:
  - Only degree-3 central atoms.
  - Exclude term when OOP parameters are missing.
  - Add three OOP permutations per valid center.
- `DefaultTorsionBondSmarts`:
  - Default SMARTS `[!$(*#*)&!D1]~[!$(*#*)&!D1]`.
  - Thread-safe `std::call_once` branch under `RDK_BUILD_THREADSAFE_SSS`.
  - Non-threadsafe static `created` branch otherwise.
- `addTorsions`:
  - Use default query or parse caller SMARTS and delete non-default query.
  - Require two-atom substructure matches.
  - Central atoms must be SP2 or SP3.
  - Enumerate non-central neighbor pairs excluding the central bond and
    three-membered-ring case where `idx4 == idx1`.
  - Add term only when torsion parameters are nonzero/successful.
- `addVdW`:
  - Optional fragment mapping when `ignoreInterfragInteractions` is true.
  - Skip interfragment pairs when requested.
  - Include only relation `RELATION_1_4` or `RELATION_1_X`.
  - Skip pairs whose conformer distance exceeds `nonBondedThresh`.
  - Add term only when VdW parameters are available.
- `addEle`:
  - Same fragment, relation, and distance gates as VdW.
  - Skip if either partial charge is zero by `isDoubleZero`.
  - Compute `chargeTerm` using dielectric constant and pass dielectric model
    and 1-4 flag.
- `addNonbonded`:
  - Combined VdW/electrostatic term path.
  - Current active branch skips true interfragment pairs by fragment mapping;
    inactive `#else` branch documents an older `cell`-based skip.
  - Add combined term when either VdW or electrostatic branch is present.
  - Respect term toggles `getMMFFVdWTerm` and `getMMFFEleTerm`.
- `constructForceField` overload without properties:
  - Build default `MMFFMolProperties`.
  - Precondition valid atom types.
  - Delegate to properties overload.
- `constructForceField` overload with properties:
  - Copy conformer coordinates by pointer into `ForceField::positions`.
  - Call `initialize`.
  - Add bond, angle, stretch-bend, OOP, torsion, and combined nonbonded terms
    only when their property toggles are enabled.

## MMFF Contrib And Helper Files To Keep In The Port Plan

The Step 675 audit focused on `Params`, `AtomTyper`, `MMFF.h`, and `Builder`,
but the builder depends on the MMFF contrib and utility files below:

- `third_party/rdkit/Code/ForceField/MMFF/AngleBend.h`
- `third_party/rdkit/Code/ForceField/MMFF/AngleBend.cpp`
- `third_party/rdkit/Code/ForceField/MMFF/BondStretch.h`
- `third_party/rdkit/Code/ForceField/MMFF/BondStretch.cpp`
- `third_party/rdkit/Code/ForceField/MMFF/Contribs.h`
- `third_party/rdkit/Code/ForceField/MMFF/Nonbonded.h`
- `third_party/rdkit/Code/ForceField/MMFF/Nonbonded.cpp`
- `third_party/rdkit/Code/ForceField/MMFF/OopBend.h`
- `third_party/rdkit/Code/ForceField/MMFF/OopBend.cpp`
- `third_party/rdkit/Code/ForceField/MMFF/StretchBend.h`
- `third_party/rdkit/Code/ForceField/MMFF/StretchBend.cpp`
- `third_party/rdkit/Code/ForceField/MMFF/TorsionAngle.h`
- `third_party/rdkit/Code/ForceField/MMFF/TorsionAngle.cpp`
- `third_party/rdkit/Code/ForceField/MMFF/AngleConstraint.h`
- `third_party/rdkit/Code/ForceField/MMFF/DistanceConstraint.h`
- `third_party/rdkit/Code/ForceField/MMFF/PositionConstraint.h`
- `third_party/rdkit/Code/ForceField/MMFF/TorsionConstraint.h`

Builder source files that must remain named in consistency checks:

- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/Builder.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/Builder.cpp`

## Global Port Hazards

- Do not implement MMFF atom typing, charge distribution, equivalence-level
  lookup, empirical parameter rules, or nonbonded inclusion with heuristics.
- Preserve `RDKIT_MMFF_PARAMS_USE_STD_MAP` branches as a conscious port choice:
  either reproduce vector-table lookup behavior or document a structured
  complexity/status reason for a different indexed representation.
- Preserve unsupported and missing-parameter behavior explicitly. Missing
  data often returns null or excludes a term instead of inventing parameters.
- Preserve MMFF94 versus MMFF94s data-set selection for OOP and torsion
  collections.
- Preserve preprocessor-controlled branches where behavior changes:
  `RDKIT_MMFF_PARAMS_USE_STD_MAP`, `RDK_BUILD_THREADSAFE_SSS`, and inactive
  empirical-rule blocks where RDKit intentionally excludes formula variants.
- Preserve mutability semantics: RDKit `MMFFMolProperties` may sanitize and set
  molecule properties while COSMolKit public APIs must maintain documented
  value semantics or use explicit operation machinery.
